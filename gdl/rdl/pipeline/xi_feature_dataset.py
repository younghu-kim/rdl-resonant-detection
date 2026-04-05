"""
=============================================================================
[Project RDL] Xi Feature Dataset — 푸리에 기저 특징 확장 래퍼
=============================================================================
스칼라 t → 다차원 푸리에 특징 벡터로 확장하여 MasterResonantNetwork에 투입.
mpmath 계산 결과를 .pt 캐시로 저장하여 재계산 방지.
"""

import os
import math
import torch
from torch.utils.data import Dataset, DataLoader
from typing import Tuple

from gdl.rdl.constants import PrecisionManager

try:
    import mpmath
    MPMATH_AVAILABLE = True
except ImportError:
    MPMATH_AVAILABLE = False


# LMFDB 알려진 비자명 영점 (처음 10개) — 기본 범위 t∈[10,50]
KNOWN_ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832
]


def compute_zeros_in_range(t_min, t_max, dps=20):
    """mpmath로 [t_min, t_max] 범위의 제타 영점을 계산하여 반환"""
    if not MPMATH_AVAILABLE:
        raise ImportError("mpmath 필요: pip install mpmath")
    mpmath.mp.dps = dps
    import math as _math
    # N(T) ~ (T/2pi)*log(T/2pi) - T/2pi 으로 시작 인덱스 추정
    if t_min > 20:
        k_start = max(1, int((t_min / (2 * _math.pi)) * _math.log(t_min / (2 * _math.pi))
                              - t_min / (2 * _math.pi)) - 2)
    else:
        k_start = 1
    zeros = []
    k = k_start
    while True:
        z = mpmath.zetazero(k)
        t_val = float(z.imag)
        if t_val > t_max:
            break
        if t_val >= t_min:
            zeros.append(round(t_val, 6))
        k += 1
    print(f"[Zeros] {len(zeros)}개 영점 발견 in [{t_min}, {t_max}] (k={k_start}~{k-1})")
    return zeros


def compute_xi(s):
    """완성 xi-함수: xi(s) = (1/2)s(s-1)pi^{-s/2}Gamma(s/2)zeta(s)"""
    half = mpmath.mpf('0.5')
    return half * s * (s - 1) * mpmath.pi ** (-s / 2) * mpmath.gamma(s / 2) * mpmath.zeta(s)


def build_xi_cache(t_min, t_max, num_points, mp_dps=50, include_zeros=True, zeros_list=None):
    """mpmath로 xi(1/2+it) 계산 후 텐서 캐시 반환"""
    if not MPMATH_AVAILABLE:
        raise ImportError("mpmath 필요: pip install mpmath")

    mpmath.mp.dps = mp_dps
    dtype = PrecisionManager.REAL_DTYPE

    t_values = torch.linspace(t_min, t_max, num_points, dtype=dtype).tolist()
    if include_zeros:
        use_zeros = zeros_list if zeros_list is not None else KNOWN_ZEROS
        for z in use_zeros:
            if t_min <= z <= t_max:
                t_values.append(z)
        t_values = sorted(set(t_values))

    N = len(t_values)
    t_tensor = torch.tensor(t_values, dtype=dtype)
    xi_real = torch.empty(N, dtype=dtype)
    xi_imag = torch.empty(N, dtype=dtype)

    print(f"[Xi Cache] {N}개 포인트 계산 시작 (dps={mp_dps})...")
    for i, t_val in enumerate(t_values):
        s = mpmath.mpf('0.5') + mpmath.mpf(str(t_val)) * mpmath.j
        xi_val = compute_xi(s)
        xi_real[i] = float(xi_val.real)
        xi_imag[i] = float(xi_val.imag)
        if (i + 1) % 100 == 0 or i == N - 1:
            print(f"  [{i+1}/{N}] t={t_val:.4f}, |xi|={math.sqrt(float(xi_val.real)**2 + float(xi_val.imag)**2):.6f}")

    return {'t': t_tensor, 'xi_real': xi_real, 'xi_imag': xi_imag}


def get_or_build_cache(cache_path, t_min, t_max, num_points, mp_dps=50, zeros_list=None):
    """캐시 파일이 있으면 로드, 없으면 계산 후 저장"""
    if os.path.exists(cache_path):
        print(f"[Xi Cache] 캐시 로드: {cache_path}")
        return torch.load(cache_path, weights_only=True)

    data = build_xi_cache(t_min, t_max, num_points, mp_dps, zeros_list=zeros_list)
    os.makedirs(os.path.dirname(cache_path) or '.', exist_ok=True)
    torch.save(data, cache_path)
    print(f"[Xi Cache] 캐시 저장: {cache_path}")
    return data


class XiFeatureDataset(Dataset):
    """
    스칼라 t를 푸리에 기저로 확장하여 in_features 차원 벡터를 반환.

    특징 벡터: [t_norm, sin(w1*t), cos(w1*t), sin(w2*t), cos(w2*t), ...]
    타겟: [Re(xi), Im(xi)]
    """
    def __init__(self, cache_data: dict, in_features: int = 64):
        super().__init__()
        self.t = cache_data['t']
        self.xi_real = cache_data['xi_real']
        self.xi_imag = cache_data['xi_imag']
        self.in_features = in_features

        # 푸리에 주파수 설계: 영점 간격 (~5-10) 에 맞춘 다중 스케일
        # 구조: [t_norm, sin(w1*t), cos(w1*t), sin(w2*t), cos(w2*t), ...]
        # in_features가 홀수면 마지막은 sin만
        self.num_freqs = (in_features - 1) // 2
        self.freqs = torch.arange(1, self.num_freqs + 1, dtype=PrecisionManager.REAL_DTYPE)
        t_range = self.t.max() - self.t.min() + 1e-8
        self.omega = 2.0 * math.pi * self.freqs / t_range

        # 정규화 통계
        self.t_mean = self.t.mean()
        self.t_std = self.t.std() + 1e-8

        # xi 진폭/위상 사전 계산
        self.xi_amp = torch.sqrt(self.xi_real**2 + self.xi_imag**2)
        self.xi_phase = torch.atan2(self.xi_imag, self.xi_real)

        # 영점 마스크
        self.is_near_zero = self.xi_amp < (self.xi_amp.median() * 0.01)

    def __len__(self):
        return len(self.t)

    def _build_features(self, t_val):
        """t 스칼라 → in_features 차원 푸리에 기저 벡터"""
        t_norm = (t_val - self.t_mean) / self.t_std
        phases = self.omega * t_val  # [num_freqs]
        sin_vals = torch.sin(phases)
        cos_vals = torch.cos(phases)

        features = torch.zeros(self.in_features, dtype=PrecisionManager.REAL_DTYPE)
        features[0] = t_norm
        # sin/cos 교대 배치
        for k in range(self.num_freqs):
            idx_sin = 1 + 2 * k
            idx_cos = 2 + 2 * k
            if idx_sin < self.in_features:
                features[idx_sin] = sin_vals[k]
            if idx_cos < self.in_features:
                features[idx_cos] = cos_vals[k]
        return features

    def __getitem__(self, idx) -> Tuple[torch.Tensor, torch.Tensor]:
        features = self._build_features(self.t[idx])

        xi_target = torch.stack([self.xi_real[idx], self.xi_imag[idx]])
        return features, xi_target

    def get_known_zero_indices(self):
        """알려진 영점에 해당하는 인덱스 반환"""
        indices = []
        for z in KNOWN_ZEROS:
            diffs = (self.t - z).abs()
            min_idx = diffs.argmin().item()
            if diffs[min_idx] < 1e-3:
                indices.append(min_idx)
        return indices

    def get_features_at_t(self, t_values):
        """임의의 t 값들에 대해 특징 벡터 배치 생성"""
        t_tensor = torch.tensor(t_values, dtype=PrecisionManager.REAL_DTYPE)
        features = torch.stack([self._build_features(t) for t in t_tensor])
        return features


def get_xi_feature_dataloaders(in_features=64, batch_size=32,
                                t_min=10.0, t_max=50.0, num_points=500,
                                val_ratio=0.2, seed=42, cache_dir='outputs',
                                zeros_list=None):
    """Xi 특징 데이터셋의 훈련/검증 DataLoader 팩토리"""
    cache_path = os.path.join(cache_dir, f'xi_cache_t{t_min}-{t_max}_n{num_points}.pt')
    cache_data = get_or_build_cache(cache_path, t_min, t_max, num_points, zeros_list=zeros_list)

    dataset = XiFeatureDataset(cache_data, in_features=in_features)

    val_size = int(len(dataset) * val_ratio)
    train_size = len(dataset) - val_size

    train_ds, val_ds = torch.utils.data.random_split(
        dataset, [train_size, val_size], generator=torch.Generator().manual_seed(seed)
    )

    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True, drop_last=True)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False, drop_last=False)

    return train_loader, val_loader, dataset


if __name__ == "__main__":
    PrecisionManager.setup_precision()
    print("\n--- [RDL] Xi Feature Dataset Test ---")

    train_loader, val_loader, dataset = get_xi_feature_dataloaders(
        in_features=64, num_points=100, t_min=12.0, t_max=16.0,
        cache_dir='/home/k0who029/Desktop/XI_function/outputs'
    )

    X, y = next(iter(train_loader))
    print(f"배치: X={X.shape}, y={y.shape}, dtype={X.dtype}")
    print(f"데이터셋 크기: {len(dataset)}")

    zeros = dataset.get_known_zero_indices()
    print(f"영점 인덱스: {zeros}")
    for idx in zeros:
        print(f"  t={dataset.t[idx]:.6f}, |xi|={dataset.xi_amp[idx]:.8f}")

    print("[OK] Xi Feature Dataset 테스트 통과")
