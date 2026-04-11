"""
=============================================================================
[Project RDL] 결정론적 vs 랜덤 푸리에 특성 비교 실험
=============================================================================
Tancik et al. (NeurIPS 2020) 방식의 랜덤 푸리에 특성 B ~ N(0, σ²)과
현재 RDL의 결정론적 등간격 ω_k = 2πk/W 를 동일 조건에서 비교.

측정: 영점 MSE, 비영점 MSE, F₂ 기반 검출, 학습 곡선
"""

import sys, os, time, math
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
from gdl.rdl.constants import PrecisionManager
from gdl.rdl.pipeline.xi_feature_dataset import (
    get_or_build_cache, compute_zeros_in_range, XiFeatureDataset
)
from gdl.rdl.layers.resonant_block import ResonantBlock
from gdl.rdl.layers.embedding import ComplexLinearEmbedding

PrecisionManager.setup_precision()

# ─── 랜덤 푸리에 데이터셋 ───
class RandomFourierDataset(Dataset):
    """Tancik 방식: B ~ N(0, σ²), features = [sin(2πBt), cos(2πBt)]"""
    def __init__(self, cache_data, in_features=128, sigma=10.0, seed=0):
        self.t = cache_data['t']
        self.xi_real = cache_data['xi_real']
        self.xi_imag = cache_data['xi_imag']
        self.in_features = in_features

        # 랜덤 B 행렬 생성 (고정 시드)
        rng = np.random.RandomState(seed)
        num_freqs = in_features // 2
        self.B = torch.tensor(
            rng.randn(num_freqs) * sigma,
            dtype=PrecisionManager.REAL_DTYPE
        )

        self.t_mean = self.t.mean()
        self.t_std = self.t.std() + 1e-8
        self.xi_amp = torch.sqrt(self.xi_real**2 + self.xi_imag**2)
        self.is_near_zero = self.xi_amp < (self.xi_amp.median() * 0.01)

    def __len__(self):
        return len(self.t)

    def __getitem__(self, idx):
        t_val = self.t[idx]
        t_norm = (t_val - self.t_mean) / self.t_std

        # Tancik: γ(x) = [sin(2πBx), cos(2πBx)]
        phases = 2.0 * math.pi * self.B * t_norm
        sin_vals = torch.sin(phases)
        cos_vals = torch.cos(phases)
        features = torch.cat([sin_vals, cos_vals])  # [in_features]

        xi_target = torch.stack([self.xi_real[idx], self.xi_imag[idx]])
        return features, xi_target


# ─── 간단한 비교 모델 (공정한 비교를 위해 동일 아키텍처) ───
class SimpleXiRegressor(nn.Module):
    """임베딩 후 3-layer MLP로 xi 예측"""
    def __init__(self, in_features, hidden=64):
        super().__init__()
        self.embed = ComplexLinearEmbedding(in_features, hidden)
        self.blocks = nn.Sequential(
            ResonantBlock(hidden, hidden),
            ResonantBlock(hidden, hidden),
        )
        self.head = nn.Linear(hidden * 2, 2, dtype=PrecisionManager.REAL_DTYPE)

    def forward(self, x):
        z = self.embed(x)
        for blk in self.blocks:
            z = blk(z)
        # 복소 → 실수 (real, imag concat)
        out = torch.cat([z.real, z.imag], dim=-1)
        return self.head(out)


def train_one(dataset, seed, epochs=100, hidden=64, batch_size=32, lr=1e-3):
    """한 설정에 대해 훈련 + 평가"""
    torch.manual_seed(seed)
    np.random.seed(seed)

    val_size = int(len(dataset) * 0.2)
    train_size = len(dataset) - val_size
    train_ds, val_ds = torch.utils.data.random_split(
        dataset, [train_size, val_size],
        generator=torch.Generator().manual_seed(seed)
    )
    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True, drop_last=True)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False)

    model = SimpleXiRegressor(dataset.in_features, hidden)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss()

    best_val = float('inf')
    for ep in range(epochs):
        model.train()
        for X, y in train_loader:
            pred = model(X)
            loss = criterion(pred, y)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        if (ep + 1) % 20 == 0 or ep == epochs - 1:
            model.eval()
            val_loss = 0; n = 0
            with torch.no_grad():
                for X, y in val_loader:
                    pred = model(X)
                    val_loss += criterion(pred, y).item() * X.size(0)
                    n += X.size(0)
            val_loss /= max(n, 1)
            best_val = min(best_val, val_loss)

    # 영점/비영점 분석
    model.eval()
    all_mse_zero = []
    all_mse_nonzero = []
    with torch.no_grad():
        for i in range(len(dataset)):
            X, y = dataset[i]
            pred = model(X.unsqueeze(0))
            mse = ((pred.squeeze() - y) ** 2).mean().item()
            if dataset.is_near_zero[i]:
                all_mse_zero.append(mse)
            else:
                all_mse_nonzero.append(mse)

    mse_zero = np.mean(all_mse_zero) if all_mse_zero else 0
    mse_nonzero = np.mean(all_mse_nonzero) if all_mse_nonzero else 1
    ratio = mse_zero / (mse_nonzero + 1e-12)

    return best_val, mse_zero, mse_nonzero, ratio


def main():
    t_min, t_max = 100.0, 200.0
    num_points = 1000
    in_features = 128
    seeds = [42, 7, 123]
    epochs = 100
    hidden = 64
    sigmas = [1.0, 10.0, 100.0]

    out_lines = []
    def log(msg):
        print(msg)
        out_lines.append(msg)

    log("=" * 72)
    log("  결정론적 vs 랜덤 푸리에 특성 비교 실험")
    log("=" * 72)
    log(f"  t∈[{t_min}, {t_max}], {num_points}점, in_features={in_features}")
    log(f"  seeds={seeds}, epochs={epochs}, hidden={hidden}")
    log(f"  랜덤 σ: {sigmas}")
    log("")

    # 캐시 준비
    zeros_list = compute_zeros_in_range(t_min, t_max)
    cache_path = os.path.expanduser(
        f'~/Desktop/gdl_unified/outputs/xi_cache_t{t_min}-{t_max}_n{num_points}.pt'
    )
    cache_data = get_or_build_cache(cache_path, t_min, t_max, num_points, zeros_list=zeros_list)

    start = time.time()
    results = {}

    # ── 1) 결정론적 (현재 RDL) ──
    log("─" * 72)
    log("  [1] 결정론적 푸리에 (ω_k = 2πk/W)")
    log("─" * 72)
    det_vals = []
    for s in seeds:
        ds = XiFeatureDataset(cache_data, in_features=in_features)
        val, mz, mnz, ratio = train_one(ds, s, epochs, hidden)
        log(f"    seed={s}: val={val:.6f}, zero_MSE={mz:.6f}, nonzero_MSE={mnz:.6f}, ratio={ratio:.4f}")
        det_vals.append((val, mz, mnz, ratio))

    det_ratio_mean = np.mean([v[3] for v in det_vals])
    det_ratio_std = np.std([v[3] for v in det_vals])
    det_val_mean = np.mean([v[0] for v in det_vals])
    log(f"  → 결정론적 앙상블: ratio={det_ratio_mean:.4f}±{det_ratio_std:.4f}, val={det_val_mean:.6f}")
    results['deterministic'] = det_vals

    # ── 2) 랜덤 푸리에 (σ sweep) ──
    for sigma in sigmas:
        log("")
        log("─" * 72)
        log(f"  [2] 랜덤 푸리에 (σ={sigma})")
        log("─" * 72)
        rand_vals = []
        for s in seeds:
            ds = RandomFourierDataset(cache_data, in_features=in_features, sigma=sigma, seed=s)
            val, mz, mnz, ratio = train_one(ds, s, epochs, hidden)
            log(f"    seed={s}: val={val:.6f}, zero_MSE={mz:.6f}, nonzero_MSE={mnz:.6f}, ratio={ratio:.4f}")
            rand_vals.append((val, mz, mnz, ratio))

        r_mean = np.mean([v[3] for v in rand_vals])
        r_std = np.std([v[3] for v in rand_vals])
        r_val = np.mean([v[0] for v in rand_vals])
        log(f"  → σ={sigma} 앙상블: ratio={r_mean:.4f}±{r_std:.4f}, val={r_val:.6f}")
        results[f'random_sigma{sigma}'] = rand_vals

    elapsed = time.time() - start

    # ── 요약 ──
    log("")
    log("=" * 72)
    log("  요약")
    log("=" * 72)
    log(f"  {'모드':<25} {'val_MSE':>10} {'zero/nonzero ratio':>20}")
    log(f"  {'-'*25} {'-'*10} {'-'*20}")

    for name, vals in results.items():
        v_mean = np.mean([v[0] for v in vals])
        r_mean = np.mean([v[3] for v in vals])
        r_std = np.std([v[3] for v in vals])
        log(f"  {name:<25} {v_mean:>10.6f} {r_mean:>10.4f}±{r_std:.4f}")

    log(f"\n  총 실행 시간: {elapsed:.1f}s")

    # 판정
    best_random = min(
        [(k, np.mean([v[0] for v in vs])) for k, vs in results.items() if k != 'deterministic'],
        key=lambda x: x[1]
    )
    det_mean = np.mean([v[0] for v in results['deterministic']])
    if best_random[1] < det_mean * 0.95:
        verdict = f"양성: 랜덤 ({best_random[0]})이 결정론적보다 {(1-best_random[1]/det_mean)*100:.1f}% 우세"
    elif best_random[1] > det_mean * 1.05:
        verdict = "음성: 결정론적 푸리에가 랜덤보다 우세"
    else:
        verdict = "중립: 유의미한 차이 없음 (5% 이내)"

    log(f"\n  판정: {verdict}")

    # 저장
    os.makedirs(os.path.expanduser('~/Desktop/gdl_unified/results'), exist_ok=True)
    out_path = os.path.expanduser('~/Desktop/gdl_unified/results/compare_fourier_modes.txt')
    with open(out_path, 'w') as f:
        f.write('\n'.join(out_lines))
    log(f"\n  결과 저장: {out_path}")


if __name__ == '__main__':
    main()
