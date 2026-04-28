"""
=============================================================================
[Project RDL] 다발 예측 실험 #1: FP = 곡률 without 모노드로미
=============================================================================
Conjecture 3 검증:
  - True zero: κ → ∞ AND monodromy = ±π
  - False positive: κ가 높지만 monodromy ≈ 0

방법:
  1. 기존 FP anatomy 파이프라인으로 FP/TP 위치 수집 (5 시드)
  2. 각 위치에서 mpmath로 모노드로미 계산
  3. TP vs FP의 모노드로미 분포 비교
  4. 이중 기준 (κ + monodromy) 적용 시 정밀도 변화 측정

예측: 이중 기준 정밀도 >> 단일 기준(κ만) 정밀도
=============================================================================
"""

import sys, os, time
import numpy as np
import torch
import torch.nn as nn
import mpmath
from scipy import stats

os.environ.setdefault("OMP_NUM_THREADS", "6")
os.environ.setdefault("MKL_NUM_THREADS", "6")
torch.set_num_threads(6)

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))
from gdl.rdl.constants import PrecisionManager
from gdl.rdl.pipeline.xi_feature_dataset import (
    get_or_build_cache, compute_zeros_in_range, XiFeatureDataset
)
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from bundle_utils import xi_func, monodromy_contour

PrecisionManager.setup_precision()

# ─── 설정 ───
T_MIN, T_MAX = 14.0, 50.0
N_POINTS = 1000
HIDDEN = 64
EPOCHS = 100
LR = 1e-3
BATCH = 32
IN_FEATURES = 128
SEEDS = [42, 7, 123, 314, 2024]
MONO_RADII = [0.1, 0.01, 0.001]  # 모노드로미 계산용 반지름

mpmath.mp.dps = 30  # 30자리 정밀도


# ─── mpmath 기반 곡률 계산 ───

def compute_curvature(t):
    """임계선 위 점 t에서 곡률 κ = |ξ'/ξ|² 계산"""
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
    xi_val = xi_func(s)
    if abs(xi_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return float('inf')
    h = mpmath.mpf('1e-15')
    xi_deriv = (xi_func(s + h) - xi_func(s - h)) / (2 * h)
    L = xi_deriv / xi_val
    return float(abs(L)**2)


# ─── F₂ 풍경 계산 & FP/TP 수집 ───

def train_and_get_predictions(dataset, seed):
    """모델 훈련 후 F₂ 풍경에서 FP/TP 위치 반환"""
    torch.manual_seed(seed); np.random.seed(seed)

    val_size = int(len(dataset) * 0.2)
    train_size = len(dataset) - val_size
    train_ds, val_ds = torch.utils.data.random_split(
        dataset, [train_size, val_size],
        generator=torch.Generator().manual_seed(seed)
    )
    train_loader = torch.utils.data.DataLoader(
        train_ds, batch_size=BATCH, shuffle=True, drop_last=True
    )

    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1,
        lambda_tgt=1.0, lambda_pqo=0.5, pqo_mode="cos2",
    )
    opt = torch.optim.Adam(model.parameters(), lr=LR)

    model.train()
    for ep in range(EPOCHS):
        for X, y in train_loader:
            X_in = X.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            opt.zero_grad(set_to_none=True)
            outputs = model(X_in)
            total_loss, _ = loss_fn(**outputs)
            if torch.isnan(total_loss) or torch.isinf(total_loss):
                continue
            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()

    # 전체 데이터에서 F₂ 계산
    model.eval()
    all_X = torch.stack([dataset[i][0] for i in range(len(dataset))])
    all_X = all_X.to(dtype=PrecisionManager.REAL_DTYPE)
    with torch.enable_grad():
        all_X.requires_grad_(True)
        out = model(all_X)
        phi = out['phi'].squeeze()
        psi = out['psi'].squeeze()
        L_G = out['L_G'].squeeze()

    with torch.no_grad():
        phi_real = phi.detach().to(dtype=PrecisionManager.REAL_DTYPE)
        rot = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
        psi_c = psi.detach().to(dtype=PrecisionManager.COMPLEX_DTYPE)
        f2 = (rot * (L_G.detach() - psi_c)).imag.mean(dim=-1)
        f2_abs = torch.abs(f2).cpu().numpy()

    t_vals = dataset.t.numpy()

    # 극소점 찾기
    threshold = np.median(f2_abs) * 0.1
    minima_idx = []
    for i in range(1, len(f2_abs) - 1):
        if f2_abs[i] < f2_abs[i-1] and f2_abs[i] < f2_abs[i+1]:
            if f2_abs[i] < threshold:
                minima_idx.append(i)

    pred_t = t_vals[minima_idx]

    # 진짜 영점
    true_zeros = compute_zeros_in_range(T_MIN, T_MAX)
    true_t = np.array([float(z) for z in true_zeros])

    # TP/FP 분류
    tp_t, fp_t = [], []
    tol = (T_MAX - T_MIN) / N_POINTS * 4  # 4 grid points tolerance
    for pt in pred_t:
        dists = np.abs(true_t - pt)
        if dists.min() < tol:
            tp_t.append(pt)
        else:
            fp_t.append(pt)

    return np.array(tp_t), np.array(fp_t), true_t


# ─── 메인 실험 ───

def main():
    out_path = os.path.expanduser(
        '~/Desktop/gdl_unified/results/bundle_prediction_fp_monodromy.txt'
    )

    print("=" * 70)
    print("다발 예측 실험 #1: FP = 곡률 without 모노드로미")
    print("=" * 70)

    # 데이터셋 준비
    cache_path = os.path.expanduser(
        f'~/Desktop/gdl_unified/outputs/cache/xi_cache_{T_MIN}_{T_MAX}_{N_POINTS}.npz'
    )
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    cache = get_or_build_cache(cache_path, T_MIN, T_MAX, N_POINTS)
    dataset = XiFeatureDataset(cache, in_features=IN_FEATURES)

    all_tp = []
    all_fp = []

    # true_t는 해석적 영점 — 시드 무관, 1회만 계산
    true_t = None

    for seed in SEEDS:
        print(f"\n--- Seed {seed} ---")
        tp_t, fp_t, tz = train_and_get_predictions(dataset, seed)
        if true_t is None:
            true_t = tz
        print(f"  TP: {len(tp_t)}, FP: {len(fp_t)}, True zeros: {len(tz)}")
        all_tp.extend(tp_t.tolist())
        all_fp.extend(fp_t.tolist())

    # 중복 제거 (±0.5 이내)
    def deduplicate(arr, tol=0.5):
        if len(arr) == 0:
            return arr
        arr = sorted(arr)
        result = [arr[0]]
        for x in arr[1:]:
            if x - result[-1] > tol:
                result.append(x)
        return np.array(result)

    tp_unique = deduplicate(all_tp)
    fp_unique = deduplicate(all_fp)

    print(f"\n합산 (중복 제거): TP {len(tp_unique)}, FP {len(fp_unique)}")

    # 모노드로미 & 곡률 계산
    print("\n모노드로미 & 곡률 계산 중...")
    print(f"  (TP는 nearest true_zero에서 계산 — 사이클 27 코드 오류 수정)")

    results = []

    print("\n  [True Positives] — monodromy at nearest true_zero")
    for i, t_det in enumerate(tp_unique[:20]):  # 최대 20개
        # 핵심 수정: 검출점이 아닌 nearest true_zero에서 monodromy 계산
        dists = np.abs(true_t - t_det)
        nearest_idx = np.argmin(dists)
        t_zero = float(true_t[nearest_idx])
        dist = dists[nearest_idx]
        # 다중 반지름으로 모노드로미 계산
        monos_by_r = {}
        for r in MONO_RADII:
            monos_by_r[r] = monodromy_contour(t_zero, radius=r)
        mono = monos_by_r[0.001]  # 가장 작은 반지름 사용
        kappa = compute_curvature(t_det)  # 곡률은 검출점에서 (신경망이 보는 값)
        results.append(('TP', t_det, mono, kappa, t_zero, dist, monos_by_r))
        print(f"    t_det={t_det:.4f} → t_zero={t_zero:.4f} (Δ={dist:.4f}): "
              f"mono(r=0.1)={monos_by_r[0.1]/np.pi:.4f}π, "
              f"mono(r=0.01)={monos_by_r[0.01]/np.pi:.4f}π, "
              f"mono(r=0.001)={monos_by_r[0.001]/np.pi:.4f}π, κ={kappa:.2e}")

    print("\n  [False Positives] — monodromy at detection point")
    for i, t in enumerate(fp_unique[:30]):  # 최대 30개
        monos_by_r = {}
        for r in MONO_RADII:
            monos_by_r[r] = monodromy_contour(t, radius=r)
        mono = monos_by_r[0.001]
        kappa = compute_curvature(t)
        results.append(('FP', t, mono, kappa, None, None, monos_by_r))
        print(f"    t={t:.4f}: "
              f"mono(r=0.1)={monos_by_r[0.1]/np.pi:.4f}π, "
              f"mono(r=0.01)={monos_by_r[0.01]/np.pi:.4f}π, "
              f"mono(r=0.001)={monos_by_r[0.001]/np.pi:.4f}π, κ={kappa:.2e}")

    # 통계
    tp_monos = [abs(r[2]) for r in results if r[0] == 'TP']
    fp_monos = [abs(r[2]) for r in results if r[0] == 'FP']
    tp_kappas = [r[3] for r in results if r[0] == 'TP']
    fp_kappas = [r[3] for r in results if r[0] == 'FP']

    print("\n" + "=" * 70)
    print("통계 요약")
    print("=" * 70)

    # 반지름별 통계
    for r in MONO_RADII:
        tp_m = [abs(res[6][r]) for res in results if res[0] == 'TP']
        fp_m = [abs(res[6][r]) for res in results if res[0] == 'FP']
        if tp_m and fp_m:
            print(f"\n  radius={r}:")
            print(f"    TP |mono|/π: mean={np.mean(tp_m)/np.pi:.4f}, std={np.std(tp_m)/np.pi:.4f}")
            print(f"    FP |mono|/π: mean={np.mean(fp_m)/np.pi:.4f}, std={np.std(fp_m)/np.pi:.4f}")
            ks_stat, ks_p = stats.ks_2samp(tp_m, fp_m)
            print(f"    KS test: stat={ks_stat:.4f}, p={ks_p:.6f}")

    if tp_monos:
        print(f"\n주 통계 (r=0.001):")
        print(f"  TP |mono|/π: mean={np.mean(tp_monos)/np.pi:.4f}, "
              f"std={np.std(tp_monos)/np.pi:.4f}")
    if fp_monos:
        print(f"  FP |mono|/π: mean={np.mean(fp_monos)/np.pi:.4f}, "
              f"std={np.std(fp_monos)/np.pi:.4f}")

    if tp_kappas:
        print(f"\nTP 곡률 κ: mean={np.mean(tp_kappas):.2e}, "
              f"median={np.median(tp_kappas):.2e}")
    if fp_kappas:
        print(f"FP 곡률 κ: mean={np.mean(fp_kappas):.2e}, "
              f"median={np.median(fp_kappas):.2e}")

    # KS 검정 (r=0.001 주 결과)
    ks_stat_main, ks_p_main = None, None
    if len(tp_monos) >= 2 and len(fp_monos) >= 2:
        ks_stat_main, ks_p_main = stats.ks_2samp(tp_monos, fp_monos)
        print(f"\nKS 검정 (TP vs FP, r=0.001): stat={ks_stat_main:.4f}, p={ks_p_main:.6f}")

    # 성공 기준 체크
    print("\n" + "=" * 70)
    print("성공 기준 체크")
    print("=" * 70)

    # 기준 1: FP |mono| < 0.3 비율 > 70%
    fp_low_mono = sum(1 for m in fp_monos if m < 0.3) / len(fp_monos) if fp_monos else 0
    crit1 = fp_low_mono > 0.7
    print(f"  FP |mono| < 0.3 비율: {fp_low_mono:.1%} ({'✅' if crit1 else '❌'} > 70%)")

    # 기준 2: TP |mono| > 2.0 비율 > 70%
    tp_high_mono = sum(1 for m in tp_monos if m > 2.0) / len(tp_monos) if tp_monos else 0
    crit2 = tp_high_mono > 0.7
    print(f"  TP |mono| > 2.0 비율: {tp_high_mono:.1%} ({'✅' if crit2 else '❌'} > 70%)")

    # 이중 기준 테스트
    mono_threshold = np.pi * 0.5  # |Δarg| > π/2 = 모노드로미 있음

    tp_pass = sum(1 for m in tp_monos if m > mono_threshold)
    fp_pass = sum(1 for m in fp_monos if m > mono_threshold)

    print(f"\n이중 기준 (κ + |mono| > π/2):")
    print(f"  TP 통과: {tp_pass}/{len(tp_monos)}")
    print(f"  FP 통과: {fp_pass}/{len(fp_monos)}")

    new_precision, old_precision = 0, 0
    if tp_pass + fp_pass > 0:
        new_precision = tp_pass / (tp_pass + fp_pass)
        old_precision = len(tp_monos) / (len(tp_monos) + len(fp_monos))
        print(f"\n  기존 정밀도 (κ만): {old_precision:.1%}")
        print(f"  새 정밀도 (κ+mono): {new_precision:.1%}")
        print(f"  개선: {new_precision/old_precision:.1f}×")

    # 기준 3: 이중 기준 정밀도 > 40%
    crit3 = new_precision > 0.4
    print(f"  이중기준 정밀도 > 40%: {'✅' if crit3 else '❌'} ({new_precision:.1%})")

    # 기준 4: KS p < 0.01
    crit4 = (ks_p_main is not None and ks_p_main < 0.01)
    print(f"  KS p < 0.01: {'✅' if crit4 else '❌'} (p={ks_p_main})")

    overall = "양성" if (crit1 and crit2 and crit3 and crit4) else \
              "약한 양성" if sum([crit1, crit2, crit3, crit4]) >= 2 else "음성"
    print(f"\n  *** 종합 판정: {overall} ({sum([crit1,crit2,crit3,crit4])}/4 기준 충족) ***")

    # 결과 저장
    with open(out_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("다발 예측 실험 #1: FP = 곡률 without 모노드로미\n")
        f.write(f"Conjecture 3 검증 — 사이클 #402\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"설정: t∈[{T_MIN},{T_MAX}], K={IN_FEATURES}, "
                f"시드 {SEEDS}, {EPOCHS} epochs\n")
        f.write(f"모노드로미 반지름: {MONO_RADII}\n")
        f.write(f"수정: TP monodromy를 nearest true_zero에서 계산, "
                f"bundle_utils.monodromy_contour 사용\n\n")

        # 다중 반지름 테이블
        header = f"{'Type':>4} {'t_det':>10} {'t_zero':>10} {'Δt':>7}"
        for r in MONO_RADII:
            header += f" {'mono_r='+str(r):>14}"
        header += f" {'κ':>12}"
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")
        for res in results:
            if res[0] == 'TP':
                line = f"{res[0]:>4} {res[1]:>10.4f} {res[4]:>10.4f} {res[5]:>7.4f}"
            else:
                line = f"{res[0]:>4} {res[1]:>10.4f} {'—':>10} {'—':>7}"
            for r in MONO_RADII:
                line += f" {res[6][r]/np.pi:>14.4f}π"
            line += f" {res[3]:>12.2e}"
            f.write(line + "\n")

        f.write(f"\n\n반지름별 통계:\n")
        for r in MONO_RADII:
            tp_m = [abs(res[6][r]) for res in results if res[0] == 'TP']
            fp_m = [abs(res[6][r]) for res in results if res[0] == 'FP']
            if tp_m and fp_m:
                ks_s, ks_p = stats.ks_2samp(tp_m, fp_m)
                f.write(f"\n  radius={r}:\n")
                f.write(f"    TP |mono|/π: {np.mean(tp_m)/np.pi:.4f} ± {np.std(tp_m)/np.pi:.4f}\n")
                f.write(f"    FP |mono|/π: {np.mean(fp_m)/np.pi:.4f} ± {np.std(fp_m)/np.pi:.4f}\n")
                f.write(f"    KS: stat={ks_s:.4f}, p={ks_p:.6f}\n")

        f.write(f"\n\n주 통계 (r=0.001):\n")
        if tp_monos:
            f.write(f"  TP |mono|/π: {np.mean(tp_monos)/np.pi:.4f} ± {np.std(tp_monos)/np.pi:.4f}\n")
        if fp_monos:
            f.write(f"  FP |mono|/π: {np.mean(fp_monos)/np.pi:.4f} ± {np.std(fp_monos)/np.pi:.4f}\n")

        f.write(f"\n성공 기준:\n")
        f.write(f"  FP |mono|<0.3 비율: {fp_low_mono:.1%} ({'PASS' if crit1 else 'FAIL'})\n")
        f.write(f"  TP |mono|>2.0 비율: {tp_high_mono:.1%} ({'PASS' if crit2 else 'FAIL'})\n")
        f.write(f"  이중기준 정밀도: {new_precision:.1%} ({'PASS' if crit3 else 'FAIL'})\n")
        f.write(f"  KS p-value: {ks_p_main} ({'PASS' if crit4 else 'FAIL'})\n")
        f.write(f"\n  이중 기준 정밀도: {new_precision:.1%} (기존 {old_precision:.1%})\n")
        f.write(f"  종합 판정: {overall} ({sum([crit1,crit2,crit3,crit4])}/4)\n")

    print(f"\n결과 저장: {out_path}")
    print("완료.")


if __name__ == '__main__':
    main()
