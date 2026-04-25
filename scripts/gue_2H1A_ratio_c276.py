#!/usr/bin/env python3
"""
C-276: GUE 랜덤 행렬 2H₁/A 비율 계산
GUE(N) 고유값에서 A, H₁, B 계산 → 2H₁/A ≈ 2/3 보편 예측 검증
+ ρ(A, gap_min) A-gap 현상 GUE 검증

결과: results/gue_2H1A_ratio_c276.txt
"""

import numpy as np
from scipy import stats
import os
import sys
import time

def generate_gue_eigvals(N, seed=None):
    """GUE(N) 행렬 생성 후 실수 고유값 반환 (정렬됨).
    H = (A + A†)/(2√N), A_{ij} ~ CN(0,1) → 고유값 ∈ [-2,2] 반원 분포.
    """
    rng = np.random.default_rng(seed)
    # 실수부/허수부 각각 N(0,1)
    real = rng.standard_normal((N, N))
    imag = rng.standard_normal((N, N))
    A = (real + 1j * imag) / np.sqrt(2.0)
    H = (A + A.conj().T) / (2.0 * np.sqrt(N))
    eigvals = np.linalg.eigvalsh(H)  # 자동 정렬
    return eigvals.real


def semicircle_cdf(x):
    """반원 분포 CDF: ρ(x)=(1/(2π))√(4-x²), x∈[-2,2].
    반환: F(x) = ∫_{-2}^{x} ρ(t) dt ∈ [0,1]
    """
    x = np.clip(x, -2.0 + 1e-12, 2.0 - 1e-12)
    return (np.arcsin(x / 2.0) + (x / 2.0) * np.sqrt(np.maximum(0.0, 1.0 - x**2 / 4.0))) / np.pi + 0.5


def unfold(lambdas):
    """GUE 반원법으로 고유값 펼치기. 평균 간격 ≈ 1."""
    return len(lambdas) * semicircle_cdf(lambdas)


def compute_metrics_single(lambdas):
    """
    단일 GUE 실현에서 A, H₁, S₁, 2H₁/A 계산 (벡터화).

    Returns dict:
      ratio_bare  : 2H₁/A_bare  (B=S₁, smooth 미보정)
      ratio_corr  : 2H₁/A_corr  (B=S₁−N*λ/2, smooth 보정)
      h1nn_ratio  : H₁^{NN}/H₁
      A_bare      : 각 고유값별 A_bare
      A_corr      : 각 고유값별 A_corr
      gap_min     : unfolded 최소 간격 (좌우 중 작은 쪽)
      gap_right   : unfolded 우측 간격
    (모두 내부 80% 고유값 기준)
    """
    N = len(lambdas)
    xi = unfold(lambdas)          # 펼친 고유값

    # ─── N×N 역차이 행렬 ───
    # diff[i,j] = λ_i - λ_j,  대각=NaN
    diff = lambdas[:, None] - lambdas[None, :]   # (N,N)
    np.fill_diagonal(diff, np.nan)

    inv1 = 1.0 / diff                  # 1/(λ_i-λ_j)
    inv2 = 1.0 / diff**2               # 1/(λ_i-λ_j)²

    S1  = np.nansum(inv1, axis=1)      # (N,)  S₁(j) = Σ_{k≠j} 1/(λ_j-λ_k)
    H1  = np.nansum(inv2, axis=1)      # (N,)  H₁(j)

    # ─── B_smooth = N·λ/2  (반원 힐버트 변환) ───
    B_smooth = N * lambdas / 2.0
    B_corr   = S1 - B_smooth

    A_bare = S1**2    + 2.0 * H1
    A_corr = B_corr**2 + 2.0 * H1

    ratio_bare = np.where(A_bare > 0, 2.0 * H1 / A_bare, np.nan)
    ratio_corr = np.where(A_corr > 0, 2.0 * H1 / A_corr, np.nan)

    # ─── H₁^{NN} (nearest neighbor only) ───
    gaps = np.diff(lambdas)            # (N-1,), 양수
    inv_g2 = 1.0 / gaps**2
    H1_nn = np.zeros(N)
    H1_nn[1:]  += inv_g2              # 왼쪽 NN 기여
    H1_nn[:-1] += inv_g2              # 오른쪽 NN 기여
    h1nn_ratio = np.where(H1 > 0, H1_nn / H1, np.nan)

    # ─── unfolded gap ───
    xi_gaps  = np.diff(xi)            # (N-1,)
    gap_left  = np.concatenate([[np.nan], xi_gaps])
    gap_right = np.concatenate([xi_gaps, [np.nan]])
    gap_min   = np.minimum(gap_left, gap_right)

    # ─── 내부 80% 선택 ───
    i0 = max(1, int(N * 0.10))        # 양 끝 10% 제외 + 경계 안전
    i1 = min(N - 1, int(N * 0.90))
    sl = slice(i0, i1)

    return dict(
        ratio_bare = ratio_bare[sl],
        ratio_corr = ratio_corr[sl],
        h1nn_ratio = h1nn_ratio[sl],
        A_bare     = A_bare[sl],
        A_corr     = A_corr[sl],
        gap_min    = gap_min[sl],
        gap_right  = gap_right[sl],
    )


def run_ensemble(N, n_ens=100, base_seed=0):
    """N×n_ens 앙상블 실행 후 통계 집계."""
    t0 = time.time()
    print(f"\n{'='*55}")
    print(f"  N = {N:4d},  ensemble = {n_ens}")
    print(f"{'='*55}")
    sys.stdout.flush()

    collectors = {k: [] for k in
                  ['ratio_bare','ratio_corr','h1nn_ratio',
                   'A_bare','A_corr','gap_min','gap_right']}

    for ens in range(n_ens):
        seed = base_seed + ens * 7919 + N
        lambdas = generate_gue_eigvals(N, seed=seed)
        d = compute_metrics_single(lambdas)
        for k in collectors:
            collectors[k].extend(d[k].tolist())

        if (ens + 1) % 20 == 0:
            print(f"  ens {ens+1:3d}/{n_ens} done  [{time.time()-t0:.1f}s]")
            sys.stdout.flush()

    # ─── numpy 변환 + NaN 필터 ───
    arrs = {k: np.array(v) for k, v in collectors.items()}
    mask_b = np.isfinite(arrs['ratio_bare']) & np.isfinite(arrs['gap_min']) & np.isfinite(arrs['A_bare'])
    mask_c = np.isfinite(arrs['ratio_corr']) & np.isfinite(arrs['gap_min']) & np.isfinite(arrs['A_corr'])

    rb = arrs['ratio_bare'][mask_b]
    rc = arrs['ratio_corr'][mask_c]
    gm_b = arrs['gap_min'][mask_b]
    gm_c = arrs['gap_min'][mask_c]
    Ab = arrs['A_bare'][mask_b]
    Ac = arrs['A_corr'][mask_c]
    h1nn = arrs['h1nn_ratio']
    h1nn = h1nn[np.isfinite(h1nn)]

    # ─── 통계 ───
    def safe_pearson(x, y):
        if len(x) < 10:
            return np.nan, np.nan
        return stats.pearsonr(x, y)

    rho_Ab_gm, p_Ab_gm = safe_pearson(Ab, gm_b)
    rho_Ac_gm, p_Ac_gm = safe_pearson(Ac, gm_c)

    res = dict(
        N=N, n_ens=n_ens,
        n_bare=len(rb), n_corr=len(rc),
        mean_bare   = float(np.mean(rb))   if len(rb) > 0 else np.nan,
        median_bare = float(np.median(rb)) if len(rb) > 0 else np.nan,
        std_bare    = float(np.std(rb))    if len(rb) > 0 else np.nan,
        mean_corr   = float(np.mean(rc))   if len(rc) > 0 else np.nan,
        median_corr = float(np.median(rc)) if len(rc) > 0 else np.nan,
        std_corr    = float(np.std(rc))    if len(rc) > 0 else np.nan,
        h1nn_mean   = float(np.mean(h1nn)) if len(h1nn) > 0 else np.nan,
        rho_Ab_gm=float(rho_Ab_gm), p_Ab_gm=float(p_Ab_gm),
        rho_Ac_gm=float(rho_Ac_gm), p_Ac_gm=float(p_Ac_gm),
        elapsed=time.time()-t0,
    )

    print(f"  2H₁/A_bare:  mean={res['mean_bare']:.4f}  median={res['median_bare']:.4f}  std={res['std_bare']:.4f}")
    print(f"  2H₁/A_corr:  mean={res['mean_corr']:.4f}  median={res['median_corr']:.4f}  std={res['std_corr']:.4f}")
    print(f"  H₁^NN/H₁:   mean={res['h1nn_mean']:.4f}")
    print(f"  ρ(A_bare, gap_min) = {rho_Ab_gm:.4f}  p={p_Ab_gm:.3e}  n={len(rb)}")
    print(f"  ρ(A_corr, gap_min) = {rho_Ac_gm:.4f}  p={p_Ac_gm:.3e}  n={len(rc)}")
    print(f"  소요: {res['elapsed']:.1f}s")
    sys.stdout.flush()

    return res


def main():
    t_total = time.time()

    # ─── 실험 설정 ───
    N_values  = [50, 100, 200, 500]
    n_ens     = 100
    out_path  = "/home/k0who029/Desktop/gdl_unified/results/gue_2H1A_ratio_c276.txt"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    print("=" * 65)
    print("C-276: GUE 2H₁/A 비율 계산 + A-gap 상관 측정")
    print(f"N={N_values}, n_ens={n_ens}")
    print(f"목표: 2H₁/A ≈ 2/3 = 0.6667 ? (L-함수 d≥4 관측: ~0.66)")
    print("=" * 65)
    sys.stdout.flush()

    all_results = {}
    for N in N_values:
        res = run_ensemble(N, n_ens=n_ens)
        all_results[N] = res

    # ─── 결과 파일 작성 ───
    lines = []
    lines.append("=" * 70)
    lines.append("C-276: GUE 랜덤 행렬 2H₁/A 비율 + A-gap 상관")
    lines.append(f"실행일시: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"총 소요: {time.time()-t_total:.1f}s")
    lines.append("=" * 70)
    lines.append("")
    lines.append("【설정】")
    lines.append(f"  N_values = {N_values}, n_ens = {n_ens}")
    lines.append(f"  GUE: H=(A+A†)/(2√N), A_ij~CN(0,1), 고유값∈[-2,2]")
    lines.append(f"  B_smooth = N·λ/2  (반원 힐버트 변환)")
    lines.append(f"  A_bare = S₁² + 2H₁,  A_corr = (S₁−B_smooth)² + 2H₁")
    lines.append(f"  내부 고유값: 양 끝 10% 제외")
    lines.append("")
    lines.append("【주요 결과: 2H₁/A 비율】")
    lines.append("")
    lines.append(f"{'N':>6}  {'bare mean':>10}  {'bare med':>9}  {'corr mean':>10}  {'corr med':>9}  {'H1NN/H1':>8}  {'n':>7}")
    lines.append("-" * 70)
    for N in N_values:
        r = all_results[N]
        lines.append(
            f"{N:>6}  {r['mean_bare']:>10.4f}  {r['median_bare']:>9.4f}  "
            f"{r['mean_corr']:>10.4f}  {r['median_corr']:>9.4f}  "
            f"{r['h1nn_mean']:>8.4f}  {r['n_bare']:>7}"
        )
    lines.append("")
    lines.append(f"  참조: L-함수 d≥4 관측값  ≈ 0.661–0.676")
    lines.append(f"  참조: 2/3 = {2/3:.4f}")
    lines.append("")
    lines.append("【A-gap 현상: ρ(A, gap_min_GUE)】")
    lines.append("")
    lines.append(f"{'N':>6}  {'ρ(A_bare,gm)':>13}  {'p':>10}  {'ρ(A_corr,gm)':>13}  {'p':>10}")
    lines.append("-" * 65)
    for N in N_values:
        r = all_results[N]
        lines.append(
            f"{N:>6}  {r['rho_Ab_gm']:>13.4f}  {r['p_Ab_gm']:>10.3e}  "
            f"{r['rho_Ac_gm']:>13.4f}  {r['p_Ac_gm']:>10.3e}"
        )
    lines.append("")
    lines.append(f"  참조: L-함수 d=6 ρ(A_L, gap_min_GUE) = -0.4225  (p=2.4e-4)")
    lines.append(f"  참조: L-함수 d=5 ρ(A_L, gap_min_GUE) = -0.485   (p=3.7e-9)")
    lines.append("")

    # ─── 판정 ───
    lines.append("【판정】")
    lines.append("")

    # 2H₁/A_bare N→∞ 추세
    bares = [all_results[N]['mean_bare'] for N in N_values]
    corrs = [all_results[N]['mean_corr'] for N in N_values]
    lines.append(f"  2H₁/A_bare (N=50→500): {bares[0]:.4f} → {bares[-1]:.4f}")
    lines.append(f"  2H₁/A_corr (N=50→500): {corrs[0]:.4f} → {corrs[-1]:.4f}")

    # bare 대 2/3
    target = 2.0/3.0
    diff_bare = abs(bares[-1] - target)
    diff_corr = abs(corrs[-1] - target)
    lines.append(f"  |2H₁/A_bare(N=500) − 2/3| = {diff_bare:.4f}")
    lines.append(f"  |2H₁/A_corr(N=500) − 2/3| = {diff_corr:.4f}")

    # A-gap in GUE
    rho_largest = all_results[N_values[-1]]['rho_Ab_gm']
    p_largest   = all_results[N_values[-1]]['p_Ab_gm']
    lines.append(f"  ρ(A_bare, gap_min) @ N=500: {rho_largest:.4f}  (p={p_largest:.3e})")

    lines.append("")
    if diff_bare < 0.02:
        lines.append("  ★★★★ GUE_BARE: 2H₁/A_bare → 2/3 확인 (L-함수 d≥4 GUE 보편성)")
    elif diff_bare < 0.05:
        lines.append(f"  ★★★  GUE_BARE: 2H₁/A_bare ≈ {bares[-1]:.4f}, 2/3에 근접 (경계)")
    else:
        lines.append(f"  ★★   GUE_BARE: 2H₁/A_bare = {bares[-1]:.4f} (2/3과 유의한 차이)")

    if diff_corr < 0.02:
        lines.append("  ★★★★ GUE_CORR: 2H₁/A_corr → 2/3 확인")
    elif diff_corr < 0.05:
        lines.append(f"  ★★★  GUE_CORR: 2H₁/A_corr ≈ {corrs[-1]:.4f}")
    else:
        lines.append(f"  ★★   GUE_CORR: 2H₁/A_corr = {corrs[-1]:.4f}")

    if p_largest < 0.05 and rho_largest < -0.1:
        lines.append(f"  ★★★★★ A-GAP: GUE에서도 ρ(A,gap_min) = {rho_largest:.4f} < 0 (유의)")
        lines.append("         → A-gap 현상은 GUE 보편 현상 (산술 구조 불필요)")
    elif p_largest < 0.05 and rho_largest < 0:
        lines.append(f"  ★★★  A-GAP: GUE ρ = {rho_largest:.4f} 약한 음의 상관 (경계)")
    else:
        lines.append(f"  ★★   A-GAP: GUE ρ = {rho_largest:.4f} (비유의 또는 양수)")
        lines.append("         → A-gap 현상은 L-함수 산술 구조의 산물")
    lines.append("")

    # 상세 N별
    lines.append("【N별 상세】")
    for N in N_values:
        r = all_results[N]
        lines.append(f"  N={N}: bare={r['mean_bare']:.4f}±{r['std_bare']:.4f}, "
                     f"corr={r['mean_corr']:.4f}±{r['std_corr']:.4f}, "
                     f"H1NN/H1={r['h1nn_mean']:.4f}, "
                     f"ρ_bare={r['rho_Ab_gm']:.4f}(p={r['p_Ab_gm']:.2e}), "
                     f"ρ_corr={r['rho_Ac_gm']:.4f}(p={r['p_Ac_gm']:.2e}), "
                     f"n={r['n_bare']}, t={r['elapsed']:.1f}s")
    lines.append("")
    lines.append("=" * 70)

    txt = "\n".join(lines)
    with open(out_path, "w") as f:
        f.write(txt)

    print("\n" + txt)
    print(f"\n✅ 결과 저장: {out_path}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
