#!/usr/bin/env python3
"""
C-277: GUE A-gap N-확장 + Spearman 재측정
비평가 Q1,Q2,Q3 동시 해소:
  Q1: N=1000,2000에서 ρ(A, gap_min)은?
  Q2: Spearman으로 측정하면 |ρ| 감소 추세 동일한가?
  Q3: edge exclusion 비율(5%,10%,15%) 변화 시 결과 안정한가?

C-276 기반. 핵심 변경:
  - N=50,100,200,500,1000,2000 (6개)
  - Pearson + Spearman 동시 계산
  - edge exclusion 감도 분석 (N=500에서 5%,10%,15%)
  - n_ens: N≤500=100, N=1000=50, N=2000=20

결과: results/gue_agap_extension_c277.txt
"""

import numpy as np
from scipy import stats
import os
import sys
import time


def generate_gue_eigvals(N, seed=None):
    """GUE(N) 행렬 생성 후 실수 고유값 반환 (정렬됨)."""
    rng = np.random.default_rng(seed)
    real = rng.standard_normal((N, N))
    imag = rng.standard_normal((N, N))
    A = (real + 1j * imag) / np.sqrt(2.0)
    H = (A + A.conj().T) / (2.0 * np.sqrt(N))
    eigvals = np.linalg.eigvalsh(H)
    return eigvals.real


def semicircle_cdf(x):
    """반원 분포 CDF."""
    x = np.clip(x, -2.0 + 1e-12, 2.0 - 1e-12)
    return (np.arcsin(x / 2.0) + (x / 2.0) * np.sqrt(np.maximum(0.0, 1.0 - x**2 / 4.0))) / np.pi + 0.5


def unfold(lambdas):
    return len(lambdas) * semicircle_cdf(lambdas)


def compute_metrics_single(lambdas, edge_frac=0.10):
    """단일 GUE 실현에서 A, H₁, gap 계산. edge_frac으로 양 끝 제외 비율 조정."""
    N = len(lambdas)
    xi = unfold(lambdas)

    diff = lambdas[:, None] - lambdas[None, :]
    np.fill_diagonal(diff, np.nan)
    inv1 = 1.0 / diff
    inv2 = 1.0 / diff**2

    S1 = np.nansum(inv1, axis=1)
    H1 = np.nansum(inv2, axis=1)

    B_smooth = N * lambdas / 2.0
    B_corr = S1 - B_smooth

    A_bare = S1**2 + 2.0 * H1
    A_corr = B_corr**2 + 2.0 * H1

    ratio_bare = np.where(A_bare > 0, 2.0 * H1 / A_bare, np.nan)

    # H₁^{NN}
    gaps = np.diff(lambdas)
    inv_g2 = 1.0 / gaps**2
    H1_nn = np.zeros(N)
    H1_nn[1:] += inv_g2
    H1_nn[:-1] += inv_g2
    h1nn_ratio = np.where(H1 > 0, H1_nn / H1, np.nan)

    # unfolded gap
    xi_gaps = np.diff(xi)
    gap_left = np.concatenate([[np.nan], xi_gaps])
    gap_right = np.concatenate([xi_gaps, [np.nan]])
    gap_min = np.minimum(gap_left, gap_right)

    # 내부 선택
    i0 = max(1, int(N * edge_frac))
    i1 = min(N - 1, int(N * (1.0 - edge_frac)))
    sl = slice(i0, i1)

    return dict(
        ratio_bare=ratio_bare[sl],
        h1nn_ratio=h1nn_ratio[sl],
        A_bare=A_bare[sl],
        A_corr=A_corr[sl],
        gap_min=gap_min[sl],
    )


def run_ensemble(N, n_ens, edge_frac=0.10, base_seed=0):
    """N×n_ens 앙상블 실행. Pearson + Spearman 동시 계산."""
    t0 = time.time()
    print(f"\n{'='*60}")
    print(f"  N = {N:5d},  ens = {n_ens},  edge_frac = {edge_frac}")
    print(f"{'='*60}")
    sys.stdout.flush()

    keys = ['ratio_bare', 'h1nn_ratio', 'A_bare', 'A_corr', 'gap_min']
    collectors = {k: [] for k in keys}

    for ens in range(n_ens):
        seed = base_seed + ens * 7919 + N
        lambdas = generate_gue_eigvals(N, seed=seed)
        d = compute_metrics_single(lambdas, edge_frac=edge_frac)
        for k in keys:
            collectors[k].extend(d[k].tolist())

        if (ens + 1) % max(1, n_ens // 5) == 0:
            print(f"  ens {ens+1:3d}/{n_ens} done  [{time.time()-t0:.1f}s]")
            sys.stdout.flush()

    arrs = {k: np.array(v) for k, v in collectors.items()}
    mask = (np.isfinite(arrs['ratio_bare']) &
            np.isfinite(arrs['gap_min']) &
            np.isfinite(arrs['A_bare']))

    rb = arrs['ratio_bare'][mask]
    Ab = arrs['A_bare'][mask]
    gm = arrs['gap_min'][mask]
    h1nn = arrs['h1nn_ratio']
    h1nn = h1nn[np.isfinite(h1nn)]

    # Pearson
    if len(Ab) > 10:
        rho_p, p_p = stats.pearsonr(Ab, gm)
    else:
        rho_p, p_p = np.nan, np.nan

    # Spearman
    if len(Ab) > 10:
        rho_s, p_s = stats.spearmanr(Ab, gm)
    else:
        rho_s, p_s = np.nan, np.nan

    elapsed = time.time() - t0
    res = dict(
        N=N, n_ens=n_ens, edge_frac=edge_frac, n=int(len(Ab)),
        mean_bare=float(np.mean(rb)),
        median_bare=float(np.median(rb)),
        h1nn_mean=float(np.mean(h1nn)),
        rho_pearson=float(rho_p), p_pearson=float(p_p),
        rho_spearman=float(rho_s), p_spearman=float(p_s),
        elapsed=elapsed,
    )

    print(f"  2H₁/A_bare: mean={res['mean_bare']:.4f}  median={res['median_bare']:.4f}")
    print(f"  H₁^NN/H₁:  mean={res['h1nn_mean']:.4f}")
    print(f"  ρ_Pearson(A,gm) = {rho_p:.4f}  p={p_p:.3e}")
    print(f"  ρ_Spearman(A,gm)= {rho_s:.4f}  p={p_s:.3e}")
    print(f"  n={len(Ab)}, 소요: {elapsed:.1f}s")
    sys.stdout.flush()

    return res


def main():
    t_total = time.time()
    out_path = "/home/k0who029/Desktop/gdl_unified/results/gue_agap_extension_c277.txt"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    print("=" * 70)
    print("C-277: GUE A-gap N-확장 + Spearman 재측정")
    print("비평가 Q1(N확장), Q2(Spearman), Q3(edge 감도) 동시 해소")
    print("=" * 70)
    sys.stdout.flush()

    # ─── Part 1: N 확장 (edge=10%) ───
    configs = [
        (50, 100), (100, 100), (200, 100), (500, 100),
        (1000, 50), (2000, 20),
    ]

    main_results = {}
    for N, n_ens in configs:
        res = run_ensemble(N, n_ens, edge_frac=0.10)
        main_results[N] = res

    # ─── Part 2: Edge 감도 분석 (N=500) ───
    edge_results = {}
    for ef in [0.05, 0.10, 0.15]:
        res = run_ensemble(500, 100, edge_frac=ef, base_seed=0)
        edge_results[ef] = res

    # ─── 결과 파일 작성 ───
    lines = []
    lines.append("=" * 70)
    lines.append("C-277: GUE A-gap N-확장 + Spearman 재측정")
    lines.append(f"실행일시: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"총 소요: {time.time()-t_total:.1f}s")
    lines.append("=" * 70)
    lines.append("")

    lines.append("【설정】")
    lines.append(f"  N = {[c[0] for c in configs]}")
    lines.append(f"  n_ens = {dict(configs)}")
    lines.append(f"  Edge exclusion = 10% (기본), 감도: 5%/10%/15% (N=500)")
    lines.append(f"  상관 측정: Pearson + Spearman 동시")
    lines.append("")

    # ─── Part 1 결과 ───
    lines.append("【Part 1: N-확장 — ρ(A_bare, gap_min) 추세】")
    lines.append("")
    lines.append(f"{'N':>6}  {'ρ_Pearson':>10}  {'p_P':>10}  {'ρ_Spearman':>11}  {'p_S':>10}  {'2H₁/A':>7}  {'H1NN/H1':>8}  {'n':>7}")
    lines.append("-" * 80)
    for N in sorted(main_results.keys()):
        r = main_results[N]
        lines.append(
            f"{N:>6}  {r['rho_pearson']:>10.4f}  {r['p_pearson']:>10.3e}  "
            f"{r['rho_spearman']:>11.4f}  {r['p_spearman']:>10.3e}  "
            f"{r['mean_bare']:>7.4f}  {r['h1nn_mean']:>8.4f}  {r['n']:>7}"
        )
    lines.append("")

    # N→∞ 추세 분석
    Ns = sorted(main_results.keys())
    rho_p_vals = [main_results[N]['rho_pearson'] for N in Ns]
    rho_s_vals = [main_results[N]['rho_spearman'] for N in Ns]

    lines.append("  N→∞ ρ_Pearson 추세:  " + " → ".join(f"{v:.4f}" for v in rho_p_vals))
    lines.append("  N→∞ ρ_Spearman 추세: " + " → ".join(f"{v:.4f}" for v in rho_s_vals))
    lines.append("")

    # 판정: N=1000,2000에서 |ρ|가 증가/유지/감소?
    rho_500_s = main_results[500]['rho_spearman']
    rho_1000_s = main_results[1000]['rho_spearman']
    rho_2000_s = main_results[2000]['rho_spearman']

    if abs(rho_2000_s) >= abs(rho_500_s) * 0.8:
        trend = "안정/증가"
        verdict_q1 = "A-gap은 N→∞에서도 존재 (GUE 보편)"
        stars_q1 = "★★★★★"
    elif abs(rho_2000_s) >= abs(rho_500_s) * 0.5:
        trend = "감쇠 (느린)"
        verdict_q1 = "A-gap 감쇠 진행 중, N→∞ 거동 불확실"
        stars_q1 = "★★★"
    else:
        trend = "급감"
        verdict_q1 = "A-gap은 유한-N 효과일 가능성"
        stars_q1 = "★★"

    lines.append(f"  Q1 판정: N=500→2000 Spearman 추세 = {trend}")
    lines.append(f"           {stars_q1} {verdict_q1}")
    lines.append("")

    # ─── Part 2 결과 ───
    lines.append("【Part 2: Edge Exclusion 감도 (N=500)】")
    lines.append("")
    lines.append(f"{'edge%':>6}  {'ρ_Pearson':>10}  {'ρ_Spearman':>11}  {'2H₁/A':>7}  {'n':>7}")
    lines.append("-" * 50)
    for ef in [0.05, 0.10, 0.15]:
        r = edge_results[ef]
        lines.append(
            f"{ef*100:>5.0f}%  {r['rho_pearson']:>10.4f}  {r['rho_spearman']:>11.4f}  "
            f"{r['mean_bare']:>7.4f}  {r['n']:>7}"
        )
    lines.append("")

    rho_5 = edge_results[0.05]['rho_spearman']
    rho_10 = edge_results[0.10]['rho_spearman']
    rho_15 = edge_results[0.15]['rho_spearman']
    span = max(abs(rho_5), abs(rho_10), abs(rho_15)) - min(abs(rho_5), abs(rho_10), abs(rho_15))

    if span < 0.05:
        lines.append(f"  Q3 판정: edge 감도 span = {span:.4f} — 안정 ✅")
    else:
        lines.append(f"  Q3 판정: edge 감도 span = {span:.4f} — 불안정 ⚠️")
    lines.append("")

    # ─── 종합 판정 ───
    lines.append("【종합 판정】")
    lines.append("")

    # Pearson vs Spearman 비교
    diff_ps = abs(main_results[500]['rho_pearson'] - main_results[500]['rho_spearman'])
    lines.append(f"  Q2: Pearson vs Spearman 차이 (N=500) = {diff_ps:.4f}")
    if diff_ps < 0.05:
        lines.append("      → 일치 — Pearson/Spearman 무관하게 A-gap 존재")
    else:
        lines.append("      → 차이 존재 — 비선형 관계 가능성")
    lines.append("")

    # 참조값
    lines.append("  참조: L-함수 ρ_Spearman 범위: -0.42 ~ -0.58 (GL(1)~GL(6))")
    lines.append(f"  GUE ρ_Spearman @ N=2000: {rho_2000_s:.4f}")
    lines.append("")

    # 최종
    if abs(rho_2000_s) > 0.15 and main_results[2000]['p_spearman'] < 0.05:
        lines.append("  ★★★★★ GUE A-gap: N=2000에서도 유의한 반상관 존재")
        lines.append("          → A-gap 현상의 GUE 보편성 강화 확인")
    elif abs(rho_2000_s) > 0.05:
        lines.append("  ★★★   GUE A-gap: N=2000에서 약한 반상관")
    else:
        lines.append("  ★★    GUE A-gap: N=2000에서 소멸 추세")
    lines.append("")

    # H₁^{NN}/H₁ 추세
    h1nn_vals = [main_results[N]['h1nn_mean'] for N in Ns]
    lines.append(f"  H₁^NN/H₁ 추세: " + " → ".join(f"{v:.4f}" for v in h1nn_vals))
    lines.append(f"  N=2000: {h1nn_vals[-1]:.4f} (2/3={2/3:.4f})")
    lines.append("")

    lines.append("【N별 상세】")
    for N in Ns:
        r = main_results[N]
        lines.append(
            f"  N={N}: bare={r['mean_bare']:.4f}, H1NN={r['h1nn_mean']:.4f}, "
            f"ρ_P={r['rho_pearson']:.4f}(p={r['p_pearson']:.2e}), "
            f"ρ_S={r['rho_spearman']:.4f}(p={r['p_spearman']:.2e}), "
            f"n={r['n']}, t={r['elapsed']:.1f}s"
        )
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
