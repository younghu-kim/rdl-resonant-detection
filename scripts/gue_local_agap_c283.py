#!/usr/bin/env python3
"""
C-283: GUE ±300 이웃 제한 + trim 재측정
목적: GUE A_bare를 ζ(s) C-282b와 동일한 조건(±300 이웃, trim 20%)으로 재측정.
     A_bare_full vs A_bare_local 비교 → 방법론 불일치 해소.

핵심 변경 (C-277 → C-283):
  1. A_bare 계산 시 전체 고유값(full) + ±N_MAX 이웃(local) 모두 출력
  2. N_MAX = 300 (ζ(s) C-282b와 동일)
  3. trim 20% (양쪽 20% 제거, 중앙 60%) — C-277의 edge 10%와 다름
  4. N = [200, 500, 1000, 2000], 앙상블 = {200:100, 500:50, 1000:30, 2000:20}
  5. Spearman만 (B-47: Pearson은 아티팩트)

성공 기준:
  - C-277 재현: ρ_S(A_full, gap_min) ≈ -0.857 ± 0.01 (N=2000)
  - 핵심 측정: ρ_S(A_local, gap_min) → 값 자체가 정보

결과: results/gue_local_agap_c283.txt
"""

import numpy as np
from scipy import stats
import os
import sys
import time

N_MAX = 300  # ±300 이웃 (ζ(s) C-282b와 동일)
TRIM_FRAC = 0.20  # 양쪽 20% 제거 → 중앙 60%


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
    """반원 분포 CDF (unfolding용)."""
    x = np.clip(x, -2.0 + 1e-12, 2.0 - 1e-12)
    return (np.arcsin(x / 2.0) + (x / 2.0) * np.sqrt(np.maximum(0.0, 1.0 - x**2 / 4.0))) / np.pi + 0.5


def unfold(lambdas):
    """반원 CDF 기반 unfolding."""
    return len(lambdas) * semicircle_cdf(lambdas)


def compute_metrics_single(lambdas):
    """단일 GUE 실현에서 A_full, A_local, gap_min 계산.

    trim 20% (양쪽 20% 제거, 중앙 60%) 적용.
    """
    N = len(lambdas)
    xi = unfold(lambdas)

    # trim 범위 (중앙 60%)
    i0 = int(N * TRIM_FRAC)
    i1 = int(N * (1.0 - TRIM_FRAC))

    # unfolded gap (gap_min = min(left, right))
    xi_gaps = np.diff(xi)

    # 결과 저장
    A_full_arr = []
    A_local_arr = []
    gap_min_arr = []

    for i in range(i0, i1):
        # gap_min 계산
        gap_left = xi[i] - xi[i - 1] if i > 0 else np.nan
        gap_right = xi[i + 1] - xi[i] if i < N - 1 else np.nan

        if np.isnan(gap_left) or np.isnan(gap_right):
            continue
        gm = min(gap_left, gap_right)
        if gm <= 0 or not np.isfinite(gm):
            continue

        # === A_full: 전체 고유값 사용 (C-277 방식) ===
        diffs_full = lambdas[i] - lambdas
        diffs_full_nz = np.concatenate([diffs_full[:i], diffs_full[i+1:]])
        # 0 방지
        mask_ok = np.abs(diffs_full_nz) > 1e-15
        if np.sum(mask_ok) < 10:
            continue
        inv1_full = 1.0 / diffs_full_nz[mask_ok]
        inv2_full = 1.0 / diffs_full_nz[mask_ok]**2
        S1_full = np.sum(inv1_full)
        H1_full = np.sum(inv2_full)
        A_full = S1_full**2 + 2.0 * H1_full

        # === A_local: ±N_MAX 이웃만 사용 (ζ(s) C-282b 방식) ===
        lo = max(0, i - N_MAX)
        hi = min(N - 1, i + N_MAX)
        # i 제외한 [lo, hi] 범위의 인덱스
        local_idx = list(range(lo, i)) + list(range(i + 1, hi + 1))
        if len(local_idx) < 10:
            continue
        diffs_local = lambdas[i] - lambdas[local_idx]
        mask_ok_l = np.abs(diffs_local) > 1e-15
        if np.sum(mask_ok_l) < 10:
            continue
        inv1_local = 1.0 / diffs_local[mask_ok_l]
        inv2_local = 1.0 / diffs_local[mask_ok_l]**2
        S1_local = np.sum(inv1_local)
        H1_local = np.sum(inv2_local)
        A_local = S1_local**2 + 2.0 * H1_local

        if not (np.isfinite(A_full) and np.isfinite(A_local) and A_full > 0 and A_local > 0):
            continue

        A_full_arr.append(A_full)
        A_local_arr.append(A_local)
        gap_min_arr.append(gm)

    return (np.array(A_full_arr), np.array(A_local_arr), np.array(gap_min_arr))


def run_ensemble(N, n_ens, base_seed=42):
    """N×n_ens 앙상블 실행."""
    t0 = time.time()
    print(f"\n{'='*60}")
    print(f"  N = {N:5d},  ens = {n_ens},  trim = {TRIM_FRAC*100:.0f}%,  N_MAX = {N_MAX}")
    print(f"{'='*60}")
    sys.stdout.flush()

    all_A_full = []
    all_A_local = []
    all_gap_min = []

    for ens in range(n_ens):
        seed = base_seed + ens * 7919 + N
        lambdas = generate_gue_eigvals(N, seed=seed)
        A_f, A_l, gm = compute_metrics_single(lambdas)
        all_A_full.extend(A_f.tolist())
        all_A_local.extend(A_l.tolist())
        all_gap_min.extend(gm.tolist())

        if (ens + 1) % max(1, n_ens // 5) == 0:
            print(f"  ens {ens+1:3d}/{n_ens} done  [{time.time()-t0:.1f}s]")
            sys.stdout.flush()

    A_f = np.array(all_A_full)
    A_l = np.array(all_A_local)
    gm = np.array(all_gap_min)

    n = len(gm)

    # Spearman 계산
    if n > 10:
        rho_full, p_full = stats.spearmanr(A_f, gm)
        rho_local, p_local = stats.spearmanr(A_l, gm)
    else:
        rho_full, p_full = np.nan, np.nan
        rho_local, p_local = np.nan, np.nan

    # 2H₁/A 비율 (A_full 기반으로 추정 — S1²+2H1에서 2H1/A)
    # 실제로는 개별 H1, S1을 저장해야 하지만, 이 비교에서는 ρ가 핵심

    # N=200일 때 ±300이면 사실상 전체
    effective_local = min(2 * N_MAX, N - 1)
    is_effectively_full = effective_local >= N - 1

    elapsed = time.time() - t0

    print(f"  n = {n}")
    print(f"  ρ_S(A_full,  gap_min) = {rho_full:.4f}  (p={p_full:.3e})")
    print(f"  ρ_S(A_local, gap_min) = {rho_local:.4f}  (p={p_local:.3e})")
    print(f"  |Δρ| = {abs(abs(rho_local) - abs(rho_full)):.4f}")
    if is_effectively_full:
        print(f"  ⚠️ N={N} < 2*N_MAX={2*N_MAX} → local ≈ full (의미 제한)")
    print(f"  소요: {elapsed:.1f}s")
    sys.stdout.flush()

    return dict(
        N=N, n_ens=n_ens, n=n,
        rho_full=float(rho_full), p_full=float(p_full),
        rho_local=float(rho_local), p_local=float(p_local),
        is_effectively_full=is_effectively_full,
        elapsed=elapsed,
        mean_gap_min=float(np.mean(gm)),
        std_gap_min=float(np.std(gm)),
        mean_A_full=float(np.mean(A_f)),
        mean_A_local=float(np.mean(A_l)),
    )


def main():
    t_total = time.time()
    out_path = "/home/k0who029/Desktop/gdl_unified/results/gue_local_agap_c283.txt"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    print("=" * 70)
    print("C-283: GUE ±300 이웃 제한 + trim 20% 재측정")
    print("목적: A_bare_full vs A_bare_local(±300) → 방법론 불일치 해소")
    print(f"N_MAX = {N_MAX}, trim = {TRIM_FRAC*100:.0f}%")
    print("=" * 70)
    sys.stdout.flush()

    configs = [
        (200, 100),
        (500, 50),
        (1000, 30),
        (2000, 20),
    ]

    results = {}
    for N, n_ens in configs:
        res = run_ensemble(N, n_ens)
        results[N] = res

    # ─── 결과 파일 작성 ───
    lines = []
    lines.append("=" * 70)
    lines.append("C-283: GUE ±300 이웃 제한 + trim 20% 재측정")
    lines.append(f"실행일시: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"총 소요: {time.time()-t_total:.1f}s")
    lines.append("=" * 70)
    lines.append("")

    lines.append("【설정】")
    lines.append(f"  N_MAX = {N_MAX} (±300 이웃, ζ(s) C-282b와 동일)")
    lines.append(f"  trim = {TRIM_FRAC*100:.0f}% (양쪽 {TRIM_FRAC*100:.0f}% 제거, 중앙 {(1-2*TRIM_FRAC)*100:.0f}%)")
    lines.append(f"  N = {[c[0] for c in configs]}")
    lines.append(f"  n_ens = {dict(configs)}")
    lines.append(f"  상관 측정: Spearman only (B-47)")
    lines.append("")

    # ─── 핵심 비교표 ───
    lines.append("【핵심 비교표: A_full vs A_local(±300)】")
    lines.append("")
    lines.append(f"{'N':>6}  {'ρ_S(A_full)':>12}  {'p_full':>10}  {'ρ_S(A_local)':>13}  {'p_local':>10}  {'|Δρ|':>7}  {'n':>7}  {'비고':>10}")
    lines.append("-" * 90)
    for N in sorted(results.keys()):
        r = results[N]
        delta = abs(abs(r['rho_local']) - abs(r['rho_full']))
        note = "≈full" if r['is_effectively_full'] else ""
        lines.append(
            f"{N:>6}  {r['rho_full']:>12.4f}  {r['p_full']:>10.3e}  "
            f"{r['rho_local']:>13.4f}  {r['p_local']:>10.3e}  "
            f"{delta:>7.4f}  {r['n']:>7}  {note:>10}"
        )
    lines.append("")

    # ─── 판정 ───
    lines.append("【판정】")
    lines.append("")

    # C-277 재현 검증
    r2000 = results[2000]
    c277_target = -0.857
    c277_diff = abs(r2000['rho_full'] - c277_target)
    if c277_diff < 0.03:
        lines.append(f"  ✅ C-277 재현: ρ_S(A_full, N=2000) = {r2000['rho_full']:.4f} ≈ -0.857 (차이 {c277_diff:.4f})")
    else:
        lines.append(f"  ⚠️ C-277 차이: ρ_S(A_full, N=2000) = {r2000['rho_full']:.4f} vs -0.857 (차이 {c277_diff:.4f})")
        lines.append(f"     (trim 방식 차이: C-277=edge 10%, C-283=trim 20% → 차이 가능)")
    lines.append("")

    # 핵심 판정: local vs full
    rho_full_2k = r2000['rho_full']
    rho_local_2k = r2000['rho_local']
    delta_2k = abs(rho_local_2k) - abs(rho_full_2k)

    lines.append(f"  핵심: ρ_S(A_local, N=2000) = {rho_local_2k:.4f}")
    lines.append(f"        ρ_S(A_full,  N=2000) = {rho_full_2k:.4f}")
    lines.append(f"        Δ = |ρ_local| - |ρ_full| = {delta_2k:+.4f}")
    lines.append("")

    if delta_2k > 0.05:
        lines.append("  → 해석: ρ_local > ρ_full → 국소 합이 gap 예측력이 더 강함")
        lines.append("    GUE에서도 동일 패턴 → ζ(s) C-282b의 -0.93은 국소 효과 포함")
    elif delta_2k > -0.02:
        lines.append("  → 해석: ρ_local ≈ ρ_full → 국소 제한 효과 미미")
        lines.append("    ζ(s) -0.93 > GUE -0.86은 L-함수 구조의 실재적 효과")
    else:
        lines.append("  → 해석: ρ_local < ρ_full → 국소 제한이 오히려 약화")
    lines.append("")

    # ζ(s)와의 비교
    lines.append("【ζ(s) C-282b와의 비교】")
    lines.append("")
    lines.append(f"  ζ(s) A_bare(±300, trim 20%): ρ_S = -0.929")
    lines.append(f"  GUE  A_local(±300, trim 20%): ρ_S = {rho_local_2k:.4f}")
    lines.append(f"  GUE  A_full (전체,  trim 20%): ρ_S = {rho_full_2k:.4f}")
    lines.append("")

    zeta_vs_local = abs(-0.929) - abs(rho_local_2k)
    zeta_vs_full = abs(-0.929) - abs(rho_full_2k)
    lines.append(f"  ζ(s) vs GUE_local: Δ = {zeta_vs_local:+.4f}")
    lines.append(f"  ζ(s) vs GUE_full:  Δ = {zeta_vs_full:+.4f}")
    lines.append("")

    if abs(rho_local_2k) > 0.90:
        lines.append("  ★ GUE_local ≈ -0.93 → 'ζ(s) > GUE' 소멸. 국소성의 효과.")
    elif abs(rho_local_2k) > 0.83 and abs(rho_local_2k) <= 0.90:
        lines.append("  ★ GUE_local ≈ -0.87~-0.90 → 국소 효과 존재 + ζ(s) 약간 강함")
    elif abs(rho_local_2k) <= 0.83:
        lines.append("  ★ GUE_local ≈ GUE_full → 국소 제한 무관. ζ(s)가 진짜 GUE 초과.")
    lines.append("")

    # N 의존성
    lines.append("【N 의존성 분석】")
    lines.append("")
    for N in sorted(results.keys()):
        r = results[N]
        lines.append(
            f"  N={N:4d}: A_full ρ={r['rho_full']:.4f}, A_local ρ={r['rho_local']:.4f}, "
            f"|Δ|={abs(abs(r['rho_local'])-abs(r['rho_full'])):.4f}, "
            f"<gap_min>={r['mean_gap_min']:.3f}±{r['std_gap_min']:.3f}, "
            f"<A_full>={r['mean_A_full']:.1f}, <A_local>={r['mean_A_local']:.1f}, "
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
