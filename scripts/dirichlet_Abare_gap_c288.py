#!/usr/bin/env python3
"""
[사이클 #288] 디리클레 L-함수 A_bare gap_min Spearman ρ_S — d=1 보편성 검증

  목적: B-54 경계 탐사 — d-의존적 감쇠의 d=1 내 보편성 확인
    Paper 4 통일표: ζ(s) δ=0.038, GL(2) δ≈0.24
    질문: Dirichlet (d=1, q>1)에서 δ는 ζ(s)와 같은가?
    YES → d-의존성 확정 (conductor 무관)
    NO → conductor/character 의존성 → 더 풍부한 구조

  설계:
    - L-함수: ζ(s) + χ₃(mod 3) + χ₄(mod 4) + χ₅(mod 5) — 모두 d=1
    - T=[5, 500], PARI lfunzeros로 영점 수집
    - A_bare = S₁² + 2H₁ (전체 영점, smooth 미보정)
    - A_L   = (S₁ - B_smooth)² + 2H₁ (Γ smooth 보정)
    - gap_min_GUE: 국소 밀도 정규화
    - Trim 20% (각 끝 10%)
    - GUE_full 기준 = -0.912
    - 동일 T=500에서 ζ(s)도 측정 → 같은 조건 비교

  체크리스트:
    [x] PARI lfuncreate(Mod(g,q)) 사용
    [x] A_bare = S1^2 + 2*H1 (같은편만, GL(d)와 동일)
    [x] B_smooth = (1/2)Im[log(q/π) + ψ((s+a)/2)]  for Dirichlet
    [x] gap_min_GUE = min(gap_r, gap_l) * d_bar
    [x] Spearman (scipy.stats.spearmanr)
    [x] python -u
    [x] trim 20% (각 끝 10%)
    [x] NaN 필터
"""

import sys
import os
import math
import time

import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 50  # t up to 500 → 50자리

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(512 * 10**6)
    pari.set_real_precision(50)
    print("cypari2 OK (512 MB, precision=50)")
except Exception as e:
    print(f"FATAL: cypari2 로드 실패: {e}")
    sys.exit(1)

# ── 설정 ──────────────────────────────────────────────────────────
T_MIN = 5.0
T_MAX = 500.0
TRIM_FRAC = 0.10  # 각 끝에서 10% 제거 (총 20%)
GUE_FULL_REF = -0.912  # Paper 4 기준

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/dirichlet_Abare_gap_c288.txt'
)

# L-함수 정의
LFUNCS = [
    {
        'name': 'ζ(s)',
        'label': 'Riemann zeta',
        'd': 1,
        'q': 1,
        'a': 0,       # even character → a=0
        'mu': [0],
        'pari_init': "lfuncreate(1)",
    },
    {
        'name': 'χ₃ (mod 3)',
        'label': 'Dirichlet mod 3, odd',
        'd': 1,
        'q': 3,
        'a': 1,       # odd character → a=1
        'mu': [1],
        'pari_init': "lfuncreate(Mod(2,3))",
    },
    {
        'name': 'χ₄ (mod 4)',
        'label': 'Dirichlet mod 4, odd',
        'd': 1,
        'q': 4,
        'a': 1,
        'mu': [1],
        'pari_init': "lfuncreate(Mod(3,4))",
    },
    {
        'name': 'χ₅ (mod 5)',
        'label': 'Dirichlet mod 5, complex',
        'd': 1,
        'q': 5,
        'a': 1,
        'mu': [1],
        'pari_init': "lfuncreate(Mod(2,5))",
    },
]


def pf(x):
    """PARI 객체를 float으로 변환"""
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeros(lf, t_max):
    """PARI lfunzeros로 L-함수 영점 수집."""
    name = lf['name']
    print(f"  [{name}] 영점 수집 T=[0, {t_max}] ...")
    t0 = time.time()

    try:
        Lobj = pari(lf['pari_init'])
        Linit = pari.lfuninit(Lobj, [0, int(t_max) + 10])
        zvec = pari.lfunzeros(Linit, t_max)
        n = int(pari.length(zvec))
        zeros = []
        for i in range(n):
            t = float(zvec[i])
            if not math.isnan(t) and t > 0.5:
                zeros.append(t)
        zeros = sorted(zeros)
        print(f"    {n}개 원시, {len(zeros)}개 유효 ({time.time()-t0:.1f}s)")
        return zeros
    except Exception as e:
        print(f"    ERROR: {e}")
        return []


def gamma_smooth_im(gamma_0, lf):
    """
    L-함수의 Γ smooth part 허수부.

    Completed L-function:
      Λ(s) = (N/π^d)^{s/2} ∏ Γ_R(s+μ_j) · L(s)

    Γ_∞'/Γ_∞(s) = Σ_j (1/2)[ψ((s+μ_j)/2) - log(π)] + (1/2)log(N)

    Im at s = 1/2 + iγ:
      B_smooth = Im{ Σ_j (1/2)ψ((1/2+μ_j+iγ)/2) } + 0  (log terms are real)

    d=1일 때:
      ζ(s): μ=0 → B = (1/2)Im[ψ(0.25 + iγ/2)]
      χ odd: μ=1 → B = (1/2)Im[ψ(0.75 + iγ/2)]
    """
    s = mpmath.mpc(0.5, gamma_0)
    total = mpmath.mpc(0)
    for mu in lf['mu']:
        total += mpmath.digamma((s + mu) / 2) - mpmath.log(mpmath.pi)
    total /= 2
    # log(N_cond)/2 contribution (real, doesn't affect Im)
    # total += mpmath.log(lf['q']) / 2  — real, omit
    return float(mpmath.im(total))


def compute_A(zeros, idx):
    """
    영점 idx에서 A_bare 계산 (전체 영점 사용, GL(2)와 동일).
    A_bare = S₁² + 2H₁
    """
    gamma_0 = zeros[idx]
    n_total = len(zeros)

    S1 = 0.0
    H1 = 0.0

    for k in range(n_total):
        if k == idx:
            continue
        dg = gamma_0 - zeros[k]
        if abs(dg) < 1e-15:
            continue
        S1 += 1.0 / dg
        H1 += 1.0 / (dg * dg)

    A_bare = S1 ** 2 + 2.0 * H1
    return S1, H1, A_bare


def analyze_one(lf, all_zeros):
    """하나의 L-함수에 대한 A_bare-gap 분석."""
    name = lf['name']
    n_all = len(all_zeros)

    if n_all < 30:
        print(f"  [{name}] 영점 부족 ({n_all}개) — 건너뜀")
        return None

    # T_MIN 이상 필터
    zeros = [z for z in all_zeros if z >= T_MIN]
    n_range = len(zeros)
    print(f"  [{name}] T≥{T_MIN}: {n_range}개")
    print(f"    t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")

    dt_vals = [zeros[i+1] - zeros[i] for i in range(len(zeros)-1)]
    print(f"    dt: min={min(dt_vals):.4f} max={max(dt_vals):.4f} mean={np.mean(dt_vals):.4f}")

    # A_bare 계산 (전체 영점에서)
    print(f"    A_bare 계산 ({n_range}개) ...")
    t0 = time.time()
    data = []
    for i in range(n_range):
        # all_zeros 내에서의 인덱스 찾기
        all_idx = all_zeros.index(zeros[i])
        S1, H1, A_bare = compute_A(all_zeros, all_idx)

        # Smooth part
        sm_im = gamma_smooth_im(zeros[i], lf)
        S1_L = S1 - sm_im
        A_L = S1_L ** 2 + 2.0 * H1

        if (math.isnan(A_bare) or math.isinf(A_bare) or A_bare <= 0
                or math.isnan(A_L) or math.isinf(A_L) or A_L <= 0):
            continue

        data.append({
            't': zeros[i],
            'S1': S1,
            'H1': H1,
            'A_bare': A_bare,
            'A_L': A_L,
            'sm_im': sm_im,
        })

        if (len(data)) % 100 == 0:
            print(f"      ...{len(data)} ({time.time()-t0:.1f}s)")

    print(f"    유효: {len(data)}/{n_range} ({time.time()-t0:.1f}s)")

    if len(data) < 20:
        print(f"    데이터 부족 — 건너뜀")
        return None

    # Gap 계산 + 내부 영점 추출
    valid = []
    for pos in range(1, len(data) - 1):
        d = data[pos]
        t_prev = data[pos - 1]['t']
        t_curr = d['t']
        t_next = data[pos + 1]['t']

        gap_r = t_next - t_curr
        gap_l = t_curr - t_prev
        if gap_r <= 0 or gap_l <= 0:
            continue

        d_bar = 2.0 / (t_next - t_prev)
        d['gap_r_gue'] = gap_r * d_bar
        d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
        valid.append(d)

    # Trim 20% (각 끝 10%)
    n_before_trim = len(valid)
    trim_n = int(n_before_trim * TRIM_FRAC)
    if trim_n > 0:
        trimmed = valid[trim_n:-trim_n]
    else:
        trimmed = valid
    n_trimmed = len(trimmed)

    print(f"    내부: {n_before_trim} → trim {2*trim_n} → {n_trimmed}")

    if n_trimmed < 15:
        print(f"    trim 후 데이터 부족 — 건너뜀")
        return None

    # Spearman 상관
    A_bare_arr = np.array([d['A_bare'] for d in trimmed])
    A_L_arr = np.array([d['A_L'] for d in trimmed])
    gm_arr = np.array([d['gap_min_gue'] for d in trimmed])

    rho_bare, p_bare = stats.spearmanr(A_bare_arr, gm_arr)
    rho_L, p_L = stats.spearmanr(A_L_arr, gm_arr)

    result = {
        'name': name,
        'q': lf['q'],
        'a': lf['a'],
        'n_all': n_all,
        'n_range': n_range,
        'n_valid': len(data),
        'n_inner': n_before_trim,
        'n_trimmed': n_trimmed,
        'rho_bare': rho_bare,
        'p_bare': p_bare,
        'rho_L': rho_L,
        'p_L': p_L,
        'delta_bare': GUE_FULL_REF - rho_bare,  # δ = |GUE| - |observed|
        'delta_L': GUE_FULL_REF - rho_L,
        'mean_A_bare': float(np.mean(A_bare_arr)),
        'mean_A_L': float(np.mean(A_L_arr)),
        'mean_gap_min': float(np.mean(gm_arr)),
        'mean_sm_im': float(np.mean([d['sm_im'] for d in trimmed])),
        't_range': (trimmed[0]['t'], trimmed[-1]['t']),
    }

    sig = lambda p: '✅' if p < 0.001 else ('⚠️' if p < 0.01 else '❌')
    print(f"\n    === {name} (q={lf['q']}) 결과 ===")
    print(f"    n_trimmed = {n_trimmed}")
    print(f"    ρ_S(A_bare, gap_min) = {rho_bare:+.4f}  p={p_bare:.3e}  {sig(p_bare)}")
    print(f"    ρ_S(A_L,    gap_min) = {rho_L:+.4f}  p={p_L:.3e}  {sig(p_L)}")
    print(f"    δ_bare = {result['delta_bare']:+.3f}  (GUE_full={GUE_FULL_REF})")
    print(f"    δ_L    = {result['delta_L']:+.3f}")

    return result


def main():
    print("=" * 70)
    print("  사이클 #288 — 디리클레 A_bare gap_min d=1 보편성 검증")
    print("=" * 70)
    print(f"  T=[{T_MIN}, {T_MAX}]  Trim={TRIM_FRAC*2:.0%}")
    print(f"  GUE_full ref = {GUE_FULL_REF}")
    print()

    results = []

    for lf in LFUNCS:
        print(f"\n{'='*60}")
        print(f"  {lf['name']} (q={lf['q']}, a={lf['a']})")
        print(f"{'='*60}")

        zeros = get_zeros(lf, T_MAX)
        if not zeros:
            continue

        r = analyze_one(lf, zeros)
        if r is not None:
            results.append(r)

    # ── 비교표 ────────────────────────────────────────────────────────
    print("\n" + "=" * 90)
    print("  d=1 보편성 비교표: A_bare gap_min Spearman ρ_S")
    print("=" * 90)
    hdr = f"{'L-함수':<20} {'q':>3} {'n':>5} {'ρ_S(A_bare)':>12} {'p':>12} {'δ_bare':>8} {'ρ_S(A_L)':>12} {'δ_L':>8} {'유의':>4}"
    print(hdr)
    print("-" * 90)

    for r in results:
        sig = '✅' if r['p_bare'] < 0.001 else ('⚠️' if r['p_bare'] < 0.01 else '❌')
        print(f"{r['name']:<20} {r['q']:>3} {r['n_trimmed']:>5} "
              f"{r['rho_bare']:>+12.4f} {r['p_bare']:>12.3e} {r['delta_bare']:>+8.3f} "
              f"{r['rho_L']:>+12.4f} {r['delta_L']:>+8.3f}  {sig}")

    # ── δ 통계 ────────────────────────────────────────────────────────
    if len(results) >= 2:
        delta_bare_all = [r['delta_bare'] for r in results]
        delta_zeta = results[0]['delta_bare'] if results[0]['q'] == 1 else None
        delta_dirichlet = [r['delta_bare'] for r in results if r['q'] > 1]

        print(f"\n  δ_bare 통계:")
        print(f"    전체 평균: {np.mean(delta_bare_all):+.3f} ± {np.std(delta_bare_all):.3f}")
        if delta_zeta is not None:
            print(f"    ζ(s): δ = {delta_zeta:+.3f}")
        if delta_dirichlet:
            print(f"    Dirichlet 평균: {np.mean(delta_dirichlet):+.3f} ± {np.std(delta_dirichlet):.3f}")
            print(f"    ζ vs Dirichlet 차이: {abs(delta_zeta - np.mean(delta_dirichlet)):.3f}")

    # Paper 4 통일표 비교
    print(f"\n  Paper 4 비교:")
    print(f"    d=1 ζ(s) [C-282, T=2000, ±300, trim20%]: ρ=-0.929, δ=0.038")
    print(f"    d=2 GL(2) [C-282c, T=500, full, trim20%]: ρ≈-0.67, δ≈0.24")
    print(f"    ⚠️ 위 수치는 T/N_MAX/trim 차이 존재. 동일 조건 비교는 이 표 내에서만 유효.")

    # ── 판정 ──────────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("  판정")
    print(f"{'='*70}")

    if len(results) >= 3:
        zeta_rho = None
        dir_rhos = []
        for r in results:
            if r['q'] == 1:
                zeta_rho = r['rho_bare']
            else:
                dir_rhos.append(r['rho_bare'])

        if zeta_rho and dir_rhos:
            max_diff = max(abs(zeta_rho - dr) for dr in dir_rhos)
            mean_dir = np.mean(dir_rhos)

            if max_diff < 0.10:
                print(f"  ★★★★★ d=1 보편성 확정: ζ({zeta_rho:+.3f}) ≈ Dirichlet({mean_dir:+.3f})")
                print(f"  최대 차이 {max_diff:.3f} < 0.10. Conductor 무관.")
                verdict = "★★★★★ 양성"
            elif max_diff < 0.15:
                print(f"  ★★★★ d=1 보편성 강 양성: ζ({zeta_rho:+.3f}) ≈ Dirichlet({mean_dir:+.3f})")
                print(f"  최대 차이 {max_diff:.3f} < 0.15. Conductor 효과 미미.")
                verdict = "★★★★ 양성"
            elif max_diff < 0.20:
                print(f"  ★★★ 조건부 양성: 차이 {max_diff:.3f}, conductor 효과 가능성")
                verdict = "★★★ 조건부"
            else:
                print(f"  ★★ 이탈: ζ({zeta_rho:+.3f}) vs Dirichlet({mean_dir:+.3f}), 차이 {max_diff:.3f}")
                print(f"  Conductor 의존성 또는 T-limited 효과.")
                verdict = "★★ 이탈"
        else:
            verdict = "데이터 부족"
    else:
        verdict = "데이터 부족"

    # ── 결과 저장 ─────────────────────────────────────────────────────
    with open(RESULT_PATH, 'w') as f:
        f.write(f"# C-288 디리클레 A_bare gap_min d=1 보편성 검증\n")
        f.write(f"# T=[{T_MIN},{T_MAX}], Trim={TRIM_FRAC*2:.0%}, GUE_full={GUE_FULL_REF}\n")
        f.write(f"# 판정: {verdict}\n\n")

        f.write("## 비교표\n")
        f.write(f"{'L-함수':<20} {'q':>3} {'n':>5} {'ρ_bare':>10} {'p':>12} {'δ_bare':>8} {'ρ_L':>10} {'δ_L':>8}\n")
        f.write("-" * 80 + "\n")
        for r in results:
            f.write(f"{r['name']:<20} {r['q']:>3} {r['n_trimmed']:>5} "
                    f"{r['rho_bare']:>+10.4f} {r['p_bare']:>12.3e} {r['delta_bare']:>+8.3f} "
                    f"{r['rho_L']:>+10.4f} {r['delta_L']:>+8.3f}\n")

        f.write(f"\n## 보조 통계\n")
        for r in results:
            f.write(f"\n### {r['name']} (q={r['q']})\n")
            f.write(f"n_all={r['n_all']}, n_range={r['n_range']}, n_valid={r['n_valid']}\n")
            f.write(f"n_inner={r['n_inner']}, n_trimmed={r['n_trimmed']}\n")
            f.write(f"t_range=[{r['t_range'][0]:.2f}, {r['t_range'][1]:.2f}]\n")
            f.write(f"mean_A_bare={r['mean_A_bare']:.4f}\n")
            f.write(f"mean_A_L={r['mean_A_L']:.4f}\n")
            f.write(f"mean_gap_min_GUE={r['mean_gap_min']:.4f}\n")
            f.write(f"mean_B_smooth={r['mean_sm_im']:.6f}\n")

        f.write(f"\n## 판정\n{verdict}\n")
        if len(results) >= 3:
            if zeta_rho and dir_rhos:
                f.write(f"ζ ρ_bare = {zeta_rho:+.4f}\n")
                f.write(f"Dirichlet mean ρ_bare = {mean_dir:+.4f} ± {np.std(dir_rhos):.4f}\n")
                f.write(f"최대 차이 = {max_diff:.4f}\n")

    print(f"\n결과 저장: {RESULT_PATH}")
    print("완료.")


if __name__ == '__main__':
    main()
