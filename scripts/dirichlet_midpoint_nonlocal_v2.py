"""
=============================================================================
[Project RDL] 디리클레 L-함수 중간점 비국소 분석 v2 — 해석적 Λ'/Λ 사용
=============================================================================
사이클 27 / 수학자 지시 2026-04-14 03:18

v1 문제:
  bare L'/L(s,χ) 사용 시 digamma ψ((s+a)/2) ~ log(t) 항이 포함
  → κ ∝ t(monotone) → ρ(κ, density)=0.97 (t-추세 혼입)
  → ρ(κ, gap)=-0.43 (nn 상쇄 기준 실패)
  즉, bare L'/L = Λ'/Λ - (1/2)log(q/π) - (1/2)ψ((s+a)/2) 이므로
  ψ 항이 t와 함께 성장하여 κ_bare가 density를 추적함.

v2 수정:
  Λ'/Λ를 해석적 공식으로 계산:
    Λ'/Λ(s,χ) = (1/2)log(q/π) + (1/2)ψ((s+a)/2) + L'(s,χ)/L(s,χ)
  여기서:
    - (1/2)log(q/π): 실수 상수, 안전
    - (1/2)ψ((s+a)/2): mpmath.digamma, O(log t), 안전
    - L'(s,χ)/L(s,χ): mpmath.dirichlet 수치미분, L은 O(1) (underflow 없음)
  → Λ 자체를 계산하지 않으므로 exp(-πt/4) underflow 완전 회피
  → ζ v2b와 동일한 철학

Hadamard 상쇄 (중간점):
  Λ'/Λ(s,χ) = B + Σ_ρ [1/(s-ρ) + ...]
  At midpoint m = (γₙ+γₙ₊₁)/2:
    1/(m-ρₙ) + 1/(m-ρₙ₊₁) = -2i/gap + 2i/gap = 0 (완전 상쇄)
  → κ_Λ = |비국소 기여|² 만 반영. Bare L'/L보다 훨씬 작고 GRH 구조를 반영.

성공 기준 (ζ v2b와 동일):
  - cap < 5%
  - (a) |ρ(κ_mid, gap)| < 0.3 (nn 상쇄 확인)
  - (b) 3개 중 2개 이상: |ρ(κ_mid, gap_next)| > 0.3, p < 0.05 → 양성

결과 파일: results/dirichlet_midpoint_nonlocal.txt (덮어쓰기)
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT     = os.path.join(BASE_DIR, "results", "dirichlet_midpoint_nonlocal.txt")
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

DPS_SCAN  = 50    # 영점 탐색용 (속도)
DPS_KAPPA = 100   # κ 계산용 (정밀도)
CAP_VALUE = 1e18
H_DIFF    = mpmath.mpf(10) ** (-50)  # dps=100 → 유효 50자리

log_lines = []
def log(msg=""):
    print(msg, flush=True)
    log_lines.append(msg)

def save_results():
    with open(OUTPUT, 'w') as f:
        f.write("\n".join(log_lines))


# ─── 지표 정의 ───────────────────────────────────────────────────────────────
CHARACTERS = {
    'χ₃ (mod 3)': {
        'chi': [0, 1, -1],
        'q': 3,
        'a': 1,    # χ₃(-1) = χ₃(2) = -1
    },
    'χ₄ (mod 4)': {
        'chi': [0, 1, 0, -1],
        'q': 4,
        'a': 1,    # χ₄(-1) = χ₄(3) = -1
    },
    'χ₅ (mod 5)': {
        'chi': [0, 1, -1, -1, 1],   # quadratic char mod 5
        'q': 5,
        'a': 0,    # χ₅(-1) = χ₅(4) = +1
    },
}


# ─── 완비 L-함수 (영점 탐색용) ────────────────────────────────────────────────

def completed_L_func(s, char_info):
    """Λ(s,χ) = (q/π)^{s/2} Γ((s+a)/2) L(s,χ) — 탐색 전용, underflow 있음"""
    q   = mpmath.mpf(char_info['q'])
    a   = mpmath.mpf(char_info['a'])
    chi = char_info['chi']
    L   = mpmath.dirichlet(s, chi)
    G   = mpmath.gamma((s + a) / 2)
    pre = mpmath.power(q / mpmath.pi, s / 2)
    return pre * G * L


# ─── 해석적 Λ'/Λ (κ 계산용) ──────────────────────────────────────────────────

def connection_Lambda_analytic(s, char_info):
    """
    Λ'/Λ(s,χ) = (1/2)log(q/π) + (1/2)ψ((s+a)/2) + L'(s,χ)/L(s,χ)

    ζ v2b의 ξ'/ξ 해석적 공식과 동일한 철학:
    - Λ 자체를 계산하지 않음 → underflow 완전 회피
    - L(s,χ) = O(1) 이므로 수치미분 dps=100에서 안전

    반환: mpmath complex (연결 형식 값)
    예외: ValueError (L 영점 근방)
    """
    mpmath.mp.dps = DPS_KAPPA

    chi = char_info['chi']
    q   = mpmath.mpf(char_info['q'])
    a   = mpmath.mpf(char_info['a'])

    # L(s,χ) 값 (중간점에서 ≠ 0이어야 함)
    L_center = mpmath.dirichlet(s, chi)
    if abs(L_center) < mpmath.mpf(10) ** (-DPS_KAPPA + 15):
        raise ValueError(f"L(s,χ) near zero at s={s}")

    # L'(s,χ)/L(s,χ) 수치미분
    L_plus  = mpmath.dirichlet(s + H_DIFF, chi)
    L_minus = mpmath.dirichlet(s - H_DIFF, chi)
    L_prime = (L_plus - L_minus) / (2 * H_DIFF)
    log_deriv_L = L_prime / L_center

    # 해석적 항들
    term_log = mpmath.log(q / mpmath.pi) / 2        # (1/2)log(q/π)
    term_psi = mpmath.digamma((s + a) / 2) / 2      # (1/2)ψ((s+a)/2)

    return term_log + term_psi + log_deriv_L


def kappa_mid(t_mid_float, char_info):
    """
    중간점 t_mid에서 κ = |Λ'/Λ(s,χ)|²
    반환: float 또는 None (cap/실패)
    """
    try:
        s = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_mid_float)))
        conn = connection_Lambda_analytic(s, char_info)
        kappa = float(abs(conn) ** 2)
        if not (kappa == kappa) or kappa > CAP_VALUE:
            return None
        return kappa
    except Exception as e:
        print(f"  WARNING: t_mid={t_mid_float:.4f}: 실패 → {e}", flush=True)
        return None


# ─── 영점 탐색 ───────────────────────────────────────────────────────────────

def find_zeros_L(char_info, t_min=1.0, t_max=200.0, n_scan=2000):
    """
    완비 Λ(s,χ)의 Re 부호 변화로 임계선 영점 탐색.
    Λ는 임계선에서 실수값 (Im≈0) → Re 부호 변화 = 영점.
    underflow 허용: local_tol = 함수 스케일 × 1e-10
    """
    chi = char_info['chi']
    ts  = np.linspace(t_min, t_max, n_scan)
    zeros = []
    fail_count = 0

    prev_re = None
    prev_t  = None

    mpmath.mp.dps = DPS_SCAN

    for t in ts:
        s = mpmath.mpc('0.5', str(t))
        try:
            val     = completed_L_func(s, char_info)
            curr_re = float(mpmath.re(val))
        except Exception as e:
            prev_re = None; prev_t = t
            continue

        if prev_re is not None and prev_re * curr_re < 0:
            try:
                def f_re(t_var, ci=char_info):
                    sv = mpmath.mpc('0.5', str(float(t_var)))
                    return mpmath.re(completed_L_func(sv, ci))

                # 적응형 tol: Λ가 대형 t에서 underflow하므로 함수 스케일 기반
                func_scale = max(abs(prev_re), abs(curr_re))
                local_tol  = max(func_scale * mpmath.mpf('1e-10'),
                                 mpmath.mpf(10) ** (-DPS_SCAN + 5))

                t_zero     = mpmath.findroot(f_re, (prev_t, t), tol=local_tol)
                t_zero_f   = float(t_zero)

                if not zeros or abs(t_zero_f - zeros[-1]) > 0.05:
                    zeros.append(t_zero_f)
            except Exception as e:
                fail_count += 1

        prev_re = curr_re
        prev_t  = t

    zeros = sorted(set(zeros))

    # 중복 제거
    filtered = []
    for z in zeros:
        if not filtered or abs(z - filtered[-1]) > 0.08:
            filtered.append(z)

    if fail_count > len(filtered) // 2 and fail_count > 5:
        print(f"  ⚠️ findroot 실패 {fail_count}회 (영점 {len(filtered)}개) — 탐색 로직 점검", flush=True)

    if len(filtered) == 0:
        print("⚠️ 영점 0개 — 탐색 로직 점검 필요", flush=True)

    mpmath.mp.dps = DPS_KAPPA  # κ 계산 전 고정밀도 복원
    return np.array(filtered)


# ─── 상관 분석 ───────────────────────────────────────────────────────────────

def partial_correlation(x, y, z):
    """편상관 ρ(x, y | z) — Spearman 기반"""
    rho_xy, _ = stats.spearmanr(x, y)
    rho_xz, _ = stats.spearmanr(x, z)
    rho_yz, _ = stats.spearmanr(y, z)
    denom = np.sqrt((1 - rho_xz**2) * (1 - rho_yz**2))
    if denom < 1e-10:
        return float('nan')
    return (rho_xy - rho_xz * rho_yz) / denom


def analyze(kappa_arr, gap_arr, gap_next_arr, density_arr, char_name):
    n = len(kappa_arr)
    rho_a, p_a = stats.spearmanr(kappa_arr, gap_arr)
    rho_b, p_b = stats.spearmanr(kappa_arr, gap_next_arr)
    rho_c, p_c = stats.spearmanr(kappa_arr, density_arr)
    rho_gg, p_gg = stats.spearmanr(gap_arr, gap_next_arr)
    pc = partial_correlation(kappa_arr, gap_next_arr, gap_arr)

    q25 = np.percentile(gap_next_arr, 25)
    q75 = np.percentile(gap_next_arr, 75)
    kl  = kappa_arr[gap_next_arr <= q25].mean()
    kh  = kappa_arr[gap_next_arr >= q75].mean()

    return {'char': char_name, 'n': n,
            'rho_a': rho_a, 'p_a': p_a,
            'rho_b': rho_b, 'p_b': p_b,
            'rho_c': rho_c, 'p_c': p_c,
            'rho_gg': rho_gg, 'p_gg': p_gg,
            'partial_corr': pc,
            'kappa_low25': kl, 'kappa_high75': kh,
            'q25': q25, 'q75': q75}


# ─── 메인 ────────────────────────────────────────────────────────────────────

log("=" * 70)
log("디리클레 L-함수 중간점 비국소 분석 v2 (해석적 Λ'/Λ)")
log(f"시작 시각: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log(f"dps 탐색={DPS_SCAN}, κ계산={DPS_KAPPA}, h=1e-50, CAP={CAP_VALUE:.0e}")
log("인프라: Λ'/Λ = (1/2)log(q/π) + (1/2)ψ((s+a)/2) + L'(s,χ)/L(s,χ)")
log("  → Λ 자체 계산 불필요, underflow 완전 회피 (ζ v2b와 동일 철학)")
log("=" * 70)
log()

all_results = []

for char_name, char_info in CHARACTERS.items():
    t0  = time.time()
    log(f"{'='*60}")
    log(f"[{char_name}] q={char_info['q']}, a={char_info['a']}")
    log(f"{'='*60}")

    # 1단계: 영점 탐색
    log(f"\n1단계: 영점 탐색 (t∈[1,200], n_scan=2000, dps={DPS_SCAN}) ...")
    zeros = find_zeros_L(char_info, t_min=1.0, t_max=200.0, n_scan=2000)
    N = len(zeros)
    log(f"  영점 수: {N}")
    if N == 0:
        log(f"  ⚠️ 영점 없음 — 스킵")
        continue
    log(f"  t 범위: {zeros[0]:.3f} ~ {zeros[-1]:.3f}")
    log(f"  연속 쌍 수: {N-1}")

    # 2단계: 중간점 κ 계산
    log(f"\n2단계: 중간점 κ = |Λ'/Λ(s,χ)|² 계산 (dps={DPS_KAPPA}) ...")
    pairs     = N - 1
    recs      = []
    cap_count = 0

    for i in range(pairs):
        g_n   = zeros[i]
        g_np1 = zeros[i + 1]
        m     = (g_n + g_np1) / 2.0
        gap   = g_np1 - g_n

        if i < pairs - 1:
            gap_next = zeros[i + 2] - zeros[i + 1]
        else:
            gap_next = None

        kappa = kappa_mid(m, char_info)

        if kappa is None:
            cap_count += 1
            continue
        if gap_next is None:
            continue  # 마지막 쌍

        # 국소 밀도
        window = 5
        i_lo = max(0, i - window)
        i_hi = min(N - 1, i + 1 + window)
        span = zeros[i_hi] - zeros[i_lo]
        density = (i_hi - i_lo) / span if span > 0 else 0.0

        recs.append({'i': i, 'm': m, 'gap': gap,
                     'gap_next': gap_next, 'density': density,
                     'kappa': kappa})

        if (i + 1) % 20 == 0:
            log(f"  [{i+1}/{pairs}] m={m:.2f}, gap={gap:.4f}, "
                f"κ={kappa:.4f}, gap_next={gap_next:.4f}")

    n_valid  = len(recs)
    cap_pct  = cap_count / pairs * 100 if pairs > 0 else 0

    log(f"\n  유효 데이터: {n_valid}개 / {pairs}쌍")
    log(f"  cap 발생: {cap_count}개 ({cap_pct:.1f}%) [목표: <5%]")
    if cap_pct > 10:
        log(f"  ⚠️ cap {cap_pct:.1f}% > 10% — 인프라 점검 필요")
    elif cap_pct > 5:
        log(f"  ⚠️ cap {cap_pct:.1f}% > 5%")
    else:
        log(f"  ✅ cap {cap_pct:.1f}% < 5%")

    if n_valid < 30:
        log(f"  ⚠️ 유효 데이터 부족 ({n_valid}) — 스킵")
        continue

    kappa_arr    = np.array([r['kappa']    for r in recs])
    gap_arr      = np.array([r['gap']      for r in recs])
    gap_next_arr = np.array([r['gap_next'] for r in recs])
    density_arr  = np.array([r['density']  for r in recs])

    log(f"\n  κ: min={kappa_arr.min():.4f}, max={kappa_arr.max():.4f}, "
        f"mean={kappa_arr.mean():.4f}, median={np.median(kappa_arr):.4f}")
    log(f"  gap: min={gap_arr.min():.4f}, max={gap_arr.max():.4f}, "
        f"mean={gap_arr.mean():.4f}")

    # 3단계: 상관 분석
    log(f"\n3단계: 상관 분석 (n={n_valid}) ...")
    res = analyze(kappa_arr, gap_arr, gap_next_arr, density_arr, char_name)
    res['cap_pct'] = cap_pct
    res['n_zeros'] = N
    res['elapsed'] = time.time() - t0
    all_results.append(res)

    log(f"  (a) ρ(κ_mid, gap)      = {res['rho_a']:+.4f}  p={res['p_a']:.3e}")
    log(f"  (b) ρ(κ_mid, gap_next) = {res['rho_b']:+.4f}  p={res['p_b']:.3e}  ← 핵심")
    log(f"  (c) ρ(κ_mid, density)  = {res['rho_c']:+.4f}  p={res['p_c']:.3e}")
    log(f"  gap 자기상관: ρ(gap, gap_next) = {res['rho_gg']:+.4f}  p={res['p_gg']:.3e}")
    log(f"  편상관 ρ(κ_mid, gap_next | gap) = {res['partial_corr']:+.4f}")
    log(f"\n  분위 비교 (gap_next 기준):")
    log(f"    하위 25% (≤{res['q25']:.3f}): κ 평균 = {res['kappa_low25']:.4f}")
    log(f"    상위 25% (≥{res['q75']:.3f}): κ 평균 = {res['kappa_high75']:.4f}")
    ratio = (res['kappa_low25'] / res['kappa_high75']
             if res['kappa_high75'] > 0 else float('nan'))
    log(f"    비율: {ratio:.2f}×")
    log(f"\n  소요: {res['elapsed']:.1f}초")
    log()

    save_results()


# ─── 최종 비교표 ─────────────────────────────────────────────────────────────
log()
log("=" * 70)
log("최종 비교표 — ζ vs 디리클레 L-함수 중간점 비국소 분석")
log("=" * 70)
log()
log(f"{'L-함수':<20} {'영점수':>6} {'cap%':>6} {'ρ(a)':>8} {'ρ(b)':>8} "
    f"{'p(b)':>12} {'편상관':>8} {'판정':>8}")
log("-" * 85)

# ζ 기준 (결과 #18)
log(f"{'ζ (결과 #18)':<20} {'234':>6} {'0.0':>6} {'-0.031':>8} {'-0.654':>8} "
    f"{'7.3e-30':>12} {'-0.666':>8} {'✅ 양성':>8}")

positive = 0
for res in all_results:
    cap_ok = res['cap_pct'] < 5
    nn_ok  = abs(res['rho_a']) < 0.3
    nl_ok  = abs(res['rho_b']) > 0.3 and res['p_b'] < 0.05

    if cap_ok and nn_ok and nl_ok:
        verdict = "✅ 양성"
        positive += 1
    elif cap_ok and nn_ok:
        verdict = "❌ 음성"
    elif not cap_ok:
        verdict = "⚠️ cap"
    else:
        verdict = "⚠️ 혼재"

    log(f"{res['char']:<20} {res['n_zeros']:>6} {res['cap_pct']:>6.1f} "
        f"{res['rho_a']:>+8.3f} {res['rho_b']:>+8.3f} {res['p_b']:>12.2e} "
        f"{res['partial_corr']:>+8.3f} {verdict:>8}")

log()
total = len(all_results)
log(f"성공 기준 달성: {positive}/{total} 개 L-함수에서 양성")
log()

if positive >= 2:
    log("★ 최종 판정: 양성 (2개 이상 L-함수에서 비국소 현상 재현)")
    log("  → 결과 #18은 ζ 고유가 아닌 L-함수 일반에서 보편적 현상")
elif positive == 0:
    log("★ 최종 판정: 음성 (3개 L-함수 모두 (a)·(b) 기준 미달)")
    log("  → 결과 #18은 ζ 고유 현상일 수 있음. 수학자 판단 필요.")
elif positive == 1:
    log(f"★ 최종 판정: 중립 ({positive}/{total} 양성)")
    log("  → 추가 L-함수 실험 필요")

log()
log(f"결과 파일: {OUTPUT}")
log(f"종료 시각: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 70)

save_results()
print(f"\n✅ 완료. 결과 저장: {OUTPUT}", flush=True)
