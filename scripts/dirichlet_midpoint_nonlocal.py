"""
=============================================================================
[Project RDL] 디리클레 L-함수 중간점 비국소 분석 — 결과 #18 일반화 검증
=============================================================================
사이클 27 / 수학자 지시 2026-04-14 03:18

목표:
  ζ에서 확립된 결과 #18 (κ_mid ↔ gap_next, ρ=-0.654)이
  디리클레 L-함수에서도 재현되는지 확인.

대상: χ mod 3, χ mod 4, χ mod 5

방법:
  - 각 L-함수에서 영점 수집 (t∈[0,200])
  - 중간점 m = (γₙ + γₙ₊₁)/2에서 κ_mid = |L'(s,χ)/L(s,χ)|²
    (완비 함수 Λ 대신 순수 L-함수 사용: underflow 없음)
  - 상관 분석 (a)(b)(c) + 편상관
  - ζ 결과와 비교표

수학적 근거 (중간점 Hadamard 상쇄):
  Dirichlet L-함수의 Hadamard 전개에서 s = 1/2 + it_mid 일 때,
  중간점 t_mid = (γₙ+γₙ₊₁)/2에서 최근접 두 영점의 기여:
    1/(s-ρₙ) + 1/(s-ρₙ₊₁)
  = -i/(t_mid-γₙ) + (-i/(t_mid-γₙ₊₁))
  = -i/(gap/2) + i/(gap/2) = 0 (완전 상쇄)
  → κ_mid = 비국소 기여만 반영. ζ와 동일한 메커니즘.

κ 계산 인프라:
  L(s,χ) = mpmath.dirichlet(s, chi) → 절대값 O(1), underflow 없음
  L'(s,χ)/L(s,χ): 수치미분 h=10^{-50}, dps=100
  (ζ v2b의 해석적 공식과 달리, L 자체가 안전하므로 수치미분 직접 사용)

성공 기준:
  - cap < 5% (각 L-함수) → 필수
  - (a) |ρ(κ_mid, gap)| < 0.3 → 필수
  - (b) 3개 중 2개 이상 |ρ(κ_mid, gap_next)| > 0.3, p < 0.05 → 양성
  - 3개 모두 |ρ| < 0.3 → 음성

결과 파일: results/dirichlet_midpoint_nonlocal.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats
from datetime import datetime

# ─── 경로 설정 ───────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT     = os.path.join(BASE_DIR, "results", "dirichlet_midpoint_nonlocal.txt")

os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# ─── 정밀도 설정 ─────────────────────────────────────────────────────────────
# 영점 탐색: dps=50 (속도), κ 계산: dps=100 (정밀도)
DPS_SCAN  = 50
DPS_KAPPA = 100
CAP_VALUE = 1e18
H_DIFF = mpmath.mpf(10) ** (-50)   # dps=100 → 유효 자릿수 충분
mpmath.mp.dps = DPS_SCAN  # 기본값: 탐색용

# ─── 로그 ────────────────────────────────────────────────────────────────────
log_lines = []
def log(msg=""):
    print(msg, flush=True)
    log_lines.append(msg)

def save_results():
    with open(OUTPUT, 'w') as f:
        f.write("\n".join(log_lines))


# ─── 지표 정의 ───────────────────────────────────────────────────────────────
# Dirichlet character `a` 파라미터: χ(-1)=+1→a=0, χ(-1)=-1→a=1

CHARACTERS = {
    'χ₃ (mod 3)': {
        'chi': [0, 1, -1],          # χ(0)=0, χ(1)=1, χ(2)=-1
        'q': 3,
        'a': 1,                      # χ(-1)=χ(2)=-1
    },
    'χ₄ (mod 4)': {
        'chi': [0, 1, 0, -1],       # χ(0)=0, χ(1)=1, χ(2)=0, χ(3)=-1
        'q': 4,
        'a': 1,                      # χ(-1)=χ(3)=-1
    },
    'χ₅ (mod 5)': {
        # 실수 지표 mod 5: χ(1)=1, χ(2)=-1, χ(3)=-1, χ(4)=1
        # (χ₅²형, quadratic character)
        'chi': [0, 1, -1, -1, 1],   # χ(0)=0, χ(1)=1, χ(2)=-1, χ(3)=-1, χ(4)=1
        'q': 5,
        'a': 0,                      # χ(-1)=χ(4)=+1 → a=0
    },
}


# ─── κ 계산 함수 ─────────────────────────────────────────────────────────────

def L_val(s, chi):
    """순수 디리클레 L-함수 L(s,χ)"""
    return mpmath.dirichlet(s, chi)


def logderiv_L(s, chi):
    """
    L'(s,χ)/L(s,χ) 계산.
    L(s,χ) = mpmath.dirichlet 자체가 O(1)이므로 수치미분 안전.
    h = 10^{-50}, dps=100 → 유효자릿수 50자리 확보.
    """
    L_center = L_val(s, chi)

    # 영점 근방 판정
    if abs(L_center) < mpmath.mpf(10) ** (-mpmath.mp.dps + 15):
        raise ValueError(f"L(s,χ) near zero at s={s}")

    L_plus  = L_val(s + H_DIFF, chi)
    L_minus = L_val(s - H_DIFF, chi)
    L_prime = (L_plus - L_minus) / (2 * H_DIFF)

    return L_prime / L_center


def kappa_mid(t_mid_float, chi):
    """
    중간점 t_mid에서 κ = |L'/L(s, χ)|²
    s = 1/2 + i*t_mid
    반환: float 또는 None (cap/실패)
    """
    try:
        s = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_mid_float)))
        ld = logderiv_L(s, chi)
        kappa = float(abs(ld) ** 2)
        if not (kappa == kappa) or kappa > CAP_VALUE:  # NaN or cap
            return None
        return kappa
    except Exception as e:
        print(f"  WARNING: t_mid={t_mid_float:.4f}: 실패 → {e}", flush=True)
        return None


# ─── 완비 L-함수 (영점 탐색용) ────────────────────────────────────────────────

def completed_L_func(s, char_info):
    """
    완비 L-함수: Λ(s, χ) = (q/π)^{s/2} Γ((s+a)/2) L(s, χ)
    Λ는 임계선에서 실수값에 가까운 Z-함수 거동 → Re(Λ) 부호변화로 영점 탐색.
    """
    q   = mpmath.mpf(char_info['q'])
    a   = mpmath.mpf(char_info['a'])
    chi = char_info['chi']
    L   = mpmath.dirichlet(s, chi)
    G   = mpmath.gamma((s + a) / 2)
    pre = mpmath.power(q / mpmath.pi, s / 2)
    return pre * G * L


# ─── 영점 탐색 ───────────────────────────────────────────────────────────────

def find_zeros_L(char_info, t_min=1.0, t_max=200.0, n_scan=2000):
    """
    임계선 σ=1/2에서 Re(Λ(s,χ)) 부호 변화 + findroot로 영점 정밀화.
    (dirichlet_bundle_verification.py 방식과 동일)
    """
    chi = char_info['chi']
    ts  = np.linspace(t_min, t_max, n_scan)
    zeros = []
    fail_count = 0

    prev_re = None
    prev_t  = None

    mpmath.mp.dps = DPS_SCAN  # 탐색 시 낮은 정밀도로 속도 확보

    for t in ts:
        s = mpmath.mpc('0.5', str(t))
        try:
            val    = completed_L_func(s, char_info)
            curr_re = float(mpmath.re(val))
        except Exception as e:
            print(f"  WARNING: completed_L 실패 t={t:.3f}: {e}", flush=True)
            prev_re = None
            prev_t  = t
            continue

        if prev_re is not None and prev_re * curr_re < 0:
            # 부호 변화 → findroot (bisection 방식, tol=1e-8로 underflow 허용)
            try:
                def f_re(t_var, ci=char_info):
                    sv = mpmath.mpc('0.5', str(float(t_var)))
                    return mpmath.re(completed_L_func(sv, ci))

                # Λ는 대형 t에서 exp(-πt/4)으로 소멸 → tol을 적응적으로 설정
                func_scale = max(abs(prev_re), abs(curr_re))
                local_tol  = max(func_scale * mpmath.mpf('1e-10'),
                                 mpmath.mpf(10) ** (-DPS_SCAN + 5))
                t_zero = mpmath.findroot(f_re, (prev_t, t), tol=local_tol)
                t_zero_f = float(t_zero)

                # Λ는 임계선에서 실수→ 부호 변화 = 영점 (대형 t에서도 타당)
                if not zeros or abs(t_zero_f - zeros[-1]) > 0.05:
                    zeros.append(t_zero_f)
            except Exception as e:
                fail_count += 1

        prev_re = curr_re
        prev_t  = t

    zeros = sorted(set(zeros))

    # 중복 제거 (0.1 이내)
    filtered = []
    for z in zeros:
        if not filtered or abs(z - filtered[-1]) > 0.1:
            filtered.append(z)

    if fail_count > len(filtered) // 2 and fail_count > 5:
        print(f"  ⚠️ findroot 실패 {fail_count}회 — 탐색 로직 점검 필요", flush=True)

    if len(filtered) == 0:
        print("⚠️ 영점 0개 — 탐색 로직 점검 필요", flush=True)

    # κ 계산 전 고정밀도로 복원
    mpmath.mp.dps = DPS_KAPPA

    return np.array(filtered)


# ─── 상관 분석 ───────────────────────────────────────────────────────────────

def partial_correlation(x, y, z):
    """
    편상관 ρ(x, y | z) — z를 통제한 후 x, y의 Spearman 상관.
    """
    rho_xy, _ = stats.spearmanr(x, y)
    rho_xz, _ = stats.spearmanr(x, z)
    rho_yz, _ = stats.spearmanr(y, z)

    denom = np.sqrt((1 - rho_xz**2) * (1 - rho_yz**2))
    if denom < 1e-10:
        return float('nan')
    return (rho_xy - rho_xz * rho_yz) / denom


def run_correlation_analysis(kappa_arr, gap_arr, gap_next_arr, density_arr, char_name):
    """상관 분석 (a)(b)(c) + 편상관"""
    n = len(kappa_arr)

    # (a) ρ(κ_mid, gap)
    rho_a, p_a = stats.spearmanr(kappa_arr, gap_arr)

    # (b) ρ(κ_mid, gap_next)
    rho_b, p_b = stats.spearmanr(kappa_arr, gap_next_arr)

    # (c) ρ(κ_mid, density)
    rho_c, p_c = stats.spearmanr(kappa_arr, density_arr)

    # gap 자기상관
    rho_gg, p_gg = stats.spearmanr(gap_arr, gap_next_arr)

    # 편상관 ρ(κ_mid, gap_next | gap)
    pc = partial_correlation(kappa_arr, gap_next_arr, gap_arr)

    # 분위 비교
    q25 = np.percentile(gap_next_arr, 25)
    q75 = np.percentile(gap_next_arr, 75)
    kappa_low  = kappa_arr[gap_next_arr <= q25].mean()
    kappa_high = kappa_arr[gap_next_arr >= q75].mean()

    return {
        'char': char_name,
        'n': n,
        'rho_a': rho_a, 'p_a': p_a,
        'rho_b': rho_b, 'p_b': p_b,
        'rho_c': rho_c, 'p_c': p_c,
        'rho_gg': rho_gg, 'p_gg': p_gg,
        'partial_corr': pc,
        'kappa_low25': kappa_low,
        'kappa_high75': kappa_high,
        'q25': q25, 'q75': q75,
    }


# ─── 메인 실험 ───────────────────────────────────────────────────────────────

log("=" * 70)
log("디리클레 L-함수 중간점 비국소 분석 — 결과 #18 일반화 검증")
log(f"시작 시각: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log(f"dps 탐색={DPS_SCAN}, κ계산={DPS_KAPPA}, h=1e-50, CAP={CAP_VALUE:.0e}")
log("인프라: L'/L 수치미분 (L(s,χ) = mpmath.dirichlet, underflow 없음)")
log("=" * 70)
log()

all_results = []

for char_name, char_info in CHARACTERS.items():
    t0 = time.time()
    chi = char_info['chi']
    log(f"{'='*60}")
    log(f"[{char_name}] q={char_info['q']}, a={char_info['a']}")
    log(f"{'='*60}")

    # ── 1단계: 영점 수집 ──────────────────────────────────────────────────
    log(f"\n1단계: 영점 탐색 (t∈[1,200], n_scan=2000, dps={DPS_SCAN}) ...")
    zeros = find_zeros_L(char_info, t_min=1.0, t_max=200.0, n_scan=2000)
    N = len(zeros)
    log(f"  영점 수: {N}")
    if N == 0:
        log(f"  ⚠️ 영점 없음 — 스킵")
        continue
    log(f"  t 범위: {zeros[0]:.3f} ~ {zeros[-1]:.3f}")
    log(f"  연속 쌍 수: {N-1}")

    # ── 2단계: 중간점 κ 계산 ─────────────────────────────────────────────
    log(f"\n2단계: 중간점 κ 계산 ...")
    pairs = N - 1
    recs = []
    cap_count = 0

    for i in range(pairs):
        g_n   = zeros[i]
        g_np1 = zeros[i + 1]
        m     = (g_n + g_np1) / 2.0
        gap   = g_np1 - g_n

        # gap_next (마지막 쌍 제외)
        if i < pairs - 1:
            gap_next = zeros[i + 2] - zeros[i + 1]
        else:
            gap_next = None  # 마지막 쌍

        kappa = kappa_mid(m, chi)

        if kappa is None:
            cap_count += 1
            continue
        if gap_next is None:
            continue  # 마지막 쌍 제외 (gap_next 없음)

        # 국소 밀도: 인근 ±5 영점의 밀도
        window = 5
        i_lo = max(0, i - window)
        i_hi = min(N - 1, i + 1 + window)
        span = zeros[i_hi] - zeros[i_lo]
        density = (i_hi - i_lo) / span if span > 0 else 0.0

        recs.append({
            'i': i, 'm': m, 'gap': gap,
            'gap_next': gap_next, 'density': density,
            'kappa': kappa,
        })

        if (i + 1) % 20 == 0:
            log(f"  [{i+1}/{pairs}] m={m:.2f}, gap={gap:.4f}, κ={kappa:.4f}, gap_next={gap_next:.4f}")

    n_valid = len(recs)
    cap_pct = cap_count / pairs * 100 if pairs > 0 else 0

    log(f"\n  유효 데이터: {n_valid}개 / {pairs}쌍")
    log(f"  cap 발생: {cap_count}개 ({cap_pct:.1f}%) [목표: <5%]")

    if cap_pct > 10:
        log(f"  ⚠️ cap {cap_pct:.1f}% > 10% — 인프라 점검 필요")
    elif cap_pct > 5:
        log(f"  ⚠️ cap {cap_pct:.1f}% > 5% — 주의")
    else:
        log(f"  ✅ cap {cap_pct:.1f}% < 5%")

    if n_valid < 30:
        log(f"  ⚠️ 유효 데이터 부족 ({n_valid} < 30) — 분석 신뢰도 낮음")
        continue

    # 기술 통계
    kappa_arr   = np.array([r['kappa']    for r in recs])
    gap_arr     = np.array([r['gap']      for r in recs])
    gap_next_arr= np.array([r['gap_next'] for r in recs])
    density_arr = np.array([r['density']  for r in recs])

    log(f"\n  κ 기술통계: min={kappa_arr.min():.4f}, max={kappa_arr.max():.4f}, "
        f"mean={kappa_arr.mean():.4f}, median={np.median(kappa_arr):.4f}")
    log(f"  gap 기술통계: min={gap_arr.min():.4f}, max={gap_arr.max():.4f}, "
        f"mean={gap_arr.mean():.4f}")

    # ── 3단계: 상관 분석 ─────────────────────────────────────────────────
    log(f"\n3단계: 상관 분석 (n={n_valid}) ...")
    res = run_correlation_analysis(kappa_arr, gap_arr, gap_next_arr, density_arr, char_name)
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
    ratio = res['kappa_low25'] / res['kappa_high75'] if res['kappa_high75'] > 0 else float('nan')
    log(f"    비율: {ratio:.2f}×")

    log(f"\n  소요: {res['elapsed']:.1f}초")
    log()

    # 중간 저장
    save_results()


# ─── 최종 비교표 ─────────────────────────────────────────────────────────────
log()
log("=" * 70)
log("최종 비교표 — ζ vs 디리클레 L-함수 중간점 비국소 분석")
log("=" * 70)
log()

# ζ 기준값 (결과 #18)
log(f"{'L-함수':<20} {'영점수':>6} {'cap%':>6} {'ρ(a)':>8} {'ρ(b)':>8} {'p(b)':>12} {'편상관':>8} {'판정':>8}")
log("-" * 85)

# ζ 기준
log(f"{'ζ (결과 #18)':<20} {'234':>6} {'0.0':>6} {'-0.031':>8} {'-0.654':>8} {'7.3e-30':>12} {'-0.666':>8} {'✅ 양성':>8}")

for res in all_results:
    rho_b = res['rho_b']
    p_b   = res['p_b']
    pc    = res['partial_corr']

    # 판정
    cap_ok  = res['cap_pct'] < 5
    nn_ok   = abs(res['rho_a']) < 0.3
    nl_ok   = abs(rho_b) > 0.3 and p_b < 0.05

    if cap_ok and nn_ok and nl_ok:
        verdict = "✅ 양성"
    elif cap_ok and nn_ok:
        verdict = "❌ 음성"
    elif not cap_ok:
        verdict = "⚠️ cap"
    else:
        verdict = "⚠️ 혼재"

    log(f"{res['char']:<20} {res['n_zeros']:>6} {res['cap_pct']:>6.1f} "
        f"{res['rho_a']:>+8.3f} {rho_b:>+8.3f} {p_b:>12.2e} "
        f"{pc:>+8.3f} {verdict:>8}")

log()

# 전체 판정
if len(all_results) >= 2:
    positive = sum(1 for r in all_results
                   if r['cap_pct'] < 5
                   and abs(r['rho_a']) < 0.3
                   and abs(r['rho_b']) > 0.3
                   and r['p_b'] < 0.05)
    total = len(all_results)

    log(f"성공 기준 달성: {positive}/{total} 개 L-함수에서 양성")
    log()
    if positive >= 2:
        log("★ 최종 판정: 양성 (2개 이상 L-함수에서 비국소 현상 재현)")
        log("  → 결과 #18은 ζ 고유가 아닌 L-함수 일반에서 보편적 현상")
    elif positive == 0:
        log("★ 최종 판정: 음성 (3개 L-함수 모두 |ρ| < 0.3)")
        log("  → 결과 #18은 ζ 고유 현상. 디리클레에서 재현 불가.")
    else:
        log(f"★ 최종 판정: 중립 ({positive}/{total} 양성)")
        log("  → 추가 L-함수 실험 필요")

log()
log(f"결과 파일: {OUTPUT}")
log(f"종료 시각: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 70)

save_results()
print(f"\n✅ 완료. 결과 저장: {OUTPUT}", flush=True)
