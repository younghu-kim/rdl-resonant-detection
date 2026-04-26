"""
=============================================================================
[Project RDL] C-344 — B-68 E비 결정인자 분석
=============================================================================
질문: E비(에너지 집중도 E(0.5)/E(0.3))를 결정하는 해석적 양은 무엇인가?

후보 변수:
  1. γ₁ (첫 영점 높이)
  2. |L(1,χ)| (임계값)
  3. |τ(χ)|/√q (정규화 Gauss sum — 원시이면 항상 1)
  4. arg(τ(χ)) (Gauss sum 위상 = root number와 관련)
  5. order(χ)
  6. parity a(χ)
  7. log(E(0.5)) 자체 (E(0.3)이 거의 불변이면 E(0.5)가 결정인자)

데이터: q=3,4,5,7,8,11의 모든 원시 비자명 지표 (20개).
E비는 기존 결과 파일에서 수동 입력 (재계산 불필요).
해석적 양은 mpmath로 직접 계산.

결과: results/b68_eratio_determinant_c344.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

from bundle_utils import completed_L

mpmath.mp.dps = 80

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 결과 파일 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
RESULT_DIR = os.path.expanduser('~/Desktop/gdl_unified/results')
RESULT_FILE = os.path.join(RESULT_DIR, 'b68_eratio_determinant_c344.txt')

out = open(RESULT_FILE, 'w')

def log(msg=''):
    print(msg, flush=True)
    out.write(msg + '\n')

log("=" * 72)
log("[Project RDL] C-344 — B-68 E비 결정인자 분석")
log(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}")
log(f"정밀도: {mpmath.mp.dps} dps")
log("=" * 72)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 지표 정의 + E비 (기존 결과에서 추출)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 이산로그 테이블
_DLOG7 = {1: 0, 2: 2, 3: 1, 4: 4, 5: 5, 6: 3}  # base g=3 mod 7
_DLOG11 = {1: 0, 2: 1, 3: 8, 4: 2, 5: 4, 6: 9, 7: 7, 8: 3, 9: 6, 10: 5}  # base g=2 mod 11

def make_chi(q, j, dlog_table, group_order):
    """일반 디리클레 지표 χ_j mod q 값 배열 생성"""
    omega = mpmath.exp(2j * mpmath.pi / group_order)
    omega_pow = [mpmath.power(omega, k) for k in range(group_order)]
    chi = [mpmath.mpf(0)] * q
    for n in range(1, q):
        if n in dlog_table:
            k = dlog_table[n]
            chi[n] = omega_pow[(j * k) % group_order]
        # else: chi[n] = 0 (gcd(n,q) > 1)
    return chi

def chi_parity(chi_vals, q):
    """a=1 if χ(-1)=-1, a=0 if χ(-1)=+1"""
    val = chi_vals[q - 1]  # χ(-1) = χ(q-1)
    if abs(val - 1) < 1e-10:
        return 0
    elif abs(val + 1) < 1e-10:
        return 1
    else:
        raise ValueError(f"χ(-1) = {val}, 실수 ±1이 아님")

def chi_order(j, group_order):
    """지표의 order = group_order / gcd(j, group_order)"""
    from math import gcd
    return group_order // gcd(j, group_order)

def gauss_sum(chi_vals, q):
    """τ(χ) = Σ_{n=0}^{q-1} χ(n) e^{2πin/q}"""
    tau = mpmath.mpc(0, 0)
    for n in range(q):
        if abs(chi_vals[n]) > 1e-30:
            tau += chi_vals[n] * mpmath.exp(2j * mpmath.pi * n / q)
    return tau

def compute_L1(chi_vals, q, a):
    """L(1, χ) = -(1/q) Σ_{a=1}^{q-1} χ(a) ψ(a/q)  (유한합, 정확)
    mpmath.dirichlet이 복소 지표에서 오작동하므로 digamma 공식 사용."""
    q_int = int(q)
    total = mpmath.mpc(0, 0)
    for n in range(1, q_int):
        cv = chi_vals[n]
        if abs(cv) > 1e-30:
            total += cv * mpmath.digamma(mpmath.mpf(n) / q_int)
    return -total / q_int

def compute_gamma1(chi_vals, q, a, t_min=1.0, t_max=200.0, n_scan=4000):
    """첫 영점 γ₁ 탐색 (임계선 위)"""
    # completed L 값의 부호변화 탐색
    char_info = {'chi': chi_vals, 'q': q, 'a': a}

    ts = np.linspace(t_min, t_max, n_scan)
    prev_re = None

    for i, t in enumerate(ts):
        s = mpmath.mpc(0.5, t)
        try:
            val = completed_L(s, char_info)
            re_val = float(mpmath.re(val))
        except:
            continue

        if prev_re is not None and prev_re * re_val < 0:
            # 부호 변화 → brent 탐색
            t_lo, t_hi = float(ts[i-1]), float(t)
            try:
                def f(tt):
                    ss = mpmath.mpc(0.5, tt)
                    return float(mpmath.re(completed_L(ss, char_info)))
                from scipy.optimize import brentq
                gamma1 = brentq(f, t_lo, t_hi, xtol=1e-12)
                return gamma1
            except:
                return (t_lo + t_hi) / 2
        prev_re = re_val

    return None

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 전 지표 데이터 수집
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# E비 데이터 (기존 결과에서 수동 추출 — 재계산 불필요)
# format: (label, q, j_or_id, order, a, E_ratio, conjugate_label)
ERATIO_DATA = [
    # q=3: φ=2, 비자명 1개
    ('chi3', 3, 1, 2, 1, 21.6, None),
    # q=4: φ=2, 비자명 1개
    ('chi4', 4, 1, 2, 1, 454.6, None),
    # q=5: 번들 검증에서 1개만 테스트 (order 4, a=1, 복소)
    ('chi5', 5, 1, 4, 1, 14.5, None),

    # q=7: φ=6, 비자명 5개 (j=1..5)
    ('chi7_1', 7, 1, 6, 1, 104.4, 'chi7_5'),
    ('chi7_2', 7, 2, 3, 0, 30.1, 'chi7_4'),
    ('chi7_3', 7, 3, 2, 1, 166.9, None),  # Legendre
    ('chi7_4', 7, 4, 3, 0, 337.8, 'chi7_2'),
    ('chi7_5', 7, 5, 6, 1, 557.3, 'chi7_1'),

    # q=8: (Z/8Z)* = Z/2×Z/2. 모든 비자명 지표 order 2 (실수).
    # χ_B: primitive cond=8, odd. χ_C: primitive cond=8, even. (자기 켤레)
    ('chi8_B', 8, 'B', 2, 1, 53.1, None),  # 원시 cond=8, 홀
    ('chi8_C', 8, 'C', 2, 0, 20.4, None),  # 원시 cond=8, 짝

    # q=11: φ=10, 비자명 9개 (j=1..9, 모두 원시)
    ('chi11_1', 11, 1, 10, 1, 75.9, 'chi11_9'),
    ('chi11_2', 11, 2, 5, 0, 364.3, 'chi11_8'),
    ('chi11_3', 11, 3, 10, 1, 6.3, 'chi11_7'),
    ('chi11_4', 11, 4, 5, 0, 7.1, 'chi11_6'),
    ('chi11_5', 11, 5, 2, 1, 60.2, None),  # Legendre
    ('chi11_6', 11, 6, 5, 0, 1023.4, 'chi11_4'),
    ('chi11_7', 11, 7, 10, 1, 25.3, 'chi11_3'),
    ('chi11_8', 11, 8, 5, 0, 326.5, 'chi11_2'),
    ('chi11_9', 11, 9, 10, 1, 111.8, 'chi11_1'),
]

log(f"\n총 지표 수: {len(ERATIO_DATA)}개")
log("")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 각 지표별 해석적 양 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("━" * 72)
log("1단계: 해석적 양 계산 (|L(1,χ)|, τ(χ), γ₁)")
log("━" * 72)

results = []

for (label, q, j_id, order, a, e_ratio, conj) in ERATIO_DATA:
    log(f"\n[{label}] q={q}, order={order}, a={a}, E비={e_ratio}×")

    # 지표 값 생성
    if q == 3:
        chi_vals = [0, 1, -1]
    elif q == 4:
        chi_vals = [0, 1, 0, -1]
    elif q == 5:
        if j_id == 1:
            chi_vals = [0, 1, 1j, -1j, -1]
        else:  # j=2 (Legendre, order 2)
            chi_vals = [0, 1, -1, -1, 1]
    elif q == 7:
        chi_vals = make_chi(7, j_id, _DLOG7, 6)
    elif q == 8:
        # q=8: (Z/8Z)* = {1,3,5,7} ≅ Z/2 × Z/2. 모든 비자명 order 2 (실수).
        # χ_A: [+1,-1,+1,-1] 유도(cond=4) — 제외
        # χ_B: [+1,+1,-1,-1] 원시(cond=8, 홀, Kronecker (8/n))
        # χ_C: [+1,-1,-1,+1] 원시(cond=8, 짝, Kronecker (-8/n))
        if j_id == 'B':
            chi_vals = [0, 1, 0, 1, 0, -1, 0, -1]
        else:
            chi_vals = [0, 1, 0, -1, 0, -1, 0, 1]
    elif q == 11:
        chi_vals = make_chi(11, j_id, _DLOG11, 10)

    chi_mpmath = [mpmath.mpc(v) for v in chi_vals]

    # |L(1,χ)|
    L1_val = compute_L1(chi_mpmath, q, a)
    L1_abs = float(abs(L1_val)) if L1_val is not None else None
    L1_arg = float(mpmath.arg(L1_val)) if L1_val is not None else None
    log(f"  |L(1,χ)| = {L1_abs:.6f}" + (f", arg = {L1_arg:.4f}" if L1_arg else ""))

    # Gauss sum τ(χ)
    tau = gauss_sum(chi_mpmath, q)
    tau_abs = float(abs(tau))
    tau_normalized = tau_abs / float(mpmath.sqrt(q))
    tau_arg = float(mpmath.arg(tau))
    log(f"  |τ(χ)| = {tau_abs:.6f}, |τ|/√q = {tau_normalized:.6f}, arg(τ) = {tau_arg:.4f}")

    # Root number w(χ) = τ(χ) / (i^a √q)
    ia = mpmath.power(1j, a)
    w_chi = tau / (ia * mpmath.sqrt(q))
    w_abs = float(abs(w_chi))
    w_arg = float(mpmath.arg(w_chi))
    log(f"  w(χ) = τ/(i^a √q): |w| = {w_abs:.6f}, arg(w) = {w_arg:.4f}")

    # γ₁ (첫 영점) — 느리므로 q=11은 skip하고 결과파일에서 추출
    gamma1 = None

    # q=7에서 결과에서 추출한 γ₁ (직접 읽기)
    gamma1_known = {
        # q=3,4: 결과 파일 부재 — 계산으로 구함
        # q=7: C-324 결과
        'chi7_1': 13.85,  # j=1, ord6
        'chi7_2': 10.74,  # j=2, ord3
        'chi7_3': 11.16,  # j=3, ord2 (Legendre, 첫 양의 영점)
        'chi7_4': 11.01,  # j=4, ord3
        'chi7_5': 12.26,  # j=5, ord6
        # q=8: C-326 결과
        'chi8_B': 12.34,  # 원시 cond=8, 홀
        'chi8_C': 10.81,  # 원시 cond=8, 짝
        # q=11: C-328 결과
        'chi11_1': 11.59,
        'chi11_2': 12.25,
        'chi11_3': 11.26,
        'chi11_4': 11.09,
        'chi11_5': 10.11,
        'chi11_6': 11.28,
        'chi11_7': 10.45,
        'chi11_8': 10.92,
        'chi11_9': 11.01,
    }

    if label in gamma1_known:
        gamma1 = gamma1_known[label]
        log(f"  γ₁ = {gamma1:.2f} (결과 파일)")
    else:
        log(f"  γ₁: 계산 중...")
        gamma1 = compute_gamma1(chi_mpmath, q, a, t_min=1.0, t_max=100.0, n_scan=2000)
        if gamma1:
            log(f"  γ₁ = {gamma1:.4f}")
        else:
            log(f"  γ₁: 탐색 실패 ([1,100] 범위)")

    results.append({
        'label': label,
        'q': q,
        'j': j_id,
        'order': order,
        'a': a,
        'e_ratio': e_ratio,
        'log_e_ratio': np.log(e_ratio),
        'L1_abs': L1_abs,
        'L1_arg': L1_arg,
        'tau_abs': tau_abs,
        'tau_norm': tau_normalized,
        'tau_arg': tau_arg,
        'w_abs': w_abs,
        'w_arg': w_arg,
        'gamma1': gamma1,
        'conj': conj,
    })

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2단계: 전체 표 출력
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("\n" + "━" * 72)
log("2단계: 전체 데이터 표")
log("━" * 72)

header = f"{'지표':<12} {'q':>3} {'ord':>4} {'a':>2} {'E비':>8} {'|L(1)|':>8} {'|τ|/√q':>7} {'arg(τ)':>7} {'arg(w)':>7} {'γ₁':>7}"
log(header)
log("-" * len(header))

for r in results:
    gamma1_str = f"{r['gamma1']:.2f}" if r['gamma1'] else "?"
    log(f"{r['label']:<12} {r['q']:>3} {r['order']:>4} {r['a']:>2} {r['e_ratio']:>8.1f} "
        f"{r['L1_abs']:>8.4f} {r['tau_norm']:>7.4f} {r['tau_arg']:>7.3f} "
        f"{r['w_arg']:>7.3f} {gamma1_str:>7}")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3단계: Spearman 상관 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("\n" + "━" * 72)
log("3단계: log(E비) vs 해석적 양 — Spearman 순위 상관")
log("━" * 72)

log_e = np.array([r['log_e_ratio'] for r in results])
n = len(results)

candidates = {
    '|L(1,χ)|': np.array([r['L1_abs'] for r in results]),
    'log|L(1,χ)|': np.log(np.array([r['L1_abs'] for r in results])),
    '1/|L(1,χ)|': 1.0 / np.array([r['L1_abs'] for r in results]),
    'arg(τ)': np.array([r['tau_arg'] for r in results]),
    '|arg(τ)|': np.abs(np.array([r['tau_arg'] for r in results])),
    'arg(w)': np.array([r['w_arg'] for r in results]),
    '|arg(w)|': np.abs(np.array([r['w_arg'] for r in results])),
    'order': np.array([float(r['order']) for r in results]),
    'parity a': np.array([float(r['a']) for r in results]),
}

# γ₁ (결측 없는 것만)
gamma1_valid = [(r['gamma1'], r['log_e_ratio']) for r in results if r['gamma1'] is not None]
if gamma1_valid:
    g1_arr = np.array([g[0] for g in gamma1_valid])
    le_arr = np.array([g[1] for g in gamma1_valid])
    rho_g1, p_g1 = stats.spearmanr(g1_arr, le_arr)
    candidates['γ₁'] = None  # 별도 처리

log(f"\n{'변수':<16} {'ρ_Spearman':>12} {'p-value':>12} {'판정':>8}")
log("-" * 52)

best_rho = 0
best_name = ""

for name, vals in candidates.items():
    if vals is None:
        # γ₁ 별도
        rho, p = rho_g1, p_g1
        log(f"{'γ₁':<16} {rho:>12.4f} {p:>12.2e} {'★' if abs(rho)>0.5 else ''}")
        if abs(rho) > abs(best_rho):
            best_rho, best_name = rho, 'γ₁'
        continue

    rho, p = stats.spearmanr(vals, log_e)
    marker = "★★★" if abs(rho) > 0.7 else "★★" if abs(rho) > 0.5 else "★" if abs(rho) > 0.3 else ""
    log(f"{name:<16} {rho:>12.4f} {p:>12.2e} {marker}")
    if abs(rho) > abs(best_rho):
        best_rho, best_name = rho, name

log(f"\n최강 상관: {best_name} (ρ={best_rho:.4f})")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4단계: q=11 내부 분석 (동일 conductor에서 변수 고립)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("\n" + "━" * 72)
log("4단계: q=11 내부 분석 (conductor 고정 → q 효과 제거)")
log("━" * 72)

chi11_results = [r for r in results if r['q'] == 11]
n11 = len(chi11_results)

log_e11 = np.array([r['log_e_ratio'] for r in chi11_results])

cands11 = {
    '|L(1,χ)|': np.array([r['L1_abs'] for r in chi11_results]),
    'log|L(1,χ)|': np.log(np.array([r['L1_abs'] for r in chi11_results])),
    '1/|L(1,χ)|': 1.0 / np.array([r['L1_abs'] for r in chi11_results]),
    'arg(τ)': np.array([r['tau_arg'] for r in chi11_results]),
    '|arg(τ)|': np.abs(np.array([r['tau_arg'] for r in chi11_results])),
    'arg(w)': np.array([r['w_arg'] for r in chi11_results]),
    '|arg(w)|': np.abs(np.array([r['w_arg'] for r in chi11_results])),
    'order': np.array([float(r['order']) for r in chi11_results]),
    'parity a': np.array([float(r['a']) for r in chi11_results]),
}

g1_11 = np.array([r['gamma1'] for r in chi11_results if r['gamma1'] is not None])
le_11 = np.array([r['log_e_ratio'] for r in chi11_results if r['gamma1'] is not None])
if len(g1_11) > 3:
    rho_g1_11, p_g1_11 = stats.spearmanr(g1_11, le_11)
    cands11['γ₁'] = None

log(f"\nq=11 (n={n11})")
log(f"{'변수':<16} {'ρ_Spearman':>12} {'p-value':>12} {'판정':>8}")
log("-" * 52)

best_rho11 = 0
best_name11 = ""

for name, vals in cands11.items():
    if vals is None:
        rho, p = rho_g1_11, p_g1_11
        log(f"{'γ₁':<16} {rho:>12.4f} {p:>12.2e} {'★' if abs(rho)>0.5 else ''}")
        if abs(rho) > abs(best_rho11):
            best_rho11, best_name11 = rho, 'γ₁'
        continue
    rho, p = stats.spearmanr(vals, log_e11)
    marker = "★★★" if abs(rho) > 0.7 else "★★" if abs(rho) > 0.5 else "★" if abs(rho) > 0.3 else ""
    log(f"{name:<16} {rho:>12.4f} {p:>12.2e} {marker}")
    if abs(rho) > abs(best_rho11):
        best_rho11, best_name11 = rho, name

log(f"\nq=11 최강 상관: {best_name11} (ρ={best_rho11:.4f})")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5단계: 켤레 쌍 E비 비율 vs 해석적 양 비율
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("\n" + "━" * 72)
log("5단계: 켤레 쌍 (χ, χ̄) E비 비율 vs |L(1,χ)|비율 / arg(τ) 차이")
log("━" * 72)

# 켤레 쌍 구성
pairs = []
seen = set()
for r in results:
    if r['conj'] and r['label'] not in seen:
        conj_r = next((x for x in results if x['label'] == r['conj']), None)
        if conj_r:
            pairs.append((r, conj_r))
            seen.add(r['label'])
            seen.add(conj_r['label'])

log(f"\n켤레 쌍 수: {len(pairs)}")
log(f"\n{'쌍':<25} {'E비(A)':>8} {'E비(B)':>8} {'비율':>7} {'|L1|A':>8} {'|L1|B':>8} {'L1비율':>7} {'argτA':>7} {'argτB':>7}")
log("-" * 100)

e_ratios_pair = []
l1_ratios_pair = []
tau_arg_diffs = []

for r_a, r_b in pairs:
    # A가 E비가 큰 쪽
    if r_a['e_ratio'] < r_b['e_ratio']:
        r_a, r_b = r_b, r_a

    e_pair_ratio = r_a['e_ratio'] / r_b['e_ratio']
    l1_pair_ratio = r_a['L1_abs'] / r_b['L1_abs'] if r_b['L1_abs'] > 0 else float('inf')
    tau_diff = abs(r_a['tau_arg'] - r_b['tau_arg'])

    e_ratios_pair.append(e_pair_ratio)
    l1_ratios_pair.append(l1_pair_ratio)
    tau_arg_diffs.append(tau_diff)

    log(f"{r_a['label']}↔{r_b['label']:<12} {r_a['e_ratio']:>8.1f} {r_b['e_ratio']:>8.1f} {e_pair_ratio:>7.1f} "
        f"{r_a['L1_abs']:>8.4f} {r_b['L1_abs']:>8.4f} {l1_pair_ratio:>7.2f} "
        f"{r_a['tau_arg']:>7.3f} {r_b['tau_arg']:>7.3f}")

if len(pairs) >= 4:
    rho_pair, p_pair = stats.spearmanr(e_ratios_pair, l1_ratios_pair)
    log(f"\nρ_S(E비 비율, |L(1)|비율) = {rho_pair:.4f} (p={p_pair:.3e})")
    rho_tau, p_tau = stats.spearmanr(e_ratios_pair, tau_arg_diffs)
    log(f"ρ_S(E비 비율, |Δarg(τ)|) = {rho_tau:.4f} (p={p_tau:.3e})")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6단계: E(0.5) vs E(0.3) 분리 — 어느 쪽이 변동의 원인인가?
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("\n" + "━" * 72)
log("6단계: 해석적 예측 — Γ인자 기반 E(0.3) 추정")
log("━" * 72)

log("""
E비 = E(0.5) / E(0.3).
E(0.5)은 영점 곡률 기여로 크게 변동.
E(0.3)은 영점에서 먼 σ에서의 배경 수준 — Γ인자로 지배.

Γ((s+a)/2) at σ=0.3 vs σ=0.5:
  a=0: Γ(0.15+it/2) / Γ(0.25+it/2) — 짝
  a=1: Γ(0.65+it/2) / Γ(0.75+it/2) — 홀

만약 E(0.3)이 a에만 의존한다면, E비 변동은 E(0.5)에서 기인.
E(0.5)는 영점 밀도와 직접 연관.
""")

# 이론적 배경 수준 비교
for a_val in [0, 1]:
    chars_a = [r for r in results if r['a'] == a_val]
    e_ratios_a = [r['e_ratio'] for r in chars_a]
    log(f"  a={a_val} ({'짝' if a_val==0 else '홀'}): n={len(chars_a)}, "
        f"E비 범위 [{min(e_ratios_a):.1f}, {max(e_ratios_a):.1f}], "
        f"변동계수 CV = {np.std(e_ratios_a)/np.mean(e_ratios_a):.2f}")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7단계: 종합 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("\n" + "━" * 72)
log("7단계: 종합 판정")
log("━" * 72)

log(f"""
[결과 요약]
  전체 (n={n}):  최강 상관 = {best_name} (ρ={best_rho:.4f})
  q=11 (n={n11}): 최강 상관 = {best_name11} (ρ={best_rho11:.4f})

[해석]
  (분석 후 작성)

[판정]: (분석 후)
""")

out.close()
log(f"\n결과 저장: {RESULT_FILE}")
print(f"\n완료. 결과: {RESULT_FILE}")
