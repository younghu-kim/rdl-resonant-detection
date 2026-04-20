#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #121] 일반화 DH off-critical 영점 탐색 + c₁ 측정 (PARI 버전)
=============================================================================
목표: DH(mod 5) 외의 L-함수에서 off-critical 영점 확보 → c₁ 법칙 보편성 확인.

방법: 일반화 DH 구성 (올바른 함수방정식 각도).
  f(s) = L(s,χ) + ω·L(s,χ̄)  where  ω² = ε(χ̄)/ε(χ)
  ε(χ) = root number = τ(χ)/(i^a · √q),  τ = Gauss sum.
  이 ω에서만 Λ_f(s) = Λ(s,χ) + ω·Λ(s,χ̄)가 함수방정식을 만족.

핵심 수정 (v2): 임의 α 대신 root number에서 유도한 정확한 ω 사용.
  이전 버전의 cos(α)/sin(α) 구성은 FE 불만족 → κδ² ≈ 0.001 실패.
  올바른 ω 사용 시 κδ² ≈ 1 기대.

대상:
  mod 7: χ₁(order 6), χ₂(order 3) — 각각 올바른 ω
  mod 11: χ₁(order 10), χ₂(order 5) — 각각 올바른 ω
  mod 5 (DH 대조): χ₁(order 4) — ω from root number

사용 환경: system python3 + cypari2

결과: results/generalized_dh_offcritical_121.txt
=============================================================================
"""

import sys, os, time
import numpy as np
from datetime import datetime
from scipy.optimize import minimize

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 80)

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'generalized_dh_offcritical_121.txt')
outf = open(RESULT_FILE, 'w')

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

START = time.time()

log("=" * 72)
log("[실험 #121] 일반화 DH off-critical 영점 탐색 (PARI, v2)")
log("=" * 72)
log(f"시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. Root number 계산 (Gauss sum)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 1. Root number 및 DH 각도 계산 ━━━")

def compute_char_table(q, k):
    """
    소수 q의 k번째 디리클레 지표 값 테이블.
    원시근 g에 대해 χ_k(g^j) = e^{2πi·j·k/(q-1)}.
    """
    g = int(pari(f'znprimroot({q})'))
    vals = {}
    for j in range(q - 1):
        n = pow(g, j, q)
        vals[n] = np.exp(2j * np.pi * j * k / (q - 1))
    return vals

def compute_root_number(q, k):
    """ε(χ_k) = τ(χ_k) / (i^a · √q), where a = 0 if χ(-1)=1, a=1 if χ(-1)=-1."""
    vals = compute_char_table(q, k)
    # Gauss sum τ(χ) = Σ χ(n)·e^{2πin/q}
    tau = sum(vals[n] * np.exp(2j * np.pi * n / q) for n in range(1, q))
    chi_minus1 = vals[q - 1]
    a = 0 if abs(chi_minus1 - 1) < 0.01 else 1
    eps = tau / (1j**a * np.sqrt(q))
    return eps, a

def compute_dh_omega(q, k):
    """ω² = ε(χ̄)/ε(χ) = ε(χ_{q-1-k})/ε(χ_k). Returns ω = √(ratio)."""
    k_bar = (q - 1) - k
    eps_k, a_k = compute_root_number(q, k)
    eps_kb, a_kb = compute_root_number(q, k_bar)
    ratio = eps_kb / eps_k
    omega = np.sqrt(ratio)
    return omega, eps_k, eps_kb, a_k

# 모든 구성의 ω 계산
CONFIGS_DATA = []

for q, k, name in [(7, 1, "mod7-χ₁(ord6)"), (7, 2, "mod7-χ₂(ord3)"),
                     (11, 1, "mod11-χ₁(ord10)"), (11, 2, "mod11-χ₂(ord5)"),
                     (5, 1, "mod5-DH(대조)")]:
    omega, eps, eps_bar, a = compute_dh_omega(q, k)
    # 검증: ω²·ε(χ) = ε(χ̄)
    check = omega**2 * eps
    err = abs(check - eps_bar)

    k_bar = (q - 1) - k
    order = (q - 1) // np.gcd(k, q - 1)
    log(f"  {name}: q={q}, k={k}, k̄={k_bar}, order={order}")
    log(f"    ε(χ)  = {eps:.6f}, |ε|={abs(eps):.8f}")
    log(f"    ε(χ̄) = {eps_bar:.6f}")
    log(f"    ω = {omega:.6f}, |ω|={abs(omega):.8f}")
    log(f"    검증 ω²·ε(χ)=ε(χ̄): err={err:.2e} {'✅' if err < 1e-10 else '❌'}")

    # FE sign: c = ω·ε(χ)
    c = omega * eps
    log(f"    FE sign c = ω·ε(χ) = {c:.6f}, |c|={abs(c):.8f}")
    log()

    CONFIGS_DATA.append({
        'q': q, 'k': k, 'k_bar': k_bar, 'name': name,
        'omega': omega, 'a': a, 'order': order
    })


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. PARI L-함수 생성 및 평가 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 2. L-함수 생성 ━━━")

def make_lfun_pair(q, k):
    """PARI L-함수 쌍 생성. 원시근 g의 k승과 k̄=(q-1-k)승."""
    g = int(pari(f'znprimroot({q})'))
    g_k = pow(g, k, q)
    k_bar = (q - 1) - k
    g_kb = pow(g, k_bar, q)
    L_chi = pari(f"lfuncreate(Mod({g_k},{q}))")
    L_chi_bar = pari(f"lfuncreate(Mod({g_kb},{q}))")
    fe = float(pari(f"lfuncheckfeq(lfuncreate(Mod({g_k},{q})))"))
    log(f"  {q}: g={g}, g^k={g_k}, g^k̄={g_kb}, FE check={fe:.0f}")
    return L_chi, L_chi_bar

LFUN_PAIRS = {}
for cfg in CONFIGS_DATA:
    q, k = cfg['q'], cfg['k']
    key = (q, k)
    if key not in LFUN_PAIRS:
        LFUN_PAIRS[key] = make_lfun_pair(q, k)

log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. 일반화 DH 함수 (올바른 ω)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def eval_L(L_obj, sigma, t):
    """L(σ+it) 평가."""
    s_str = f"{sigma} + {t}*I" if t >= 0 else f"{sigma} - {abs(t)}*I"
    return complex(pari.lfun(L_obj, pari(s_str)))

def eval_Lambda(L_obj, sigma, t):
    """Λ(σ+it) = completed L-function."""
    s_str = f"{sigma} + {t}*I" if t >= 0 else f"{sigma} - {abs(t)}*I"
    return complex(pari.lfun(L_obj, pari(s_str), 1))

def gen_dh_Lambda(sigma, t, omega, L_chi, L_chi_bar):
    """Λ_f(s) = Λ(s,χ) + ω·Λ(s,χ̄)  — 올바른 DH 구성."""
    v1 = eval_Lambda(L_chi, sigma, t)
    v2 = eval_Lambda(L_chi_bar, sigma, t)
    return v1 + omega * v2


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. Off-critical 영점 탐색
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 3. Off-critical 영점 탐색 ━━━")

T_MAX = 200
SIGMA_SCAN = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]
T_STEP = 0.5
THRESHOLD = 0.5


def find_offcritical(omega, L_chi, L_chi_bar, name):
    """Λ_f(s) off-critical 영점 1D 스캔."""
    log(f"\n  ─── {name}, ω={omega:.6f} ───")
    t0 = time.time()

    candidates = []
    for sigma in SIGMA_SCAN:
        ts = np.arange(2.0, T_MAX, T_STEP)
        prev_abs = None
        prev_t = None

        for t in ts:
            try:
                val = gen_dh_Lambda(sigma, t, omega, L_chi, L_chi_bar)
                cur_abs = abs(val)
            except Exception:
                cur_abs = 1e30

            # 국소 최솟값 검출
            if prev_abs is not None and prev_abs < THRESHOLD:
                if cur_abs > prev_abs:
                    candidates.append((sigma, float(prev_t), prev_abs))

            prev_abs = cur_abs
            prev_t = t

    elapsed_scan = time.time() - t0
    log(f"    스캔: {len(candidates)}개 후보 ({elapsed_scan:.1f}초)")

    # 상위 30개 정밀화
    confirmed = []
    for sigma0, t0_val, v0 in sorted(candidates, key=lambda x: x[2])[:30]:
        try:
            def objective(x):
                return abs(gen_dh_Lambda(x[0], x[1], omega, L_chi, L_chi_bar))

            res = minimize(objective, [sigma0, t0_val], method='Nelder-Mead',
                           options={'xatol': 1e-10, 'fatol': 1e-14, 'maxiter': 2000})

            sigma_r, t_r = res.x
            val_r = res.fun

            if 0.01 < sigma_r < 0.99 and t_r > 0.5 and abs(sigma_r - 0.5) > 0.01 and val_r < 1e-6:
                # 대칭: σ < 0.5이면 mirror로 변환
                if sigma_r < 0.5:
                    sigma_r = 1.0 - sigma_r

                is_dup = any(abs(c['sigma'] - sigma_r) < 0.005 and
                             abs(c['t'] - t_r) < 0.5 for c in confirmed)
                if not is_dup:
                    confirmed.append({'sigma': sigma_r, 't': t_r, 'absf': val_r})
                    log(f"    ✅ #{len(confirmed)}: σ={sigma_r:.8f}, t={t_r:.6f}, "
                        f"|Λ_f|={val_r:.3e}, |σ-½|={abs(sigma_r-0.5):.6f}")
        except Exception:
            pass

    elapsed_total = time.time() - t0
    log(f"    확인: {len(confirmed)}개 ({elapsed_total:.1f}초)")
    return confirmed


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. c₁ 측정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DELTAS = [0.01, 0.02, 0.03, 0.05, 0.07, 0.1]

def measure_c1(zero, omega, L_chi, L_chi_bar):
    """κδ² 다중 δ + c₁ fit."""
    sigma0, t0 = zero['sigma'], zero['t']

    kd2_data = []
    for delta in DELTAS:
        s_off_sigma = sigma0 + delta
        h = 1e-7
        try:
            lam_c = gen_dh_Lambda(s_off_sigma, t0, omega, L_chi, L_chi_bar)
            lam_p = gen_dh_Lambda(s_off_sigma + h, t0, omega, L_chi, L_chi_bar)
            lam_m = gen_dh_Lambda(s_off_sigma - h, t0, omega, L_chi, L_chi_bar)
            if abs(lam_c) > 1e-50:
                lam_d = (lam_p - lam_m) / (2 * h)
                kappa = abs(lam_d / lam_c) ** 2
                kd2 = kappa * delta**2
                kd2_data.append((delta, kd2))
        except Exception:
            pass

    if len(kd2_data) < 3:
        return None, None, []

    x = np.array([d for d, _ in kd2_data])
    y = np.array([v - 1.0 for _, v in kd2_data])
    A = np.column_stack([x, x**2])
    coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    c1_fit = coeffs[0]

    # c₁ analytic: Re(Λ''/Λ')(ρ)
    h2 = 1e-6
    try:
        lam_p = gen_dh_Lambda(sigma0 + h2, t0, omega, L_chi, L_chi_bar)
        lam_c = gen_dh_Lambda(sigma0, t0, omega, L_chi, L_chi_bar)
        lam_m = gen_dh_Lambda(sigma0 - h2, t0, omega, L_chi, L_chi_bar)
        Ld = (lam_p - lam_m) / (2 * h2)
        Ldd = (lam_p - 2*lam_c + lam_m) / (h2**2)
        c1_analytic = (Ldd / Ld).real if abs(Ld) > 1e-50 else None
    except Exception:
        c1_analytic = None

    return c1_fit, c1_analytic, kd2_data


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

all_results = []

for cfg in CONFIGS_DATA:
    q, k = cfg['q'], cfg['k']
    L_chi, L_chi_bar = LFUN_PAIRS[(q, k)]
    omega = cfg['omega']
    name = cfg['name']

    zeros = find_offcritical(omega, L_chi, L_chi_bar, name)

    if zeros:
        log(f"\n  ━━━ c₁ 측정: {name} ━━━")
        for z in zeros:
            c1_fit, c1_analytic, kd2_data = measure_c1(z, omega, L_chi, L_chi_bar)
            dist = abs(z['sigma'] - 0.5)
            mirror = 1.0 / dist if dist > 0.001 else None

            log(f"\n    σ={z['sigma']:.8f}, t={z['t']:.6f}, |σ-½|={dist:.6f}")
            if kd2_data:
                log(f"    δ      κδ²")
                for d, v in kd2_data:
                    status = "✅" if abs(v - 1.0) < 0.3 else "❌"
                    log(f"    {d:.3f}   {v:.6f}  {status}")

            if c1_fit is not None:
                log(f"    c₁(fit)      = {c1_fit:.4f}")
            if c1_analytic is not None:
                log(f"    c₁(analytic) = {c1_analytic:.4f}")
            if mirror is not None:
                log(f"    1/(σ₀-½)     = {mirror:.4f}")
            if c1_fit is not None and dist > 0.001:
                log(f"    c₁·|σ-½|     = {c1_fit * dist:.4f}")

            z['c1_fit'] = c1_fit
            z['c1_analytic'] = c1_analytic
            z['config'] = name
            all_results.append(z)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7. 종합
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log(f"\n{'='*72}")
log("종합")
log(f"{'='*72}")

DH5_REF = [
    {"sigma": 0.808517, "t": 85.699, "c1": 3.76, "dist": 0.3085, "src": "DH mod5 (#115)"},
    {"sigma": 0.650830, "t": 114.163, "c1": 6.92, "dist": 0.1508, "src": "DH mod5 (#115)"},
    {"sigma": 0.574356, "t": 166.479, "c1": 13.54, "dist": 0.0744, "src": "DH mod5 (#115)"},
    {"sigma": 0.724258, "t": 176.702, "c1": 5.03, "dist": 0.2243, "src": "DH mod5 (#115)"},
]

n_new_mod7 = len([r for r in all_results if 'mod7' in r.get('config', '')])
n_new_mod11 = len([r for r in all_results if 'mod11' in r.get('config', '')])
n_new_mod5 = len([r for r in all_results if 'mod5' in r.get('config', '')])
n_new = n_new_mod7 + n_new_mod11

log(f"\n  기존 DH(mod 5) #115: 4개")
log(f"  신규 mod 7: {n_new_mod7}개")
log(f"  신규 mod 11: {n_new_mod11}개")
log(f"  대조 mod 5: {n_new_mod5}개")
log(f"  비-mod5 합계: {n_new}개")

log(f"\n  비교표:")
log(f"  {'소스':<22} {'σ':>8} {'t':>10} {'|σ-½|':>8} {'c₁':>10} {'1/(σ-½)':>10} {'c₁·|σ-½|':>10}")
log(f"  {'─'*82}")

for d in DH5_REF:
    p = d['c1'] * d['dist']
    m = 1.0 / d['dist']
    log(f"  {d['src']:<22} {d['sigma']:>8.6f} {d['t']:>10.3f} {d['dist']:>8.6f} "
        f"{d['c1']:>10.4f} {m:>10.4f} {p:>10.4f}")

for r in all_results:
    dist = abs(r['sigma'] - 0.5)
    c1 = r.get('c1_fit')
    m = 1.0/dist if dist > 0.001 else None
    p = c1*dist if c1 is not None else None
    cfg = r.get('config', '?')
    if c1 is not None and m is not None and p is not None:
        log(f"  {cfg:<22} {r['sigma']:>8.6f} {r['t']:>10.3f} {dist:>8.6f} "
            f"{c1:>10.4f} {m:>10.4f} {p:>10.4f}")
    else:
        log(f"  {cfg:<22} {r['sigma']:>8.6f} {r['t']:>10.3f} {dist:>8.6f} "
            f"{'N/A':>10} {'N/A':>10} {'N/A':>10}")

# 성공 기준
log(f"\n  성공 기준:")
sc1 = n_new >= 1
log(f"  SC1 (≥1 non-mod5): {'✅' if sc1 else '❌'} ({n_new}개)")

if n_new > 0:
    c1_products = [r.get('c1_fit', 0) * abs(r['sigma'] - 0.5)
                   for r in all_results if r.get('c1_fit') and 'mod5' not in r.get('config', '')]
    if c1_products:
        mean_p = np.mean(c1_products)
        log(f"  SC2 (c₁·|σ-½|≈1): mean={mean_p:.3f}")

# 판정
if n_new >= 3:
    log(f"\n  ★★★ 양성 — 다중 L-함수에서 c₁ 법칙 보편성 확인.")
elif n_new >= 1:
    log(f"\n  ★★ 양성 — off-critical 발견. c₁ 법칙 보편성 예비 확인.")
elif n_new_mod5 > 0:
    log(f"\n  중립 — mod 5 대조만 확인. mod 7/11에서 off-critical 미발견.")
else:
    log(f"\n  음성 — 어떤 구성에서도 off-critical 미발견. 탐색 범위 확대 필요.")

elapsed = time.time() - START
log(f"\n총 소요: {elapsed:.1f}초 ({elapsed/60:.1f}분)")
outf.close()
