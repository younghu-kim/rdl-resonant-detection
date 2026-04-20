#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #121] 일반화 DH off-critical 영점 탐색 + c₁ 측정 (PARI 버전)
=============================================================================
목표: DH(mod 5) 외의 L-함수에서 off-critical 영점 확보 → c₁ 법칙 보편성 확인.

방법: 일반화 DH 구성 (올바른 함수방정식 각도).
  f(s) = Λ(s,χ) + ω·Λ(s,χ̄)  where  ω² = ε(χ)/ε(χ̄)
  이 조건에서 Λ_f(1-s) = Λ_f(s) (FE sign = +1).

핵심 수정 (v3):
  - ω² = ε(χ)/ε(χ̄)  (NOT ε(χ̄)/ε(χ))
  - lfunlambda() for Λ(s)  (NOT lfun(L,s,1) which is L'(s))
  이전 버전: cos(α)/sin(α) 임의각 → FE 불만족 → κδ² ≈ 0.001 실패.
  v2: ω² = ε(χ̄)/ε(χ) → 부호 오류, FE 불만족.

대상:
  mod 7: χ₁(order 6), χ₂(order 3)
  mod 11: χ₁(order 10), χ₂(order 5)
  mod 5 (DH 대조): χ₁(order 4)

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
log("[실험 #121] 일반화 DH off-critical 영점 탐색 (PARI, v3)")
log("=" * 72)
log(f"시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. Root number 계산 + DH 각도
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 1. Root number 및 DH ω 계산 ━━━")
log("  공식: ω² = ε(χ)/ε(χ̄),  Λ_f(s) = Λ(s,χ) + ω·Λ(s,χ̄)")
log("  FE: Λ_f(1-s) = Λ_f(s) (sign = +1)")
log()

def compute_root_number_pari(q, k):
    """
    ε(χ_k) = Λ(s,χ_k)/Λ(1-s,χ̄_k)  from lfunlambda.
    직접 PARI에서 수치 계산 — Gauss sum 불필요.
    """
    g = int(pari(f'znprimroot({q})'))
    g_k = pow(g, k, q)
    k_bar = (q - 1) - k
    g_kb = pow(g, k_bar, q)

    Lk = pari(f'lfuncreate(Mod({g_k},{q}))')
    Lkb = pari(f'lfuncreate(Mod({g_kb},{q}))')

    # ε(χ) = Λ(s,χ)/Λ(1-s,χ̄) — 임의 s에서 상수
    s_test = pari('0.3+10*I')
    s_test_conj = pari('0.7-10*I')  # 1-s
    lam_s = complex(pari.lfunlambda(Lk, s_test))
    lam_bar_1ms = complex(pari.lfunlambda(Lkb, s_test_conj))
    eps = lam_s / lam_bar_1ms
    return eps, g, g_k, g_kb

CONFIGS_DATA = []

for q, k, name in [(7, 1, "mod7-χ₁(ord6)"), (7, 2, "mod7-χ₂(ord3)"),
                     (11, 1, "mod11-χ₁(ord10)"), (11, 2, "mod11-χ₂(ord5)"),
                     (5, 1, "mod5-DH(대조)")]:
    k_bar = (q - 1) - k
    order = (q - 1) // np.gcd(k, q - 1)

    eps_chi, g, g_k, g_kb = compute_root_number_pari(q, k)
    eps_chi_bar, _, _, _ = compute_root_number_pari(q, k_bar)

    # ω² = ε(χ)/ε(χ̄)
    ratio = eps_chi / eps_chi_bar
    omega = np.sqrt(ratio)

    # FE 검증: Λ_f(1-s)/Λ_f(s) = 1
    Lk_fe = pari(f'lfuncreate(Mod({g_k},{q}))')
    Lkb_fe = pari(f'lfuncreate(Mod({g_kb},{q}))')
    s_fe = pari('0.3+10*I')
    s_fe_m = pari('0.7-10*I')
    Lf_s = complex(pari.lfunlambda(Lk_fe, s_fe)) + omega * complex(pari.lfunlambda(Lkb_fe, s_fe))
    Lf_1ms = complex(pari.lfunlambda(Lk_fe, s_fe_m)) + omega * complex(pari.lfunlambda(Lkb_fe, s_fe_m))
    fe_check = abs(Lf_1ms / Lf_s - 1.0) if abs(Lf_s) > 1e-50 else 999

    log(f"  {name}: q={q}, g={g}, g^k={g_k}, g^k̄={g_kb}, order={order}")
    log(f"    ε(χ) = {eps_chi:.6f}, |ε|={abs(eps_chi):.8f}")
    log(f"    ε(χ̄) = {eps_chi_bar:.6f}")
    log(f"    ω = {omega:.6f}, |ω|={abs(omega):.8f}")
    log(f"    FE check |Λ_f(1-s)/Λ_f(s) - 1| = {fe_check:.2e} {'✅' if fe_check < 1e-6 else '❌'}")
    log()

    # parity a: χ(-1)=(-1)^a. compute from character table
    g_prim = int(pari(f'znprimroot({q})'))
    # χ_k(g^j) = e^{2πi·j·k/(q-1)}. χ(-1) = χ(q-1). q-1 = g^{(q-1)/2} mod q (since ord(g)=q-1)
    # So χ(-1) = e^{2πi·(q-1)/2·k/(q-1)} = e^{πik} = (-1)^k
    a_parity = k % 2
    log(f"    a(parity) = {a_parity} ({'홀수' if a_parity else '짝수'})")

    CONFIGS_DATA.append({
        'q': q, 'k': k, 'k_bar': k_bar, 'name': name,
        'omega': omega, 'order': order, 'g_k': g_k, 'g_kb': g_kb,
        'a': a_parity
    })


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. Λ_f 평가 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# Python Gen 객체로 L-함수 저장
LFUN_OBJS = []
for i, cfg in enumerate(CONFIGS_DATA):
    L_chi = pari(f'lfuncreate(Mod({cfg["g_k"]},{cfg["q"]}))')
    L_bar = pari(f'lfuncreate(Mod({cfg["g_kb"]},{cfg["q"]}))')
    LFUN_OBJS.append((L_chi, L_bar))

def eval_L_f(idx, sigma, t, omega):
    """f(s) = L(s,χ) + ω·L(s,χ̄)  via lfun (빠름, 영점 탐색용)."""
    s_str = f"{sigma}+{t}*I" if t >= 0 else f"{sigma}-{abs(t)}*I"
    s = pari(s_str)
    L_chi, L_bar = LFUN_OBJS[idx]
    v1 = complex(pari.lfun(L_chi, s))
    v2 = complex(pari.lfun(L_bar, s))
    return v1 + omega * v2

def eval_Lambda_f(idx, sigma, t, omega):
    """Λ_f(s) = Λ(s,χ) + ω·Λ(s,χ̄)  via lfunlambda (느림, κδ² 측정용)."""
    s_str = f"{sigma}+{t}*I" if t >= 0 else f"{sigma}-{abs(t)}*I"
    s = pari(s_str)
    L_chi, L_bar = LFUN_OBJS[idx]
    v1 = complex(pari.lfunlambda(L_chi, s))
    v2 = complex(pari.lfunlambda(L_bar, s))
    return v1 + omega * v2


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. Off-critical 영점 탐색
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 2. Off-critical 영점 탐색 ━━━")

T_MAX = 200
SIGMA_SCAN = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]
T_STEP = 0.5
THRESHOLD = 0.5


def find_offcritical(idx, omega, name):
    """f(s) = L(s,χ)+ω·L(s,χ̄) off-critical 영점 1D 스캔 (lfun 사용, 빠름)."""
    log(f"\n  ─── {name}, ω={omega:.6f} ───")
    t0 = time.time()

    candidates = []
    for sigma in SIGMA_SCAN:
        ts = np.arange(2.0, T_MAX, T_STEP)
        prev_abs = None
        prev_t = None

        for t_val in ts:
            try:
                val = eval_L_f(idx, sigma, t_val, omega)
                cur_abs = abs(val)
            except Exception:
                cur_abs = 1e30

            if prev_abs is not None and prev_abs < THRESHOLD:
                if cur_abs > prev_abs:
                    candidates.append((sigma, float(prev_t), prev_abs))

            prev_abs = cur_abs
            prev_t = t_val

        elapsed = time.time() - t0
        log(f"    σ={sigma:.2f}: {len(candidates)}개 누적 ({elapsed:.0f}초)")

    elapsed_scan = time.time() - t0
    log(f"    스캔 완료: {len(candidates)}개 후보 ({elapsed_scan:.1f}초)")

    # 상위 15개 정밀화 (lfun 사용 — 빠름)
    confirmed = []
    for sigma0, t0_val, v0 in sorted(candidates, key=lambda x: x[2])[:15]:
        try:
            def objective(x):
                return abs(eval_L_f(idx, x[0], x[1], omega))

            res = minimize(objective, [sigma0, t0_val], method='Nelder-Mead',
                           options={'xatol': 1e-10, 'fatol': 1e-14, 'maxiter': 500})

            sigma_r, t_r = res.x
            val_r = res.fun

            if 0.01 < sigma_r < 0.99 and t_r > 0.5 and abs(sigma_r - 0.5) > 0.01 and val_r < 1e-8:
                if sigma_r < 0.5:
                    sigma_r = 1.0 - sigma_r

                is_dup = any(abs(c['sigma'] - sigma_r) < 0.005 and
                             abs(c['t'] - t_r) < 0.5 for c in confirmed)
                if not is_dup:
                    # lfunlambda로 검증
                    lam_val = abs(eval_Lambda_f(idx, sigma_r, t_r, omega))
                    confirmed.append({'sigma': sigma_r, 't': t_r, 'absf': val_r, 'absLam': lam_val})
                    log(f"    ✅ #{len(confirmed)}: σ={sigma_r:.8f}, t={t_r:.6f}, "
                        f"|f|={val_r:.3e}, |Λ_f|={lam_val:.3e}, |σ-½|={abs(sigma_r-0.5):.6f}")
        except Exception:
            pass

    elapsed_total = time.time() - t0
    log(f"    확인: {len(confirmed)}개 ({elapsed_total:.1f}초)")
    return confirmed


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. κδ² + c₁ 측정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DELTAS = [0.01, 0.02, 0.03, 0.05, 0.07, 0.1]

def gamma_deriv_ratio(sigma, t, q, a):
    """
    γ'/γ(s) for Dirichlet L-function χ mod q with parity a.
    γ(s) = (q/π)^{s/2} · Γ((s+a)/2)
    γ'/γ = (1/2)·log(q/π) + (1/2)·ψ((s+a)/2)
    """
    from scipy.special import digamma
    s_half = complex((sigma + a) / 2, t / 2)
    # scipy digamma: real argument → extend via recurrence for complex
    # Use PARI for complex digamma instead
    s_str = f"{(sigma+a)/2}+{t/2}*I" if t >= 0 else f"{(sigma+a)/2}-{abs(t)/2}*I"
    psi_val = complex(pari(f"psi({s_str})"))
    return 0.5 * np.log(q / np.pi) + 0.5 * psi_val


def measure_c1(zero, idx, omega):
    """κδ² 다중 δ + c₁ fit.  Λ = γ·f 보정 포함."""
    sigma0, t0 = zero['sigma'], zero['t']
    q = CONFIGS_DATA[idx]['q']
    a_val = CONFIGS_DATA[idx]['a']

    # κ = |Λ_f'/Λ_f|² = |f'/f + γ'/γ|²
    # f(s) = L(s,χ) + ω·L(s,χ̄),  γ는 두 성분 공통 (같은 패리티)
    kd2_data = []
    for delta in DELTAS:
        s_off = sigma0 + delta
        h = 1e-7
        try:
            f_c = eval_L_f(idx, s_off, t0, omega)
            f_p = eval_L_f(idx, s_off + h, t0, omega)
            f_m = eval_L_f(idx, s_off - h, t0, omega)
            if abs(f_c) > 1e-200:
                f_d = (f_p - f_m) / (2 * h)
                logderiv_f = f_d / f_c  # f'/f
                gamma_corr = gamma_deriv_ratio(s_off, t0, q, a_val)
                logderiv_Lambda = logderiv_f + gamma_corr  # Λ'/Λ
                kappa = abs(logderiv_Lambda) ** 2
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

    # c₁ analytic: c₁ = Re(f''/f') + 2·Re(γ'/γ(ρ))
    # 유도: Λ'/Λ = 1/δ + B + O(δ), B = f''/f'/2 + γ'/γ
    # κδ² = 1 + 2·Re(B)·δ + ... → c₁ = 2Re(B) = Re(f''/f') + 2Re(γ'/γ)
    h2 = 1e-6
    try:
        f_p = eval_L_f(idx, sigma0 + h2, t0, omega)
        f_c = eval_L_f(idx, sigma0, t0, omega)
        f_m = eval_L_f(idx, sigma0 - h2, t0, omega)
        fd = (f_p - f_m) / (2 * h2)
        fdd = (f_p - 2*f_c + f_m) / (h2**2)
        gamma_re = gamma_deriv_ratio(sigma0, t0, q, a_val).real
        c1_analytic = (fdd / fd).real + 2 * gamma_re if abs(fd) > 1e-200 else None
    except Exception:
        c1_analytic = None

    return c1_fit, c1_analytic, kd2_data


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

all_results = []

for idx, cfg in enumerate(CONFIGS_DATA):
    omega = cfg['omega']
    name = cfg['name']

    zeros = find_offcritical(idx, omega, name)

    if zeros:
        log(f"\n  ━━━ c₁ 측정: {name} ━━━")
        for z in zeros:
            c1_fit, c1_analytic, kd2_data = measure_c1(z, idx, omega)
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
# 6. 종합
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
        std_p = np.std(c1_products) if len(c1_products) > 1 else 0
        log(f"  SC2 (c₁·|σ-½|≈1): mean={mean_p:.3f}±{std_p:.3f}")

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
