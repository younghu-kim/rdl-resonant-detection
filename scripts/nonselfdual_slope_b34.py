#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 B-34] 비자기쌍대 L-함수의 κδ² slope 검증
=============================================================================
목적:
  Slope Universality Theorem (#214)의 sharpness 검증.
  정리 조건: FE + RC (자기쌍대). RC 없이(비자기쌍대) slope=2가 깨지는지 확인.

이론적 예측:
  - 자기쌍대 (RC 성립): FE+RC ⟹ Re(c₀)=0 ⟹ δ¹ 소멸 ⟹ slope=2
  - 비자기쌍대 (RC 불성립): Re(c₀)≠0 가능 ⟹ δ¹ 잔존 ⟹ slope=1
    (κδ²-1 ≈ 2·Re(c₀)·δ for small δ)

실험 설계:
  A. 비자기쌍대 (실험군):
     - χ mod 5, chi=[1]: χ(2)=i, order 4, 복소
     - χ mod 7, chi=[1]: χ(3)=e^(πi/3), order 6, 복소
     - χ mod 8, chi=[1,1]: 비원시적 → mod 8 generator order 2 → 실수...
       → 대신 χ mod 13, chi=[1]: order 12, 복소
  B. 자기쌍대 (대조군):
     - χ₋₃ (Kronecker, 실수): slope=2 확인됨 (#214)
     - χ₋₇ (Kronecker, 실수): slope=2 확인됨 (#214)

각 L-함수에서 5영점, 총 25영점 측정.

결과: results/nonselfdual_slope_b34.txt
=============================================================================
"""
import sys, os, time
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'nonselfdual_slope_b34.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 B-34] 비자기쌍대 L-함수의 κδ² slope 검증")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PARI 초기화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import cypari2
gp = cypari2.Pari()
gp.allocatemem(1000 * 1024 * 1024)
gp("default(realprecision, 150)")
log(f"{T()} PARI OK, realprecision=150")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 이론 서술
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("PREDICTION (Non-Self-Dual Slope Deviation)")
log("=" * 72)
log()
log("For self-dual Λ(s): FE + RC ⟹ Re(c₀) = 0 ⟹ slope = 2 (Thm #214)")
log()
log("For non-self-dual Λ(s,χ) with χ ≠ χ̄:")
log("  FE: Λ(s,χ) = ε·Λ(1-s,χ̄)  → Λ'/Λ(s,χ) = -Λ'/Λ(1-s,χ̄)")
log("  RC violated: Λ(s̄,χ) = Λ̄(s,χ̄) ≠ Λ̄(s,χ) in general")
log()
log("  Laurent at ρ=½+iγ (zero of χ):")
log("    Λ'/Λ(½+δ+iγ,χ) = 1/δ + c₀ + c₁δ + ...")
log("  FE gives: c₀(χ,ρ) = -c₀(χ̄,ρ̄)")
log("  But c₀(χ̄,ρ̄) is a Laurent coeff of a DIFFERENT L-function at")
log("  a DIFFERENT point → no constraint on Re(c₀(χ,ρ))")
log()
log("  If Re(c₀) ≠ 0:")
log("    κδ² - 1 = 2·Re(c₀)·δ + O(δ²)")
log("    log(κδ²-1) ≈ log(2|Re(c₀)|) + 1·log(δ)")
log("    ⟹ slope = 1 (NOT 2)")
log()
log("  If Re(c₀) = 0 (despite RC violation):")
log("    slope = 2 → some hidden symmetry at play")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# L-함수 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

LFUNCS = [
    # 비자기쌍대 (실험군)
    {
        "name": "L(s,χ₅ complex)",
        "setup": "G5 = znstar(5, 1); Lcur = lfuncreate([G5, [1]])",
        "self_dual": False,
        "q": 5,
        "chi_label": "χ mod 5, order 4 (complex)",
    },
    {
        "name": "L(s,χ₇ complex)",
        "setup": "G7 = znstar(7, 1); Lcur = lfuncreate([G7, [1]])",
        "self_dual": False,
        "q": 7,
        "chi_label": "χ mod 7, order 6 (complex)",
    },
    {
        "name": "L(s,χ₁₃ complex)",
        "setup": "G13 = znstar(13, 1); Lcur = lfuncreate([G13, [1]])",
        "self_dual": False,
        "q": 13,
        "chi_label": "χ mod 13, order 12 (complex)",
    },
    # 자기쌍대 (대조군 — #214에서 slope=2 확인됨)
    {
        "name": "L(s,χ₋₃)",
        "setup": "Lcur = lfuncreate(-3)",
        "self_dual": True,
        "q": 3,
        "chi_label": "Kronecker (-3) (real, self-dual)",
    },
    {
        "name": "L(s,χ₋₇)",
        "setup": "Lcur = lfuncreate(-7)",
        "self_dual": True,
        "q": 7,
        "chi_label": "Kronecker (-7) (real, self-dual)",
    },
]

N_ZEROS_PER = 5
C0_DELTAS = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01])
SLOPE_DELTAS = np.array([0.01, 0.015, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2])
T_MAX = 50.0

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def get_log_deriv(sigma, t0):
    """Λ'/Λ(σ+it₀) via lfunlambda"""
    gp(f"s_eval = {sigma:.20f} + I*{t0:.20f}")
    gp("Lv = lfunlambda(Linit, s_eval)")
    gp("dLv = lfunlambda(Linit, s_eval, 1)")
    gp("r_val = dLv / Lv")
    return complex(float(gp("real(r_val)")), float(gp("imag(r_val)")))

def extract_c0_symmetric(t0):
    """c₀ 추출 — (f(½+δ) + f(½-δ))/2 ≈ c₀ + O(δ²)"""
    c0_estimates = []
    for d in C0_DELTAS:
        f_plus = get_log_deriv(0.5 + d, t0)
        f_minus = get_log_deriv(0.5 - d, t0)
        c0_est = (f_plus + f_minus) / 2.0
        c0_estimates.append(c0_est)
    c0_val = np.mean(c0_estimates[:3])
    c0_spread_re = np.std([v.real for v in c0_estimates[:3]])
    return c0_val, c0_spread_re

def extract_c1_symmetric(t0):
    """c₁ 추출 — [(f(½+δ) - f(½-δ))/2 - 1/δ] / δ ≈ c₁ + O(δ²)"""
    c1_estimates = []
    for d in C0_DELTAS:
        f_plus = get_log_deriv(0.5 + d, t0)
        f_minus = get_log_deriv(0.5 - d, t0)
        antisym = (f_plus - f_minus) / 2.0
        c1_est = (antisym - 1.0/d) / d
        c1_estimates.append(c1_est)
    c1_val = np.mean(c1_estimates[:3])
    return c1_val

def measure_slope(t0):
    """κδ² log-log slope 측정 — 두 가지 모델 피팅"""
    kd2_vals = []
    for d in SLOPE_DELTAS:
        ratio = get_log_deriv(0.5 + d, t0)
        kappa = abs(ratio)**2
        kd2_vals.append(kappa * d**2)

    # slope=2 모델: log(κδ²-1) vs log(δ) 전체 피팅
    valid_d = []
    valid_lkd2 = []
    for d, kd2 in zip(SLOPE_DELTAS, kd2_vals):
        if kd2 > 1.001:
            valid_d.append(np.log(d))
            valid_lkd2.append(np.log(kd2 - 1.0))

    if len(valid_d) < 4:
        return None, None, None, None, None

    # 전체 범위 피팅
    coeffs_full = np.polyfit(valid_d, valid_lkd2, 1)
    slope_full = coeffs_full[0]
    pred_full = np.polyval(coeffs_full, valid_d)
    ss_res = np.sum((np.array(valid_lkd2) - pred_full)**2)
    ss_tot = np.sum((np.array(valid_lkd2) - np.mean(valid_lkd2))**2)
    r2_full = 1.0 - ss_res/ss_tot if ss_tot > 1e-30 else 0.0

    # 소-δ 범위 피팅 (δ ≤ 0.05만): slope=1 vs slope=2 구분에 더 민감
    small_d = []
    small_lkd2 = []
    for d, kd2 in zip(SLOPE_DELTAS, kd2_vals):
        if d <= 0.05 and kd2 > 1.001:
            small_d.append(np.log(d))
            small_lkd2.append(np.log(kd2 - 1.0))

    slope_small = None
    r2_small = None
    if len(small_d) >= 3:
        coeffs_small = np.polyfit(small_d, small_lkd2, 1)
        slope_small = coeffs_small[0]
        pred_small = np.polyval(coeffs_small, small_d)
        ss_res_s = np.sum((np.array(small_lkd2) - pred_small)**2)
        ss_tot_s = np.sum((np.array(small_lkd2) - np.mean(small_lkd2))**2)
        r2_small = 1.0 - ss_res_s/ss_tot_s if ss_tot_s > 1e-30 else 0.0

    A_meas = np.exp(coeffs_full[1]) if slope_full else None
    return slope_full, r2_full, slope_small, r2_small, A_meas

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("수치 실험")
log("=" * 72)
log()

all_results = []

for lf in LFUNCS:
    log(f"{T()} ━━━ {lf['name']} ━━━")
    log(f"  {lf['chi_label']}")
    log(f"  self-dual: {lf['self_dual']}")

    try:
        gp(lf["setup"])
    except Exception as e:
        log(f"  ⚠ lfuncreate 실패: {e}")
        log()
        continue

    fe = float(gp("lfuncheckfeq(Lcur)"))
    log(f"  FE: {fe:.0f} digits")

    gp(f"Linit = lfuninit(Lcur, [1.0, {T_MAX + 20}])")

    gp(f"zvec = lfunzeros(Linit, {T_MAX})")
    n_found = int(gp("length(zvec)"))
    n_use = min(n_found, N_ZEROS_PER)
    zeros = [float(gp(f"zvec[{i}]")) for i in range(1, n_use + 1)]
    log(f"  영점 {n_found}개, 사용: {[f'{z:.4f}' for z in zeros]}")
    log()

    for idx, t0 in enumerate(zeros):
        log(f"  [영점 #{idx+1}] t₀ = {t0:.6f}")

        # c₀
        c0, c0_spread = extract_c0_symmetric(t0)
        log(f"    c₀ = ({c0.real:.6e}) + ({c0.imag:.6f})i")
        log(f"    |Re(c₀)| = {abs(c0.real):.6e}  [self-dual 이론: 0]")

        # c₁
        c1 = extract_c1_symmetric(t0)
        log(f"    c₁ = ({c1.real:.6f}) + ({c1.imag:.6e})i")
        log(f"    |Im(c₁)| = {abs(c1.imag):.6e}  [self-dual 이론: 0]")

        # A (이론 예측 — Re(c₀)=0이면)
        A_pred = c0.imag**2 + 2*c1.real
        log(f"    A_pred = Im(c₀)²+2Re(c₁) = {c0.imag**2:.4f}+{2*c1.real:.4f} = {A_pred:.4f}")

        # slope
        slope_full, r2_full, slope_small, r2_small, A_meas = measure_slope(t0)
        if slope_full is not None:
            log(f"    slope(전체) = {slope_full:.4f}  R² = {r2_full:.8f}")
            if slope_small is not None:
                log(f"    slope(δ≤0.05) = {slope_small:.4f}  R² = {r2_small:.8f}")
            log(f"    A_meas = {A_meas:.4f}")
            if A_pred > 0:
                log(f"    A 오차 = {abs(A_pred - A_meas)/A_pred*100:.1f}%")
        else:
            log(f"    slope: 측정 불가 (유효 데이터 부족)")

        # 이론 판정
        re_c0_sig = abs(c0.real) > max(c0_spread * 3, 1e-10)
        if lf['self_dual']:
            expect = "Re(c₀)=0, slope=2"
            if abs(c0.real) < 1e-6 and slope_full and abs(slope_full - 2.0) < 0.05:
                verdict = "★★★ PASS (대조군 확인)"
            else:
                verdict = "✗ FAIL (대조군 이상)"
        else:
            expect = "Re(c₀)≠0 가능, slope≠2 가능"
            if re_c0_sig:
                if slope_small and abs(slope_small - 1.0) < 0.3:
                    verdict = f"★★★ slope≈1 확인 (Re(c₀)={c0.real:.4e})"
                elif slope_full and abs(slope_full - 2.0) > 0.1:
                    verdict = f"★★ slope≠2 확인 (slope={slope_full:.3f})"
                else:
                    verdict = f"△ Re(c₀)≠0이나 slope 모호"
            else:
                if slope_full and abs(slope_full - 2.0) < 0.05:
                    verdict = "△ 비자기쌍대이나 Re(c₀)≈0 — 추가 대칭?"
                else:
                    verdict = f"? 예상 밖 패턴 (Re(c₀)≈0, slope={slope_full:.3f})"

        log(f"    기대: {expect}")
        log(f"    판정: {verdict}")
        log()

        all_results.append({
            'name': lf['name'],
            'self_dual': lf['self_dual'],
            't0': t0,
            'Re_c0': c0.real,
            'Im_c0': c0.imag,
            'Re_c1': c1.real,
            'Im_c1': c1.imag,
            'slope_full': slope_full,
            'r2_full': r2_full,
            'slope_small': slope_small,
            'r2_small': r2_small,
            'A_pred': A_pred,
            'A_meas': A_meas,
        })

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 종합 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("종합 분석")
log("=" * 72)
log()

# 자기쌍대 vs 비자기쌍대 통계
sd_results = [r for r in all_results if r['self_dual']]
nsd_results = [r for r in all_results if not r['self_dual']]

if sd_results:
    sd_re_c0 = [abs(r['Re_c0']) for r in sd_results]
    sd_slopes = [r['slope_full'] for r in sd_results if r['slope_full'] is not None]
    log("━━━ 자기쌍대 (대조군) ━━━")
    log(f"  N = {len(sd_results)} 영점")
    log(f"  max|Re(c₀)| = {max(sd_re_c0):.3e}")
    log(f"  mean|Re(c₀)| = {np.mean(sd_re_c0):.3e}")
    if sd_slopes:
        log(f"  mean slope = {np.mean(sd_slopes):.6f} ± {np.std(sd_slopes):.6f}")
    log(f"  판정: {'★★★ PASS (slope=2 대조 확인)' if sd_slopes and abs(np.mean(sd_slopes)-2.0)<0.01 else '확인 필요'}")
    log()

if nsd_results:
    nsd_re_c0 = [abs(r['Re_c0']) for r in nsd_results]
    nsd_slopes_full = [r['slope_full'] for r in nsd_results if r['slope_full'] is not None]
    nsd_slopes_small = [r['slope_small'] for r in nsd_results if r['slope_small'] is not None]
    log("━━━ 비자기쌍대 (실험군) ━━━")
    log(f"  N = {len(nsd_results)} 영점")
    log(f"  max|Re(c₀)| = {max(nsd_re_c0):.6e}")
    log(f"  mean|Re(c₀)| = {np.mean(nsd_re_c0):.6e}")
    log(f"  min|Re(c₀)| = {min(nsd_re_c0):.6e}")
    if nsd_slopes_full:
        log(f"  mean slope(전체) = {np.mean(nsd_slopes_full):.4f} ± {np.std(nsd_slopes_full):.4f}")
    if nsd_slopes_small:
        log(f"  mean slope(δ≤0.05) = {np.mean(nsd_slopes_small):.4f} ± {np.std(nsd_slopes_small):.4f}")
    log()

    # 핵심 질문: Re(c₀)가 실제로 ≠ 0인가?
    n_sig = sum(1 for r in nsd_results if abs(r['Re_c0']) > 1e-6)
    log(f"  |Re(c₀)| > 1e-6: {n_sig}/{len(nsd_results)} 영점")
    log()

    if n_sig > 0:
        log("  결론: RC 조건이 진정 필요. 비자기쌍대에서 Re(c₀)≠0 확인.")
        if nsd_slopes_full and abs(np.mean(nsd_slopes_full) - 2.0) > 0.1:
            log(f"  slope = {np.mean(nsd_slopes_full):.3f} ≠ 2 → Theorem sharpness 입증!")
        elif nsd_slopes_small and abs(np.mean(nsd_slopes_small) - 1.0) < 0.3:
            log(f"  slope(소-δ) ≈ {np.mean(nsd_slopes_small):.3f} → 1에 가까움 (δ¹ 잔존)")
    else:
        log("  결론: 비자기쌍대에서도 Re(c₀)≈0 — 추가 대칭 존재 가능?")
        log("  (임계선 위의 대칭이 RC보다 더 근본적일 수 있음)")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 상세 결과표
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("━━━ 상세 결과표 ━━━")
hdr = f"{'L-함수':25s} {'SD':3s} {'t₀':>8s} {'Re(c₀)':>12s} {'Im(c₀)':>10s} {'Re(c₁)':>9s} {'Im(c₁)':>10s} {'slope':>7s} {'slope_sm':>9s} {'A_pred':>7s} {'A_meas':>7s}"
log(hdr)
log("-" * len(hdr))
for r in all_results:
    sd_str = "✓" if r['self_dual'] else "✗"
    sl_str = f"{r['slope_full']:.4f}" if r['slope_full'] is not None else "N/A"
    sl_sm_str = f"{r['slope_small']:.4f}" if r['slope_small'] is not None else "N/A"
    am_str = f"{r['A_meas']:.4f}" if r['A_meas'] is not None else "N/A"
    log(f"{r['name']:25s} {sd_str:3s} {r['t0']:8.3f} {r['Re_c0']:12.6e} {r['Im_c0']:10.4f} {r['Re_c1']:9.4f} {r['Im_c1']:10.3e} {sl_str:>7s} {sl_sm_str:>9s} {r['A_pred']:7.4f} {am_str:>7s}")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 최종 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("최종 판정")
log("=" * 72)
log()

if nsd_results:
    nsd_re_c0_vals = [abs(r['Re_c0']) for r in nsd_results]
    nsd_slopes_all = [r['slope_full'] for r in nsd_results if r['slope_full'] is not None]

    if np.mean(nsd_re_c0_vals) > 1e-6 and nsd_slopes_all:
        mean_slope = np.mean(nsd_slopes_all)
        if abs(mean_slope - 2.0) > 0.1:
            log("★★★ 강양성 — Theorem sharpness 입증")
            log(f"  비자기쌍대: mean|Re(c₀)| = {np.mean(nsd_re_c0_vals):.4e}, slope = {mean_slope:.3f}")
            log("  RC(자기쌍대) 조건은 Slope Universality Theorem에 필수.")
            log("  비자기쌍대에서 δ¹ 항이 잔존하여 slope≠2.")
        elif abs(mean_slope - 1.0) < 0.3:
            log("★★★ 강양성 — slope≈1 확인 (δ¹ 지배)")
            log(f"  비자기쌍대: Re(c₀)≠0 → κδ²-1 ~ 2Re(c₀)δ → slope=1")
        else:
            log(f"★★ 양성 — Re(c₀)≠0 확인, slope={mean_slope:.3f} (2와 구별됨)")
    elif np.mean(nsd_re_c0_vals) <= 1e-6:
        log("△ 중립/음성 — 비자기쌍대에서도 Re(c₀)≈0")
        log("  RC 조건은 slope=2에 불필요할 수 있음. 추가 조사 필요.")
        log("  가능한 설명: 임계선 위 대칭, 또는 수치 정밀도 한계.")
    else:
        log("△ 중립 — 결과 모호. 추가 분석 필요.")

log()
log(f"총 소요시간: {time.time()-START:.1f}초")
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

outf.close()
