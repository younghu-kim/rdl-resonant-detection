#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 C-240] Laurent 계수 패리티 정리 수치 검증
=============================================================================
목적:
  Theorem (Laurent coefficient parity) 검증:
  FE + SR → cₙ = (-1)^(n+1) c̄ₙ at zeros on the critical line.
  - n 짝수: Re(cₙ) = 0 (순허수)
  - n 홀수: Im(cₙ) = 0 (순실수)

  기존: Re(c₀)=0, Im(c₁)=0 확인됨.
  **새 예측**: Re(c₂)=0, Im(c₃)=0 — 최초 수치 검증.

방법:
  Λ'/Λ(ρ+δ) = 1/δ + c₀ + c₁δ + c₂δ² + c₃δ³ + ...
  대칭/반대칭 분리 후 δ² 다항식 피팅으로 c₀-c₃동시 추출.

대상: ζ(s) 5영점, L(s,11a1) 3영점, L(s,χ₋₃) 3영점 = 총 11영점
결과: results/laurent_parity_c240.txt
=============================================================================
"""
import sys, os, time
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'laurent_parity_c240.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 C-240] Laurent 계수 패리티 정리 수치 검증")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()
log("정리: FE + SR → cₙ = (-1)^(n+1) c̄ₙ")
log("  n 짝수 → Re(cₙ) = 0 (순허수)")
log("  n 홀수 → Im(cₙ) = 0 (순실수)")
log("  기존 확인: Re(c₀)=0, Im(c₁)=0")
log("  **새 예측**: Re(c₂)=0, Im(c₃)=0")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PARI 초기화 — 200자리 정밀도 (c₂, c₃ 안정성)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
import cypari2
gp = cypari2.Pari()
gp.allocatemem(1000 * 1024 * 1024)
gp("default(realprecision, 200)")
log(f"{T()} PARI OK, realprecision=200")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# L-함수 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
LFUNCS = [
    {
        "name": "ζ(s)",
        "setup": "Lcur = lfuncreate(1)",
        "degree": 1,
        "n_zeros": 5,
        "t_max": 50.0,
    },
    {
        "name": "L(s,11a1)",
        "setup": 'E11 = ellinit([0,-1,1,-10,-20]); Lcur = lfuncreate(E11)',
        "degree": 2,
        "n_zeros": 3,
        "t_max": 20.0,
    },
    {
        "name": "L(s,χ₋₃)",
        "setup": 'Lcur = lfuncreate(Mod(2,3))',
        "degree": 1,
        "n_zeros": 3,
        "t_max": 30.0,
    },
]

# δ 범위: 0.001~0.05 (18개 — 다항식 피팅에 충분)
DELTAS = np.array([0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.007,
                    0.01, 0.012, 0.015, 0.02, 0.025, 0.03, 0.035,
                    0.04, 0.045, 0.05, 0.06])

SIGMA_C = 0.5  # 기본값, L-함수별 갱신

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def get_log_deriv(sigma, t0):
    """Λ'/Λ(σ+it₀) via lfunlambda"""
    gp(f"s_eval = {sigma:.20f} + I*{t0:.20f}")
    gp("Lv = lfunlambda(Lcur, s_eval)")
    gp("dLv = lfunlambda(Lcur, s_eval, 1)")
    gp("r_val = dLv / Lv")
    return complex(float(gp("real(r_val)")), float(gp("imag(r_val)")))


def extract_laurent_coefficients(t0, deltas):
    """
    c₀, c₁, c₂, c₃ 동시 추출 — δ² 다항식 피팅.

    대칭 부분: sym(δ) = [f(+δ)+f(-δ)]/2 = c₀ + c₂δ² + c₄δ⁴ + ...
    → δ²의 다항식으로 피팅 → c₀ (절편), c₂ (1차 계수)

    반대칭 부분: anti(δ) = [(f(+δ)-f(-δ))/2 - 1/δ]/δ = c₁ + c₃δ² + c₅δ⁴ + ...
    → δ²의 다항식으로 피팅 → c₁ (절편), c₃ (1차 계수)
    """
    sym_vals = []
    anti_vals = []
    good_deltas = []

    for d in deltas:
        try:
            f_plus = get_log_deriv(SIGMA_C + d, t0)
            f_minus = get_log_deriv(SIGMA_C - d, t0)

            sym = (f_plus + f_minus) / 2.0
            anti_raw = (f_plus - f_minus) / 2.0
            anti = (anti_raw - 1.0/d) / d  # c₁ + c₃δ² + ...

            sym_vals.append(sym)
            anti_vals.append(anti)
            good_deltas.append(d)
        except Exception as e:
            log(f"    ⚠️ δ={d}: {e}")

    if len(good_deltas) < 6:
        return None, None, None, None, False

    d_arr = np.array(good_deltas)
    d2 = d_arr ** 2

    # ━━ 대칭 부분: sym(δ) = c₀ + c₂δ² + c₄δ⁴ ━━
    # 2차 다항식 피팅 in δ²: sym = a + b*δ² + c*δ⁴
    sym_re = np.array([s.real for s in sym_vals])
    sym_im = np.array([s.imag for s in sym_vals])

    # Re(sym) = Re(c₀) + Re(c₂)*δ² + Re(c₄)*δ⁴
    p_sym_re = np.polyfit(d2, sym_re, 2)  # [c₄_coeff, c₂_coeff, c₀_coeff]
    c0_re = p_sym_re[2]
    c2_re = p_sym_re[1]

    # Im(sym) = Im(c₀) + Im(c₂)*δ² + Im(c₄)*δ⁴
    p_sym_im = np.polyfit(d2, sym_im, 2)
    c0_im = p_sym_im[2]
    c2_im = p_sym_im[1]

    c0 = complex(c0_re, c0_im)
    c2 = complex(c2_re, c2_im)

    # ━━ 반대칭 부분: anti(δ) = c₁ + c₃δ² + c₅δ⁴ ━━
    anti_re = np.array([a.real for a in anti_vals])
    anti_im = np.array([a.imag for a in anti_vals])

    p_anti_re = np.polyfit(d2, anti_re, 2)
    c1_re = p_anti_re[2]
    c3_re = p_anti_re[1]

    p_anti_im = np.polyfit(d2, anti_im, 2)
    c1_im = p_anti_im[2]
    c3_im = p_anti_im[1]

    c1 = complex(c1_re, c1_im)
    c3 = complex(c3_re, c3_im)

    # 피팅 품질: R² 확인
    sym_pred_re = np.polyval(p_sym_re, d2)
    ss_res = np.sum((sym_re - sym_pred_re)**2)
    ss_tot = np.sum((sym_re - np.mean(sym_re))**2)
    r2_sym = 1 - ss_res/ss_tot if ss_tot > 1e-30 else 0

    anti_pred_re = np.polyval(p_anti_re, d2)
    ss_res_a = np.sum((anti_re - anti_pred_re)**2)
    ss_tot_a = np.sum((anti_re - np.mean(anti_re))**2)
    r2_anti = 1 - ss_res_a/ss_tot_a if ss_tot_a > 1e-30 else 0

    return c0, c1, c2, c3, True


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

all_results = []
total_pass = 0
total_tested = 0

for lf in LFUNCS:
    log(f"{T()} ━━━ {lf['name']} (d={lf['degree']}) ━━━")

    try:
        setup_cmd = lf['setup'].strip().replace('\n', '; ')
        gp(setup_cmd)

        k_val = float(str(gp("Lcur[4]")))
        SIGMA_C = k_val / 2.0
        log(f"  k = {k_val}, σ_c = {SIGMA_C}")

        gp(f"zeros = lfunzeros(Lcur, {lf['t_max']:.1f})")
        n_zeros_found = int(gp("length(zeros)"))

        if n_zeros_found == 0:
            log(f"  ⚠️ 영점 없음")
            continue

        zeros_to_use = min(lf['n_zeros'], n_zeros_found)
        zero_list = []
        for j in range(1, zeros_to_use + 1):
            zt = float(gp(f"zeros[{j}]"))
            zero_list.append(zt)

        log(f"  영점 {n_zeros_found}개, 사용: {[f'{z:.4f}' for z in zero_list]}")
        log()

        for idx, t0 in enumerate(zero_list, 1):
            log(f"  [영점 #{idx}] t₀ = {t0:.6f}")

            c0, c1, c2, c3, ok = extract_laurent_coefficients(t0, DELTAS)

            if not ok:
                log(f"    ⚠️ 추출 실패")
                continue

            total_tested += 1

            # 패리티 판정
            log(f"    c₀ = ({c0.real:.6e}) + ({c0.imag:.6f})i")
            log(f"    |Re(c₀)| = {abs(c0.real):.3e}  ← 예측: 0  {'✓' if abs(c0.real) < 1e-6 else '✗'}")
            log(f"    c₁ = ({c1.real:.6f}) + ({c1.imag:.6e})i")
            log(f"    |Im(c₁)| = {abs(c1.imag):.3e}  ← 예측: 0  {'✓' if abs(c1.imag) < 1e-6 else '✗'}")
            log(f"    c₂ = ({c2.real:.6e}) + ({c2.imag:.6f})i")
            log(f"    |Re(c₂)| = {abs(c2.real):.3e}  ← **새 예측: 0**  {'✓' if abs(c2.real) < 1e-3 else '✗'}")
            log(f"    c₃ = ({c3.real:.6f}) + ({c3.imag:.6e})i")
            log(f"    |Im(c₃)| = {abs(c3.imag):.3e}  ← **새 예측: 0**  {'✓' if abs(c3.imag) < 1e-3 else '✗'}")

            # 상대 오차 계산 (0이 아닌 부분 기준)
            # c₀: |Re(c₀)| / |Im(c₀)| (c₀는 순허수)
            rel_c0 = abs(c0.real) / abs(c0.imag) if abs(c0.imag) > 1e-30 else float('inf')
            # c₁: |Im(c₁)| / |Re(c₁)| (c₁는 순실수)
            rel_c1 = abs(c1.imag) / abs(c1.real) if abs(c1.real) > 1e-30 else float('inf')
            # c₂: |Re(c₂)| / |Im(c₂)| (c₂는 순허수)
            rel_c2 = abs(c2.real) / abs(c2.imag) if abs(c2.imag) > 1e-30 else float('inf')
            # c₃: |Im(c₃)| / |Re(c₃)| (c₃는 순실수)
            rel_c3 = abs(c3.imag) / abs(c3.real) if abs(c3.real) > 1e-30 else float('inf')

            log(f"    상대 오차: c₀={rel_c0:.2e}, c₁={rel_c1:.2e}, c₂={rel_c2:.2e}, c₃={rel_c3:.2e}")

            # 판정: 상대 오차 < 1e-3 이면 PASS (0.1% 이내)
            # c₂, c₃는 δ⁴ 피팅 오차가 있으므로 1% 기준
            c0_pass = rel_c0 < 1e-3
            c1_pass = rel_c1 < 1e-3
            c2_pass = rel_c2 < 1e-2
            c3_pass = rel_c3 < 1e-2
            all_pass = c0_pass and c1_pass and c2_pass and c3_pass

            verdict = "★★★ PASS" if all_pass else "✗ FAIL"
            log(f"    판정: {verdict}")
            log()

            if all_pass:
                total_pass += 1

            all_results.append({
                'lfunc': lf['name'],
                'degree': lf['degree'],
                't0': t0,
                'c0': c0, 'c1': c1, 'c2': c2, 'c3': c3,
                'rel_c0': rel_c0, 'rel_c1': rel_c1,
                'rel_c2': rel_c2, 'rel_c3': rel_c3,
                'pass': all_pass,
            })

    except Exception as e:
        log(f"  ❌ 에러: {e}")
        import traceback
        traceback.print_exc()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 통합 요약
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log()
log("=" * 72)
log("통합 요약")
log("=" * 72)
log()

if all_results:
    # 계수별 통계
    c0_rels = [r['rel_c0'] for r in all_results]
    c1_rels = [r['rel_c1'] for r in all_results]
    c2_rels = [r['rel_c2'] for r in all_results]
    c3_rels = [r['rel_c3'] for r in all_results]

    log("계수별 |vanishing ratio| 통계:")
    log(f"  |Re(c₀)|/|Im(c₀)|: mean={np.mean(c0_rels):.2e}, max={np.max(c0_rels):.2e}  [기존 확인]")
    log(f"  |Im(c₁)|/|Re(c₁)|: mean={np.mean(c1_rels):.2e}, max={np.max(c1_rels):.2e}  [기존 확인]")
    log(f"  |Re(c₂)|/|Im(c₂)|: mean={np.mean(c2_rels):.2e}, max={np.max(c2_rels):.2e}  [**새 검증**]")
    log(f"  |Im(c₃)|/|Re(c₃)|: mean={np.mean(c3_rels):.2e}, max={np.max(c3_rels):.2e}  [**새 검증**]")
    log()

    # 영점별 테이블
    log("영점별 결과:")
    log(f"{'L-func':<15} {'d':>2} {'t₀':>10} {'|Re(c₀)|':>10} {'|Im(c₁)|':>10} {'|Re(c₂)|':>10} {'|Im(c₃)|':>10} {'판정':>6}")
    log("-" * 80)
    for r in all_results:
        log(f"{r['lfunc']:<15} {r['degree']:>2} {r['t0']:>10.4f} "
            f"{abs(r['c0'].real):>10.2e} {abs(r['c1'].imag):>10.2e} "
            f"{abs(r['c2'].real):>10.2e} {abs(r['c3'].imag):>10.2e} "
            f"{'PASS' if r['pass'] else 'FAIL':>6}")
    log()

    log(f"총 판정: {total_pass}/{total_tested} PASS")
    log()

    if total_pass == total_tested:
        log("★★★★ 정리 수치 확립: cₙ = (-1)^(n+1) c̄ₙ, n=0,1,2,3 모두 검증됨")
        log("  Re(c₀)=0, Im(c₁)=0: 기존 결과 재확인")
        log("  Re(c₂)=0, Im(c₃)=0: **새로운 수치 검증 (최초)**")
    else:
        log("⚠️ 일부 실패 — 피팅 정밀도 점검 필요")
else:
    log("결과 없음")

log()
log(f"총 소요: {time.time()-START:.1f}s")
log("=" * 72)
outf.close()
