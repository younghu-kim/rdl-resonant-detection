#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 C-241] 비자기쌍대 Laurent 패리티 sharpness 테스트
=============================================================================
목적:
  C-240에서 확립된 Laurent 패리티 정리의 sharpness 검증.
  정리: FE + SR → cₙ = (-1)^(n+1) c̄ₙ  (자기쌍대 L-함수에서)

  비자기쌍대(non-self-dual) L-함수에서는 SR이 성립하지 않으므로
  패리티가 **깨져야** 한다. 이를 확인하면 정리가 sharp (self-dual이 필요조건).

이론적 배경:
  자기쌍대: Λ(s,χ) = ε·Λ(1-s,χ)  AND  Λ(s̄,χ) = Λ̄(s,χ)
    → FE+SR → cₙ(ρ) = (-1)^(n+1) c̄ₙ(ρ)
    → n 짝수: Re(cₙ)=0, n 홀수: Im(cₙ)=0

  비자기쌍대: Λ(s,χ) = ε·Λ(1-s,χ̄)  where χ̄≠χ
    → SR: Λ(s̄,χ) = Λ̄(s,χ̄) ≠ Λ̄(s,χ)
    → FE는 χ와 χ̄의 Laurent 계수를 연결하지만 같은 L-함수 내 제약 없음
    → Re(c₀)≠0, Im(c₁)≠0 가능 → 패리티 깨짐

실험군: χ₅ (mod 5, order 4), χ₇ (mod 7, order 6) — 각 3영점
대조군: χ₋₃ (Kronecker, self-dual) — 3영점
결과: results/nonselfdual_parity_c241.txt
=============================================================================
"""
import sys, os, time
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'nonselfdual_parity_c241.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 C-241] 비자기쌍대 Laurent 패리티 sharpness 테스트")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()
log("정리 (C-240 확립): FE + SR → cₙ = (-1)^(n+1) c̄ₙ")
log("  자기쌍대: n 짝수 → Re(cₙ)=0, n 홀수 → Im(cₙ)=0")
log("  비자기쌍대: SR 실패 → 패리티 깨짐 예측")
log()
log("성공 기준:")
log("  대조군 (χ₋₃): |ratio| < 1e-5 → PASS (패리티 성립)")
log("  실험군 (χ₅, χ₇): |ratio| > 0.01 → BROKEN (패리티 깨짐 = sharpness)")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PARI 초기화
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
    # 비자기쌍대 (실험군)
    {
        "name": "L(s,χ₅)",
        "setup": "G5 = znstar(5, 1); Lcur = lfuncreate([G5, [1]])",
        "self_dual": False,
        "chi_label": "χ mod 5, order 4 (complex)",
        "n_zeros": 3,
        "t_max": 50.0,
    },
    {
        "name": "L(s,χ₇)",
        "setup": "G7 = znstar(7, 1); Lcur = lfuncreate([G7, [1]])",
        "self_dual": False,
        "chi_label": "χ mod 7, order 6 (complex)",
        "n_zeros": 3,
        "t_max": 50.0,
    },
    # 자기쌍대 (대조군)
    {
        "name": "L(s,χ₋₃)",
        "setup": "Lcur = lfuncreate(Mod(2,3))",
        "self_dual": True,
        "chi_label": "Kronecker (-3) (real, self-dual)",
        "n_zeros": 3,
        "t_max": 30.0,
    },
]

# δ 범위: C-240과 동일 (18개, 다항식 피팅용)
DELTAS = np.array([0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.007,
                    0.01, 0.012, 0.015, 0.02, 0.025, 0.03, 0.035,
                    0.04, 0.045, 0.05, 0.06])

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def get_log_deriv(sigma, t0):
    """Λ'/Λ(σ+it₀) via lfunlambda — Linit 사용"""
    gp(f"s_eval = {sigma:.20f} + I*{t0:.20f}")
    gp("Lv = lfunlambda(Linit, s_eval)")
    gp("dLv = lfunlambda(Linit, s_eval, 1)")
    gp("r_val = dLv / Lv")
    return complex(float(gp("real(r_val)")), float(gp("imag(r_val)")))


def extract_laurent_coefficients(t0, deltas, sigma_c=0.5):
    """
    c₀, c₁, c₂, c₃ 동시 추출 — δ² 다항식 피팅 (C-240과 동일).
    """
    sym_vals = []
    anti_vals = []
    good_deltas = []

    for d in deltas:
        try:
            f_plus = get_log_deriv(sigma_c + d, t0)
            f_minus = get_log_deriv(sigma_c - d, t0)

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

    # 대칭 부분: sym = c₀ + c₂δ² + c₄δ⁴
    sym_re = np.array([s.real for s in sym_vals])
    sym_im = np.array([s.imag for s in sym_vals])

    p_sym_re = np.polyfit(d2, sym_re, 2)
    c0_re = p_sym_re[2]
    c2_re = p_sym_re[1]

    p_sym_im = np.polyfit(d2, sym_im, 2)
    c0_im = p_sym_im[2]
    c2_im = p_sym_im[1]

    c0 = complex(c0_re, c0_im)
    c2 = complex(c2_re, c2_im)

    # 반대칭 부분: anti = c₁ + c₃δ² + c₅δ⁴
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

    return c0, c1, c2, c3, True


def parity_ratios(c0, c1, c2, c3):
    """패리티 비율 계산."""
    def ratio(vanishing, nonvanishing):
        if abs(nonvanishing) < 1e-30:
            return float('inf')
        return abs(vanishing) / abs(nonvanishing)

    r0 = ratio(c0.real, c0.imag)   # Re(c₀)/Im(c₀)
    r1 = ratio(c1.imag, c1.real)   # Im(c₁)/Re(c₁)
    r2 = ratio(c2.real, c2.imag)   # Re(c₂)/Im(c₂)
    r3 = ratio(c3.imag, c3.real)   # Im(c₃)/Re(c₃)
    return r0, r1, r2, r3


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

all_results = []
sd_results = []    # 자기쌍대
nsd_results = []   # 비자기쌍대

for lf in LFUNCS:
    log(f"{T()} ━━━ {lf['name']} ({lf['chi_label']}) ━━━")
    log(f"  자기쌍대: {'예' if lf['self_dual'] else '아니오'}")

    try:
        setup_cmd = lf['setup'].strip().replace('\n', '; ')
        gp(setup_cmd)

        # k (weight) 추출
        k_val = float(str(gp("Lcur[4]")))
        sigma_c = k_val / 2.0
        log(f"  k = {k_val}, σ_c = {sigma_c}")

        # lfuninit
        gp(f"Linit = lfuninit(Lcur, [1.0, {lf['t_max'] + 20}])")

        # 영점 탐색
        gp(f"zeros = lfunzeros(Linit, {lf['t_max']:.1f})")
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

            c0, c1, c2, c3, ok = extract_laurent_coefficients(t0, DELTAS, sigma_c)

            if not ok:
                log(f"    ⚠️ 추출 실패")
                continue

            r0, r1, r2, r3 = parity_ratios(c0, c1, c2, c3)

            log(f"    c₀ = ({c0.real:.6e}) + ({c0.imag:.6e})i")
            log(f"    |Re(c₀)|/|Im(c₀)| = {r0:.3e}  {'← 0 예측' if lf['self_dual'] else '← O(1) 예측'}")
            log(f"    c₁ = ({c1.real:.6e}) + ({c1.imag:.6e})i")
            log(f"    |Im(c₁)|/|Re(c₁)| = {r1:.3e}  {'← 0 예측' if lf['self_dual'] else '← O(1) 예측'}")
            log(f"    c₂ = ({c2.real:.6e}) + ({c2.imag:.6e})i")
            log(f"    |Re(c₂)|/|Im(c₂)| = {r2:.3e}  {'← 0 예측' if lf['self_dual'] else '← O(1) 예측'}")
            log(f"    c₃ = ({c3.real:.6e}) + ({c3.imag:.6e})i")
            log(f"    |Im(c₃)|/|Re(c₃)| = {r3:.3e}  {'← 0 예측' if lf['self_dual'] else '← O(1) 예측'}")

            # 판정
            if lf['self_dual']:
                c0_ok = r0 < 1e-3
                c1_ok = r1 < 1e-3
                c2_ok = r2 < 1e-2
                c3_ok = r3 < 1e-2
                all_pass = c0_ok and c1_ok and c2_ok and c3_ok
                verdict = "PASS (패리티 성립)" if all_pass else "FAIL (대조군 이상!)"
            else:
                any_broken = (r0 > 0.01) or (r1 > 0.01) or (r2 > 0.01) or (r3 > 0.01)
                verdict = "BROKEN (패리티 깨짐 = sharpness)" if any_broken else "INTACT (패리티 유지 = 예상 밖)"

            log(f"    판정: {verdict}")
            log()

            result = {
                'lfunc': lf['name'],
                'self_dual': lf['self_dual'],
                'chi_label': lf['chi_label'],
                't0': t0,
                'c0': c0, 'c1': c1, 'c2': c2, 'c3': c3,
                'r0': r0, 'r1': r1, 'r2': r2, 'r3': r3,
                'verdict': verdict,
            }
            all_results.append(result)
            if lf['self_dual']:
                sd_results.append(result)
            else:
                nsd_results.append(result)

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
    # 대조군
    log("━━ 대조군: 자기쌍대 (패리티 성립 예측) ━━")
    if sd_results:
        log(f"{'L-func':<15} {'t₀':>10} {'|Re(c₀)|/|Im|':>15} {'|Im(c₁)|/|Re|':>15} {'|Re(c₂)|/|Im|':>15} {'|Im(c₃)|/|Re|':>15} {'판정':>8}")
        log("-" * 100)
        for r in sd_results:
            log(f"{r['lfunc']:<15} {r['t0']:>10.4f} {r['r0']:>15.2e} {r['r1']:>15.2e} {r['r2']:>15.2e} {r['r3']:>15.2e} {'PASS' if 'PASS' in r['verdict'] else 'FAIL':>8}")
        sd_pass = sum(1 for r in sd_results if 'PASS' in r['verdict'])
        log(f"  대조군: {sd_pass}/{len(sd_results)} PASS")
    log()

    # 실험군
    log("━━ 실험군: 비자기쌍대 (패리티 깨짐 예측) ━━")
    if nsd_results:
        log(f"{'L-func':<15} {'t₀':>10} {'|Re(c₀)|/|Im|':>15} {'|Im(c₁)|/|Re|':>15} {'|Re(c₂)|/|Im|':>15} {'|Im(c₃)|/|Re|':>15} {'판정':>8}")
        log("-" * 100)
        for r in nsd_results:
            log(f"{r['lfunc']:<15} {r['t0']:>10.4f} {r['r0']:>15.2e} {r['r1']:>15.2e} {r['r2']:>15.2e} {r['r3']:>15.2e} {'BROKEN' if 'BROKEN' in r['verdict'] else 'INTACT':>8}")
        nsd_broken = sum(1 for r in nsd_results if 'BROKEN' in r['verdict'])
        log(f"  실험군: {nsd_broken}/{len(nsd_results)} BROKEN")
    log()

    # 대비 통계
    log("━━ 대비 통계 ━━")
    if sd_results and nsd_results:
        for label, idx in [("c₀ (Re/Im)", 'r0'), ("c₁ (Im/Re)", 'r1'), ("c₂ (Re/Im)", 'r2'), ("c₃ (Im/Re)", 'r3')]:
            sd_vals = [r[idx] for r in sd_results]
            nsd_vals = [r[idx] for r in nsd_results]
            log(f"  {label}:")
            log(f"    자기쌍대:   mean={np.mean(sd_vals):.2e}, max={np.max(sd_vals):.2e}")
            log(f"    비자기쌍대: mean={np.mean(nsd_vals):.2e}, max={np.max(nsd_vals):.2e}")
            if np.mean(sd_vals) > 1e-30:
                log(f"    비율 (nsd/sd): {np.mean(nsd_vals)/np.mean(sd_vals):.1f}×")
            log()

    # 최종 판정
    log("━━ 최종 판정 ━━")
    sd_all_pass = all('PASS' in r['verdict'] for r in sd_results) if sd_results else False
    nsd_all_broken = all('BROKEN' in r['verdict'] for r in nsd_results) if nsd_results else False

    if sd_all_pass and nsd_all_broken:
        log("★★★★ SHARPNESS 확립:")
        log("  - 자기쌍대: 패리티 성립 (FE+SR → cₙ=(-1)^(n+1)c̄ₙ)")
        log("  - 비자기쌍대: 패리티 깨짐 (SR 없음 → 제약 해소)")
        log("  → 정리가 SHARP: self-dual 조건은 필요충분")
    elif sd_all_pass and not nsd_all_broken:
        log("⚠️ 예상 밖: 비자기쌍대에서도 패리티 유지")
        log("  → FE만으로 충분? 이론적 재검토 필요")
    elif not sd_all_pass:
        log("❌ 대조군 실패: 자기쌍대에서 패리티 깨짐")
        log("  → 수치 정밀도 또는 코드 버그 점검 필요")
    else:
        log("⚠️ 혼합 결과: 개별 분석 필요")

else:
    log("결과 없음")

log()
log(f"총 소요: {time.time()-START:.1f}s")
log("=" * 72)
outf.close()
