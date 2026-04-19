#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #98b — B-05 σ-FAIL S(σ,d) — Part B 재개 (GL(5))
=============================================================================
상황: #98 Part A(GL(4)) 완료, Part B(GL(5)) σ=2.1(offset=-0.4) 1개만 완료 후 종료.
      Part B 나머지 6개 σ + Part C 통합 분석 실행.

Part A 확정 결과 (hardcoded):
  GL(4) sym³(11a1): S=45 for ALL σ → ratio=1.000 (완전 평탄)

Part B 부분 결과 (σ=2.1, offset=-0.4 → S=60 이미 확인):
  GL(5) sym⁴(11a1), center=2.5
  나머지 offsets: -0.2, -0.1, 0.0, +0.1, +0.2, +0.4
=============================================================================
"""

import sys, os, time
import numpy as np

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results",
    "sigma_fail_quantification_98b.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    lines.append(msg)
    print(msg, flush=True)

def flush_file():
    with open(OUTFILE, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")

log("=" * 72)
log("[Project RDL] 결과 #98b — B-05 σ-FAIL S(σ,d) Part B 재개 (GL(5))")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()
log("Part A GL(4) 결과 (기확보): S=45 for ALL σ → ratio=1.000")
log("Part B 재개: GL(5) sym⁴(11a1), σ=2.1(S=60) 기확인, 나머지 6개 σ 실행")
log()
flush_file()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PARI 초기화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)
gp("default(realprecision, 50)")
log("PARI 초기화: 2GB 메모리, realprecision=50")
log()
flush_file()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 유틸리티
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def count_sign_changes(vals):
    signs = np.sign(vals)
    signs = signs[signs != 0]
    return int(np.sum(np.diff(signs) != 0))


def sigma_sweep_pari_v2(L_name, Linit_name, center, sigma_offsets, t_min, t_max, t_step=0.1):
    t_vals = np.arange(t_min, t_max, t_step)
    n_pts = len(t_vals)
    log(f"\n  σ-스윕: {L_name}")
    log(f"  center={center}, t∈[{t_min:.1f}, {t_max:.1f}], {n_pts}점, Δt={t_step}")

    results = {}
    for offset in sigma_offsets:
        sigma = center + offset
        t0 = time.time()

        re_vals = []
        for t in t_vals:
            try:
                cmd = f"real(lfunlambda({Linit_name}, {sigma:.6f} + {t:.6f}*I))"
                val = float(gp(cmd))
                re_vals.append(val)
            except Exception as e:
                print(f"WARNING: σ={sigma:.1f}, t={t:.2f}: {e}", flush=True)
                re_vals.append(0.0)

        re_arr = np.array(re_vals)
        n_sc = count_sign_changes(re_arr)
        elapsed = time.time() - t0
        results[sigma] = n_sc
        log(f"    σ={sigma:.1f} (offset={offset:+.1f}): S={n_sc} 부호변환 ({elapsed:.1f}s)")
        flush_file()

    return results


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Part A 확정 결과 (GL(4) sym³(11a1))
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CENTER_4 = 2.0
OFFSETS = [-0.4, -0.2, -0.1, 0.0, +0.1, +0.2, +0.4]

# 기확보 결과
results_gl4 = {
    1.6: 45,  # offset=-0.4
    1.8: 45,  # offset=-0.2
    1.9: 45,  # offset=-0.1
    2.0: 45,  # center
    2.1: 45,  # offset=+0.1
    2.2: 45,  # offset=+0.2
    2.4: 45,  # offset=+0.4
}
s_center_4 = results_gl4[CENTER_4]

log("=" * 72)
log("[Part A] GL(4) sym³(11a1) — 기확보 결과 로드")
log("=" * 72)
log(f"  center={CENTER_4}, S(center)={s_center_4}")
log("  S(σ)/S(center) 프로파일:")
for offset in OFFSETS:
    sigma = CENTER_4 + offset
    s_val = results_gl4[sigma]
    ratio = s_val / s_center_4
    log(f"    σ={sigma:.1f} (offset={offset:+.1f}): S={s_val}, ratio={ratio:.3f}")
log()
flush_file()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Part B GL(5) sym⁴(11a1) — center=2.5, degree 5
# σ=2.1 (offset=-0.4) → S=60 이미 확인
# 나머지 6개 σ 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("=" * 72)
log("[Part B] GL(5) sym⁴(11a1) — σ-스윕 (재개)")
log("=" * 72)
log("  L-함수: lfunsympow(ellinit([0,-1,1,0,0]), 4)")
log("  degree=5, w=1, center=2.5")
log("  기확인: σ=2.1 (offset=-0.4) → S=60")
log()
flush_file()

t_start_b = time.time()

gp("E11 = ellinit([0,-1,1,0,0])")
gp("L5_11a1 = lfunsympow(E11, 4)")
log("  lfunsympow(E11, 4) 생성 완료")

gp("L5init = lfuninit(L5_11a1, [30, 0])")
log("  lfuninit 완료 (T=30)")

zeros_str = str(gp("lfunzeros(L5_11a1, 30)"))
zeros_5 = [float(x) for x in zeros_str.strip("[]").split(", ") if x.strip()]
log(f"  영점 {len(zeros_5)}개 발견 (t₁={zeros_5[0]:.4f})")
flush_file()

CENTER_5 = 2.5
T_MIN_5 = max(zeros_5[0] + 0.5, 2.0)
T_MAX_5 = min(zeros_5[0] + 50.0, 30.0)

# 기확인 σ=2.1 (offset=-0.4) 제외, 나머지 6개만 실행
OFFSETS_REMAINING = [-0.2, -0.1, 0.0, +0.1, +0.2, +0.4]
log(f"  남은 σ offsets: {OFFSETS_REMAINING}")
log(f"  t∈[{T_MIN_5:.1f}, {T_MAX_5:.1f}]")
log()
flush_file()

# 나머지 6개 σ 실행
results_gl5_new = sigma_sweep_pari_v2(
    "sym⁴(11a1)", "L5init", CENTER_5, OFFSETS_REMAINING,
    T_MIN_5, T_MAX_5, t_step=0.1
)

# 기확인 결과와 합산
results_gl5 = {CENTER_5 - 0.4: 60}  # σ=2.1 (offset=-0.4)
results_gl5.update(results_gl5_new)

log(f"\n  [Part B 완료] {time.time()-t_start_b:.0f}s")

s_center_5 = results_gl5.get(CENTER_5, 1)
log(f"\n  S(σ)/S(center) 프로파일 (center={CENTER_5}):")
for offset in OFFSETS:
    sigma = CENTER_5 + offset
    s_val = results_gl5.get(sigma, 0)
    ratio = s_val / max(s_center_5, 1)
    already = " ← 기확인" if offset == -0.4 else ""
    log(f"    σ={sigma:.1f} (offset={offset:+.1f}): S={s_val}, ratio={ratio:.3f}{already}")
flush_file()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Part C 5-degree 통합 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log()
log("=" * 72)
log("[Part C] 5-degree 통합 분석")
log("=" * 72)
log()

log("기존 데이터 요약:")
log()
log("  d=1 GL(1):")
log("    ζ(s):  S(σ)=3 at all σ → σ-무관 (FAIL)")
log("    χ₃:   S(0.5)=21, S(off)=8 → PASS (ratio=0.38)")
log("    χ₄:   S(0.5)=9, S(1.5)=10 → FAIL")
log("    χ₅:   S(σ)=10 at all σ → FAIL")
log("    χ₇:   S(σ)=12 at all σ → FAIL")
log("    → GL(1) 평균 패턴: 대부분 σ-무관, χ₃만 PASS")
log()
log("  d=2 GL(2):")
log("    11a1: S(σ)=17 at all σ → 완전 평탄 (FAIL)")
log("    37a1: S(σ)=23 at all σ → 완전 평탄 (FAIL)")
log("    Δ:    S(σ)∈{8,9} at σ → 거의 평탄 (FAIL)")
log()
log("  d=3 GL(3):")
gl3_data = {0.3: 31, 0.5: 27, 0.7: 31, 0.9: 33}
s_crit_3 = gl3_data[0.5]
log(f"    sym²(11a1): S(σ) = {gl3_data}")
log(f"    center=0.5, S(center)={s_crit_3}")
log(f"    패턴: 비단조, off-center에서 증가 (FAIL)")
log()
log("  d=4 GL(4):")
log(f"    sym³(11a1): S(σ) = {results_gl4}")
log(f"    center={CENTER_4}, 완전 평탄 (ratio=1.000 all)")
log()
log("  d=5 GL(5):")
log(f"    sym⁴(11a1): S(σ) = {results_gl5}")
log(f"    center={CENTER_5}")
log()
flush_file()

# σ-비대칭 지표 A(d) 계산
log("=" * 72)
log("[Part C-2] σ-비대칭 지표 A(d) = max(S_off)/S_crit - 1")
log("=" * 72)
log()

asymmetry = {}

a1_pass = 8/21 - 1   # χ₃: off_max/crit - 1 < 0 → PASS
a1_fail = 3/3 - 1    # ζ: 0.0 → FAIL
log(f"  d=1 GL(1) χ₃ (PASS): A = {a1_pass:.3f}")
log(f"  d=1 GL(1) ζ  (FAIL): A = {a1_fail:.3f}")
asymmetry[1] = a1_pass  # χ₃ 대표

a2 = 17/17 - 1
log(f"  d=2 GL(2) 11a1: A = {a2:.3f}")
asymmetry[2] = a2

a3 = max(gl3_data.values()) / gl3_data[0.5] - 1
log(f"  d=3 GL(3) sym²(11a1): A = {a3:.3f}")
asymmetry[3] = a3

# d=4
off_max_4 = max(v for k, v in results_gl4.items() if k != CENTER_4)
a4 = off_max_4 / s_center_4 - 1
log(f"  d=4 GL(4) sym³(11a1): A = {a4:.3f}")
asymmetry[4] = a4

# d=5
if s_center_5 > 0:
    off_max_5 = max(v for k, v in results_gl5.items() if k != CENTER_5)
    a5 = off_max_5 / s_center_5 - 1
    log(f"  d=5 GL(5) sym⁴(11a1): A = {a5:.3f}")
    asymmetry[5] = a5
else:
    log(f"  d=5 GL(5): S(center)=0, 비대칭 지표 계산 불가")
    a5 = float('nan')
    asymmetry[5] = float('nan')

log()

# 5-degree 통합 요약표
log("=" * 72)
log("[Part C-3] 5-degree 통합 요약표")
log("=" * 72)
log()
log(f"  {'d':>3} | {'L-함수':<20} | {'S(crit)':>8} | {'S(off) range':>14} | {'A(d)':>8} | {'판정':>6}")
log(f"  {'---':>3}-+-{'-'*20}-+-{'-'*8}-+-{'-'*14}-+-{'-'*8}-+-{'-'*6}")

log(f"  {'1':>3} | {'χ₃ (GL(1))':.<20} | {'21':>8} | {'8':>14} | {a1_pass:>8.3f} | {'PASS':>6}")
log(f"  {'1':>3} | {'ζ (GL(1))':.<20} | {'3':>8} | {'3':>14} | {a1_fail:>8.3f} | {'FAIL':>6}")
log(f"  {'2':>3} | {'11a1 (GL(2))':.<20} | {'17':>8} | {'17':>14} | {a2:>8.3f} | {'FAIL':>6}")
log(f"  {'3':>3} | {'sym²(11a1) (GL(3))':.<20} | {s_crit_3:>8} | {f'{min(gl3_data.values())}-{max(gl3_data.values())}':>14} | {a3:>8.3f} | {'FAIL':>6}")

off_vals_4 = [v for k, v in sorted(results_gl4.items()) if k != CENTER_4]
log(f"  {'4':>3} | {'sym³(11a1) (GL(4))':.<20} | {s_center_4:>8} | {f'{min(off_vals_4)}-{max(off_vals_4)}':>14} | {a4:>8.3f} | {'PASS' if a4 < 0 else 'FAIL':>6}")

if s_center_5 > 0:
    off_vals_5 = [v for k, v in sorted(results_gl5.items()) if k != CENTER_5]
    s5_min = min(off_vals_5) if off_vals_5 else s_center_5
    s5_max = max(off_vals_5) if off_vals_5 else s_center_5
    log(f"  {'5':>3} | {'sym⁴(11a1) (GL(5))':.<20} | {s_center_5:>8} | {f'{s5_min}-{s5_max}':>14} | {a5:>8.3f} | {'PASS' if a5 < 0 else 'FAIL':>6}")

log()
flush_file()

# σ-오프셋별 비율 비교표
log("=" * 72)
log("[Part C-4] σ-오프셋별 S(σ)/S(center) 비교")
log("=" * 72)
log()
log(f"  {'offset':>8} | {'d=1(ζ)':>8} | {'d=2(11a1)':>10} | {'d=3(sym²)':>10} | {'d=4(sym³)':>10} | {'d=5(sym⁴)':>10}")
log(f"  {'-'*8}-+-{'-'*8}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}")

gl3_norm = {}
for sig, val in gl3_data.items():
    gl3_norm[sig - 0.5] = val / max(s_crit_3, 1)

for offset in OFFSETS:
    d1_r = 1.0  # ζ
    d2_r = 1.0  # 11a1
    d3_r = gl3_norm.get(offset, float('nan'))
    sig4 = CENTER_4 + offset
    d4_r = results_gl4.get(sig4, 0) / max(s_center_4, 1)
    sig5 = CENTER_5 + offset
    d5_r = results_gl5.get(sig5, 0) / max(s_center_5, 1) if s_center_5 > 0 else float('nan')

    d3_str = f"{d3_r:.3f}" if not np.isnan(d3_r) else "N/A"
    d4_str = f"{d4_r:.3f}"
    d5_str = f"{d5_r:.3f}" if not np.isnan(d5_r) else "N/A"

    log(f"  {offset:>+8.1f} | {d1_r:>8.3f} | {d2_r:>10.3f} | {d3_str:>10} | {d4_str:>10} | {d5_str:>10}")

log()
flush_file()

# A(d) vs degree 피팅
log("=" * 72)
log("[Part C-5] A(d) vs degree 피팅")
log("=" * 72)
log()

degrees = []
A_vals = []
for d in [1, 2, 3, 4, 5]:
    if d in asymmetry and not np.isnan(asymmetry[d]):
        degrees.append(d)
        A_vals.append(asymmetry[d])

degrees = np.array(degrees, dtype=float)
A_vals = np.array(A_vals)

log(f"  피팅 데이터: {len(degrees)}점")
for d, a in zip(degrees, A_vals):
    log(f"    d={int(d)}: A={a:.4f}")
log()

if len(degrees) >= 3:
    coeffs_lin = np.polyfit(degrees, A_vals, 1)
    A_pred_lin = np.polyval(coeffs_lin, degrees)
    ss_res = np.sum((A_vals - A_pred_lin)**2)
    ss_tot = np.sum((A_vals - np.mean(A_vals))**2)
    R2_lin = 1 - ss_res / max(ss_tot, 1e-30)
    log(f"  선형 피팅: A(d) = {coeffs_lin[0]:.4f}·d + {coeffs_lin[1]:.4f}")
    log(f"    R² = {R2_lin:.4f}")
    log()

    if len(degrees) >= 4:
        coeffs_quad = np.polyfit(degrees, A_vals, 2)
        A_pred_quad = np.polyval(coeffs_quad, degrees)
        ss_res_q = np.sum((A_vals - A_pred_quad)**2)
        R2_quad = 1 - ss_res_q / max(ss_tot, 1e-30)
        log(f"  2차 피팅: A(d) = {coeffs_quad[0]:.4f}·d² + {coeffs_quad[1]:.4f}·d + {coeffs_quad[2]:.4f}")
        log(f"    R² = {R2_quad:.4f}")
        log()

# 최종 판정
log("=" * 72)
log("[결론] B-05 σ-FAIL S(σ,d) 정량화 (5-degree 통합)")
log("=" * 72)
log()
log("  1. σ-프로파일 패턴 (degree별):")
log("     d=1: 대부분 평탄 (ζ,χ₄,χ₅,χ₇), χ₃만 center 극대 (PASS)")
log("     d=2: 완전 평탄 — σ 무관")
log(f"     d=3: center에서 극소, off-center 증가 (FAIL, 비단조)")

# d=4 판정
if abs(a4) < 0.05:
    pattern_4 = "완전 평탄 (σ-포화 체제)"
elif a4 < -0.05:
    pattern_4 = "center 극대 (PASS 패턴)"
else:
    pattern_4 = f"비대칭 A={a4:.3f} (FAIL)"
log(f"     d=4: {pattern_4}")

# d=5 판정
if not np.isnan(a5):
    if abs(a5) < 0.05:
        pattern_5 = "완전 평탄 (σ-포화 체제 — d=4 동일)"
    elif a5 < -0.05:
        pattern_5 = "center 극대 (PASS 패턴)"
    else:
        pattern_5 = f"비대칭 A={a5:.3f}"
    log(f"     d=5: {pattern_5}")

log()
log("  2. σ-비대칭 지표 A(d) 요약:")
log("     A < 0 → center가 극대 (σ-유일성 PASS)")
log("     A ≈ 0 → 완전 평탄 (σ-유일성 퇴화)")
log("     A > 0 → center에서 극소 (σ-유일성 역전)")
for d in [1, 2, 3, 4, 5]:
    if d in asymmetry and not np.isnan(asymmetry[d]):
        a = asymmetry[d]
        status = "PASS" if a < -0.05 else ("평탄" if abs(a) < 0.05 else "FAIL")
        log(f"     d={d}: A = {a:.4f} ({status})")
log()
log("  3. 경계 B-05 갱신:")
if 5 in asymmetry and not np.isnan(asymmetry[5]):
    if abs(a5) < 0.05:
        log("     d=4,5 모두 S=const → ★★ d≥4 포화 경계 확립")
        log("     σ-유일성은 d=1~3 (불규칙 영점)에서만 판별력 보유")
        log("     d≥4 (정현파적 Hardy Z, 균일 영점)에서는 퇴화")
    else:
        log(f"     d=5 비단조 (A={a5:.3f}) → 포화 경계 재검토 필요")
log()
log("  4. 논문 반영 (#99 지시에 따라):")
log("     - σ-유일성 섹션: d=4,5 S(σ)=const 표 추가")
log("     - 'conditional validity (d=1~3) + saturation (d≥4)' 기술")
log("     - Summary Table #98 행: B-05 σ-FAIL 경계")

total_time = time.time() - t_start_b
log()
log(f"Part B+C 소요 시간: {total_time:.0f}s ({total_time/60:.1f}분)")
flush_file()

log()
log(f"결과 저장: {OUTFILE}")
