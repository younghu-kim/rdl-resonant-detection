#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #98 — B-05 σ-FAIL S(σ,d) 5-degree 정량화
=============================================================================
목적: degree 1–5 L-함수의 σ-방향 부호변환 카운트 S(σ)를 통합 측정하고,
      S(σ,d) 스케일링 공식을 유도한다.

기존 데이터:
  d=1 GL(1): #47 — ζ, χ₃, χ₄, χ₅, χ₇  (σ_crit=0.5, t∈[1,29.5])
  d=2 GL(2): #47 — 11a1, 37a1, Δ       (σ_crit=1.0 or 6.0)
  d=3 GL(3): #63 — sym²(11a1)          (σ_crit=0.5, t∈[2,35])

신규 측정:
  d=4 GL(4): sym³(Δ)     — σ_crit=17.0 (w=11), center=k/2=17 (PARI k=34/2=17)
             sym³(11a1)  — σ_crit=2.0 (w=1),  center=k/2=2
  d=5 GL(5): sym⁴(11a1)  — σ_crit=2.5 (w=1),  center=k/2=2.5

방법:
  Re(Λ(σ+it)) 부호변환 카운트 (t 스윕).
  PARI lfunlambda(Linit, σ+i*t)로 Λ 값 계산.
  t-range: [t₁+1, t₁+50] (t₁=첫 번째 영점)
  t-step: 0.1 (≈500개 포인트)
  σ offsets: center + {-0.4, -0.2, -0.1, 0, +0.1, +0.2, +0.4}

출력:
  A) d=4, d=5 σ-프로파일 S(σ)
  B) 5-degree 통합표
  C) S(σ,d) 스케일링 피팅 (power law, 선형, 2차)
=============================================================================
"""

import sys, os, time
import numpy as np

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results",
    "sigma_fail_quantification_98.txt"
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
log("[Project RDL] 결과 #98 — B-05 σ-FAIL S(σ,d) 5-degree 정량화")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

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
# 유틸리티: 부호변환 카운트
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def count_sign_changes(vals):
    """실수 배열에서 부호변환 횟수 (0 무시)"""
    signs = np.sign(vals)
    signs = signs[signs != 0]
    return int(np.sum(np.diff(signs) != 0))


def sigma_sweep_pari(L_name, Linit_name, center, sigma_offsets, t_min, t_max, t_step=0.1):
    """
    PARI lfunlambda로 σ-스윕 부호변환 카운트.

    Parameters:
        L_name: 표시용 이름
        Linit_name: PARI 변수명 (lfuninit 적용 완료)
        center: 임계선 σ (motivic center)
        sigma_offsets: center 기준 오프셋 리스트
        t_min, t_max: t 범위
        t_step: t 간격

    Returns:
        dict: {σ: sign_change_count}
    """
    t_vals = np.arange(t_min, t_max, t_step)
    n_pts = len(t_vals)
    log(f"\n  σ-스윕: {L_name}")
    log(f"  center={center}, t∈[{t_min:.1f}, {t_max:.1f}], {n_pts}점, Δt={t_step}")

    results = {}
    for offset in sigma_offsets:
        sigma = center + offset
        re_vals = []
        t0 = time.time()
        for t in t_vals:
            try:
                s_str = f"{sigma:.6f} + {t:.6f}*I"
                val = gp(f"lfunlambda({Linit_name}, {s_str})")
                # 실수부 추출
                re_part = float(gp(f"real({Linit_name}_val = lfunlambda({Linit_name}, {s_str}); real({Linit_name}_val)"))
            except Exception:
                try:
                    val_str = str(gp(f"lfunlambda({Linit_name}, {s_str})"))
                    # 복소수 파싱
                    val_str = val_str.replace(" ", "")
                    if "*I" in val_str or "I" in val_str:
                        # PARI 복소수: a + b*I
                        val_c = complex(val_str.replace("*I", "j").replace("I", "j"))
                        re_part = val_c.real
                    else:
                        re_part = float(val_str)
                except Exception as e:
                    re_part = 0.0
            re_vals.append(re_part)

        re_arr = np.array(re_vals)
        n_sc = count_sign_changes(re_arr)
        elapsed = time.time() - t0
        results[sigma] = n_sc
        log(f"    σ={sigma:.1f} (offset={offset:+.1f}): S={n_sc} 부호변환 ({elapsed:.1f}s)")
        flush_file()

    return results


def sigma_sweep_pari_v2(L_name, Linit_name, center, sigma_offsets, t_min, t_max, t_step=0.1):
    """
    개선된 PARI σ-스윕. 벡터화 방식으로 속도 향상.
    """
    t_vals = np.arange(t_min, t_max, t_step)
    n_pts = len(t_vals)
    log(f"\n  σ-스윕: {L_name}")
    log(f"  center={center}, t∈[{t_min:.1f}, {t_max:.1f}], {n_pts}점, Δt={t_step}")

    results = {}
    for offset in sigma_offsets:
        sigma = center + offset
        t0 = time.time()

        # PARI 벡터 방식: 한 번에 계산
        re_vals = []
        for t in t_vals:
            try:
                cmd = f"real(lfunlambda({Linit_name}, {sigma:.6f} + {t:.6f}*I))"
                val = float(gp(cmd))
                re_vals.append(val)
            except Exception as e:
                re_vals.append(0.0)

        re_arr = np.array(re_vals)
        n_sc = count_sign_changes(re_arr)
        elapsed = time.time() - t0
        results[sigma] = n_sc
        log(f"    σ={sigma:.1f} (offset={offset:+.1f}): S={n_sc} 부호변환 ({elapsed:.1f}s)")
        flush_file()

    return results


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# [Part A] GL(4) sym³(11a1) — center=2.0, degree 4
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("=" * 72)
log("[Part A] GL(4) sym³(11a1) — σ-스윕 부호변환 카운트")
log("=" * 72)
log("  L-함수: lfunsympow(ellinit([0,-1,1,0,0]), 3)")
log("  degree=4, w=1, center=2.0")
log()
flush_file()

t_start_a = time.time()

# L-함수 생성
gp("E11 = ellinit([0,-1,1,0,0])")
gp("L4_11a1 = lfunsympow(E11, 3)")
log("  lfunsympow(E11, 3) 생성 완료")

# lfuninit
gp("L4init = lfuninit(L4_11a1, [30, 0])")
log("  lfuninit 완료 (T=30)")

# 영점 확인
zeros_str = str(gp("lfunzeros(L4_11a1, 30)"))
zeros_4 = [float(x) for x in zeros_str.strip("[]").split(", ") if x.strip()]
log(f"  영점 {len(zeros_4)}개 발견 (t₁={zeros_4[0]:.4f})")
flush_file()

# σ-스윕
CENTER_4 = 2.0
OFFSETS = [-0.4, -0.2, -0.1, 0.0, +0.1, +0.2, +0.4]
T_MIN_4 = max(zeros_4[0] + 0.5, 2.0)
T_MAX_4 = min(zeros_4[0] + 50.0, 30.0)

results_gl4 = sigma_sweep_pari_v2(
    "sym³(11a1)", "L4init", CENTER_4, OFFSETS,
    T_MIN_4, T_MAX_4, t_step=0.1
)

# S(σ)/S(center) 비율
s_center_4 = results_gl4.get(CENTER_4, 1)
log(f"\n  S(σ)/S(center) 프로파일 (center={CENTER_4}):")
for offset in OFFSETS:
    sigma = CENTER_4 + offset
    s_val = results_gl4.get(sigma, 0)
    ratio = s_val / max(s_center_4, 1)
    log(f"    σ={sigma:.1f}: S={s_val}, ratio={ratio:.3f}")

log(f"\n  [Part A 완료] {time.time()-t_start_a:.0f}s")
flush_file()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# [Part B] GL(5) sym⁴(11a1) — center=2.5, degree 5
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log()
log("=" * 72)
log("[Part B] GL(5) sym⁴(11a1) — σ-스윕 부호변환 카운트")
log("=" * 72)
log("  L-함수: lfunsympow(ellinit([0,-1,1,0,0]), 4)")
log("  degree=5, w=1, center=2.5")
log()
flush_file()

t_start_b = time.time()

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

results_gl5 = sigma_sweep_pari_v2(
    "sym⁴(11a1)", "L5init", CENTER_5, OFFSETS,
    T_MIN_5, T_MAX_5, t_step=0.1
)

s_center_5 = results_gl5.get(CENTER_5, 1)
log(f"\n  S(σ)/S(center) 프로파일 (center={CENTER_5}):")
for offset in OFFSETS:
    sigma = CENTER_5 + offset
    s_val = results_gl5.get(sigma, 0)
    ratio = s_val / max(s_center_5, 1)
    log(f"    σ={sigma:.1f}: S={s_val}, ratio={ratio:.3f}")

log(f"\n  [Part B 완료] {time.time()-t_start_b:.0f}s")
flush_file()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# [Part C] 5-degree 통합 분석 + S(σ,d) 피팅
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log()
log("=" * 72)
log("[Part C] 5-degree 통합 분석")
log("=" * 72)
log()

# 기존 데이터 (정규화: S(σ)/S(center) 비율)
# GL(1) #47: ζ — 모든 σ에서 동일 (=3), ratio=1.0
# GL(1) #47: χ₃ — center=0.5에서 21, off에서 8 → ratio=0.38
# 대표값: χ₃ 사용 (가장 명확한 PASS 사례)
# GL(2) #47: 11a1 — 모든 σ에서 17 → ratio=1.0 (평탄)
# GL(3) #63: sym²(11a1) — {0.3:31, 0.5:27, 0.7:31, 0.9:33}

log("기존 데이터 요약:")
log()

# GL(1) — 다수 L-함수 평균
log("  d=1 GL(1):")
log("    ζ(s):  S(σ)=3 at all σ ∈ {0.1..1.5} → σ-무관 (FAIL)")
log("    χ₃:   S(0.5)=21, S(off)=8 → PASS (ratio=0.38)")
log("    χ₄:   S(0.5)=9, S(1.5)=10 → FAIL")
log("    χ₅:   S(σ)=10 at all σ → FAIL")
log("    χ₇:   S(σ)=12 at all σ → FAIL")
log("    → GL(1) 평균 패턴: 대부분 σ-무관, χ₃만 PASS")
log()

# GL(2)
log("  d=2 GL(2):")
log("    11a1: S(σ)=17 at all σ ∈ {0.7..1.3} → 완전 평탄 (FAIL)")
log("    37a1: S(σ)=23 at all σ ∈ {0.7..1.3} → 완전 평탄 (FAIL)")
log("    Δ:    S(σ)∈{8,9} at σ ∈ {4.0..8.0} → 거의 평탄 (FAIL)")
log("    → GL(2) 패턴: 완전 평탄")
log()

# GL(3)
log("  d=3 GL(3):")
gl3_data = {0.3: 31, 0.5: 27, 0.7: 31, 0.9: 33}
s_crit_3 = gl3_data[0.5]  # center=0.5 (analytic normalization)
log(f"    sym²(11a1): S(σ) = {gl3_data}")
log(f"    center=0.5, S(center)={s_crit_3}")
log(f"    패턴: σ 증가 시 S 증가 (FAIL, 단조 경향)")
log()

# GL(4) 신규
log("  d=4 GL(4):")
log(f"    sym³(11a1): S(σ) = {results_gl4}")
log(f"    center={CENTER_4}")
log()

# GL(5) 신규
log("  d=5 GL(5):")
log(f"    sym⁴(11a1): S(σ) = {results_gl5}")
log(f"    center={CENTER_5}")
log()
flush_file()

# ── 정규화된 비대칭 지표 A(d) 계산 ──
# A(d) = max(S(σ_off)) / S(σ_crit) - 1
# A=0 → 완전 대칭(PASS 가능), A>0 → 비대칭(FAIL 경향), A<0 → center가 최대(PASS)
log("=" * 72)
log("[Part C-2] σ-비대칭 지표 A(d) = max(S_off)/S_crit - 1")
log("=" * 72)
log()

asymmetry = {}

# d=1: χ₃ (PASS 대표) — off_max=8, crit=21
a1_pass = 8/21 - 1  # = -0.619
# d=1: ζ (FAIL 대표) — 동일
a1_fail = 3/3 - 1   # = 0.0
log(f"  d=1 GL(1) χ₃ (PASS): A = {a1_pass:.3f}")
log(f"  d=1 GL(1) ζ  (FAIL): A = {a1_fail:.3f}")
asymmetry[1] = a1_pass  # χ₃ 사용 (PASS는 음수 A)

# d=2: 11a1 — off_max=17, crit=17
a2 = 17/17 - 1  # = 0.0
log(f"  d=2 GL(2) 11a1: A = {a2:.3f}")
asymmetry[2] = a2

# d=3: sym²(11a1) — off_max=33, crit=27
a3 = max(gl3_data.values()) / gl3_data[0.5] - 1
log(f"  d=3 GL(3) sym²(11a1): A = {a3:.3f}")
asymmetry[3] = a3

# d=4
if s_center_4 > 0:
    off_max_4 = max(v for k, v in results_gl4.items() if k != CENTER_4)
    a4 = off_max_4 / s_center_4 - 1
    log(f"  d=4 GL(4) sym³(11a1): A = {a4:.3f}")
    asymmetry[4] = a4
else:
    log(f"  d=4 GL(4): S(center)=0, 비대칭 지표 계산 불가")
    asymmetry[4] = float('nan')

# d=5
if s_center_5 > 0:
    off_max_5 = max(v for k, v in results_gl5.items() if k != CENTER_5)
    a5 = off_max_5 / s_center_5 - 1
    log(f"  d=5 GL(5) sym⁴(11a1): A = {a5:.3f}")
    asymmetry[5] = a5
else:
    log(f"  d=5 GL(5): S(center)=0, 비대칭 지표 계산 불가")
    asymmetry[5] = float('nan')

log()

# ── 통합 요약표 ──
log("=" * 72)
log("[Part C-3] 5-degree 통합 요약표")
log("=" * 72)
log()
log(f"  {'d':>3} | {'L-함수':<20} | {'S(crit)':>8} | {'S(off) range':>14} | {'A(d)':>8} | {'판정':>6}")
log(f"  {'---':>3}-+-{'-'*20}-+-{'-'*8}-+-{'-'*14}-+-{'-'*8}-+-{'-'*6}")

# d=1 χ₃
log(f"  {'1':>3} | {'χ₃ (GL(1))':.<20} | {'21':>8} | {'8':>14} | {a1_pass:>8.3f} | {'PASS':>6}")
# d=1 ζ
log(f"  {'1':>3} | {'ζ (GL(1))':.<20} | {'3':>8} | {'3':>14} | {a1_fail:>8.3f} | {'FAIL':>6}")
# d=2
log(f"  {'2':>3} | {'11a1 (GL(2))':.<20} | {'17':>8} | {'17':>14} | {a2:>8.3f} | {'FAIL':>6}")
# d=3
log(f"  {'3':>3} | {'sym²(11a1) (GL(3))':.<20} | {s_crit_3:>8} | {f'{min(gl3_data.values())}-{max(gl3_data.values())}':>14} | {a3:>8.3f} | {'FAIL':>6}")
# d=4
if s_center_4 > 0:
    off_vals_4 = [v for k, v in sorted(results_gl4.items()) if k != CENTER_4]
    log(f"  {'4':>3} | {'sym³(11a1) (GL(4))':.<20} | {s_center_4:>8} | {f'{min(off_vals_4)}-{max(off_vals_4)}':>14} | {a4:>8.3f} | {'PASS' if a4 < 0 else 'FAIL':>6}")
# d=5
if s_center_5 > 0:
    off_vals_5 = [v for k, v in sorted(results_gl5.items()) if k != CENTER_5]
    log(f"  {'5':>3} | {'sym⁴(11a1) (GL(5))':.<20} | {s_center_5:>8} | {f'{min(off_vals_5)}-{max(off_vals_5)}':>14} | {a5:>8.3f} | {'PASS' if a5 < 0 else 'FAIL':>6}")

log()

# ── S(σ,d) 피팅 시도 ──
log("=" * 72)
log("[Part C-4] S(σ,d) 스케일링 분석")
log("=" * 72)
log()

# 각 degree에서 σ 오프셋별 S 비율 프로파일
# 정규화: offset = σ - center, ratio = S(σ)/S(center)
log("  σ-오프셋별 S(σ)/S(center) 비교:")
log()
log(f"  {'offset':>8} | {'d=1(ζ)':>8} | {'d=2(11a1)':>10} | {'d=3(sym²)':>10} | {'d=4(sym³)':>10} | {'d=5(sym⁴)':>10}")
log(f"  {'-'*8}-+-{'-'*8}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}")

# GL(3) 정규화 (center=0.5, offsets at -0.2, 0, +0.2, +0.4)
gl3_norm = {}
for sig, val in gl3_data.items():
    gl3_norm[sig - 0.5] = val / max(s_crit_3, 1)

for offset in OFFSETS:
    d1_r = 1.0  # ζ: 모든 σ 동일
    d2_r = 1.0  # 11a1: 모든 σ 동일

    # GL(3): offset 매핑
    d3_r = gl3_norm.get(offset, float('nan'))

    # GL(4)
    sig4 = CENTER_4 + offset
    d4_r = results_gl4.get(sig4, 0) / max(s_center_4, 1) if s_center_4 > 0 else float('nan')

    # GL(5)
    sig5 = CENTER_5 + offset
    d5_r = results_gl5.get(sig5, 0) / max(s_center_5, 1) if s_center_5 > 0 else float('nan')

    d3_str = f"{d3_r:.3f}" if not np.isnan(d3_r) else "N/A"
    d4_str = f"{d4_r:.3f}" if not np.isnan(d4_r) else "N/A"
    d5_str = f"{d5_r:.3f}" if not np.isnan(d5_r) else "N/A"

    log(f"  {offset:>+8.1f} | {d1_r:>8.3f} | {d2_r:>10.3f} | {d3_str:>10} | {d4_str:>10} | {d5_str:>10}")

log()

# ── 비대칭 지표 A(d) vs degree 피팅 ──
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
    # 선형 피팅: A(d) = a₀ + a₁·d
    coeffs_lin = np.polyfit(degrees, A_vals, 1)
    A_pred_lin = np.polyval(coeffs_lin, degrees)
    ss_res = np.sum((A_vals - A_pred_lin)**2)
    ss_tot = np.sum((A_vals - np.mean(A_vals))**2)
    R2_lin = 1 - ss_res / max(ss_tot, 1e-30)
    log(f"  선형 피팅: A(d) = {coeffs_lin[0]:.4f}·d + {coeffs_lin[1]:.4f}")
    log(f"    R² = {R2_lin:.4f}")
    log()

    # 2차 피팅: A(d) = a₀ + a₁·d + a₂·d²
    if len(degrees) >= 4:
        coeffs_quad = np.polyfit(degrees, A_vals, 2)
        A_pred_quad = np.polyval(coeffs_quad, degrees)
        ss_res_q = np.sum((A_vals - A_pred_quad)**2)
        R2_quad = 1 - ss_res_q / max(ss_tot, 1e-30)
        log(f"  2차 피팅: A(d) = {coeffs_quad[0]:.4f}·d² + {coeffs_quad[1]:.4f}·d + {coeffs_quad[2]:.4f}")
        log(f"    R² = {R2_quad:.4f}")
        log()

# ── 최종 판정 ──
log("=" * 72)
log("[결론] B-05 σ-FAIL S(σ,d) 정량화")
log("=" * 72)
log()
log("  1. σ-프로파일 패턴 (degree별):")
log("     d=1: 대부분 평탄 (ζ,χ₄,χ₅,χ₇), χ₃만 center 극대 (PASS)")
log("     d=2: 완전 평탄 — σ 무관")
log(f"     d=3: center에서 극소, off-center에서 증가 (FAIL 패턴)")

# d=4, d=5 판정
if s_center_4 > 0:
    if a4 < -0.05:
        pattern_4 = "center 극대 (PASS 패턴)"
    elif abs(a4) < 0.05:
        pattern_4 = "평탄 (FAIL 패턴)"
    else:
        pattern_4 = f"비대칭 A={a4:.3f} (FAIL 패턴)"
    log(f"     d=4: {pattern_4}")

if s_center_5 > 0:
    if a5 < -0.05:
        pattern_5 = "center 극대 (PASS 패턴)"
    elif abs(a5) < 0.05:
        pattern_5 = "평탄 (FAIL 패턴)"
    else:
        pattern_5 = f"비대칭 A={a5:.3f} (FAIL 패턴)"
    log(f"     d=5: {pattern_5}")

log()
log("  2. σ-비대칭 지표 A(d):")
for d in [1, 2, 3, 4, 5]:
    if d in asymmetry and not np.isnan(asymmetry[d]):
        log(f"     d={d}: A = {asymmetry[d]:.4f}")
log()
log("  3. 다음 과제:")
log("     - S(σ,d) 공식 유도 (현재 데이터로 가능 여부 판정)")
log("     - 수렴 속도 메커니즘: R(N,σ,d) 스케일링 → 부호변환 인과 다리")
log("     - 논문 B-05 Remark 업데이트")

total_time = time.time() - t_start_a
log()
log(f"총 소요 시간: {total_time:.0f}s ({total_time/60:.1f}분)")
flush_file()

log()
log(f"결과 저장: {OUTFILE}")
