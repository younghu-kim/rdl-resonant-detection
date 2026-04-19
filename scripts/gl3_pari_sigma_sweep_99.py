#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #99 — GL(3) sym²(11a1) PARI σ-sweep 부호변환 재측정
=============================================================================
목적: #98에서 d=3만 비대칭(A=0.222)이었으나, 기존 #63 측정은
      AFE + Δt=0.8 (41점) → PARI + Δt=0.1 (280점)과 조건 불균일.
      동일 PARI lfunlambda + Δt=0.1 조건으로 재측정하여
      d=3 비대칭이 실제인지 AFE 아티팩트인지 판별.

기존 #63 결과:
  sym²(11a1), center=0.5 (analytic normalization)
  S(σ) = {0.3: 31, 0.5: 27, 0.7: 31, 0.9: 33}
  A(d=3) = 33/27 - 1 = 0.222

본 실험:
  PARI lfunsympow(E, 2) → lfuninit → lfunlambda
  center = k/2 = 1.5 (PARI motivic normalization, k=3)
  σ offsets: [-0.4, -0.2, -0.1, 0.0, +0.1, +0.2, +0.4]
  t-range: 동적 (영점 기반), Δt=0.1
  → #98과 동일 방법/해상도
=============================================================================
"""

import sys, os, time
import numpy as np

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results",
    "gl3_pari_sigma_sweep_99.txt"
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
log("[Project RDL] 결과 #99 — GL(3) sym²(11a1) PARI σ-sweep 재측정")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ─── PARI 초기화 ────────────────────────────────────────────────────────
import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)
gp("default(realprecision, 50)")
log("PARI 초기화: 2GB 메모리, realprecision=50")
log()
flush_file()

# ─── L-함수 생성 ────────────────────────────────────────────────────────
log("=" * 72)
log("[1] GL(3) sym²(11a1) L-함수 생성")
log("=" * 72)
t0 = time.time()

gp("E11 = ellinit([0,-1,1,0,0])")
gp("L3 = lfunsympow(E11, 2)")
log("  lfunsympow(E11, 2) 생성 완료")

# lfuninit — 충분히 넓은 범위
gp("L3init = lfuninit(L3, [50, 0])")
log("  lfuninit 완료 (T=50)")

# 영점 확인
zeros_str = str(gp("lfunzeros(L3, 50)"))
zeros = [float(x) for x in zeros_str.strip("[]").split(", ") if x.strip()]
log(f"  영점 {len(zeros)}개 발견")
if len(zeros) > 0:
    log(f"  t₁={zeros[0]:.4f}, t_last={zeros[-1]:.4f}")

# center 확인
try:
    k_val = str(gp("lfunparams(L3)[3]"))
    log(f"  PARI k={k_val}")
except:
    pass

CENTER = 1.5  # sym²(11a1): degree 3, w=2, k=3, center=k/2=1.5
log(f"  center = {CENTER} (motivic, PARI k=3)")
log(f"  생성 시간: {time.time()-t0:.1f}s")
log()
flush_file()

# ─── σ-sweep ────────────────────────────────────────────────────────────
log("=" * 72)
log("[2] σ-sweep 부호변환 카운트 (PARI lfunlambda, Δt=0.1)")
log("=" * 72)
log()

OFFSETS = [-0.4, -0.2, -0.1, 0.0, +0.1, +0.2, +0.4]
T_MIN = max(zeros[0] + 0.5, 2.0)
T_MAX = min(zeros[0] + 50.0, 50.0)  # lfuninit T=50까지

t_vals = np.arange(T_MIN, T_MAX, 0.1)
n_pts = len(t_vals)
log(f"  center={CENTER}, t∈[{T_MIN:.1f}, {T_MAX:.1f}], {n_pts}점, Δt=0.1")
log(f"  σ offsets: {OFFSETS}")
log()
flush_file()

results = {}
for offset in OFFSETS:
    sigma = CENTER + offset
    t_start = time.time()
    re_vals = []

    for t in t_vals:
        try:
            cmd = f"real(lfunlambda(L3init, {sigma:.6f} + {t:.6f}*I))"
            val = float(gp(cmd))
            re_vals.append(val)
        except Exception:
            re_vals.append(0.0)

    re_arr = np.array(re_vals)
    # 부호변환 카운트
    signs = np.sign(re_arr)
    signs = signs[signs != 0]
    n_sc = int(np.sum(np.diff(signs) != 0))

    elapsed = time.time() - t_start
    results[sigma] = n_sc
    log(f"  σ={sigma:.1f} (offset={offset:+.1f}): S={n_sc} 부호변환 ({elapsed:.1f}s)")
    flush_file()

# ─── 분석 ───────────────────────────────────────────────────────────────
log()
log("=" * 72)
log("[3] 분석 — #63 AFE 결과와 비교")
log("=" * 72)
log()

s_center = results.get(CENTER, 1)
log(f"  S(σ)/S(center) 프로파일 (center={CENTER}):")
for offset in OFFSETS:
    sigma = CENTER + offset
    s_val = results.get(sigma, 0)
    ratio = s_val / max(s_center, 1)
    log(f"    σ={sigma:.1f}: S={s_val}, ratio={ratio:.3f}")

log()

# 비대칭 지표
if s_center > 0:
    off_max = max(v for k, v in results.items() if k != CENTER)
    A_new = off_max / s_center - 1
    log(f"  비대칭 지표 A(d=3) [PARI] = {A_new:.4f}")
else:
    A_new = float('nan')
    log(f"  비대칭 지표 계산 불가 (S(center)=0)")

log()

# #63 AFE 결과 (analytic normalization, center=0.5)
# σ = {0.3: 31, 0.5: 27, 0.7: 31, 0.9: 33}
A_old = 33/27 - 1  # = 0.222
log(f"  비교:")
log(f"    #63 AFE (Δt=0.8, 41점, center=0.5): A = {A_old:.4f}")
log(f"    #99 PARI (Δt=0.1, {n_pts}점, center={CENTER}): A = {A_new:.4f}")
log()

if abs(A_new) < 0.05:
    verdict = "평탄 (d=3 비대칭은 AFE 아티팩트)"
    log(f"  ★ 판정: {verdict}")
    log(f"    → #63의 비대칭(A=0.222)은 AFE 방법 + 낮은 해상도의 산물")
    log(f"    → d=3도 d=2,4,5와 동일하게 평탄")
    log(f"    → 결론: ALL d≥2에서 S(σ)≈const")
elif A_new > 0.05:
    verdict = f"비대칭 재현 (A={A_new:.4f})"
    log(f"  ★ 판정: {verdict}")
    log(f"    → d=3 비대칭은 실제 현상 (AFE 아티팩트 아님)")
    log(f"    → d=3 고유 σ-의존성 확립")
    log(f"    → d=2,4,5 평탄 vs d=3 비대칭 → degree 3 특수성")
elif A_new < -0.05:
    verdict = f"역전 (center 극대, A={A_new:.4f})"
    log(f"  ★ 판정: {verdict}")
    log(f"    → d=3에서 center가 극대 (PASS 패턴)")
    log(f"    → #63과 완전 불일치 → 측정 방법 의존성 강함")

log()

# ─── #98과 통합 ─────────────────────────────────────────────────────────
log("=" * 72)
log("[4] 5-degree 통합 (PARI 통일 측정 기준)")
log("=" * 72)
log()

log(f"  {'d':>3} | {'L-함수':<20} | {'S(crit)':>8} | {'S range':>10} | {'A(d)':>8} | {'판정':>10}")
log(f"  {'---':>3}-+-{'-'*20}-+-{'-'*8}-+-{'-'*10}-+-{'-'*8}-+-{'-'*10}")
log(f"  {'1':>3} | {'χ₃ (GL(1))':.<20} | {'21':>8} | {'8-21':>10} | {'-0.619':>8} | {'PASS':>10}")
log(f"  {'2':>3} | {'11a1 (GL(2))':.<20} | {'17':>8} | {'17':>10} | {'0.000':>8} | {'평탄':>10}")
log(f"  {'3':>3} | {'sym²(11a1) (GL(3))':.<20} | {s_center:>8} | {f'{min(results.values())}-{max(results.values())}':>10} | {A_new:>8.4f} | {'평탄' if abs(A_new)<0.05 else '비대칭':>10}")
log(f"  {'4':>3} | {'sym³(11a1) (GL(4))':.<20} | {'45':>8} | {'45':>10} | {'0.000':>8} | {'평탄':>10}")
log(f"  {'5':>3} | {'sym⁴(11a1) (GL(5))':.<20} | {'59':>8} | {'59-60':>10} | {'0.017':>8} | {'평탄':>10}")

log()

total = time.time() - t0
log(f"총 소요 시간: {total:.0f}s ({total/60:.1f}분)")
flush_file()
log(f"결과 저장: {OUTFILE}")
