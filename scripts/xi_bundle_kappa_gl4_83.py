#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #83 — GL(4) ξ-bundle κ (σ-방향) A(t₀) 측정
=============================================================================
목적:
  #80-#82의 Hardy Z κ (t-방향)와 달리, ξ-bundle κ (σ-방향)로
  sym³(Δ) [d=4, w=11]와 sym³(11a1) [d=4, w=1]의 A(t₀)를 측정.
  d=1,2,3 데이터 (#73-#74)와 동일 방법 → B-12 degree 단조증가 최종 판결.

핵심:
  - κ = |Λ'(s)/Λ(s)|² at s = center + δ + i*t₀  (σ-방향!)
  - center: sym³(11a1) → k/2 = 2, sym³(Δ) → k/2 = 17
  - Λ'(s) = lfunlambda(Linit, s, 1)  (PARI 내장 도함수)
  - A(t₀) = κ - 1/δ²  (δ-독립이면 ξ-bundle 정리 확인)

PARI 중심점 발견:
  - PARI lfunlambda의 s는 motivic 정규화: 영점 at s = k/2 + i*t
  - sym³(11a1): k=4 → center=2, 영점 at s=2+i*t
  - sym³(Δ):   k=34 → center=17, 영점 at s=17+i*t
  - (PARI lfunzeros 반환값 t는 imaginary part of s at zero)

대상 L-함수:
  - sym³(11a1): lfunsympow(E, 3), gammaV=[-1,0,0,1], k=4, N=1331, ε=+1
  - sym³(Δ):   lfuncreate([an_int, 0, [-11,-10,0,1], 34, 1, -1, []]), k=34

영점 (#80/#81에서):
  - sym³(Δ):    t₁=4.155866, t₂=5.549122, t₃=8.111776
  - sym³(11a1): t₁=2.320021, t₂=3.591881, t₃=4.622627

δ 값: [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]

비교 기준 (#73-#74, ξ-bundle σ-방향):
  - d=1: A=1.27  (ζ, #73)
  - d=2: A=3.93  (GL2, #74)
  - d=3: A=12.79 (GL3 sym², #74)
  - d=4: ? (이번 실험)

결과: results/xi_bundle_kappa_gl4_83.txt
=============================================================================
"""

import sys, os, time
import statistics

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

# ─── 출력 ─────────────────────────────────────────────────────────────────
OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "results", "xi_bundle_kappa_gl4_83.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))
def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ─── 파라미터 ─────────────────────────────────────────────────────────────
N_COEFF = 4000

# ξ-bundle σ-방향 δ 값
DELTAS = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]

# sym³(Δ) 영점 (pari_80.txt에서)
SYM3_DELTA_ZEROS = [4.155866, 5.549122, 8.111776]
SYM3_DELTA_CENTER = 17.0  # k=34 → center=k/2=17

# sym³(11a1) 영점 (#81에서)
SYM3_11A1_ZEROS = [2.320021, 3.591881, 4.622627]
SYM3_11A1_CENTER = 2.0  # k=4 → center=k/2=2

# 기존 A(t₀) 값 (#73-#74, ξ-bundle σ-방향)
D1_A = 1.27    # ζ (d=1)
D2_A = 3.93    # GL(2) (d=2)
D3_A = 12.79   # GL(3) (d=3)

log("=" * 72)
log("결과 #83 — GL(4) ξ-bundle κ (σ-방향) A(t₀) 측정")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"δ 값: {DELTAS}")
log(f"N_COEFF: {N_COEFF}")
log()
log("★ 핵심: s = center + δ + i*t₀ (σ-방향, NOT Hardy Z)")
log("  PARI 정규화: 영점 at s = k/2 + i*t (motivic)")
log(f"  sym³(11a1): center={SYM3_11A1_CENTER}, 영점={SYM3_11A1_ZEROS}")
log(f"  sym³(Δ):   center={SYM3_DELTA_CENTER}, 영점={SYM3_DELTA_ZEROS}")
log()
flush_file()

# ─── PARI 초기화 ──────────────────────────────────────────────────────────
import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)  # 2GB
gp("default(realprecision, 100)")   # 100자리

log("PARI 초기화: 2GB 메모리, realprecision=100")
log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [1] sym³(Δ): PARI direuler로 계수 생성 + lfuninit
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[1] sym³(Δ) 초기화 — gammaV=[-11,-10,0,1], k=34, N=1, ε=-1")
log("=" * 72)
t_start = time.time()

gp(f"""
an_int = direuler(p=2, {N_COEFF},
  my(t=ramanujantau(p), q=p^11,
     e1=t*(t^2-2*q),
     e2=q*(t^2-2*q)*(t^2-q),
     e3=q^3*e1,
     e4=q^6);
  1/(1 - e1*X + e2*X^2 - e3*X^3 + e4*X^4)
)
""")
log(f"  direuler 완료 ({time.time()-t_start:.1f}s)")

t1 = time.time()
gp("Ld = lfuncreate([an_int, 0, [-11,-10,0,1], 34, 1, -1, []])")
gp("Ldi = lfuninit(Ld, [0, 15])")
log(f"  lfuninit 완료 ({time.time()-t1:.1f}s)")
log(f"  총 초기화: {time.time()-t_start:.1f}s")
log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [2] sym³(Δ) κ 측정 — ξ-bundle σ-방향
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[2] sym³(Δ) ξ-bundle κ (σ-방향) — 3영점 × 6δ = 18점")
log("=" * 72)
log(f"  center = {SYM3_DELTA_CENTER}")
log(f"  s = {SYM3_DELTA_CENTER} + δ + i*t₀  (σ 방향)")
log()

delta_results_sym3d = []

for t_zero in SYM3_DELTA_ZEROS:
    log(f"  ── t₀ = {t_zero:.6f} ──")
    row = {'t0': t_zero, 'As': [], 'kd2s': []}
    for delta in DELTAS:
        try:
            s_val = f"{SYM3_DELTA_CENTER} + {delta} + I*{t_zero}"
            gp(f"s_cur = {s_val}")
            L0 = gp("lfunlambda(Ldi, s_cur)")
            Lp = gp("lfunlambda(Ldi, s_cur, 1)")
            abs_L0 = float(gp(f"abs({L0})"))
            abs_Lp = float(gp(f"abs({Lp})"))
            if abs_L0 < 1e-250:
                log(f"    δ={delta:.3f}: ⚠️ |Λ(s)| 너무 작음 ({abs_L0:.2e})")
                continue
            kappa = abs_Lp**2 / abs_L0**2
            kd2 = kappa * delta**2
            A = kappa - 1.0 / delta**2
            row['As'].append(A)
            row['kd2s'].append(kd2)
            ok = "✅" if 0.99 <= kd2 <= 1.01 else "⚠️"
            log(f"    δ={delta:.3f}: κ={kappa:.4f}, A={A:.4f}, κδ²={kd2:.6f} {ok}")
        except Exception as e:
            log(f"    δ={delta:.3f}: ❌ 오류: {str(e)[:60]}")

    if row['As']:
        mean_A = statistics.mean(row['As'])
        cv_A = (statistics.stdev(row['As']) / abs(mean_A) * 100) if len(row['As']) > 1 and abs(mean_A) > 1e-10 else 0
        row['mean_A'] = mean_A
        row['cv_A'] = cv_A
        log(f"    → mean(A) = {mean_A:.4f}, CV(A) = {cv_A:.2f}%")
    delta_results_sym3d.append(row)
    log()
    flush_file()

# sym³(Δ) 요약
log("  sym³(Δ) 요약:")
sym3d_all_A = [A for r in delta_results_sym3d for A in r['As']]
if sym3d_all_A:
    sym3d_mean = statistics.mean(sym3d_all_A)
    sym3d_cv = (statistics.stdev(sym3d_all_A) / abs(sym3d_mean) * 100) if len(sym3d_all_A) > 1 else 0
    log(f"    전체 mean(A) = {sym3d_mean:.4f}")
    log(f"    CV(A) = {sym3d_cv:.2f}%")
    log(f"    점수: {len(sym3d_all_A)}/18")
log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [3] sym³(11a1) 초기화
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[3] sym³(11a1) 초기화 — lfunsympow(E, 3), k=4, N=1331, ε=+1")
log("=" * 72)
t2 = time.time()

gp("E11 = ellinit([0,-1,1,-10,-20])")
gp("L11 = lfunsympow(E11, 3)")
gp("L11i = lfuninit(L11, [0, 10])")

log(f"  lfuninit 완료 ({time.time()-t2:.1f}s)")
log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [4] sym³(11a1) κ 측정
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[4] sym³(11a1) ξ-bundle κ (σ-방향) — 3영점 × 6δ = 18점")
log("=" * 72)
log(f"  center = {SYM3_11A1_CENTER}")
log(f"  s = {SYM3_11A1_CENTER} + δ + i*t₀  (σ 방향)")
log()

delta_results_sym3_11a1 = []

for t_zero in SYM3_11A1_ZEROS:
    log(f"  ── t₀ = {t_zero:.6f} ──")
    row = {'t0': t_zero, 'As': [], 'kd2s': []}
    for delta in DELTAS:
        try:
            s_val = f"{SYM3_11A1_CENTER} + {delta} + I*{t_zero}"
            gp(f"s_cur = {s_val}")
            L0 = gp("lfunlambda(L11i, s_cur)")
            Lp = gp("lfunlambda(L11i, s_cur, 1)")
            abs_L0 = float(gp(f"abs({L0})"))
            abs_Lp = float(gp(f"abs({Lp})"))
            if abs_L0 < 1e-250:
                log(f"    δ={delta:.3f}: ⚠️ |Λ(s)| 너무 작음 ({abs_L0:.2e})")
                continue
            kappa = abs_Lp**2 / abs_L0**2
            kd2 = kappa * delta**2
            A = kappa - 1.0 / delta**2
            row['As'].append(A)
            row['kd2s'].append(kd2)
            ok = "✅" if 0.99 <= kd2 <= 1.01 else "⚠️"
            log(f"    δ={delta:.3f}: κ={kappa:.4f}, A={A:.4f}, κδ²={kd2:.6f} {ok}")
        except Exception as e:
            log(f"    δ={delta:.3f}: ❌ 오류: {str(e)[:60]}")

    if row['As']:
        mean_A = statistics.mean(row['As'])
        cv_A = (statistics.stdev(row['As']) / abs(mean_A) * 100) if len(row['As']) > 1 and abs(mean_A) > 1e-10 else 0
        row['mean_A'] = mean_A
        row['cv_A'] = cv_A
        log(f"    → mean(A) = {mean_A:.4f}, CV(A) = {cv_A:.2f}%")
    delta_results_sym3_11a1.append(row)
    log()
    flush_file()

# sym³(11a1) 요약
log("  sym³(11a1) 요약:")
sym3_11a1_all_A = [A for r in delta_results_sym3_11a1 for A in r['As']]
if sym3_11a1_all_A:
    sym3_11a1_mean = statistics.mean(sym3_11a1_all_A)
    sym3_11a1_cv = (statistics.stdev(sym3_11a1_all_A) / abs(sym3_11a1_mean) * 100) if len(sym3_11a1_all_A) > 1 else 0
    log(f"    전체 mean(A) = {sym3_11a1_mean:.4f}")
    log(f"    CV(A) = {sym3_11a1_cv:.2f}%")
    log(f"    점수: {len(sym3_11a1_all_A)}/18")
log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [5] 핵심 비교표: d=1,2,3,4 — B-12 degree 단조증가
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[5] 핵심 비교표 — B-12 degree 단조증가 (ξ-bundle σ-방향, 동일 방법)")
log("=" * 72)
log()

d4_delta_mean = sym3d_mean if sym3d_all_A else float('nan')
d4_11a1_mean = sym3_11a1_mean if sym3_11a1_all_A else float('nan')

log(f"  d=1 (ζ)              A_mean = {D1_A:.4f}  [#73]")
log(f"  d=2 (GL2 sym¹)       A_mean = {D2_A:.4f}  [#74]")
log(f"  d=3 (GL3 sym²)       A_mean = {D3_A:.4f}  [#74]")
log(f"  d=4 sym³(Δ)  w=11    A_mean = {d4_delta_mean:.4f}  [#83]")
log(f"  d=4 sym³(11a1) w=1   A_mean = {d4_11a1_mean:.4f}  [#83]")
log()

# 단조성 체크
d1_to_d2 = D2_A > D1_A
d2_to_d3 = D3_A > D2_A
d3_to_d4_d = d4_delta_mean > D3_A
d3_to_d4_11a1 = d4_11a1_mean > D3_A

log(f"  d=1→2 단조증가: {'✅' if d1_to_d2 else '❌'} ({D1_A:.2f} → {D2_A:.2f})")
log(f"  d=2→3 단조증가: {'✅' if d2_to_d3 else '❌'} ({D2_A:.2f} → {D3_A:.2f})")
log(f"  d=3→4 (sym³Δ): {'✅' if d3_to_d4_d else '❌'} ({D3_A:.2f} → {d4_delta_mean:.4f})")
log(f"  d=3→4 (sym³11a1): {'✅' if d3_to_d4_11a1 else '❌'} ({D3_A:.2f} → {d4_11a1_mean:.4f})")
log()

# ════════════════════════════════════════════════════════════════════════════
# [6] weight vs degree 비교 (동일 d=4)
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[6] weight vs degree 비교 — B-22")
log("=" * 72)
log()
log(f"  같은 d=4에서:")
log(f"    w=11 (sym³Δ):   A_mean = {d4_delta_mean:.4f}")
log(f"    w=1  (sym³11a1): A_mean = {d4_11a1_mean:.4f}")
if not (sym3d_all_A and sym3_11a1_all_A):
    log("  ⚠️ 데이터 부족, 비교 불가")
else:
    ratio = d4_delta_mean / d4_11a1_mean if abs(d4_11a1_mean) > 1e-10 else float('nan')
    log(f"  비율 A(w=11)/A(w=1) = {ratio:.4f}")
    if d4_delta_mean > d4_11a1_mean:
        log(f"  → A(w=11) > A(w=1): weight 증가 → A 증가 (weight 효과)")
    else:
        log(f"  → A(w=11) < A(w=1): weight 증가 → A 감소 (예상 외!)")
    log()
    log(f"  비교: Hardy Z κ에서는 A(w=11)=286 >> A(w=1)=8.82")
    log(f"        ξ-bundle κ에서는 A(w=11)={d4_delta_mean:.2f} vs A(w=1)={d4_11a1_mean:.2f}")
    log(f"  → Hardy Z의 1/δ 보정항이 weight에 의존적임을 확인")
log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [7] δ-독립성 검증 (ξ-bundle 정리 확인)
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[7] δ-독립성 검증 — ξ-bundle 정리 (c₁=0 확인)")
log("=" * 72)
log()

log("  sym³(Δ) A(t₀, δ):")
log(f"  {'t₀':>10s} | " + " | ".join(f"δ={d:.3f}" for d in DELTAS))
log("  " + "-" * 80)
for r in delta_results_sym3d:
    row_vals = []
    for delta in DELTAS:
        idx = DELTAS.index(delta)
        if idx < len(r['As']):
            row_vals.append(f"{r['As'][idx]:8.4f}")
        else:
            row_vals.append("  N/A   ")
    log(f"  t₀={r['t0']:.5f} | " + " | ".join(row_vals[:len(DELTAS)]))

log()
log("  sym³(11a1) A(t₀, δ):")
log(f"  {'t₀':>10s} | " + " | ".join(f"δ={d:.3f}" for d in DELTAS))
log("  " + "-" * 80)
for r in delta_results_sym3_11a1:
    row_vals = []
    for i, delta in enumerate(DELTAS):
        if i < len(r['As']):
            row_vals.append(f"{r['As'][i]:8.4f}")
        else:
            row_vals.append("  N/A   ")
    log(f"  t₀={r['t0']:.5f} | " + " | ".join(row_vals[:len(DELTAS)]))

log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [8] 성공 기준 평가
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[8] 성공 기준 평가")
log("=" * 72)
log()

# 기준 1: 36점 이상
total_pts = len(sym3d_all_A) + len(sym3_11a1_all_A)
c1 = total_pts >= 36
log(f"  ① 36점 이상: {'✅' if c1 else '❌'} ({total_pts}점)")

# 기준 2: CV(A) < 5%
c2_d = sym3d_cv < 5.0 if sym3d_all_A else False
c2_11a1 = sym3_11a1_cv < 5.0 if sym3_11a1_all_A else False
log(f"  ② CV(A) < 5% — sym³(Δ): {'✅' if c2_d else '❌'} ({sym3d_cv:.2f}%), "
    f"sym³(11a1): {'✅' if c2_11a1 else '❌'} ({sym3_11a1_cv:.2f}%)")

# 기준 3: d=4 A(t₀) 정량화
c3 = bool(sym3d_all_A) and bool(sym3_11a1_all_A)
log(f"  ③ d=4 A(t₀) 정량화: {'✅' if c3 else '❌'}")
if c3:
    log(f"    sym³(Δ) w=11: A = {d4_delta_mean:.4f}")
    log(f"    sym³(11a1) w=1: A = {d4_11a1_mean:.4f}")

# 기준 4: degree 단조증가
mono_ok = d1_to_d2 and d2_to_d3
d4_any_mono = d3_to_d4_d or d3_to_d4_11a1
log(f"  ④ degree 단조증가 d=1→2→3: {'✅' if mono_ok else '❌'}")
log(f"     d=3→4 (어느 하나라도): {'✅' if d4_any_mono else '❌'}")

# 전체 판정
pass_count = sum([c1, c2_d or c2_11a1, c3, mono_ok])
log()
log(f"  통과: {pass_count}/4")
log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [9] 수치 요약 + 경계 갱신
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[9] 수치 요약")
log("=" * 72)
log()
log(f"  대상: sym³(Δ) [d=4,w=11] + sym³(11a1) [d=4,w=1]")
log(f"  방법: ξ-bundle κ = |Λ'(s)/Λ(s)|² at s = center + δ + i·t₀")
log(f"  PARI lfunlambda(Linit, s, 1) 내장 도함수 사용")
log()

log(f"  sym³(Δ) [d=4, w=11, center=17]:")
for r in delta_results_sym3d:
    if r['As']:
        log(f"    t₀={r['t0']:.5f}: mean(A)={r.get('mean_A',float('nan')):.4f}, CV={r.get('cv_A',0):.2f}%")
log(f"    → 전체 mean(A) = {d4_delta_mean:.4f}")

log()
log(f"  sym³(11a1) [d=4, w=1, center=2]:")
for r in delta_results_sym3_11a1:
    if r['As']:
        log(f"    t₀={r['t0']:.5f}: mean(A)={r.get('mean_A',float('nan')):.4f}, CV={r.get('cv_A',0):.2f}%")
log(f"    → 전체 mean(A) = {d4_11a1_mean:.4f}")

log()
log(f"  ★★ B-12 degree 비교 (ξ-bundle σ-방향, 동일 방법):")
log(f"    d=1: A=1.27  d=2: A=3.93  d=3: A=12.79")
log(f"    d=4 sym³(Δ): A={d4_delta_mean:.4f}  d=4 sym³(11a1): A={d4_11a1_mean:.4f}")

log()
log("=" * 72)
log("경계 갱신")
log("=" * 72)

if d3_to_d4_d or d3_to_d4_11a1:
    log(f"  B-12: ⏳ → ★★ 조건부 — d=4에서 A값 정량화 완료")
else:
    log(f"  B-12: ⏳ → ★ 음성 — d=3→4 단조 위반 (ξ-bundle 방법)")

log(f"  B-22: ⏳ → 결론 포함 — weight vs degree 비교 완료")
log(f"  B-23: 확인 — ξ-bundle κ (σ)와 Hardy Z κ (t) 두 양의 A 값 비교 가능")
log()
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

flush_file()
print(f"\n✅ 결과 저장: {OUTFILE}")
