#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #94b — Hadamard A(t₀)=B²+2H₁ GL(3) sym²(11a1) 재검증
=============================================================================
#94 결과 분석:
  - Re(c₀)≈0: 8/8 ✅
  - Hadamard 수렴 방향: N=5→10→64 단조 ✅ (N=64에서 16-25% 부족)
  - EM 보정 후: 0/8 통과 ❌ (20-55% 과보정)

버그 진단:
  - 기존 #94 스크립트 GL(3) 밀도 공식: N'(T) = (3/π)*log(T/(2π))
  - 올바른 GL(3) 공식:                  N'(T) = (3/(2π))*log(T/(2π))
  - 비율: 2배 과대추정 → B_tail, H₁_tail 2배 과보정

해석학적 검증:
  GL(n) Weyl 밀도 통일 공식:
    GL(1): N'(T) = 1/(2π) * log(T/(2π))
    GL(2): N'(T) = 2/(2π) * log(T/(2π)) = (1/π)*log(T/(2π))
    GL(3): N'(T) = 3/(2π) * log(T/(2π))
  → 계수 d/(2π), 기존 GL(3) 스크립트는 d/π = 2×(d/(2π))를 사용 (오류)

수정 및 추가 개선:
  [1] EM 공식 수정: 3/(2π) 사용
  [2] 영점 확보: T≤150 (N≈200+)  — 수학자 "대안" 권고 반영
  [3] raw Hadamard (EM 없이) 수렴도 확인: N=[50,100,200,N_all]

예측 결과:
  수정된 EM으로 N=64에서도 8/8 <1% 달성 가능
  N=200으로 raw Hadamard 수렴 확인 추가

결과: results/hadamard_gl3_universality_94b.txt
=============================================================================
"""

import sys
import os
import time
import math
import numpy as np

OUTFILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/hadamard_gl3_universality_94b.txt')
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []

def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))

def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

log("=" * 72)
log("결과 #94b — Hadamard A(t₀)=B²+2H₁ GL(3) sym²(11a1) 재검증")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()
log("[수정] EM 밀도 공식 버그 수정: (3/π)→(3/(2π)), 영점 N≥200")
log()
flush_file()

# ━━━━━━━━━━━ PARI 초기화 ━━━━━━━━━━━
log("[Init] PARI 초기화")
import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)
gp("default(realprecision, 100)")
log("  2GB 메모리, realprecision=100 ✅")
log()
flush_file()

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
CENTER       = 1.5
DELTA_DIRECT = 0.01
N_TEST       = 10
T_ZEROS_MAX  = 150.0   # 수학자 권고: N≥200 (T≈150이면 ~200개 기대)

log(f"파라미터: center={CENTER}, δ={DELTA_DIRECT}, N_test={N_TEST}")
log(f"영점 탐색: t ∈ [0, {T_ZEROS_MAX}]  (이전: T=50)")
log()
flush_file()

# ━━━━━━━━━━━ EM 꼬리 보정 (GL(3) 수정버전) ━━━━━━━━━━━
try:
    from scipy import integrate as _sci_int

    def tail_gl3_fixed(t0_val, t_last, mode):
        """
        GL(3) 올바른 Weyl 꼬리 보정:
        N'(T) = (3/(2π)) * log(T/(2π))   ← 수정: 3/π → 3/(2π)

        B_tail  = ∫_{t_last}^∞ [2t₀/(t²-t₀²)] * N'(t) dt
        H₁_tail = ∫_{t_last}^∞ [2/t²]           * N'(t) dt
        """
        c = 3.0 / (2.0 * math.pi)   # 올바른 계수: 3/(2π) ≈ 0.477
        c2pi = 2.0 * math.pi
        if mode == 'B':
            fn = lambda t: 2.0*t0_val/(t**2 - t0_val**2) * c * math.log(t/c2pi) if t > t0_val+1e-6 else 0.0
        else:  # H1
            fn = lambda t: 2.0/t**2 * c * math.log(t/c2pi)
        v, _ = _sci_int.quad(fn, t_last, 1e7, limit=500, epsabs=1e-12)
        return v

    def tail_gl3_buggy(t0_val, t_last, mode):
        """
        #94 기존 버그 버전: N'(T) = (3/π)*log(T/(2π)) — 2× 과대
        비교용으로만 사용
        """
        c = 3.0 / math.pi   # 버그: 3/π (2배 과대)
        c2pi = 2.0 * math.pi
        if mode == 'B':
            fn = lambda t: 2.0*t0_val/(t**2 - t0_val**2) * c * math.log(t/c2pi) if t > t0_val+1e-6 else 0.0
        else:
            fn = lambda t: 2.0/t**2 * c * math.log(t/c2pi)
        v, _ = _sci_int.quad(fn, t_last, 1e7, limit=500, epsabs=1e-12)
        return v

    HAS_TAIL = True
    log("[EM tail] scipy 사용 가능 ✅")
    log("  [수정] 밀도 공식: N'(T) = 3/(2π) * log(T/(2π))")
    log("  비교용 [버그]: N'(T) = 3/π * log(T/(2π))  (이전 #94)")
    log()
    log("  GL(n) Weyl 밀도 통일 공식 검증:")
    log("    GL(1): 1/(2π) * log = 0.159 * log")
    log("    GL(2): 2/(2π) * log = 0.318 * log  [#92에서 확인]")
    log("    GL(3): 3/(2π) * log = 0.477 * log  [이번 수정]")
    log("  #94 기존: 3/π * log = 0.955 * log  ← 2× 과대 → 이번 수정")
except ImportError:
    HAS_TAIL = False
    def tail_gl3_fixed(t0_val, t_last, mode): return 0.0
    def tail_gl3_buggy(t0_val, t_last, mode): return 0.0
    log("[EM tail] scipy 없음")
log()
flush_file()

# ━━━━━━━━━━━ Step 1: sym²(11a1) 영점 확보 (N≥200 목표) ━━━━━━━━━━━
log("=" * 72)
log(f"[Step 1] PARI — sym²(11a1) 영점 확보 (T≤{T_ZEROS_MAX})")
log("=" * 72)
t0_s = time.time()

gp("E11 = ellinit([0,-1,1,-10,-20])")   # 11a1
gp("L11 = lfunsympow(E11, 2)")

k11_raw = str(gp("L11[4]"))
try:
    k11 = int(round(float(k11_raw)))
except Exception:
    k11 = 3
log(f"  L11[4] = {k11_raw} → k={k11}, center={k11/2.0:.1f}")

# lfuninit with extended range
log(f"  lfuninit([0, {T_ZEROS_MAX}]) 시작...")
t_init = time.time()
gp(f"L11i = lfuninit(L11, [0, {T_ZEROS_MAX}])")
log(f"  lfuninit 완료 ({time.time()-t_init:.1f}s)")

# lfunzeros
log(f"  lfunzeros(L11i, {T_ZEROS_MAX}) 시작...")
t_zv = time.time()
gp(f"zvec = lfunzeros(L11i, {T_ZEROS_MAX})")
n_zeros_raw = int(gp("length(zvec)"))
log(f"  lfunzeros 완료: {n_zeros_raw}개 ({time.time()-t_zv:.1f}s)")

all_zeros_list = []
for i in range(1, n_zeros_raw + 1):
    try:
        z = float(gp(f"zvec[{i}]"))
        all_zeros_list.append(z)
    except Exception as e:
        log(f"  WARNING 영점 {i}: {e}")

all_zeros = np.array(sorted(all_zeros_list))
N_ALL = len(all_zeros)
t_last = float(all_zeros[-1]) if N_ALL > 0 else 0.0

log(f"  최종: N={N_ALL}개, t ∈ [{all_zeros[0]:.4f}, {t_last:.4f}]")
log(f"  처음 5개: {[f'{z:.4f}' for z in all_zeros[:5]]}")
log(f"  밀도: {N_ALL/(t_last-all_zeros[0]):.3f} zeros/unit")

if N_ALL < 10:
    log("⚠️ 영점 10개 미만 — 중단")
    flush_file()
    sys.exit(1)

n_prev = 64  # 이전 N=64 대비
log(f"  이전 #94: N=64 (T=50). 이번: N={N_ALL} (T={T_ZEROS_MAX}). 증가: {N_ALL}×")
log(f"  전체 소요: {time.time()-t0_s:.1f}s")
log()
flush_file()

# ━━━━━━━━━━━ Hadamard 함수 ━━━━━━━━━━━
def hadamard_B_H1(t0, tzeros):
    """
    B(t₀) = -1/(2t₀) + Σ_{n≠0}[-1/(t₀-t_n) - 1/(t₀+t_n)]
    H₁(t₀) = 1/(4t₀²) + Σ_{n≠0}[1/(t₀-t_n)² + 1/(t₀+t_n)²]
    """
    mask = np.abs(tzeros - t0) > 1e-6
    tv   = tzeros[mask]
    d    = t0 - tv
    s    = t0 + tv
    ok   = (np.abs(d) > 1e-12) & (np.abs(s) > 1e-12)
    d, s = d[ok], s[ok]
    B    = -1.0/(2.0*t0) + np.sum(-1.0/d - 1.0/s)
    H1   = 1.0/(4.0*t0**2) + np.sum(1.0/d**2 + 1.0/s**2)
    return B, H1, B**2 + 2.0*H1

# ━━━━━━━━━━━ Step 2: 테스트 영점별 계산 ━━━━━━━━━━━
log("=" * 72)
log(f"[Step 2] 테스트 영점별 A(t₀) 계산")
log("=" * 72)

test_zeros = all_zeros[:N_TEST]
N_HAD_LIST = sorted(set([50, 100, N_ALL] + ([200] if N_ALL >= 200 else [])))

log(f"  테스트 영점: t = {test_zeros[0]:.4f} ~ {test_zeros[-1]:.4f}")
log(f"  Hadamard N 목록: {N_HAD_LIST}")
log(f"  EM 꼬리: t_last={t_last:.2f}")
log()
flush_file()

results = []
for idx, t0_val in enumerate(test_zeros):
    log(f"  ── 영점 #{idx+1}: t₀={t0_val:.6f} ──")
    tz = time.time()

    # ── A_direct: PARI lfunlambda ──
    sigma_near = CENTER + DELTA_DIRECT
    try:
        gp(f"snear = {sigma_near:.10f} + I*{t0_val:.10f}")
        gp(f"Lnear = lfunlambda(L11i, snear)")
        gp(f"Lpnear = lfunlambda(L11i, snear, 1)")
        abs_L0 = float(gp("abs(Lnear)"))
        abs_Lp = float(gp("abs(Lpnear)"))
        if abs_L0 < 1e-200:
            log(f"    ⚠️ |Λ| 너무 작음 — skip")
            results.append({'t': t0_val, 'skip': True})
            log(); flush_file()
            continue
        kappa  = (abs_Lp / abs_L0)**2
        A_dir  = float(kappa - 1.0/DELTA_DIRECT**2)
    except Exception as e:
        log(f"    WARNING A_direct: {e}")
        results.append({'t': t0_val, 'skip': True})
        log(); flush_file()
        continue

    log(f"    A_direct = {A_dir:.6f}  ({time.time()-tz:.1f}s)")

    # ── Re(c₀) 근사: Re(Λ'/Λ(center+δ+it₀)) - 1/δ ──
    try:
        re_conn  = float(gp("real(Lpnear / Lnear)"))
        Re_c0_approx = re_conn - 1.0/DELTA_DIRECT
    except Exception:
        Re_c0_approx = float('nan')
    log(f"    Re(c₀)≈ {Re_c0_approx:.4e}  {'✅<0.05' if np.isfinite(Re_c0_approx) and abs(Re_c0_approx)<0.05 else '⚠️'}")

    # ── Hadamard 수렴 추세 ──
    A_had_N = {}
    for Nv in N_HAD_LIST:
        Nv_eff = min(Nv, N_ALL)
        B_n, H1_n, A_n = hadamard_B_H1(t0_val, all_zeros[:Nv_eff])
        A_had_N[Nv] = (B_n, H1_n, A_n)
    log(f"    Had: " + ", ".join([f"N={Nv}→{A_had_N[Nv][2]:.4f}" for Nv in N_HAD_LIST]))

    # ── EM 꼬리 보정 (수정버전) ──
    B_had, H1_had, A_had_max = A_had_N[N_HAD_LIST[-1]]
    Bt_fix = H1t_fix = 0.0
    Bt_bug = H1t_bug = 0.0
    if HAS_TAIL and t_last > t0_val + 1.0:
        try:
            Bt_fix  = tail_gl3_fixed(t0_val, t_last, 'B')
            H1t_fix = tail_gl3_fixed(t0_val, t_last, 'H1')
            Bt_bug  = tail_gl3_buggy(t0_val, t_last, 'B')
            H1t_bug = tail_gl3_buggy(t0_val, t_last, 'H1')
        except Exception as e:
            log(f"    WARNING EM: {e}")

    A_corr_fix = (B_had + Bt_fix)**2 + 2.0*(H1_had + H1t_fix)
    A_corr_bug = (B_had + Bt_bug)**2 + 2.0*(H1_had + H1t_bug)
    err_fix    = abs(A_corr_fix - A_dir)/(abs(A_dir) + 1e-12)
    err_bug    = abs(A_corr_bug - A_dir)/(abs(A_dir) + 1e-12)
    err_raw    = abs(A_had_max  - A_dir)/(abs(A_dir) + 1e-12)   # EM 없이

    log(f"    EM [수정] B_tail={Bt_fix:.5f}, H₁_tail={H1t_fix:.5f}")
    log(f"    A_corr[수정] = {A_corr_fix:.6f}  err={err_fix*100:.2f}%  "
        f"{'✓<1%' if err_fix<0.01 else ('✓<2%' if err_fix<0.02 else '✗')}")
    log(f"    EM [버그]  B_tail={Bt_bug:.5f}, H₁_tail={H1t_bug:.5f}")
    log(f"    A_corr[버그]  = {A_corr_bug:.6f}  err={err_bug*100:.2f}%  "
        f"(#94 재현)")
    log(f"    Raw Hadamard  = {A_had_max:.6f}  err={err_raw*100:.2f}%  "
        f"(EM 없이)")
    log(f"    소요: {time.time()-tz:.1f}s")

    results.append({
        't': t0_val, 'A_direct': A_dir,
        'B_had': B_had, 'H1_had': H1_had, 'A_had': A_had_max,
        'Bt_fix': Bt_fix, 'H1t_fix': H1t_fix,
        'Bt_bug': Bt_bug, 'H1t_bug': H1t_bug,
        'A_corr_fix': A_corr_fix, 'err_fix': err_fix,
        'A_corr_bug': A_corr_bug, 'err_bug': err_bug,
        'err_raw': err_raw,
        'Re_c0': Re_c0_approx, 'A_had_N': A_had_N,
        'skip': False
    })
    log(); flush_file()

# ━━━━━━━━━━━ Step 3: N 수렴 추세 (EM 없이) ━━━━━━━━━━━
log("=" * 72)
log("[Step 3] N 수렴 추세: raw Hadamard (EM 없이)")
log("=" * 72)
valid = [r for r in results if not r['skip']]
for idx in [0, 4, min(9, len(valid)-1)]:
    if idx >= len(valid): continue
    r = valid[idx]
    log(f"  t₀={r['t']:.4f}: A_direct={r['A_direct']:.4f}")
    for Nv in sorted(r['A_had_N'].keys()):
        An = r['A_had_N'][Nv][2]
        dp = abs(An - r['A_direct'])/(abs(r['A_direct'])+1e-12)*100
        log(f"    N={Nv:>4}: A_had={An:.4f}, err={dp:.1f}%")
    log(f"    +EM[수정]: A_corr={r['A_corr_fix']:.4f}, err={r['err_fix']*100:.1f}%")
    log(f"    +EM[버그]:  A_corr={r['A_corr_bug']:.4f}, err={r['err_bug']*100:.1f}%")
log()
flush_file()

# ━━━━━━━━━━━ Step 4: B-20 패턴 ━━━━━━━━━━━
log("=" * 72)
log("[Step 4] B-20 패턴 (근접 영점 ↔ 큰 A 상관)")
log("=" * 72)
if valid:
    Apairs = sorted([(r['A_direct'], r) for r in valid], reverse=True)
    for lab, r in [("최대 A", Apairs[0][1]), ("최소 A", Apairs[-1][1])]:
        t_x = r['t']
        tv   = all_zeros[np.abs(all_zeros - t_x) > 1e-6]
        near = tv[np.abs(tv - t_x) < 5.0]
        H1n  = float(np.sum(1.0/(t_x-near)**2 + 1.0/(t_x+near)**2)) if len(near) > 0 else 0.0
        log(f"  {lab}: t₀={t_x:.4f}, A_direct={r['A_direct']:.4f}, "
            f"근접(Δt<5) {len(near)}개, H₁_near={H1n:.4f}")
        for tn in near[:5]:  # 최대 5개 표시
            log(f"    t_n={tn:.4f}, Δt={abs(t_x-tn):.4f}")
log()
flush_file()

# ━━━━━━━━━━━ 최종 판정 ━━━━━━━━━━━
log("=" * 72)
log("최종 판정 — Hadamard GL(3) 수정 EM")
log("=" * 72)
log()

nv    = len(valid)
n_fix_1pct = sum(1 for r in valid if r['err_fix'] < 0.01)
n_fix_2pct = sum(1 for r in valid if r['err_fix'] < 0.02)
n_bug_2pct = sum(1 for r in valid if r['err_bug'] < 0.02)
n_raw_2pct = sum(1 for r in valid if r['err_raw'] < 0.02)
n_rec0     = sum(1 for r in valid
                 if np.isfinite(r.get('Re_c0', float('nan')))
                 and abs(r['Re_c0']) < 0.05)

log(f"{'#':>3} {'t₀':>10} {'A_direct':>10} {'A_corr[수정]':>14} {'err_fix%':>9} {'err_raw%':>9}  pass")
log("-" * 75)
for i, r in enumerate(results):
    if r['skip']:
        log(f"  {i+1:>2} {r['t']:>10.4f}  SKIP")
        continue
    ok = '✓<1%' if r['err_fix'] < 0.01 else ('✓<2%' if r['err_fix'] < 0.02 else '✗')
    log(f"  {i+1:>2} {r['t']:>10.4f} {r['A_direct']:>10.4f} {r['A_corr_fix']:>14.4f} "
        f"{r['err_fix']*100:>9.2f} {r['err_raw']*100:>9.2f}  {ok}")
log()

log(f"판정 요약 (N={N_ALL}개 영점, T≤{T_ZEROS_MAX}):")
log(f"  [수정 EM] <1%:  {n_fix_1pct}/{nv}개  <2%: {n_fix_2pct}/{nv}개")
log(f"  [버그  EM] <2%:  {n_bug_2pct}/{nv}개  (이전 #94 재현)")
log(f"  [raw  Had] <2%:  {n_raw_2pct}/{nv}개  (EM 없이)")
log(f"  Re(c₀)<0.05: {n_rec0}/{nv}개  (FE 대칭성)")
log()

# 성공 기준
ok_strict  = (nv >= 10) and (n_fix_1pct >= 7)
ok_crit    = (nv >= 10) and (n_fix_2pct >= 7)
ok_rec0    = (n_rec0 >= nv - 1)

log("성공 기준:")
log(f"  영점 10개 이상: {'✅' if nv >= 10 else '❌'}")
log(f"  수정EM <1%:  {n_fix_1pct}/10  {'✅' if n_fix_1pct>=7 else '❌'} (기준 ≥7/10)")
log(f"  수정EM <2%:  {n_fix_2pct}/10  {'✅' if n_fix_2pct>=7 else '❌'}")
log(f"  Re(c₀)<0.05: {n_rec0}/{nv}  {'✅' if ok_rec0 else '⚠️'}")
log()

if ok_strict and ok_rec0:
    log("  ★★★ 강양성 — GL(3) Hadamard A=B²+2H₁ 1% 기준 통과 (수정 EM)")
    log("  → degree-보편성 d=1(★★★), d=2(★★), d=3(★★★) 수립")
    log(f"  → 버그 수정 핵심: EM 밀도 계수 3/π → 3/(2π) (d/(2π) 통일 공식)")
elif ok_crit and ok_rec0:
    log("  ★★ 조건부 양성 — GL(3) Hadamard 2% 기준 통과 (수정 EM)")
    log("  → degree-보편성 d=1(★★★), d=2(★★), d=3(★★) 수립")
elif ok_crit:
    log("  ★ 약양성 — 2% 기준 통과, Re(c₀) 일부 큼")
else:
    log("  ⚠️ 미결 — 추가 분석 필요")

log()
log("도표 (degree 비교, 수정 후):")
log(f"  d=1 GL(1) ζ       : 13/13 <0.003%  ★★★  (N=2000, #90)")
log(f"  d=2 GL(2) L(s,Δ)  : 7/10  <1.3%    ★★   (N=86, #92)")
log(f"  d=3 GL(3) sym²    : {n_fix_1pct}/10  <1%, {n_fix_2pct}/10 <2%  (N={N_ALL}, #94b)")
log()
log("[EM 버그 분석]")
log("  기존 #94: N'(T) = (3/π)*log(T/(2π)) — GL(n) 공식의 분모를 π로 잘못 설정")
log("  올바른:   N'(T) = (3/(2π))*log(T/(2π)) — d/(2π) 통일 공식")
log("  영향: B_tail 2× 과대추정 → A_corr > A_direct (20-55% 과보정)")
log(f"  수정 후: err ≈ {'<1%' if ok_strict else '<2%' if ok_crit else '>2%'}")
log()
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"결과: {OUTFILE}")
log("=" * 72)
flush_file()
