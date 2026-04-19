#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #94c — Hadamard GL(3) sym²(11a1): Re(c₀) 보정 재검증
=============================================================================
#94b 분석:
  - 수정 EM (3/(2π)): 2/10 <2% (개선되었지만 부족)
  - 주요 잔존 오차: A_direct 자체의 Re(c₀)/δ 오염

수학적 진단:
  A_direct = κ(ρ₀+δ) - 1/δ²
           = 2·Re(c₀)/δ + |c₀|² + 2·Re(c₁) + O(δ)
           = 2·Re(c₀)/δ + B² + 2H₁                [이론값]

  Re(c₀) ≠ 0 (수치적으로 ~1-3e-3) 이면:
    A_direct ≠ B² + 2H₁  — 2·Re(c₀)/δ 오염 ≈ 0.2-0.6 (A의 2-4%)

  수정: A_true = A_direct - 2·Re(c₀)_rich/δ  ≈ B² + 2H₁
  → Hadamard+EM vs A_true 비교: 오차 대폭 감소 예상

접근법:
  [단계 1] Richardson으로 Re(c₀) 정밀 추출
  [단계 2] A_true = A_direct - 2·Re(c₀)/δ
  [단계 3] A_corr (N=271, 수정 EM) vs A_true 비교

이론:
  - A_true = A_direct - 2Re(c₀)/δ ≈ B² + 2H₁  (진짜 A)
  - A_corr = Hadamard+EM ≈ B² + 2H₁           (Hadamard 추정)
  - 두 추정이 일치하면 A = B² + 2H₁ 수치 확인

결과: results/hadamard_gl3_universality_94c.txt
=============================================================================
"""

import sys
import os
import time
import math
import numpy as np

OUTFILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/hadamard_gl3_universality_94c.txt')
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))
def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

log("=" * 72)
log("결과 #94c — Hadamard GL(3) sym²(11a1): Re(c₀) 보정 재검증")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()
log("핵심 수정: A_true = A_direct - 2·Re(c₀)/δ  (Re(c₀) Richardson 보정)")
log()
flush_file()

# ━━━━━━━━━━━ PARI 초기화 ━━━━━━━━━━━
import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)
gp("default(realprecision, 100)")
log("[Init] PARI: 2GB, prec=100 ✅")
log()
flush_file()

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
CENTER       = 1.5
DELTA_DIRECT = 0.01
N_TEST       = 10
T_ZEROS_MAX  = 150.0

log(f"파라미터: center={CENTER}, δ={DELTA_DIRECT}, N_test={N_TEST}")
log()
flush_file()

# ━━━━━━━━━━━ EM 꼬리 보정 (GL(3) 수정 공식) ━━━━━━━━━━━
try:
    from scipy import integrate as _sci_int

    def tail_gl3(t0_val, t_last, mode):
        """N'(T) = 3/(2π) * log(T/(2π))"""
        c = 3.0 / (2.0 * math.pi)
        c2pi = 2.0 * math.pi
        if mode == 'B':
            fn = lambda t: 2.0*t0_val/(t**2 - t0_val**2) * c * math.log(t/c2pi) if t > t0_val+1e-6 else 0.0
        else:
            fn = lambda t: 2.0/t**2 * c * math.log(t/c2pi)
        v, _ = _sci_int.quad(fn, t_last, 1e7, limit=500, epsabs=1e-12)
        return v

    HAS_TAIL = True
    log("[EM tail] scipy ✅, 밀도: 3/(2π)*log(T/(2π))")
except ImportError:
    HAS_TAIL = False
    def tail_gl3(t0_val, t_last, mode): return 0.0
    log("[EM tail] scipy 없음")
log()
flush_file()

# ━━━━━━━━━━━ Step 1: 영점 확보 ━━━━━━━━━━━
log("=" * 72)
log(f"[Step 1] sym²(11a1) 영점 확보 (T≤{T_ZEROS_MAX})")
log("=" * 72)
t0_s = time.time()

gp("E11 = ellinit([0,-1,1,-10,-20])")
gp("L11 = lfunsympow(E11, 2)")
log(f"  lfuninit([0, {T_ZEROS_MAX}]) ...")
gp(f"L11i = lfuninit(L11, [0, {T_ZEROS_MAX}])")
gp(f"zvec = lfunzeros(L11i, {T_ZEROS_MAX})")
n_zeros_raw = int(gp("length(zvec)"))

all_zeros_list = []
for i in range(1, n_zeros_raw + 1):
    try:
        all_zeros_list.append(float(gp(f"zvec[{i}]")))
    except Exception:
        pass

all_zeros = np.array(sorted(all_zeros_list))
N_ALL = len(all_zeros)
t_last = float(all_zeros[-1]) if N_ALL > 0 else 0.0

log(f"  N={N_ALL}개, t ∈ [{all_zeros[0]:.4f}, {t_last:.4f}]")
log(f"  소요: {time.time()-t0_s:.1f}s")

if N_ALL < 10:
    log("⚠️ 영점 10개 미만 — 중단"); flush_file(); sys.exit(1)
log()
flush_file()

# ━━━━━━━━━━━ Hadamard 함수 ━━━━━━━━━━━
def hadamard_B_H1(t0, tzeros):
    mask = np.abs(tzeros - t0) > 1e-6
    tv   = tzeros[mask]
    d    = t0 - tv
    s    = t0 + tv
    ok   = (np.abs(d) > 1e-12) & (np.abs(s) > 1e-12)
    d, s = d[ok], s[ok]
    B    = -1.0/(2.0*t0) + np.sum(-1.0/d - 1.0/s)
    H1   = 1.0/(4.0*t0**2) + np.sum(1.0/d**2 + 1.0/s**2)
    return B, H1, B**2 + 2.0*H1

# ━━━━━━━━━━━ Step 2: Richardson Re(c₀) 계산 함수 ━━━━━━━━━━━
def compute_Re_c0_richardson(t0_val, delta=DELTA_DIRECT):
    """
    c₀ = [Λ'/Λ(ρ₀+ε) - 1/ε] as ε→0  (Richardson 외삽)
    Re(c₀)가 0에 가까운지 확인 (FE 대칭성)
    Return: A_direct, Re_c0, 성공 여부
    """
    rho0_re = CENTER
    rho0_im = t0_val

    # A_direct: δ=0.01 사용
    try:
        gp(f"snear = {CENTER+delta:.10f} + I*{t0_val:.10f}")
        gp(f"Lnear = lfunlambda(L11i, snear)")
        gp(f"Lpnear = lfunlambda(L11i, snear, 1)")
        abs_L0 = float(gp("abs(Lnear)"))
        abs_Lp = float(gp("abs(Lpnear)"))
        if abs_L0 < 1e-200:
            return float('nan'), float('nan'), False
        kappa = (abs_Lp/abs_L0)**2
        A_dir = float(kappa - 1.0/delta**2)
    except Exception as e:
        log(f"    WARNING A_direct: {e}")
        return float('nan'), float('nan'), False

    # Richardson Re(c₀): ε=[1e-3, 5e-4, 2.5e-4, 1.25e-4]
    eps_list = [1e-3, 5e-4, 2.5e-4, 1.25e-4]
    c0_raw = []
    for eps_f in eps_list:
        try:
            gp(f"srich = {CENTER:.6f} + {eps_f:.10f} + I*{t0_val:.10f}")
            gp(f"Lrich = lfunlambda(L11i, srich)")
            gp(f"Lprich = lfunlambda(L11i, srich, 1)")
            abs_Lr = float(gp("abs(Lrich)"))
            if abs_Lr < 1e-200:
                c0_raw.append(float('nan'))
                continue
            conn_re = float(gp("real(Lprich / Lrich)"))
            conn_im = float(gp("imag(Lprich / Lrich)"))
            # c₀_raw = Λ'/Λ(ρ₀+ε) - 1/ε  (ε is real, so 1/ε subtracts from Re part)
            c0_raw_val = conn_re - 1.0/eps_f + 1j*(conn_im)
            c0_raw.append(c0_raw_val)
        except Exception as e:
            log(f"    WARNING Richardson eps={eps_f}: {e}")
            c0_raw.append(float('nan'))

    # Richardson 외삽 (같은 구조: #92 참조)
    valid_c0 = [(i, v) for i, v in enumerate(c0_raw) if np.isfinite(np.real(v))]
    if len(valid_c0) < 3:
        return A_dir, float('nan'), True

    c0_arr = np.array([v for _, v in valid_c0])
    r1 = [(4*c0_arr[i+1] - c0_arr[i])/3 for i in range(len(c0_arr)-1)]
    if len(r1) >= 2:
        r2 = [(4*r1[i+1] - r1[i])/3 for i in range(len(r1)-1)]
        c0_rich = r2[0] if len(r2) == 1 else (4*r2[-1] - r2[-2])/3 if len(r2) >= 2 else r1[0]
    else:
        c0_rich = r1[0]

    Re_c0 = float(np.real(c0_rich))
    return A_dir, Re_c0, True

# ━━━━━━━━━━━ Step 3: 테스트 영점별 계산 ━━━━━━━━━━━
log("=" * 72)
log(f"[Step 3] 테스트 영점별 계산 (Re(c₀) Richardson + Hadamard+EM)")
log("=" * 72)

test_zeros = all_zeros[:N_TEST]
log(f"  테스트 영점: t = {test_zeros[0]:.4f} ~ {test_zeros[-1]:.4f}")
log(f"  N_Had = N_ALL = {N_ALL}개")
log()
flush_file()

results = []
for idx, t0_val in enumerate(test_zeros):
    log(f"  ── 영점 #{idx+1}: t₀={t0_val:.6f} ──")
    tz = time.time()

    # A_direct + Re(c₀)
    A_dir, Re_c0, success = compute_Re_c0_richardson(t0_val)
    log(f"    A_direct = {A_dir:.6f}")
    log(f"    Re(c₀)_rich = {Re_c0:.4e}  {'✅' if np.isfinite(Re_c0) and abs(Re_c0)<0.01 else '⚠️'}")

    if not np.isfinite(A_dir):
        log(f"    ⚠️ skip")
        results.append({'t': t0_val, 'skip': True})
        log(); flush_file()
        continue

    # A_true = A_direct - 2*Re(c₀)/δ  (Re(c₀) 보정)
    if np.isfinite(Re_c0):
        A_true = A_dir - 2.0 * Re_c0 / DELTA_DIRECT
        correction = 2.0 * Re_c0 / DELTA_DIRECT
    else:
        A_true = A_dir   # 보정 불가 시 A_direct 그대로 사용
        correction = 0.0
    log(f"    A_true = {A_true:.6f}  (보정량: {correction:.4f} = 2·Re(c₀)/δ)")

    # Hadamard sum (all N_ALL zeros)
    B_had, H1_had, A_had = hadamard_B_H1(t0_val, all_zeros)
    log(f"    B_had = {B_had:.6f}, H₁_had = {H1_had:.6f}")
    log(f"    A_had (N={N_ALL}) = {A_had:.6f}")

    # EM 꼬리 보정 (수정 공식)
    Bt = H1t = 0.0
    if HAS_TAIL and t_last > t0_val + 1.0:
        try:
            Bt  = tail_gl3(t0_val, t_last, 'B')
            H1t = tail_gl3(t0_val, t_last, 'H1')
        except Exception as e:
            log(f"    WARNING EM: {e}")

    A_corr = (B_had + Bt)**2 + 2.0*(H1_had + H1t)
    log(f"    EM 보정: B_tail={Bt:.5f}, H₁_tail={H1t:.5f}")
    log(f"    A_corr = {A_corr:.6f}")

    # 비교 1: A_corr vs A_direct (이전 방식)
    err_old = abs(A_corr - A_dir) / (abs(A_dir) + 1e-12)
    # 비교 2: A_corr vs A_true (Re(c₀) 보정 후)
    err_new = abs(A_corr - A_true) / (abs(A_true) + 1e-12)

    log(f"    err vs A_direct: {err_old*100:.2f}%  (이전 기준)")
    log(f"    err vs A_true:   {err_new*100:.2f}%  (Re(c₀) 보정 후) "
        f"{'✓<1%' if err_new<0.01 else ('✓<2%' if err_new<0.02 else '✗')}")
    log(f"    소요: {time.time()-tz:.1f}s")

    results.append({
        't': t0_val, 'A_direct': A_dir, 'A_true': A_true,
        'Re_c0': Re_c0, 'correction': correction,
        'B_had': B_had, 'H1_had': H1_had, 'A_had': A_had,
        'Bt': Bt, 'H1t': H1t, 'A_corr': A_corr,
        'err_old': err_old, 'err_new': err_new,
        'skip': False
    })
    log(); flush_file()

# ━━━━━━━━━━━ Step 4: B-20 패턴 ━━━━━━━━━━━
log("=" * 72)
log("[Step 4] B-20 패턴 (근접 영점 ↔ 큰 A 상관)")
log("=" * 72)
valid = [r for r in results if not r['skip']]
if valid:
    Apairs = sorted([(r['A_true'], r) for r in valid], reverse=True)
    for lab, r in [("최대 A_true", Apairs[0][1]), ("최소 A_true", Apairs[-1][1])]:
        t_x = r['t']
        tv   = all_zeros[np.abs(all_zeros - t_x) > 1e-6]
        near = tv[np.abs(tv - t_x) < 5.0]
        H1n  = float(np.sum(1.0/(t_x-near)**2 + 1.0/(t_x+near)**2)) if len(near) > 0 else 0.0
        log(f"  {lab}: t₀={t_x:.4f}, A_true={r['A_true']:.4f}, "
            f"근접 {len(near)}개, H₁_near={H1n:.4f}")
        for tn in near[:5]:
            log(f"    t_n={tn:.4f}, Δt={abs(t_x-tn):.4f}")
log()
flush_file()

# ━━━━━━━━━━━ 최종 판정 ━━━━━━━━━━━
log("=" * 72)
log("최종 판정 — Hadamard GL(3) Re(c₀) 보정 재검증")
log("=" * 72)
log()

nv      = len(valid)
n_old_1 = sum(1 for r in valid if r['err_old'] < 0.01)
n_old_2 = sum(1 for r in valid if r['err_old'] < 0.02)
n_new_1 = sum(1 for r in valid if r['err_new'] < 0.01)
n_new_2 = sum(1 for r in valid if r['err_new'] < 0.02)
n_rec0  = sum(1 for r in valid
              if np.isfinite(r.get('Re_c0', float('nan'))) and abs(r['Re_c0']) < 0.01)

log(f"{'#':>3} {'t₀':>10} {'A_direct':>10} {'A_true':>10} {'A_corr':>10} {'err_old%':>9} {'err_new%':>9}  pass")
log("-" * 80)
for i, r in enumerate(results):
    if r['skip']:
        log(f"  {i+1:>2} {r['t']:>10.4f}  SKIP")
        continue
    ok = '✓<1%' if r['err_new'] < 0.01 else ('✓<2%' if r['err_new'] < 0.02 else '✗')
    log(f"  {i+1:>2} {r['t']:>10.4f} {r['A_direct']:>10.4f} {r['A_true']:>10.4f} "
        f"{r['A_corr']:>10.4f} {r['err_old']*100:>9.2f} {r['err_new']*100:>9.2f}  {ok}")
log()

log(f"판정 요약 (N={N_ALL}개 영점, T≤{T_ZEROS_MAX}):")
log(f"  [A_direct 기준, 수정 EM] <1%: {n_old_1}/{nv}  <2%: {n_old_2}/{nv}")
log(f"  [A_true  기준, 수정 EM] <1%: {n_new_1}/{nv}  <2%: {n_new_2}/{nv}")
log(f"  Re(c₀)<0.01: {n_rec0}/{nv}  (FE 대칭성)")
log()

ok_new_strict = (nv >= 10) and (n_new_1 >= 7)
ok_new_crit   = (nv >= 10) and (n_new_2 >= 7)
ok_rec0       = (n_rec0 >= nv - 2)

log("성공 기준 (A_true 기준):")
log(f"  영점 10개: {'✅' if nv >= 10 else '❌'}")
log(f"  <1% : {n_new_1}/10  {'✅' if n_new_1>=7 else '❌'}")
log(f"  <2% : {n_new_2}/10  {'✅' if n_new_2>=7 else '❌'}")
log(f"  Re(c₀)<0.01: {n_rec0}/{nv}  {'✅' if ok_rec0 else '⚠️'}")
log()

if ok_new_strict and ok_rec0:
    log("  ★★★ 강양성 — A=B²+2H₁ GL(3) 엄격 1% 기준 통과 (Re(c₀) 보정 후)")
    log("  → degree-보편성 d=1(★★★), d=2(★★), d=3(★★★)")
elif ok_new_crit and ok_rec0:
    log("  ★★ 조건부 양성 — A=B²+2H₁ GL(3) 2% 기준 통과 (Re(c₀) 보정 후)")
    log("  → degree-보편성 d=1(★★★), d=2(★★), d=3(★★)")
elif ok_new_crit:
    log("  ★ 약양성 — 2% 기준 통과, Re(c₀) 일부 큼")
else:
    log("  ⚠️ 미결 — 추가 분석 필요")

log()
log("[Re(c₀) 보정 방법론]")
log("  A = κ(ρ₀+δ) - 1/δ² = 2·Re(c₀)/δ + B² + 2H₁ + O(δ)")
log("  Re(c₀) ≠ 0 (수치적): 2·Re(c₀)/δ ≈ 0.2-0.6 (A의 2-5% 오염)")
log("  A_true = A_direct - 2·Re(c₀)_rich/δ ≈ B² + 2H₁")
log("  이는 A의 두 독립 추정: (1) Richardson, (2) Hadamard+EM")
log("  두 추정이 일치 → A = B² + 2H₁ 수치 확인")
log()
log("도표 (degree 비교, Re(c₀) 보정 포함):")
log(f"  d=1 GL(1) ζ      : 13/13 <0.003%  ★★★  (N=2000, #90)")
log(f"  d=2 GL(2) L(s,Δ) : 7/10  <1.3%    ★★   (N=86, #92)")
log(f"  d=3 GL(3) sym²   : {n_new_1}/10 <1%, {n_new_2}/10 <2%  (N={N_ALL}, #94c, Re(c₀) 보정)")
log()
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"결과: {OUTFILE}")
log("=" * 72)
flush_file()
