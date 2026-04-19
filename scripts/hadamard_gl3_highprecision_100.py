#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #100 — Hadamard GL(3) sym²(11a1): 고정밀 재검증 (T=500)
=============================================================================
#94c 분석:
  - N=271 영점 (T≤150), 3/10 <2% — EM tail 기여 3-5%가 주요 오차 원인
  - t₀=6~12 영역에서 tail 기여가 크기 때문에 수렴 부족

#100 목적:
  - T=150 → T=500: N~1000 영점 확보
  - EM tail 기여 ~3× 감소 예상 (tail ∝ 1/N)
  - 동일 10개 테스트 영점(t₁~t₁₀)에서 #94c와 직접 비교
  - 성공 기준: 7/10 이상 err<2% 또는 5/10 이상 err<1%

수학적 접근:
  A_direct = κ(ρ₀+δ) - 1/δ²
           = 2·Re(c₀)/δ + B² + 2H₁ + O(δ)

  Re(c₀) Richardson 보정: A_true = A_direct - 2·Re(c₀)/δ
  Hadamard+EM: A_corr = (B_had + B_tail)² + 2(H₁_had + H₁_tail)

  → A_corr vs A_true: 두 독립 추정 비교

EM 꼬리 밀도 (GL(3) = 수정):
  N'(T) = d/(2π)·log(T/(2π)) = 3/(2π)·log(T/(2π))

파라미터:
  - center=1.5 (motivic, k=3)
  - δ=0.01 (DELTA_DIRECT)
  - T=500.0
  - PARI 메모리: 3GB (T=500에서 증가)
  - 동일 10개 테스트 영점: #94c t₁~t₁₀ (3.8993~13.3353)

결과: results/hadamard_gl3_highprecision_100.txt
=============================================================================
"""

import sys
import os
import time
import math
import numpy as np

OUTFILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/hadamard_gl3_highprecision_100.txt')
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))
def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

log("=" * 72)
log("결과 #100 — Hadamard GL(3) sym²(11a1): 고정밀 재검증 (T=500)")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()
log("목적: #94c (N=271, 3/10<2%) → T=500(N~1000)으로 EM tail ~3× 감소")
log("핵심 수정: A_true = A_direct - 2·Re(c₀)/δ  (Re(c₀) Richardson 보정)")
log()
flush_file()

# ━━━━━━━━━━━ PARI 초기화 ━━━━━━━━━━━
import cypari2
gp = cypari2.Pari()
gp.allocatemem(3000 * 1024 * 1024)   # T=500: 3GB (94c는 2GB)
gp("default(realprecision, 100)")
log("[Init] PARI: 3GB, prec=100 ✅")
log()
flush_file()

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
CENTER       = 1.5
DELTA_DIRECT = 0.01
N_TEST       = 10
T_ZEROS_MAX  = 500.0   # #94c: 150.0 → #100: 500.0

log(f"파라미터: center={CENTER}, δ={DELTA_DIRECT}, N_test={N_TEST}, T_max={T_ZEROS_MAX}")
log(f"  (#94c 비교: T=150 → T=500, N예상: 271 → ~1000)")
log()
flush_file()

# ━━━━━━━━━━━ EM 꼬리 보정 (GL(3) 수정 공식) ━━━━━━━━━━━
try:
    from scipy import integrate as _sci_int

    def tail_gl3(t0_val, t_last, mode):
        """N'(T) = d/(2π)*log(T/(2π)) = 3/(2π)*log(T/(2π))"""
        c = 3.0 / (2.0 * math.pi)   # d=3
        c2pi = 2.0 * math.pi
        if mode == 'B':
            fn = lambda t: 2.0*t0_val/(t**2 - t0_val**2) * c * math.log(t/c2pi) if t > t0_val+1e-6 else 0.0
        else:
            fn = lambda t: 2.0/t**2 * c * math.log(t/c2pi)
        v, _ = _sci_int.quad(fn, t_last, 1e7, limit=500, epsabs=1e-12)
        return v

    HAS_TAIL = True
    log("[EM tail] scipy ✅, 밀도: d/(2π)*log(T/(2π)) = 3/(2π)*log(T/(2π))")
except ImportError:
    HAS_TAIL = False
    def tail_gl3(t0_val, t_last, mode): return 0.0
    log("[EM tail] scipy 없음 — tail 보정 없이 진행")
log()
flush_file()

# ━━━━━━━━━━━ Step 1: 영점 확보 ━━━━━━━━━━━
log("=" * 72)
log(f"[Step 1] sym²(11a1) 영점 확보 (T≤{T_ZEROS_MAX})")
log(f"  예상 소요: lfuninit ~3-5분, lfunzeros 추가 ~수분")
log("=" * 72)
t0_s = time.time()

gp("E11 = ellinit([0,-1,1,-10,-20])")
gp("L11 = lfunsympow(E11, 2)")
log(f"  lfuninit([0, {T_ZEROS_MAX}]) 시작...")
flush_file()

gp(f"L11i = lfuninit(L11, [0, {T_ZEROS_MAX}])")
log(f"  lfuninit 완료 ({time.time()-t0_s:.1f}s). 영점 계산 시작...")
flush_file()

gp(f"zvec = lfunzeros(L11i, {T_ZEROS_MAX})")
n_zeros_raw = int(gp("length(zvec)"))
log(f"  lfunzeros 완료. raw 영점 수: {n_zeros_raw}")
flush_file()

all_zeros_list = []
for i in range(1, n_zeros_raw + 1):
    try:
        all_zeros_list.append(float(gp(f"zvec[{i}]")))
    except Exception as e:
        log(f"    WARNING 영점 #{i} 읽기 실패: {e}")

all_zeros = np.array(sorted(all_zeros_list))
N_ALL = len(all_zeros)
t_last = float(all_zeros[-1]) if N_ALL > 0 else 0.0

log(f"  N={N_ALL}개, t ∈ [{all_zeros[0]:.4f}, {t_last:.4f}]")
log(f"  소요: {time.time()-t0_s:.1f}s")

if N_ALL < 10:
    log("⚠️ 영점 10개 미만 — 중단"); flush_file(); sys.exit(1)
if N_ALL < 200:
    log(f"⚠️ N={N_ALL} < 200 — 예상보다 적음. 계속 진행.")
log(f"  → #94c(N=271) 대비 {N_ALL/271:.1f}× 증가")
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
            # c₀_raw = Λ'/Λ(ρ₀+ε) - 1/ε
            c0_raw_val = conn_re - 1.0/eps_f + 1j*(conn_im)
            c0_raw.append(c0_raw_val)
        except Exception as e:
            log(f"    WARNING Richardson eps={eps_f}: {e}")
            c0_raw.append(float('nan'))

    # Richardson 외삽
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
log(f"  전체 N={N_ALL}개 영점 사용 (T≤{T_ZEROS_MAX})")
log("=" * 72)

# 동일 10개 테스트 영점 (#94c와 같은 t₁~t₁₀)
test_zeros = all_zeros[:N_TEST]
log(f"  테스트 영점: t = {test_zeros[0]:.4f} ~ {test_zeros[-1]:.4f}")
log(f"  N_Had = N_ALL = {N_ALL}개  (#94c 대비 {N_ALL/271:.1f}×)")
log()
flush_file()

# #94c 기준값 (비교용)
ref_err_94c = {
    1: (3.8993, 1.82),
    2: (4.7346, 1.77),
    3: (6.1895, 3.05),
    4: (7.3120, 3.34),
    5: (8.6501, 4.30),
    6: (10.1282, 3.89),
    7: (10.9368, 3.80),
    8: (12.0706, 3.79),
    9: (12.7446, 3.05),
   10: (13.3353, 1.65),
}
log("  [#94c 기준 err_new% 참고]:")
for k, (t_ref, e_ref) in ref_err_94c.items():
    log(f"    #{k}: t₀={t_ref:.4f}, err={e_ref:.2f}%")
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
        log(f"    ⚠️ skip (A_dir NaN)")
        results.append({'t': t0_val, 'skip': True})
        log(); flush_file()
        continue

    # A_true = A_direct - 2*Re(c₀)/δ
    if np.isfinite(Re_c0):
        A_true = A_dir - 2.0 * Re_c0 / DELTA_DIRECT
        correction = 2.0 * Re_c0 / DELTA_DIRECT
    else:
        A_true = A_dir
        correction = 0.0
    log(f"    A_true = {A_true:.6f}  (보정량: {correction:.4f} = 2·Re(c₀)/δ)")

    # Hadamard sum (all N_ALL zeros)
    B_had, H1_had, A_had = hadamard_B_H1(t0_val, all_zeros)
    log(f"    B_had = {B_had:.6f}, H₁_had = {H1_had:.6f}")
    log(f"    A_had (N={N_ALL}) = {A_had:.6f}")

    # EM 꼬리 보정
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

    # 비교
    err_old = abs(A_corr - A_dir) / (abs(A_dir) + 1e-12)
    err_new = abs(A_corr - A_true) / (abs(A_true) + 1e-12)

    # #94c 대비 개선량
    ref_err = ref_err_94c.get(idx+1, (t0_val, float('nan')))[1]
    delta_err = err_new*100 - ref_err if not math.isnan(ref_err) else float('nan')
    delta_str = f"(vs #94c: {ref_err:.2f}% → Δ={delta_err:+.2f}%)" if not math.isnan(ref_err) else ""

    log(f"    err vs A_direct: {err_old*100:.2f}%  (이전 기준)")
    log(f"    err vs A_true:   {err_new*100:.2f}%  (Re(c₀) 보정 후) "
        f"{'✓<1%' if err_new<0.01 else ('✓<2%' if err_new<0.02 else '✗')}  {delta_str}")
    log(f"    소요: {time.time()-tz:.1f}s")

    results.append({
        't': t0_val, 'A_direct': A_dir, 'A_true': A_true,
        'Re_c0': Re_c0, 'correction': correction,
        'B_had': B_had, 'H1_had': H1_had, 'A_had': A_had,
        'Bt': Bt, 'H1t': H1t, 'A_corr': A_corr,
        'err_old': err_old, 'err_new': err_new,
        'ref_err_94c': ref_err,
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
log("최종 판정 — Hadamard GL(3) 고정밀 재검증 (#100)")
log("=" * 72)
log()

nv      = len(valid)
n_old_1 = sum(1 for r in valid if r['err_old'] < 0.01)
n_old_2 = sum(1 for r in valid if r['err_old'] < 0.02)
n_new_1 = sum(1 for r in valid if r['err_new'] < 0.01)
n_new_2 = sum(1 for r in valid if r['err_new'] < 0.02)
n_rec0  = sum(1 for r in valid
              if np.isfinite(r.get('Re_c0', float('nan'))) and abs(r['Re_c0']) < 0.01)

log(f"{'#':>3} {'t₀':>10} {'A_direct':>10} {'A_true':>10} {'A_corr':>10} {'err_old%':>9} {'err_new%':>9} {'#94c%':>7}  pass")
log("-" * 90)
for i, r in enumerate(results):
    if r['skip']:
        log(f"  {i+1:>2} {r['t']:>10.4f}  SKIP")
        continue
    ok = '✓<1%' if r['err_new'] < 0.01 else ('✓<2%' if r['err_new'] < 0.02 else '✗')
    ref_str = f"{r['ref_err_94c']:.2f}" if not math.isnan(r.get('ref_err_94c', float('nan'))) else "  -"
    log(f"  {i+1:>2} {r['t']:>10.4f} {r['A_direct']:>10.4f} {r['A_true']:>10.4f} "
        f"{r['A_corr']:>10.4f} {r['err_old']*100:>9.2f} {r['err_new']*100:>9.2f} {ref_str:>7}  {ok}")
log()

log(f"판정 요약 (N={N_ALL}개 영점, T≤{T_ZEROS_MAX}):")
log(f"  [A_direct 기준, EM] <1%: {n_old_1}/{nv}  <2%: {n_old_2}/{nv}")
log(f"  [A_true  기준, EM] <1%: {n_new_1}/{nv}  <2%: {n_new_2}/{nv}")
log(f"  Re(c₀)<0.01: {n_rec0}/{nv}  (FE 대칭성)")
log()

# #94c 비교
if valid:
    n_94c_2 = sum(1 for r in valid
                  if not math.isnan(r.get('ref_err_94c', float('nan'))) and r['ref_err_94c'] < 2.0)
    n_improved = sum(1 for r in valid
                     if not math.isnan(r.get('ref_err_94c', float('nan')))
                     and r['err_new']*100 < r['ref_err_94c'])
    log(f"  #94c(N=271) 비교: <2%: {n_94c_2}/10 → #100(N={N_ALL}): {n_new_2}/{nv}")
    log(f"  개선된 영점: {n_improved}/{nv}개")
log()

ok_strict  = (nv >= 10) and (n_new_1 >= 7)
ok_crit    = (nv >= 10) and (n_new_2 >= 7)
ok_medium  = (nv >= 10) and (n_new_1 >= 5)
ok_rec0    = (n_rec0 >= nv - 2)

log("성공 기준 (A_true 기준):")
log(f"  영점 10개: {'✅' if nv >= 10 else '❌'} ({nv}개)")
log(f"  <1% 7개 이상: {n_new_1}/10  {'✅' if n_new_1>=7 else '❌'}")
log(f"  <2% 7개 이상: {n_new_2}/10  {'✅' if n_new_2>=7 else '❌'}")
log(f"  <1% 5개 이상: {n_new_1}/10  {'✅' if n_new_1>=5 else '❌'}")
log(f"  Re(c₀)<0.01: {n_rec0}/{nv}  {'✅' if ok_rec0 else '⚠️'}")
log()

if ok_strict and ok_rec0:
    verdict = "★★★ 강양성 — A=B²+2H₁ GL(3) 엄격 1% 기준 통과 (T=500 고정밀)"
    verdict2 = "  → degree-보편성: d=1(★★★), d=2(★★), d=3(★★★)"
elif ok_crit and ok_rec0:
    verdict = "★★ 조건부 양성 — A=B²+2H₁ GL(3) 2% 기준 통과 (T=500 고정밀)"
    verdict2 = "  → degree-보편성: d=1(★★★), d=2(★★), d=3(★★)"
elif ok_medium and ok_rec0:
    verdict = "★ 약양성 — 1% 기준 5/10 통과, 수렴 부분적"
    verdict2 = "  → d=3 Hadamard 정밀도 한계 존재"
elif n_new_2 > 3:
    verdict = f"⚠️ 부분 — <2%: {n_new_2}/10 (#94c 대비 {'개선' if n_new_2 > 3 else '유사'})"
    verdict2 = "  → d=3 구조적 한계로 기록"
else:
    verdict = f"❌ 미결 — <2%: {n_new_2}/10 (#94c: 3/10 대비 {'개선' if n_new_2 > 3 else '미개선'})"
    verdict2 = "  → EM tail이 아닌 구조적 원인 분석 필요"

log(f"  {verdict}")
log(verdict2)
log()

log("[Re(c₀) 보정 방법론]")
log("  A = κ(ρ₀+δ) - 1/δ² = 2·Re(c₀)/δ + B² + 2H₁ + O(δ)")
log("  Re(c₀) ≠ 0 (수치적): 2·Re(c₀)/δ ≈ 0.2-0.6 (A의 2-5% 오염)")
log("  A_true = A_direct - 2·Re(c₀)_rich/δ ≈ B² + 2H₁")
log("  이는 A의 두 독립 추정: (1) Richardson, (2) Hadamard+EM")
log()
log("[EM tail 기여 분석]")
log(f"  #94c: t_last={149.75:.2f}, tail 기여 3-5%")
log(f"  #100: t_last={t_last:.2f}, tail 기여 ~{3*149.75/max(t_last,1):.1f}% 예상 (∝1/N)")
log()

log("도표 (degree 비교, T=500 고정밀 포함):")
log(f"  d=1 GL(1) ζ       : 13/13 <0.003%  ★★★  (N=2000, #90)")
log(f"  d=2 GL(2) L(s,Δ)  :  7/10 <1.3%    ★★   (N=86, #92)")
log(f"  d=3 GL(3) sym²    :")
log(f"    #94c (N=271)    :  3/10 <2%       ★    (T=150)")
log(f"    #100 (N={N_ALL:3d})    :  {n_new_1}/10 <1%, {n_new_2}/10 <2%  (T={T_ZEROS_MAX:.0f})")
log()
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"결과: {OUTFILE}")
log("=" * 72)
flush_file()
