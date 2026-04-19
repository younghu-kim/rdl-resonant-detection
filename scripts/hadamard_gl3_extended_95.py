#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #95 — Hadamard GL(3) sym²(11a1) 확장 재검증 (N≥200)
=============================================================================
동기:
  #94에서 N=64 영점으로 raw Hadamard 오차 16-33%, EM 과보정 20-55%.
  → T_MAX=150으로 200+개 영점 확보, EM 없이 raw 수렴 확인.

방법:
  [1] PARI lfunsympow(E,2) → lfunzeros(Ls2, 150) → N≥200 영점
  [2] 테스트 영점 12개 (t ∈ [10, 40] — 양쪽 가장자리 제외)
  [3] Hadamard paired sum at N=20,50,100,150,200,ALL → 수렴 추세
  [4] raw Hadamard만으로 err<2% 달성 가능한지 확인
  [5] EM은 참조용으로만 계산 (주 판정 기준 아님)

성공 기준:
  - raw Hadamard (EM 없이): err<5% 비율 ≥ 7/12
  - N 증가 시 단조 수렴 확인
  - Re(c₀)≈0: 전부 <0.01
  - GL(1) 0.003% → GL(2) 0.8% → GL(3) ?% 의 degree 스케일링 패턴

결과: results/hadamard_gl3_extended_95.txt
=============================================================================
"""

import sys, os, time, math
import numpy as np

OUTFILE = os.path.expanduser('~/Desktop/gdl_unified/results/hadamard_gl3_extended_95.txt')
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))

def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

log("=" * 80)
log("결과 #95 — Hadamard GL(3) sym²(11a1) 확장 재검증 (N≥200, EM-free)")
log("=" * 80)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"동기: #94에서 N=64 부족 (raw 16-33%, EM 과보정). T_MAX=150으로 확장.")
log(f"이전: GL(1) 13/13 0.003%, GL(2) 7/10 <1.3%, GL(3) #94 0/8 미결")
log()
flush_file()

# ─── PARI 초기화 ──────────────────────────────────────────────────────────
import cypari2
gp = cypari2.Pari()
gp.allocatemem(4000 * 1024 * 1024)  # 4GB (T_MAX=150 대비)
gp("default(realprecision, 100)")
log("PARI: 4GB, realprecision=100")
log()
flush_file()

# ─── 파라미터 ─────────────────────────────────────────────────────────────
DELTA = 0.01
CENTER = 1.5  # k=3, sym²(E) — 절대 2.0 금지
T_MAX = 150.0
N_TEST = 12

# ─── sym²(11a1) 초기화 ────────────────────────────────────────────────────
log("=" * 80)
log("[1] sym²(11a1) 초기화 + 영점 수집 (T_MAX=150)")
log("=" * 80)
t0_total = time.time()

gp("E = ellinit([0,-1,1,-10,-20])")
gp("Ls2 = lfunsympow(E, 2)")

k_val = int(round(float(str(gp("Ls2[4]")))))
center = k_val / 2.0
assert abs(center - CENTER) < 1e-9, f"center={center} != 1.5"
log(f"  k={k_val}, center={center} ✓")

# lfuninit — T_MAX=150
log(f"  lfuninit([0, {T_MAX}]) ...")
t_init = time.time()
gp(f"Linit = lfuninit(Ls2, [0, {T_MAX}])")
log(f"  lfuninit 완료 ({time.time()-t_init:.1f}s)")

# 영점 수집
log(f"  lfunzeros(Ls2, {T_MAX}) ...")
t_zz = time.time()
gp(f"zvec = lfunzeros(Ls2, {T_MAX})")
n_zeros = int(gp("length(zvec)"))
log(f"  {n_zeros}개 영점 수집 ({time.time()-t_zz:.1f}s)")

all_zeros = np.array([float(gp(f"zvec[{i}]")) for i in range(1, n_zeros + 1)])
all_zeros = np.sort(all_zeros)
log(f"  범위: t ∈ [{all_zeros[0]:.4f}, {all_zeros[-1]:.4f}]")
log(f"  밀도: {n_zeros/(all_zeros[-1]-all_zeros[0]):.3f} zeros/unit")
log()
flush_file()

# ─── 테스트 영점 선택 ────────────────────────────────────────────────────
# t ∈ [10, 100] 범위에서 균등 샘플링 (가장자리 제외)
mask = (all_zeros >= 10.0) & (all_zeros <= 100.0)
cands = all_zeros[mask]
if len(cands) >= N_TEST:
    step = len(cands) // N_TEST
    test_zeros = cands[::step][:N_TEST]
else:
    test_zeros = cands[:N_TEST]

log(f"[2] 테스트 영점 {len(test_zeros)}개 (t ∈ [10, 100] 균등)")
for i, tz in enumerate(test_zeros):
    log(f"  #{i+1}: t₀ = {tz:.6f}")
log()
flush_file()

# ─── Hadamard N 체크포인트 ───────────────────────────────────────────────
N_CHECKS = sorted(set([20, 50, 100, 150, min(200, n_zeros), n_zeros]))
N_CHECKS = [n for n in N_CHECKS if n <= n_zeros]
log(f"[3] Hadamard N 체크포인트: {N_CHECKS}")
log()
flush_file()

# ─── 핵심 함수들 ─────────────────────────────────────────────────────────
def compute_A_direct(t0_val, delta=DELTA):
    """A_direct = κ - 1/δ² via PARI lfunlambda"""
    sigma = center + delta
    gp(f"scur = {sigma:.12f} + I*{t0_val:.10f}")
    gp("L0v = lfunlambda(Linit, scur)")
    gp("L1v = lfunlambda(Linit, scur, 1)")
    abs_L0 = float(gp("abs(L0v)"))
    if abs_L0 < 1e-150:
        return float('nan'), float('nan')
    kappa = (float(gp("abs(L1v)")) / abs_L0) ** 2
    return kappa - 1.0 / delta**2, kappa

def compute_Re_c0(t0_val):
    """Richardson 외삽으로 Re(c₀) 추정"""
    eps_list = [5e-3, 2.5e-3, 1.25e-3]
    c0_raw = []
    for eps in eps_list:
        sigma = center + eps
        gp(f"sr = {sigma:.12f} + I*{t0_val:.10f}")
        gp("L0r = lfunlambda(Linit, sr)")
        gp("L1r = lfunlambda(Linit, sr, 1)")
        abs_L0 = float(gp("abs(L0r)"))
        if abs_L0 < 1e-150:
            return float('nan')
        re_conn = float(gp("real(L1r / L0r)"))
        c0_raw.append(re_conn - 1.0 / eps)
    # Richardson 2단계
    r1 = [(4*c0_raw[i+1] - c0_raw[i]) / 3.0 for i in range(2)]
    return (4*r1[1] - r1[0]) / 3.0

def hadamard_B_H1(t0, tzeros):
    """Hadamard paired sum: B, H₁, A=B²+2H₁"""
    mask = np.abs(tzeros - t0) > 1e-6
    tv = tzeros[mask]
    d = t0 - tv
    s = t0 + tv
    ok = (np.abs(d) > 1e-12) & (np.abs(s) > 1e-12)
    d, s = d[ok], s[ok]
    B = -1.0 / (2.0 * t0) + (-np.sum(1.0/d + 1.0/s))
    H1 = 1.0 / (4.0 * t0**2) + np.sum(1.0/d**2 + 1.0/s**2)
    return float(B), float(H1), float(B**2 + 2.0*H1)

# ─── [4] 메인 계산 루프 ──────────────────────────────────────────────────
log("=" * 80)
log(f"[4] 테스트 영점별 A(t₀) 계산 (δ={DELTA})")
log("=" * 80)
log()
flush_file()

results = []
for idx, t0_val in enumerate(test_zeros):
    log(f"── #{idx+1}: t₀ = {t0_val:.6f} ──")
    tz = time.time()

    # A_direct
    A_dir, kappa = compute_A_direct(t0_val)
    if not np.isfinite(A_dir):
        log(f"  A_direct 실패 — SKIP")
        results.append({'t': t0_val, 'skip': True})
        log()
        flush_file()
        continue

    log(f"  A_direct = {A_dir:.4f} (κ={kappa:.4f})")

    # Re(c₀)
    re_c0 = compute_Re_c0(t0_val)
    log(f"  Re(c₀) = {re_c0:.4e}")

    # Hadamard at multiple N
    had_data = {}
    for N in N_CHECKS:
        B_n, H1_n, A_n = hadamard_B_H1(t0_val, all_zeros[:N])
        err_n = abs(A_n - A_dir) / abs(A_dir) if abs(A_dir) > 1e-10 else float('nan')
        had_data[N] = {'B': B_n, 'H1': H1_n, 'A': A_n, 'err': err_n}

    # 수렴 추세 출력
    for N in N_CHECKS:
        d = had_data[N]
        log(f"  N={N:>4}: A_had={d['A']:>10.4f}  err={d['err']*100:>7.2f}%")

    # 최대 N에서의 결과
    Nmax = N_CHECKS[-1]
    best = had_data[Nmax]

    log(f"  ({time.time()-tz:.1f}s)")
    log()

    results.append({
        't': t0_val,
        'A_direct': A_dir,
        'kappa': kappa,
        're_c0': re_c0,
        'had_data': had_data,
        'A_best': best['A'],
        'err_best': best['err'],
        'skip': False,
    })
    flush_file()

# ─── [5] 결과 요약 ───────────────────────────────────────────────────────
log("=" * 80)
log("결과 요약")
log("=" * 80)
log()

valid = [r for r in results if not r['skip']]
nv = len(valid)

# 다양한 기준으로 통과율
thresholds = [0.02, 0.05, 0.10]
for th in thresholds:
    n_pass = sum(1 for r in valid if r['err_best'] < th)
    log(f"  err < {th*100:.0f}%: {n_pass}/{nv}")
log()

# 요약 테이블
Nmax = N_CHECKS[-1]
log(f"{'#':>3} {'t₀':>10} {'A_direct':>10} {'A_had':>10} {'err%':>8} {'Re(c₀)':>12}")
log("-" * 60)
for i, r in enumerate(results):
    if r['skip']:
        log(f"  {i+1:>2} {r['t']:>10.4f}  SKIP")
        continue
    log(f"  {i+1:>2} {r['t']:>10.4f} {r['A_direct']:>10.4f} {r['A_best']:>10.4f} "
        f"{r['err_best']*100:>8.2f} {r['re_c0']:>12.4e}")
log()

# Re(c₀) 통계
re_c0_vals = [abs(r['re_c0']) for r in valid if np.isfinite(r['re_c0'])]
if re_c0_vals:
    log(f"Re(c₀) 통계: max={max(re_c0_vals):.4e}, mean={np.mean(re_c0_vals):.4e}")
    log(f"Re(c₀)≈0: {sum(1 for v in re_c0_vals if v < 0.01)}/{len(re_c0_vals)} (<0.01)")
log()

# ─── [6] 수렴 추세 분석 ──────────────────────────────────────────────────
log("=" * 80)
log("N 수렴 추세 (EM-free, raw Hadamard)")
log("=" * 80)
log()

# N별 평균 오차
for N in N_CHECKS:
    errs = [r['had_data'][N]['err'] for r in valid if N in r['had_data']]
    if errs:
        log(f"  N={N:>4}: mean_err={np.mean(errs)*100:.2f}%, "
            f"median={np.median(errs)*100:.2f}%, "
            f"min={min(errs)*100:.2f}%, max={max(errs)*100:.2f}%")
log()

# 수렴 속도 추정: err ~ C/N^α
if len(N_CHECKS) >= 3 and len(valid) >= 3:
    log("수렴 속도 추정 (err ~ C·N^{-α}):")
    for r in valid[:3]:
        Ns = np.array([N for N in N_CHECKS if N >= 20])
        errs = np.array([r['had_data'][N]['err'] for N in Ns])
        mask_ok = errs > 0
        if mask_ok.sum() >= 2:
            # log-log 선형 회귀
            logN = np.log(Ns[mask_ok])
            logE = np.log(errs[mask_ok])
            if len(logN) >= 2:
                slope, intercept = np.polyfit(logN, logE, 1)
                alpha = -slope
                C = np.exp(intercept)
                # N needed for 2%
                if alpha > 0:
                    N_2pct = (C / 0.02) ** (1.0 / alpha)
                    N_5pct = (C / 0.05) ** (1.0 / alpha)
                else:
                    N_2pct = N_5pct = float('inf')
                log(f"  t₀={r['t']:.2f}: α={alpha:.2f}, C={C:.3f} → "
                    f"N(2%)≈{N_2pct:.0f}, N(5%)≈{N_5pct:.0f}")
    log()
flush_file()

# ─── [7] 최종 판정 ───────────────────────────────────────────────────────
log("=" * 80)
log("최종 판정 — Hadamard A=B²+2H₁ GL(3) 확장 검증")
log("=" * 80)
log()

n_pass_2 = sum(1 for r in valid if r['err_best'] < 0.02)
n_pass_5 = sum(1 for r in valid if r['err_best'] < 0.05)
n_pass_10 = sum(1 for r in valid if r['err_best'] < 0.10)

log(f"  영점 수: {n_zeros}개 (T_MAX={T_MAX})")
log(f"  테스트: {nv}개 유효")
log(f"  err<2%: {n_pass_2}/{nv},  err<5%: {n_pass_5}/{nv},  err<10%: {n_pass_10}/{nv}")
log()
log(f"  degree-보편 정리 수치 현황:")
log(f"    GL(1) ζ(s):     13/13, mean 0.003%  ★★★")
log(f"    GL(2) L(s,Δ):    7/10 <1.3%         ★★")
log(f"    GL(3) sym²(11a1): {n_pass_5}/{nv} <5% (N={n_zeros})")
log()

if n_pass_2 >= 7:
    verdict = "★★★ 강양성 — GL(3) err<2% 확립. degree-보편 정리 완전 지지."
elif n_pass_5 >= 7:
    verdict = "★★ 양성 — GL(3) err<5%. degree 증가에 따른 수렴 지연은 정상 범위."
elif n_pass_10 >= 7:
    verdict = "★ 조건부 양성 — 수렴 방향 확인. err<10% 다수. 영점 추가 확보 시 개선 기대."
elif n_pass_10 >= 3:
    verdict = "⚠️ 미결 — 수렴 추세 존재하나 N 부족. 추가 영점(N≥500) 필요."
else:
    verdict = "❌ 음성 — 수렴 불확실. 이론 재검토 또는 계산 방법 변경 필요."

log(f"  판정: {verdict}")
log()
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"총 소요: {time.time()-t0_total:.1f}s")
log(f"결과: {OUTFILE}")
log("=" * 80)
flush_file()
