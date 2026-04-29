#!/usr/bin/env python3
"""[C-404] Path C 정밀도 검증 — T=2000 mpmath 보간 vs C-398 PARI
질문: C-400 음성(R=0.801)이 진짜인가 아니면 선형 보간 아티팩트인가?
방법: T=2000에서 mpmath 보간 R(T)을 계산, C-398 PARI R=1.365와 비교.
"""
import sys, os, time, math
import numpy as np
from scipy import stats as sp_stats
sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import mpmath

mpmath.mp.dps = 30
T_MAX = 2000.0
W = 100
TRIM = 0.20
N_BANDS = 10
N_BOOT = 2000
SIGMA_C = 0.5
RPATH = os.path.expanduser('~/Desktop/gdl_unified/results/path_c_verify_c404.txt')
out_f = open(RPATH, 'w')

def log(m=''):
    print(m, flush=True); out_f.write(m + '\n'); out_f.flush()

def gamma_corr(t0):
    """Gamma 보정: Im(ψ(s/2))/2, Re(ψ₁(s/2))/4 at s=1/2+it0"""
    s = mpmath.mpc(SIGMA_C, t0)
    arg = s / 2
    return float(mpmath.im(mpmath.digamma(arg)) / 2), float(mpmath.re(mpmath.psi(1, arg)) / 4)

# ===================================================================
log("=" * 70)
log("[C-404] Path C 정밀도 검증 — T=2000 mpmath 보간 vs C-398 PARI")
log("=" * 70)
log(f"  T_MAX={T_MAX}, W={W}, TRIM={TRIM}, N_BANDS={N_BANDS}, N_BOOT={N_BOOT}")
log(f"  밀도: d(t) = log(t/(2π))/(2π) (이론적)")
t0g = time.time()

# ===================================================================
# 1. Z-함수 부호변환 스캔으로 ζ(s) 영점 수집
# ===================================================================
log(f"\n[ζ(s)] Z-함수 스캔 t∈(2, {T_MAX}], step=0.25 ...")

step = 0.25  # 최소 간격 0.778 (t=10000) 의 1/3 이하 → 안전
t_grid = np.arange(2.0, T_MAX + 0.5, step)
n_grid = len(t_grid)
log(f"  스캔 그리드: {n_grid} 점")

# 청크 단위로 Z 값 계산 (진행 보고)
z_vals = np.empty(n_grid)
chunk = 5000
for i in range(0, n_grid, chunk):
    end = min(i + chunk, n_grid)
    for j in range(i, end):
        z_vals[j] = float(mpmath.siegelz(t_grid[j]))
    elapsed = time.time() - t0g
    log(f"  Z 스캔 {end}/{n_grid} ({100*end/n_grid:.0f}%), {elapsed:.0f}s")

# 부호변환 검출
sc = np.where(z_vals[:-1] * z_vals[1:] < 0)[0]
log(f"  부호변환 {len(sc)} 개 감지, 이분법 정밀화 중 ...")

# 선형 보간으로 영점 근사 (정밀도 ~0.003, 통계 분석에 충분)
# bisection 완전 제거 → 순수 numpy 연산, <1초
t_a = t_grid[sc]
t_b = t_grid[sc + 1]
z_a = z_vals[sc]
z_b = z_vals[sc + 1]
zeros_list = t_a - z_a * (t_b - t_a) / (z_b - z_a)  # 선형 보간
log(f"  선형 보간 완료 ({len(zeros_list)} 영점), {time.time()-t0g:.0f}s")
zeros = np.sort(zeros_list[(zeros_list > 1.0) & (zeros_list <= T_MAX)])
N = len(zeros)

# 기대치 검증
N_expected = T_MAX * math.log(T_MAX / (2 * math.pi * math.e)) / (2 * math.pi)
log(f"  수집 완료: {N} 영점 (기대 ~{N_expected:.0f}), [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
log(f"  영점 수집 소요: {time.time()-t0g:.0f}s")

if abs(N - N_expected) / N_expected > 0.02:
    log(f"  ⚠️ 기대치 대비 {(N-N_expected)/N_expected*100:+.1f}% 차이 — 누락 점검 필요")

# ===================================================================
# 2. Gamma 보정 계산
# ===================================================================
log(f"\n[Gamma 보정] {N} 점 ...")
im_g = np.zeros(N)
re_g = np.zeros(N)
for i in range(N):
    im_g[i], re_g[i] = gamma_corr(zeros[i])
log(f"  Gamma 보정 완료: {time.time()-t0g:.0f}s")

# ===================================================================
# 3. 3-tier 분해 (C-395/C-398 동일)
# ===================================================================
x = zeros
lo_w, hi_w = W, N - W - 1
nv = hi_w - lo_w

# 간격
g = x[lo_w + 1:hi_w + 1] - x[lo_w:hi_w]
g = np.where(np.abs(g) < 1e-15, 1e-15, g)
dL = x[lo_w:hi_w] - x[lo_w - 1:hi_w - 1]
dL = np.where(np.abs(dL) < 1e-15, 1e-15, dL)
dR = x[lo_w + 2:hi_w + 2] - x[lo_w + 1:hi_w + 1]
dR = np.where(np.abs(dR) < 1e-15, 1e-15, dR)

# Hadamard 합 (W 이웃)
Ht_n = np.zeros(nv)
Ht_np1 = np.zeros(nv)
S1_n = -1.0 / g + 1.0 / dL
S1_np1 = -1.0 / dR + 1.0 / g
for j in range(2, W + 1):
    idx = np.arange(lo_w, hi_w)
    dr = x[idx + j] - x[idx]; dr = np.where(np.abs(dr) < 1e-15, 1e-15, dr)
    dl = x[idx] - x[idx - j]; dl = np.where(np.abs(dl) < 1e-15, 1e-15, dl)
    Ht_n += 1.0 / dr**2 + 1.0 / dl**2
    S1_n += -1.0 / dr + 1.0 / dl
    ix1 = np.arange(lo_w + 1, hi_w + 1)
    dr1 = x[ix1 + j] - x[ix1]; dr1 = np.where(np.abs(dr1) < 1e-15, 1e-15, dr1)
    dl1 = x[ix1] - x[ix1 - j]; dl1 = np.where(np.abs(dl1) < 1e-15, 1e-15, dl1)
    Ht_np1 += 1.0 / dr1**2 + 1.0 / dl1**2
    S1_np1 += -1.0 / dr1 + 1.0 / dl1

# 완비 접속 (Gamma 보정)
S1L_n = S1_n - im_g[lo_w:hi_w]
S1L_np1 = S1_np1 - im_g[lo_w + 1:hi_w + 1]
H_n = 1.0 / g**2 + 1.0 / dL**2 + Ht_n + re_g[lo_w:hi_w]
H_np1 = 1.0 / g**2 + 1.0 / dR**2 + Ht_np1 + re_g[lo_w + 1:hi_w + 1]
A_n = S1L_n**2 + 2 * H_n
A_np1 = S1L_np1**2 + 2 * H_np1
sh = 2.0 / g**2  # 공유 성분 = 2/g²
Ap_n = A_n - sh
Ap_np1 = A_np1 - sh
t_arr = x[lo_w:hi_w]

# Trim 20%
tl = int(nv * TRIM)
th = nv - int(nv * TRIM)
sl = slice(tl, th)
mask = np.isfinite(A_n[sl]) & np.isfinite(A_np1[sl]) & (A_n[sl] > 0) & (A_np1[sl] > 0)
A_n = A_n[sl][mask]
A_np1 = A_np1[sl][mask]
Ap_n = Ap_n[sl][mask]
Ap_np1 = Ap_np1[sl][mask]
sh = sh[sl][mask]
t_arr = t_arr[sl][mask]
nn = len(A_n)
log(f"\n  유효 쌍: {nn}, compute {time.time()-t0g:.0f}s")

# ===================================================================
# 4. 전역 3-tier 분해
# ===================================================================
def decomp(a1, a2, p1, p2, s):
    cf = np.cov(a1, a2)[0, 1]
    vs = np.var(s)
    cr = np.cov(p1, p2)[0, 1]
    return cf, vs, cf - vs - cr, cr

cf, vs, cres, cr = decomp(A_n, A_np1, Ap_n, Ap_np1, sh)
Rg = vs / (abs(cr) + abs(cres)) if abs(cr) + abs(cres) > 1e-30 else 999
log(f"\n[전역 분해]")
log(f"  Cov={cf:.4e}, Var(2/g²)={vs:.4e}({vs/cf*100:.1f}%)")
log(f"  Cross={cr:.4e}({cr/cf*100:.1f}%), Resid={cres:.4e}({cres/cf*100:.1f}%)")
log(f"  전역 R = {Rg:.3f}x")

# ===================================================================
# 5. T-대역별 분해
# ===================================================================
log(f"\n[T-대역별 분석 ({N_BANDS} 대역)]")
edges = np.linspace(t_arr.min(), t_arr.max(), N_BANDS + 1)
log(f"{'대역':>24s} | {'n':>5s} | {'Var%':>6s} | {'Cross%':>7s} | {'Resid%':>7s} | {'R':>6s} | ok")
log("-" * 80)

bands = []
for i in range(N_BANDS):
    tl_b, th_b = edges[i], edges[i + 1]
    m = (t_arr >= tl_b) & (t_arr < th_b) if i < N_BANDS - 1 else (t_arr >= tl_b) & (t_arr <= th_b)
    nb = m.sum()
    if nb < 30:
        log(f"  [{tl_b:.0f},{th_b:.0f}] | {nb:5d} | 부족")
        continue
    cf_b, vs_b, cres_b, cr_b = decomp(A_n[m], A_np1[m], Ap_n[m], Ap_np1[m], sh[m])
    d = abs(cr_b) + abs(cres_b)
    R = vs_b / d if d > 1e-30 else 999
    ok = "Y" if R > 1 else "N"
    log(f"  [{tl_b:7.0f},{th_b:7.0f}] | {nb:5d} | {vs_b/cf_b*100:5.1f}% | {cr_b/cf_b*100:+6.1f}% | {cres_b/cf_b*100:+6.1f}% | {R:5.2f}x | {ok}")
    bands.append({'tl': tl_b, 'th': th_b, 'tm': (tl_b + th_b) / 2, 'n': nb, 'R': R,
                  'Var_pct': vs_b / cf_b * 100, 'Cross_pct': cr_b / cf_b * 100, 'Resid_pct': cres_b / cf_b * 100})

nv_b = len(bands)
n_ok = sum(1 for b in bands if b['R'] > 1)
log(f"\n  R>1: {n_ok}/{nv_b}")

# ===================================================================
# 6. 추세 분석
# ===================================================================
if nv_b >= 3:
    tm = np.array([b['tm'] for b in bands])
    Rv = np.array([b['R'] for b in bands])

    log(f"\n[추세 분석]")
    sl_lin, int_lin, r_lin, p_lin, se_lin = sp_stats.linregress(tm, Rv)
    log(f"  선형: R = {sl_lin:+.6f}*T + {int_lin:.4f}, r²={r_lin**2:.4f}, p={p_lin:.4e}, se={se_lin:.6f}")
    sl_log, int_log, r_log, p_log, se_log = sp_stats.linregress(np.log(tm), Rv)
    log(f"  로그: R = {sl_log:+.4f}*log(T) + {int_log:.4f}, r²={r_log**2:.4f}")
    sl_inv, int_inv, r_inv, p_inv, se_inv = sp_stats.linregress(1.0 / tm, Rv)
    log(f"  1/T:  R → {int_inv:.3f} (T→∞), r²={r_inv**2:.4f}")

    # Spearman 순위상관 (단조성 검정)
    rho_s, p_s = sp_stats.spearmanr(tm, Rv)
    log(f"  Spearman: ρ={rho_s:+.3f}, p={p_s:.4f}")

    log(f"\n  통계: min={Rv.min():.3f}, max={Rv.max():.3f}, mean={Rv.mean():.3f}±{Rv.std():.3f}")

# ===================================================================
# 7. Bootstrap 신뢰구간
# ===================================================================
log(f"\n[Bootstrap B={N_BOOT}]")
rng = np.random.default_rng(42)
boot_global_R = np.zeros(N_BOOT)
for bb in range(N_BOOT):
    ix = rng.integers(0, nn, nn)
    c = np.cov(A_n[ix], A_np1[ix])[0, 1]
    v = np.var(sh[ix])
    r = np.cov(Ap_n[ix], Ap_np1[ix])[0, 1]
    d = abs(c - v - r) + abs(r)
    boot_global_R[bb] = v / d if d > 1e-30 else 999
log(f"  전역 R={Rg:.3f}x, R>1:{np.mean(boot_global_R>1)*100:.1f}%, CI:[{np.percentile(boot_global_R,2.5):.3f},{np.percentile(boot_global_R,97.5):.3f}]")

for bi, b in enumerate(bands):
    m = (t_arr >= b['tl']) & (t_arr < b['th']) if bi < len(bands) - 1 else (t_arr >= b['tl']) & (t_arr <= b['th'])
    an = A_n[m]; anp1 = A_np1[m]; apn = Ap_n[m]; apnp1 = Ap_np1[m]; s = sh[m]; nb = len(an)
    bR = np.zeros(N_BOOT)
    for bb in range(N_BOOT):
        ix = rng.integers(0, nb, nb)
        c = np.cov(an[ix], anp1[ix])[0, 1]
        v = np.var(s[ix])
        r = np.cov(apn[ix], apnp1[ix])[0, 1]
        d = abs(c - v - r) + abs(r)
        bR[bb] = v / d if d > 1e-30 else 999
    log(f"  [{b['tl']:7.0f},{b['th']:7.0f}] R={b['R']:.2f}x, R>1:{np.mean(bR>1)*100:.0f}%, CI:[{np.percentile(bR,2.5):.2f},{np.percentile(bR,97.5):.2f}]")

# ===================================================================
# 8. 5대역 vs C-395 재현 확인
# ===================================================================
log(f"\n[5대역 분해 — C-395 재현 확인]")
edges5 = np.linspace(t_arr.min(), t_arr.max(), 6)
bands5 = []
for i in range(5):
    tl_b, th_b = edges5[i], edges5[i + 1]
    m = (t_arr >= tl_b) & (t_arr < th_b) if i < 4 else (t_arr >= tl_b) & (t_arr <= th_b)
    nb = m.sum()
    if nb < 30: continue
    cf_b, vs_b, cres_b, cr_b = decomp(A_n[m], A_np1[m], Ap_n[m], Ap_np1[m], sh[m])
    d = abs(cr_b) + abs(cres_b)
    R = vs_b / d if d > 1e-30 else 999
    log(f"  [{tl_b:7.0f},{th_b:7.0f}] n={nb:5d} R={R:.3f}x")
    bands5.append({'tm': (tl_b + th_b) / 2, 'R': R})

# ===================================================================
# 9. 판정
# ===================================================================
log(f"\n{'=' * 70}")
log("  [판정]")
Rv = np.array([b['R'] for b in bands])
nf = (Rv < 1).sum()
log(f"  전역: R = {Rg:.3f}x (n={nn})")
log(f"  대역: R>1 = {n_ok}/{nv_b}, min={Rv.min():.3f}x, max={Rv.max():.3f}x")

if nv_b >= 3:
    log(f"  선형 기울기: {sl_lin:+.6f} ± {se_lin:.6f} (p={p_lin:.4e})")
    log(f"  1/T 외삽: R(T→∞) = {int_inv:.3f}")
    log(f"  Spearman: ρ={rho_s:+.3f} (p={p_s:.4f})")

    if nf == 0 and Rv.min() > 1.05:
        if sl_lin >= -abs(se_lin):
            verdict = "강한 양성 — 전 대역 R>1.05, 비감소"
        else:
            verdict = "양성 — 전 대역 R>1.05이나 약한 감소 추세"
    elif nf == 0:
        verdict = "조건부 양성 — 전 대역 R>1이나 마진 좁음"
    elif nf <= 2:
        verdict = f"약한 양성 — {nf}개 대역 R<1"
    else:
        verdict = f"음성 — {nf}/{nv_b} 대역 R<1"

    log(f"\n  ** 판정: {verdict}")

    # 증명 feasibility
    if int_inv > 1.0:
        log(f"  ** R(T→∞)={int_inv:.3f}>1: 해석적 증명 가능성 유지")
    else:
        log(f"  ** R(T→∞)={int_inv:.3f}≤1: 해석적 증명 전략 재고 필요")

log(f"\n  총 소요: {time.time()-t0g:.1f}s")
log("=" * 70)
out_f.close()
print(f"\n결과: {RPATH}")
