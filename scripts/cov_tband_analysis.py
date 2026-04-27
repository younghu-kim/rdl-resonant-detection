#!/usr/bin/env python3
"""[C-398] Path C R(T) T-대역별 추세 (T=2000, 10대역)"""
import sys, os, time, math
import numpy as np
from scipy import stats as sp_stats
sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import mpmath; mpmath.mp.dps = 50
import cypari2; pari = cypari2.Pari()
pari.allocatemem(1024 * 10**6); pari.set_real_precision(80)

T_MAX, W, TRIM = 2000.0, 100, 0.20
N_BANDS, N_BOOT = 10, 2000
SIGMA_C = 0.5
RPATH = os.path.expanduser('~/Desktop/gdl_unified/results/path_c_asymptotic_c398.txt')
out_f = open(RPATH, 'w')
def log(m=''):
    print(m, flush=True); out_f.write(m+'\n'); out_f.flush()

def gamma_corr(t0):
    s = mpmath.mpc(SIGMA_C, t0); arg = (s + 0) / 2
    return float(mpmath.im(mpmath.digamma(arg))/2), float(mpmath.re(mpmath.psi(1, arg))/4)

log("="*70)
log("[C-398] Path C R(T) T-대역별 추세 (T=2000, 10대역)")
log("="*70)
t0g = time.time()

pari('L=lfuncreate(1)')
pari(f'Li=lfuninit(L,[0,{int(T_MAX)+5}])')
pari(f'zv=lfunzeros(Li,{T_MAX})')
n_z = int(str(pari('#zv')))
raw = str(pari('zv')).strip('[]')
zeros = sorted([float(p) for p in raw.replace(',',' ').split() if float(p)>0.5])
zeros = np.array(zeros)
log(f"  {len(zeros)}개 영점, [{zeros[0]:.1f}, {zeros[-1]:.1f}], {time.time()-t0g:.0f}s")

N = len(zeros)
im_g = np.zeros(N); re_g = np.zeros(N)
for i in range(N):
    im_g[i], re_g[i] = gamma_corr(zeros[i])
log(f"  Gamma 보정 {time.time()-t0g:.0f}s")

x = zeros; lo_w, hi_w = W, N-W-1; nv = hi_w-lo_w
g = x[lo_w+1:hi_w+1]-x[lo_w:hi_w]; g = np.where(np.abs(g)<1e-15,1e-15,g)
dL = x[lo_w:hi_w]-x[lo_w-1:hi_w-1]; dL = np.where(np.abs(dL)<1e-15,1e-15,dL)
dR = x[lo_w+2:hi_w+2]-x[lo_w+1:hi_w+1]; dR = np.where(np.abs(dR)<1e-15,1e-15,dR)

Ht_n = np.zeros(nv); Ht_np1 = np.zeros(nv)
S1_n = -1.0/g+1.0/dL; S1_np1 = -1.0/dR+1.0/g
for j in range(2, W+1):
    idx = np.arange(lo_w, hi_w)
    dr = x[idx+j]-x[idx]; dr = np.where(np.abs(dr)<1e-15,1e-15,dr)
    dl = x[idx]-x[idx-j]; dl = np.where(np.abs(dl)<1e-15,1e-15,dl)
    Ht_n += 1.0/dr**2+1.0/dl**2; S1_n += -1.0/dr+1.0/dl
    ix1 = np.arange(lo_w+1, hi_w+1)
    dr1 = x[ix1+j]-x[ix1]; dr1 = np.where(np.abs(dr1)<1e-15,1e-15,dr1)
    dl1 = x[ix1]-x[ix1-j]; dl1 = np.where(np.abs(dl1)<1e-15,1e-15,dl1)
    Ht_np1 += 1.0/dr1**2+1.0/dl1**2; S1_np1 += -1.0/dr1+1.0/dl1

S1L_n = S1_n-im_g[lo_w:hi_w]; S1L_np1 = S1_np1-im_g[lo_w+1:hi_w+1]
H_n = 1.0/g**2+1.0/dL**2+Ht_n+re_g[lo_w:hi_w]
H_np1 = 1.0/g**2+1.0/dR**2+Ht_np1+re_g[lo_w+1:hi_w+1]
A_n = S1L_n**2+2*H_n; A_np1 = S1L_np1**2+2*H_np1
sh = 2.0/g**2; Ap_n = A_n-sh; Ap_np1 = A_np1-sh
t_arr = x[lo_w:hi_w]

tl = int(nv*TRIM); th = nv-int(nv*TRIM); sl = slice(tl,th)
mask = np.isfinite(A_n[sl])&np.isfinite(A_np1[sl])&(A_n[sl]>0)&(A_np1[sl]>0)
A_n=A_n[sl][mask]; A_np1=A_np1[sl][mask]
Ap_n=Ap_n[sl][mask]; Ap_np1=Ap_np1[sl][mask]
sh=sh[sl][mask]; t_arr=t_arr[sl][mask]
nn=len(A_n)
log(f"  {nn} 쌍, compute {time.time()-t0g:.0f}s")

def decomp(a1,a2,p1,p2,s):
    cf=np.cov(a1,a2)[0,1]; vs=np.var(s); cr=np.cov(p1,p2)[0,1]
    return cf,vs,cf-vs-cr,cr

cf,vs,cr,cres = decomp(A_n,A_np1,Ap_n,Ap_np1,sh)
Rg = vs/(abs(cr)+abs(cres)) if abs(cr)+abs(cres)>1e-30 else 999
log(f"\n  전체: Cov={cf:.4e}, Var(sh)={vs:.4e}({vs/cf*100:.1f}%), Cross={cr:.4e}({cr/cf*100:.1f}%), Resid={cres:.4e}({cres/cf*100:.1f}%)")
log(f"  전역 R = {Rg:.3f}x")

edges = np.linspace(t_arr.min(), t_arr.max(), N_BANDS+1)
log(f"\n{'대역':>24s} | {'n':>5s} | {'Var%':>6s} | {'Cross%':>7s} | {'Resid%':>7s} | {'R':>6s} | ok")
log("-"*80)
bands = []
for i in range(N_BANDS):
    tl_b, th_b = edges[i], edges[i+1]
    m = (t_arr>=tl_b)&(t_arr<th_b) if i<N_BANDS-1 else (t_arr>=tl_b)&(t_arr<=th_b)
    nb = m.sum()
    if nb<20: log(f"  [{tl_b:.0f},{th_b:.0f}] | {nb:5d} | 부족"); continue
    cf_b,vs_b,cr_b,cres_b = decomp(A_n[m],A_np1[m],Ap_n[m],Ap_np1[m],sh[m])
    d = abs(cr_b)+abs(cres_b)
    R = vs_b/d if d>1e-30 else 999
    ok = "Y" if R>1 else "N"
    log(f"  [{tl_b:7.0f},{th_b:7.0f}] | {nb:5d} | {vs_b/cf_b*100:5.1f}% | {cr_b/cf_b*100:+6.1f}% | {cres_b/cf_b*100:+6.1f}% | {R:5.2f}x | {ok}")
    bands.append({'tl':tl_b,'th':th_b,'tm':(tl_b+th_b)/2,'n':nb,'R':R})

nv_b=len(bands); n_ok=sum(1 for b in bands if b['R']>1)
log(f"\n  R>1: {n_ok}/{nv_b}")

if nv_b>=3:
    tm=np.array([b['tm'] for b in bands]); Rv=np.array([b['R'] for b in bands])
    sl_lin,int_lin,r_lin,p_lin,se_lin = sp_stats.linregress(tm,Rv)
    log(f"\n  선형: R = {sl_lin:+.6f}*T + {int_lin:.4f}, r2={r_lin**2:.4f}, p={p_lin:.4e}")
    sl_log,int_log,r_log,p_log,se_log = sp_stats.linregress(np.log(tm),Rv)
    log(f"  로그: R = {sl_log:+.4f}*log(T) + {int_log:.4f}, r2={r_log**2:.4f}")
    sl_inv,int_inv,r_inv,p_inv,se_inv = sp_stats.linregress(1.0/tm,Rv)
    log(f"  1/T:  R -> {int_inv:.3f} (T->inf), r2={r_inv**2:.4f}")
    log(f"  통계: min={Rv.min():.3f}, max={Rv.max():.3f}, mean={Rv.mean():.3f}, std={Rv.std():.3f}")

# Bootstrap
log(f"\n  [Bootstrap B={N_BOOT}]")
rng = np.random.default_rng(42)
for bi,b in enumerate(bands):
    m = (t_arr>=b['tl'])&(t_arr<b['th']) if bi<len(bands)-1 else (t_arr>=b['tl'])&(t_arr<=b['th'])
    an=A_n[m];anp1=A_np1[m];apn=Ap_n[m];apnp1=Ap_np1[m];s=sh[m]; nb=len(an)
    bR=np.zeros(N_BOOT)
    for bb in range(N_BOOT):
        ix=rng.integers(0,nb,nb)
        c=np.cov(an[ix],anp1[ix])[0,1]; v=np.var(s[ix]); r=np.cov(apn[ix],apnp1[ix])[0,1]
        d=abs(c-v-r)+abs(r); bR[bb]=v/d if d>1e-30 else 999
    log(f"  [{b['tl']:7.0f},{b['th']:7.0f}] R={b['R']:.2f}x, R>1:{np.mean(bR>1)*100:.0f}%, CI:[{np.percentile(bR,2.5):.2f},{np.percentile(bR,97.5):.2f}]")

log(f"\n{'='*70}")
log("  [판정]")
Rv=np.array([b['R'] for b in bands]); nf=(Rv<1).sum()
log(f"  R>1: {n_ok}/{nv_b}, min={Rv.min():.3f}x, 전역={Rg:.3f}x")
if nv_b>=3:
    log(f"  기울기: {sl_lin:+.6f}+/-{se_lin:.6f} (p={p_lin:.4e}), T->inf: {int_inv:.3f}")
    if nf==0 and Rv.min()>1.05:
        v="강한 양성" if sl_lin>=-abs(se_lin) else "양성"
    elif nf==0: v="조건부 양성 (마진 좁음)"
    elif nf<=2: v=f"약한 양성 ({nf}개 R<1)"
    else: v=f"음성 ({nf}개 R<1)"
    log(f"\n  ** 판정: {v}")
    log(f"  ** T=2000 범위. T>2000 확장은 PARI 최적화 후 별도 실행.")
log(f"\n총 소요: {time.time()-t0g:.1f}s")
log("="*70)
out_f.close()
print(f"\n결과: {RPATH}")
