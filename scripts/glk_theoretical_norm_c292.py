#!/usr/bin/env python3
"""
=============================================================================
[사이클 #292] GL(3)~GL(6) 이론적 밀도 정규화 재측정
=============================================================================

목적:
  C-291에서 발견: GL(2) 경험적 d_bar vs 이론적 d_bar 차이가 ρ를 -0.93→-0.66으로 약화.
  GL(3)~GL(6)도 이론적 d_bar 사용 시 ρ≈-0.93이 되는가?

  이론적 밀도 추정: W-zero 윈도우 스무딩
    d̄_smooth(t_i) = W / (zeros[i+W//2] - zeros[i-W//2])
  이것은 지역 밀도 평균으로, 개별 갭과 상관이 없음.

설계:
  - Sym^k(11a1) k=2,3,4,5 (GL(3)~GL(6))
  - Sym^2(37a1), Sym^3(37a1) (도도체 N 독립성 확인)
  - T=[3,35], N_MAX=300, trim 20%
  - 각 곡선: smooth d_bar vs empirical d_bar ρ_S 모두 출력
  - C-280/C-275 기존 경험적 ρ 비교

결과: results/glk_theoretical_norm_c292.txt
"""

import sys, os, math, time
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 30  # T<35 → 30자리 충분

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(1024**3)  # 1 GB
    pari.set_real_precision(100)
    print("cypari2 OK (1 GB, precision=100)", flush=True)
except Exception as e:
    print(f"FATAL: {e}")
    sys.exit(1)

# ── 설정 ──────────────────────────────────────────────────────────────────
T_MAX = 35.0
T_MIN = 3.0
N_MAX = 300
TRIM_FRAC = 0.20
W_SMOOTH = 40   # 스무딩 윈도우 (zeros 개수)

# C-291 기준값
REF_ZETA_THEO  = -0.929   # ζ(s) 이론적 d_bar
REF_ZETA_EMP   = -0.636   # ζ(s) 경험적 d_bar
REF_GL2_THEO   = -0.934   # GL(2) 이론적 d_bar (평균)
REF_GL2_EMP    = -0.669   # GL(2) 경험적 d_bar (평균)
GUE_LOCAL_REF  = -0.967   # GUE local (C-283)
GUE_FULL_REF   = -0.912   # GUE full (C-283)

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/glk_theoretical_norm_c292.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)
out_f = open(RESULT_PATH, 'w')


def log(msg=''):
    print(msg, flush=True)
    out_f.write(str(msg) + '\n')
    out_f.flush()


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


# ── 영점 수집 ─────────────────────────────────────────────────────────────

def get_zeros_sympow(E_coeffs, k, t_max):
    """PARI lfunsympow(E, k)로 Sym^k(E) 영점 수집."""
    t0 = time.time()
    pari(f'E_c292 = ellinit({E_coeffs})')
    pari(f'L_c292 = lfunsympow(E_c292, {k})')
    pari(f'Li_c292 = lfuninit(L_c292, [0, {t_max + 2}])')
    pari(f'zv_c292 = lfunzeros(Li_c292, {t_max})')
    n = int(str(pari('#zv_c292')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv_c292[{i}]'))
        if not math.isnan(t) and t > 0.5:
            zeros.append(t)
    zeros = sorted(zeros)
    log(f"  {len(zeros)}개 영점, t∈[{zeros[0]:.3f},{zeros[-1]:.3f}]  ({time.time()-t0:.1f}s)")
    if len(zeros) == 0:
        log("  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
    return zeros


# ── A_bare 계산 ───────────────────────────────────────────────────────────

def compute_A_bare(all_zeros, idx, n_max=N_MAX):
    """A_bare = S1^2 + 2*H1 (±n_max, same-side)"""
    gamma_0 = all_zeros[idx]
    n_total = len(all_zeros)
    S1 = 0.0
    H1 = 0.0
    for k in range(max(0, idx - n_max), min(n_total, idx + n_max + 1)):
        if k == idx:
            continue
        dg = gamma_0 - all_zeros[k]
        if abs(dg) < 1e-15:
            continue
        S1 += 1.0 / dg
        H1 += 1.0 / (dg * dg)
    return S1 * S1 + 2.0 * H1


# ── 스무딩 밀도 ───────────────────────────────────────────────────────────

def smooth_density(zeros, idx, W=W_SMOOTH):
    """
    윈도우 W 스무딩 밀도:
    d̄(t_idx) = W_actual / (zeros[idx+W//2] - zeros[idx-W//2])
    이것은 W개 영점에 걸친 평균 밀도로, 개별 갭과 상관이 없음.
    """
    n = len(zeros)
    lo = max(0, idx - W // 2)
    hi = min(n - 1, idx + W // 2)
    n_span = hi - lo
    if n_span < 4 or abs(zeros[hi] - zeros[lo]) < 1e-10:
        return float('nan')
    return n_span / (zeros[hi] - zeros[lo])


# ── 분석 루틴 ─────────────────────────────────────────────────────────────

def analyze_curve(label, all_zeros, ref_emp_rho=None):
    """
    trim 20% → gap + smooth density 계산 → Spearman ρ (smooth vs empirical)
    """
    n_all = len(all_zeros)
    log(f"\n  [{label}] n_zeros={n_all}")

    if n_all < 30:
        log(f"  영점 부족 ({n_all} < 30)")
        return None

    # A_bare 계산
    data = []
    n_fail = 0
    for i in range(n_all):
        try:
            A = compute_A_bare(all_zeros, i)
            if math.isnan(A) or math.isinf(A) or A <= 0:
                n_fail += 1
                continue
            data.append({'idx': i, 't': all_zeros[i], 'A_bare': A})
        except Exception as e:
            n_fail += 1
            if n_fail <= 3:
                log(f"    WARNING i={i}: {e}")

    n_data = len(data)
    log(f"  A_bare 계산: {n_data}/{n_all} 유효, 실패={n_fail}")

    if n_data < 20:
        log(f"  데이터 부족")
        return None

    # trim 20%
    lo = int(n_data * TRIM_FRAC)
    hi = int(n_data * (1.0 - TRIM_FRAC))
    trimmed = data[lo:hi]
    n_trim = len(trimmed)

    if n_trim < 10:
        log(f"  trim 후 부족 ({n_trim})")
        return None

    # gap + 밀도 계산
    valid_smooth = []
    valid_emp = []
    n_skip = 0

    for pos in range(len(trimmed)):
        d = trimmed[pos]
        idx_global = d['idx']
        t_n = d['t']

        # 갭: 전체 zeros에서 직전/직후
        if idx_global <= 0 or idx_global >= n_all - 1:
            n_skip += 1
            continue
        t_prev = all_zeros[idx_global - 1]
        t_next = all_zeros[idx_global + 1]
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        if gap_r <= 0 or gap_l <= 0:
            n_skip += 1
            continue
        gap_min = min(gap_r, gap_l)

        # 스무딩 밀도 (이론적)
        d_bar_smooth = smooth_density(all_zeros, idx_global, W=W_SMOOTH)
        # 경험적 밀도
        d_bar_emp = 2.0 / (t_next - t_prev)

        if math.isnan(d_bar_smooth) or d_bar_smooth <= 0:
            n_skip += 1
            continue

        valid_smooth.append({
            'A_bare': d['A_bare'],
            'gap_min_gue_smooth': gap_min * d_bar_smooth,
            'gap_min_gue_emp': gap_min * d_bar_emp,
        })

    n_final = len(valid_smooth)
    log(f"  trim={n_trim}, 유효={n_final}, 스킵={n_skip}")

    if n_final < 10:
        log(f"  최종 데이터 부족")
        return None

    A_arr = np.array([d['A_bare'] for d in valid_smooth])
    gm_smooth = np.array([d['gap_min_gue_smooth'] for d in valid_smooth])
    gm_emp = np.array([d['gap_min_gue_emp'] for d in valid_smooth])

    # NaN 필터
    mask = np.isfinite(A_arr) & np.isfinite(gm_smooth) & np.isfinite(gm_emp) & (A_arr > 0)
    A_arr  = A_arr[mask]
    gm_smooth = gm_smooth[mask]
    gm_emp = gm_emp[mask]
    n_final = len(A_arr)

    if n_final < 5:
        log(f"  NaN 필터 후 {n_final} < 5")
        return None

    rho_smooth, p_smooth = stats.spearmanr(A_arr, gm_smooth)
    rho_emp,    p_emp    = stats.spearmanr(A_arr, gm_emp)
    se = 1.0 / math.sqrt(n_final - 3) if n_final > 3 else float('nan')

    sig_s = '✅' if p_smooth < 0.001 else ('⚠️' if p_smooth < 0.05 else '❌')
    sig_e = '✅' if p_emp    < 0.001 else ('⚠️' if p_emp    < 0.05 else '❌')
    log(f"  ρ_S(A_bare, gap_min_smooth) = {rho_smooth:+.6f}  p={p_smooth:.3e}  SE={se:.4f}  {sig_s}")
    log(f"  ρ_S(A_bare, gap_min_emp)    = {rho_emp:+.6f}  p={p_emp:.3e}  SE={se:.4f}  {sig_e}")
    log(f"  Δρ (smooth-emp) = {rho_smooth - rho_emp:+.4f}  (기대: ≈{-0.27:+.2f})")
    if ref_emp_rho is not None:
        log(f"  참조 (경험적 C-280/etc): {ref_emp_rho}")

    return {
        'label': label,
        'n_zeros': n_all,
        'n_final': n_final,
        'rho_smooth': float(rho_smooth),
        'rho_emp': float(rho_emp),
        'p_smooth': float(p_smooth),
        'p_emp': float(p_emp),
        'se': float(se),
    }


# ===== 메인 ===================================================================
t_start = time.time()

log("=" * 70)
log("  [C-292] GL(3)~GL(6) 이론적 밀도 정규화 재측정")
log("=" * 70)
log(f"  T=[{T_MIN},{T_MAX}], N_MAX={N_MAX}, trim {TRIM_FRAC*100:.0f}%")
log(f"  smooth 윈도우 W={W_SMOOTH} zeros")
log(f"  참조: ζ(s) 이론적={REF_ZETA_THEO}, GL(2) 이론적≈{REF_GL2_THEO}")
log()

# 타원곡선 계수
E11A1 = '[0,-1,1,-10,-20]'
E37A1 = '[0,0,1,-1,0]'  # 37a1 (C-285에서 검증됨)

all_results = []

# ── GL(3) Sym²(11a1) ──────────────────────────────────────────────────────
log("━" * 70)
log("  GL(3): Sym²(11a1)  — lfunsympow(E, 2)")
log("━" * 70)
try:
    z_gl3_11 = get_zeros_sympow(E11A1, 2, T_MAX)
    if z_gl3_11:
        r = analyze_curve("GL(3) Sym²(11a1)", z_gl3_11, ref_emp_rho=-0.492)
        if r:
            r['degree'] = 3
            all_results.append(r)
except Exception as e:
    log(f"  ⚠️ GL(3) Sym²(11a1) 실패: {e}")
    import traceback; traceback.print_exc()

log()

# ── GL(3) Sym²(37a1) ──────────────────────────────────────────────────────
log("━" * 70)
log("  GL(3): Sym²(37a1)  — lfunsympow(E, 2)")
log("━" * 70)
try:
    z_gl3_37 = get_zeros_sympow(E37A1, 2, T_MAX)
    if z_gl3_37:
        r = analyze_curve("GL(3) Sym²(37a1)", z_gl3_37, ref_emp_rho=-0.474)
        if r:
            r['degree'] = 3
            all_results.append(r)
except Exception as e:
    log(f"  ⚠️ GL(3) Sym²(37a1) 실패: {e}")
    import traceback; traceback.print_exc()

log()

# ── GL(4) Sym³(11a1) ──────────────────────────────────────────────────────
log("━" * 70)
log("  GL(4): Sym³(11a1)  — lfunsympow(E, 3)")
log("━" * 70)
try:
    z_gl4_11 = get_zeros_sympow(E11A1, 3, T_MAX)
    if z_gl4_11:
        r = analyze_curve("GL(4) Sym³(11a1)", z_gl4_11, ref_emp_rho=-0.520)
        if r:
            r['degree'] = 4
            all_results.append(r)
except Exception as e:
    log(f"  ⚠️ GL(4) Sym³(11a1) 실패: {e}")
    import traceback; traceback.print_exc()

log()

# ── GL(4) Sym³(37a1) ──────────────────────────────────────────────────────
log("━" * 70)
log("  GL(4): Sym³(37a1)  — lfunsympow(E, 3)")
log("━" * 70)
try:
    z_gl4_37 = get_zeros_sympow(E37A1, 3, T_MAX)
    if z_gl4_37:
        r = analyze_curve("GL(4) Sym³(37a1)", z_gl4_37, ref_emp_rho=-0.514)
        if r:
            r['degree'] = 4
            all_results.append(r)
except Exception as e:
    log(f"  ⚠️ GL(4) Sym³(37a1) 실패: {e}")
    import traceback; traceback.print_exc()

log()

# ── GL(5) Sym⁴(11a1) ──────────────────────────────────────────────────────
log("━" * 70)
log("  GL(5): Sym⁴(11a1)  — lfunsympow(E, 4)")
log("━" * 70)
try:
    z_gl5_11 = get_zeros_sympow(E11A1, 4, T_MAX)
    if z_gl5_11:
        r = analyze_curve("GL(5) Sym⁴(11a1)", z_gl5_11, ref_emp_rho=-0.485)
        if r:
            r['degree'] = 5
            all_results.append(r)
except Exception as e:
    log(f"  ⚠️ GL(5) Sym⁴(11a1) 실패: {e}")
    import traceback; traceback.print_exc()

log()

# ── GL(6) Sym⁵(11a1) ──────────────────────────────────────────────────────
log("━" * 70)
log("  GL(6): Sym⁵(11a1)  — lfunsympow(E, 5)")
log("━" * 70)
try:
    z_gl6_11 = get_zeros_sympow(E11A1, 5, T_MAX)
    if z_gl6_11:
        r = analyze_curve("GL(6) Sym⁵(11a1)", z_gl6_11, ref_emp_rho=-0.423)
        if r:
            r['degree'] = 6
            all_results.append(r)
except Exception as e:
    log(f"  ⚠️ GL(6) Sym⁵(11a1) 실패: {e}")
    import traceback; traceback.print_exc()

log()

# ── 최종 비교표 ───────────────────────────────────────────────────────────
log("=" * 70)
log("  [최종 비교표] 이론적 vs 경험적 d_bar")
log("=" * 70)
log()
log(f"  {'L-함수':<22} {'d':>2} {'ρ_smooth':>10} {'ρ_emp':>10} {'Δρ':>8} {'n':>5}")
log("  " + "-" * 60)

# ζ(s) 및 GL(2) 참조 추가
log(f"  {'ζ(s) [C-282/C-290]':<22} {'1':>2} {REF_ZETA_THEO:>+10.4f} {REF_ZETA_EMP:>+10.4f} {REF_ZETA_THEO-REF_ZETA_EMP:>+8.4f} {'—':>5}")
log(f"  {'GL(2) 11a1 [C-291]':<22} {'2':>2} {REF_GL2_THEO:>+10.4f} {REF_GL2_EMP:>+10.4f} {REF_GL2_THEO-REF_GL2_EMP:>+8.4f} {'—':>5}")

for r in all_results:
    delta = r['rho_smooth'] - r['rho_emp']
    log(f"  {r['label']:<22} {r['degree']:>2} {r['rho_smooth']:>+10.4f} {r['rho_emp']:>+10.4f} {delta:>+8.4f} {r['n_final']:>5}")

log()

# GUE 기준
log(f"  {'GUE_local [C-283]':<22} {'—':>2} {'(≈-0.967)':>10} {GUE_LOCAL_REF:>+10.4f} {'—':>8}")
log()

# ── 판정 ──────────────────────────────────────────────────────────────────
log("=" * 70)
log("  [판정]")
log("=" * 70)

if all_results:
    rho_smooth_vals = [r['rho_smooth'] for r in all_results]
    rho_emp_vals = [r['rho_emp'] for r in all_results]
    avg_smooth = np.mean(rho_smooth_vals)
    avg_emp = np.mean(rho_emp_vals)

    log(f"\n  GL(3)~GL(6) 스무딩 정규화 평균 ρ = {avg_smooth:+.4f}")
    log(f"  GL(3)~GL(6) 경험적 정규화 평균 ρ   = {avg_emp:+.4f}")
    log(f"  ζ(s) 이론적 ρ                     = {REF_ZETA_THEO:+.4f}")
    log(f"  GL(2) 이론적 ρ 평균               = {REF_GL2_THEO:+.4f}")
    log()

    # 이론적(스무딩) 정규화에서 d-감쇠 확인
    delta_smooth = abs(REF_ZETA_THEO) - abs(avg_smooth)
    log(f"  δ(이론적): |ζ(s)| - |GL(3-6)_smooth| = {abs(REF_ZETA_THEO):.4f} - {abs(avg_smooth):.4f} = {delta_smooth:+.4f}")

    if abs(delta_smooth) < 0.05:
        log(f"\n  ★★★★★ 보편성 확인: 이론적 정규화 시 GL(1)~GL(6) 모두 ρ≈-0.93")
        log(f"         d-의존 감쇠는 정규화 아티팩트. Paper 4 핵심 수정 필요.")
    elif delta_smooth > 0.05:
        log(f"\n  ★★★★ 부분 d-감쇠 실재: δ={delta_smooth:.4f} > 0.05")
        log(f"         이론적 정규화 후에도 GL(3-6)가 ζ(s)보다 약함.")
    else:
        log(f"\n  ★★★ GL(3-6)가 ζ(s)보다 강한 상관 (δ<0). 스무딩 과대 보정 가능성.")

    # 경험적 ρ 개선 확인
    delta_emp = abs(avg_emp) - abs(avg_smooth)
    log()
    log(f"  경험적 ρ 대비 스무딩 ρ 변화: {delta_emp:+.4f}")
    if delta_emp < -0.10:
        log(f"  → 이론적 정규화가 경험적보다 ρ 크게 강화 (Δ={abs(delta_emp):.3f})")
    elif delta_emp < -0.05:
        log(f"  → 이론적 정규화가 경험적보다 ρ 중등 강화 (Δ={abs(delta_emp):.3f})")
    else:
        log(f"  → 정규화 방법 차이 미미 (GL(3-6)는 GL(2)와 다른 양상)")

else:
    log("  결과 없음")

elapsed = time.time() - t_start
log()
log(f"  총 소요: {elapsed:.1f}s ({elapsed/60:.1f}분)")

out_f.close()
log(f"\n결과 저장: {RESULT_PATH}")
log("[C-292 완료]")
