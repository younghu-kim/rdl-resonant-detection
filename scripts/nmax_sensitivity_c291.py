#!/usr/bin/env python3
"""
=============================================================================
[사이클 #291] N_MAX 감도 검증 — Paper 4 통일표 무결성 결정적 실험
=============================================================================

목적:
  Paper 4 통일표에서
    - d=1 ζ(s): N_MAX=300 (C-282) → ρ=-0.929
    - d=2 11a1: N_MAX=200 (C-282c) → ρ=-0.658
  이 0.271 차이가 N_MAX 아티팩트인가, 진정한 d-의존 감쇠인가?

실험 설계:
  - GL(2) 11a1: T_MAX=2000, trim 중앙 60%, N_MAX=200 AND N_MAX=300
  - ζ(s):       T_MAX=2000, trim 중앙 60%, N_MAX=200 AND N_MAX=300
  - 고정: d_bar = 2/(t_next-t_prev) (local GUE, C-290 방식으로 통일)
  - 모든 조건 동일, N_MAX만 변경

판정 기준 (수학자 지시):
  - GL(2) |ρ(N300) - ρ(N200)| > 0.15: N_MAX 감도 확인
  - GL(2) N_MAX=300에서 ρ > -0.85: 음성 (d-감쇠 아티팩트 확정) → 긴급 수정
  - GL(2) N_MAX=300에서 ρ < -0.75: 양성 (N_MAX 효과 있으나 d-감쇠 실재)
  - GL(2) N_MAX=300에서 ρ ∈ [-0.85, -0.75]: 중립 (부분 아티팩트)

체크리스트:
  [x] A_bare = S1^2 + 2*H1 (±n_max, same-side only)
  [x] gap_min_GUE = min(gap_r, gap_l) * (2/(t_next-t_prev))
  [x] trim 20% 양쪽 (중앙 60%)
  [x] Spearman ρ_S (scipy.stats.spearmanr)
  [x] NaN/Inf 체크
  [x] python -u
  [x] GL(2) T=2000 실패 시 T=1000 폴백
  [x] conductor 검증 (ellglobalred)
  [x] mpmath dps=80 (t≤2000)

결과: results/nmax_sensitivity_c291.txt
=============================================================================
"""

import sys
import os
import math
import time

import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 80  # t ≤ 2000에서 충분

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(2 * 1024**3)  # 2 GB (GL(2) T=2000 대비)
    pari.set_real_precision(80)
    print("cypari2 OK (2 GB, precision=80)", flush=True)
except Exception as e:
    print(f"FATAL: cypari2 로드 실패: {e}")
    sys.exit(1)

# ── 설정 ──────────────────────────────────────────────────────────────────
T_MAX_ZETA = 2000.0
T_MAX_GL2  = 2000.0
T_MAX_GL2_FALLBACK = 1000.0

CENTER_ZETA = 0.5
CENTER_GL2  = 1.0
MU_ZETA     = [0]
MU_GL2      = [0, 1]
N_COND_ZETA = 1
N_COND_GL2  = 11

TRIM_FRAC   = 0.20   # 양쪽 20% 제외 → 중앙 60%
N_MAX_LIST  = [200, 300]

# C-282c 기준값 (GL(2) N_MAX=200, T=500, trim 20%)
C282C_11A1_N200 = -0.658
# C-282 기준값 (ζ(s) N_MAX=300, T=2000, trim 20%, 이론적 d_bar)
C282_ZETA_N300  = -0.929
# C-290 기준값 (ζ(s) N_MAX=200, T=2000, trim 20%, local d_bar)
C290_ZETA_N200  = -0.636
# GUE 기준값
GUE_LOCAL_REF = -0.967  # ζ(s) 동일조건 GUE (C-283)
GUE_FULL_REF  = -0.912  # GL(2) 동일조건 GUE (C-283)

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/nmax_sensitivity_c291.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

out_lines = []


def log(msg=''):
    print(msg, flush=True)
    out_lines.append(str(msg))


def save():
    with open(RESULT_PATH, 'w') as f:
        f.write('\n'.join(out_lines) + '\n')
    log(f"  저장: {RESULT_PATH}")


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


# ── 영점 수집 ─────────────────────────────────────────────────────────────

def get_zeros_zeta(t_max):
    """PARI lfunzeros로 ζ(s) 영점 수집."""
    log(f"[ζ(s) 영점] lfuninit([0, {int(t_max)+10}]) ...")
    t0 = time.time()
    pari('L_z91 = lfuncreate(1)')
    pari(f'Li_z91 = lfuninit(L_z91, [0, {int(t_max) + 10}])')
    pari(f'zv_z91 = lfunzeros(Li_z91, {t_max})')
    n = int(str(pari('#zv_z91')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv_z91[{i}]'))
        if not math.isnan(t) and t > 0.5:
            zeros.append(t)
    zeros = sorted(zeros)
    log(f"  {len(zeros)}개, t ∈ [{zeros[0]:.4f}, {zeros[-1]:.4f}]  ({time.time()-t0:.1f}s)")
    return zeros


def get_zeros_gl2(coeffs, t_max, n_cond_expected):
    """PARI lfunzeros로 GL(2) 타원곡선 영점 수집. conductor 검증 포함."""
    log(f"[GL(2) 영점] 11a1 lfuninit([0, {int(t_max)+10}]) ...")
    t0 = time.time()
    pari(f'E_gl2c291 = ellinit({coeffs})')
    # conductor 검증
    try:
        N_actual = int(str(pari('ellglobalred(E_gl2c291)[1]')))
        log(f"  conductor 검증: {N_actual} (기대 {n_cond_expected})")
        if N_actual != n_cond_expected:
            log(f"  ⚠️ conductor 불일치! a-invariants 점검 필요")
    except Exception as e:
        log(f"  WARNING: conductor 검증 실패: {e}")
    pari(f'Li_gl2c291 = lfuninit(E_gl2c291, [0, {int(t_max) + 10}])')
    pari(f'zv_gl2c291 = lfunzeros(Li_gl2c291, {t_max})')
    n = int(str(pari('#zv_gl2c291')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv_gl2c291[{i}]'))
        if not math.isnan(t) and t > 0.5:
            zeros.append(t)
    zeros = sorted(zeros)
    log(f"  {len(zeros)}개, t ∈ [{zeros[0]:.4f}, {zeros[-1]:.4f}]  ({time.time()-t0:.1f}s)")
    if len(zeros) == 0:
        log("  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
    return zeros


# ── A_bare 계산 ───────────────────────────────────────────────────────────

def compute_A_bare(all_zeros, idx, n_max):
    """
    A_bare(γ₀) = S₁² + 2H₁
    S₁ = Σ_{|k-idx|≤n_max, k≠idx} 1/(γ₀ - γₖ)  (same-side only)
    H₁ = Σ_{|k-idx|≤n_max, k≠idx} 1/(γ₀ - γₖ)²
    """
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
    A_bare = S1 * S1 + 2.0 * H1
    return A_bare


# ── 분석 루틴 ─────────────────────────────────────────────────────────────

def run_nmax_analysis(label, all_zeros, t_max):
    """
    전체 영점 → central 60% trim → gap_min_GUE 계산
    → N_MAX=200, 300 각각 A_bare → Spearman ρ_S
    """
    n_all = len(all_zeros)
    log()
    log(f"  [{label}] n_zeros={n_all}, T_MAX={t_max}")
    log(f"  trim: 중앙 60% (양쪽 20%)")

    if n_all < 50:
        log("  FATAL: 영점 부족 (<50)")
        return None

    # ── trim 인덱스 ────────────────────────────────────────────────
    lo = int(n_all * TRIM_FRAC)
    hi = int(n_all * (1.0 - TRIM_FRAC))
    log(f"  lo={lo}, hi={hi}, n_trim_range={hi-lo}")

    # ── gap 계산 (전체 zeros 사용, 인접 이웃 참조) ─────────────────
    gapped = []
    for idx in range(lo, hi):
        # 경계 보호: 전체 zeros에서 idx-1, idx+1이 존재해야 함
        if idx <= 0 or idx >= n_all - 1:
            continue
        t_n    = all_zeros[idx]
        t_prev = all_zeros[idx - 1]
        t_next = all_zeros[idx + 1]
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        if gap_r <= 0 or gap_l <= 0:
            continue
        # local GUE 정규화: d_bar = 2/(t_next - t_prev)
        d_bar = 2.0 / (t_next - t_prev)
        gap_min_gue = min(gap_r, gap_l) * d_bar
        gapped.append({'idx': idx, 't': t_n, 'gap_min_gue': gap_min_gue})

    n_gapped = len(gapped)
    log(f"  gap 계산 완료: {n_gapped}개")
    if n_gapped < 30:
        log("  FATAL: 내부 영점 부족 (<30)")
        return None
    log(f"  t ∈ [{gapped[0]['t']:.3f}, {gapped[-1]['t']:.3f}]")

    # ── N_MAX별 A_bare 계산 + Spearman ────────────────────────────
    results = {}
    for n_max in N_MAX_LIST:
        log(f"\n  [N_MAX={n_max}] A_bare 계산 ({n_gapped}개)...")
        t0 = time.time()

        A_arr = []
        gm_arr = []
        n_fail = 0

        for entry in gapped:
            idx = entry['idx']
            try:
                A = compute_A_bare(all_zeros, idx, n_max)
                if math.isnan(A) or math.isinf(A) or A <= 0:
                    n_fail += 1
                    continue
                A_arr.append(A)
                gm_arr.append(entry['gap_min_gue'])
            except Exception as e:
                n_fail += 1
                if n_fail <= 3:
                    log(f"    WARNING idx={idx}: {e}")

        elapsed = time.time() - t0
        n_valid = len(A_arr)
        log(f"    완료: {n_valid}/{n_gapped} 유효, 실패={n_fail}, ({elapsed:.1f}s)")

        if n_fail > n_gapped // 2:
            log(f"    ⚠️ 실패율 {n_fail}/{n_gapped} — 절반 초과, 결과 신뢰도 낮음")

        if n_valid < 20:
            log(f"    ⚠️ 유효 데이터 부족 ({n_valid}), 건너뜀")
            continue

        A_np  = np.array(A_arr)
        gm_np = np.array(gm_arr)

        # NaN/Inf 최종 필터
        mask = np.isfinite(A_np) & np.isfinite(gm_np) & (A_np > 0)
        A_np  = A_np[mask]
        gm_np = gm_np[mask]
        n_final = len(A_np)

        if n_final < 10:
            log(f"    ⚠️ NaN 필터 후 {n_final} < 10")
            continue

        rho_S, p_val = stats.spearmanr(A_np, gm_np)
        se = 1.0 / math.sqrt(n_final - 3) if n_final > 3 else float('nan')

        sig = '✅' if p_val < 0.001 else ('⚠️' if p_val < 0.05 else '❌')
        log(f"    ρ_S(A_bare, gap_min_GUE) = {rho_S:+.6f}  p={p_val:.3e}  "
            f"SE={se:.4f}  n={n_final}  {sig}")

        results[n_max] = {
            'rho_S': float(rho_S),
            'p_val': float(p_val),
            'se': float(se),
            'n_valid': n_final,
        }

    # N_MAX 감도 요약
    if 200 in results and 300 in results:
        delta = abs(results[300]['rho_S'] - results[200]['rho_S'])
        log(f"\n  N_MAX 감도: |ρ(N300) - ρ(N200)| = {delta:.4f}  "
            f"({'감도 확인 ⚠️' if delta > 0.15 else '감도 미미'})")

    return results


# ===== 메인 ===================================================================
t_start = time.time()

log("=" * 70)
log("  사이클 #291 — N_MAX 감도 검증 (Paper 4 무결성)")
log("=" * 70)
log(f"  N_MAX 목록: {N_MAX_LIST}")
log(f"  trim: 중앙 60% (양쪽 {int(TRIM_FRAC*100)}% 제거)")
log(f"  GUE 정규화: local d_bar = 2/(t_next-t_prev)")
log(f"  참조: C-282 ζ(s) N_MAX=300 → ρ={C282_ZETA_N300}")
log(f"  참조: C-290 ζ(s) N_MAX=200 → ρ={C290_ZETA_N200}")
log(f"  참조: C-282c GL(2) 11a1 N_MAX=200, T=500 → ρ={C282C_11A1_N200}")
log()

all_results = {}  # {'zeta': {200: {...}, 300: {...}}, 'gl2_11a1': {...}}

# ── 1. ζ(s) ──────────────────────────────────────────────────────────────
log("━" * 70)
log("  [1] ζ(s)  — CENTER=0.5, mu=[0], N_cond=1")
log("━" * 70)
try:
    zeros_zeta = get_zeros_zeta(T_MAX_ZETA)
    r_zeta = run_nmax_analysis("ζ(s)", zeros_zeta, T_MAX_ZETA)
    if r_zeta:
        all_results['zeta'] = r_zeta
except Exception as e:
    log(f"  ⚠️ ζ(s) 실패: {e}")
    import traceback; traceback.print_exc()

log()

# ── 2. GL(2) 11a1 ─────────────────────────────────────────────────────────
log("━" * 70)
log("  [2] GL(2) 11a1  — CENTER=1.0, mu=[0,1], N_cond=11")
log("━" * 70)
gl2_zeros = None
gl2_t_max_used = None

try:
    gl2_zeros = get_zeros_gl2('[0,-1,1,-10,-20]', T_MAX_GL2, N_COND_GL2)
    gl2_t_max_used = T_MAX_GL2
except Exception as e:
    log(f"  ⚠️ GL(2) T={T_MAX_GL2} 실패: {e}")
    import traceback; traceback.print_exc()
    log(f"  → T={T_MAX_GL2_FALLBACK} 폴백 시도...")
    try:
        gl2_zeros = get_zeros_gl2('[0,-1,1,-10,-20]', T_MAX_GL2_FALLBACK, N_COND_GL2)
        gl2_t_max_used = T_MAX_GL2_FALLBACK
        log(f"  폴백 성공 (T={T_MAX_GL2_FALLBACK})")
    except Exception as e2:
        log(f"  ⚠️ 폴백도 실패: {e2}")

if gl2_zeros and len(gl2_zeros) >= 50:
    try:
        r_gl2 = run_nmax_analysis(f"GL(2) 11a1 [T={gl2_t_max_used}]",
                                   gl2_zeros, gl2_t_max_used)
        if r_gl2:
            all_results['gl2_11a1'] = r_gl2
    except Exception as e:
        log(f"  ⚠️ GL(2) 분석 실패: {e}")
        import traceback; traceback.print_exc()
elif gl2_zeros is not None:
    log(f"  ⚠️ GL(2) 영점 {len(gl2_zeros)}개 — 분석 최소치 미달")

# ── 3. 최종 비교표 ────────────────────────────────────────────────────────
log()
log("=" * 70)
log("  [최종 비교표] N_MAX 감도 — Paper 4 핵심 결과")
log("=" * 70)
log()

header = f"  {'L-함수':<18} {'N_MAX':>6} {'ρ_S':>10} {'SE':>7} {'n':>7}"
log(header)
log("  " + "-" * 50)

summary_rows = []
for key, label, gue_ref in [
        ('zeta',      'ζ(s)',        GUE_LOCAL_REF),
        ('gl2_11a1',  'GL(2) 11a1',  GUE_FULL_REF)]:
    res = all_results.get(key, {})
    if not res:
        log(f"  {label}: 결과 없음")
        continue
    for nm in N_MAX_LIST:
        r = res.get(nm, {})
        rho = r.get('rho_S', float('nan'))
        se  = r.get('se',    float('nan'))
        nv  = r.get('n_valid', 0)
        row = f"  {label:<18} {nm:>6} {rho:>+10.6f} {se:>7.4f} {nv:>7}"
        log(row)
        summary_rows.append((key, label, nm, rho, se, nv, gue_ref))

log()

# δ_arith = |GUE_ref| - |ρ| (how much weaker than GUE)
log(f"  {'L-함수':<18} {'N_MAX':>6} {'δ_arith':>10}  (|GUE_ref| - |ρ|)")
log("  " + "-" * 40)
for (key, label, nm, rho, se, nv, gue_ref) in summary_rows:
    if not math.isnan(rho):
        delta = abs(gue_ref) - abs(rho)
        log(f"  {label:<18} {nm:>6} {delta:>+10.4f}")

log()

# ── 4. 판정 ───────────────────────────────────────────────────────────────
log("=" * 70)
log("  [판정]")
log("=" * 70)

# ── 4a. ζ(s) 재현 검증 ───
log()
log("  ζ(s) N_MAX 감도 (C-282/C-290 재현 확인):")
if 'zeta' in all_results:
    rz = all_results['zeta']
    rho_z200 = rz.get(200, {}).get('rho_S', float('nan'))
    rho_z300 = rz.get(300, {}).get('rho_S', float('nan'))
    log(f"    N_MAX=200 → ρ = {rho_z200:+.4f}  (C-290 참조: {C290_ZETA_N200})")
    log(f"    N_MAX=300 → ρ = {rho_z300:+.4f}  (C-282 참조: {C282_ZETA_N300}, 이론적 d_bar 사용)")
    if not (math.isnan(rho_z200) or math.isnan(rho_z300)):
        dz = abs(rho_z300 - rho_z200)
        log(f"    |Δρ_ζ(N200→N300)| = {dz:.4f}")
        if dz > 0.20:
            log(f"    → ζ(s) N_MAX 감도 확인 (Δ={dz:.4f}) ✅")
        elif dz > 0.10:
            log(f"    → ζ(s) N_MAX 감도 중등 (Δ={dz:.4f}) ⚠️")
        else:
            log(f"    → ζ(s) N_MAX 감도 미미 (Δ={dz:.4f})")
            log(f"    → 주의: C-282와 C-290의 차이는 d_bar 공식(이론 vs 로컬) 때문일 수 있음")
else:
    log("    결과 없음")

# ── 4b. GL(2) 핵심 판정 ───
log()
log("  GL(2) 11a1 N_MAX 감도 (핵심 판정):")
if 'gl2_11a1' in all_results:
    rg = all_results['gl2_11a1']
    rho_g200 = rg.get(200, {}).get('rho_S', float('nan'))
    rho_g300 = rg.get(300, {}).get('rho_S', float('nan'))
    log(f"    N_MAX=200 → ρ = {rho_g200:+.4f}  (C-282c 참조: {C282C_11A1_N200}, T=500)")
    log(f"    N_MAX=300 → ρ = {rho_g300:+.4f}  (신규 측정)")
    log(f"    GUE_full 기준: {GUE_FULL_REF}")

    if not (math.isnan(rho_g200) or math.isnan(rho_g300)):
        dg = abs(rho_g300 - rho_g200)
        log(f"    |Δρ_GL2(N200→N300)| = {dg:.4f}")

        # 수학자 판정 기준 적용
        log()
        if dg > 0.15:
            log(f"    ★ N_MAX 감도 확인: |Δ|={dg:.4f} > 0.15")

        if rho_g300 > -0.75:
            # ρ=-0.65처럼 GUE에서 여전히 멀면 → d-감쇠 실재 증거
            log(f"    판정: 중립/양성 방향")
            log(f"      ρ(N_MAX=300)={rho_g300:.3f} > -0.75")
            log(f"      GL(2)는 N_MAX=300에서도 GUE({GUE_FULL_REF})보다 여전히 약함")
            delta_arith_n300 = abs(GUE_FULL_REF) - abs(rho_g300)
            log(f"      δ_arith(N_MAX=300) = {delta_arith_n300:+.4f}")
            if dg > 0.15:
                log(f"      N_MAX 효과 존재하나 d=2 감쇠 방향 유지")
            else:
                log(f"      N_MAX 효과 미미 → d-감쇠는 N_MAX와 무관")
        elif rho_g300 < -0.85:
            # ρ=-0.90처럼 GUE에 근접 → 아티팩트 가능성
            log(f"    판정: 음성 (아티팩트 위험)")
            log(f"      ρ(N_MAX=300)={rho_g300:.3f} < -0.85")
            log(f"      GL(2) N_MAX=300이 ζ(s) 수준에 근접 → d-감쇠가 N_MAX 아티팩트일 수 있음")
            log(f"      → Paper 4 통일표 방법론 재검토 필요!")
        else:
            # ρ ∈ [-0.85, -0.75]
            log(f"    판정: 중립 (부분 아티팩트)")
            log(f"      ρ(N_MAX=300)={rho_g300:.3f} ∈ [-0.85, -0.75]")
            log(f"      N_MAX 효과와 d-감쇠 혼재. 추가 분석 필요.")

        # d-감쇠 보존 여부
        log()
        log(f"  d-감쇠 비교 (같은 N_MAX=300 기준):")
        if 'zeta' in all_results and 300 in all_results['zeta']:
            rho_z300 = all_results['zeta'][300].get('rho_S', float('nan'))
            if not math.isnan(rho_z300) and not math.isnan(rho_g300):
                d_gap = abs(rho_z300) - abs(rho_g300)
                log(f"    |ρ_ζ(N300)| - |ρ_GL2(N300)| = {abs(rho_z300):.4f} - {abs(rho_g300):.4f} = {d_gap:+.4f}")
                if d_gap > 0.05:
                    log(f"    ζ(s)가 GL(2)보다 N_MAX=300에서도 여전히 강함 → d-감쇠 방향 유지")
                else:
                    log(f"    N_MAX=300에서 ζ(s)와 GL(2) 거의 동일 → d-감쇠는 N_MAX 아티팩트")
else:
    log("    결과 없음")

# ── 5. 논문 파급 ──────────────────────────────────────────────────────────
log()
log("=" * 70)
log("  [Paper 4 파급]")
log("=" * 70)
if 'gl2_11a1' in all_results and 300 in all_results['gl2_11a1']:
    rho_g300 = all_results['gl2_11a1'][300].get('rho_S', float('nan'))
    if not math.isnan(rho_g300):
        if rho_g300 > -0.75:
            log("  → 현재 Paper 4 통일표 유효. N_MAX=300 재측정 결과를 주석으로 추가.")
        elif rho_g300 < -0.85:
            log("  → 통일표의 d=1 vs d=2 비교가 N_MAX 아티팩트. 긴급 수정 필요.")
            log("    모든 L-함수에 대해 동일한 N_MAX로 재측정 필요.")
        else:
            log("  → 부분 아티팩트. N_MAX을 명시하고 감도 분석 결과를 포함할 것.")

elapsed = time.time() - t_start
log()
log(f"  총 소요: {elapsed:.1f}s ({elapsed/60:.1f}분)")
log()

# 저장
save()

log("[C-291 완료]")
