#!/usr/bin/env python3
"""
=============================================================================
[사이클 #270] gap_min t-안정성 + Cauchy 저-t 재현 — B-46/B-45 판정
=============================================================================

목적:
  C-268에서 발견된 gap_right 상관의 t-의존성(B-46)이
  gap_min에서도 발생하는지 확인.
  Paper 4 Theorem 1 수치적 유효범위 결정.

측정 항목:
  PART A: T=2000 데이터, 8 t-bin별 gap_min/gap_right 상관
    - ρ(A_L,  gap_min_GUE) per t-bin → t-안정성 확인
    - ρ(A_Λ,  gap_min_GUE) per t-bin → 동일
    - ρ(A_L,  gap_right_GUE) per t-bin → C-268과 재확인
    - ρ(A_Λ,  gap_right_GUE) per t-bin → 동일
    - 교차표: t-bin × gap_metric × method 2×2×8

  PART B: T=300 Cauchy 저-t 재현
    - mpmath.zetazero + xi 코시 윤곽으로 A_Λ 직접 계산
    - ρ(A_Λ, gap_right_GUE) → C-256(-0.578) vs C-268 bin1(-0.686) 비교
    - B-45 해소 가능성 탐색

핵심 가설:
  B-46: gap_right ρ는 t와 함께 감쇠. gap_min ρ는 t-안정적.
  B-45: Cauchy 저-t ρ ≈ C-256(-0.578) → 방법론 동일 조건 재현

성공 기준:
  ★★★★ gap_min ρ가 8 bin에서 ±0.10 이내 안정
  ★★★  gap_min ρ도 t-감쇠 (음성이나 가치)
  ★★★  Cauchy low-t ρ ≈ −0.55~−0.65
  ★★   Cauchy low-t ρ ≠ C-256 (방법론 재검토)

결과: results/gap_stability_c270.txt
=============================================================================
"""

import sys, os, time
import numpy as np
from scipy import stats
import mpmath

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import cypari2
from bundle_utils import xi_func

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CENTER     = 0.5
T_MAX_A    = 2000      # Part A: T=2000 (C-268과 동일)
T_MAX_B    = 300       # Part B: T=300 (저-t 재현)
TRIM_FRAC  = 0.20      # 양쪽 20% 제외 → 중앙 60%
N_BINS     = 8         # t-bin 수

# Cauchy 파라미터 (C-256과 동일)
CONTOUR_RADIUS = mpmath.mpf('0.01')
N_CONTOUR_PTS  = 128
NU_DIFF_H      = mpmath.mpf(1) / mpmath.mpf(10**18)

# 참조값
RHO_GUE        = -0.50    # GUE 이론
RHO_C256       = -0.5898  # C-256 Cauchy ρ(A, gap_right_GUE)
RHO_C268_FULL  = -0.4553  # C-268 Hadamard ρ(A_Λ, gap_right_GUE) 전체
RHO_C268_BIN1  = -0.686   # C-268 bin 1 (t̄=626) ρ(A_Λ, gap_right_GUE)

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gap_stability_c270.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 로그
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

_log_buf = []

def log(msg=''):
    print(msg, flush=True)
    _log_buf.append(str(msg))

def save():
    with open(RESULT_PATH, 'w') as f:
        f.write('\n'.join(_log_buf))
    log(f"[저장] {RESULT_PATH}")

sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Cauchy 적분 — C-256 방법 (xi_func 기반)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_A_cauchy(t_zero):
    """
    Cauchy 적분으로 c₀, c₁ 계산 (C-256 방법 재사용).

    ξ'/ξ(ρ + u) = 1/u + c₀ + c₁u + ...
    g(u) = ξ'/ξ(ρ+u) - 1/u  (정칙 부분)
    c₀ = (1/N) Σ g(u_k),  c₁ = (1/N) Σ g(u_k)/u_k

    A = Im(c₀)² + 2Re(c₁)  (Cor 4.2)
    반환: (A, c0, c1) — float, complex, complex
    """
    rho = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_zero)))
    r   = CONTOUR_RADIUS
    N   = N_CONTOUR_PTS
    h   = NU_DIFF_H

    c0_sum = mpmath.mpc(0, 0)
    c1_sum = mpmath.mpc(0, 0)

    for k in range(N):
        theta = 2 * mpmath.pi * k / N
        u     = r * mpmath.exp(mpmath.mpc(0, theta))
        s     = rho + u

        xi_val = xi_func(s)
        xi_p   = xi_func(s + h)
        xi_m   = xi_func(s - h)
        L_val  = (xi_p - xi_m) / (2 * h * xi_val)

        g = L_val - 1 / u

        c0_sum += g
        c1_sum += g / u

    c0 = c0_sum / N
    c1 = c1_sum / N

    c0_f = complex(c0)
    c1_f = complex(c1)

    A = c0_f.imag**2 + 2 * c1_f.real
    return A, c0_f, c1_f


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    log("=" * 80)
    log("[C-270] gap_min t-안정성 + Cauchy 저-t 재현")
    log("=" * 80)
    log(f"  PART A: T={T_MAX_A}, trim={TRIM_FRAC*2*100:.0f}% 제외, N_BINS={N_BINS}")
    log(f"  PART B: T={T_MAX_B}, Cauchy 윤곽 방법 (C-256 재현)")
    log(f"  참조: ρ_GUE={RHO_GUE}, ρ_C256={RHO_C256}, ρ_C268_full={RHO_C268_FULL}")
    log()

    # ══════════════════════════════════════════════════════════════════════
    # PART A: T=2000 zero-sum + Hadamard, t-bin별 gap_min/gap_right 분석
    # ══════════════════════════════════════════════════════════════════════

    log("=" * 80)
    log("PART A: T=2000 zero-sum + Hadamard — t-bin 교차표")
    log("=" * 80)
    log()

    # ──────────────────────────────────────────────────────────────────────
    # A-STEP 1: ζ(s) 영점 수집 (PARI, T=2000)
    # ──────────────────────────────────────────────────────────────────────
    log("[A-STEP 1] ζ(s) 영점 수집 (PARI, T=2000)...")
    pari = cypari2.Pari()
    pari.allocatemem(256 * 10**6)
    pari.set_real_precision(38)
    mpmath.mp.dps = 30

    pari('Li_z = lfuninit(lfuncreate(1), [0, 2100])')
    pari('zv = lfunzeros(Li_z, 2000)')
    n_z = int(str(pari('#zv')))

    zeros = []
    for i in range(1, n_z + 1):
        t = float(str(pari(f'zv[{i}]')).replace(' E', 'e'))
        if t > 5:
            zeros.append(t)
    zeros = sorted(zeros)

    log(f"  {len(zeros)}개 영점, t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    log(f"  소요: {time.time()-t_start:.1f}초")
    log()

    # ──────────────────────────────────────────────────────────────────────
    # A-STEP 2: A_L / A_Λ 병렬 계산 (C-268과 동일 방법)
    # ──────────────────────────────────────────────────────────────────────
    log("[A-STEP 2] A_L / A_Λ 계산 (zero-sum + Hadamard)...")
    log(f"  K = 전체({len(zeros)}개)")
    log()

    data = []
    n_fail = 0

    for i in range(len(zeros)):
        g0 = zeros[i]

        try:
            # 같은 부호 영점 합산
            S1_same = 0.0
            H1_same = 0.0
            for k in range(len(zeros)):
                if k == i:
                    continue
                diff = g0 - zeros[k]
                S1_same += 1.0 / diff
                H1_same += 1.0 / diff**2

            # 켤레 영점 합산
            S1_conj = sum(1.0 / (g0 + zeros[k]) for k in range(len(zeros)))
            H1_conj = sum(1.0 / (g0 + zeros[k])**2 for k in range(len(zeros)))

            # A_L: pure zero-sum (Gamma 미포함)
            S1_L = S1_same + S1_conj
            H1_L = H1_same + H1_conj
            A_L  = S1_L**2 + 2.0 * H1_L

            # Gamma 보정
            s_val      = mpmath.mpc(CENTER, g0)
            gamma_S    = -mpmath.log(mpmath.pi) / 2 + mpmath.digamma(s_val / 2) / 2
            im_gamma_S = float(mpmath.im(gamma_S))
            psi1_val   = mpmath.psi(1, s_val / 2)
            re_gamma_H = float(mpmath.re(psi1_val)) / 4.0

            # A_Λ: completed (Gamma 포함)
            S1_Lambda = S1_L - im_gamma_S
            H1_Lambda = H1_L + re_gamma_H
            A_Lambda  = S1_Lambda**2 + 2.0 * H1_Lambda

            if A_L <= 0 or A_Lambda <= 0:
                n_fail += 1
                continue

            data.append({
                't':        g0,
                'A_L':      A_L,
                'A_Lambda': A_Lambda,
            })

        except Exception as e:
            n_fail += 1
            if n_fail <= 5:
                print(f"WARNING: i={i}, t={g0:.3f}: {e}", flush=True)

        if (i + 1) % 300 == 0:
            log(f"  {i+1}/{len(zeros)} 완료... ({time.time()-t_start:.0f}s, "
                f"유효={len(data)}, 실패={n_fail})")

    log(f"  계산 완료: 유효={len(data)}, 실패={n_fail}, "
        f"소요={time.time()-t_start:.1f}초")

    if len(data) < 100:
        log("⚠️ 유효 데이터 100개 미만 — 중단")
        save()
        return

    if n_fail > len(zeros) // 2:
        log(f"⚠️ 실패율 {n_fail}/{len(zeros)} — 절반 초과, 신뢰도 낮음")

    log()

    # ──────────────────────────────────────────────────────────────────────
    # A-STEP 3: 중앙 60% 선택 + gap 계산
    # ──────────────────────────────────────────────────────────────────────
    log("[A-STEP 3] 중앙 60% 선택 + GUE 정규화 gap 계산...")

    N = len(data)
    lo = int(N * TRIM_FRAC)
    hi = int(N * (1.0 - TRIM_FRAC))

    valid = []
    for idx in range(lo, hi):
        if idx <= 0 or idx >= N - 1:
            continue
        d = data[idx]
        gap_r = data[idx+1]['t'] - d['t']
        gap_l = d['t'] - data[idx-1]['t']
        if gap_r <= 0 or gap_l <= 0:
            continue
        d_bar = np.log(d['t'] / (2.0 * np.pi)) / (2.0 * np.pi)
        d['gap_r_gue']   = gap_r * d_bar
        d['gap_l_gue']   = gap_l * d_bar
        d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
        d['d_bar']       = d_bar
        valid.append(d)

    log(f"  내부 영점: {len(valid)} (중앙 {(1.0-2*TRIM_FRAC)*100:.0f}%)")
    log()

    if len(valid) < 50:
        log("⚠️ 내부 영점 50개 미만 — 분석 중단")
        save()
        return

    # 배열 추출
    n_v    = len(valid)
    t_arr  = np.array([d['t']          for d in valid])
    AL_arr = np.array([d['A_L']        for d in valid])
    ALa_arr= np.array([d['A_Lambda']   for d in valid])
    gm_arr = np.array([d['gap_min_gue'] for d in valid])
    gr_arr = np.array([d['gap_r_gue']   for d in valid])

    # ──────────────────────────────────────────────────────────────────────
    # A-STEP 4: 전체 상관 (4가지)
    # ──────────────────────────────────────────────────────────────────────
    log("[A-STEP 4] 전체 상관 분석 (n={})".format(n_v))
    log()
    log(f"  {'항목':<38} {'ρ':>8}  {'p':>10}  {'판정'}")
    log("  " + "-" * 68)

    corr_global = {}
    for name, x, y in [
        ("A_L,   gap_min_GUE",   AL_arr,  gm_arr),
        ("A_L,   gap_right_GUE", AL_arr,  gr_arr),
        ("A_Λ,   gap_min_GUE",   ALa_arr, gm_arr),
        ("A_Λ,   gap_right_GUE", ALa_arr, gr_arr),
    ]:
        try:
            r, p = stats.spearmanr(x, y)
        except Exception as e:
            log(f"  WARNING: {name} — {e}")
            continue
        corr_global[name] = (r, p)
        log(f"  ρ({name:<36}) = {r:+.4f}  (p={p:.3e})  {sig(p)}")

    log()

    # ──────────────────────────────────────────────────────────────────────
    # A-STEP 5: t-bin × gap_metric 교차표 [핵심]
    # ──────────────────────────────────────────────────────────────────────
    log("[A-STEP 5] t-bin × gap_metric 교차표 (핵심: B-46 판정)")
    log()
    log(f"  빈당 예상 n ≈ {n_v // N_BINS}")
    log()

    t_quantiles = np.percentile(t_arr, np.linspace(0, 100, N_BINS + 1))

    # 헤더 출력
    log(f"  {'bin':>4}  {'t̄':>8}  "
        f"{'ρ(A_L,gm)':>11}  {'ρ(A_Λ,gm)':>11}  "
        f"{'ρ(A_L,gr)':>11}  {'ρ(A_Λ,gr)':>11}  {'n':>5}")
    log("  " + "-" * 80)

    bin_results = []
    for b in range(N_BINS):
        if b < N_BINS - 1:
            mask = (t_arr >= t_quantiles[b]) & (t_arr < t_quantiles[b+1])
        else:
            mask = (t_arr >= t_quantiles[b]) & (t_arr <= t_quantiles[b+1])

        n_bin = int(np.sum(mask))
        if n_bin < 5:
            log(f"  {b+1:>4}  빈 n={n_bin} < 5 — 건너뜀")
            continue

        t_mean = float(np.mean(t_arr[mask]))

        def corr_bin(x, y, n_min=5):
            if np.sum(mask) < n_min:
                return np.nan, np.nan
            try:
                return stats.spearmanr(x[mask], y[mask])
            except Exception as e:
                print(f"WARNING corr_bin: {e}", flush=True)
                return np.nan, np.nan

        r_AL_gm, p_AL_gm   = corr_bin(AL_arr, gm_arr)
        r_ALa_gm, p_ALa_gm = corr_bin(ALa_arr, gm_arr)
        r_AL_gr, p_AL_gr   = corr_bin(AL_arr, gr_arr)
        r_ALa_gr, p_ALa_gr = corr_bin(ALa_arr, gr_arr)

        bin_results.append({
            'b': b+1, 't_mean': t_mean, 'n': n_bin,
            'r_AL_gm': r_AL_gm, 'p_AL_gm': p_AL_gm,
            'r_ALa_gm': r_ALa_gm, 'p_ALa_gm': p_ALa_gm,
            'r_AL_gr': r_AL_gr, 'p_AL_gr': p_AL_gr,
            'r_ALa_gr': r_ALa_gr, 'p_ALa_gr': p_ALa_gr,
        })

        log(f"  {b+1:>4}  {t_mean:>8.1f}  "
            f"{r_AL_gm:>+11.4f}  {r_ALa_gm:>+11.4f}  "
            f"{r_AL_gr:>+11.4f}  {r_ALa_gr:>+11.4f}  {n_bin:>5}")

    log()

    # ── Kendall τ로 단조성 검증 ──
    if len(bin_results) >= 4:
        t_seq     = [x['t_mean']  for x in bin_results]
        seq_AL_gm = [x['r_AL_gm']  for x in bin_results]
        seq_ALa_gm= [x['r_ALa_gm'] for x in bin_results]
        seq_AL_gr = [x['r_AL_gr']  for x in bin_results]
        seq_ALa_gr= [x['r_ALa_gr'] for x in bin_results]

        def kendall_test(seq_name, seq):
            vals = [v for v in seq if not np.isnan(v)]
            ts   = [t_seq[i] for i, v in enumerate(seq) if not np.isnan(v)]
            if len(vals) < 4:
                log(f"  Kendall τ({seq_name}): N/A")
                return np.nan, np.nan
            tau, p = stats.kendalltau(ts, vals)
            log(f"  Kendall τ({seq_name:<30}) = {tau:+.4f}  (p={p:.3e})  {sig(p)}")
            return tau, p

        log("[A-STEP 5.2] t-bin 단조성 — Kendall τ")
        log()
        tau_AL_gm, p_AL_gm_kt   = kendall_test("A_L,   gap_min",  seq_AL_gm)
        tau_ALa_gm, p_ALa_gm_kt = kendall_test("A_Λ,   gap_min",  seq_ALa_gm)
        tau_AL_gr,  p_AL_gr_kt  = kendall_test("A_L,   gap_right", seq_AL_gr)
        tau_ALa_gr, p_ALa_gr_kt = kendall_test("A_Λ,   gap_right", seq_ALa_gr)
        log()

        # gap_min 안정성 판정
        gm_range_AL  = max(seq_AL_gm)  - min(seq_AL_gm)  if seq_AL_gm else np.nan
        gm_range_ALa = max(seq_ALa_gm) - min(seq_ALa_gm) if seq_ALa_gm else np.nan

        log("[A-STEP 5.3] gap_min 안정성 판정")
        log()
        log(f"  ρ(A_L,  gap_min) 범위: [{min(seq_AL_gm):+.4f}, {max(seq_AL_gm):+.4f}]  "
            f"span={gm_range_AL:.4f}")
        log(f"  ρ(A_Λ,  gap_min) 범위: [{min(seq_ALa_gm):+.4f}, {max(seq_ALa_gm):+.4f}]  "
            f"span={gm_range_ALa:.4f}")
        log()

        # ★★★★ 기준: ±0.10 이내 (span < 0.20)
        stable_AL_gm  = gm_range_AL  < 0.20
        stable_ALa_gm = gm_range_ALa < 0.20
        decay_ALa_gr  = abs(tau_ALa_gr) > 0.5 if not np.isnan(tau_ALa_gr) else False

        log(f"  [B-46 판정] gap_min 안정 (span < 0.20):")
        log(f"    A_L  gap_min: {'✅ 안정' if stable_AL_gm  else '⚠️ 불안정'}  "
            f"(span={gm_range_AL:.4f})")
        log(f"    A_Λ  gap_min: {'✅ 안정' if stable_ALa_gm else '⚠️ 불안정'}  "
            f"(span={gm_range_ALa:.4f})")
        log(f"    A_Λ  gap_right: {'✅ 감쇠 확인' if decay_ALa_gr else '❌ 감쇠 없음'}  "
            f"(τ={tau_ALa_gr:+.4f})")
        log()

        if stable_AL_gm and stable_ALa_gm:
            log("  ★★★★ B-46 양성: gap_min은 t-안정적 → Paper 4 Thm 1 수치 유효범위 확정")
            log("        gap_right만 t-의존적 → gap_min이 정준(canonical) 메트릭")
        elif not stable_ALa_gm:
            log("  ★★★  B-46 음성: gap_min도 t-감쇠 → 비보편성 발견 (더 큰 논문)")
        else:
            log("  ★★★  B-46 부분: A_L 안정 / A_Λ 불안정 → Gamma 오염 심각")
        log()

    else:
        log("⚠️ 빈 수 < 4 — Kendall τ 분석 불가")
        tau_AL_gm = tau_ALa_gm = tau_AL_gr = tau_ALa_gr = np.nan

    save()

    # ══════════════════════════════════════════════════════════════════════
    # PART B: T=300 Cauchy 저-t 재현 (C-256 조건 복원)
    # ══════════════════════════════════════════════════════════════════════

    log()
    log("=" * 80)
    log("PART B: Cauchy T=300 저-t 재현 — B-45 해소 시도")
    log("=" * 80)
    log()

    # ──────────────────────────────────────────────────────────────────────
    # B-STEP 1: T=300 이하 영점 수집 (mpmath.zetazero)
    # ──────────────────────────────────────────────────────────────────────
    log("[B-STEP 1] T=300 이하 영점 수집 (mpmath.zetazero)...")
    mpmath.mp.dps = 60

    zeros_b = []
    n_b = 1
    while True:
        try:
            t = float(mpmath.zetazero(n_b).imag)
        except Exception as e:
            log(f"  WARNING: zetazero({n_b}): {e}")
            break
        if t > T_MAX_B:
            break
        zeros_b.append(t)
        n_b += 1

    zeros_b = np.array(zeros_b)
    n_raw = len(zeros_b)
    log(f"  {n_raw}개 영점, t ∈ [{zeros_b[0]:.2f}, {zeros_b[-1]:.2f}]")
    log(f"  소요: {time.time()-t_start:.1f}초")
    log()

    if n_raw < 20:
        log("⚠️ 영점 20개 미만 — Part B 중단")
        save()
        return

    # ──────────────────────────────────────────────────────────────────────
    # B-STEP 2: Cauchy 코시 윤곽으로 A_Λ 계산
    # ──────────────────────────────────────────────────────────────────────
    log("[B-STEP 2] Cauchy 윤곽으로 A_Λ 계산 (C-256 방법)...")
    log(f"  반경={CONTOUR_RADIUS}, N_pts={N_CONTOUR_PTS}, dps={mpmath.mp.dps}")
    log()

    A_b_list  = []
    c0_b_list = []
    t_b_list  = []
    n_fail_b  = 0

    for j in range(n_raw):
        t_n = zeros_b[j]
        if (j + 1) % 20 == 0 or j == 0:
            log(f"  영점 #{j+1}/{n_raw}  t={t_n:.2f}  경과={time.time()-t_start:.0f}s")

        try:
            A, c0, c1 = compute_A_cauchy(t_n)
            if np.isnan(A) or np.isinf(A):
                n_fail_b += 1
                if n_fail_b <= 5:
                    log(f"  ⚠️ 영점 #{j+1} A=NaN/Inf — 건너뜀")
                continue
            A_b_list.append(A)
            c0_b_list.append(c0)
            t_b_list.append(t_n)
        except Exception as e:
            n_fail_b += 1
            if n_fail_b <= 5:
                log(f"  WARNING: 영점 #{j+1} t={t_n:.2f}: {e}")

    log()
    log(f"  계산 완료: 유효={len(A_b_list)}/{n_raw}, 실패={n_fail_b}, "
        f"소요={time.time()-t_start:.1f}초")
    log()

    if len(A_b_list) < 20:
        log("⚠️ 유효 데이터 20개 미만 — Part B 분석 중단")
        save()
        return

    if n_fail_b > n_raw // 2:
        log(f"⚠️ 실패율 {n_fail_b}/{n_raw} — 절반 초과, 신뢰도 낮음")

    # ──────────────────────────────────────────────────────────────────────
    # B-STEP 3: gap 계산 + 상관 분석
    # ──────────────────────────────────────────────────────────────────────
    log("[B-STEP 3] gap 계산 + 상관 분석...")

    A_b  = np.array(A_b_list)
    t_b  = np.array(t_b_list)
    n_b_valid = len(A_b)

    # 유효 인덱스에서 내부 영점만 (edge 제외)
    inner_b = []
    for k in range(len(t_b)):
        # t_b[k]가 t_b 배열에서 위치 찾기 (zeros_b 기준 인접)
        pos = np.searchsorted(zeros_b, t_b[k])
        if pos <= 0 or pos >= len(zeros_b) - 1:
            continue
        gap_r = zeros_b[pos + 1] - zeros_b[pos]
        gap_l = zeros_b[pos] - zeros_b[pos - 1]
        if gap_r <= 0 or gap_l <= 0:
            continue
        d_bar = np.log(t_b[k] / (2.0 * np.pi)) / (2.0 * np.pi)
        inner_b.append({
            't': t_b[k],
            'A': A_b[k],
            'gap_r_gue':   gap_r * d_bar,
            'gap_l_gue':   gap_l * d_bar,
            'gap_min_gue': min(gap_r, gap_l) * d_bar,
        })

    log(f"  내부 영점: {len(inner_b)} (edge 제외)")
    log()

    if len(inner_b) < 10:
        log("⚠️ 내부 영점 10개 미만 — 상관 분석 중단")
        save()
        return

    A_ib  = np.array([d['A']          for d in inner_b])
    gr_ib = np.array([d['gap_r_gue']  for d in inner_b])
    gm_ib = np.array([d['gap_min_gue'] for d in inner_b])
    t_ib  = np.array([d['t']          for d in inner_b])

    n_ib = len(A_ib)

    log(f"  n = {n_ib}, t ∈ [{t_ib.min():.1f}, {t_ib.max():.1f}]")
    log()
    log(f"  {'항목':<38} {'ρ':>8}  {'p':>10}  {'판정'}")
    log("  " + "-" * 65)

    for name, x, y in [
        ("A_Λ(Cauchy), gap_right_GUE", A_ib, gr_ib),
        ("A_Λ(Cauchy), gap_min_GUE",   A_ib, gm_ib),
    ]:
        try:
            r, p = stats.spearmanr(x, y)
            log(f"  ρ({name:<36}) = {r:+.4f}  (p={p:.3e})  {sig(p)}")
        except Exception as e:
            log(f"  WARNING: {name} — {e}")

    log()

    # ──────────────────────────────────────────────────────────────────────
    # B-STEP 4: C-256 / C-268 비교
    # ──────────────────────────────────────────────────────────────────────
    log("[B-STEP 4] Cauchy 방법 간 비교 — B-45 판정")
    log()

    try:
        rho_b_gr, p_b_gr = stats.spearmanr(A_ib, gr_ib)
        rho_b_gm, p_b_gm = stats.spearmanr(A_ib, gm_ib)
    except Exception as e:
        log(f"  WARNING: {e}")
        rho_b_gr = rho_b_gm = np.nan
        p_b_gr = p_b_gm = np.nan

    log(f"  C-256 Cauchy (T~396, n=198):    ρ(A, gap_right_GUE) = {RHO_C256:+.4f}")
    log(f"  C-270 Cauchy (T≤300, n={n_ib}):  ρ(A, gap_right_GUE) = {rho_b_gr:+.4f}  "
        f"(p={p_b_gr:.3e})  {sig(p_b_gr) if not np.isnan(rho_b_gr) else ''}")
    log(f"  C-268 bin-1  (t̄=626, Hadamard): ρ(A_Λ, gap_right_GUE) = {RHO_C268_BIN1:+.4f}")
    log()

    if not np.isnan(rho_b_gr):
        diff_vs_c256 = abs(rho_b_gr - RHO_C256)
        diff_vs_c268b1 = abs(rho_b_gr - RHO_C268_BIN1)
        log(f"  |ρ_C270 - ρ_C256|    = {diff_vs_c256:.4f}")
        log(f"  |ρ_C270 - ρ_C268b1|  = {diff_vs_c268b1:.4f}")
        log()

        if diff_vs_c256 < 0.10:
            log("  ★★★ [B-45 해소 기반 확보] C-270 ≈ C-256")
            log("       방법 동일 조건 재현 성공 → 차이는 t-범위 차이로 설명됨")
        elif diff_vs_c256 < 0.15:
            log("  ★★  [B-45 부분 해소] C-270 ≈ C-256 (10~15% 차이)")
        else:
            log(f"  ⚠️  [B-45 미해소] 차이 {diff_vs_c256:.4f} — 방법론 재검토 필요")

        log()

        # t-의존성 설명 가능 여부
        log("  [t-의존성으로 B-45 설명 가능 여부]")
        # C-256: t up to ~400, avg t ≈ 200
        # C-268 full: t up to 2000, avg t ≈ 1000
        # C-268 bin1: t̄ ≈ 626
        # C-270: t up to 300, avg t ≈ 150-200
        log(f"  C-270 t̄ ≈ {t_ib.mean():.1f}  (가장 낮은 t)")
        log(f"  C-256 t̄ ≈ 200 (N_ZEROS=200, 저-t)")
        log(f"  C-268 bin1 t̄ ≈ 626")
        log(f"  C-268 전체 t̄ ≈ 1000")
        log()
        log("  가설: ρ(t낮은) ≈ −0.65~−0.70, ρ(t높은) ≈ −0.12")
        log(f"  C-270 관측: {rho_b_gr:+.4f}")
        if rho_b_gr < -0.50:
            log("  → 저-t에서 강한 음상관 확인 ✅ t-의존성 가설 지지")
        else:
            log("  → 저-t에서도 약한 상관 — 추가 분석 필요")

    log()
    save()

    # ══════════════════════════════════════════════════════════════════════
    # 최종 요약
    # ══════════════════════════════════════════════════════════════════════

    log()
    log("=" * 80)
    log("[C-270 최종 요약] — B-46 / B-45 판정")
    log("=" * 80)
    log()

    # gap_min 안정성 요약
    if len(bin_results) >= 4:
        seq_AL_gm_v  = [x['r_AL_gm']  for x in bin_results if not np.isnan(x['r_AL_gm'])]
        seq_ALa_gm_v = [x['r_ALa_gm'] for x in bin_results if not np.isnan(x['r_ALa_gm'])]
        seq_ALa_gr_v = [x['r_ALa_gr'] for x in bin_results if not np.isnan(x['r_ALa_gr'])]

        log("  [A. gap_min t-안정성 — B-46 핵심]")
        log()
        log(f"  t-bin별 ρ(A_Λ, gap_min_GUE) 요약:")
        for br in bin_results:
            stable_flag = '≈' if not np.isnan(br['r_ALa_gm']) else '?'
            log(f"    bin {br['b']} (t̄={br['t_mean']:6.1f}): "
                f"ρ_gm={br['r_ALa_gm']:+.4f}  ρ_gr={br['r_ALa_gr']:+.4f}  "
                f"n={br['n']}")
        log()

        gm_range = max(seq_ALa_gm_v) - min(seq_ALa_gm_v) if seq_ALa_gm_v else np.nan
        gr_range = max(seq_ALa_gr_v) - min(seq_ALa_gr_v) if seq_ALa_gr_v else np.nan
        log(f"  gap_min span: {gm_range:.4f}  (안정 기준 < 0.20)")
        log(f"  gap_right span: {gr_range:.4f}  (C-268에서 확인된 감쇠)")
        log()

        if gm_range < 0.20:
            log("  B-46 판정: ★★★★ gap_min은 t-안정적")
            log("    → Paper 4 Thm 1 (A ≥ 4/gap_min²)의 수치 유효범위 전 t에서 확정")
            log("    → gap_right만 Gamma-오염 / 높이 의존적")
            log("    → gap_min이 정준(canonical) 메트릭으로 확정")
        elif gm_range < 0.40:
            log("  B-46 판정: ★★★ gap_min도 중등 t-의존성")
            log("    → Paper 4 서사에서 t-범위 제한 명시 필요")
        else:
            log("  B-46 판정: ★★★ gap_min도 강한 t-의존성 (비보편성 발견)")
            log("    → Paper 4 서사 근본 수정 필요")
    else:
        log("  B-46: 빈 수 부족으로 판정 불가")
        gm_range = np.nan

    log()
    log("  [B. Cauchy 저-t 재현 — B-45]")
    log()
    log(f"  C-256 (Cauchy, t≤396, n=198):  ρ = {RHO_C256:+.4f}")
    if not np.isnan(rho_b_gr):
        log(f"  C-270 (Cauchy, t≤300, n={n_ib}):  ρ = {rho_b_gr:+.4f}")
        diff = abs(rho_b_gr - RHO_C256)
        log(f"  차이: {diff:.4f}")
        if diff < 0.10:
            log("  B-45 해소 기반: ✅ Cauchy 방법 재현 성공")
            log("         B-45 차이(-0.12)는 Hadamard vs Cauchy 방법론 차이에 기인")
        else:
            log(f"  B-45: 차이 {diff:.4f} — 방법론 외 요인 존재 가능성")
    else:
        log("  C-270 Cauchy: 계산 실패")

    log()
    log(f"  총 소요시간: {time.time()-t_start:.1f}초")
    log()

    save()
    log(f"[완료] 결과: {RESULT_PATH}")


if __name__ == '__main__':
    main()
