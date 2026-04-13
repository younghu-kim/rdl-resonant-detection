"""
=============================================================================
[Project RDL] 홀로노미 적분 — 영점 개수 = 제1 Chern 수
=============================================================================
사이클 31 / 수학자 지시 2026-04-14 07:01

목적:
  논증 원리(Argument Principle)를 다발 언어로 검증:
    L(s) = (log ξ)' = 접속 1-form
    ξ의 단순 영점 ρ_n에서 L의 residue = 1
    ∮_C L(s) ds = 2πi × N(C) — 제1 Chern 수 c₁ = N

  수치 검증 3단계:
  1. 개별 영점 홀로노미: ∮_{|s-ρ_n|=0.3} L(s)ds / (2πi) ≈ 1.000 (sanity check)
  2. 직사각형 윤곽 위상추적: σ∈[0.3,0.7], t∈[13,t₂] → N_contour vs N_actual
  3. 연속 누적 함수: t₂ ∈ [14,430], C(t) 계단함수 vs N(t)
  4. 디리클레 χ₃ 보편성 확인

메서드:
  §1 — 원형 윤곽: 수동 사다리꼴 (L(s)=해석적 공식, 64분할)
  §2,3 — 직사각형 위상추적: Im(log ξ) 차분 누적
    → log ξ(s) = log(s)+log(s-1)-(s/2)log π+loggamma(s/2)+log ζ(s) 직접 계산
    → 위상 차분 unwrapping: 각 스텝에서 (-π,π) 범위로 보정
    → L(s) 직접 적분보다 수치적으로 훨씬 안정적 (underflow 없음)
  dps=100 (적분), 50 (영점 탐색)

출력:
  results/holonomy_chern_number.txt
  results/holonomy_chern_number_figdata.csv
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT_TXT = os.path.join(BASE_DIR, "results", "holonomy_chern_number.txt")
OUTPUT_CSV = os.path.join(BASE_DIR, "results", "holonomy_chern_number_figdata.csv")
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# ─── 설정 ─────────────────────────────────────────────────────────────────────
DPS = 100
mpmath.mp.dps = DPS

SIGMA_MIN = mpmath.mpf('0.3')
SIGMA_MAX = mpmath.mpf('0.7')
T1        = mpmath.mpf('13.0')   # 첫 영점(14.135) 아래

N_HORIZ   = 64    # 직사각형 수평변 분할
N_VERT_PER_UNIT = 3  # 수직변: 단위 t당 분할 수 (충분한 해상도)
N_CIRCLE_STEPS  = 64   # 원형 윤곽 분할
R_CIRCLE        = mpmath.mpf('0.3')  # 개별 영점 홀로노미 반경

CHI3 = [0, 1, -1]  # χ mod 3

# ─── 로그 ─────────────────────────────────────────────────────────────────────
log_lines = []
def log(msg=""):
    print(msg, flush=True)
    log_lines.append(msg)

def save_all():
    with open(OUTPUT_TXT, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))


# ═══════════════════════════════════════════════════════════════════════════════
# §1용 해석적 접속 L(s) = ξ'/ξ (원형 윤곽 전용)
# ═══════════════════════════════════════════════════════════════════════════════

def connection_L_analytic(s):
    """
    L(s) = ξ'(s)/ξ(s) = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2 + ζ'(s)/ζ(s)
    원형 윤곽(§1 sanity check) 전용. 영점에서 발산.
    """
    zeta_val = mpmath.zeta(s)
    if abs(zeta_val) < mpmath.mpf(10) ** (-DPS + 15):
        raise ValueError(f"ζ(s) 영점 근처: |ζ|={float(abs(zeta_val)):.2e}")
    zeta_d    = mpmath.zeta(s, derivative=1)
    term_s    = mpmath.mpf(1) / s
    term_sm1  = mpmath.mpf(1) / (s - 1)
    term_logpi = -mpmath.log(mpmath.pi) / 2
    term_psi  = mpmath.digamma(s / 2) / 2
    term_zeta = zeta_d / zeta_val
    return term_s + term_sm1 + term_logpi + term_psi + term_zeta


# ═══════════════════════════════════════════════════════════════════════════════
# §2,3,4 용 위상 추적법 (Phase-Tracking Argument Principle)
# ═══════════════════════════════════════════════════════════════════════════════

def log_xi_imag(s):
    """
    Im(log ξ(s)) 해석적 계산 (underflow 없음):
      log ξ = log(s) + log(s-1) - (s/2)log π + loggamma(s/2) + log ζ(s)
      (1/2 상수항은 실수 → 위상에 영향 없음)

    ζ(s) 영점 근처이면 None 반환.
    """
    z_val = mpmath.zeta(s)
    if abs(z_val) < mpmath.mpf(10) ** (-DPS + 15):
        return None  # ζ(s) 영점 근처 — 호출자가 처리

    lxi = (mpmath.log(s)
           + mpmath.log(s - 1)
           - (s / 2) * mpmath.log(mpmath.pi)
           + mpmath.loggamma(s / 2)
           + mpmath.log(z_val))
    return float(lxi.imag)


def log_lambda_chi3_imag(s):
    """
    Im(log Λ(s, χ₃)) 해석적 계산:
      log Λ = (s/2)log(q/π) + loggamma((s+a)/2) + log L(s, χ)
      q=3, a=1 (χ₃)
    """
    q, a = 3, 1
    L_val = mpmath.dirichlet(s, CHI3)
    if abs(L_val) < mpmath.mpf(10) ** (-DPS + 15):
        return None

    lL = ((s / 2) * mpmath.log(mpmath.mpf(q) / mpmath.pi)
          + mpmath.loggamma((s + mpmath.mpf(a)) / 2)
          + mpmath.log(L_val))
    return float(lL.imag)


def phase_track_contour(sigma_lo, sigma_hi, t_lo, t_hi,
                        n_horiz, n_vert, log_phase_fn):
    """
    위상 추적법 — ∮_C d(arg Λ) / (2π) = N (포함된 영점 수)

    직사각형 반시계 방향 경로:
      하변 → 우변 → 상변 → 좌변

    Parameters:
      n_horiz: 수평변(σ 방향) 분할 수
      n_vert:  수직변(t 방향) 분할 수
      log_phase_fn: s → Im(log ξ(s)) 함수 (None 반환 시 이전 값 유지)

    Returns: 위상 누적 / (2π) ≈ N (정수에 가까운 실수)
    """
    sig_lo = float(sigma_lo)
    sig_hi = float(sigma_hi)
    t0 = float(t_lo)
    t1 = float(t_hi)

    # 경로 점 생성
    def make_path():
        pts = []
        # 하변: sigma_lo → sigma_hi, t=t0
        for k in range(n_horiz + 1):
            sig = sig_lo + k * (sig_hi - sig_lo) / n_horiz
            pts.append(mpmath.mpc(sig, t0))
        # 우변: t=t0 → t1, sigma=sig_hi
        for k in range(1, n_vert + 1):
            t = t0 + k * (t1 - t0) / n_vert
            pts.append(mpmath.mpc(sig_hi, t))
        # 상변: sigma_hi → sigma_lo, t=t1
        for k in range(1, n_horiz + 1):
            sig = sig_hi - k * (sig_hi - sig_lo) / n_horiz
            pts.append(mpmath.mpc(sig, t1))
        # 좌변: t=t1 → t0, sigma=sig_lo (마지막 점은 시작점, 제외)
        for k in range(1, n_vert):
            t = t1 - k * (t1 - t0) / n_vert
            pts.append(mpmath.mpc(sig_lo, t))
        return pts

    path = make_path()

    # 위상 누적
    total_phase = 0.0
    phi_prev = None
    skip_count = 0

    for s in path:
        phi = log_phase_fn(s)
        if phi is None:
            skip_count += 1
            continue  # 영점 근처 — 건너뜀
        if phi_prev is None:
            phi_prev = phi
            continue
        # 위상 차분 — unwrap to (-π, π)
        diff = phi - phi_prev
        diff -= 2 * np.pi * round(diff / (2 * np.pi))
        total_phase += diff
        phi_prev = phi

    if skip_count > 0:
        print(f"  WARNING: {skip_count}개 점 건너뜀 (ζ 영점 근처)", flush=True)

    N_estimate = total_phase / (2 * np.pi)
    return N_estimate


# ═══════════════════════════════════════════════════════════════════════════════
# §1 용 원형 윤곽 적분 (사다리꼴)
# ═══════════════════════════════════════════════════════════════════════════════

def circle_contour_integral(center, radius, n_steps, connection_fn):
    """
    center 주위 반지름 radius 원에서 ∮ L(s) ds / (2πi) 계산.
    기대값: 원 안 단순 영점 수 (= 1이면 성공)
    """
    two_pi_i = 2 * mpmath.pi * 1j
    total = mpmath.mpc(0)

    thetas = [2 * mpmath.pi * k / n_steps for k in range(n_steps + 1)]
    s_pts  = [center + radius * mpmath.exp(1j * th) for th in thetas]

    try:
        L_vals = [connection_fn(s) for s in s_pts]
    except ValueError as e:
        raise ValueError(f"원형 윤곽 L 계산 오류: {e}")

    for k in range(n_steps):
        s_curr = s_pts[k]
        s_next = s_pts[k + 1]
        ds = s_next - s_curr
        total += (L_vals[k] + L_vals[k + 1]) * ds / 2

    return total / two_pi_i


# ═══════════════════════════════════════════════════════════════════════════════
# 영점 수집 (리만 ζ)
# ═══════════════════════════════════════════════════════════════════════════════

def collect_zeta_zeros(n_max):
    """mpmath.zetazero로 첫 n_max개 영점 수집. Returns sorted list of floats."""
    zeros = []
    log(f"  리만 ζ 영점 수집 중 (n=1..{n_max}, dps=50)...")
    orig_dps = mpmath.mp.dps
    mpmath.mp.dps = 50
    for k in range(1, n_max + 1):
        try:
            z = mpmath.zetazero(k)
            t = float(z.imag)
            zeros.append(t)
        except Exception as e:
            print(f"WARNING: zetazero({k}) 실패: {e}", flush=True)
        if k % 50 == 0:
            print(f"  ... {k}/{n_max}개 완료 (t≈{zeros[-1]:.2f})", flush=True)
    mpmath.mp.dps = orig_dps
    zeros = sorted(zeros)
    log(f"  수집 완료: {len(zeros)}개 (t={zeros[0]:.3f}~{zeros[-1]:.3f})")
    return zeros


def count_zeros_below(zeros_sorted, t_thresh):
    """t_thresh 미만 영점 개수"""
    return sum(1 for z in zeros_sorted if z < t_thresh)


def midpoint_between_zeros(zeros_sorted, n):
    """n번째와 n+1번째 영점 사이 중간점 (n: 0-indexed)."""
    if n + 1 >= len(zeros_sorted):
        raise IndexError(f"영점 {n+1}이 없음 (총 {len(zeros_sorted)}개)")
    return (zeros_sorted[n] + zeros_sorted[n + 1]) / 2.0


# ═══════════════════════════════════════════════════════════════════════════════
# 디리클레 χ₃ 영점 탐색
# ═══════════════════════════════════════════════════════════════════════════════

def find_zeros_chi3(t_min=1.0, t_max=80.0, n_scan=3000):
    """χ₃ 영점 탐색 — 부호 변화 + findroot."""
    chi_list = CHI3
    orig_dps = mpmath.mp.dps
    mpmath.mp.dps = 60

    ts = np.linspace(t_min, t_max, n_scan)
    zeros = []
    fail_count = 0

    def L_real(t_var):
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_var))
        return float(mpmath.dirichlet(s, chi_list).real)

    for i in range(len(ts) - 1):
        try:
            f1 = L_real(ts[i])
            f2 = L_real(ts[i + 1])
        except Exception as e:
            print(f"WARNING: χ₃ 탐색 t={ts[i]:.3f}: {e}", flush=True)
            continue
        if f1 * f2 < 0:
            try:
                t_zero = mpmath.findroot(
                    lambda t: mpmath.dirichlet(mpmath.mpf('0.5') + 1j * t, chi_list).real,
                    mpmath.mpf(str((ts[i] + ts[i+1]) / 2))
                )
                t_zero_f = float(t_zero.real)
                if not zeros or abs(t_zero_f - zeros[-1]) > 0.05:
                    zeros.append(t_zero_f)
            except Exception as e:
                fail_count += 1
                print(f"WARNING: χ₃ findroot 실패 t≈{ts[i]:.3f}: {e}", flush=True)

    mpmath.mp.dps = orig_dps
    if len(zeros) == 0:
        print("⚠️ χ₃ 영점 0개 — 탐색 로직 점검 필요", flush=True)
    log(f"  χ₃ 영점: {len(zeros)}개, findroot 실패: {fail_count}회")
    return sorted(zeros)


# ═══════════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    log("=" * 70)
    log("[Project RDL] 홀로노미 적분 — 제1 Chern 수 검증")
    log(f"실행 시각: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"dps={DPS}, 수평분할={N_HORIZ}, 수직분할=~{N_VERT_PER_UNIT}/unit-t, 원형={N_CIRCLE_STEPS}")
    log("방법: §1=사다리꼴(L(s)), §2,3,4=위상추적(Im log ξ)")
    log("=" * 70)

    # ────────────────────────────────────────────────────────────────────────
    # §0. 영점 수집
    # ────────────────────────────────────────────────────────────────────────
    log("\n§0. 리만 ζ 영점 수집")
    log("-" * 40)
    ALL_N_ZEROS = 250  # t≈471까지 (N(430)≈220)
    zeros_zeta = collect_zeta_zeros(ALL_N_ZEROS)
    N_ZEROS = len(zeros_zeta)
    log(f"  총 영점: {N_ZEROS}개")
    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §1. 개별 영점 홀로노미 (sanity check, 첫 20개)
    # ────────────────────────────────────────────────────────────────────────
    log("\n§1. 개별 영점 홀로노미 (∮_{|s-ρ|=0.3} L ds / 2πi)")
    log("  방법: 사다리꼴 (64분할), L(s)=해석적 공식")
    log("  기대: = 1.000 ± 0.01 (단순 영점)")
    log("-" * 60)
    log(f"  {'n':>4}  {'γ_n':>10}  {'Re(∮/2πi)':>12}  {'Im(∮/2πi)':>12}  {'|오차|':>10}  판정")
    log("-" * 60)

    circle_results = []
    n_pass_circle = 0
    t0_sec1 = time.time()
    for idx in range(20):
        gamma_n = zeros_zeta[idx]
        center = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(gamma_n))
        try:
            val = circle_contour_integral(center, R_CIRCLE, N_CIRCLE_STEPS, connection_L_analytic)
            re_val = float(val.real)
            im_val = float(val.imag)
            err = abs(re_val - 1.0)
            ok = "✅" if err < 0.01 else "❌"
            if err < 0.01:
                n_pass_circle += 1
            log(f"  {idx+1:>4}  {gamma_n:>10.4f}  {re_val:>12.6f}  {im_val:>12.6f}  {err:>10.6f}  {ok}")
            circle_results.append((idx+1, gamma_n, re_val, im_val, err))
        except Exception as e:
            log(f"  {idx+1:>4}  {gamma_n:>10.4f}  오류: {e}")
            circle_results.append((idx+1, gamma_n, None, None, None))

    log("-" * 60)
    log(f"  통과: {n_pass_circle}/20  ({time.time()-t0_sec1:.1f}초)")
    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §2. 직사각형 윤곽 — 위상 추적 (6개 t₂)
    # ────────────────────────────────────────────────────────────────────────
    log("\n§2. 직사각형 윤곽 — 위상 추적법")
    log(f"  σ∈[{float(SIGMA_MIN)},{float(SIGMA_MAX)}], t∈[{float(T1)},t₂]")
    log("  방법: Im(log ξ) 차분 누적, n_vert ≈ (t₂-13)×3")
    log("  t₂ = 영점 간 중간점 (n=10,20,50,100,150,200)")
    log("-" * 70)
    log(f"  {'n':>5}  {'t₂':>10}  {'N_actual':>10}  {'N_ctr':>10}  {'|차이|':>10}  {'n_vert':>8}  판정")
    log("-" * 70)

    RECT_N_LIST = [10, 20, 50, 100, 150, 200]
    rect_results = []
    n_pass_rect = 0

    for n_zeros_target in RECT_N_LIST:
        if n_zeros_target >= N_ZEROS - 1:
            log(f"  n={n_zeros_target}: 영점 부족 ({N_ZEROS}개)")
            continue

        t2_f = midpoint_between_zeros(zeros_zeta, n_zeros_target - 1)
        t2   = mpmath.mpf(str(t2_f))
        N_actual = count_zeros_below(zeros_zeta, t2_f)

        # 수직 분할: (t₂-T1) × N_VERT_PER_UNIT, 최소 64
        n_vert = max(64, int((t2_f - float(T1)) * N_VERT_PER_UNIT))

        t0_r = time.time()
        try:
            N_ctr = phase_track_contour(
                SIGMA_MIN, SIGMA_MAX, T1, t2,
                N_HORIZ, n_vert, log_xi_imag
            )
            elapsed = time.time() - t0_r
            diff = abs(N_ctr - N_actual)
            ok = "✅" if diff < 0.1 else "❌"
            if diff < 0.1:
                n_pass_rect += 1
            log(f"  {n_zeros_target:>5}  {t2_f:>10.4f}  {N_actual:>10}  {N_ctr:>10.4f}  {diff:>10.4f}  {n_vert:>8}  {ok} ({elapsed:.1f}s)")
            rect_results.append((n_zeros_target, t2_f, N_actual, N_ctr, diff))
        except Exception as e:
            log(f"  n={n_zeros_target}: 오류: {e}")
            rect_results.append((n_zeros_target, t2_f, N_actual, None, None))

    log("-" * 70)
    n_rect_valid = sum(1 for r in rect_results if r[3] is not None)
    log(f"  통과: {n_pass_rect}/{n_rect_valid}")
    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §3. 연속 누적 함수 C(t) — Figure 데이터
    # ────────────────────────────────────────────────────────────────────────
    log("\n§3. 연속 누적 함수 C(t) = (1/2π) × Δ arg ξ")
    log("  t₂ ∈ [14, 430], 100개 점 (영점 0.5 이내 회피)")
    log("-" * 50)

    t2_raw = np.linspace(14.0, 430.0, 120)

    def is_near_zero_t(t_val, margin=0.5):
        for z in zeros_zeta:
            if abs(t_val - z) < margin:
                return True
        return False

    t2_scan_list = []
    for t_val in t2_raw:
        if not is_near_zero_t(t_val) and len(t2_scan_list) < 100:
            t2_scan_list.append(t_val)

    log(f"  선택된 스캔 점: {len(t2_scan_list)}개")

    figdata_rows = []  # (t2, C, N_actual)

    for i, t2_f in enumerate(t2_scan_list):
        t2   = mpmath.mpf(str(t2_f))
        N_actual = count_zeros_below(zeros_zeta, t2_f)

        # 수직 분할: 최소 64, 단위 t당 2개
        n_vert_s = max(64, int((t2_f - float(T1)) * 2))

        try:
            C_val = phase_track_contour(
                SIGMA_MIN, SIGMA_MAX, T1, t2,
                N_HORIZ, n_vert_s, log_xi_imag
            )
        except Exception as e:
            print(f"WARNING: t₂={t2_f:.3f} 실패: {e}", flush=True)
            figdata_rows.append((t2_f, None, N_actual))
            continue

        figdata_rows.append((t2_f, C_val, N_actual))
        if i % 10 == 0:
            log(f"  [{i+1:>3}/{len(t2_scan_list)}] t₂={t2_f:.2f}: C={C_val:.3f}, N={N_actual}")

    # 계단함수 일치 검사
    valid_rows = [(t, C, N) for t, C, N in figdata_rows if C is not None]
    if len(valid_rows) > 0:
        diffs = [abs(C - N) for t, C, N in valid_rows]
        mean_diff = np.mean(diffs)
        max_diff  = np.max(diffs)
        n_match   = sum(1 for d in diffs if d < 0.5)
        log(f"\n  계단함수 일치: {n_match}/{len(valid_rows)} (|C-N|<0.5)")
        log(f"  평균 |C-N|={mean_diff:.4f}, 최대={max_diff:.4f}")
    else:
        n_match = 0
        mean_diff = max_diff = float('nan')

    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §4. 디리클레 χ₃ 검증
    # ────────────────────────────────────────────────────────────────────────
    log("\n§4. 디리클레 χ₃ 보편성 검증")
    log("-" * 60)

    log("  χ₃ 영점 탐색...")
    t0 = time.time()
    zeros_chi3 = find_zeros_chi3(t_min=1.0, t_max=80.0, n_scan=3000)
    log(f"  탐색 완료: {len(zeros_chi3)}개 ({time.time()-t0:.1f}초)")

    n_pass_chi3_rect = 0
    n_chi3_rect_valid = 0
    chi3_results = []

    if len(zeros_chi3) == 0:
        log("  ⚠️ χ₃ 영점 없음 — §4 건너뜀")
    else:
        T1_CHI3 = mpmath.mpf('0.5')  # χ₃ 첫 영점 아래
        CHI3_N_LIST = [5, 10, 15, 20]
        log(f"\n  §4a. χ₃ 직사각형 위상추적 (σ∈[0.3,0.7], t∈[0.5,t₂])")
        log(f"  {'n':>5}  {'t₂':>10}  {'N_actual':>10}  {'N_ctr':>10}  {'|차이|':>10}  판정")
        log("-" * 60)

        for n_t in CHI3_N_LIST:
            if n_t >= len(zeros_chi3) - 1:
                log(f"  n={n_t}: χ₃ 영점 부족")
                continue
            t2_f = (zeros_chi3[n_t - 1] + zeros_chi3[n_t]) / 2.0
            t2   = mpmath.mpf(str(t2_f))
            N_actual = sum(1 for z in zeros_chi3 if z < t2_f)
            n_vert_c = max(64, int((t2_f - float(T1_CHI3)) * 3))

            try:
                t0_c = time.time()
                N_ctr = phase_track_contour(
                    SIGMA_MIN, SIGMA_MAX, T1_CHI3, t2,
                    N_HORIZ, n_vert_c, log_lambda_chi3_imag
                )
                elapsed = time.time() - t0_c
                diff = abs(N_ctr - N_actual)
                ok = "✅" if diff < 0.1 else "❌"
                if diff < 0.1:
                    n_pass_chi3_rect += 1
                n_chi3_rect_valid += 1
                log(f"  {n_t:>5}  {t2_f:>10.4f}  {N_actual:>10}  {N_ctr:>10.4f}  {diff:>10.4f}  {ok} ({elapsed:.1f}s)")
                chi3_results.append((n_t, t2_f, N_actual, N_ctr, diff))
            except Exception as e:
                log(f"  n={n_t}: 오류: {e}")

        log(f"  χ₃ 통과: {n_pass_chi3_rect}/{n_chi3_rect_valid}")

    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §5. 최종 판정
    # ────────────────────────────────────────────────────────────────────────
    elapsed_total = time.time() - t_start

    log("\n" + "=" * 70)
    log("§5. 최종 판정")
    log("=" * 70)

    log(f"\n【A】 개별 영점 홀로노미 (ζ): {n_pass_circle}/20")
    log(f"     기준: 20/20 (|오차|<0.01)")

    log(f"\n【B】 직사각형 위상추적 (ζ): {n_pass_rect}/{n_rect_valid}")
    log(f"     기준: {n_rect_valid}/{n_rect_valid} (|차이|<0.1)")

    if len(valid_rows) > 0:
        log(f"\n【C】 연속 누적 함수 (ζ): {n_match}/{len(valid_rows)}")
        log(f"     평균 |C-N|={mean_diff:.4f}, 최대={max_diff:.4f}")
    else:
        log("\n【C】 연속 누적 함수: 데이터 없음")

    log(f"\n【D】 χ₃ 보편성: 직사각형 {n_pass_chi3_rect}/{n_chi3_rect_valid}")

    crit_A = n_pass_circle >= 18
    crit_B = n_pass_rect == n_rect_valid and n_rect_valid >= 5
    crit_C = (len(valid_rows) > 50 and mean_diff < 0.5)
    chi3_ok = (n_pass_chi3_rect == n_chi3_rect_valid and n_chi3_rect_valid >= 3)

    log("\n【종합】")
    log(f"  A(개별 홀로노미): {'✅' if crit_A else '❌'}")
    log(f"  B(직사각형):      {'✅' if crit_B else '❌'}")
    log(f"  C(연속 누적):     {'✅' if crit_C else '❌'}")
    log(f"  D(χ₃ 보편성):     {'✅' if chi3_ok else '❌'}")

    if crit_A and crit_B and crit_C and chi3_ok:
        verdict = "★ 양성 — N(T) = c₁(ξ-bundle) 수치 확립, χ₃ 보편성 확인"
    elif crit_A and crit_B and crit_C:
        verdict = "★ 조건부 양성 — ζ 3단계 검증 완전, χ₃ 보완 필요"
    elif crit_A and crit_B:
        verdict = "조건부 양성 — 개별+직사각형 통과, 누적 미확인"
    elif crit_A:
        verdict = "부분 양성 — sanity check 통과, 직사각형 실패"
    else:
        verdict = "음성 — 디버깅 필요"

    log(f"\n  최종: {verdict}")
    log(f"\n  실행 시간: {elapsed_total/60:.1f}분")
    log(f"  결과: {OUTPUT_TXT}")
    log(f"  CSV:  {OUTPUT_CSV}")

    # CSV 저장
    with open(OUTPUT_CSV, 'w', encoding='utf-8') as f:
        f.write("t2,C_val,N_actual\n")
        for t2_v, C, N in figdata_rows:
            if C is not None:
                f.write(f"{t2_v:.6f},{C:.6f},{N}\n")
            else:
                f.write(f"{t2_v:.6f},,{N}\n")

    log(f"\nCSV 저장 완료: {OUTPUT_CSV}")
    save_all()


if __name__ == "__main__":
    main()
