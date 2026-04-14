"""
=============================================================================
[Project RDL] Gauss-Bonnet 면적분 검증 — ∫∫ F₂ dσ∧dt = 2π·N
=============================================================================
사이클 34 / 수학자 지시 2026-04-14 10:01

목적:
  Chern 수 노선(결과 #22-23, 윤곽 적분)의 면적분 버전을 독립 확인.
  곡률 2-form F₂의 영역 적분이 2π × (영점 개수)임을 수치 검증.

수학적 배경:
  U(1) 다발에서 제1 Chern 수 = (1/2π) ∫∫_R F₂ dσ∧dt = N (영점 개수)

  F₂의 정의 — 격자 플라켓 곡률 (Lattice Plaquette Curvature):
    θ(σ,t) = arg(ξ(σ+it)) = Im(log ξ(σ+it))
    각 격자 셀 (i,j)의 곡률 = 셀 경계를 따른 위상 감김 (wrapped)
    Σ(모든 셀 곡률) = 2π × N  (Stokes 정리의 이산 버전)

  이는 격자 게이지 이론의 표준 방법: ξ의 영점이 위상 와동(vortex)이며,
  각 와동은 정확히 2π의 자속(flux)을 기여.

  참고: ∂_σ Im(L) 단독은 σ₁↔σ₂ 대칭 직사각형에서 0을 주므로 사용 불가.
  대신 전체 곡률 밀도 — arg(ξ)의 라플라시안 — 을 격자 위에서 계산.
  Stokes: ∫∫ Δ(log|ξ|) = ∮ ∇(log|ξ|)·n̂ dl = ∮ Im(L ds) = 2πN (논증 원리)

메서드:
  §A — 단일 영점: ζ 첫 영점 γ₁=14.1347 포함 영역
  §B — 다영점 누적: k=1,5,10,20 영점 포함 영역
  §C — 빈 영역: 영점 불포함 구간
  §D — 디리클레 확장: χ₅=(·/5) 짝수 지표

  격자 해상도: Nσ=64 (σ방향), Nt는 영역 높이에 비례 (최소 64, 단위당 10)
  dps=100

출력:
  results/gauss_bonnet_area_integral.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT_TXT = os.path.join(BASE_DIR, "results", "gauss_bonnet_area_integral.txt")
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# ─── 설정 ─────────────────────────────────────────────────────────────────────
DPS = 100
mpmath.mp.dps = DPS

SIGMA_MIN = 0.3
SIGMA_MAX = 0.7
N_SIGMA   = 64    # σ 방향 격자 점 수

# χ₅ = (·/5) Legendre symbol, 짝수 (a=0)
CHI5_LIST = [0, 1, -1, -1, 1]
CHI5_Q    = 5

# ─── 로그 ─────────────────────────────────────────────────────────────────────
log_lines = []
def log(msg=""):
    print(msg, flush=True)
    log_lines.append(str(msg))

def save_all():
    with open(OUTPUT_TXT, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))


# ═══════════════════════════════════════════════════════════════════════════════
# 위상 함수 — θ(σ,t) = Im(log ξ(σ+it))
# ═══════════════════════════════════════════════════════════════════════════════

def log_xi_imag(sigma, t):
    """
    ζ 함수용: θ(σ,t) = Im(log ξ(σ+it))
    log ξ(s) = log(s) + log(s-1) - (s/2)log(π) + loggamma(s/2) + log(ζ(s))
    ζ(s) 영점 근처이면 None 반환.
    """
    s = mpmath.mpc(sigma, t)
    z_val = mpmath.zeta(s)
    if abs(z_val) < mpmath.mpf(10) ** (-DPS + 15):
        return None
    lxi = (mpmath.log(s)
           + mpmath.log(s - 1)
           - (s / 2) * mpmath.log(mpmath.pi)
           + mpmath.loggamma(s / 2)
           + mpmath.log(z_val))
    return float(lxi.imag)


def log_lambda_chi5_imag(sigma, t):
    """
    χ₅=(·/5) 용: θ(σ,t) = Im(log Λ(σ+it, χ₅))
    Λ(s,χ) = (q/π)^(s/2) Γ(s/2) L(s,χ)  [a=0, 짝수]
    L(s,χ) 영점 근처이면 None 반환.
    """
    s = mpmath.mpc(sigma, t)
    L_val = mpmath.dirichlet(s, CHI5_LIST)
    if abs(L_val) < mpmath.mpf(10) ** (-DPS + 15):
        return None
    lL = ((s / 2) * mpmath.log(mpmath.mpf(CHI5_Q) / mpmath.pi)
          + mpmath.loggamma(s / 2)
          + mpmath.log(L_val))
    return float(lL.imag)


# ═══════════════════════════════════════════════════════════════════════════════
# 격자 플라켓 곡률 (Lattice Plaquette Curvature)
# ═══════════════════════════════════════════════════════════════════════════════

def wrap_phase(dphi):
    """위상 차분을 [-π, π) 범위로 감김."""
    return dphi - 2.0 * np.pi * np.round(dphi / (2.0 * np.pi))


def compute_plaquette_flux(theta_grid):
    """
    2D 위상 격자에서 각 플라켓(셀)의 자속 계산.

    θ_{i,j}가 격자 위의 위상 값일 때, 플라켓 (i,j)의 자속:
      F = wrap(θ_{i+1,j}-θ_{i,j}) + wrap(θ_{i+1,j+1}-θ_{i+1,j})
        + wrap(θ_{i,j+1}-θ_{i+1,j+1}) + wrap(θ_{i,j}-θ_{i,j+1})

    각 영점(와동)은 ≈±2π, 빈 셀은 ≈0.

    Parameters:
        theta_grid: (Nσ, Nt) 배열, θ(σ_i, t_j) 값

    Returns:
        flux_grid: (Nσ-1, Nt-1) 배열, 각 플라켓의 자속
        total_flux: 전체 자속 합
    """
    Ns, Nt = theta_grid.shape

    # 수평/수직 위상 차분 (wrapped)
    # d_sigma[i,j] = θ[i+1,j] - θ[i,j]  (σ 방향)
    # d_t[i,j]     = θ[i,j+1] - θ[i,j]  (t 방향)
    d_sigma = wrap_phase(theta_grid[1:, :] - theta_grid[:-1, :])    # (Ns-1, Nt)
    d_t     = wrap_phase(theta_grid[:, 1:] - theta_grid[:, :-1])    # (Ns, Nt-1)

    # 플라켓 자속: 반시계 경로를 따른 감김
    # F[i,j] = d_sigma[i,j] + d_t[i+1,j] - d_sigma[i,j+1] - d_t[i,j]
    flux = (d_sigma[:, :-1]      # 하변: σ 방향 →
            + d_t[1:, :]         # 우변: t 방향 ↑
            - d_sigma[:, 1:]     # 상변: σ 방향 ← (부호 반전)
            - d_t[:-1, :])       # 좌변: t 방향 ↓ (부호 반전)

    total_flux = float(np.sum(flux))
    return flux, total_flux


def compute_area_integral(phase_fn, sigma_min, sigma_max, t_min, t_max,
                          n_sigma=64, n_t=None, label=""):
    """
    2D 영역에서 격자 플라켓 곡률의 면적분 계산.

    Parameters:
        phase_fn: (sigma, t) → Im(log f(σ+it)) 함수
        sigma_min, sigma_max, t_min, t_max: 영역 경계
        n_sigma: σ 방향 격자 점 수
        n_t: t 방향 격자 점 수 (None이면 높이에 비례 자동 설정)
        label: 로그용 라벨

    Returns:
        total_flux: 전체 자속 (≈ 2π × N)
        flux_grid: 플라켓별 자속 배열
        sigma_arr, t_arr: 격자 좌표
    """
    height = t_max - t_min
    if n_t is None:
        n_t = max(64, int(height * 10))  # 단위당 최소 10개 점

    sigma_arr = np.linspace(sigma_min, sigma_max, n_sigma)
    t_arr     = np.linspace(t_min, t_max, n_t)

    log(f"  격자: Nσ={n_sigma}, Nt={n_t}, 총 {n_sigma*n_t}점")
    log(f"  영역: σ∈[{sigma_min},{sigma_max}], t∈[{t_min:.4f},{t_max:.4f}]")

    # 위상 격자 계산
    theta_grid = np.zeros((n_sigma, n_t))
    skip_count = 0

    t0 = time.time()
    for i, sig in enumerate(sigma_arr):
        for j, t in enumerate(t_arr):
            val = phase_fn(float(sig), float(t))
            if val is None:
                skip_count += 1
                # 영점 근처: 이웃 값으로 보간 (간단히 이전 값 사용)
                if j > 0:
                    theta_grid[i, j] = theta_grid[i, j-1]
                elif i > 0:
                    theta_grid[i, j] = theta_grid[i-1, j]
                else:
                    theta_grid[i, j] = 0.0
            else:
                theta_grid[i, j] = val

    elapsed = time.time() - t0
    log(f"  계산 시간: {elapsed:.1f}초, 건너뜀: {skip_count}점")

    # 플라켓 자속 계산
    flux_grid, total_flux = compute_plaquette_flux(theta_grid)

    N_est = total_flux / (2 * np.pi)
    log(f"  전체 자속: {total_flux:.6f}")
    log(f"  N = 자속/(2π) = {N_est:.6f}")

    return total_flux, flux_grid, sigma_arr, t_arr


# ═══════════════════════════════════════════════════════════════════════════════
# 영점 탐색 (mpmath 기본 제공)
# ═══════════════════════════════════════════════════════════════════════════════

def get_zeta_zeros(n_max):
    """ζ 영점 γ₁...γ_{n_max} 반환"""
    zeros = []
    for n in range(1, n_max + 1):
        z = mpmath.zetazero(n)
        zeros.append(float(z.imag))
    return zeros


def find_zeros_chi5(t_min=1.0, t_max=60.0, n_scan=6000):
    """χ₅=(·/5) 영점 탐색: |L(1/2+it)| 최소화 방식 (보편적 방법)"""
    ts = np.linspace(t_min, t_max, n_scan)
    half = mpmath.mpf('0.5')

    # |L| 프로파일 스캔
    abs_vals = []
    for t in ts:
        s = half + 1j * mpmath.mpf(str(t))
        val = mpmath.dirichlet(s, CHI5_LIST)
        abs_vals.append(float(abs(val)))
    abs_vals = np.array(abs_vals)

    # 극소점 찾기
    zeros = []
    for i in range(1, len(abs_vals) - 1):
        if abs_vals[i] < abs_vals[i-1] and abs_vals[i] < abs_vals[i+1]:
            if abs_vals[i] < 0.1:  # 후보 필터
                t_guess = ts[i]
                try:
                    def f_abs(t_var):
                        sv = half + 1j * mpmath.mpf(t_var)
                        return abs(mpmath.dirichlet(sv, CHI5_LIST))
                    # findroot에 |L| 사용 (미분 가능하지 않으므로 secant)
                    # 대신 Re(L)과 Im(L)의 동시 0을 찾거나 minimize
                    # 간단히: findroot로 Re와 Im을 동시에
                    def f_re(t_var):
                        sv = half + 1j * mpmath.mpf(t_var)
                        return mpmath.re(mpmath.dirichlet(sv, CHI5_LIST))
                    def f_im(t_var):
                        sv = half + 1j * mpmath.mpf(t_var)
                        return mpmath.im(mpmath.dirichlet(sv, CHI5_LIST))

                    # |L| < threshold인 극소점 → findroot Re 또는 Im
                    t_refined = float(mpmath.findroot(f_re, t_guess))
                    sv_check = half + 1j * mpmath.mpf(str(t_refined))
                    if abs(mpmath.dirichlet(sv_check, CHI5_LIST)) < 1e-6:
                        if not zeros or abs(t_refined - zeros[-1]) > 0.1:
                            zeros.append(t_refined)
                            continue

                    t_refined = float(mpmath.findroot(f_im, t_guess))
                    sv_check = half + 1j * mpmath.mpf(str(t_refined))
                    if abs(mpmath.dirichlet(sv_check, CHI5_LIST)) < 1e-6:
                        if not zeros or abs(t_refined - zeros[-1]) > 0.1:
                            zeros.append(t_refined)
                except Exception as e:
                    pass  # findroot 실패 — 다음 후보

    zeros.sort()
    return zeros


# ═══════════════════════════════════════════════════════════════════════════════
# 보조: 플라켓 자속 공간 분포 분석
# ═══════════════════════════════════════════════════════════════════════════════

def analyze_flux_distribution(flux_grid, sigma_arr, t_arr, zeros_in_region, label=""):
    """
    영점 위치와 플라켓 자속의 공간 분포를 비교.
    각 영점 근방의 자속 기여, 빈 영역의 자속 잔차 등을 보고.
    """
    hs = sigma_arr[1] - sigma_arr[0]
    ht = t_arr[1] - t_arr[0]

    log(f"  --- 자속 공간 분포 ({label}) ---")

    if len(zeros_in_region) == 0:
        log(f"  영점 없음 → 전체 자속이 ≈0이어야 함")
        return

    # 각 영점 근방(±1.0)의 자속 기여
    for k, gamma in enumerate(zeros_in_region[:5]):  # 최대 5개만 상세 출력
        # γ 근방의 t 인덱스 범위
        t_lo = gamma - 1.0
        t_hi = gamma + 1.0
        j_lo = max(0, np.searchsorted(t_arr, t_lo) - 1)
        j_hi = min(flux_grid.shape[1], np.searchsorted(t_arr, t_hi))

        local_flux = float(np.sum(flux_grid[:, j_lo:j_hi]))
        log(f"  영점 #{k+1} (γ={gamma:.4f}): 근방 자속 = {local_flux:.4f} "
            f"(≈{local_flux/(2*np.pi):.4f} × 2π)")

    if len(zeros_in_region) > 5:
        log(f"  ... (영점 {len(zeros_in_region)}개 중 5개만 출력)")


# ═══════════════════════════════════════════════════════════════════════════════
# 메인 실험
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    log("=" * 77)
    log(f"[Project RDL] Gauss-Bonnet 면적분 검증")
    log(f"∫∫_R F₂ dσ∧dt = 2π × N (격자 플라켓 곡률)")
    log(f"실행: {now}")
    log("=" * 77)
    log()
    log("방법: 격자 플라켓 곡률 (Lattice Plaquette Curvature)")
    log("  θ(σ,t) = Im(log ξ(σ+it))를 2D 격자 위에서 계산")
    log("  각 격자 셀의 위상 감김(wrapped sum) = 2π × (셀 내 영점 수)")
    log("  전체 합 = 2π × N (Stokes 정리 ≡ Gauss-Bonnet)")
    log()
    log(f"설정: dps={DPS}, Nσ={N_SIGMA}, σ∈[{SIGMA_MIN},{SIGMA_MAX}]")
    log()

    # ─── ζ 영점 준비 ─────────────────────────────────────────────────────────
    log("영점 준비...")
    zeta_zeros = get_zeta_zeros(25)  # 25개면 t≈100까지 충분
    log(f"  ζ 영점 {len(zeta_zeros)}개 로드 (γ₁={zeta_zeros[0]:.4f} ~ γ₂₅={zeta_zeros[-1]:.4f})")
    log()

    results = {}

    # ═══════════════════════════════════════════════════════════════════════════
    # §A: 단일 영점 (ζ 첫 영점 γ₁=14.1347)
    # ═══════════════════════════════════════════════════════════════════════════
    log("═" * 77)
    log("§A: 단일 영점 — ζ 첫 영점 γ₁=14.1347 포함 영역")
    log("═" * 77)
    log()
    log("R₁ = [0.3, 0.7] × [13.0, 15.5]")
    log("기대: ∫∫ F₂ ≈ 2π = 6.2832 (영점 1개)")
    log()

    t_min_A = 13.0
    t_max_A = 15.5

    flux_A, fgrid_A, sig_A, t_A = compute_area_integral(
        log_xi_imag, SIGMA_MIN, SIGMA_MAX, t_min_A, t_max_A,
        n_sigma=N_SIGMA, label="§A"
    )

    N_A = flux_A / (2 * np.pi)
    err_A = abs(N_A - 1.0)
    pass_A = err_A < 0.05
    log(f"  → N = {N_A:.6f}, 기대=1, |오차|={err_A:.6f}")
    log(f"  → 판정: {'✅ PASS' if pass_A else '❌ FAIL'} (기준: |오차| < 0.05)")

    # 공간 분포
    zeros_in_A = [g for g in zeta_zeros if t_min_A < g < t_max_A]
    analyze_flux_distribution(fgrid_A, sig_A, t_A, zeros_in_A, "§A")

    results['A'] = {'N': N_A, 'expected': 1, 'err': err_A, 'pass': pass_A}
    log()
    save_all()

    # ═══════════════════════════════════════════════════════════════════════════
    # §B: 다영점 누적 (k=1, 5, 10, 20)
    # ═══════════════════════════════════════════════════════════════════════════
    log("═" * 77)
    log("§B: 다영점 누적 — ∫∫ F₂ ≈ 2πk (k=1,5,10,20)")
    log("═" * 77)
    log()

    t_min_B = 0.5  # 첫 영점(14.13) 아래에서 시작
    results_B = []

    for k in [1, 5, 10, 20]:
        log(f"--- k={k} ---")

        # t₂: k번째와 k+1번째 영점의 중간점
        if k < len(zeta_zeros):
            t_max_B = (zeta_zeros[k-1] + zeta_zeros[k]) / 2.0
        else:
            t_max_B = zeta_zeros[k-1] + 1.0

        log(f"  R = [0.3,0.7] × [{t_min_B}, {t_max_B:.4f}]")
        log(f"  기대: ∫∫ F₂ ≈ 2π×{k} = {2*np.pi*k:.4f}")

        flux_Bk, fgrid_Bk, sig_Bk, t_Bk = compute_area_integral(
            log_xi_imag, SIGMA_MIN, SIGMA_MAX, t_min_B, t_max_B,
            n_sigma=N_SIGMA, label=f"§B k={k}"
        )

        N_Bk = flux_Bk / (2 * np.pi)
        err_Bk = abs(N_Bk - k) / k
        pass_Bk = err_Bk < 0.05

        log(f"  → N = {N_Bk:.6f}, 기대={k}, 상대오차={err_Bk:.6f}")
        log(f"  → 판정: {'✅ PASS' if pass_Bk else '❌ FAIL'}")

        results_B.append({'k': k, 'N': N_Bk, 'err': err_Bk, 'pass': pass_Bk,
                          't_max': t_max_B})
        log()
        save_all()

    results['B'] = results_B

    # ═══════════════════════════════════════════════════════════════════════════
    # §C: 빈 영역 (영점 불포함)
    # ═══════════════════════════════════════════════════════════════════════════
    log("═" * 77)
    log("§C: 빈 영역 — 영점 불포함 구간")
    log("═" * 77)
    log()

    # γ₁=14.1347, γ₂=21.0220 사이에 영점 없음
    t_min_C = zeta_zeros[0] + 0.5   # ≈ 14.63
    t_max_C = zeta_zeros[1] - 0.5   # ≈ 20.52

    log(f"  R_empty = [0.3,0.7] × [{t_min_C:.4f}, {t_max_C:.4f}]")
    log(f"  γ₁={zeta_zeros[0]:.4f}, γ₂={zeta_zeros[1]:.4f}")
    log(f"  기대: ∫∫ F₂ ≈ 0 (영점 없음)")
    log()

    flux_C, fgrid_C, sig_C, t_C = compute_area_integral(
        log_xi_imag, SIGMA_MIN, SIGMA_MAX, t_min_C, t_max_C,
        n_sigma=N_SIGMA, label="§C"
    )

    N_C = flux_C / (2 * np.pi)
    pass_C = abs(flux_C) < 0.5

    log(f"  → 전체 자속 = {flux_C:.6f}, N = {N_C:.6f}")
    log(f"  → 판정: {'✅ PASS' if pass_C else '❌ FAIL'} (기준: |자속| < 0.5)")

    results['C'] = {'flux': flux_C, 'N': N_C, 'pass': pass_C}
    log()
    save_all()

    # ═══════════════════════════════════════════════════════════════════════════
    # §D: 디리클레 확장 — χ₅=(·/5), 짝수 지표
    # ═══════════════════════════════════════════════════════════════════════════
    log("═" * 77)
    log("§D: 디리클레 확장 — χ₅=(·/5) mod 5, 짝수 지표 (a=0)")
    log("═" * 77)
    log()

    # χ₅ 영점 탐색
    log("χ₅ 영점 탐색 (|L| 최소화 방식)...")
    chi5_zeros = find_zeros_chi5(t_min=1.0, t_max=60.0, n_scan=6000)
    log(f"  {len(chi5_zeros)}개 영점 발견")
    if len(chi5_zeros) >= 5:
        log(f"  첫 5개: {[f'{z:.4f}' for z in chi5_zeros[:5]]}")

    if len(chi5_zeros) < 5:
        log("  ⚠️ χ₅ 영점 5개 미만 — §D 건너뜀")
        results['D'] = {'pass': False, 'reason': '영점 부족'}
    else:
        results_D = []

        # §D-1: 단일 영점
        log()
        log("--- §D-1: χ₅ 단일 영점 ---")
        gamma1_chi5 = chi5_zeros[0]
        t_min_D1 = gamma1_chi5 - 2.0
        t_max_D1 = gamma1_chi5 + 2.0

        log(f"  R = [0.3,0.7] × [{t_min_D1:.4f}, {t_max_D1:.4f}]")
        log(f"  χ₅ 첫 영점: γ₁={gamma1_chi5:.4f}")
        log(f"  기대: ∫∫ F₂ ≈ 2π")

        flux_D1, _, _, _ = compute_area_integral(
            log_lambda_chi5_imag, SIGMA_MIN, SIGMA_MAX, t_min_D1, t_max_D1,
            n_sigma=N_SIGMA, label="§D-1"
        )

        N_D1 = flux_D1 / (2 * np.pi)
        err_D1 = abs(N_D1 - 1.0)
        pass_D1 = err_D1 < 0.05
        log(f"  → N = {N_D1:.6f}, |오차|={err_D1:.6f}")
        log(f"  → 판정: {'✅ PASS' if pass_D1 else '❌ FAIL'}")
        results_D.append({'test': 'D-1 (단일)', 'N': N_D1, 'expected': 1,
                          'err': err_D1, 'pass': pass_D1})
        log()
        save_all()

        # §D-2: 5영점 영역
        log("--- §D-2: χ₅ 5영점 영역 ---")
        t_min_D2 = chi5_zeros[0] - 2.0
        if len(chi5_zeros) > 5:
            t_max_D2 = (chi5_zeros[4] + chi5_zeros[5]) / 2.0
        else:
            t_max_D2 = chi5_zeros[4] + 1.0

        log(f"  R = [0.3,0.7] × [{t_min_D2:.4f}, {t_max_D2:.4f}]")
        log(f"  기대: ∫∫ F₂ ≈ 2π×5 = {10*np.pi:.4f}")

        flux_D2, _, _, _ = compute_area_integral(
            log_lambda_chi5_imag, SIGMA_MIN, SIGMA_MAX, t_min_D2, t_max_D2,
            n_sigma=N_SIGMA, label="§D-2"
        )

        N_D2 = flux_D2 / (2 * np.pi)
        err_D2 = abs(N_D2 - 5.0) / 5.0
        pass_D2 = err_D2 < 0.05
        log(f"  → N = {N_D2:.6f}, 기대=5, 상대오차={err_D2:.6f}")
        log(f"  → 판정: {'✅ PASS' if pass_D2 else '❌ FAIL'}")
        results_D.append({'test': 'D-2 (5영점)', 'N': N_D2, 'expected': 5,
                          'err': err_D2, 'pass': pass_D2})

        results['D'] = results_D

    log()
    save_all()

    # ═══════════════════════════════════════════════════════════════════════════
    # 종합
    # ═══════════════════════════════════════════════════════════════════════════
    log()
    log("═" * 77)
    log("종합 결과")
    log("═" * 77)
    log()

    # 요약 테이블
    log("┌──────────┬───────────┬────────┬──────────┬────────┐")
    log("│ 테스트   │ ∫∫F₂/(2π) │ 기대 N │ |오차|   │ 판정   │")
    log("├──────────┼───────────┼────────┼──────────┼────────┤")

    # A
    r = results['A']
    log(f"│ §A 단일  │ {r['N']:9.4f} │ {r['expected']:6d} │ {r['err']:.6f} │ {'✅' if r['pass'] else '❌'}     │")

    # B
    for rb in results['B']:
        log(f"│ §B k={rb['k']:>2d}  │ {rb['N']:9.4f} │ {rb['k']:6d} │ {rb['err']:.6f} │ {'✅' if rb['pass'] else '❌'}     │")

    # C
    r = results['C']
    log(f"│ §C 빈영역│ {r['N']:9.4f} │      0 │ {abs(r['N']):.6f} │ {'✅' if r['pass'] else '❌'}     │")

    # D
    if isinstance(results.get('D'), list):
        for rd in results['D']:
            exp = rd['expected']
            log(f"│ {rd['test']:8s} │ {rd['N']:9.4f} │ {exp:6d} │ {rd['err']:.6f} │ {'✅' if rd['pass'] else '❌'}     │")

    log("└──────────┴───────────┴────────┴──────────┴────────┘")
    log()

    # 통과율
    all_passes = [results['A']['pass']]
    all_passes += [rb['pass'] for rb in results['B']]
    all_passes += [results['C']['pass']]
    if isinstance(results.get('D'), list):
        all_passes += [rd['pass'] for rd in results['D']]

    n_pass = sum(all_passes)
    n_total = len(all_passes)

    log(f"통과: {n_pass}/{n_total}")
    log()

    # 물리적 해석
    log("─── 물리적 해석 ───")
    log()
    log("격자 플라켓 곡률은 U(1) 다발의 곡률 2-form을 이산화한 것.")
    log("각 격자 셀의 위상 감김(plaquette flux)은:")
    log("  - 영점(위상 와동) 포함 시: ≈ ±2π")
    log("  - 영점 미포함 시: ≈ 0")
    log("전체 합 = ∮_∂R d(arg ξ) = 2πN  (Stokes 정리)")
    log()
    log("이는 결과 #22-23 (윤곽 적분 ∮ L ds = 2πiN)의")
    log("면적분 등가물: ∫∫_R F₂ dσ∧dt = 2πN (Gauss-Bonnet)")
    log()

    if n_pass == n_total:
        log("★ 판정: 양성 (완전) — Gauss-Bonnet for ξ-bundle 확립")
    elif n_pass >= n_total - 1:
        log("★ 판정: 조건부 양성 — 대부분 통과, 일부 정밀도 개선 필요")
    else:
        log("★ 판정: 부분 양성/음성 — 수학자 재검토 필요")

    elapsed_total = time.time() - t_start
    log()
    log(f"총 실행 시간: {elapsed_total:.1f}초 ({elapsed_total/60:.1f}분)")
    log(f"dps={DPS}, Nσ={N_SIGMA}")

    save_all()
    log()
    log(f"결과 저장: {OUTPUT_TXT}")


if __name__ == "__main__":
    main()
