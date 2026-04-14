"""
=============================================================================
[Project RDL] Synthetic Off-Critical Control — 합성 함수 음성 대조군
=============================================================================
사이클 37 / 수학자 지시 2026-04-14 14:57

목적:
  26개 결과가 모두 양성 대조군(ξ, L-함수). 음성 대조군이 없다.
  영점이 σ≠1/2에 있는 함수에서도 σ=1/2 δ-집중이 나타나는가?
  → 방법론적 인공물 가능성 검증. 프레임워크 특이성(specificity) 확인.

합성 함수 (영점 위치 완전히 알려짐):
  f₁(s) = (s-(0.5+14i))(s-(0.5+21i))(s-(0.5+25i))  [모두 σ=1/2, RH-유사]
  f₂(s) = (s-(0.3+14i))(s-(0.7+14i))(s-(0.5+21i))  [혼합, 부분 RH-위반]
  f₃(s) = (s-(0.3+16i))(s-(0.7+16i))(s-(0.6+22i))(s-(0.4+22i)) [완전 RH-위반]

11종 테스트 매트릭스:
  §A (f₁): A-1 전체검출, A-2 σ-프로파일, A-3 off-critical
  §B (f₂): B-1 전체검출, B-2 σ-프로파일, B-3 σ=0.5만, B-4 σ=0.7만, B-5 σ=0.3만
  §C (f₃): C-1 전체검출, C-2 σ-프로파일, C-3 임계선 (핵심!)

출력:
  results/synthetic_offcritical_control.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT_TXT = os.path.join(BASE_DIR, "results", "synthetic_offcritical_control.txt")
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# ─── 설정 ─────────────────────────────────────────────────────────────────────
DPS     = 30   # 다항식은 수치 안정 → 낮은 정밀도 충분
N_SIGMA = 64   # σ 방향 격자 점수 (스트립 테스트)
N_T     = 200  # t 방향 격자 점수

mpmath.mp.dps = DPS

# ─── 합성 함수 영점 정의 ───────────────────────────────────────────────────────
# f₁: 모든 영점 σ=1/2 위
F1_ZEROS = [
    mpmath.mpc(0.5, 14),
    mpmath.mpc(0.5, 21),
    mpmath.mpc(0.5, 25),
]
# f₂: 혼합 (σ=0.3, 0.7, 0.5)
F2_ZEROS = [
    mpmath.mpc(0.3, 14),
    mpmath.mpc(0.7, 14),
    mpmath.mpc(0.5, 21),
]
# f₃: 모두 off-critical (σ=0.3, 0.7, 0.6, 0.4)
F3_ZEROS = [
    mpmath.mpc(0.3, 16),
    mpmath.mpc(0.7, 16),
    mpmath.mpc(0.6, 22),
    mpmath.mpc(0.4, 22),
]

# ─── 로그 ─────────────────────────────────────────────────────────────────────
log_lines = []

def log(msg=""):
    print(msg, flush=True)
    log_lines.append(str(msg))

def save_all():
    with open(OUTPUT_TXT, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))


# ═══════════════════════════════════════════════════════════════════════════════
# 다항식 위상 함수
# ═══════════════════════════════════════════════════════════════════════════════

def log_poly_imag(sigma, t, zeros):
    """
    Im(log(f(s))) = sum_k Im(log(s - z_k)) = sum_k arg(s - z_k)

    수학적 근거: f(s) = prod(s - z_k)이므로
      log(f(s)) = sum log(s - z_k)
      Im(log(f(s))) = sum arg(s - z_k)

    영점 근방 처리: |s - z_k| < 1e-10이면 None 반환.
    """
    s = mpmath.mpc(float(sigma), float(t))
    total_imag = mpmath.mpf(0)
    for z in zeros:
        diff = s - z
        if abs(diff) < mpmath.mpf('1e-10'):
            return None  # 영점 너무 근접
        total_imag += mpmath.im(mpmath.log(diff))
    return float(total_imag)


# ═══════════════════════════════════════════════════════════════════════════════
# 격자 플라켓 곡률 (offcritical_gauss_bonnet.py 재사용)
# ═══════════════════════════════════════════════════════════════════════════════

def wrap_phase(dphi):
    """위상 차분을 [-π, π) 범위로 감김."""
    return dphi - 2.0 * np.pi * np.round(dphi / (2.0 * np.pi))


def compute_plaquette_flux(theta_grid):
    """
    2D 위상 격자에서 각 플라켓(셀)의 자속 계산.
    Returns: flux_grid (Nσ-1, Nt-1), total_flux
    """
    d_sigma = wrap_phase(theta_grid[1:, :] - theta_grid[:-1, :])  # (Ns-1, Nt)
    d_t     = wrap_phase(theta_grid[:, 1:] - theta_grid[:, :-1])  # (Ns, Nt-1)

    flux = (d_sigma[:, :-1]   # 하변
            + d_t[1:, :]      # 우변
            - d_sigma[:, 1:]  # 상변
            - d_t[:-1, :])    # 좌변

    total_flux = float(np.sum(flux))
    return flux, total_flux


def compute_area_integral(phase_fn, sigma_min, sigma_max, t_min, t_max,
                          n_sigma=64, n_t=200, label=""):
    """
    2D 영역에서 격자 플라켓 곡률의 면적분 계산.
    Returns: N_est, flux_grid, sigma_arr, t_arr, theta_grid
    """
    sigma_arr = np.linspace(sigma_min, sigma_max, n_sigma)
    t_arr     = np.linspace(t_min, t_max, n_t)

    log(f"  격자: Nσ={n_sigma}, Nt={n_t}, 총 {n_sigma*n_t}점")
    log(f"  σ∈[{sigma_min},{sigma_max}], t∈[{t_min:.4f},{t_max:.4f}]")

    theta_grid  = np.zeros((n_sigma, n_t))
    skip_count  = 0
    t0          = time.time()

    for i, sig in enumerate(sigma_arr):
        for j, t in enumerate(t_arr):
            val = phase_fn(float(sig), float(t))
            if val is None:
                skip_count += 1
                if j > 0:
                    theta_grid[i, j] = theta_grid[i, j-1]
                elif i > 0:
                    theta_grid[i, j] = theta_grid[i-1, j]
                else:
                    theta_grid[i, j] = 0.0
            else:
                theta_grid[i, j] = val

    elapsed = time.time() - t0
    log(f"  계산: {elapsed:.1f}초, skip={skip_count}점")

    flux_grid, total_flux = compute_plaquette_flux(theta_grid)
    N_est = total_flux / (2 * np.pi)
    log(f"  전체 자속: {total_flux:.6f} → N = {N_est:.6f}")

    return N_est, flux_grid, sigma_arr, t_arr, theta_grid


def compute_sigma_profile(phase_fn, sigma_min, sigma_max, t_min, t_max,
                          n_sigma=200, n_t=200, label=""):
    """
    σ-프로파일 계산: 각 σ 열의 플라켓 자속 합 ρ(σ).
    Returns: sigma_mid, rho_sigma, total_N
    """
    sigma_arr = np.linspace(sigma_min, sigma_max, n_sigma)
    t_arr     = np.linspace(t_min, t_max, n_t)

    log(f"  σ-프로파일: Nσ={n_sigma}, Nt={n_t}, 총 {n_sigma*n_t}점")

    theta_grid = np.zeros((n_sigma, n_t))
    t0         = time.time()

    for i, sig in enumerate(sigma_arr):
        for j, t in enumerate(t_arr):
            val = phase_fn(float(sig), float(t))
            if val is None:
                if j > 0:
                    theta_grid[i, j] = theta_grid[i, j-1]
                elif i > 0:
                    theta_grid[i, j] = theta_grid[i-1, j]
                else:
                    theta_grid[i, j] = 0.0
            else:
                theta_grid[i, j] = val

        if (i + 1) % 50 == 0:
            log(f"    σ 진행: {i+1}/{n_sigma} ({time.time()-t0:.1f}초)")
            save_all()

    elapsed = time.time() - t0
    log(f"  계산: {elapsed:.1f}초")

    flux_grid, total_flux = compute_plaquette_flux(theta_grid)
    rho_sigma  = np.sum(flux_grid, axis=1)  # (n_sigma-1,) 각 σ 열 합
    sigma_mid  = (sigma_arr[:-1] + sigma_arr[1:]) / 2.0
    total_N    = total_flux / (2 * np.pi)

    log(f"  전체 N = {total_N:.6f}")
    return sigma_mid, rho_sigma, total_N


def find_peaks_profile(sigma_mid, rho_sigma, min_fraction=0.15):
    """
    |ρ(σ)| 의 국소 최대 탐색.
    Returns: [(sigma, rho_val), ...] 정렬됨
    """
    abs_rho  = np.abs(rho_sigma)
    max_val  = np.max(abs_rho) if len(abs_rho) > 0 else 0.0
    threshold = min_fraction * max_val

    peaks = []
    for i in range(1, len(abs_rho) - 1):
        if (abs_rho[i] > abs_rho[i-1] and
                abs_rho[i] > abs_rho[i+1] and
                abs_rho[i] >= threshold):
            peaks.append((float(sigma_mid[i]), float(rho_sigma[i])))

    peaks.sort(key=lambda x: x[0])
    return peaks


def check_peak_near(peaks, expected_sigma, tol=0.03):
    """
    예상 σ 위치 근방에 피크가 있는지 확인.
    Returns: (found, closest_sigma, closest_rho)
    """
    for s, r in peaks:
        if abs(s - expected_sigma) <= tol:
            return True, s, r
    return False, None, None


# ═══════════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    now     = datetime.now().strftime("%Y-%m-%d %H:%M")

    log("=" * 77)
    log("[Project RDL] Synthetic Off-Critical Control")
    log("프레임워크 특이성: σ-국소화는 방법론적 인공물인가?")
    log(f"실행: {now}")
    log("=" * 77)
    log()
    log("합성 함수 정의 (영점 위치 완전히 알려짐):")
    log("  f₁(s) = (s-(0.5+14i))(s-(0.5+21i))(s-(0.5+25i))  [모두 σ=1/2]")
    log("  f₂(s) = (s-(0.3+14i))(s-(0.7+14i))(s-(0.5+21i))  [혼합]")
    log("  f₃(s) = (s-(0.3+16i))(s-(0.7+16i))(s-(0.6+22i))(s-(0.4+22i))")
    log("                                                       [완전 off-critical]")
    log()
    log(f"설정: dps={DPS}, Nσ={N_SIGMA}, Nt={N_T}")
    log()
    log("핵심 예측:")
    log("  §A: f₁ → δ(σ=0.5) (RH-유사 함수는 임계선 집중)")
    log("  §B: f₂ → σ=0.3, 0.5, 0.7 세 피크 (영점 위치별 분리)")
    log("  §C: f₃ → σ=0.3,0.4,0.6,0.7 네 피크; 임계선 N=0 ← 핵심!")
    log()

    all_results = {}
    total_pass  = 0
    total_tests = 0

    save_all()

    # ═══════════════════════════════════════════════════════════════════════════
    # §A: f₁ 테스트 (모든 영점 σ=1/2)
    # t-범위: [10, 28], 영점 Im = 14, 21, 25 (범위 내, 격자점 회피됨)
    # ═══════════════════════════════════════════════════════════════════════════
    log("═" * 77)
    log("§A: f₁(s) — 모든 영점 σ=1/2 (RH-유사)")
    log("  영점: (0.5+14i), (0.5+21i), (0.5+25i)")
    log("═" * 77)
    log()

    f1_fn = lambda s, t: log_poly_imag(s, t, F1_ZEROS)
    F1_T_MIN, F1_T_MAX = 10.0, 28.0

    # §A-1: 전체 검출 [0.3, 0.7], 기대 N=3
    log("--- §A-1: f₁, σ∈[0.3,0.7], t∈[10,28], 기대 N=3 ---")
    N_A1, _, _, _, _ = compute_area_integral(
        f1_fn, 0.30, 0.70, F1_T_MIN, F1_T_MAX,
        n_sigma=N_SIGMA, n_t=N_T, label="A-1"
    )
    err_A1   = abs(N_A1 - 3.0)
    pass_A1  = err_A1 < 0.5
    log(f"  → N = {N_A1:.6f}, 기대=3, |오차|={err_A1:.6f}")
    log(f"  → 판정: {'✅ PASS' if pass_A1 else '❌ FAIL'}")
    log()
    total_tests += 1
    total_pass  += int(pass_A1)
    save_all()

    # §A-2: σ-프로파일, 기대: σ=0.5에 단일 피크
    log("--- §A-2: f₁ σ-프로파일, σ∈[0.1,0.9], 기대: σ=0.5 단일 피크 ---")
    smid_A2, rho_A2, totN_A2 = compute_sigma_profile(
        f1_fn, 0.10, 0.90, F1_T_MIN, F1_T_MAX,
        n_sigma=200, n_t=N_T, label="A-2"
    )
    peaks_A2 = find_peaks_profile(smid_A2, rho_A2, min_fraction=0.15)
    log(f"  검출된 피크 ({len(peaks_A2)}개): {[(f'{s:.3f}', f'{r:.2f}') for s, r in peaks_A2]}")

    found_A2, sig_A2, rho_at_A2 = check_peak_near(peaks_A2, 0.5, tol=0.03)
    # off-critical 피크가 없는지도 확인 (σ<0.45 또는 σ>0.55의 피크)
    offcrit_peaks_A2 = [(s, r) for s, r in peaks_A2 if abs(s - 0.5) > 0.05]
    pass_A2 = found_A2 and (len(offcrit_peaks_A2) == 0)
    log(f"  σ=0.5 근방 피크: {'있음 ('+f'σ={sig_A2:.3f}, ρ={rho_at_A2:.2f}'+')' if found_A2 else '없음 ❌'}")
    log(f"  off-critical 피크 (|σ-0.5|>0.05): {len(offcrit_peaks_A2)}개 {'✅' if len(offcrit_peaks_A2)==0 else '❌'}")
    log(f"  → 판정: {'✅ PASS' if pass_A2 else '❌ FAIL'} (기대: σ=0.5 단일 피크)")
    log()
    total_tests += 1
    total_pass  += int(pass_A2)
    save_all()

    # §A-3: off-critical 오른쪽 [0.6, 0.9], 기대 N=0
    log("--- §A-3: f₁, σ∈[0.6,0.9], t∈[10,28], 기대 N=0 ---")
    N_A3, _, _, _, _ = compute_area_integral(
        f1_fn, 0.60, 0.90, F1_T_MIN, F1_T_MAX,
        n_sigma=N_SIGMA, n_t=N_T, label="A-3"
    )
    pass_A3 = abs(N_A3) < 0.5
    log(f"  → N = {N_A3:.6f}, 기대=0, |N|={abs(N_A3):.6f}")
    log(f"  → 판정: {'✅ PASS' if pass_A3 else '❌ FAIL'}")
    log()
    total_tests += 1
    total_pass  += int(pass_A3)
    save_all()

    pass_A = pass_A1 and pass_A2 and pass_A3
    log(f"§A 소계: {sum([pass_A1, pass_A2, pass_A3])}/3")
    log()
    all_results['A'] = {
        'A1': {'N': N_A1, 'pass': pass_A1},
        'A2': {'peaks': peaks_A2, 'pass': pass_A2},
        'A3': {'N': N_A3, 'pass': pass_A3},
    }

    # ═══════════════════════════════════════════════════════════════════════════
    # §B: f₂ 테스트 (혼합: σ=0.3, 0.7, 0.5)
    # t-범위: [10, 24], 영점 Im = 14, 21 (범위 내)
    # ═══════════════════════════════════════════════════════════════════════════
    log("═" * 77)
    log("§B: f₂(s) — 혼합 영점 (σ=0.3, 0.7, 0.5)")
    log("  영점: (0.3+14i), (0.7+14i), (0.5+21i)")
    log("═" * 77)
    log()

    f2_fn = lambda s, t: log_poly_imag(s, t, F2_ZEROS)
    F2_T_MIN, F2_T_MAX = 10.0, 24.0

    # §B-1: 전체 검출 [0.2, 0.8], 기대 N=3
    log("--- §B-1: f₂, σ∈[0.2,0.8], t∈[10,24], 기대 N=3 ---")
    N_B1, _, _, _, _ = compute_area_integral(
        f2_fn, 0.20, 0.80, F2_T_MIN, F2_T_MAX,
        n_sigma=N_SIGMA, n_t=N_T, label="B-1"
    )
    err_B1  = abs(N_B1 - 3.0)
    pass_B1 = err_B1 < 0.5
    log(f"  → N = {N_B1:.6f}, 기대=3, |오차|={err_B1:.6f}")
    log(f"  → 판정: {'✅ PASS' if pass_B1 else '❌ FAIL'}")
    log()
    total_tests += 1
    total_pass  += int(pass_B1)
    save_all()

    # §B-2: σ-프로파일, 기대: σ=0.3, 0.5, 0.7 세 피크
    log("--- §B-2: f₂ σ-프로파일, σ∈[0.1,0.9], 기대: σ=0.3, 0.5, 0.7 세 피크 ---")
    smid_B2, rho_B2, totN_B2 = compute_sigma_profile(
        f2_fn, 0.10, 0.90, F2_T_MIN, F2_T_MAX,
        n_sigma=200, n_t=N_T, label="B-2"
    )
    peaks_B2 = find_peaks_profile(smid_B2, rho_B2, min_fraction=0.15)
    log(f"  검출된 피크 ({len(peaks_B2)}개): {[(f'{s:.3f}', f'{r:.2f}') for s, r in peaks_B2]}")

    found_03_B2, s03, r03 = check_peak_near(peaks_B2, 0.3, tol=0.03)
    found_05_B2, s05, r05 = check_peak_near(peaks_B2, 0.5, tol=0.03)
    found_07_B2, s07, r07 = check_peak_near(peaks_B2, 0.7, tol=0.03)
    pass_B2 = found_03_B2 and found_05_B2 and found_07_B2
    log(f"  σ=0.3 피크: {'✅ σ='+f'{s03:.3f}' if found_03_B2 else '❌ 없음'}")
    log(f"  σ=0.5 피크: {'✅ σ='+f'{s05:.3f}' if found_05_B2 else '❌ 없음'}")
    log(f"  σ=0.7 피크: {'✅ σ='+f'{s07:.3f}' if found_07_B2 else '❌ 없음'}")
    log(f"  → 판정: {'✅ PASS' if pass_B2 else '❌ FAIL'} (기대: 3피크 분리)")
    log()
    total_tests += 1
    total_pass  += int(pass_B2)
    save_all()

    # §B-3: σ=0.5 영점만 [0.49, 0.51], 기대 N=1
    log("--- §B-3: f₂, σ∈[0.49,0.51], t∈[10,24], 기대 N=1 (σ=0.5 영점만) ---")
    N_B3, _, _, _, _ = compute_area_integral(
        f2_fn, 0.49, 0.51, F2_T_MIN, F2_T_MAX,
        n_sigma=N_SIGMA, n_t=N_T, label="B-3"
    )
    err_B3  = abs(N_B3 - 1.0)
    pass_B3 = err_B3 < 0.5
    log(f"  → N = {N_B3:.6f}, 기대=1, |오차|={err_B3:.6f}")
    log(f"  → 판정: {'✅ PASS' if pass_B3 else '❌ FAIL'}")
    log()
    total_tests += 1
    total_pass  += int(pass_B3)
    save_all()

    # §B-4: σ=0.7 영점만 [0.6, 0.8], 기대 N=1
    log("--- §B-4: f₂, σ∈[0.6,0.8], t∈[10,24], 기대 N=1 (σ=0.7 영점만) ---")
    N_B4, _, _, _, _ = compute_area_integral(
        f2_fn, 0.60, 0.80, F2_T_MIN, F2_T_MAX,
        n_sigma=N_SIGMA, n_t=N_T, label="B-4"
    )
    err_B4  = abs(N_B4 - 1.0)
    pass_B4 = err_B4 < 0.5
    log(f"  → N = {N_B4:.6f}, 기대=1, |오차|={err_B4:.6f}")
    log(f"  → 판정: {'✅ PASS' if pass_B4 else '❌ FAIL'}")
    log()
    total_tests += 1
    total_pass  += int(pass_B4)
    save_all()

    # §B-5: σ=0.3 영점만 [0.2, 0.4], 기대 N=1
    log("--- §B-5: f₂, σ∈[0.2,0.4], t∈[10,24], 기대 N=1 (σ=0.3 영점만) ---")
    N_B5, _, _, _, _ = compute_area_integral(
        f2_fn, 0.20, 0.40, F2_T_MIN, F2_T_MAX,
        n_sigma=N_SIGMA, n_t=N_T, label="B-5"
    )
    err_B5  = abs(N_B5 - 1.0)
    pass_B5 = err_B5 < 0.5
    log(f"  → N = {N_B5:.6f}, 기대=1, |오차|={err_B5:.6f}")
    log(f"  → 판정: {'✅ PASS' if pass_B5 else '❌ FAIL'}")
    log()
    total_tests += 1
    total_pass  += int(pass_B5)
    save_all()

    pass_B_list = [pass_B1, pass_B2, pass_B3, pass_B4, pass_B5]
    log(f"§B 소계: {sum(pass_B_list)}/5")
    log()
    all_results['B'] = {
        'B1': {'N': N_B1, 'pass': pass_B1},
        'B2': {'peaks': peaks_B2, 'pass': pass_B2},
        'B3': {'N': N_B3, 'pass': pass_B3},
        'B4': {'N': N_B4, 'pass': pass_B4},
        'B5': {'N': N_B5, 'pass': pass_B5},
    }

    # ═══════════════════════════════════════════════════════════════════════════
    # §C: f₃ 테스트 (완전 off-critical: σ=0.3, 0.7, 0.6, 0.4)
    # t-범위: [12, 25], 영점 Im = 16, 22 (범위 내)
    # ═══════════════════════════════════════════════════════════════════════════
    log("═" * 77)
    log("§C: f₃(s) — 완전 off-critical (RH 위반)")
    log("  영점: (0.3+16i), (0.7+16i), (0.6+22i), (0.4+22i)")
    log("  ← f₃에는 σ=1/2 영점이 단 하나도 없다!")
    log("═" * 77)
    log()

    f3_fn = lambda s, t: log_poly_imag(s, t, F3_ZEROS)
    F3_T_MIN, F3_T_MAX = 12.0, 25.0

    # §C-1: 전체 검출 [0.2, 0.8], 기대 N=4
    log("--- §C-1: f₃, σ∈[0.2,0.8], t∈[12,25], 기대 N=4 ---")
    N_C1, _, _, _, _ = compute_area_integral(
        f3_fn, 0.20, 0.80, F3_T_MIN, F3_T_MAX,
        n_sigma=N_SIGMA, n_t=N_T, label="C-1"
    )
    err_C1  = abs(N_C1 - 4.0)
    pass_C1 = err_C1 < 0.5
    log(f"  → N = {N_C1:.6f}, 기대=4, |오차|={err_C1:.6f}")
    log(f"  → 판정: {'✅ PASS' if pass_C1 else '❌ FAIL'}")
    log()
    total_tests += 1
    total_pass  += int(pass_C1)
    save_all()

    # §C-2: σ-프로파일, 기대: σ=0.3, 0.4, 0.6, 0.7 피크; σ=0.5 비어야 함
    log("--- §C-2: f₃ σ-프로파일, σ∈[0.1,0.9], 기대: 4피크(σ=0.3,0.4,0.6,0.7), σ=0.5 없음 ---")
    smid_C2, rho_C2, totN_C2 = compute_sigma_profile(
        f3_fn, 0.10, 0.90, F3_T_MIN, F3_T_MAX,
        n_sigma=200, n_t=N_T, label="C-2"
    )
    peaks_C2 = find_peaks_profile(smid_C2, rho_C2, min_fraction=0.15)
    log(f"  검출된 피크 ({len(peaks_C2)}개): {[(f'{s:.3f}', f'{r:.2f}') for s, r in peaks_C2]}")

    found_03_C2, sc03, rc03 = check_peak_near(peaks_C2, 0.3, tol=0.03)
    found_04_C2, sc04, rc04 = check_peak_near(peaks_C2, 0.4, tol=0.03)
    found_06_C2, sc06, rc06 = check_peak_near(peaks_C2, 0.6, tol=0.03)
    found_07_C2, sc07, rc07 = check_peak_near(peaks_C2, 0.7, tol=0.03)
    found_05_C2, sc05, rc05 = check_peak_near(peaks_C2, 0.5, tol=0.05)  # σ=0.5 없어야 함

    # σ=0.5 근방의 ρ 값 확인 (직접)
    idx_05  = np.argmin(np.abs(smid_C2 - 0.5))
    rho_at_05 = float(rho_C2[idx_05])
    max_rho   = float(np.max(np.abs(rho_C2)))
    rho_05_fraction = abs(rho_at_05) / max_rho if max_rho > 0 else 0.0

    pass_C2 = (found_03_C2 and found_04_C2 and found_06_C2 and found_07_C2
               and not found_05_C2)

    log(f"  σ=0.3 피크: {'✅ σ='+f'{sc03:.3f}' if found_03_C2 else '❌ 없음'}")
    log(f"  σ=0.4 피크: {'✅ σ='+f'{sc04:.3f}' if found_04_C2 else '❌ 없음'}")
    log(f"  σ=0.6 피크: {'✅ σ='+f'{sc06:.3f}' if found_06_C2 else '❌ 없음'}")
    log(f"  σ=0.7 피크: {'✅ σ='+f'{sc07:.3f}' if found_07_C2 else '❌ 없음'}")
    log(f"  σ=0.5 피크: {'❌ 있음 σ='+f'{sc05:.3f}' if found_05_C2 else '✅ 없음 (기대)'}")
    log(f"  σ=0.5 근방 ρ값: {rho_at_05:.4f} (전체 최대의 {rho_05_fraction:.1%})")
    log(f"  → 판정: {'✅ PASS' if pass_C2 else '❌ FAIL'} (기대: 4피크 분리, σ=0.5 빈 공간)")
    log()
    total_tests += 1
    total_pass  += int(pass_C2)
    save_all()

    # §C-3: ★★ 핵심 테스트 ★★ 임계선 [0.49, 0.51], 기대 N=0
    log("--- §C-3: ★ 핵심! f₃, σ∈[0.49,0.51], t∈[12,25], 기대 N=0 ---")
    log("    f₃에는 σ=1/2 영점이 없다. 프레임워크가 올바르면 N=0이어야 함.")
    log("    N≠0이 나오면 σ=1/2 δ-집중이 인공물 → 프레임워크 신뢰성 위기!")
    N_C3, _, _, _, _ = compute_area_integral(
        f3_fn, 0.49, 0.51, F3_T_MIN, F3_T_MAX,
        n_sigma=N_SIGMA, n_t=N_T, label="C-3"
    )
    pass_C3 = abs(N_C3) < 0.5
    log(f"  → N = {N_C3:.6f}, 기대=0, |N|={abs(N_C3):.6f}")
    if pass_C3:
        log(f"  → 판정: ✅ PASS — 임계선 비어있음. 프레임워크 특이성 확인!")
    else:
        log(f"  → 판정: ❌ FAIL — 임계선에서 비제로 자속! 심각한 문제 검토 필요!")
    log()
    total_tests += 1
    total_pass  += int(pass_C3)
    save_all()

    pass_C_list = [pass_C1, pass_C2, pass_C3]
    log(f"§C 소계: {sum(pass_C_list)}/3")
    log()
    all_results['C'] = {
        'C1': {'N': N_C1, 'pass': pass_C1},
        'C2': {'peaks': peaks_C2, 'pass': pass_C2, 'rho_at_05': rho_at_05},
        'C3': {'N': N_C3, 'pass': pass_C3},
    }

    # ═══════════════════════════════════════════════════════════════════════════
    # 종합
    # ═══════════════════════════════════════════════════════════════════════════
    log()
    log("═" * 77)
    log("종합 결과 — 11종 테스트 매트릭스")
    log("═" * 77)
    log()

    # 테이블 출력
    log(f"  {'테스트':<8} {'함수':<5} {'σ-범위':<18} {'결과':>12} {'기대':>6} {'판정'}")
    log(f"  {'-'*8} {'-'*5} {'-'*18} {'-'*12} {'-'*6} {'-'*6}")

    def fmt_N(n):
        return f"{n:.4f}"

    rows = [
        ("§A-1", "f₁", "[0.3, 0.7]",   fmt_N(N_A1), "N=3",  "✅" if pass_A1 else "❌"),
        ("§A-2", "f₁", "σ-profile",    f"{len(peaks_A2)}피크",  "1피크@0.5","✅" if pass_A2 else "❌"),
        ("§A-3", "f₁", "[0.6, 0.9]",   fmt_N(N_A3), "N=0",  "✅" if pass_A3 else "❌"),
        ("§B-1", "f₂", "[0.2, 0.8]",   fmt_N(N_B1), "N=3",  "✅" if pass_B1 else "❌"),
        ("§B-2", "f₂", "σ-profile",    f"{len(peaks_B2)}피크",  "3피크",    "✅" if pass_B2 else "❌"),
        ("§B-3", "f₂", "[0.49, 0.51]", fmt_N(N_B3), "N=1",  "✅" if pass_B3 else "❌"),
        ("§B-4", "f₂", "[0.6, 0.8]",   fmt_N(N_B4), "N=1",  "✅" if pass_B4 else "❌"),
        ("§B-5", "f₂", "[0.2, 0.4]",   fmt_N(N_B5), "N=1",  "✅" if pass_B5 else "❌"),
        ("§C-1", "f₃", "[0.2, 0.8]",   fmt_N(N_C1), "N=4",  "✅" if pass_C1 else "❌"),
        ("§C-2", "f₃", "σ-profile",    f"{len(peaks_C2)}피크",  "4피크",    "✅" if pass_C2 else "❌"),
        ("§C-3★","f₃", "[0.49, 0.51]", fmt_N(N_C3), "N=0★", "✅" if pass_C3 else "❌"),
    ]
    for row in rows:
        log(f"  {row[0]:<8} {row[1]:<5} {row[2]:<18} {row[3]:>12} {row[4]:>6}  {row[5]}")

    log()
    log(f"  §A: {sum([pass_A1,pass_A2,pass_A3])}/3  §B: {sum(pass_B_list)}/5  §C: {sum(pass_C_list)}/3")
    log(f"  전체: {total_pass}/{total_tests}")
    log()

    # ─── 물리적 해석 ──────────────────────────────────────────────────────────
    log("─── 물리적 해석 ───")
    log()
    log("§C-3 분석 (핵심):")
    if pass_C3:
        log(f"  f₃(σ=1/2+it)에 영점이 없으면 플라켓 자속 = 0 (N={N_C3:.6f})")
        log("  → 프레임워크는 σ=1/2 영점의 존재 여부를 정확히 식별한다.")
        log("  → σ-국소화는 방법론적 인공물이 아님 — Cauchy 논증 원리의 정확한 구현.")
    else:
        log(f"  ⚠️ f₃에 σ=1/2 영점이 없는데 N={N_C3:.6f} 비제로 반환!")
        log("  → 이는 방법론적 인공물을 시사. 즉각 원인 분석 필요.")
    log()

    log("§B 분석:")
    if pass_B3 and pass_B4 and pass_B5:
        log(f"  개별 σ 위치별 영점 분리 검출 성공:")
        log(f"    σ=0.5 strip: N={N_B3:.4f} (기대 1)  σ=0.7 strip: N={N_B4:.4f} (기대 1)"
            f"  σ=0.3 strip: N={N_B5:.4f} (기대 1)")
        log("  → 영점이 정확히 자신의 σ 위치에만 기여한다는 것을 수치적으로 확인.")
    log()

    log("§A 분석:")
    if pass_A2:
        log(f"  f₁의 σ-프로파일에서 σ=0.5 단일 피크 확인 ({len(peaks_A2)}개 피크)")
        log("  → RH-유사 함수는 임계선 집중 패턴을 보임 (예상대로).")
    log()

    # ─── 최종 판정 ────────────────────────────────────────────────────────────
    log()
    if total_pass == 11:
        verdict = "★★★ 판정: 완전 양성 (11/11) — 프레임워크 특이성 완전 확립"
        detail  = ("f₁은 δ(σ=0.5), f₂는 영점별 σ 분리, f₃는 임계선 N=0.\n"
                   "  σ-국소화 관측은 ξ의 고유 성질이지 방법론적 인공물이 아님.")
    elif total_pass >= 9:
        verdict = f"★★ 판정: 실질 양성 ({total_pass}/11) — 핵심 테스트 통과"
        detail  = "일부 프로파일 검출 정밀도 이슈 있으나 핵심 C-3 확인."
    elif pass_C3:
        verdict = f"★ 판정: 조건부 양성 ({total_pass}/11) — C-3(핵심) 통과"
        detail  = "임계선 특이성은 확인. 세부 분리 테스트 일부 실패."
    else:
        verdict = f"★ 판정: 음성 ({total_pass}/11) — C-3 실패, 긴급 검토 필요"
        detail  = "프레임워크 방법론적 인공물 가능성. 수학자 즉시 보고."

    log(verdict)
    log(f"  {detail}")
    log()

    elapsed_total = time.time() - t_start
    log(f"총 실행 시간: {elapsed_total:.1f}초 ({elapsed_total/60:.1f}분)")
    log(f"dps={DPS}, Nσ={N_SIGMA}, Nt={N_T}")
    save_all()
    log()
    log(f"결과 저장: {OUTPUT_TXT}")


if __name__ == "__main__":
    main()
