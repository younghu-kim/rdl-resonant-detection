"""
=============================================================================
[Project RDL] Off-critical Gauss-Bonnet — 위상적 전하의 σ-국소화
=============================================================================
사이클 35 / 수학자 지시 2026-04-14 11:15

목적:
  Gauss-Bonnet(#24)이 "영점 포함 영역의 면적분 = 2πN"을 확인했다.
  모든 테스트가 σ∈[0.3,0.7]으로 임계선을 걸치는 영역이었다.
  이번 질문: **위상적 전하가 정확히 σ=1/2에 집중되어 있는가?**

수학적 의미:
  RH가 참이면 모든 영점이 σ=1/2 위에 있으므로:
  - σ=1/2를 포함하지 않는 영역 → ∫∫F₂ = 0 (영점 없음)
  - σ=1/2를 아무리 좁게 감싸도 → ∫∫F₂ = 2πN (모든 영점 포함)
  이것은 RH의 기하학적 동치:
  "ξ-다발의 위상적 전하가 σ=1/2에 δ-함수적으로 집중"

설계:
  §A — σ-strip 변조: 7가지 σ 구간에서 t∈[0.5,78.24] (20개 영점)
  §B — σ-프로파일: ρ(σ) = σ열의 플라켓 자속 합, σ∈[0.1,0.9] 200등분
  §C — 디리클레 확장: χ₅에서 좁은 strip + off-critical

출력:
  results/offcritical_gauss_bonnet.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT_TXT = os.path.join(BASE_DIR, "results", "offcritical_gauss_bonnet.txt")
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# ─── 설정 ─────────────────────────────────────────────────────────────────────
DPS = 100
mpmath.mp.dps = DPS

N_SIGMA_BASE = 64   # σ 방향 기본 격자 점 수 (A3 극좁은 strip은 증가)
T_MIN_GLOBAL = 0.5
T_MAX_GLOBAL = 78.24  # ζ 20번째 영점(γ₂₀=77.1445) 이후 중간점 근방

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
# 위상 함수
# ═══════════════════════════════════════════════════════════════════════════════

def log_xi_imag(sigma, t):
    """θ(σ,t) = Im(log ξ(σ+it)) for ζ"""
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
    """θ(σ,t) = Im(log Λ(σ+it, χ₅)) for χ₅=(·/5), 짝수 (a=0)"""
    s = mpmath.mpc(sigma, t)
    L_val = mpmath.dirichlet(s, CHI5_LIST)
    if abs(L_val) < mpmath.mpf(10) ** (-DPS + 15):
        return None
    lL = ((s / 2) * mpmath.log(mpmath.mpf(CHI5_Q) / mpmath.pi)
          + mpmath.loggamma(s / 2)
          + mpmath.log(L_val))
    return float(lL.imag)


# ═══════════════════════════════════════════════════════════════════════════════
# 격자 플라켓 곡률 (gauss_bonnet_area_integral.py에서 재사용)
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
                          n_sigma=64, n_t=None, label=""):
    """
    2D 영역에서 격자 플라켓 곡률의 면적분 계산.
    Returns: total_flux, flux_grid, sigma_arr, t_arr
    """
    height = t_max - t_min
    if n_t is None:
        n_t = max(64, int(height * 8))  # 단위당 최소 8개 점

    sigma_arr = np.linspace(sigma_min, sigma_max, n_sigma)
    t_arr     = np.linspace(t_min, t_max, n_t)

    log(f"  격자: Nσ={n_sigma}, Nt={n_t}, 총 {n_sigma*n_t}점")
    log(f"  σ∈[{sigma_min},{sigma_max}], t∈[{t_min:.6f},{t_max:.4f}]")

    theta_grid = np.zeros((n_sigma, n_t))
    skip_count = 0

    t0 = time.time()
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

    return total_flux, flux_grid, sigma_arr, t_arr, theta_grid


# ═══════════════════════════════════════════════════════════════════════════════
# χ₅ 영점 탐색
# ═══════════════════════════════════════════════════════════════════════════════

def find_zeros_chi5(t_min=1.0, t_max=60.0, n_scan=6000):
    """χ₅=(·/5) 영점 탐색: |L(1/2+it)| 최소화 방식"""
    ts = np.linspace(t_min, t_max, n_scan)
    half = mpmath.mpf('0.5')

    abs_vals = []
    for t in ts:
        s = half + 1j * mpmath.mpf(str(t))
        val = mpmath.dirichlet(s, CHI5_LIST)
        abs_vals.append(float(abs(val)))
    abs_vals = np.array(abs_vals)

    zeros = []
    fail_count = 0
    for i in range(1, len(abs_vals) - 1):
        if abs_vals[i] < abs_vals[i-1] and abs_vals[i] < abs_vals[i+1]:
            if abs_vals[i] < 0.1:
                t_guess = ts[i]
                found = False
                for f_name, f_fn in [
                    ('Re', lambda tv: mpmath.re(mpmath.dirichlet(
                        half + 1j*mpmath.mpf(tv), CHI5_LIST))),
                    ('Im', lambda tv: mpmath.im(mpmath.dirichlet(
                        half + 1j*mpmath.mpf(tv), CHI5_LIST))),
                ]:
                    try:
                        t_ref = float(mpmath.findroot(f_fn, t_guess))
                        sv = half + 1j * mpmath.mpf(str(t_ref))
                        if abs(mpmath.dirichlet(sv, CHI5_LIST)) < 1e-6:
                            if not zeros or abs(t_ref - zeros[-1]) > 0.1:
                                zeros.append(t_ref)
                                found = True
                                break
                    except Exception as e:
                        pass
                if not found:
                    fail_count += 1

    zeros.sort()
    log(f"  χ₅ 영점 {len(zeros)}개 발견, findroot 실패 {fail_count}회")
    if len(zeros) > 0:
        log(f"  첫 5개: {[f'{z:.4f}' for z in zeros[:5]]}")
    return zeros


# ═══════════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    log("=" * 77)
    log("[Project RDL] Off-critical Gauss-Bonnet — σ-국소화 검증")
    log("위상적 전하가 σ=1/2에 δ-함수적으로 집중되어 있는가?")
    log(f"실행: {now}")
    log("=" * 77)
    log()
    log("RH의 기하학적 동치 테스트:")
    log("  σ=1/2를 포함하지 않는 strip → ∫∫F₂ = 0 (오프-임계, 영점 없음)")
    log("  σ=1/2를 극히 좁게 감싸도  → ∫∫F₂ = 2πN (임계선 포함, 전부 포착)")
    log()
    log(f"설정: dps={DPS}, t∈[{T_MIN_GLOBAL},{T_MAX_GLOBAL}] (ζ 영점 20개)")
    log()

    # ζ 영점 준비 (20개)
    log("ζ 영점 준비...")
    zeta_zeros_20 = [float(mpmath.zetazero(n).imag) for n in range(1, 21)]
    log(f"  γ₁={zeta_zeros_20[0]:.4f}, γ₂₀={zeta_zeros_20[19]:.4f}")
    # t_max: γ₂₀과 γ₂₁의 중간점
    gamma21 = float(mpmath.zetazero(21).imag)
    T_MAX = (zeta_zeros_20[19] + gamma21) / 2.0
    log(f"  T_MAX = (γ₂₀+γ₂₁)/2 = {T_MAX:.4f}")
    log()

    results = {}
    save_all()

    # ═══════════════════════════════════════════════════════════════════════════
    # §A: σ-strip 변조 (핵심 테스트)
    # ═══════════════════════════════════════════════════════════════════════════
    log("═" * 77)
    log("§A: σ-strip 변조 — 위상적 전하의 σ 위치 의존성")
    log("═" * 77)
    log()
    log("테스트 설계:")
    log("  A1 [0.3, 0.7]    기대=20  baseline (#24 재현)")
    log("  A2 [0.49, 0.51]  기대=20  좁은 strip (폭 0.02)")
    log("  A3 [0.499, 0.501] 기대=20  극좁은 strip (폭 0.002)")
    log("  A4 [0.6, 0.9]    기대=0   off-critical 오른쪽")
    log("  A5 [0.1, 0.4]    기대=0   off-critical 왼쪽")
    log("  A6 [0.3, 0.499]  기대=0   반-strip (임계선 미포함)")
    log("  A7 [0.501, 0.7]  기대=0   반-strip (임계선 미포함)")
    log()

    strip_tests = [
        ("A1", 0.30,  0.70,  20, "baseline"),
        ("A2", 0.49,  0.51,  20, "좁은 strip (폭 0.02)"),
        ("A3", 0.499, 0.501, 20, "극좁은 strip (폭 0.002, Nσ=128)"),
        ("A4", 0.60,  0.90,   0, "off-critical 오른쪽"),
        ("A5", 0.10,  0.40,   0, "off-critical 왼쪽"),
        ("A6", 0.30,  0.499,  0, "반-strip 왼쪽 (0.5 미포함)"),
        ("A7", 0.501, 0.70,   0, "반-strip 오른쪽 (0.5 미포함)"),
    ]

    results_A = []
    for test_id, sig_lo, sig_hi, expected_N, desc in strip_tests:
        log(f"--- {test_id}: σ∈[{sig_lo},{sig_hi}] 기대={expected_N} ({desc}) ---")

        # A3는 극좁은 strip → Nσ=128으로 증가
        n_sig = 128 if test_id == "A3" else N_SIGMA_BASE

        flux, fgrid, sig_arr, t_arr, theta = compute_area_integral(
            log_xi_imag, sig_lo, sig_hi, T_MIN_GLOBAL, T_MAX,
            n_sigma=n_sig, label=test_id
        )

        N_calc = flux / (2 * np.pi)

        if expected_N == 0:
            err = abs(N_calc)
            pass_crit = err < 0.5
            log(f"  → N = {N_calc:.6f}, 기대=0, |N|={err:.6f}")
            log(f"  → 판정: {'✅ PASS' if pass_crit else '❌ FAIL'} (기준: |N| < 0.5)")
        else:
            err = abs(N_calc - expected_N)
            pass_crit = err < 1.0  # 기대 N=20에서 1개 이내
            log(f"  → N = {N_calc:.6f}, 기대={expected_N}, |오차|={err:.6f}")
            log(f"  → 판정: {'✅ PASS' if pass_crit else '❌ FAIL'} (기준: |오차| < 1.0)")

        results_A.append({
            'id': test_id, 'sig_lo': sig_lo, 'sig_hi': sig_hi,
            'N': N_calc, 'expected': expected_N, 'err': err,
            'pass': pass_crit, 'desc': desc
        })
        log()
        save_all()

    results['A'] = results_A

    # ═══════════════════════════════════════════════════════════════════════════
    # §B: σ-프로파일 ρ(σ)
    # ═══════════════════════════════════════════════════════════════════════════
    log("═" * 77)
    log("§B: σ-프로파일 — ρ(σ) = σ열의 플라켓 자속 합")
    log("σ∈[0.1,0.9] 200등분, σ=1/2 근방에 sharp peak 예상")
    log("═" * 77)
    log()

    # 전략: 전체 격자 σ∈[0.1,0.9], t∈[T_MIN_GLOBAL, T_MAX]를 한 번에 계산 후
    # 각 σ 열(σ 방향 슬라이스)의 자속 합을 집계
    N_SIGMA_B = 200  # σ 방향 200점
    N_T_B = max(128, int((T_MAX - T_MIN_GLOBAL) * 6))  # 적당한 해상도
    SIG_MIN_B = 0.1
    SIG_MAX_B = 0.9

    log(f"  전체 격자: Nσ={N_SIGMA_B}, Nt={N_T_B}, σ∈[{SIG_MIN_B},{SIG_MAX_B}]")
    log(f"  총 {N_SIGMA_B * N_T_B}점 계산")
    log()

    sigma_B = np.linspace(SIG_MIN_B, SIG_MAX_B, N_SIGMA_B)
    t_B     = np.linspace(T_MIN_GLOBAL, T_MAX, N_T_B)

    log("  위상 격자 계산 중...")
    theta_B = np.zeros((N_SIGMA_B, N_T_B))
    skip_B = 0
    t0_B = time.time()

    for i, sig in enumerate(sigma_B):
        for j, t in enumerate(t_B):
            val = log_xi_imag(float(sig), float(t))
            if val is None:
                skip_B += 1
                if j > 0:
                    theta_B[i, j] = theta_B[i, j-1]
                elif i > 0:
                    theta_B[i, j] = theta_B[i-1, j]
                else:
                    theta_B[i, j] = 0.0
            else:
                theta_B[i, j] = val
        if (i + 1) % 20 == 0:
            log(f"    σ 진행: {i+1}/{N_SIGMA_B} ({time.time()-t0_B:.1f}초)")
            save_all()

    elapsed_B = time.time() - t0_B
    log(f"  완료: {elapsed_B:.1f}초, skip={skip_B}")

    # 플라켓 자속 계산
    flux_B_full, total_B = compute_plaquette_flux(theta_B)
    log(f"  전체 자속: {total_B:.4f}, N={total_B/(2*np.pi):.4f}")
    log()

    # σ-프로파일: 각 σ 열의 자속 합 (flux_B_full shape: (Nσ-1, Nt-1))
    # 각 σ 인덱스 i (0..N_SIGMA_B-2)의 자속 = flux_B_full[i, :].sum()
    rho_sigma = np.sum(flux_B_full, axis=1)  # (N_SIGMA_B-1,)
    sigma_mid = (sigma_B[:-1] + sigma_B[1:]) / 2.0  # 플라켓 중심 σ

    # σ=0.5 근방 peak 분석
    idx_peak = np.argmax(np.abs(rho_sigma))
    sigma_peak = sigma_mid[idx_peak]
    rho_peak = rho_sigma[idx_peak]
    log(f"  ρ(σ) 최대 위치: σ={sigma_peak:.4f}, ρ={rho_peak:.4f}")

    # FWHM 계산 (반최대폭)
    half_max = abs(rho_peak) / 2.0
    above_half = np.where(np.abs(rho_sigma) >= half_max)[0]
    if len(above_half) >= 2:
        fwhm = sigma_mid[above_half[-1]] - sigma_mid[above_half[0]]
    else:
        fwhm = 0.0
    log(f"  FWHM (반최대폭): {fwhm:.4f}")
    log(f"  격자 해상도 δσ: {sigma_B[1]-sigma_B[0]:.4f}")
    log(f"  FWHM/δσ ≈ {fwhm/(sigma_B[1]-sigma_B[0]):.1f}")

    pass_B = fwhm < 0.01
    log(f"  → 판정: {'✅ PASS' if pass_B else '❌ FAIL'} (기준: FWHM < 0.01)")
    log()

    # 대표값 출력 (10점마다)
    log("  ρ(σ) 프로파일 (10점 간격):")
    log(f"  {'σ':>8}  {'ρ(σ)':>12}")
    for k in range(0, len(sigma_mid), max(1, len(sigma_mid)//20)):
        log(f"  {sigma_mid[k]:8.4f}  {rho_sigma[k]:12.4f}")

    # CSV 저장
    csv_path = os.path.join(BASE_DIR, "results", "offcritical_gauss_bonnet_profile.csv")
    np.savetxt(csv_path, np.column_stack([sigma_mid, rho_sigma]),
               delimiter=',', header='sigma,rho_sigma', comments='')
    log(f"\n  σ-프로파일 CSV 저장: {csv_path}")

    results['B'] = {'fwhm': fwhm, 'sigma_peak': sigma_peak, 'rho_peak': rho_peak,
                    'total_N': total_B/(2*np.pi), 'pass': pass_B}
    log()
    save_all()

    # ═══════════════════════════════════════════════════════════════════════════
    # §C: 디리클레 확장 (χ₅)
    # ═══════════════════════════════════════════════════════════════════════════
    log("═" * 77)
    log("§C: 디리클레 확장 — χ₅=(·/5) off-critical 검증")
    log("═" * 77)
    log()

    # χ₅ 영점 탐색
    log("χ₅ 영점 탐색...")
    chi5_zeros = find_zeros_chi5(t_min=1.0, t_max=50.0, n_scan=5000)

    if len(chi5_zeros) < 5:
        log("  ⚠️ χ₅ 영점 5개 미만 → §C 건너뜀")
        results['C'] = {'pass': False, 'reason': '영점 부족'}
    else:
        # χ₅ t_max: 5번째와 6번째 영점의 중간점
        if len(chi5_zeros) > 5:
            chi5_t_max = (chi5_zeros[4] + chi5_zeros[5]) / 2.0
        else:
            chi5_t_max = chi5_zeros[4] + 1.0
        chi5_t_min = T_MIN_GLOBAL

        log(f"  χ₅ 5영점 범위: t∈[{chi5_t_min},{chi5_t_max:.4f}]")
        log(f"  영점들: {[f'{z:.4f}' for z in chi5_zeros[:6]]}")
        log()

        results_C = []

        # C1: 좁은 strip [0.49, 0.51] — 기대 N=5
        log("--- C1: χ₅ 좁은 strip [0.49, 0.51], 기대=5 ---")
        flux_C1, _, _, _, _ = compute_area_integral(
            log_lambda_chi5_imag, 0.49, 0.51, chi5_t_min, chi5_t_max,
            n_sigma=N_SIGMA_BASE, label="C1"
        )
        N_C1 = flux_C1 / (2 * np.pi)
        err_C1 = abs(N_C1 - 5.0)
        pass_C1 = err_C1 < 1.0
        log(f"  → N = {N_C1:.6f}, 기대=5, |오차|={err_C1:.6f}")
        log(f"  → 판정: {'✅ PASS' if pass_C1 else '❌ FAIL'}")
        results_C.append({'id': 'C1', 'desc': '좁은 strip [0.49,0.51]',
                          'N': N_C1, 'expected': 5, 'err': err_C1, 'pass': pass_C1})
        log()
        save_all()

        # C2: off-critical [0.6, 0.9] — 기대 N=0
        log("--- C2: χ₅ off-critical [0.6, 0.9], 기대=0 ---")
        flux_C2, _, _, _, _ = compute_area_integral(
            log_lambda_chi5_imag, 0.60, 0.90, chi5_t_min, chi5_t_max,
            n_sigma=N_SIGMA_BASE, label="C2"
        )
        N_C2 = flux_C2 / (2 * np.pi)
        pass_C2 = abs(N_C2) < 0.5
        log(f"  → N = {N_C2:.6f}, 기대=0, |N|={abs(N_C2):.6f}")
        log(f"  → 판정: {'✅ PASS' if pass_C2 else '❌ FAIL'}")
        results_C.append({'id': 'C2', 'desc': 'off-critical [0.6,0.9]',
                          'N': N_C2, 'expected': 0, 'err': abs(N_C2), 'pass': pass_C2})

        results['C'] = results_C
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
    log("§A: σ-strip 변조")
    log(f"  {'테스트':<6} {'σ-범위':<18} {'N 계산':>10} {'기대':>6} {'|오차|':>10} {'판정'}")
    log(f"  {'-'*6} {'-'*18} {'-'*10} {'-'*6} {'-'*10} {'-'*6}")
    n_pass_A = 0
    for r in results['A']:
        mark = '✅' if r['pass'] else '❌'
        sig_range = f"[{r['sig_lo']},{r['sig_hi']}]"
        log(f"  {r['id']:<6} {sig_range:<18} {r['N']:10.4f} {r['expected']:6d} {r['err']:10.6f} {mark}")
        if r['pass']:
            n_pass_A += 1
    log(f"  → A 통과: {n_pass_A}/7")
    log()

    log("§B: σ-프로파일")
    rb = results['B']
    log(f"  피크 위치: σ={rb['sigma_peak']:.4f} (기대: 0.5)")
    log(f"  FWHM: {rb['fwhm']:.4f} (기준: < 0.01)")
    log(f"  전체 N: {rb['total_N']:.4f} (기대: ~20)")
    log(f"  → B 판정: {'✅ PASS' if rb['pass'] else '❌ FAIL (격자 해상도 의존)'}")
    log()

    log("§C: 디리클레 확장 (χ₅)")
    n_pass_C = 0
    if isinstance(results.get('C'), list):
        for r in results['C']:
            mark = '✅' if r['pass'] else '❌'
            log(f"  {r['id']} ({r['desc']}): N={r['N']:.4f}, 기대={r['expected']}, {mark}")
            if r['pass']:
                n_pass_C += 1
        log(f"  → C 통과: {n_pass_C}/2")
    else:
        log(f"  → C 건너뜀: {results.get('C', {}).get('reason', '알 수 없음')}")
    log()

    # 전체 통과율
    all_pass = [r['pass'] for r in results['A']]
    all_pass.append(results['B']['pass'])
    if isinstance(results.get('C'), list):
        all_pass += [r['pass'] for r in results['C']]

    n_total = len(all_pass)
    n_total_pass = sum(all_pass)
    log(f"전체 통과: {n_total_pass}/{n_total}")
    log()

    # 물리적 해석
    log("─── 물리적 해석 ───")
    log()
    log("σ-국소화 결과 분석:")
    r_A2 = next(r for r in results['A'] if r['id'] == 'A2')
    r_A3 = next(r for r in results['A'] if r['id'] == 'A3')
    r_A4 = next(r for r in results['A'] if r['id'] == 'A4')
    r_A5 = next(r for r in results['A'] if r['id'] == 'A5')
    r_A6 = next(r for r in results['A'] if r['id'] == 'A6')
    r_A7 = next(r for r in results['A'] if r['id'] == 'A7')

    log(f"  좁은 strip A2 (폭 0.02): N={r_A2['N']:.4f} → {'임계선 집중 확인' if r_A2['pass'] else '이상'}")
    log(f"  극좁은 strip A3 (폭 0.002): N={r_A3['N']:.4f} → {'δ-함수적 집중' if r_A3['pass'] else '격자 해상도 한계'}")
    log(f"  off-critical A4 (σ∈[0.6,0.9]): N={r_A4['N']:.4f} → {'RH 기하학적 확인' if r_A4['pass'] else '반례 가능성!'}")
    log(f"  off-critical A5 (σ∈[0.1,0.4]): N={r_A5['N']:.4f} → {'RH 기하학적 확인' if r_A5['pass'] else '반례 가능성!'}")
    log(f"  반-strip A6 (σ∈[0.3,0.499]): N={r_A6['N']:.4f} → {'임계선 필수 확인' if r_A6['pass'] else '이상'}")
    log(f"  반-strip A7 (σ∈[0.501,0.7]): N={r_A7['N']:.4f} → {'임계선 필수 확인' if r_A7['pass'] else '이상'}")
    log()

    # 최종 판정
    a_critical_pass = all(r['pass'] for r in results['A'])
    a_offcrit_pass = all(r['pass'] for r in results['A'] if r['expected'] == 0)
    a_narrow_pass = all(r['pass'] for r in results['A'] if r['expected'] == 20)

    if a_critical_pass and results['B']['pass'] and n_pass_C >= 2:
        verdict = "★ 판정: 양성 (완전) — σ=1/2 δ-함수적 집중 확립"
    elif a_critical_pass and results['B']['pass']:
        verdict = "★ 판정: 양성 — A+B 완전 통과 (C는 부분)"
    elif a_offcrit_pass and a_narrow_pass:
        verdict = "★ 판정: 조건부 양성 — off-critical=0 + 임계선 포착 확인"
    elif a_offcrit_pass:
        verdict = "★ 판정: 부분 양성 — off-critical 0 확인, 좁은 strip 정밀도 부족"
    else:
        verdict = "★ 판정: 검토 필요 — off-critical에서 비제로 자속 감지"

    log(verdict)
    log()

    elapsed_total = time.time() - t_start
    log(f"총 실행 시간: {elapsed_total:.1f}초 ({elapsed_total/60:.1f}분)")
    log(f"dps={DPS}")

    save_all()
    log()
    log(f"결과 저장: {OUTPUT_TXT}")


if __name__ == "__main__":
    main()
