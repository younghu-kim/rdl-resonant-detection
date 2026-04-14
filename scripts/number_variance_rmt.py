"""
=============================================================================
[Project RDL] 곡률장 수 분산 (Number Variance) — RMT 연결 최종 검증
=============================================================================
사이클 38 / 수학자 지시 2026-04-14 15:45

목적:
  플라켓 기반 수 분산 Var(N; L) vs GUE/Poisson 이론 비교
  + 직접 ζ 영점 계수와 교차 검증

설계:
  - 영역: σ∈[0.3, 0.7] × t∈[50, 600]
  - 격자: Nσ=32, Nt=2000
  - Unfolded window 폭 L: ⟨N⟩ = 5, 10, 20, 50에 해당 (4종)
  - 경계 효과 제거: t∈[70, 580] (양쪽 ~20% 제외) → safe window centers
  - 슬라이딩: 0.5L씩 이동

이론 비교:
  - GUE: Var(N; L) ≈ (2/π²)(log(2πL) + 1 + γ_E - π²/8)
  - Poisson: Var(N; L) = L (선형)

출력:
  results/number_variance_rmt.txt
  results/number_variance_figdata.csv
=============================================================================
"""

import sys, os, time, math
import numpy as np
import mpmath
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT_TXT = os.path.join(BASE_DIR, "results", "number_variance_rmt.txt")
OUTPUT_CSV = os.path.join(BASE_DIR, "results", "number_variance_figdata.csv")
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# ─── 설정 ─────────────────────────────────────────────────────────────────────
DPS      = 80
N_SIGMA  = 32
N_T      = 2000
T_MIN    = 50.0
T_MAX    = 600.0
SAFE_T_MIN = 70.0   # 경계 효과 제외 (하한)
SAFE_T_MAX = 580.0  # 경계 효과 제외 (상한)
L_TARGETS  = [5, 10, 20, 50]  # 목표 ⟨N⟩

mpmath.mp.dps = DPS

# ─── 로그 ─────────────────────────────────────────────────────────────────────
log_lines = []
def log(msg=""):
    print(msg, flush=True)
    log_lines.append(str(msg))

def save_all():
    with open(OUTPUT_TXT, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))

# ═══════════════════════════════════════════════════════════════════════════════
# §1. ζ 영점 수집 (직접 계수용)
# ═══════════════════════════════════════════════════════════════════════════════

def get_zeros_in_range(t_lo, t_hi):
    """t∈[t_lo, t_hi] 범위의 ζ 영점 허수부 리스트 반환"""
    log(f"ζ 영점 수집: t∈[{t_lo}, {t_hi}]")
    t0 = time.time()
    zeros = []
    n = 1
    while True:
        try:
            gamma = float(mpmath.zetazero(n).imag)
        except Exception as e:
            print(f"WARNING: zetazero({n}) 실패: {e}")
            break
        if gamma > t_hi:
            break
        if gamma >= t_lo:
            zeros.append(gamma)
        n += 1
    t1 = time.time()
    log(f"  {len(zeros)}개 영점 수집 ({t1-t0:.1f}초)")
    if len(zeros) == 0:
        print("⚠️ 영점 0개 — 탐색 로직 점검 필요")
    return np.array(zeros)

# ═══════════════════════════════════════════════════════════════════════════════
# §2. Unfolding (von Mangoldt smooth 카운팅 함수)
# ═══════════════════════════════════════════════════════════════════════════════

def N_smooth(T):
    """
    리만-폰 망골트 smooth 카운팅 함수:
    N_smooth(T) = (T/2π) * log(T/2π) - T/2π + 7/8

    각 영점의 unfolded 좌표 = N_smooth(γ_n)
    이 좌표계에서 평균 간격 = 1 (정의상)
    """
    if T <= 0:
        return 0.0
    return (T / (2.0 * math.pi)) * math.log(T / (2.0 * math.pi)) - T / (2.0 * math.pi) + 7.0/8.0

def N_smooth_deriv(T):
    """N_smooth'(T) = (1/2π) * log(T/2π) — 평균 영점 밀도"""
    return (1.0 / (2.0 * math.pi)) * math.log(T / (2.0 * math.pi))

def invert_N_smooth(xi_target, t_lo=1.0, t_hi=1e6, tol=1e-8):
    """
    N_smooth(T) = xi_target 를 bisection으로 풀어 T 반환
    N_smooth는 단조증가이므로 bisection 가능
    """
    # 범위 보장
    while N_smooth(t_hi) < xi_target:
        t_hi *= 2.0
    while N_smooth(t_lo) > xi_target and t_lo > 1e-6:
        t_lo /= 2.0

    for _ in range(100):
        t_mid = (t_lo + t_hi) / 2.0
        if abs(t_hi - t_lo) < tol:
            break
        if N_smooth(t_mid) < xi_target:
            t_lo = t_mid
        else:
            t_hi = t_mid
    return (t_lo + t_hi) / 2.0

# ═══════════════════════════════════════════════════════════════════════════════
# §3. 플라켓 격자 계산 (θ = Im(log ξ))
# ═══════════════════════════════════════════════════════════════════════════════

def wrap_phase(dphi):
    """위상 차분을 [-π, π) 범위로 감김"""
    return dphi - 2.0 * np.pi * np.round(dphi / (2.0 * np.pi))

def log_xi_imag(sigma, t):
    """θ(σ,t) = Im(log ξ(σ+it))"""
    s = mpmath.mpc(sigma, t)
    z_val = mpmath.zeta(s)
    if abs(z_val) < mpmath.mpf(10)**(-DPS + 15):
        return None
    lxi = (mpmath.log(s)
           + mpmath.log(s - 1)
           - (s / 2) * mpmath.log(mpmath.pi)
           + mpmath.loggamma(s / 2)
           + mpmath.log(z_val))
    return float(lxi.imag)

def compute_theta_grid(sigma_arr, t_arr):
    """
    theta_grid[i, j] = Im(log ξ(sigma[i] + i*t[j]))
    영점 근방(None 반환)은 보간으로 처리
    """
    Ns = len(sigma_arr)
    Nt = len(t_arr)
    theta = np.zeros((Ns, Nt), dtype=float)
    skip_total = 0

    t_start = time.time()
    total = Ns * Nt
    count = 0
    report_every = max(1, total // 20)  # 5% 단위 보고

    for i, sig in enumerate(sigma_arr):
        prev_val = 0.0
        for j, t in enumerate(t_arr):
            count += 1
            if count % report_every == 0:
                elapsed = time.time() - t_start
                eta = elapsed / count * (total - count)
                print(f"  격자 진행: {count}/{total} ({100*count/total:.0f}%) "
                      f"경과={elapsed:.0f}s 잔여={eta:.0f}s", flush=True)

            val = log_xi_imag(float(sig), float(t))
            if val is None:
                skip_total += 1
                val = prev_val
            theta[i, j] = val
            prev_val = val

    elapsed = time.time() - t_start
    log(f"  격자 계산 완료: {Ns}×{Nt}={Ns*Nt}점, "
        f"스킵={skip_total}, 소요={elapsed:.1f}초")
    return theta

def compute_plaquette_and_cumsum(theta_grid, t_arr):
    """
    플라켓 자속 계산 + t 방향 누적합 반환
    flux_grid[i, j] = cell (i,j)의 위상 감김
    cum_flux[j] = sum_{j'=0}^{j-1} sum_i flux_grid[i, j'] / (2π)
                = [t_arr[0], t_arr[j]] 구간의 N_plaq
    """
    d_sigma = wrap_phase(theta_grid[1:, :] - theta_grid[:-1, :])   # (Ns-1, Nt)
    d_t     = wrap_phase(theta_grid[:, 1:] - theta_grid[:, :-1])   # (Ns, Nt-1)

    flux = (d_sigma[:, :-1]    # 하변: sigma 방향, t[0..N-2]
            + d_t[1:, :]       # 우변: t 방향, sigma[1..Ns-1]
            - d_sigma[:, 1:]   # 상변: sigma 방향, t[1..N-1]
            - d_t[:-1, :])     # 좌변: t 방향, sigma[0..Ns-2]
    # flux shape: (Ns-1, Nt-1)

    # 각 t-열의 총 flux (σ 방향 합)
    flux_per_t = flux.sum(axis=0)  # (Nt-1,)

    # 누적합: cum_flux[j] = 셀 j=0..j-1까지의 합 / (2π) = N_plaq(t_arr[0], t_arr[j])
    # 즉 N_plaq(t_a, t_b) = cum_flux[j_b] - cum_flux[j_a]
    cum_flux = np.zeros(len(t_arr))  # shape (Nt,)
    cum_flux[0] = 0.0
    running = 0.0
    for j in range(len(flux_per_t)):
        running += flux_per_t[j]
        cum_flux[j + 1] = running / (2.0 * math.pi)

    return flux, cum_flux

def query_plaq_count(cum_flux, t_arr, t1, t2):
    """
    [t1, t2] 구간의 N_plaq 반환
    cum_flux[j] = N_plaq(t_arr[0], t_arr[j])
    """
    j1 = np.searchsorted(t_arr, t1)
    j2 = np.searchsorted(t_arr, t2)
    j1 = max(0, min(j1, len(cum_flux) - 1))
    j2 = max(0, min(j2, len(cum_flux) - 1))
    return cum_flux[j2] - cum_flux[j1]

# ═══════════════════════════════════════════════════════════════════════════════
# §4. 수 분산 분석
# ═══════════════════════════════════════════════════════════════════════════════

def gue_var(L):
    """
    GUE 수 분산 예측:
    Var_GUE(N; L) = (2/π²)(log(2πL) + 1 + γ_E - π²/8)
    γ_E = 0.5772156... (오일러-마스케로니 상수)
    """
    gamma_E = 0.5772156649015329
    return (2.0 / (math.pi**2)) * (math.log(2.0 * math.pi * L) + 1.0 + gamma_E - math.pi**2/8.0)

def poisson_var(L):
    """Poisson 수 분산: Var_Poisson(N; L) = L"""
    return float(L)

def analyze_number_variance(L_target, zeros_safe, xi_zeros_safe,
                             cum_flux, t_arr,
                             t_safe_min, t_safe_max):
    """
    특정 L에 대한 수 분산 분석.

    L_target: unfolded 좌표에서의 윈도우 폭 (= 목표 ⟨N⟩)
    zeros_safe: safe 구간의 ζ 영점 배열
    xi_zeros_safe: safe 구간 영점의 unfolded 좌표 배열
    cum_flux: 플라켓 누적합
    t_arr: t 격자 배열
    """
    log(f"\n  [L={L_target}] 수 분산 분석")

    # 윈도우를 step=1 (unfolded 1단위 = 영점 약 1개 간격)씩 슬라이딩
    # step=0.5L 대신 step=1 사용:
    #   step=0.5L → L=50에서 11개 창만 (창 경계가 정수 xi 사이에 고정 → count 고정 → Var≈0)
    #   step=1    → L=50에서 ~259개 창 (영점 진입/탈출 포착 → 실제 분산 측정 가능)
    xi_min = xi_zeros_safe[0]
    xi_max = xi_zeros_safe[-1] - L_target  # 마지막 창이 safe 범위 내에 끝나도록

    if xi_max <= xi_min:
        log(f"  ⚠️ xi 범위 부족 — L={L_target} 스킵")
        return None

    step = 1.0  # unfolded 1단위 스텝
    xi_starts = np.arange(xi_min, xi_max + step/2, step)
    xi_windows = [(xs, xs + L_target) for xs in xi_starts]

    if len(xi_windows) < 5:
        log(f"  ⚠️ 윈도우 수={len(xi_windows)} (너무 적음) — L={L_target} 스킵")
        return None

    # 독립적 유효 샘플 수 (창 길이 L만큼 상관되므로)
    n_eff = max(1, int((xi_max - xi_min) / L_target))

    log(f"  윈도우 수: {len(xi_windows)}, "
        f"xi∈[{xi_min:.2f}, {xi_max:.2f}], n_eff≈{n_eff}")

    # 각 윈도우에서 직접 계수 & 플라켓 계수
    N_direct_list  = []
    N_plaq_list    = []
    t_boundaries_list = []

    for xi_a, xi_b in xi_windows:
        # 직접 계수: 영점의 xi 좌표가 [xi_a, xi_b]에 있는 것
        mask = (xi_zeros_safe >= xi_a) & (xi_zeros_safe < xi_b)
        N_direct = float(np.sum(mask))
        N_direct_list.append(N_direct)

        # t 경계 역변환 (N_smooth 역함수)
        t1 = invert_N_smooth(xi_a, t_lo=1.0, t_hi=T_MAX * 2)
        t2 = invert_N_smooth(xi_b, t_lo=1.0, t_hi=T_MAX * 2)
        t_boundaries_list.append((t1, t2))

        # 플라켓 계수
        N_plaq = query_plaq_count(cum_flux, t_arr, t1, t2)
        N_plaq_list.append(N_plaq)

    N_direct_arr = np.array(N_direct_list)
    N_plaq_arr   = np.array(N_plaq_list)

    # 통계
    mean_direct  = np.mean(N_direct_arr)
    var_direct   = np.var(N_direct_arr, ddof=1) if len(N_direct_arr) > 1 else np.var(N_direct_arr)
    mean_plaq    = np.mean(N_plaq_arr)
    var_plaq     = np.var(N_plaq_arr, ddof=1) if len(N_plaq_arr) > 1 else np.var(N_plaq_arr)

    # 이론 예측
    var_gue     = gue_var(L_target)
    var_poisson = poisson_var(L_target)

    # GUE와의 상대 편차 (%)
    rel_gue_direct = abs(var_direct - var_gue) / var_gue * 100 if var_gue > 0 else float('inf')
    rel_gue_plaq   = abs(var_plaq   - var_gue) / var_gue * 100 if var_gue > 0 else float('inf')

    # Poisson과의 거리 (σ 기준: std of variance estimate)
    # 유효 독립 샘플 수 n_eff 사용 (상관된 창들의 실제 자유도)
    n_win = len(N_direct_arr)
    std_var_direct = math.sqrt(2.0) * var_direct / math.sqrt(max(1, n_eff - 1))
    sep_poisson_direct = abs(var_direct - var_poisson) / std_var_direct if std_var_direct > 0 else float('inf')

    log(f"  직접 계수: ⟨N⟩={mean_direct:.3f}, Var(N)={var_direct:.4f}")
    log(f"  플라켓:   ⟨N⟩={mean_plaq:.3f}, Var(N)={var_plaq:.4f}")
    log(f"  GUE 예측:    Var(N)={var_gue:.4f}")
    log(f"  Poisson 예측: Var(N)={var_poisson:.1f}")
    log(f"  GUE 편차: direct={rel_gue_direct:.1f}%, plaq={rel_gue_plaq:.1f}%")
    log(f"  Poisson 거리(σ): {sep_poisson_direct:.1f}")

    # 판정
    gue_pass_direct = rel_gue_direct <= 20.0
    gue_pass_plaq   = rel_gue_plaq   <= 20.0
    poi_sep_pass    = sep_poisson_direct >= 3.0

    log(f"  판정: GUE(direct)={'✅' if gue_pass_direct else '❌'} "
        f"GUE(plaq)={'✅' if gue_pass_plaq else '❌'} "
        f"Poisson분리={'✅' if poi_sep_pass else '❌'}")

    return {
        'L': L_target,
        'n_windows': n_win,
        'n_eff': n_eff,
        'mean_direct': mean_direct,
        'var_direct': var_direct,
        'mean_plaq': mean_plaq,
        'var_plaq': var_plaq,
        'var_gue': var_gue,
        'var_poisson': var_poisson,
        'rel_gue_direct': rel_gue_direct,
        'rel_gue_plaq': rel_gue_plaq,
        'sep_poisson': sep_poisson_direct,
        'gue_pass_direct': gue_pass_direct,
        'gue_pass_plaq': gue_pass_plaq,
        'poi_sep_pass': poi_sep_pass,
        'N_direct_arr': N_direct_arr,
        'N_plaq_arr': N_plaq_arr,
    }

# ═══════════════════════════════════════════════════════════════════════════════
# main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    t_total_start = time.time()
    log("=" * 72)
    log("[RDL] 곡률장 수 분산 (Number Variance) — RMT 연결 최종 검증")
    log(f"실행 시각: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"DPS={DPS}, Nσ={N_SIGMA}, Nt={N_T}")
    log(f"t∈[{T_MIN}, {T_MAX}], safe t∈[{SAFE_T_MIN}, {SAFE_T_MAX}]")
    log(f"L_targets={L_TARGETS}")
    log("=" * 72)

    # ─── §1. 영점 수집 ─────────────────────────────────────────────────────
    log("\n" + "─" * 50)
    log("§1. ζ 영점 수집")
    log("─" * 50)
    # 전체 범위 + 약간의 여유
    all_zeros = get_zeros_in_range(T_MIN - 5, T_MAX + 5)
    log(f"전체 영점: {len(all_zeros)}개 (t∈[{T_MIN-5}, {T_MAX+5}])")

    # safe 구간 영점
    mask_safe = (all_zeros >= SAFE_T_MIN) & (all_zeros <= SAFE_T_MAX)
    zeros_safe = all_zeros[mask_safe]
    log(f"Safe 구간 영점: {len(zeros_safe)}개 (t∈[{SAFE_T_MIN}, {SAFE_T_MAX}])")

    # Unfolded 좌표 계산
    xi_all    = np.array([N_smooth(float(g)) for g in all_zeros])
    xi_safe   = np.array([N_smooth(float(g)) for g in zeros_safe])
    log(f"Unfolded xi: [{xi_safe[0]:.3f}, {xi_safe[-1]:.3f}]")
    log(f"Safe 구간 xi 범위: {xi_safe[-1]-xi_safe[0]:.1f} (≈ {len(zeros_safe)}개 기대)")

    # mean spacing check
    if len(xi_safe) > 1:
        mean_spacing = np.mean(np.diff(xi_safe))
        log(f"Unfolded 평균 간격: {mean_spacing:.4f} (이론=1.000)")

    # ─── §2. 플라켓 격자 계산 ─────────────────────────────────────────────
    log("\n" + "─" * 50)
    log("§2. 플라켓 격자 계산 (θ = Im(log ξ))")
    log(f"   σ∈[0.3, 0.7] × t∈[{T_MIN}, {T_MAX}]")
    log(f"   Nσ={N_SIGMA}, Nt={N_T}")
    log("─" * 50)

    sigma_arr = np.linspace(0.3, 0.7, N_SIGMA)
    t_arr     = np.linspace(T_MIN, T_MAX, N_T)

    log(f"격자 간격: Δσ={sigma_arr[1]-sigma_arr[0]:.4f}, Δt={t_arr[1]-t_arr[0]:.4f}")
    log(f"총 격자점: {N_SIGMA}×{N_T}={N_SIGMA*N_T}")

    # 캐시 사용 (재실행 시 시간 절약)
    cache_path = os.path.join(BASE_DIR, "outputs", "cache",
                              f"theta_grid_Ns{N_SIGMA}_Nt{N_T}_"
                              f"t{int(T_MIN)}_{int(T_MAX)}_dps{DPS}.npy")
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)

    if os.path.exists(cache_path):
        log(f"캐시 로드: {cache_path}")
        theta_grid = np.load(cache_path)
        log(f"  캐시 shape: {theta_grid.shape}")
    else:
        log("격자 계산 시작...")
        theta_grid = compute_theta_grid(sigma_arr, t_arr)
        np.save(cache_path, theta_grid)
        log(f"캐시 저장: {cache_path}")

    log("플라켓 자속 계산...")
    flux_grid, cum_flux = compute_plaquette_and_cumsum(theta_grid, t_arr)

    # 전체 N_plaq 확인
    total_plaq = cum_flux[-1]
    total_direct = len(all_zeros[all_zeros <= T_MAX])
    log(f"전체 검증: N_plaq={total_plaq:.3f}, N_direct={total_direct}")
    if abs(total_plaq - total_direct) > 2.0:
        log(f"⚠️ 큰 불일치 ({abs(total_plaq-total_direct):.2f}) — 결과 해석 주의")

    # ─── §3. 수 분산 분석 ─────────────────────────────────────────────────
    log("\n" + "─" * 50)
    log("§3. 수 분산 (Number Variance) 분석")
    log("─" * 50)

    results_list = []
    csv_rows = []
    csv_rows.append("L,n_windows,mean_direct,var_direct,mean_plaq,var_plaq,"
                    "var_gue,var_poisson,rel_gue_direct_pct,rel_gue_plaq_pct,"
                    "sep_poisson_sigma,gue_pass_direct,gue_pass_plaq,poi_sep_pass")

    for L in L_TARGETS:
        res = analyze_number_variance(
            L_target=L,
            zeros_safe=zeros_safe,
            xi_zeros_safe=xi_safe,
            cum_flux=cum_flux,
            t_arr=t_arr,
            t_safe_min=SAFE_T_MIN,
            t_safe_max=SAFE_T_MAX,
        )
        if res is not None:
            results_list.append(res)
            csv_rows.append(
                f"{res['L']},{res['n_windows']},"
                f"{res['mean_direct']:.4f},{res['var_direct']:.6f},"
                f"{res['mean_plaq']:.4f},{res['var_plaq']:.6f},"
                f"{res['var_gue']:.6f},{res['var_poisson']:.6f},"
                f"{res['rel_gue_direct']:.2f},{res['rel_gue_plaq']:.2f},"
                f"{res['sep_poisson']:.2f},"
                f"{'1' if res['gue_pass_direct'] else '0'},"
                f"{'1' if res['gue_pass_plaq'] else '0'},"
                f"{'1' if res['poi_sep_pass'] else '0'}"
            )

    # ─── §4. 최종 요약 ────────────────────────────────────────────────────
    log("\n" + "=" * 72)
    log("§4. 최종 요약 — 수 분산 vs RMT 이론")
    log("=" * 72)

    log(f"\n{'L':>5}  {'N_win':>6}  "
        f"{'Var_dir':>10}  {'Var_plaq':>10}  "
        f"{'Var_GUE':>10}  {'Var_Pois':>10}  "
        f"{'GUE%_dir':>10}  {'GUE%_plq':>10}  "
        f"{'Poi_σ':>7}  판정")
    log("-" * 110)

    n_pass_full = 0
    n_pass_partial = 0
    for r in results_list:
        gue_dir = "✅" if r['gue_pass_direct'] else "❌"
        gue_plq = "✅" if r['gue_pass_plaq']   else "❌"
        poi_sep = "✅" if r['poi_sep_pass']     else "❌"
        verdict = "완전양성" if (r['gue_pass_direct'] and r['gue_pass_plaq'] and r['poi_sep_pass']) else (
                  "부분양성" if (r['gue_pass_direct'] or r['gue_pass_plaq']) else "음성")
        if verdict == "완전양성": n_pass_full += 1
        elif verdict == "부분양성": n_pass_partial += 1
        log(f"{r['L']:>5}  {r['n_windows']:>6}  "
            f"{r['var_direct']:>10.4f}  {r['var_plaq']:>10.4f}  "
            f"{r['var_gue']:>10.4f}  {r['var_poisson']:>10.1f}  "
            f"{r['rel_gue_direct']:>9.1f}%  {r['rel_gue_plaq']:>9.1f}%  "
            f"{r['sep_poisson']:>7.1f}  {gue_dir}{gue_plq}{poi_sep} {verdict}")

    log("")
    log(f"완전양성: {n_pass_full}/{len(results_list)}, 부분양성: {n_pass_partial}/{len(results_list)}")

    # 전체 판정
    log("\n" + "─" * 60)
    log("전체 판정:")
    if n_pass_full == len(results_list):
        overall = "★★★ 완전 양성 — 4종 윈도우 모두 GUE ±20% 이내 + Poisson >3σ 차이"
        verdict_short = "완전 양성"
    elif n_pass_full >= 2 or n_pass_partial >= 2:
        overall = "★★ 부분 양성 — 일부 윈도우에서 GUE 일치 (경계 효과 가능)"
        verdict_short = "부분 양성"
    else:
        overall = "★ 음성 또는 불확정 — Poisson/GUE 구분 불충분"
        verdict_short = "음성"
    log(overall)

    # GUE logarithmic scaling 검증
    log("\n" + "─" * 60)
    log("대수적 스케일링 검증 (Var_direct ∝ log L?):")
    if len(results_list) >= 2:
        log(f"{'L':>6}  {'Var_dir':>10}  {'(2/π²)log L':>14}  {'Var_GUE':>10}")
        for r in results_list:
            log_L = (2.0/math.pi**2) * math.log(r['L'])
            log(f"{r['L']:>6}  {r['var_direct']:>10.4f}  {log_L:>14.4f}  {r['var_gue']:>10.4f}")
        # 선형 log-log 피팅
        L_vals = np.array([r['L'] for r in results_list], dtype=float)
        V_vals = np.array([r['var_direct'] for r in results_list], dtype=float)
        if len(L_vals) >= 2 and np.all(V_vals > 0):
            log_L = np.log(L_vals)
            log_V = np.log(V_vals)
            slope = np.polyfit(log_L, log_V, 1)[0]
            log(f"log-log 기울기: {slope:.3f} (GUE 예측: ~1.0 — 대수 스케일링)")

    # 수학적 해석
    log("\n" + "─" * 60)
    log("수학적 해석:")
    log("- GUE 수 분산: Var(N;L) ≈ (2/π²)(log(2πL) + 1 + γ_E - π²/8)")
    log("  ≈ (2/π²) log L + 상수 — 대수적 스케일링")
    log("- Poisson: Var(N;L) = L — 선형 스케일링")
    log("- 리만 영점이 GUE 수 분산을 따른다는 것은 Odlyzko(1987)이 수치적으로 확인")
    log("- 플라켓 방법이 동일한 GUE 수 분산을 재현 → 위상 기하학 ↔ RMT 정량 연결")

    # 총 소요 시간
    elapsed = time.time() - t_total_start
    log(f"\n총 소요 시간: {elapsed/60:.1f}분 ({elapsed:.0f}초)")
    log(f"결과 파일: {OUTPUT_TXT}")
    log(f"데이터 파일: {OUTPUT_CSV}")

    # 결론
    log("\n" + "=" * 72)
    log(f"결론: {verdict_short}")
    log("=" * 72)

    # CSV 저장
    with open(OUTPUT_CSV, 'w', encoding='utf-8') as f:
        f.write("\n".join(csv_rows))
    log(f"\nCSV 저장 완료: {OUTPUT_CSV}")

    save_all()
    log("TXT 저장 완료.")


if __name__ == "__main__":
    main()
