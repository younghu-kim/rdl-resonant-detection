"""
=============================================================================
[Project RDL] 三重 위상적 일관성 표 (Triple Topological Consistency Table)
=============================================================================
사이클 36 / 수학자 지시 2026-04-14 12:50

목적:
  4본주 프레임워크(곡률, Chern, GB, σ-국소화)의 내적 정합성 최종 검증.
  동일 영역에서 3가지 독립 위상 불변량 계산법이 동일한 정수를 주는지 직접 비교.

  | N  | t2 (approx) | Method A: Direct | Method B: Chern ∮ | Method C: GB ∫∫ |
  |----|--------------|-----------------|--------------------|-----------------|
  | 10 | ~51.4        | 10              | ?                  | ?               |
  | 20 | ~78.2        | 20              | ?                  | ?               |
  | 50 | ~144.6       | 50              | ?                  | ?               |
  |100 | ~237.2       | 100             | ?                  | ?               |

  추가: σ-국소화 고높이 확인 (N=100)
    - 좁은 strip [0.49,0.51] → N=100 기대
    - off-critical [0.6,0.9] → N=0 기대

성공 기준:
  - Direct Count: 정확한 정수 (정의상 exact)
  - Chern: |N_ctr - N| < 0.05 (기존 오프셋 ~0.012)
  - GB: |N_GB - N| < 0.001 (위상적 양자화, 기존 오차 0.000000)
  - σ-국소화 (N=100): |N_narrow - 100| < 0.001, |N_off| < 0.001

출력:
  results/triple_consistency_table.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT_TXT = os.path.join(BASE_DIR, "results", "triple_consistency_table.txt")
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# ─── 설정 ─────────────────────────────────────────────────────────────────────
DPS       = 100
mpmath.mp.dps = DPS

SIGMA_MIN = 0.3
SIGMA_MAX = 0.7
T1        = 13.0   # 첫 영점(14.135) 아래

N_SIGMA   = 64     # σ 방향 격자 점 수 (GB)
N_HORIZ   = 64     # 수평변 분할 (Chern)
N_VERT_PER_UNIT = 3  # 수직변: 단위 t당 분할 수 (Chern)

# ─── 로그 ─────────────────────────────────────────────────────────────────────
log_lines = []
def log(msg=""):
    print(str(msg), flush=True)
    log_lines.append(str(msg))

def save_all():
    with open(OUTPUT_TXT, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))


# ═══════════════════════════════════════════════════════════════════════════════
# §1. 위상 함수 — θ(σ,t) = Im(log ξ(σ+it))
# (holonomy_chern_number.py + gauss_bonnet_area_integral.py 공통 로직)
# ═══════════════════════════════════════════════════════════════════════════════

def log_xi_imag_scalar(s):
    """
    mpmath.mpc s = σ+it 를 받아 Im(log ξ(s)) 반환.
    log ξ(s) = log(s) + log(s-1) - (s/2)log π + loggamma(s/2) + log ζ(s)
    ζ(s) 영점 근처이면 None 반환.
    """
    z_val = mpmath.zeta(s)
    if abs(z_val) < mpmath.mpf(10) ** (-DPS + 15):
        return None
    lxi = (mpmath.log(s)
           + mpmath.log(s - 1)
           - (s / 2) * mpmath.log(mpmath.pi)
           + mpmath.loggamma(s / 2)
           + mpmath.log(z_val))
    return float(lxi.imag)


def log_xi_imag_float(sigma, t):
    """(float sigma, float t) → Im(log ξ) | None. GB 격자용."""
    s = mpmath.mpc(sigma, t)
    return log_xi_imag_scalar(s)


# ═══════════════════════════════════════════════════════════════════════════════
# §2. Method B — Chern 위상 추적 (Rectangle Contour Phase Tracking)
# ═══════════════════════════════════════════════════════════════════════════════

def phase_track_contour(sigma_lo, sigma_hi, t_lo, t_hi, n_horiz, n_vert):
    """
    직사각형 반시계 방향 위상 추적:
      하변(σ 방향 →) + 우변(t 방향 ↑) + 상변(σ 방향 ←) + 좌변(t 방향 ↓)
    Im(log ξ)의 차분을 누적하여 총 위상 변화 / (2π) = N_Chern 추정.

    t₂는 반드시 영점 사이 중간점으로 설정해야 함 (경로가 영점을 지나지 않도록).
    """
    sig_lo = float(sigma_lo)
    sig_hi = float(sigma_hi)
    t0 = float(t_lo)
    t1 = float(t_hi)

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
        # 좌변: t=t1 → t0, sigma=sig_lo (마지막 점 제외)
        for k in range(1, n_vert):
            t = t1 - k * (t1 - t0) / n_vert
            pts.append(mpmath.mpc(sig_lo, t))
        return pts

    path = make_path()
    total_phase = 0.0
    phi_prev = None
    skip_count = 0

    for s in path:
        phi = log_xi_imag_scalar(s)
        if phi is None:
            skip_count += 1
            continue
        if phi_prev is None:
            phi_prev = phi
            continue
        diff = phi - phi_prev
        diff -= 2 * np.pi * round(diff / (2 * np.pi))
        total_phase += diff
        phi_prev = phi

    if skip_count > 0:
        log(f"  WARNING: {skip_count}개 점 건너뜀 (ξ 영점 근처)")

    N_ctr = total_phase / (2 * np.pi)
    return N_ctr


# ═══════════════════════════════════════════════════════════════════════════════
# §3. Method C — Gauss-Bonnet 격자 플라켓 면적분
# ═══════════════════════════════════════════════════════════════════════════════

def wrap_phase(dphi):
    """위상 차분을 [-π, π) 범위로 감김."""
    return dphi - 2.0 * np.pi * np.round(dphi / (2.0 * np.pi))


def compute_plaquette_flux(theta_grid):
    """
    2D 위상 격자에서 각 플라켓(셀)의 자속 계산.
    F[i,j] = wrap(Δσ[i,j]) + wrap(Δt[i+1,j]) - wrap(Δσ[i,j+1]) - wrap(Δt[i,j])
    반환: flux_grid (Nσ-1, Nt-1), total_flux
    """
    d_sigma = wrap_phase(theta_grid[1:, :] - theta_grid[:-1, :])   # (Ns-1, Nt)
    d_t     = wrap_phase(theta_grid[:, 1:] - theta_grid[:, :-1])   # (Ns, Nt-1)
    flux = (d_sigma[:, :-1]
            + d_t[1:, :]
            - d_sigma[:, 1:]
            - d_t[:-1, :])
    total_flux = float(np.sum(flux))
    return flux, total_flux


def compute_area_integral(phase_fn_2d, sigma_min, sigma_max, t_min, t_max,
                          n_sigma=64, n_t=None, label=""):
    """
    2D 영역의 격자 플라켓 곡률 면적분.
    phase_fn_2d: (sigma: float, t: float) → float|None
    반환: N_GB = total_flux / (2π)
    """
    height = t_max - t_min
    if n_t is None:
        n_t = max(64, int(height * 10))

    sigma_arr = np.linspace(sigma_min, sigma_max, n_sigma)
    t_arr     = np.linspace(t_min, t_max, n_t)

    log(f"  격자: Nσ={n_sigma}, Nt={n_t}, 총 {n_sigma*n_t}점")

    theta_grid = np.zeros((n_sigma, n_t))
    skip_count = 0
    t0_c = time.time()

    for i, sig in enumerate(sigma_arr):
        for j, t_val in enumerate(t_arr):
            val = phase_fn_2d(float(sig), float(t_val))
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

        if (i + 1) % 16 == 0:
            elapsed = time.time() - t0_c
            log(f"  ... σ {i+1}/{n_sigma} 완료 ({elapsed:.0f}초 경과)")

    elapsed = time.time() - t0_c
    log(f"  격자 계산 완료: {elapsed:.1f}초, 건너뜀: {skip_count}점")

    flux_grid, total_flux = compute_plaquette_flux(theta_grid)
    N_gb = total_flux / (2 * np.pi)
    log(f"  전체 자속 / (2π) = {N_gb:.6f}")

    return N_gb


# ═══════════════════════════════════════════════════════════════════════════════
# §4. 영점 수집
# ═══════════════════════════════════════════════════════════════════════════════

def collect_zeta_zeros(n_max):
    """mpmath.zetazero로 첫 n_max개 영점 수집."""
    zeros = []
    log(f"  ζ 영점 수집 중 (n=1..{n_max}, dps=50 임시)...")
    orig_dps = mpmath.mp.dps
    mpmath.mp.dps = 50
    for k in range(1, n_max + 1):
        try:
            z = mpmath.zetazero(k)
            zeros.append(float(z.imag))
        except Exception as e:
            print(f"WARNING: zetazero({k}) 실패: {e}", flush=True)
        if k % 25 == 0:
            print(f"  ... {k}/{n_max}개 완료 (t≈{zeros[-1]:.2f})", flush=True)
    mpmath.mp.dps = orig_dps
    zeros = sorted(zeros)
    log(f"  수집 완료: {len(zeros)}개, t범위 [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    return zeros


def midpoint_between_zeros(zeros, n):
    """n번째와 n+1번째 영점 사이 중간점 (n: 0-indexed)."""
    return (zeros[n] + zeros[n + 1]) / 2.0


def count_zeros_below(zeros, t_thresh):
    return sum(1 for z in zeros if z < t_thresh)


# ═══════════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    log("=" * 77)
    log("[Project RDL] 三重 위상적 일관성 표 (Triple Topological Consistency Table)")
    log(f"실행: {now} | dps={DPS}")
    log("=" * 77)
    log()
    log("방법 A: Direct Count (mpmath.zetazero — 정의상 정확)")
    log("방법 B: Chern 윤곽적분 (Im(log ξ) 위상 추적, σ∈[0.3,0.7])")
    log("방법 C: Gauss-Bonnet 면적분 (격자 플라켓 자속, σ∈[0.3,0.7])")
    log("공통: t₁=13.0, t₂=γ_N과 γ_{N+1} 사이 중간점 (경로 영점 회피)")
    log()

    # ─── 영점 수집 (101개 필요: N=100일 때 γ_100, γ_101 모두 필요) ─────────
    log("§0. ζ 영점 수집 (101개)")
    log("-" * 50)
    zeros = collect_zeta_zeros(101)
    save_all()

    # ─── Triple Consistency Table ──────────────────────────────────────────
    log()
    log("§1. Triple Consistency Table")
    log("=" * 77)
    log(f"  {'N':>5}  {'t2':>10}  {'A(Direct)':>12}  {'B(Chern)':>12}  {'C(GB)':>12}  "
        f"{'|B-N|':>8}  {'|C-N|':>8}  판정")
    log("=" * 77)

    N_LIST = [10, 20, 50, 100]
    table_rows = []

    # 성공 기준
    CHERN_TOL = 0.05   # |N_Chern - N| < 0.05
    GB_TOL    = 0.001  # |N_GB - N| < 0.001

    for N_target in N_LIST:
        log()
        log(f"━━━ N = {N_target} ━━━")

        # t2: γ_N과 γ_{N+1} 사이 중간점 (0-indexed: N_target-1 과 N_target)
        t2 = midpoint_between_zeros(zeros, N_target - 1)
        N_direct = count_zeros_below(zeros, t2)  # = N_target (정의상)

        log(f"  t2 = {t2:.4f} (γ_{N_target}={zeros[N_target-1]:.4f} ~ γ_{N_target+1}={zeros[N_target]:.4f})")
        log(f"  Direct Count = {N_direct}")

        # Method B: Chern 위상 추적
        log(f"  [B] Chern 위상 추적 시작...")
        n_vert_b = max(64, int((t2 - T1) * N_VERT_PER_UNIT))
        t0_b = time.time()
        N_chern = phase_track_contour(
            SIGMA_MIN, SIGMA_MAX, T1, t2,
            N_HORIZ, n_vert_b
        )
        elapsed_b = time.time() - t0_b
        log(f"  [B] N_Chern = {N_chern:.6f}  ({elapsed_b:.1f}초)")
        save_all()

        # Method C: Gauss-Bonnet 면적분
        log(f"  [C] GB 면적분 시작...")
        t0_c = time.time()
        N_gb = compute_area_integral(
            log_xi_imag_float,
            SIGMA_MIN, SIGMA_MAX,
            0.5, t2,   # t_min=0.5 (t>1 영역, 첫 영점 훨씬 아래)
            n_sigma=N_SIGMA,
            n_t=None,  # 자동 (높이에 비례)
            label=f"N={N_target}"
        )
        elapsed_c = time.time() - t0_c
        log(f"  [C] N_GB = {N_gb:.6f}  ({elapsed_c:.1f}초)")
        save_all()

        # 판정
        err_b = abs(N_chern - N_direct)
        err_c = abs(N_gb - N_direct)
        ok_b = "✅" if err_b < CHERN_TOL else "❌"
        ok_c = "✅" if err_c < GB_TOL else "❌"
        row_ok = "✅" if (err_b < CHERN_TOL and err_c < GB_TOL) else "❌"

        log(f"  판정: Chern {ok_b} (|오차|={err_b:.6f}), GB {ok_c} (|오차|={err_c:.6f})")
        table_rows.append((N_target, t2, N_direct, N_chern, N_gb, err_b, err_c, ok_b, ok_c, row_ok))

        # 중간 요약 행 출력
        log()
        log(f"  {'N':>5}  {'t2':>10}  {'A':>12}  {'B':>12}  {'C':>12}  "
            f"{'|B-A|':>8}  {'|C-A|':>8}  판정")
        log(f"  {N_target:>5}  {t2:>10.4f}  {N_direct:>12}  {N_chern:>12.6f}  {N_gb:>12.6f}  "
            f"{err_b:>8.6f}  {err_c:>8.6f}  {row_ok}")

    save_all()

    # ─── 최종 표 출력 ─────────────────────────────────────────────────────
    log()
    log("=" * 77)
    log("【 最終 Triple Consistency Table 】")
    log("=" * 77)
    log(f"  {'N':>5}  {'t2':>10}  {'A(Direct)':>10}  {'B(Chern)':>12}  {'C(GB)':>12}  "
        f"{'|B-A|':>10}  {'|C-A|':>10}  {'Chern':>6}  {'GB':>4}  종합")
    log("-" * 100)
    for row in table_rows:
        N, t2, A, B, C, eb, ec, okb, okc, rowok = row
        log(f"  {N:>5}  {t2:>10.4f}  {A:>10}  {B:>12.6f}  {C:>12.6f}  "
            f"{eb:>10.6f}  {ec:>10.6f}  {okb:>6}  {okc:>4}  {rowok}")
    log("-" * 100)

    n_pass_chern = sum(1 for r in table_rows if r[7] == "✅")
    n_pass_gb    = sum(1 for r in table_rows if r[8] == "✅")
    n_pass_total = sum(1 for r in table_rows if r[9] == "✅")
    log(f"  Chern 통과: {n_pass_chern}/4,  GB 통과: {n_pass_gb}/4,  종합: {n_pass_total}/4")

    # ─── σ-국소화 고높이 확인 (N=100) ─────────────────────────────────────
    log()
    log("§2. σ-국소화 고높이 확인 (N=100)")
    log("=" * 60)

    # N=100의 t2
    t2_100 = midpoint_between_zeros(zeros, 99)  # γ_100 ~ γ_101 중간
    log(f"  t2 = {t2_100:.4f} (N=100용)")
    log()

    # §2-A: 좁은 strip [0.49, 0.51] → N=100 기대
    log("  §2-A: 좁은 strip σ∈[0.49, 0.51] → N=100 기대")
    t0_2a = time.time()
    N_gb_narrow = compute_area_integral(
        log_xi_imag_float,
        0.49, 0.51,
        0.5, t2_100,
        n_sigma=64,   # 좁은 범위에서도 64점
        n_t=None,
        label="N100-narrow"
    )
    log(f"  §2-A 결과: N = {N_gb_narrow:.6f}, 기대=100, |오차|={abs(N_gb_narrow-100):.6f}  "
        f"{'✅' if abs(N_gb_narrow-100) < GB_TOL else '❌'}  ({time.time()-t0_2a:.1f}초)")
    save_all()

    # §2-B: off-critical σ∈[0.6, 0.9] → N=0 기대
    log()
    log("  §2-B: off-critical σ∈[0.6, 0.9] → N=0 기대")
    t0_2b = time.time()
    N_gb_off = compute_area_integral(
        log_xi_imag_float,
        0.6, 0.9,
        0.5, t2_100,
        n_sigma=64,
        n_t=None,
        label="N100-offcrit"
    )
    log(f"  §2-B 결과: N = {N_gb_off:.6f}, 기대=0, |오차|={abs(N_gb_off):.6f}  "
        f"{'✅' if abs(N_gb_off) < GB_TOL else '❌'}  ({time.time()-t0_2b:.1f}초)")
    save_all()

    # ─── 최종 판정 ─────────────────────────────────────────────────────────
    log()
    log("=" * 77)
    log("【 최종 판정 】")
    log("=" * 77)

    # 성공 집계
    all_pass_triple = (n_pass_total == 4)
    pass_sigma_narrow = abs(N_gb_narrow - 100) < GB_TOL
    pass_sigma_off    = abs(N_gb_off) < GB_TOL

    # 방법별 정밀도 비교 (기대 오프셋 기록)
    log()
    log("방법별 정밀도 비교:")
    log(f"  Direct Count (A): 항상 정확 (정의상)")
    if table_rows:
        offsets_b = [r[5] for r in table_rows]
        offsets_c = [r[6] for r in table_rows]
        log(f"  Chern (B):        오차 범위 {min(offsets_b):.6f} ~ {max(offsets_b):.6f} (체계적 오프셋 ~0.012)")
        log(f"  GB (C):           오차 범위 {min(offsets_c):.6f} ~ {max(offsets_c):.6f} (위상적 양자화)")
    log()
    log(f"  Chern 체계적 오프셋: GB보다 ~{max(offsets_b)/max(offsets_c+[1e-9]):.0f}배 큰 오차 (수평변 수치 이산화)")
    log(f"  → 두 방법의 정밀도 차이 자체가 위상적 불변량의 본질 반영")
    log()

    log("Triple Consistency Table 결과:")
    log(f"  Chern {n_pass_chern}/4 {'✅' if n_pass_chern==4 else '❌'}")
    log(f"  GB    {n_pass_gb}/4 {'✅' if n_pass_gb==4 else '❌'}")
    log(f"  종합  {n_pass_total}/4 {'✅' if all_pass_triple else '❌'}")
    log()
    log("σ-국소화 (N=100):")
    log(f"  §2-A (좁은 strip): {'✅' if pass_sigma_narrow else '❌'} N={N_gb_narrow:.6f}")
    log(f"  §2-B (off-critical): {'✅' if pass_sigma_off else '❌'} N={N_gb_off:.6f}")
    log()

    if all_pass_triple and pass_sigma_narrow and pass_sigma_off:
        verdict = "★ 완전 양성 (6/6)"
        log(f"최종: {verdict}")
        log()
        log("해석: 세 독립적 위상 불변량 계산법(Direct, Chern, GB)이 N=10,20,50,100에서 모두 일치.")
        log("      σ-국소화도 N=100에서 확인됨 (좁은 strip: 100, off-critical: 0).")
        log("      프레임워크의 내적 정합성 최종 확인.")
    elif n_pass_total >= 3 and (pass_sigma_narrow or pass_sigma_off):
        verdict = f"★ 부분 양성 ({n_pass_total}/4 triple + σ-일부)"
        log(f"최종: {verdict}")
    else:
        verdict = "❌ 불일치 발생 — 재검토 필요"
        log(f"최종: {verdict}")

    total_elapsed = time.time() - t_start
    log()
    log(f"총 실행 시간: {total_elapsed:.1f}초 ({total_elapsed/60:.1f}분)")
    log(f"결과 파일: {OUTPUT_TXT}")
    log("=" * 77)

    save_all()
    print(f"\n[완료] {verdict} — 결과: {OUTPUT_TXT}", flush=True)


if __name__ == "__main__":
    main()
