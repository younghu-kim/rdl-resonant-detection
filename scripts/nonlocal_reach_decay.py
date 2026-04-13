"""
=============================================================================
[Project RDL] 비국소 도달 범위 측정 — κ_mid ↔ gap_{n+k} 상관 감쇠
=============================================================================
사이클 30 / 수학자 지시 2026-04-14 06:08

목적:
  결과 #18에서 ρ(κ_mid, gap_next) = -0.654 (k=1) 확인.
  k=0..10 (전방) 및 k=1..10 (후방)으로 확장하여 비국소 도달 범위 측정.
  특성 길이(correlation length) ξ 정의 시도.
  ζ + χ₃(+ χ₄, χ₅) 보편성 비교.

메서드:
  1. 리만 ζ: 235개 영점 (t∈[0,450]), 해석적 log(ξ)' 공식 (v2b)
  2. 디리클레 χ₃/χ₄/χ₅: find_zeros_L 탐색 (t∈[1,200])
     해석적 Λ'/Λ 공식 (v2)
  3. 전방 k=0..10: Spearman ρ(κ_n, gap_{n+k})
  4. 후방 k=1..10: Spearman ρ(κ_n, gap_{n-k})
  5. 편상관: ρ(κ, gap_{n+k} | gap_n) — gap_n 기여 제거
  6. 감쇠 피팅: |ρ(k)| ~ A·exp(-k/ξ) vs A·k^{-α}
  7. 전방/후방 비교 (방향성 검증)

출력:
  results/nonlocal_reach_decay.txt
  results/nonlocal_reach_decay_figdata.csv
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats
from scipy.optimize import curve_fit
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT_TXT = os.path.join(BASE_DIR, "results", "nonlocal_reach_decay.txt")
OUTPUT_CSV = os.path.join(BASE_DIR, "results", "nonlocal_reach_decay_figdata.csv")
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# ─── 설정 ────────────────────────────────────────────────────────────────────
DPS_ZETA   = 100   # ζ κ 계산
DPS_DIR_K  = 100   # 디리클레 κ 계산
DPS_DIR_S  = 50    # 디리클레 영점 탐색
H_DIFF     = mpmath.mpf(10) ** (-50)   # dps=100 → 유효 50자리
CAP_VALUE  = 1e18
K_MAX      = 10    # 전방/후방 최대 k

# ─── 로그 ─────────────────────────────────────────────────────────────────────
log_lines = []
def log(msg=""):
    print(msg, flush=True)
    log_lines.append(msg)

def save_all():
    with open(OUTPUT_TXT, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))


# ═══════════════════════════════════════════════════════════════════════════════
# 리만 ζ — 해석적 log(ξ)' 공식 (v2b 재사용)
# ═══════════════════════════════════════════════════════════════════════════════

def connection_xi_analytic(s):
    """
    L(s) = ξ'(s)/ξ(s) = 1/s + 1/(s-1) - log(π)/2 + (1/2)ψ(s/2) + ζ'(s)/ζ(s)
    ξ 자체 계산 불필요 → underflow 완전 회피
    """
    zeta_val = mpmath.zeta(s)
    if abs(zeta_val) < mpmath.mpf(10) ** (-mpmath.mp.dps + 15):
        raise ValueError(f"ζ(s) near zero: |ζ|={float(abs(zeta_val)):.2e}")
    zeta_deriv = mpmath.diff(mpmath.zeta, s)
    term_zeta  = zeta_deriv / zeta_val
    term_s     = mpmath.mpf(1) / s
    term_sm1   = mpmath.mpf(1) / (s - 1)
    term_log   = -mpmath.log(mpmath.pi) / 2
    term_psi   = mpmath.digamma(s / 2) / 2
    return term_s + term_sm1 + term_log + term_psi + term_zeta


def kappa_zeta(t_float):
    """중간점 t에서 κ_ζ = |ξ'/ξ|². None=cap/실패"""
    try:
        mpmath.mp.dps = DPS_ZETA
        s = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_float)))
        L = connection_xi_analytic(s)
        k = float(abs(L) ** 2)
        if not (k == k) or k > CAP_VALUE:
            return None
        return k
    except Exception as e:
        print(f"  WARNING ζ κ: t={t_float:.4f} → {e}", flush=True)
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# 디리클레 L-함수 — 해석적 Λ'/Λ 공식 (v2 재사용)
# ═══════════════════════════════════════════════════════════════════════════════

def connection_Lambda_analytic(s, char_info):
    """
    Λ'/Λ(s,χ) = (1/2)log(q/π) + (1/2)ψ((s+a)/2) + L'(s,χ)/L(s,χ)
    Λ 자체 계산 불필요 → underflow 완전 회피
    """
    mpmath.mp.dps = DPS_DIR_K
    chi = char_info['chi']
    q   = mpmath.mpf(char_info['q'])
    a   = mpmath.mpf(char_info['a'])
    L_c = mpmath.dirichlet(s, chi)
    if abs(L_c) < mpmath.mpf(10) ** (-DPS_DIR_K + 15):
        raise ValueError(f"L(s,χ) near zero")
    L_plus  = mpmath.dirichlet(s + H_DIFF, chi)
    L_minus = mpmath.dirichlet(s - H_DIFF, chi)
    log_deriv_L = (L_plus - L_minus) / (2 * H_DIFF * L_c)
    term_log = mpmath.log(q / mpmath.pi) / 2
    term_psi = mpmath.digamma((s + a) / 2) / 2
    return term_log + term_psi + log_deriv_L


def kappa_dirichlet(t_float, char_info):
    """중간점 t에서 κ_Λ = |Λ'/Λ|². None=cap/실패"""
    try:
        s = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_float)))
        conn = connection_Lambda_analytic(s, char_info)
        k = float(abs(conn) ** 2)
        if not (k == k) or k > CAP_VALUE:
            return None
        return k
    except Exception as e:
        print(f"  WARNING Λ κ: t={t_float:.4f} → {e}", flush=True)
        return None


def completed_L_func(s, char_info):
    """Λ(s,χ) = (q/π)^{s/2} Γ((s+a)/2) L(s,χ) — 영점 탐색 전용"""
    q   = mpmath.mpf(char_info['q'])
    a   = mpmath.mpf(char_info['a'])
    chi = char_info['chi']
    L   = mpmath.dirichlet(s, chi)
    G   = mpmath.gamma((s + a) / 2)
    pre = mpmath.power(q / mpmath.pi, s / 2)
    return pre * G * L


def find_zeros_L(char_info, t_min=1.0, t_max=200.0, n_scan=3000):
    """
    완비 Λ(s,χ)의 Re 부호 변화로 임계선 영점 탐색 (v2와 동일 로직).
    적응형 tol로 underflow 대응.
    """
    ts = np.linspace(t_min, t_max, n_scan)
    zeros = []
    fail_count = 0
    mpmath.mp.dps = DPS_DIR_S

    prev_re = None
    prev_t  = None

    for t in ts:
        s = mpmath.mpc('0.5', str(t))
        try:
            val     = completed_L_func(s, char_info)
            curr_re = float(mpmath.re(val))
        except Exception:
            prev_re = None; prev_t = t
            continue

        if prev_re is not None and prev_re * curr_re < 0:
            try:
                def f_re(t_var, ci=char_info):
                    sv = mpmath.mpc('0.5', str(float(t_var)))
                    return mpmath.re(completed_L_func(sv, ci))
                func_scale = max(abs(prev_re), abs(curr_re))
                local_tol  = max(func_scale * mpmath.mpf('1e-10'),
                                 mpmath.mpf(10) ** (-DPS_DIR_S + 5))
                t_zero     = mpmath.findroot(f_re, (prev_t, t), tol=local_tol)
                t_zero_f   = float(t_zero)
                if not zeros or abs(t_zero_f - zeros[-1]) > 0.05:
                    zeros.append(t_zero_f)
            except Exception as e:
                fail_count += 1

        prev_re = curr_re
        prev_t  = t

    return sorted(zeros), fail_count


# ═══════════════════════════════════════════════════════════════════════════════
# 핵심 분석 함수
# ═══════════════════════════════════════════════════════════════════════════════

def partial_spearman(x, y, z):
    """
    편 Spearman ρ(x, y | z) — z 통제 후 x와 y의 상관.
    rank 기반: r_xy - r_xz * r_yz / sqrt((1-r_xz²)(1-r_yz²))
    """
    rx = stats.rankdata(x).astype(float)
    ry = stats.rankdata(y).astype(float)
    rz = stats.rankdata(z).astype(float)
    r_xy, _ = stats.pearsonr(rx, ry)
    r_xz, _ = stats.pearsonr(rx, rz)
    r_yz, _ = stats.pearsonr(ry, rz)
    denom = np.sqrt(max(1e-12, (1 - r_xz**2) * (1 - r_yz**2)))
    return (r_xy - r_xz * r_yz) / denom


def compute_reach_analysis(zeros_list, kappa_func, label, char_info=None):
    """
    영점 배열 + κ 함수로 k=0..K_MAX 전방/후방 상관 감쇠 계산.

    Returns:
        dict {
          'label': str,
          'N': int,
          'n_valid': int,
          'cap_pct': float,
          'forward': [(k, rho, p, partial_rho, n_used), ...],  k=0..K_MAX
          'backward': [(k, rho, p, partial_rho, n_used), ...], k=1..K_MAX
        }
    """
    zeros = np.array(sorted(zeros_list))
    N = len(zeros)
    pairs = N - 1

    log(f"\n{'─'*70}")
    log(f"  L-함수: {label}")
    log(f"  영점 수: {N}  (쌍 수: {pairs})")

    if N < 30:
        log(f"  ⚠️ 영점 {N}개 < 30 — 분석 불가")
        return None

    # ── 중간점 κ 계산 ─────────────────────────────────────────────────────────
    log(f"  κ 계산 중 (dps={DPS_ZETA if char_info is None else DPS_DIR_K}) ...")
    kappa_list = []
    cap_count  = 0
    for i in range(pairs):
        m = (zeros[i] + zeros[i+1]) / 2.0
        if char_info is None:
            k_val = kappa_func(m)
        else:
            k_val = kappa_func(m, char_info)

        if k_val is None:
            cap_count += 1
            kappa_list.append(None)
        else:
            kappa_list.append(k_val)

        if (i + 1) % 50 == 0:
            kstr = f"{k_val:.4f}" if k_val is not None else "cap"
            log(f"    [{i+1}/{pairs}] m={m:.2f}, κ={kstr}")

    n_valid = sum(1 for k in kappa_list if k is not None)
    cap_pct = cap_count / pairs * 100
    log(f"  유효 κ: {n_valid}/{pairs}  cap {cap_pct:.1f}%")

    if cap_pct > 10:
        log(f"  ⚠️ cap {cap_pct:.1f}% > 10%")
    elif cap_pct > 5:
        log(f"  ⚠️ cap {cap_pct:.1f}% > 5%")
    else:
        log(f"  ✅ cap {cap_pct:.1f}% < 5%")

    if n_valid < 30:
        log(f"  ⚠️ 유효 데이터 {n_valid} < 30 — 분석 불가")
        return None

    # gaps (i → gap_i = zeros[i+1] - zeros[i], 길이 = pairs)
    gaps = zeros[1:] - zeros[:-1]   # shape (pairs,)

    # ── 전방 상관 k=0..K_MAX ─────────────────────────────────────────────────
    forward_results = []
    log(f"\n  [전방 상관] k=0..{K_MAX}")
    log(f"  {'k':>3}  {'ρ':>8}  {'p':>10}  {'편상관':>8}  {'n':>5}")
    log(f"  {'─'*3}  {'─'*8}  {'─'*10}  {'─'*8}  {'─'*5}")

    for k in range(K_MAX + 1):
        # gap_{n+k} = zeros[n+k+1] - zeros[n+k]
        # index n+k+1 ≤ pairs-1 → n ≤ pairs-k-2
        # n runs 0..pairs-1, but also κ_n must be valid
        valid_indices = []
        for i in range(pairs - k - 1 + 1):   # i ≤ pairs - k - 1
            if i + k + 1 <= pairs - 1 and kappa_list[i] is not None:
                valid_indices.append(i)

        # gap_{n+k}: for i in valid_indices → gaps[i+k]
        if len(valid_indices) < 20:
            log(f"  {k:>3}  {'n/a':>8}  {'n/a':>10}  {'n/a':>8}  {len(valid_indices):>5}  (데이터 부족)")
            forward_results.append((k, np.nan, np.nan, np.nan, len(valid_indices)))
            continue

        kappa_k = np.array([kappa_list[i] for i in valid_indices])
        gap_k   = np.array([gaps[i + k]   for i in valid_indices])
        gap_0   = np.array([gaps[i]        for i in valid_indices])

        rho, p   = stats.spearmanr(kappa_k, gap_k)
        if k == 0:
            partial = np.nan   # k=0에서 gap_{n+0} = gap_n → 편상관 trivial
        else:
            partial = partial_spearman(kappa_k, gap_k, gap_0)

        n_used = len(valid_indices)
        log(f"  {k:>3}  {rho:>8.4f}  {p:>10.3e}  {partial:>8.4f}  {n_used:>5}")
        forward_results.append((k, rho, p, partial, n_used))

    # ── 후방 상관 k=1..K_MAX ─────────────────────────────────────────────────
    backward_results = []
    log(f"\n  [후방 상관] k=1..{K_MAX}")
    log(f"  {'k':>3}  {'ρ':>8}  {'p':>10}  {'편상관':>8}  {'n':>5}")
    log(f"  {'─'*3}  {'─'*8}  {'─'*10}  {'─'*8}  {'─'*5}")

    for k in range(1, K_MAX + 1):
        # gap_{n-k} = zeros[n-k+1] - zeros[n-k]
        # n-k ≥ 0 → n ≥ k
        valid_indices = []
        for i in range(k, pairs):
            if kappa_list[i] is not None:
                valid_indices.append(i)

        if len(valid_indices) < 20:
            log(f"  {k:>3}  {'n/a':>8}  {'n/a':>10}  {'n/a':>8}  {len(valid_indices):>5}  (데이터 부족)")
            backward_results.append((k, np.nan, np.nan, np.nan, len(valid_indices)))
            continue

        kappa_k  = np.array([kappa_list[i] for i in valid_indices])
        gap_bk   = np.array([gaps[i - k]   for i in valid_indices])
        gap_0    = np.array([gaps[i]        for i in valid_indices])

        rho, p   = stats.spearmanr(kappa_k, gap_bk)
        partial  = partial_spearman(kappa_k, gap_bk, gap_0)

        n_used = len(valid_indices)
        log(f"  {k:>3}  {rho:>8.4f}  {p:>10.3e}  {partial:>8.4f}  {n_used:>5}")
        backward_results.append((k, rho, p, partial, n_used))

    return {
        'label': label,
        'N': N,
        'n_valid': n_valid,
        'cap_pct': cap_pct,
        'forward': forward_results,
        'backward': backward_results,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 감쇠 피팅
# ═══════════════════════════════════════════════════════════════════════════════

def fit_decay(k_vals, rho_vals, label):
    """
    |ρ(k)| vs k 감쇠 피팅: 지수 vs 멱법칙.
    k=0 포함 (k=0의 |ρ|는 nn 상쇄 확인값, 기준점으로 사용).
    유효값(nan 아님)만 사용.
    """
    k_arr = np.array(k_vals, dtype=float)
    r_arr = np.abs(np.array(rho_vals, dtype=float))

    # k=0 제외하고 k=1..K_MAX만 피팅 (k=0은 nan에 해당하는 경우 있음)
    # 유효 데이터 필터
    mask = (~np.isnan(r_arr)) & (k_arr >= 1)
    k_fit = k_arr[mask]
    r_fit = r_arr[mask]

    if len(k_fit) < 4:
        log(f"\n  [{label}] 피팅 데이터 부족 ({len(k_fit)}개 < 4)")
        return None

    results_fit = {}

    # 지수 피팅: A * exp(-k/xi)
    try:
        def exp_decay(k, A, xi):
            return A * np.exp(-k / xi)
        p0_exp = [r_fit[0], 2.0]
        popt_e, _ = curve_fit(exp_decay, k_fit, r_fit, p0=p0_exp,
                               bounds=([0, 0.1], [2.0, 100.0]), maxfev=5000)
        r_pred_e = exp_decay(k_fit, *popt_e)
        ss_res = np.sum((r_fit - r_pred_e)**2)
        ss_tot = np.sum((r_fit - r_fit.mean())**2)
        r2_e = 1 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
        results_fit['exp'] = {'A': popt_e[0], 'xi': popt_e[1], 'R2': r2_e}
        log(f"\n  [{label}] 지수 피팅: |ρ| ~ {popt_e[0]:.4f}·exp(-k/{popt_e[1]:.3f}),  R²={r2_e:.4f}")
    except Exception as e:
        log(f"\n  [{label}] 지수 피팅 실패: {e}")
        results_fit['exp'] = None

    # 멱법칙 피팅: A * k^{-alpha}
    try:
        def power_law(k, A, alpha):
            return A * k ** (-alpha)
        p0_pow = [r_fit[0], 1.0]
        popt_p, _ = curve_fit(power_law, k_fit, r_fit, p0=p0_pow,
                               bounds=([0, 0.01], [2.0, 10.0]), maxfev=5000)
        r_pred_p = power_law(k_fit, *popt_p)
        ss_res = np.sum((r_fit - r_pred_p)**2)
        ss_tot = np.sum((r_fit - r_fit.mean())**2)
        r2_p = 1 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
        results_fit['pow'] = {'A': popt_p[0], 'alpha': popt_p[1], 'R2': r2_p}
        log(f"  [{label}] 멱법칙 피팅: |ρ| ~ {popt_p[0]:.4f}·k^{{-{popt_p[1]:.3f}}},  R²={r2_p:.4f}")
    except Exception as e:
        log(f"  [{label}] 멱법칙 피팅 실패: {e}")
        results_fit['pow'] = None

    # 어느 모델이 더 나은가?
    if results_fit.get('exp') and results_fit.get('pow'):
        r2_e = results_fit['exp']['R2']
        r2_p = results_fit['pow']['R2']
        better = '지수' if r2_e >= r2_p else '멱법칙'
        log(f"  [{label}] 우세 모델: {better}  (지수 R²={r2_e:.4f}, 멱법칙 R²={r2_p:.4f})")

        # 특성 길이 판정 (지수 모델 기준)
        if better == '지수':
            xi = results_fit['exp']['xi']
            if xi > 2.0:
                log(f"  [{label}] ★ 특성 길이 ξ={xi:.3f} > 2 → 장거리 상관 확인")
            elif xi > 1.0:
                log(f"  [{label}] △ 특성 길이 ξ={xi:.3f} > 1 → 단거리 비국소성")
            else:
                log(f"  [{label}] ξ={xi:.3f} ≤ 1 → nearest-neighbor 고유")

    return results_fit


# ═══════════════════════════════════════════════════════════════════════════════
# 메인 실행
# ═══════════════════════════════════════════════════════════════════════════════

log("=" * 70)
log("비국소 도달 범위 측정 — κ_mid ↔ gap_{n+k} 상관 감쇠")
log(f"시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log(f"dps(ζ)={DPS_ZETA}, dps(Λ 탐색)={DPS_DIR_S}, dps(Λ κ)={DPS_DIR_K}")
log(f"k 범위: 전방 0..{K_MAX}, 후방 1..{K_MAX}")
log("=" * 70)

all_results = []
csv_rows = []   # (func_label, direction, k, rho, p, partial_rho, n_used)


# ────────────────────────────────────────────────────────────────────────────
# 1. 리만 ζ
# ────────────────────────────────────────────────────────────────────────────
log("\n" + "=" * 70)
log("1. 리만 ζ — 영점 수집 (k=1..235, t∈[0,450])")
log("=" * 70)

mpmath.mp.dps = DPS_ZETA
zeta_zeros = []
for k in range(1, 236):
    try:
        z = mpmath.zetazero(k)
        t = float(z.imag)
        zeta_zeros.append(t)
    except Exception as e:
        print(f"WARNING: zetazero({k}) 실패: {e}", flush=True)

zeta_zeros = sorted(zeta_zeros)
log(f"  수집된 영점 수: {len(zeta_zeros)}")
log(f"  t 범위: {zeta_zeros[0]:.3f} ~ {zeta_zeros[-1]:.3f}")

if len(zeta_zeros) == 0:
    print("⚠️ 영점 0개 — 탐색 로직 점검 필요")
    sys.exit(1)

res_zeta = compute_reach_analysis(zeta_zeros, kappa_zeta, "리만 ζ", char_info=None)
if res_zeta:
    all_results.append(res_zeta)
    for row in res_zeta['forward']:
        csv_rows.append(("ζ", "forward", *row))
    for row in res_zeta['backward']:
        csv_rows.append(("ζ", "backward", *row))

    # 감쇠 피팅 (전방 k=1..K_MAX)
    log("\n  [ζ] 전방 감쇠 피팅:")
    fwd_k   = [r[0] for r in res_zeta['forward']]
    fwd_rho = [r[1] for r in res_zeta['forward']]
    fit_zeta_fwd = fit_decay(fwd_k, fwd_rho, "ζ 전방")

    # 전방/후방 비교
    log("\n  [ζ] 전방 vs 후방 |ρ| 비교:")
    log(f"  {'k':>3}  {'|ρ_fwd|':>9}  {'|ρ_bwd|':>9}  {'비대칭?':>10}")
    for k in range(1, K_MAX + 1):
        fwd_entry = next((r for r in res_zeta['forward']  if r[0] == k), None)
        bwd_entry = next((r for r in res_zeta['backward'] if r[0] == k), None)
        if fwd_entry and bwd_entry:
            rf = abs(fwd_entry[1]) if not np.isnan(fwd_entry[1]) else float('nan')
            rb = abs(bwd_entry[1]) if not np.isnan(bwd_entry[1]) else float('nan')
            asym = "비대칭" if not np.isnan(rf) and not np.isnan(rb) and abs(rf - rb) > 0.1 else "대칭"
            log(f"  {k:>3}  {rf:>9.4f}  {rb:>9.4f}  {asym:>10}")


# ────────────────────────────────────────────────────────────────────────────
# 2. 디리클레 L-함수
# ────────────────────────────────────────────────────────────────────────────

DIRICHLET_CHARS = {
    'χ₃ (mod 3)': {'chi': [0, 1, -1],        'q': 3, 'a': 1},
    'χ₄ (mod 4)': {'chi': [0, 1, 0, -1],     'q': 4, 'a': 1},
    'χ₅ (mod 5)': {'chi': [0, 1, -1, -1, 1], 'q': 5, 'a': 0},
}

for char_label, char_info in DIRICHLET_CHARS.items():
    log("\n" + "=" * 70)
    log(f"2. 디리클레 {char_label}")
    log("=" * 70)

    log(f"  영점 탐색: t∈[1,200], n_scan=3000, dps={DPS_DIR_S} ...")
    t0 = time.time()
    dir_zeros, fail_cnt = find_zeros_L(char_info, t_min=1.0, t_max=200.0, n_scan=3000)
    log(f"  탐색 완료: {len(dir_zeros)}개 영점, 실패={fail_cnt}, {time.time()-t0:.1f}초")

    if len(dir_zeros) == 0:
        print(f"⚠️ {char_label} 영점 0개 — 탐색 로직 점검 필요")
        continue

    res_dir = compute_reach_analysis(dir_zeros, kappa_dirichlet, char_label, char_info=char_info)
    if res_dir:
        all_results.append(res_dir)
        for row in res_dir['forward']:
            csv_rows.append((char_label, "forward", *row))
        for row in res_dir['backward']:
            csv_rows.append((char_label, "backward", *row))

        # 감쇠 피팅 (전방)
        log(f"\n  [{char_label}] 전방 감쇠 피팅:")
        fwd_k   = [r[0] for r in res_dir['forward']]
        fwd_rho = [r[1] for r in res_dir['forward']]
        fit_decay(fwd_k, fwd_rho, f"{char_label} 전방")

        # 전방/후방 비교
        log(f"\n  [{char_label}] 전방 vs 후방 |ρ| 비교:")
        log(f"  {'k':>3}  {'|ρ_fwd|':>9}  {'|ρ_bwd|':>9}  {'비대칭?':>10}")
        for k in range(1, K_MAX + 1):
            fwd_entry = next((r for r in res_dir['forward']  if r[0] == k), None)
            bwd_entry = next((r for r in res_dir['backward'] if r[0] == k), None)
            if fwd_entry and bwd_entry:
                rf = abs(fwd_entry[1]) if not np.isnan(fwd_entry[1]) else float('nan')
                rb = abs(bwd_entry[1]) if not np.isnan(bwd_entry[1]) else float('nan')
                asym = "비대칭" if not np.isnan(rf) and not np.isnan(rb) and abs(rf - rb) > 0.1 else "대칭"
                log(f"  {k:>3}  {rf:>9.4f}  {rb:>9.4f}  {asym:>10}")

    save_all()


# ────────────────────────────────────────────────────────────────────────────
# 3. 최종 종합 판정
# ────────────────────────────────────────────────────────────────────────────
log("\n" + "=" * 70)
log("최종 종합 판정")
log("=" * 70)

log(f"\n  분석된 L-함수: {len(all_results)}개")

# k=1에서 #18 재현 확인
log("\n  k=1 재현 (기존 #18 ρ(ζ)=-0.654):")
for res in all_results:
    fwd_k1 = next((r for r in res['forward'] if r[0] == 1), None)
    if fwd_k1 and not np.isnan(fwd_k1[1]):
        reprod = "✅ 재현" if abs(fwd_k1[1]) > 0.5 else ("△ 약함" if abs(fwd_k1[1]) > 0.3 else "❌ 재현 실패")
        log(f"    {res['label']:20s}: ρ(k=1)={fwd_k1[1]:.4f}  {reprod}")

# 감쇠 패턴 요약
log("\n  감쇠 패턴 요약:")
for res in all_results:
    k_vals  = [r[0] for r in res['forward'] if r[0] >= 1 and not np.isnan(r[1])]
    rho_abs = [abs(r[1]) for r in res['forward'] if r[0] >= 1 and not np.isnan(r[1])]
    if len(rho_abs) >= 3:
        monotone = all(rho_abs[i] >= rho_abs[i+1] for i in range(len(rho_abs)-1))
        k2_sig   = rho_abs[1] > 0.15 if len(rho_abs) > 1 else False
        k3_sig   = rho_abs[2] > 0.10 if len(rho_abs) > 2 else False
        pattern  = "단조감소" if monotone else "진동"
        reach    = "장거리(k≥3)" if k3_sig else ("중거리(k=2)" if k2_sig else "단거리(k=1)")
        log(f"    {res['label']:20s}: 패턴={pattern}, 도달={reach}, |ρ(k=1..3)|={[f'{x:.3f}' for x in rho_abs[:3]]}")

# 전방/후방 대칭성
log("\n  전방/후방 대칭성 (k=1 기준):")
for res in all_results:
    fwd1 = next((r for r in res['forward']  if r[0] == 1), None)
    bwd1 = next((r for r in res['backward'] if r[0] == 1), None)
    if fwd1 and bwd1 and not np.isnan(fwd1[1]) and not np.isnan(bwd1[1]):
        rf, rb = abs(fwd1[1]), abs(bwd1[1])
        sym = "대칭 (시간반전 대칭 시사)" if abs(rf - rb) < 0.05 else f"비대칭 (|ρ_fwd|-|ρ_bwd|={rf-rb:+.3f})"
        log(f"    {res['label']:20s}: |ρ_fwd|={rf:.4f}, |ρ_bwd|={rb:.4f} → {sym}")

# 성공 기준 달성 여부
log("\n  성공 기준 평가:")
success_k1  = any(
    next((r for r in res['forward'] if r[0] == 1), None) is not None and
    abs(next((r for r in res['forward'] if r[0] == 1), (0,0))[1]) > 0.5
    for res in all_results if res['label'] == '리만 ζ'
)
success_k2  = any(
    next((r for r in res['forward'] if r[0] == 2), None) is not None and
    abs(next((r for r in res['forward'] if r[0] == 2), (0,0))[1]) > 0.15
    for res in all_results
)
log(f"  ✅ k=1에서 |ρ| > 0.5 (기존 #18 재현): {'✅ PASS' if success_k1 else '❌ FAIL'}")
log(f"  k=2 이상에서 감쇠 패턴 관측: {'✅ 있음' if success_k2 else '❌ 없음 (k=1 고유)'}")

if success_k1:
    if success_k2:
        verdict = "★ 양성 — 비국소 도달 범위 k≥2 확인. 특성 길이 ξ 정의 가능."
    else:
        verdict = "△ 조건부 양성 — k=1 재현. k≥2에서 빠른 감쇠 (nearest-neighbor+1 고유)."
else:
    verdict = "❌ k=1 재현 실패 — 인프라 점검 필요"

log(f"\n  최종 판정: {verdict}")
log()

# ── CSV 저장 ──────────────────────────────────────────────────────────────────
csv_header = "func,direction,k,rho,p,partial_rho,n_used\n"
csv_body   = ""
for row in csv_rows:
    func, direction, k, rho, p, partial_rho, n_used = row
    csv_body += f"{func},{direction},{k},{rho:.6f},{p:.4e},{partial_rho:.6f},{n_used}\n"

with open(OUTPUT_CSV, 'w', encoding='utf-8') as f:
    f.write(csv_header + csv_body)
log(f"\n  CSV 저장: {OUTPUT_CSV}")

save_all()
log(f"\n결과 저장: {OUTPUT_TXT}")
log(f"완료: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
