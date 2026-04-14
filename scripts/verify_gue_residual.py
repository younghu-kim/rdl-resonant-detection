#!/usr/bin/env python3
"""
독립 검증: 결과 #34 (Δκ 잔차 vs GUE NNS 상관)
================================================
목적: curvature_gue_residual.py / results/curvature_gue_residual.txt 결과를
      완전히 독립적으로 spot-check 및 통계 재현.

검증 항목:
  1. spot-check 3개 데이터 포인트 (N=203, 5253, 16767)
  2. OLS R²=0.0152 재현
  3. Pearson r(잔차, NNS_unfolded) = -0.4735 재현
  4. 2변수 R²=0.9185 (Lorentzian² 모델) 재현
  5. Lorentzian² 공식 차원·이론 점검
"""

import sys, os
import numpy as np
import mpmath
from scipy import stats as scipy_stats

# ─────────────────────────────────────────────────────
# 정밀도 설정
# ─────────────────────────────────────────────────────
def get_dps(t):
    if t > 30000: return 120
    if t > 10000: return 100
    if t > 5000:  return 80
    return max(50, int(30 + t / 200))


# ─────────────────────────────────────────────────────
# 수식 구현 (원본과 독립)
# ─────────────────────────────────────────────────────
def G_component(s):
    """G(s) = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2"""
    return (mpmath.mpf(1)/s
            + mpmath.mpf(1)/(s - 1)
            - mpmath.log(mpmath.pi)/2
            + mpmath.digamma(s/2)/2)


def zeta_log_deriv_cd(s, h=None):
    """ζ'/ζ(s) — h=1e-6 중앙 차분 (원본과 동일 접근)"""
    if h is None:
        h = mpmath.mpf('1e-6')
    z = mpmath.zeta(s)
    if abs(z) < mpmath.mpf(10)**(-(mpmath.mp.dps - 10)):
        return mpmath.mpc(1e8, 0)
    zd = (mpmath.zeta(s + h) - mpmath.zeta(s - h)) / (2 * h)
    return zd / z


def xi_log_deriv(s):
    """ξ'/ξ(s) = G(s) + ζ'/ζ(s)"""
    return G_component(s) + zeta_log_deriv_cd(s)


def compute_all(t, delta=0.03):
    """
    t : 영점 허수부
    delta : δ 오프셋 (0.03)
    반환: dict with t, s, xi_ld, kappa, kappa_base, delta_kappa, R, Re_R, Im_R
    """
    dps = get_dps(t)
    with mpmath.workdps(dps):
        s = mpmath.mpc(0.5 + delta, t)
        xi_ld = xi_log_deriv(s)
        kappa = float(abs(xi_ld)**2)
        kappa_base = 1.0 / delta**2
        delta_kappa = kappa - kappa_base
        R = xi_ld - mpmath.mpf(1) / delta   # R = ξ'/ξ - 1/δ
        Re_R = float(mpmath.re(R))
        Im_R = float(mpmath.im(R))
    return dict(t=t, dps=dps, xi_ld=xi_ld, kappa=kappa,
                kappa_base=kappa_base, delta_kappa=delta_kappa,
                R=R, Re_R=Re_R, Im_R=Im_R)


def get_zero(n, dps=100):
    """n번째 리만 영점 허수부"""
    with mpmath.workdps(dps):
        z = mpmath.zetazero(n)
        return float(mpmath.im(z))


# ─────────────────────────────────────────────────────
# 보고된 값 (results/curvature_gue_residual.txt)
# ─────────────────────────────────────────────────────
REPORTED = {
    203:   dict(t=401.8392,   dk=6.3285,    Re_R=0.0558,  Im_R=1.6138,
                NNS_unf=0.6768, NNS_min=1.0227, gap_L=1.8541, gap_R=1.0227),
    5253:  dict(t=5682.0992,  dk=83.4610,   Re_R=0.9047,  Im_R=4.7255,
                NNS_unf=0.2047, NNS_min=0.1890, gap_L=0.8576, gap_R=0.1890),
    16767: dict(t=15471.5508, dk=568.0846,  Re_R=5.5320,  Im_R=12.9879,
                NNS_unf=0.0843, NNS_min=0.0678, gap_L=1.3427, gap_R=0.0678),
}

# ─────────────────────────────────────────────────────
# 파트 1: Spot-check 3 data points
# ─────────────────────────────────────────────────────
def part1_spotcheck():
    print()
    print("=" * 70)
    print("파트 1: Spot-check — N=203, 5253, 16767")
    print("=" * 70)

    DELTA = 0.03
    results_spot = {}

    for N in [203, 5253, 16767]:
        rep = REPORTED[N]
        print(f"\n  ─── N = {N} ───")

        # 1a. zetazero
        dps = get_dps(rep['t'])
        t_n = get_zero(N, dps=dps)
        t_err = abs(t_n - rep['t'])
        print(f"  t_n 계산값    : {t_n:.6f}")
        print(f"  t_n 보고값    : {rep['t']:.4f}")
        print(f"  |오차|        : {t_err:.2e}   {'OK' if t_err < 0.01 else 'FAIL'}")

        # 1b. s = 1/2 + δ + i*t_n
        s_val = complex(0.5 + DELTA, t_n)
        print(f"  s             : {0.5+DELTA:.4f} + {t_n:.4f}i")

        # 1c. G(s) + ζ'/ζ(s) → ξ'/ξ
        res = compute_all(t_n, delta=DELTA)
        xi_ld_complex = complex(mpmath.re(res['xi_ld']), mpmath.im(res['xi_ld']))

        print(f"  ξ'/ξ(s)      : {xi_ld_complex.real:.6f} + {xi_ld_complex.imag:.6f}i")

        # 1d. κ, Δκ, R, Im[R]
        kappa      = res['kappa']
        dk         = res['delta_kappa']
        Re_R       = res['Re_R']
        Im_R       = res['Im_R']

        dk_err  = abs(dk  - rep['dk'])
        ReR_err = abs(Re_R - rep['Re_R'])
        ImR_err = abs(Im_R - rep['Im_R'])

        print(f"  κ             : {kappa:.6f}")
        print(f"  1/δ²          : {1/DELTA**2:.4f}")
        print(f"  Δκ 계산값     : {dk:.4f}   보고값: {rep['dk']:.4f}   오차: {dk_err:.4f}  {'OK' if dk_err < 0.5 else 'FAIL'}")
        print(f"  Re[R] 계산값  : {Re_R:.4f}   보고값: {rep['Re_R']:.4f}   오차: {ReR_err:.4f}  {'OK' if ReR_err < 0.01 else 'FAIL'}")
        print(f"  Im[R] 계산값  : {Im_R:.4f}   보고값: {rep['Im_R']:.4f}   오차: {ImR_err:.4f}  {'OK' if ImR_err < 0.02 else 'FAIL'}")

        # 1e. NNS 계산
        t_prev = get_zero(N - 1, dps=min(dps, 80))
        t_next = get_zero(N + 1, dps=min(dps, 80))
        gap_L = t_n - t_prev
        gap_R = t_next - t_n
        NNS_min = min(gap_L, gap_R)
        density = np.log(t_n / (2 * np.pi)) / (2 * np.pi)
        NNS_unf = NNS_min * density

        nns_err = abs(NNS_unf - rep['NNS_unf'])
        print(f"  gap_L         : {gap_L:.4f}   보고값: {rep['gap_L']:.4f}")
        print(f"  gap_R         : {gap_R:.4f}   보고값: {rep['gap_R']:.4f}")
        print(f"  NNS_min       : {NNS_min:.4f}   보고값: {rep['NNS_min']:.4f}")
        print(f"  density D(t)  : {density:.6f}")
        print(f"  NNS_unfolded  : {NNS_unf:.4f}   보고값: {rep['NNS_unf']:.4f}   오차: {nns_err:.4f}  {'OK' if nns_err < 0.005 else 'FAIL'}")

        results_spot[N] = dict(t=t_n, dk=dk, Re_R=Re_R, Im_R=Im_R,
                                NNS_min=NNS_min, NNS_unf=NNS_unf,
                                gap_L=gap_L, gap_R=gap_R)

    return results_spot


# ─────────────────────────────────────────────────────
# 파트 2: 통계 재현 (OLS R², Pearson r, 2변수 R²)
# ─────────────────────────────────────────────────────

def parse_result_file():
    """결과 파일에서 99개 데이터 파싱"""
    fpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '..', 'results', 'curvature_gue_residual.txt')
    Ns, ts, dks, nns_mins, nns_unfs = [], [], [], [], []
    Re_Rs, Im_Rs = [], []

    with open(fpath, 'r') as f:
        for line in f:
            line = line.strip()
            # 데이터 행: 숫자로 시작하는 고정폭 행
            # 형식: "   1    203    401.8392    6.3285    1.8541    1.0227    0.6768    0.0558    1.6138  ✅"
            parts = line.split()
            if len(parts) >= 9 and parts[0].isdigit():
                try:
                    idx   = int(parts[0])
                    N     = int(parts[1])
                    t     = float(parts[2])
                    dk    = float(parts[3])
                    gL    = float(parts[4])
                    gR    = float(parts[5])
                    nns_u = float(parts[6])
                    Re_R  = float(parts[7])
                    Im_R  = float(parts[8])
                    # NNS_min = min(gL, gR)
                    nns_m = min(gL, gR)
                    Ns.append(N); ts.append(t); dks.append(dk)
                    nns_mins.append(nns_m); nns_unfs.append(nns_u)
                    Re_Rs.append(Re_R); Im_Rs.append(Im_R)
                except (ValueError, IndexError):
                    pass

    return (np.array(Ns), np.array(ts), np.array(dks),
            np.array(nns_mins), np.array(nns_unfs),
            np.array(Re_Rs), np.array(Im_Rs))


def part2_statistics():
    print()
    print("=" * 70)
    print("파트 2: 통계 재현 (결과 파일 데이터 사용)")
    print("=" * 70)

    Ns, ts, dks, nns_mins, nns_unfs, Re_Rs, Im_Rs = parse_result_file()
    print(f"  파싱된 데이터 포인트 수: {len(Ns)}")

    # ── 2a. 1변수 OLS: Δκ ~ a·log(t/2π) + b ──
    log_t = np.log(ts / (2 * np.pi))
    X1 = np.column_stack([log_t, np.ones(len(log_t))])
    # OLS: β = (X'X)^{-1} X'y
    beta1, res1, rank1, sv1 = np.linalg.lstsq(X1, dks, rcond=None)
    a1, b1 = beta1
    dk_pred1 = X1 @ beta1
    ss_res1 = np.sum((dks - dk_pred1)**2)
    ss_tot  = np.sum((dks - np.mean(dks))**2)
    R2_1var = 1 - ss_res1 / ss_tot
    residuals = dks - dk_pred1

    print(f"\n  1변수 OLS: Δκ ~ a·log(t/2π) + b")
    print(f"    a      = {a1:.4f}   보고값: 9.6846")
    print(f"    b      = {b1:.4f}   보고값: -45.1474")
    print(f"    R²     = {R2_1var:.4f}   보고값: 0.0152   {'OK' if abs(R2_1var-0.0152)<0.001 else 'FAIL'}")

    # ── 2b. Pearson r(residuals, NNS_unfolded) ──
    r_pearson, p_pearson = scipy_stats.pearsonr(residuals, nns_unfs)
    print(f"\n  res_n vs NNS_unfolded:")
    print(f"    Pearson r  = {r_pearson:.4f}   보고값: -0.4735   {'OK' if abs(r_pearson-(-0.4735))<0.002 else 'FAIL'}")
    print(f"    p-value    = {p_pearson:.4e}   보고값: 7.4368e-07")

    r_spearman, p_spearman = scipy_stats.spearmanr(residuals, nns_unfs)
    print(f"    Spearman ρ = {r_spearman:.4f}   보고값: -0.7387   {'OK' if abs(r_spearman-(-0.7387))<0.002 else 'FAIL'}")

    # ── 2c. 2변수 R² — Lorentzian² 모델 ──
    # f₁ = δ²/(δ² + NNS_min²)²
    DELTA = 0.03
    f1 = DELTA**2 / (DELTA**2 + nns_mins**2)**2
    X2 = np.column_stack([log_t, f1, np.ones(len(log_t))])
    beta2, res2, rank2, sv2 = np.linalg.lstsq(X2, dks, rcond=None)
    dk_pred2 = X2 @ beta2
    ss_res2 = np.sum((dks - dk_pred2)**2)
    R2_2var = 1 - ss_res2 / ss_tot

    print(f"\n  2변수 OLS: Δκ ~ a·log(t/2π) + b + c·f₁(NNS_min)")
    print(f"  f₁ = δ²/(δ²+NNS_min²)²  (Lorentzian²)")
    print(f"    R²_2var    = {R2_2var:.4f}   보고값: 0.9185   {'OK' if abs(R2_2var-0.9185)<0.002 else 'FAIL'}")
    print(f"    ΔR²        = {R2_2var - R2_1var:.4f}   보고값: +0.9033")

    # F-통계량
    n = len(dks)
    p_full = 3
    p_restr = 2
    F_stat = ((ss_res1 - ss_res2) / (p_full - p_restr)) / (ss_res2 / (n - p_full))
    print(f"    F-통계량   = {F_stat:.3f}   보고값: 1064.239   {'OK' if abs(F_stat-1064.239)/1064.239<0.01 else 'CHECK'}")

    # res vs Lorentzian² = f₁ (원본 스크립트 line 281: (delta/(delta²+NNS_min²))² 계산)
    # 결과 파일 레이블이 "1/NNS_min²"로 잘못 표기되었으나 실제로는 f₁ Lorentzian²임
    # Pearson r은 스케일 불변이므로 (delta/(delta²+NNS²))² 와 (delta²/(delta²+NNS²))² 동일
    lorentz2_partB = (DELTA / (DELTA**2 + nns_mins**2))**2
    r_lor, p_lor = scipy_stats.pearsonr(residuals, lorentz2_partB)
    print(f"\n  res_n vs Lorentzian² (원본 Part B 실제 계산값):")
    print(f"    Pearson r  = {r_lor:.4f}   보고값: 0.9534   {'OK' if abs(r_lor-0.9534)<0.001 else 'FAIL'}")
    print(f"  [주의] 결과 파일 레이블 '1/NNS_min²'는 오기 — 실제로는 Lorentzian² f₁")
    print(f"         문자 그대로 1/NNS_min² 계산 시 r={scipy_stats.pearsonr(residuals,1/nns_mins**2)[0]:.4f} (다름)")

    return dict(R2_1var=R2_1var, r_pearson=r_pearson, p_pearson=p_pearson,
                R2_2var=R2_2var, F_stat=F_stat,
                Ns=Ns, ts=ts, dks=dks, nns_mins=nns_mins, nns_unfs=nns_unfs,
                residuals=residuals)


# ─────────────────────────────────────────────────────
# 파트 3: Lorentzian² 공식 이론 점검
# ─────────────────────────────────────────────────────

def part3_lorentzian_check():
    print()
    print("=" * 70)
    print("파트 3: Lorentzian² 공식 점검")
    print("=" * 70)
    DELTA = 0.03

    # 공식: f₁ = δ²/(δ² + NNS_min²)²
    # 차원 점검: δ와 NNS_min은 둘 다 무차원 (unfolded 좌표 또는 실수 간격)
    # 여기서 NNS_min은 원 간격 (원 단위, Hz 아님) — 차원 일치.
    print(f"\n  f₁ = δ²/(δ²+NNS_min²)²,  δ={DELTA}")
    print(f"  차원 분석:")
    print(f"    분자: δ² = {DELTA**2:.6f}  [무차원²]")
    print(f"    분모: (δ²+NNS_min²)² = [...] [무차원⁴]")
    print(f"    f₁ 단위: [무차원²/무차원⁴] = [무차원^-2]")
    print(f"    → f₁은 순수 수치 (무차원) — 단, 분모가 4승이므로 단위 상쇄 확인 필요")
    print(f"    주의: δ는 복소 오프셋(σ 방향), NNS_min은 t 방향 간격")
    print(f"          두 양이 같은 단위(실수 길이)로 취급 — 이론적 정당화는 별도 필요")

    print(f"\n  특수값 계산:")
    for nns in [0.01, 0.03, 0.06, 0.1, 0.2, 0.5, 1.0]:
        f1 = DELTA**2 / (DELTA**2 + nns**2)**2
        print(f"    NNS_min={nns:.3f}:  f₁={f1:.4f}")

    # NNS_min = δ일 때 최대? 미분으로 확인
    # df/d(nns) = δ² * (-2)*(δ²+nns²) * 2*nns / (δ²+nns²)^4 = -4δ²*nns/(δ²+nns²)^3
    # 음수 → nns 증가시 f₁ 단조 감소 (nns>0 범위에서)
    # nns→0 일 때 f₁→1/δ² = 1111.1
    f1_limit_0 = 1.0 / DELTA**2
    f1_at_delta = DELTA**2 / (2*DELTA**2)**2
    print(f"\n  극한·특이점:")
    print(f"    NNS_min→0:  f₁→ 1/δ² = {f1_limit_0:.4f}   (영점이 거의 겹칠 때 발산)")
    print(f"    NNS_min=δ:  f₁ = {f1_at_delta:.6f}")
    print(f"    NNS_min→∞: f₁→ 0     (영점 간격이 크면 기여 없음)")
    print(f"  → 단조감소, NNS_min≈0(영점 군집)에서 Δκ 스파이크 설명 가능")
    print(f"  → 이론: ξ'/ξ 극점(유사-영점)이 NNS_min≈0 근방에서 s와 가까워져 κ 폭발")
    print(f"  → Lorentzian² 형태는 Breit-Wigner/Cauchy 분포의 제곱으로 표준적")

    # 이론 근거: ξ'/ξ의 Hadamard 전개에서 각 영점 ρ_k의 기여 ~ 1/(s-ρ_k)
    # δ=0.03 오프셋에서 s=σ+it_n, ρ_k=1/2+it_k
    # 가장 가까운 영점(NNS_min~gap)의 기여: 1/(δ + i*(t_n-t_k))
    # |기여|² ~ 1/(δ²+(t_n-t_k)²) = Lorentzian
    # κ = |합|² 에서 교차항 무시 시 ~ 1/(δ²+NNS²)²? → 정확히는 아님
    # 실제로는 κ=|ξ'/ξ|² 전체이고 f₁은 회귀 기저함수로 선택된 것
    print(f"\n  이론적 해석:")
    print(f"    Hadamard 전개: ξ'/ξ(s) ~ Σ_ρ [1/(s-ρ) + 1/(s-ρ̄)]")
    print(f"    가장 가까운 영점(ρ_k=1/2+it_k) 기여: ~1/(δ+i·Δt), |기여|=1/√(δ²+Δt²)")
    print(f"    Δκ ~ κ - 1/δ² ≈ 주도항 기여 + 교차항")
    print(f"    f₁=δ²/(δ²+NNS²)² 는 Lorentzian의 제곱 → 피팅 기저함수로 타당")
    print(f"    단, NNS_min을 Δt 대리변수로 사용하는 근사 포함")


# ─────────────────────────────────────────────────────
# 파트 4: 종합 판정
# ─────────────────────────────────────────────────────

def part4_summary(spot_ok, stat_res):
    print()
    print("=" * 70)
    print("파트 4: 종합 판정")
    print("=" * 70)

    r2_1 = stat_res['R2_1var']
    r_p  = stat_res['r_pearson']
    r2_2 = stat_res['R2_2var']

    checks = [
        ("Spot-check 통과",                           spot_ok),
        ("R²_1var ≈ 0.0152 (±0.001)",                abs(r2_1 - 0.0152) < 0.001),
        ("Pearson r ≈ -0.4735 (±0.002)",             abs(r_p - (-0.4735)) < 0.002),
        ("R²_2var ≈ 0.9185 (±0.002)",                abs(r2_2 - 0.9185) < 0.002),
    ]

    all_pass = True
    for name, ok in checks:
        sym = "PASS" if ok else "FAIL"
        print(f"  [{sym}] {name}")
        if not ok:
            all_pass = False

    print()
    if all_pass:
        print("  ★ 전체 검증 통과: 결과 #34 독립 재현 성공")
    else:
        print("  ✗ 일부 항목 불일치 — 상세 내용 확인 필요")

    print(f"\n  실측값 요약:")
    print(f"    R²_1var  = {r2_1:.4f}  (보고: 0.0152)")
    print(f"    r_pearson = {r_p:.4f}  (보고: -0.4735)")
    print(f"    R²_2var  = {r2_2:.4f}  (보고: 0.9185)")


# ─────────────────────────────────────────────────────
# 메인
# ─────────────────────────────────────────────────────

if __name__ == '__main__':
    print("독립 검증: 결과 #34 Δκ 잔차 vs GUE NNS 상관")
    print(f"mpmath version: {mpmath.__version__}")

    # 파트 1: Spot-check
    spot_results = part1_spotcheck()

    # Spot-check 통과 여부 판단
    spot_ok = True
    for N, rep in REPORTED.items():
        if N not in spot_results:
            spot_ok = False
            continue
        calc = spot_results[N]
        if abs(calc['dk'] - rep['dk']) > 1.0:      # Δκ 오차 1 이내
            spot_ok = False
        if abs(calc['Im_R'] - rep['Im_R']) > 0.05:  # Im[R] 오차 0.05 이내
            spot_ok = False
        if abs(calc['NNS_unf'] - rep['NNS_unf']) > 0.01:
            spot_ok = False

    # 파트 2: 통계 재현
    stat_res = part2_statistics()

    # 파트 3: Lorentzian² 이론 점검
    part3_lorentzian_check()

    # 파트 4: 종합 판정
    part4_summary(spot_ok, stat_res)
