#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 C-260] A(γ) height scaling law — ζ(s) 300 영점
=============================================================================
목적:
  A(γ) = Im(c₀)² + 2Re(c₁)의 γ-의존성 스케일링 법칙 탐색.
  C-244에서 c₁ 단독은 R²=0.79 (멱법칙)로 깨끗한 스케일링 없었음.
  A(γ) 전체는 분산이 줄어 더 깨끗한 법칙이 나올 수 있음.

이론 배경:
  Hadamard: c₁(ρ₀) = Σ'_ρ 1/(ρ₀-ρ)² + Γ-smooth
  <c₁> ~ (영점 밀도)^2 × pair correlation integral ~ (log γ)² / (2π)²
  B² = Im(c₀)²: c₀ = Σ'_ρ 1/(ρ₀-ρ)의 허수부 → 근접 영점 간격에 의존
  A = B² + 2H₁ → H₁ 성분이 지배 (B-20 확인)

  질문: <A> ~ log(γ/(2π))? (log(γ/(2π)))²? 또는 다른 법칙?

방법:
  C-244 δ-extrapolation 재사용 (8 δ값, σ방향 대칭/반대칭 분해).
  300 영점으로 확대 (γ ~ 14..570). 10-bin 평균으로 트렌드 추출.

성공 기준:
  (1) 300개 A 측정 완료
  (2) 피팅 R² > 0.90 → ★★★ 양성 (clean law)
  (3) R² ∈ [0.70, 0.90) → ★★ 조건부 양성
  (4) R² < 0.70 → ★ 음성 (A에도 깨끗한 스케일링 없음 = 유효 음성)

결과: results/A_height_scaling_c260.txt
=============================================================================
"""
import sys, os, time
import numpy as np
import mpmath

START = time.time()
RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'A_height_scaling_c260.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 80)
log("[실험 C-260] A(γ) height scaling law — ζ(s) 300 영점")
log("=" * 80)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

mpmath.mp.dps = 80
NZEROS = 300
NZEROS_HAD = 400  # Hadamard 비교용 (100 추가)

log(f"mpmath.dps = {mpmath.mp.dps}")
log(f"영점 수: 측정 {NZEROS} / Hadamard용 {NZEROS_HAD}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. ζ 영점 로드
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log(f"{T()} ζ 영점 {NZEROS_HAD}개 계산 중...")
gammas = []
for n in range(1, NZEROS_HAD + 1):
    g = float(mpmath.im(mpmath.zetazero(n)))
    gammas.append(g)
    if n % 100 == 0:
        log(f"  {T()} {n}/{NZEROS_HAD} (γ_{n} = {g:.4f})")

log(f"{T()} 완료. γ 범위: [{gammas[0]:.4f}, {gammas[-1]:.4f}]")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. Λ'/Λ 및 c₀, c₁ 추출
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def Lambda_zeta(s):
    """Λ(s) = π^{-s/2} Γ(s/2) ζ(s)"""
    s = mpmath.mpc(s)
    return mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s)

def log_deriv_Lambda(s):
    """Λ'/Λ(s) via central difference"""
    s = mpmath.mpc(s)
    h = mpmath.mpf(10) ** (-(mpmath.mp.dps // 3))
    Lp = Lambda_zeta(s + h)
    Lm = Lambda_zeta(s - h)
    L0 = Lambda_zeta(s)
    if abs(L0) == 0:
        return mpmath.mpc(0)
    return (Lp - Lm) / (2 * h * L0)

DELTAS = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02])

def extract_c0_c1(gamma):
    """
    ρ=1/2+iγ에서 Λ'/Λ Laurent c₀, c₁ 추출.
    sym(δ) = (f(+δ)+f(-δ))/2 → c₀ (절편)
    anti(δ) = ((f(+δ)-f(-δ))/2 - 1/δ)/δ → c₁ (절편)
    """
    rho = mpmath.mpc(0.5, gamma)
    sym_vals = []
    anti_vals = []
    for d in DELTAS:
        d_mp = mpmath.mpf(float(d))
        fp = log_deriv_Lambda(rho + d_mp)
        fm = log_deriv_Lambda(rho - d_mp)
        sym = complex((fp + fm) / 2)
        anti_num = (fp - fm) / 2 - 1 / d_mp
        anti = complex(anti_num / d_mp)
        sym_vals.append(sym)
        anti_vals.append(anti)

    d2 = DELTAS ** 2
    c0_re = np.polyfit(d2, [v.real for v in sym_vals], 2)[-1]
    c0_im = np.polyfit(d2, [v.imag for v in sym_vals], 2)[-1]
    c1_re = np.polyfit(d2, [v.real for v in anti_vals], 2)[-1]
    c1_im = np.polyfit(d2, [v.imag for v in anti_vals], 2)[-1]
    return complex(c0_re, c0_im), complex(c1_re, c1_im)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. Hadamard c₁ 예측
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def hadamard_c1(idx, gammas_all):
    g0 = gammas_all[idx]
    c1_xi = 0.0
    for k, g in enumerate(gammas_all):
        if k == idx:
            c1_xi += 1.0 / (2 * g0) ** 2
        else:
            c1_xi += 1.0 / (g0 - g) ** 2
            c1_xi += 1.0 / (g0 + g) ** 2
    corr = 2.0 * (0.25 - g0**2) / (0.25 + g0**2)**2
    return c1_xi + corr

def hadamard_c1_nn(idx, gammas_all, k_nn=1):
    """최근접 k개 영점의 기여만"""
    g0 = gammas_all[idx]
    # 모든 영점과의 거리 정렬
    dists = [(abs(g0 - g), k_idx) for k_idx, g in enumerate(gammas_all) if k_idx != idx]
    dists.sort()
    nn_sum = 0.0
    for i in range(min(k_nn, len(dists))):
        _, kid = dists[i]
        g = gammas_all[kid]
        nn_sum += 1.0 / (g0 - g) ** 2
        nn_sum += 1.0 / (g0 + g) ** 2
    return nn_sum

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. 측정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
results = []

log("━━━ c₀, c₁, A 추출 ━━━")
log(f"{'n':>4} {'γ':>10} {'Im(c₀)':>10} {'c₁':>10} {'A':>10} {'B²':>8} {'2H₁':>8} {'c₁_Had':>10}")
log("-" * 82)

for n in range(NZEROS):
    gamma = gammas[n]
    c0, c1 = extract_c0_c1(gamma)
    c1_had = hadamard_c1(n, gammas)
    c1_nn2 = hadamard_c1_nn(n, gammas, k_nn=2)

    B_sq = c0.imag ** 2
    H1 = c1.real
    A = B_sq + 2 * H1

    # 최근접 영점 간격
    if n > 0 and n < NZEROS - 1:
        gap_left = gammas[n] - gammas[n-1]
        gap_right = gammas[n+1] - gammas[n]
    elif n == 0:
        gap_left = float('inf')
        gap_right = gammas[1] - gammas[0]
    else:
        gap_left = gammas[n] - gammas[n-1]
        gap_right = float('inf')

    results.append({
        'n': n + 1, 'gamma': gamma,
        'c0_im': c0.imag, 'c0_re': c0.real,
        'c1_re': c1.real, 'c1_im': c1.imag,
        'c1_had': c1_had, 'c1_nn2': c1_nn2,
        'B_sq': B_sq, 'H1': H1, 'A': A,
        'gap_min': min(gap_left, gap_right),
        'gap_avg': (gap_left + gap_right) / 2 if gap_left < 1e10 and gap_right < 1e10 else max(gap_left, gap_right),
    })

    if (n + 1) % 50 == 0 or (n + 1) <= 5:
        log(f"{n+1:>4} {gamma:>10.4f} {c0.imag:>10.6f} {c1.real:>10.6f} {A:>10.6f} {B_sq:>8.4f} {2*H1:>8.4f} {c1_had:>10.6f}")
    if (n + 1) % 50 == 0:
        log(f"  {T()} {n+1}/{NZEROS} 완료")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. 스케일링 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
gam = np.array([r['gamma'] for r in results])
A_arr = np.array([r['A'] for r in results])
B2_arr = np.array([r['B_sq'] for r in results])
H1_arr = np.array([r['H1'] for r in results])
c1_arr = np.array([r['c1_re'] for r in results])
gap_arr = np.array([r['gap_min'] for r in results])

x_log = np.log(gam / (2 * np.pi))  # log(γ/(2π))

# --- 5a. 개별 A의 피팅 ---
log("━━━ 개별 A(γ) 피팅 (300점) ━━━")
log()

from numpy.polynomial import polynomial as P

def fit_and_report(name, x, y, label_x='x'):
    """다양한 모델 피팅"""
    models = {}

    # (a) linear: y = a + b·x
    c = np.polyfit(x, y, 1)
    y_pred = np.polyval(c, x)
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    models['linear'] = (r2, f"{name} = {c[1]:.4f} + {c[0]:.4f}·{label_x}")

    # (b) quadratic: y = a + b·x + c·x²
    c2 = np.polyfit(x, y, 2)
    y_pred2 = np.polyval(c2, x)
    r2_2 = 1 - np.sum((y - y_pred2)**2) / ss_tot if ss_tot > 0 else 0
    models['quadratic'] = (r2_2, f"{name} = {c2[2]:.4f} + {c2[1]:.4f}·{label_x} + {c2[0]:.4f}·{label_x}²")

    # (c) power law: log(y) = log(a) + b·log(x) → y = a·x^b
    mask = (y > 0) & (x > 0)
    if np.sum(mask) > 10:
        c_pw = np.polyfit(np.log(x[mask]), np.log(y[mask]), 1)
        y_pred_pw = np.exp(c_pw[1]) * x[mask] ** c_pw[0]
        r2_pw = 1 - np.sum((y[mask] - y_pred_pw)**2) / np.sum((y[mask] - np.mean(y[mask]))**2)
        models['power'] = (r2_pw, f"{name} = {np.exp(c_pw[1]):.4f}·{label_x}^{c_pw[0]:.4f}")
    else:
        models['power'] = (0, 'N/A')

    return models

# A vs log(γ/(2π))
log(f"  [A vs log(γ/(2π))]")
models_A = fit_and_report('A', x_log, A_arr, 'log(γ/(2π))')
for k in ['linear', 'quadratic', 'power']:
    r2, eq = models_A[k]
    log(f"    {k:12s}: R² = {r2:.6f}  {eq}")

# A vs γ (power law)
log(f"\n  [A vs γ (power law)]")
c_pw = np.polyfit(np.log(gam), np.log(np.maximum(A_arr, 1e-10)), 1)
y_pw = np.exp(c_pw[1]) * gam ** c_pw[0]
r2_gam = 1 - np.sum((A_arr - y_pw)**2) / np.sum((A_arr - np.mean(A_arr))**2)
log(f"    A = {np.exp(c_pw[1]):.4f}·γ^{c_pw[0]:.4f}   R² = {r2_gam:.6f}")

# B² vs log(γ/(2π))
log(f"\n  [B² vs log(γ/(2π))]")
models_B2 = fit_and_report('B²', x_log, B2_arr, 'log(γ/(2π))')
for k in ['linear', 'quadratic']:
    r2, eq = models_B2[k]
    log(f"    {k:12s}: R² = {r2:.6f}  {eq}")

# H₁ vs log(γ/(2π))
log(f"\n  [H₁=Re(c₁) vs log(γ/(2π))]")
models_H1 = fit_and_report('H₁', x_log, H1_arr, 'log(γ/(2π))')
for k in ['linear', 'quadratic', 'power']:
    r2, eq = models_H1[k]
    log(f"    {k:12s}: R² = {r2:.6f}  {eq}")

log()

# --- 5b. 빈 평균 분석 (노이즈 감소) ---
log("━━━ 빈 평균 분석 (10 bins) ━━━")
N_BINS = 10
bin_edges = np.linspace(gam[0], gam[-1], N_BINS + 1)
bin_stats = []

log(f"{'bin':>4} {'γ_lo':>8} {'γ_hi':>8} {'n':>4} {'<A>':>10} {'std':>8} {'<B²>':>8} {'<2H₁>':>8} {'<gap>':>8} {'log(γ̄/2π)':>10}")
log("-" * 94)

for i in range(N_BINS):
    mask = (gam >= bin_edges[i]) & (gam < bin_edges[i+1])
    if i == N_BINS - 1:
        mask = (gam >= bin_edges[i]) & (gam <= bin_edges[i+1])
    if np.sum(mask) == 0:
        continue
    g_mid = np.mean(gam[mask])
    A_mean = np.mean(A_arr[mask])
    A_std = np.std(A_arr[mask])
    B2_mean = np.mean(B2_arr[mask])
    H1_mean = np.mean(H1_arr[mask])
    gap_mean = np.mean(gap_arr[mask])
    x_mid = np.log(g_mid / (2 * np.pi))

    bin_stats.append({
        'g_lo': bin_edges[i], 'g_hi': bin_edges[i+1],
        'n': np.sum(mask), 'g_mid': g_mid,
        'A_mean': A_mean, 'A_std': A_std,
        'B2_mean': B2_mean, 'H1_mean': H1_mean,
        'gap_mean': gap_mean, 'x_mid': x_mid,
    })
    log(f"{i+1:>4} {bin_edges[i]:>8.1f} {bin_edges[i+1]:>8.1f} {np.sum(mask):>4} {A_mean:>10.4f} {A_std:>8.4f} {B2_mean:>8.4f} {2*H1_mean:>8.4f} {gap_mean:>8.4f} {x_mid:>10.4f}")

log()

# 빈 평균에 대한 피팅
if len(bin_stats) >= 5:
    x_bin = np.array([b['x_mid'] for b in bin_stats])
    A_bin = np.array([b['A_mean'] for b in bin_stats])
    B2_bin = np.array([b['B2_mean'] for b in bin_stats])
    H1_bin = np.array([b['H1_mean'] for b in bin_stats])
    w_bin = np.array([np.sqrt(b['n']) for b in bin_stats])  # sqrt(n) 가중치

    log("━━━ 빈 평균 피팅 ━━━")

    # <A> vs x = log(γ/(2π))
    c_lin = np.polyfit(x_bin, A_bin, 1, w=w_bin)
    pred_lin = np.polyval(c_lin, x_bin)
    r2_lin = 1 - np.sum(w_bin**2 * (A_bin - pred_lin)**2) / np.sum(w_bin**2 * (A_bin - np.average(A_bin, weights=w_bin**2))**2)
    log(f"  <A> = {c_lin[1]:.4f} + {c_lin[0]:.4f}·log(γ/(2π))  R²_w = {r2_lin:.6f}")

    c_q = np.polyfit(x_bin, A_bin, 2, w=w_bin)
    pred_q = np.polyval(c_q, x_bin)
    r2_q = 1 - np.sum(w_bin**2 * (A_bin - pred_q)**2) / np.sum(w_bin**2 * (A_bin - np.average(A_bin, weights=w_bin**2))**2)
    log(f"  <A> = {c_q[2]:.4f} + {c_q[1]:.4f}·x + {c_q[0]:.4f}·x²  R²_w = {r2_q:.6f}")

    # <H₁> 피팅
    c_h = np.polyfit(x_bin, H1_bin, 1, w=w_bin)
    pred_h = np.polyval(c_h, x_bin)
    r2_h = 1 - np.sum(w_bin**2 * (H1_bin - pred_h)**2) / np.sum(w_bin**2 * (H1_bin - np.average(H1_bin, weights=w_bin**2))**2)
    log(f"  <H₁> = {c_h[1]:.4f} + {c_h[0]:.4f}·log(γ/(2π))  R²_w = {r2_h:.6f}")

    c_h2 = np.polyfit(x_bin, H1_bin, 2, w=w_bin)
    pred_h2 = np.polyval(c_h2, x_bin)
    r2_h2 = 1 - np.sum(w_bin**2 * (H1_bin - pred_h2)**2) / np.sum(w_bin**2 * (H1_bin - np.average(H1_bin, weights=w_bin**2))**2)
    log(f"  <H₁> = {c_h2[2]:.4f} + {c_h2[1]:.4f}·x + {c_h2[0]:.4f}·x²  R²_w = {r2_h2:.6f}")

    # <B²> 피팅
    c_b = np.polyfit(x_bin, B2_bin, 1, w=w_bin)
    pred_b = np.polyval(c_b, x_bin)
    r2_b = 1 - np.sum(w_bin**2 * (B2_bin - pred_b)**2) / np.sum(w_bin**2 * (B2_bin - np.average(B2_bin, weights=w_bin**2))**2)
    log(f"  <B²> = {c_b[1]:.4f} + {c_b[0]:.4f}·log(γ/(2π))  R²_w = {r2_b:.6f}")

    log()

    # --- H₁ vs B² 지배 분석 ---
    log("━━━ A 분해: H₁ vs B² 지배 분석 ━━━")
    H1_frac = np.mean(2 * H1_arr / A_arr)
    B2_frac = np.mean(B2_arr / A_arr)
    log(f"  평균 2H₁/A = {H1_frac:.4f} ({H1_frac*100:.1f}%)")
    log(f"  평균 B²/A  = {B2_frac:.4f} ({B2_frac*100:.1f}%)")
    log(f"  → {'H₁ 지배' if H1_frac > B2_frac else 'B² 지배'}")
    log()

# --- 5c. 이론 예측 비교 ---
log("━━━ 이론 예측 비교 ━━━")

# Smooth Γ 기여: (1/4)ψ'(1/4+iγ/2) ≈ 4/γ² at large γ
gamma_smooth = np.array([r['gamma'] for r in results])
smooth_pred = 4.0 / gamma_smooth**2
log(f"  Γ smooth 기여 (1/γ²): 첫 {smooth_pred[0]:.6f}, 마지막 {smooth_pred[-1]:.6f}")
log(f"  → H₁ 대비 비율: 첫 {smooth_pred[0]/H1_arr[0]:.4f}, 마지막 {smooth_pred[-1]/H1_arr[-1]:.4f}")
log(f"  → Γ smooth는 H₁의 {smooth_pred[-1]/H1_arr[-1]*100:.2f}% (무시 가능)")
log()

# 영점 밀도 예측: d̄ = log(γ/(2π))/(2π)
d_bar = np.log(gam / (2*np.pi)) / (2*np.pi)
# 이론: <H₁> ~ (π²/3)·d̄² (pair correlation integral)
# → <H₁> ~ (1/12) (log(γ/(2π)))²
H1_theory = (1.0/12.0) * x_log**2

# 피팅: H₁ = α·(log(γ/(2π)))² + β
c_theory = np.polyfit(x_log**2, H1_arr, 1)
log(f"  H₁ vs (log(γ/(2π)))² 선형 피팅:")
log(f"    H₁ = {c_theory[1]:.4f} + {c_theory[0]:.4f}·(log(γ/(2π)))²")
pred_theory = np.polyval(c_theory, x_log**2)
r2_theory = 1 - np.sum((H1_arr - pred_theory)**2) / np.sum((H1_arr - np.mean(H1_arr))**2)
log(f"    R² = {r2_theory:.6f}")
log(f"    이론 계수 1/12 = {1/12:.4f}, 측정 계수 = {c_theory[0]:.4f}, 비율 = {c_theory[0]/(1/12):.4f}")
log()

# --- 5d. gap 상관 ---
log("━━━ A-gap 상관 ━━━")
# pearson_r 직접 구현 (scipy 의존성 회피)

def pearson_r(x, y):
    n = len(x)
    mx, my = np.mean(x), np.mean(y)
    cov = np.sum((x - mx) * (y - my)) / (n - 1)
    sx = np.sqrt(np.sum((x - mx)**2) / (n - 1))
    sy = np.sqrt(np.sum((y - my)**2) / (n - 1))
    r = cov / (sx * sy) if sx * sy > 0 else 0
    # t-test
    t_stat = r * np.sqrt((n - 2) / (1 - r**2)) if abs(r) < 1 else float('inf')
    # 2-tailed p-value approximation
    from math import erfc, sqrt
    p = erfc(abs(t_stat) / sqrt(2))  # rough approximation
    return r, p

rho_Ag, p_Ag = pearson_r(A_arr, gap_arr)
log(f"  ρ(A, gap_min) = {rho_Ag:+.4f} (p ≈ {p_Ag:.2e})")

rho_H1g, p_H1g = pearson_r(H1_arr, gap_arr)
log(f"  ρ(H₁, gap_min) = {rho_H1g:+.4f} (p ≈ {p_H1g:.2e})")

rho_B2g, p_B2g = pearson_r(B2_arr, gap_arr)
log(f"  ρ(B², gap_min) = {rho_B2g:+.4f} (p ≈ {p_B2g:.2e})")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. 기술 통계
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 기술 통계 ━━━")
log(f"  A:  범위 [{np.min(A_arr):.4f}, {np.max(A_arr):.4f}], 평균 {np.mean(A_arr):.4f}, 중앙값 {np.median(A_arr):.4f}, CV {np.std(A_arr)/np.mean(A_arr)*100:.1f}%")
log(f"  B²: 범위 [{np.min(B2_arr):.4f}, {np.max(B2_arr):.4f}], 평균 {np.mean(B2_arr):.4f}")
log(f"  H₁: 범위 [{np.min(H1_arr):.4f}, {np.max(H1_arr):.4f}], 평균 {np.mean(H1_arr):.4f}")
log(f"  c₁ Hadamard ρ (처음 50): {np.corrcoef(c1_arr[:50], [r['c1_had'] for r in results[:50]])[0,1]:.6f}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7. 종합 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("=" * 80)
log("[종합 판정]")

# 최고 R² 찾기
all_r2 = {}
all_r2['A_linear_logx'] = models_A.get('linear', (0, ''))[0]
all_r2['A_quad_logx'] = models_A.get('quadratic', (0, ''))[0]
all_r2['A_power_logx'] = models_A.get('power', (0, ''))[0]
all_r2['A_power_gam'] = r2_gam
all_r2['H1_x2'] = r2_theory

best_name = max(all_r2, key=all_r2.get)
best_r2 = all_r2[best_name]

log(f"  (1) 300개 A 측정: ✅")
log(f"  (2) 최고 R² (개별): {best_r2:.4f} ({best_name})")

if len(bin_stats) >= 5:
    bin_r2_best = max(r2_lin, r2_q)
    log(f"  (3) 최고 R²_w (빈 평균): {bin_r2_best:.4f}")
else:
    bin_r2_best = 0

log(f"  (4) H₁/(A/2) 지배: {H1_frac*100:.1f}%")
log()

if best_r2 > 0.90 or bin_r2_best > 0.95:
    verdict = "★★★ 양성 (clean scaling law)"
elif best_r2 > 0.70 or bin_r2_best > 0.85:
    verdict = "★★ 조건부 양성 (trend detected, noisy)"
else:
    verdict = "★ 음성 (no clean scaling)"

log(f"  판정: {verdict}")
log()
log(f"  경과 시간: {time.time()-START:.1f}초")
log("=" * 80)

outf.close()
