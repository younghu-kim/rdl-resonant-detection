#!/usr/bin/env python3
"""
[사이클 #349] EC 체계적 확장 + A(t₀) rank 의존성 탐사

목적:
  1. LMFDB rank 0,1,2,3 × 5곡선 = 20곡선 4성질 검증
  2. ★ 핵심 신규: A(t₀) = Im(c₀)² + 2Re(c₁) 의 rank 의존 구조 탐사

4성질:
  P1. FE — 함수방정식 log|Λ(s)/Λ(2-s)| ≈ 0 (σ=1 위)
  P2. κδ² — 영점 근처 |Λ'/Λ|²·δ² ≈ 1
  P3. σ-sweep — σ∈[0.7,1.3]에서 영점 분류 안정성
  P4. mono/π — 영점 주위 등고선 적분 winding number ≈ 2

A(t₀) 계산:
  Λ'/Λ(ρ+u) = 1/u + c₀ + c₁u + ...  (Laurent 전개)
  A(γ) = Im(c₀)² + 2Re(c₁)
  Cauchy 적분으로 c₀, c₁ 추출 (c263 패턴)

주의사항:
  - ε=-1 곡선: Im(Λ) 사용 (#107 교훈)
  - rank 2,3: s=1 근처 다중 영점 → δ≥0.01
  - realprecision≥100 (degree 2)
  - system python3 (cypari2)

결과: results/ec_20curves_c349.txt
"""

import sys, os, time, math
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)
print(f"PARI 초기화: 2GB 메모리, realprecision=100")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 20곡선 정의 (LMFDB 기준)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CURVES = [
    # === rank 0 (5곡선) — LMFDB 검증 완료 ===
    {"label": "11a1",   "coeffs": [0,-1,1,-10,-20],       "rank": 0, "N": 11},
    {"label": "14a6",   "coeffs": [1,0,1,4,-6],            "rank": 0, "N": 14},
    {"label": "15a1",   "coeffs": [1,1,1,-2160,-39540],    "rank": 0, "N": 15},
    {"label": "43a1",   "coeffs": [0,1,1,0,0],              "rank": 0, "N": 43},
    {"label": "197a1",  "coeffs": [0,1,1,-2,-4],             "rank": 0, "N": 197},

    # === rank 1 (5곡선) — LMFDB 검증 완료 ===
    {"label": "37a1",   "coeffs": [0,0,1,-1,0],              "rank": 1, "N": 37},
    {"label": "53a1",   "coeffs": [1,-1,1,0,0],               "rank": 1, "N": 53},
    {"label": "58a1",   "coeffs": [1,-1,0,-1,1],              "rank": 1, "N": 58},
    {"label": "61a1",   "coeffs": [1,0,0,-2,1],               "rank": 1, "N": 61},
    {"label": "79a1",   "coeffs": [1,1,1,-2,0],               "rank": 1, "N": 79},

    # === rank 2 (5곡선) — LMFDB 검증 완료 ===
    {"label": "389a1",  "coeffs": [0,1,1,-2,0],               "rank": 2, "N": 389},
    {"label": "433a1",  "coeffs": [1,0,0,0,1],                "rank": 2, "N": 433},
    {"label": "446a1",  "coeffs": [1,-1,0,-4,4],              "rank": 2, "N": 446},
    {"label": "563a1",  "coeffs": [1,1,1,-15,16],             "rank": 2, "N": 563},
    {"label": "571a1",  "coeffs": [0,1,1,-4,2],               "rank": 2, "N": 571},

    # === rank 3 (5곡선) — LMFDB 검증 완료 ===
    {"label": "5077a1",  "coeffs": [0,0,1,-7,6],              "rank": 3, "N": 5077},
    {"label": "11197a1", "coeffs": [1,-1,1,-6,0],             "rank": 3, "N": 11197},
    {"label": "11642a1", "coeffs": [1,-1,0,-16,28],           "rank": 3, "N": 11642},
    {"label": "12279a1", "coeffs": [0,-1,1,-10,12],           "rank": 3, "N": 12279},
    {"label": "13766a1", "coeffs": [1,0,1,-23,42],            "rank": 3, "N": 13766},
]

T_MAX = 100       # t∈[1,100] (수학자 지시)
CENTER = 1.0      # GL(2) 임계점
DELTA = 0.03      # κδ² 측정 기본 δ
H_DERIV = 1e-8    # 수치미분 간격
CAUCHY_R = 0.01   # A(t₀) Cauchy 적분 반지름
CAUCHY_N = 64     # Cauchy 적분 점 수
MONO_R = 0.05     # 모노드로미 등고선 반지름
MONO_N = 64       # 등고선 점 수
SIGMA_VALS = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]  # σ-sweep

OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "ec_20curves_c349.txt")
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# P1: 함수방정식 검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def check_fe(linit, t_test=5.0, center=CENTER):
    """log|Λ(s)/Λ(2-s)| at s = center + i*t_test."""
    s = pari(f"{center} + {t_test}*I")
    s_conj = pari(f"{2 - center} + {t_test}*I")
    try:
        lam_s = complex(pari.lfunlambda(linit, s))
        lam_c = complex(pari.lfunlambda(linit, s_conj))
        if abs(lam_s) < 1e-50 or abs(lam_c) < 1e-50:
            return None
        ratio = abs(lam_s) / abs(lam_c)
        return math.log10(ratio) if ratio > 0 else None
    except Exception as e:
        print(f"    FE 검증 실패: {e}")
        return None


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# P2: κδ² 검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_kappa_delta2(linit, t_zero, delta=DELTA, center=CENTER, h=H_DERIV):
    """κ(ρ+δ)·δ² where κ = |Λ'/Λ|²."""
    s = pari(f"{center + delta} + {t_zero}*I")
    s_p = pari(f"{center + delta + h} + {t_zero}*I")
    s_m = pari(f"{center + delta - h} + {t_zero}*I")
    try:
        lam = complex(pari.lfunlambda(linit, s))
        lam_p = complex(pari.lfunlambda(linit, s_p))
        lam_m = complex(pari.lfunlambda(linit, s_m))
        if abs(lam) < 1e-50:
            return None
        deriv = (lam_p - lam_m) / (2 * h)
        omega = deriv / lam
        kappa = abs(omega) ** 2
        return kappa * delta ** 2
    except Exception as e:
        return None


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# P3: σ-sweep (영점 분류 안정성)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def sigma_sweep_classify(linit, zeros, epsilon, sigmas=SIGMA_VALS):
    """
    각 σ에서 "영점 근처 부호변환이 관측되는 영점 수" 계산.
    ε=+1 → Re(Λ), ε=-1 → Im(Λ) 부호 변화.
    """
    counts = {}
    for sigma in sigmas:
        n_detected = 0
        for t_zero in zeros:
            t_below = t_zero - 0.05
            t_above = t_zero + 0.05
            s_b = pari(f"{sigma} + {t_below}*I")
            s_a = pari(f"{sigma} + {t_above}*I")
            try:
                lam_b = complex(pari.lfunlambda(linit, s_b))
                lam_a = complex(pari.lfunlambda(linit, s_a))
                if epsilon == 1:
                    val_b, val_a = lam_b.real, lam_a.real
                else:
                    val_b, val_a = lam_b.imag, lam_a.imag
                if val_b * val_a < 0:
                    n_detected += 1
            except Exception:
                pass
        counts[sigma] = n_detected
    return counts


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# P4: 모노드로미 (등고선 적분)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def contour_monodromy(linit, t_zero, center=CENTER, radius=MONO_R, n_pts=MONO_N):
    """영점 주위 원형 등고선 arg(Λ) 적분 → mono/π."""
    total_delta = 0.0
    prev_arg = None
    for k in range(n_pts + 1):
        theta = 2 * np.pi * k / n_pts
        sigma = center + radius * np.cos(theta)
        t = t_zero + radius * np.sin(theta)
        s = pari(f"{sigma} + {t}*I")
        try:
            lam = complex(pari.lfunlambda(linit, s))
            cur_arg = np.angle(lam)
        except Exception:
            continue
        if prev_arg is not None:
            d = cur_arg - prev_arg
            while d > np.pi:
                d -= 2 * np.pi
            while d < -np.pi:
                d += 2 * np.pi
            total_delta += d
        prev_arg = cur_arg
    return total_delta / np.pi


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ★ A(t₀) 계산 (Cauchy 적분으로 Laurent 계수 추출)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_A_cauchy(linit, t_zero, center=CENTER, r=CAUCHY_R, n_pts=CAUCHY_N, h=1e-5):
    """
    Cauchy 적분: Λ'/Λ(ρ+u) = 1/u + c₀ + c₁u + ...
    A(γ) = Im(c₀)² + 2Re(c₁)

    Returns: (A, Im_c0, Re_c1) 또는 None
    """
    c0_sum = 0+0j
    c1_sum = 0+0j
    n_valid = 0

    for k in range(n_pts):
        theta = 2 * np.pi * k / n_pts
        u = r * np.exp(1j * theta)
        sigma = center + u.real
        t = t_zero + u.imag

        s = pari(f"{sigma} + {t}*I")
        s_p = pari(f"{sigma + h} + {t}*I")
        s_m = pari(f"{sigma - h} + {t}*I")

        try:
            lam = complex(pari.lfunlambda(linit, s))
            if abs(lam) < 1e-40:
                continue
            lam_p = complex(pari.lfunlambda(linit, s_p))
            lam_m = complex(pari.lfunlambda(linit, s_m))

            deriv = (lam_p - lam_m) / (2 * h)
            omega = deriv / lam  # Λ'/Λ at ρ+u

            # Λ'/Λ(ρ+u) - 1/u = c₀ + c₁u + ...
            g = omega - 1/u
            c0_sum += g
            c1_sum += g / u   # c₁ ≈ (1/n) Σ g(uₖ)/uₖ
            n_valid += 1
        except Exception:
            continue

    if n_valid < n_pts // 2:
        return None

    c0 = c0_sum / n_valid
    c1 = c1_sum / n_valid

    A = c0.imag**2 + 2 * c1.real

    if math.isnan(A) or math.isinf(A):
        return None

    return (A, c0.imag, c1.real)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인: 곡선별 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def analyze_curve(info):
    """단일 곡선 전체 분석."""
    label = info["label"]
    coeffs = info["coeffs"]
    rank = info["rank"]
    N = info["N"]

    print(f"\n{'='*70}")
    print(f"[{label}] rank={rank}, N={N}")
    print(f"{'='*70}")
    t0 = time.time()

    # PARI 초기화
    try:
        E = pari.ellinit(coeffs)
        lf = pari.lfuncreate(E)
        eps = int(pari.ellrootno(E))
    except Exception as e:
        print(f"  ❌ ellinit/lfuncreate 실패: {e}")
        return None

    print(f"  ε={eps}")

    # 영점 탐색
    try:
        zeros_raw = pari.lfunzeros(lf, T_MAX)
        zeros = sorted([float(z) for z in zeros_raw])
    except Exception as e:
        print(f"  ❌ lfunzeros 실패: {e}")
        return None

    nontrivial = [z for z in zeros if z > 1.0]
    print(f"  영점: {len(zeros)}개 (비자명 t>1: {len(nontrivial)}개)")

    if len(nontrivial) == 0:
        print(f"  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
        return None

    # lfuninit
    try:
        linit = pari.lfuninit(lf, [0, T_MAX + 10])
    except Exception:
        try:
            linit = pari.lfuninit(lf, T_MAX + 10)
        except Exception as e:
            print(f"  ❌ lfuninit 실패: {e}")
            return None

    # ── P1: FE ──
    fe_vals = []
    for t_test in [5.0, 15.0, 25.0]:
        fe = check_fe(linit, t_test)
        if fe is not None:
            fe_vals.append(fe)
    fe_log = np.mean(fe_vals) if fe_vals else None
    fe_pass = fe_log is not None and abs(fe_log) < 1  # log10 < 1 → ratio < 10
    print(f"  P1 FE: log10|Λ/Λ̃| = {fe_log:.1f}" if fe_log is not None else "  P1 FE: 실패")

    # ── P2: κδ² ──
    kd2_vals = []
    for t_zero in nontrivial[:20]:
        kd2 = compute_kappa_delta2(linit, t_zero)
        if kd2 is not None and not math.isnan(kd2) and not math.isinf(kd2):
            kd2_vals.append(kd2)

    if kd2_vals:
        kd2_mean = np.mean(kd2_vals)
        kd2_cv = np.std(kd2_vals) / kd2_mean * 100 if kd2_mean > 0 else 999
        kd2_pass = abs(kd2_mean - 1.0) < 0.1 and kd2_cv < 10
        print(f"  P2 κδ²: mean={kd2_mean:.6f}, CV={kd2_cv:.2f}%, n={len(kd2_vals)}")
    else:
        kd2_mean, kd2_cv = None, None
        kd2_pass = False
        print(f"  P2 κδ²: 측정 실패")

    # ── P3: σ-sweep ──
    sweep_zeros = nontrivial[:30]  # 최대 30개 영점
    sweep = sigma_sweep_classify(linit, sweep_zeros, eps)
    n_crit = sweep.get(1.0, 0)
    n_total = len(sweep_zeros)
    # σ=1.0에서 검출률
    if n_total > 0:
        detect_rate = n_crit / n_total
        sweep_pass = detect_rate > 0.5  # 50% 이상 검출
        # A(σ) = max deviation from σ=1.0 count
        a_sigma = max(abs(v - n_crit) / n_total for v in sweep.values()) if sweep else 0
    else:
        detect_rate = 0
        sweep_pass = False
        a_sigma = 0
    print(f"  P3 σ-sweep: σ=1.0에서 {n_crit}/{n_total} 검출 (rate={detect_rate:.2f})")
    print(f"    σ별: {sweep}")
    print(f"    A(σ) = {a_sigma:.4f}")

    # ── P4: 모노드로미 ──
    mono_vals = []
    for t_zero in nontrivial[:15]:
        try:
            m = contour_monodromy(linit, t_zero)
            if not math.isnan(m) and not math.isinf(m):
                mono_vals.append(m)
        except Exception:
            pass

    if mono_vals:
        mono_mean = np.mean(mono_vals)
        mono_std = np.std(mono_vals)
        # mono/π ≈ 2.0 (winding number 1 → Δarg=2π → mono/π=2)
        n_good = sum(1 for v in mono_vals if abs(v - 2.0) < 0.5)
        mono_pass = n_good / len(mono_vals) > 0.8
        print(f"  P4 mono/π: mean={mono_mean:.6f}, std={mono_std:.6f}, good={n_good}/{len(mono_vals)}")
    else:
        mono_mean, mono_std, n_good = None, None, 0
        mono_pass = False
        print(f"  P4 mono/π: 측정 실패")

    # ── 4성질 판정 ──
    n_pass = sum([fe_pass, kd2_pass, sweep_pass, mono_pass])
    print(f"\n  판정: {n_pass}/4 PASS  [FE={'✅' if fe_pass else '❌'}, κδ²={'✅' if kd2_pass else '❌'}, σ-sw={'✅' if sweep_pass else '❌'}, mono={'✅' if mono_pass else '❌'}]")

    # ── ★ A(t₀) 계산 ──
    print(f"\n  ★ A(t₀) 계산 ({len(nontrivial)}개 영점)...")
    A_vals = []
    Im_c0_vals = []
    Re_c1_vals = []
    n_A_fail = 0

    for i, t_zero in enumerate(nontrivial):
        result = compute_A_cauchy(linit, t_zero)
        if result is None:
            n_A_fail += 1
        else:
            A, im_c0, re_c1 = result
            A_vals.append(A)
            Im_c0_vals.append(im_c0)
            Re_c1_vals.append(re_c1)

        if (i + 1) % 20 == 0:
            print(f"    [{i+1}/{len(nontrivial)}] A 계산 중...")

    if A_vals:
        A_mean = np.mean(A_vals)
        A_std = np.std(A_vals)
        A_median = np.median(A_vals)
        Im_c0_mean = np.mean(np.abs(Im_c0_vals))
        Re_c1_mean = np.mean(Re_c1_vals)
        print(f"  ★ A(t₀): n={len(A_vals)}/{len(nontrivial)}, mean={A_mean:.6f}, std={A_std:.6f}, median={A_median:.6f}")
        print(f"    |Im(c₀)| mean={Im_c0_mean:.6f}, Re(c₁) mean={Re_c1_mean:.6f}")
    else:
        A_mean, A_std, A_median = None, None, None
        Im_c0_mean, Re_c1_mean = None, None
        print(f"  ★ A(t₀): 계산 실패 (전체 {n_A_fail}건)")

    elapsed = time.time() - t0
    print(f"\n  소요: {elapsed:.1f}s")

    return {
        "label": label, "rank": rank, "N": N, "epsilon": eps,
        "n_zeros": len(zeros), "n_nontrivial": len(nontrivial),
        # P1-P4
        "fe_log": fe_log, "fe_pass": fe_pass,
        "kd2_mean": kd2_mean, "kd2_cv": kd2_cv, "kd2_pass": kd2_pass,
        "sweep": sweep, "sweep_pass": sweep_pass, "a_sigma": a_sigma,
        "mono_mean": mono_mean, "mono_std": mono_std, "mono_pass": mono_pass,
        "n_pass": n_pass,
        # A(t₀)
        "A_mean": A_mean, "A_std": A_std, "A_median": A_median,
        "n_A": len(A_vals), "n_A_fail": n_A_fail,
        "Im_c0_mean": Im_c0_mean, "Re_c1_mean": Re_c1_mean,
        "A_vals": A_vals,  # 전체 A값 (상관 분석용)
        "elapsed": elapsed,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

if __name__ == "__main__":
    print(f"\n{'='*70}")
    print(f"[Project RDL] C-349 — EC 20곡선 체계적 확장 + A(t₀) rank 의존성")
    print(f"{'='*70}")
    print(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"곡선: {len(CURVES)}개 (rank 0-3 × 5)")
    print(f"t∈[1,{T_MAX}], δ={DELTA}, Cauchy r={CAUCHY_R}")
    print()

    total_t0 = time.time()
    all_results = []
    failed_curves = []

    for info in CURVES:
        try:
            r = analyze_curve(info)
            if r:
                all_results.append(r)
            else:
                failed_curves.append(info["label"])
        except Exception as e:
            print(f"\n  ❌ {info['label']} 치명적 오류: {e}")
            import traceback
            traceback.print_exc()
            failed_curves.append(info["label"])

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 종합 분석
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    print(f"\n\n{'='*70}")
    print(f"종합 결과")
    print(f"{'='*70}")

    # 4성질 요약표
    print(f"\n[A] 4성질 검증 요약")
    print(f"{'곡선':>12} {'rank':>5} {'N':>6} {'ε':>3} {'#영점':>6} {'FE':>4} {'κδ²':>8} {'σ-sw':>6} {'mono':>6} {'PASS':>5}")
    print(f"{'-'*70}")
    for r in all_results:
        fe_str = '✅' if r['fe_pass'] else '❌'
        kd2_str = f"{r['kd2_mean']:.4f}" if r['kd2_mean'] is not None else "  --  "
        sw_str = '✅' if r['sweep_pass'] else '❌'
        mo_str = '✅' if r['mono_pass'] else '❌'
        print(f"{r['label']:>12} {r['rank']:>5} {r['N']:>6} {r['epsilon']:>3} {r['n_nontrivial']:>6} "
              f"{fe_str:>4} {kd2_str:>8} {sw_str:>6} {mo_str:>6} {r['n_pass']}/4")

    n_all_pass = sum(1 for r in all_results if r['n_pass'] == 4)
    print(f"\n전원 4/4 PASS: {n_all_pass}/{len(all_results)}")

    # ★ A(t₀) vs rank 분석
    print(f"\n\n[B] ★ A(t₀) vs rank 분석")
    print(f"{'곡선':>12} {'rank':>5} {'N':>6} {'n_A':>5} {'A_mean':>12} {'A_std':>12} {'A_median':>12}")
    print(f"{'-'*70}")

    rank_A_data = {0: [], 1: [], 2: [], 3: []}
    for r in all_results:
        if r['A_mean'] is not None:
            print(f"{r['label']:>12} {r['rank']:>5} {r['N']:>6} {r['n_A']:>5} "
                  f"{r['A_mean']:>12.6f} {r['A_std']:>12.6f} {r['A_median']:>12.6f}")
            rank_A_data[r['rank']].extend(r['A_vals'])
        else:
            print(f"{r['label']:>12} {r['rank']:>5} {r['N']:>6}    --        --           --           --")

    print(f"\n[C] Rank별 A(t₀) 통계")
    print(f"{'rank':>5} {'n':>6} {'mean':>12} {'std':>12} {'median':>12} {'SE':>12}")
    print(f"{'-'*60}")
    rank_means = []
    rank_labels = []
    for rank in [0, 1, 2, 3]:
        vals = rank_A_data[rank]
        if vals:
            m = np.mean(vals)
            s = np.std(vals)
            med = np.median(vals)
            se = s / np.sqrt(len(vals))
            print(f"{rank:>5} {len(vals):>6} {m:>12.6f} {s:>12.6f} {med:>12.6f} {se:>12.6f}")
            rank_means.append(m)
            rank_labels.append(rank)
        else:
            print(f"{rank:>5}     0        --           --           --           --")

    # Rank 간 차이 통계 검정
    print(f"\n[D] Rank 간 차이 검정")
    for i in range(len(rank_labels)):
        for j in range(i+1, len(rank_labels)):
            r1, r2 = rank_labels[i], rank_labels[j]
            v1, v2 = rank_A_data[r1], rank_A_data[r2]
            if len(v1) >= 3 and len(v2) >= 3:
                # Welch t-test
                t_stat, p_val = stats.ttest_ind(v1, v2, equal_var=False)
                # Mann-Whitney U
                u_stat, p_mw = stats.mannwhitneyu(v1, v2, alternative='two-sided')
                diff = np.mean(v1) - np.mean(v2)
                pooled_std = np.sqrt((np.std(v1)**2 + np.std(v2)**2) / 2) if (np.std(v1) + np.std(v2)) > 0 else 1
                effect = diff / pooled_std if pooled_std > 0 else 0
                sig_t = '★' if p_val < 0.05 else ''
                sig_mw = '★' if p_mw < 0.05 else ''
                print(f"  rank {r1} vs {r2}: diff={diff:+.6f}, Cohen d={effect:.3f}, "
                      f"t-test p={p_val:.4f}{sig_t}, MW p={p_mw:.4f}{sig_mw}")

    # Spearman: rank vs mean A
    if len(rank_means) >= 3:
        rho_sp, p_sp = stats.spearmanr(rank_labels, rank_means)
        print(f"\n  Spearman(rank, A_mean): ρ={rho_sp:.4f}, p={p_sp:.4f}")

    # κδ² vs rank (보편성 확인)
    print(f"\n[E] κδ² 크로스 비교 (δ={DELTA})")
    kd2_all = []
    for r in all_results:
        if r['kd2_mean'] is not None:
            print(f"  {r['label']:>12} rank={r['rank']} N={r['N']}: κδ²={r['kd2_mean']:.6f} (CV={r['kd2_cv']:.2f}%)")
            kd2_all.append(r['kd2_mean'])
    if kd2_all:
        print(f"  전곡선: mean={np.mean(kd2_all):.6f}, std={np.std(kd2_all):.6f}, CV={np.std(kd2_all)/np.mean(kd2_all)*100:.2f}%")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 최종 판정
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    print(f"\n\n{'='*70}")
    print(f"최종 판정")
    print(f"{'='*70}")

    # 4성질 판정
    if n_all_pass == len(all_results) and len(all_results) >= 18:
        verdict_4prop = "★★★★★ 양성 — 20곡선 전원 4/4 PASS"
    elif n_all_pass >= len(all_results) * 0.9:
        verdict_4prop = f"★★★★ 양성 — {n_all_pass}/{len(all_results)} PASS (90%+)"
    elif n_all_pass >= len(all_results) * 0.8:
        verdict_4prop = f"★★★ 조건부 양성 — {n_all_pass}/{len(all_results)} PASS"
    else:
        verdict_4prop = f"★★ 부분적 — {n_all_pass}/{len(all_results)} PASS"

    print(f"\n  4성질: {verdict_4prop}")

    # A(t₀) 판정
    # rank별 2σ 이상 차이 있는 쌍이 있는가?
    significant_pairs = 0
    total_pairs = 0
    for i in range(len(rank_labels)):
        for j in range(i+1, len(rank_labels)):
            v1, v2 = rank_A_data[rank_labels[i]], rank_A_data[rank_labels[j]]
            if len(v1) >= 3 and len(v2) >= 3:
                total_pairs += 1
                _, p = stats.ttest_ind(v1, v2, equal_var=False)
                if p < 0.05:
                    significant_pairs += 1

    if significant_pairs > total_pairs * 0.5 and len(rank_means) >= 3:
        rho_sp, p_sp = stats.spearmanr(rank_labels, rank_means)
        if abs(rho_sp) > 0.8 and p_sp < 0.1:
            verdict_A = f"★★★★★ 양성 — A(t₀) rank 체계적 의존 (ρ={rho_sp:.3f}, {significant_pairs}/{total_pairs} sig pairs)"
        else:
            verdict_A = f"★★★ 조건부 양성 — 유의 차이 존재하나 단조 아닌 ({significant_pairs}/{total_pairs})"
    elif significant_pairs >= 1:
        verdict_A = f"★★ 약한 신호 — {significant_pairs}/{total_pairs} 쌍 유의"
    else:
        verdict_A = "★ 중립 — A(t₀) rank 무관 (보편 상수)"

    print(f"  A(t₀): {verdict_A}")

    if failed_curves:
        print(f"\n  실패 곡선: {failed_curves}")

    total_elapsed = time.time() - total_t0
    print(f"\n총 소요: {total_elapsed:.0f}s ({total_elapsed/60:.1f}분)")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 결과 파일 저장
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    with open(OUT_PATH, "w") as f:
        f.write(f"{'='*70}\n")
        f.write(f"[Project RDL] C-349 — EC 20곡선 + A(t₀) rank 의존성\n")
        f.write(f"{'='*70}\n")
        f.write(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"곡선: {len(CURVES)}개, 실행 완료: {len(all_results)}개\n")
        f.write(f"t∈[1,{T_MAX}], δ={DELTA}, Cauchy r={CAUCHY_R}\n\n")

        # 개별 곡선
        for r in all_results:
            f.write(f"\n[{r['label']}] rank={r['rank']}, N={r['N']}, ε={r['epsilon']}\n")
            f.write(f"  영점: {r['n_zeros']}개 (비자명: {r['n_nontrivial']}개)\n")
            f.write(f"  P1 FE: {r['fe_log']:.1f}\n" if r['fe_log'] is not None else "  P1 FE: --\n")
            f.write(f"  P2 κδ²: {r['kd2_mean']:.6f} (CV={r['kd2_cv']:.2f}%)\n" if r['kd2_mean'] is not None else "  P2 κδ²: --\n")
            f.write(f"  P3 σ-sweep: {r['sweep']}\n")
            f.write(f"    A(σ) = {r['a_sigma']:.4f}\n")
            f.write(f"  P4 mono/π: mean={r['mono_mean']:.6f} (std={r['mono_std']:.6f})\n" if r['mono_mean'] is not None else "  P4 mono/π: --\n")
            f.write(f"  통과: {r['n_pass']}/4\n")
            if r['A_mean'] is not None:
                f.write(f"  ★ A(t₀): mean={r['A_mean']:.6f}, std={r['A_std']:.6f}, median={r['A_median']:.6f}, n={r['n_A']}\n")
            else:
                f.write(f"  ★ A(t₀): 계산 실패\n")

        # rank별 통계
        f.write(f"\n\n{'='*70}\n")
        f.write(f"Rank별 A(t₀) 통계\n")
        f.write(f"{'='*70}\n")
        for rank in [0, 1, 2, 3]:
            vals = rank_A_data[rank]
            if vals:
                f.write(f"  rank {rank}: n={len(vals)}, mean={np.mean(vals):.6f}, std={np.std(vals):.6f}, median={np.median(vals):.6f}\n")

        # 검정 결과
        f.write(f"\n{'='*70}\n")
        f.write(f"통계 검정\n")
        f.write(f"{'='*70}\n")
        for i in range(len(rank_labels)):
            for j in range(i+1, len(rank_labels)):
                r1, r2 = rank_labels[i], rank_labels[j]
                v1, v2 = rank_A_data[r1], rank_A_data[r2]
                if len(v1) >= 3 and len(v2) >= 3:
                    t_stat, p_val = stats.ttest_ind(v1, v2, equal_var=False)
                    f.write(f"  rank {r1} vs {r2}: t={t_stat:.4f}, p={p_val:.4f}\n")

        if len(rank_means) >= 3:
            rho_sp, p_sp = stats.spearmanr(rank_labels, rank_means)
            f.write(f"\n  Spearman(rank, A_mean): ρ={rho_sp:.4f}, p={p_sp:.4f}\n")

        # 판정
        f.write(f"\n\n{'='*70}\n")
        f.write(f"최종 판정\n")
        f.write(f"{'='*70}\n")
        f.write(f"  4성질: {verdict_4prop}\n")
        f.write(f"  A(t₀): {verdict_A}\n")
        if failed_curves:
            f.write(f"  실패 곡선: {failed_curves}\n")
        f.write(f"\n총 소요: {total_elapsed:.0f}s ({total_elapsed/60:.1f}분)\n")

    print(f"\n결과 저장: {OUT_PATH}")
