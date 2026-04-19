#!/usr/bin/env python3
"""
결과 #111 — 타원곡선 κ_near conductor-스케일링

질문: κ_near의 O(1) 보정항 Δκ = κ_near - 1/δ² 가 conductor N에 어떻게 의존하는가?
#55 (GL(3)): 6개 곡선에서 κ_near가 conductor-독립 (α=0.003, R²=0.0002).
#109: 타원곡선 8개에서 κ_near CV=0.47%.

이번: rank=0 타원곡선 10개 (N=11→10099) + rank=1 타원곡선 6개 → Δκ vs log(N) 스케일링.

PARI/GP via cypari2.
"""

import time, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)
print(f"PARI 초기화 완료\n")

# rank 0 곡선 (ε=+1), conductor 증가순
CURVES_R0 = [
    {"label": "11a1",    "coeffs": [0,-1,1,-10,-20],     "N": 11},
    {"label": "14a1",    "coeffs": [1,0,1,4,-6],          "N": 14},
    {"label": "37b1",    "coeffs": [0,1,1,-23,-50],        "N": 37},  # rank 0
    {"label": "46a1",    "coeffs": [1,-1,0,-2,-1],         "N": 46},
    {"label": "99d1",    "coeffs": [0,0,1,-1,-1],          "N": 99},
    {"label": "197a1",   "coeffs": [0,1,1,-2,-4],          "N": 197},
    {"label": "571a1",   "coeffs": [0,-1,1,-929,-10595],   "N": 571},
    {"label": "997c1",   "coeffs": [0,0,1,1,-1],           "N": 997},
    {"label": "2999b1",  "coeffs": [0,1,1,4,-4],           "N": 2999},
    {"label": "10099a1", "coeffs": [0,0,1,-4,4],           "N": 10099},
]

# rank 1 곡선 (ε=-1)
CURVES_R1 = [
    {"label": "37a1",    "coeffs": [0,0,1,-1,0],           "N": 37},
    {"label": "53a1",    "coeffs": [1,-1,1,0,0],            "N": 53},
    {"label": "79a1",    "coeffs": [1,1,1,-2,0],            "N": 79},
    {"label": "389a1_r1","coeffs": [1,0,0,-4,4],            "N": 389},  # rank 1 curve near N=389
    {"label": "997b1",   "coeffs": [0,-1,1,-2,2],           "N": 997},
    {"label": "5077a1",  "coeffs": [0,0,1,-7,6],            "N": 5077},
]

T_MAX = 30
CENTER = 1.0
DELTA = 0.03  # 고정 δ (논문 기준점)
H_DERIV = 1e-8

OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "ec_kappa_conductor_111.txt")
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)


def compute_kappa(linit, s, h=H_DERIV):
    """κ(s) = |Λ'(s)/Λ(s)|²"""
    lam_s = complex(pari.lfunlambda(linit, pari(f"{s.real} + {s.imag}*I")))
    lam_p = complex(pari.lfunlambda(linit, pari(f"{s.real + h} + {s.imag}*I")))
    lam_m = complex(pari.lfunlambda(linit, pari(f"{s.real - h} + {s.imag}*I")))
    lam_d = (lam_p - lam_m) / (2 * h)
    if abs(lam_s) < 1e-50:
        return None
    return abs(lam_d / lam_s) ** 2


def measure_kappa_near(label, coeffs, N, rank_label):
    """κ_near 측정 + Δκ 계산"""
    E = pari.ellinit(coeffs)
    lf = pari.lfuncreate(E)
    eps = int(pari.ellrootno(E))

    zeros = [float(z) for z in pari.lfunzeros(lf, T_MAX)]
    nontrivial = [z for z in zeros if z > 1.0]

    if not nontrivial:
        return None

    try:
        linit = pari.lfuninit(lf, [0, T_MAX])
    except Exception:
        linit = pari.lfuninit(lf, T_MAX)

    kappas = []
    for t_z in nontrivial[:15]:
        s = complex(CENTER + DELTA, t_z)
        try:
            k = compute_kappa(linit, s)
            if k and k > 0:
                kappas.append(k)
        except Exception:
            pass

    if not kappas:
        return None

    med = np.median(kappas)
    kd2 = med * DELTA**2
    delta_k = med - 1.0 / DELTA**2  # O(1) 보정항

    print(f"  {label:>12} N={N:<6} ε={eps:+d} n={len(kappas):>3}  κ_near={med:>10.2f}  κ·δ²={kd2:.5f}  Δκ={delta_k:>8.2f}")

    return {
        "label": label, "N": N, "epsilon": eps, "rank": rank_label,
        "n": len(kappas), "kappa_near": med, "kd2": kd2, "delta_k": delta_k,
    }


if __name__ == "__main__":
    print(f"{'='*70}")
    print(f"[Project RDL] 결과 #111 — 타원곡선 κ conductor-스케일링")
    print(f"{'='*70}")
    print(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"δ={DELTA}, κ_near = |Λ'/Λ|² at (center+δ+it_zero)\n")

    total_t0 = time.time()

    # rank 0
    print(f"\n[rank 0 곡선]")
    print(f"  {'label':>12} {'N':>7} {'ε':>3} {'n':>5}  {'κ_near':>12}  {'κ·δ²':>10}  {'Δκ':>10}")
    print(f"  {'-'*70}")
    r0_results = []
    for info in CURVES_R0:
        try:
            r = measure_kappa_near(info["label"], info["coeffs"], info["N"], "r0")
            if r:
                r0_results.append(r)
        except Exception as e:
            print(f"  {info['label']:>12} ❌ {e}")

    # rank 1
    print(f"\n[rank 1 곡선]")
    print(f"  {'label':>12} {'N':>7} {'ε':>3} {'n':>5}  {'κ_near':>12}  {'κ·δ²':>10}  {'Δκ':>10}")
    print(f"  {'-'*70}")
    r1_results = []
    for info in CURVES_R1:
        try:
            r = measure_kappa_near(info["label"], info["coeffs"], info["N"], "r1")
            if r:
                r1_results.append(r)
        except Exception as e:
            print(f"  {info['label']:>12} ❌ {e}")

    all_results = r0_results + r1_results

    # 스케일링 분석
    print(f"\n{'='*70}")
    print(f"스케일링 분석: Δκ vs log(N)")
    print(f"{'='*70}")

    for rank_label, results in [("rank 0", r0_results), ("rank 1", r1_results)]:
        if len(results) < 3:
            print(f"\n  [{rank_label}] 데이터 부족 ({len(results)}개)")
            continue

        log_N = np.array([np.log(r["N"]) for r in results])
        delta_k = np.array([r["delta_k"] for r in results])

        # 선형 회귀: Δκ = a * log(N) + b
        A = np.vstack([log_N, np.ones(len(log_N))]).T
        slope, intercept = np.linalg.lstsq(A, delta_k, rcond=None)[0]
        residuals = delta_k - (slope * log_N + intercept)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((delta_k - np.mean(delta_k))**2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        print(f"\n  [{rank_label}] n={len(results)}")
        print(f"    Δκ = {slope:.4f} · log(N) + {intercept:.4f}")
        print(f"    R² = {r_squared:.4f}")
        print(f"    Δκ range: [{min(delta_k):.2f}, {max(delta_k):.2f}]")
        print(f"    mean Δκ = {np.mean(delta_k):.2f}, std = {np.std(delta_k):.2f}, CV = {np.std(delta_k)/np.mean(delta_k)*100:.2f}%")

        if r_squared < 0.5:
            print(f"    → ✅ conductor-독립 (R² < 0.5)")
        else:
            print(f"    → ⚠️ conductor 의존 가능 (R² ≥ 0.5)")

    total_elapsed = time.time() - total_t0
    print(f"\n총 소요: {total_elapsed:.0f}s")

    # 저장
    with open(OUT_PATH, "w") as f:
        f.write(f"{'='*70}\n")
        f.write(f"[Project RDL] 결과 #111 — 타원곡선 κ conductor-스케일링\n")
        f.write(f"{'='*70}\n")
        f.write(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"δ={DELTA}\n\n")

        for r in all_results:
            f.write(f"{r['label']} N={r['N']} ε={r['epsilon']} rank={r['rank']}: κ_near={r['kappa_near']:.2f}, Δκ={r['delta_k']:.2f}, κ·δ²={r['kd2']:.5f}\n")

        f.write(f"\n총 소요: {total_elapsed:.0f}s\n")

    print(f"\n결과 저장: {OUT_PATH}")
