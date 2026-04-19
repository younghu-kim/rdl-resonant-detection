#!/usr/bin/env python3
"""
결과 #109 — 타원곡선 κ_near (논문 정의) 검증

논문 정의: κ(s) = |Λ'(s)/Λ(s)|²
영점 ρ = center + it₀ 근처에서 κ(ρ + δ) ≈ 1/δ² (단순영점).

테스트: 8개 타원곡선 (rank 0–3) × 5 δ값 × 다수 영점
질문: κ·δ² ≈ 1이 rank/conductor 무관하게 성립하는가?

PARI/GP: lfunlambda + 수치미분으로 Λ'/Λ 계산.
"""

import time, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)
print(f"PARI 초기화: 2GB 메모리, realprecision=100\n")

CURVES = [
    {"label": "11a1",   "coeffs": [0,-1,1,-10,-20],  "rank": 0, "N": 11},
    {"label": "43a1",   "coeffs": [0,1,1,0,0],        "rank": 0, "N": 43},
    {"label": "197a1",  "coeffs": [0,1,1,-2,-4],       "rank": 0, "N": 197},
    {"label": "37a1",   "coeffs": [0,0,1,-1,0],        "rank": 1, "N": 37},
    {"label": "53a1",   "coeffs": [1,-1,1,0,0],         "rank": 1, "N": 53},
    {"label": "79a1",   "coeffs": [1,1,1,-2,0],         "rank": 1, "N": 79},
    {"label": "389a1",  "coeffs": [0,1,1,-2,0],         "rank": 2, "N": 389},
    {"label": "5077a1", "coeffs": [0,0,1,-7,6],         "rank": 3, "N": 5077},
]

T_MAX = 30
CENTER = 1.0  # weight 2 critical point
DELTAS = [0.01, 0.02, 0.03, 0.05, 0.10]
H_DERIV = 1e-8  # 수치미분 간격

OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "ec_kappa_near_109.txt")
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

def compute_kappa(linit, s, h=H_DERIV):
    """
    κ(s) = |Λ'(s)/Λ(s)|²
    Λ' 는 수치미분: (Λ(s+h) - Λ(s-h)) / (2h)
    """
    lam_s = complex(pari.lfunlambda(linit, pari(f"{s.real} + {s.imag}*I")))
    lam_plus = complex(pari.lfunlambda(linit, pari(f"{s.real + h} + {s.imag}*I")))
    lam_minus = complex(pari.lfunlambda(linit, pari(f"{s.real - h} + {s.imag}*I")))

    lam_deriv = (lam_plus - lam_minus) / (2 * h)

    if abs(lam_s) < 1e-50:
        return None  # 영점 위에서는 발산

    omega = lam_deriv / lam_s
    kappa = abs(omega) ** 2
    return kappa


def measure_curve(info):
    label = info["label"]
    coeffs = info["coeffs"]
    rank = info["rank"]
    N = info["N"]

    print(f"\n{'='*60}")
    print(f"[{label}] rank={rank}, N={N}")
    print(f"{'='*60}")
    t0_time = time.time()

    E = pari.ellinit(coeffs)
    lf = pari.lfuncreate(E)
    eps = int(pari.ellrootno(E))
    print(f"  ε={eps}")

    # 영점 탐색
    zeros_raw = pari.lfunzeros(lf, T_MAX)
    zeros = [float(z) for z in zeros_raw]
    # 비자명 영점만 (t > 1.0으로 center zero 제외)
    nontrivial = [z for z in zeros if z > 1.0]
    print(f"  영점: {len(zeros)}개 (비자명 t>1: {len(nontrivial)}개)")

    if not nontrivial:
        print(f"  ⚠️ 비자명 영점 없음, 스킵")
        return None

    # lfuninit
    try:
        linit = pari.lfuninit(lf, [0, T_MAX])
    except Exception:
        linit = pari.lfuninit(lf, T_MAX)

    # 각 영점에서 κ 측정
    results_by_delta = {d: [] for d in DELTAS}

    for t_zero in nontrivial[:10]:  # 최대 10개 영점
        for delta in DELTAS:
            s = complex(CENTER + delta, t_zero)
            try:
                kappa = compute_kappa(linit, s)
                if kappa is not None and kappa > 0:
                    kd2 = kappa * delta**2
                    results_by_delta[delta].append({
                        "t": t_zero,
                        "kappa": kappa,
                        "kd2": kd2,
                    })
            except Exception:
                pass

    # 요약 출력
    print(f"\n  δ별 κ·δ² 통계 (비자명 영점):")
    print(f"  {'δ':>8} {'1/δ²':>10} {'κ_near(med)':>14} {'κ·δ²(med)':>12} {'CV':>8} {'n':>4}")
    print(f"  {'-'*60}")

    summary = {}
    for delta in DELTAS:
        data = results_by_delta[delta]
        if data:
            kappas = [d["kappa"] for d in data]
            kd2s = [d["kd2"] for d in data]
            med_k = np.median(kappas)
            med_kd2 = np.median(kd2s)
            cv = np.std(kappas) / np.mean(kappas) * 100 if np.mean(kappas) > 0 else 0
            print(f"  {delta:>8.2f} {1/delta**2:>10.0f} {med_k:>14.2f} {med_kd2:>12.5f} {cv:>7.2f}% {len(data):>4}")
            summary[delta] = {"median_kappa": med_k, "median_kd2": med_kd2, "cv": cv, "n": len(data)}
        else:
            print(f"  {delta:>8.2f} {'--':>10} {'--':>14} {'--':>12} {'--':>8} {0:>4}")

    elapsed = time.time() - t0_time
    print(f"\n  소요: {elapsed:.1f}s")

    return {
        "label": label, "rank": rank, "N": N, "epsilon": eps,
        "n_zeros_total": len(zeros), "n_nontrivial": len(nontrivial),
        "summary": summary, "elapsed": elapsed,
    }


if __name__ == "__main__":
    print(f"{'='*60}")
    print(f"[Project RDL] 결과 #109 — 타원곡선 κ_near (논문 정의)")
    print(f"{'='*60}")
    print(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"κ(s) = |Λ'(s)/Λ(s)|², 영점 근처 δ-오프셋 측정")
    print(f"곡선: {len(CURVES)}개, δ: {DELTAS}")

    total_t0 = time.time()
    all_results = []

    for info in CURVES:
        try:
            r = measure_curve(info)
            if r:
                all_results.append(r)
        except Exception as e:
            print(f"\n  ❌ {info['label']} 실패: {e}")
            import traceback
            traceback.print_exc()

    # ── 크로스 비교 ──
    print(f"\n{'='*60}")
    print(f"크로스 비교: δ=0.03 기준 κ·δ² 비교")
    print(f"{'='*60}")

    print(f"\n  {'곡선':>10} {'rank':>5} {'N':>6} {'ε':>3} {'κ_near':>12} {'κ·δ²':>10} {'CV':>8}")
    print(f"  {'-'*60}")
    for r in all_results:
        s = r["summary"].get(0.03)
        if s:
            print(f"  {r['label']:>10} {r['rank']:>5} {r['N']:>6} {r['epsilon']:>3} {s['median_kappa']:>12.2f} {s['median_kd2']:>10.5f} {s['cv']:>7.2f}%")

    # κ·δ² 일관성 판정
    kd2_vals = []
    for r in all_results:
        s = r["summary"].get(0.03)
        if s:
            kd2_vals.append(s["median_kd2"])

    if kd2_vals:
        mean_kd2 = np.mean(kd2_vals)
        std_kd2 = np.std(kd2_vals)
        cv_kd2 = std_kd2 / mean_kd2 * 100 if mean_kd2 > 0 else 0
        print(f"\n  전곡선 κ·δ² 통계: mean={mean_kd2:.5f}, std={std_kd2:.5f}, CV={cv_kd2:.2f}%")
        print(f"  판정: {'✅ 보편 성립 (CV < 5%)' if cv_kd2 < 5 else '❌ 비보편 (CV ≥ 5%)'}")

    total_elapsed = time.time() - total_t0
    print(f"\n총 소요: {total_elapsed:.0f}s ({total_elapsed/60:.1f}분)")

    # 파일 저장
    with open(OUT_PATH, "w") as f:
        f.write(f"{'='*60}\n")
        f.write(f"[Project RDL] 결과 #109 — 타원곡선 κ_near (논문 정의)\n")
        f.write(f"{'='*60}\n")
        f.write(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"κ(s) = |Λ'(s)/Λ(s)|², δ-오프셋 측정\n\n")

        for r in all_results:
            f.write(f"\n[{r['label']}] rank={r['rank']}, N={r['N']}, ε={r['epsilon']}\n")
            f.write(f"  영점: {r['n_zeros_total']}개 (비자명: {r['n_nontrivial']}개)\n")
            for delta in DELTAS:
                s = r["summary"].get(delta)
                if s:
                    f.write(f"  δ={delta:.2f}: κ_near={s['median_kappa']:.2f}, κ·δ²={s['median_kd2']:.5f}, CV={s['cv']:.2f}%, n={s['n']}\n")

        f.write(f"\n{'='*60}\n")
        f.write(f"크로스 비교 (δ=0.03)\n")
        f.write(f"{'='*60}\n")
        for r in all_results:
            s = r["summary"].get(0.03)
            if s:
                f.write(f"  {r['label']} rank={r['rank']} N={r['N']}: κ·δ²={s['median_kd2']:.5f}\n")

        if kd2_vals:
            f.write(f"\n전곡선: mean={mean_kd2:.5f}, std={std_kd2:.5f}, CV={cv_kd2:.2f}%\n")

        f.write(f"\n총 소요: {total_elapsed:.0f}s\n")

    print(f"\n결과 저장: {OUT_PATH}")
