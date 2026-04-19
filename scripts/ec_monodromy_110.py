#!/usr/bin/env python3
"""
결과 #110 — 타원곡선 모노드로미 정밀 측정

논문 기준: 각 영점에서 Δarg(Λ) = ±π (simple zero → phase jump π).
GL(1)–GL(5): mono/π = 2.0000 (argument principle, 등고선 적분 기준).
여기서는 임계선을 따라 영점 전후 위상 차이 측정.

AP-17 준수: log-space 계산으로 overflow 방지.
AP-18 준수: ε=-1에서도 arg(Λ) 측정은 유효 (arg는 Re/Im 무관).

PARI/GP via cypari2.
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
CENTER = 1.0
EPS_OFFSET = 0.05  # 영점 전후 측정 거리

OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "ec_monodromy_110.txt")
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)


def measure_monodromy(linit, t_zero, center, eps=EPS_OFFSET):
    """
    영점 t_zero 전후 Λ 위상 변화 측정.
    Returns: Δarg/π (expected: ±1 for simple zero)
    """
    s_before = pari(f"{center} + {t_zero - eps}*I")
    s_after = pari(f"{center} + {t_zero + eps}*I")

    lam_before = complex(pari.lfunlambda(linit, s_before))
    lam_after = complex(pari.lfunlambda(linit, s_after))

    arg_before = np.angle(lam_before)
    arg_after = np.angle(lam_after)

    # branch cut 보정 (AP-07)
    delta_arg = arg_after - arg_before
    while delta_arg > np.pi:
        delta_arg -= 2 * np.pi
    while delta_arg < -np.pi:
        delta_arg += 2 * np.pi

    return delta_arg / np.pi


def measure_contour_monodromy(linit, t_zero, center, radius=0.05, n_points=64):
    """
    영점 주위 작은 원형 등고선에서 arg(Λ) 적분.
    Returns: winding number (expected: 1 for simple zero → Δarg = 2π → mono/π = 2)
    """
    total_delta = 0.0
    prev_arg = None

    for k in range(n_points + 1):
        theta = 2 * np.pi * k / n_points
        sigma = center + radius * np.cos(theta)
        t = t_zero + radius * np.sin(theta)
        s = pari(f"{sigma} + {t}*I")

        try:
            lam = complex(pari.lfunlambda(linit, s))
            current_arg = np.angle(lam)
        except Exception:
            continue

        if prev_arg is not None:
            d = current_arg - prev_arg
            while d > np.pi:
                d -= 2 * np.pi
            while d < -np.pi:
                d += 2 * np.pi
            total_delta += d

        prev_arg = current_arg

    return total_delta / np.pi


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

    # 영점
    zeros_raw = pari.lfunzeros(lf, T_MAX)
    zeros = [float(z) for z in zeros_raw]
    nontrivial = [z for z in zeros if z > 1.0]
    print(f"  영점: {len(zeros)}개 (비자명 t>1: {len(nontrivial)}개)")

    if not nontrivial:
        print(f"  ⚠️ 비자명 영점 없음")
        return None

    try:
        linit = pari.lfuninit(lf, [0, T_MAX])
    except Exception:
        linit = pari.lfuninit(lf, T_MAX)

    # 각 영점에서 측정
    linear_results = []
    contour_results = []

    for t_zero in nontrivial[:15]:  # 최대 15개
        # 1. 임계선 위상 점프
        try:
            mono_lin = measure_monodromy(linit, t_zero, CENTER)
            linear_results.append({"t": t_zero, "mono_pi": mono_lin})
        except Exception:
            pass

        # 2. 등고선 적분 (mono/π = 2 기대)
        try:
            mono_cont = measure_contour_monodromy(linit, t_zero, CENTER)
            contour_results.append({"t": t_zero, "mono_pi": mono_cont})
        except Exception:
            pass

    # 결과 출력
    print(f"\n  [A] 임계선 위상 점프 (Δarg/π):")
    if linear_results:
        vals = [r["mono_pi"] for r in linear_results]
        print(f"    n={len(vals)}, mean={np.mean(vals):.6f}, std={np.std(vals):.6f}")
        print(f"    |Δarg/π| mean={np.mean(np.abs(vals)):.6f}")
        for r in linear_results[:5]:
            print(f"      t={r['t']:.4f}: Δarg/π = {r['mono_pi']:+.6f}")
        if len(linear_results) > 5:
            print(f"      ... ({len(linear_results) - 5}개 더)")

    print(f"\n  [B] 등고선 적분 (winding number × 2):")
    if contour_results:
        vals = [r["mono_pi"] for r in contour_results]
        print(f"    n={len(vals)}, mean={np.mean(vals):.6f}, std={np.std(vals):.6f}")
        for r in contour_results[:5]:
            print(f"      t={r['t']:.4f}: mono/π = {r['mono_pi']:.6f}")
        if len(contour_results) > 5:
            print(f"      ... ({len(contour_results) - 5}개 더)")

    elapsed = time.time() - t0_time
    print(f"\n  소요: {elapsed:.1f}s")

    return {
        "label": label, "rank": rank, "N": N, "epsilon": eps,
        "linear": linear_results, "contour": contour_results,
        "elapsed": elapsed,
    }


if __name__ == "__main__":
    print(f"{'='*60}")
    print(f"[Project RDL] 결과 #110 — 타원곡선 모노드로미 정밀 측정")
    print(f"{'='*60}")
    print(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"임계선 위상 점프 + 등고선 적분 (argument principle)")

    total_t0 = time.time()
    all_results = []

    for info in CURVES:
        try:
            r = measure_curve(info)
            if r:
                all_results.append(r)
        except Exception as e:
            print(f"\n  ❌ {info['label']} 실패: {e}")

    # 크로스 비교
    print(f"\n{'='*60}")
    print(f"크로스 비교")
    print(f"{'='*60}")

    print(f"\n  {'곡선':>10} {'rank':>5} {'ε':>3} {'|Δarg/π| (linear)':>20} {'mono/π (contour)':>20}")
    print(f"  {'-'*65}")
    for r in all_results:
        lin_vals = [abs(x["mono_pi"]) for x in r["linear"]]
        cont_vals = [x["mono_pi"] for x in r["contour"]]
        lin_mean = np.mean(lin_vals) if lin_vals else 0
        cont_mean = np.mean(cont_vals) if cont_vals else 0
        print(f"  {r['label']:>10} {r['rank']:>5} {r['epsilon']:>3} {lin_mean:>20.6f} {cont_mean:>20.6f}")

    total_elapsed = time.time() - total_t0
    print(f"\n총 소요: {total_elapsed:.0f}s ({total_elapsed/60:.1f}분)")

    # 저장
    with open(OUT_PATH, "w") as f:
        f.write(f"{'='*60}\n")
        f.write(f"[Project RDL] 결과 #110 — 타원곡선 모노드로미 정밀 측정\n")
        f.write(f"{'='*60}\n")
        f.write(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        for r in all_results:
            f.write(f"\n[{r['label']}] rank={r['rank']}, N={r['N']}, ε={r['epsilon']}\n")
            lin_vals = [abs(x["mono_pi"]) for x in r["linear"]]
            cont_vals = [x["mono_pi"] for x in r["contour"]]
            if lin_vals:
                f.write(f"  Linear |Δarg/π|: mean={np.mean(lin_vals):.6f}, std={np.std(lin_vals):.6f}, n={len(lin_vals)}\n")
            if cont_vals:
                f.write(f"  Contour mono/π: mean={np.mean(cont_vals):.6f}, std={np.std(cont_vals):.6f}, n={len(cont_vals)}\n")
            for x in r["linear"][:5]:
                f.write(f"    t={x['t']:.4f}: linear={x['mono_pi']:+.6f}\n")
            for x in r["contour"][:5]:
                f.write(f"    t={x['t']:.4f}: contour={x['mono_pi']:.6f}\n")

        f.write(f"\n총 소요: {total_elapsed:.0f}s\n")

    print(f"\n결과 저장: {OUT_PATH}")
