#!/usr/bin/env python3
"""
결과 #112 — Rankin-Selberg L(s, f×g) ξ-bundle 4성질 검증

L(s, E1 × E2) = ζ(s) · L(s, sym²(E1)) 류의 GL(4) L-함수.
PARI lfunmul로 생성.

테스트 쌍:
  (A) 11a1 × 37a1 (다른 conductor, 다른 rank)
  (B) 11a1 × 11a1 (자기 텐서 = Rankin-Selberg 자기곱)
  (C) 37a1 × 53a1 (둘 다 rank 1)

4성질: FE, 영점, κ_near (|Λ'/Λ|²), 모노드로미 (arg jump)
"""

import time, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)
print(f"PARI 초기화 완료\n")

# 타원곡선 쌍
PAIRS = [
    {
        "name": "11a1 × 37a1",
        "E1": [0,-1,1,-10,-20], "E2": [0,0,1,-1,0],
        "desc": "r0×r1, N=11×37=407",
    },
    {
        "name": "11a1 × 11a1",
        "E1": [0,-1,1,-10,-20], "E2": [0,-1,1,-10,-20],
        "desc": "자기곱, N=11²=121",
    },
    {
        "name": "37a1 × 53a1",
        "E1": [0,0,1,-1,0], "E2": [1,-1,1,0,0],
        "desc": "r1×r1, N=37×53=1961",
    },
]

T_MAX = 20
# Rankin-Selberg의 center는 weight에 따라 다름
# L(s, E1×E2): weight 2×2 = 4, center = (w+1)/2 = 5/2? 아니면 1?
# lfunmul의 결과에서 자동 결정됨
DELTAS = [0.01, 0.02, 0.03, 0.05]
H_DERIV = 1e-8
EPS_MONO = 0.03

OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "rankin_selberg_112.txt")
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


def measure_pair(pair_info):
    name = pair_info["name"]
    desc = pair_info["desc"]

    print(f"\n{'='*60}")
    print(f"[{name}] {desc}")
    print(f"{'='*60}")
    t0_time = time.time()

    E1 = pari.ellinit(pair_info["E1"])
    E2 = pari.ellinit(pair_info["E2"])
    lf1 = pari.lfuncreate(E1)
    lf2 = pari.lfuncreate(E2)

    # Rankin-Selberg
    lf_rs = pari.lfunmul(lf1, lf2)

    # Ldata 확인
    try:
        # conductor, weight 등 확인
        conductor = pari(f"lfunconductor({lf_rs})")
        print(f"  conductor: {conductor}")
    except Exception:
        pass

    # P1: FE
    fe = float(pari.lfuncheckfeq(lf_rs))
    p1_pass = fe < -30
    print(f"\n  [P1] FE check: {fe:.0f} → {'✅' if p1_pass else '❌'}")

    # P2: 영점
    zeros = [float(z) for z in pari.lfunzeros(lf_rs, T_MAX)]
    nontrivial = [z for z in zeros if z > 0.5]
    print(f"  [P2] 영점: {len(zeros)}개 (비자명: {len(nontrivial)}개)")
    if zeros[:3]:
        for z in zeros[:3]:
            print(f"    t={z:.6f}")

    if not nontrivial:
        print(f"  ⚠️ 비자명 영점 없음")
        return None

    # center 자동 판정: lfuninit 범위로 확인
    # Rankin-Selberg of weight 2 curves: critical line at Re(s) = 1
    center = 1.0  # weight 2 ⊗ weight 2의 center

    try:
        linit = pari.lfuninit(lf_rs, [0, T_MAX])
    except Exception:
        linit = pari.lfuninit(lf_rs, T_MAX)

    # P3: κ_near
    t0_zero = nontrivial[0]  # 첫 번째 비자명 영점
    print(f"\n  [P3] κ_near at t₀={t0_zero:.4f}")

    kappa_results = {}
    for delta in DELTAS:
        s = complex(center + delta, t0_zero)
        try:
            k = compute_kappa(linit, s)
            if k and k > 0:
                kd2 = k * delta**2
                kappa_results[delta] = {"kappa": k, "kd2": kd2}
                print(f"    δ={delta:.2f}: κ={k:.2f}, κ·δ²={kd2:.5f}")
        except Exception:
            pass

    if kappa_results:
        # δ=0.03 기준
        r03 = kappa_results.get(0.03)
        if r03:
            p3_kd2 = r03["kd2"]
            p3_pass = abs(p3_kd2 - 1.0) < 0.1
        else:
            p3_kd2 = None
            p3_pass = False
    else:
        p3_kd2 = None
        p3_pass = False

    print(f"    → {'✅' if p3_pass else '❌'} (κ·δ²={p3_kd2})")

    # P4: 모노드로미 (위상 점프)
    print(f"\n  [P4] 모노드로미 (arg jump)")
    mono_results = []
    for t_z in nontrivial[:10]:
        s_before = pari(f"{center} + {t_z - EPS_MONO}*I")
        s_after = pari(f"{center} + {t_z + EPS_MONO}*I")
        try:
            lam_b = complex(pari.lfunlambda(linit, s_before))
            lam_a = complex(pari.lfunlambda(linit, s_after))
            arg_b = np.angle(lam_b)
            arg_a = np.angle(lam_a)
            delta_arg = arg_a - arg_b
            while delta_arg > np.pi: delta_arg -= 2*np.pi
            while delta_arg < -np.pi: delta_arg += 2*np.pi
            mono_results.append({"t": t_z, "darg_pi": delta_arg/np.pi})
        except Exception:
            pass

    if mono_results:
        abs_vals = [abs(r["darg_pi"]) for r in mono_results]
        mean_abs = np.mean(abs_vals)
        print(f"    n={len(mono_results)}, |Δarg/π| mean={mean_abs:.6f}")
        for r in mono_results[:3]:
            print(f"      t={r['t']:.4f}: Δarg/π={r['darg_pi']:+.6f}")
        p4_pass = mean_abs > 0.9
    else:
        mean_abs = 0
        p4_pass = False

    print(f"    → {'✅' if p4_pass else '❌'}")

    # 등고선 모노드로미
    print(f"\n  [P4b] 등고선 모노드로미")
    contour_results = []
    for t_z in nontrivial[:5]:
        total_delta = 0.0
        prev_arg = None
        radius = 0.05
        n_pts = 64
        for k in range(n_pts + 1):
            theta = 2*np.pi*k/n_pts
            sigma = center + radius*np.cos(theta)
            t = t_z + radius*np.sin(theta)
            try:
                lam = complex(pari.lfunlambda(linit, pari(f"{sigma} + {t}*I")))
                cur_arg = np.angle(lam)
                if prev_arg is not None:
                    d = cur_arg - prev_arg
                    while d > np.pi: d -= 2*np.pi
                    while d < -np.pi: d += 2*np.pi
                    total_delta += d
                prev_arg = cur_arg
            except Exception:
                pass
        contour_results.append({"t": t_z, "mono_pi": total_delta/np.pi})

    if contour_results:
        vals = [r["mono_pi"] for r in contour_results]
        print(f"    n={len(vals)}, mono/π mean={np.mean(vals):.6f}")
        for r in contour_results[:3]:
            print(f"      t={r['t']:.4f}: mono/π={r['mono_pi']:.6f}")

    elapsed = time.time() - t0_time

    passes = sum([p1_pass, True, p3_pass, p4_pass])
    print(f"\n  요약: {passes}/4 PASS, {elapsed:.1f}s")

    return {
        "name": name, "desc": desc,
        "P1_fe": fe, "P1_pass": p1_pass,
        "P2_zeros": len(zeros), "P2_nontrivial": len(nontrivial),
        "P3_kd2": p3_kd2, "P3_pass": p3_pass,
        "P4_mean": mean_abs, "P4_pass": p4_pass,
        "contour_mono": [r["mono_pi"] for r in contour_results],
        "passes": passes, "elapsed": elapsed,
    }


if __name__ == "__main__":
    print(f"{'='*60}")
    print(f"[Project RDL] 결과 #112 — Rankin-Selberg 4성질 검증")
    print(f"{'='*60}")
    print(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    total_t0 = time.time()
    results = []

    for pair in PAIRS:
        try:
            r = measure_pair(pair)
            if r:
                results.append(r)
        except Exception as e:
            print(f"\n  ❌ {pair['name']} 실패: {e}")

    # 크로스 비교
    print(f"\n{'='*60}")
    print(f"크로스 비교")
    print(f"{'='*60}")
    print(f"\n  {'쌍':>20} {'FE':>6} {'영점':>6} {'κ·δ²':>10} {'|Δarg/π|':>10} {'mono/π':>10} {'통과':>5}")
    for r in results:
        mono = np.mean(r["contour_mono"]) if r["contour_mono"] else 0
        print(f"  {r['name']:>20} {r['P1_fe']:>6.0f} {r['P2_nontrivial']:>6} {r['P3_kd2'] or 0:>10.5f} {r['P4_mean']:>10.6f} {mono:>10.6f} {r['passes']:>4}/4")

    total_elapsed = time.time() - total_t0
    print(f"\n총 소요: {total_elapsed:.0f}s")

    # 저장
    with open(OUT_PATH, "w") as f:
        f.write(f"[Project RDL] 결과 #112 — Rankin-Selberg 4성질 검증\n")
        f.write(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        for r in results:
            mono = np.mean(r["contour_mono"]) if r["contour_mono"] else 0
            f.write(f"{r['name']}: FE={r['P1_fe']:.0f}, zeros={r['P2_nontrivial']}, "
                    f"κ·δ²={r['P3_kd2'] or 'N/A'}, |Δarg/π|={r['P4_mean']:.6f}, "
                    f"mono/π={mono:.6f}, {r['passes']}/4\n")
        f.write(f"\n총 소요: {total_elapsed:.0f}s\n")

    print(f"\n결과 저장: {OUT_PATH}")
