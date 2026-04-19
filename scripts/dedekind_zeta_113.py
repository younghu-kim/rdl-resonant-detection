#!/usr/bin/env python3
"""
결과 #113 — Dedekind ζ 함수 4성질 검증

ζ_K(s) for imaginary quadratic fields K = Q(√D).
ζ_K(s) = ζ(s) · L(s, χ_D) — 오일러 곱 있음 (degree 2).

테스트: D = -3, -4, -7, -8, -11, -23, -163 (class number 1,1,1,1,1,3,1)

PARI nfinit + lfuncreate.
"""

import time, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)
print(f"PARI 초기화 완료\n")

FIELDS = [
    {"D": -3,   "poly": "x^2 + x + 1",   "name": "Q(√-3)",   "h": 1},
    {"D": -4,   "poly": "x^2 + 1",        "name": "Q(i)",      "h": 1},
    {"D": -7,   "poly": "x^2 + x + 2",   "name": "Q(√-7)",   "h": 1},
    {"D": -8,   "poly": "x^2 + 2",        "name": "Q(√-2)",   "h": 1},
    {"D": -11,  "poly": "x^2 + x + 3",   "name": "Q(√-11)",  "h": 1},
    {"D": -23,  "poly": "x^2 + x + 6",   "name": "Q(√-23)",  "h": 3},
    {"D": -163, "poly": "x^2 + x + 41",  "name": "Q(√-163)", "h": 1},
]

T_MAX = 30
CENTER = 0.5  # Dedekind ζ center = 1/2
DELTAS = [0.01, 0.02, 0.03, 0.05]
H_DERIV = 1e-8
EPS_MONO = 0.03

OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "dedekind_zeta_113.txt")
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)


def compute_kappa(linit, s, h=H_DERIV):
    lam_s = complex(pari.lfunlambda(linit, pari(f"{s.real} + {s.imag}*I")))
    lam_p = complex(pari.lfunlambda(linit, pari(f"{s.real + h} + {s.imag}*I")))
    lam_m = complex(pari.lfunlambda(linit, pari(f"{s.real - h} + {s.imag}*I")))
    lam_d = (lam_p - lam_m) / (2 * h)
    if abs(lam_s) < 1e-50:
        return None
    return abs(lam_d / lam_s) ** 2


def measure_field(field_info):
    name = field_info["name"]
    D = field_info["D"]
    h = field_info["h"]

    print(f"\n{'='*60}")
    print(f"[{name}] D={D}, h={h}")
    print(f"{'='*60}")
    t0_time = time.time()

    K = pari.nfinit(pari(field_info["poly"]))
    lf = pari.lfuncreate(K)

    # P1: FE
    fe = float(pari.lfuncheckfeq(lf))
    p1_pass = fe < -30
    print(f"  [P1] FE: {fe:.0f} → {'✅' if p1_pass else '❌'}")

    # P2: 영점
    zeros = [float(z) for z in pari.lfunzeros(lf, T_MAX)]
    nontrivial = [z for z in zeros if z > 0.5]
    print(f"  [P2] 영점: {len(zeros)}개 (비자명: {len(nontrivial)}개)")

    if not nontrivial:
        print(f"  ⚠️ 비자명 영점 없음")
        elapsed = time.time() - t0_time
        return {"name": name, "D": D, "h": h, "fe": fe, "zeros": len(zeros), "passes": 1, "elapsed": elapsed}

    try:
        linit = pari.lfuninit(lf, [0, T_MAX])
    except Exception:
        linit = pari.lfuninit(lf, T_MAX)

    # P3: κ_near
    t0_zero = nontrivial[0]
    kd2_03 = None
    for delta in DELTAS:
        s = complex(CENTER + delta, t0_zero)
        try:
            k = compute_kappa(linit, s)
            if k and k > 0:
                kd2 = k * delta**2
                if delta == 0.03:
                    kd2_03 = kd2
        except Exception:
            pass

    p3_pass = kd2_03 is not None and abs(kd2_03 - 1.0) < 0.1
    print(f"  [P3] κ·δ²(0.03) = {kd2_03:.5f}" if kd2_03 else "  [P3] 측정 실패")
    print(f"    → {'✅' if p3_pass else '❌'}")

    # P4: 모노드로미
    mono_vals = []
    for t_z in nontrivial[:10]:
        try:
            lam_b = complex(pari.lfunlambda(linit, pari(f"{CENTER} + {t_z - EPS_MONO}*I")))
            lam_a = complex(pari.lfunlambda(linit, pari(f"{CENTER} + {t_z + EPS_MONO}*I")))
            da = np.angle(lam_a) - np.angle(lam_b)
            while da > np.pi: da -= 2*np.pi
            while da < -np.pi: da += 2*np.pi
            mono_vals.append(abs(da/np.pi))
        except Exception:
            pass

    if mono_vals:
        mean_mono = np.mean(mono_vals)
        p4_pass = mean_mono > 0.9
        print(f"  [P4] |Δarg/π| mean={mean_mono:.6f}, n={len(mono_vals)}")
    else:
        mean_mono = 0
        p4_pass = False

    # 등고선
    contour_vals = []
    for t_z in nontrivial[:5]:
        total = 0.0
        prev = None
        for k in range(65):
            theta = 2*np.pi*k/64
            sig = CENTER + 0.05*np.cos(theta)
            t = t_z + 0.05*np.sin(theta)
            try:
                lam = complex(pari.lfunlambda(linit, pari(f"{sig} + {t}*I")))
                cur = np.angle(lam)
                if prev is not None:
                    d = cur - prev
                    while d > np.pi: d -= 2*np.pi
                    while d < -np.pi: d += 2*np.pi
                    total += d
                prev = cur
            except Exception:
                pass
        contour_vals.append(total/np.pi)

    if contour_vals:
        mean_cont = np.mean(contour_vals)
        print(f"  [P4b] mono/π mean={mean_cont:.6f}")
    else:
        mean_cont = 0

    passes = sum([p1_pass, True, p3_pass, p4_pass])
    elapsed = time.time() - t0_time
    print(f"  통과: {passes}/4, {elapsed:.1f}s")

    return {
        "name": name, "D": D, "h": h,
        "fe": fe, "zeros": len(nontrivial),
        "kd2": kd2_03, "mono_linear": mean_mono, "mono_contour": mean_cont,
        "passes": passes, "elapsed": elapsed,
    }


if __name__ == "__main__":
    print(f"{'='*60}")
    print(f"[Project RDL] 결과 #113 — Dedekind ζ 4성질 검증")
    print(f"{'='*60}")
    print(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    total_t0 = time.time()
    results = []

    for field in FIELDS:
        try:
            r = measure_field(field)
            results.append(r)
        except Exception as e:
            print(f"\n  ❌ {field['name']} 실패: {e}")

    # 크로스 비교
    print(f"\n{'='*60}")
    print(f"크로스 비교")
    print(f"{'='*60}")
    print(f"\n  {'수체':>12} {'D':>5} {'h':>3} {'FE':>6} {'영점':>5} {'κ·δ²':>10} {'|Δarg/π|':>10} {'mono/π':>10} {'통과':>5}")
    for r in results:
        print(f"  {r['name']:>12} {r['D']:>5} {r['h']:>3} {r['fe']:>6.0f} {r['zeros']:>5} "
              f"{r.get('kd2', 0) or 0:>10.5f} {r.get('mono_linear', 0):>10.6f} "
              f"{r.get('mono_contour', 0):>10.6f} {r['passes']:>4}/4")

    total_elapsed = time.time() - total_t0
    print(f"\n총 소요: {total_elapsed:.0f}s")

    with open(OUT_PATH, "w") as f:
        f.write(f"[Project RDL] 결과 #113 — Dedekind ζ 4성질 검증\n")
        f.write(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        for r in results:
            f.write(f"{r['name']} D={r['D']} h={r['h']}: FE={r['fe']:.0f}, zeros={r['zeros']}, "
                    f"κ·δ²={r.get('kd2', 'N/A')}, |Δarg/π|={r.get('mono_linear', 0):.6f}, "
                    f"mono/π={r.get('mono_contour', 0):.6f}, {r['passes']}/4\n")
        f.write(f"\n총 소요: {total_elapsed:.0f}s\n")

    print(f"\n결과 저장: {OUT_PATH}")
