#!/usr/bin/env python3
"""
결과 #114 — Epstein ζ (non-Euler) 4성질 검증

핵심 질문: ξ-bundle 4성질이 오일러 곱 없는 L-함수에서도 성립하는가?

Epstein ζ: Z_Q(s) = Σ' Q(m,n)^{-s}  (Q = 양의 정부호 이차형식)
- class number h=1이면 Z_Q = c·ζ_K(s) (Euler product 있음)
- h>1이면 개별 형식의 Z_Q는 Euler product 없음 (non-multiplicative)
- 모든 경우 함수방정식 존재

대조군: [1,0;0,1] (disc=-4, h=1, Euler)
실험군: [1,0;0,5] (disc=-20, h=2, non-Euler)
        [2,1;1,3] (disc=-23, h=3, non-Euler)
        [1,0;0,6] (disc=-24, h=2, non-Euler)
"""

import time, sys, os
import numpy as np

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)
print("PARI 초기화 완료\n")

FORMS = [
    {"M": "[1,0;0,1]", "name": "x²+y²",        "disc": -4,  "h": 1, "euler": True,
     "note": "= 4·ζ(s)·L(s,χ₋₄)"},
    {"M": "[1,0;0,5]", "name": "x²+5y²",       "disc": -20, "h": 2, "euler": False,
     "note": "class group Z/2Z, non-Euler"},
    {"M": "[2,1;1,3]", "name": "2x²+xy+3y²",   "disc": -23, "h": 3, "euler": False,
     "note": "class group Z/3Z, non-Euler"},
    {"M": "[1,0;0,6]", "name": "x²+6y²",       "disc": -24, "h": 2, "euler": False,
     "note": "class group Z/2Z, non-Euler"},
]

T_MAX = 30
CENTER = 0.5  # Epstein ζ center = 1/2
DELTAS = [0.01, 0.02, 0.03, 0.05]
H_DERIV = 1e-8
EPS_MONO = 0.03
CONTOUR_RADIUS = 0.05
CONTOUR_STEPS = 64

OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "epstein_non_euler_114.txt")
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)


def lfunlambda_safe(linit, s_re, s_im):
    """안전한 lfunlambda 호출. 복소수 반환."""
    s_str = f"{s_re} + {s_im}*I" if s_im >= 0 else f"{s_re} - {abs(s_im)}*I"
    val = pari.lfunlambda(linit, pari(s_str))
    return complex(val)


def compute_kappa(linit, s_re, s_im, h=H_DERIV):
    """κ = |Λ'/Λ|² 수치 미분."""
    lam_s = lfunlambda_safe(linit, s_re, s_im)
    lam_p = lfunlambda_safe(linit, s_re + h, s_im)
    lam_m = lfunlambda_safe(linit, s_re - h, s_im)
    lam_d = (lam_p - lam_m) / (2 * h)
    if abs(lam_s) < 1e-50:
        return None
    return abs(lam_d / lam_s) ** 2


def measure_form(form_info):
    """하나의 이차형식에 대해 4성질 측정."""
    M_str = form_info["M"]
    name = form_info["name"]
    disc = form_info["disc"]
    h = form_info["h"]
    is_euler = form_info["euler"]

    print(f"\n{'='*60}")
    print(f"[{name}] disc={disc}, h={h}, Euler={is_euler}")
    print(f"  {form_info['note']}")
    print(f"{'='*60}")
    t0_time = time.time()

    # L-함수 생성
    lf_str = f"lfunqf({M_str})"
    lf = pari(lf_str)

    # P1: 함수방정식 확인
    fe = float(pari(f"lfuncheckfeq({lf_str})"))
    p1_pass = fe < -30
    print(f"  [P1] FE: {fe:.0f} → {'✅' if p1_pass else '❌'}")

    # P2: 영점
    zeros = [float(z) for z in pari(f"lfunzeros({lf_str}, {T_MAX})")]
    nontrivial = [z for z in zeros if z > 0.5]
    print(f"  [P2] 영점: {len(nontrivial)}개")

    if not nontrivial:
        elapsed = time.time() - t0_time
        return {"name": name, "disc": disc, "h": h, "euler": is_euler,
                "fe": fe, "zeros": 0, "kd2": None, "mono_linear": 0,
                "mono_contour": 0, "passes": 1, "elapsed": elapsed}

    # lfuninit (wide range to avoid insufficient init warning)
    try:
        linit = pari(f"lfuninit({lf_str}, [0, {T_MAX + 5}])")
    except Exception:
        linit = pari(f"lfuninit({lf_str}, {T_MAX + 5})")

    # P3: κ·δ² ≈ 1
    kd2_list = []
    kd2_report = None
    for t_z in nontrivial[:10]:
        for delta in DELTAS:
            try:
                k = compute_kappa(linit, CENTER + delta, t_z)
                if k and k > 0:
                    kd2 = k * delta**2
                    if delta == 0.03:
                        kd2_list.append(kd2)
            except Exception:
                pass

    if kd2_list:
        kd2_mean = np.mean(kd2_list)
        kd2_report = kd2_mean
        p3_pass = abs(kd2_mean - 1.0) < 0.1
        print(f"  [P3] κ·δ²(δ=0.03) = {kd2_mean:.6f} (n={len(kd2_list)})")
    else:
        p3_pass = False
        print(f"  [P3] 측정 실패")

    print(f"    → {'✅' if p3_pass else '❌'}")

    # P4a: 모노드로미 (선형 — 임계선 위 위상 점프)
    mono_vals = []
    for t_z in nontrivial[:15]:
        try:
            lam_b = lfunlambda_safe(linit, CENTER, t_z - EPS_MONO)
            lam_a = lfunlambda_safe(linit, CENTER, t_z + EPS_MONO)
            da = np.angle(lam_a) - np.angle(lam_b)
            while da > np.pi: da -= 2*np.pi
            while da < -np.pi: da += 2*np.pi
            mono_vals.append(abs(da / np.pi))
        except Exception:
            pass

    if mono_vals:
        mean_mono_linear = np.mean(mono_vals)
        p4a_pass = mean_mono_linear > 0.9
        print(f"  [P4a] |Δarg/π| mean={mean_mono_linear:.6f}, n={len(mono_vals)}")
    else:
        mean_mono_linear = 0
        p4a_pass = False

    # P4b: 모노드로미 (폐곡선 적분)
    contour_vals = []
    for t_z in nontrivial[:10]:
        total = 0.0
        prev_angle = None
        for k in range(CONTOUR_STEPS + 1):
            theta = 2 * np.pi * k / CONTOUR_STEPS
            sig = CENTER + CONTOUR_RADIUS * np.cos(theta)
            t = t_z + CONTOUR_RADIUS * np.sin(theta)
            try:
                lam = lfunlambda_safe(linit, sig, t)
                cur_angle = np.angle(lam)
                if prev_angle is not None:
                    d = cur_angle - prev_angle
                    while d > np.pi: d -= 2*np.pi
                    while d < -np.pi: d += 2*np.pi
                    total += d
                prev_angle = cur_angle
            except Exception:
                pass
        contour_vals.append(total / np.pi)

    if contour_vals:
        mean_mono_contour = np.mean(contour_vals)
        p4b_pass = abs(mean_mono_contour - 2.0) < 0.2
        print(f"  [P4b] mono/π mean={mean_mono_contour:.6f}, n={len(contour_vals)}")
    else:
        mean_mono_contour = 0
        p4b_pass = False

    p4_pass = p4a_pass and p4b_pass
    print(f"  [P4] → {'✅' if p4_pass else '❌'}")

    passes = sum([p1_pass, True, p3_pass, p4_pass])
    elapsed = time.time() - t0_time
    print(f"  통과: {passes}/4, {elapsed:.1f}s")

    return {
        "name": name, "disc": disc, "h": h, "euler": is_euler,
        "fe": fe, "zeros": len(nontrivial),
        "kd2": kd2_report,
        "mono_linear": mean_mono_linear, "mono_contour": mean_mono_contour,
        "passes": passes, "elapsed": elapsed,
        "zero_list": nontrivial[:5],  # 처음 5개 영점 기록
    }


def verify_non_euler(form_info):
    """비-Euler 검증: 개별 Epstein ζ 영점이 Dedekind ζ 영점과 다른지 확인."""
    M_str = form_info["M"]
    disc = form_info["disc"]

    epstein_zeros = [float(z) for z in pari(f"lfunzeros(lfunqf({M_str}), 20)")]

    # Dedekind ζ: disc에 대응하는 수체
    # disc = -4n (fundamental discriminant 아닐 수 있음)
    # 간단한 비교용으로만 사용
    abs_d = abs(disc)
    poly_map = {
        4: "x^2+1", 20: "x^2+5", 23: "x^2+x+6", 24: "x^2+6",
    }
    poly = poly_map.get(abs_d)
    if not poly:
        return None

    try:
        K = pari.nfinit(pari(poly))
        Lk = pari.lfuncreate(K)
        ded_zeros = [float(z) for z in pari.lfunzeros(Lk, 20)]
    except Exception:
        return None

    # 영점 비교
    matched = 0
    for ez in epstein_zeros[:5]:
        for dz in ded_zeros:
            if abs(ez - dz) < 0.01:
                matched += 1
                break

    return {
        "epstein_count": len(epstein_zeros),
        "dedekind_count": len(ded_zeros),
        "matched": matched,
        "first_epstein": epstein_zeros[:3] if epstein_zeros else [],
        "first_dedekind": ded_zeros[:3] if ded_zeros else [],
    }


if __name__ == "__main__":
    print(f"{'='*60}")
    print(f"[Project RDL] 결과 #114 — Epstein ζ (non-Euler) 4성질 검증")
    print(f"{'='*60}")
    print(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"T_MAX={T_MAX}, δ={DELTAS}, contour_r={CONTOUR_RADIUS}")

    total_t0 = time.time()
    results = []

    for form in FORMS:
        try:
            r = measure_form(form)
            results.append(r)
        except Exception as e:
            print(f"\n  ❌ {form['name']} 실패: {e}")
            import traceback; traceback.print_exc()

    # 비-Euler 검증
    print(f"\n{'='*60}")
    print(f"비-Euler 검증: Epstein vs Dedekind 영점 비교")
    print(f"{'='*60}")
    for form in FORMS:
        if not form["euler"]:
            comp = verify_non_euler(form)
            if comp:
                print(f"\n  [{form['name']}] disc={form['disc']}")
                print(f"    Epstein 영점 {comp['epstein_count']}개 vs Dedekind 영점 {comp['dedekind_count']}개")
                print(f"    처음 3개: Ep={[f'{z:.3f}' for z in comp['first_epstein']]}")
                print(f"              Ded={[f'{z:.3f}' for z in comp['first_dedekind']]}")
                print(f"    일치: {comp['matched']}/5 → {'❌ Euler' if comp['matched'] >= 4 else '✅ non-Euler 확인'}")

    # 크로스 비교
    print(f"\n{'='*60}")
    print(f"크로스 비교")
    print(f"{'='*60}")
    print(f"\n  {'형식':>15} {'disc':>5} {'h':>3} {'Euler':>6} {'FE':>6} {'영점':>5} {'κ·δ²':>10} {'|Δarg/π|':>10} {'mono/π':>10} {'통과':>5}")
    for r in results:
        euler_str = "Y" if r["euler"] else "N"
        kd2_str = f"{r['kd2']:.6f}" if r['kd2'] else "N/A"
        print(f"  {r['name']:>15} {r['disc']:>5} {r['h']:>3} {euler_str:>6} {r['fe']:>6.0f} "
              f"{r['zeros']:>5} {kd2_str:>10} {r['mono_linear']:>10.6f} "
              f"{r['mono_contour']:>10.6f} {r['passes']:>4}/4")

    total_elapsed = time.time() - total_t0
    print(f"\n총 소요: {total_elapsed:.0f}s")

    # 판정
    euler_results = [r for r in results if r["euler"]]
    non_euler_results = [r for r in results if not r["euler"]]

    euler_pass_all = all(r["passes"] == 4 for r in euler_results)
    non_euler_pass_all = all(r["passes"] == 4 for r in non_euler_results)
    non_euler_pass_any = any(r["passes"] == 4 for r in non_euler_results)

    print(f"\n{'='*60}")
    print(f"판정")
    print(f"{'='*60}")
    print(f"  대조군 (Euler): {'4/4 ALL PASS' if euler_pass_all else 'FAIL 존재'}")
    print(f"  실험군 (non-Euler): {len([r for r in non_euler_results if r['passes']==4])}/{len(non_euler_results)} PASS")

    if non_euler_pass_all:
        print(f"\n  ★★★ 강양성: 오일러 곱 없는 L-함수에서도 ξ-bundle 4성질 성립!")
        print(f"      → 프레임워크는 함수방정식에만 의존, 오일러 곱은 불필요")
    elif non_euler_pass_any:
        print(f"\n  ★★ 조건부 양성: 일부 non-Euler 형식에서 4성질 성립")
    else:
        print(f"\n  ★ 음성: non-Euler 형식에서 4성질 실패")
        print(f"      → 오일러 곱이 프레임워크의 필요 조건일 가능성")

    # 결과 저장
    with open(OUT_PATH, "w") as f:
        f.write(f"[Project RDL] 결과 #114 — Epstein ζ (non-Euler) 4성질 검증\n")
        f.write(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"T_MAX={T_MAX}, δ={DELTAS}, contour_r={CONTOUR_RADIUS}\n\n")

        for r in results:
            euler_str = "Euler" if r["euler"] else "non-Euler"
            kd2_str = f"{r['kd2']:.6f}" if r['kd2'] else "N/A"
            f.write(f"{r['name']} disc={r['disc']} h={r['h']} ({euler_str}): "
                    f"FE={r['fe']:.0f}, zeros={r['zeros']}, "
                    f"κ·δ²={kd2_str}, |Δarg/π|={r['mono_linear']:.6f}, "
                    f"mono/π={r['mono_contour']:.6f}, {r['passes']}/4\n")

        f.write(f"\n판정: {'★★★ 강양성' if non_euler_pass_all else '★★ 조건부' if non_euler_pass_any else '★ 음성'}\n")
        f.write(f"총 소요: {total_elapsed:.0f}s\n")

    print(f"\n결과 저장: {OUT_PATH}")
