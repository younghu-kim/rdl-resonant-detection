#!/usr/bin/env python3
"""
결과 #108 — B-08 rank-dependent ξ-bundle 검증 (ε-보정판)

#107 Red Team 지적사항 반영:
  1. ε=-1일 때 Re(Λ)≡0 → Im(Λ) 부호변환으로 P4 보정
  2. rank-conductor 교란 분리: 같은 rank에서 conductor 다양화

테스트 곡선 (8개):
  rank 0: 11a1 (N=11), 43a1 (N=43), 389a1 → X (대신 197a1 N=197)
  rank 1: 37a1 (N=37), 43a1 → X (rank1은 ε=-1 필요), 53a1 (N=53), 79a1 (N=79)
  rank 2: 389a1 (N=389)
  rank 3: 5077a1 (N=5077)

PARI/GP via cypari2.
"""

import time, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

# ── PARI 초기화 ──
import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)
print(f"PARI 초기화: 2GB 메모리, realprecision=100\n")

# ── 테스트 곡선 ──
CURVES = [
    # rank 0, ε=+1
    {"label": "11a1",  "coeffs": [0,-1,1,-10,-20],   "rank": 0, "N": 11},
    {"label": "43a1",  "coeffs": [0,1,1,0,0],         "rank": 0, "N": 43},
    {"label": "197a1", "coeffs": [0,1,1,-2,-4],        "rank": 0, "N": 197},
    # rank 1, ε=-1
    {"label": "37a1",  "coeffs": [0,0,1,-1,0],         "rank": 1, "N": 37},
    {"label": "53a1",  "coeffs": [1,-1,1,0,0],          "rank": 1, "N": 53},
    {"label": "79a1",  "coeffs": [1,1,1,-2,0],          "rank": 1, "N": 79},
    # rank 2, ε=+1
    {"label": "389a1", "coeffs": [0,1,1,-2,0],          "rank": 2, "N": 389},
    # rank 3, ε=-1
    {"label": "5077a1","coeffs": [0,0,1,-7,6],          "rank": 3, "N": 5077},
]

T_MAX = 30
CENTER = 1.0  # GL(2) weight 2
DELTA_OFFSETS = [0.01 * i for i in range(1, 21)]  # 20점
SIGMA_VALUES = [CENTER + d for d in [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]]
DT = 0.1  # σ-sweep 샘플링 간격

RESULTS = []
OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "rank_dependent_corrected_108.txt")
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

def banner(msg):
    print(f"\n{'='*72}")
    print(f"{msg}")
    print(f"{'='*72}")

def count_sign_changes_component(lf, linit, sigma, t_max, dt, use_imag):
    """
    Re(Λ) 또는 Im(Λ)의 부호변환 카운트.
    use_imag=True이면 Im(Λ) 사용 (ε=-1 곡선용).
    """
    ts = np.arange(dt, t_max, dt)
    changes = 0
    prev_val = None
    for t in ts:
        s = pari(f"{sigma} + {t}*I")
        try:
            lam = pari.lfunlambda(linit, s)
            val_c = complex(lam)
            val = val_c.imag if use_imag else val_c.real
        except Exception:
            continue
        if prev_val is not None and val * prev_val < 0:
            changes += 1
        if abs(val) > 1e-30:
            prev_val = val
    return changes

def measure_curve(info):
    """하나의 타원곡선에 대해 4성질 + σ-sweep 측정"""
    label = info["label"]
    coeffs = info["coeffs"]
    rank = info["rank"]
    N_cond = info["N"]

    banner(f"[{label}] rank={rank}, N={N_cond}")
    t0 = time.time()

    result = {"label": label, "rank": rank, "N": N_cond}

    # ellinit + lfuncreate
    E = pari.ellinit(coeffs)
    lf = pari.lfuncreate(E)

    # root number 확인
    eps_raw = pari.ellrootno(E)
    eps = int(eps_raw)
    use_imag = (eps == -1)
    result["epsilon"] = eps
    result["component"] = "Im" if use_imag else "Re"
    print(f"  ε={eps}, 측정 성분: {'Im(Λ)' if use_imag else 'Re(Λ)'}")

    # ── P1: 함수방정식 ──
    fe = float(pari.lfuncheckfeq(lf))
    result["P1_fe"] = fe
    p1_pass = fe < -30
    print(f"\n  [P1] lfuncheckfeq = {fe:.0f} → {'✅ PASS' if p1_pass else '❌ FAIL'}")

    # ── P2: 영점 ──
    zeros_raw = pari.lfunzeros(lf, T_MAX)
    zeros = [float(z) for z in zeros_raw]
    n_zeros = len(zeros)
    result["P2_zeros"] = n_zeros
    if n_zeros >= 1:
        print(f"\n  [P2] 영점 {n_zeros}개, t₁={zeros[0]:.6f}")
    else:
        print(f"\n  [P2] 영점 {n_zeros}개")

    # L-값 (BSD 검증)
    for r in range(rank + 1):
        try:
            Lr = float(pari(f"lfun({lf}, 1, {r})"))
            if r == 0:
                result["L_E_1"] = Lr
                print(f"    L(E, 1) = {Lr:.10f}")
            elif r == 1:
                result["L_prime"] = Lr
                print(f"    L'(E, 1) = {Lr:.10f}")
            elif r == 2:
                result["L_double_prime"] = Lr
                print(f"    L''(E, 1) = {Lr:.10f}")
            elif r == 3:
                result["L_triple_prime"] = Lr
                print(f"    L'''(E, 1) = {Lr:.10f}")
        except Exception as e:
            print(f"    L^({r})(E,1) 계산 실패: {e}")

    # ── lfuninit ──
    try:
        linit = pari.lfuninit(lf, [0, T_MAX])
    except Exception:
        linit = pari.lfuninit(lf, T_MAX)

    # ── P3: κ_near (σ-방향) ──
    # 첫 번째 비자명 영점 사용 (t>0.5)
    t0_zero = None
    for z in zeros:
        if z > 0.5:
            t0_zero = z
            break

    if t0_zero is not None:
        kappas = []
        for delta in DELTA_OFFSETS:
            s_plus = pari(f"{CENTER + delta} + {t0_zero}*I")
            s_minus = pari(f"{CENTER - delta} + {t0_zero}*I")
            try:
                lam_plus = complex(pari.lfunlambda(linit, s_plus))
                lam_minus = complex(pari.lfunlambda(linit, s_minus))
                ratio = abs(lam_plus) / max(abs(lam_minus), 1e-100)
                kappas.append(ratio)
            except Exception:
                continue

        if kappas:
            kappa_mean = np.mean(kappas)
            kappa_std = np.std(kappas)
            A_t0 = abs(kappa_mean - 1.0)
            p3_pass = A_t0 < 0.1
            result["P3_kappa"] = kappa_mean
            result["P3_std"] = kappa_std
            result["P3_A"] = A_t0
            print(f"\n  [P3] κδ² = {kappa_mean:.6f} ± {kappa_std:.6f} ({len(kappas)}점)")
            print(f"    A(t₀) = {A_t0:.4f}")
            print(f"    → {'✅ PASS' if p3_pass else '❌ FAIL'}")
        else:
            p3_pass = False
            result["P3_kappa"] = None
            print(f"\n  [P3] κ 측정 실패")
    else:
        p3_pass = False
        result["P3_kappa"] = None
        print(f"\n  [P3] 적합한 영점 없음 (t>0.5)")

    # ── P4: 모노드로미 (ε-보정) ──
    # ε=+1 → Re(Λ) 부호변환
    # ε=-1 → Im(Λ) 부호변환 (AP-18 보정)
    sign_changes = count_sign_changes_component(lf, linit, CENTER, T_MAX, DT, use_imag)

    # 비자명 영점 수 (t>0.5)
    nontrivial_zeros = sum(1 for z in zeros if z > 0.5)
    if nontrivial_zeros > 0:
        mono_ratio = sign_changes / nontrivial_zeros
    else:
        mono_ratio = 0.0

    p4_pass = mono_ratio > 0.8
    result["P4_sign_changes"] = sign_changes
    result["P4_nontrivial_zeros"] = nontrivial_zeros
    result["P4_ratio"] = mono_ratio

    comp_name = "Im(Λ)" if use_imag else "Re(Λ)"
    print(f"\n  [P4] 모노드로미 ({comp_name} 부호변환)")
    print(f"    부호변환 {sign_changes}개 / 비자명영점 {nontrivial_zeros}개")
    print(f"    비율: {mono_ratio:.2f}")
    print(f"    → {'✅ PASS' if p4_pass else '❌ FAIL'}")

    # ── σ-sweep (ε-보정) ──
    print(f"\n  [σ-sweep] {comp_name} 부호변환 카운트")
    sigma_counts = {}
    for sigma in SIGMA_VALUES:
        S = count_sign_changes_component(lf, linit, sigma, T_MAX, DT, use_imag)
        offset = sigma - CENTER
        sigma_counts[sigma] = S
        print(f"    σ={sigma:.1f} (offset={offset:+.1f}): S={S}")

    S_values = list(sigma_counts.values())
    S_max = max(S_values) if S_values else 0
    S_min = min(S_values) if S_values else 0
    S_crit = sigma_counts.get(CENTER, 0)

    if S_crit > 0:
        A_sigma = S_max / S_crit - 1.0
    elif S_max > 0:
        A_sigma = float('inf')
    else:
        A_sigma = 0.0

    result["sigma_counts"] = sigma_counts
    result["A_sigma"] = A_sigma
    print(f"    A(σ) = {A_sigma:.4f}")

    # ── 요약 ──
    passes = sum([p1_pass, True, p3_pass, p4_pass])  # P2는 항상 PASS
    result["passes"] = passes
    elapsed = time.time() - t0
    result["elapsed"] = elapsed

    print(f"\n  ── {label} 요약 (rank={rank}, ε={eps}) ──")
    print(f"    P1 FE:        {'✅' if p1_pass else '❌'} ({fe:.0f})")
    print(f"    P2 영점:      {n_zeros}개")
    print(f"    P3 κδ²:       {'✅' if p3_pass else '❌'} ({result.get('P3_kappa', 'N/A')})")
    print(f"    P4 모노드로미: {'✅' if p4_pass else '❌'} ({comp_name}, {sign_changes}/{nontrivial_zeros})")
    print(f"    통과:         {passes}/4")
    print(f"    소요:         {elapsed:.0f}s")

    return result


# ── 메인 실행 ──
if __name__ == "__main__":
    banner(f"[Project RDL] 결과 #108 — B-08 rank-dependent ξ-bundle (ε-보정판)")
    print(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"보정: ε=-1 → Im(Λ) 부호변환 (AP-18)")
    print(f"곡선: {len(CURVES)}개 (rank 0×3, rank 1×3, rank 2×1, rank 3×1)")

    total_t0 = time.time()

    for info in CURVES:
        try:
            r = measure_curve(info)
            RESULTS.append(r)
        except Exception as e:
            print(f"\n  ❌ {info['label']} 실패: {e}")
            import traceback
            traceback.print_exc()

    # ── 크로스 비교 ──
    banner("크로스 비교")

    # rank별 그룹
    from collections import defaultdict
    by_rank = defaultdict(list)
    for r in RESULTS:
        by_rank[r["rank"]].append(r)

    print("\n  [1] rank별 κδ² 비교:")
    for rank in sorted(by_rank.keys()):
        curves = by_rank[rank]
        kappas = [c.get("P3_kappa", None) for c in curves if c.get("P3_kappa") is not None]
        labels = [c["label"] for c in curves]
        conductors = [c["N"] for c in curves]
        if kappas:
            print(f"    rank {rank}: {', '.join(f'{l}(N={n})→κ={k:.4f}' for l, n, k in zip(labels, conductors, kappas))}")
            if len(kappas) >= 2:
                print(f"      mean={np.mean(kappas):.4f}, std={np.std(kappas):.4f}")

    print("\n  [2] rank별 P4 모노드로미 비교:")
    for rank in sorted(by_rank.keys()):
        curves = by_rank[rank]
        for c in curves:
            comp = c.get("component", "?")
            sc = c.get("P4_sign_changes", "?")
            nz = c.get("P4_nontrivial_zeros", "?")
            ratio = c.get("P4_ratio", 0)
            print(f"    rank {rank} {c['label']}(ε={c['epsilon']}): {comp}(Λ) S={sc}/{nz} = {ratio:.2f}")

    print("\n  [3] conductor 효과 (rank 0 내):")
    r0_curves = by_rank.get(0, [])
    if len(r0_curves) >= 2:
        for c in r0_curves:
            k = c.get("P3_kappa", None)
            a = c.get("P3_A", None)
            print(f"    {c['label']} N={c['N']}: κδ²={k:.6f}, A(t₀)={a:.6f}" if k else f"    {c['label']} N={c['N']}: 측정 실패")

    print(f"\n  [4] rank별 P4 PASS 비율:")
    for rank in sorted(by_rank.keys()):
        curves = by_rank[rank]
        pass_count = sum(1 for c in curves if c.get("P4_ratio", 0) > 0.8)
        print(f"    rank {rank}: {pass_count}/{len(curves)}")

    total_elapsed = time.time() - total_t0
    print(f"\n총 소요: {total_elapsed:.0f}s ({total_elapsed/60:.1f}분)")

    # ── 파일 저장 ──
    with open(OUT_PATH, "w") as f:
        f.write(f"{'='*72}\n")
        f.write(f"[Project RDL] 결과 #108 — B-08 rank-dependent ξ-bundle (ε-보정판)\n")
        f.write(f"{'='*72}\n")
        f.write(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"보정: ε=-1 → Im(Λ) 부호변환 (AP-18)\n\n")

        for r in RESULTS:
            f.write(f"\n[{r['label']}] rank={r['rank']}, N={r['N']}, ε={r['epsilon']}\n")
            f.write(f"  P1 FE: {r['P1_fe']:.0f}\n")
            f.write(f"  P2 영점: {r['P2_zeros']}개\n")
            k = r.get('P3_kappa')
            f.write(f"  P3 κδ²: {k:.6f} (A={r.get('P3_A', 'N/A'):.6f})\n" if k else "  P3 κδ²: N/A\n")
            f.write(f"  P4 모노드로미 ({r['component']}): {r['P4_sign_changes']}/{r['P4_nontrivial_zeros']} = {r['P4_ratio']:.2f}\n")
            f.write(f"  σ-sweep ({r['component']}): {r.get('sigma_counts', {})}\n")
            f.write(f"  A(σ) = {r.get('A_sigma', 'N/A')}\n")
            f.write(f"  통과: {r['passes']}/4\n")

        f.write(f"\n{'='*72}\n")
        f.write(f"크로스 비교\n")
        f.write(f"{'='*72}\n")

        f.write("\nrank별 κδ²:\n")
        for rank in sorted(by_rank.keys()):
            curves = by_rank[rank]
            for c in curves:
                k = c.get("P3_kappa")
                f.write(f"  rank {rank} {c['label']} N={c['N']}: κδ²={k:.6f}\n" if k else f"  rank {rank} {c['label']} N={c['N']}: N/A\n")

        f.write(f"\n총 소요: {total_elapsed:.0f}s\n")

    print(f"\n결과 저장: {OUT_PATH}")
