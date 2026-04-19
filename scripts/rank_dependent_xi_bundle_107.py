#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #107 — B-08 rank-dependent ξ-bundle 검증
=============================================================================
목적: rank 0/1/2/3 타원곡선 L-함수에서 ξ-bundle 4성질이 rank에 의존하는지 측정.
      BSD 추측과의 구조적 연결 가설 (B-08)의 첫 실험적 검증.

타원곡선:
  rank 0: 11a1  [0,-1,1,-10,-20]  N=11,  ε=+1
  rank 1: 37a1  [0,0,1,-1,0]      N=37,  ε=-1
  rank 2: 389a1 [0,1,1,-2,0]      N=389, ε=+1
  rank 3: 5077a1 [0,0,1,-7,6]     N=5077,ε=-1

4성질:
  P1 (FE):       lfuncheckfeq ≤ -30
  P2 (영점):     lfunzeros(L, T) 확인
  P3 (κ_near):   κδ² ≈ 1 (δ=0.001, σ-방향)
  P4 (모노드로미): Hardy Z 부호변환 (단순 영점)

ξ-bundle A(t₀):
  A(t₀) = κδ² - 1 at σ = center + δ + i*t₀

rank-dependent 관측 목표:
  - rank ↑ → A(t₀) 변화?
  - rank ↑ → κ_near 변화?
  - rank ↑ → σ-유일성 패턴 변화?
  - s=1 영점 다중도 = rank 인지 확인 (BSD 핵심)
=============================================================================
"""

import sys, os, time
import numpy as np

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results",
    "rank_dependent_xi_bundle_107.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    lines.append(msg)
    print(msg, flush=True)

def flush_file():
    with open(OUTFILE, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")

log("=" * 72)
log("[Project RDL] 결과 #107 — B-08 rank-dependent ξ-bundle 검증")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ─── PARI 초기화 ────────────────────────────────────────────────────────
import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)
gp("default(realprecision, 100)")
log("PARI 초기화: 2GB 메모리, realprecision=100")
log()
flush_file()

# ─── 타원곡선 정의 ──────────────────────────────────────────────────────
CURVES = [
    {"name": "11a1",   "ainvs": "[0,-1,1,-10,-20]", "N": 11,   "rank": 0, "eps": +1},
    {"name": "37a1",   "ainvs": "[0,0,1,-1,0]",     "N": 37,   "rank": 1, "eps": -1},
    {"name": "389a1",  "ainvs": "[0,1,1,-2,0]",     "N": 389,  "rank": 2, "eps": +1},
    {"name": "5077a1", "ainvs": "[0,0,1,-7,6]",     "N": 5077, "rank": 3, "eps": -1},
]

CENTER = 1.0  # GL(2) weight 2: center = k/2 = 1
T_MAX = 30    # 영점 탐색 범위
XI_DELTAS = [0.0005, 0.001, 0.005, 0.01]
SIGMA_OFFSETS = [-0.3, -0.2, -0.1, 0.0, +0.1, +0.2, +0.3]
N_ZEROS_KAPPA = 5   # κ_near 측정용 영점 수
N_ZEROS_MONO = 20   # 모노드로미 측정용 영점 수

# ─── 공통 측정 함수 ─────────────────────────────────────────────────────

def measure_curve(curve):
    """타원곡선 L-함수의 4성질 + A(t₀) + σ-sweep 측정"""
    name = curve["name"]
    ainvs = curve["ainvs"]
    rank = curve["rank"]

    log(f"\n{'='*72}")
    log(f"[{name}] rank={rank}, N={curve['N']}, ε={curve['eps']}")
    log(f"{'='*72}")
    t_start = time.time()

    # L-함수 생성
    var_name = f"E_{name.replace('.','')}"
    L_name = f"L_{name.replace('.','')}"
    Linit_name = f"Li_{name.replace('.','')}"

    gp(f"{var_name} = ellinit({ainvs})")
    gp(f"{L_name} = lfuncreate({var_name})")
    log(f"  ellinit + lfuncreate 완료")

    results = {"name": name, "rank": rank, "N": curve["N"]}

    # ── P1: FE 검증 ──
    log(f"\n  [P1] 함수방정식 검증")
    try:
        fe = float(gp(f"lfuncheckfeq({L_name})"))
        results["fe"] = fe
        fe_pass = fe <= -30
        log(f"    lfuncheckfeq = {fe:.0f} → {'✅ PASS' if fe_pass else '❌ FAIL'}")
    except Exception as e:
        log(f"    FE 에러: {e}")
        results["fe"] = 0
        fe_pass = False
    flush_file()

    # ── P2: 영점 탐색 ──
    log(f"\n  [P2] 영점 탐색 (T={T_MAX})")
    try:
        zeros_str = str(gp(f"lfunzeros({L_name}, {T_MAX})"))
        zeros = [float(x) for x in zeros_str.strip("[]").split(", ") if x.strip()]
        results["n_zeros"] = len(zeros)
        results["zeros"] = zeros
        log(f"    영점 {len(zeros)}개 발견")
        if len(zeros) > 0:
            log(f"    t₁={zeros[0]:.6f}")
            if len(zeros) > 1:
                log(f"    t₂={zeros[1]:.6f}")

        # rank와 s=1 영점 관계 확인
        if rank >= 1:
            # s=1 (t=0)에서의 L값 확인
            try:
                L_at_1 = float(gp(f"real(lfun({L_name}, 1))"))
                log(f"    L(E, 1) = {L_at_1:.10f} (rank={rank} → {'≈0 예상' if rank>=1 else '≠0 예상'})")
                results["L_at_1"] = L_at_1
            except:
                pass

            # s=1에서의 도함수
            if rank >= 1:
                try:
                    Lp_at_1 = float(gp(f"real(lfun({L_name}, 1, 1))"))  # 1차 도함수
                    log(f"    L'(E, 1) = {Lp_at_1:.10f}")
                    results["Lp_at_1"] = Lp_at_1
                except:
                    pass
            if rank >= 2:
                try:
                    Lpp_at_1 = float(gp(f"real(lfun({L_name}, 1, 2))"))  # 2차 도함수
                    log(f"    L''(E, 1) = {Lpp_at_1:.10f}")
                    results["Lpp_at_1"] = Lpp_at_1
                except:
                    pass
            if rank >= 3:
                try:
                    Lppp_at_1 = float(gp(f"real(lfun({L_name}, 1, 3))"))
                    log(f"    L'''(E, 1) = {Lppp_at_1:.10f}")
                    results["Lppp_at_1"] = Lppp_at_1
                except:
                    pass
        else:
            try:
                L_at_1 = float(gp(f"real(lfun({L_name}, 1))"))
                log(f"    L(E, 1) = {L_at_1:.10f} (rank=0 → ≠0 예상)")
                results["L_at_1"] = L_at_1
            except:
                pass

    except Exception as e:
        log(f"    영점 에러: {e}")
        zeros = []
        results["n_zeros"] = 0
        results["zeros"] = []
    flush_file()

    if len(zeros) < 2:
        log(f"    ⚠️ 영점 부족 — P3/P4/A(t₀) 측정 제한")

    # ── lfuninit ──
    try:
        gp(f"{Linit_name} = lfuninit({L_name}, [{T_MAX}, 0])")
        log(f"\n  lfuninit 완료 (T={T_MAX})")
    except Exception as e:
        log(f"\n  lfuninit 에러: {e}")
        results["elapsed"] = time.time() - t_start
        flush_file()
        return results

    # ── P3: κ_near (ξ-bundle σ-방향) ──
    log(f"\n  [P3] κ_near (σ-방향, ξ-bundle)")
    test_zeros = zeros[:min(N_ZEROS_KAPPA, len(zeros))]
    kappa_delta2_list = []
    A_list = []

    for t0 in test_zeros:
        for delta in XI_DELTAS:
            sigma = CENTER + delta
            try:
                # Λ(s) and Λ(s+h) for numerical derivative
                h = delta * 0.01  # small step for derivative
                s1 = f"{sigma:.10f} + {t0:.10f}*I"
                s2 = f"{sigma+h:.10f} + {t0:.10f}*I"

                L1_str = str(gp(f"lfunlambda({Linit_name}, {s1})"))
                L2_str = str(gp(f"lfunlambda({Linit_name}, {s2})"))

                # 복소수 파싱
                def parse_complex(s):
                    s = s.strip()
                    if '*I' in s or 'I' in s:
                        return complex(s.replace('*I', 'j').replace('I', 'j').replace(' ', ''))
                    return complex(float(s), 0)

                L1 = parse_complex(L1_str)
                L2 = parse_complex(L2_str)

                if abs(L1) > 1e-50:
                    dLdS = (L2 - L1) / h
                    kappa = abs(dLdS / L1) ** 2
                    kd2 = kappa * delta**2
                    A_val = kd2 - 1
                    kappa_delta2_list.append(kd2)
                    A_list.append(A_val)
            except:
                pass

    if kappa_delta2_list:
        mean_kd2 = np.mean(kappa_delta2_list)
        std_kd2 = np.std(kappa_delta2_list)
        mean_A = np.mean(A_list)
        results["mean_kd2"] = mean_kd2
        results["std_kd2"] = std_kd2
        results["mean_A"] = mean_A
        results["n_kappa_pts"] = len(kappa_delta2_list)
        kd2_pass = abs(mean_kd2 - 1.0) < 0.01
        log(f"    κδ² = {mean_kd2:.6f} ± {std_kd2:.6f} ({len(kappa_delta2_list)}점)")
        log(f"    A(t₀) = {mean_A:.4f}")
        log(f"    → {'✅ PASS' if kd2_pass else '❌ FAIL'}")
    else:
        log(f"    κ 측정 실패")
        kd2_pass = False
    flush_file()

    # ── P4: 모노드로미 (Hardy Z 부호변환) ──
    log(f"\n  [P4] 모노드로미 (Hardy Z 부호변환)")
    mono_zeros = zeros[:min(N_ZEROS_MONO, len(zeros))]

    if len(mono_zeros) >= 2:
        # Hardy Z: t 방향 sign changes
        t_scan_min = mono_zeros[0] - 0.5
        t_scan_max = mono_zeros[-1] + 0.5
        dt = 0.02
        t_scan = np.arange(max(0.5, t_scan_min), t_scan_max, dt)

        re_vals = []
        for t in t_scan:
            try:
                val = float(gp(f"real(lfunlambda({Linit_name}, {CENTER:.6f} + {t:.6f}*I))"))
                re_vals.append(val)
            except:
                re_vals.append(0.0)

        re_arr = np.array(re_vals)
        signs = np.sign(re_arr)
        signs = signs[signs != 0]
        n_sign_changes = int(np.sum(np.diff(signs) != 0))

        # 단순 영점 비율: sign changes / zeros
        if len(mono_zeros) > 0:
            simple_ratio = n_sign_changes / len(mono_zeros)
        else:
            simple_ratio = 0

        results["n_sign_changes"] = n_sign_changes
        results["simple_ratio"] = simple_ratio
        mono_pass = simple_ratio >= 0.8
        log(f"    부호변환 {n_sign_changes}개 / 영점 {len(mono_zeros)}개")
        log(f"    단순 영점 비율: {simple_ratio:.2f}")
        log(f"    → {'✅ PASS' if mono_pass else '❌ FAIL'}")
    else:
        log(f"    영점 부족으로 측정 불가")
        mono_pass = False
    flush_file()

    # ── σ-sweep (부호변환 카운트) ──
    log(f"\n  [σ-sweep] 부호변환 카운트")
    if len(zeros) >= 3:
        t_sw_min = zeros[0] - 0.5
        t_sw_max = min(zeros[-1] + 0.5, T_MAX)
        t_sweep = np.arange(max(0.5, t_sw_min), t_sw_max, 0.1)

        sigma_results = {}
        for offset in SIGMA_OFFSETS:
            sigma = CENTER + offset
            re_vals = []
            for t in t_sweep:
                try:
                    val = float(gp(f"real(lfunlambda({Linit_name}, {sigma:.6f} + {t:.6f}*I))"))
                    re_vals.append(val)
                except:
                    re_vals.append(0.0)

            re_arr = np.array(re_vals)
            sgn = np.sign(re_arr)
            sgn = sgn[sgn != 0]
            n_sc = int(np.sum(np.diff(sgn) != 0))
            sigma_results[sigma] = n_sc
            log(f"    σ={sigma:.1f} (offset={offset:+.1f}): S={n_sc}")

        s_center = sigma_results.get(CENTER, 1)
        if s_center > 0:
            off_max = max(v for k, v in sigma_results.items() if k != CENTER)
            A_sigma = off_max / s_center - 1
            results["sigma_A"] = A_sigma
            log(f"    A(σ) = {A_sigma:.4f}")
        results["sigma_sweep"] = sigma_results
    else:
        log(f"    영점 부족으로 σ-sweep 불가")
    flush_file()

    # ── 요약 ──
    elapsed = time.time() - t_start
    results["elapsed"] = elapsed

    n_pass = sum([fe_pass, len(zeros)>=5, kd2_pass, mono_pass])
    results["n_pass"] = n_pass

    log(f"\n  ── {name} 요약 (rank={rank}) ──")
    log(f"    P1 FE:        {'✅' if fe_pass else '❌'} ({results.get('fe', 'N/A')})")
    log(f"    P2 영점:      {len(zeros)}개")
    log(f"    P3 κδ²:       {'✅' if kd2_pass else '❌'} ({results.get('mean_kd2', 'N/A')})")
    log(f"    P4 모노드로미: {'✅' if mono_pass else '❌'}")
    log(f"    A(t₀):        {results.get('mean_A', 'N/A')}")
    log(f"    통과:         {n_pass}/4")
    log(f"    소요:         {elapsed:.0f}s")
    flush_file()

    return results


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

all_results = []
for curve in CURVES:
    try:
        res = measure_curve(curve)
        all_results.append(res)
    except Exception as e:
        log(f"\n  ❌ {curve['name']} 전체 에러: {e}")
        all_results.append({"name": curve["name"], "rank": curve["rank"], "error": str(e)})
    flush_file()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 통합 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"\n{'='*72}")
log("[통합 분석] rank-dependent ξ-bundle 비교")
log(f"{'='*72}")
log()

# 요약표
log(f"  {'rank':>4} | {'곡선':<8} | {'N':>5} | {'FE':>5} | {'영점':>4} | {'κδ²':>10} | {'A(t₀)':>8} | {'mono':>5} | {'통과':>4}")
log(f"  {'-'*4}-+-{'-'*8}-+-{'-'*5}-+-{'-'*5}-+-{'-'*4}-+-{'-'*10}-+-{'-'*8}-+-{'-'*5}-+-{'-'*4}")

for res in all_results:
    if "error" in res:
        log(f"  {res['rank']:>4} | {res['name']:<8} | {'ERR':>5} | {'ERR':>5} | {'ERR':>4} | {'ERR':>10} | {'ERR':>8} | {'ERR':>5} | {'ERR':>4}")
        continue

    fe_str = f"{res.get('fe', 0):.0f}"
    kd2_str = f"{res.get('mean_kd2', 0):.6f}" if 'mean_kd2' in res else "N/A"
    A_str = f"{res.get('mean_A', 0):.4f}" if 'mean_A' in res else "N/A"
    mono_str = f"{res.get('simple_ratio', 0):.2f}" if 'simple_ratio' in res else "N/A"

    log(f"  {res['rank']:>4} | {res['name']:<8} | {res['N']:>5} | {fe_str:>5} | {res.get('n_zeros', 0):>4} | {kd2_str:>10} | {A_str:>8} | {mono_str:>5} | {res.get('n_pass', 0):>4}/4")

log()

# rank-dependent 패턴 분석
log("  rank-dependent 패턴:")
for res in all_results:
    if "error" in res:
        continue
    log(f"    rank {res['rank']} ({res['name']}): ", end="")

    # L(E,1) 관련
    if "L_at_1" in res:
        log(f"L(1)={res['L_at_1']:.6f}", end="")
    if "Lp_at_1" in res:
        log(f", L'(1)={res['Lp_at_1']:.6f}", end="")
    if "Lpp_at_1" in res:
        log(f", L''(1)={res['Lpp_at_1']:.6f}", end="")
    if "Lppp_at_1" in res:
        log(f", L'''(1)={res['Lppp_at_1']:.6f}", end="")
    log()

log()

# A(t₀) vs rank
log("  A(t₀) vs rank:")
for res in all_results:
    if "mean_A" in res:
        log(f"    rank {res['rank']}: A = {res['mean_A']:.4f}")

log()

# σ-비대칭 vs rank
log("  σ-비대칭 A(σ) vs rank:")
for res in all_results:
    if "sigma_A" in res:
        log(f"    rank {res['rank']}: A(σ) = {res['sigma_A']:.4f}")

log()

# ── 결론 ──
log(f"{'='*72}")
log("[결론] B-08 rank-dependent ξ-bundle")
log(f"{'='*72}")
log()
log("  1. 4성질 PASS/FAIL: rank에 따른 차이?")
log("  2. A(t₀): rank 증가에 따른 체계적 변화?")
log("  3. L(E,1) 도함수: BSD ord_{s=1} = rank 확인?")
log("  4. σ-sweep: rank에 따른 패턴 변화?")

total_time = sum(r.get("elapsed", 0) for r in all_results)
log(f"\n총 소요 시간: {total_time:.0f}s ({total_time/60:.1f}분)")
flush_file()
log(f"\n결과 저장: {OUTFILE}")
