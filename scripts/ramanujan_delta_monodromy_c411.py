"""
=============================================================================
[C-411] Ramanujan Δ 함수 (weight 12, level 1) 모노드로미 TP/FP 분리 검증
=============================================================================
목표: 모노드로미 TP/FP 분리가 weight 12에서도 성립하는지 확인.
  - 기존 검증: weight 1 (ζ, χ₅), weight 2 (11a1, 37a1)
  - Ramanujan Δ는 weight 12 유일 정규화 cuspform, level 1
  - 임계선: σ = (k-1)/2 = 11/2 = 5.5 (기존 σ=1/2, σ=1과 완전히 다름)
  - "low weight에서만 성립한다"는 비판 봉쇄

핵심 수학:
  - Λ(s, Δ) = (2π)^{-s} Γ(s) L(s, Δ)
  - 대칭: Λ(s) = Λ(12-s)
  - root number ε = 1 (자기 쌍대)
  - 임계선: Re(s) = 6 (FE 대칭 중심)... 아니,
    weight k modular form의 FE: Λ(s) = (-1)^{k/2} Λ(k-s)
    Δ: k=12, Λ(s) = Λ(12-s), 대칭축 s=6, 임계선 Re(s)=6
    BUT PARI 정규화: 영점은 Im축 위 (critical strip 0 < Re(s) < 1로 정규화)
    PARI lfunlambda는 analytic normalization 사용.

  실제 PARI 동작: lfuncreate([v, 0, [0,1], 12, 1, 1])
  - weight=12, conductor=1
  - 영점의 Re(s) = (weight-1)/2 = 5.5? NO.
  - PARI는 "analytic" L-function으로 정규화: 영점이 Re(s)=1/2 위에 놓임
  - lfunzeros가 반환하는 것은 Im(s) (Re(s)=1/2 가정)
  - lfunlambda(lf, s)에서 s = 1/2 + i*t로 평가

  확인 필요: PARI가 s=1/2+it를 기대하는지, s=6+it를 기대하는지.

방법: PARI cypari2 (C-408 패턴)
  1. ramanujantau 계수 500개로 L-함수 생성
  2. lfunzeros로 영점 수집 (T≤50)
  3. FP 생성 (인접 영점 중점 + 랜덤)
  4. 폐곡선 모노드로미 (center는 PARI 정규화에 맞춤)

성공 기준:
  - TP |mono|/π > 1.5 비율 ≥ 70%
  - FP |mono|/π < 0.3 비율 ≥ 70%
  - KS p < 0.01
  - 이중기준 정밀도 > 40%
  - 4/4 → 강한 양성
=============================================================================
"""

import sys, os, time
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)
print(f"PARI 초기화: 2GB 메모리, realprecision=100")

# ─── 설정 ───
T_MAX = 50.0
N_COEFFS = 500  # ramanujantau 계수 수
MONO_RADII = [0.1, 0.01, 0.001]  # 수학자 지시: 3종
MONO_NSTEPS = 128  # 수학자 지시: 128
SEEDS = [42, 7, 123, 314, 2024]

DELTA_INFO = {
    'label': 'Ramanujan_Delta',
    'weight': 12,
    'level': 1,
    'root_number': 1,  # ε = +1 (자기 쌍대)
}

OUT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/ramanujan_delta_monodromy_c411.txt'
)


def init_delta_lfun():
    """PARI Ramanujan Delta L-함수 초기화 (ramanujantau 계수)"""
    lf = pari(f'lfuncreate([vector({N_COEFFS}, n, ramanujantau(n)), 0, [0, 1], 12, 1, 1])')
    return lf


def delta_lambda(lf, sigma, t):
    """PARI로 Λ(s) 계산. s = sigma + i*t"""
    s = pari(f"{sigma} + {t}*I")
    try:
        val = complex(pari.lfunlambda(lf, s))
        return val
    except Exception:
        return None


def detect_critical_sigma(lf, zeros):
    """
    PARI 정규화에서 실제 임계선 σ 탐지.
    영점에서 |Λ(σ+it)|가 최소가 되는 σ를 찾는다.
    후보: 0.5 (analytic normalization) vs 5.5 (arithmetic) vs 6.0 (대칭 중심)
    """
    print("\n[임계선 탐지] 영점에서 |Λ| 최소화하는 σ 탐색...")
    candidates = [0.5, 1.0, 5.5, 6.0]

    test_zeros = zeros[:min(3, len(zeros))]
    best_sigma = None
    best_avg = float('inf')

    for sigma in candidates:
        vals = []
        for t in test_zeros:
            v = delta_lambda(lf, sigma, t)
            if v is not None:
                vals.append(abs(v))
        if vals:
            avg = np.mean(vals)
            print(f"  σ={sigma}: avg|Λ| = {avg:.2e}  ({[f'{v:.2e}' for v in vals]})")
            if avg < best_avg:
                best_avg = avg
                best_sigma = sigma

    print(f"  → 최적 σ = {best_sigma} (avg|Λ| = {best_avg:.2e})")
    return best_sigma


def monodromy_contour_delta(lf, t_zero, radius=0.1, n_steps=128, center=None):
    """
    Ramanujan Delta L-함수의 점 s=center+i*t_zero 주위 폐곡선 모노드로미.
    arg(Λ(s)) 누적 변화를 계산.
    """
    if center is None:
        raise ValueError("center must be specified")

    total_delta = 0.0
    prev_arg = None
    n_valid = 0

    for k in range(n_steps + 1):
        theta = 2 * np.pi * k / n_steps
        sigma = center + radius * np.cos(theta)
        t = t_zero + radius * np.sin(theta)

        val = delta_lambda(lf, sigma, t)
        if val is None or abs(val) < 1e-300:
            continue

        curr_arg = np.angle(val)

        if prev_arg is not None:
            delta = curr_arg - prev_arg
            while delta > np.pi:
                delta -= 2 * np.pi
            while delta < -np.pi:
                delta += 2 * np.pi
            total_delta += delta

        prev_arg = curr_arg
        n_valid += 1

    return total_delta


def curvature_delta(lf, t, center, delta_h=0.01):
    """Ramanujan Delta L-함수의 곡률 κ = |Λ'/Λ|² 근사 (차분)"""
    s_center = delta_lambda(lf, center + delta_h, t)
    s_plus = delta_lambda(lf, center + delta_h, t + delta_h)
    s_minus = delta_lambda(lf, center + delta_h, t - delta_h)
    if s_center is None or s_plus is None or s_minus is None:
        return 0.0
    if abs(s_center) < 1e-300:
        return float('inf')
    deriv = (s_plus - s_minus) / (2 * delta_h * 1j)
    conn = deriv / s_center
    return abs(conn) ** 2


def generate_fp_points(true_zeros, t_min, t_max, seed, n_random=10):
    """FP 후보 생성: 인접 영점 중점 + 랜덤 비영점"""
    rng = np.random.RandomState(seed)
    fps = []

    # 1. 인접 영점 중점
    if len(true_zeros) >= 2:
        for i in range(len(true_zeros) - 1):
            mid = (true_zeros[i] + true_zeros[i + 1]) / 2.0
            fps.append(mid)

    # 2. 랜덤 (영점에서 >1.0 거리)
    min_dist = 1.0
    added = 0
    attempts = 0
    while added < n_random and attempts < n_random * 20:
        t_rand = rng.uniform(t_min, t_max)
        if len(true_zeros) > 0:
            dists = np.abs(np.array(true_zeros) - t_rand)
            if dists.min() < min_dist:
                attempts += 1
                continue
        if fps and min(abs(t_rand - f) for f in fps) < 0.5:
            attempts += 1
            continue
        fps.append(t_rand)
        added += 1
        attempts += 1

    return np.array(sorted(fps))


def main():
    START = time.time()

    print("=" * 70)
    print("[C-411] Ramanujan Δ 함수 모노드로미 TP/FP 검증")
    print(f"  weight: {DELTA_INFO['weight']}, level: {DELTA_INFO['level']}")
    print(f"  root number ε: {DELTA_INFO['root_number']}")
    print(f"  T_MAX: {T_MAX}, n_steps={MONO_NSTEPS}, n_coeffs={N_COEFFS}")
    print(f"  핵심: weight 12 → σ=5.5 또는 PARI 정규화 σ=? (자동 탐지)")
    print("=" * 70)

    # ─── 1. L-함수 초기화 ───
    print("\n[1단계] Ramanujan Delta L-함수 초기화...")
    lf = init_delta_lfun()
    print(f"  lfuncreate 성공 (ramanujantau, {N_COEFFS}개 계수)")
    print(f"  tau(1..5) = {[int(pari(f'ramanujantau({n})')) for n in range(1, 6)]}")

    # ─── 2. 영점 수집 ───
    print(f"\n[2단계] 영점 수집 (T ≤ {T_MAX})...")
    zeros_raw = pari.lfunzeros(lf, T_MAX)
    true_zeros = [float(z) for z in zeros_raw if float(z) > 1.0]
    print(f"  발견 영점 수: {len(true_zeros)}")

    if len(true_zeros) < 5:
        print("⚠️ 영점 부족. T_MAX=100으로 확장...")
        zeros_raw = pari.lfunzeros(lf, 100.0)
        true_zeros = [float(z) for z in zeros_raw if float(z) > 1.0]
        print(f"  확장 후 영점 수: {len(true_zeros)}")
        if len(true_zeros) < 5:
            print("❌ 영점 부족. 중단.")
            return

    print(f"  영점 목록 (처음 10개):")
    for i, z in enumerate(true_zeros[:10]):
        print(f"    [{i+1}] t = {z:.6f}")
    if len(true_zeros) > 10:
        print(f"    ... ({len(true_zeros)}개 총)")

    # ─── 3. 임계선 자동 탐지 ───
    critical_sigma = detect_critical_sigma(lf, true_zeros)
    print(f"\n  사용할 임계선: σ = {critical_sigma}")

    # ─── 3.5. 건전성 검사 ───
    print("\n[3.5단계] Λ 건전성 검사...")
    print(f"  Functional equation: Λ(s) = Λ({DELTA_INFO['weight']}-s), ε = {DELTA_INFO['root_number']}")

    for z in true_zeros[:5]:
        val = delta_lambda(lf, critical_sigma, z)
        if val is not None:
            print(f"  Λ({critical_sigma}+i·{z:.4f}) = {val.real:.2e} + {val.imag:.2e}i  |Λ| = {abs(val):.2e}")

    # FE 대칭 확인: Λ(s) vs Λ(k-s)
    print("\n  FE 대칭 검증 (비영점에서):")
    for t_test in [15.0, 25.0]:
        s_val = delta_lambda(lf, critical_sigma, t_test)
        k = DELTA_INFO['weight']
        # PARI 정규화: s → k-s는 PARI가 내부적으로 처리
        # analytic normalization에서는 s → 1-s
        s_conj = delta_lambda(lf, 1.0 - critical_sigma, t_test)
        if s_val is not None and s_conj is not None:
            ratio = abs(s_val) / abs(s_conj) if abs(s_conj) > 1e-300 else float('inf')
            print(f"  |Λ({critical_sigma}+{t_test}i)| = {abs(s_val):.6e}, "
                  f"|Λ({1.0-critical_sigma}+{t_test}i)| = {abs(s_conj):.6e}, "
                  f"ratio = {ratio:.4f}")

    # ─── 4. FP 생성 ───
    print("\n[4단계] FP 후보 생성...")
    all_fp = []
    for seed in SEEDS:
        fps = generate_fp_points(true_zeros, 2.0, T_MAX, seed, n_random=10)
        all_fp.extend(fps.tolist())

    # 중복 제거 (±0.5)
    all_fp = sorted(set(all_fp))
    fp_dedup = [all_fp[0]] if all_fp else []
    for f in all_fp[1:]:
        if f - fp_dedup[-1] > 0.5:
            fp_dedup.append(f)
    fp_unique = np.array(fp_dedup)

    # 영점과 겹치는 FP 제거
    fp_clean = []
    for f in fp_unique:
        if np.min(np.abs(np.array(true_zeros) - f)) > 0.5:
            fp_clean.append(f)
    fp_unique = np.array(fp_clean)
    print(f"  FP 후보: {len(fp_unique)}개 (영점에서 >0.5 거리)")

    # ─── 5. 모노드로미 계산 ───
    print("\n[5단계] 모노드로미 계산...")
    results = []

    n_tp = min(20, len(true_zeros))
    print(f"\n  [True Positives] — {n_tp}개 영점 (center=σ={critical_sigma})")
    for i in range(n_tp):
        t = true_zeros[i]
        monos_by_r = {}
        for r in MONO_RADII:
            try:
                monos_by_r[r] = monodromy_contour_delta(
                    lf, t, radius=r, n_steps=MONO_NSTEPS, center=critical_sigma)
            except Exception as e:
                print(f"    WARNING: monodromy r={r} t={t:.4f} 실패: {e}")
                monos_by_r[r] = 0.0

        kappa = curvature_delta(lf, t, critical_sigma)
        mono_main = monos_by_r[MONO_RADII[0]]
        results.append(('TP', t, mono_main, kappa, monos_by_r))
        parts = [f"mono(r={r})={monos_by_r[r]/np.pi:.4f}π" for r in MONO_RADII]
        print(f"    t={t:.4f}: {', '.join(parts)}, κ={kappa:.2e}")

    n_fp = min(30, len(fp_unique))
    print(f"\n  [False Positives] — {n_fp}개 비영점")
    for i in range(n_fp):
        t = float(fp_unique[i])
        monos_by_r = {}
        for r in MONO_RADII:
            try:
                monos_by_r[r] = monodromy_contour_delta(
                    lf, t, radius=r, n_steps=MONO_NSTEPS, center=critical_sigma)
            except Exception as e:
                print(f"    WARNING: monodromy r={r} t={t:.4f} 실패: {e}")
                monos_by_r[r] = 0.0

        kappa = curvature_delta(lf, t, critical_sigma)
        mono_main = monos_by_r[MONO_RADII[0]]
        results.append(('FP', t, mono_main, kappa, monos_by_r))
        parts = [f"mono(r={r})={monos_by_r[r]/np.pi:.4f}π" for r in MONO_RADII]
        print(f"    t={t:.4f}: {', '.join(parts)}, κ={kappa:.2e}")

    # ─── 6. 통계 ───
    print("\n" + "=" * 70)
    print("통계 요약")
    print("=" * 70)

    tp_results = [r for r in results if r[0] == 'TP']
    fp_results = [r for r in results if r[0] == 'FP']

    if not tp_results or not fp_results:
        print("❌ TP 또는 FP가 0개 — 통계 불가")
        return

    for r in MONO_RADII:
        tp_m = [abs(res[4][r]) for res in tp_results]
        fp_m = [abs(res[4][r]) for res in fp_results]
        if tp_m and fp_m:
            ks_stat, ks_p = stats.ks_2samp(tp_m, fp_m)
            print(f"\n  radius={r}:")
            print(f"    TP |mono|/π: mean={np.mean(tp_m)/np.pi:.4f}, std={np.std(tp_m)/np.pi:.4f}")
            print(f"    FP |mono|/π: mean={np.mean(fp_m)/np.pi:.4f}, std={np.std(fp_m)/np.pi:.4f}")
            print(f"    KS test: stat={ks_stat:.4f}, p={ks_p:.6e}")

    # 주 통계 (r=0.1)
    tp_monos = [abs(r[2]) for r in tp_results]
    fp_monos = [abs(r[2]) for r in fp_results]

    print(f"\n주 통계 (r={MONO_RADII[0]}, n_steps={MONO_NSTEPS}):")
    print(f"  TP |mono|/π: mean={np.mean(tp_monos)/np.pi:.4f}, std={np.std(tp_monos)/np.pi:.4f}")
    print(f"  FP |mono|/π: mean={np.mean(fp_monos)/np.pi:.4f}, std={np.std(fp_monos)/np.pi:.4f}")

    ks_stat_main, ks_p_main = stats.ks_2samp(tp_monos, fp_monos)
    print(f"  KS test: stat={ks_stat_main:.4f}, p={ks_p_main:.6e}")

    # ─── 7. 성공 기준 ───
    print("\n" + "=" * 70)
    print("성공 기준 체크")
    print("=" * 70)

    tp_high = sum(1 for m in tp_monos if m / np.pi > 1.5) / len(tp_monos)
    crit1 = tp_high >= 0.7
    print(f"  TP |mono|/π > 1.5 비율: {tp_high:.1%} ({'✅' if crit1 else '❌'} ≥ 70%)")

    fp_low = sum(1 for m in fp_monos if m / np.pi < 0.3) / len(fp_monos)
    crit2 = fp_low >= 0.7
    print(f"  FP |mono|/π < 0.3 비율: {fp_low:.1%} ({'✅' if crit2 else '❌'} ≥ 70%)")

    crit3 = ks_p_main < 0.01
    print(f"  KS p < 0.01: {'✅' if crit3 else '❌'} (p={ks_p_main:.6e})")

    mono_threshold = np.pi * 0.5
    tp_pass = sum(1 for m in tp_monos if m > mono_threshold)
    fp_pass = sum(1 for m in fp_monos if m > mono_threshold)
    if tp_pass + fp_pass > 0:
        dual_precision = tp_pass / (tp_pass + fp_pass)
    else:
        dual_precision = 0.0
    crit4 = dual_precision > 0.4 if (tp_pass + fp_pass) > 0 else False

    print(f"\n  이중 기준 (|mono| > π/2):")
    print(f"    TP 통과: {tp_pass}/{len(tp_monos)}")
    print(f"    FP 통과: {fp_pass}/{len(fp_monos)}")
    print(f"    정밀도: {dual_precision:.1%} ({'✅' if crit4 else '❌'} > 40%)")

    n_pass = sum([crit1, crit2, crit3, crit4])
    if n_pass >= 4:
        overall = "강한 양성"
    elif n_pass >= 3:
        overall = "양성"
    elif n_pass >= 2:
        overall = "약한 양성"
    else:
        overall = "음성"

    print(f"\n  *** 종합 판정: {overall} ({n_pass}/4 기준 충족) ***")

    # ─── 8. 결과 저장 ───
    elapsed = time.time() - START
    print(f"\n소요: {elapsed:.1f}초 ({elapsed/60:.1f}분)")

    with open(OUT_PATH, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("[C-411] Ramanujan Δ 함수 모노드로미 TP/FP 검증\n")
        f.write(f"weight: {DELTA_INFO['weight']}\n")
        f.write(f"level: {DELTA_INFO['level']}\n")
        f.write(f"critical_line: sigma={critical_sigma}\n")
        f.write(f"functional_eq: Lambda(s) = Lambda({DELTA_INFO['weight']}-s)\n")
        f.write(f"root_number: +{DELTA_INFO['root_number']}\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write("=" * 70 + "\n\n")

        f.write("설정:\n")
        f.write(f"  T_MAX: {T_MAX}\n")
        f.write(f"  PARI realprecision: 100\n")
        f.write(f"  ramanujantau 계수: {N_COEFFS}\n")
        f.write(f"  n_steps: {MONO_NSTEPS}\n")
        f.write(f"  반지름: {MONO_RADII}\n")
        f.write(f"  FP 시드: {SEEDS}\n")
        f.write(f"  root number ε: +{DELTA_INFO['root_number']}\n\n")

        f.write(f"영점 수: {len(true_zeros)}\n")
        f.write(f"  목록: {[f'{z:.6f}' for z in true_zeros]}\n\n")

        # 결과 테이블
        header = f"{'Type':>4} {'t':>12}"
        for r in MONO_RADII:
            header += f" {'mono_r=' + str(r):>16}"
        header += f" {'κ':>12}"
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")
        for res in results:
            line = f"{res[0]:>4} {res[1]:>12.6f}"
            for r in MONO_RADII:
                line += f" {res[4][r] / np.pi:>15.4f}π"
            line += f" {res[3]:>12.2e}"
            f.write(line + "\n")

        f.write(f"\n\n반지름별 통계:\n")
        for r in MONO_RADII:
            tp_m = [abs(res[4][r]) for res in tp_results]
            fp_m = [abs(res[4][r]) for res in fp_results]
            if tp_m and fp_m:
                ks_s, ks_p = stats.ks_2samp(tp_m, fp_m)
                f.write(f"\n  radius={r}:\n")
                f.write(f"    TP |mono|/π: {np.mean(tp_m) / np.pi:.4f} ± {np.std(tp_m) / np.pi:.4f}\n")
                f.write(f"    FP |mono|/π: {np.mean(fp_m) / np.pi:.4f} ± {np.std(fp_m) / np.pi:.4f}\n")
                f.write(f"    KS: stat={ks_s:.4f}, p={ks_p:.6e}\n")

        f.write(f"\n\n성공 기준:\n")
        f.write(f"  TP |mono|/π > 1.5 비율: {tp_high:.1%} ({'PASS' if crit1 else 'FAIL'})\n")
        f.write(f"  FP |mono|/π < 0.3 비율: {fp_low:.1%} ({'PASS' if crit2 else 'FAIL'})\n")
        f.write(f"  KS p-value: {ks_p_main:.6e} ({'PASS' if crit3 else 'FAIL'})\n")
        f.write(f"  이중기준 정밀도: {dual_precision:.1%} ({'PASS' if crit4 else 'FAIL'})\n")
        f.write(f"\n  종합 판정: {overall} ({n_pass}/4)\n")

        # 누적 통계
        f.write(f"\n누적 (C-403+406+407+408+409+410+411):\n")
        f.write(f"  L-함수: ζ + χ₅_even + χ₅_odd_complex + 11a1 + 37a1 + Ramanujan_Δ\n")
        f.write(f"  degree: 1-2, weight: 1-12, rank 0-1\n")
        f.write(f"  이전: 105 TP / 174 FP (dedicated, 5 L-함수)\n")
        f.write(f"  이번: {len(tp_results)} TP / {len(fp_results)} FP\n")
        f.write(f"  합산: {105 + len(tp_results)} TP / {174 + len(fp_results)} FP (6 L-함수)\n")

        f.write(f"\n소요: {elapsed:.1f}초\n")

    print(f"\n결과 저장: {OUT_PATH}")
    print("완료.")


if __name__ == '__main__':
    main()
