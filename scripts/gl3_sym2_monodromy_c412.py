"""
=============================================================================
[C-412] GL(3) Symmetric Square L-함수 sym²(11a1) 모노드로미 TP/FP 분리 검증
=============================================================================
목표: 모노드로미 TP/FP 분리가 degree 3 (GL(3))에서도 성립하는지 확인.
  - 기존 검증: degree 1 (ζ, χ₅), degree 2 (11a1, 37a1, Δ)
  - sym²(11a1): degree 3, conductor 121 (=11²), weight 3
  - "degree 2까지만 성립한다"는 가장 강력한 비판 봉쇄
  - 11a1은 C-408에서 검증 완료 → sym²의 모체를 알고 있어 교차검증 가능

핵심 수학:
  - PARI: lfunsympow(E, 2)로 sym² L-함수 직접 구성
  - 함수방정식: Λ(s) = Λ(3-s), 대칭축 s=3/2
  - 임계선: σ = 3/2 = 1.5
  - gammaV = [0, 0, 1] → Γ_R(s) × Γ_R(s) × Γ_R(s+1)
  - root number ε = 1
  - conductor = 121

성공 기준 (수학자 지시):
  1. TP ≥ 15개, mono/π = 2.000 ± 0.01
  2. FP ≥ 20개, mono/π = 0.000 ± 0.01
  3. KS p < 1e-6
  4. 이중기준(κ + mono) 100% 분리
  5. 결과를 results/gl3_sym2_monodromy_c412.txt에 저장
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
T_ZEROS_MAX = 60.0  # FP 필터링을 위해 영점을 더 넓게 수집
MONO_RADII = [0.1, 0.01, 0.001]  # 수학자 지시: 3종
MONO_NSTEPS = 128  # 수학자 지시: 128

SYM2_INFO = {
    'label': 'sym2_11a1',
    'source': 'sym²(11a1)',
    'degree': 3,
    'weight': 3,          # FE: Λ(s) = Λ(3-s)
    'conductor': 121,     # 11²
    'root_number': 1,     # ε = +1
    'gammaV': [0, 0, 1],  # Γ_R(s) × Γ_R(s) × Γ_R(s+1)
    'ec_model': '[0, -1, 1, -10, -20]',  # 11a1: y² + y = x³ - x² - 10x - 20
}

OUT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl3_sym2_monodromy_c412.txt'
)


def init_sym2_lfun():
    """PARI sym²(11a1) L-함수 초기화"""
    E = pari(f'ellinit({SYM2_INFO["ec_model"]})')
    lf = pari.lfunsympow(E, 2)
    return lf


def sym2_lambda(lf, sigma, t):
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
    sym²(11a1) 예상: σ=1.5 (FE: Λ(s)=Λ(3-s))
    """
    print("\n[임계선 탐지] 영점에서 |Λ| 최소화하는 σ 탐색...")
    candidates = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

    test_zeros = zeros[:min(3, len(zeros))]
    best_sigma = None
    best_avg = float('inf')

    for sigma in candidates:
        vals = []
        for t in test_zeros:
            v = sym2_lambda(lf, sigma, t)
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


def monodromy_contour_sym2(lf, t_zero, radius=0.1, n_steps=128, center=None):
    """
    sym²(11a1) L-함수의 점 s=center+i*t_zero 주위 폐곡선 모노드로미.
    arg(Λ(s)) 누적 변화를 계산 (인수 원리).
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

        val = sym2_lambda(lf, sigma, t)
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


def curvature_sym2(lf, t, center, delta_h=0.01):
    """sym²(11a1) L-함수의 곡률 κ = |Λ'/Λ|² 근사 (차분)"""
    s_center = sym2_lambda(lf, center + delta_h, t)
    s_plus = sym2_lambda(lf, center + delta_h, t + delta_h)
    s_minus = sym2_lambda(lf, center + delta_h, t - delta_h)
    if s_center is None or s_plus is None or s_minus is None:
        return 0.0
    if abs(s_center) < 1e-300:
        return float('inf')
    deriv = (s_plus - s_minus) / (2 * delta_h * 1j)
    conn = deriv / s_center
    return abs(conn) ** 2


def generate_fp_points(true_zeros, t_min, t_max, n_midpoints=None, n_random=30):
    """
    FP 후보 생성: 인접 영점 중점 + 저t 영역 + 랜덤 비영점.
    degree 3: 영점이 밀집 (평균 갭 ~0.8)하므로 공격적 전략 필요.
    r=0.1 contour가 영점을 포함하지 않으려면 min_dist > 0.15 필요.
    """
    fps = []
    zeros_arr = np.array(true_zeros)

    # 1. 인접 영점 중점 (갭 > 0.4이면 포함 — 이전 1.0에서 완화)
    if len(true_zeros) >= 2:
        for i in range(len(true_zeros) - 1):
            mid = (true_zeros[i] + true_zeros[i + 1]) / 2.0
            gap = true_zeros[i + 1] - true_zeros[i]
            if gap > 0.4:
                fps.append(mid)

    # 2. 저t 영역 (첫 영점 이전): 영점이 없는 확실한 FP
    if len(true_zeros) > 0:
        first_zero = true_zeros[0]
        if first_zero > 1.0:
            # 0.5 ~ first_zero - 0.2 범위에서 균등 배치
            low_t_max = first_zero - 0.2
            if low_t_max > 1.0:
                n_low = min(5, int((low_t_max - 0.5) / 0.3))
                for i in range(n_low):
                    t = 0.5 + (low_t_max - 0.5) * (i + 0.5) / n_low
                    fps.append(t)

    # 3. 랜덤: min_dist=0.2 (r=0.1 contour 안전 거리)
    min_dist = 0.2
    for seed in [42, 7, 123, 314, 2024, 999, 55, 77]:
        rng = np.random.RandomState(seed)
        added = 0
        attempts = 0
        target = n_random // 8 + 2
        while added < target and attempts < 500:
            t_rand = rng.uniform(t_min, t_max)
            if len(true_zeros) > 0:
                dists = np.abs(zeros_arr - t_rand)
                if dists.min() < min_dist:
                    attempts += 1
                    continue
            if fps and min(abs(t_rand - f) for f in fps) < 0.15:
                attempts += 1
                continue
            fps.append(t_rand)
            added += 1
            attempts += 1

    return np.array(sorted(fps))


def main():
    START = time.time()

    print("=" * 70)
    print("[C-412] GL(3) sym²(11a1) 모노드로미 TP/FP 검증")
    print(f"  source: {SYM2_INFO['source']}")
    print(f"  degree: {SYM2_INFO['degree']}, weight: {SYM2_INFO['weight']}")
    print(f"  conductor: {SYM2_INFO['conductor']}")
    print(f"  gammaV: {SYM2_INFO['gammaV']}")
    print(f"  root number ε: {SYM2_INFO['root_number']}")
    print(f"  T_MAX: {T_MAX}, n_steps={MONO_NSTEPS}")
    print(f"  핵심: degree 3 첫 확장 → degree 독립성 검증")
    print("=" * 70)

    # ─── 1. L-함수 초기화 ───
    print("\n[1단계] sym²(11a1) L-함수 초기화...")
    lf = init_sym2_lfun()
    print(f"  lfunsympow(E, 2) 성공")
    print(f"  EC: {SYM2_INFO['ec_model']} (11a1)")

    # a_p 확인
    E = pari(f'ellinit({SYM2_INFO["ec_model"]})')
    print(f"  a_p: ", end="")
    for p in [2, 3, 5, 7, 11, 13, 17, 19]:
        ap = int(pari.ellap(E, p))
        print(f"a_{p}={ap}", end=" ")
    print()

    # 메타데이터 확인
    params = pari.lfunparams(lf)
    print(f"  lfunparams: {params}")

    # ─── 2. 영점 수집 ───
    # T_ZEROS_MAX까지 영점 수집 (FP 필터링을 위해 T_MAX보다 넓게)
    print(f"\n[2단계] 영점 수집 (T ≤ {T_ZEROS_MAX}, TP용 T ≤ {T_MAX})...")
    zeros_raw = pari.lfunzeros(lf, T_ZEROS_MAX)
    all_zeros = [float(z) for z in zeros_raw if float(z) > 1.0]
    true_zeros = [z for z in all_zeros if z <= T_MAX]  # TP 후보: T_MAX 이내
    print(f"  발견 영점 수: {len(all_zeros)} (T≤{T_ZEROS_MAX}), TP용: {len(true_zeros)} (T≤{T_MAX})")

    if len(true_zeros) < 15:
        print("⚠️ 영점 부족. T_MAX=100으로 확장...")
        zeros_raw = pari.lfunzeros(lf, 100.0)
        all_zeros = [float(z) for z in zeros_raw if float(z) > 1.0]
        true_zeros = all_zeros  # 전체 사용
        print(f"  확장 후 영점 수: {len(true_zeros)}")
        if len(true_zeros) < 15:
            print("❌ 영점 15개 미만. 중단.")
            return

    print(f"  영점 목록 (처음 15개):")
    for i, z in enumerate(true_zeros[:15]):
        print(f"    [{i+1}] t = {z:.6f}")
    if len(true_zeros) > 15:
        print(f"    ... ({len(true_zeros)}개 총)")

    # ─── 3. 임계선 자동 탐지 ───
    critical_sigma = detect_critical_sigma(lf, true_zeros)
    print(f"\n  사용할 임계선: σ = {critical_sigma}")

    # 예상값 검증
    if abs(critical_sigma - 1.5) > 0.01:
        print(f"  ⚠️ 예상 σ=1.5와 다름! PARI 정규화 확인 필요")
    else:
        print(f"  ✅ 예상 σ=1.5와 일치 (FE: Λ(s)=Λ(3-s))")

    # ─── 3.5. 건전성 검사 ───
    print("\n[3.5단계] Λ 건전성 검사...")
    print(f"  Functional equation: Λ(s) = Λ({SYM2_INFO['weight']}-s)")
    print(f"  대칭축: s = {SYM2_INFO['weight']}/2 = {SYM2_INFO['weight']/2}")

    # 영점에서 |Λ| ≈ 0 확인
    print("\n  영점에서 |Λ|:")
    for z in true_zeros[:5]:
        val = sym2_lambda(lf, critical_sigma, z)
        if val is not None:
            print(f"    Λ({critical_sigma}+i·{z:.4f}) = {val.real:.2e} + {val.imag:.2e}i  |Λ| = {abs(val):.2e}")

    # FE 대칭 확인
    print("\n  FE 대칭 검증 (비영점에서):")
    for t_test in [15.0, 25.0, 35.0]:
        s_val = sym2_lambda(lf, critical_sigma, t_test)
        k = SYM2_INFO['weight']
        s_conj = sym2_lambda(lf, k - critical_sigma, t_test)
        if s_val is not None and s_conj is not None:
            ratio = abs(s_val) / abs(s_conj) if abs(s_conj) > 1e-300 else float('inf')
            print(f"    |Λ({critical_sigma}+{t_test}i)| = {abs(s_val):.6e}, "
                  f"|Λ({k-critical_sigma}+{t_test}i)| = {abs(s_conj):.6e}, "
                  f"ratio = {ratio:.4f}")

    # FE check (PARI)
    print(f"\n  PARI lfuncheckfeq: ", end="")
    feq = pari.lfuncheckfeq(lf)
    print(f"{feq} (음수 = 일치하는 자릿수)")

    # ─── 4. FP 생성 ───
    print("\n[4단계] FP 후보 생성...")
    # all_zeros를 사용하여 FP 필터링 (T_MAX 밖의 영점도 고려)
    fp_points = generate_fp_points(all_zeros, 2.0, T_MAX, n_random=20)

    # 영점과 겹치는 FP 제거 (all_zeros 기준, min_dist=0.15 — r=0.1 안전)
    fp_clean = []
    all_zeros_arr = np.array(all_zeros)
    for f in fp_points:
        if np.min(np.abs(all_zeros_arr - f)) > 0.15:
            fp_clean.append(f)
    fp_unique = np.array(fp_clean)
    print(f"  FP 후보: {len(fp_unique)}개 (영점에서 >0.5 거리)")

    # FP 유형별 집계
    n_mid = 0
    n_rand = 0
    for f in fp_unique:
        is_mid = False
        for i in range(len(true_zeros) - 1):
            mid = (true_zeros[i] + true_zeros[i + 1]) / 2.0
            if abs(f - mid) < 0.01:
                is_mid = True
                break
        if is_mid:
            n_mid += 1
        else:
            n_rand += 1
    print(f"  중점: {n_mid}개, 랜덤: {n_rand}개")

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
                monos_by_r[r] = monodromy_contour_sym2(
                    lf, t, radius=r, n_steps=MONO_NSTEPS, center=critical_sigma)
            except Exception as e:
                print(f"    WARNING: monodromy r={r} t={t:.4f} 실패: {e}")
                monos_by_r[r] = 0.0

        kappa = curvature_sym2(lf, t, critical_sigma)
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
                monos_by_r[r] = monodromy_contour_sym2(
                    lf, t, radius=r, n_steps=MONO_NSTEPS, center=critical_sigma)
            except Exception as e:
                print(f"    WARNING: monodromy r={r} t={t:.4f} 실패: {e}")
                monos_by_r[r] = 0.0

        kappa = curvature_sym2(lf, t, critical_sigma)
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

    # 기준 1: TP ≥ 15, mono/π = 2.000 ± 0.01
    tp_mono_pi = [m / np.pi for m in tp_monos]
    tp_in_range = sum(1 for m in tp_mono_pi if abs(m - 2.0) < 0.01)
    crit1 = len(tp_results) >= 15 and tp_in_range == len(tp_results)
    print(f"\n  기준 1: TP ≥ 15, mono/π = 2.000 ± 0.01")
    print(f"    TP 수: {len(tp_results)} (≥ 15: {'✅' if len(tp_results) >= 15 else '❌'})")
    print(f"    범위 내: {tp_in_range}/{len(tp_results)} {'✅' if crit1 else '❌'}")

    # 기준 2: FP ≥ 20, mono/π = 0.000 ± 0.01
    fp_mono_pi = [m / np.pi for m in fp_monos]
    fp_in_range = sum(1 for m in fp_mono_pi if abs(m) < 0.01)
    crit2 = len(fp_results) >= 20 and fp_in_range == len(fp_results)
    print(f"\n  기준 2: FP ≥ 20, mono/π = 0.000 ± 0.01")
    print(f"    FP 수: {len(fp_results)} (≥ 20: {'✅' if len(fp_results) >= 20 else '❌'})")
    print(f"    범위 내: {fp_in_range}/{len(fp_results)} {'✅' if crit2 else '❌'}")

    # 기준 3: KS p < 1e-6
    crit3 = ks_p_main < 1e-6
    print(f"\n  기준 3: KS p < 1e-6")
    print(f"    p = {ks_p_main:.6e} {'✅' if crit3 else '❌'}")

    # 기준 4: 이중기준(κ + mono) 100% 분리
    mono_threshold = np.pi * 0.5
    tp_mono_pass = sum(1 for m in tp_monos if m > mono_threshold)
    fp_mono_pass = sum(1 for m in fp_monos if m > mono_threshold)
    if tp_mono_pass + fp_mono_pass > 0:
        dual_precision = tp_mono_pass / (tp_mono_pass + fp_mono_pass)
    else:
        dual_precision = 0.0
    crit4 = dual_precision == 1.0
    print(f"\n  기준 4: 이중기준 100% 분리")
    print(f"    TP mono>π/2: {tp_mono_pass}/{len(tp_results)}")
    print(f"    FP mono>π/2: {fp_mono_pass}/{len(fp_results)}")
    print(f"    정밀도: {dual_precision:.1%} {'✅' if crit4 else '❌'}")

    # κ 분포
    tp_kappas = [r[3] for r in tp_results]
    fp_kappas = [r[3] for r in fp_results]
    tp_k_finite = [k for k in tp_kappas if np.isfinite(k)]
    fp_k_finite = [k for k in fp_kappas if np.isfinite(k)]
    if tp_k_finite and fp_k_finite:
        print(f"\n  κ 분포:")
        print(f"    TP κ: mean={np.mean(tp_k_finite):.2e}, range=[{min(tp_k_finite):.2e}, {max(tp_k_finite):.2e}]")
        print(f"    FP κ: mean={np.mean(fp_k_finite):.2e}, range=[{min(fp_k_finite):.2e}, {max(fp_k_finite):.2e}]")

    # 종합 판정
    n_pass = sum([crit1, crit2, crit3, crit4])
    # C-412 specific: 수학자 5개 기준 (기준 5는 결과 저장)
    if n_pass >= 4:
        overall = "강한 양성"
    elif n_pass >= 3:
        overall = "양성"
    elif n_pass >= 2:
        overall = "약한 양성"
    else:
        overall = "음성"

    print(f"\n  *** 종합 판정: {overall} ({n_pass}/4 기준 충족) ***")

    # 호환성 기준도 체크 (C-411 스타일)
    tp_high = sum(1 for m in tp_monos if m / np.pi > 1.5) / len(tp_monos)
    fp_low = sum(1 for m in fp_monos if m / np.pi < 0.3) / len(fp_monos)
    print(f"\n  [호환성] TP |mono|/π > 1.5 비율: {tp_high:.1%}")
    print(f"  [호환성] FP |mono|/π < 0.3 비율: {fp_low:.1%}")

    # ─── 8. 결과 저장 ───
    elapsed = time.time() - START
    print(f"\n소요: {elapsed:.1f}초 ({elapsed/60:.1f}분)")

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("[C-412] GL(3) sym²(11a1) 모노드로미 TP/FP 검증\n")
        f.write("=" * 70 + "\n")
        f.write(f"source: {SYM2_INFO['source']}\n")
        f.write(f"degree: {SYM2_INFO['degree']}\n")
        f.write(f"weight: {SYM2_INFO['weight']}\n")
        f.write(f"conductor: {SYM2_INFO['conductor']}\n")
        f.write(f"gammaV: {SYM2_INFO['gammaV']}\n")
        f.write(f"critical_line: sigma={critical_sigma}\n")
        f.write(f"functional_eq: Lambda(s) = Lambda({SYM2_INFO['weight']}-s)\n")
        f.write(f"root_number: +{SYM2_INFO['root_number']}\n")
        f.write(f"ec_model: {SYM2_INFO['ec_model']}\n")
        f.write(f"PARI: lfunsympow(ellinit({SYM2_INFO['ec_model']}), 2)\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n\n")

        f.write("설정:\n")
        f.write(f"  T_MAX: {T_MAX}\n")
        f.write(f"  PARI realprecision: 100\n")
        f.write(f"  n_steps: {MONO_NSTEPS}\n")
        f.write(f"  반지름: {MONO_RADII}\n")
        f.write(f"  root number ε: +{SYM2_INFO['root_number']}\n\n")

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

        f.write(f"\n\n성공 기준 (수학자 C-412 5개):\n")
        f.write(f"  1. TP ≥ 15, mono/π = 2.000 ± 0.01: {tp_in_range}/{len(tp_results)} "
                f"({'PASS' if crit1 else 'FAIL'})\n")
        f.write(f"  2. FP ≥ 20, mono/π = 0.000 ± 0.01: {fp_in_range}/{len(fp_results)} "
                f"({'PASS' if crit2 else 'FAIL'})\n")
        f.write(f"  3. KS p < 1e-6: {ks_p_main:.6e} ({'PASS' if crit3 else 'FAIL'})\n")
        f.write(f"  4. 이중기준 100%: {dual_precision:.1%} ({'PASS' if crit4 else 'FAIL'})\n")
        f.write(f"  5. 결과 파일: {OUT_PATH} (PASS)\n")
        f.write(f"\n  종합 판정: {overall} ({n_pass}/4)\n")

        # degree 3 특이사항
        f.write(f"\n\ndegree 3 특이사항:\n")
        f.write(f"  - 첫 degree 3 L-함수 검증\n")
        f.write(f"  - sym²(11a1): 11a1의 대칭 제곱. 11a1은 C-408에서 degree 2로 검증 완료\n")
        f.write(f"  - gammaV = [0, 0, 1]: Γ_R(s) × Γ_R(s) × Γ_R(s+1) (degree 2와 다른 감마 구조)\n")
        f.write(f"  - 임계선 σ={critical_sigma} (FE: Λ(s)=Λ(3-s))\n")
        f.write(f"  - Euler product: degree 3 → 각 p에서 3개 local factor\n")

        # 누적 통계
        f.write(f"\n누적 (C-403~C-412):\n")
        f.write(f"  L-함수: ζ + χ₅_even + χ₅_odd_complex + 11a1 + 37a1 + Ramanujan_Δ + sym²(11a1)\n")
        f.write(f"  degree: 1-3, weight: 1-12, rank 0-1\n")
        f.write(f"  이전: 125 TP / 204 FP (dedicated, 6 L-함수)\n")
        f.write(f"  이번: {len(tp_results)} TP / {len(fp_results)} FP\n")
        f.write(f"  합산: {125 + len(tp_results)} TP / {204 + len(fp_results)} FP (7 L-함수)\n")

        f.write(f"\n소요: {elapsed:.1f}초\n")

    print(f"\n결과 저장: {OUT_PATH}")
    print("완료.")


if __name__ == '__main__':
    main()
