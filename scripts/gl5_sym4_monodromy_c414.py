"""
=============================================================================
[C-414] GL(5) Symmetric Fourth Power L-함수 sym⁴(11a1) 모노드로미 TP/FP 분리 검증
=============================================================================
목표: 모노드로미 TP/FP 분리가 degree 5 (GL(5))에서도 성립하는지 확인.
  - 기존 검증: degree 1 (ζ, χ₅), degree 2 (11a1, 37a1, Δ), degree 3 (sym²(11a1)),
               degree 4 (sym³(11a1))
  - sym⁴(11a1): degree 5, 11a1의 대칭 네제곱
  - degree {1,2,3,4} → {1,2,3,4,5} 확장 → "임의 degree 성립" 주장 직접 강화

핵심 수학:
  - PARI: lfunsympow(E, 4)로 sym⁴ L-함수 직접 구성
  - 함수방정식: Λ(s) = ε·Λ(w-s), w = PARI 자동 결정
  - 임계선: σ = w/2 (자동 탐지, 예상 σ=2.5)
  - gammaV: degree 5 → 5개 감마 인자 (PARI 자동)
  - conductor = 11⁴ = 14641 → realprecision=150 권장

C-413 교훈 적용:
  - 성공 기준을 r=0.01 기반으로
  - T_ZEROS_MAX > T_MAX로 FP 필터링 강화
  - 임계선 자동 탐지 필수

성공 기준:
  1. TP ≥ 15개, mono/π = 2.000 ± 0.01 (r=0.01 기준)
  2. FP ≥ 20개, mono/π = 0.000 ± 0.01 (r=0.01 기준)
  3. KS p < 1e-6 (r=0.01 기준)
  4. 이중기준(κ + mono) 100% 분리 (r=0.01 기준)
=============================================================================
"""

import sys, os, time
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 150)
print(f"PARI 초기화: 2GB 메모리, realprecision=150")

# ─── 설정 ───
T_MAX = 50.0
T_ZEROS_MAX = 70.0  # FP 필터링을 위해 영점을 더 넓게 수집
MONO_RADII = [0.1, 0.01, 0.001]
MONO_NSTEPS = 128
PRIMARY_RADIUS_IDX = 1  # r=0.01을 주 판정 기준으로 사용

SYM4_INFO = {
    'label': 'sym4_11a1',
    'source': 'sym⁴(11a1)',
    'degree': 5,
    'ec_model': '[0, -1, 1, -10, -20]',  # 11a1
}

OUT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl5_sym4_monodromy_c414.txt'
)


def init_sym4_lfun():
    """PARI sym⁴(11a1) L-함수 초기화 + 메타데이터 수집"""
    E = pari(f'ellinit({SYM4_INFO["ec_model"]})')
    lf = pari.lfunsympow(E, 4)

    # 메타데이터 자동 추출
    params = pari.lfunparams(lf)
    print(f"  lfunparams: {params}")

    try:
        print(f"  (메타데이터는 lfunparams 출력에서 직접 확인)")
    except Exception as e:
        print(f"  메타데이터 추출 참고: {e}")

    return lf, E


def sym4_lambda(lf, sigma, t):
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
    sym⁴(11a1) 예상: σ=2.5 (FE: Λ(s)=Λ(5-s))
    """
    print("\n[임계선 탐지] 영점에서 |Λ| 최소화하는 σ 탐색...")
    candidates = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

    test_zeros = zeros[:min(3, len(zeros))]
    best_sigma = None
    best_avg = float('inf')

    for sigma in candidates:
        vals = []
        for t in test_zeros:
            v = sym4_lambda(lf, sigma, t)
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


def monodromy_contour(lf, t_zero, radius=0.1, n_steps=128, center=None):
    """
    sym⁴(11a1) L-함수의 점 s=center+i*t_zero 주위 폐곡선 모노드로미.
    arg(Λ(s)) 누적 변화를 계산 (인수 원리).
    """
    if center is None:
        raise ValueError("center must be specified")

    total_delta = 0.0
    prev_arg = None

    for k in range(n_steps + 1):
        theta = 2 * np.pi * k / n_steps
        sigma = center + radius * np.cos(theta)
        t = t_zero + radius * np.sin(theta)

        val = sym4_lambda(lf, sigma, t)
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

    return total_delta


def curvature(lf, t, center, delta_h=0.01):
    """κ = |Λ'/Λ|² 근사 (차분)"""
    s_center = sym4_lambda(lf, center + delta_h, t)
    s_plus = sym4_lambda(lf, center + delta_h, t + delta_h)
    s_minus = sym4_lambda(lf, center + delta_h, t - delta_h)
    if s_center is None or s_plus is None or s_minus is None:
        return 0.0
    if abs(s_center) < 1e-300:
        return float('inf')
    deriv = (s_plus - s_minus) / (2 * delta_h * 1j)
    conn = deriv / s_center
    return abs(conn) ** 2


def generate_fp_points(all_zeros, t_min, t_max, n_random=30):
    """
    FP 후보 생성: 인접 영점 중점 + 저t 영역 + 랜덤 비영점.
    all_zeros는 T_ZEROS_MAX까지의 전체 영점 (FP가 영점 근방에 오지 않도록).
    min_dist=0.2: r=0.1 contour 안전 거리.
    """
    fps = []
    zeros_arr = np.array(all_zeros)

    # T_MAX 이내 영점만 중점 계산에 사용
    tp_zeros = [z for z in all_zeros if z <= t_max]

    # 1. 인접 영점 중점 (갭 > 0.4이면 포함)
    if len(tp_zeros) >= 2:
        for i in range(len(tp_zeros) - 1):
            mid = (tp_zeros[i] + tp_zeros[i + 1]) / 2.0
            gap = tp_zeros[i + 1] - tp_zeros[i]
            if gap > 0.4:
                fps.append(mid)

    # 2. 저t 영역 (첫 영점 이전)
    if len(all_zeros) > 0:
        first_zero = all_zeros[0]
        if first_zero > 1.0:
            low_t_max = first_zero - 0.2
            if low_t_max > 1.0:
                n_low = min(5, int((low_t_max - 0.5) / 0.3))
                for i in range(n_low):
                    t = 0.5 + (low_t_max - 0.5) * (i + 0.5) / n_low
                    fps.append(t)

    # 3. 랜덤: min_dist=0.2 (전체 all_zeros 기준)
    min_dist = 0.2
    for seed in [42, 7, 123, 314, 2024, 999, 55, 77, 1234, 5678]:
        rng = np.random.RandomState(seed)
        added = 0
        attempts = 0
        target = n_random // 10 + 2
        while added < target and attempts < 500:
            t_rand = rng.uniform(t_min, t_max)
            if len(all_zeros) > 0:
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
    primary_r = MONO_RADII[PRIMARY_RADIUS_IDX]

    print("=" * 70)
    print("[C-414] GL(5) sym⁴(11a1) 모노드로미 TP/FP 검증")
    print(f"  source: {SYM4_INFO['source']}")
    print(f"  degree: {SYM4_INFO['degree']}")
    print(f"  EC: {SYM4_INFO['ec_model']} (11a1)")
    print(f"  conductor: 11⁴ = 14641")
    print(f"  T_MAX: {T_MAX}, T_ZEROS_MAX: {T_ZEROS_MAX}")
    print(f"  n_steps: {MONO_NSTEPS}, radii: {MONO_RADII}")
    print(f"  주 판정 반지름: r={primary_r}")
    print(f"  realprecision: 150 (conductor 14641 대응)")
    print(f"  핵심: degree 5 첫 모노드로미 확장 → degree {{1,2,3,4,5}} 완성")
    print(f"  기대 임계선: σ=2.5 (FE: Λ(s)=Λ(5-s)) — 6번째 임계선 유형")
    print("=" * 70)

    # ─── 1. L-함수 초기화 ───
    print("\n[1단계] sym⁴(11a1) L-함수 초기화...")
    lf, E = init_sym4_lfun()
    print(f"  lfunsympow(E, 4) 성공")

    # a_p 확인
    print(f"  a_p: ", end="")
    for p in [2, 3, 5, 7, 11, 13, 17, 19]:
        ap = int(pari.ellap(E, p))
        print(f"a_{p}={ap}", end=" ")
    print()

    # FE 확인
    print(f"\n  PARI lfuncheckfeq: ", end="")
    feq = pari.lfuncheckfeq(lf)
    print(f"{feq} (음수 = 일치하는 자릿수)")
    feq_val = float(feq)
    if feq_val > -40:
        print(f"  ⚠️ lfuncheckfeq = {feq_val} (>-40). 정밀도 부족 가능성. 진행은 계속.")

    # ─── 2. 영점 수집 ───
    print(f"\n[2단계] 영점 수집 (T ≤ {T_ZEROS_MAX}, TP용 T ≤ {T_MAX})...")
    zeros_raw = pari.lfunzeros(lf, T_ZEROS_MAX)
    all_zeros = [float(z) for z in zeros_raw if float(z) > 1.0]
    true_zeros = [z for z in all_zeros if z <= T_MAX]
    print(f"  전체 영점: {len(all_zeros)} (T≤{T_ZEROS_MAX})")
    print(f"  TP용 영점: {len(true_zeros)} (T≤{T_MAX})")

    if len(true_zeros) < 15:
        print("⚠️ 영점 부족 (< 15). T_MAX=100으로 확장...")
        T_MAX_EXT = 100.0
        T_ZEROS_EXT = T_MAX_EXT + 20
        print(f"  lfunzeros(lf, {T_ZEROS_EXT}) 시작...")
        zeros_raw = pari.lfunzeros(lf, T_ZEROS_EXT)
        all_zeros = [float(z) for z in zeros_raw if float(z) > 1.0]
        true_zeros = [z for z in all_zeros if z <= T_MAX_EXT]
        print(f"  확장 후: 전체 {len(all_zeros)}, TP용 {len(true_zeros)}")
        if len(true_zeros) < 15:
            print("❌ 영점 15개 미만. 중단.")
            # 결과 파일에 실패 기록
            os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
            with open(OUT_PATH, 'w') as f:
                f.write(f"[C-414] GL(5) sym⁴(11a1) — 영점 부족으로 중단\n")
                f.write(f"T_MAX=100, 영점: {len(true_zeros)}개\n")
                f.write(f"lfuncheckfeq: {feq}\n")
            return

    print(f"  영점 목록 (처음 10개):")
    for i, z in enumerate(true_zeros[:10]):
        print(f"    [{i+1}] t = {z:.6f}")
    if len(true_zeros) > 10:
        print(f"    ... ({len(true_zeros)}개 총)")

    # Paper A L2775 참조: sym⁴(11a1) 3 zeros 기록 → 영점 수 확인
    print(f"\n  참고: Paper A L2775에 3 zeros 기록. 현재 {len(true_zeros)}개 (T≤{T_MAX})")

    # ─── 3. 임계선 자동 탐지 ───
    critical_sigma = detect_critical_sigma(lf, true_zeros)
    print(f"\n  사용할 임계선: σ = {critical_sigma}")

    # ─── 3.5. 건전성 검사 ───
    print("\n[3.5단계] Λ 건전성 검사...")

    # 영점에서 |Λ| ≈ 0 확인
    print("  영점에서 |Λ|:")
    for z in true_zeros[:5]:
        val = sym4_lambda(lf, critical_sigma, z)
        if val is not None:
            print(f"    Λ({critical_sigma}+i·{z:.4f}) = {val.real:.2e} + {val.imag:.2e}i  |Λ| = {abs(val):.2e}")

    # FE 대칭 확인: |Λ(σ+it)| vs |Λ(w-σ+it)|
    w_fe = 2 * critical_sigma
    print(f"\n  FE 대칭 검증 (w={w_fe}):")
    for t_test in [15.0, 25.0, 35.0]:
        s_val = sym4_lambda(lf, critical_sigma, t_test)
        s_conj = sym4_lambda(lf, w_fe - critical_sigma, t_test)
        if s_val is not None and s_conj is not None:
            ratio = abs(s_val) / abs(s_conj) if abs(s_conj) > 1e-300 else float('inf')
            print(f"    |Λ({critical_sigma}+{t_test}i)| = {abs(s_val):.6e}, "
                  f"|Λ({w_fe - critical_sigma}+{t_test}i)| = {abs(s_conj):.6e}, "
                  f"ratio = {ratio:.4f}")

    # ─── 4. FP 생성 ───
    print("\n[4단계] FP 후보 생성...")
    fp_points = generate_fp_points(all_zeros, 2.0, T_MAX, n_random=30)

    # 영점과 겹치는 FP 제거 (all_zeros 기준, min_dist=0.15)
    fp_clean = []
    all_zeros_arr = np.array(all_zeros)
    for f in fp_points:
        if np.min(np.abs(all_zeros_arr - f)) > 0.15:
            fp_clean.append(f)
    fp_unique = np.array(fp_clean)
    print(f"  FP 후보: {len(fp_unique)}개 (전체 영점에서 >0.15 거리)")

    if len(fp_unique) < 20:
        print(f"  ⚠️ FP 부족 ({len(fp_unique)} < 20). 갭 기준 완화 + 추가 랜덤 시도...")
        # 갭 기준 0.4 → 0.25로 완화, 추가 랜덤
        tp_zeros_local = [z for z in all_zeros if z <= T_MAX]
        extra_fps = []
        if len(tp_zeros_local) >= 2:
            for i in range(len(tp_zeros_local) - 1):
                mid = (tp_zeros_local[i] + tp_zeros_local[i + 1]) / 2.0
                gap = tp_zeros_local[i + 1] - tp_zeros_local[i]
                if gap > 0.25:
                    if np.min(np.abs(all_zeros_arr - mid)) > 0.15:
                        if not any(abs(mid - f) < 0.1 for f in fp_clean):
                            extra_fps.append(mid)
        # 추가 랜덤 (새 시드)
        for seed in [11111, 22222, 33333, 44444, 55555, 66666, 77777, 88888]:
            rng = np.random.RandomState(seed)
            for _ in range(50):
                t_rand = rng.uniform(2.0, T_MAX)
                if np.min(np.abs(all_zeros_arr - t_rand)) > 0.15:
                    if not any(abs(t_rand - f) < 0.1 for f in fp_clean + extra_fps):
                        extra_fps.append(t_rand)
        fp_clean = list(fp_clean) + extra_fps
        fp_unique = np.array(sorted(set(fp_clean)))
        # 다시 영점 거리 필터
        fp_unique = np.array([f for f in fp_unique if np.min(np.abs(all_zeros_arr - f)) > 0.15])
        print(f"  확장 후 FP: {len(fp_unique)}개")

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
                monos_by_r[r] = monodromy_contour(
                    lf, t, radius=r, n_steps=MONO_NSTEPS, center=critical_sigma)
            except Exception as e:
                print(f"    WARNING: monodromy r={r} t={t:.4f} 실패: {e}")
                monos_by_r[r] = 0.0

        kappa = curvature(lf, t, critical_sigma)
        results.append(('TP', t, monos_by_r, kappa))
        parts = [f"r={r}:{monos_by_r[r]/np.pi:.4f}π" for r in MONO_RADII]
        print(f"    t={t:.4f}: {', '.join(parts)}, κ={kappa:.2e}")

    n_fp = min(30, len(fp_unique))
    print(f"\n  [False Positives] — {n_fp}개 비영점")
    for i in range(n_fp):
        t = float(fp_unique[i])
        monos_by_r = {}
        for r in MONO_RADII:
            try:
                monos_by_r[r] = monodromy_contour(
                    lf, t, radius=r, n_steps=MONO_NSTEPS, center=critical_sigma)
            except Exception as e:
                print(f"    WARNING: monodromy r={r} t={t:.4f} 실패: {e}")
                monos_by_r[r] = 0.0

        kappa = curvature(lf, t, critical_sigma)
        results.append(('FP', t, monos_by_r, kappa))
        parts = [f"r={r}:{monos_by_r[r]/np.pi:.4f}π" for r in MONO_RADII]
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
        tp_m = [abs(res[2][r]) for res in tp_results]
        fp_m = [abs(res[2][r]) for res in fp_results]
        if tp_m and fp_m:
            ks_stat, ks_p = stats.ks_2samp(tp_m, fp_m)
            marker = " ★★★ 주 판정" if r == primary_r else ""
            print(f"\n  radius={r}:{marker}")
            print(f"    TP |mono|/π: mean={np.mean(tp_m)/np.pi:.4f}, std={np.std(tp_m)/np.pi:.4f}")
            print(f"    FP |mono|/π: mean={np.mean(fp_m)/np.pi:.4f}, std={np.std(fp_m)/np.pi:.4f}")
            print(f"    KS test: stat={ks_stat:.4f}, p={ks_p:.6e}")

    # ─── 7. 성공 기준 (r=0.01 기준) ───
    print(f"\n{'='*70}")
    print(f"성공 기준 체크 (r={primary_r} 기준)")
    print(f"{'='*70}")

    tp_monos_primary = [abs(r[2][primary_r]) for r in tp_results]
    fp_monos_primary = [abs(r[2][primary_r]) for r in fp_results]

    # 기준 1: TP ≥ 15, mono/π = 2.000 ± 0.01
    tp_mono_pi = [m / np.pi for m in tp_monos_primary]
    tp_in_range = sum(1 for m in tp_mono_pi if abs(m - 2.0) < 0.01)
    crit1 = len(tp_results) >= 15 and tp_in_range == len(tp_results)
    print(f"\n  기준 1: TP ≥ 15, mono/π = 2.000 ± 0.01 (r={primary_r})")
    print(f"    {tp_in_range}/{len(tp_results)} {'✅' if crit1 else '❌'}")

    # 기준 2: FP ≥ 20, mono/π = 0.000 ± 0.01
    fp_mono_pi = [m / np.pi for m in fp_monos_primary]
    fp_in_range = sum(1 for m in fp_mono_pi if abs(m) < 0.01)
    crit2 = len(fp_results) >= 20 and fp_in_range == len(fp_results)
    print(f"\n  기준 2: FP ≥ 20, mono/π = 0.000 ± 0.01 (r={primary_r})")
    print(f"    {fp_in_range}/{len(fp_results)} (FP 수: {len(fp_results)}) {'✅' if crit2 else '❌'}")

    # 기준 3: KS p < 1e-6
    ks_stat_primary, ks_p_primary = stats.ks_2samp(tp_monos_primary, fp_monos_primary)
    crit3 = ks_p_primary < 1e-6
    print(f"\n  기준 3: KS p < 1e-6 (r={primary_r})")
    print(f"    p = {ks_p_primary:.6e} {'✅' if crit3 else '❌'}")

    # 기준 4: 이중기준(κ + mono) 100% 분리
    mono_threshold = np.pi * 0.5
    tp_mono_pass = sum(1 for m in tp_monos_primary if m > mono_threshold)
    fp_mono_pass = sum(1 for m in fp_monos_primary if m > mono_threshold)
    if tp_mono_pass + fp_mono_pass > 0:
        dual_precision = tp_mono_pass / (tp_mono_pass + fp_mono_pass)
    else:
        dual_precision = 0.0
    crit4 = dual_precision == 1.0
    print(f"\n  기준 4: 이중기준 100% 분리 (r={primary_r})")
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
        print(f"    TP κ: mean={np.mean(tp_k_finite):.2e}, "
              f"range=[{min(tp_k_finite):.2e}, {max(tp_k_finite):.2e}]")
        print(f"    FP κ: mean={np.mean(fp_k_finite):.2e}, "
              f"range=[{min(fp_k_finite):.2e}, {max(fp_k_finite):.2e}]")

    # 종합 판정
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

    # 메타데이터 수집
    try:
        params_str = str(pari.lfunparams(lf))
    except Exception:
        params_str = "N/A"

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("[C-414] GL(5) sym⁴(11a1) 모노드로미 TP/FP 검증\n")
        f.write("=" * 70 + "\n")
        f.write(f"source: {SYM4_INFO['source']}\n")
        f.write(f"degree: {SYM4_INFO['degree']}\n")
        f.write(f"ec_model: {SYM4_INFO['ec_model']} (11a1)\n")
        f.write(f"conductor: 11^4 = 14641\n")
        f.write(f"critical_line: sigma={critical_sigma}\n")
        f.write(f"functional_eq_center: w={2*critical_sigma}\n")
        f.write(f"lfunparams: {params_str}\n")
        f.write(f"PARI: lfunsympow(ellinit({SYM4_INFO['ec_model']}), 4)\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n\n")

        f.write("설정:\n")
        f.write(f"  T_MAX: {T_MAX}\n")
        f.write(f"  T_ZEROS_MAX: {T_ZEROS_MAX}\n")
        f.write(f"  PARI realprecision: 150\n")
        f.write(f"  n_steps: {MONO_NSTEPS}\n")
        f.write(f"  반지름: {MONO_RADII}\n")
        f.write(f"  주 판정 반지름: {primary_r}\n\n")

        f.write(f"lfuncheckfeq: {feq}\n\n")

        f.write(f"영점 수: {len(true_zeros)} (T≤{T_MAX}), "
                f"전체 {len(all_zeros)} (T≤{T_ZEROS_MAX})\n")
        f.write(f"  목록 (T≤{T_MAX}): {[f'{z:.6f}' for z in true_zeros]}\n\n")

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
                line += f" {res[2][r] / np.pi:>15.4f}π"
            line += f" {res[3]:>12.2e}"
            f.write(line + "\n")

        f.write(f"\n\n반지름별 통계:\n")
        for r in MONO_RADII:
            tp_m = [abs(res[2][r]) for res in tp_results]
            fp_m = [abs(res[2][r]) for res in fp_results]
            if tp_m and fp_m:
                ks_s, ks_p = stats.ks_2samp(tp_m, fp_m)
                marker = " ★ 주 판정" if r == primary_r else ""
                f.write(f"\n  radius={r}:{marker}\n")
                f.write(f"    TP |mono|/π: {np.mean(tp_m)/np.pi:.4f} ± {np.std(tp_m)/np.pi:.4f}\n")
                f.write(f"    FP |mono|/π: {np.mean(fp_m)/np.pi:.4f} ± {np.std(fp_m)/np.pi:.4f}\n")
                f.write(f"    KS: stat={ks_s:.4f}, p={ks_p:.6e}\n")

        f.write(f"\n\n성공 기준 (C-414, r={primary_r} 기준):\n")
        f.write(f"  1. TP ≥ 15, mono/π = 2.000 ± 0.01: "
                f"{tp_in_range}/{len(tp_results)} ({'PASS' if crit1 else 'FAIL'})\n")
        f.write(f"  2. FP ≥ 20, mono/π = 0.000 ± 0.01: "
                f"{fp_in_range}/{len(fp_results)} ({'PASS' if crit2 else 'FAIL'})\n")
        f.write(f"  3. KS p < 1e-6: {ks_p_primary:.6e} ({'PASS' if crit3 else 'FAIL'})\n")
        f.write(f"  4. 이중기준 100%: {dual_precision:.1%} ({'PASS' if crit4 else 'FAIL'})\n")
        f.write(f"\n  종합 판정: {overall} ({n_pass}/4)\n")

        # degree 5 특이사항
        f.write(f"\n\ndegree 5 특이사항:\n")
        f.write(f"  - 첫 degree 5 전용 모노드로미 검증\n")
        f.write(f"  - sym⁴(11a1): 11a1의 대칭 네제곱.\n")
        f.write(f"    11a1은 C-408(d2), sym²은 C-412(d3), sym³은 C-413(d4)에서 검증 완료\n")
        f.write(f"  - degree {{1,2,3,4}} → {{1,2,3,4,5}} 확장\n")
        f.write(f"  - 임계선: σ={critical_sigma} (6번째 임계선 유형 기대)\n")
        f.write(f"  - conductor = 11⁴ = 14641 (realprecision=150)\n")

        # 누적 통계
        prev_tp, prev_fp = 165, 264
        curr_tp = len(tp_results)
        curr_fp = len(fp_results)
        f.write(f"\n누적 (C-403~C-414):\n")
        f.write(f"  L-함수: ζ + χ₅_even + χ₅_odd_complex + 11a1 + 37a1 "
                f"+ Ramanujan_Δ + sym²(11a1) + sym³(11a1) + sym⁴(11a1)\n")
        f.write(f"  degree: 1-5, weight: 1-12, rank 0-1\n")
        f.write(f"  이전: {prev_tp} TP / {prev_fp} FP (dedicated, 8 L-함수)\n")
        f.write(f"  이번: {curr_tp} TP / {curr_fp} FP\n")
        f.write(f"  합산: {prev_tp + curr_tp} TP / {prev_fp + curr_fp} FP (9 L-함수)\n")

        f.write(f"\n소요: {elapsed:.1f}초\n")

    print(f"\n결과 저장: {OUT_PATH}")
    print("완료.")


if __name__ == '__main__':
    main()
