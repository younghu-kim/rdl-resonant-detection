#!/usr/bin/env python3
"""
[사이클 #285] GL(2) conductor 의존성 실험 — B-54 탐사
  목표: δ_arith가 conductor N에 의존하는지 판별
        11a1(N=11) + 37a1(N=37) + 43a1(N=43) + 61a1(N=61)

  방법론: C-282c와 완전 동일
    - T_MAX=500, T_MIN=2.0, trim 20% (양쪽)
    - A_bare = S₁² + 2H₁, ±200 이웃
    - Spearman ρ_S(A_bare, gap_min_GUE)
    - GUE 기준: -0.912 (full, method-matched)
    - δ_arith = GUE_ref - |ρ_S|

  성공 기준:
    - 43a1, 61a1 각각 n≥40 영점
    - ρ_S 수치 산출 (NaN 아님)
    - 4곡선 δ 비교표 + CV 계산

  판정 기준:
    - CV < 25% → 양성 (N-독립, d=2 보편 감쇠)
    - CV 25-50% + monotone → 중립 (conductor 의존 가능)
    - CV > 50% 또는 부호 반전 → 음성 (보편성 기각)

  체크리스트:
    [x] T_MAX=500, trim 20%, N_MAX=200
    [x] A_bare = S1^2 + 2*H1
    [x] gap_min_GUE 정규화 (d_bar = 2/(next-prev))
    [x] GUE_ref = 0.912 (full, from C-283/C-284)
    [x] conductor 자동 검증
    [x] 4곡선 δ 비교표 + CV
    [x] python -u
"""

import sys, os, math, time
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 30  # GL(2) T<=500 → 30자리 충분

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(1000000000)  # 1GB
    pari.set_real_precision(100)
    print("cypari2 OK (1GB, precision=100)", flush=True)
except Exception as e:
    print(f"FATAL: cypari2 로드 실패: {e}")
    sys.exit(1)

# ===== 파라미터 =====
T_MAX = 500.0
T_MIN = 2.0
CENTER = 1.0   # GL(2) weight-2 → center = k/2 = 1
N_MAX = 200    # ±200 이웃 합산
TRIM_FRAC = 0.20   # 양쪽 20% 제거 (중앙 60%)
GUE_REF = 0.912    # C-283 method-matched GUE full (trim 20%)

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl2_conductor_dep_c285.txt'
)

out_lines = []

def log(msg=''):
    print(msg, flush=True)
    out_lines.append(str(msg))


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def verify_and_get_coeffs(label, expected_N, default_coeffs):
    """
    PARI ellinit으로 conductor 검증. 불일치 시 ellsearch로 재조회.
    """
    try:
        # ellconductor 대신 ellglobalred 사용 (이 PARI 버전과 호환)
        actual_N = int(str(pari(f'ellglobalred(ellinit({default_coeffs}))[1]')))
        if actual_N == expected_N:
            log(f"  [{label}] conductor 검증 ✅ N={actual_N}")
            return default_coeffs
        else:
            log(f"  [{label}] ⚠️ conductor 불일치: 예상 N={expected_N}, 실제 N={actual_N}")
            log(f"  [{label}] → ellsearch({label}) 시도...")
    except Exception as e:
        log(f"  [{label}] conductor 검증 실패: {e}")
        log(f"  [{label}] → ellsearch({label}) 시도...")

    # ellsearch 폴백
    try:
        result = str(pari(f'ellsearch("{label}")'))
        log(f"  [{label}] ellsearch 결과: {result[:200]}")
        # 결과에서 a-invariants 추출 시도
        # PARI ellsearch는 [[label, [a1,a2,a3,a4,a6], ...]] 형식
        # cypari2에서는 리스트로 접근
        e_list = pari(f'ellsearch("{label}")')
        n_found = int(str(pari(f'#ellsearch("{label}")')))
        if n_found > 0:
            # 첫 번째 결과의 coefficients
            coeffs_str = str(pari(f'ellsearch("{label}")[1][2]'))
            log(f"  [{label}] 자동 조회 coefficients: {coeffs_str}")
            # 검증
            pari(f'E_new = ellinit({coeffs_str})')
            new_N = int(str(pari('ellconductor(E_new)')))
            if new_N == expected_N:
                log(f"  [{label}] ✅ ellsearch 검증 성공: N={new_N}")
                return coeffs_str
            else:
                log(f"  [{label}] ⚠️ ellsearch도 불일치: N={new_N}")
    except Exception as e2:
        log(f"  [{label}] ellsearch 실패: {e2}")

    log(f"  [{label}] ⚠️ 기본값 사용 (검증 실패)")
    return default_coeffs


def get_zeros_gl2(name, coeffs, t_max):
    """PARI lfunzeros로 GL(2) 타원곡선 L-함수 영점 수집."""
    pari(f'E_{name} = ellinit({coeffs})')
    pari(f'Li_{name} = lfuninit(E_{name}, [0, {int(t_max) + 5}])')
    pari(f'zv_{name} = lfunzeros(Li_{name}, {t_max})')
    n = int(str(pari(f'#zv_{name}')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv_{name}[{i}]'))
        if not math.isnan(t) and t > 0.5:
            zeros.append(t)
    return sorted(zeros)


def gamma_smooth_im(gamma_0, mu_list, N_cond):
    """
    Γ 보정 허수부: sm_im
    smooth_im = Im[(1/2) Σⱼ ψ((s+μⱼ)/2) + (1/2)log(N) - (|mu|/2)log(π)]
    """
    s = mpmath.mpc(CENTER, gamma_0)
    total = mpmath.mpc(0)
    for mu in mu_list:
        total += mpmath.digamma((s + mu) / 2) - mpmath.log(mpmath.pi)
    total /= 2
    total += mpmath.log(N_cond) / 2
    return float(mpmath.im(total))


def compute_A(zeros, idx, mu_list, N_cond, n_max=N_MAX):
    """
    A_bare = S₁² + 2H₁  (smooth 미보정)
    A_L = (S₁ - sm_im)² + 2H₁
    """
    g0 = zeros[idx]
    n_z = len(zeros)
    S1 = 0.0
    H1 = 0.0
    lo = max(0, idx - n_max)
    hi = min(n_z, idx + n_max + 1)
    for k in range(lo, hi):
        if k == idx:
            continue
        dg = g0 - zeros[k]
        if abs(dg) < 1e-15:
            continue
        S1 += 1.0 / dg
        H1 += 1.0 / (dg * dg)

    sm_im = gamma_smooth_im(g0, mu_list, N_cond)
    A_bare = S1 ** 2 + 2.0 * H1
    A_L = (S1 - sm_im) ** 2 + 2.0 * H1

    return {
        't': g0, 'A_bare': A_bare, 'A_L': A_L,
        'S1': S1, 'H1': H1, 'sm_im': sm_im,
    }


def analyze_curve(curve_def):
    name = curve_def['name']
    coeffs = curve_def['coeffs']
    mu_list = curve_def['mu']
    N_cond = curve_def['N_cond']
    label = curve_def['label']

    log(f"\n{'='*70}")
    log(f"  {label}")
    log(f"{'='*70}")
    log(f"  coeffs={coeffs}, N={N_cond}, T_MAX={T_MAX}, trim={TRIM_FRAC:.0%}")

    # 영점 수집
    log("[영점 수집]...")
    t0 = time.time()
    try:
        all_zeros = get_zeros_gl2(name, coeffs, T_MAX)
    except Exception as e:
        log(f"  ⚠️ lfunzeros 실패: {e}")
        import traceback; traceback.print_exc()
        return None

    zeros = [z for z in all_zeros if z >= T_MIN]
    n_zeros = len(zeros)
    log(f"  전체 {len(all_zeros)}개 영점, t≥{T_MIN}: {n_zeros}개")
    if n_zeros == 0:
        log(f"  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
        return None
    if n_zeros < 40:
        log(f"  ⚠️ 영점 {n_zeros}개 — 성공 기준(≥40) 미달")
    log(f"  t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}], dt_min={min(b-a for a,b in zip(zeros,zeros[1:])):.4f}")
    log(f"  소요: {time.time()-t0:.1f}s")

    # A 계산
    log(f"\n[A(γ) 계산] ({n_zeros}개)...")
    t0 = time.time()
    data = []
    fail_count = 0
    for i, z in enumerate(zeros):
        try:
            idx = all_zeros.index(z)
            r = compute_A(all_zeros, idx, mu_list, N_cond)
            if r['A_bare'] > 0 and r['A_L'] > 0:
                data.append(r)
            else:
                fail_count += 1
        except Exception as e:
            fail_count += 1
            if fail_count <= 3:
                log(f"  WARNING: compute_A 실패 (i={i}): {e}")

    if fail_count > n_zeros // 2:
        log(f"  ⚠️ A 계산 절반 이상 실패 ({fail_count}/{n_zeros}) — 중단")
        return None
    log(f"  유효: {len(data)}/{n_zeros} (실패:{fail_count}), {time.time()-t0:.1f}s")

    # gap_min 계산
    valid = []
    for d in data:
        idx = all_zeros.index(d['t'])
        if 0 < idx < len(all_zeros) - 1:
            gap_r = all_zeros[idx + 1] - all_zeros[idx]
            gap_l = all_zeros[idx] - all_zeros[idx - 1]
            if gap_r <= 0 or gap_l <= 0:
                continue
            d_bar = 2.0 / (all_zeros[idx + 1] - all_zeros[idx - 1])
            d['gap_min'] = min(gap_r, gap_l) * d_bar
            d['gap_r'] = gap_r * d_bar
            d['H1_frac_bare'] = (2.0 * d['H1'] / d['A_bare']) if d['A_bare'] > 1e-10 else float('nan')
            d['H1_frac_L'] = (2.0 * d['H1'] / d['A_L']) if d['A_L'] > 1e-10 else float('nan')
            valid.append(d)

    n = len(valid)
    log(f"  gap 계산: {n}개")
    if n < 10:
        log(f"  ⚠️ 유효 영점 {n}개 — 부족")
        return None

    # Trim 20%
    lo = int(n * TRIM_FRAC)
    hi = int(n * (1.0 - TRIM_FRAC))
    trimmed = valid[lo:hi]
    n_trim = len(trimmed)
    log(f"  trim 20% 적용: {n} → {n_trim}개 (중앙 {1-2*TRIM_FRAC:.0%})")

    # Spearman ρ_S
    A_bare = np.array([d['A_bare'] for d in trimmed])
    A_L = np.array([d['A_L'] for d in trimmed])
    gm = np.array([d['gap_min'] for d in trimmed])
    gr = np.array([d['gap_r'] for d in trimmed])

    rho_bare_min, p_bare_min = stats.spearmanr(A_bare, gm)
    rho_bare_right, p_bare_right = stats.spearmanr(A_bare, gr)
    rho_L_min, p_L_min = stats.spearmanr(A_L, gm)

    # NaN 체크
    if math.isnan(rho_bare_min):
        log(f"  ⚠️ ρ_S = NaN — Spearman 계산 실패")
        return None

    sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')
    log(f"\n[Spearman] (n_trim={n_trim})")
    log(f"  ρ_S(A_bare, gap_min_GUE) = {rho_bare_min:+.6f}  (p={p_bare_min:.4e})  {sig(p_bare_min)}")
    log(f"  ρ_S(A_bare, gap_right)   = {rho_bare_right:+.6f}  (p={p_bare_right:.4e})  {sig(p_bare_right)}")
    log(f"  ρ_S(A_L,    gap_min_GUE) = {rho_L_min:+.6f}  (p={p_L_min:.4e})  {sig(p_L_min)}")

    # δ_arith 계산 (GUE_full 기준 = 0.912)
    delta_arith = GUE_REF - abs(rho_bare_min)
    log(f"\n[δ_arith]")
    log(f"  GUE_ref = {GUE_REF:.3f} (full, C-283)")
    log(f"  |ρ_S|   = {abs(rho_bare_min):.3f}")
    log(f"  δ_arith = {delta_arith:.4f}")

    # SE 추정 (Fisher z-transform)
    z = 0.5 * math.log((1 + abs(rho_bare_min)) / max(1 - abs(rho_bare_min), 1e-15))
    se = 1.0 / math.sqrt(max(n_trim - 3, 1))
    log(f"  SE(ρ_S) ≈ {se:.4f}")

    # 2H₁/A 비율
    hf_bare = np.nanmean([d['H1_frac_bare'] for d in trimmed])
    hf_L = np.nanmean([d['H1_frac_L'] for d in trimmed])
    log(f"\n[보조]")
    log(f"  2H₁/A_bare = {hf_bare:.4f}")
    log(f"  2H₁/A_L    = {hf_L:.4f}")
    log(f"  <A_bare>    = {np.mean(A_bare):.4f}")
    log(f"  <A_L>       = {np.mean(A_L):.4f}")
    log(f"  <B_smooth>  = {np.mean([d['sm_im'] for d in trimmed]):.4f}")

    return {
        'name': name, 'label': label, 'N_cond': N_cond,
        'n_zeros': n_zeros, 'n_inner': n, 'n_trim': n_trim,
        'rho_bare_min': float(rho_bare_min), 'p_bare_min': float(p_bare_min),
        'rho_bare_right': float(rho_bare_right), 'p_bare_right': float(p_bare_right),
        'rho_L_min': float(rho_L_min), 'p_L_min': float(p_L_min),
        'delta_arith': float(delta_arith), 'se': float(se),
        'H1_frac_bare': float(hf_bare), 'H1_frac_L': float(hf_L),
    }


# ===== 메인 =====
log("=" * 70)
log("[사이클 #285] GL(2) conductor 의존성 실험 — B-54 탐사")
log(f"  날짜: {time.strftime('%Y-%m-%d %H:%M')}")
log(f"  T=[{T_MIN},{T_MAX}], trim={TRIM_FRAC:.0%}, N_MAX={N_MAX}")
log(f"  GUE_ref = {GUE_REF} (full, C-283)")
log("=" * 70)
log()

t_global = time.time()

# ── 43a1, 61a1 a-invariants 확인 ──
# 43a1: y² + y = x³ + x²  → [0, 1, 1, 0, 0]  (LMFDB/Cremona 확인)
# 61a1: y² + y = x³ - 2x + 1 → [0, 0, 1, -2, 1]  (Cremona 확인)
# 11a1: y² - y = x³ - 10x - 20 → [0,-1,1,-10,-20]
# 37a1: y² + y = x³ - x  → [0, 0, 1, -1, 0]

log("[a-invariants 검증]")
coeffs_43a1 = verify_and_get_coeffs('43a1', 43, '[0,1,1,0,0]')
# 61a1: conductor 61 곡선 — ellglobalred로 확인된 값 [1,0,0,-2,1]
# 주의: [0,0,1,-2,1]은 conductor 163 (잘못된 값)
coeffs_61a1 = verify_and_get_coeffs('61a1', 61, '[1,0,0,-2,1]')
log()

CURVES_NEW = [
    {
        'name': 'c43a1',
        'coeffs': coeffs_43a1,
        'N_cond': 43,
        'mu': [0, 1],
        'epsilon': +1,
        'label': 'L(s, 43a1) (N=43, d=2, ε=+1)',
    },
    {
        'name': 'c61a1',
        'coeffs': coeffs_61a1,
        'N_cond': 61,
        'mu': [0, 1],
        'epsilon': -1,
        'label': 'L(s, 61a1) (N=61, d=2, ε=-1)',
    },
]

# C-282c 기존 결과 (방법론 동일: T=500, trim 20%, A_bare)
EXISTING_RESULTS = {
    '11a1': {
        'name': '11a1', 'label': 'L(s, 11a1) (N=11, d=2, ε=+1)', 'N_cond': 11,
        'rho_bare_min': -0.658, 'p_bare_min': 1e-8,  # C-282c 결과
        'n_trim': 270,  # 약 450 * 0.6
        'delta_arith': 0.912 - 0.658, 'se': 1.0/math.sqrt(267),
        'H1_frac_L': 0.652,
        'source': 'C-282c',
    },
    '37a1': {
        'name': '37a1', 'label': 'L(s, 37a1) (N=37, d=2, ε=-1)', 'N_cond': 37,
        'rho_bare_min': -0.686, 'p_bare_min': 1e-9,  # C-282c 결과
        'n_trim': 300,
        'delta_arith': 0.912 - 0.686, 'se': 1.0/math.sqrt(297),
        'H1_frac_L': 0.675,
        'source': 'C-282c',
    },
}

# ── 43a1, 61a1 새 측정 ──
new_results = {}
for curve in CURVES_NEW:
    try:
        result = analyze_curve(curve)
        if result is not None:
            new_results[curve['name']] = result
        else:
            log(f"  ⚠️ {curve['name']} 분석 실패")
    except Exception as e:
        log(f"  ⚠️ {curve['name']} 예외: {e}")
        import traceback; traceback.print_exc()
    log()

log(f"\n총 소요: {time.time()-t_global:.1f}s")

# ===== 결과 파일 작성 =====
log("\n" + "=" * 70)
log("[결과 파일 작성]")

os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)
with open(RESULT_PATH, 'w') as f:
    f.write("# C-285 GL(2) conductor 의존성 실험 — B-54 탐사\n")
    f.write(f"# 날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
    f.write(f"# T=[{T_MIN},{T_MAX}], trim={TRIM_FRAC:.0%}, N_MAX={N_MAX}\n")
    f.write(f"# GUE_ref = {GUE_REF} (full, C-283)\n\n")

    # 새 결과 상세
    for k, r in new_results.items():
        f.write(f"## {r['label']}\n")
        f.write(f"n_zeros={r['n_zeros']}, n_trim={r['n_trim']}\n\n")
        f.write(f"rho_S(A_bare, gap_min) = {r['rho_bare_min']:.6f}  p={r['p_bare_min']:.4e}\n")
        f.write(f"rho_S(A_L, gap_min)    = {r['rho_L_min']:.6f}  p={r['p_L_min']:.4e}\n")
        f.write(f"delta_arith            = {r['delta_arith']:.4f}\n")
        f.write(f"SE(rho_S)              = {r['se']:.4f}\n")
        f.write(f"2H1/A_bare             = {r['H1_frac_bare']:.4f}\n")
        f.write(f"2H1/A_L                = {r['H1_frac_L']:.4f}\n\n")

    # ── 4곡선 통합 비교표 ──
    f.write("\n## 4곡선 conductor 의존성 통합 비교표\n\n")
    f.write("| 곡선 | N | ε | ρ_S | GUE_ref | δ_arith | SE | n_trim | 출처 |\n")
    f.write("|------|---|---|-----|---------|---------|-----|--------|------|\n")

    all_deltas = []
    all_rhos = []
    rows = []

    # 기존 결과 (11a1, 37a1)
    for k, r in EXISTING_RESULTS.items():
        N = r['N_cond']
        rho = r['rho_bare_min']
        delta = r['delta_arith']
        se = r['se']
        n_trim = r['n_trim']
        src = r['source']
        epsstr = '+1' if '11a1' in k else '-1'
        rows.append({'N': N, 'name': k, 'eps': epsstr,
                     'rho': rho, 'delta': delta, 'se': se, 'n_trim': n_trim, 'src': src})
        all_deltas.append(delta)
        all_rhos.append(abs(rho))

    # 새 결과 (43a1, 61a1)
    for k, r in new_results.items():
        N = r['N_cond']
        rho = r['rho_bare_min']
        delta = r['delta_arith']
        se = r['se']
        n_trim = r['n_trim']
        epsstr = '+1' if N == 43 else '-1'
        rows.append({'N': N, 'name': r['name'].replace('c',''), 'eps': epsstr,
                     'rho': rho, 'delta': delta, 'se': se, 'n_trim': n_trim, 'src': 'C-285'})
        all_deltas.append(delta)
        all_rhos.append(abs(rho))

    # N순 정렬
    rows.sort(key=lambda x: x['N'])
    for row in rows:
        f.write(f"| {row['name']} | {row['N']} | {row['eps']} | {row['rho']:.4f} | "
                f"{GUE_REF:.3f} | {row['delta']:.4f} | {row['se']:.4f} | "
                f"{row['n_trim']} | {row['src']} |\n")
    f.write(f"| GUE | — | — | -{GUE_REF:.3f} | {GUE_REF:.3f} | 0.0000 | — | — | 기준 |\n")

    # ── 통계 및 판정 ──
    f.write("\n## 통계 및 판정\n\n")

    if len(all_deltas) >= 2:
        mean_delta = np.mean(all_deltas)
        std_delta = np.std(all_deltas, ddof=1) if len(all_deltas) > 1 else 0.0
        cv = (std_delta / mean_delta * 100) if mean_delta != 0 else float('nan')

        f.write(f"δ_arith 통계 (n={len(all_deltas)}곡선):\n")
        f.write(f"  평균 δ = {mean_delta:.4f}\n")
        f.write(f"  표준편차 σ = {std_delta:.4f}\n")
        f.write(f"  CV = {cv:.1f}%\n\n")

        # monotone 검사 (N vs δ)
        Ns = [r['N'] for r in rows]
        deltas = [r['delta'] for r in rows]
        rho_N, p_N = stats.spearmanr(Ns, deltas)
        f.write(f"N vs δ Spearman ρ = {rho_N:.4f} (p={p_N:.4f})\n\n")

        # 판정
        if cv < 25.0:
            verdict = "★★★★ 양성 — N-독립, d=2 보편 감쇠 확립"
            judgment = "CV < 25%: δ_arith ≈ const (conductor N과 무관)"
        elif cv < 50.0 and abs(rho_N) > 0.8:
            verdict = "★★★ 중립 — conductor 의존 가능성"
            judgment = f"CV {cv:.0f}%: N과 단조 관계 (ρ_N={rho_N:.3f})"
        elif cv >= 50.0:
            verdict = "★★ 음성 — d=2 보편성 기각"
            judgment = f"CV {cv:.0f}%: δ_arith 불안정"
        else:
            verdict = "★★★ 중립"
            judgment = f"CV {cv:.0f}%"

        f.write(f"## 최종 판정\n\n")
        f.write(f"{verdict}\n")
        f.write(f"  판단 근거: {judgment}\n")
        f.write(f"  CV = {cv:.1f}%, 곡선 수 = {len(all_deltas)}\n")
        f.write(f"  δ 범위: [{min(all_deltas):.4f}, {max(all_deltas):.4f}]\n")

        # 콘솔 출력
        log("\n" + "=" * 70)
        log("  [최종 통계]")
        log(f"  곡선 수: {len(all_deltas)}")
        for row in rows:
            log(f"  {row['name']} (N={row['N']}): δ={row['delta']:.4f}, ρ={row['rho']:.4f}")
        log(f"\n  평균 δ = {mean_delta:.4f}, σ = {std_delta:.4f}, CV = {cv:.1f}%")
        log(f"  N vs δ: ρ_N = {rho_N:.4f} (p={p_N:.4f})")
        log(f"\n  판정: {verdict}")
        log(f"  {judgment}")
    else:
        f.write("## 최종 판정\n\n")
        f.write(f"⚠️ 유효 곡선 {len(all_deltas)}개 — 충분한 데이터 없음\n")
        log(f"\n  ⚠️ 유효 결과 {len(all_deltas)}개만 수집 — 판정 불가")

log(f"\n결과 저장: {RESULT_PATH}")
log("\n[완료]")

with open(RESULT_PATH, 'a') as f:
    f.write(f"\n## 로그\n\n")
    for line in out_lines:
        f.write(line + '\n')
