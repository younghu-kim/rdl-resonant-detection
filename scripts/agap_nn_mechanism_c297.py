#!/usr/bin/env python3
"""
[C-297] A-gap 반상관 메커니즘 분해: 최근접 이웃(NN) vs 원거리(Rest)

핵심 아이디어:
  A_bare = S1² + 2·H1  where  S1 = Σ_j 1/(t_n - t_{n±j}), H1 = Σ_j 1/(t_n - t_{n±j})²

  j=1(최근접 이웃) 기여를 A_NN, j≥2 기여를 A_rest로 분리.
  만약 A-gap 반상관이 NN 항에서 발생한다면, 이는 Hadamard 곱의 직접적 귀결.

측정:
  1. A_NN / A_bare 비율 분포 (NN 기여율)
  2. ρ(A_bare, gap), ρ(A_NN, gap), ρ(A_rest, gap) — 반상관의 원천 식별
  3. ρ(H1_NN, gap), ρ(S1_NN², gap) — H1 vs S1² 중 어느 것이 주도적?
  4. 8/gap² vs A_bare 상관 (단순 모델 예측력)
  5. W=5,20,100,500에서 분해 안정성
  6. ζ(s) vs GUE 비교

체크리스트:
  [x] 이론적 d_bar = log(t/(2π))/(2π) for ζ(s)
  [x] GUE: 반원 CDF unfolding
  [x] trim 20%
  [x] python -u
  [x] Spearman only
"""

import sys, os, time
import numpy as np
from scipy import stats

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(1024 * 10**6)
    pari.set_real_precision(80)
    print("cypari2 OK", flush=True)
except Exception as e:
    print(f"FATAL: {e}", flush=True)
    sys.exit(1)

# ── 설정 ──────────────────────────────────────────────────────────
T_MAX = 2000.0
TRIM_FRAC = 0.20
W_LIST = [5, 20, 100, 500]
GUE_N = 2000
GUE_ENS = 30
GUE_SEED = 42

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/agap_nn_mechanism_c297.txt'
)

out_f = open(RESULT_PATH, 'w')
def log(msg=''):
    print(msg, flush=True)
    out_f.write(msg + '\n')
    out_f.flush()


# ── ζ(s) 영점 ────────────────────────────────────────────────────
def get_zeta_zeros():
    log(f"[ζ(s)] 영점 수집 t∈(0, {T_MAX}] ...")
    t0 = time.time()
    pari('L_zeta = lfuncreate(1)')
    pari(f'Li_zeta = lfuninit(L_zeta, [0, {int(T_MAX) + 5}])')
    pari(f'zv_zeta = lfunzeros(Li_zeta, {T_MAX})')
    n = int(str(pari('#zv_zeta')))
    zeros = []
    for i in range(1, n + 1):
        s = str(pari(f'zv_zeta[{i}]')).strip().replace(' E', 'e').replace('E ', 'e')
        try:
            t = float(s)
            if t > 0.5:
                zeros.append(t)
        except ValueError:
            pass
    zeros = np.array(sorted(zeros))
    log(f"  {len(zeros)}개 영점, t∈[{zeros[0]:.3f}, {zeros[-1]:.3f}], {time.time()-t0:.1f}s")
    return zeros


# ── GUE ──────────────────────────────────────────────────────────
def generate_gue_eigvals(N, seed=None):
    rng = np.random.default_rng(seed)
    real = rng.standard_normal((N, N))
    imag = rng.standard_normal((N, N))
    A = (real + 1j * imag) / np.sqrt(2.0)
    H = (A + A.conj().T) / (2.0 * np.sqrt(N))
    return np.linalg.eigvalsh(H).real


def semicircle_cdf(x):
    x = np.clip(x, -2.0 + 1e-12, 2.0 - 1e-12)
    return (np.arcsin(x / 2.0) + (x / 2.0) * np.sqrt(
        np.maximum(0.0, 1.0 - x**2 / 4.0))) / np.pi + 0.5


def unfold(lambdas):
    return len(lambdas) * semicircle_cdf(lambdas)


# ── 분해된 A_bare 계산 ──────────────────────────────────────────
def compute_A_decomposed(x, W):
    """A_bare를 NN(j=1)과 Rest(j≥2)로 분해.

    Returns: dict with keys:
      A_bare, A_NN, A_rest,
      S1_NN, S1_rest, S1,
      H1_NN, H1_rest, H1,
      gap_min, gap_l, gap_r,
      indices
    """
    n = len(x)
    if n < 2 * W + 1:
        return None

    valid_lo = W
    valid_hi = n - W
    n_valid = valid_hi - valid_lo
    if n_valid < 1:
        return None

    # NN (j=1) 기여
    diff_r1 = x[valid_lo + 1 : valid_hi + 1] - x[valid_lo : valid_hi]
    diff_l1 = x[valid_lo : valid_hi] - x[valid_lo - 1 : valid_hi - 1]
    diff_r1 = np.where(np.abs(diff_r1) < 1e-15, 1e-15, diff_r1)
    diff_l1 = np.where(np.abs(diff_l1) < 1e-15, 1e-15, diff_l1)

    S1_NN = -1.0 / diff_r1 + 1.0 / diff_l1
    H1_NN = 1.0 / diff_r1**2 + 1.0 / diff_l1**2

    # Rest (j≥2) 기여
    S1_rest = np.zeros(n_valid)
    H1_rest = np.zeros(n_valid)

    for j in range(2, W + 1):
        diff_r = x[valid_lo + j : valid_hi + j] - x[valid_lo : valid_hi]
        diff_l = x[valid_lo : valid_hi] - x[valid_lo - j : valid_hi - j]
        diff_r = np.where(np.abs(diff_r) < 1e-15, 1e-15, diff_r)
        diff_l = np.where(np.abs(diff_l) < 1e-15, 1e-15, diff_l)
        S1_rest += -1.0 / diff_r + 1.0 / diff_l
        H1_rest += 1.0 / diff_r**2 + 1.0 / diff_l**2

    S1 = S1_NN + S1_rest
    H1 = H1_NN + H1_rest
    A_bare = S1**2 + 2.0 * H1
    A_NN = S1_NN**2 + 2.0 * H1_NN
    A_rest = A_bare - A_NN  # 교차항 + rest 자체

    # gap
    gap_l = diff_l1.copy()
    gap_r = diff_r1.copy()
    gap_min = np.minimum(gap_l, gap_r)

    indices = np.arange(valid_lo, valid_hi)

    return {
        'A_bare': A_bare, 'A_NN': A_NN, 'A_rest': A_rest,
        'S1': S1, 'S1_NN': S1_NN, 'S1_rest': S1_rest,
        'H1': H1, 'H1_NN': H1_NN, 'H1_rest': H1_rest,
        'gap_min': gap_min, 'gap_l': gap_l, 'gap_r': gap_r,
        'indices': indices,
    }


def apply_trim_and_density(data, x, density_func=None):
    """밀도 정규화 + trim 적용."""
    A_bare = data['A_bare'].copy()
    A_NN = data['A_NN'].copy()
    A_rest = data['A_rest'].copy()
    H1_NN = data['H1_NN'].copy()
    S1_NN_sq = data['S1_NN']**2
    gap_min = data['gap_min'].copy()
    indices = data['indices']

    # 밀도 정규화
    if density_func is not None:
        d_bar = np.array([density_func(xi) for xi in x[indices]])
        gap_min = gap_min * d_bar

    # trim
    n_pts = len(A_bare)
    trim_lo = int(n_pts * TRIM_FRAC)
    trim_hi = n_pts - int(n_pts * TRIM_FRAC)
    if trim_hi - trim_lo < 20:
        return None

    sl = slice(trim_lo, trim_hi)
    result = {
        'A_bare': A_bare[sl], 'A_NN': A_NN[sl], 'A_rest': A_rest[sl],
        'H1_NN': H1_NN[sl], 'S1_NN_sq': S1_NN_sq[sl],
        'gap_min': gap_min[sl],
    }

    # NaN/inf 필터
    mask = np.ones(len(result['A_bare']), dtype=bool)
    for k in result:
        mask &= np.isfinite(result[k])
    mask &= (result['A_bare'] > 0) & (result['gap_min'] > 0)

    return {k: v[mask] for k, v in result.items()}


def analyze_decomposition(data, label, W):
    """분해 분석 결과 출력."""
    if data is None or len(data['A_bare']) < 20:
        log(f"  [{label}] W={W}: 데이터 부족")
        return None

    A_bare = data['A_bare']
    A_NN = data['A_NN']
    A_rest = data['A_rest']
    H1_NN = data['H1_NN']
    S1_NN_sq = data['S1_NN_sq']
    gap = data['gap_min']
    n = len(A_bare)

    # 1. NN 기여율
    nn_ratio = A_NN / A_bare
    nn_ratio_valid = nn_ratio[np.isfinite(nn_ratio) & (nn_ratio > 0) & (nn_ratio < 10)]

    # 2. 상관 계수
    rho_full, p_full = stats.spearmanr(A_bare, gap)
    rho_nn, p_nn = stats.spearmanr(A_NN, gap)
    rho_rest, p_rest = stats.spearmanr(A_rest, gap)
    rho_h1nn, p_h1nn = stats.spearmanr(H1_NN, gap)
    rho_s1nn, p_s1nn = stats.spearmanr(S1_NN_sq, gap)

    # 3. 단순 모델: 2/gap² ≈ H1_NN ≈ A_bare?
    simple_model = 2.0 / gap**2
    rho_simple, p_simple = stats.spearmanr(simple_model, A_bare)

    # 4. H1_NN이 H1 전체의 몇 %인지
    H1_total = data.get('H1_NN', H1_NN)  # rest가 없으면 NN만
    h1_ratio = np.mean(H1_NN) / np.mean(H1_NN + 1e-30)  # 전체는 나중에

    log(f"  [{label}] W={W}: n={n}")
    log(f"    NN기여율: mean={np.mean(nn_ratio_valid):.4f}, "
        f"median={np.median(nn_ratio_valid):.4f}, "
        f"std={np.std(nn_ratio_valid):.4f}")
    log(f"    ρ(A_bare, gap) = {rho_full:.4f}  (p={p_full:.2e})")
    log(f"    ρ(A_NN,   gap) = {rho_nn:.4f}  (p={p_nn:.2e})")
    log(f"    ρ(A_rest, gap) = {rho_rest:.4f}  (p={p_rest:.2e})")
    log(f"    ρ(H1_NN,  gap) = {rho_h1nn:.4f}  (p={p_h1nn:.2e})")
    log(f"    ρ(S1_NN², gap) = {rho_s1nn:.4f}  (p={p_s1nn:.2e})")
    log(f"    ρ(2/gap², A_bare) = {rho_simple:.4f}  (단순모델 예측력)")

    return {
        'n': n,
        'nn_ratio_mean': np.mean(nn_ratio_valid),
        'nn_ratio_median': np.median(nn_ratio_valid),
        'rho_full': rho_full, 'p_full': p_full,
        'rho_nn': rho_nn, 'p_nn': p_nn,
        'rho_rest': rho_rest, 'p_rest': p_rest,
        'rho_h1nn': rho_h1nn,
        'rho_s1nn': rho_s1nn,
        'rho_simple': rho_simple,
    }


# ── ζ(s) 밀도 ───────────────────────────────────────────────────
def zeta_density(t):
    if t < 2.0:
        return 1.0
    return np.log(t / (2.0 * np.pi)) / (2.0 * np.pi)


# ══════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════
log("=" * 70)
log("[C-297] A-gap 분해: 최근접 이웃(NN) 메커니즘 정량화")
log("=" * 70)

t_start = time.time()

# ── ζ(s) ─────────────────────────────────────────────────────────
zeros_zeta = get_zeta_zeros()

log()
log("=" * 70)
log("  [Phase 1] ζ(s) — W별 A_bare 분해")
log("=" * 70)

zeta_results = {}
for W in W_LIST:
    log()
    data = compute_A_decomposed(zeros_zeta, W)
    if data is None:
        log(f"  W={W}: 데이터 부족, 스킵")
        continue
    processed = apply_trim_and_density(data, zeros_zeta, zeta_density)
    r = analyze_decomposition(processed, "ζ", W)
    if r:
        zeta_results[W] = r

# ── GUE ──────────────────────────────────────────────────────────
log()
log("=" * 70)
log(f"  [Phase 2] GUE — {GUE_ENS}앙상블 × N={GUE_N}")
log("=" * 70)

gue_results = {}
for W in W_LIST:
    log()
    log(f"--- GUE W={W} ---")
    all_rho_full = []
    all_rho_nn = []
    all_rho_rest = []
    all_nn_ratio = []

    for e in range(GUE_ENS):
        eigvals = generate_gue_eigvals(GUE_N, seed=GUE_SEED + e)
        unf = unfold(eigvals)
        # 중앙 80% 사용
        n_eig = len(unf)
        lo = int(n_eig * 0.1)
        hi = int(n_eig * 0.9)
        unf_core = np.sort(unf[lo:hi])

        if len(unf_core) < 2 * W + 10:
            continue

        data = compute_A_decomposed(unf_core, W)
        if data is None:
            continue
        processed = apply_trim_and_density(data, unf_core, density_func=None)
        if processed is None or len(processed['A_bare']) < 20:
            continue

        A_b = processed['A_bare']
        A_n = processed['A_NN']
        A_r = processed['A_rest']
        gap = processed['gap_min']

        r_f, _ = stats.spearmanr(A_b, gap)
        r_n, _ = stats.spearmanr(A_n, gap)
        r_r, _ = stats.spearmanr(A_r, gap)
        nn_rat = np.mean(A_n / A_b)

        all_rho_full.append(r_f)
        all_rho_nn.append(r_n)
        all_rho_rest.append(r_r)
        all_nn_ratio.append(nn_rat)

    if len(all_rho_full) > 0:
        r_f_m, r_f_s = np.mean(all_rho_full), np.std(all_rho_full)
        r_n_m, r_n_s = np.mean(all_rho_nn), np.std(all_rho_nn)
        r_r_m, r_r_s = np.mean(all_rho_rest), np.std(all_rho_rest)
        nn_m = np.mean(all_nn_ratio)
        log(f"  [GUE] W={W}: {len(all_rho_full)}/{GUE_ENS}앙상블")
        log(f"    NN기여율: mean={nn_m:.4f}")
        log(f"    ρ(A_bare, gap) = {r_f_m:.4f} ± {r_f_s:.4f}")
        log(f"    ρ(A_NN,   gap) = {r_n_m:.4f} ± {r_n_s:.4f}")
        log(f"    ρ(A_rest, gap) = {r_r_m:.4f} ± {r_r_s:.4f}")
        gue_results[W] = {
            'nn_ratio_mean': nn_m,
            'rho_full': r_f_m, 'rho_full_std': r_f_s,
            'rho_nn': r_n_m, 'rho_nn_std': r_n_s,
            'rho_rest': r_r_m, 'rho_rest_std': r_r_s,
        }

# ── 비교표 ───────────────────────────────────────────────────────
log()
log("=" * 70)
log("  [Phase 3] ζ(s) vs GUE 비교표")
log("=" * 70)
log()
log(f"{'W':>6}  {'NN%_ζ':>8}  {'NN%_GUE':>8}  "
    f"{'ρ(full)_ζ':>10}  {'ρ(NN)_ζ':>10}  {'ρ(rest)_ζ':>10}  "
    f"{'ρ(full)_G':>10}  {'ρ(NN)_G':>10}  {'ρ(rest)_G':>10}")
log("─" * 100)

for W in W_LIST:
    z = zeta_results.get(W, {})
    g = gue_results.get(W, {})
    log(f"{W:6d}  "
        f"{z.get('nn_ratio_mean',float('nan'))*100:7.1f}%  "
        f"{g.get('nn_ratio_mean',float('nan'))*100:7.1f}%  "
        f"{z.get('rho_full',float('nan')):10.4f}  "
        f"{z.get('rho_nn',float('nan')):10.4f}  "
        f"{z.get('rho_rest',float('nan')):10.4f}  "
        f"{g.get('rho_full',float('nan')):10.4f}  "
        f"{g.get('rho_nn',float('nan')):10.4f}  "
        f"{g.get('rho_rest',float('nan')):10.4f}")

# ── 해석 ─────────────────────────────────────────────────────────
log()
log("=" * 70)
log("  [Phase 4] 해석 + 정리 후보")
log("=" * 70)

# NN 기여율이 높고 ρ(rest,gap)≈0이면 → 반상관은 NN 메커니즘
if zeta_results:
    W_ref = max(zeta_results.keys())
    z_ref = zeta_results[W_ref]
    nn_pct = z_ref['nn_ratio_mean'] * 100
    rho_rest = z_ref['rho_rest']
    rho_full = z_ref['rho_full']
    rho_nn = z_ref['rho_nn']

    log()
    log(f"  참조: W={W_ref}")
    log(f"  NN 기여율: {nn_pct:.1f}%")
    log(f"  ρ(A_bare, gap) = {rho_full:.4f}")
    log(f"  ρ(A_NN, gap)   = {rho_nn:.4f}")
    log(f"  ρ(A_rest, gap) = {rho_rest:.4f}")
    log()

    # 판정
    if abs(rho_rest) < 0.3 and abs(rho_nn) > 0.7:
        log("  ★★★★★ A-gap 반상관은 최근접 이웃 메커니즘에 의해 지배됨")
        log("  → Hadamard 곱의 1/(t_n - t_{n±1})² 항이 상관을 발생시킴")
        log("  → 원거리 항 A_rest는 gap과 약상관 (정보 추가 미미)")
        verdict = "양성"
    elif abs(rho_rest) > 0.5:
        log("  ★★ 원거리 항도 유의미한 상관 기여 → NN만으로 설명 불충분")
        verdict = "취약"
    else:
        log("  ★★★ 혼합: NN 지배적이나 원거리도 일부 기여")
        verdict = "조건부 양성"

    log()
    log(f"  판정: {verdict}")

    # 정리 후보
    log()
    log("  [정리 후보]")
    log("  Proposition (Nearest-Neighbor Dominance):")
    log("    Let {γ_n} be zeros of ξ(s) on the critical line.")
    log("    Define A_bare(n) = S1(n)² + 2·H1(n) with window W,")
    log("    and decompose A_bare = A_NN + A_cross + A_rest.")
    log("    Then:")
    log("    (a) A_NN/A_bare → known fraction as W→∞")
    log("    (b) ρ_S(A_NN, gap) ≈ ρ_S(A_bare, gap)")
    log("    (c) ρ_S(A_rest, gap) ≈ 0")
    log("    implying the A-gap anti-correlation is a direct consequence")
    log("    of the Hadamard product structure, driven by the nearest-neighbor")
    log("    1/(γ_n - γ_{n±1})² terms.")

# ── W-의존 추세 ──────────────────────────────────────────────────
if len(zeta_results) >= 3:
    log()
    log("  [W-의존 추세]")
    W_vals = sorted(zeta_results.keys())
    nn_pcts = [zeta_results[w]['nn_ratio_mean'] * 100 for w in W_vals]
    rho_rests = [zeta_results[w]['rho_rest'] for w in W_vals]
    for i, w in enumerate(W_vals):
        log(f"    W={w:4d}: NN%={nn_pcts[i]:.1f}%, ρ(rest,gap)={rho_rests[i]:.4f}")

    if nn_pcts[0] > nn_pcts[-1]:
        log("  → NN 기여율은 W 증가에 따라 감소 (원거리 항 축적)")
        log("  → 그러나 반상관 강도(ρ_full)는 W에 약의존 (C-295/C-296 확인)")

log()
elapsed = time.time() - t_start
log(f"총 소요: {elapsed:.1f}s")
log("=" * 70)

out_f.close()
print(f"\n결과 저장: {RESULT_PATH}")
