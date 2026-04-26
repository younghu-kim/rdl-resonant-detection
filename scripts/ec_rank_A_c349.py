#!/usr/bin/env python3
"""
=============================================================================
[사이클 #349] EC rank vs A(t₀) — Paper B 핵심 실험
=============================================================================

목적:
  타원곡선 L-함수에서 ξ-다발 곡률 보정항 A(t₀)가 BSD rank에 의존하는지 검증.
  양성이면: "ξ-다발 곡률이 산술적 rank를 감지" — RH↔BSD 교차점.
  음성이면: A(t₀)는 순수 해석적 양 (rank 무관).

방법:
  1. rank 0, 1, 2, 3 각 5곡선 (총 20곡선)
  2. PARI lfunzeros로 영점 탐색, Cauchy 윤곽으로 c₀, c₁ 추출
  3. A(t₀) = Im(c₀)² + 2·Re(c₁)
  4. rank별 A(t₀) 분포 비교: Kruskal-Wallis + 쌍별 Mann-Whitney

체크리스트:
  [x] Cauchy 적분: lfunlambda (GL(2) 완성 함수)
  [x] 영점: lfunzeros → s=center+it (center = k/2 = 1)
  [x] python -u (버퍼링 방지)
  [x] NaN/Inf 체크
  [x] 결과 파일에 설정·곡선별 결과·통계·판정 포함
  [x] ε=-1 곡선: Im(Λ) 기반 부호변환 → 영점 찾기 OK (PARI 내부 처리)
  [x] rank≥2: center 영점(t=0)은 제외 → A(t₀) 측정에 부적합 (다중 영점)

결과: results/ec_rank_A_c349.txt
=============================================================================
"""

import sys, os, time, math
import numpy as np
from scipy import stats

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')

import cypari2
pari = cypari2.Pari()
pari.set_real_precision(50)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

T_MIN = 0.5       # 영점 탐색 하한 (t=0 근방 제외)
T_MAX = 60.0      # 영점 탐색 상한
CAUCHY_R = 0.01   # Cauchy 적분 반지름
N_PTS = 64        # Cauchy 적분 점 수
H_DEC = '0.00001' # 수치 미분 간격

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/ec_rank_A_c349.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 20곡선 정의 — LMFDB 표준 라벨
# rank 0: ε=+1, L(E,1) ≠ 0
# rank 1: ε=-1, L(E,1) = 0 (단순)
# rank 2: ε=+1, L(E,1) = 0 (이중)
# rank 3: ε=-1, L(E,1) = 0 (삼중)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CURVES = [
    # --- Rank 0 (PARI 검증 완료) ---
    {'name': '11a1',   'coeffs': '[0,-1,1,-10,-20]',   'N': 11,   'rank': 0},
    {'name': '14a1',   'coeffs': '[1,0,1,4,-6]',       'N': 14,   'rank': 0},
    {'name': '15a1',   'coeffs': '[1,1,1,-10,-10]',    'N': 15,   'rank': 0},
    {'name': '19a1',   'coeffs': '[0,1,1,-9,-15]',     'N': 19,   'rank': 0},
    {'name': '20a1',   'coeffs': '[0,1,0,4,4]',        'N': 20,   'rank': 0},
    # --- Rank 1 (PARI 검증 완료) ---
    {'name': '37a1',   'coeffs': '[0,0,1,-1,0]',       'N': 37,   'rank': 1},
    {'name': '43a1',   'coeffs': '[0,1,1,0,0]',        'N': 43,   'rank': 1},
    {'name': '53a1',   'coeffs': '[1,-1,1,0,0]',       'N': 53,   'rank': 1},
    {'name': '57a',    'coeffs': '[0,-1,1,-2,2]',      'N': 57,   'rank': 1},
    {'name': '58a1',   'coeffs': '[1,-1,0,-1,1]',      'N': 58,   'rank': 1},
    # --- Rank 2 (PARI ellrank=[2,2] 확인) ---
    {'name': 'r2_389', 'coeffs': '[0,1,1,-2,0]',       'N': 389,  'rank': 2},
    {'name': 'r2_571', 'coeffs': '[0,-2,1,-3,6]',      'N': 571,  'rank': 2},
    {'name': 'r2_664', 'coeffs': '[0,0,0,-7,10]',      'N': 664,  'rank': 2},
    {'name': 'r2_681', 'coeffs': '[0,-1,1,0,2]',       'N': 681,  'rank': 2},
    {'name': 'r2_709', 'coeffs': '[0,-1,1,-2,0]',      'N': 709,  'rank': 2},
    # --- Rank 3 (PARI ellrank=[3,3] 확인) ---
    {'name': 'r3_5077',  'coeffs': '[0,0,1,-7,6]',     'N': 5077,  'rank': 3},
    {'name': 'r3_12279', 'coeffs': '[0,2,1,-9,2]',     'N': 12279, 'rank': 3},
    {'name': 'r3_16811', 'coeffs': '[0,-3,1,2,6]',     'N': 16811, 'rank': 3},
]

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 유틸
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

LOG_LINES = []

def log(msg):
    print(msg, flush=True)
    LOG_LINES.append(msg)

def save_results():
    with open(RESULT_PATH, 'w') as f:
        f.write('\n'.join(LOG_LINES))
    print(f"\n[저장] {RESULT_PATH}", flush=True)

def pari_to_float(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    if s in ('None', ''):
        return float('nan')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeros(coeffs_str, t_min, t_max, lf_var):
    """PARI lfunzeros로 영점 목록 반환."""
    try:
        pari(f'__zz = lfunzeros({lf_var}, [{t_min}, {t_max}])')
        n = int(pari_to_float(pari(f'#__zz')))
        zeros = []
        for i in range(1, n + 1):
            t = pari_to_float(pari(f'__zz[{i}]'))
            if not math.isnan(t) and t > 0.3:
                zeros.append(t)
        return sorted(zeros)
    except Exception as e:
        log(f"  ⚠️ lfunzeros 실패: {e}")
        return []


def compute_A_cauchy(t_zero, lf_var, r=CAUCHY_R, n_pts=N_PTS, h_dec=H_DEC):
    """
    Cauchy 윤곽으로 Laurent 계수 c₀, c₁ → A(γ) = Im(c₀)² + 2Re(c₁).
    center = 1 (weight 2 EC).
    """
    script = (
        f'rho=1+{t_zero:.15f}*I; c0s=0+0.*I; c1s=0+0.*I; cnt=0; '
        f'for(k=0,{n_pts-1}, '
        f'th=2*Pi*k/{n_pts}; uu={r}*exp(th*I); ss=rho+uu; '
        f'Lv=lfunlambda({lf_var},ss); '
        f'if(abs(Lv)<1e-20,next); '
        f'Lp=lfunlambda({lf_var},ss+{h_dec}); '
        f'Lm=lfunlambda({lf_var},ss-{h_dec}); '
        f'LpL=(Lp-Lm)/(2*{h_dec}*Lv); '
        f'gg=LpL-1/uu; c0s=c0s+gg; c1s=c1s+gg/uu; cnt=cnt+1); '
        f'[c0s/cnt, c1s/cnt, cnt]'
    )
    try:
        res = pari(script)
        c0_re = pari_to_float(pari(f'real(({res})[1])'))
        c0_im = pari_to_float(pari(f'imag(({res})[1])'))
        c1_re = pari_to_float(pari(f'real(({res})[2])'))
        n_cnt = pari_to_float(pari(f'({res})[3]'))
        if any(math.isnan(v) for v in [c0_re, c0_im, c1_re]):
            return None
        if n_cnt < n_pts // 2:
            return None
        A = c0_im**2 + 2 * c1_re
        if math.isnan(A) or math.isinf(A):
            return None
        return {'A': A, 'Re_c0': c0_re, 'Im_c0': c0_im, 'Re_c1': c1_re}
    except Exception:
        return None


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 78)
log("[C-349] EC rank vs A(t₀) — Paper B 핵심 실험")
log("=" * 78)
log(f"  곡선: {len(CURVES)}개 (rank 0×5, rank 1×5, rank 2×5, rank 3×3)")
log(f"  T=[{T_MIN},{T_MAX}], Cauchy_R={CAUCHY_R}, N_PTS={N_PTS}")
log(f"  중심: s = 1 + it (weight 2, center = k/2 = 1)")
log("")

# Step 0: 곡선 초기화 + rank 검증
log("[Step 0] PARI L-함수 초기화 + rank 검증...")
t0_global = time.time()

valid_curves = []
for curve in CURVES:
    name = curve['name']
    coeffs = curve['coeffs']
    claimed_rank = curve['rank']
    try:
        pari(f'E_{name} = ellinit({coeffs})')
        actual_rank = int(pari_to_float(pari(f'ellrank(E_{name})[1]')))
        pari(f'L_{name} = lfuncreate(E_{name})')
        fe_check = pari_to_float(pari(f'lfuncheckfeq(L_{name})'))
        log(f"  {name:12s} N={curve['N']:5d} rank={claimed_rank} "
            f"ellrank≥{actual_rank} FE={fe_check:.0f}자리 ✅")
        if actual_rank < claimed_rank:
            log(f"    ⚠️ ellrank 하한({actual_rank}) < claimed({claimed_rank}). 진행하되 주의.")
        curve['lf_var'] = f'L_{name}'
        curve['fe_digits'] = fe_check
        valid_curves.append(curve)
    except Exception as e:
        log(f"  {name:12s} ❌ 실패: {e}")

log(f"\n  유효 곡선: {len(valid_curves)}/{len(CURVES)}")
log(f"  초기화 소요: {time.time() - t0_global:.1f}s")
log("")

# Step 1: 영점 탐색 + A(t₀) 계산
log("=" * 78)
log("[Step 1] 영점 탐색 + A(t₀) 계산")
log("=" * 78)

all_results = {}  # name → {'rank', 'N', 'n_zeros', 'A_values', 'A_mean', 'A_std'}

for curve in valid_curves:
    name = curve['name']
    lf = curve['lf_var']
    rank = curve['rank']
    N_cond = curve['N']

    log(f"\n--- {name} (rank {rank}, N={N_cond}) ---")
    t0_curve = time.time()

    # 영점 탐색
    zeros = get_zeros(curve['coeffs'], T_MIN, T_MAX, lf)
    if len(zeros) < 3:
        log(f"  영점 {len(zeros)}개 — 부족, 스킵")
        all_results[name] = {'rank': rank, 'N': N_cond, 'n_zeros': len(zeros),
                             'A_values': [], 'A_mean': float('nan'), 'A_std': float('nan')}
        continue

    log(f"  영점 {len(zeros)}개  t∈[{zeros[0]:.3f}, {zeros[-1]:.3f}]")

    # A(t₀) 계산
    A_vals = []
    Re_c0_list = []
    n_fail = 0

    for i, tz in enumerate(zeros):
        result = compute_A_cauchy(tz, lf)
        if result is None:
            n_fail += 1
            continue
        A_vals.append(result['A'])
        Re_c0_list.append(result['Re_c0'])

    elapsed = time.time() - t0_curve

    if len(A_vals) == 0:
        log(f"  ⚠️ A 계산 전부 실패 ({n_fail}건)")
        all_results[name] = {'rank': rank, 'N': N_cond, 'n_zeros': len(zeros),
                             'A_values': [], 'A_mean': float('nan'), 'A_std': float('nan')}
        continue

    A_arr = np.array(A_vals)
    A_mean = np.mean(A_arr)
    A_std = np.std(A_arr)
    A_med = np.median(A_arr)
    Re_c0_max = max(abs(x) for x in Re_c0_list) if Re_c0_list else float('nan')

    log(f"  A(t₀): mean={A_mean:.4f} ± {A_std:.4f}, median={A_med:.4f}")
    log(f"  |Re(c₀)|_max={Re_c0_max:.2e} (slope=2 검증)")
    log(f"  유효={len(A_vals)}, 실패={n_fail}, 소요={elapsed:.1f}s")

    all_results[name] = {
        'rank': rank, 'N': N_cond, 'n_zeros': len(zeros),
        'A_values': A_vals, 'A_mean': A_mean, 'A_std': A_std,
        'A_median': A_med, 'Re_c0_max': Re_c0_max,
    }

log("")

# Step 2: rank별 통계 비교
log("=" * 78)
log("[Step 2] rank별 A(t₀) 비교")
log("=" * 78)

rank_A = {0: [], 1: [], 2: [], 3: []}  # rank → [모든 A 값]
rank_means = {0: [], 1: [], 2: [], 3: []}  # rank → [곡선별 mean A]

for name, res in all_results.items():
    r = res['rank']
    if res['A_values']:
        rank_A[r].extend(res['A_values'])
        rank_means[r].append(res['A_mean'])

log("\n  [개별 곡선 요약]")
log(f"  {'곡선':12s} {'rank':>4s} {'N':>6s} {'n':>4s} {'A_mean':>10s} {'A_std':>10s} {'A_med':>10s}")
log(f"  {'-'*62}")
for name, res in sorted(all_results.items(), key=lambda x: (x[1]['rank'], x[1]['N'])):
    if res['A_values']:
        log(f"  {name:12s} {res['rank']:4d} {res['N']:6d} {len(res['A_values']):4d} "
            f"{res['A_mean']:10.4f} {res['A_std']:10.4f} {res.get('A_median', float('nan')):10.4f}")
    else:
        log(f"  {name:12s} {res['rank']:4d} {res['N']:6d} {res['n_zeros']:4d}   (실패)")

log(f"\n  [rank별 통합]")
log(f"  {'rank':>4s} {'n_curves':>8s} {'n_zeros':>8s} {'A_mean':>10s} {'A_std':>10s} {'A_med':>10s}")
log(f"  {'-'*52}")
for r in [0, 1, 2, 3]:
    if rank_A[r]:
        arr = np.array(rank_A[r])
        log(f"  {r:4d} {len(rank_means[r]):8d} {len(arr):8d} "
            f"{np.mean(arr):10.4f} {np.std(arr):10.4f} {np.median(arr):10.4f}")
    else:
        log(f"  {r:4d}        0        0       N/A       N/A       N/A")

# Step 3: 통계 검정
log(f"\n  [통계 검정]")

# Kruskal-Wallis (전체)
groups = [np.array(rank_A[r]) for r in [0, 1, 2, 3] if len(rank_A[r]) > 0]
if len(groups) >= 2:
    H_stat, H_p = stats.kruskal(*groups)
    log(f"  Kruskal-Wallis: H={H_stat:.4f}, p={H_p:.4e}")
    if H_p < 0.05:
        log(f"  → ★ rank간 A(t₀) 분포 유의차 있음 (p<0.05)")
    else:
        log(f"  → rank간 A(t₀) 분포 유의차 없음 (p≥0.05)")
else:
    log(f"  Kruskal-Wallis: 비교 가능 그룹 부족")
    H_p = float('nan')

# 쌍별 Mann-Whitney
log(f"\n  [쌍별 Mann-Whitney U]")
log(f"  {'rank_i':>6s} {'rank_j':>6s} {'U':>10s} {'p':>12s} {'A_i-A_j':>10s}")
log(f"  {'-'*48}")
for i in range(4):
    for j in range(i + 1, 4):
        if len(rank_A[i]) >= 3 and len(rank_A[j]) >= 3:
            U, p = stats.mannwhitneyu(rank_A[i], rank_A[j], alternative='two-sided')
            diff = np.mean(rank_A[i]) - np.mean(rank_A[j])
            sig = "★" if p < 0.05 else " "
            log(f"  {i:6d} {j:6d} {U:10.1f} {p:12.4e} {diff:10.4f} {sig}")

# Spearman: rank vs A_mean (곡선 단위)
all_ranks_flat = []
all_A_flat = []
for name, res in all_results.items():
    if res['A_values']:
        all_ranks_flat.append(res['rank'])
        all_A_flat.append(res['A_mean'])

if len(all_ranks_flat) >= 5:
    rho_s, p_s = stats.spearmanr(all_ranks_flat, all_A_flat)
    log(f"\n  Spearman ρ(rank, A_mean) = {rho_s:.4f} (p={p_s:.4e})")
    if p_s < 0.05:
        log(f"  → ★ rank와 A_mean 유의한 상관")
    else:
        log(f"  → rank와 A_mean 유의한 상관 없음")

# Spearman: conductor vs A_mean (교란 변수)
all_N_flat = [all_results[n]['N'] for n in all_results if all_results[n]['A_values']]
if len(all_N_flat) >= 5:
    rho_N, p_N = stats.spearmanr(all_N_flat, all_A_flat)
    log(f"  Spearman ρ(N, A_mean) = {rho_N:.4f} (p={p_N:.4e})")
    log(f"  → {'⚠️ conductor 교란 의심' if p_N < 0.05 else 'conductor 교란 없음'}")

# 편상관: rank vs A_mean | N (conductor 통제)
if len(all_ranks_flat) >= 7:
    from scipy.stats import spearmanr
    # rank vs A | N
    rho_rA = spearmanr(all_ranks_flat, all_A_flat)[0]
    rho_rN = spearmanr(all_ranks_flat, all_N_flat)[0]
    rho_AN = spearmanr(all_A_flat, all_N_flat)[0]
    # partial correlation
    numer = rho_rA - rho_rN * rho_AN
    denom = math.sqrt((1 - rho_rN**2) * (1 - rho_AN**2))
    if denom > 1e-10:
        rho_partial = numer / denom
        log(f"  편상관 ρ(rank, A | N) = {rho_partial:.4f}")

# Step 4: 최종 판정
log("")
log("=" * 78)
log("[C-349 최종 판정]")
log("=" * 78)

if not math.isnan(H_p):
    if H_p < 0.01:
        log("  ★★★★ 강양성: rank간 A(t₀) 분포 유의차 (p<0.01)")
        log("  → ξ-다발 곡률이 BSD rank를 감지. Paper B 핵심 결과.")
    elif H_p < 0.05:
        log("  ★★★ 양성: rank간 A(t₀) 분포 유의차 (p<0.05)")
        log("  → 추가 곡선으로 확인 필요.")
    elif H_p < 0.10:
        log("  ★★ 조건부 양성: 약한 신호 (p<0.10)")
        log("  → 곡선 수 증가 또는 T 범위 확대 필요.")
    else:
        log("  ★ 음성: rank간 A(t₀) 유의차 없음 (p≥0.10)")
        log("  → A(t₀)는 순수 해석적 양. rank 무관.")
        log("  → 음성이어도 가치: '다발 곡률은 산술을 보지 못한다' 확인.")

log(f"\n  총 소요: {time.time() - t0_global:.1f}s")
save_results()
