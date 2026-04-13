"""
κ 잔차 분석 — 곡률의 비국소 영점 정보 검증 (결과 #18 후보)
사이클 25 / 수학자 지시 2026-04-14 01:43

핵심 질문: κ_residual = κ_total - κ_nn 이 비국소(non-local) 영점 정보를 인코딩하는가?

실험 설계:
  1. mpmath.zetazero(k) k=1..235 (t∈[0,450])
  2. 중간점 m = (γₙ + γₙ₊₁)/2에서:
     - κ_total = |ξ'/ξ|² (bundle_utils.curvature_zeta)
     - κ_nn = 1/(m-γₙ)² + 1/(m-γₙ₊₁)² (Hadamard 최근접 2항)
     - κ_residual = κ_total - κ_nn
  3. 상관 분석 (a)~(e)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import mpmath
import numpy as np
from scipy import stats

from bundle_utils import curvature_zeta, xi_func

# ─── 정밀도 설정 ───────────────────────────────────────────────────────────
mpmath.mp.dps = 50  # t < 450이므로 50 충분

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           "results", "curvature_residual_nonlocal.txt")

# ─── 로그 출력 헬퍼 ────────────────────────────────────────────────────────
log_lines = []
def log(msg=""):
    print(msg, flush=True)
    log_lines.append(msg)

# ─── 1단계: 영점 수집 k=1..235 ─────────────────────────────────────────────
log("=" * 65)
log("κ 잔차 분석 — 비국소 영점 정보 검증")
log("=" * 65)
log()
log("1단계: 영점 수집 (k=1..235, t∈[0,450]) ...")

zeros = []
for k in range(1, 236):
    try:
        z = mpmath.zetazero(k)
        t = float(z.imag)
        zeros.append(t)
    except Exception as e:
        print(f"WARNING: k={k} 영점 계산 실패: {e}")

zeros = sorted(zeros)
N = len(zeros)
log(f"  수집된 영점 수: {N}")
log(f"  t 범위: {zeros[0]:.3f} ~ {zeros[-1]:.3f}")

# 연속 쌍: 234개 (인덱스 0..233)
# gap[i] = zeros[i+1] - zeros[i]
# gap_next[i] = zeros[i+2] - zeros[i+1]
# 중간점 m[i] = (zeros[i] + zeros[i+1]) / 2

pairs = N - 1  # 234쌍
log(f"  연속 쌍 수: {pairs}")
log()

# ─── 2단계: 중간점에서 κ_total, κ_nn, κ_residual 계산 ──────────────────────
log("2단계: 각 중간점에서 κ 계산 ...")
log(f"  (dps={mpmath.mp.dps})")
log()

results = []
cap_count = 0
nan_count = 0

for i in range(pairs):
    g_n   = zeros[i]
    g_np1 = zeros[i+1]
    m     = (g_n + g_np1) / 2.0
    gap   = g_np1 - g_n

    # gap_next: 마지막 쌍(i=233)은 None
    if i < pairs - 1:
        gap_next = zeros[i+2] - zeros[i+1]
    else:
        gap_next = None

    # κ_total via bundle_utils.curvature_zeta
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(m))
    try:
        kappa_total = curvature_zeta(s)
    except Exception as e:
        print(f"WARNING: i={i}, m={m:.4f}: curvature_zeta 실패 → {e}")
        nan_count += 1
        continue

    if not np.isfinite(kappa_total) or kappa_total > 1e18:
        # 중간점은 영점에서 gap/2 이상 떨어져 있어 cap이 발생해선 안 됨
        # 발생 시 경고
        print(f"WARNING: i={i}, m={m:.4f}: κ_total={kappa_total:.3e} — 비정상 (gap={gap:.4f})")
        cap_count += 1
        continue

    # κ_nn = 1/(m-γₙ)² + 1/(m-γₙ₊₁)²
    d_n   = m - g_n    # = gap/2
    d_np1 = g_np1 - m  # = gap/2
    kappa_nn = 1.0 / (d_n**2) + 1.0 / (d_np1**2)

    kappa_residual = kappa_total - kappa_nn
    nn_ratio = kappa_nn / kappa_total if kappa_total > 0 else float('nan')

    results.append({
        'i': i,
        'm': m,
        'g_n': g_n,
        'g_np1': g_np1,
        'gap': gap,
        'gap_next': gap_next,
        'kappa_total': kappa_total,
        'kappa_nn': kappa_nn,
        'kappa_residual': kappa_residual,
        'nn_ratio': nn_ratio,
    })

    if (i + 1) % 50 == 0:
        log(f"  {i+1}/{pairs} 완료 ... (최근: m={m:.2f}, κ_total={kappa_total:.3e}, "
            f"κ_nn={kappa_nn:.3e}, ratio={nn_ratio:.3f})")

log()
log(f"  유효 데이터: {len(results)}개 / {pairs}쌍")
log(f"  cap 발생: {cap_count}개, NaN/실패: {nan_count}개")

if len(results) < 100:
    print("⚠️ 유효 데이터 100개 미만 — 분석 신뢰도 낮음. 스크립트 점검 필요.")

# ─── 데이터 배열 준비 ────────────────────────────────────────────────────────
log()
log("3단계: 상관 분석 ...")

# gap_next가 없는 마지막 쌍 제외 (상관 (c)용)
res_full = results
res_hasnext = [r for r in results if r['gap_next'] is not None]

kappa_total_arr   = np.array([r['kappa_total']   for r in res_full])
kappa_nn_arr      = np.array([r['kappa_nn']       for r in res_full])
kappa_res_arr     = np.array([r['kappa_residual'] for r in res_full])
nn_ratio_arr      = np.array([r['nn_ratio']       for r in res_full])
gap_arr           = np.array([r['gap']            for r in res_full])
m_arr             = np.array([r['m']              for r in res_full])

kappa_res_hn      = np.array([r['kappa_residual'] for r in res_hasnext])
gap_arr_hn        = np.array([r['gap']            for r in res_hasnext])
gap_next_arr      = np.array([r['gap_next']       for r in res_hasnext])

# 국소 밀도: |γ - m| < 5 내 영점 수
zeros_np = np.array(zeros)
local_density = np.array([
    int(np.sum(np.abs(zeros_np - r['m']) < 5.0))
    for r in res_full
])

log()
log("─" * 65)
log("(a) Control: Spearman ρ(κ_total, gap)")
rho_a, p_a = stats.spearmanr(kappa_total_arr, gap_arr)
log(f"    ρ = {rho_a:.4f}, p = {p_a:.3e}  (n={len(gap_arr)})")
log(f"    → 기준 |ρ| > 0.8: {'✅ PASS' if abs(rho_a) > 0.8 else '❌ FAIL'}")
if abs(rho_a) <= 0.8:
    log("    ⚠️ Control 실패 → 스크립트 버그 가능성 점검 필요")

log()
log("(b) Spearman ρ(κ_residual, gap)")
rho_b, p_b = stats.spearmanr(kappa_res_arr, gap_arr)
log(f"    ρ = {rho_b:.4f}, p = {p_b:.3e}  (n={len(gap_arr)})")
log(f"    → 해석: 잔차가 현재 간격과 얼마나 상관?")

log()
log("(c) [핵심] Spearman ρ(κ_residual, gap_next) — 비국소 정보 검증")
rho_c, p_c = stats.spearmanr(kappa_res_hn, gap_next_arr)
log(f"    ρ = {rho_c:.4f}, p = {p_c:.3e}  (n={len(gap_next_arr)})")
log(f"    → 기준 |ρ| > 0.3 AND p < 0.01: "
    f"{'✅ 양성 (비국소 정보 존재)' if abs(rho_c) > 0.3 and p_c < 0.01 else '❌ 음성 (nn 지배, 비국소 없음)'}")

log()
log("(d) Spearman ρ(κ_residual, local_density)")
rho_d, p_d = stats.spearmanr(kappa_res_arr, local_density)
log(f"    ρ = {rho_d:.4f}, p = {p_d:.3e}  (n={len(gap_arr)})")
log(f"    → 기준 |ρ| > 0.4: {'✅ 보조 양성 (밀도 추적)' if abs(rho_d) > 0.4 else '음성'}")

log()
log("(e) κ_nn / κ_total 지배도 분석")
# 이상치 있는 경우 제외 (비율이 음수이거나 >2인 경우)
valid_ratio = nn_ratio_arr[np.isfinite(nn_ratio_arr) & (nn_ratio_arr > 0)]
log(f"    n(유효 비율) = {len(valid_ratio)}")
log(f"    평균 κ_nn/κ_total = {np.mean(valid_ratio):.4f}")
log(f"    표준편차         = {np.std(valid_ratio):.4f}")
log(f"    중앙값           = {np.median(valid_ratio):.4f}")
log(f"    최솟값/최댓값    = {np.min(valid_ratio):.4f} / {np.max(valid_ratio):.4f}")
pct_nn_dom = np.mean(valid_ratio >= 0.95) * 100
pct_nonlocal_dom = np.mean(valid_ratio < 0.8) * 100
log(f"    nn_ratio >= 0.95 (nn 지배): {pct_nn_dom:.1f}%")
log(f"    nn_ratio <  0.80 (비국소 지배): {pct_nonlocal_dom:.1f}%")

mean_ratio = np.mean(valid_ratio)
log(f"    → 기준 평균 > 0.95: {'✅ nn 지배 확인' if mean_ratio > 0.95 else '❌'}")

# ─── 비국소 지배 점 분석 ────────────────────────────────────────────────────
log()
log("─" * 65)
log("4단계: 비국소 지배 점 특성 분석 (nn_ratio < 0.8)")

nonlocal_mask = np.array([
    (r['nn_ratio'] < 0.8 and np.isfinite(r['nn_ratio']) and r['nn_ratio'] > 0)
    for r in res_full
])
n_nonlocal = np.sum(nonlocal_mask)
log(f"  비국소 지배 점 수: {n_nonlocal}")
if n_nonlocal > 0:
    nl_gaps = gap_arr[nonlocal_mask]
    nl_m    = m_arr[nonlocal_mask]
    nl_res  = kappa_res_arr[nonlocal_mask]
    nl_ratio = nn_ratio_arr[nonlocal_mask]
    log(f"  t 범위: {nl_m.min():.1f} ~ {nl_m.max():.1f}")
    log(f"  평균 간격: {nl_gaps.mean():.4f} vs 전체 평균 {gap_arr.mean():.4f}")
    log(f"  평균 nn_ratio: {nl_ratio.mean():.4f}")
    log(f"  평균 κ_residual: {nl_res.mean():.3e}")

# ─── 기술 통계 ──────────────────────────────────────────────────────────────
log()
log("─" * 65)
log("기술 통계:")
log(f"  gap 범위: {gap_arr.min():.4f} ~ {gap_arr.max():.4f}, 평균 {gap_arr.mean():.4f}")
log(f"  κ_total: mean={kappa_total_arr.mean():.3e}, median={np.median(kappa_total_arr):.3e}")
log(f"  κ_nn:    mean={kappa_nn_arr.mean():.3e}, median={np.median(kappa_nn_arr):.3e}")
log(f"  κ_residual: mean={kappa_res_arr.mean():.3e}, median={np.median(kappa_res_arr):.3e}")
log(f"  κ_residual < 0 비율: {np.mean(kappa_res_arr < 0)*100:.1f}%  (Hadamard 고차항 영향)")

# ─── 최종 판정 ──────────────────────────────────────────────────────────────
log()
log("=" * 65)
log("최종 판정")
log("=" * 65)

ctrl_pass  = abs(rho_a) > 0.8
key_pass   = abs(rho_c) > 0.3 and p_c < 0.01
aux_pass   = abs(rho_d) > 0.4
nn_dom     = mean_ratio > 0.95

log(f"  (a) Control |ρ(κ_total, gap)| > 0.8:          {'✅ PASS' if ctrl_pass else '❌ FAIL'} (ρ={rho_a:.3f})")
log(f"  (c) [핵심] |ρ(κ_res, gap_next)| > 0.3 & p<0.01: {'✅ PASS' if key_pass else '❌ FAIL'} (ρ={rho_c:.3f}, p={p_c:.3e})")
log(f"  (d) |ρ(κ_res, density)| > 0.4:                {'✅ PASS' if aux_pass else '❌ FAIL'} (ρ={rho_d:.3f})")
log(f"  (e) 평균 nn_ratio > 0.95:                     {'✅ nn지배' if nn_dom else '❌'} ({mean_ratio:.3f})")
log()

if not ctrl_pass:
    verdict = "⚠️ INVALID — Control 실패 (스크립트 버그 점검 필요)"
elif key_pass:
    verdict = "★ 양성 — κ residual이 비국소(다음 간격) 정보를 인코딩함"
    if aux_pass:
        verdict += " + 밀도 추적 보조 확인"
else:
    if nn_dom:
        verdict = "음성 — κ ≈ nearest-neighbor Hadamard (nn 지배, 비국소 정보 없음)"
    else:
        verdict = "음성/불명 — nn_ratio < 0.95이나 gap_next 상관 없음"

log(f"  판정: {verdict}")
log()
log(f"  [결과 #18 후보] {'양성 등록' if key_pass and ctrl_pass else '음성 등록'}")

# ─── 파일 저장 ──────────────────────────────────────────────────────────────
with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
    f.write('\n'.join(log_lines))

log()
log(f"결과 저장: {OUTPUT_FILE}")
log("완료.")
