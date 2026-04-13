"""
중간점 곡률 비국소 분석 v2 — 결과 #18 후보 재실험
사이클 26 / 수학자 지시 2026-04-14 02:27

핵심 질문: 중간점 κ_total은 비국소(nearest-neighbor 이상의) 영점 정보를 인코딩하는가?

설계 근거:
  - 중간점 m=(γₙ+γₙ₊₁)/2에서 두 최근접 영점의 Hadamard L(s) 기여가 정확히 상쇄.
  - 따라서 κ_total = |ξ'/ξ|² 자체가 "비국소 기여만"을 측정하는 프로브.
  - κ_nn을 빼지 않아도 된다 — 자연이 이미 빼 주었음.

인프라 수정 (사이클 25 cap 82% 해결):
  - bundle_utils.connection_zeta의 h=10^{-20} 대신 mpmath.diff(xi_func, s) 사용.
  - mpmath.diff(): Richardson 외삽법, 자동 스텝 크기 → h 의존성 없음.
  - dps=100: t≤450에서 cap < 5% 목표.

실험 설계:
  1. zetazero(k), k=1..235 → 234쌍
  2. 각 중간점 m에서 κ = |mpmath.diff(xi_func, s) / xi_func(s)|², dps=100
  3. 상관 (a)~(d):
     (a) ρ(κ_mid, gap)         — 새니티: nn 상쇄 → |ρ|<0.3 예상
     (b) ρ(κ_mid, gap_next)    — 핵심: 비국소 정보 검증
     (c) ρ(κ_mid, density)     — 밀도 추적 보조
     (d) ρ(κ_mid, d_second)    — 2차 이웃 지배 여부
     (e) 기술 통계: cap 비율 (<5% 목표)

주의: bundle_utils.py 수정 금지. 로컬 κ 계산 함수 사용.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import mpmath
import numpy as np
from scipy import stats

# bundle_utils에서 xi_func만 import (curvature_zeta 사용 안 함)
from bundle_utils import xi_func

# ─── 정밀도 설정 ────────────────────────────────────────────────────────────
mpmath.mp.dps = 100  # dps=100: t≤450에서 cap 없어야 함. cap>10%이면 dps=150 권고.
CAP_VALUE = 1e18

OUTPUT_FILE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "results", "curvature_midpoint_nonlocal_v2.txt"
)

# ─── 로그 헬퍼 ──────────────────────────────────────────────────────────────
log_lines = []
def log(msg=""):
    print(msg, flush=True)
    log_lines.append(msg)

# ─── 로컬 κ 계산 함수 (bundle_utils.curvature_zeta 미사용) ──────────────────
def kappa_local_mpmath_diff(m_float):
    """
    중간점 m에서 κ = |ξ'(s)/ξ(s)|² 계산.
    mpmath.diff() 사용 → Richardson 외삽, h 의존성 없음.
    dps=100 상태에서 호출해야 함.

    반환: float 또는 None(cap/실패)
    """
    try:
        s = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(m_float)))
        xi_val = xi_func(s)
        # 영점 근처 판정 (중간점은 영점에서 gap/2 떨어져 있으므로 일반적으로 아님)
        if abs(xi_val) < mpmath.mpf(10) ** (-mpmath.mp.dps + 10):
            print(f"  WARNING: m={m_float:.4f}: xi near-zero ({float(abs(xi_val)):.2e})")
            return None
        # mpmath.diff(): 내부적으로 Richardson 외삽법으로 도함수 계산
        xi_deriv = mpmath.diff(xi_func, s)
        L = xi_deriv / xi_val
        kappa = float(abs(L) ** 2)
        if not (kappa == kappa) or kappa > CAP_VALUE:  # NaN or cap
            return None
        return kappa
    except Exception as e:
        print(f"  WARNING: m={m_float:.4f}: kappa 계산 실패 → {e}")
        return None


# ─── 1단계: 영점 수집 k=1..235 ──────────────────────────────────────────────
log("=" * 70)
log("중간점 곡률 비국소 분석 v2 (결과 #18 후보)")
log(f"dps={mpmath.mp.dps}, CAP={CAP_VALUE:.0e}, mpmath.diff() 사용")
log("=" * 70)
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
if N == 0:
    print("⚠️ 영점 0개 — 탐색 로직 점검 필요")
    sys.exit(1)

pairs = N - 1  # 234쌍
log(f"  연속 쌍 수: {pairs}")
log()


# ─── 2단계: 중간점 κ_total 계산 ──────────────────────────────────────────────
log("2단계: 중간점 κ 계산 (mpmath.diff, dps=100) ...")
log()

results = []
cap_count = 0

for i in range(pairs):
    g_n   = zeros[i]
    g_np1 = zeros[i + 1]
    m     = (g_n + g_np1) / 2.0
    gap   = g_np1 - g_n

    # gap_next: 마지막 쌍(i=233)은 None
    gap_next = zeros[i + 2] - zeros[i + 1] if i < pairs - 1 else None

    # 2차 이웃 거리: d_second = min(m - γₙ₋₁, γₙ₊₂ - m) (경계 처리)
    if i > 0 and i < pairs - 1:
        d_prev = m - zeros[i - 1]    # m - γₙ₋₁
        d_next = zeros[i + 2] - m    # γₙ₊₂ - m
        d_second = min(d_prev, d_next)
    elif i == 0:
        d_second = zeros[i + 2] - m  # 앞쪽 없음
    else:
        d_second = m - zeros[i - 1]  # 뒤쪽 없음

    # κ 계산
    kappa = kappa_local_mpmath_diff(m)

    if kappa is None:
        cap_count += 1
        # 통계를 위해 cap 기록
        continue

    results.append({
        'i': i,
        'm': m,
        'g_n': g_n,
        'g_np1': g_np1,
        'gap': gap,
        'gap_next': gap_next,
        'd_second': d_second,
        'kappa': kappa,
    })

    if (i + 1) % 30 == 0:
        log(f"  [{i+1}/{pairs}] m={m:.2f}, gap={gap:.4f}, κ={kappa:.3e}")

log()
n_valid = len(results)
cap_pct = cap_count / pairs * 100
log(f"  유효 데이터: {n_valid}개 / {pairs}쌍")
log(f"  cap 발생: {cap_count}개 ({cap_pct:.1f}%) [목표: <5%]")

if cap_pct > 10:
    print(f"⚠️ cap 비율 {cap_pct:.1f}% > 10% — dps=150 상향 권고")
elif cap_pct > 5:
    print(f"⚠️ cap 비율 {cap_pct:.1f}% > 5% — 주의")
else:
    log(f"  ✅ cap 비율 {cap_pct:.1f}% < 5% — 인프라 수정 성공")

if n_valid < 100:
    print(f"⚠️ 유효 데이터 {n_valid}개 < 100 — 분석 신뢰도 낮음")
    sys.exit(1)


# ─── 3단계: 상관 분석 ────────────────────────────────────────────────────────
log()
log("3단계: 상관 분석 ...")
log()

kappa_arr   = np.array([r['kappa']    for r in results])
gap_arr     = np.array([r['gap']      for r in results])
m_arr       = np.array([r['m']        for r in results])
d2_arr      = np.array([r['d_second'] for r in results])

# gap_next가 있는 쌍만 (핵심 분석용)
res_hn       = [r for r in results if r['gap_next'] is not None]
kappa_hn     = np.array([r['kappa']    for r in res_hn])
gap_next_arr = np.array([r['gap_next'] for r in res_hn])
gap_hn       = np.array([r['gap']      for r in res_hn])

# 국소 밀도: |γ - m| < 5 내 영점 수
zeros_np = np.array(zeros)
local_density = np.array([
    int(np.sum(np.abs(zeros_np - r['m']) < 5.0))
    for r in results
])

# ─── (a) 새니티: ρ(κ_mid, gap) — nn 상쇄 확인 ──────────────────────────────
log("─" * 70)
log("(a) 새니티: Spearman ρ(κ_midpoint, gap)")
log("    [nn 상쇄 → 약한 상관 예상 → |ρ| < 0.3이면 상쇄 확인]")
rho_a, p_a = stats.spearmanr(kappa_arr, gap_arr)
log(f"    ρ = {rho_a:.4f}, p = {p_a:.3e}  (n={len(gap_arr)})")
sanity_ok = abs(rho_a) < 0.3
log(f"    → 기준 |ρ| < 0.3: {'✅ nn 상쇄 확인' if sanity_ok else f'⚠️ |ρ|={abs(rho_a):.3f} ≥ 0.3 (상쇄 예상보다 강한 상관)'}")

# ─── (b) 핵심: ρ(κ_mid, gap_next) — 비국소 정보 ────────────────────────────
log()
log("(b) [핵심] Spearman ρ(κ_midpoint, gap_next)")
log("    [비국소 정보 검증: |ρ| > 0.3, p < 0.01 → 양성]")
rho_b, p_b = stats.spearmanr(kappa_hn, gap_next_arr)
log(f"    ρ = {rho_b:.4f}, p = {p_b:.3e}  (n={len(gap_next_arr)})")
key_pos = abs(rho_b) > 0.3 and p_b < 0.01
log(f"    → 기준 |ρ|>0.3 AND p<0.01: {'✅ 양성 (비국소 정보 존재)' if key_pos else '❌ 음성 (비국소 정보 없음)'}")

# ─── (c) 밀도: ρ(κ_mid, local_density) ─────────────────────────────────────
log()
log("(c) Spearman ρ(κ_midpoint, local_density)")
log("    [local_density = #{γ : |γ-m| < 5}]")
rho_c, p_c = stats.spearmanr(kappa_arr, local_density)
log(f"    ρ = {rho_c:.4f}, p = {p_c:.3e}  (n={len(kappa_arr)})")
aux_density = abs(rho_c) > 0.4
log(f"    → 기준 |ρ| > 0.4: {'✅ 보조 양성 (밀도 추적)' if aux_density else '음성'}")

# ─── (d) 2차 이웃: ρ(κ_mid, d_second) ──────────────────────────────────────
log()
log("(d) Spearman ρ(κ_midpoint, d_second)")
log("    [d_second = min(m-γₙ₋₁, γₙ₊₂-m) — 2차 이웃이 κ를 지배하는지]")
rho_d, p_d = stats.spearmanr(kappa_arr, d2_arr)
log(f"    ρ = {rho_d:.4f}, p = {p_d:.3e}  (n={len(kappa_arr)})")
aux_second = abs(rho_d) > 0.5
log(f"    → 기준 |ρ| > 0.5: {'✅ 2차 이웃 지배 (메커니즘 해명)' if aux_second else '해당 없음'}")

# ─── (e) 기술 통계 ──────────────────────────────────────────────────────────
log()
log("─" * 70)
log("(e) 기술 통계:")
log(f"  cap 비율: {cap_pct:.1f}%  {'✅' if cap_pct < 5 else ('⚠️' if cap_pct < 10 else '❌')}")
log(f"  κ_midpoint 범위: {kappa_arr.min():.3e} ~ {kappa_arr.max():.3e}")
log(f"  κ_midpoint mean={kappa_arr.mean():.3e}, median={np.median(kappa_arr):.3e}, std={kappa_arr.std():.3e}")
log(f"  gap 범위: {gap_arr.min():.4f} ~ {gap_arr.max():.4f}, mean={gap_arr.mean():.4f}")
log(f"  d_second 범위: {d2_arr.min():.4f} ~ {d2_arr.max():.4f}, mean={d2_arr.mean():.4f}")
log(f"  local_density: mean={local_density.mean():.2f}, range={local_density.min()}~{local_density.max()}")

# κ 분포: 상위 5% 극단값 확인
q95 = np.percentile(kappa_arr, 95)
q50 = np.percentile(kappa_arr, 50)
log(f"  κ P50={q50:.3e}, P95={q95:.3e}")

# cap 분포 확인: t 구간별 cap 비율
log()
log("  t 구간별 cap 현황:")
thresholds = [100, 200, 300, 400, 450]
# 각 영점 쌍의 중간점 t 범위
all_midpoints = [(zeros[i] + zeros[i+1])/2 for i in range(pairs)]
for idx_t, t_thresh in enumerate(thresholds):
    if idx_t == 0:
        t_lo = 0
    else:
        t_lo = thresholds[idx_t - 1]
    t_hi = t_thresh
    total_in = sum(1 for m2 in all_midpoints if t_lo <= m2 < t_hi)
    valid_in = sum(1 for r in results if t_lo <= r['m'] < t_hi)
    cap_in = total_in - valid_in
    if total_in > 0:
        log(f"    t∈[{t_lo},{t_hi}): 총{total_in}쌍, 유효{valid_in}, cap{cap_in} ({cap_in/total_in*100:.1f}%)")


# ─── 최종 판정 ──────────────────────────────────────────────────────────────
log()
log("=" * 70)
log("최종 판정")
log("=" * 70)
log()

infra_ok   = cap_pct < 5.0
sanity_ok2 = abs(rho_a) < 0.3  # nn 상쇄 확인

log(f"  [인프라] cap 비율 < 5%:                  {'✅ PASS' if infra_ok else '❌ FAIL'} ({cap_pct:.1f}%)")
log(f"  [실험유효성] |ρ(κ_mid, gap)| < 0.3:      {'✅ nn 상쇄 확인' if sanity_ok2 else '⚠️ 예상보다 강한 상관'} (|ρ|={abs(rho_a):.3f})")
log(f"  [핵심] |ρ(κ_mid, gap_next)| > 0.3 & p<0.01: {'✅ 양성' if key_pos else '❌ 음성'} (ρ={rho_b:.3f}, p={p_b:.3e})")
log(f"  [보조] |ρ(κ_mid, density)| > 0.4:        {'✅ 보조 양성' if aux_density else '음성'} (ρ={rho_c:.3f})")
log(f"  [메커니즘] |ρ(κ_mid, d_second)| > 0.5:   {'✅ 2차 이웃 지배' if aux_second else '해당 없음'} (ρ={rho_d:.3f})")
log()

if not infra_ok:
    verdict = "❌ INFRA FAIL — cap > 5%, dps=150으로 재시도 필요"
elif key_pos and aux_density:
    verdict = "★ 강한 양성 — κ_midpoint가 비국소 영점 구조(gap_next + 밀도)를 인코딩"
elif key_pos:
    verdict = "★ 양성 — κ_midpoint가 비국소 정보(gap_next)를 인코딩"
    if aux_second:
        verdict += " + 2차 이웃 지배 메커니즘 확인"
elif aux_density:
    verdict = "△ 보조 양성 — gap_next 기준 미달, 그러나 밀도 추적 확인"
else:
    verdict = "음성 — κ_midpoint는 gap_next/density와 무관, 비국소 정보 없음"

log(f"  최종 판정: {verdict}")
log()

# 결과 #18 등록 여부
if key_pos or aux_density:
    log("  → 결과 #18: 양성 등록 가능")
else:
    log("  → 결과 #18: 음성 등록 (κ_midpoint는 비국소 정보 없음)")

log()
log("  [수학자 판단 자료]")
log("  - 양성 서술: 'κ_midpoint (at nn 상쇄점)는 nearest-neighbor 이상의 비국소")
log("    영점 구조를 인코딩 → 프레임워크의 비자명한 예측력'")
log("  - 음성 서술: 'κ_midpoint는 gap_next/density와 무관 → κ ≈ local density의")
log("    단순 반영이며, nn 상쇄점에서 비국소 정보는 검출되지 않음'")
log("  - 결과 #16 (ρ=0.835)는 어느 결과든 '실용적 가치'로 남음")

# ─── 결과 저장 ──────────────────────────────────────────────────────────────
os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
    f.write('\n'.join(log_lines))

log()
log(f"결과 저장: {OUTPUT_FILE}")
log("완료.")
