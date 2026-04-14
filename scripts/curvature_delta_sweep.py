#!/usr/bin/env python3
"""
결과 #35: δ-sweep — κ 차선도항의 δ-의존성 검증
=================================================
목적: 사이클 49 수학자 지시
  - #34의 99개 영점 그대로 사용
  - δ = {0.01, 0.02, 0.03, 0.05, 0.1} (5개)
  - 각 δ에서: κ(δ) = |ξ'/ξ(1/2+δ+it)|², Δκ(δ) = κ-1/δ²
  - NNS_min은 δ-독립 (영점 간격 고정)

파트 A: 영점 수집 + NNS 계산 (1회)
파트 B: δ별 1변수 회귀 (log(t/2π)) → a(δ), a(δ)·δ 상수성
파트 C: δ별 2변수 회귀 (Lorentzian²) → R²_2var(δ) > 0.88?
파트 D: 통합 회귀 (495 = 99×5 데이터) → R² > 0.90?
파트 E: δ=0.01 수치 정밀도 진단

이론:
  f₁ = δ²/(δ²+NNS_min²)²  (Lorentzian²)
  a(δ)·δ ≈ 상수 ≈ 0.048  (이론 예측)
  통합 모델: Δκ ~ (a/δ)·log(t/2π) + b/δ² + c·δ²/(δ²+NNS²)²
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats as scipy_stats

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.path.join(BASE_DIR, '..', 'results', 'curvature_delta_sweep.txt')
os.makedirs(os.path.join(BASE_DIR, '..', 'results'), exist_ok=True)

_log_buf = []
def log(msg=''):
    print(msg, flush=True)
    _log_buf.append(str(msg))

def save():
    with open(RESULTS_PATH, 'w') as f:
        f.write('\n'.join(_log_buf))

# ═══════════════════════════════════════════════════════════════════════════
# 상수
# ═══════════════════════════════════════════════════════════════════════════

DELTAS = [0.01, 0.02, 0.03, 0.05, 0.1]

# #34와 동일한 99개 영점 인덱스
RAW_INDICES = np.linspace(1, 20000, 100).astype(int)
ZERO_INDICES = [int(n) for n in RAW_INDICES if n >= 5]

print(f"총 영점 수: {len(ZERO_INDICES)}개  (5~20000)", flush=True)
print(f"δ 값: {DELTAS}", flush=True)
print(f"총 계산: {len(ZERO_INDICES) * len(DELTAS)}회", flush=True)


def get_dps(t, delta=None):
    """정밀도 설정 — δ=0.01은 최소 80"""
    if t > 30000: base = 120
    elif t > 10000: base = 100
    elif t > 5000:  base = 80
    else: base = max(50, int(30 + t / 200))
    # δ=0.01: 수치 안정성 위해 최소 80
    if delta is not None and delta <= 0.01:
        base = max(base, 80)
    return base


# ═══════════════════════════════════════════════════════════════════════════
# 해석적 계산 함수 (bundle_utils 함정 회피 — 직접 구현)
# ═══════════════════════════════════════════════════════════════════════════

def G_component(s):
    """
    G(s) = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2
    ξ'/ξ(s)의 감마+극 기여 (복소수 반환)
    """
    return (mpmath.mpf(1)/s
            + mpmath.mpf(1)/(s - 1)
            - mpmath.log(mpmath.pi)/2
            + mpmath.digamma(s/2)/2)


def zeta_log_deriv(s):
    """
    ζ'(s)/ζ(s): h=1e-6 중앙 차분 (mpmath.diff 금지)
    """
    z = mpmath.zeta(s)
    if abs(z) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e8, 0)
    h = mpmath.mpf('1e-6')
    zd = (mpmath.zeta(s + h) - mpmath.zeta(s - h)) / (2 * h)
    return zd / z


def xi_log_deriv(s):
    """ξ'/ξ(s) = G(s) + ζ'/ζ(s) — 해석적 공식"""
    return G_component(s) + zeta_log_deriv(s)


def get_zetazero_t(n, dps=50):
    """n번째 리만 제타 영점의 허수부 반환"""
    with mpmath.workdps(dps):
        z = mpmath.zetazero(n)
        return float(mpmath.im(z))


def compute_kappa(t, delta):
    """
    κ = |ξ'/ξ(s)|², Δκ = κ - 1/δ²
    반환: (kappa, delta_kappa)
    """
    dps = get_dps(t, delta)
    with mpmath.workdps(dps):
        s = mpmath.mpc(0.5 + delta, t)
        xi_ld = xi_log_deriv(s)
        kappa = float(abs(xi_ld)**2)
    delta_kappa = kappa - 1.0 / delta**2
    return kappa, delta_kappa


# ═══════════════════════════════════════════════════════════════════════════
# 파트 A: 영점 수집 + NNS 계산
# ═══════════════════════════════════════════════════════════════════════════

def part_a_collect():
    log("=" * 70)
    log("파트 A: 영점 수집 + NNS 계산")
    log("=" * 70)
    log(f"  영점 수: {len(ZERO_INDICES)}개  인덱스 {ZERO_INDICES[0]}~{ZERO_INDICES[-1]}")
    log("")

    t_start = time.time()

    # 영점 좌표 수집 (이웃 영점 포함)
    # NNS_min 계산을 위해 각 n에서 t_{n-1}, t_n, t_{n+1} 필요
    # 효율: 모든 필요 인덱스를 한 번에 계산
    needed = set()
    for n in ZERO_INDICES:
        needed.add(n)
        if n > 1: needed.add(n - 1)
        needed.add(n + 1)
    needed = sorted(needed)

    log(f"  수집 영점 수 (이웃 포함): {len(needed)}개")

    t_cache = {}
    fail_count = 0
    for i, n in enumerate(needed):
        try:
            dps_n = get_dps(1000)  # 영점 위치 계산 — 표준 dps
            t_cache[n] = get_zetazero_t(n, dps=dps_n)
        except Exception as e:
            print(f"WARNING: zetazero({n}) 실패 — {e}", flush=True)
            fail_count += 1
            t_cache[n] = None
        if (i + 1) % 50 == 0:
            t_disp = f"{t_cache[n]:.2f}" if t_cache[n] is not None else "FAIL"
            log(f"  [{i+1}/{len(needed)}] n={n}, t={t_disp}")
            save()

    if fail_count > len(needed) // 2:
        log(f"⚠️ 영점 수집 실패율 과다: {fail_count}/{len(needed)} — 탐색 로직 점검 필요")
        return None

    # NNS_min 계산 (δ-독립)
    records = []
    for n in ZERO_INDICES:
        tn = t_cache.get(n)
        t_prev = t_cache.get(n - 1) if n > 1 else None
        t_next = t_cache.get(n + 1)
        if tn is None:
            print(f"WARNING: t_{n} 없음 — 건너뜀", flush=True)
            continue

        gaps = []
        if t_prev is not None:
            gaps.append(tn - t_prev)
        if t_next is not None:
            gaps.append(t_next - tn)

        if len(gaps) == 0:
            print(f"WARNING: n={n} 이웃 영점 없음 — 건너뜀", flush=True)
            continue

        nns_min = min(gaps)

        # NNS_unfolded = NNS_min × log(t/(2π))/(2π)  (Montgomery-Odlyzko 밀도)
        density = np.log(tn / (2 * np.pi)) / (2 * np.pi)
        nns_unfolded = nns_min * density

        records.append({
            'n': n,
            't': tn,
            'nns_min': nns_min,
            'nns_unfolded': nns_unfolded,
        })

    log("")
    log(f"  유효 영점: {len(records)}개 / {len(ZERO_INDICES)}개")
    log(f"  t 범위: {records[0]['t']:.2f} ~ {records[-1]['t']:.2f}")
    log(f"  NNS_min 범위: {min(r['nns_min'] for r in records):.4f} ~ {max(r['nns_min'] for r in records):.4f}")
    log(f"  NNS_unf 범위: {min(r['nns_unfolded'] for r in records):.4f} ~ {max(r['nns_unfolded'] for r in records):.4f}")
    log(f"  소요 시간: {time.time()-t_start:.1f}초")

    return records


# ═══════════════════════════════════════════════════════════════════════════
# κ 계산 — 모든 (영점, δ) 쌍
# ═══════════════════════════════════════════════════════════════════════════

def compute_all_kappa(records):
    """
    99개 영점 × 5 δ = 495회 κ 계산
    반환: data[δ] = list of (n, t, nns_min, nns_unf, kappa, delta_kappa)
    """
    log("=" * 70)
    log("κ 계산 — 모든 (영점, δ) 쌍")
    log("=" * 70)

    t_start = time.time()
    data = {d: [] for d in DELTAS}
    total = len(records) * len(DELTAS)
    done = 0

    for rec in records:
        n = rec['n']
        t = rec['t']
        nns_min = rec['nns_min']
        nns_unf = rec['nns_unfolded']

        for delta in DELTAS:
            try:
                kappa, dk = compute_kappa(t, delta)
                if not (np.isfinite(kappa) and np.isfinite(dk)):
                    print(f"WARNING: n={n}, δ={delta} → NaN/Inf, 건너뜀", flush=True)
                    continue
                data[delta].append({
                    'n': n, 't': t,
                    'nns_min': nns_min, 'nns_unfolded': nns_unf,
                    'kappa': kappa, 'delta_kappa': dk,
                    'log_t2pi': np.log(t / (2 * np.pi)),
                    'lorentzian2': delta**2 / (delta**2 + nns_min**2)**2,
                })
            except Exception as e:
                print(f"WARNING: n={n}, δ={delta} → {e}", flush=True)
            done += 1
            if done % 50 == 0:
                elapsed = time.time() - t_start
                eta = elapsed / done * (total - done)
                log(f"  [{done}/{total}] n={n}, δ={delta:.2f} → Δκ={dk:.2f}  ETA={eta:.0f}초")
                save()

    log("")
    for delta in DELTAS:
        log(f"  δ={delta:.2f}: {len(data[delta])}개 유효 / {len(records)}개")
    log(f"  총 소요: {time.time()-t_start:.1f}초")
    return data


# ═══════════════════════════════════════════════════════════════════════════
# 파트 B: δ별 1변수 회귀
# ═══════════════════════════════════════════════════════════════════════════

def part_b_1var(data):
    log("=" * 70)
    log("파트 B: δ별 1변수 회귀  Δκ ~ a(δ)·log(t/2π) + b(δ)")
    log("이론 예측: a(δ)·δ ≈ 상수 ≈ 0.048")
    log("=" * 70)

    results = {}
    log(f"{'δ':>8} {'N':>5} {'a(δ)':>12} {'b(δ)':>12} {'R²':>8} {'a·δ':>10}")
    log("-" * 65)
    for delta in DELTAS:
        recs = data[delta]
        if len(recs) < 10:
            log(f"  δ={delta:.2f}: 데이터 부족 ({len(recs)}개) — 건너뜀")
            continue
        dk = np.array([r['delta_kappa'] for r in recs])
        log_t = np.array([r['log_t2pi'] for r in recs])

        slope, intercept, r, p, se = scipy_stats.linregress(log_t, dk)
        r2 = r**2
        a_delta = slope * delta

        log(f"  {delta:>6.2f}  {len(recs):>5d}  {slope:>12.4f}  {intercept:>12.4f}  {r2:>8.4f}  {a_delta:>10.4f}")
        results[delta] = {'slope': slope, 'intercept': intercept, 'r2': r2, 'a_delta': a_delta, 'n': len(recs)}

    # a(δ)·δ 일관성 검증
    log("")
    a_deltas = [results[d]['a_delta'] for d in DELTAS if d in results]
    if len(a_deltas) >= 2:
        mean_ad = np.mean(a_deltas)
        std_ad = np.std(a_deltas)
        cv = std_ad / abs(mean_ad) * 100 if abs(mean_ad) > 1e-10 else 999
        log(f"  a(δ)·δ 통계: 평균={mean_ad:.4f}, std={std_ad:.4f}, CV={cv:.1f}%")
        log(f"  이론 예측: ~0.048")
        if cv < 20:
            log(f"  ✅ a(δ)·δ 변동 < 20% → 이론 예측 확인 (δ-독립)")
        else:
            log(f"  ⚠️ a(δ)·δ 변동 >= 20% → δ-의존성 발견")

    return results


# ═══════════════════════════════════════════════════════════════════════════
# 파트 C: δ별 2변수 회귀 (Lorentzian²)
# ═══════════════════════════════════════════════════════════════════════════

def part_c_2var(data):
    log("=" * 70)
    log("파트 C: δ별 2변수 회귀  Δκ ~ a·log(t/2π) + b + c·f₁(NNS)")
    log("f₁ = δ²/(δ²+NNS_min²)²")
    log("성공 기준: R²_2var > 0.88 (강한 양성: 모든 δ)")
    log("=" * 70)

    results = {}
    log(f"{'δ':>8} {'N':>5} {'R²_1var':>9} {'R²_2var':>9} {'ΔR²':>8} {'c(δ)':>12} {'판정':>10}")
    log("-" * 70)

    for delta in DELTAS:
        recs = data[delta]
        if len(recs) < 10:
            log(f"  δ={delta:.2f}: 데이터 부족 — 건너뜀")
            continue

        dk = np.array([r['delta_kappa'] for r in recs])
        log_t = np.array([r['log_t2pi'] for r in recs])
        f1 = np.array([r['lorentzian2'] for r in recs])

        # 1변수 OLS: dk ~ a·log_t + b
        X1 = np.column_stack([log_t, np.ones(len(log_t))])
        coef1, res1, _, _ = np.linalg.lstsq(X1, dk, rcond=None)
        ss_res1 = np.sum((dk - X1 @ coef1)**2)
        ss_tot = np.sum((dk - np.mean(dk))**2)
        r2_1var = 1 - ss_res1 / ss_tot if ss_tot > 0 else 0

        # 2변수 OLS: dk ~ a·log_t + b + c·f1
        X2 = np.column_stack([log_t, np.ones(len(log_t)), f1])
        coef2, _, _, _ = np.linalg.lstsq(X2, dk, rcond=None)
        ss_res2 = np.sum((dk - X2 @ coef2)**2)
        r2_2var = 1 - ss_res2 / ss_tot if ss_tot > 0 else 0

        delta_r2 = r2_2var - r2_1var
        c_val = coef2[2]

        judgment = "✅ 강한양성" if r2_2var > 0.88 else ("✅ 양성" if r2_2var > 0.80 else "⚠️ 약함")
        log(f"  {delta:>6.2f}  {len(recs):>5d}  {r2_1var:>9.4f}  {r2_2var:>9.4f}  {delta_r2:>8.4f}  {c_val:>12.4f}  {judgment}")
        results[delta] = {
            'r2_1var': r2_1var, 'r2_2var': r2_2var,
            'delta_r2': delta_r2, 'c': c_val,
            'coef2': coef2, 'n': len(recs)
        }

    # 요약 판정
    log("")
    r2_vals = [results[d]['r2_2var'] for d in DELTAS if d in results]
    if len(r2_vals) >= 4:
        count_pass = sum(1 for v in r2_vals if v > 0.88)
        count_good = sum(1 for v in r2_vals if v > 0.80)
        log(f"  R²_2var > 0.88: {count_pass}/{len(r2_vals)}개 δ")
        log(f"  R²_2var > 0.80: {count_good}/{len(r2_vals)}개 δ")
        if count_pass == len(r2_vals):
            log("  ★★ 강한 양성: 모든 δ에서 R² > 0.88 — Lorentzian² 보편성 확인")
        elif count_pass >= 4:
            log("  ★ 양성: 4/5 이상 δ에서 R² > 0.88")
        elif count_good >= 4:
            log("  ★ 조건부 양성: 4/5 이상 δ에서 R² > 0.80")
        else:
            log("  ⚠️ 음성: 3개 이상 δ에서 R² < 0.80")

    return results


# ═══════════════════════════════════════════════════════════════════════════
# 파트 D: 통합 회귀 (5δ × 99 = 495 데이터)
# ═══════════════════════════════════════════════════════════════════════════

def part_d_unified(data):
    log("=" * 70)
    log("파트 D: 통합 회귀 (495 = 99×5 데이터)")
    log("모델: Δκ ~ (a/δ)·log(t/2π) + b/δ² + c·δ²/(δ²+NNS_min²)²")
    log("자유 파라미터 3개. 성공 기준: 통합 R² > 0.90")
    log("=" * 70)

    rows = []
    for delta in DELTAS:
        for r in data[delta]:
            rows.append({
                'dk': r['delta_kappa'],
                'log_t_over_delta': r['log_t2pi'] / delta,
                'inv_delta2': 1.0 / delta**2,
                'lorentzian2': r['lorentzian2'],
                'delta': delta, 't': r['t'], 'n': r['n'],
            })

    if len(rows) < 50:
        log("⚠️ 데이터 부족 — 파트 D 건너뜀")
        return None

    log(f"  총 데이터: {len(rows)}개 (목표: 495)")

    dk = np.array([r['dk'] for r in rows])
    X_log = np.array([r['log_t_over_delta'] for r in rows])  # log(t/2π)/δ
    X_inv = np.array([r['inv_delta2'] for r in rows])         # 1/δ²
    X_f1 = np.array([r['lorentzian2'] for r in rows])         # δ²/(δ²+NNS²)²

    X = np.column_stack([X_log, X_inv, X_f1])
    coef, _, _, _ = np.linalg.lstsq(X, dk, rcond=None)
    a_unified, b_unified, c_unified = coef

    dk_pred = X @ coef
    ss_res = np.sum((dk - dk_pred)**2)
    ss_tot = np.sum((dk - np.mean(dk))**2)
    r2_unified = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    log("")
    log(f"  회귀 계수:")
    log(f"    a (log(t/2π)/δ 앞) = {a_unified:.4f}  (이론 예측: ~0.048)")
    log(f"    b (1/δ² 앞)        = {b_unified:.4f}  (예상: ~0, 잔차 보정)")
    log(f"    c (Lorentzian² 앞) = {c_unified:.4f}")
    log("")
    log(f"  통합 R² = {r2_unified:.4f}")
    log(f"  성공 기준 (R² > 0.90): {'✅ 달성' if r2_unified > 0.90 else ('✅ 근접 (>0.85)' if r2_unified > 0.85 else '⚠️ 미달')}")

    # δ별 예측 품질 분해
    log("")
    log(f"{'δ':>8} {'N':>5} {'RMSE':>10} {'R² (부분)':>12}")
    log("-" * 40)
    for delta in DELTAS:
        idx = [i for i, r in enumerate(rows) if r['delta'] == delta]
        if len(idx) == 0:
            continue
        dk_d = dk[idx]
        pred_d = dk_pred[idx]
        rmse = np.sqrt(np.mean((dk_d - pred_d)**2))
        ss_r_d = np.sum((dk_d - pred_d)**2)
        ss_t_d = np.sum((dk_d - np.mean(dk_d))**2)
        r2_d = 1 - ss_r_d / ss_t_d if ss_t_d > 0 else 0
        log(f"  {delta:>6.2f}  {len(idx):>5d}  {rmse:>10.4f}  {r2_d:>12.4f}")

    return {
        'r2': r2_unified,
        'a': a_unified, 'b': b_unified, 'c': c_unified,
        'n_total': len(rows)
    }


# ═══════════════════════════════════════════════════════════════════════════
# 파트 E: δ=0.01 수치 정밀도 진단
# ═══════════════════════════════════════════════════════════════════════════

def part_e_precision(data):
    log("=" * 70)
    log("파트 E: δ=0.01 수치 정밀도 진단")
    log("  1/δ²=10000. Δκ = κ - 10000 → 유효 자릿수 감소 위험.")
    log("  dps=80으로 상향하여 계산. R²가 급락하면 수치 아티팩트 보고.")
    log("=" * 70)

    d001 = data.get(0.01, [])
    d003 = data.get(0.03, [])

    if len(d001) < 5:
        log("⚠️ δ=0.01 데이터 부족 — 진단 불가")
        return

    # Δκ 통계
    dk001 = np.array([r['delta_kappa'] for r in d001])
    dk003 = np.array([r['delta_kappa'] for r in d003]) if d003 else None

    log(f"  δ=0.01  N={len(dk001)}  Δκ: 평균={np.mean(dk001):.2f}, std={np.std(dk001):.2f}")
    log(f"  δ=0.01  Δκ 범위: [{np.min(dk001):.2f}, {np.max(dk001):.2f}]")
    if dk003 is not None:
        log(f"  δ=0.03  N={len(dk003)}  Δκ: 평균={np.mean(dk003):.2f}, std={np.std(dk003):.2f}")

    # 이론 스케일링 검증: Δκ(0.01)/Δκ(0.03) ≈ (0.03/0.01)² × (수정 인자)
    # 1변수 log(t) 항: a(δ)·log(t)/δ → δ=0.01은 δ=0.03 대비 3배
    # 그러므로 Δκ_mean(0.01) / Δκ_mean(0.03) ≈ 3 (대략)
    if dk003 is not None and len(dk003) > 0:
        ratio = np.mean(dk001) / np.mean(dk003)
        log(f"  Δκ 평균 비율 (0.01/0.03) = {ratio:.2f}  (이론 예측: ~3)")

    # 수치 노이즈 진단: Δκ < 0 인 비율 (물리적으로 불가 — 아티팩트 신호)
    neg_count = np.sum(dk001 < 0)
    log(f"  Δκ < 0 개수: {neg_count}/{len(dk001)} {'⚠️ 수치 아티팩트 의심' if neg_count > 5 else '✅ 정상'}")

    # R² vs δ=0.03 비교
    log_t001 = np.array([r['log_t2pi'] for r in d001])
    f1_001 = np.array([r['lorentzian2'] for r in d001])
    X2 = np.column_stack([log_t001, np.ones(len(log_t001)), f1_001])
    coef2, _, _, _ = np.linalg.lstsq(X2, dk001, rcond=None)
    ss_res = np.sum((dk001 - X2 @ coef2)**2)
    ss_tot = np.sum((dk001 - np.mean(dk001))**2)
    r2_001 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    log(f"  δ=0.01 R²_2var = {r2_001:.4f}")

    if r2_001 < 0.70:
        log("  ⚠️ δ=0.01 R² 급락 → 수치 아티팩트 의심. 1/δ²=10000 빼기 과정에서 정밀도 손실.")
    else:
        log("  ✅ δ=0.01 R² 정상 — 수치 안정")


# ═══════════════════════════════════════════════════════════════════════════
# 종합 판정
# ═══════════════════════════════════════════════════════════════════════════

def final_judgment(part_b, part_c, part_d):
    log("=" * 70)
    log("★ 종합 판정")
    log("=" * 70)

    r2_2var_vals = [part_c.get(d, {}).get('r2_2var', 0) for d in DELTAS if d in part_c]
    count_088 = sum(1 for v in r2_2var_vals if v > 0.88)
    count_080 = sum(1 for v in r2_2var_vals if v > 0.80)

    # a(δ)·δ 변동
    a_deltas = [part_b.get(d, {}).get('a_delta', None) for d in DELTAS if d in part_b]
    a_deltas = [v for v in a_deltas if v is not None]
    if a_deltas:
        cv_ad = np.std(a_deltas) / abs(np.mean(a_deltas)) * 100 if np.mean(a_deltas) != 0 else 999
    else:
        cv_ad = 999

    r2_unified = part_d['r2'] if part_d else 0

    log("")
    log(f"  성공 기준 체크:")
    log(f"    5개 δ 전부 R²_2var > 0.88: {count_088}/5  {'✅' if count_088 == 5 else '❌'}")
    log(f"    통합 R² > 0.90: {r2_unified:.4f}  {'✅' if r2_unified > 0.90 else '❌'}")
    log(f"    a(δ)·δ 변동 < 20%: {cv_ad:.1f}%  {'✅' if cv_ad < 20 else '❌'}")
    log("")

    if count_088 == 5 and r2_unified > 0.90 and cv_ad < 20:
        verdict = "★★ 강한 양성: 5개 δ 전부 R² > 0.88, 통합 R² > 0.90, a(δ)·δ 변동 <20%"
    elif count_088 >= 4 and r2_unified > 0.85:
        verdict = "★ 양성: 4/5 δ에서 R² > 0.88, 통합 R² > 0.85"
    elif count_080 >= 4:
        verdict = "★ 조건부 양성: 4/5 이상 δ에서 R² > 0.80"
    else:
        verdict = "⚠️ 음성: 3개 이상 δ에서 R² < 0.80 → δ-특이적 아티팩트"

    log(f"  {verdict}")

    log("")
    log(f"  핵심 발견:")
    if a_deltas and cv_ad < 30:
        log(f"    - a(δ)·δ ≈ {np.mean(a_deltas):.4f} ± {np.std(a_deltas):.4f} → {'δ-독립 확인' if cv_ad < 20 else '약한 δ-의존'}")
    if part_d:
        log(f"    - 통합 계수: a={part_d['a']:.4f}, b={part_d['b']:.4f}, c={part_d['c']:.4f}")
    log(f"    - Lorentzian² 모델 보편성: {count_088}/5 δ에서 R² > 0.88")

    return verdict


# ═══════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════

def main():
    t_global = time.time()
    log("결과 #35: δ-sweep — κ 차선도항의 δ-의존성 검증")
    log(f"실행 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"δ 값: {DELTAS}")
    log(f"영점 수: {len(ZERO_INDICES)}개")
    log(f"총 계산 횟수: {len(ZERO_INDICES) * len(DELTAS)}회")
    log("")

    # 파트 A: 영점 + NNS 수집
    records = part_a_collect()
    if records is None or len(records) == 0:
        log("⚠️ 파트 A 실패 — 종료")
        save()
        return
    save()

    # κ 계산 — 모든 (영점, δ) 쌍
    data = compute_all_kappa(records)
    save()

    # 파트 B: δ별 1변수 회귀
    part_b = part_b_1var(data)
    save()

    # 파트 C: δ별 2변수 회귀
    part_c = part_c_2var(data)
    save()

    # 파트 D: 통합 회귀
    part_d = part_d_unified(data)
    save()

    # 파트 E: 수치 정밀도 진단
    part_e_precision(data)
    save()

    # 종합 판정
    verdict = final_judgment(part_b, part_c, part_d)

    log("")
    log(f"총 소요 시간: {time.time()-t_global:.1f}초")
    log("")
    log(f"결과 파일: results/curvature_delta_sweep.txt")
    save()
    print("✅ 완료", flush=True)


if __name__ == '__main__':
    main()
