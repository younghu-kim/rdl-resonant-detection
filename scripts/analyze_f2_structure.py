#!/usr/bin/env python3
"""
|F₂| 잔차와 영점 미세구조의 상관관계 분석.

검증 대상:
1. |F₂| vs |Z'(t_k)|  — 영점에서의 미분 (영점의 "깊이")
2. |F₂| vs s_k         — 정규화 간격 (인접 영점과의 거리)
3. |F₂| vs min|Z(τ)|   — 영점 사이 극값 (Lehmer 지표)

만약 깔끔한 함수 관계가 있다면:
→ 게이지 동역학이 영점 미세구조를 인코딩하는 메커니즘 규명 가능
"""

import os
import sys
import json
import glob
import numpy as np

# mpmath for Z'(t) computation
import mpmath
mpmath.mp.dps = 30

OVERNIGHT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "outputs", "overnight")
ANALYSIS_DIR = os.path.join(os.path.dirname(OVERNIGHT_DIR), "analysis")
os.makedirs(ANALYSIS_DIR, exist_ok=True)


def hardy_Z(t):
    """Hardy Z-function: Z(t) = e^{iθ(t)} ζ(1/2+it), 실수값."""
    return float(mpmath.siegelz(t))


def hardy_Z_derivative(t, h=1e-6):
    """Z'(t) by central difference."""
    return (hardy_Z(t + h) - hardy_Z(t - h)) / (2 * h)


def hardy_Z_between_min(t1, t2, n_sample=20):
    """두 영점 사이에서 |Z(τ)|의 최솟값 (Lehmer 측정)."""
    ts = np.linspace(float(t1), float(t2), n_sample + 2)[1:-1]
    vals = [abs(hardy_Z(float(tau))) for tau in ts]
    return min(vals) if vals else float('inf')


def load_all_results():
    """모든 result JSON에서 영점 목록과 |F₂|를 로드."""
    files = sorted(glob.glob(os.path.join(OVERNIGHT_DIR, "result_t*.json")))
    all_zeros = []
    all_f2 = []

    for fpath in files:
        with open(fpath) as f:
            data = json.load(f)

        zeros = data.get("zeros_list", [])
        f2s = data.get("final_F2", [])

        if len(zeros) == len(f2s):
            all_zeros.extend(zeros)
            all_f2.extend(f2s)
        else:
            print(f"[경고] {os.path.basename(fpath)}: zeros={len(zeros)}, F2={len(f2s)} 불일치, 스킵")

    return np.array(all_zeros), np.array(all_f2)


def compute_correlations(zeros, f2_vals, sample_size=500):
    """영점 샘플에서 Z'(t_k), 간격, 극값을 계산하고 상관관계 분석."""

    n = len(zeros)
    if n < 3:
        print("영점 부족")
        return

    # 샘플링 (전부 계산하면 너무 오래 걸림)
    if n > sample_size:
        indices = np.sort(np.random.default_rng(42).choice(n, sample_size, replace=False))
    else:
        indices = np.arange(n)

    # 정규화 간격
    spacings = np.diff(zeros)
    # 평균 간격으로 정규화
    mean_spacing_local = []
    for i in indices:
        # 로컬 평균 간격 (주변 10개)
        lo = max(0, i - 5)
        hi = min(len(spacings), i + 5)
        if lo < hi:
            mean_spacing_local.append(np.mean(spacings[lo:hi]))
        else:
            mean_spacing_local.append(1.0)

    results = []
    print(f"\n  {len(indices)}개 영점 분석 중 (Z' 계산 포함)...\n")

    for count, idx in enumerate(indices):
        t_k = zeros[idx]
        f2_k = f2_vals[idx]

        # Z'(t_k)
        Zp = hardy_Z_derivative(t_k)

        # 정규화 간격 (왼쪽)
        if idx > 0:
            s_left = spacings[idx - 1] / mean_spacing_local[count]
        else:
            s_left = float('nan')

        # 정규화 간격 (오른쪽)
        if idx < len(spacings):
            s_right = spacings[idx] / mean_spacing_local[count]
        else:
            s_right = float('nan')

        s_min = np.nanmin([s_left, s_right])

        # 진행 표시
        if (count + 1) % 50 == 0:
            print(f"    [{count+1}/{len(indices)}] t={t_k:.2f}")

        results.append({
            "t": t_k,
            "F2": f2_k,
            "Z_prime": Zp,
            "abs_Z_prime": abs(Zp),
            "s_min": s_min,
            "s_left": s_left,
            "s_right": s_right,
        })

    return results


def analyze_and_report(results):
    """상관관계 계산 및 보고."""

    f2 = np.array([r["F2"] for r in results])
    abs_zp = np.array([r["abs_Z_prime"] for r in results])
    s_min = np.array([r["s_min"] for r in results])

    # NaN 제거
    valid = ~np.isnan(s_min)
    f2_v = f2[valid]
    abs_zp_v = abs_zp[valid]
    s_min_v = s_min[valid]

    # Pearson 상관계수
    corr_f2_zp = np.corrcoef(f2_v, abs_zp_v)[0, 1]
    corr_f2_s = np.corrcoef(f2_v, s_min_v)[0, 1]
    corr_f2_inv_zp = np.corrcoef(f2_v, 1.0 / (abs_zp_v + 1e-15))[0, 1]

    # log-log 상관 (power law 탐색)
    log_f2 = np.log(f2_v + 1e-15)
    log_zp = np.log(abs_zp_v + 1e-15)
    log_s = np.log(s_min_v + 1e-15)
    corr_log_f2_log_zp = np.corrcoef(log_f2, log_zp)[0, 1]
    corr_log_f2_log_s = np.corrcoef(log_f2, log_s)[0, 1]

    # 최악 영점 (|F₂| 상위 10개)
    worst_idx = np.argsort(f2_v)[-10:][::-1]

    report = []
    report.append("=" * 72)
    report.append("  |F₂| 잔차 vs 영점 미세구조 — 상관관계 분석")
    report.append("=" * 72)
    report.append(f"\n  분석 영점 수: {len(f2_v)}")
    report.append(f"  |F₂| 범위: [{f2_v.min():.6f}, {f2_v.max():.6f}]")
    report.append(f"  |Z'| 범위: [{abs_zp_v.min():.6f}, {abs_zp_v.max():.6f}]")
    report.append(f"  s_min 범위: [{s_min_v.min():.6f}, {s_min_v.max():.6f}]")

    report.append(f"\n  ── Pearson 상관계수 ──")
    report.append(f"  corr(|F₂|, |Z'|)      = {corr_f2_zp:+.4f}")
    report.append(f"  corr(|F₂|, 1/|Z'|)    = {corr_f2_inv_zp:+.4f}")
    report.append(f"  corr(|F₂|, s_min)      = {corr_f2_s:+.4f}")

    report.append(f"\n  ── log-log 상관 (멱법칙 탐색) ──")
    report.append(f"  corr(ln|F₂|, ln|Z'|)   = {corr_log_f2_log_zp:+.4f}")
    report.append(f"  corr(ln|F₂|, ln s_min)  = {corr_log_f2_log_s:+.4f}")

    # 해석
    report.append(f"\n  ── 해석 ──")
    if abs(corr_f2_inv_zp) > 0.5:
        report.append(f"  ★ |F₂| ∝ 1/|Z'|  상관 {corr_f2_inv_zp:+.4f}")
        report.append(f"    → 영점의 미분이 작을수록(얕은 교차) F₂ 잔차가 큼")
        report.append(f"    → 게이지 동역학이 Z'(t_k)를 직접 인코딩")
    if abs(corr_f2_s) > 0.3:
        report.append(f"  ★ |F₂| vs s_min  상관 {corr_f2_s:+.4f}")
        report.append(f"    → 인접 영점 간격이 좁을수록 F₂ 잔차가 큼")
        report.append(f"    → GUE 간격 통계와 직접 연결")
    if abs(corr_log_f2_log_zp) > 0.5:
        # 기울기 추정 (power law exponent)
        slope = np.polyfit(log_zp, log_f2, 1)[0]
        report.append(f"  ★ 멱법칙: |F₂| ∝ |Z'|^{slope:.3f}")

    report.append(f"\n  ── |F₂| 상위 10 영점 ──")
    report.append(f"  {'t':>12} {'|F₂|':>10} {'|Z_prime|':>12} {'s_min':>10}")
    report.append(f"  {'-'*12} {'-'*10} {'-'*12} {'-'*10}")
    for i in worst_idx:
        r = results[int(np.where(valid)[0][i])]
        report.append(f"  {r['t']:12.4f} {r['F2']:10.6f} {r['abs_Z_prime']:12.6f} {r['s_min']:10.6f}")

    report_text = "\n".join(report)
    print(report_text)

    # 저장
    out_path = os.path.join(ANALYSIS_DIR, "f2_structure_correlation.txt")
    with open(out_path, "w") as f:
        f.write(report_text)
    print(f"\n  → 저장: {out_path}")

    # 상세 데이터도 JSON으로
    json_path = os.path.join(ANALYSIS_DIR, "f2_structure_data.json")
    with open(json_path, "w") as f:
        json.dump({
            "correlations": {
                "pearson_f2_abszp": corr_f2_zp,
                "pearson_f2_inv_abszp": corr_f2_inv_zp,
                "pearson_f2_smin": corr_f2_s,
                "loglog_f2_abszp": corr_log_f2_log_zp,
                "loglog_f2_smin": corr_log_f2_log_s,
            },
            "data": results,
        }, f, indent=2)
    print(f"  → 저장: {json_path}")


def main():
    print("\n  영점 데이터 로드 중...")
    zeros, f2_vals = load_all_results()
    print(f"  총 영점: {len(zeros)}개")

    print("\n  Z'(t_k), 간격 계산 시작...")
    results = compute_correlations(zeros, f2_vals, sample_size=500)

    if results:
        analyze_and_report(results)


if __name__ == "__main__":
    main()
