#!/usr/bin/env python3
"""
3단계: |F₂| 자기상관 함수 → 클러스터 상관 길이 측정.
"""
import os, json, glob
import numpy as np

OVERNIGHT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "outputs", "overnight")
ANALYSIS_DIR = os.path.join(os.path.dirname(OVERNIGHT_DIR), "analysis")

def load_all_sequential():
    """전체 영점을 t 순서로 로드."""
    files = sorted(glob.glob(os.path.join(OVERNIGHT_DIR, "result_t*.json")))
    all_data = []
    for fp in files:
        with open(fp) as f:
            d = json.load(f)
        zl, fl = d.get("zeros_list",[]), d.get("final_F2",[])
        if len(zl)==len(fl):
            for t, f2 in zip(zl, fl):
                all_data.append((t, f2))
    all_data.sort(key=lambda x: x[0])
    return np.array([x[0] for x in all_data]), np.array([x[1] for x in all_data])

def autocorrelation(x, max_lag):
    """자기상관 함수."""
    x_centered = x - x.mean()
    var = np.var(x)
    if var == 0:
        return np.zeros(max_lag)
    result = []
    for lag in range(1, max_lag + 1):
        c = np.mean(x_centered[:-lag] * x_centered[lag:]) / var
        result.append(c)
    return np.array(result)

def main():
    zeros, f2 = load_all_sequential()
    n = len(zeros)
    spacings = np.diff(zeros)
    mean_spacing = np.mean(spacings)
    
    lines = []
    lines.append("=" * 72)
    lines.append("  3단계: |F₂| 자기상관 — 클러스터 상관 길이")
    lines.append("=" * 72)
    lines.append(f"\n  총 영점: {n}개")
    lines.append(f"  평균 간격: {mean_spacing:.6f}")
    
    # 전체 자기상관
    max_lag = min(200, n // 5)
    ac = autocorrelation(f2, max_lag)
    
    # 상관 길이: C(Δ) < 1/e 되는 첫 Δ
    corr_length = max_lag  # default
    for i, c in enumerate(ac):
        if c < 1/np.e:
            corr_length = i + 1
            break
    
    # C(Δ) < 0 되는 첫 Δ
    zero_crossing = max_lag
    for i, c in enumerate(ac):
        if c < 0:
            zero_crossing = i + 1
            break
    
    lines.append(f"\n  ── 전체 자기상관 C(Δ) ──")
    lines.append(f"  상관 길이 (C < 1/e): Δ = {corr_length} 영점 ≈ {corr_length * mean_spacing:.2f} t-단위")
    lines.append(f"  영교차 (C < 0):      Δ = {zero_crossing} 영점 ≈ {zero_crossing * mean_spacing:.2f} t-단위")
    
    lines.append(f"\n  Δ(영점)  Δ(t-단위)      C(Δ)")
    lines.append(f"  {'-'*8} {'-'*10} {'-'*10}")
    display_lags = [1,2,3,5,7,10,15,20,30,50,70,100,150,200]
    for lag in display_lags:
        if lag <= max_lag:
            t_dist = lag * mean_spacing
            bar = "█" * int(max(0, ac[lag-1]) * 40)
            lines.append(f"  {lag:>8} {t_dist:>10.2f} {ac[lag-1]:>+10.4f}  {bar}")
    
    # 높이별 분리
    lines.append(f"\n  ── 높이별 상관 길이 ──")
    ranges = [(10, 1000), (1000, 5000), (5000, 10000), (10000, 13000)]
    for lo, hi in ranges:
        mask = (zeros >= lo) & (zeros < hi)
        if mask.sum() < 50:
            continue
        f2_sub = f2[mask]
        sp_sub = np.diff(zeros[mask])
        ms = np.mean(sp_sub) if len(sp_sub) > 0 else mean_spacing
        
        sub_max_lag = min(100, len(f2_sub) // 5)
        if sub_max_lag < 5:
            continue
        ac_sub = autocorrelation(f2_sub, sub_max_lag)
        
        cl = sub_max_lag
        for i, c in enumerate(ac_sub):
            if c < 1/np.e:
                cl = i + 1
                break
        
        lines.append(f"  t∈[{lo},{hi}]: n={mask.sum()}, 상관길이={cl}영점 ≈ {cl*ms:.1f}t, C(1)={ac_sub[0]:+.4f}, C(5)={ac_sub[min(4,len(ac_sub)-1)]:+.4f}")
    
    # 알려진 스케일과 비교
    lines.append(f"\n  ── 알려진 스케일과 비교 ──")
    t_mid = np.median(zeros[zeros >= 10000])
    scale_avg = 2 * np.pi / np.log(t_mid / (2 * np.pi))  # 평균 영점 간격
    lines.append(f"  t={t_mid:.0f}에서 평균 영점 간격: 2π/ln(t/2π) = {scale_avg:.4f}")
    lines.append(f"  상관 길이: {corr_length * mean_spacing:.2f} t-단위")
    lines.append(f"  비율: 상관길이 / 평균간격 = {corr_length * mean_spacing / scale_avg:.1f} 간격")
    
    # 2π/ln(p) for small primes
    primes = [2, 3, 5, 7, 11, 13]
    lines.append(f"\n  소수 스케일 비교:")
    for p in primes:
        ps = np.log(p) / (2 * np.pi)
        lines.append(f"    p={p}: ln(p)/2π = {ps:.4f},  상관길이/이 스케일 = {corr_length * mean_spacing / ps:.1f}")
    
    report = "\n".join(lines)
    print(report)
    
    out = os.path.join(ANALYSIS_DIR, "f2_correlation_length.txt")
    with open(out, "w") as f:
        f.write(report)
    print(f"\n  → 저장: {out}")

if __name__ == "__main__":
    main()
