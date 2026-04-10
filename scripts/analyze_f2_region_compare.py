#!/usr/bin/env python3
"""
2단계: 고-F₂ 구간 vs 저-F₂ 구간 차이 특성화.
"""
import os, json, glob
import numpy as np
from scipy.stats import mannwhitneyu, ks_2samp
import mpmath
mpmath.mp.dps = 25

OVERNIGHT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "outputs", "overnight")
ANALYSIS_DIR = os.path.join(os.path.dirname(OVERNIGHT_DIR), "analysis")

def Z(t): return float(mpmath.siegelz(t))
def S_func(t):
    z = mpmath.zeta(0.5 + 1j*mpmath.mpf(t))
    return float(mpmath.arg(z)) / float(mpmath.pi)

def load_fine_regions():
    files = sorted(glob.glob(os.path.join(OVERNIGHT_DIR, "result_t*.json")))
    regions = []
    for fp in files:
        with open(fp) as f:
            d = json.load(f)
        zl, fl = d.get("zeros_list",[]), d.get("final_F2",[])
        width = d["t_max"] - d["t_min"]
        if len(zl)==len(fl) and len(zl)>0 and width <= 200 and d["t_min"] >= 10000:
            regions.append({
                "t_min": d["t_min"], "t_max": d["t_max"],
                "zeros": np.array(zl), "f2s": np.array(fl),
                "f2_mean": np.mean(fl),
            })
    return sorted(regions, key=lambda x: x["t_min"])

def main():
    regions = load_fine_regions()
    f2_means = [r["f2_mean"] for r in regions]
    mean_val = np.mean(f2_means)
    std_val = np.std(f2_means)
    
    high_thresh = mean_val + std_val
    low_thresh = mean_val - std_val * 0.5
    
    high_r = [r for r in regions if r["f2_mean"] > high_thresh]
    low_r = [r for r in regions if r["f2_mean"] < low_thresh]
    
    # 각 구간에서 특성 추출
    def extract_features(r, sample_S=5):
        zeros = r["zeros"]
        spacings = np.diff(zeros)
        norm_spacings = spacings / spacings.mean() if len(spacings) > 0 else spacings
        
        feat = {}
        feat["density"] = len(zeros) / (r["t_max"] - r["t_min"])
        feat["spacing_mean"] = spacings.mean() if len(spacings) > 0 else 0
        feat["spacing_std"] = spacings.std() if len(spacings) > 1 else 0
        feat["spacing_cv"] = feat["spacing_std"] / (feat["spacing_mean"] + 1e-15)
        
        # GUE 편차: 정규화 간격의 분산 (GUE 이론값 ≈ 0.178)
        feat["norm_spacing_var"] = np.var(norm_spacings) if len(norm_spacings) > 1 else 0
        
        # 간격 교대성 (lag-1 자기상관)
        if len(norm_spacings) > 2:
            feat["spacing_ac1"] = np.corrcoef(norm_spacings[:-1], norm_spacings[1:])[0,1]
        else:
            feat["spacing_ac1"] = 0
        
        # S(t) 샘플 (구간의 시작, 중간, 끝)
        t_samples = np.linspace(r["t_min"]+1, r["t_max"]-1, sample_S)
        s_vals = [S_func(float(t)) for t in t_samples]
        feat["S_mean"] = np.mean(s_vals)
        feat["S_std"] = np.std(s_vals)
        feat["S_range"] = max(s_vals) - min(s_vals)
        feat["S_abs_mean"] = np.mean(np.abs(s_vals))
        
        # Z(t) 영점 사이 극값 (처음 5쌍)
        extrema = []
        for j in range(min(5, len(spacings))):
            mid_t = (zeros[j] + zeros[j+1]) / 2
            extrema.append(abs(Z(float(mid_t))))
        feat["Z_extrema_mean"] = np.mean(extrema) if extrema else 0
        feat["Z_extrema_min"] = np.min(extrema) if extrema else 0
        
        # 소수 간격: ln(p)/2π에 가까운 간격 비율
        primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]
        prime_spacings = [np.log(p)/(2*np.pi) for p in primes]
        prime_hits = 0
        for s in spacings:
            for ps in prime_spacings:
                if abs(s - ps) < 0.05:
                    prime_hits += 1
                    break
        feat["prime_spacing_frac"] = prime_hits / (len(spacings) + 1)
        
        return feat
    
    lines = []
    lines.append("=" * 80)
    lines.append("  2단계: 고-F₂ vs 저-F₂ 구간 특성 비교")
    lines.append("=" * 80)
    lines.append(f"\n  고-F₂ 구간: {len(high_r)}개 (F₂ mean > {high_thresh:.6f})")
    lines.append(f"  저-F₂ 구간: {len(low_r)}개 (F₂ mean < {low_thresh:.6f})")
    
    print("  고-F₂ 구간 특성 추출...")
    high_feats = [extract_features(r) for r in high_r]
    print("  저-F₂ 구간 특성 추출...")
    low_feats = [extract_features(r) for r in low_r]
    
    keys = ["density", "spacing_mean", "spacing_cv", "norm_spacing_var",
            "spacing_ac1", "S_mean", "S_std", "S_range", "S_abs_mean",
            "Z_extrema_mean", "Z_extrema_min", "prime_spacing_frac"]
    
    lines.append(f"\n  {'특성':<22} {'고-F₂ mean':>12} {'저-F₂ mean':>12} {'차이':>8} {'p-value':>10}")
    lines.append(f"  {'-'*22} {'-'*12} {'-'*12} {'-'*8} {'-'*10}")
    
    for key in keys:
        h_vals = [f[key] for f in high_feats]
        l_vals = [f[key] for f in low_feats]
        h_mean = np.mean(h_vals)
        l_mean = np.mean(l_vals)
        diff_pct = (h_mean - l_mean) / (abs(l_mean) + 1e-15) * 100
        
        if len(h_vals) >= 2 and len(l_vals) >= 2:
            try:
                _, p = mannwhitneyu(h_vals, l_vals, alternative='two-sided')
            except:
                p = 1.0
        else:
            p = 1.0
        
        marker = " ★" if p < 0.1 else ""
        lines.append(f"  {key:<22} {h_mean:>12.6f} {l_mean:>12.6f} {diff_pct:>+7.1f}% {p:>10.4f}{marker}")
    
    # 상세 구간 목록
    lines.append(f"\n  ── 고-F₂ 구간 상세 ──")
    for r, f in zip(high_r, high_feats):
        lines.append(f"  [{int(r['t_min'])},{int(r['t_max'])}] F₂={r['f2_mean']:.6f} S_mean={f['S_mean']:+.4f} Z_ext={f['Z_extrema_mean']:.3f} prime={f['prime_spacing_frac']:.3f}")
    
    lines.append(f"\n  ── 저-F₂ 구간 상세 ──")
    for r, f in zip(low_r, low_feats):
        lines.append(f"  [{int(r['t_min'])},{int(r['t_max'])}] F₂={r['f2_mean']:.6f} S_mean={f['S_mean']:+.4f} Z_ext={f['Z_extrema_mean']:.3f} prime={f['prime_spacing_frac']:.3f}")
    
    report = "\n".join(lines)
    print(report)
    
    out = os.path.join(ANALYSIS_DIR, "f2_region_comparison.txt")
    with open(out, "w") as f:
        f.write(report)
    print(f"\n  → 저장: {out}")

if __name__ == "__main__":
    main()
