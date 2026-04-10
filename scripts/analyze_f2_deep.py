#!/usr/bin/env python3
"""
|F₂| vs 고차 구조 상관분석.
Z'와 간격에 상관 없음을 확인했으므로, 더 깊은 구조를 탐색.
"""
import os, sys, json, glob
import numpy as np
from scipy import stats
import mpmath
mpmath.mp.dps = 30

OVERNIGHT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "outputs", "overnight")
ANALYSIS_DIR = os.path.join(os.path.dirname(OVERNIGHT_DIR), "analysis")
os.makedirs(ANALYSIS_DIR, exist_ok=True)

def Z(t):
    return float(mpmath.siegelz(t))

def Z_deriv(t, n=1, h=1e-5):
    """n차 미분 (중심차분)."""
    if n == 1:
        return (Z(t+h) - Z(t-h)) / (2*h)
    elif n == 2:
        return (Z(t+h) - 2*Z(t) + Z(t-h)) / (h**2)
    elif n == 3:
        return (Z(t+2*h) - 2*Z(t+h) + 2*Z(t-h) - Z(t-2*h)) / (2*h**3)

def S_function(t):
    """S(t) = (1/pi) * arg zeta(1/2+it) — 위상 편차."""
    z = mpmath.zeta(0.5 + 1j*mpmath.mpf(t))
    return float(mpmath.arg(z)) / float(mpmath.pi)

def theta_deriv(t):
    """Riemann-Siegel theta'(t) = 위상 속도."""
    h = 1e-5
    def theta(s):
        return float(mpmath.siegeltheta(s))
    return (theta(t+h) - theta(t-h)) / (2*h)

def load_all():
    files = sorted(glob.glob(os.path.join(OVERNIGHT_DIR, "result_t*.json")))
    zeros, f2s = [], []
    for fp in files:
        with open(fp) as f:
            d = json.load(f)
        zl, fl = d.get("zeros_list",[]), d.get("final_F2",[])
        if len(zl)==len(fl):
            zeros.extend(zl); f2s.extend(fl)
    return np.array(zeros), np.array(f2s)

def main():
    print("  데이터 로드...")
    zeros, f2 = load_all()
    n = len(zeros)
    print(f"  총 영점: {n}개")
    
    # 300개 샘플
    rng = np.random.default_rng(42)
    idx = np.sort(rng.choice(n, min(300, n), replace=False))
    
    spacings = np.diff(zeros)
    
    # 계산
    data = []
    print(f"  {len(idx)}개 영점 분석 (Z'', Z''', S(t), theta')...")
    for count, i in enumerate(idx):
        t_k = zeros[i]
        f2_k = f2[i]
        
        zp1 = Z_deriv(t_k, 1)
        zp2 = Z_deriv(t_k, 2)
        zp3 = Z_deriv(t_k, 3)
        st = S_function(t_k)
        td = theta_deriv(t_k)
        
        # 로컬 영점 밀도 (±2 이내 영점 수)
        mask = np.abs(zeros - t_k) < 2.0
        local_density = mask.sum() - 1  # 자기 자신 제외
        
        # 간격 비율 (s_left / s_right) — 비대칭도
        s_left = spacings[i-1] if i > 0 else np.nan
        s_right = spacings[i] if i < len(spacings) else np.nan
        asym = s_left / s_right if (s_right and s_right > 0) else np.nan
        
        # Z의 곡률 비: |Z''| / |Z'| (영점의 "날카로움")
        sharpness = abs(zp2) / (abs(zp1) + 1e-15)
        
        data.append({
            "t": t_k, "F2": f2_k,
            "Z_p1": zp1, "Z_p2": zp2, "Z_p3": zp3,
            "S_t": st, "theta_p": td,
            "local_density": local_density,
            "asym": asym, "sharpness": sharpness,
        })
        
        if (count+1) % 50 == 0:
            print(f"    [{count+1}/{len(idx)}] t={t_k:.1f}")
    
    # 상관 분석
    f2_arr = np.array([d["F2"] for d in data])
    
    quantities = {
        "|Z'(t_k)|":      np.array([abs(d["Z_p1"]) for d in data]),
        "|Z''(t_k)|":     np.array([abs(d["Z_p2"]) for d in data]),
        "|Z'''(t_k)|":    np.array([abs(d["Z_p3"]) for d in data]),
        "Z''/ Z' (날카로움)": np.array([d["sharpness"] for d in data]),
        "S(t_k)":         np.array([d["S_t"] for d in data]),
        "|S(t_k)|":       np.array([abs(d["S_t"]) for d in data]),
        "theta'(t_k)":    np.array([d["theta_p"] for d in data]),
        "로컬밀도(±2)":    np.array([d["local_density"] for d in data], dtype=float),
        "간격비(비대칭)":   np.array([d["asym"] for d in data]),
    }
    
    lines = []
    lines.append("=" * 72)
    lines.append("  |F₂| 잔차 vs 고차 영점 구조 — 심층 상관분석")
    lines.append("=" * 72)
    lines.append(f"\n  분석 영점: {len(f2_arr)}개\n")
    lines.append(f"  {'관측량':<22} {'Pearson r':>10} {'p-value':>12} {'Spearman ρ':>12} {'p-value':>12}")
    lines.append(f"  {'-'*22} {'-'*10} {'-'*12} {'-'*12} {'-'*12}")
    
    for name, vals in quantities.items():
        valid = ~np.isnan(vals) & ~np.isnan(f2_arr)
        if valid.sum() < 10:
            continue
        fv, vv = f2_arr[valid], vals[valid]
        pr, pp = stats.pearsonr(fv, vv)
        sr, sp = stats.spearmanr(fv, vv)
        marker = " ★" if abs(sr) > 0.15 else ""
        lines.append(f"  {name:<22} {pr:>+10.4f} {pp:>12.2e} {sr:>+12.4f} {sp:>12.2e}{marker}")
    
    # 상위 10 영점 상세
    worst_idx = np.argsort(f2_arr)[-10:][::-1]
    lines.append(f"\n  ── |F₂| 상위 10 영점 상세 ──")
    lines.append(f"  {'t':>11} {'|F₂|':>8} {'|Z_p|':>8} {'|Z_pp|':>8} {'sharp':>8} {'S(t)':>8} {'dens':>5}")
    lines.append(f"  {'-'*11} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*5}")
    for i in worst_idx:
        d = data[i]
        lines.append(f"  {d['t']:11.2f} {d['F2']:8.4f} {abs(d['Z_p1']):8.3f} {abs(d['Z_p2']):8.3f} {d['sharpness']:8.3f} {d['S_t']:8.4f} {d['local_density']:5.0f}")
    
    report = "\n".join(lines)
    print(report)
    
    out = os.path.join(ANALYSIS_DIR, "f2_deep_correlation.txt")
    with open(out, "w") as f:
        f.write(report)
    print(f"\n  → 저장: {out}")
    
    # JSON
    jout = os.path.join(ANALYSIS_DIR, "f2_deep_data.json")
    with open(jout, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  → 저장: {jout}")

if __name__ == "__main__":
    main()
