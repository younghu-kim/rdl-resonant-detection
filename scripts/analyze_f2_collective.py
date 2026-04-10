#!/usr/bin/env python3
"""
|F₂| vs 다중 영점 집단 배치 패턴 분석.
국소 미분(Z', Z'', Z''')으로 23%만 설명 → 집단적/위상적 구조 탐색.

테스트:
1. k-근접 이웃 간격 패턴 (2, 3, 5, 7개 이웃의 간격 분산/엔트로피)
2. 로컬 간격 시퀀스의 자기상관 (이웃 간격들이 교대하는 패턴?)
3. ∫|Z(τ)|dτ in ±δ (영점 주변 Z함수의 "총 면적" — 얕은 교차 측정)
4. Z의 로컬 푸리에 구조 (영점 근방 주파수 내용)
5. 영점의 "이탈도" — 평균 위치에서 얼마나 벗어났는지 (N(t) 잔차)
6. 이웃 영점들의 |F₂| 상관 (F₂가 클러스터링 되는지)
"""
import os, sys, json, glob
import numpy as np
from scipy.stats import spearmanr, pearsonr
import mpmath
mpmath.mp.dps = 25

OVERNIGHT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "outputs", "overnight")
ANALYSIS_DIR = os.path.join(os.path.dirname(OVERNIGHT_DIR), "analysis")

def Z(t): return float(mpmath.siegelz(t))

def Z_area(t_k, delta=0.5, n_pts=30):
    """영점 주변 ±delta에서 ∫|Z(τ)|dτ."""
    ts = np.linspace(t_k - delta, t_k + delta, n_pts)
    vals = [abs(Z(float(tau))) for tau in ts]
    return np.trapezoid(vals, ts)

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
    zeros, f2 = load_all()
    n = len(zeros)
    spacings = np.diff(zeros)
    print(f"  총 영점: {n}개")
    
    # 양쪽에 충분한 이웃이 있는 영점만
    margin = 10
    candidates = np.arange(margin, n - margin)
    rng = np.random.default_rng(77)
    idx = np.sort(rng.choice(candidates, min(400, len(candidates)), replace=False))
    
    print(f"  {len(idx)}개 영점 집단 구조 분석...")
    
    data = []
    for count, i in enumerate(idx):
        t_k = zeros[i]
        f2_k = f2[i]
        
        # --- 1. k-근접 이웃 간격 통계 ---
        # 좌우 k개씩의 간격
        features = {}
        for k in [2, 3, 5, 7]:
            left_spacings = spacings[max(0,i-k):i]
            right_spacings = spacings[i:min(len(spacings),i+k)]
            local_s = np.concatenate([left_spacings, right_spacings])
            if len(local_s) > 1:
                # 로컬 평균으로 정규화
                norm_s = local_s / local_s.mean()
                features[f"s_var_{k}"] = np.var(norm_s)          # 간격 분산 (균일 vs 불균일)
                features[f"s_range_{k}"] = norm_s.max() - norm_s.min()  # 간격 범위
                features[f"s_skew_{k}"] = float(np.mean((norm_s - 1)**3))  # 비대칭도
        
        # --- 2. 간격 교대 패턴 (자기상관) ---
        local_sp = spacings[max(0,i-5):min(len(spacings),i+5)]
        if len(local_sp) > 2:
            norm_sp = local_sp / local_sp.mean()
            # lag-1 자기상관: 교대 패턴이면 음수
            ac1 = np.corrcoef(norm_sp[:-1], norm_sp[1:])[0,1] if len(norm_sp) > 2 else 0
            features["spacing_autocorr"] = ac1
        
        # --- 3. 좌우 비대칭 (간격 패턴) ---
        left_mean = spacings[max(0,i-3):i].mean() if i > 0 else 0
        right_mean = spacings[i:min(len(spacings),i+3)].mean() if i < len(spacings) else 0
        features["lr_ratio"] = left_mean / (right_mean + 1e-15)
        features["lr_diff"] = abs(left_mean - right_mean) / (0.5*(left_mean+right_mean) + 1e-15)
        
        # --- 4. 이웃 |F₂| (클러스터링) ---
        neighbor_f2 = []
        for di in [-3, -2, -1, 1, 2, 3]:
            if 0 <= i+di < n:
                neighbor_f2.append(f2[i+di])
        features["neighbor_f2_mean"] = np.mean(neighbor_f2)
        features["neighbor_f2_max"] = np.max(neighbor_f2)
        
        # --- 5. 영점 이탈도: t_k vs 평균 위치 ---
        # N(t) ≈ (t/2π)ln(t/2πe) + 7/8
        expected_t_count = (t_k/(2*np.pi)) * np.log(t_k/(2*np.pi*np.e)) + 7/8
        actual_rank = i  # 대략적 순위
        features["rank_deviation"] = actual_rank - expected_t_count
        
        # --- 6. ∫|Z|dτ (영점 주변 면적, 50개만 — 느림) ---
        if count < 50:
            area = Z_area(t_k, delta=0.3, n_pts=20)
            features["Z_area_03"] = area
        
        features["t"] = t_k
        features["F2"] = f2_k
        data.append(features)
        
        if (count+1) % 100 == 0:
            print(f"    [{count+1}/{len(idx)}]")
    
    # --- 상관 분석 ---
    f2_arr = np.array([d["F2"] for d in data])
    
    test_keys = [
        "s_var_2", "s_var_3", "s_var_5", "s_var_7",
        "s_range_2", "s_range_3", "s_range_5", "s_range_7",
        "s_skew_3", "s_skew_5",
        "spacing_autocorr",
        "lr_ratio", "lr_diff",
        "neighbor_f2_mean", "neighbor_f2_max",
        "rank_deviation",
    ]
    
    lines = []
    lines.append("=" * 72)
    lines.append("  |F₂| vs 다중 영점 집단 배치 — 상관분석")
    lines.append("=" * 72)
    lines.append(f"\n  분석 영점: {len(f2_arr)}개\n")
    lines.append(f"  {'관측량':<24} {'Pearson r':>10} {'p-value':>12} {'Spearman ρ':>12} {'p-value':>12}")
    lines.append(f"  {'-'*24} {'-'*10} {'-'*12} {'-'*12} {'-'*12}")
    
    for key in test_keys:
        vals = np.array([d.get(key, np.nan) for d in data])
        valid = ~np.isnan(vals) & ~np.isnan(f2_arr) & ~np.isinf(vals)
        if valid.sum() < 20:
            continue
        pr, pp = pearsonr(f2_arr[valid], vals[valid])
        sr, sp = spearmanr(f2_arr[valid], vals[valid])
        marker = " ★★" if abs(sr) > 0.3 else (" ★" if abs(sr) > 0.15 else "")
        lines.append(f"  {key:<24} {pr:>+10.4f} {pp:>12.2e} {sr:>+12.4f} {sp:>12.2e}{marker}")
    
    # Z_area (50개만)
    area_data = [(d["F2"], d["Z_area_03"]) for d in data if "Z_area_03" in d]
    if len(area_data) > 10:
        af2, aarea = zip(*area_data)
        af2, aarea = np.array(af2), np.array(aarea)
        pr, pp = pearsonr(af2, aarea)
        sr, sp = spearmanr(af2, aarea)
        marker = " ★★" if abs(sr) > 0.3 else (" ★" if abs(sr) > 0.15 else "")
        lines.append(f"  {'Z_area(±0.3)':<24} {pr:>+10.4f} {pp:>12.2e} {sr:>+12.4f} {sp:>12.2e}{marker}")
    
    # 이웃 F₂ 클러스터링 상세
    lines.append(f"\n  ── 이웃 |F₂| 클러스터링 상세 ──")
    nf2_mean = np.array([d["neighbor_f2_mean"] for d in data])
    sr_n, sp_n = spearmanr(f2_arr, nf2_mean)
    lines.append(f"  자신의 |F₂| vs 이웃평균 |F₂|: Spearman={sr_n:+.4f} (p={sp_n:.2e})")
    if abs(sr_n) > 0.15:
        lines.append(f"  → |F₂|가 공간적으로 클러스터링됨! 국소 현상이 아닌 영역 현상")
    
    # F₂ 상위 20개의 간격 패턴
    worst20 = np.argsort(f2_arr)[-20:][::-1]
    lines.append(f"\n  ── |F₂| 상위 20개의 집단 구조 ──")
    lines.append(f"  {'t':>11} {'|F₂|':>8} {'s_var5':>8} {'autocorr':>9} {'nbr_F2':>8} {'lr_diff':>8}")
    lines.append(f"  {'-'*11} {'-'*8} {'-'*8} {'-'*9} {'-'*8} {'-'*8}")
    for i in worst20:
        d = data[i]
        lines.append(f"  {d['t']:11.2f} {d['F2']:8.4f} {d.get('s_var_5',0):8.4f} {d.get('spacing_autocorr',0):9.4f} {d['neighbor_f2_mean']:8.4f} {d.get('lr_diff',0):8.4f}")
    
    report = "\n".join(lines)
    print(report)
    
    out = os.path.join(ANALYSIS_DIR, "f2_collective.txt")
    with open(out, "w") as f:
        f.write(report)
    print(f"\n  → 저장: {out}")

if __name__ == "__main__":
    main()
