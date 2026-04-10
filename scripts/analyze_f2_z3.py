#!/usr/bin/env python3
"""
|F₂| vs Z'''(t_k) 심층 분석.
1. 함수형 피팅: |F₂| = a·|Z'''|^b + c
2. 높이별 분리: 낮은 t vs 높은 t에서 관계가 동일한지
3. Z''' 제거 후 잔차: 다른 숨은 변수가 있는지
4. |Z'''|이 작은 영점들의 공통 특징
"""
import os, sys, json, glob
import numpy as np
import mpmath
mpmath.mp.dps = 30

OVERNIGHT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "outputs", "overnight")
ANALYSIS_DIR = os.path.join(os.path.dirname(OVERNIGHT_DIR), "analysis")

def Z(t): return float(mpmath.siegelz(t))
def Z3(t, h=1e-5):
    return (Z(t+2*h) - 2*Z(t+h) + 2*Z(t-h) - Z(t-2*h)) / (2*h**3)
def Z1(t, h=1e-5):
    return (Z(t+h) - Z(t-h)) / (2*h)
def Z2(t, h=1e-5):
    return (Z(t+h) - 2*Z(t) + Z(t-h)) / (h**2)

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
    print(f"  총 영점: {n}개")
    
    rng = np.random.default_rng(123)
    idx = np.sort(rng.choice(n, min(400, n), replace=False))
    
    spacings = np.diff(zeros)
    
    print(f"  {len(idx)}개 영점 Z''', Z'', Z' 계산...")
    ts, f2s, z3s, z2s, z1s = [], [], [], [], []
    for count, i in enumerate(idx):
        t_k = zeros[i]
        ts.append(t_k)
        f2s.append(f2[i])
        z3s.append(Z3(t_k))
        z2s.append(Z2(t_k))
        z1s.append(Z1(t_k))
        if (count+1) % 100 == 0:
            print(f"    [{count+1}/{len(idx)}]")
    
    ts = np.array(ts)
    f2s = np.array(f2s)
    z3s = np.array(z3s)
    z2s = np.array(z2s)
    z1s = np.array(z1s)
    abs_z3 = np.abs(z3s)
    abs_z2 = np.abs(z2s)
    abs_z1 = np.abs(z1s)
    
    lines = []
    lines.append("=" * 72)
    lines.append("  |F₂| vs Z'''(t_k) — 심층 분석")
    lines.append("=" * 72)
    
    # 1. 함수형 탐색: log-log
    valid = (abs_z3 > 0) & (f2s > 0)
    lf2 = np.log(f2s[valid])
    lz3 = np.log(abs_z3[valid])
    slope, intercept = np.polyfit(lz3, lf2, 1)
    a_fit = np.exp(intercept)
    
    lines.append(f"\n  ── 1. 멱법칙 피팅 ──")
    lines.append(f"  |F₂| ≈ {a_fit:.6f} · |Z'''|^({slope:.4f})")
    lines.append(f"  음의 지수 → |Z'''| 작을수록 |F₂| 큼 확인")
    
    # 잔차
    predicted_lf2 = slope * lz3 + intercept
    residual = lf2 - predicted_lf2
    r2 = 1 - np.var(residual) / np.var(lf2)
    lines.append(f"  R² = {r2:.4f}  (설명력)")
    
    # 2. 높이별 분리
    lines.append(f"\n  ── 2. 높이별 |F₂|-Z''' 관계 ──")
    for lo, hi in [(10, 1000), (1000, 5000), (5000, 10000), (10000, 12000)]:
        mask = (ts >= lo) & (ts < hi) & valid
        if mask.sum() < 10:
            continue
        from scipy.stats import spearmanr
        sr, sp = spearmanr(f2s[mask], abs_z3[mask])
        sl, _ = np.polyfit(np.log(abs_z3[mask]), np.log(f2s[mask]), 1) if mask.sum() > 2 else (0, 0)
        lines.append(f"  t∈[{lo},{hi}]: n={mask.sum():3d}, Spearman={sr:+.3f} (p={sp:.1e}), 기울기={sl:.3f}")
    
    # 3. Z''' 제거 후 잔차 vs 다른 변수
    lines.append(f"\n  ── 3. Z''' 효과 제거 후 잔차 분석 ──")
    resid_f2 = residual  # log scale residual after Z''' regression
    
    from scipy.stats import spearmanr
    # 잔차 vs |Z'|
    sr1, sp1 = spearmanr(resid_f2, np.log(abs_z1[valid]+1e-15))
    lines.append(f"  잔차 vs ln|Z'|:  Spearman={sr1:+.4f} (p={sp1:.2e})")
    # 잔차 vs |Z''|
    sr2, sp2 = spearmanr(resid_f2, np.log(abs_z2[valid]+1e-15))
    lines.append(f"  잔차 vs ln|Z''|: Spearman={sr2:+.4f} (p={sp2:.2e})")
    # 잔차 vs t (높이 자체)
    sr3, sp3 = spearmanr(resid_f2, ts[valid])
    lines.append(f"  잔차 vs t:       Spearman={sr3:+.4f} (p={sp3:.2e})")
    # 잔차 vs ln(t)
    sr4, sp4 = spearmanr(resid_f2, np.log(ts[valid]))
    lines.append(f"  잔차 vs ln(t):   Spearman={sr4:+.4f} (p={sp4:.2e})")
    
    # 4. |Z'''| 가장 작은 10개 (= |F₂| 가장 큰 후보)
    small_z3_idx = np.argsort(abs_z3)[:15]
    lines.append(f"\n  ── 4. |Z'''| 최소 15개 영점 (F₂ 난점 후보) ──")
    lines.append(f"  {'t':>11} {'|F₂|':>8} {'|Z_p|':>8} {'|Z_pp|':>9} {'|Z_ppp|':>10}")
    lines.append(f"  {'-'*11} {'-'*8} {'-'*8} {'-'*9} {'-'*10}")
    for i in small_z3_idx:
        lines.append(f"  {ts[i]:11.2f} {f2s[i]:8.4f} {abs_z1[i]:8.3f} {abs_z2[i]:9.3f} {abs_z3[i]:10.3f}")
    
    # 5. 부호 분석: Z''' > 0 vs Z''' < 0
    pos = z3s > 0
    neg = z3s < 0
    lines.append(f"\n  ── 5. Z''' 부호별 |F₂| ──")
    lines.append(f"  Z''' > 0: n={pos.sum()}, mean|F₂|={f2s[pos].mean():.6f}, max={f2s[pos].max():.6f}")
    lines.append(f"  Z''' < 0: n={neg.sum()}, mean|F₂|={f2s[neg].mean():.6f}, max={f2s[neg].max():.6f}")
    
    # 6. 핵심 수치
    lines.append(f"\n  ── 6. 요약 ──")
    lines.append(f"  |F₂| ∝ |Z'''|^({slope:.4f})")
    lines.append(f"  R² = {r2:.4f}")
    lines.append(f"  Z'''가 |F₂| 분산의 {r2*100:.1f}%를 설명")
    lines.append(f"  나머지 {(1-r2)*100:.1f}%는 미지 — 추가 탐색 필요")
    
    report = "\n".join(lines)
    print(report)
    
    out = os.path.join(ANALYSIS_DIR, "f2_z3_deep.txt")
    with open(out, "w") as f:
        f.write(report)
    print(f"\n  → 저장: {out}")

if __name__ == "__main__":
    main()
