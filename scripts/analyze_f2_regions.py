#!/usr/bin/env python3
"""
1단계: 임계선의 "영역" 정의
100 단위 구간별 평균 |F₂|를 시계열로 나열, 고/저 구간 분류.
"""
import os, json, glob
import numpy as np

OVERNIGHT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "outputs", "overnight")
ANALYSIS_DIR = os.path.join(os.path.dirname(OVERNIGHT_DIR), "analysis")
os.makedirs(ANALYSIS_DIR, exist_ok=True)

def load_all():
    """모든 result JSON 로드, 구간별 + 전체 영점별."""
    files = sorted(glob.glob(os.path.join(OVERNIGHT_DIR, "result_t*.json")))
    regions = []
    all_zeros, all_f2 = [], []
    for fp in files:
        with open(fp) as f:
            d = json.load(f)
        zl = d.get("zeros_list", [])
        fl = d.get("final_F2", [])
        if len(zl) != len(fl) or len(zl) == 0:
            continue
        regions.append({
            "t_min": d["t_min"], "t_max": d["t_max"],
            "n_zeros": len(zl),
            "f2_mean": np.mean(fl),
            "f2_median": np.median(fl),
            "f2_max": np.max(fl),
            "f2_std": np.std(fl),
            "zeros": zl, "f2s": fl,
        })
        all_zeros.extend(zl)
        all_f2.extend(fl)
    return regions, np.array(all_zeros), np.array(all_f2)

def main():
    regions, all_zeros, all_f2 = load_all()
    
    # 100 단위 구간만 (t >= 10000)
    fine_regions = [r for r in regions if (r["t_max"] - r["t_min"]) <= 200]
    wide_regions = [r for r in regions if (r["t_max"] - r["t_min"]) > 200]
    
    lines = []
    lines.append("=" * 80)
    lines.append("  1단계: 임계선 영역 분류 — |F₂| 지형도")
    lines.append("=" * 80)
    
    # 전체 통계
    global_mean = all_f2.mean()
    global_std = all_f2.std()
    lines.append(f"\n  전체 영점: {len(all_f2)}개")
    lines.append(f"  전체 |F₂|: mean={global_mean:.6f}, std={global_std:.6f}")
    
    # 넓은 구간 요약
    lines.append(f"\n  ── 넓은 구간 (t < 10000) ──")
    lines.append(f"  {'범위':<20} {'영점':>6} {'F₂ mean':>10} {'F₂ max':>10}")
    for r in sorted(wide_regions, key=lambda x: x["t_min"]):
        label = f"[{int(r['t_min'])},{int(r['t_max'])}]"
        lines.append(f"  {label:<20} {r['n_zeros']:>6} {r['f2_mean']:>10.6f} {r['f2_max']:>10.6f}")
    
    # 100 단위 구간 시계열
    if fine_regions:
        fine_sorted = sorted(fine_regions, key=lambda x: x["t_min"])
        f2_means = [r["f2_mean"] for r in fine_sorted]
        f2_medians = [r["f2_median"] for r in fine_sorted]
        
        fine_mean = np.mean(f2_means)
        fine_std = np.std(f2_means)
        threshold_high = fine_mean + fine_std
        threshold_low = fine_mean - fine_std * 0.5
        
        lines.append(f"\n  ── 100 단위 구간 (t ≥ 10000): {len(fine_sorted)}개 ──")
        lines.append(f"  구간 F₂ 평균의 평균: {fine_mean:.6f}")
        lines.append(f"  구간 F₂ 평균의 표준편차: {fine_std:.6f}")
        lines.append(f"  고-F₂ 임계: > {threshold_high:.6f}")
        lines.append(f"  저-F₂ 임계: < {threshold_low:.6f}")
        
        lines.append(f"\n  {'구간':<18} {'영점':>5} {'F₂ mean':>9} {'F₂ max':>9} {'분류':>6}")
        lines.append(f"  {'-'*18} {'-'*5} {'-'*9} {'-'*9} {'-'*6}")
        
        high_regions = []
        low_regions = []
        mid_regions = []
        
        for r in fine_sorted:
            label = f"[{int(r['t_min'])},{int(r['t_max'])}]"
            if r["f2_mean"] > threshold_high:
                cls = "▲ 고"
                high_regions.append(r)
            elif r["f2_mean"] < threshold_low:
                cls = "▽ 저"
                low_regions.append(r)
            else:
                cls = "─ 중"
                mid_regions.append(r)
            
            # 막대 그래프
            bar_len = int(r["f2_mean"] / fine_mean * 20)
            bar = "█" * min(bar_len, 40)
            lines.append(f"  {label:<18} {r['n_zeros']:>5} {r['f2_mean']:>9.6f} {r['f2_max']:>9.6f} {cls}  {bar}")
        
        lines.append(f"\n  ── 분류 요약 ──")
        lines.append(f"  고-F₂ 구간: {len(high_regions)}개")
        lines.append(f"  중-F₂ 구간: {len(mid_regions)}개")
        lines.append(f"  저-F₂ 구간: {len(low_regions)}개")
        
        # 연속성 분석: 고-F₂가 연속으로 나타나는 패턴
        lines.append(f"\n  ── 구간 경계 분석: 연속 패턴 ──")
        labels = []
        for r in fine_sorted:
            if r["f2_mean"] > threshold_high:
                labels.append("H")
            elif r["f2_mean"] < threshold_low:
                labels.append("L")
            else:
                labels.append("M")
        
        pattern = "".join(labels)
        lines.append(f"  패턴: {pattern}")
        
        # 연속 구간 길이
        runs = []
        current = pattern[0]
        count = 1
        for c in pattern[1:]:
            if c == current:
                count += 1
            else:
                runs.append((current, count))
                current = c
                count = 1
        runs.append((current, count))
        
        h_runs = [length for label, length in runs if label == "H"]
        l_runs = [length for label, length in runs if label == "L"]
        lines.append(f"  고-F₂ 연속 길이: {h_runs} (평균 {np.mean(h_runs):.1f})" if h_runs else "  고-F₂ 연속: 없음")
        lines.append(f"  저-F₂ 연속 길이: {l_runs} (평균 {np.mean(l_runs):.1f})" if l_runs else "  저-F₂ 연속: 없음")
        
        # 전이 빈도
        transitions = sum(1 for i in range(len(pattern)-1) if pattern[i] != pattern[i+1])
        lines.append(f"  전이 횟수: {transitions} / {len(pattern)-1} ({transitions/(len(pattern)-1)*100:.0f}%)")
        if transitions / (len(pattern)-1) > 0.7:
            lines.append(f"  → 급변 패턴 (구간 간 빠르게 전환)")
        else:
            lines.append(f"  → 점진 패턴 (같은 성격이 이어짐)")
    
    report = "\n".join(lines)
    print(report)
    
    out = os.path.join(ANALYSIS_DIR, "f2_regions.txt")
    with open(out, "w") as f:
        f.write(report)
    print(f"\n  → 저장: {out}")

if __name__ == "__main__":
    main()
