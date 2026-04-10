#!/usr/bin/env python3
"""
analyze_zero_107.py — t=107.1686 국소 검증 (N1a).

t=173.41 과 동일한 분석을 t=107.17 에 적용. 가설: Lehmer-유사 압축 쌍의
멤버이면 peak 상대크기 < 0.2 가 나와야 한다.

출력: outputs/analysis/zero_107_analysis.{json,txt}
"""
import json
from pathlib import Path
import numpy as np
import torch
from scipy.signal import hilbert

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "outputs" / "analysis"
OVERNIGHT = ROOT / "outputs" / "overnight"
OUT_JSON = ANALYSIS / "zero_107_analysis.json"
OUT_TXT = ANALYSIS / "zero_107_analysis.txt"

T_FOCUS = 107.1686


def main():
    lines = []
    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    log("=" * 72)
    log(f"  analyze_zero_107 — t={T_FOCUS} Lehmer 가설 검증 (N1a)")
    log("=" * 72)

    cache = torch.load(OVERNIGHT / "xi_cache_t100-200_n500.pt",
                       weights_only=True)
    t = cache["t"].numpy()
    xr = cache["xi_real"].numpy()
    xr_norm = xr / np.abs(xr).max()
    z = hilbert(xr_norm)

    with open(ANALYSIS / "phase_isolation_experiment.json") as f:
        iso = json.load(f)
    zeros = np.array(iso["zeros"])
    seeds = sorted(iso["results"].keys())

    with open(ANALYSIS / "angle_at_zeros_inventory.json") as f:
        inv = json.load(f)
    cstd_per_zero = np.array(inv["D1_seed_consistency"]["circular_std_per_zero"])
    xi_prime = np.array(inv["D2_xi_prime_sign"]["xi_prime_at_zero"])

    idx = int(np.argmin(np.abs(zeros - T_FOCUS)))
    log(f"\n[1] 기본 정보")
    log(f"    target={T_FOCUS}, zero_idx={idx}, 실제값={zeros[idx]:.6f}")

    # 이웃 간격
    all_gaps = np.diff(zeros)
    gap_prev = zeros[idx] - zeros[idx-1] if idx > 0 else float('nan')
    gap_next = zeros[idx+1] - zeros[idx] if idx < len(zeros)-1 else float('nan')
    pct_prev = float((all_gaps < gap_prev).mean() * 100)
    pct_next = float((all_gaps < gap_next).mean() * 100)
    log(f"    이전 3개 zero: {zeros[max(0,idx-3):idx].tolist()}")
    log(f"    현재: {zeros[idx]:.6f}")
    log(f"    이후 3개 zero: {zeros[idx+1:idx+4].tolist()}")
    log(f"    gap_prev = {gap_prev:.4f}  (하위 {pct_prev:.0f}%)")
    log(f"    gap_next = {gap_next:.4f}  (하위 {pct_next:.0f}%)")

    # ------------------------------------------------------------------
    # peak 분석 (173.41 분석과 동일 방법)
    # ------------------------------------------------------------------
    log(f"\n[2] 국소 |ξ_real| peak 분석")
    def peak_between(t_lo, t_hi):
        m = (t > t_lo) & (t < t_hi)
        if not m.any():
            return 0.0, None
        sub = np.abs(xr[m])
        t_sub = t[m]
        j = int(np.argmax(sub))
        return float(sub[j]), float(t_sub[j])

    peaks_all = []
    for i in range(len(zeros)-1):
        p, _ = peak_between(zeros[i], zeros[i+1])
        peaks_all.append(p)
    peaks_all = np.array(peaks_all)

    for i in range(max(0, idx-2), min(len(zeros)-1, idx+3)):
        p, tp = peak_between(zeros[i], zeros[i+1])
        mark = " <-- focus" if i in (idx-1, idx) else ""
        log(f"    [{zeros[i]:.4f}, {zeros[i+1]:.4f}] peak |ξ|={p:.3e}  @t={tp:.3f}{mark}")

    def local_relative(i):
        if i < 0 or i >= len(peaks_all):
            return float('nan')
        p = peaks_all[i]
        lo = max(0, i-2); hi = min(len(peaks_all), i+3)
        neighbors = np.delete(peaks_all[lo:hi], i-lo)
        return p / neighbors.mean() if len(neighbors) else float('nan')

    rel_prev = local_relative(idx-1)
    rel_next = local_relative(idx)
    log(f"\n    peak 상대크기 (주변 이웃 평균 대비):")
    log(f"      [{zeros[idx-1]:.4f}, {zeros[idx]:.4f}] 상대 = {rel_prev:.3f}")
    log(f"      [{zeros[idx]:.4f}, {zeros[idx+1]:.4f}] 상대 = {rel_next:.3f}")
    lehmer_hit = (rel_prev < 0.2) or (rel_next < 0.2)
    log(f"    Lehmer 가설(<0.2): {'확인 ✓' if lehmer_hit else '기각 ✗'}")

    # ------------------------------------------------------------------
    # 네트워크 각 및 ξ'
    # ------------------------------------------------------------------
    log(f"\n[3] 네트워크 arg Z_out 시드 간 분포")
    a = [iso["results"][s]["phi_z_circ"][idx] for s in seeds]
    log(f"    arg Z_out: {[f'{x:+.3f}' for x in a]}  cstd={cstd_per_zero[idx]:.3f}")
    log(f"    ξ'({T_FOCUS})={xi_prime[idx]:+.3e}")

    # 이웃 영점 네트워크 비교
    log(f"\n    이웃 영점 비교:")
    for j in [idx-1, idx, idx+1]:
        if 0 <= j < len(zeros):
            aa = [iso["results"][s]["phi_z_circ"][j] for s in seeds]
            mark = " <-- focus" if j == idx else ""
            log(f"      t={zeros[j]:.4f}  cstd={cstd_per_zero[j]:.3f}  "
                f"arg={[f'{x:+.2f}' for x in aa]}{mark}")

    # ------------------------------------------------------------------
    # 저장
    # ------------------------------------------------------------------
    out = {
        "target_t": T_FOCUS,
        "target_idx": idx,
        "neighbors_before": zeros[max(0,idx-3):idx].tolist(),
        "neighbors_after": zeros[idx+1:idx+4].tolist(),
        "gap_prev": float(gap_prev),
        "gap_next": float(gap_next),
        "gap_percentile_prev": pct_prev,
        "gap_percentile_next": pct_next,
        "peak_rel_before": float(rel_prev),
        "peak_rel_after": float(rel_next),
        "lehmer_hit": bool(lehmer_hit),
        "angles_per_seed": a,
        "cstd": float(cstd_per_zero[idx]),
        "xi_prime": float(xi_prime[idx]),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")
    log(f"저장: {OUT_TXT}")


if __name__ == "__main__":
    main()
