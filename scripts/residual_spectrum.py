#!/usr/bin/env python3
"""
residual_spectrum.py — N4c. 앙상블 잔차/분산의 구조 분석.

N2 가 알려준 것:
  - 어려움 = 시드 간 크기 분산 (σ) 큼, 부호는 일치
  - 앙상블 이득 0 → 시스템적 bias

이 스크립트는 기존 ensemble_hard_zeros.json 을 사용해:
  1. per-seed residual = F2_seed - F2_ensemble 의 영점 인덱스 공간 FFT
  2. |F₂|-mean vs 로컬 구조 (gap, log t) 상관
  3. 어려운 영점 군집이 있는가 (인접 인덱스 연관)
  4. 시드 간 cross-correlation (어떤 두 시드가 서로 가까운 표현을 배웠나)

학습 없음. 수 초.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "outputs" / "analysis"
OVERNIGHT = ROOT / "outputs" / "overnight"
OUT_JSON = ANALYSIS / "residual_spectrum.json"
OUT_TXT = ANALYSIS / "residual_spectrum.txt"


def main():
    lines = []
    def log(s=""):
        print(s, flush=True); lines.append(s)

    log("=" * 72)
    log("  residual_spectrum — N4c")
    log("=" * 72)

    with open(ANALYSIS / "ensemble_hard_zeros.json") as f:
        ens = json.load(f)
    seeds = list(ens["per_seed_F2"].keys())
    zeros = np.array(ens["zeros_list"])
    N = len(zeros)
    F2_stack = np.stack([np.array(ens["per_seed_F2"][s]) for s in seeds], 0)  # [K,N]
    F2_ens = np.array(ens["ensemble_F2"])
    abs_mean = np.array(ens["per_zero_absF2_mean"])
    abs_std = np.array(ens["per_zero_absF2_std"])

    log(f"  K={len(seeds)} seeds, N={N} zeros")

    # ---------- 1. 잔차 = F2_seed - F2_ensemble ----------
    R = F2_stack - F2_ens[None, :]  # [K,N]
    log("\n[1] 잔차 스케일")
    log(f"    mean |R| per seed: {[f'{float(np.abs(R[k]).mean()):.4f}' for k in range(len(seeds))]}")
    log(f"    std per zero: mean={R.std(axis=0).mean():.4f}  max={R.std(axis=0).max():.4f}")

    # ---------- 2. 영점 인덱스 공간 FFT ----------
    log("\n[2] 시드 평균 잔차 |R|.mean(axis=0) 의 인덱스-공간 FFT")
    r_profile = np.abs(R).mean(axis=0)  # [N]
    r_profile_dm = r_profile - r_profile.mean()
    F = np.fft.rfft(r_profile_dm)
    P = np.abs(F) ** 2
    freqs = np.fft.rfftfreq(N, d=1.0)  # 인덱스 당 주파수
    log(f"    DC 제거 후 power top5 주파수:")
    top = np.argsort(-P)[:5]
    for k in top:
        # 주파수 → 인덱스 주기
        period = 1.0 / freqs[k] if freqs[k] > 0 else float('inf')
        log(f"      freq={freqs[k]:.4f}  period≈{period:.2f} zeros  P={P[k]:.4e}")

    # |F₂|-mean 프로파일 FFT (난이도 자체의 주기성)
    log("\n    난이도 |F₂|-mean 프로파일 FFT:")
    am_dm = abs_mean - abs_mean.mean()
    Fm = np.fft.rfft(am_dm)
    Pm = np.abs(Fm) ** 2
    top_m = np.argsort(-Pm)[:5]
    for k in top_m:
        period = 1.0 / freqs[k] if freqs[k] > 0 else float('inf')
        log(f"      freq={freqs[k]:.4f}  period≈{period:.2f} zeros  P={Pm[k]:.4e}")

    # ---------- 3. 로컬 구조 상관 ----------
    log("\n[3] 난이도 |F₂|-mean 과 로컬 구조 상관")
    gaps = np.diff(zeros)
    # 각 zero k 의 최소 이웃 간격 (양쪽 중 작은 것)
    gap_min = np.zeros(N)
    gap_min[0] = gaps[0]; gap_min[-1] = gaps[-1]
    gap_min[1:-1] = np.minimum(gaps[:-1], gaps[1:])
    log_t = np.log(zeros)

    for name, x in [("log_t", log_t), ("gap_min", gap_min)]:
        r_am = float(np.corrcoef(x, abs_mean)[0, 1])
        r_std = float(np.corrcoef(x, abs_std)[0, 1])
        log(f"    corr({name:>8s}, |F₂|-mean) = {r_am:+.4f}   "
            f"corr({name:>8s}, std) = {r_std:+.4f}")

    # ---------- 4. 어려운 영점 군집 ----------
    log("\n[4] 어려운 영점 군집 (|F₂|-mean 상위 10)")
    hard_idx = np.argsort(-abs_mean)[:10]
    log(f"    top 10 인덱스 (정렬): {sorted(hard_idx.tolist())}")
    log(f"    해당 t: {sorted(zeros[hard_idx].tolist())}")
    # 연속 인덱스 군집 탐지
    sorted_idx = sorted(hard_idx.tolist())
    runs = []
    cur = [sorted_idx[0]]
    for i in sorted_idx[1:]:
        if i - cur[-1] <= 2:
            cur.append(i)
        else:
            if len(cur) >= 2:
                runs.append(cur)
            cur = [i]
    if len(cur) >= 2:
        runs.append(cur)
    log(f"    인접 군집 (|Δidx|≤2): {runs}")

    # ---------- 5. 시드 간 cross-correlation ----------
    log("\n[5] 시드 간 F₂ 벡터 상관 (5×5)")
    for i, s1 in enumerate(seeds):
        row = []
        for j, s2 in enumerate(seeds):
            c = float(np.corrcoef(F2_stack[i], F2_stack[j])[0, 1])
            row.append(f"{c:+.3f}")
        log(f"    {s1:>5s}:  " + " ".join(row))

    # 평균 off-diag corr (1-다양성 지표)
    K = len(seeds)
    off = []
    for i in range(K):
        for j in range(K):
            if i != j:
                off.append(float(np.corrcoef(F2_stack[i], F2_stack[j])[0, 1]))
    log(f"\n    평균 off-diag corr: {np.mean(off):+.4f}  "
        f"(→ 1 이면 모두 같은 예측, → 0 이면 다양, <0 이면 반대)")

    # ---------- 6. hard-focus 영점의 시드별 위치 ----------
    log("\n[6] t=173.41, 107.17 시드별 |F₂| 랭크")
    for t_h in [173.4115, 107.1686]:
        idx = int(np.argmin(np.abs(zeros - t_h)))
        log(f"  t={zeros[idx]:.4f} (idx={idx})")
        for k, s in enumerate(seeds):
            rank = int((np.abs(F2_stack[k]) > abs(F2_stack[k, idx])).sum()) + 1
            log(f"    seed {s}: |F₂|={abs(F2_stack[k, idx]):.4f}  "
                f"난이도 랭크={rank}/{N}")

    # ---------- 저장 ----------
    out = {
        "seeds": seeds,
        "residual_shape": list(R.shape),
        "residual_mean_abs": float(np.abs(R).mean()),
        "corr_log_t_abs_mean": float(np.corrcoef(log_t, abs_mean)[0, 1]),
        "corr_gap_min_abs_mean": float(np.corrcoef(gap_min, abs_mean)[0, 1]),
        "hard_idx_top10": sorted(hard_idx.tolist()),
        "hard_clusters": runs,
        "seed_off_diag_corr_mean": float(np.mean(off)),
        "top_freqs_abs_residual": [
            {"freq": float(freqs[k]), "P": float(P[k])} for k in top],
        "top_freqs_abs_mean_difficulty": [
            {"freq": float(freqs[k]), "P": float(Pm[k])} for k in top_m],
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")


if __name__ == "__main__":
    main()
