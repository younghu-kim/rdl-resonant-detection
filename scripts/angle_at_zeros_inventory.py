#!/usr/bin/env python3
"""
angle_at_zeros_inventory.py — "영점의 각" 데이터 수준 인벤토리 (학습 없음).

네 개 관측량을 한 번에 계산하여 outputs/analysis/angle_at_zeros_inventory.{json,txt}
에 저장한다. 모든 계산은 기존 outputs/ 데이터만 사용하며 수 초면 끝난다.

(D1) 고립실험 arg Z_out 의 시드 간 일관성 (per-영점 원형 표준편차)
(D2) 직접 계산된 영점 접근각 arg ξ'(t_k) ∈ {+1,-1}  (xi_real 유한차분)
(D3) 해석신호 z = xi_real + i H[xi_real] 의 arg 를 영점에서 평가 (Prop 8.1)
(D4) 49개 overnight 밴드의 F₂ 난이도 요약
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import torch
from scipy.signal import hilbert
from scipy.stats import circstd

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "outputs" / "analysis"
OVERNIGHT = ROOT / "outputs" / "overnight"
OUT_JSON = ANALYSIS / "angle_at_zeros_inventory.json"
OUT_TXT = ANALYSIS / "angle_at_zeros_inventory.txt"


def circ_dist(a, b):
    d = a - b
    return np.arctan2(np.sin(d), np.cos(d))


def circ_std_axis(angles):
    # angles shape [S, N]; compute per-column circular std (axis=0)
    s = np.sin(angles).mean(axis=0)
    c = np.cos(angles).mean(axis=0)
    R = np.sqrt(s * s + c * c)
    R = np.clip(R, 1e-12, 1.0)
    return np.sqrt(-2.0 * np.log(R))


def main():
    lines = []
    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    log("=" * 72)
    log("  angle_at_zeros_inventory — 학습 없는 데이터 수준 분석")
    log("=" * 72)

    # ------------------------------------------------------------------
    # 0. 영점 목록 및 ξ 캐시 로드 (t∈[100,200])
    # ------------------------------------------------------------------
    iso_path = ANALYSIS / "phase_isolation_experiment.json"
    with open(iso_path) as f:
        iso = json.load(f)
    zeros = np.array(iso["zeros"], dtype=np.float64)
    nonzeros = np.array(iso["nonzeros"], dtype=np.float64)
    log(f"\n[0] 로드: {len(zeros)} zeros, {len(nonzeros)} nonzeros in t∈[100,200]")

    cache = torch.load(OVERNIGHT / "xi_cache_t100-200_n500.pt",
                       weights_only=True)
    t = cache["t"].numpy()
    xr = cache["xi_real"].numpy()
    xi = cache["xi_imag"].numpy()
    log(f"    xi cache: {len(t)} samples, xi_imag/xi_real max ratio "
        f"= {np.abs(xi).max()/np.abs(xr).max():.2e}")

    # ------------------------------------------------------------------
    # (D1) 고립실험 arg Z_out 의 시드 간 일관성
    # ------------------------------------------------------------------
    log("\n" + "=" * 72)
    log("  (D1) arg Z_out(t_k) 시드 간 일관성 — 신호 존재 증명")
    log("=" * 72)
    seeds = sorted(iso["results"].keys())
    phi_mat = np.array([iso["results"][s]["phi_z_circ"] for s in seeds])  # [S, N]
    log(f"    shape: seeds={seeds}, per-영점 arg 행렬 {phi_mat.shape}")

    cstd = circ_std_axis(phi_mat)  # per-zero circular std across seeds
    log(f"    per-영점 circular std (평균 0.0=완전일치, π/√3≈1.81=균일):")
    log(f"        mean={cstd.mean():.4f}  median={np.median(cstd):.4f}  "
        f"max={cstd.max():.4f}  min={cstd.min():.4f}")

    # 가장 일관된/비일관된 영점 상위 5개
    order = np.argsort(cstd)
    log("\n    가장 시드-일관적 상위 5 영점:")
    for idx in order[:5]:
        log(f"        t_{idx:2d}={zeros[idx]:.4f}  cstd={cstd[idx]:.4f}  "
            f"angles={[f'{phi_mat[i,idx]:+.3f}' for i in range(len(seeds))]}")
    log("    가장 시드-비일관적 상위 5 영점:")
    for idx in order[-5:]:
        log(f"        t_{idx:2d}={zeros[idx]:.4f}  cstd={cstd[idx]:.4f}  "
            f"angles={[f'{phi_mat[i,idx]:+.3f}' for i in range(len(seeds))]}")

    # 비영점과 비교
    phi_nz_mat = np.array([iso["results"][s]["phi_nz_circ"] for s in seeds])
    cstd_nz = circ_std_axis(phi_nz_mat)
    log(f"\n    비교: 비영점 cstd 평균={cstd_nz.mean():.4f}  "
        f"max={cstd_nz.max():.4f}")
    log(f"    → 영점 cstd mean({cstd.mean():.3f}) vs 비영점({cstd_nz.mean():.3f}): "
        f"{'영점이 더 일관' if cstd.mean() < cstd_nz.mean() else '영점이 더 노이지'}")

    # ------------------------------------------------------------------
    # (D2) 직접 계산된 영점 접근각 arg ξ'(t_k)
    # ------------------------------------------------------------------
    log("\n" + "=" * 72)
    log("  (D2) arg ξ'(t_k) — xi_real 유한차분 접근각 (확정값)")
    log("=" * 72)
    dt = np.diff(t).mean()
    dxr = np.gradient(xr, t)  # same grid

    # 각 영점에서 ξ'(t_k) 값 (선형 보간)
    xi_prime_at_zero = np.interp(zeros, t, dxr)
    sign_xi_prime = np.sign(xi_prime_at_zero)
    n_pos = int((sign_xi_prime > 0).sum())
    n_neg = int((sign_xi_prime < 0).sum())
    log(f"    영점 50개에서 sign(ξ'(t_k)): +={n_pos}, -={n_neg}")
    log(f"    ξ'(t_k) 절대값 범위: "
        f"[{np.abs(xi_prime_at_zero).min():.2e}, "
        f"{np.abs(xi_prime_at_zero).max():.2e}]")
    log(f"    |ξ'| 의 t-의존성 (Lehmer 후보 = |ξ'| 가 작은 영점):")
    order_lehmer = np.argsort(np.abs(xi_prime_at_zero))
    log("      |ξ'| 최소 상위 5 (Lehmer 후보):")
    for idx in order_lehmer[:5]:
        log(f"        t_{idx:2d}={zeros[idx]:.4f}  |ξ'|={np.abs(xi_prime_at_zero[idx]):.3e}  "
            f"sign={int(sign_xi_prime[idx])}")
    log("      |ξ'| 최대 상위 3:")
    for idx in order_lehmer[-3:]:
        log(f"        t_{idx:2d}={zeros[idx]:.4f}  |ξ'|={np.abs(xi_prime_at_zero[idx]):.3e}")

    # 인접 영점 간 sign alternation 검증: 단순 영점이면 부호가 교차해야 함
    alternations = int((np.diff(sign_xi_prime) != 0).sum())
    log(f"\n    부호 교차수 (49개 인접쌍 중): {alternations}  "
        f"(모두 simple zero 이면 49)")

    # ------------------------------------------------------------------
    # (D3) 해석신호 z(t) = xi_real + i H[xi_real] 의 arg 를 영점에서 평가
    # ------------------------------------------------------------------
    log("\n" + "=" * 72)
    log("  (D3) arg(analytic signal)(t_k) — Prop 8.1 직접 검증")
    log("=" * 72)
    # xi_real 스케일이 너무 작아서 수치 안정성 위해 정규화
    xr_norm = xr / np.abs(xr).max()
    z_analytic = hilbert(xr_norm)  # complex analytic signal
    arg_z = np.angle(z_analytic)

    arg_z_at_zero = np.interp(zeros, t, np.unwrap(arg_z))
    # wrap back to (-π, π]
    arg_z_at_zero_wrapped = np.arctan2(np.sin(arg_z_at_zero),
                                        np.cos(arg_z_at_zero))
    # |arg - π/2| 과 |arg + π/2| 중 작은 값 → ±π/2 포화도
    dev_plus = np.abs(circ_dist(arg_z_at_zero_wrapped, np.pi/2))
    dev_minus = np.abs(circ_dist(arg_z_at_zero_wrapped, -np.pi/2))
    dev_half_pi = np.minimum(dev_plus, dev_minus)
    sign_half_pi = np.where(dev_plus < dev_minus, +1, -1)
    log(f"    |arg z - (±π/2)| 통계 (0 = Prop 8.1 정확 만족):")
    log(f"        mean={dev_half_pi.mean():.4f}  median={np.median(dev_half_pi):.4f}  "
        f"max={dev_half_pi.max():.4f}")
    log(f"    sign 분포 (+π/2 vs -π/2): "
        f"+={int((sign_half_pi>0).sum())}  -={int((sign_half_pi<0).sum())}")

    # (D2) sign(ξ') 와 (D3) sign(arg z) 의 대응
    agree = int((sign_xi_prime == sign_half_pi).sum())
    log(f"\n    sign(ξ'(t_k)) vs sign(arg z(t_k)-근사): "
        f"일치 {agree}/50  "
        f"({'강한 상관' if agree >= 45 or agree <= 5 else '약한 상관'})")

    # ------------------------------------------------------------------
    # (D1') 시드간 일관성이 높은 각과 (D2) 접근각 부호의 대응
    # ------------------------------------------------------------------
    # 각 영점의 3시드 circular-mean arg Z_out
    phi_mean_per_zero = np.arctan2(
        np.sin(phi_mat).mean(axis=0), np.cos(phi_mat).mean(axis=0)
    )
    sign_phi_net = np.sign(phi_mean_per_zero)
    agree_net = int((sign_phi_net == sign_xi_prime).sum())
    log(f"\n    sign(mean arg Z_out) vs sign(ξ'(t_k)): "
        f"일치 {agree_net}/50  (랜덤 기대값 25)")

    # ------------------------------------------------------------------
    # (D4) 49개 overnight 밴드 F₂ 난이도 요약
    # ------------------------------------------------------------------
    log("\n" + "=" * 72)
    log("  (D4) overnight 밴드별 F₂ 난이도 — detection rate 및 F2_max 의 t-의존성")
    log("=" * 72)
    def _f(x):
        if x is None:
            return None
        if isinstance(x, str) and "/" in x:
            a, b = x.split("/")
            try:
                return float(a) / float(b)
            except (TypeError, ValueError, ZeroDivisionError):
                return None
        try:
            return float(x)
        except (TypeError, ValueError):
            return None
    def _i(x):
        try:
            return int(x)
        except (TypeError, ValueError):
            return None
    band_files = sorted(OVERNIGHT.glob("result_t*.json"))
    bands = []
    for bf in band_files:
        with open(bf) as f:
            d = json.load(f)
        bands.append({
            "name": bf.stem,
            "t_min": _f(d.get("t_min")),
            "t_max": _f(d.get("t_max")),
            "n_zeros": _i(d.get("n_zeros")),
            "detection_rate": _f(d.get("detection_rate")),
            "F2_mean": _f(d.get("F2_mean")),
            "F2_max": _f(d.get("F2_max")),
            "cos2_mean": _f(d.get("cos2_mean")),
            "cos2_max": _f(d.get("cos2_max")),
            "best_val_loss": _f(d.get("best_val_loss")),
            "gauge_transition_epoch": _i(d.get("gauge_transition_epoch")),
            "hardest_zeros": d.get("hardest_zeros"),
            "easiest_zeros": d.get("easiest_zeros"),
        })
    bands.sort(key=lambda b: (b["t_min"] if b["t_min"] is not None else 0))
    log(f"    {len(bands)} 밴드 로드")
    def _row(b):
        tm = b['t_min'] if b['t_min'] is not None else 0
        tM = b['t_max'] if b['t_max'] is not None else 0
        nz = b['n_zeros'] if b['n_zeros'] is not None else 0
        dr = (b['detection_rate'] or 0) * 100
        fm = b['F2_mean'] or 0.0
        fM = b['F2_max'] or 0.0
        return (f"    [{tm:>7.0f}, {tM:>7.0f}] {nz:>4d} "
                f"{dr:>6.1f}% {fm:>10.4f} {fM:>10.4f}")
    log(f"\n    {'range':>22s} {'n':>4s} {'det%':>7s} {'F2_mean':>10s} {'F2_max':>10s}")
    for b in bands[:5]:
        log(_row(b))
    log("    ...")
    for b in bands[-5:]:
        log(_row(b))

    # 전체 집계
    det_rates = [b["detection_rate"] for b in bands if b["detection_rate"] is not None]
    f2_means = [b["F2_mean"] for b in bands if b["F2_mean"] is not None]
    log(f"\n    overall detection rate: mean={np.mean(det_rates)*100:.2f}%, "
        f"min={np.min(det_rates)*100:.1f}%, max={np.max(det_rates)*100:.1f}%")
    log(f"    F2_mean over all bands: "
        f"min={np.min(f2_means):.4f}, max={np.max(f2_means):.4f}, "
        f"mean={np.mean(f2_means):.4f}")

    total_zeros_scanned = sum(b["n_zeros"] or 0 for b in bands)
    log(f"    총 스캔된 영점: {total_zeros_scanned}")

    # hardest zero 집계 — t 의 전역 분포
    all_hardest = []
    for b in bands:
        if b["hardest_zeros"]:
            for item in b["hardest_zeros"]:
                if isinstance(item, (list, tuple)) and len(item) >= 2:
                    all_hardest.append((float(item[0]), float(item[1]),
                                        b["t_min"], b["t_max"]))
    all_hardest.sort(key=lambda x: -x[1])  # F2 큰 순
    log(f"\n    전체 밴드 통합 가장 어려운 영점 상위 10 (F₂ 큰 순):")
    for t_k, f2, tm, tM in all_hardest[:10]:
        log(f"        t={t_k:>10.4f}  F₂={f2:.4f}  band=[{tm:.0f},{tM:.0f}]")

    # ------------------------------------------------------------------
    # 저장
    # ------------------------------------------------------------------
    out = {
        "D1_seed_consistency": {
            "seeds": seeds,
            "zeros_t": zeros.tolist(),
            "phi_z_per_seed": phi_mat.tolist(),
            "phi_nz_per_seed": phi_nz_mat.tolist(),
            "circular_std_per_zero": cstd.tolist(),
            "circular_std_per_nonzero": cstd_nz.tolist(),
            "mean_cstd_zero": float(cstd.mean()),
            "mean_cstd_nonzero": float(cstd_nz.mean()),
        },
        "D2_xi_prime_sign": {
            "zeros_t": zeros.tolist(),
            "xi_prime_at_zero": xi_prime_at_zero.tolist(),
            "sign_xi_prime": sign_xi_prime.astype(int).tolist(),
            "sign_alternations": int(alternations),
            "lehmer_candidates_idx": order_lehmer[:5].tolist(),
        },
        "D3_analytic_signal": {
            "zeros_t": zeros.tolist(),
            "arg_z_at_zero": arg_z_at_zero_wrapped.tolist(),
            "deviation_from_half_pi": dev_half_pi.tolist(),
            "sign_half_pi": sign_half_pi.astype(int).tolist(),
            "mean_deviation": float(dev_half_pi.mean()),
            "agree_with_xi_prime_sign": int(agree),
        },
        "D1_vs_D2_correspondence": {
            "mean_arg_Zout_per_zero": phi_mean_per_zero.tolist(),
            "sign_mean_arg_Zout": sign_phi_net.astype(int).tolist(),
            "agree_with_xi_prime_sign": int(agree_net),
        },
        "D4_overnight_bands": bands,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")
    log(f"저장: {OUT_TXT}")


if __name__ == "__main__":
    main()
