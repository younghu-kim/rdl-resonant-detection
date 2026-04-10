#!/usr/bin/env python3
"""
difficulty_regression.py — F₂ 난이도 다중회귀 요인 분석 (N4).

목적: 영점 난이도(final_F2, 부호 무시한 절댓값)의 주 요인을 분리.
- N1/N1a/N1c/N1c+ 는 Lehmer 압축을 단일 원인으로 가정했는데 marginal.
- 여기서 여러 후보 요인을 동시에 회귀에 넣어 각자의 **부분 기여도**를 본다.

데이터:
  A) 전역: 51 밴드 × 전 영점 = 15,422 점
     features: log_t, gap_prev, gap_next, gap_min, gap_asym,
               g_prev_norm, g_next_norm, min_g_norm (Riemann-von Mangoldt)
     label: |final_F2|
  B) 풍부 (t∈[100,200], 50점): 위 전부 + cstd, |ξ'|, peak_rel
     label: |final_F2| 그리고 F₂ 자체

출력: outputs/analysis/difficulty_regression.{json,txt}
"""
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "outputs" / "analysis"
OVERNIGHT = ROOT / "outputs" / "overnight"
OUT_JSON = ANALYSIS / "difficulty_regression.json"
OUT_TXT = ANALYSIS / "difficulty_regression.txt"


def normalized_gap(g_left, g_right_midpoint_t):
    """Riemann-von Mangoldt 정규화: gap * density. 평균=1."""
    density = math.log(g_right_midpoint_t / (2 * math.pi)) / (2 * math.pi)
    return g_left * density


def collect_global():
    """모든 밴드의 영점을 features + label 로 수집."""
    result_files = sorted(OVERNIGHT.glob("result_t*.json"))
    rows = []
    for rf in result_files:
        with open(rf) as f:
            d = json.load(f)
        zeros = np.array(d["zeros_list"])
        f2 = np.array(d["final_F2"])
        if len(zeros) != len(f2) or len(zeros) < 3:
            continue
        gaps = np.diff(zeros)  # len = n-1
        # 각 영점 k (1..n-2) 에 대해 gap_prev=gaps[k-1], gap_next=gaps[k]
        for k in range(1, len(zeros) - 1):
            t = zeros[k]
            gp = gaps[k-1]
            gn = gaps[k]
            g_min = min(gp, gn)
            g_asym = abs(gp - gn) / (gp + gn + 1e-12)
            # 정규화 (밀도 logt/2π)
            density = math.log(t / (2 * math.pi)) / (2 * math.pi)
            gp_norm = gp * density
            gn_norm = gn * density
            g_min_norm = g_min * density
            rows.append({
                "band": rf.stem.replace("result_", ""),
                "t": float(t),
                "log_t": float(math.log(t)),
                "gap_prev": float(gp),
                "gap_next": float(gn),
                "gap_min": float(g_min),
                "gap_asym": float(g_asym),
                "g_prev_norm": float(gp_norm),
                "g_next_norm": float(gn_norm),
                "g_min_norm": float(g_min_norm),
                "F2": float(f2[k]),
                "abs_F2": float(abs(f2[k])),
            })
    return rows


def standardize(X):
    """각 열 표준화."""
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=1)
    sd[sd < 1e-12] = 1.0
    return (X - mu) / sd, mu, sd


def ols(X, y):
    """OLS, returns β, σ², t-stat, R²."""
    n, p = X.shape
    # 1 열 추가
    X1 = np.column_stack([np.ones(n), X])
    XtX = X1.T @ X1
    XtX_inv = np.linalg.pinv(XtX)
    beta = XtX_inv @ X1.T @ y
    y_hat = X1 @ beta
    resid = y - y_hat
    sse = float(resid @ resid)
    sst = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - sse / sst if sst > 0 else 0.0
    r2_adj = 1 - (1 - r2) * (n - 1) / (n - p - 1) if n > p + 1 else r2
    sigma2 = sse / max(1, n - p - 1)
    se = np.sqrt(np.diag(XtX_inv) * sigma2)
    t_stat = beta / np.maximum(se, 1e-12)
    return beta, se, t_stat, r2, r2_adj


def vif(X):
    """각 feature 의 VIF (분산 팽창 계수)."""
    n, p = X.shape
    out = []
    for i in range(p):
        y = X[:, i]
        Xi = np.delete(X, i, axis=1)
        try:
            _, _, _, r2, _ = ols(Xi, y)
            out.append(1 / max(1e-12, 1 - r2))
        except Exception:
            out.append(float('nan'))
    return out


def main():
    lines = []
    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    log("=" * 72)
    log("  difficulty_regression — F₂ 난이도 다중회귀 (N4)")
    log("=" * 72)

    # ------------------------------------------------------------------
    # PART A: 전역 회귀
    # ------------------------------------------------------------------
    log("\n[A] 전역 회귀 (51 밴드, 전 영점)")
    rows = collect_global()
    log(f"    데이터: N={len(rows)}")
    feat_names = ["log_t", "gap_prev", "gap_next", "gap_min", "gap_asym",
                  "g_prev_norm", "g_next_norm", "g_min_norm"]
    X = np.array([[r[f] for f in feat_names] for r in rows])
    y = np.array([r["abs_F2"] for r in rows])

    log(f"    label: |final_F2|  mean={y.mean():.4e}  std={y.std():.4e}  max={y.max():.4e}")

    # 표준화 후 OLS
    Xz, mu, sd = standardize(X)
    beta, se, tstat, r2, r2_adj = ols(Xz, y)

    log("\n    표준화 회귀계수 (β, t-stat):")
    log(f"      {'feature':>14s} {'β':>12s} {'se':>10s} {'t':>8s} {'|t|>3?':>8s}")
    log(f"      {'(intercept)':>14s} {beta[0]:>12.4e} {se[0]:>10.4e} {tstat[0]:>8.2f}")
    for i, fn in enumerate(feat_names):
        mark = "***" if abs(tstat[i+1]) > 3 else "**" if abs(tstat[i+1]) > 2 else ""
        log(f"      {fn:>14s} {beta[i+1]:>12.4e} {se[i+1]:>10.4e} "
            f"{tstat[i+1]:>8.2f} {mark:>8s}")
    log(f"\n    R² = {r2:.4f}   adj R² = {r2_adj:.4f}")

    # VIF
    vifs = vif(Xz)
    log("\n    VIF (다중공선성, >10 경고):")
    for fn, v in zip(feat_names, vifs):
        mark = " ⚠" if v > 10 else ""
        log(f"      {fn:>14s}  VIF = {v:>8.2f}{mark}")

    # 단순 상관 (각 feature 단독)
    log("\n    단변량 상관 (표준화 y vs 각 X 열):")
    yz = (y - y.mean()) / (y.std() + 1e-12)
    for i, fn in enumerate(feat_names):
        r = float((Xz[:, i] * yz).mean())
        log(f"      {fn:>14s}  corr = {r:+.4f}")

    # 축소 회귀: 정규화 변수만 (다중공선성 제거 후보)
    log("\n    축소 모델 A1 (log_t + g_min_norm 만):")
    keep = [0, 7]  # log_t, g_min_norm
    Xz_small = Xz[:, keep]
    beta_s, se_s, t_s, r2_s, _ = ols(Xz_small, y)
    log(f"      β[log_t]      = {beta_s[1]:+.4e}  t={t_s[1]:.2f}")
    log(f"      β[g_min_norm] = {beta_s[2]:+.4e}  t={t_s[2]:.2f}")
    log(f"      R² = {r2_s:.4f}")

    log("\n    축소 모델 A2 (log_t 만):")
    beta_s2, se_s2, t_s2, r2_s2, _ = ols(Xz[:, [0]], y)
    log(f"      β[log_t] = {beta_s2[1]:+.4e}  t={t_s2[1]:.2f}")
    log(f"      R² = {r2_s2:.4f}")

    global_result = {
        "n": len(rows),
        "features": feat_names,
        "beta": beta.tolist(),
        "se": se.tolist(),
        "t_stat": tstat.tolist(),
        "r2": r2, "r2_adj": r2_adj,
        "vif": vifs,
        "reduced_logt_gmin_r2": r2_s,
        "reduced_logt_only_r2": r2_s2,
    }

    # ------------------------------------------------------------------
    # PART B: 풍부 feature 회귀 (t∈[100,200])
    # ------------------------------------------------------------------
    log("\n" + "=" * 72)
    log("[B] 풍부 회귀 (t∈[100,200], 50 영점)")

    try:
        with open(ANALYSIS / "angle_at_zeros_inventory.json") as f:
            inv = json.load(f)
        with open(ANALYSIS / "phase_isolation_experiment.json") as f:
            iso = json.load(f)
        with open(OVERNIGHT / "result_t100-200.json") as f:
            r100 = json.load(f)
        zeros_b = np.array(iso["zeros"])
        f2_b = np.array(r100["final_F2"])
        cstd = np.array(inv["D1_seed_consistency"]["circular_std_per_zero"])
        xi_p = np.array(inv["D2_xi_prime_sign"]["xi_prime_at_zero"])
        log(f"    영점 수: {len(zeros_b)}  (cstd len={len(cstd)}, xi_p len={len(xi_p)})")

        # peak_rel 은 xi_cache 에서 다시 계산
        import torch
        cache = torch.load(OVERNIGHT / "xi_cache_t100-200_n500.pt",
                           weights_only=True)
        t_arr = cache["t"].numpy()
        xr = cache["xi_real"].numpy()

        def peak_between(t_lo, t_hi):
            m = (t_arr > t_lo) & (t_arr < t_hi)
            return float(np.abs(xr[m]).max()) if m.any() else 0.0

        peaks = np.array([peak_between(zeros_b[i], zeros_b[i+1])
                          for i in range(len(zeros_b) - 1)])

        def peak_rel(i):
            if i < 0 or i >= len(peaks):
                return float('nan')
            lo, hi = max(0, i - 2), min(len(peaks), i + 3)
            nb = np.delete(peaks[lo:hi], i - lo)
            m = nb.mean() if len(nb) > 0 else 0
            return float(peaks[i] / m) if m > 0 else float('nan')

        # 각 영점 k (1..48) feature
        rows_b = []
        gaps_b = np.diff(zeros_b)
        for k in range(1, len(zeros_b) - 1):
            rel_left = peak_rel(k - 1)
            rel_right = peak_rel(k)
            rel_min = min(r for r in (rel_left, rel_right) if not math.isnan(r))
            rows_b.append({
                "t": float(zeros_b[k]),
                "log_t": float(math.log(zeros_b[k])),
                "gap_min": float(min(gaps_b[k-1], gaps_b[k])),
                "cstd": float(cstd[k]),
                "abs_xi_prime": float(abs(xi_p[k])),
                "peak_rel_min": float(rel_min),
                "abs_F2": float(abs(f2_b[k])),
            })

        fnB = ["log_t", "gap_min", "cstd", "abs_xi_prime", "peak_rel_min"]
        Xb = np.array([[r[f] for f in fnB] for r in rows_b])
        yb = np.array([r["abs_F2"] for r in rows_b])
        Xbz, _, _ = standardize(Xb)
        beta_b, se_b, tb, r2b, r2b_adj = ols(Xbz, yb)

        log(f"    N={len(rows_b)}   label |F₂| mean={yb.mean():.4f}")
        log(f"\n    {'feature':>14s} {'β':>12s} {'t':>8s}")
        log(f"    {'(intercept)':>14s} {beta_b[0]:>12.4e} {tb[0]:>8.2f}")
        for i, fn in enumerate(fnB):
            mark = "***" if abs(tb[i+1]) > 3 else "**" if abs(tb[i+1]) > 2 else ""
            log(f"    {fn:>14s} {beta_b[i+1]:>12.4e} {tb[i+1]:>8.2f} {mark}")
        log(f"\n    R² = {r2b:.4f}  adj R² = {r2b_adj:.4f}")

        # 단변량 상관
        log("\n    단변량 상관:")
        ybz = (yb - yb.mean()) / (yb.std() + 1e-12)
        for i, fn in enumerate(fnB):
            r = float((Xbz[:, i] * ybz).mean())
            log(f"      {fn:>14s}  corr = {r:+.4f}")

        # 잔차 구조: 상위 5 잔차 영점
        X1 = np.column_stack([np.ones(len(rows_b)), Xbz])
        yhat = X1 @ beta_b
        resid = yb - yhat
        top_resid_idx = np.argsort(-np.abs(resid))[:5]
        log("\n    최대 잔차 5개 영점 (모델이 예측 못한 난이도):")
        for i in top_resid_idx:
            r = rows_b[i]
            log(f"      t={r['t']:.4f}  |F₂|={r['abs_F2']:.4f}  "
                f"pred={yhat[i]:.4f}  resid={resid[i]:+.4f}  "
                f"cstd={r['cstd']:.3f}  peak_rel={r['peak_rel_min']:.3f}")

        rich_result = {
            "n": len(rows_b),
            "features": fnB,
            "beta": beta_b.tolist(),
            "t_stat": tb.tolist(),
            "r2": r2b, "r2_adj": r2b_adj,
            "residual_top5_idx": [int(i) for i in top_resid_idx],
            "residual_top5_rows": [rows_b[int(i)] for i in top_resid_idx],
        }
    except Exception as e:
        log(f"    ⚠ 풍부 회귀 실패: {e}")
        import traceback; traceback.print_exc()
        rich_result = {"error": str(e)}

    # ------------------------------------------------------------------
    # 저장
    # ------------------------------------------------------------------
    out = {
        "global": global_result,
        "rich_t100_200": rich_result,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")
    log(f"저장: {OUT_TXT}")


if __name__ == "__main__":
    main()
