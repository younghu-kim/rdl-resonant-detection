#!/usr/bin/env python3
"""
lehmer_mpmath_scan.py — 고-t 영점에 대한 mpmath 기반 Lehmer 서명 검증 (N1c+).

문제: xi_cache_*.pt 파일은 밴드당 n=500 고정 샘플 → t>500 에서는 dt > 이웃 영점
간격이 되어 저-t 스캔(lehmer_global_scan.py)이 4개 밴드만 처리 가능했음.

해결: mpmath 로 고-t 영점 근방을 직접 조회하여, **정규화 간격**
      g_n = (γ_{n+1}-γ_n) * log(γ_n/(2π)) / (2π)  [평균 = 1]
      을 Lehmer 지표로 사용. g 가 작을수록 Lehmer-유사 압축.

프로토콜:
  for each result_t*.json 밴드 B:
    hardest = 상위 3 F₂-hardest 영점 t_h
    for t_h:
      n = nzeros(t_h) 로 인덱스 추정
      γ_{n-1}, γ_n, γ_{n+1} 조회 → 2 gap → min g 기록
    control = 해당 밴드 내 랜덤 3 영점 → 같은 방식
  hardest 의 min_g 분포 vs control 분포 비교.

출력: outputs/analysis/lehmer_mpmath_scan.{json,txt}
"""
import json
import math
import random
import re
import time
from pathlib import Path

import mpmath

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "outputs" / "analysis"
OVERNIGHT = ROOT / "outputs" / "overnight"
OUT_JSON = ANALYSIS / "lehmer_mpmath_scan.json"
OUT_TXT = ANALYSIS / "lehmer_mpmath_scan.txt"

mpmath.mp.dps = 15
random.seed(42)

HARDEST_PER_BAND = 3
CONTROL_PER_BAND = 3


def normalized_gap(gamma_left, gamma_right):
    """Riemann-von Mangoldt 정규화 간격. 평균 = 1."""
    dg = gamma_right - gamma_left
    t_mid = (gamma_left + gamma_right) / 2
    density = math.log(t_mid / (2 * math.pi)) / (2 * math.pi)  # 평균 영점 밀도
    return dg * density


def parse_band(name):
    """'result_t10000-10100.json' -> (10000.0, 10100.0)"""
    m = re.match(r"result_t([\d\.]+)-([\d\.]+)\.json", name)
    if not m:
        return None
    return float(m.group(1)), float(m.group(2))


def nearest_n(t_target):
    """t 근처의 Riemann 영점 인덱스 추정."""
    return int(mpmath.nzeros(t_target))


def triplet_gaps(n_center, cache):
    """γ_{n-1}, γ_n, γ_{n+1} 에서 두 간격 계산. 캐시 재사용."""
    out = []
    for k in (n_center - 1, n_center, n_center + 1):
        if k not in cache:
            cache[k] = float(mpmath.zetazero(k).imag)
        out.append(cache[k])
    g_left = normalized_gap(out[0], out[1])
    g_right = normalized_gap(out[1], out[2])
    return out[1], g_left, g_right


def main():
    lines = []
    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    log("=" * 72)
    log("  lehmer_mpmath_scan — 고-t Lehmer 서명 (정규화 간격)")
    log("=" * 72)

    result_files = sorted(OVERNIGHT.glob("result_t*.json"))
    log(f"  밴드 수: {len(result_files)}")
    log(f"  밴드당 hardest: {HARDEST_PER_BAND},  control: {CONTROL_PER_BAND}")
    log(f"  정규화 간격 g: 평균=1, g<0.5 = Lehmer-유사")

    hardest_records = []
    control_records = []
    band_summaries = []

    t_start = time.time()
    for bi, rf in enumerate(result_files):
        band = parse_band(rf.name)
        if band is None:
            continue
        t_lo, t_hi = band
        with open(rf) as f:
            rd = json.load(f)
        hz = rd.get("hardest_zeros", [])
        hardest_t = [float(x[0]) for x in hz[:HARDEST_PER_BAND]
                     if isinstance(x, list) and len(x) >= 1]
        if not hardest_t:
            continue

        band_cache = {}  # n -> γ_n
        band_hard = []
        band_ctrl = []

        # 1. hardest 영점 분석
        for th in hardest_t:
            try:
                n = nearest_n(th)
                gamma, gL, gR = triplet_gaps(n, band_cache)
                # 실제 th 가 γ_n 과 맞는지 확인 (틀리면 한 단계 조정)
                if abs(gamma - th) > 1.0:
                    # 인덱스 보정 시도
                    for dn in (-1, 1, -2, 2):
                        if n + dn not in band_cache:
                            band_cache[n + dn] = float(mpmath.zetazero(n + dn).imag)
                        if abs(band_cache[n + dn] - th) < abs(gamma - th):
                            n = n + dn
                            gamma, gL, gR = triplet_gaps(n, band_cache)
                            break
                min_g = min(gL, gR)
                rec = {
                    "band": rf.name, "t_requested": th, "n": n,
                    "gamma": gamma, "match_err": abs(gamma - th),
                    "g_left": gL, "g_right": gR, "min_g": min_g,
                }
                band_hard.append(rec)
                hardest_records.append(rec)
            except Exception as e:
                log(f"  ⚠ {rf.name} t={th}: {e}")

        # 2. control: 밴드 내 랜덤 t 선택 → 근처 영점의 min_g
        # 밴드의 영점 개수를 추정한 뒤 무작위 인덱스 선택
        try:
            n_lo = int(mpmath.nzeros(t_lo))
            n_hi = int(mpmath.nzeros(t_hi))
            if n_hi - n_lo >= CONTROL_PER_BAND + 2:
                chosen_idx = random.sample(range(n_lo + 1, n_hi - 1),
                                           min(CONTROL_PER_BAND, n_hi - n_lo - 2))
                for n in chosen_idx:
                    try:
                        gamma, gL, gR = triplet_gaps(n, band_cache)
                        min_g = min(gL, gR)
                        rec = {
                            "band": rf.name, "n": n, "gamma": gamma,
                            "g_left": gL, "g_right": gR, "min_g": min_g,
                        }
                        band_ctrl.append(rec)
                        control_records.append(rec)
                    except Exception as e:
                        log(f"  ⚠ ctrl {rf.name} n={n}: {e}")
        except Exception as e:
            log(f"  ⚠ {rf.name} control sampling: {e}")

        band_summaries.append({
            "band": rf.name,
            "n_hard": len(band_hard),
            "n_ctrl": len(band_ctrl),
            "hard_min_g_mean": (sum(r["min_g"] for r in band_hard) / len(band_hard)
                                if band_hard else None),
            "ctrl_min_g_mean": (sum(r["min_g"] for r in band_ctrl) / len(band_ctrl)
                                if band_ctrl else None),
        })

        elapsed = time.time() - t_start
        log(f"  [{bi+1:>2}/{len(result_files)}] {rf.name:<30} "
            f"hard_min_g={band_summaries[-1]['hard_min_g_mean']}  "
            f"ctrl_min_g={band_summaries[-1]['ctrl_min_g_mean']}  "
            f"({elapsed:.1f}s)")

    # ------------------------------------------------------------------
    # 전역 통계
    # ------------------------------------------------------------------
    log("\n" + "=" * 72)
    log("  전역 집계")
    log("=" * 72)

    def stats(records, label):
        if not records:
            log(f"  {label}: (비어있음)")
            return None
        gs = [r["min_g"] for r in records]
        gs_sorted = sorted(gs)
        mean = sum(gs) / len(gs)
        median = gs_sorted[len(gs) // 2]
        frac_lt_05 = sum(1 for g in gs if g < 0.5) / len(gs)
        frac_lt_03 = sum(1 for g in gs if g < 0.3) / len(gs)
        log(f"  {label}:  N={len(gs)}  mean={mean:.4f}  median={median:.4f}")
        log(f"    P(min_g<0.5)={frac_lt_05*100:.1f}%   "
            f"P(min_g<0.3)={frac_lt_03*100:.1f}%")
        return {"n": len(gs), "mean": mean, "median": median,
                "frac_lt_05": frac_lt_05, "frac_lt_03": frac_lt_03,
                "values": gs}

    log("")
    h_stats = stats(hardest_records, "F₂-hardest 영점")
    c_stats = stats(control_records, "랜덤 대조군")

    if h_stats and c_stats:
        log("")
        log(f"  hardest/control mean 비율: "
            f"{h_stats['mean']/c_stats['mean']:.3f}  (<1 이면 압축)")
        # 간단한 Welch-유사 z
        import statistics
        try:
            sh = statistics.stdev(h_stats["values"])
            sc = statistics.stdev(c_stats["values"])
            se = math.sqrt(sh**2 / h_stats["n"] + sc**2 / c_stats["n"])
            if se > 0:
                z = (h_stats["mean"] - c_stats["mean"]) / se
                log(f"  Welch z = {z:+.2f}  "
                    f"({'hardest 가 유의하게 작음' if z < -2 else 'hardest 가 작음' if z < -1 else '유의하지 않음'})")
        except Exception:
            pass

        # 저-압축(min_g<0.5) 비율 오즈비
        if c_stats["frac_lt_05"] > 0:
            odds = h_stats["frac_lt_05"] / c_stats["frac_lt_05"]
            log(f"  P(min_g<0.5) 비율 = {odds:.2f}×  "
                f"(hardest 가 Lehmer-유사일 오즈)")

    # 저장
    out = {
        "config": {
            "hardest_per_band": HARDEST_PER_BAND,
            "control_per_band": CONTROL_PER_BAND,
            "n_bands": len(result_files),
            "mpmath_dps": mpmath.mp.dps,
        },
        "hardest_records": hardest_records,
        "control_records": control_records,
        "band_summaries": band_summaries,
        "global": {
            "hardest_stats": {k: v for k, v in (h_stats or {}).items()
                              if k != "values"},
            "control_stats": {k: v for k, v in (c_stats or {}).items()
                              if k != "values"},
        } if h_stats and c_stats else {},
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")
    log(f"저장: {OUT_TXT}")
    log(f"총 시간: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
