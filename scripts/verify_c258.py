"""
verify_c258.py — C-258 B-42 편상관 분석 독립 수치 검증
사이클 #258 결과의 핵심 수치를 C-256 원데이터에서 재현한다.
"""

import re
import sys
import math

# mpmath 가용성 확인
try:
    import mpmath
    MPMATH_OK = True
except ImportError:
    MPMATH_OK = False
    print("[경고] mpmath 미설치 — κ_mid 3점 계산 건너뜀")

try:
    from scipy import stats as spstats
    SCIPY_OK = True
except ImportError:
    SCIPY_OK = False
    print("[경고] scipy 미설치 — Spearman 수동 계산으로 대체")

DATA_PATH = "/home/k0who029/Desktop/gdl_unified/results/A_gap_correlation_c256.txt"

# ──────────────────────────────────────────
# 1. C-256 데이터 파싱 (처음 20개 영점)
# ──────────────────────────────────────────
print("=" * 70)
print("§1  C-256 데이터 파싱 (처음 20개 영점)")
print("=" * 70)

rows = []
with open(DATA_PATH, "r") as fh:
    for line in fh:
        line = line.strip()
        # 숫자로 시작하는 데이터 행 파싱
        m = re.match(
            r"^\s*(\d+)\s+([\d.]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([\d.e+]+)\s+([\d.e+nan]+)",
            line,
        )
        if m:
            idx   = int(m.group(1))
            t     = float(m.group(2))
            re_c0 = float(m.group(3))
            im_c0 = float(m.group(4))
            re_c1 = float(m.group(5))
            im_c1 = float(m.group(6))
            A     = float(m.group(7))
            gr_str = m.group(8)
            gap_right = float("nan") if gr_str == "nan" else float(gr_str)
            rows.append(dict(idx=idx, t=t, re_c0=re_c0, im_c0=im_c0,
                             re_c1=re_c1, im_c1=im_c1, A=A, gap_right=gap_right))
        if len(rows) >= 20:
            break

print(f"파싱된 영점 수: {len(rows)}")
print(f"{'#':>3}  {'t':>10}  {'im_c0':>10}  {'re_c1':>10}  {'A':>10}  {'gap_right':>10}")
for r in rows:
    print(f"{r['idx']:>3}  {r['t']:>10.4f}  {r['im_c0']:>10.6f}  "
          f"{r['re_c1']:>10.6f}  {r['A']:>10.6f}  {r['gap_right']:>10.4f}")

# ──────────────────────────────────────────
# 2. gap_left 계산 (n≥1, 0-indexed)
# ──────────────────────────────────────────
print()
print("=" * 70)
print("§2  gap_left = t[n] - t[n-1] 계산")
print("=" * 70)

for i in range(1, len(rows)):
    rows[i]["gap_left"] = rows[i]["t"] - rows[i - 1]["t"]
rows[0]["gap_left"] = float("nan")

print(f"{'#':>3}  {'gap_left':>10}  {'gap_right':>10}")
for r in rows:
    gl = rows.index(r)  # 실제로 enumerate 쓰자
for i, r in enumerate(rows):
    gl = r.get("gap_left", float("nan"))
    print(f"{r['idx']:>3}  {gl:>10.4f}  {r['gap_right']:>10.4f}")

# ──────────────────────────────────────────
# 3. Im(c₀)_NN = 1/gap_right - 1/gap_left 계산
# ──────────────────────────────────────────
print()
print("=" * 70)
print("§3  Im(c₀)_NN = 1/gap_right - 1/gap_left  (n≥2, 1-indexed)")
print("=" * 70)

# n≥2 (1-indexed) = i≥1 (0-indexed), gap_right 유효
for i, r in enumerate(rows):
    if i >= 1 and not math.isnan(r.get("gap_left", float("nan"))) and not math.isnan(r["gap_right"]):
        r["im_c0_NN"] = 1.0 / r["gap_right"] - 1.0 / r["gap_left"]
    else:
        r["im_c0_NN"] = float("nan")

print(f"{'#':>3}  {'Im(c0)_NN':>12}  {'im_c0 (원본)':>14}  {'im_c0_bg':>12}")
for r in rows:
    nn = r.get("im_c0_NN", float("nan"))
    if not math.isnan(nn):
        bg = r["im_c0"] - nn
        print(f"{r['idx']:>3}  {nn:>12.6f}  {r['im_c0']:>14.6f}  {bg:>12.6f}")
    else:
        print(f"{r['idx']:>3}  {'nan':>12}  {r['im_c0']:>14.6f}  {'nan':>12}")

# ──────────────────────────────────────────
# 4. A_NN, A_bg, A_cross 분해 및 검증
# ──────────────────────────────────────────
print()
print("=" * 70)
print("§4  A 분해: A_NN + A_bg + A_cross ≈ A 검증")
print("=" * 70)

valid = []
for r in rows:
    nn = r.get("im_c0_NN", float("nan"))
    if math.isnan(nn):
        continue
    im_c0_bg = r["im_c0"] - nn
    A_NN    = nn ** 2
    A_bg    = im_c0_bg ** 2 + 2.0 * r["re_c1"]
    A_cross = 2.0 * nn * im_c0_bg
    A_recon = A_NN + A_bg + A_cross
    residual = abs(A_recon - r["A"])
    valid.append(dict(idx=r["idx"], A=r["A"], A_NN=A_NN, A_bg=A_bg,
                      A_cross=A_cross, A_recon=A_recon, residual=residual))
    r["A_NN"] = A_NN
    r["A_bg"] = A_bg
    r["A_cross"] = A_cross

print(f"{'#':>3}  {'A_orig':>9}  {'A_NN':>9}  {'A_bg':>9}  {'A_cross':>9}  "
      f"{'A_recon':>9}  {'잔차':>10}")
for v in valid:
    print(f"{v['idx']:>3}  {v['A']:>9.6f}  {v['A_NN']:>9.6f}  {v['A_bg']:>9.6f}  "
          f"{v['A_cross']:>9.6f}  {v['A_recon']:>9.6f}  {v['residual']:>10.2e}")

rms_resid = math.sqrt(sum(v["residual"]**2 for v in valid) / len(valid))
print(f"\nRMS 잔차: {rms_resid:.4e}  (C-258 보고값: 9.31e-07)")

A_NN_vals    = [v["A_NN"]    for v in valid]
A_bg_vals    = [v["A_bg"]    for v in valid]
A_cross_vals = [v["A_cross"] for v in valid]
A_orig_vals  = [v["A"]       for v in valid]

total_mean = (sum(A_NN_vals) + sum(A_bg_vals) + sum(A_cross_vals)) / len(valid)
frac_NN    = sum(A_NN_vals)    / len(valid) / total_mean * 100
frac_bg    = sum(A_bg_vals)    / len(valid) / total_mean * 100
frac_cross = sum(A_cross_vals) / len(valid) / total_mean * 100

print(f"\n[20개 부분집합 분산 기여율 (평균 기반)]")
print(f"  A_NN:    {frac_NN:.1f}%   (전체 C-258 보고값: 2.6%)")
print(f"  A_bg:    {frac_bg:.1f}%   (전체 C-258 보고값: 50.8%)")
print(f"  A_cross: {frac_cross:.1f}%   (전체 C-258 보고값: 10.9%)")

# 절대 기여율 재계산 (A_orig 평균 대비)
A_orig_mean = sum(A_orig_vals) / len(A_orig_vals)
frac_NN_abs    = sum(A_NN_vals)    / len(valid) / A_orig_mean * 100
frac_bg_abs    = sum(A_bg_vals)    / len(valid) / A_orig_mean * 100
frac_cross_abs = sum(A_cross_vals) / len(valid) / A_orig_mean * 100
print(f"\n[A_orig 평균({A_orig_mean:.4f}) 대비 기여율]")
print(f"  A_NN:    {frac_NN_abs:.1f}%")
print(f"  A_bg:    {frac_bg_abs:.1f}%")
print(f"  A_cross: {frac_cross_abs:.1f}%")

# ──────────────────────────────────────────
# 5. Spearman 편상관 방향성 확인 (20개 부분집합)
# ──────────────────────────────────────────
print()
print("=" * 70)
print("§5  편상관 방향성 확인 (20개 부분집합)")
print("=" * 70)

# κ_mid 계산에는 전체 rows 필요 — 20개 안에서 인접 쌍 구성
# κ_mid 대신 gap 기반 근사: κ_mid ≈ 1/(gap²) 로 단순화하지 않고
# 실제 kappa_mid 없이 ρ(A_bg_n, A_bg_n+1) 만 계산

# A_bg 인접 쌍 상관
A_bg_n  = [v["A_bg"] for v in valid[:-1]]
A_bg_n1 = [v["A_bg"] for v in valid[1:]]

def spearman_r(x, y):
    """수동 Spearman rank 상관계수"""
    n = len(x)
    def rank(lst):
        sorted_pairs = sorted(enumerate(lst), key=lambda t: t[1])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j < n - 1 and sorted_pairs[j+1][1] == sorted_pairs[j][1]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[sorted_pairs[k][0]] = avg_rank
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    mean_rx = sum(rx) / n
    mean_ry = sum(ry) / n
    num = sum((rx[i] - mean_rx) * (ry[i] - mean_ry) for i in range(n))
    den = math.sqrt(sum((rx[i] - mean_rx)**2 for i in range(n)) *
                    sum((ry[i] - mean_ry)**2 for i in range(n)))
    return num / den if den != 0 else 0.0

if SCIPY_OK:
    rho_bg_adj, p_bg_adj = spstats.spearmanr(A_bg_n, A_bg_n1)
    print(f"ρ(A_bg_n, A_bg_n+1) [scipy] = {rho_bg_adj:+.6f}  (p={p_bg_adj:.3e})")
else:
    rho_bg_adj = spearman_r(A_bg_n, A_bg_n1)
    print(f"ρ(A_bg_n, A_bg_n+1) [수동]  = {rho_bg_adj:+.6f}")

print(f"방향 확인: {'양수 ✓' if rho_bg_adj > 0 else '음수 ✗'}  "
      f"(C-258 보고값: +0.580143)")

# A_orig 인접 상관
A_n  = [v["A"] for v in valid[:-1]]
A_n1 = [v["A"] for v in valid[1:]]
if SCIPY_OK:
    rho_A_adj, p_A_adj = spstats.spearmanr(A_n, A_n1)
    print(f"ρ(A_n, A_n+1) 원시 [scipy]  = {rho_A_adj:+.6f}  (p={p_A_adj:.3e})")
else:
    rho_A_adj = spearman_r(A_n, A_n1)
    print(f"ρ(A_n, A_n+1) 원시 [수동]   = {rho_A_adj:+.6f}")
print(f"방향 확인: {'양수 ✓' if rho_A_adj > 0 else '음수 ✗'}  "
      f"(C-256 보고값: +0.400891)")

# kappa_mid 근사 (gap 중간값 기반, mpmath 없이)
# κ_mid 대리변수: A_avg = (A[i] + A[i+1])/2
A_avg_subset = [(valid[i]["A"] + valid[i+1]["A"]) / 2.0 for i in range(len(valid) - 1)]
# kappa_mid 대리: 1/gap_avg² 로 근사 — 정확 계산은 §6에서
kappa_proxy = []
for i in range(len(valid) - 1):
    r_i   = rows[valid[i]["idx"] - 1]    # 1-indexed → 0-indexed
    r_i1  = rows[valid[i+1]["idx"] - 1]
    if not math.isnan(r_i["gap_right"]):
        kappa_proxy.append(1.0 / r_i["gap_right"]**2)
    else:
        kappa_proxy.append(float("nan"))

# nan 제거
pairs_kA = [(k, a) for k, a in zip(kappa_proxy, A_avg_subset) if not math.isnan(k)]
if pairs_kA:
    kp, aa = zip(*pairs_kA)
    if SCIPY_OK:
        rho_kA, p_kA = spstats.spearmanr(kp, aa)
        print(f"\nρ(κ_proxy, A_avg) [scipy]   = {rho_kA:+.6f}  (p={p_kA:.3e})")
    else:
        rho_kA = spearman_r(list(kp), list(aa))
        print(f"\nρ(κ_proxy, A_avg) [수동]    = {rho_kA:+.6f}")
    print(f"방향 확인: {'양수 ✓' if rho_kA > 0 else '음수 ✗'}  "
          f"(C-258 보고값: +0.350739 ~ +0.436167)")

# ──────────────────────────────────────────
# 6. κ_mid 3점 독립 계산 (mpmath, dps=40)
# ──────────────────────────────────────────
print()
print("=" * 70)
print("§6  κ_mid 3점 독립 계산 (mpmath dps=40)")
print("=" * 70)

if not MPMATH_OK:
    print("[건너뜀] mpmath 미설치")
else:
    mpmath.mp.dps = 40

    def kappa_at(t_val):
        """
        κ = |ξ'/ξ(s)|² 에서  ξ'/ξ = (log ξ)'
        log ξ(s) = log ξ(1-s)  (함수방정식)
        ξ'(s)/ξ(s) = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2 + ζ'(s)/ζ(s)
        ψ = digamma
        """
        s = mpmath.mpc("0.5", str(t_val))
        # 각 항 계산
        term_s    = 1 / s
        term_sm1  = 1 / (s - 1)
        term_logpi = -mpmath.log(mpmath.pi) / 2
        term_psi  = mpmath.digamma(s / 2) / 2
        # ζ'/ζ: mpmath.zeta 는 수치적으로 불안정할 수 있으나
        # 크리티컬 라인에서 zeta 가 비영점이면 정상 계산
        try:
            zv  = mpmath.zeta(s)
            zdv = mpmath.diff(mpmath.zeta, s)
            term_zeta = zdv / zv
        except (ZeroDivisionError, ValueError):
            return float("nan")
        xi_log_deriv = term_s + term_sm1 + term_logpi + term_psi + term_zeta
        kappa = abs(xi_log_deriv) ** 2
        return float(kappa)

    # 3개 중간점: (t[1]+t[2])/2, (t[5]+t[6])/2, (t[10]+t[11])/2  (0-indexed)
    t_pts = [
        (rows[1]["t"] + rows[2]["t"]) / 2.0,   # 쌍 (2,3)
        (rows[5]["t"] + rows[6]["t"]) / 2.0,   # 쌍 (6,7)
        (rows[10]["t"] + rows[11]["t"]) / 2.0, # 쌍 (11,12)
    ]
    labels = ["(t[2]+t[3])/2", "(t[6]+t[7])/2", "(t[11]+t[12])/2"]

    print(f"{'중간점':>18}  {'t_mid':>10}  {'κ_mid':>14}")
    for lbl, t_mid in zip(labels, t_pts):
        kv = kappa_at(t_mid)
        print(f"{lbl:>18}  {t_mid:>10.4f}  {kv:>14.6f}")

    # 비교: C-256 κ_mid (ρ에서 mean=0.6502, std=0.3731)
    print(f"\n참고 — C-256/258 전체 κ_mid 통계: mean≈0.6502, std≈0.3731")
    print(f"3점이 합리적 범위(~0.3~1.5)에 있으면 ✓")

# ──────────────────────────────────────────
# 최종 요약
# ──────────────────────────────────────────
print()
print("=" * 70)
print("종합 검증 요약")
print("=" * 70)
print(f"[1] RMS 잔차 |A_recon - A|:  {rms_resid:.4e}  (보고값 9.31e-07, 20개 부분집합)")
print(f"    → 재구성 정밀도: {'✓ 양호' if rms_resid < 1e-4 else '주의'}")
print(f"[2] A_bg/A_orig 평균 기여율: {frac_bg_abs:.1f}%  (전체 보고값 50.8%)")
print(f"    → 방향 일치: {'✓' if frac_bg_abs > 30 else '주의'}")
print(f"[3] A_NN/A_orig 평균 기여율: {frac_NN_abs:.1f}%  (전체 보고값 2.6%)")
print(f"    → 방향 일치: {'✓' if frac_NN_abs < 20 else '주의'}")
print(f"[4] ρ(A_bg_n, A_bg_n+1):    {rho_bg_adj:+.4f}  (보고값 +0.5801)")
print(f"    → 양수 방향: {'✓' if rho_bg_adj > 0 else '✗ 불일치'}")
print(f"[5] ρ(A_n, A_n+1) 원시:      {rho_A_adj:+.4f}  (보고값 +0.4009)")
print(f"    → 양수 방향: {'✓' if rho_A_adj > 0 else '✗ 불일치'}")
if pairs_kA:
    print(f"[6] ρ(κ_proxy, A_avg):       {rho_kA:+.4f}  (보고값 +0.35~+0.44)")
    print(f"    → 양수 방향: {'✓' if rho_kA > 0 else '✗ 불일치'}")
