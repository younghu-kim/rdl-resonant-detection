"""
독립 검증: 결과 #36b 핵심 수치 재현
검증 항목:
  1. 영점 개수 (4901개)
  2. t>600 필터 후 영점 수 (4659개)
  3. r̄ (0.6161)
  4. E[r̃]_GUE = 0.6027 (수치적분)
  5. KS p-값 (scipy.stats.kstest)
"""

import numpy as np
from scipy import stats

CACHE = "/home/k0who029/Desktop/gdl_unified/outputs/cache/spacing_ratio_zeros_n100_5000.npy"

# ── 1. 영점 로드 ──
t_zeros = np.load(CACHE)
print(f"[1] 영점 개수: {len(t_zeros)}  (기대: 4901, {'OK' if len(t_zeros)==4901 else 'FAIL'})")
print(f"    t 범위: {t_zeros[0]:.2f} ~ {t_zeros[-1]:.2f}")

# ── 2. t>600 필터 ──
t600 = t_zeros[t_zeros > 600]
print(f"[2] t>600 영점 수: {len(t600)}  (기대: 4659, {'OK' if len(t600)==4659 else 'FAIL'})")

# ── 3. spacing ratio 계산 ──
def compute_ratios(t_sub):
    spacings = np.diff(t_sub)
    r_list = []
    for i in range(len(spacings) - 1):
        s1, s2 = spacings[i], spacings[i+1]
        if s1 > 0 and s2 > 0:
            r_list.append(min(s1, s2) / max(s1, s2))
    return np.array(r_list)

r600 = compute_ratios(t600)
r_mean = r600.mean()
print(f"[3] r̄ (t>600): {r_mean:.4f}  (기대: 0.6161, {'OK' if abs(r_mean-0.6161)<0.0005 else 'FAIL'})")
print(f"    N_r: {len(r600)}  (기대: 4657)")

# ── 4. E[r̃]_GUE 수치적분 ──
def p_gue(r):
    C = 81.0 * np.sqrt(3) / (4.0 * np.pi)
    return C * (r + r**2)**2 / (1.0 + r + r**2)**4

x = np.linspace(0, 1, 1_000_001)
# p_gue의 [0,1] 적분 ≈ 0.5 → r̃ 정의역에서 밀도는 2*p_gue
E_gue = np.trapz(2.0 * x * p_gue(x), x)
print(f"[4] E[r̃]_GUE (수치적분): {E_gue:.6f}  (기대: 0.6027, {'OK' if abs(E_gue-0.6027)<0.0001 else 'FAIL'})")

# ── 5. KS 검정 (GUE CDF 직접 구성) ──
def build_cdf(pdf_func, n=50000):
    r_g = np.linspace(0, 1, n+1)
    p_v = pdf_func(r_g)
    dr = 1.0 / n
    cdf = np.zeros(n+1)
    for i in range(1, n+1):
        cdf[i] = cdf[i-1] + 0.5*(p_v[i-1]+p_v[i])*dr
    cdf /= cdf[-1]
    return r_g, cdf

r_grid, cdf_gue = build_cdf(p_gue)
D, p_val = stats.kstest(r600, lambda x: np.interp(x, r_grid, cdf_gue))
print(f"[5] KS 검정 (t>600 vs GUE):")
print(f"    D = {D:.5f}  (기대: 0.02993)")
print(f"    p = {p_val:.4e}  (기대: 4.66e-04)")
p_match = abs(p_val - 4.6634e-4) / 4.6634e-4 < 0.10  # 10% 이내
print(f"    p-값 재현: {'OK' if p_match else 'FAIL'} (10% 허용오차)")

# ── 요약 ──
print()
print("=== 독립 검증 요약 ===")
checks = [
    ("영점 개수 4901", len(t_zeros)==4901),
    ("t>600 영점 4659", len(t600)==4659),
    ("r̄ 0.6161", abs(r_mean-0.6161)<0.0005),
    ("E[r̃]_GUE 0.6027", abs(E_gue-0.6027)<0.0001),
    ("KS p 재현", p_match),
]
for name, ok in checks:
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")

all_pass = all(ok for _, ok in checks)
print(f"\n최종: {'전체 PASS' if all_pass else '일부 FAIL'}")
