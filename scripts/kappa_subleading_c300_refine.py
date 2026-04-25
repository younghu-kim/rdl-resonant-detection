"""
[C-300 정제] κδ² Taylor 계수의 전체 Hadamard 공식 검증
초기 음성 원인: c₂ ≈ -2H₁+const 가설이 S₁² 지배항을 무시.
전체 공식: c₂ = (S₁ - β)² + α² - 2H₁  (α=Re(C₀), β=Im(C₀))
"""
import numpy as np
from scipy import stats

# C-300 데이터 직접 입력
zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
         37.586178, 40.918719, 43.327073, 48.005151, 49.773832]
c1 = [-1.1599, -1.5117, -1.1735, -1.9148, -0.9342,
      -1.6809, -1.7201, -0.9026, -2.4057, -0.8572]
c2 = [0.2150, 0.3298, 0.0556, 0.4454, -0.2971,
      0.2806, 0.1060, -0.4051, 0.6092, -0.7676]
S1 = [-1.2238, -1.2906, -1.0637, -1.3593, -0.8364,
      -1.1506, -1.1285, -0.6913, -1.3853, -0.5915]
H1 = [0.0509, 0.1240, 0.1509, 0.2646, 0.2822,
      0.2387, 0.3604, 0.3359, 0.4989, 0.5400]

c1, c2, S1, H1 = map(np.array, [c1, c2, S1, H1])
N = len(zeros)

print("=" * 60)
print("[C-300 정제] κδ² 전체 Hadamard 공식 검증")
print("=" * 60)

# ━━━ 1. c₁ = 2(S₁ - β) 검증 ━━━
print("\n[1] c₁ = 2(S₁ - β) 모델")
# c₁ = 2S₁ - 2β → c₁ = a + 2·S₁ (a = -2β)
# 기존 자유 기울기 적합: b=1.52
# 이론 기울기 = 2 고정 적합:
beta_from_c1 = (2.0 * S1 - c1).mean() / 2.0  # β = mean((2S₁-c₁)/2)
c1_pred = 2.0 * (S1 - beta_from_c1)
c1_res = c1 - c1_pred
print(f"  β 추정 (c₁에서): {beta_from_c1:.4f}")
print(f"  c₁ 잔차 σ: {c1_res.std():.4f}")
print(f"  c₁ 잔차 max|e|: {np.max(np.abs(c1_res)):.4f}")

# 자유 적합 vs 고정 기울기 비교
A_free = np.column_stack([np.ones(N), S1])
cf = np.linalg.lstsq(A_free, c1, rcond=None)[0]
c1_pred_free = A_free @ cf
rss_free = np.sum((c1 - c1_pred_free)**2)
rss_fixed = np.sum(c1_res**2)
print(f"  자유 기울기: a={cf[0]:.4f}, b={cf[1]:.4f} (RSS={rss_free:.4f})")
print(f"  고정 기울기 2: RSS={rss_fixed:.4f}")
# F-test (추가 자유도 1, n-2=8)
if rss_free > 1e-10:
    F_stat = (rss_fixed - rss_free) / (rss_free / (N-2))
    print(f"  F-test (기울기=2 제약): F={F_stat:.2f}")

# ━━━ 2. c₂ = (S₁-β)² + α² - 2H₁ 검증 ━━━
print(f"\n[2] c₂ = (S₁-β)² + α² - 2H₁ 모델")

# 모델: c₂ = S₁² - 2β·S₁ + β² + α² - 2H₁
# = (α²+β²) - 2β·S₁ + S₁² - 2H₁
# 3개 예측변수: 1, S₁, S₁²  (H₁ 계수 = -2 고정)
# c₂ + 2H₁ = a + b·S₁ + S₁²  (S₁² 계수 = 1 고정 이론)

target = c2 + 2.0 * H1  # c₂+2H₁ = α²+β² - 2β·S₁ + S₁²

# 모델 A: 자유 적합 c₂+2H₁ = a + b·S₁ + c·S₁²
A3 = np.column_stack([np.ones(N), S1, S1**2])
cf3 = np.linalg.lstsq(A3, target, rcond=None)[0]
pred3 = A3 @ cf3
ss_res3 = np.sum((target - pred3)**2)
ss_tot = np.sum((target - target.mean())**2)
r2_3 = 1.0 - ss_res3/ss_tot if ss_tot > 0 else 0
print(f"  자유 적합: a={cf3[0]:.4f}, b={cf3[1]:.4f}, c(S₁²)={cf3[2]:.4f}")
print(f"  이론 예측: c(S₁²)=1.0, b=-2β, a=α²+β²")
print(f"  R² = {r2_3:.6f}")

# β 추정 (c₂에서)
beta_from_c2 = -cf3[1] / 2.0
alpha_sq = cf3[0] - beta_from_c2**2
print(f"  β 추정 (c₂에서): {beta_from_c2:.4f}")
print(f"  α² 추정: {alpha_sq:.4f}")
print(f"  β 일관성 (c₁에서 {beta_from_c1:.4f} vs c₂에서 {beta_from_c2:.4f}): "
      f"Δβ = {abs(beta_from_c1 - beta_from_c2):.4f}")

# 모델 B: 이론 기울기 고정 c₂+2H₁ = a + b·S₁ + 1·S₁²
# = a + b·S₁ + S₁²
A2 = np.column_stack([np.ones(N), S1])
target_B = target - S1**2  # c₂+2H₁-S₁² = a + b·S₁
cf2 = np.linalg.lstsq(A2, target_B, rcond=None)[0]
pred_B = A2 @ cf2 + S1**2
ss_res_B = np.sum((target - pred_B)**2)
r2_B = 1.0 - ss_res_B/ss_tot if ss_tot > 0 else 0
print(f"\n  고정 c(S₁²)=1: a={cf2[0]:.4f}, b={cf2[1]:.4f}, R²={r2_B:.6f}")

# 모델 C: 완전 이론 (β from c₁ 사용, 자유 파라미터 = α² 만)
print(f"\n  완전 이론 (β={beta_from_c1:.4f} from c₁):")
target_C = target - S1**2 + 2*beta_from_c1*S1 - beta_from_c1**2
# = α² (상수여야)
alpha_sq_C = target_C.mean()
print(f"  α² 각 영점: {target_C}")
print(f"  평균 α² = {alpha_sq_C:.4f}, σ = {target_C.std():.4f}, "
      f"CV = {target_C.std()/abs(alpha_sq_C)*100:.1f}%")

# ━━━ 3. 전체 모델 정확도 ━━━
print(f"\n[3] 전체 예측 c₂ = (S₁-β)² + α² - 2H₁")
beta_avg = (beta_from_c1 + beta_from_c2) / 2
alpha_sq_avg = max(alpha_sq_C, 0)
c2_pred = (S1 - beta_avg)**2 + alpha_sq_avg - 2*H1
c2_err = c2 - c2_pred
print(f"  β = {beta_avg:.4f}, α² = {alpha_sq_avg:.4f}")
for i in range(N):
    print(f"  γ_{i+1}: c₂={c2[i]:+.4f}, pred={c2_pred[i]:+.4f}, err={c2_err[i]:+.4f}")
print(f"  MAE = {np.mean(np.abs(c2_err)):.4f}")
print(f"  max|err| = {np.max(np.abs(c2_err)):.4f}")
r_pred, p_pred = stats.pearsonr(c2, c2_pred)
print(f"  ρ(c₂, c₂_pred) = {r_pred:.4f}, p = {p_pred:.2e}")

# 모델 R²
ss_res_full = np.sum(c2_err**2)
ss_tot_c2 = np.sum((c2 - c2.mean())**2)
r2_full = 1.0 - ss_res_full/ss_tot_c2 if ss_tot_c2 > 0 else 0
print(f"  R² (전체 모델) = {r2_full:.4f}")

# ━━━ 4. 종합 판정 ━━━
print(f"\n{'='*60}")
print("[4] 종합 판정")

# 판정 기준
test_c1 = abs(beta_from_c1 - beta_from_c2) < 0.3  # β 일관성
test_S1sq = abs(cf3[2] - 1.0) < 0.3  # S₁² 계수 ≈ 1
test_r2 = r2_3 > 0.90  # 3변수 모델 R²
test_pred = r_pred > 0.9  # 전체 예측 상관

print(f"  β 일관성: {test_c1} (Δβ={abs(beta_from_c1-beta_from_c2):.3f})")
print(f"  S₁² 계수 ≈ 1: {test_S1sq} (={cf3[2]:.3f})")
print(f"  3변수 R²>0.9: {test_r2} (={r2_3:.4f})")
print(f"  전체 예측 ρ>0.9: {test_pred} (={r_pred:.4f})")

n_pass = sum([test_c1, test_S1sq, test_r2, test_pred])
if n_pass >= 4:
    v = "★★★★★ 강양성 — 전체 Hadamard 공식 확정"
elif n_pass >= 3:
    v = "★★★★ 양성 — Hadamard 공식 대부분 확인"
elif n_pass >= 2:
    v = "★★★ 중립 — 부분 확인, 추가 검증 필요"
else:
    v = "★★ 음성"
print(f"\n  판정: {v}")
print(f"  통과: {n_pass}/4")
