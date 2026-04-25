#!/usr/bin/env python3
"""
[C-265 보충] S₁² 기여 포함 GUE A-gap 상관

L-함수에서 A = S₁² + 2H₁. GUE 실수 고유값에서 S₁ 유사량:
  S₁ ≡ Σ_{k≠n} 1/(E_n - E_k)  (실수)
  A_full = S₁² + 2H₁

이것이 gap_right 상관을 -0.50 → -0.57로 이동시키는지 확인.
"""

import numpy as np
from scipy import stats
import time

np.random.seed(42)
N_MATRIX = 200
N_TRIALS = 500
BULK_FRAC = 0.3
t0 = time.time()

print("=" * 60)
print("[C-265 보충] S₁² 포함 A = S₁² + 2H₁ GUE 상관")
print("=" * 60)

all_A_2H1 = []
all_A_full = []  # S₁² + 2H₁
all_S1_sq = []
all_gap_min = []
all_gap_right = []

for trial in range(N_TRIALS):
    H = np.random.randn(N_MATRIX, N_MATRIX) + 1j * np.random.randn(N_MATRIX, N_MATRIX)
    H = (H + H.conj().T) / (2 * np.sqrt(2 * N_MATRIX))
    evals = np.sort(np.linalg.eigvalsh(H))

    n = len(evals)
    lo = int(n * (0.5 - BULK_FRAC / 2))
    hi = int(n * (0.5 + BULK_FRAC / 2))
    bulk = evals[lo:hi]
    spacings = np.diff(bulk)
    mean_sp = np.mean(spacings)
    if mean_sp <= 0:
        continue
    s_norm = spacings / mean_sp

    m = len(bulk)
    # K=10 이웃으로 S₁, H₁ 계산 (속도와 정밀도 균형)
    K = 10
    for i in range(K, m - K):
        S1 = 0.0
        H1 = 0.0
        for k in range(max(0, i - K), min(m, i + K + 1)):
            if k == i:
                continue
            d = bulk[i] - bulk[k]
            if abs(d) > 0:
                S1 += 1.0 / (d / mean_sp)    # 정규화
                H1 += 1.0 / (d / mean_sp)**2

        s_L = s_norm[i - 1] if i > 0 else 999
        s_R = s_norm[i] if i < len(s_norm) else 999

        all_A_2H1.append(2 * H1)
        all_A_full.append(S1**2 + 2 * H1)
        all_S1_sq.append(S1**2)
        all_gap_min.append(min(s_L, s_R))
        all_gap_right.append(s_R)

    if (trial + 1) % 100 == 0:
        print(f"  [{time.time()-t0:.0f}s] {trial+1}/{N_TRIALS}")

A_2H1 = np.array(all_A_2H1)
A_full = np.array(all_A_full)
S1_sq = np.array(all_S1_sq)
gm = np.array(all_gap_min)
gr = np.array(all_gap_right)

n = len(A_full)
print(f"\n  데이터: {n}개 포인트")

# 기여도
frac_H1 = np.mean(2 * (A_full - S1_sq) / (2 * A_full))
frac_S1 = np.mean(S1_sq / A_full)
print(f"\n[1] 기여도 분해")
print(f"  <2H₁ / A> = {1 - frac_S1:.3f} ({(1-frac_S1)*100:.1f}%)")
print(f"  <S₁² / A> = {frac_S1:.3f} ({frac_S1*100:.1f}%)")

# 상관
print(f"\n[2] Spearman 상관")
r1, _ = stats.spearmanr(A_2H1, gm)
r2, _ = stats.spearmanr(A_2H1, gr)
r3, _ = stats.spearmanr(A_full, gm)
r4, _ = stats.spearmanr(A_full, gr)
r5, _ = stats.spearmanr(S1_sq, gm)
r6, _ = stats.spearmanr(S1_sq, gr)

print(f"  ρ(2H₁,     gap_min)  = {r1:+.4f}")
print(f"  ρ(2H₁,     gap_right)= {r2:+.4f}")
print(f"  ρ(S₁²+2H₁, gap_min) = {r3:+.4f}")
print(f"  ρ(S₁²+2H₁, gap_right)= {r4:+.4f}")
print(f"  ρ(S₁²,     gap_min)  = {r5:+.4f}")
print(f"  ρ(S₁²,     gap_right)= {r6:+.4f}")

# L-함수 비교
print(f"\n[3] 관측 비교")
print(f"  GUE A=S₁²+2H₁ vs gap_right: {r4:+.4f}")
print(f"  관측 A_Λ        vs gap_right: -0.578 (GL(1)+GL(2) 평균)")
print(f"  차이: {abs(r4 - (-0.578)):.4f}")
print()
print(f"  GUE A=S₁²+2H₁ vs gap_min: {r3:+.4f}")
print(f"  관측 A_L        vs gap_min: -0.455 (GL(3) 평균)")
print(f"  차이: {abs(r3 - (-0.455)):.4f}")

# 판정
delta_right = abs(r4 - (-0.578))
delta_min = abs(r3 - (-0.455))
print(f"\n[4] 판정")
if delta_right < 0.05:
    print(f"  gap_right: ★★★ 일치 (Δ={delta_right:.4f} < 0.05)")
elif delta_right < 0.10:
    print(f"  gap_right: ★★ 근사 일치 (Δ={delta_right:.4f} < 0.10)")
else:
    print(f"  gap_right: ⚠️ 불일치 (Δ={delta_right:.4f})")

if delta_min < 0.05:
    print(f"  gap_min:   ★★★ 일치 (Δ={delta_min:.4f} < 0.05)")
elif delta_min < 0.10:
    print(f"  gap_min:   ★★ 근사 일치 (Δ={delta_min:.4f} < 0.10)")
else:
    print(f"  gap_min:   ⚠️ 불일치 (Δ={delta_min:.4f})")

print(f"\n  총 소요: {time.time()-t0:.1f}초")
