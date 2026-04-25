#!/usr/bin/env python3
"""
[사이클 #265] GUE 이론 예측: A-gap 상관의 보편값

목적:
  Obs 3 (A-gap anti-correlation ρ ≈ -0.57)가 GUE 통계의 필연적 귀결인지 검증.

  논증:
    1. Hadamard 분해: A(γ) = Im(c₀)² + 2Re(c₁)
    2. NN 근사: A ≈ 2H₁ ≈ 2d̄²(1/s_L² + 1/s_R²), s=정규화 간격
    3. GUE 고유값의 간격 분포에서 ρ(A_NN, gap_min) 계산
    4. 관측 ρ ≈ -0.57과 비교

  만약 GUE 시뮬레이션이 ρ ≈ -0.57을 재현하면:
    → A-gap anti-correlation = Hadamard structure + GUE universality (QED)
    → Paper 4 중심 정리의 이론적 근거 확보

방법:
  1. GUE 랜덤 행렬 (N×N) 반복 생성
  2. 벌크 고유값 추출 + 국소 언폴딩
  3. NN Hadamard 근사 A_NN 계산
  4. Full Hadamard (K-이웃) 근사 A_full 계산
  5. Spearman ρ(A, gap_min), ρ(A, gap_right) 측정
  6. K 의존성: NN only vs 5-이웃 vs 10-이웃 vs 전체
"""

import numpy as np
from scipy import stats
import time

np.random.seed(42)

print("=" * 70)
print("[사이클 #265] GUE 이론 예측: A-gap Spearman 상관")
print("=" * 70)

# ── 설정 ──
N_MATRIX = 300       # 행렬 크기
N_TRIALS = 1000      # 반복 횟수
BULK_FRAC = 0.3      # 벌크 중앙 비율 (양쪽 가장자리 제외)
K_NEIGHBORS = [1, 2, 5, 10, 50]  # Hadamard 이웃 수

t0 = time.time()

# ── GUE 샘플링 ──
print(f"\n[1] GUE 행렬 생성 (N={N_MATRIX}, {N_TRIALS}회)...")

# 각 K에 대해 수집
results_by_K = {K: {'A': [], 'gap_min': [], 'gap_right': [], 'gap_left': []}
                for K in K_NEIGHBORS}
results_by_K['full'] = {'A': [], 'gap_min': [], 'gap_right': [], 'gap_left': []}

for trial in range(N_TRIALS):
    # GUE 행렬 생성
    H = np.random.randn(N_MATRIX, N_MATRIX) + 1j * np.random.randn(N_MATRIX, N_MATRIX)
    H = (H + H.conj().T) / (2 * np.sqrt(2 * N_MATRIX))  # Wigner 정규화
    evals = np.sort(np.linalg.eigvalsh(H))

    # 벌크 추출 (중앙 BULK_FRAC)
    n = len(evals)
    lo = int(n * (0.5 - BULK_FRAC / 2))
    hi = int(n * (0.5 + BULK_FRAC / 2))
    bulk = evals[lo:hi]

    # 국소 언폴딩: 간격 / 평균 간격
    spacings = np.diff(bulk)
    mean_sp = np.mean(spacings)
    if mean_sp <= 0:
        continue
    s = spacings / mean_sp  # 정규화 간격

    # 각 "영점" i에 대해 (i=1..len(s)-1, 양쪽 간격 필요)
    for i in range(1, len(s)):
        s_L = s[i - 1]  # 왼쪽 간격
        s_R = s[i]  if i < len(s) else None
        if s_R is None or s_L < 1e-10 or s_R < 1e-10:
            continue

        g_min = min(s_L, s_R)
        g_right = s_R
        g_left = s_L

        # Full Hadamard: Σ_{k≠i} 1/(E_i - E_k)² (정규화 단위)
        # 인덱스 i는 bulk[i]에 해당, 간격 s 인덱스와는 1 차이
        # bulk 인덱스에서 i+lo가 원래 고유값 인덱스
        # 대신 정규화 거리로 계산

        # K-이웃 Hadamard
        for K in K_NEIGHBORS:
            H1_K = 0.0
            for k in range(max(0, i - K), min(len(s) + 1, i + K + 1)):
                if k == i:
                    continue
                # 정규화 거리 (간격 누적)
                if k < i:
                    dist = sum(s[k:i])
                else:
                    dist = sum(s[i:k])
                if dist > 0:
                    H1_K += 1.0 / dist**2
            A_K = 2 * H1_K  # 2H₁ (S₁² 무시 — 실수 고유값이므로)
            results_by_K[K]['A'].append(A_K)
            results_by_K[K]['gap_min'].append(g_min)
            results_by_K[K]['gap_right'].append(g_right)
            results_by_K[K]['gap_left'].append(g_left)

        # Full Hadamard (전체 벌크)
        H1_full = 0.0
        for k in range(len(s) + 1):
            if k == i:
                continue
            if k < i:
                dist = sum(s[k:i])
            else:
                dist = sum(s[i:k])
            if dist > 0:
                H1_full += 1.0 / dist**2
        A_full = 2 * H1_full
        results_by_K['full']['A'].append(A_full)
        results_by_K['full']['gap_min'].append(g_min)
        results_by_K['full']['gap_right'].append(g_right)
        results_by_K['full']['gap_left'].append(g_left)

    if (trial + 1) % 200 == 0:
        elapsed = time.time() - t0
        print(f"  [{elapsed:.0f}s] {trial+1}/{N_TRIALS} 행렬 완료")

elapsed = time.time() - t0
n_total = len(results_by_K[1]['A'])
print(f"  완료: {n_total}개 데이터 포인트, {elapsed:.1f}초")

# ── Spearman 상관 계산 ──
print(f"\n[2] Spearman 상관 측정")
print("-" * 70)
print(f"{'K':>6} | {'n':>8} | {'ρ(A,gap_min)':>14} | {'ρ(A,gap_right)':>14} | {'ρ(A,gap_left)':>14}")
print("-" * 70)

for K_label in K_NEIGHBORS + ['full']:
    data = results_by_K[K_label]
    A_arr = np.array(data['A'])
    gm_arr = np.array(data['gap_min'])
    gr_arr = np.array(data['gap_right'])
    gl_arr = np.array(data['gap_left'])

    n = len(A_arr)
    if n < 10:
        continue

    rho_min, p_min = stats.spearmanr(A_arr, gm_arr)
    rho_right, p_right = stats.spearmanr(A_arr, gr_arr)
    rho_left, p_left = stats.spearmanr(A_arr, gl_arr)

    print(f"{'K='+str(K_label):>6} | {n:>8} | {rho_min:>+.4f} (p={p_min:.1e}) | "
          f"{rho_right:>+.4f} (p={p_right:.1e}) | {rho_left:>+.4f} (p={p_left:.1e})")

# ── NN 기여도 분석 ──
print(f"\n[3] NN 기여도 분석 (Full 대비)")
data_full = results_by_K['full']
data_nn = results_by_K[1]
A_full_arr = np.array(data_full['A'])
A_nn_arr = np.array(data_nn['A'])

if len(A_full_arr) > 0 and len(A_nn_arr) > 0:
    nn_frac = np.mean(A_nn_arr / A_full_arr)
    rho_nn_full, _ = stats.spearmanr(A_nn_arr, A_full_arr)
    print(f"  <A_NN / A_full> = {nn_frac:.4f} ({nn_frac*100:.1f}%)")
    print(f"  ρ(A_NN, A_full) = {rho_nn_full:.4f}")

# ── 부트스트랩 신뢰구간 ──
print(f"\n[4] 부트스트랩 95% CI (K=1, gap_min)")
data = results_by_K[1]
A_arr = np.array(data['A'])
gm_arr = np.array(data['gap_min'])
n_boot = 2000
rho_boot = []
rng = np.random.default_rng(123)
for _ in range(n_boot):
    idx = rng.choice(len(A_arr), size=len(A_arr), replace=True)
    r, _ = stats.spearmanr(A_arr[idx], gm_arr[idx])
    rho_boot.append(r)
rho_boot = np.sort(rho_boot)
ci_lo, ci_hi = rho_boot[int(0.025 * n_boot)], rho_boot[int(0.975 * n_boot)]
rho_mean = np.mean(rho_boot)
print(f"  ρ = {rho_mean:.4f}, 95% CI = [{ci_lo:.4f}, {ci_hi:.4f}]")

# ── 비교: 관측 ρ 값 ──
print(f"\n[5] 관측 vs GUE 예측 비교")
print("-" * 70)
observed = {
    'ζ(s)': -0.5898,
    'χ₃ mod 3': -0.5733,
    'χ₄ mod 4': -0.5530,
    'χ₅ mod 5': -0.6334,
    '11a1 GL(2)': -0.5688,
    '37a1 GL(2)': -0.5500,
    'Sym²(11a1) GL(3) gap_min': -0.42,
    'Sym²(37a1) GL(3) gap_min': -0.49,
}

# gap_right 관측과 비교
gue_rho_min_nn = stats.spearmanr(np.array(results_by_K[1]['A']),
                                  np.array(results_by_K[1]['gap_min']))[0]
gue_rho_right_nn = stats.spearmanr(np.array(results_by_K[1]['A']),
                                    np.array(results_by_K[1]['gap_right']))[0]
gue_rho_min_full = stats.spearmanr(np.array(results_by_K['full']['A']),
                                    np.array(results_by_K['full']['gap_min']))[0]
gue_rho_right_full = stats.spearmanr(np.array(results_by_K['full']['A']),
                                      np.array(results_by_K['full']['gap_right']))[0]

print(f"  GUE 예측 (NN,  gap_min):  ρ = {gue_rho_min_nn:+.4f}")
print(f"  GUE 예측 (NN,  gap_right): ρ = {gue_rho_right_nn:+.4f}")
print(f"  GUE 예측 (full, gap_min):  ρ = {gue_rho_min_full:+.4f}")
print(f"  GUE 예측 (full, gap_right): ρ = {gue_rho_right_full:+.4f}")
print()

obs_gap_right = [v for k, v in observed.items() if 'gap_min' not in k]
obs_gap_min = [v for k, v in observed.items() if 'gap_min' in k]
mean_obs_right = np.mean(obs_gap_right)
mean_obs_min = np.mean(obs_gap_min) if obs_gap_min else None

print(f"  관측 평균 (gap_right, GL(1)+GL(2)):  ρ = {mean_obs_right:+.4f}")
if mean_obs_min is not None:
    print(f"  관측 평균 (gap_min,  GL(3)):           ρ = {mean_obs_min:+.4f}")

# ── S₁² 효과 추정 (L-함수 고유) ──
print(f"\n[6] L-함수 고유 효과: S₁² 기여")
print("  GUE 실수 고유값에서 S₁=Σ 1/(E_n-E_k)는 자명하지 않음")
print("  L-함수에서 c₀ = ib (순허수) → S₁² = Im(c₀)² = b²")
print("  GUE에서 유사량: (Σ 1/(E_n-E_k))² — 이것도 측정")

# S₁ 유사량 측정 (full Hadamard)
S1_data = []
for trial_idx in range(min(200, N_TRIALS)):
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
    s = spacings / mean_sp

    for i in range(1, min(len(s), 10)):  # 처음 10개만 샘플
        S1 = 0.0
        H1 = 0.0
        for k in range(len(s) + 1):
            if k == i:
                continue
            if k < i:
                dist = sum(s[k:i])
            else:
                dist = sum(s[i:k])
            if dist > 0:
                S1 += 1.0 / dist * (1 if k > i else -1)
                H1 += 1.0 / dist**2
        S1_data.append({'S1_sq': S1**2, 'H1': H1, 'A': S1**2 + 2*H1,
                        'gap_min': min(s[i-1], s[i] if i < len(s) else s[i-1])})

if S1_data:
    S1_sq_arr = np.array([d['S1_sq'] for d in S1_data])
    H1_arr = np.array([d['H1'] for d in S1_data])
    A_with_S1 = np.array([d['A'] for d in S1_data])
    gm_S1 = np.array([d['gap_min'] for d in S1_data])

    frac_H1 = np.mean(2*H1_arr / A_with_S1)
    frac_S1 = np.mean(S1_sq_arr / A_with_S1)
    rho_with_S1, _ = stats.spearmanr(A_with_S1, gm_S1)

    print(f"  <2H₁/A> = {frac_H1:.3f} ({frac_H1*100:.1f}%)")
    print(f"  <S₁²/A> = {frac_S1:.3f} ({frac_S1*100:.1f}%)")
    print(f"  ρ(A=S₁²+2H₁, gap_min) = {rho_with_S1:+.4f}")

# ── 인접 쌍 상관 ──
print(f"\n[7] 인접 쌍 상관: A_n vs A_{n+1}")
data = results_by_K['full']
A_arr = np.array(data['A'])
# 같은 행렬 내 인접 쌍만 의미 있지만, 데이터가 섞여있으므로 근사
# 실제로는 연속 인덱스가 같은 행렬 내 인접을 의미
# 1개 행렬당 ~BULK_FRAC*N_MATRIX 개 포인트
pts_per_matrix = int(BULK_FRAC * N_MATRIX) - 2
if pts_per_matrix > 2 and len(A_arr) > pts_per_matrix:
    adj_corrs = []
    for start in range(0, len(A_arr) - pts_per_matrix, pts_per_matrix):
        chunk = A_arr[start:start + pts_per_matrix]
        if len(chunk) > 3:
            r, _ = stats.spearmanr(chunk[:-1], chunk[1:])
            if not np.isnan(r):
                adj_corrs.append(r)
    if adj_corrs:
        print(f"  ρ(A_n, A_{{n+1}}) = {np.mean(adj_corrs):+.4f} ± {np.std(adj_corrs):.4f}")
        print(f"  관측: ρ ≈ +0.28 ~ +0.42")

# ── 결론 ──
print(f"\n{'=' * 70}")
print("[결론]")
diff_right = abs(gue_rho_right_full - mean_obs_right)
match_right = "일치" if diff_right < 0.05 else ("근사" if diff_right < 0.10 else "불일치")
print(f"  GUE ρ(A_full, gap_right) = {gue_rho_right_full:+.4f}  vs  관측 {mean_obs_right:+.4f}  → {match_right} (Δ={diff_right:.4f})")

if mean_obs_min is not None:
    diff_min = abs(gue_rho_min_full - mean_obs_min)
    match_min = "일치" if diff_min < 0.05 else ("근사" if diff_min < 0.10 else "불일치")
    print(f"  GUE ρ(A_full, gap_min)  = {gue_rho_min_full:+.4f}  vs  관측 {mean_obs_min:+.4f}  → {match_min} (Δ={diff_min:.4f})")

print()
total = time.time() - t0
print(f"  총 소요: {total:.1f}초")
print(f"{'=' * 70}")
