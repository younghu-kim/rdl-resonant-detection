#!/usr/bin/env python3
"""
사이클 #258: B-42 편상관 분석 + A(γ) NN 분해
C-256 데이터 재활용 (200영점, Laurent 계수 + gap_right 기존재)
새로 계산: gap_left (t차분), κ_mid (mpmath 중간점 ξ'/ξ 평가)
"""

import sys, os, time, math
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from scipy import stats
import mpmath

# ────────────────────────────────────────────────────────────────────────
# 설정
# ────────────────────────────────────────────────────────────────────────
DATA_FILE   = os.path.join(os.path.dirname(__file__), '..', 'results', 'A_gap_correlation_c256.txt')
OUTPUT_FILE = os.path.join(os.path.dirname(__file__), '..', 'results', 'B42_partial_correlation_c258.txt')
mpmath.mp.dps = 40   # κ_mid는 고정밀도 불필요, 40 dps 충분

t0_global = time.time()

# ────────────────────────────────────────────────────────────────────────
# 1. C-256 데이터 파싱
# ────────────────────────────────────────────────────────────────────────
print("[1/5] C-256 데이터 파싱 ...", flush=True)

rows = []
with open(DATA_FILE, 'r') as f:
    for line in f:
        line = line.strip()
        # 데이터 행: "   1    14.1347   ..."  (숫자로 시작)
        parts = line.split()
        if len(parts) >= 8:
            try:
                idx     = int(parts[0])
                t_val   = float(parts[1])
                re_c0   = float(parts[2])
                im_c0   = float(parts[3])
                re_c1   = float(parts[4])
                im_c1   = float(parts[5])
                A_val   = float(parts[6])
                g_right = parts[7]
                if g_right == 'nan':
                    g_right = float('nan')
                else:
                    g_right = float(g_right)
                rows.append((idx, t_val, re_c0, im_c0, re_c1, im_c1, A_val, g_right))
            except (ValueError, IndexError):
                pass

rows.sort(key=lambda r: r[0])
N = len(rows)
print(f"  파싱된 영점 수: {N}", flush=True)

t_arr      = np.array([r[1] for r in rows])
re_c0_arr  = np.array([r[2] for r in rows])
im_c0_arr  = np.array([r[3] for r in rows])
re_c1_arr  = np.array([r[4] for r in rows])
im_c1_arr  = np.array([r[5] for r in rows])
A_arr      = np.array([r[6] for r in rows])
gap_r_arr  = np.array([r[7] for r in rows])   # 마지막 nan 포함

# gap_left 계산: gap_left[n] = t[n] - t[n-1], n=1..N-1 (0-indexed: 1..N-1)
gap_l_arr  = np.full(N, float('nan'))
gap_l_arr[1:] = t_arr[1:] - t_arr[:-1]

print(f"  gap_right: {np.sum(np.isfinite(gap_r_arr))}개 유효")
print(f"  gap_left:  {np.sum(np.isfinite(gap_l_arr))}개 유효")

# ────────────────────────────────────────────────────────────────────────
# 2. κ_mid 계산 (ξ'/ξ 중간점 mpmath 평가)
# ────────────────────────────────────────────────────────────────────────
print("[2/5] κ_mid 계산 (199쌍 mpmath) ...", flush=True)

def xi_func(s):
    """ξ(s) = ½s(s-1)π^{-s/2}Γ(s/2)ζ(s)"""
    s = mpmath.mpc(s)
    return mpmath.fp.matrix  # 사용 안 함 — 아래에서 직접 계산

def xi_log_deriv(s):
    """log ξ'(s)/ξ(s) = d/ds [log ξ(s)] via mpmath"""
    s = mpmath.mpc(s)
    # ξ(s) = ½·s·(s-1)·π^{-s/2}·Γ(s/2)·ζ(s)
    # log ξ = log(½) + log(s) + log(s-1) - (s/2)·log(π) + log Γ(s/2) + log ζ(s)
    # d/ds log ξ = 1/s + 1/(s-1) - log(π)/2 + (1/2)·ψ(s/2) + ζ'(s)/ζ(s)
    psi_val  = mpmath.digamma(s / 2) / 2
    zeta_val = mpmath.zeta(s)
    if abs(zeta_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 5):
        return None   # ζ 영점 근방 (없어야 하지만 안전장치)
    dzeta = mpmath.diff(mpmath.zeta, s)
    return (1/s + 1/(s-1)
            - mpmath.log(mpmath.pi) / 2
            + psi_val
            + dzeta / zeta_val)

kappa_mid = np.full(N-1, float('nan'))
fail_count = 0

for i in range(N-1):
    if not np.isfinite(gap_r_arr[i]):
        continue
    t_mid = (t_arr[i] + t_arr[i+1]) / 2.0
    s_mid = mpmath.mpc(0.5, t_mid)
    try:
        ld = xi_log_deriv(s_mid)
        if ld is None:
            fail_count += 1
            print(f"  WARNING: ξ'/ξ 실패 i={i} t_mid={t_mid:.4f}", flush=True)
            continue
        kappa_mid[i] = float(abs(ld)**2)
    except Exception as e:
        fail_count += 1
        print(f"  WARNING: κ_mid 계산 실패 i={i} t_mid={t_mid:.4f}: {e}", flush=True)

valid_kappa = np.sum(np.isfinite(kappa_mid))
print(f"  κ_mid: {valid_kappa}/{N-1}개 유효, 실패={fail_count}", flush=True)
print(f"  κ_mid: mean={np.nanmean(kappa_mid):.4f}, std={np.nanstd(kappa_mid):.4f}", flush=True)
elapsed = time.time() - t0_global
print(f"  소요: {elapsed:.1f}초", flush=True)

if fail_count > (N-1) // 2:
    print("⚠️ κ_mid 실패율 50% 초과 — 중단", flush=True)
    sys.exit(1)

# ────────────────────────────────────────────────────────────────────────
# 3. NN 분해
# ────────────────────────────────────────────────────────────────────────
print("[3/5] NN 분해 ...", flush=True)

# Im(c₀)_NN = 1/gap_right - 1/gap_left  (Hadamard 최근접 기여)
# n≥2만 유효 (0-indexed: i≥1, 즉 gap_left[i] 유효한 경우)
im_c0_NN  = np.full(N, float('nan'))
valid_mask = np.isfinite(gap_r_arr) & np.isfinite(gap_l_arr)  # n≥2 조건
im_c0_NN[valid_mask] = (1.0 / gap_r_arr[valid_mask]
                        - 1.0 / gap_l_arr[valid_mask])

im_c0_bg  = im_c0_arr - im_c0_NN   # Im(c₀) - NN항

# A 성분 분해
A_NN    = im_c0_NN ** 2
A_bg    = im_c0_bg ** 2 + 2 * re_c1_arr
A_cross = 2 * im_c0_NN * im_c0_bg

# 분산 기여율 (valid 영점만)
vmask = valid_mask & np.isfinite(A_NN) & np.isfinite(A_bg) & np.isfinite(A_cross)
A_total_check = A_NN[vmask] + A_bg[vmask] + A_cross[vmask]
A_actual      = A_arr[vmask]
resid_rms = np.sqrt(np.mean((A_total_check - A_actual)**2))

var_NN    = np.var(A_NN[vmask])
var_bg    = np.var(A_bg[vmask])
var_cross = np.var(A_cross[vmask])
var_total = np.var(A_actual)

print(f"  NN 분해 유효: {np.sum(vmask)}개")
print(f"  재구성 RMS 잔차: {resid_rms:.2e}")
print(f"  분산 기여율 — A_NN: {var_NN/var_total*100:.1f}%  "
      f"A_bg: {var_bg/var_total*100:.1f}%  "
      f"A_cross: {var_cross/var_total*100:.1f}%")

# ────────────────────────────────────────────────────────────────────────
# 4. 편상관 (B-42) + 짝 패턴 분리
# ────────────────────────────────────────────────────────────────────────
print("[4/5] 편상관 + 짝 패턴 분리 ...", flush=True)

def partial_corr_spearman(x, y, z):
    """
    Spearman 편상관 ρ(x, y | z) 계산
    방법: z에 대한 선형 잔차(rank 공간)
    """
    # rank 변환
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    # z에 대한 rx, ry 회귀 잔차
    slope_xz, intercept_xz, _, _, _ = stats.linregress(rz, rx)
    slope_yz, intercept_yz, _, _, _ = stats.linregress(rz, ry)
    resid_x = rx - (slope_xz * rz + intercept_xz)
    resid_y = ry - (slope_yz * rz + intercept_yz)
    rho, p = stats.pearsonr(resid_x, resid_y)
    return rho, p

def partial_corr_spearman_2z(x, y, z1, z2):
    """
    ρ(x, y | z1, z2) — z1, z2 두 변수 동시 통제
    """
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz1 = stats.rankdata(z1)
    rz2 = stats.rankdata(z2)
    Z = np.column_stack([rz1, rz2, np.ones(len(x))])
    # 최소제곱 잔차
    resid_x = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    resid_y = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    rho, p = stats.pearsonr(resid_x, resid_y)
    return rho, p

# ── B-42 편상관용 쌍별 데이터 준비 ──
# 쌍 i: (i, i+1), i=1..N-2 (i≥1 → gap_left 유효, i+1 < N-1 → gap_right 유효)
pair_idx = []
for i in range(1, N-1):
    if (np.isfinite(kappa_mid[i]) and
        np.isfinite(A_NN[i]) and np.isfinite(A_NN[i+1]) and
        np.isfinite(gap_r_arr[i]) and np.isfinite(gap_l_arr[i])):
        pair_idx.append(i)

pair_idx = np.array(pair_idx, dtype=int)
n_pairs = len(pair_idx)
print(f"  B-42 편상관 유효 쌍: {n_pairs}개", flush=True)

kap   = kappa_mid[pair_idx]
A_avg = (A_arr[pair_idx] + A_arr[pair_idx + 1]) / 2.0
gap   = gap_r_arr[pair_idx]
A_NN_avg = (A_NN[pair_idx] + A_NN[pair_idx + 1]) / 2.0

# 비조건부 상관 (재현)
rho_kap_A,   p_kap_A   = stats.spearmanr(kap, A_avg)
rho_kap_gap, p_kap_gap = stats.spearmanr(kap, gap)

# 편상관
rho_kap_A_cond_gap,   p1 = partial_corr_spearman(kap, A_avg, gap)
rho_kap_A_cond_ANN,   p2 = partial_corr_spearman(kap, A_avg, A_NN_avg)
rho_kap_gap_cond_A,   p3 = partial_corr_spearman(kap, gap, A_avg)

print(f"  ρ(κ, A_avg)           = {rho_kap_A:+.4f}  (p={p_kap_A:.3e})")
print(f"  ρ(κ, gap)             = {rho_kap_gap:+.4f}  (p={p_kap_gap:.3e})")
print(f"  ρ(κ, A_avg | gap)     = {rho_kap_A_cond_gap:+.4f}  (p={p1:.3e})")
print(f"  ρ(κ, A_avg | A_NN)    = {rho_kap_A_cond_ANN:+.4f}  (p={p2:.3e})")
print(f"  ρ(κ, gap | A_avg)     = {rho_kap_gap_cond_A:+.4f}  (p={p3:.3e})")

# ── 짝 패턴 분리: A_bg 인접 상관 ──
bg_idx = []
for i in range(N-1):
    if np.isfinite(A_bg[i]) and np.isfinite(A_bg[i+1]):
        bg_idx.append(i)

A_bg_n   = A_bg[bg_idx]
A_bg_n1  = A_bg[[i+1 for i in bg_idx]]
rho_bg, p_bg = stats.spearmanr(A_bg_n, A_bg_n1)

# 원시 A_n, A_{n+1} (재현)
rho_A_raw, p_A_raw = stats.spearmanr(A_arr[:-1], A_arr[1:])
print(f"  ρ(A_n, A_n+1) 원시         = {rho_A_raw:+.4f}  (p={p_A_raw:.3e})")
print(f"  ρ(A_bg_n, A_bg_n+1) NN제거  = {rho_bg:+.4f}  (p={p_bg:.3e})")

# ── gap 비대칭 정량: ρ(A, A_NN) vs ρ(A, 2c₁) ──
vmask2 = valid_mask & np.isfinite(A_NN)
rho_A_ANN, p_A_ANN = stats.spearmanr(A_arr[vmask2], A_NN[vmask2])
rho_A_2c1, p_A_2c1 = stats.spearmanr(A_arr[vmask2], 2 * re_c1_arr[vmask2])
print(f"  ρ(A, A_NN)  = {rho_A_ANN:+.4f}  (p={p_A_ANN:.3e})")
print(f"  ρ(A, 2c₁)   = {rho_A_2c1:+.4f}  (p={p_A_2c1:.3e})")

# ────────────────────────────────────────────────────────────────────────
# 5. 결과 저장
# ────────────────────────────────────────────────────────────────────────
print("[5/5] 결과 저장 ...", flush=True)
elapsed_total = time.time() - t0_global

def fmt_p(p):
    if p < 0.01:
        return f"{p:.3e}  유의(p<0.01)"
    elif p < 0.05:
        return f"{p:.3e}  유의(p<0.05)"
    else:
        return f"{p:.3e}  비유의"

# 판정 로직
if abs(rho_kap_A_cond_gap) > 0.3:
    verdict = "★★★ 양성 — A가 gap-독립적으로 κ_mid 설명 (비자명)"
elif abs(rho_kap_A_cond_gap) > 0.1:
    verdict = "⚠️ 중간 — 약한 gap-독립 상관, 추가분석 필요"
else:
    verdict = "음성 — gap 통제 시 A-κ 상관 소멸 (가성상관)"

bonus = "★ 보너스 양성" if rho_bg > 0 and p_bg < 0.05 else "보너스 음성"

lines = []
lines.append("=" * 80)
lines.append("[사이클 #258] B-42 편상관 분석 + A(γ) NN 분해")
lines.append(f"생성: {time.strftime('%Y-%m-%d %H:%M:%S')}")
lines.append(f"소요: {elapsed_total:.1f}초  N={N}  dps={mpmath.mp.dps}")
lines.append("=" * 80)

lines.append("")
lines.append("설정:")
lines.append(f"  C-256 데이터: {DATA_FILE}")
lines.append(f"  κ_mid 유효: {valid_kappa}/{N-1}")
lines.append(f"  B-42 편상관 쌍: {n_pairs}")

lines.append("")
lines.append("=" * 80)
lines.append("§1 NN 분해 결과")
lines.append("=" * 80)
lines.append("")
lines.append("  Im(c₀) = Im(c₀)_NN + Im(c₀)_bg")
lines.append("  Im(c₀)_NN = 1/gap_right - 1/gap_left  [Hadamard NN 기여]")
lines.append("")
lines.append(f"  유효 영점 (n≥2): {np.sum(vmask)}개")
lines.append(f"  재구성 RMS 잔차 |A_NN+A_bg+A_cross - A|: {resid_rms:.2e}")
lines.append("")
lines.append("  분산 기여율:")
lines.append(f"    A_NN    = Im(c₀)_NN²          : {var_NN/var_total*100:.1f}%")
lines.append(f"    A_bg    = Im(c₀)_bg²+2Re(c₁) : {var_bg/var_total*100:.1f}%")
lines.append(f"    A_cross = 2·Im(c₀)_NN·Im(c₀)_bg: {var_cross/var_total*100:.1f}%")
lines.append(f"    합계 (교차항 포함): ~100%")
lines.append("")
lines.append("  성분 통계 (유효 영점):")
lines.append(f"    A_NN:   mean={np.nanmean(A_NN):.4f}, std={np.nanstd(A_NN):.4f}")
lines.append(f"    A_bg:   mean={np.nanmean(A_bg):.4f}, std={np.nanstd(A_bg):.4f}")
lines.append(f"    A_cross: mean={np.nanmean(A_cross):.4f}, std={np.nanstd(A_cross):.4f}")

lines.append("")
lines.append("=" * 80)
lines.append("§2 gap 비대칭 정량")
lines.append("=" * 80)
lines.append("")
lines.append("  가설: A ~ Im(c₀)_NN² (gap_right 지배)이면 ρ(A, A_NN) >> ρ(A, 2c₁)")
lines.append("")
lines.append(f"  ρ(A, A_NN)  = {rho_A_ANN:+.6f}  ({fmt_p(p_A_ANN)})")
lines.append(f"  ρ(A, 2Re(c₁)) = {rho_A_2c1:+.6f}  ({fmt_p(p_A_2c1)})")
lines.append("")
if abs(rho_A_ANN) > abs(rho_A_2c1) + 0.1:
    lines.append("  → NN항이 A 분산을 지배: Hadamard gap_right 기여 확인")
elif abs(rho_A_2c1) > abs(rho_A_ANN) + 0.1:
    lines.append("  → c₁항이 A 분산을 지배: 배경항 중요")
else:
    lines.append("  → NN항과 c₁항 기여 비슷")

lines.append("")
lines.append("=" * 80)
lines.append("§3 B-42 편상관 분석")
lines.append("=" * 80)
lines.append("")
lines.append(f"  쌍 수: {n_pairs}  (i=1..N-2, κ_mid·A_NN 유효)")
lines.append("")
lines.append("  비조건부 상관 (재현):")
lines.append(f"    ρ(κ_mid, A_avg)  = {rho_kap_A:+.6f}  ({fmt_p(p_kap_A)})")
lines.append(f"    ρ(κ_mid, gap)    = {rho_kap_gap:+.6f}  ({fmt_p(p_kap_gap)})")
lines.append("")
lines.append("  편상관:")
lines.append(f"    ρ(κ_mid, A_avg | gap)  = {rho_kap_A_cond_gap:+.6f}  ({fmt_p(p1)})")
lines.append(f"    ρ(κ_mid, A_avg | A_NN) = {rho_kap_A_cond_ANN:+.6f}  ({fmt_p(p2)})")
lines.append(f"    ρ(κ_mid, gap | A_avg)  = {rho_kap_gap_cond_A:+.6f}  ({fmt_p(p3)})")
lines.append("")
lines.append(f"  판정: {verdict}")
lines.append("")
if abs(rho_kap_A_cond_gap) > 0.3:
    lines.append("  해석: gap을 통제해도 κ_mid ~ A_avg 상관이 유지됨")
    lines.append("        → A(γ)는 gap의 대리변수가 아닌 독립 정보를 κ_mid와 공유")
    lines.append("        → B-42: κ_mid ↔ A 경로가 gap-독립적으로 실재")
elif abs(rho_kap_A_cond_gap) < 0.1:
    lines.append("  해석: gap 통제 시 κ_mid ~ A_avg 상관 소멸")
    lines.append("        → C-256의 B-42 상관은 gap을 통한 가성상관")
    lines.append("        → B-42 닫힘")
else:
    lines.append("  해석: 부분적 gap-독립 상관. 추가 분석 필요.")

lines.append("")
lines.append("=" * 80)
lines.append("§4 짝 패턴 분리 (NN 효과 제거)")
lines.append("=" * 80)
lines.append("")
lines.append("  원시 A vs NN항 제거 후 A_bg의 인접 상관 비교:")
lines.append(f"    ρ(A_n, A_n+1) 원시          = {rho_A_raw:+.6f}  ({fmt_p(p_A_raw)})")
lines.append(f"    ρ(A_bg_n, A_bg_n+1) NN제거  = {rho_bg:+.6f}  ({fmt_p(p_bg)})")
lines.append("")
if rho_bg > 0 and p_bg < 0.05:
    lines.append(f"  {bonus}: NN 제거 후에도 인접 A_bg 상관 양성")
    lines.append("  → 짝 패턴이 NN 공유 효과 너머 비자명 구조를 가짐")
    if rho_bg < rho_A_raw - 0.05:
        frac_nn = (rho_A_raw - rho_bg) / rho_A_raw * 100
        lines.append(f"  → 원시 상관의 ~{frac_nn:.0f}%는 NN 공유 효과, 나머지는 비자명")
else:
    lines.append(f"  {bonus}: NN 제거 후 A_bg 인접 상관 소멸/음성")
    lines.append("  → 짝 패턴은 NN 공유 효과(자명)로 설명 가능")

lines.append("")
lines.append("=" * 80)
lines.append("종합 판정")
lines.append("=" * 80)
lines.append("")
lines.append(f"  주판정: {verdict}")
lines.append(f"  보너스: {bonus}")
lines.append("")
lines.append("  B-42 경계 상태:")
if abs(rho_kap_A_cond_gap) > 0.3:
    lines.append("    ⏳→✅ 부분 해결: A가 gap-독립적으로 κ_mid와 상관")
    lines.append("    완전 해결을 위해: A_NN vs A_bg 중 어떤 성분이 κ_mid와 연결되는지 추가 필요")
elif abs(rho_kap_A_cond_gap) < 0.1:
    lines.append("    ⏳→✅ 닫힘: gap 가성상관으로 판명")
else:
    lines.append("    ⏳ 열린 문제: 부분적 증거, 추가 분석 필요")

lines.append("")
lines.append(f"  소요: {elapsed_total:.1f}초")

with open(OUTPUT_FILE, 'w') as f:
    f.write("\n".join(lines) + "\n")

print(f"\n결과 저장: {OUTPUT_FILE}")
print(f"총 소요: {elapsed_total:.1f}초")
print("\n[완료]", flush=True)
