"""
=============================================================================
[Project RDL] 디리클레 off-critical 해부
=============================================================================
목표: 고 conductor (mod 7,8,11)에서 타σ 위상점프=8이 체계적으로 출현하는 원인 규명

비교 대상:
  (a) ξ(s,χ) = Λ(s,χ) = (q/π)^{s/2} Γ((s+a)/2) L(s,χ)  — 완비 함수
  (b) L(s,χ) = mpmath.dirichlet(s, chi)                    — 순수 L 함수 (Gamma 제거)

방법: ζ off-critical (bundle_prediction_2_offcritical.py)과 동일
  - branch cut 보정 있음 (±π wrap)
  - threshold = π/2
  - t∈[10,40], 2000점

실험:
  1. χ₇(mod 7): σ ∈ {0.10, 0.15, ..., 0.90} (17점) 전체 프로파일
     - Λ(s,χ₇) vs L(s,χ₇) 위상점프 프로파일 비교
  2. χ₈, χ₁₁: σ=0.5 + 타σ 최대점프 σ에서만 확인 (효율)

진단 추가:
  - 구 방법 (threshold=2.5, no branch cut) 재현: 왜 8이 나왔는지 확인

성공 기준:
  - Γ 제거 후 타σ 점프=0 → 양성 (Gamma 아티팩트, Conj 1은 L 레벨에서 성립)
  - Γ 제거 후에도 점프 잔존하되 감소 → 중립 (부분적 설명)
  - Γ 제거 후에도 타σ 점프=8 그대로 → 음성 (본질적 한계)

결과 파일: results/dirichlet_offcritical_anatomy.txt
=============================================================================
"""

import sys, os, time, cmath
import numpy as np
import mpmath
from datetime import datetime

# bundle_utils import
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))
from bundle_utils import completed_L

mpmath.mp.dps = 50  # t<50이므로 30 충분, 안전 여유치 50

T_MIN, T_MAX = 10.0, 40.0
N_POINTS = 2000

print(f"[시작] {datetime.now().strftime('%H:%M:%S')}  "
      f"정밀도={mpmath.mp.dps}dps  t∈[{T_MIN},{T_MAX}]  {N_POINTS}점",
      flush=True)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 고 conductor 지표 정의 (dirichlet_high_conductor.py와 동일)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# mod 7: primitive root = 3, order 6 character
w6 = cmath.exp(2j * cmath.pi / 6)
_chi7 = [0, 1, w6**2, w6, w6**4, w6**5, w6**3]  # 인덱스: 3^k mod 7 → 1→0,3→1,2→2,6→3,4→4,5→5

# mod 8: real character (Kronecker (-1/n)), χ(-1)=+1 → a=0
_chi8 = [0, 1, 0, -1, 0, -1, 0, 1]

# mod 11: primitive root = 2, order 10 character
w10 = cmath.exp(2j * cmath.pi / 10)
_inds11 = {1:0, 2:1, 4:2, 8:3, 5:4, 10:5, 9:6, 7:7, 3:8, 6:9}
_chi11 = [0] + [w10**_inds11[k] for k in range(1, 11)]

CHARS_HIGH = {
    'χ₇ (mod 7)': {'chi': _chi7, 'q': 7, 'a': 1, 'type': 'complex'},
    'χ₈ (mod 8)': {'chi': _chi8, 'q': 8, 'a': 0, 'type': 'real'},
    'χ₁₁ (mod 11)': {'chi': _chi11, 'q': 11, 'a': 1, 'type': 'complex'},
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 함수 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _near_zero(val):
    """dps 기반 영점 판정"""
    return abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10)


def pure_L(s, char_info):
    """순수 L-함수: L(s,χ) — Gamma/prefactor 없음"""
    return mpmath.dirichlet(s, char_info['chi'])


def count_jumps_new(sigma, char_info, fn, t_min=T_MIN, t_max=T_MAX, n_points=N_POINTS):
    """
    ζ off-critical과 동일한 방법:
    - branch cut 보정 있음 (±π wrap)
    - threshold = π/2
    fn: callable (s, char_info) → 복소수
    """
    ts = np.linspace(t_min, t_max, n_points)
    jumps = 0
    prev_arg = None

    for t in ts:
        s = mpmath.mpc(sigma, t)
        try:
            val = fn(s, char_info)
        except Exception as e:
            print(f"  WARNING fn({sigma},{t:.2f}): {e}", flush=True)
            prev_arg = None
            continue

        if _near_zero(val):
            prev_arg = None
            continue

        curr_arg = float(mpmath.arg(val))

        if prev_arg is not None:
            delta = curr_arg - prev_arg
            # branch cut 보정
            while delta > np.pi:
                delta -= 2 * np.pi
            while delta < -np.pi:
                delta += 2 * np.pi
            if abs(delta) > np.pi / 2:
                jumps += 1

        prev_arg = curr_arg

    return jumps


def count_jumps_old(sigma, char_info, fn, t_min=T_MIN, t_max=T_MAX, n_points=1500,
                    threshold=2.5):
    """
    구 방법 (dirichlet_high_conductor.py와 동일):
    - branch cut 보정 없음
    - threshold=2.5
    - 1500점
    """
    ts = np.linspace(t_min, t_max, n_points)
    jumps = 0
    prev_phase = None

    for t in ts:
        s = mpmath.mpc(sigma, t)
        try:
            val = fn(s, char_info)
        except Exception as e:
            print(f"  WARNING fn_old({sigma},{t:.2f}): {e}", flush=True)
            continue

        if abs(val) > 0:
            try:
                phase = float(mpmath.arg(val))
            except Exception:
                continue
            if prev_phase is not None:
                diff = abs(phase - prev_phase)
                if diff > threshold:
                    jumps += 1
            prev_phase = phase

    return jumps


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실험 1: χ₇ 전체 프로파일 (Λ vs L, 새 방법)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

print("\n" + "="*70, flush=True)
print("실험 1: χ₇(mod 7) σ 프로파일 — Λ(s,χ) vs L(s,χ) [새 방법]", flush=True)
print("="*70, flush=True)

ci7 = CHARS_HIGH['χ₇ (mod 7)']
sigmas_17 = np.round(np.arange(0.10, 0.91, 0.05), 2)  # 17점: 0.10,0.15,...,0.90
print(f"σ 값: {sigmas_17}", flush=True)
print(f"총 {len(sigmas_17)}개 σ × 2함수 × {N_POINTS}점 = {len(sigmas_17)*2*N_POINTS}점 평가\n", flush=True)

jumps_lambda_7 = {}  # Λ(s,χ₇)
jumps_pure_7   = {}  # L(s,χ₇)

t0_exp1 = time.time()

for sigma in sigmas_17:
    marker = "★" if abs(sigma - 0.5) < 0.01 else " "

    j_lam = count_jumps_new(sigma, ci7, completed_L)
    j_pur = count_jumps_new(sigma, ci7, pure_L)

    jumps_lambda_7[sigma] = j_lam
    jumps_pure_7[sigma]   = j_pur

    print(f"  {marker} σ={sigma:.2f}:  Λ={j_lam:3d}  L={j_pur:3d}", flush=True)

elapsed1 = time.time() - t0_exp1
print(f"\n[실험 1 소요: {elapsed1:.0f}초]", flush=True)

# 요약
j_half_lam7 = jumps_lambda_7[0.50]
j_half_pur7 = jumps_pure_7[0.50]
other_lam7  = {s: v for s, v in jumps_lambda_7.items() if abs(s - 0.5) > 0.01}
other_pur7  = {s: v for s, v in jumps_pure_7.items()   if abs(s - 0.5) > 0.01}
max_other_lam7 = max(other_lam7.values()) if other_lam7 else 0
max_other_pur7 = max(other_pur7.values()) if other_pur7 else 0
argmax_lam7 = max(other_lam7, key=other_lam7.get) if other_lam7 else None
argmax_pur7 = max(other_pur7, key=other_pur7.get) if other_pur7 else None

print(f"\n▶ χ₇ 요약 [새 방법]:", flush=True)
print(f"  Λ(s,χ₇): σ=0.5→{j_half_lam7}  타σ 최대={max_other_lam7}  (σ={argmax_lam7})", flush=True)
print(f"  L(s,χ₇): σ=0.5→{j_half_pur7}  타σ 최대={max_other_pur7}  (σ={argmax_pur7})", flush=True)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 진단: 구 방법 재현 (σ=0.5 + 최대점프 σ에서만)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

print("\n" + "="*70, flush=True)
print("진단: 구 방법 (threshold=2.5, no branch cut) 재현 — χ₇", flush=True)
print("="*70, flush=True)

diag_sigmas = [0.50]
if argmax_lam7 is not None:
    diag_sigmas.append(float(argmax_lam7))
diag_sigmas = sorted(set(diag_sigmas))

for sigma in diag_sigmas:
    j_old_lam = count_jumps_old(sigma, ci7, completed_L)
    j_old_pur = count_jumps_old(sigma, ci7, pure_L)
    print(f"  σ={sigma:.2f}: 구방법 Λ={j_old_lam}  L={j_old_pur}", flush=True)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실험 2: χ₈, χ₁₁ 확인 (σ=0.5 + χ₇에서 최대점프가 나타난 σ)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

print("\n" + "="*70, flush=True)
print("실험 2: χ₈, χ₁₁ 확인 [새 방법]", flush=True)
print("="*70, flush=True)

# 확인할 σ: 0.5 + χ₇ Λ에서 최대 점프가 나타난 σ
check_sigmas = sorted(set([0.50] + ([float(argmax_lam7)] if argmax_lam7 is not None else [])))
print(f"확인 σ: {check_sigmas}", flush=True)

exp2_results = {}
for char_name in ['χ₈ (mod 8)', 'χ₁₁ (mod 11)']:
    ci = CHARS_HIGH[char_name]
    print(f"\n  [{char_name}]", flush=True)
    r = {}
    for sigma in check_sigmas:
        j_lam = count_jumps_new(sigma, ci, completed_L)
        j_pur = count_jumps_new(sigma, ci, pure_L)
        r[sigma] = {'lambda': j_lam, 'pure': j_pur}
        marker = "★" if abs(sigma - 0.5) < 0.01 else " "
        print(f"    {marker} σ={sigma:.2f}:  Λ={j_lam:3d}  L={j_pur:3d}", flush=True)
    exp2_results[char_name] = r


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 판정 및 결과 저장
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

print("\n" + "="*70, flush=True)
print("최종 판정", flush=True)
print("="*70, flush=True)

# Gamma 영향 평가
gamma_explains = (max_other_lam7 > 0 and max_other_pur7 == 0)
gamma_partial  = (max_other_lam7 > 0 and 0 < max_other_pur7 < max_other_lam7)
gamma_none     = (max_other_pur7 >= max_other_lam7)

if gamma_explains:
    verdict = "양성: Gamma 아티팩트 규명 → Conj 1은 L 레벨에서 성립"
elif gamma_partial:
    verdict = f"중립: Gamma 부분 설명 (Λ={max_other_lam7} → L={max_other_pur7}), 추가 분석 필요"
else:
    verdict = f"음성: Gamma 제거 후에도 타σ 점프 잔존 (Λ={max_other_lam7}, L={max_other_pur7})"

print(f"\n  χ₇ Λ 타σ 최대: {max_other_lam7}", flush=True)
print(f"  χ₇ L 타σ 최대: {max_other_pur7}", flush=True)
print(f"  판정: {verdict}", flush=True)

# 결과 파일 저장
out_path = os.path.expanduser(
    '~/Desktop/gdl_unified/results/dirichlet_offcritical_anatomy.txt')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

with open(out_path, 'w') as f:
    f.write("=" * 72 + "\n")
    f.write("디리클레 off-critical 해부: Λ(s,χ) vs L(s,χ) 위상점프 비교\n")
    f.write(f"날짜: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    f.write(f"정밀도: {mpmath.mp.dps} dps\n")
    f.write(f"구간: t∈[{T_MIN},{T_MAX}], {N_POINTS}점/σ\n")
    f.write(f"방법: branch cut 보정 있음, threshold=π/2 (ζ off-critical 동일)\n")
    f.write("=" * 72 + "\n\n")

    f.write("실험 1: χ₇(mod 7) 전체 σ 프로파일\n")
    f.write(f"{'σ':>6}  {'Λ(s,χ₇)':>10}  {'L(s,χ₇)':>10}  비고\n")
    f.write("-" * 40 + "\n")
    for sigma in sigmas_17:
        marker = "★" if abs(sigma - 0.5) < 0.01 else " "
        f.write(f"{marker}{sigma:>5.2f}  {jumps_lambda_7[sigma]:>10}  "
                f"{jumps_pure_7[sigma]:>10}\n")

    f.write(f"\n▶ 요약:\n")
    f.write(f"  Λ: σ=0.5→{j_half_lam7},  타σ 최대={max_other_lam7}  (σ={argmax_lam7})\n")
    f.write(f"  L: σ=0.5→{j_half_pur7},  타σ 최대={max_other_pur7}  (σ={argmax_pur7})\n")

    f.write("\n" + "=" * 72 + "\n")
    f.write("진단: 구 방법 재현 (threshold=2.5, no branch cut) — χ₇\n")
    f.write(f"{'σ':>6}  {'구방법 Λ':>10}  {'구방법 L':>10}\n")
    f.write("-" * 35 + "\n")
    for sigma in diag_sigmas:
        j_old_lam = count_jumps_old(sigma, ci7, completed_L)
        j_old_pur = count_jumps_old(sigma, ci7, pure_L)
        f.write(f"{sigma:>6.2f}  {j_old_lam:>10}  {j_old_pur:>10}\n")

    f.write("\n" + "=" * 72 + "\n")
    f.write("실험 2: χ₈(mod 8), χ₁₁(mod 11) 확인\n")
    f.write(f"확인 σ: {check_sigmas}\n\n")
    for char_name, r in exp2_results.items():
        f.write(f"{char_name}:\n")
        for sigma, vals in r.items():
            marker = "★" if abs(sigma - 0.5) < 0.01 else " "
            f.write(f"  {marker}σ={sigma:.2f}:  Λ={vals['lambda']:3d}  L={vals['pure']:3d}\n")
        f.write("\n")

    f.write("=" * 72 + "\n")
    f.write("최종 판정\n")
    f.write("=" * 72 + "\n")
    f.write(f"\n{verdict}\n\n")
    f.write(f"  χ₇ Λ(s,χ) 타σ 최대: {max_other_lam7}  (기존 high_conductor 8)\n")
    f.write(f"  χ₇ L(s,χ) 타σ 최대: {max_other_pur7}\n")
    f.write(f"  σ=0.5 비교: Λ={j_half_lam7}, L={j_half_pur7}\n\n")
    f.write(f"소요 시간: 실험1 {elapsed1:.0f}초\n")

print(f"\n결과 저장: {out_path}", flush=True)
print(f"[완료] {datetime.now().strftime('%H:%M:%S')}", flush=True)
