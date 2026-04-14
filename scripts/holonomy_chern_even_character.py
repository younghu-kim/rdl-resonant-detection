"""
=============================================================================
[Project RDL] Chern 수 — 짝수 디리클레 지표 검증 (결과 #23 후보)
=============================================================================
사이클 33 / 수학자 지시 2026-04-14 09:05

목적:
  결과 #22 (홀수 지표 χ₃, χ₄)에 이어 짝수 디리클레 지표에서 Chern 수 검증.

  핵심:
  1. Legendre symbol (·/5): χ₅ = [0,1,-1,-1,1] mod 5
     - χ₅(-1) = χ₅(4) = 1 → 짝수 지표 (a=0)
     - Γ 인자: Γ(s/2) (a=0용)
     - L(1/2+it, χ₅)는 실수 → Re(L)=0 방식 유효
  2. 방식 A (Re(L(1/2+it))=0 부호변화) vs 방식 B (|L| 최소화) 교차 검증
     - 짝수 지표에서 두 방식이 같은 영점 목록 → 방법론 보편성 확립
  3. 직사각형 위상추적: n=5,10,15,20
  4. (선택) mod 8 짝수 지표 [0,1,0,-1,0,-1,0,1]: n=5,10

수학적 근거:
  - 짝수 실수 원시 지표: L(s,χ)가 s=1/2에서 대칭 → L(1/2+it)∈ℝ
  - Λ(s,χ) = (q/π)^(s/2) Γ(s/2) L(s,χ)  [a=0의 경우]
    log Λ = (s/2)log(q/π) + loggamma(s/2) + log L
  - 홀수 지표(a=1): Γ((s+1)/2) 사용. 짝수(a=0): Γ(s/2) 사용.
  - 차이: 체크리스트 "짝수 디리클레 지표: a=0, Γ(s/2) 사용"

파라미터:
  dps=100, 수평분할=64, σ∈[0.3,0.7], t∈[1,60], 스캔=6000점

출력:
  results/holonomy_chern_even_character.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT_TXT = os.path.join(BASE_DIR, "results", "holonomy_chern_even_character.txt")
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# ─── 설정 ─────────────────────────────────────────────────────────────────────
DPS = 100
mpmath.mp.dps = DPS

SIGMA_MIN = mpmath.mpf('0.3')
SIGMA_MAX = mpmath.mpf('0.7')
N_HORIZ   = 64

# 짝수 원시 지표 정의
# (·/5) = Legendre symbol mod 5: QR={1,4}, NR={2,3}
CHI5_LIST = [0, 1, -1, -1, 1]   # chi(0)=0, chi(1)=1, chi(2)=-1, chi(3)=-1, chi(4)=1
CHI5_Q, CHI5_A = 5, 0            # a=0 → 짝수, Γ(s/2)

# (선택) mod 8 짝수 지표: chi(1)=1, chi(3)=-1, chi(5)=-1, chi(7)=1
CHI8_LIST = [0, 1, 0, -1, 0, -1, 0, 1]  # chi(-1)=chi(7)=1 → 짝수
CHI8_Q, CHI8_A = 8, 0            # a=0 → 짝수, Γ(s/2)

# ─── 로그 ─────────────────────────────────────────────────────────────────────
log_lines = []
def log(msg=""):
    print(msg, flush=True)
    log_lines.append(msg)

def save_all():
    with open(OUTPUT_TXT, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))


# ═══════════════════════════════════════════════════════════════════════════════
# Λ(s,χ) 위상 계산 (짝수 지표: a=0, Γ(s/2))
# ═══════════════════════════════════════════════════════════════════════════════

def log_lambda_imag_even(s, chi_list, q):
    """
    Im(log Λ(s, χ)) — 짝수 지표 (a=0) 전용.
    Λ(s,χ) = (q/π)^(s/2) Γ(s/2) L(s,χ)
    log Λ = (s/2)log(q/π) + loggamma(s/2) + log L(s,χ)

    L(s,χ) 영점 근처이면 None 반환.
    """
    L_val = mpmath.dirichlet(s, chi_list)
    if abs(L_val) < mpmath.mpf(10) ** (-DPS + 15):
        return None
    lL = (
        (s / 2) * mpmath.log(mpmath.mpf(q) / mpmath.pi)
        + mpmath.loggamma(s / 2)
        + mpmath.log(L_val)
    )
    return float(lL.imag)


def log_lambda_imag_odd(s, chi_list, q):
    """
    Im(log Λ(s, χ)) — 홀수 지표 (a=1) 전용. 비교용으로 구현.
    Λ(s,χ) = (q/π)^(s/2) Γ((s+1)/2) L(s,χ)
    """
    L_val = mpmath.dirichlet(s, chi_list)
    if abs(L_val) < mpmath.mpf(10) ** (-DPS + 15):
        return None
    lL = (
        (s / 2) * mpmath.log(mpmath.mpf(q) / mpmath.pi)
        + mpmath.loggamma((s + mpmath.mpf(1)) / 2)
        + mpmath.log(L_val)
    )
    return float(lL.imag)


# ═══════════════════════════════════════════════════════════════════════════════
# 위상 추적법 (직사각형 윤곽 적분)
# ═══════════════════════════════════════════════════════════════════════════════

def phase_track_contour(sigma_lo, sigma_hi, t_lo, t_hi,
                        n_horiz, n_vert, log_phase_fn):
    """
    위상 추적법 — ∮_C d(arg Λ) / (2π) = N (포함된 영점 수)
    직사각형 반시계 방향: 하변 → 우변 → 상변 → 좌변
    """
    sig_lo = float(sigma_lo)
    sig_hi = float(sigma_hi)
    t0 = float(t_lo)
    t1 = float(t_hi)

    def make_path():
        pts = []
        for k in range(n_horiz + 1):
            sig = sig_lo + k * (sig_hi - sig_lo) / n_horiz
            pts.append(mpmath.mpc(sig, t0))
        for k in range(1, n_vert + 1):
            t = t0 + k * (t1 - t0) / n_vert
            pts.append(mpmath.mpc(sig_hi, t))
        for k in range(1, n_horiz + 1):
            sig = sig_hi - k * (sig_hi - sig_lo) / n_horiz
            pts.append(mpmath.mpc(sig, t1))
        for k in range(1, n_vert):
            t = t1 - k * (t1 - t0) / n_vert
            pts.append(mpmath.mpc(sig_lo, t))
        return pts

    path = make_path()
    total_phase = 0.0
    phi_prev = None
    skip_count = 0

    for s in path:
        phi = log_phase_fn(s)
        if phi is None:
            skip_count += 1
            continue
        if phi_prev is None:
            phi_prev = phi
            continue
        diff = phi - phi_prev
        diff -= 2 * np.pi * round(diff / (2 * np.pi))
        total_phase += diff
        phi_prev = phi

    if skip_count > 0:
        print(f"  WARNING: {skip_count}개 점 건너뜀 (L 영점 근처)", flush=True)

    return total_phase / (2 * np.pi)


# ═══════════════════════════════════════════════════════════════════════════════
# 방식 A: Re(L(1/2+it))=0 부호변화 영점 탐색 (짝수 지표 전용)
# ═══════════════════════════════════════════════════════════════════════════════

def find_zeros_method_A(chi_list, t_min=1.0, t_max=60.0, n_scan=6000,
                        label="χ"):
    """
    방식 A: Re(L(1/2+it, χ))=0 부호변화 탐색.
    짝수 원시 지표에서 L(1/2+it)∈ℝ → 이 방식 유효.
    홀수 지표에서는 사용 금지.

    Returns: 영점 t 값 목록 (오름차순)
    """
    orig_dps = mpmath.mp.dps
    mpmath.mp.dps = 60

    ts = np.linspace(t_min, t_max, n_scan)
    re_vals = []

    log(f"  {label} 방식A — Re(L(1/2+it)) 스캔 ({n_scan}점, t∈[{t_min:.1f},{t_max:.1f}])...")
    t0 = time.time()

    for t_val in ts:
        try:
            s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_val))
            lv = mpmath.dirichlet(s, chi_list)
            re_vals.append(float(lv.real))
        except Exception as e:
            re_vals.append(0.0)
            print(f"  WARNING: {label} Re(L) t={t_val:.3f}: {e}", flush=True)

    re_vals = np.array(re_vals)
    log(f"  스캔 완료: {time.time()-t0:.1f}초")

    # 부호변화 찾기
    sign_changes = []
    for i in range(len(re_vals) - 1):
        if re_vals[i] * re_vals[i+1] < 0:  # 부호변화
            # 이분법으로 정밀화
            try:
                t_a, t_b = ts[i], ts[i+1]
                t_root = mpmath.findroot(
                    lambda t_var: mpmath.re(mpmath.dirichlet(
                        mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var), chi_list)),
                    (mpmath.mpf(str(t_a)), mpmath.mpf(str(t_b))),
                    solver='illinois'
                )
                t_zero = float(t_root.real)
                # 검증: |L| 작아야 함
                sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero))
                residual = float(abs(mpmath.dirichlet(sv, chi_list)))
                if residual < 0.01:
                    sign_changes.append(t_zero)
                else:
                    print(f"  WARNING: {label} 방식A 위부호변화 t≈{t_zero:.4f}, |L|={residual:.4e}", flush=True)
            except Exception as e:
                # 단순 중점 사용
                t_mid = (ts[i] + ts[i+1]) / 2
                sign_changes.append(t_mid)
                print(f"  WARNING: {label} 방식A findroot 실패 t≈{t_mid:.4f}: {e}", flush=True)

    # 중복 제거 (가까운 점 합치기)
    zeros = []
    for z in sorted(sign_changes):
        if not zeros or abs(z - zeros[-1]) > 0.05:
            zeros.append(z)

    mpmath.mp.dps = orig_dps

    if len(zeros) == 0:
        print(f"⚠️ {label} 방식A 영점 0개 — 탐색 로직 점검 필요", flush=True)

    log(f"  {label} 방식A 영점: {len(zeros)}개")
    if zeros:
        log(f"  첫 10개: {[f'{z:.4f}' for z in zeros[:10]]}")

    return sorted(zeros)


# ═══════════════════════════════════════════════════════════════════════════════
# 방식 B: |L| 최소화 (범용, 사이클 32에서 검증)
# ═══════════════════════════════════════════════════════════════════════════════

def find_zeros_method_B(chi_list, t_min=1.0, t_max=60.0, n_scan=6000,
                        abs_tol=0.15, min_sep=0.05, label="χ"):
    """
    방식 B: |L(1/2+it, χ)| 극소점에서 영점 탐색.
    짝수/홀수 지표 모두에서 사용 가능한 보편적 방식.
    """
    orig_dps = mpmath.mp.dps
    mpmath.mp.dps = 60

    ts = np.linspace(t_min, t_max, n_scan)
    L_abs_vals = []

    log(f"  {label} 방식B — |L(1/2+it)| 스캔 ({n_scan}점, t∈[{t_min:.1f},{t_max:.1f}])...")
    t0 = time.time()

    for t_val in ts:
        try:
            s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_val))
            lv = abs(mpmath.dirichlet(s, chi_list))
            L_abs_vals.append(float(lv))
        except Exception as e:
            L_abs_vals.append(1e10)
            print(f"  WARNING: {label} |L| t={t_val:.3f}: {e}", flush=True)

    L_abs_vals = np.array(L_abs_vals)
    log(f"  스캔 완료: {time.time()-t0:.1f}초, |L| min={L_abs_vals.min():.6f}")

    # 극소점 발견
    minima_idx = []
    for i in range(1, len(L_abs_vals) - 1):
        if (L_abs_vals[i] < L_abs_vals[i-1] and
                L_abs_vals[i] < L_abs_vals[i+1] and
                L_abs_vals[i] < abs_tol):
            minima_idx.append(i)

    log(f"  극소점 후보: {len(minima_idx)}개 (|L|<{abs_tol})")

    zeros = []
    fail_count = 0

    for idx in minima_idx:
        t_guess = ts[idx]
        try:
            t_result = mpmath.findroot(
                lambda t_var: abs(mpmath.dirichlet(
                    mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var), chi_list)),
                mpmath.mpf(str(t_guess))
            )
            t_zero = float(t_result.real)

            # 수렴 검증
            sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero))
            residual = float(abs(mpmath.dirichlet(sv, chi_list)))
            if residual > 1e-6:
                fail_count += 1
                continue

            # 범위 체크
            if t_zero < t_min or t_zero > t_max:
                continue

            # 중복 제거
            if not zeros or abs(t_zero - zeros[-1]) > min_sep:
                zeros.append(t_zero)
            else:
                sv_prev = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(zeros[-1]))
                res_prev = float(abs(mpmath.dirichlet(sv_prev, chi_list)))
                if residual < res_prev:
                    zeros[-1] = t_zero

        except Exception as e:
            fail_count += 1
            print(f"  WARNING: {label} 방식B findroot 실패 t≈{t_guess:.3f}: {e}", flush=True)

    mpmath.mp.dps = orig_dps
    zeros = sorted(zeros)

    if len(zeros) == 0:
        print(f"⚠️ {label} 방식B 영점 0개 — 탐색 로직 점검 필요", flush=True)

    log(f"  {label} 방식B 영점: {len(zeros)}개, findroot 실패: {fail_count}회")
    if zeros:
        log(f"  첫 10개: {[f'{z:.4f}' for z in zeros[:10]]}")

    return sorted(zeros)


# ═══════════════════════════════════════════════════════════════════════════════
# 방식 A ↔ B 교차 검증
# ═══════════════════════════════════════════════════════════════════════════════

def cross_validate_zeros(zeros_A, zeros_B, tol=0.05, label="χ"):
    """
    두 방식의 영점 목록 교차 비교.
    각 방식A 영점에 대해 방식B에서 tol 이내 짝이 있으면 일치.

    Returns: (n_match_AB, n_only_A, n_only_B, match_rate)
    """
    log(f"\n  {label} 방식A↔B 교차 검증 (허용 오차: {tol})")

    matched_A = set()
    matched_B = set()

    for i, za in enumerate(zeros_A):
        for j, zb in enumerate(zeros_B):
            if j not in matched_B and abs(za - zb) < tol:
                matched_A.add(i)
                matched_B.add(j)
                break

    n_match = len(matched_A)
    n_only_A = len(zeros_A) - n_match
    n_only_B = len(zeros_B) - n_match
    n_total = max(len(zeros_A), len(zeros_B))
    match_rate = n_match / n_total * 100 if n_total > 0 else 0.0

    log(f"  방식A: {len(zeros_A)}개, 방식B: {len(zeros_B)}개")
    log(f"  일치: {n_match}쌍, A전용: {n_only_A}개, B전용: {n_only_B}개")
    log(f"  일치율: {match_rate:.1f}%")

    # 불일치 영점 표시
    if n_only_A > 0:
        only_A_zeros = [zeros_A[i] for i in range(len(zeros_A)) if i not in matched_A]
        log(f"  A전용 영점: {[f'{z:.4f}' for z in only_A_zeros]}")
    if n_only_B > 0:
        only_B_zeros = [zeros_B[j] for j in range(len(zeros_B)) if j not in matched_B]
        log(f"  B전용 영점: {[f'{z:.4f}' for z in only_B_zeros]}")

    # 일치 쌍 목록 (상세)
    log(f"  {'방식A':>12}  {'방식B':>12}  {'|차이|':>8}")
    count = 0
    for i, za in enumerate(zeros_A):
        if i in matched_A:
            # 짝 찾기
            best_j, best_diff = -1, 999.
            for j, zb in enumerate(zeros_B):
                if abs(za - zb) < tol and abs(za - zb) < best_diff:
                    best_j, best_diff = j, abs(za - zb)
            if best_j >= 0:
                log(f"  {za:>12.6f}  {zeros_B[best_j]:>12.6f}  {best_diff:>8.6f}")
                count += 1
                if count >= 15:
                    log(f"  ... (이하 생략)")
                    break

    return n_match, n_only_A, n_only_B, match_rate


# ═══════════════════════════════════════════════════════════════════════════════
# 직사각형 위상추적 검증
# ═══════════════════════════════════════════════════════════════════════════════

def run_rectangle_verification(zeros_list, chi_list, q, a, label,
                                t1_start=0.5, n_list=None):
    """
    직사각형 위상추적 검증.
    n_list의 각 n에 대해 t₂ = n번째와 n+1번째 영점 사이 중간점.

    a: 0 (짝수 지표) 또는 1 (홀수 지표)
    """
    if n_list is None:
        n_list = [5, 10, 15, 20]

    T1 = mpmath.mpf(str(t1_start))

    if a == 0:
        def phase_fn(s):
            return log_lambda_imag_even(s, chi_list, q)
    else:
        def phase_fn(s):
            return log_lambda_imag_odd(s, chi_list, q)

    log(f"\n  {label} 직사각형 위상추적 (σ∈[0.3,0.7], t∈[{t1_start:.1f},t₂], a={a})")
    log(f"  {'n':>5}  {'t₂':>10}  {'N_actual':>10}  {'N_ctr':>10}  {'|차이|':>10}  {'n_vert':>8}  판정")
    log("-" * 72)

    results = []
    n_pass = 0
    n_valid = 0

    for n_t in n_list:
        if n_t >= len(zeros_list):
            log(f"  n={n_t}: 영점 부족 ({len(zeros_list)}개)")
            continue

        # t₂ = n번째 영점과 n+1번째 영점 사이 중간점
        t2_f = (zeros_list[n_t - 1] + zeros_list[n_t]) / 2.0
        t2 = mpmath.mpf(str(t2_f))
        N_actual = sum(1 for z in zeros_list if z < t2_f)
        n_vert = max(64, int((t2_f - t1_start) * 3))

        t0_r = time.time()
        try:
            N_ctr = phase_track_contour(
                SIGMA_MIN, SIGMA_MAX, T1, t2,
                N_HORIZ, n_vert, phase_fn
            )
            elapsed = time.time() - t0_r
            diff = abs(N_ctr - N_actual)
            ok = "✅" if diff < 0.1 else "❌"
            if diff < 0.1:
                n_pass += 1
            n_valid += 1
            log(f"  {n_t:>5}  {t2_f:>10.4f}  {N_actual:>10}  {N_ctr:>10.4f}  {diff:>10.4f}  {n_vert:>8}  {ok} ({elapsed:.1f}s)")
            results.append((n_t, t2_f, N_actual, N_ctr, diff))
        except Exception as e:
            log(f"  n={n_t}: 오류: {e}")
            n_valid += 1

    log(f"  {label} 통과: {n_pass}/{n_valid}")
    return results, n_pass, n_valid


# ═══════════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    log("=" * 72)
    log("[Project RDL] Chern 수 — 짝수 디리클레 지표 검증 (결과 #23 후보)")
    log(f"실행 시각: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"dps={DPS}, 수평분할={N_HORIZ}")
    log("=" * 72)

    log("\n【수학적 배경】")
    log("  짝수 원시 지표: χ(-1)=+1 → a=0 파라미터")
    log("  완전화: Λ(s,χ) = (q/π)^(s/2) Γ(s/2) L(s,χ)  [홀수: Γ((s+1)/2)]")
    log("  실수 짝수 지표: L(1/2+it, χ)∈ℝ → Re(L)=0 방식 유효")
    log("  교차 검증: 방식A(Re=0)와 방식B(|L|최소화)가 같은 영점 → 보편성 확립")
    log("")
    log("【지표】")
    log("  χ₅ = (·/5) = [0,1,-1,-1,1] mod 5  (conductor=5, χ₅(-1)=χ₅(4)=1, 짝수)")
    log("  χ₈ = [0,1,0,-1,0,-1,0,1] mod 8   (conductor=8, χ₈(-1)=χ₈(7)=1, 짝수)")
    log("  cf. 사이클 32: χ₃=[0,1,-1] a=1(홀수), χ₄=[0,1,0,-1] a=1(홀수)")

    # ────────────────────────────────────────────────────────────────────────
    # §1. χ₅ = (·/5) 영점 탐색 — 방식 A
    # ────────────────────────────────────────────────────────────────────────
    log("\n" + "=" * 72)
    log("§1. χ₅ = (·/5) mod 5 — 방식 A: Re(L(1/2+it))=0 부호변화")
    log("=" * 72)

    zeros_chi5_A = find_zeros_method_A(
        CHI5_LIST, t_min=1.0, t_max=60.0, n_scan=6000, label="χ₅"
    )

    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §2. χ₅ 영점 탐색 — 방식 B
    # ────────────────────────────────────────────────────────────────────────
    log("\n" + "=" * 72)
    log("§2. χ₅ = (·/5) mod 5 — 방식 B: |L(1/2+it)| 최소화")
    log("=" * 72)

    zeros_chi5_B = find_zeros_method_B(
        CHI5_LIST, t_min=1.0, t_max=60.0, n_scan=6000,
        abs_tol=0.2, min_sep=0.05, label="χ₅"
    )

    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §3. χ₅ 방식 A ↔ B 교차 검증
    # ────────────────────────────────────────────────────────────────────────
    log("\n" + "=" * 72)
    log("§3. χ₅ 방식 A ↔ B 교차 검증")
    log("=" * 72)

    n_match_chi5, n_only_A_chi5, n_only_B_chi5, match_rate_chi5 = cross_validate_zeros(
        zeros_chi5_A, zeros_chi5_B, tol=0.05, label="χ₅"
    )

    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §4. χ₅ 직사각형 위상추적 검증 (방식 B의 영점 목록 사용)
    # ────────────────────────────────────────────────────────────────────────
    log("\n" + "=" * 72)
    log("§4. χ₅ 직사각형 위상추적 검증 (방식B 영점 기준, n=5,10,15,20)")
    log("=" * 72)

    # 방식B 영점이 더 신뢰할 수 있음 (|L| 직접 최소화)
    # 두 방식 모두 충분한 영점 있으면 방식B 우선 사용
    zeros_chi5_phase = zeros_chi5_B if len(zeros_chi5_B) >= 21 else zeros_chi5_A

    if len(zeros_chi5_phase) < 5:
        log("  ⚠️ χ₅ 영점 5개 미만 — 위상추적 건너뜀.")
        chi5_results = []
        n_pass_chi5, n_valid_chi5 = 0, 0
    else:
        chi5_results, n_pass_chi5, n_valid_chi5 = run_rectangle_verification(
            zeros_chi5_phase, CHI5_LIST, CHI5_Q, CHI5_A,
            label="χ₅", t1_start=0.5, n_list=[5, 10, 15, 20]
        )

    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §5. (선택) χ₈ = mod 8 짝수 지표
    # ────────────────────────────────────────────────────────────────────────
    log("\n" + "=" * 72)
    log("§5. (선택) χ₈ = [0,1,0,-1,0,-1,0,1] mod 8 — n=5,10만 검증")
    log("=" * 72)

    log("\n  χ₈ 방식A 영점 탐색...")
    zeros_chi8_A = find_zeros_method_A(
        CHI8_LIST, t_min=1.0, t_max=60.0, n_scan=6000, label="χ₈"
    )
    log("\n  χ₈ 방식B 영점 탐색...")
    zeros_chi8_B = find_zeros_method_B(
        CHI8_LIST, t_min=1.0, t_max=60.0, n_scan=6000,
        abs_tol=0.2, min_sep=0.05, label="χ₈"
    )
    save_all()

    log("\n  χ₈ 방식A↔B 교차 검증:")
    n_match_chi8, n_only_A_chi8, n_only_B_chi8, match_rate_chi8 = cross_validate_zeros(
        zeros_chi8_A, zeros_chi8_B, tol=0.05, label="χ₈"
    )
    save_all()

    zeros_chi8_phase = zeros_chi8_B if len(zeros_chi8_B) >= 11 else zeros_chi8_A

    if len(zeros_chi8_phase) < 5:
        log("  ⚠️ χ₈ 영점 5개 미만 — 위상추적 건너뜀.")
        chi8_results = []
        n_pass_chi8, n_valid_chi8 = 0, 0
    else:
        chi8_results, n_pass_chi8, n_valid_chi8 = run_rectangle_verification(
            zeros_chi8_phase, CHI8_LIST, CHI8_Q, CHI8_A,
            label="χ₈", t1_start=0.5, n_list=[5, 10]
        )

    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §6. 영점 목록 상세
    # ────────────────────────────────────────────────────────────────────────
    log("\n" + "=" * 72)
    log("§6. 영점 목록 상세")
    log("=" * 72)

    log(f"\n  χ₅ 방식A 영점 (총 {len(zeros_chi5_A)}개, t∈[1,60]):")
    for i, z in enumerate(zeros_chi5_A[:20]):
        log(f"    {i+1:>3}. t = {z:.6f}")
    if len(zeros_chi5_A) > 20:
        log(f"    ... ({len(zeros_chi5_A)-20}개 더)")

    log(f"\n  χ₅ 방식B 영점 (총 {len(zeros_chi5_B)}개, t∈[1,60]):")
    for i, z in enumerate(zeros_chi5_B[:20]):
        log(f"    {i+1:>3}. t = {z:.6f}")
    if len(zeros_chi5_B) > 20:
        log(f"    ... ({len(zeros_chi5_B)-20}개 더)")

    log(f"\n  χ₈ 방식A 영점 (총 {len(zeros_chi8_A)}개):")
    for i, z in enumerate(zeros_chi8_A[:15]):
        log(f"    {i+1:>3}. t = {z:.6f}")

    log(f"\n  χ₈ 방식B 영점 (총 {len(zeros_chi8_B)}개):")
    for i, z in enumerate(zeros_chi8_B[:15]):
        log(f"    {i+1:>3}. t = {z:.6f}")

    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §7. 최종 판정
    # ────────────────────────────────────────────────────────────────────────
    elapsed_total = time.time() - t_start

    log("\n" + "=" * 72)
    log("§7. 최종 판정")
    log("=" * 72)

    log(f"\n【χ₅ 방식A↔B 일치율】: {match_rate_chi5:.1f}% ({n_match_chi5}쌍)")
    log(f"【χ₈ 방식A↔B 일치율】: {match_rate_chi8:.1f}% ({n_match_chi8}쌍)")

    log(f"\n【χ₅ 직사각형 위상추적】: {n_pass_chi5}/{n_valid_chi5}")
    for r in chi5_results:
        n_t, t2_f, N_act, N_ctr, diff = r
        ok = "✅" if diff < 0.1 else "❌"
        log(f"  n={n_t:>3}: t₂={t2_f:.4f}, N_act={N_act}, N_ctr={N_ctr:.4f}, |차이|={diff:.4f} {ok}")

    log(f"\n【χ₈ 직사각형 위상추적 (선택)】: {n_pass_chi8}/{n_valid_chi8}")
    for r in chi8_results:
        n_t, t2_f, N_act, N_ctr, diff = r
        ok = "✅" if diff < 0.1 else "❌"
        log(f"  n={n_t:>3}: t₂={t2_f:.4f}, N_act={N_act}, N_ctr={N_ctr:.4f}, |차이|={diff:.4f} {ok}")

    # 성공 기준: χ₅ 4/4 + 방식A↔B 100% 일치
    chi5_phase_ok = (n_pass_chi5 == n_valid_chi5 and n_valid_chi5 >= 4)
    match_ok = match_rate_chi5 >= 99.0
    chi8_phase_ok = (n_pass_chi8 == n_valid_chi8 and n_valid_chi8 >= 2)

    if chi5_phase_ok and match_ok:
        if chi8_phase_ok:
            verdict = ("★ 완전 양성 — χ₅(4/4) + 방식A↔B 100% + χ₈ 추가 확인.\n"
                       "  짝수/홀수 지표 보편성 완전 확립. 결과 #23 등록.")
        else:
            verdict = ("★ 양성 — χ₅(4/4) + 방식A↔B 100% 일치.\n"
                       "  짝수 지표 Chern 수 확립. 결과 #23 등록.")
    elif chi5_phase_ok and not match_ok:
        verdict = ("부분 양성 — 위상추적 통과하나 방식A≠B.\n"
                   "  → 영점 탐색 방법론 문제: 디버깅 필요.")
    elif not chi5_phase_ok and match_ok:
        verdict = ("부분 양성 — 방식A↔B 일치하나 위상추적 실패.\n"
                   "  → χ₅ 짝수 지표 특수성 or 경로 설계 문제.")
    elif n_valid_chi5 == 0:
        verdict = "미결 — 영점 탐색 실패 (결과 없음)"
    else:
        verdict = "음성 — 위상추적 실패 + 방식 불일치. 디버깅 필요."

    log(f"\n  최종: {verdict}")
    log(f"\n  성공 기준:")
    log(f"    - χ₅ 위상추적 4/4: {'✅' if chi5_phase_ok else '❌'} ({n_pass_chi5}/{n_valid_chi5})")
    log(f"    - 방식A↔B 일치 100%: {'✅' if match_ok else '❌'} ({match_rate_chi5:.1f}%)")
    log(f"    - χ₈ (선택): {'✅' if chi8_phase_ok else '부분/미달'} ({n_pass_chi8}/{n_valid_chi8})")
    log(f"\n  실행 시간: {elapsed_total/60:.1f}분")
    log(f"  결과 파일: {OUTPUT_TXT}")

    log("\n" + "=" * 72)
    log("【수학자 참고 — 새로 확인된 사항】")
    log("=" * 72)
    log("  1. 짝수 지표 a=0 → Λ 감마 인자 = Γ(s/2) (홀수: Γ((s+1)/2))")
    log("  2. 짝수 실수 지표: L(1/2+it)∈ℝ → Re(L)=0 방식 유효")
    log("  3. 방식A(Re=0)↔방식B(|L|) 교차 결과 위에 기록")
    log(f"  4. χ₅ 첫 영점: A={zeros_chi5_A[0]:.4f}  B={zeros_chi5_B[0]:.4f}" if zeros_chi5_A and zeros_chi5_B else "  4. 영점 탐색 실패")
    log(f"  5. χ₈ 첫 영점: A={zeros_chi8_A[0]:.4f}  B={zeros_chi8_B[0]:.4f}" if zeros_chi8_A and zeros_chi8_B else "  5. χ₈ 영점 탐색 실패")

    save_all()
    log("\n완료!")


if __name__ == "__main__":
    main()
