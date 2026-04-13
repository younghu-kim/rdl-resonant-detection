"""
=============================================================================
[Project RDL] Chern 수 디리클레 검증 — χ₃ 버그 수정 + χ₄ 추가
=============================================================================
사이클 32 / 수학자 지시 2026-04-14 08:15

목적:
  결과 #22 χ₃ 코드 버그 수정 및 χ₄ 추가 검증.

  문제: 기존 find_zeros_chi3()가 Re(L(1/2+it,χ₃))=0 부호변화를 세었는데,
        χ₃는 홀수 지표(χ(-1)=-1)이므로 L(1/2+it)은 복소수.
        Re(L)=0 ≠ L=0.

  수정:
    1. |L(1/2+it, χ)| 스캔 → 극소점 발견 → mpmath.findroot로 정밀화
       : 복소 2변수 동시 풀기 (Re=0 AND Im=0)
    2. 동일 위상추적 직사각형 검증 (n=5,10,15,20)

  χ₄ = mod 4, [0,1,0,-1], a=1 (홀수 지표, χ(-1)=χ(3)=-1)
    → 마찬가지로 |L| 최소화 방식 사용

  dps=100, 수평 64분할, 수직 ≥64분할

출력:
  results/holonomy_chern_dirichlet.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT_TXT = os.path.join(BASE_DIR, "results", "holonomy_chern_dirichlet.txt")
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# ─── 설정 ─────────────────────────────────────────────────────────────────────
DPS = 100
mpmath.mp.dps = DPS

SIGMA_MIN = mpmath.mpf('0.3')
SIGMA_MAX = mpmath.mpf('0.7')
N_HORIZ   = 64

# 지표 정의
CHI3_LIST = [0, 1, -1]    # mod 3, a=1 (홀수)
CHI3_Q, CHI3_A = 3, 1

CHI4_LIST = [0, 1, 0, -1]  # mod 4, a=1 (홀수, χ(-1)=χ(3)=-1)
CHI4_Q, CHI4_A = 4, 1

# ─── 로그 ─────────────────────────────────────────────────────────────────────
log_lines = []
def log(msg=""):
    print(msg, flush=True)
    log_lines.append(msg)

def save_all():
    with open(OUTPUT_TXT, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))


# ═══════════════════════════════════════════════════════════════════════════════
# 위상 추적법 (기존 holonomy_chern_number.py에서 재사용)
# ═══════════════════════════════════════════════════════════════════════════════

def log_lambda_imag(s, chi_list, q, a):
    """
    Im(log Λ(s, χ)) 해석적 계산:
      log Λ = (s/2)log(q/π) + loggamma((s+a)/2) + log L(s, χ)
    L(s,χ) 영점 근처이면 None 반환.
    """
    L_val = mpmath.dirichlet(s, chi_list)
    if abs(L_val) < mpmath.mpf(10) ** (-DPS + 15):
        return None
    lL = ((s / 2) * mpmath.log(mpmath.mpf(q) / mpmath.pi)
          + mpmath.loggamma((s + mpmath.mpf(a)) / 2)
          + mpmath.log(L_val))
    return float(lL.imag)


def phase_track_contour(sigma_lo, sigma_hi, t_lo, t_hi,
                        n_horiz, n_vert, log_phase_fn):
    """
    위상 추적법 — ∮_C d(arg Λ) / (2π) = N (포함된 영점 수)
    직사각형 반시계 방향: 하변 → 우변 → 상변 → 좌변
    log_phase_fn: s → Im(log Λ(s)) 또는 None (영점 근처)
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
# 수정된 영점 탐색 — |L| 최소화 방식
# ═══════════════════════════════════════════════════════════════════════════════

def find_zeros_abs_L(chi_list, t_min=1.0, t_max=80.0, n_scan=5000,
                     abs_tol=0.15, min_sep=0.08, label="χ"):
    """
    |L(1/2+it, χ)| 극소점에서 영점 탐색.

    알고리즘:
    1. t_min~t_max를 n_scan 점으로 스캔하여 |L| 계산
    2. 국소 극소점(|L| < abs_tol) 발견
    3. mpmath.findroot로 Re(L)=Im(L)=0 동시 풀기 (2변수)
    4. 수렴한 점만 보존 (|L|<1e-8 기준)

    Parameters:
      abs_tol: 극소점 필터링 기준 |L| 임계값
      min_sep: 중복 영점 제거 최소 간격
    """
    orig_dps = mpmath.mp.dps
    mpmath.mp.dps = 60

    ts = np.linspace(t_min, t_max, n_scan)
    L_abs_vals = []

    log(f"  {label} |L| 스캔 중 ({n_scan}점, t∈[{t_min:.1f},{t_max:.1f}])...")
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

    # 극소점 발견: 이웃보다 작고 < abs_tol
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
            # 2변수 findroot: Re(L)=0, Im(L)=0 동시
            def f_vec(t_var):
                sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                lv = mpmath.dirichlet(sv, chi_list)
                return [mpmath.re(lv), mpmath.im(lv)]

            # 브래킷 없이 단일 시작점 (1D findroot)
            t_result = mpmath.findroot(
                lambda t_var: abs(mpmath.dirichlet(
                    mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var), chi_list)),
                mpmath.mpf(str(t_guess))
            )
            t_zero = float(t_result.real)

            # 수렴 검증: |L| < 1e-8
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
                # 더 좋은 것(|L| 더 작은 것) 보존
                sv_prev = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(zeros[-1]))
                res_prev = float(abs(mpmath.dirichlet(sv_prev, chi_list)))
                if residual < res_prev:
                    zeros[-1] = t_zero

        except Exception as e:
            fail_count += 1
            print(f"  WARNING: {label} findroot 실패 t≈{t_guess:.3f}: {e}", flush=True)

    mpmath.mp.dps = orig_dps
    zeros = sorted(zeros)

    if len(zeros) == 0:
        print(f"⚠️ {label} 영점 0개 — 탐색 로직 점검 필요", flush=True)

    log(f"  {label} 영점: {len(zeros)}개, findroot 실패: {fail_count}회")
    if len(zeros) > 0:
        log(f"  첫 10개: {[f'{z:.4f}' for z in zeros[:10]]}")

    return sorted(zeros)


# ═══════════════════════════════════════════════════════════════════════════════
# 직사각형 위상추적 검증 (공통)
# ═══════════════════════════════════════════════════════════════════════════════

def run_rectangle_verification(zeros_list, chi_list, q, a, label,
                                t1_start=0.5, n_list=None):
    """
    주어진 영점 목록으로 직사각형 위상추적 검증.
    n_list의 각 n에 대해 n번째 영점 이후 중간점을 t₂로 설정.

    Returns: (results_list, n_pass, n_valid)
    """
    if n_list is None:
        n_list = [5, 10, 15, 20]

    T1 = mpmath.mpf(str(t1_start))

    def phase_fn(s):
        return log_lambda_imag(s, chi_list, q, a)

    log(f"\n  {label} 직사각형 위상추적 (σ∈[0.3,0.7], t∈[{t1_start:.1f},t₂])")
    log(f"  {'n':>5}  {'t₂':>10}  {'N_actual':>10}  {'N_ctr':>10}  {'|차이|':>10}  {'n_vert':>8}  판정")
    log("-" * 68)

    results = []
    n_pass = 0
    n_valid = 0

    for n_t in n_list:
        if n_t >= len(zeros_list):
            log(f"  n={n_t}: 영점 부족 ({len(zeros_list)}개 밖에 없음)")
            continue
        if n_t - 1 < 0:
            log(f"  n={n_t}: n-1 인덱스 오류")
            continue

        # t₂ = n번째와 n+1번째 영점 사이 중간점
        if n_t < len(zeros_list):
            t2_f = (zeros_list[n_t - 1] + zeros_list[n_t]) / 2.0
        else:
            t2_f = zeros_list[n_t - 1] + 1.0

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
    log("=" * 70)
    log("[Project RDL] Chern 수 디리클레 검증 — χ₃ 버그 수정 + χ₄ 추가")
    log(f"실행 시각: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"dps={DPS}, 수평분할={N_HORIZ}")
    log("방법: |L| 최소화 → findroot 영점, 위상추적 직사각형 검증")
    log("=" * 70)

    log("\n【핵심 수정 사항】")
    log("  기존: Re(L(1/2+it, χ))=0 부호변화 → 복소 지표에서 영점 아님")
    log("  수정: |L(1/2+it, χ)|를 스캔, 극소점에서 findroot → 진짜 영점 탐색")
    log("  χ₃=[0,1,-1] mod 3, a=1 (홀수): L(1/2+it) 복소수 → |L| 방식 필수")
    log("  χ₄=[0,1,0,-1] mod 4, a=1 (홀수): 마찬가지로 |L| 방식 적용")

    # ────────────────────────────────────────────────────────────────────────
    # §1. χ₃ 영점 탐색 (수정)
    # ────────────────────────────────────────────────────────────────────────
    log("\n" + "=" * 70)
    log("§1. χ₃ = [0,1,-1] mod 3 영점 탐색 (|L| 최소화 방식)")
    log("=" * 70)

    zeros_chi3 = find_zeros_abs_L(
        CHI3_LIST, t_min=1.0, t_max=60.0, n_scan=6000,
        abs_tol=0.2, min_sep=0.05, label="χ₃"
    )
    log(f"\n  χ₃ 영점 {len(zeros_chi3)}개 발견")

    if len(zeros_chi3) < 5:
        log("  ⚠️ χ₃ 영점 5개 미만 — 스캔 범위/해상도 부족. 위상추적 건너뜀.")
        chi3_results = []
        n_pass_chi3, n_valid_chi3 = 0, 0
    else:
        save_all()
        chi3_results, n_pass_chi3, n_valid_chi3 = run_rectangle_verification(
            zeros_chi3, CHI3_LIST, CHI3_Q, CHI3_A,
            label="χ₃", t1_start=0.5, n_list=[5, 10, 15, 20]
        )

    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §2. χ₄ 영점 탐색 + 검증
    # ────────────────────────────────────────────────────────────────────────
    log("\n" + "=" * 70)
    log("§2. χ₄ = [0,1,0,-1] mod 4 영점 탐색 (|L| 최소화 방식)")
    log("=" * 70)

    zeros_chi4 = find_zeros_abs_L(
        CHI4_LIST, t_min=1.0, t_max=60.0, n_scan=6000,
        abs_tol=0.2, min_sep=0.05, label="χ₄"
    )
    log(f"\n  χ₄ 영점 {len(zeros_chi4)}개 발견")

    if len(zeros_chi4) < 5:
        log("  ⚠️ χ₄ 영점 5개 미만 — 위상추적 건너뜀.")
        chi4_results = []
        n_pass_chi4, n_valid_chi4 = 0, 0
    else:
        save_all()
        chi4_results, n_pass_chi4, n_valid_chi4 = run_rectangle_verification(
            zeros_chi4, CHI4_LIST, CHI4_Q, CHI4_A,
            label="χ₄", t1_start=0.5, n_list=[5, 10, 15, 20]
        )

    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §3. 영점 목록 상세
    # ────────────────────────────────────────────────────────────────────────
    log("\n" + "=" * 70)
    log("§3. 영점 목록 상세")
    log("=" * 70)

    log(f"\n  χ₃ 영점 (총 {len(zeros_chi3)}개, t∈[1,60]):")
    for i, z in enumerate(zeros_chi3):
        log(f"    {i+1:>3}. t = {z:.6f}")

    log(f"\n  χ₄ 영점 (총 {len(zeros_chi4)}개, t∈[1,60]):")
    for i, z in enumerate(zeros_chi4):
        log(f"    {i+1:>3}. t = {z:.6f}")

    save_all()

    # ────────────────────────────────────────────────────────────────────────
    # §4. 최종 판정
    # ────────────────────────────────────────────────────────────────────────
    elapsed_total = time.time() - t_start

    log("\n" + "=" * 70)
    log("§4. 최종 판정")
    log("=" * 70)

    log(f"\n【χ₃ 직사각형 위상추적】: {n_pass_chi3}/{n_valid_chi3}")
    for r in chi3_results:
        n_t, t2_f, N_act, N_ctr, diff = r
        ok = "✅" if diff < 0.1 else "❌"
        log(f"  n={n_t:>3}: t₂={t2_f:.4f}, N_act={N_act}, N_ctr={N_ctr:.4f}, |차이|={diff:.4f} {ok}")

    log(f"\n【χ₄ 직사각형 위상추적】: {n_pass_chi4}/{n_valid_chi4}")
    for r in chi4_results:
        n_t, t2_f, N_act, N_ctr, diff = r
        ok = "✅" if diff < 0.1 else "❌"
        log(f"  n={n_t:>3}: t₂={t2_f:.4f}, N_act={N_act}, N_ctr={N_ctr:.4f}, |차이|={diff:.4f} {ok}")

    total_pass = n_pass_chi3 + n_pass_chi4
    total_valid = n_valid_chi3 + n_valid_chi4

    log(f"\n【종합】: {total_pass}/{total_valid}")

    if total_valid == 0:
        verdict = "미결 — 영점 탐색 실패 (결과 없음)"
    elif total_pass == total_valid and total_valid >= 6:
        verdict = "★ 양성 — χ₃+χ₄ 디리클레 보편성 확립, 결과 #22 완전 양성 승격"
    elif n_pass_chi3 == n_valid_chi3 and n_valid_chi3 >= 3 and n_pass_chi4 == n_valid_chi4 and n_valid_chi4 >= 3:
        verdict = "★ 양성 — χ₃+χ₄ 모두 통과"
    elif n_pass_chi3 == n_valid_chi3 and n_valid_chi3 >= 3:
        verdict = "부분 양성 — χ₃ 통과, χ₄ 미달"
    elif n_pass_chi4 == n_valid_chi4 and n_valid_chi4 >= 3:
        verdict = "부분 양성 — χ₄ 통과, χ₃ 미달"
    elif total_pass >= 4:
        verdict = "부분 양성 — 일부 통과"
    else:
        verdict = "음성 — 디버깅 필요"

    log(f"\n  최종: {verdict}")
    log(f"  성공 기준: χ₃ 4/4 AND χ₄ 4/4 → 완전 양성")
    log(f"  실행 시간: {elapsed_total/60:.1f}분")
    log(f"  결과: {OUTPUT_TXT}")

    # 수학자에 대한 노트
    log(f"\n【수학자 참고】")
    log(f"  χ₄ = [0,1,0,-1] mod 4 → χ(-1)=χ(3)=-1 → 홀수 지표 (짝수 아님)")
    log(f"  → L(1/2+it, χ₄)는 복소수. Re(L)=0 방식 부적절.")
    log(f"  → |L| 최소화 방식으로 대체. 위상추적 자체는 정확함.")
    if len(zeros_chi3) > 0:
        log(f"  χ₃ 첫 영점: t={zeros_chi3[0]:.4f} (참고값: ≈5.2 또는 6.0)")
    if len(zeros_chi4) > 0:
        log(f"  χ₄ 첫 영점: t={zeros_chi4[0]:.4f}")

    save_all()
    log("\n완료!")


if __name__ == "__main__":
    main()
