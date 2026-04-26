"""
=============================================================================
[Project RDL] C-328 — χ mod 11 (φ(11)=10) 비자명 지표 4성질 완전 검증
=============================================================================
목표:
  - q=11 (소수 conductor)의 9개 비자명 지표에 대해 4성질 검증
  - 모든 지표가 원시(primitive) — 11은 소수이므로 유도 지표 없음
  - order 2 (Legendre), order 5 (복소), order 10 (복소) 분류별 집계
  - Im-스캔 아티팩트: Legendre(order 2)만 spurious 필터 적용
  - 결과 저장: results/dirichlet_chi11_c328.txt

q=11 지표 구조:
  - φ(11) = 10, 비자명 지표 9개 (모두 원시)
  - (Z/11Z)* ≅ Z/10 (순환군), 생성원 g=2
  - 이산로그 (base 2) mod 11:
      n:  1  2  3  4  5  6  7  8  9  10
      k:  0  1  8  2  4  9  7  3  6   5
  - χ_j(n) = ω^(j * k_n mod 10), ω = e^(2πi/10)
  - order(χ_j) = 10 / gcd(j, 10)
  - 홀짝성: a=1 if χ(-1)=-1, χ_j(-1) = ω^(5j) = (-1)^j
      j=홀수 → a=1(홀), j=짝수 → a=0(짝)

  지표 분류:
    order 2:  χ_5 (Legendre symbol, 실수 ±1)
    order 5:  χ_2, χ_4, χ_6, χ_8 (복소 5th root)
    order 10: χ_1, χ_3, χ_7, χ_9 (복소 10th root)

  켤레쌍: χ_j와 χ_{10-j}가 켤레 → 같은 영점 집합
    (χ_1, χ_9), (χ_2, χ_8), (χ_3, χ_7), (χ_4, χ_6)

4성질:
  (a) σ-유일성: E(σ=0.5)/E(σ=0.3) 비율 (임계선 집중도)
  (b) κδ² 정규화: κ(σ=0.5+δ, t=γ)·δ² ≈ constant (단순 극 검증)
  (c) monodromy ±π: 영점 통과 시 위상 점프 ≈ ±π
  (d) detect율: κ 피크 → 영점 확인 비율 ≥ 90%

주요 변경 (C-328 vs C-326):
  - q: 8 → 11 (소수, 유도 지표 없음)
  - 지표 수: 3 → 9 (복소 포함)
  - 복소 지표: order 5, 10 → mpmath 복소수 연산 (5th/10th root of unity)
  - Im-스캔: χ_5 (Legendre, order 2) 이외 복소 지표는 자동 활성화
  - spurious 필터: χ_5에만 범위 필터 적용 (C-326과 동일 패턴)

=============================================================================
"""

import sys, os, time, math
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

# bundle_utils import — 직접 구현 금지
from bundle_utils import (
    completed_L, connection_dirichlet, curvature_dirichlet,
    find_zeros_dirichlet, evaluate_predictions,
)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 정밀도 설정 — t>100 포함 구간 → dps≥80 (체크리스트)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
mpmath.mp.dps = 80

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# q=11 지표 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# (Z/11Z)* 이산로그 (base g=2):
#   2^0=1, 2^1=2, 2^2=4, 2^3=8, 2^4=5, 2^5=10,
#   2^6=9, 2^7=7, 2^8=3, 2^9=6
#
# n  : 1  2  3  4  5  6  7  8  9  10
# k  : 0  1  8  2  4  9  7  3  6   5
#
# χ_j(n) = ω^(j*k_n mod 10), ω = e^(2πi/10)
# χ_j(0) = 0 (gcd(0,11)=11>1)

# 이산로그 테이블 (base 2, mod 11)
_DLOG11 = {1: 0, 2: 1, 3: 8, 4: 2, 5: 4, 6: 9, 7: 7, 8: 3, 9: 6, 10: 5}

# ω = e^(2πi/10) 10th root of unity — mpmath 고정밀
_OMEGA10 = mpmath.exp(2j * mpmath.pi / 10)

# ω 거듭제곱 미리 계산 (exponent 0..9)
_OMEGA_POW = [mpmath.power(_OMEGA10, k) for k in range(10)]


def _make_chi11(j):
    """
    j ∈ {1,..,9}: χ_j mod 11 지표 값 배열 (길이 11).
    chi[0] = 0 (gcd(0,11)>1), chi[n] = ω^(j*k_n mod 10) for n=1..10
    """
    chi = [mpmath.mpf(0)] * 11
    for n in range(1, 11):
        k = _DLOG11[n]
        exp_idx = (j * k) % 10
        chi[n] = _OMEGA_POW[exp_idx]
    return chi


def _order(j):
    """χ_j의 order = 10 / gcd(j, 10)"""
    from math import gcd
    return 10 // gcd(j, 10)


def _parity(j):
    """
    a=1 if χ_j(-1)=-1 (홀수), a=0 if χ_j(-1)=+1 (짝수).
    χ_j(-1) = χ_j(10) = ω^(j*5 mod 10) = ω^(5j mod 10)
    ω^5 = e^(iπ) = -1, ω^0 = 1
    → (-1)^j: j 홀수 → -1 (a=1), j 짝수 → +1 (a=0)
    """
    return 1 if j % 2 == 1 else 0


def _is_real_chi(j):
    """order 2 (Legendre symbol)만 실수"""
    return _order(j) == 2


def _conj_pair(j):
    """켤레 쌍: χ_j ↔ χ_{10-j}"""
    pair = 10 - j
    if pair == 10:
        pair = 0  # trivial
    return pair if pair != j else None


# 9개 비자명 지표 구성 (j=1..9)
CHARACTERS_MOD11 = {}
for _j in range(1, 10):
    _ord = _order(_j)
    _a   = _parity(_j)
    _real = _is_real_chi(_j)
    _conj = _conj_pair(_j)

    # 실수 여부 판별용: order 2 → 값이 ±1 (mpf), 나머지 → mpc
    _chi = _make_chi11(_j)

    # order 2이면 정수값 사용 (is_real_char 판별 호환)
    if _real:
        _chi_clean = []
        for v in _chi:
            re_v = float(mpmath.re(v))
            im_v = float(mpmath.im(v))
            if abs(im_v) < 1e-10:
                _chi_clean.append(mpmath.mpf(round(re_v)))
            else:
                _chi_clean.append(v)
    else:
        _chi_clean = _chi

    _prim_label = f"chi11_{_j}"
    _conj_str = f"chi11_{_conj}" if _conj else "self"
    CHARACTERS_MOD11[_prim_label] = {
        'chi':        _chi_clean,
        'q':          11,
        'a':          _a,
        'conductor':  11,   # 소수 → 모두 원시
        'primitive':  True,
        'order':      _ord,
        'label':      f"chi11_{_j} (order={_ord})",
        'description': f"χ_{_j} mod 11, order {_ord}, {'홀(a=1)' if _a else '짝(a=0)'}",
        'values':     f"j={_j}",
        'conj_of':    _conj_str,
        'is_legendre': _real,
    }

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파라미터
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T_MIN = 10.0
T_MAX = 200.0
N_SCAN = 4000          # 영점 탐색 스캔 수 (범위 190 / 4000 ≈ 0.047 간격)
DELTA = 0.03           # kd2 오프셋
EPS_MONO = 0.005       # monodromy eps
ENERGY_N_T = 100       # 에너지 적분 격자 수
T_ENERGY_MAX = 100.0   # 에너지 계산 상한 (속도 최적화)
DETECT_THRESH_FACTOR = 5.0
MAX_ZEROS_KD2 = 30     # kd2 계산 최대 영점 수
MAX_ZEROS_MONO = 20    # monodromy 계산 최대 수
T_DETECT_MAX = 150.0   # detect율 계산 상한


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 유틸리티 함수 (C-326과 동일 구조)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def monodromy_at_zero(char_info, t_zero, eps=EPS_MONO):
    """단순 영점 위상 점프 Δarg(Lambda)"""
    s_plus  = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero + eps))
    s_minus = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero - eps))
    try:
        val_plus  = completed_L(s_plus,  char_info)
        val_minus = completed_L(s_minus, char_info)
        arg_plus  = float(mpmath.arg(val_plus))
        arg_minus = float(mpmath.arg(val_minus))
        delta = arg_plus - arg_minus
        while delta >  math.pi: delta -= 2 * math.pi
        while delta < -math.pi: delta += 2 * math.pi
        return delta
    except Exception as e:
        print(f"  WARNING monodromy t={t_zero:.4f}: {e}", flush=True)
        return float('nan')


def kappa_delta2_normalization(char_info, zeros, delta=DELTA):
    """kd2 = |A(s)|^2 * delta^2 at s = 0.5+delta + i*gamma"""
    kd2_vals = []
    for gamma in zeros:
        s = mpmath.mpf(str(0.5 + delta)) + 1j * mpmath.mpf(str(gamma))
        try:
            kappa = float(abs(connection_dirichlet(s, char_info))**2)
            if not math.isfinite(kappa):
                print(f"  WARNING kd2 NaN/Inf t={gamma:.4f}", flush=True)
                continue
            kd2_vals.append(kappa * delta**2)
        except Exception as e:
            print(f"  WARNING kd2 t={gamma:.4f}: {e}", flush=True)
    return np.array(kd2_vals)


def energy_ratio(char_info, t_min=T_MIN, t_max=T_ENERGY_MAX, n_t=ENERGY_N_T):
    """E(sigma) = integral kappa(sigma+it) dt, numpy 2.0 trapezoid"""
    ts = np.linspace(t_min, t_max, n_t)
    dt = (t_max - t_min) / (n_t - 1)
    energies = {}
    for sigma in [0.3, 0.5, 0.7]:
        vals = []
        for t in ts:
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            try:
                kappa = float(abs(connection_dirichlet(s, char_info))**2)
                vals.append(min(kappa, 1e6) if math.isfinite(kappa) else 1e6)
            except Exception:
                vals.append(0.0)
        energies[sigma] = np.trapezoid(vals, dx=dt)
    e05 = energies[0.5]
    e03 = energies[0.3]
    e07 = energies[0.7]
    r03 = e05 / e03 if e03 > 0 else float('inf')
    r07 = e05 / e07 if e07 > 0 else float('inf')
    return r03, r07, e05, e03, e07


def detect_rate(char_info, zeros, t_detect_max=T_DETECT_MAX, n_scan=2000, tol=0.5):
    """kappa 피크 -> 영점 매칭율 (detect 범위: T_MIN ~ t_detect_max)"""
    zeros_in_range = [z for z in zeros if T_MIN <= z <= t_detect_max]
    ts = np.linspace(T_MIN, t_detect_max, n_scan)
    kappas = []
    for t in ts:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        try:
            k = float(abs(connection_dirichlet(s, char_info))**2)
            kappas.append(min(k, 1e8) if math.isfinite(k) else 1e8)
        except Exception:
            kappas.append(0.0)

    kappas = np.array(kappas)
    threshold = np.median(kappas) * DETECT_THRESH_FACTOR
    peaks = []
    for i in range(1, len(ts) - 1):
        if kappas[i] > threshold and kappas[i] > kappas[i-1] and kappas[i] > kappas[i+1]:
            peaks.append(ts[i])

    if len(peaks) == 0:
        print(f"  WARNING: kappa 피크 0개 — 임계값 낮춤", flush=True)
        threshold = np.median(kappas) * 2.0
        for i in range(1, len(ts) - 1):
            if kappas[i] > threshold and kappas[i] > kappas[i-1] and kappas[i] > kappas[i+1]:
                peaks.append(ts[i])

    peak_arr = np.array(peaks, dtype=float)
    zero_arr = np.array(zeros_in_range, dtype=float)

    if len(zero_arr) == 0:
        return 0.0, 0.0, 0.0, len(peaks), 0

    precision, recall, f1 = evaluate_predictions(peak_arr, zero_arr, tolerance=tol)
    return precision, recall, f1, len(peaks), len(zeros_in_range)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

OUTPUT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/dirichlet_chi11_c328.txt'
)


def main():
    t_start_all = time.time()
    lines = []

    def emit(s=""):
        print(s, flush=True)
        lines.append(s)

    emit("=" * 72)
    emit("[Project RDL] C-328 --- chi mod 11 비자명 지표 4성질 검증")
    emit(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}")
    emit(f"정밀도: {mpmath.mp.dps} dps (t>100 포함)")
    emit(f"t 범위: [{T_MIN}, {T_MAX}] (영점 탐색), 에너지: [{T_MIN}, {T_ENERGY_MAX}]")
    emit(f"n_scan={N_SCAN}, delta={DELTA}, eps_mono={EPS_MONO}")
    emit(f"총 지표 수: {len(CHARACTERS_MOD11)}개 (chi11_1 ~ chi11_9)")
    emit("=" * 72)
    emit()

    emit("q=11 지표 구조:")
    emit("  (Z/11Z)* ≅ Z/10 (순환군), 생성원 g=2")
    emit("  모든 지표 원시(primitive) — 11은 소수")
    emit("  ω = e^(2πi/10), χ_j(n) = ω^(j * log_2(n) mod 10)")
    emit()
    emit("  이산로그 (base 2 mod 11):")
    emit("    n: 1  2  3  4  5  6  7  8  9  10")
    emit("    k: 0  1  8  2  4  9  7  3  6   5")
    emit()
    emit(f"  {'지표':<12} {'order':>6} {'a':>3} {'켤레':>10} {'Im-스캔':>8}")
    emit("  " + "-" * 45)
    for key, ci in CHARACTERS_MOD11.items():
        scan_str = "비활성(실수)" if ci['is_legendre'] else "활성(복소)"
        emit(f"  {ci['label']:<12} {ci['order']:>6} {ci['a']:>3} {ci['conj_of']:>10} {scan_str:>12}")
    emit()

    emit("Im-스캔 정책 (C-328):")
    emit("  - Legendre (order 2, χ_5): 실수 → Im-스캔 비활성 + 범위 필터 (C-326 패턴)")
    emit("  - order 5, 10 복소 지표: Im-스캔 활성 (bundle_utils 자동 판별)")
    emit("  - 모든 지표: 범위 필터 [T_MIN, T_MAX] 적용")
    emit()

    summary = {}

    for key, char_info in CHARACTERS_MOD11.items():
        emit("=" * 72)
        ord_str = f"order={char_info['order']}"
        parity_str = "홀(a=1)" if char_info['a'] else "짝(a=0)"
        emit(f"  [{char_info['label']}]  {ord_str}, {parity_str}, 켤레={char_info['conj_of']}")
        emit(f"  {char_info['description']}")
        leg_str = "★ Legendre symbol (실수, spurious 필터 적용)" if char_info['is_legendre'] else ""
        if leg_str:
            emit(f"  {leg_str}")
        emit("=" * 72)
        t_start = time.time()

        # ── 파트 A: 영점 탐색 ──
        emit(f"\n[파트 A] 영점 탐색 (t in [{T_MIN},{T_MAX}], n_scan={N_SCAN})")
        zeros_raw = find_zeros_dirichlet(char_info, T_MIN, T_MAX, N_SCAN)

        # 범위 필터: spurious 영점 제거 (모든 지표에 적용)
        zeros_in_range = [float(t) for t in zeros_raw if T_MIN <= float(t) <= T_MAX]
        zeros = sorted(set(zeros_in_range))
        n_zeros_raw = len(zeros_raw)
        n_zeros = len(zeros)
        n_removed = n_zeros_raw - n_zeros

        emit(f"  원시 탐색 결과: {n_zeros_raw}개")
        if n_removed > 0:
            removed_list = [float(z) for z in zeros_raw if not (T_MIN <= float(z) <= T_MAX)]
            emit(f"  ★ 범위 필터 제거 (spurious): {n_removed}개 = {[f'{z:.3f}' for z in removed_list]}")
        emit(f"  최종 영점 수: {n_zeros}개")

        if n_zeros == 0:
            emit("  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
            summary[key] = {'n_zeros': 0, 'status': 'SKIP', 'label': char_info['label'],
                            'order': char_info['order'], 'primitive': True,
                            'is_legendre': char_info['is_legendre']}
            continue
        elif n_zeros < 10:
            emit(f"  ⚠️ 영점 {n_zeros}개 (너무 적음)")

        # 첫 10개 출력
        for i, z in enumerate(zeros[:10], 1):
            emit(f"    #{i:3d}: t = {z:.8f}")
        if n_zeros > 10:
            emit(f"    ... (이하 {n_zeros - 10}개 생략, 마지막: t={zeros[-1]:.4f})")
        emit(f"  성공 기준 (>=50개): {'PASS' if n_zeros >= 50 else ('WARN' if n_zeros >= 30 else 'FAIL')}")

        # ── 파트 B: kd2 정규화 ──
        emit(f"\n[파트 B] kd2 정규화 (delta={DELTA}, 최대 {MAX_ZEROS_KD2}개 영점)")
        kd2_vals = kappa_delta2_normalization(char_info, zeros[:MAX_ZEROS_KD2])
        if len(kd2_vals) > 0:
            kd2_median = float(np.median(kd2_vals))
            kd2_mean   = float(np.mean(kd2_vals))
            kd2_std    = float(np.std(kd2_vals))
            emit(f"  kd2 중앙값: {kd2_median:.6f}")
            emit(f"  kd2 평균:   {kd2_mean:.6f} ± {kd2_std:.6f}")
            emit(f"  이론값 (단순 극): ~1.0")
            kd2_pass = (0.5 < kd2_median < 2.0)
            emit(f"  성공 기준 (0.5 < kd2 < 2.0): {'PASS' if kd2_pass else 'FAIL'}")
        else:
            kd2_median = float('nan')
            kd2_pass = False
            emit("  ⚠️ kd2 계산 실패")

        # ── 파트 C: Monodromy ──
        emit(f"\n[파트 C] Monodromy (eps={EPS_MONO}, 첫 {MAX_ZEROS_MONO}개 영점)")
        emit(f"  {'t':>12} {'mono/pi':>10} {'|dev|':>10} {'판정':>6}")
        emit(f"  {'-'*44}")

        mono_vals = []
        mono_devs = []
        for tz in zeros[:MAX_ZEROS_MONO]:
            mono = monodromy_at_zero(char_info, tz)
            if not math.isnan(mono):
                dev = abs(abs(mono) - math.pi)
                mono_vals.append(mono)
                mono_devs.append(dev)
                pass_str = "PASS" if dev < 0.01 else "FAIL"
                sign = "+" if mono > 0 else "-"
                emit(f"  {tz:>12.4f} {sign}{abs(mono)/math.pi:>9.5f} {dev:>10.6f} {pass_str:>6}")

        if len(mono_devs) > 0:
            mean_dev  = float(np.mean(mono_devs))
            pass_rate = sum(1 for d in mono_devs if d < 0.01) / len(mono_devs) * 100
            mono_pass = (mean_dev < 0.01)
            emit(f"\n  평균 |dev|: {mean_dev:.6f}")
            emit(f"  통과율 (|dev|<0.01): {pass_rate:.1f}%")
            emit(f"  성공 기준 (mean_dev < 0.01): {'PASS' if mono_pass else 'FAIL'}")
        else:
            mean_dev  = float('nan')
            pass_rate = 0.0
            mono_pass = False

        # ── 파트 D: Detect율 ──
        emit(f"\n[파트 D] Detect율 (kappa 피크 -> 영점 매칭, tol=0.5, t<={T_DETECT_MAX})")
        precision, recall, f1, n_peaks, n_zeros_detect = detect_rate(char_info, zeros)
        detect_pct = precision * 100
        emit(f"  kappa 피크 수: {n_peaks}")
        emit(f"  대상 영점 수: {n_zeros_detect} (t<={T_DETECT_MAX})")
        emit(f"  precision={precision:.3f}, recall={recall:.3f}, F1={f1:.3f}")
        emit(f"  detect율(precision): {detect_pct:.1f}%")
        detect_pass = (detect_pct >= 90.0)
        emit(f"  성공 기준 (>=90%): {'PASS' if detect_pass else 'FAIL'}")

        # ── 파트 E: sigma-유일성 ──
        emit(f"\n[파트 E] sigma-유일성 (에너지 집중도, t in [{T_MIN},{T_ENERGY_MAX}])")
        r03, r07, e05, e03, e07 = energy_ratio(char_info)
        emit(f"  E(sigma=0.5) = {e05:.4e}")
        emit(f"  E(sigma=0.3) = {e03:.4e}")
        emit(f"  E(sigma=0.7) = {e07:.4e}")
        emit(f"  E(0.5)/E(0.3) = {r03:.1f}x")
        emit(f"  E(0.5)/E(0.7) = {r07:.1f}x")
        energy_pass = (min(r03, r07) >= 5.0)
        emit(f"  성공 기준 (ratio>=5x): {'PASS' if energy_pass else 'FAIL'}")

        t_elapsed = time.time() - t_start
        emit(f"\n  소요 시간: {t_elapsed:.1f}초 ({t_elapsed/60:.1f}분)")

        summary[key] = {
            'label':          char_info['label'],
            'order':          char_info['order'],
            'a':              char_info['a'],
            'primitive':      True,
            'is_legendre':    char_info['is_legendre'],
            'conj_of':        char_info['conj_of'],
            'n_zeros':        n_zeros,
            'n_removed':      n_removed,
            'kd2_median':     kd2_median,
            'kd2_pass':       kd2_pass,
            'mean_mono_dev':  mean_dev,
            'mono_pass_rate': pass_rate,
            'mono_pass':      mono_pass,
            'detect_pct':     detect_pct,
            'detect_pass':    detect_pass,
            'energy_r03':     r03,
            'energy_r07':     r07,
            'energy_pass':    energy_pass,
            'status':         'OK',
        }

    # ─── 종합 요약 ───────────────────────────────────────────────────────
    emit()
    emit("=" * 72)
    emit("종합 요약 --- chi mod 11 (9개 비자명 지표)")
    emit("=" * 72)
    emit()

    hdr = f"  {'지표':<14} {'ord':>4} {'a':>2} {'영점':>5} {'kd2med':>8} {'monodev':>9} {'detect%':>8} {'E비':>7}"
    emit(hdr)
    emit("  " + "-" * 60)

    all_4pass  = 0
    all_ok     = 0
    by_order   = {2: [], 5: [], 10: []}

    for key, res in summary.items():
        if res.get('status') == 'SKIP':
            emit(f"  {res['label']:<14} SKIP")
            continue
        all_ok += 1
        nz   = res['n_zeros']
        kd2s = f"{res['kd2_median']:.4f}" if not math.isnan(res['kd2_median']) else "N/A"
        mdev = f"{res['mean_mono_dev']:.5f}" if not math.isnan(res['mean_mono_dev']) else "N/A"
        det  = f"{res['detect_pct']:.1f}%"
        er   = f"{res['energy_r03']:.1f}x"
        p4   = (res['kd2_pass'] and res['mono_pass'] and res['detect_pass'] and res['energy_pass'])
        if p4:
            all_4pass += 1
        leg_mark = "★" if res['is_legendre'] else " "
        st4  = "ALL-PASS" if p4 else "PARTIAL"
        emit(f"  {res['label']:<14}{leg_mark} {res['order']:>3} {res['a']:>2} {nz:>5}"
             f" {kd2s:>8} {mdev:>9} {det:>8} {er:>7}  {st4}")
        by_order[res['order']].append((key, res, p4))

    emit()
    emit("  ★ = Legendre symbol (order 2, 실수, spurious 필터 적용)")
    emit()

    # order별 집계
    emit("  [order별 집계]")
    for ord_val in [2, 5, 10]:
        items = by_order[ord_val]
        if not items:
            continue
        n_total = len(items)
        n_pass  = sum(1 for _, _, p4 in items if p4)
        avg_zeros = np.mean([r['n_zeros'] for _, r, _ in items]) if items else 0
        avg_kd2   = np.mean([r['kd2_median'] for _, r, _ in items
                             if not math.isnan(r['kd2_median'])]) if items else float('nan')
        avg_erat  = np.mean([r['energy_r03'] for _, r, _ in items]) if items else 0
        emit(f"    order {ord_val:>2}: {n_pass}/{n_total} ALL-PASS"
             f"  avg_zeros={avg_zeros:.0f}  avg_kd2={avg_kd2:.4f}  avg_E비={avg_erat:.1f}x")

    emit()
    emit("  [4성질 상세]")
    for key, res in summary.items():
        if res.get('status') != 'OK':
            continue
        emit(f"  {res['label']}  (켤레={res['conj_of']})")
        kd2s = f"{res['kd2_median']:.4f}" if not math.isnan(res['kd2_median']) else "N/A"
        emit(f"    (a) sigma-유일성: {'PASS' if res['energy_pass'] else 'FAIL'}  E비={res['energy_r03']:.1f}x")
        emit(f"    (b) kd2 ~ 1:     {'PASS' if res['kd2_pass']    else 'FAIL'}  kd2중앙={kd2s}")
        mdev_str = f"{res['mean_mono_dev']:.5f}" if not math.isnan(res['mean_mono_dev']) else "N/A"
        emit(f"    (c) mono = +-pi: {'PASS' if res['mono_pass']   else 'FAIL'}  "
             f"|dev|avg={mdev_str}  pass_rate={res['mono_pass_rate']:.1f}%")
        emit(f"    (d) detect>=90%: {'PASS' if res['detect_pass'] else 'FAIL'}  {res['detect_pct']:.1f}%")

    emit()
    emit("  [B-68 경계 탐사 — E비 conductor 의존성]")
    emit("  q: 3→4→5→7→8→11 E비 추세 (E(0.5)/E(0.3)):")
    emit("  (이번 결과에서 order별 E비를 확인 — 고conductor에서 감소 추세?)")

    # 스퓨리어스 통계
    emit()
    n_spurious = sum(r.get('n_removed', 0) for r in summary.values() if r.get('status') == 'OK')
    n_leg_spur = sum(r.get('n_removed', 0) for r in summary.values()
                     if r.get('status') == 'OK' and r.get('is_legendre'))
    n_cx_spur  = n_spurious - n_leg_spur
    emit(f"  Im-스캔/spurious 통계: 전체 제거 {n_spurious}개")
    emit(f"    - Legendre(order 2):  {n_leg_spur}개 제거")
    emit(f"    - 비-Legendre(복소):  {n_cx_spur}개 제거")
    if n_cx_spur > 0:
        emit(f"    ⚠️  비-Legendre에서 spurious {n_cx_spur}개 — 원인 분석 필요")
    else:
        emit(f"    ✅ 비-Legendre에서 spurious 0개 (정상)")

    emit()
    emit("  [누적 디리클레 검증 현황]")
    emit(f"  {'q':>4} {'phi(q)':>7} {'비자명':>7} {'conductor유형':>16} {'ALL PASS':>12}")
    emit(f"  {'-'*52}")
    emit(f"  {'3':>4} {'2':>7} {'1':>7} {'소수':>16} {'1/1':>12}  [기존]")
    emit(f"  {'4':>4} {'2':>7} {'1':>7} {'합성(cond=4)':>16} {'1/1':>12}  [기존]")
    emit(f"  {'5':>4} {'4':>7} {'2':>7} {'소수':>16} {'2/2':>12}  [기존]")
    emit(f"  {'7':>4} {'6':>7} {'5':>7} {'소수':>16} {'4/5(arti)':>12}  [C-324]")
    emit(f"  {'8':>4} {'4':>7} {'3':>7} {'혼합(4,8)':>16} {'2/2+1조건부':>12}  [C-326]")
    emit(f"  {'11':>4} {'10':>7} {'9':>7} {'소수':>16} {str(all_4pass)+'/'+str(all_ok):>12}  [C-328]")

    total_time = time.time() - t_start_all
    emit()
    emit(f"총 소요 시간: {total_time:.1f}초 ({total_time/60:.1f}분)")
    emit("=" * 72)

    # 결과 저장
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(OUTPUT_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"\n결과 저장 완료: {OUTPUT_PATH}", flush=True)


if __name__ == '__main__':
    main()
