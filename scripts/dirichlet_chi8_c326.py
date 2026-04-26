"""
=============================================================================
[Project RDL] C-326 — χ mod 8 (φ(8)=4) 비자명 지표 4성질 완전 검증
=============================================================================
목표:
  - q=8 (합성수 conductor)의 3개 비자명 지표에 대해 4성질 검증
  - 원시(primitive) vs 유도(induced) 지표 명시 구분
  - Im-스캔 아티팩트 수정: find_zeros_dirichlet 호출 후 범위 필터 적용
  - 결과 저장: results/dirichlet_chi8_c326.txt

q=8 지표 구조:
  - φ(8) = 4, 비자명 지표 3개
  - (Z/8Z)* = {1,3,5,7} ≅ Z/2 × Z/2 (Klein four-group)
  - 모든 지표가 실수 (order 1 or 2) → Im-스캔 자동 비활성화
  - χ₁: conductor=4 (유도, induced from χ mod 4)
  - χ₂: conductor=8 (원시, primitive)
  - χ₃: conductor=8 (원시, primitive, Kronecker symbol (2/n))

4성질:
  (a) σ-유일성: E(σ=0.5)/E(σ=0.3) 비율 (임계선 집중도)
  (b) κδ² 정규화: κ(σ=0.5+δ, t=γ)·δ² ≈ constant (단순 극 검증)
  (c) monodromy ±π: 영점 통과 시 위상 점프 ≈ ±π
  (d) detect율: κ 피크 → 영점 확인 비율 ≥ 90%

주요 변경 (C-326):
  - t 범위: [10, 200] (C-324 [10,70] 대비 확대, 영점 ≥30개 목표)
  - dps=80 (t>100 구간 포함 → 체크리스트 준수)
  - 범위 필터: find_zeros_dirichlet 후 [T_MIN, T_MAX] 외 영점 제거
  - 에너지 적분: [10, 100] (속도 최적화, 범위 기재)

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
# q=8 지표 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# (Z/8Z)* = {1, 3, 5, 7} ≅ Z/2 × Z/2 (Klein four-group)
# 동형사상 기저: 3 → (1,0), 5 → (0,1)  [3²≡5²≡1 mod 8, 3·5≡7 mod 8]
#
# χ₁ = χ_{10}: χ(1)=1, χ(3)=-1, χ(5)=1, χ(7)=-1
#   depends only on n mod 4 → conductor=4 (induced from χ mod 4)
#   χ(-1) = χ(7) = -1 → 홀수(a=1)
#
# χ₂ = χ_{01}: χ(1)=1, χ(3)=1, χ(5)=-1, χ(7)=-1
#   primitive mod 8  (conductor=8)
#   χ(-1) = χ(7) = -1 → 홀수(a=1)
#
# χ₃ = χ_{11}: χ(1)=1, χ(3)=-1, χ(5)=-1, χ(7)=1
#   primitive mod 8  (conductor=8), Kronecker (2/n)
#   χ(-1) = χ(7) = 1 → 짝수(a=0)

def _make_chi8_from_table(values_1357):
    """
    n=0..7에 대한 chi 배열 생성.
    values_1357: [χ(1), χ(3), χ(5), χ(7)] (정수 ±1)
    반환: length-8 list, gcd(n,8)>1이면 0
    """
    chi = [mpmath.mpf(0)] * 8
    for n, v in zip([1, 3, 5, 7], values_1357):
        chi[n] = mpmath.mpf(v)
    return chi


CHARACTERS_MOD8 = {
    'chi8_1': {
        'chi': _make_chi8_from_table([1, -1, 1, -1]),
        'q': 8,
        'a': 1,
        'conductor': 4,
        'primitive': False,
        'label': 'chi8_1 (induced, cond=4)',
        'description': 'induced from chi mod 4',
        'values': 'chi(1)=+1, chi(3)=-1, chi(5)=+1, chi(7)=-1',
    },
    'chi8_2': {
        'chi': _make_chi8_from_table([1, 1, -1, -1]),
        'q': 8,
        'a': 1,
        'conductor': 8,
        'primitive': True,
        'label': 'chi8_2 (primitive, cond=8)',
        'description': 'primitive mod 8',
        'values': 'chi(1)=+1, chi(3)=+1, chi(5)=-1, chi(7)=-1',
    },
    'chi8_3': {
        'chi': _make_chi8_from_table([1, -1, -1, 1]),
        'q': 8,
        'a': 0,
        'conductor': 8,
        'primitive': True,
        'label': 'chi8_3 (primitive, cond=8)',
        'description': 'primitive mod 8, Kronecker symbol (2/n)',
        'values': 'chi(1)=+1, chi(3)=-1, chi(5)=-1, chi(7)=+1',
    },
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
# 유틸리티 함수
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
    '~/Desktop/gdl_unified/results/dirichlet_chi8_c326.txt'
)


def main():
    t_start_all = time.time()
    lines = []

    def emit(s=""):
        print(s, flush=True)
        lines.append(s)

    emit("=" * 72)
    emit("[Project RDL] C-326 --- chi mod 8 비자명 지표 4성질 검증")
    emit(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}")
    emit(f"정밀도: {mpmath.mp.dps} dps (t>100 포함)")
    emit(f"t 범위: [{T_MIN}, {T_MAX}] (영점 탐색), 에너지: [{T_MIN}, {T_ENERGY_MAX}]")
    emit(f"n_scan={N_SCAN}, delta={DELTA}, eps_mono={EPS_MONO}")
    emit(f"총 지표 수: {len(CHARACTERS_MOD8)}개 (chi1,chi2,chi3)")
    emit("=" * 72)
    emit()

    emit("q=8 지표 구조:")
    emit("  (Z/8Z)* = {1,3,5,7} = Z/2 x Z/2 (Klein four-group)")
    emit("  모든 지표가 실수 (order <= 2) --- Im-스캔 자동 비활성화")
    emit()
    for key, ci in CHARACTERS_MOD8.items():
        prim_str = "원시(primitive)" if ci['primitive'] else "유도(induced)"
        parity = "홀수(a=1)" if ci['a'] == 1 else "짝수(a=0)"
        emit(f"  {ci['label']}")
        emit(f"    conductor={ci['conductor']}, {prim_str}, {parity}")
        emit(f"    {ci['values']}")
    emit()

    emit("Im-스캔 수정 (C-326 핵심):")
    emit("  find_zeros_dirichlet 호출 후 범위 필터:")
    emit(f"    zeros = [t for t in zeros_raw if {T_MIN} <= t <= {T_MAX}]")
    emit("  q=8 모든 지표는 실수 -> is_real_char=True -> Im-스캔 이미 비활성화")
    emit("  범위 필터는 findroot bracket 이탈 방지 추가 안전장치")
    emit()

    summary = {}

    for key, char_info in CHARACTERS_MOD8.items():
        emit("=" * 72)
        prim_str = "원시" if char_info['primitive'] else "유도(induced, cond=4)"
        emit(f"  [{char_info['label']}]  a={char_info['a']}, {prim_str}")
        emit(f"  {char_info['values']}")
        emit("=" * 72)
        t_start = time.time()

        # ── 파트 A: 영점 탐색 ──
        emit(f"\n[파트 A] 영점 탐색 (t in [{T_MIN},{T_MAX}], n_scan={N_SCAN})")
        zeros_raw = find_zeros_dirichlet(char_info, T_MIN, T_MAX, N_SCAN)

        # ★ C-326 핵심 수정: 범위 필터로 spurious 영점 제거
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
            emit("  ⚠️ 영점 0개 --- 이 지표는 건너뜀")
            summary[key] = {'n_zeros': 0, 'status': 'SKIP', 'label': char_info['label'],
                            'conductor': char_info['conductor'], 'primitive': char_info['primitive']}
            continue
        elif n_zeros < 10:
            emit(f"  ⚠️ 영점 {n_zeros}개 (너무 적음)")

        # 첫 10개 출력
        for i, z in enumerate(zeros[:10], 1):
            emit(f"    #{i:3d}: t = {z:.8f}")
        if n_zeros > 10:
            emit(f"    ... (이하 {n_zeros - 10}개 생략, 마지막: t={zeros[-1]:.4f})")
        emit(f"  성공 기준 (>=30): {'PASS' if n_zeros >= 30 else 'WARN'}")

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
            'label':         char_info['label'],
            'conductor':     char_info['conductor'],
            'primitive':     char_info['primitive'],
            'a':             char_info['a'],
            'n_zeros':       n_zeros,
            'n_removed':     n_removed,
            'kd2_median':    kd2_median,
            'kd2_pass':      kd2_pass,
            'mean_mono_dev': mean_dev,
            'mono_pass_rate': pass_rate,
            'mono_pass':     mono_pass,
            'detect_pct':    detect_pct,
            'detect_pass':   detect_pass,
            'energy_r03':    r03,
            'energy_r07':    r07,
            'energy_pass':   energy_pass,
            'status':        'OK',
        }

    # ─── 종합 요약 ───────────────────────────────────────────────────────
    emit()
    emit("=" * 72)
    emit("종합 요약 --- chi mod 8 (3개 비자명 지표)")
    emit("=" * 72)
    emit()

    hdr = f"  {'지표':<28} {'cond':>5} {'원시':>6} {'영점':>5} {'kd2med':>8} {'monodev':>9} {'detect%':>8} {'E비':>7}"
    emit(hdr)
    emit("  " + "-" * 80)

    all_4pass = 0
    all_ok    = 0
    for key, res in summary.items():
        if res.get('status') == 'SKIP':
            emit(f"  {res['label']:<28} SKIP")
            continue
        all_ok += 1
        prim_sym = "Y" if res['primitive'] else "N"
        nz   = res['n_zeros']
        kd2s = f"{res['kd2_median']:.4f}" if not math.isnan(res['kd2_median']) else "N/A"
        mdev = f"{res['mean_mono_dev']:.5f}" if not math.isnan(res['mean_mono_dev']) else "N/A"
        det  = f"{res['detect_pct']:.1f}%"
        er   = f"{res['energy_r03']:.1f}x"
        p4   = (res['kd2_pass'] and res['mono_pass'] and res['detect_pass'] and res['energy_pass'])
        if p4: all_4pass += 1
        st4  = "ALL-PASS" if p4 else "PARTIAL"
        emit(f"  {res['label']:<28} {res['conductor']:>5} {prim_sym:>6} {nz:>5}"
             f" {kd2s:>8} {mdev:>9} {det:>8} {er:>7}  {st4}")

    emit()
    emit("  [4성질 상세]")
    for key, res in summary.items():
        if res.get('status') != 'OK': continue
        prim = "원시" if res['primitive'] else "유도(cond=4)"
        emit(f"  {res['label']} [{prim}]")
        kd2s = f"{res['kd2_median']:.4f}" if not math.isnan(res['kd2_median']) else "N/A"
        emit(f"    (a) sigma-유일성: {'PASS' if res['energy_pass'] else 'FAIL'}  E비={res['energy_r03']:.1f}x")
        emit(f"    (b) kd2 ~ 1:     {'PASS' if res['kd2_pass']    else 'FAIL'}  kd2중앙={kd2s}")
        mdev_str = f"{res['mean_mono_dev']:.5f}" if not math.isnan(res['mean_mono_dev']) else "N/A"
        emit(f"    (c) mono = +-pi: {'PASS' if res['mono_pass']   else 'FAIL'}  |dev|avg={mdev_str}")
        emit(f"    (d) detect>=90%: {'PASS' if res['detect_pass'] else 'FAIL'}  {res['detect_pct']:.1f}%")

    emit()
    emit("  [원시/유도 비교 고찰]")
    emit("  * 유도 지표 chi8_1 (cond=4): L(s,chi_ind) = L(s,chi_prim_mod4) x (오일러 인자)")
    emit("    영점 위치 = chi mod 4 원시 지표의 영점. 프레임워크 동일 적용 여부 확인.")
    emit("  * 원시 지표 chi8_2,3 (cond=8): 완전히 새로운 conductor에서 4성질 검증.")
    emit("  * 모든 지표 실수(real) -> monodromy eps-기반 방식 적용 (Im-스캔 불필요)")

    emit()
    emit("  [기존 결과 비교]")
    emit(f"  {'q':>4} {'phi(q)':>7} {'비자명':>7} {'conductor유형':>16} {'ALL PASS':>12}")
    emit(f"  {'-'*52}")
    emit(f"  {'3':>4} {'2':>7} {'1':>7} {'소수':>16} {'1/1':>12}  [기존]")
    emit(f"  {'4':>4} {'2':>7} {'1':>7} {'합성(cond=4)':>16} {'1/1':>12}  [기존]")
    emit(f"  {'5':>4} {'4':>7} {'2':>7} {'소수':>16} {'2/2':>12}  [기존]")
    emit(f"  {'7':>4} {'6':>7} {'5':>7} {'소수':>16} {'4/5(arti)':>12}  [C-324]")
    emit(f"  {'8':>4} {'4':>7} {'3':>7} {'혼합(4,8)':>16} {str(all_4pass)+'/'+str(all_ok):>12}  [C-326]")

    emit()
    n_spurious = sum(r.get('n_removed', 0) for r in summary.values() if r.get('status') == 'OK')
    emit(f"  Im-스캔 수정 효과: spurious 제거 {n_spurious}개 (전체)")

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
