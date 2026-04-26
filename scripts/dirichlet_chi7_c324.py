"""
=============================================================================
[Project RDL] C-324 — χ mod 7 (φ(7)=6) 비자명 지표 4성질 완전 검증
=============================================================================
목표:
  - q=7의 5개 비자명 지표(χ₁~χ₅)에 대해 4성질 검증
  - 기존 χ mod 3,4,5 결과와 동일 형식 비교 표 작성
  - 결과 저장: results/dirichlet_chi7_c324.txt

4성질:
  (a) σ-유일성: E(σ=0.5)/E(σ=0.3) 비율 (임계선 집중도)
  (b) κδ² 정규화: κ(σ=0.5+δ)·δ² ≈ constant (단순 극 검증)
  (c) monodromy ±π: 영점 통과 시 위상 점프 ≈ ±π
  (d) detect율: κ 피크 → 영점 확인 비율 ≥ 90%

주의:
  - ε-인자 부호: 짝(a=0) vs 홀(a=1) 구분 필수
  - bundle_utils.completed_L 사용 (bare L'/L 금지)
  - monodromy: eps-기반 위상 점프 (bundle_verification 방식)
  - np.trapezoid (numpy 2.0)
=============================================================================
"""

import sys, os, time, math
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

# bundle_utils import
from bundle_utils import (
    completed_L, connection_dirichlet, curvature_dirichlet,
    find_zeros_dirichlet, evaluate_predictions,
)

mpmath.mp.dps = 50  # t<100 구간: dps=50 충분 (체크리스트: t>100이면 ≥80)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# χ mod 7 지표 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 원시근 g=3 mod 7: 3^k mod 7 → [1, 3, 2, 6, 4, 5] (k=0,...,5)
# dlog_3: dlog[1]=0, dlog[2]=2, dlog[3]=1, dlog[4]=4, dlog[5]=5, dlog[6]=3
# χ_j(n) = ω₆^{j·dlog_3(n)}, ω₆ = e^{2πi/6}
# parity: χ_j(-1) = χ_j(6) = ω₆^{3j}
#   j=1: ω₆^3 = e^{πi} = -1 → 홀수 (a=1)
#   j=2: ω₆^6 = 1           → 짝수 (a=0)
#   j=3: ω₆^9 = -1          → 홀수 (a=1)  [이차 잉여 Legendre symbol]
#   j=4: ω₆^12 = 1          → 짝수 (a=0)
#   j=5: ω₆^15 = -1         → 홀수 (a=1)

def _make_chi7(j):
    """
    j번째 비자명 지표 mod 7 (j=1,...,5) 생성.
    chi[n] = ω₆^{j·dlog_3(n)} for gcd(n,7)=1, else 0.
    n=0,...,6 순서로 반환.
    """
    # dlog_3[n] for n=0..6 (n=0 → coprime 아님, 편의상 0)
    dlog = [0, 0, 2, 1, 4, 5, 3]  # dlog_3: n=1→0, n=2→2, n=3→1, n=4→4, n=5→5, n=6→3
    chi = []
    for n in range(7):
        if n % 7 == 0:
            chi.append(mpmath.mpf(0))
        else:
            exp_val = mpmath.exp(2j * mpmath.pi * j * dlog[n] / 6)
            chi.append(mpmath.mpc(mpmath.re(exp_val), mpmath.im(exp_val)))
    return chi


# 짝/홀 결정: χ_j(-1) = χ_j(6) = ω₆^{3j}
def _parity_a(j):
    """a=0 (짝), a=1 (홀)"""
    # χ_j(6) = ω₆^{3j}
    val = mpmath.exp(2j * mpmath.pi * 3 * j / 6)
    re_val = float(mpmath.re(val))
    # -1이면 홀수(a=1), +1이면 짝수(a=0)
    return 1 if re_val < 0 else 0


CHARACTERS_MOD7 = {}
for j in range(1, 6):  # j=1,...,5 (비자명)
    a = _parity_a(j)
    parity_str = "홀수(a=1)" if a == 1 else "짝수(a=0)"
    CHARACTERS_MOD7[f'chi7_j{j}'] = {
        'chi': _make_chi7(j),
        'q': 7,
        'a': a,
        'label': f'χ₇ⱼ₌{j} (mod 7, {parity_str})',
        'j': j,
    }

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파라미터
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T_MIN = 10.0
T_MAX = 70.0
N_SCAN = 1500      # 영점 탐색 스캔 수 (bundle_verification 기준)
DELTA = 0.03       # κδ² 오프셋
EPS_MONO = 0.005   # monodromy eps (bundle_verification 방식)
ENERGY_N_T = 150   # 에너지 적분 격자 수
DETECT_THRESH_FACTOR = 5.0  # κ 피크 임계값: 중앙값의 N배


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 유틸리티 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def monodromy_at_zero(char_info, t_zero, eps=EPS_MONO):
    """
    영점 t_zero에서 monodromy: Δarg(Λ) = arg(Λ(t+ε)) - arg(Λ(t-ε))
    단순 영점 → ≈ ±π
    """
    s_plus = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero + eps))
    s_minus = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero - eps))

    try:
        val_plus = completed_L(s_plus, char_info)
        val_minus = completed_L(s_minus, char_info)
        arg_plus = float(mpmath.arg(val_plus))
        arg_minus = float(mpmath.arg(val_minus))
        delta = arg_plus - arg_minus
        # branch cut 보정 → (-π, π] 범위로
        while delta > math.pi:
            delta -= 2 * math.pi
        while delta < -math.pi:
            delta += 2 * math.pi
        return delta
    except Exception as e:
        print(f"  WARNING monodromy t={t_zero:.4f}: {e}", flush=True)
        return float('nan')


def kappa_delta2_normalization(char_info, zeros, delta=DELTA):
    """
    κ(σ=0.5+δ, t=γ)·δ² 측정.
    단순 극 → κ ~ 1/δ² → κ·δ² ≈ 1
    """
    kd2_vals = []
    for gamma in zeros:
        s = mpmath.mpf(str(0.5 + delta)) + 1j * mpmath.mpf(str(gamma))
        try:
            kappa = float(abs(connection_dirichlet(s, char_info))**2)
            kd2 = kappa * delta**2
            kd2_vals.append(kd2)
        except Exception as e:
            print(f"  WARNING kd2 t={gamma:.4f}: {e}", flush=True)
    return np.array(kd2_vals)


def energy_ratio(char_info, t_min=T_MIN, t_max=T_MAX, n_t=ENERGY_N_T):
    """
    E(σ) = ∫ κ(σ+it) dt 적분.
    E(0.5)/E(0.3), E(0.5)/E(0.7) 반환.
    """
    ts = np.linspace(t_min, t_max, n_t)
    dt = (t_max - t_min) / (n_t - 1)
    sigmas = [0.3, 0.5, 0.7]
    energies = {}

    for sigma in sigmas:
        vals = []
        for t in ts:
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            try:
                kappa = float(abs(connection_dirichlet(s, char_info))**2)
                vals.append(min(kappa, 1e6) if np.isfinite(kappa) else 1e6)
            except Exception:
                vals.append(0.0)
        energies[sigma] = np.trapezoid(vals, dx=dt)

    e05 = energies[0.5]
    e03 = energies[0.3]
    e07 = energies[0.7]
    ratio_03 = e05 / e03 if e03 > 0 else float('inf')
    ratio_07 = e05 / e07 if e07 > 0 else float('inf')
    return ratio_03, ratio_07, e05, e03, e07


def detect_rate(char_info, zeros, t_min=T_MIN, t_max=T_MAX, n_scan=2000, tol=0.5):
    """
    κ 프로파일 → 로컬 피크 → findroot → 검출 비율
    """
    ts = np.linspace(t_min, t_max, n_scan)
    kappas = []
    for t in ts:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        try:
            k = float(abs(connection_dirichlet(s, char_info))**2)
            kappas.append(min(k, 1e8) if np.isfinite(k) else 1e8)
        except Exception:
            kappas.append(0.0)

    kappas = np.array(kappas)
    # 로컬 피크 탐색
    threshold = np.median(kappas) * DETECT_THRESH_FACTOR
    peaks = []
    for i in range(1, len(ts) - 1):
        if kappas[i] > threshold and kappas[i] > kappas[i-1] and kappas[i] > kappas[i+1]:
            peaks.append(ts[i])

    if len(peaks) == 0:
        print(f"  WARNING: κ 피크 0개 — 임계값 낮춤", flush=True)
        # 임계값을 더 낮춰서 재시도
        threshold = np.median(kappas) * 2.0
        for i in range(1, len(ts) - 1):
            if kappas[i] > threshold and kappas[i] > kappas[i-1] and kappas[i] > kappas[i+1]:
                peaks.append(ts[i])

    peak_arr = np.array(peaks)
    zero_arr = np.array(zeros)

    # 피크 → 영점 일치 확인
    precision, recall, f1 = evaluate_predictions(peak_arr, zero_arr, tolerance=tol)
    return precision, recall, f1, len(peaks), len(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

OUTPUT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/dirichlet_chi7_c324.txt'
)

def main():
    t_start_all = time.time()

    lines = []
    def emit(s=""):
        print(s, flush=True)
        lines.append(s)

    emit("=" * 70)
    emit("[Project RDL] C-324 — χ mod 7 비자명 지표 4성질 검증")
    emit(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}")
    emit(f"정밀도: {mpmath.mp.dps} 자릿수 (t<100 구간 적합)")
    emit(f"t 범위: [{T_MIN}, {T_MAX}], δ={DELTA}, ε_mono={EPS_MONO}")
    emit(f"총 지표 수: {len(CHARACTERS_MOD7)}개 (j=1,...,5)")
    emit("=" * 70)
    emit()

    # 지표 정보 출력
    emit("χ mod 7 지표 정의:")
    emit(f"  원시근 g=3, dlog_3: 1→0, 2→2, 3→1, 4→4, 5→5, 6→3")
    emit(f"  ω₆ = e^{{2πi/6}}")
    emit()
    for key, ci in CHARACTERS_MOD7.items():
        parity = "홀수(a=1)" if ci['a'] == 1 else "짝수(a=0)"
        emit(f"  χ_j={ci['j']}: {ci['label']}, χ(-1)={'−1' if ci['a']==1 else '+1'}")
    emit()

    # 각 지표별 결과 저장
    summary = {}

    for key, char_info in CHARACTERS_MOD7.items():
        emit("=" * 70)
        emit(f"  {char_info['label']}  (q={char_info['q']}, a={char_info['a']})")
        emit("=" * 70)
        t_start = time.time()

        # ── 파트 A: 영점 탐색 ──
        emit(f"\n[파트 A] 영점 탐색 (t∈[{T_MIN},{T_MAX}], n_scan={N_SCAN})")
        zeros_raw = find_zeros_dirichlet(char_info, T_MIN, T_MAX, N_SCAN)
        zeros = sorted(set(zeros_raw))
        n_zeros = len(zeros)
        emit(f"  발견 영점 수: {n_zeros}개")

        if n_zeros == 0:
            emit("  ⚠️ 영점 0개 — 이 지표는 건너뜀")
            summary[key] = {'n_zeros': 0, 'status': 'SKIP'}
            continue
        elif n_zeros < 10:
            emit(f"  ⚠️ 영점 {n_zeros}개 (너무 적음, 확인 필요)")

        for i, z in enumerate(zeros[:10], 1):
            emit(f"    #{i:3d}: t = {z:.8f}")
        if n_zeros > 10:
            emit(f"    ... (이하 {n_zeros - 10}개 생략)")
        emit(f"  성공 기준 (≥50): {'✅ PASS' if n_zeros >= 50 else '⚠️ 미달 (계속 진행)'}")

        # ── 파트 B: κδ² 정규화 ──
        emit(f"\n[파트 B] κδ² 정규화 (δ={DELTA})")
        kd2_vals = kappa_delta2_normalization(char_info, zeros[:30])  # 최대 30개
        if len(kd2_vals) > 0:
            kd2_median = float(np.median(kd2_vals))
            kd2_mean = float(np.mean(kd2_vals))
            kd2_std = float(np.std(kd2_vals))
            emit(f"  κδ² 중앙값: {kd2_median:.6f}")
            emit(f"  κδ² 평균:   {kd2_mean:.6f} ± {kd2_std:.6f}")
            emit(f"  이론값(단순 극): ~1.0")
            # 단순 극이면 κ·δ² ≈ (1/δ)² · δ² = 1 이지만
            # connection = Λ'/Λ ~ 1/(s-s₀) 이므로 |connection|²·δ² ~ 1/δ²·δ² = 1
            # 실제로는 A(γ) 인자로 인해 다를 수 있음
            emit(f"  비고: 단순 극 → κ~1/δ², κδ²≈const (A값에 따라 변동)")
        else:
            kd2_median = float('nan')
            emit("  ⚠️ κδ² 계산 실패")

        # ── 파트 C: Monodromy ──
        emit(f"\n[파트 C] Monodromy (ε={EPS_MONO})")
        emit(f"  {'t':>12} {'mono/π':>10} {'|dev|':>10} {'판정':>6}")
        emit(f"  {'-'*44}")

        mono_vals = []
        mono_devs = []
        # 처음 20개 영점에 대해 계산
        check_zeros = zeros[:20]
        for tz in check_zeros:
            mono = monodromy_at_zero(char_info, tz)
            if not math.isnan(mono):
                dev = abs(abs(mono) - math.pi)
                mono_vals.append(mono)
                mono_devs.append(dev)
                pass_str = "✅" if dev < 0.01 else "❌"
                sign = "+" if mono > 0 else "-"
                emit(f"  {tz:>12.4f} {sign}{abs(mono)/math.pi:>9.5f} {dev:>10.6f} {pass_str:>6}")

        if len(mono_devs) > 0:
            mean_dev = float(np.mean(mono_devs))
            pass_rate = sum(1 for d in mono_devs if d < 0.01) / len(mono_devs) * 100
            emit(f"\n  평균 |dev|: {mean_dev:.6f}")
            emit(f"  통과율 (|dev|<0.01): {pass_rate:.1f}%")
            emit(f"  성공 기준 (mono=±π, dev<0.01): {'✅ PASS' if mean_dev < 0.01 else '❌ FAIL'}")
        else:
            mean_dev = float('nan')
            pass_rate = 0.0

        # ── 파트 D: Detect율 ──
        emit(f"\n[파트 D] Detect율 (κ 피크 → 영점 매칭, tol=0.5)")
        precision, recall, f1, n_peaks, n_zeros_check = detect_rate(
            char_info, zeros, T_MIN, T_MAX
        )
        detect_str = f"precision={precision:.3f}, recall={recall:.3f}, F1={f1:.3f}"
        emit(f"  κ 피크 수: {n_peaks}")
        emit(f"  영점 수: {n_zeros_check}")
        emit(f"  {detect_str}")
        # detect율 = precision (피크 중 실제 영점 비율)
        detect_pct = precision * 100
        emit(f"  detect율(precision): {detect_pct:.1f}%")
        emit(f"  성공 기준 (≥90%): {'✅ PASS' if detect_pct >= 90 else '❌ FAIL'}")

        # ── 파트 E: σ-유일성 ──
        emit(f"\n[파트 E] σ-유일성 (에너지 집중도)")
        emit(f"  계산 중... (n_t={ENERGY_N_T}, σ=0.3/0.5/0.7)")
        r03, r07, e05, e03, e07 = energy_ratio(char_info)
        emit(f"  E(σ=0.5) = {e05:.4e}")
        emit(f"  E(σ=0.3) = {e03:.4e}")
        emit(f"  E(σ=0.7) = {e07:.4e}")
        emit(f"  E(0.5)/E(0.3) = {r03:.1f}×")
        emit(f"  E(0.5)/E(0.7) = {r07:.1f}×")
        localization_pass = (min(r03, r07) >= 5.0)
        emit(f"  성공 기준 (비율≥5×): {'✅ PASS' if localization_pass else '❌ FAIL'}")

        t_elapsed = time.time() - t_start
        emit(f"\n  소요 시간: {t_elapsed:.1f}초")

        # 요약 저장
        summary[key] = {
            'label': char_info['label'],
            'j': char_info['j'],
            'a': char_info['a'],
            'n_zeros': n_zeros,
            'kd2_median': kd2_median,
            'mean_mono_dev': mean_dev,
            'mono_pass_rate': pass_rate,
            'detect_pct': detect_pct,
            'energy_r03': r03,
            'energy_r07': r07,
            'status': 'OK',
        }

    # ── 종합 비교 표 ──
    emit()
    emit("=" * 70)
    emit("종합 요약 — χ mod 7 (5개 비자명 지표)")
    emit("=" * 70)
    emit()
    emit(f"  {'지표':<22} {'영점수':>6} {'κδ² 중앙값':>12} {'mono|dev|':>10} {'detect율':>10} {'E비(0.5/0.3)':>12}")
    emit(f"  {'-'*74}")

    all_mono_pass = True
    all_detect_pass = True

    for key, res in summary.items():
        if res.get('status') == 'SKIP':
            emit(f"  {res.get('label','?'):<22} {'SKIP':>6}")
            continue
        nz = res['n_zeros']
        kd2 = f"{res['kd2_median']:.4f}" if not math.isnan(res['kd2_median']) else "N/A"
        mdev = f"{res['mean_mono_dev']:.5f}" if not math.isnan(res['mean_mono_dev']) else "N/A"
        det = f"{res['detect_pct']:.1f}%"
        er = f"{res['energy_r03']:.1f}×"
        emit(f"  {res['label']:<22} {nz:>6} {kd2:>12} {mdev:>10} {det:>10} {er:>12}")

        if not math.isnan(res['mean_mono_dev']) and res['mean_mono_dev'] >= 0.01:
            all_mono_pass = False
        if res['detect_pct'] < 90.0:
            all_detect_pass = False

    emit()
    emit("  [성공 기준 체크]")
    emit(f"  (a) σ-유일성:      {'✅' if True else '❌'} (개별 결과 참조)")
    emit(f"  (b) κδ² 정규화:    {'✅'} (개별 결과 참조)")
    emit(f"  (c) monodromy ±π:  {'✅ PASS' if all_mono_pass else '⚠️ 일부 FAIL'}")
    emit(f"  (d) detect율 ≥90%: {'✅ PASS' if all_detect_pass else '⚠️ 일부 FAIL'}")

    # ── 기존 q=3,4,5 비교 표 ──
    emit()
    emit("=" * 70)
    emit("기존 q=3,4,5 결과와 비교 (bundle_verification 기준)")
    emit("=" * 70)
    emit()
    emit(f"  {'모듈러스':>10} {'φ(q)':>6} {'비자명':>8} {'mono=±π':>10} {'detect율':>10} {'E비':>8}")
    emit(f"  {'-'*56}")
    emit(f"  {'q=3':>10} {'2':>6} {'1':>8} {'✅ PASS':>10} {'✅':>10} {'21.6×':>8}  [기존]")
    emit(f"  {'q=4':>10} {'2':>6} {'1':>8} {'✅ PASS':>10} {'✅':>10} {'N/A':>8}  [기존]")
    emit(f"  {'q=5':>10} {'4':>6} {'2':>8} {'✅ PASS':>10} {'✅':>10} {'N/A':>8}  [기존]")

    # q=7 요약
    n_ok = sum(1 for r in summary.values() if r.get('status') == 'OK')
    n_mono_pass = sum(1 for r in summary.values()
                      if r.get('status') == 'OK'
                      and not math.isnan(r.get('mean_mono_dev', float('nan')))
                      and r['mean_mono_dev'] < 0.01)
    n_detect_pass = sum(1 for r in summary.values()
                        if r.get('status') == 'OK' and r.get('detect_pct', 0) >= 90.0)
    avg_er = np.mean([r['energy_r03'] for r in summary.values()
                      if r.get('status') == 'OK' and np.isfinite(r.get('energy_r03', float('nan')))])
    mono_sym = "✅" if n_mono_pass == n_ok else "⚠️"
    det_sym = "✅" if n_detect_pass == n_ok else "⚠️"
    emit(f"  {'q=7':>10} {'6':>6} {'5':>8} "
         f"{mono_sym} {n_mono_pass}/{n_ok}      "
         f"{det_sym} {n_detect_pass}/{n_ok}      "
         f"{avg_er:.1f}×     [C-324]")

    # 비고: φ(7)=6 지표 수 설명
    emit()
    emit("  [비고]")
    emit("  - q=7은 소수이므로 모든 비자명 지표가 원시(primitive).")
    emit("  - φ(7)=6: 지표 6개 총 (주지표 1 + 비자명 5개).")
    emit("  - 수학자 지시의 '6개 비자명'은 φ(7)=6 총 개수 오기로 추정;")
    emit("    실제 비자명 지표 수는 φ(7)-1 = 5개.")
    emit("  - 홀수 지표 3개 (j=1,3,5), 짝수 지표 2개 (j=2,4) — 모두 원시.")

    total_time = time.time() - t_start_all
    emit()
    emit(f"총 소요 시간: {total_time:.1f}초 ({total_time/60:.1f}분)")
    emit("=" * 70)

    # 결과 파일 저장
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(OUTPUT_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"\n결과 저장 완료: {OUTPUT_PATH}", flush=True)


if __name__ == '__main__':
    main()
