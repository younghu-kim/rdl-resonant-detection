"""
=============================================================================
[Project RDL] 디리클레 L-함수 다발 검증: ξ-다발 프레임워크의 일반화
=============================================================================
리만 제타 함수에서 검증된 ξ-다발 구조가 디리클레 L-함수로 일반화되는지 확인.

대상 지표(character):
  - χ mod 3: 가장 단순한 비자명 지표
  - χ mod 4: 실수 지표 (크로네커 기호)
  - χ mod 5: 복소 지표

검증 항목:
  1. σ별 위상 점프 수 → σ=1/2에서만 점프 발생
  2. 곡률 집중도: 영점 근방 vs 일반 점
  3. 모노드로미 = ±π 양자화
  4. 에너지 프로파일 E(σ) → σ=1/2에서 피크

수학적 객체:
  - 완비 L-함수: Λ(s, χ) = (q/π)^{s/2} Γ((s+a)/2) L(s, χ)
    여기서 a=0 (χ(-1)=1), a=1 (χ(-1)=-1), q=도체(conductor)
  - 접속: L(s) = Λ'/Λ (로그 미분)
  - 곡률: κ = |L|²
  - 모노드로미: 영점 주변 Δarg(Λ)
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))

mpmath.mp.dps = 80  # 큰 t에서 정밀도 확보


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 디리클레 지표 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 지표값 리스트: chi[k] = χ(k), k=0,1,...,q-1
# mpmath.dirichlet(s, chi) 형식에 맞춤

CHARACTERS = {
    'χ mod 3': {
        'chi': [0, 1, -1],           # χ(0)=0, χ(1)=1, χ(2)=-1
        'q': 3,
        'a': 1,                       # χ(-1) = χ(2) = -1, 즉 a=1
        'label': 'χ₃ (mod 3, 비자명)',
    },
    'χ mod 4': {
        'chi': [0, 1, 0, -1],        # χ(0)=0, χ(1)=1, χ(2)=0, χ(3)=-1
        'q': 4,
        'a': 1,                       # χ(-1) = χ(3) = -1, 즉 a=1
        'label': 'χ₄ (mod 4, 실수 지표)',
    },
    'χ mod 5': {
        # 원시 지표 mod 5: χ(1)=1, χ(2)=i, χ(3)=-i, χ(4)=-1
        'chi': [0, 1, 1j, -1j, -1],
        'q': 5,
        'a': 0,                       # χ(-1) = χ(4) = -1... 재확인
        'label': 'χ₅ (mod 5, 복소 지표)',
    },
}

# χ mod 5 패리티 보정: χ(-1) = χ(4) = -1이므로 a=1
CHARACTERS['χ mod 5']['a'] = 1


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def completed_L(s, char_info):
    """
    완비 L-함수: Λ(s, χ) = (q/π)^{s/2} Γ((s+a)/2) L(s, χ)
    """
    q = mpmath.mpf(char_info['q'])
    a = mpmath.mpf(char_info['a'])
    chi = char_info['chi']

    L_val = mpmath.dirichlet(s, chi)
    gamma_val = mpmath.gamma((s + a) / 2)
    prefactor = mpmath.power(q / mpmath.pi, s / 2)

    return prefactor * gamma_val * L_val


def connection_L(s, char_info):
    """접속 L(s) = Λ'/Λ (수치 미분)"""
    h = mpmath.mpf(1) / mpmath.mpf(10**20)
    Lambda_val = completed_L(s, char_info)
    Lambda_abs = abs(Lambda_val)

    # 영점 근방 → 매우 높은 곡률
    if Lambda_abs < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e10, 0)

    Lambda_d = (completed_L(s + h, char_info) - completed_L(s - h, char_info)) / (2 * h)
    return Lambda_d / Lambda_val


def curvature_at(s, char_info):
    """곡률 κ = |L|² = |Λ'/Λ|²"""
    L = connection_L(s, char_info)
    return float(abs(L)**2)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 영점 탐색: 임계선 위 부호 변화로 영점 위치 찾기
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def find_zeros_dirichlet(char_info, t_min=10.0, t_max=40.0, n_scan=2000):
    """
    임계선 σ=1/2 위에서 Re(Λ)의 부호 변화를 탐색하여 영점 위치를 찾는다.
    mpmath.findroot로 정밀화.
    """
    chi = char_info['chi']
    ts = np.linspace(t_min, t_max, n_scan)
    zeros = []

    prev_re = None
    prev_t = None

    for t in ts:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        val = completed_L(s, char_info)
        curr_re = mpmath.re(val)

        if prev_re is not None and prev_re * curr_re < 0:
            # 부호 변화 → 영점 근사
            try:
                def f_real(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    return mpmath.re(completed_L(sv, char_info))

                t_zero = mpmath.findroot(f_real, (prev_t, t))
                t_zero_f = float(t_zero)
                # 중복 방지 (0.1 이내)
                if not zeros or abs(t_zero_f - zeros[-1]) > 0.1:
                    zeros.append(t_zero_f)
            except Exception:
                pass

        prev_re = curr_re
        prev_t = t

    # Im(Λ) 부호 변화도 확인 (일부 영점은 Re만으로 못 잡음)
    prev_im = None
    prev_t = None

    for t in ts:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        val = completed_L(s, char_info)
        curr_im = mpmath.im(val)

        if prev_im is not None and prev_im * curr_im < 0:
            try:
                def f_imag(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    return mpmath.im(completed_L(sv, char_info))

                t_zero = mpmath.findroot(f_imag, (prev_t, t))
                t_zero_f = float(t_zero)
                # |Λ| 자체가 작은지 확인 (실제 영점인지)
                sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero_f))
                if abs(completed_L(sv, char_info)) < mpmath.mpf('1e-10'):
                    if not any(abs(t_zero_f - z) < 0.1 for z in zeros):
                        zeros.append(t_zero_f)
            except Exception:
                pass

        prev_im = curr_im
        prev_t = t

    zeros.sort()
    return np.array(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실험 1: σ별 위상 점프 수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def count_phase_jumps(char_info, sigma, t_min=10.0, t_max=40.0, n_points=1500):
    """σ 고정, t를 스캔하며 |Δarg(Λ)| > π/2인 점프 수 세기"""
    ts = np.linspace(t_min, t_max, n_points)
    jumps = 0
    prev_arg = None

    for t in ts:
        s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
        val = completed_L(s, char_info)

        if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
            prev_arg = None
            continue

        curr_arg = float(mpmath.arg(val))

        if prev_arg is not None:
            delta = curr_arg - prev_arg
            while delta > np.pi:
                delta -= 2 * np.pi
            while delta < -np.pi:
                delta += 2 * np.pi

            if abs(delta) > np.pi / 2:
                jumps += 1

        prev_arg = curr_arg

    return jumps


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실험 2: 곡률 집중도
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def curvature_concentration(char_info, zeros, t_min=10.0, t_max=40.0, n_generic=100):
    """영점 근방 곡률 vs 일반 점 곡률 비교"""
    # 영점 근방 곡률 (각 영점에서 ±0.01 떨어진 점)
    zero_kappas = []
    for tz in zeros:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(tz + 0.01))
        try:
            k = curvature_at(s, char_info)
            if np.isfinite(k):
                zero_kappas.append(k)
        except Exception:
            pass

    # 일반 점 곡률 (영점에서 멀리 떨어진 점)
    generic_kappas = []
    rng = np.random.RandomState(42)
    attempts = 0
    while len(generic_kappas) < n_generic and attempts < n_generic * 3:
        t = rng.uniform(t_min, t_max)
        # 영점에서 최소 1.0 이상 떨어진 점만
        if len(zeros) > 0 and np.min(np.abs(zeros - t)) < 1.0:
            attempts += 1
            continue
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        try:
            k = curvature_at(s, char_info)
            if np.isfinite(k):
                generic_kappas.append(k)
        except Exception:
            pass
        attempts += 1

    zero_med = np.median(zero_kappas) if zero_kappas else 0
    gen_med = np.median(generic_kappas) if generic_kappas else 1
    ratio = zero_med / gen_med if gen_med > 0 else float('inf')

    return zero_med, gen_med, ratio


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실험 3: 모노드로미 양자화 검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def monodromy_at_zero(char_info, t_zero, eps=0.005):
    """영점 t_zero에서 모노드로미: Δarg(Λ) = arg(Λ(t+ε)) - arg(Λ(t-ε))"""
    s_plus = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero + eps))
    s_minus = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero - eps))

    val_plus = completed_L(s_plus, char_info)
    val_minus = completed_L(s_minus, char_info)

    arg_plus = float(mpmath.arg(val_plus))
    arg_minus = float(mpmath.arg(val_minus))

    delta = arg_plus - arg_minus
    while delta > np.pi:
        delta -= 2 * np.pi
    while delta < -np.pi:
        delta += 2 * np.pi

    return delta


def monodromy_quantization(char_info, zeros):
    """각 영점에서 모노드로미가 ±π에 얼마나 가까운지 측정"""
    results = []
    for tz in zeros:
        mono = monodromy_at_zero(char_info, tz)
        deviation = abs(abs(mono) - np.pi)  # |mono| - π 로부터의 거리
        results.append((tz, mono, deviation))
    return results


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실험 4: 에너지 σ-프로파일
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def energy_profile(char_info, sigmas, t_min=12.0, t_max=38.0, n_t=150):
    """E(σ) = ∫|L(σ+it)|² dt (사다리꼴 적분), L = Λ'/Λ"""
    ts = np.linspace(t_min, t_max, n_t)
    dt = (t_max - t_min) / (n_t - 1)
    energies = []

    for sigma in sigmas:
        vals = []
        for t in ts:
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            try:
                k = curvature_at(s, char_info)
                if np.isfinite(k):
                    vals.append(k)
                else:
                    vals.append(1e6)
            except Exception:
                vals.append(0.0)
        E = np.trapezoid(vals, dx=dt)
        energies.append(E)

    return np.array(energies)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def run_verification(char_name, char_info, out_file):
    """단일 지표에 대한 전체 검증 실행"""
    label = char_info['label']
    q = char_info['q']
    a = char_info['a']

    print(f"\n{'='*70}")
    print(f"  {label}  (q={q}, a={a})")
    print(f"{'='*70}")

    # ── 영점 탐색 ──
    print(f"\n[0] 영점 탐색: t∈[10, 40]")
    t0 = time.time()
    zeros = find_zeros_dirichlet(char_info, 10.0, 40.0, n_scan=2000)
    elapsed = time.time() - t0
    print(f"  발견된 영점: {len(zeros)}개 ({elapsed:.1f}초)")
    for i, tz in enumerate(zeros):
        print(f"    #{i+1}: t = {tz:.6f}")

    out_file.write(f"\n{'='*70}\n")
    out_file.write(f"  {label}  (q={q}, a={a})\n")
    out_file.write(f"{'='*70}\n")
    out_file.write(f"\n영점 ({len(zeros)}개, t∈[10,40]):\n")
    for i, tz in enumerate(zeros):
        out_file.write(f"  #{i+1}: t = {tz:.6f}\n")

    # ── 실험 1: σ별 위상 점프 ──
    print(f"\n[1] σ별 위상 점프 수 (t∈[10,40], 1500점)")
    sigmas_jump = [0.1, 0.2, 0.3, 0.4, 0.49, 0.5, 0.51, 0.6, 0.7, 0.8, 0.9]
    jump_results = {}
    for sigma in sigmas_jump:
        n_jumps = count_phase_jumps(char_info, sigma)
        jump_results[sigma] = n_jumps
        marker = "★" if sigma == 0.5 else " "
        print(f"  {marker} σ={sigma:.2f}: {n_jumps} jumps")

    out_file.write(f"\n[1] 위상 점프 수 (t∈[10,40], 1500점)\n")
    out_file.write(f"{'σ':>8} {'jumps':>8}\n")
    out_file.write("-" * 20 + "\n")
    for sigma in sigmas_jump:
        marker = "★" if sigma == 0.5 else " "
        out_file.write(f"{marker}{sigma:>7.2f} {jump_results[sigma]:>8}\n")

    # ── 실험 2: 곡률 집중도 ──
    print(f"\n[2] 곡률 집중도: 영점 근방 vs 일반 점")
    if len(zeros) > 0:
        zero_med, gen_med, ratio = curvature_concentration(char_info, zeros, n_generic=50)
        print(f"  영점 근방 중앙값 κ: {zero_med:.2e}")
        print(f"  일반 점 중앙값 κ:   {gen_med:.2e}")
        print(f"  집중도 비율:        {ratio:.1f}×")

        out_file.write(f"\n[2] 곡률 집중도\n")
        out_file.write(f"  영점 근방 median κ: {zero_med:.2e}\n")
        out_file.write(f"  일반 점 median κ:   {gen_med:.2e}\n")
        out_file.write(f"  비율:               {ratio:.1f}×\n")
    else:
        print("  (영점 미발견, 건너뜀)")
        out_file.write(f"\n[2] 곡률 집중도: 영점 미발견\n")

    # ── 실험 3: 모노드로미 양자화 ──
    print(f"\n[3] 모노드로미 양자화 검증")
    if len(zeros) > 0:
        mono_results = monodromy_quantization(char_info, zeros)
        deviations = [r[2] for r in mono_results]
        mean_dev = np.mean(deviations)

        out_file.write(f"\n[3] 모노드로미 양자화\n")
        out_file.write(f"{'t':>12} {'mono/π':>12} {'|dev|':>12}\n")
        out_file.write("-" * 40 + "\n")

        for tz, mono, dev in mono_results:
            sign = "+" if mono > 0 else "-"
            print(f"    t={tz:.4f}: Δarg = {sign}{abs(mono)/np.pi:.4f}π  "
                  f"(|π로부터 편차| = {dev:.4f})")
            out_file.write(f"{tz:>12.4f} {mono/np.pi:>+12.4f} {dev:>12.4f}\n")

        print(f"  평균 |π로부터 편차|: {mean_dev:.4f}")
        out_file.write(f"\n  평균 |π로부터 편차|: {mean_dev:.4f}\n")
    else:
        print("  (영점 미발견, 건너뜀)")
        out_file.write(f"\n[3] 모노드로미: 영점 미발견\n")

    # ── 실험 4: 에너지 σ-프로파일 ──
    print(f"\n[4] 에너지 σ-프로파일 (t∈[12,38], 150점)")
    sigmas_energy = np.linspace(0.1, 0.9, 17)
    energies = energy_profile(char_info, sigmas_energy, n_t=150)

    idx_peak = np.argmax(energies)
    # σ=0.5과 σ=0.3에 가장 가까운 인덱스
    idx_05 = np.argmin(np.abs(sigmas_energy - 0.5))
    idx_03 = np.argmin(np.abs(sigmas_energy - 0.3))
    e_ratio = energies[idx_05] / energies[idx_03] if energies[idx_03] > 0 else float('inf')

    print(f"  에너지 피크: σ={sigmas_energy[idx_peak]:.3f}, E={energies[idx_peak]:.2e}")
    print(f"  E(0.5)/E(0.3) = {e_ratio:.1f}×")

    out_file.write(f"\n[4] 에너지 σ-프로파일 (t∈[12,38])\n")
    out_file.write(f"{'σ':>8} {'E':>14}\n")
    out_file.write("-" * 25 + "\n")
    for sigma, E in zip(sigmas_energy, energies):
        marker = "★" if abs(sigma - 0.5) < 0.02 else " "
        out_file.write(f"{marker}{sigma:>7.3f} {E:>14.4e}\n")
    out_file.write(f"\n  에너지 피크: σ={sigmas_energy[idx_peak]:.3f}\n")
    out_file.write(f"  E(0.5)/E(0.3) = {e_ratio:.1f}×\n")

    # 결과 사전 반환 (cross_comparison용)
    return {
        'n_zeros': len(zeros),
        'zeros': zeros,
        'jump_results': jump_results,
        'curvature_ratio': ratio if len(zeros) > 0 else None,
        'mono_mean_dev': mean_dev if len(zeros) > 0 else None,
        'energy_ratio': e_ratio,
        'energy_peak_sigma': sigmas_energy[idx_peak],
    }


def main():
    out_path = os.path.expanduser(
        '~/Desktop/gdl_unified/results/dirichlet_bundle_verification.txt'
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    print("=" * 70)
    print("디리클레 L-함수 다발 검증: ξ-다발 프레임워크 일반화")
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}")
    print("=" * 70)

    with open(out_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("디리클레 L-함수 다발 검증: ξ-다발 프레임워크 일반화\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"정밀도: {mpmath.mp.dps} 자릿수\n")
        f.write("=" * 70 + "\n")

        all_results = {}
        for char_name, char_info in CHARACTERS.items():
            t0 = time.time()
            result = run_verification(char_name, char_info, f)
            elapsed = time.time() - t0
            all_results[char_name] = result
            f.write(f"\n  (소요 시간: {elapsed:.0f}초)\n")
            print(f"\n  [소요 시간: {elapsed:.0f}초]")

        # ── 종합 요약 ──
        f.write(f"\n\n{'='*70}\n")
        f.write("종합 요약\n")
        f.write(f"{'='*70}\n\n")

        f.write(f"{'지표':<20} {'영점수':>8} {'κ비율':>10} {'mono편차':>10} "
                f"{'E비율':>10} {'피크σ':>8}\n")
        f.write("-" * 70 + "\n")

        print(f"\n\n{'='*70}")
        print("종합 요약")
        print(f"{'='*70}")
        print(f"{'지표':<20} {'영점수':>8} {'κ비율':>10} {'mono편차':>10} "
              f"{'E비율':>10} {'피크σ':>8}")
        print("-" * 70)

        for char_name, r in all_results.items():
            kappa_str = f"{r['curvature_ratio']:.1f}×" if r['curvature_ratio'] else "N/A"
            mono_str = f"{r['mono_mean_dev']:.4f}" if r['mono_mean_dev'] is not None else "N/A"
            line = (f"{char_name:<20} {r['n_zeros']:>8} {kappa_str:>10} "
                    f"{mono_str:>10} {r['energy_ratio']:>9.1f}× "
                    f"{r['energy_peak_sigma']:>7.3f}")
            f.write(line + "\n")
            print(line)

    print(f"\n결과 저장: {out_path}")
    print("완료.")


if __name__ == '__main__':
    main()
