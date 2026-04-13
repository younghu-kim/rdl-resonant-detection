"""
=============================================================================
[Project RDL] 디리클레 L-함수 교차 비교: ζ vs χ₃ vs χ₄ vs χ₅
=============================================================================
리만 제타 함수와 3개 디리클레 L-함수의 ξ-다발 구조를 동일 기준으로 비교.

비교 항목 (t∈[10, 40]):
  1. 영점 수
  2. 곡률 집중도: 영점 근방 중앙값 / 일반 점 중앙값
  3. 모노드로미 양자화 정확도: |Δarg|가 ±π에 얼마나 가까운지
  4. 에너지 집중도: E(σ=0.5) / E(σ=0.3)
  5. 위상 점프 유일성: σ=0.5 점프 수 vs 타 σ 점프 수

출력: 깔끔한 비교 표
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))

mpmath.mp.dps = 80

T_MIN, T_MAX = 10.0, 40.0


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 함수 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# ── 디리클레 L-함수용 ──

CHARACTERS = {
    'χ₃ (mod 3)': {
        'chi': [0, 1, -1], 'q': 3, 'a': 1,
    },
    'χ₄ (mod 4)': {
        'chi': [0, 1, 0, -1], 'q': 4, 'a': 1,
    },
    'χ₅ (mod 5)': {
        'chi': [0, 1, 1j, -1j, -1], 'q': 5, 'a': 1,
    },
}


def completed_L_dirichlet(s, char_info):
    """완비 L-함수: Λ(s, χ) = (q/π)^{s/2} Γ((s+a)/2) L(s, χ)"""
    q = mpmath.mpf(char_info['q'])
    a = mpmath.mpf(char_info['a'])
    L_val = mpmath.dirichlet(s, char_info['chi'])
    gamma_val = mpmath.gamma((s + a) / 2)
    prefactor = mpmath.power(q / mpmath.pi, s / 2)
    return prefactor * gamma_val * L_val


def connection_dirichlet(s, char_info):
    """접속 = Λ'/Λ"""
    h = mpmath.mpf(1) / mpmath.mpf(10**20)
    val = completed_L_dirichlet(s, char_info)
    if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e10, 0)
    d = (completed_L_dirichlet(s + h, char_info) - completed_L_dirichlet(s - h, char_info)) / (2 * h)
    return d / val


# ── 리만 제타용 (기존 방식) ──

def xi_zeta(s):
    """ξ(s) = (1/2) s(s-1) π^{-s/2} Γ(s/2) ζ(s)"""
    half = mpmath.mpf('0.5')
    return half * s * (s - 1) * mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s)


def connection_zeta(s):
    """ξ'/ξ for Riemann zeta"""
    h = mpmath.mpf(1) / mpmath.mpf(10**20)
    val = xi_zeta(s)
    if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e10, 0)
    d = (xi_zeta(s + h) - xi_zeta(s - h)) / (2 * h)
    return d / val


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 통합 인터페이스
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

class LFunctionWrapper:
    """리만 ζ와 디리클레 L-함수를 동일 인터페이스로 다루는 래퍼"""

    def __init__(self, name, char_info=None):
        self.name = name
        self.char_info = char_info
        self.is_zeta = (char_info is None)

    def completed(self, s):
        if self.is_zeta:
            return xi_zeta(s)
        return completed_L_dirichlet(s, self.char_info)

    def connection(self, s):
        if self.is_zeta:
            return connection_zeta(s)
        return connection_dirichlet(s, self.char_info)

    def curvature(self, s):
        L = self.connection(s)
        return float(abs(L)**2)

    def find_zeros(self, t_min, t_max):
        if self.is_zeta:
            return self._find_zeros_zeta(t_min, t_max)
        return self._find_zeros_dirichlet(t_min, t_max)

    def _find_zeros_zeta(self, t_min, t_max):
        """mpmath.zetazero로 리만 영점 구하기"""
        zeros = []
        n = 1
        while True:
            t = float(mpmath.zetazero(n).imag)
            if t > t_max:
                break
            if t >= t_min:
                zeros.append(t)
            n += 1
        return np.array(zeros)

    def _find_zeros_dirichlet(self, t_min, t_max, n_scan=2000):
        """부호 변화 + findroot로 디리클레 영점 탐색"""
        ts = np.linspace(t_min, t_max, n_scan)
        zeros = []

        prev_re, prev_t = None, None
        for t in ts:
            s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
            curr_re = mpmath.re(self.completed(s))
            if prev_re is not None and prev_re * curr_re < 0:
                try:
                    def f_real(t_var, _self=self):
                        sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                        return mpmath.re(_self.completed(sv))
                    tz = float(mpmath.findroot(f_real, (prev_t, t)))
                    if not zeros or abs(tz - zeros[-1]) > 0.1:
                        zeros.append(tz)
                except Exception:
                    pass
            prev_re, prev_t = curr_re, t

        prev_im, prev_t = None, None
        for t in ts:
            s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
            curr_im = mpmath.im(self.completed(s))
            if prev_im is not None and prev_im * curr_im < 0:
                try:
                    def f_imag(t_var, _self=self):
                        sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                        return mpmath.im(_self.completed(sv))
                    tz = float(mpmath.findroot(f_imag, (prev_t, t)))
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(tz))
                    if abs(self.completed(sv)) < mpmath.mpf('1e-10'):
                        if not any(abs(tz - z) < 0.1 for z in zeros):
                            zeros.append(tz)
                except Exception:
                    pass
            prev_im, prev_t = curr_im, t

        zeros.sort()
        return np.array(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 측정 함수들
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def measure_curvature_concentration(lf, zeros, n_generic=50):
    """곡률 집중도: 영점 근방 중앙값 / 일반 점 중앙값"""
    if len(zeros) == 0:
        return None

    # 영점 근방
    zero_k = []
    for tz in zeros:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(tz + 0.01))
        try:
            k = lf.curvature(s)
            if np.isfinite(k):
                zero_k.append(k)
        except Exception:
            pass

    # 일반 점
    gen_k = []
    rng = np.random.RandomState(42)
    attempts = 0
    while len(gen_k) < n_generic and attempts < n_generic * 3:
        t = rng.uniform(T_MIN, T_MAX)
        if len(zeros) > 0 and np.min(np.abs(zeros - t)) < 1.0:
            attempts += 1
            continue
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        try:
            k = lf.curvature(s)
            if np.isfinite(k):
                gen_k.append(k)
        except Exception:
            pass
        attempts += 1

    if not zero_k or not gen_k:
        return None

    return np.median(zero_k) / np.median(gen_k)


def measure_monodromy_accuracy(lf, zeros, eps=0.005):
    """모노드로미 양자화 정확도: |Δarg|가 π에 얼마나 가까운지 (평균 편차)"""
    if len(zeros) == 0:
        return None, None

    deviations = []
    abs_monos = []
    for tz in zeros:
        s_plus = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(tz + eps))
        s_minus = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(tz - eps))

        val_plus = lf.completed(s_plus)
        val_minus = lf.completed(s_minus)

        delta = float(mpmath.arg(val_plus)) - float(mpmath.arg(val_minus))
        while delta > np.pi:
            delta -= 2 * np.pi
        while delta < -np.pi:
            delta += 2 * np.pi

        deviations.append(abs(abs(delta) - np.pi))
        abs_monos.append(abs(delta))

    return np.mean(deviations), np.mean(abs_monos)


def measure_energy_concentration(lf, t_min=12.0, t_max=38.0, n_t=100):
    """에너지 집중도: E(0.5) / E(0.3)"""
    ts = np.linspace(t_min, t_max, n_t)
    dt = (t_max - t_min) / (n_t - 1)

    def compute_E(sigma):
        vals = []
        for t in ts:
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            try:
                k = lf.curvature(s)
                vals.append(k if np.isfinite(k) else 1e6)
            except Exception:
                vals.append(0.0)
        return np.trapezoid(vals, dx=dt)

    E_05 = compute_E(0.5)
    E_03 = compute_E(0.3)

    return E_05 / E_03 if E_03 > 0 else float('inf')


def measure_phase_jump_uniqueness(lf, t_min=10.0, t_max=40.0, n_points=1000):
    """위상 점프 유일성: σ=0.5 점프 수 vs 나머지 σ 최대 점프 수"""
    def count_jumps(sigma):
        ts = np.linspace(t_min, t_max, n_points)
        jumps = 0
        prev_arg = None
        for t in ts:
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            val = lf.completed(s)
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

    j_05 = count_jumps(0.5)
    other_sigmas = [0.1, 0.3, 0.7, 0.9]
    j_others = [count_jumps(s) for s in other_sigmas]
    j_max_other = max(j_others) if j_others else 0

    return j_05, j_max_other


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    out_path = os.path.expanduser(
        '~/Desktop/gdl_unified/outputs/analysis/dirichlet_cross_comparison.txt'
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    print("=" * 70)
    print("교차 비교: Riemann ζ vs χ₃ vs χ₄ vs χ₅")
    print(f"구간: t∈[{T_MIN}, {T_MAX}]")
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}")
    print("=" * 70)

    # L-함수 래퍼 생성
    functions = {
        'Riemann ζ': LFunctionWrapper('Riemann ζ'),
    }
    for name, info in CHARACTERS.items():
        functions[name] = LFunctionWrapper(name, info)

    results = {}

    for fname, lf in functions.items():
        print(f"\n{'─'*50}")
        print(f"  {fname} 분석 중...")
        print(f"{'─'*50}")
        t0 = time.time()

        # 1. 영점 수
        print("  [1/5] 영점 탐색...")
        zeros = lf.find_zeros(T_MIN, T_MAX)
        print(f"        영점 {len(zeros)}개 발견")

        # 2. 곡률 집중도
        print("  [2/5] 곡률 집중도 측정...")
        kappa_ratio = measure_curvature_concentration(lf, zeros)
        print(f"        κ 비율: {kappa_ratio:.1f}×" if kappa_ratio else "        κ 비율: N/A")

        # 3. 모노드로미 양자화
        print("  [3/5] 모노드로미 양자화 검증...")
        mono_dev, mono_mean = measure_monodromy_accuracy(lf, zeros)
        if mono_dev is not None:
            print(f"        평균 |mono|/π = {mono_mean/np.pi:.4f}, 편차 = {mono_dev:.4f}")
        else:
            print("        N/A")

        # 4. 에너지 집중도
        print("  [4/5] 에너지 집중도 측정...")
        e_ratio = measure_energy_concentration(lf)
        print(f"        E(0.5)/E(0.3) = {e_ratio:.1f}×")

        # 5. 위상 점프 유일성
        print("  [5/5] 위상 점프 유일성 검사...")
        j_05, j_max_other = measure_phase_jump_uniqueness(lf)
        print(f"        σ=0.5: {j_05} jumps, 타 σ 최대: {j_max_other} jumps")

        elapsed = time.time() - t0
        print(f"  [소요 시간: {elapsed:.0f}초]")

        results[fname] = {
            'n_zeros': len(zeros),
            'kappa_ratio': kappa_ratio,
            'mono_dev': mono_dev,
            'mono_mean': mono_mean,
            'e_ratio': e_ratio,
            'j_05': j_05,
            'j_max_other': j_max_other,
        }

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 비교 표 출력
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    print(f"\n\n{'='*80}")
    print("교차 비교 결과표")
    print(f"{'='*80}\n")

    # 헤더
    names = list(results.keys())
    col_w = 14
    header = f"{'측정 항목':<25}"
    for n in names:
        header += f"{n:>{col_w}}"
    print(header)
    print("-" * (25 + col_w * len(names)))

    # 1. 영점 수
    row = f"{'영점 수 [10,40]':<25}"
    for n in names:
        row += f"{results[n]['n_zeros']:>{col_w}}"
    print(row)

    # 2. 곡률 집중도
    row = f"{'κ 집중도 (×)':<25}"
    for n in names:
        v = results[n]['kappa_ratio']
        row += f"{f'{v:.1f}':>{col_w}}" if v else f"{'N/A':>{col_w}}"
    print(row)

    # 3. 모노드로미 편차
    row = f"{'mono 편차 (|π-|Δarg||)':<25}"
    for n in names:
        v = results[n]['mono_dev']
        row += f"{f'{v:.4f}':>{col_w}}" if v is not None else f"{'N/A':>{col_w}}"
    print(row)

    # 4. 에너지 집중도
    row = f"{'E(0.5)/E(0.3) (×)':<25}"
    for n in names:
        val = results[n]['e_ratio']
        row += f"{f'{val:.1f}':>{col_w}}"
    print(row)

    # 5. 위상 점프
    row = f"{'σ=0.5 점프 수':<25}"
    for n in names:
        row += f"{results[n]['j_05']:>{col_w}}"
    print(row)

    row = f"{'타 σ 최대 점프 수':<25}"
    for n in names:
        row += f"{results[n]['j_max_other']:>{col_w}}"
    print(row)

    # 판정
    row = f"{'다발 구조 확인':<25}"
    for n in names:
        r = results[n]
        ok = (r['j_05'] > 0 and r['j_max_other'] == 0 and
              r['e_ratio'] > 10 and
              (r['mono_dev'] is None or r['mono_dev'] < 0.5))
        row += f"{'YES':>{col_w}}" if ok else f"{'PARTIAL':>{col_w}}"
    print(row)

    # ── 파일 저장 ──
    with open(out_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("교차 비교: Riemann ζ vs χ₃ vs χ₄ vs χ₅\n")
        f.write(f"구간: t∈[{T_MIN}, {T_MAX}]\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"정밀도: {mpmath.mp.dps} 자릿수\n")
        f.write("=" * 80 + "\n\n")

        # 표
        f.write(header + "\n")
        f.write("-" * (25 + col_w * len(names)) + "\n")

        # 각 행 다시 작성
        rows_data = [
            ('영점 수 [10,40]', lambda n: f"{results[n]['n_zeros']}"),
            ('κ 집중도 (×)', lambda n: f"{results[n]['kappa_ratio']:.1f}" if results[n]['kappa_ratio'] else "N/A"),
            ('mono 편차', lambda n: f"{results[n]['mono_dev']:.4f}" if results[n]['mono_dev'] is not None else "N/A"),
            ('|mono|/π 평균', lambda n: f"{results[n]['mono_mean']/np.pi:.4f}" if results[n]['mono_mean'] is not None else "N/A"),
            ('E(0.5)/E(0.3) (×)', lambda n: f"{results[n]['e_ratio']:.1f}"),
            ('σ=0.5 점프 수', lambda n: f"{results[n]['j_05']}"),
            ('타 σ 최대 점프', lambda n: f"{results[n]['j_max_other']}"),
        ]

        for label, fmt_fn in rows_data:
            row = f"{label:<25}"
            for n in names:
                row += f"{fmt_fn(n):>{col_w}}"
            f.write(row + "\n")

        f.write("\n" + "-" * (25 + col_w * len(names)) + "\n")

        # 판정
        row = f"{'다발 구조 확인':<25}"
        for n in names:
            r = results[n]
            ok = (r['j_05'] > 0 and r['j_max_other'] == 0 and
                  r['e_ratio'] > 10 and
                  (r['mono_dev'] is None or r['mono_dev'] < 0.5))
            row += f"{'YES':>{col_w}}" if ok else f"{'PARTIAL':>{col_w}}"
        f.write(row + "\n")

        # 해석
        f.write(f"\n\n{'='*80}\n")
        f.write("해석\n")
        f.write(f"{'='*80}\n\n")
        f.write("모든 L-함수에서 다음이 확인되면 ξ-다발 프레임워크의 보편성이 입증됨:\n")
        f.write("  1. 위상 점프가 σ=1/2에서만 발생 (타 σ에서 0)\n")
        f.write("  2. 곡률이 영점 근방에서 수십~수백 배 집중\n")
        f.write("  3. 모노드로미가 ±π로 양자화 (편차 < 0.5)\n")
        f.write("  4. 에너지 E(σ)가 σ=1/2에서 뚜렷한 피크 (비율 > 10×)\n\n")
        f.write("이는 RDL의 핵심 추측 — 영점 = 모노드로미, 곡률 = 검출 신호 —이\n")
        f.write("리만 제타를 넘어 일반 L-함수 전체에 적용됨을 시사한다.\n")

    print(f"\n결과 저장: {out_path}")
    print("완료.")


if __name__ == '__main__':
    main()
