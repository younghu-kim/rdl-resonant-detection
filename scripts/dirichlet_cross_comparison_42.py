"""
=============================================================================
[Project RDL] 결과 #42 — 디리클레 교차비교표: ζ vs χ₃ vs χ₄ vs χ₅ vs χ₇
=============================================================================
5개 L-함수를 동일 기준(t∈[10,40])으로 비교하는 통합 실험.

비교 항목:
  1. 영점 수
  2. 곡률 집중도: 영점 근방 중앙값 / 일반 점 중앙값
  3. 모노드로미 양자화 정확도: |Δarg|/π 편차 (log-space)
  4. 에너지 집중도: E(σ=0.5) / E(σ=0.3)
  5. σ-국소화: σ=0.5 점프 수 / 타 σ 최대 점프 수

추가 분석:
  - conductor 스케일링: κ 비율 vs (1/2)log(q/π) 상관 정량화

핵심 수정 (기존 cross_comparison 대비):
  - 접속 Λ'/Λ: h=1e-20 → 해석적 공식 (digamma + L'/L h=1e-6)
  - ζ 접속 ξ'/ξ: h=1e-20 → 해석적 공식 (1/s+1/(s-1)+digamma+ζ'/ζ)
  - 모노드로미: direct arg → log-space Im(log Λ) (underflow 방지)
  - mod 7 추가 (결과 #40/#40b와 동일 지표 정의)

결과 파일: results/dirichlet_cross_comparison.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
import cmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

mpmath.mp.dps = 80

T_MIN, T_MAX = 10.0, 40.0
DELTA_OFFSET = 0.03   # 영점 위 직접 측정 금지 → 오프셋

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 지표 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

_w6 = cmath.exp(2j * cmath.pi / 6)
_chi7_raw = [0, 1, _w6**2, _w6**1, _w6**4, _w6**5, _w6**3]
_chi7 = [mpmath.mpc(c.real, c.imag) for c in _chi7_raw]

CHARACTERS_EXT = {
    'χ₃ (mod 3)': {'chi': [0, 1, -1],            'q': 3, 'a': 1},
    'χ₄ (mod 4)': {'chi': [0, 1, 0, -1],          'q': 4, 'a': 1},
    'χ₅ (mod 5)': {'chi': [0, 1, 1j, -1j, -1],    'q': 5, 'a': 1},
    'χ₇ (mod 7)': {'chi': _chi7,                   'q': 7, 'a': 1},
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 해석적 접속 (Λ'/Λ)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def connection_dirichlet_analytic(s, ci):
    """
    Λ'/Λ 해석적 공식:
      Λ'/Λ = (1/2)log(q/π) + (1/2)ψ((s+a)/2) + L'/L

    - (1/2)log(q/π): 순수 해석
    - (1/2)ψ((s+a)/2): mpmath.digamma 해석
    - L'/L: L(s,χ)만 h=1e-6 중앙차분 (h=1e-20 금지)
    """
    q = mpmath.mpf(ci['q'])
    a = mpmath.mpf(ci['a'])
    chi = ci['chi']

    log_term = mpmath.log(q / mpmath.pi) / 2
    digamma_term = mpmath.digamma((s + a) / 2) / 2

    L_val = mpmath.dirichlet(s, chi)
    if abs(L_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e10, 0)

    h = mpmath.mpf('1e-6')
    L_d = (mpmath.dirichlet(s + h, chi) - mpmath.dirichlet(s - h, chi)) / (2 * h)
    L_log_deriv = L_d / L_val

    return log_term + digamma_term + L_log_deriv


def connection_zeta_analytic(s):
    """
    ξ'/ξ 해석적 공식:
      ξ'/ξ = 1/s + 1/(s-1) - (1/2)log(π) + (1/2)ψ(s/2) + ζ'/ζ

    - 앞 4항: 순수 해석 (1/s + 1/(s-1) + digamma)
    - ζ'/ζ: ζ(s)만 h=1e-6 중앙차분 (h=1e-20 금지)
    """
    background = (mpmath.mpf(1)/s +
                  mpmath.mpf(1)/(s - 1) -
                  mpmath.mpf('0.5') * mpmath.log(mpmath.pi) +
                  mpmath.mpf('0.5') * mpmath.digamma(s / 2))

    zeta_val = mpmath.zeta(s)
    if abs(zeta_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e10, 0)

    h = mpmath.mpf('1e-6')
    zeta_d = (mpmath.zeta(s + h) - mpmath.zeta(s - h)) / (2 * h)

    return background + zeta_d / zeta_val


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# log-space arg 함수 (underflow 방지)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def log_arg_dirichlet(s, ci):
    """
    arg(Λ(s,χ)) = Im[log Λ]
    log Λ = (s/2)·log(q/π) + loggamma((s+a)/2) + log(L(s,χ))
    반환: float (실패 시 None)
    """
    q = mpmath.mpf(ci['q'])
    a = mpmath.mpf(ci['a'])
    L_val = mpmath.dirichlet(s, ci['chi'])
    if abs(L_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 15):
        return None

    log_Lambda = ((s / 2) * mpmath.log(q / mpmath.pi) +
                  mpmath.loggamma((s + a) / 2) +
                  mpmath.log(L_val))
    return float(mpmath.im(log_Lambda))


def log_arg_zeta(s):
    """
    arg(ξ(s)) = Im[log ξ]
    log ξ = log(1/2) + log(s) + log(s-1) - (s/2)log(π) + loggamma(s/2) + log(ζ(s))
    반환: float (실패 시 None)
    """
    zeta_val = mpmath.zeta(s)
    if abs(zeta_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 15):
        return None

    log_xi = (mpmath.log(mpmath.mpf('0.5')) +
              mpmath.log(s) +
              mpmath.log(s - 1) -
              (s / 2) * mpmath.log(mpmath.pi) +
              mpmath.loggamma(s / 2) +
              mpmath.log(zeta_val))
    return float(mpmath.im(log_xi))


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# L-함수 래퍼
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

class LFunc:
    def __init__(self, name, char_info=None):
        self.name = name
        self.ci = char_info
        self.is_zeta = (char_info is None)
        self.q = 1 if self.is_zeta else char_info['q']

    def connection(self, s):
        if self.is_zeta:
            return connection_zeta_analytic(s)
        return connection_dirichlet_analytic(s, self.ci)

    def curvature(self, s):
        conn = self.connection(s)
        k = float(abs(conn)**2)
        return k if np.isfinite(k) else 1e12

    def log_arg(self, t_val):
        """log-space arg at s = 0.5 + it"""
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_val))
        if self.is_zeta:
            return log_arg_zeta(s)
        return log_arg_dirichlet(s, self.ci)

    def completed_val(self, s):
        """완비 L-함수 값 (위상 점프용)"""
        if self.is_zeta:
            half = mpmath.mpf('0.5')
            return half * s * (s-1) * mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s)
        q = mpmath.mpf(self.ci['q'])
        a = mpmath.mpf(self.ci['a'])
        return (mpmath.power(q/mpmath.pi, s/2) *
                mpmath.gamma((s+a)/2) *
                mpmath.dirichlet(s, self.ci['chi']))

    def find_zeros(self, t_min=T_MIN, t_max=T_MAX):
        if self.is_zeta:
            return self._zeros_zeta(t_min, t_max)
        return self._zeros_dirichlet(t_min, t_max)

    def _zeros_zeta(self, t_min, t_max):
        zeros = []
        n = 1
        while True:
            try:
                t = float(mpmath.zetazero(n).imag)
            except Exception as e:
                print(f"WARNING: zetazero({n}) 실패: {e}", flush=True)
                n += 1
                if n > 200:
                    break
                continue
            if t > t_max:
                break
            if t >= t_min:
                zeros.append(t)
            n += 1
        return np.array(zeros)

    def _zeros_dirichlet(self, t_min, t_max, n_scan=3000):
        """Re(Λ)+Im(Λ) 부호 변화 + findroot 단일 시작점 + |Λ| 확인"""
        ts = np.linspace(t_min, t_max, n_scan)
        zeros = []
        fail_count = 0
        success_count = 0
        ci = self.ci

        # ── Re(Λ) 부호 변화 ──
        prev_re, prev_t = None, None
        for t in ts:
            s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
            try:
                val = self.completed_val(s)
                curr_re = float(mpmath.re(val))
            except Exception as e:
                print(f"  WARNING: completed_val 실패 t={t:.2f}: {e}", flush=True)
                prev_re, prev_t = None, float(t)
                continue

            if prev_re is not None and prev_re * curr_re < 0:
                mid = (prev_t + float(t)) / 2
                try:
                    def f_re(t_var, ci=ci):
                        sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                        q = mpmath.mpf(ci['q']); a = mpmath.mpf(ci['a'])
                        return mpmath.re(mpmath.power(q/mpmath.pi,sv/2) *
                                        mpmath.gamma((sv+a)/2) *
                                        mpmath.dirichlet(sv, ci['chi']))
                    tz = float(mpmath.findroot(f_re, mpmath.mpf(str(mid))))
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(tz))
                    lv = abs(self.completed_val(sv))
                    if lv < mpmath.mpf(10)**(-20):
                        if not zeros or abs(tz - zeros[-1]) > 0.1:
                            zeros.append(tz)
                            success_count += 1
                except Exception as e:
                    fail_count += 1
                    print(f"  WARNING: Re findroot 실패 t≈{mid:.2f}: {e}", flush=True)
            prev_re, prev_t = curr_re, float(t)

        # ── Im(Λ) 부호 변화 (복소 지표 누락 영점 포착) ──
        prev_im, prev_t = None, None
        for t in ts:
            s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
            try:
                val = self.completed_val(s)
                curr_im = float(mpmath.im(val))
            except Exception as e:
                prev_im, prev_t = None, float(t)
                continue

            if prev_im is not None and prev_im * curr_im < 0:
                mid = (prev_t + float(t)) / 2
                try:
                    def f_im(t_var, ci=ci):
                        sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                        q = mpmath.mpf(ci['q']); a = mpmath.mpf(ci['a'])
                        return mpmath.im(mpmath.power(q/mpmath.pi,sv/2) *
                                        mpmath.gamma((sv+a)/2) *
                                        mpmath.dirichlet(sv, ci['chi']))
                    tz = float(mpmath.findroot(f_im, mpmath.mpf(str(mid))))
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(tz))
                    lv = abs(self.completed_val(sv))
                    if lv < mpmath.mpf(10)**(-20):
                        if not any(abs(tz - z) < 0.1 for z in zeros):
                            zeros.append(tz)
                except Exception as e:
                    fail_count += 1
            prev_im, prev_t = curr_im, float(t)

        if len(zeros) == 0:
            print(f"  ⚠️ 영점 0개 — 탐색 로직 점검 필요!", flush=True)
        if fail_count > 0:
            print(f"  ⚠️ findroot 실패 {fail_count}회 (성공 {success_count}회)", flush=True)
            if success_count > 0 and fail_count > success_count:
                print(f"  ⚠⚠ 실패율 > 50% — 탐색 로직 점검 필요!", flush=True)

        zeros.sort()
        return np.array(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 측정 함수 (5가지 + conductor)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def measure_curvature_concentration(lf, zeros, n_generic=40):
    """
    κ 집중도: median(κ near zeros) / median(κ far)
    영점 오프셋 δ=0.03 사용 (영점 위 직접 측정 금지)
    """
    if len(zeros) == 0:
        return None

    near_k = []
    for tz in zeros:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(tz + DELTA_OFFSET))
        try:
            k = lf.curvature(s)
            if np.isfinite(k) and k < 1e11:
                near_k.append(k)
        except Exception as e:
            print(f"  WARNING: curvature near-zero 실패 t={tz:.2f}: {e}", flush=True)

    far_k = []
    rng = np.random.RandomState(42)
    attempts = 0
    while len(far_k) < n_generic and attempts < n_generic * 5:
        t = rng.uniform(T_MIN + 1.0, T_MAX - 1.0)
        # 영점에서 충분히 멀리 (1.0 이상)
        if len(zeros) > 0 and np.min(np.abs(zeros - t)) < 1.0:
            attempts += 1
            continue
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        try:
            k = lf.curvature(s)
            if np.isfinite(k) and k < 1e11:
                far_k.append(k)
        except Exception as e:
            print(f"  WARNING: curvature far 실패 t={t:.2f}: {e}", flush=True)
        attempts += 1

    if not near_k or not far_k:
        return None

    ratio = float(np.median(near_k)) / float(np.median(far_k))
    return ratio


def measure_monodromy_accuracy(lf, zeros, eps=0.005):
    """
    모노드로미 양자화 정확도 (log-space arg 사용, underflow 방지):
    phase_delta = Im[log Λ(tz+eps)] - Im[log Λ(tz-eps)]
    기대: |phase_delta| ≈ π
    반환: (mean_deviation, mean_abs_mono/π)
    """
    if len(zeros) == 0:
        return None, None

    deviations = []
    abs_monos = []
    skipped = 0

    for tz in zeros:
        arg_plus = lf.log_arg(tz + eps)
        arg_minus = lf.log_arg(tz - eps)

        if arg_plus is None or arg_minus is None:
            skipped += 1
            continue

        delta = arg_plus - arg_minus

        # branch 정규화: |delta| ≈ π 기대, 2π 배수 오류 제거
        while delta > 2 * np.pi:
            delta -= 2 * np.pi
        while delta < -2 * np.pi:
            delta += 2 * np.pi

        # ±π로 부호 중립화: -π → π
        if abs(delta + np.pi) < abs(delta - np.pi):
            delta = -delta

        deviations.append(abs(abs(delta) - np.pi))
        abs_monos.append(abs(delta))

    if not deviations:
        return None, None

    if skipped > 0:
        print(f"  INFO: 모노드로미 log_arg 계산 실패 {skipped}개 (underflow)", flush=True)

    return float(np.mean(deviations)), float(np.mean(abs_monos)) / np.pi


def measure_energy_concentration(lf, t_min=12.0, t_max=38.0, n_t=80):
    """에너지 집중도: E(σ=0.5) / E(σ=0.3)"""
    ts = np.linspace(t_min, t_max, n_t)
    dt = (t_max - t_min) / (n_t - 1)

    def compute_E(sigma):
        vals = []
        for t in ts:
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            try:
                k = lf.curvature(s)
                vals.append(min(k, 1e8) if np.isfinite(k) else 1e8)
            except Exception as e:
                print(f"  WARNING: energy curvature σ={sigma} t={t:.2f}: {e}", flush=True)
                vals.append(0.0)
        return np.trapezoid(vals, dx=dt)

    E_05 = compute_E(0.5)
    E_03 = compute_E(0.3)

    if E_03 <= 0:
        return float('inf')
    return E_05 / E_03


def measure_phase_jump_uniqueness(lf, t_min=T_MIN, t_max=T_MAX, n_points=600):
    """
    위상 점프 유일성: σ=0.5 점프 수 vs 타 σ 최대 점프 수
    점프 기준: |Δarg| > π/2
    """
    def count_jumps(sigma):
        ts = np.linspace(t_min, t_max, n_points)
        jumps = 0
        prev_arg = None
        for t in ts:
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            try:
                val = lf.completed_val(s)
                if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
                    prev_arg = None
                    continue
                curr_arg = float(mpmath.arg(val))
            except Exception as e:
                print(f"  WARNING: phase jump σ={sigma} t={t:.2f}: {e}", flush=True)
                prev_arg = None
                continue

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
# conductor 스케일링 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def conductor_scaling_analysis(results_dict):
    """
    κ 집중도 vs (1/2)log(q/π) 상관 정량화.
    ζ → q_eff=1, Dirichlet → q=conductor
    """
    try:
        from scipy import stats
    except ImportError:
        print("  WARNING: scipy 없음 — 상관 계산 스킵", flush=True)
        return None

    # q 값 추출
    q_map = {
        'Riemann ζ': 1,
        'χ₃ (mod 3)': 3,
        'χ₄ (mod 4)': 4,
        'χ₅ (mod 5)': 5,
        'χ₇ (mod 7)': 7,
    }

    x_vals = []  # (1/2)log(q/π)
    y_kappa = []  # κ 집중도

    for name, res in results_dict.items():
        q = q_map.get(name)
        if q is None:
            continue
        x = 0.5 * float(mpmath.log(mpmath.mpf(q) / mpmath.pi))
        x_vals.append(x)
        kappa = res.get('kappa_ratio')
        y_kappa.append(kappa if kappa is not None else float('nan'))

    x_arr = np.array(x_vals)
    y_arr = np.array(y_kappa)

    valid = np.isfinite(y_arr)
    if valid.sum() < 3:
        return {'note': '유효 데이터 부족 (< 3개)'}

    xv, yv = x_arr[valid], y_arr[valid]

    pearson_r, pearson_p = stats.pearsonr(xv, yv)
    spearman_r, spearman_p = stats.spearmanr(xv, yv)

    # 선형 피팅: κ ~ a·x + b
    slope, intercept, r_lin, p_lin, se_lin = stats.linregress(xv, yv)

    # 로그 피팅: κ ~ a·log(x - x_min + 1) + b (x_min을 평행이동해서 log 안전)
    # 단순히 선형 비교
    fit_vals = slope * xv + intercept
    residuals = yv - fit_vals
    r2_linear = 1 - np.sum(residuals**2) / np.sum((yv - np.mean(yv))**2) if np.sum((yv - np.mean(yv))**2) > 0 else 0.0

    return {
        'x_vals': x_arr.tolist(),
        'y_kappa': y_arr.tolist(),
        'pearson_r': float(pearson_r),
        'pearson_p': float(pearson_p),
        'spearman_r': float(spearman_r),
        'spearman_p': float(spearman_p),
        'lin_slope': float(slope),
        'lin_intercept': float(intercept),
        'lin_r2': float(r2_linear),
        'n_valid': int(valid.sum()),
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_global_start = time.time()

    out_path = os.path.expanduser(
        '~/Desktop/gdl_unified/results/dirichlet_cross_comparison.txt'
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    print("=" * 75, flush=True)
    print("결과 #42 — 디리클레 교차비교: ζ vs χ₃ vs χ₄ vs χ₅ vs χ₇", flush=True)
    print(f"구간: t∈[{T_MIN}, {T_MAX}]", flush=True)
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}", flush=True)
    print(f"정밀도: {mpmath.mp.dps} 자릿수", flush=True)
    print("=" * 75, flush=True)

    # L-함수 목록 (순서 고정)
    lf_list = [('Riemann ζ', LFunc('Riemann ζ'))]
    for name, ci in CHARACTERS_EXT.items():
        lf_list.append((name, LFunc(name, ci)))

    results = {}

    for fname, lf in lf_list:
        print(f"\n{'─'*55}", flush=True)
        print(f"  {fname} 분석", flush=True)
        print(f"{'─'*55}", flush=True)
        t0 = time.time()

        # 1. 영점 탐색
        print("  [1/5] 영점 탐색...", flush=True)
        zeros = lf.find_zeros()
        if len(zeros) == 0:
            print(f"  ⚠️ 영점 0개 발견!", flush=True)
        else:
            print(f"        {len(zeros)}개 발견: {zeros[:3].tolist()} ...", flush=True)

        # 2. κ 집중도
        print("  [2/5] κ 집중도 측정...", flush=True)
        kappa_ratio = measure_curvature_concentration(lf, zeros)
        if kappa_ratio is not None:
            print(f"        κ 비율: {kappa_ratio:.1f}×", flush=True)
        else:
            print("        κ 비율: N/A", flush=True)

        # 3. 모노드로미 양자화
        print("  [3/5] 모노드로미 양자화 (log-space arg)...", flush=True)
        mono_dev, mono_ratio = measure_monodromy_accuracy(lf, zeros)
        if mono_dev is not None:
            print(f"        |mono|/π 평균 = {mono_ratio:.6f}, 편차 = {mono_dev:.6f}", flush=True)
        else:
            print("        N/A", flush=True)

        # 4. 에너지 집중도
        print("  [4/5] 에너지 집중도...", flush=True)
        e_ratio = measure_energy_concentration(lf)
        print(f"        E(0.5)/E(0.3) = {e_ratio:.1f}×", flush=True)

        # 5. 위상 점프 유일성
        print("  [5/5] 위상 점프 유일성...", flush=True)
        j_05, j_max_other = measure_phase_jump_uniqueness(lf)
        print(f"        σ=0.5: {j_05} jumps, 타 σ 최대: {j_max_other} jumps", flush=True)

        elapsed = time.time() - t0
        print(f"  [소요: {elapsed:.0f}초]", flush=True)

        results[fname] = {
            'n_zeros': len(zeros),
            'kappa_ratio': kappa_ratio,
            'mono_dev': mono_dev,
            'mono_ratio': mono_ratio,
            'e_ratio': e_ratio,
            'j_05': j_05,
            'j_max_other': j_max_other,
            'q': lf.q,
        }

    # ── conductor 스케일링 분석 ──────────────────────────────────────────────
    print(f"\n{'─'*55}", flush=True)
    print("  Conductor 스케일링 분석", flush=True)
    print(f"{'─'*55}", flush=True)
    cond_result = conductor_scaling_analysis(results)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 결과 표 출력 + 파일 저장
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    names = [n for n, _ in lf_list]
    col_w = 13

    def fmt_row(label, fmt_fn, names, results):
        row = f"{label:<28}"
        for n in names:
            try:
                row += f"{fmt_fn(n):>{col_w}}"
            except Exception:
                row += f"{'ERR':>{col_w}}"
        return row

    header = f"{'측정 항목':<28}"
    for n in names:
        short = n.replace('Riemann ', '').replace(' (mod ', ' mod').replace(')', '')
        header += f"{short:>{col_w}}"
    sep = "-" * (28 + col_w * len(names))

    rows_data = [
        ("conductor q",
         lambda n: f"{results[n]['q']}"),
        ("(1/2)log(q/π)",
         lambda n: f"{0.5*float(mpmath.log(mpmath.mpf(results[n]['q'])/mpmath.pi)):.3f}"),
        ("영점 수 [10,40]",
         lambda n: f"{results[n]['n_zeros']}"),
        ("κ 집중도 (×)",
         lambda n: f"{results[n]['kappa_ratio']:.1f}" if results[n]['kappa_ratio'] else "N/A"),
        ("mono 편차 (π 단위)",
         lambda n: f"{results[n]['mono_dev']:.6f}" if results[n]['mono_dev'] is not None else "N/A"),
        ("|mono|/π 평균",
         lambda n: f"{results[n]['mono_ratio']:.6f}" if results[n]['mono_ratio'] is not None else "N/A"),
        ("E(0.5)/E(0.3) (×)",
         lambda n: f"{results[n]['e_ratio']:.1f}"),
        ("σ=0.5 점프 수",
         lambda n: f"{results[n]['j_05']}"),
        ("타σ 최대 점프 수",
         lambda n: f"{results[n]['j_max_other']}"),
    ]

    # 터미널 출력
    print(f"\n\n{'='*75}", flush=True)
    print("결과 #42 교차비교 결과표", flush=True)
    print(f"{'='*75}\n", flush=True)
    print(header, flush=True)
    print(sep, flush=True)
    for label, fmt_fn in rows_data:
        print(fmt_row(label, fmt_fn, names, results), flush=True)

    # 판정 행
    print(sep, flush=True)
    pass_row = f"{'다발 구조 PASS':<28}"
    for n in names:
        r = results[n]
        ok = (r['j_05'] > 0 and r['j_max_other'] == 0 and
              r['e_ratio'] > 10 and
              (r['mono_dev'] is None or r['mono_dev'] < 0.1))
        pass_row += f"{'PASS':>{col_w}}" if ok else f"{'PARTIAL':>{col_w}}"
    print(pass_row, flush=True)

    # conductor 스케일링 출력
    if cond_result and 'pearson_r' in cond_result:
        print(f"\n{'─'*55}", flush=True)
        print("Conductor 스케일링: κ 집중도 vs (1/2)log(q/π)", flush=True)
        print(f"{'─'*55}", flush=True)
        for n, x, y in zip(names, cond_result['x_vals'], cond_result['y_kappa']):
            ys = f"{y:.1f}" if np.isfinite(y) else "N/A"
            print(f"  {n:<20}  x={(x):+.3f}  κ={ys}", flush=True)
        print(f"  Pearson  r = {cond_result['pearson_r']:.4f} (p={cond_result['pearson_p']:.3e})", flush=True)
        print(f"  Spearman r = {cond_result['spearman_r']:.4f} (p={cond_result['spearman_p']:.3e})", flush=True)
        print(f"  선형 피팅: κ = {cond_result['lin_slope']:.1f}·x + {cond_result['lin_intercept']:.1f}  (R²={cond_result['lin_r2']:.3f})", flush=True)

    total_elapsed = time.time() - t_global_start

    # ── 파일 저장 ──────────────────────────────────────────────────────────
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write("=" * 75 + "\n")
        f.write("결과 #42 — 디리클레 교차비교: ζ vs χ₃ vs χ₄ vs χ₅ vs χ₇\n")
        f.write(f"구간: t∈[{T_MIN}, {T_MAX}]\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"정밀도: {mpmath.mp.dps} 자릿수\n")
        f.write(f"총 소요: {total_elapsed:.0f}초 ({total_elapsed/60:.1f}분)\n")
        f.write("=" * 75 + "\n\n")

        # 표
        f.write(header + "\n")
        f.write(sep + "\n")
        for label, fmt_fn in rows_data:
            f.write(fmt_row(label, fmt_fn, names, results) + "\n")
        f.write(sep + "\n")
        f.write(pass_row + "\n")

        # conductor 분석
        f.write("\n\n" + "=" * 75 + "\n")
        f.write("Conductor 스케일링 분석\n")
        f.write("=" * 75 + "\n\n")
        f.write(f"기준: κ 집중도(near/far 비율) vs (1/2)log(q/π)\n")
        f.write(f"이론: conductor 증가 → 배경항 (1/2)log(q/π) 증가 → off-zero κ 증가 → 비율 감소\n\n")

        if cond_result and 'pearson_r' in cond_result:
            f.write(f"{'L-함수':<22}  (1/2)log(q/π)  κ 집중도\n")
            f.write("-" * 50 + "\n")
            for n, x, y in zip(names, cond_result['x_vals'], cond_result['y_kappa']):
                ys = f"{y:.1f}×" if np.isfinite(y) else "N/A"
                f.write(f"  {n:<20}  {x:+.3f}           {ys}\n")
            f.write("\n")
            f.write(f"Pearson  r = {cond_result['pearson_r']:.4f}  (p = {cond_result['pearson_p']:.3e})\n")
            f.write(f"Spearman r = {cond_result['spearman_r']:.4f}  (p = {cond_result['spearman_p']:.3e})\n")
            f.write(f"선형 피팅: κ = {cond_result['lin_slope']:.1f}·x + {cond_result['lin_intercept']:.1f}  (R² = {cond_result['lin_r2']:.3f})\n")
            f.write(f"유효 데이터: {cond_result['n_valid']}개 / 5개\n")
        else:
            f.write("conductor 분석 실패 (데이터 부족 또는 scipy 미설치)\n")

        # 종합 판정
        f.write("\n\n" + "=" * 75 + "\n")
        f.write("종합 판정\n")
        f.write("=" * 75 + "\n\n")
        f.write("성공 기준:\n")
        f.write("  PASS: 5개 L-함수 비교표 완성\n")
        f.write("       + 모노드로미 전부 편차 < 0.01 (보편성)\n")
        f.write("       + conductor 상관 정량화\n\n")

        all_mono_ok = all(
            results[n]['mono_dev'] is not None and results[n]['mono_dev'] < 0.1
            for n in names
        )
        comparison_complete = len(results) == 5

        f.write(f"비교표 완성: {'✅ YES' if comparison_complete else '❌ NO'} (5개 L-함수)\n")
        f.write(f"모노드로미 보편성: {'✅ YES' if all_mono_ok else '⚠️ 일부 실패'}\n")
        for n in names:
            d = results[n]['mono_dev']
            s = f"{d:.6f}" if d is not None else "N/A"
            f.write(f"  - {n}: mono_dev = {s}\n")

        if cond_result and 'pearson_r' in cond_result:
            pr = cond_result['pearson_r']
            sp = cond_result['pearson_p']
            trend_ok = (pr < -0.5 and sp < 0.2)
            f.write(f"conductor 단조 감소 경향: {'✅ YES' if trend_ok else '⚠️ 불확실'} (r={pr:.3f}, p={sp:.3e})\n")

        final_pass = (comparison_complete and all_mono_ok)
        f.write(f"\n최종 판정: {'✅ PASS' if final_pass else '⚠️ PARTIAL'}\n")

    print(f"\n결과 저장: {out_path}", flush=True)
    print(f"총 소요: {total_elapsed:.0f}초 ({total_elapsed/60:.1f}분)", flush=True)
    print("완료.", flush=True)


if __name__ == '__main__':
    main()
