"""
=============================================================================
[Project RDL] 다발 기하학 공용 유틸리티
=============================================================================
모든 다발/디리클레 실험 스크립트가 공유하는 핵심 함수.
새 스크립트 작성 시 이 모듈을 import하여 사용.

사용법:
    from bundle_utils import (
        xi_func, connection_zeta, curvature_zeta, monodromy_contour,
        find_zeros_zeta, find_zeros_dirichlet,
        completed_L, connection_dirichlet, curvature_dirichlet,
        evaluate_predictions, CHARACTERS,
    )

주의사항 (이 모듈이 해결하는 반복 실수들):
    - 모노드로미: eps 차분 금지 → 폐곡선 적분 (monodromy_contour)
    - 영점 판정: 1e-40 절대값 금지 → dps 기반 상대 판정
    - mpmath dps: t>100이면 80 이상 필요
    - numpy 2.0: trapezoid 사용
    - 인덱싱: dtype=int 강제
=============================================================================
"""

import numpy as np
import mpmath

# dps 설정은 호출자가 하되, 기본값 제공
if mpmath.mp.dps < 50:
    mpmath.mp.dps = 80


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 디리클레 지표 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CHARACTERS = {
    'chi_mod_3': {
        'chi': [0, 1, -1],
        'q': 3, 'a': 1,
        'label': 'χ₃ (mod 3, 비자명)',
    },
    'chi_mod_4': {
        'chi': [0, 1, 0, -1],
        'q': 4, 'a': 1,
        'label': 'χ₄ (mod 4, 실수)',
    },
    'chi_mod_5': {
        'chi': [0, 1, 1j, -1j, -1],
        'q': 5, 'a': 1,
        'label': 'χ₅ (mod 5, 복소)',
    },
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 리만 제타용 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _near_zero(val):
    """dps 기반 상대적 영점 판정 (1e-40 같은 절대값 금지)"""
    return abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10)


def xi_func(s):
    """ξ(s) = (1/2) s(s-1) π^(-s/2) Γ(s/2) ζ(s)"""
    half = mpmath.mpf('0.5')
    return half * s * (s - 1) * mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s)


def connection_zeta(s):
    """접속 L(s) = ξ'/ξ (수치 미분)"""
    h = mpmath.mpf(1) / mpmath.mpf(10**20)
    xi_val = xi_func(s)
    if _near_zero(xi_val):
        return mpmath.mpc(1e10, 0)
    xi_d = (xi_func(s + h) - xi_func(s - h)) / (2 * h)
    return xi_d / xi_val


def curvature_zeta(s):
    """곡률 κ = |L|² = |ξ'/ξ|²"""
    L = connection_zeta(s)
    return float(abs(L)**2)


def curvature_at_t(t):
    """임계선 위 점 t에서의 곡률 (편의 함수)"""
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
    return curvature_zeta(s)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 디리클레 L-함수용 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def completed_L(s, char_info):
    """완비 L-함수: Λ(s, χ) = (q/π)^{s/2} Γ((s+a)/2) L(s, χ)"""
    q = mpmath.mpf(char_info['q'])
    a = mpmath.mpf(char_info['a'])
    chi = char_info['chi']
    L_val = mpmath.dirichlet(s, chi)
    gamma_val = mpmath.gamma((s + a) / 2)
    prefactor = mpmath.power(q / mpmath.pi, s / 2)
    return prefactor * gamma_val * L_val


def connection_dirichlet(s, char_info):
    """디리클레 접속 L(s) = Λ'/Λ"""
    h = mpmath.mpf(1) / mpmath.mpf(10**20)
    Lambda_val = completed_L(s, char_info)
    if _near_zero(Lambda_val):
        return mpmath.mpc(1e10, 0)
    Lambda_d = (completed_L(s + h, char_info) - completed_L(s - h, char_info)) / (2 * h)
    return Lambda_d / Lambda_val


def curvature_dirichlet(s, char_info):
    """디리클레 곡률 κ = |Λ'/Λ|²"""
    L = connection_dirichlet(s, char_info)
    return float(abs(L)**2)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 모노드로미 (폐곡선 적분 — eps 차분 금지)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def monodromy_contour(t, radius=0.5, n_steps=64, func='zeta', char_info=None):
    """
    점 s=1/2+it 주위 반지름 radius 원에서 모노드로미 계산.
    영점이 원 안에 있으면 ≈±2π, 없으면 ≈0.

    Parameters:
        t: 임계선 위 좌표
        radius: 폐곡선 반지름
        n_steps: 적분 단계 수
        func: 'zeta' 또는 'dirichlet'
        char_info: func='dirichlet'일 때 지표 정보 dict
    """
    s_center = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))

    total_delta = 0.0
    prev_arg = None

    for k in range(n_steps + 1):
        theta = 2 * np.pi * k / n_steps
        s = s_center + radius * mpmath.exp(1j * theta)

        if func == 'zeta':
            val = xi_func(s)
        else:
            val = completed_L(s, char_info)

        if _near_zero(val):
            continue

        curr_arg = float(mpmath.arg(val))

        if prev_arg is not None:
            delta = curr_arg - prev_arg
            while delta > np.pi:
                delta -= 2 * np.pi
            while delta < -np.pi:
                delta += 2 * np.pi
            total_delta += delta

        prev_arg = curr_arg

    return total_delta


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 영점 탐색
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def find_zeros_zeta(t_min, t_max):
    """mpmath로 구간 내 리만 제타 영점 반환"""
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


def find_zeros_dirichlet(char_info, t_min=10.0, t_max=40.0, n_scan=2000):
    """
    임계선 위에서 Re(Λ)/Im(Λ) 부호 변화로 영점 탐색 + findroot 정밀화.
    검증됨: χ mod 3/4/5에서 각각 12-13개 영점 정확히 발견 (verification 결과).
    """
    ts = np.linspace(t_min, t_max, n_scan)
    zeros = []

    # Re(Λ) 부호 변화
    prev_re, prev_t = None, None
    for t in ts:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        val = completed_L(s, char_info)
        curr_re = mpmath.re(val)

        if prev_re is not None and prev_re * curr_re < 0:
            try:
                def f_real(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    return mpmath.re(completed_L(sv, char_info))
                t_zero = float(mpmath.findroot(f_real, (prev_t, t)))
                # |Λ| 검증: Re(Λ)=0 교차가 진짜 영점인지 확인
                # (실수 지표에서 Λ이 순허수일 때 Re≈0 항상 → 가짜 교차 방지)
                sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero))
                if abs(completed_L(sv, char_info)) < mpmath.mpf('1e-10'):
                    if not zeros or abs(t_zero - zeros[-1]) > 0.1:
                        zeros.append(t_zero)
            except Exception:
                pass
        prev_re, prev_t = curr_re, t

    # Im(Λ) 부호 변화 (일부 영점은 Re만으로 못 잡음)
    # 실수 지표에서는 Im(Λ)≈0 상시 → 스퓨리어스 부호변화 혼입. 비활성화.
    is_real_char = all(isinstance(v, (int, float)) or (hasattr(v, 'imag') and v.imag == 0)
                       for v in char_info['chi'])
    if is_real_char:
        zeros.sort()
        return np.array(zeros)

    prev_im, prev_t = None, None
    for t in ts:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        val = completed_L(s, char_info)
        curr_im = mpmath.im(val)

        if prev_im is not None and prev_im * curr_im < 0:
            try:
                def f_imag(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    return mpmath.im(completed_L(sv, char_info))
                t_zero = float(mpmath.findroot(f_imag, (prev_t, t)))
                sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero))
                if abs(completed_L(sv, char_info)) < mpmath.mpf('1e-10'):
                    if not any(abs(t_zero - z) < 0.1 for z in zeros):
                        zeros.append(t_zero)
            except Exception:
                pass
        prev_im, prev_t = curr_im, t

    zeros.sort()
    return np.array(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 곡률 프로파일 & 피크 탐색
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_curvature_profile(t_min, t_max, n_points, func='zeta', char_info=None):
    """t 구간에서 곡률+위상 프로파일 계산"""
    ts = np.linspace(t_min, t_max, n_points)
    kappas = np.zeros(n_points)
    args = np.zeros(n_points)

    for i, t in enumerate(ts):
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))

        if func == 'zeta':
            val = xi_func(s)
        else:
            val = completed_L(s, char_info)

        if _near_zero(val):
            kappas[i] = 1e10
            args[i] = 0
            continue

        if func == 'zeta':
            L = connection_zeta(s)
        else:
            L = connection_dirichlet(s, char_info)

        kappas[i] = float(abs(L)**2)
        args[i] = float(mpmath.arg(val))

    return ts, kappas, args


def find_curvature_peaks(ts, kappas, min_prominence=5.0):
    """곡률 극대점 찾기 (dtype=int 강제)"""
    peaks = []
    median_k = np.median(kappas[kappas < 1e9])  # capped 값 제외
    if median_k == 0:
        median_k = 1.0
    for i in range(1, len(kappas) - 1):
        if kappas[i] > kappas[i-1] and kappas[i] > kappas[i+1]:
            if kappas[i] > min_prominence * median_k:
                peaks.append(i)
    return np.array(peaks, dtype=int)  # 항상 int!


def compute_monodromy_profile(ts, args):
    """위상 배열에서 모노드로미 프로파일 (Δarg) 계산"""
    monos = np.zeros(len(ts))
    for i in range(1, len(ts)):
        delta = args[i] - args[i-1]
        while delta > np.pi:
            delta -= 2 * np.pi
        while delta < -np.pi:
            delta += 2 * np.pi
        monos[i] = delta
    return monos


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 위상 점프 수 세기
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def count_phase_jumps(sigma, t_min=10.0, t_max=50.0, n_points=2000,
                      func='zeta', char_info=None):
    """σ 고정, t를 스캔하며 |Δarg| > π/2인 점프 수 세기"""
    ts = np.linspace(t_min, t_max, n_points)
    jumps = 0
    prev_arg = None

    for t in ts:
        s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
        if func == 'zeta':
            val = xi_func(s)
        else:
            val = completed_L(s, char_info)

        if _near_zero(val):
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
# 에너지 프로파일
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def energy_profile(sigmas, t_min=13.0, t_max=34.0, n_t=200,
                   func='zeta', char_info=None):
    """E(σ) = ∫|L(σ+it)|² dt (사다리꼴 적분)"""
    ts = np.linspace(t_min, t_max, n_t)
    dt = (t_max - t_min) / (n_t - 1)
    energies = []

    for sigma in sigmas:
        vals = []
        for t in ts:
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            try:
                if func == 'zeta':
                    kappa = curvature_zeta(s)
                else:
                    kappa = curvature_dirichlet(s, char_info)
                vals.append(min(kappa, 1e6) if np.isfinite(kappa) else 1e6)
            except Exception:
                vals.append(0.0)
        E = np.trapezoid(vals, dx=dt)  # numpy 2.0
        energies.append(E)

    return np.array(energies)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 예측 평가
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def evaluate_predictions(pred_t, true_t, tolerance=0.5):
    """precision, recall, F1 계산"""
    if len(pred_t) == 0 or len(true_t) == 0:
        return 0.0, 0.0, 0.0

    tp = 0
    matched = set()
    for pt in pred_t:
        dists = np.abs(true_t - pt)
        best = np.argmin(dists)
        if dists[best] < tolerance and best not in matched:
            tp += 1
            matched.add(best)

    precision = tp / len(pred_t) if len(pred_t) > 0 else 0
    recall = tp / len(true_t) if len(true_t) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return precision, recall, f1
