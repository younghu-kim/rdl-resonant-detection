"""
=============================================================================
[C-351] GL(2) 다발 성질 검증 — 타원곡선 11a1
=============================================================================
목적: GL(1) Dirichlet에서 확립된 4가지 다발 성질이 GL(2) (타원곡선)에도 성립하는지.
대상: 11a1 (y² + y = x³ - x² - 10x - 20, conductor N=11, ε=+1)

방법: 불완전 감마 함수 기반 급속 수렴 공식 (lcalc 불필요)
  Λ(s) = Σ a_n [g(s, 2πn/√N) + ε g(2-s, 2πn/√N)]
  g(s, y) = y^{-s} Γ(s, y)

정규화: 클래식 (임계선 Re(s)=1, 함수방정식 Λ(s)=εΛ(2-s))
4성질: κδ², mono, detect, E비 (σ=1.0 vs σ=0.7 비교)
=============================================================================
"""

import numpy as np
import mpmath
import time
import sys

mpmath.mp.dps = 80

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 11a1 Dirichlet 계수 (Legendre 기호 O(p) 알고리즘)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_ap(p, N=11):
    """11a1: y² + y = x³ - x² - 10x - 20 의 a_p (Legendre 기호 O(p))"""
    if p == N:
        return 1
    if p == 2:
        count = sum(1 for x in range(2) for y in range(2)
                    if (y*y + y - x*x*x + x*x + 10*x + 20) % 2 == 0)
        return 2 + 1 - (count + 1)
    ap = 0
    neg10 = (-10) % p
    neg20 = (-20) % p
    for x in range(p):
        rhs = (pow(x, 3, p) - pow(x, 2, p) + neg10 * x + neg20) % p
        disc = (4 * rhs + 1) % p
        if disc == 0:
            pass
        else:
            leg = pow(disc, (p - 1) // 2, p)
            if leg == 1:
                ap -= 1
            else:
                ap += 1
    return ap


def compute_an_table(max_n, N=11):
    """a_n 테이블 (1-indexed, a[0]=0)"""
    a = [0] * (max_n + 1)
    a[1] = 1

    is_prime = [True] * (max_n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(max_n**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, max_n+1, i):
                is_prime[j] = False

    primes = [p for p in range(2, max_n+1) if is_prime[p]]
    print(f"a_p 계산 중 (소수 {len(primes)}개)...", flush=True)

    ap_cache = {}
    for p in primes:
        ap_cache[p] = compute_ap(p, N)
    print(f"  a_2={ap_cache[2]}, a_3={ap_cache[3]}, a_5={ap_cache[5]}, a_7={ap_cache[7]}", flush=True)

    # 소수 거듭제곱
    for p in primes:
        a[p] = ap_cache[p]
        pk = p
        k = 1
        while pk * p <= max_n:
            pk *= p
            k += 1
            if p == N:
                a[pk] = ap_cache[p] ** k
            else:
                a[pk] = ap_cache[p] * a[pk // p] - p * a[pk // (p*p)]

    # 합성수 (곱셈적)
    for n in range(2, max_n + 1):
        if a[n] != 0 or is_prime[n]:
            continue
        temp = n
        result = 1
        for p in primes:
            if p * p > temp:
                break
            if temp % p == 0:
                pe = 1
                while temp % p == 0:
                    pe *= p
                    temp //= p
                result *= a[pe]
        if temp > 1:
            result *= a[temp]
        a[n] = result

    return a


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# GL(2) 완비 L-함수 (불완전 감마 급속 수렴)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

class GL2CompletedL:
    """
    Λ(s) = Σ a_n [g(s, y_n) + ε g(2-s, y_n)]
    g(s, y) = y^{-s} Γ(s, y)   (upper incomplete gamma)
    y_n = 2πn / √N

    임계선: Re(s) = 1 (클래식 정규화)
    함수방정식: Λ(s) = ε Λ(2-s)
    """

    def __init__(self, an_table, N=11, epsilon=1):
        self.an = an_table
        self.N = N
        self.epsilon = epsilon
        self.sqrtN = mpmath.sqrt(mpmath.mpf(N))
        self.two_pi = 2 * mpmath.pi
        self.max_n = len(an_table) - 1

    def _g(self, s, y):
        """g(s, y) = y^{-s} Γ(s, y)"""
        return mpmath.power(y, -s) * mpmath.gammainc(s, y)

    def Lambda(self, s, n_terms=None):
        """완비 L-함수 Λ(s)"""
        if n_terms is None:
            # 자동: |s|에 비례하여 항 수 결정
            n_terms = max(50, int(abs(s) * float(self.sqrtN) / float(self.two_pi)) + 20)
        n_terms = min(n_terms, self.max_n)

        result = mpmath.mpf(0)
        s2 = 2 - s  # 함수방정식 대칭점
        for n in range(1, n_terms + 1):
            if self.an[n] == 0:
                continue
            y = self.two_pi * n / self.sqrtN
            term = mpmath.mpf(self.an[n]) * (self._g(s, y) + self.epsilon * self._g(s2, y))
            result += term
        return result

    def connection(self, s, n_terms=None):
        """접속 Λ'/Λ (수치 미분, h=10^{-20})"""
        h = mpmath.mpf(1) / mpmath.mpf(10**20)
        val = self.Lambda(s, n_terms)
        if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
            return mpmath.mpc(1e10, 0)
        val_d = (self.Lambda(s + h, n_terms) - self.Lambda(s - h, n_terms)) / (2 * h)
        return val_d / val

    def curvature(self, s, n_terms=None):
        """곡률 κ = |Λ'/Λ|²"""
        L = self.connection(s, n_terms)
        return float(abs(L)**2)

    def monodromy_phase_jump(self, t, eps=0.005, n_terms=None):
        """
        임계선 위상 점프 Δarg(Λ) = arg(Λ(σ_c+i(t+ε))) - arg(Λ(σ_c+i(t-ε)))
        단순 영점에서 ±π (논문/디리클레 스크립트와 동일 정의)
        σ_c = 1 (GL(2) 클래식 임계선)
        """
        s_plus = mpmath.mpf(1) + 1j * mpmath.mpf(str(t + eps))
        s_minus = mpmath.mpf(1) + 1j * mpmath.mpf(str(t - eps))
        val_plus = self.Lambda(s_plus, n_terms)
        val_minus = self.Lambda(s_minus, n_terms)
        arg_plus = float(mpmath.arg(val_plus))
        arg_minus = float(mpmath.arg(val_minus))
        delta = arg_plus - arg_minus
        import math
        while delta > math.pi:
            delta -= 2 * math.pi
        while delta < -math.pi:
            delta += 2 * math.pi
        return delta


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 영점 탐색 (Re(Λ)/Im(Λ) 부호 변화 + findroot)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def find_zeros_gl2(gl2, t_min=2.0, t_max=100.0, n_scan=4000):
    """임계선 Re(s)=1 위에서 영점 탐색"""
    ts = np.linspace(t_min, t_max, n_scan)
    zeros = []

    # Re(Λ) 부호 변화
    prev_re, prev_t = None, None
    print(f"  영점 스캔 t∈[{t_min},{t_max}] ({n_scan}점)...", flush=True)
    for i, t in enumerate(ts):
        s = mpmath.mpf(1) + 1j * mpmath.mpf(str(t))
        val = gl2.Lambda(s)
        curr_re = float(mpmath.re(val))

        if prev_re is not None and prev_re * curr_re < 0:
            try:
                def f_real(t_var):
                    sv = mpmath.mpf(1) + 1j * mpmath.mpf(t_var)
                    return mpmath.re(gl2.Lambda(sv))
                t_zero = float(mpmath.findroot(f_real, (prev_t, t)))
                # |Λ| 검증
                sv = mpmath.mpf(1) + 1j * mpmath.mpf(str(t_zero))
                if abs(gl2.Lambda(sv)) < mpmath.mpf('1e-5'):
                    if not zeros or abs(t_zero - zeros[-1]) > 0.1:
                        zeros.append(t_zero)
                        print(f"    영점 #{len(zeros)}: t={t_zero:.6f}", flush=True)
            except Exception:
                pass
        prev_re, prev_t = curr_re, t

        if (i+1) % 1000 == 0:
            print(f"    스캔 진행: {i+1}/{n_scan} (영점 {len(zeros)}개)", flush=True)

    # Im(Λ) 부호 변화 (추가 영점)
    prev_im, prev_t = None, None
    for t in ts:
        s = mpmath.mpf(1) + 1j * mpmath.mpf(str(t))
        val = gl2.Lambda(s)
        curr_im = float(mpmath.im(val))

        if prev_im is not None and prev_im * curr_im < 0:
            try:
                def f_imag(t_var):
                    sv = mpmath.mpf(1) + 1j * mpmath.mpf(t_var)
                    return mpmath.im(gl2.Lambda(sv))
                t_zero = float(mpmath.findroot(f_imag, (prev_t, t)))
                sv = mpmath.mpf(1) + 1j * mpmath.mpf(str(t_zero))
                if abs(gl2.Lambda(sv)) < mpmath.mpf('1e-5'):
                    if not any(abs(t_zero - z) < 0.1 for z in zeros):
                        zeros.append(t_zero)
                        print(f"    영점 (Im) #{len(zeros)}: t={t_zero:.6f}", flush=True)
            except Exception:
                pass
        prev_im, prev_t = curr_im, t

    zeros.sort()
    return np.array(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4성질 검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def verify_kd2(gl2, zeros, delta=0.03):
    """[성질 1] κδ² ≈ 1"""
    results = []
    for t_z in zeros:
        s = mpmath.mpf(1) + mpmath.mpf(str(delta)) + 1j * mpmath.mpf(str(t_z))
        kappa = gl2.curvature(s)
        kd2 = kappa * delta**2
        results.append(kd2)
    return np.array(results)


def verify_mono(gl2, zeros, eps=0.005):
    """[성질 2] 임계선 위상 점프 ≈ ±π"""
    results = []
    for t_z in zeros:
        m = gl2.monodromy_phase_jump(t_z, eps=eps)
        results.append(m)
    return np.array(results)


def verify_detect(gl2, zeros, t_min, t_max, n_points=3000):
    """[성질 3] 곡률 피크 → 영점 검출"""
    ts = np.linspace(t_min, t_max, n_points)
    kappas = np.zeros(n_points)
    print(f"    곡률 프로파일 계산 ({n_points}점)...", flush=True)
    for i, t in enumerate(ts):
        s = mpmath.mpf(1) + 1j * mpmath.mpf(str(t))
        kappas[i] = gl2.curvature(s)
        if (i+1) % 500 == 0:
            print(f"      {i+1}/{n_points}", flush=True)

    # 극대점 탐색
    peaks = []
    for i in range(1, len(kappas) - 1):
        if kappas[i] > kappas[i-1] and kappas[i] > kappas[i+1]:
            peaks.append(ts[i])

    # 매칭
    tol = 1.0
    eligible = [z for z in zeros if t_min + tol < z < t_max - tol]
    detected = sum(1 for z in eligible if len(peaks) > 0 and
                   min(abs(np.array(peaks) - z)) < tol)
    recall = detected / len(eligible) if eligible else 0
    return recall, detected, len(eligible)


def verify_eratio(gl2, zeros, n_per_zero=10):
    """[성질 4] E(σ=1.0)/E(σ=0.7) 에너지 비율"""
    results = []
    for t_z in zeros[:n_per_zero]:
        t_pts = np.linspace(t_z - 2.0, t_z + 2.0, 21)
        E_crit = 0.0  # σ=1.0
        E_off = 0.0   # σ=0.7
        for t in t_pts:
            s_crit = mpmath.mpf(1) + 1j * mpmath.mpf(str(t))
            s_off = mpmath.mpf('0.7') + 1j * mpmath.mpf(str(t))
            E_crit += gl2.curvature(s_crit)
            E_off += gl2.curvature(s_off)
        if E_off > 0:
            results.append(E_crit / E_off)
    return np.array(results)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    T_START = time.time()
    RESULT_FILE = "results/gl2_bundle_verification_c351.txt"

    print("=" * 70, flush=True)
    print("[C-351] GL(2) 다발 성질 검증 — 11a1", flush=True)
    print(f"  클래식 정규화: 임계선 Re(s)=1", flush=True)
    print(f"  불완전 감마 급속 수렴 공식", flush=True)
    print("=" * 70, flush=True)

    # ── 1. a_n 계수 ──
    print("\n[1/5] a_n 계수 계산...", flush=True)
    t0 = time.time()
    NUM_COEFFS = 500  # 불완전 감마는 적은 항으로 수렴
    an = compute_an_table(NUM_COEFFS)
    print(f"  {time.time()-t0:.1f}초", flush=True)

    # LMFDB 참조값 검증
    ref = {2: -2, 3: -1, 5: 1, 7: -2, 11: 1, 13: 4, 17: -2, 19: 0, 23: -1, 29: 0, 31: 7, 37: 3}
    all_ok = all(an[p] == ref[p] for p in ref)
    print(f"  a_p 검증: {'✅ ALL MATCH' if all_ok else '❌ MISMATCH'}", flush=True)
    if not all_ok:
        for p in ref:
            if an[p] != ref[p]:
                print(f"    a_{p}: got {an[p]}, expected {ref[p]}", flush=True)
        return

    # ── 2. GL(2) L-함수 설정 ──
    print("\n[2/5] GL(2) L-함수 설정...", flush=True)
    gl2 = GL2CompletedL(an, N=11, epsilon=1)

    # 함수방정식 검증
    s_test = mpmath.mpf('1.3') + 5j
    v1 = gl2.Lambda(s_test)
    v2 = gl2.Lambda(2 - s_test)
    fe_err = float(abs(v1 - v2) / abs(v1))
    print(f"  함수방정식 Λ(s)=εΛ(2-s) 오차: {fe_err:.2e}", flush=True)
    assert fe_err < 1e-10, f"함수방정식 불만족: {fe_err}"
    print(f"  ✅ 함수방정식 검증 통과", flush=True)

    # ── 3. 영점 탐색 ──
    print("\n[3/5] 영점 탐색 (t∈[2,100])...", flush=True)
    t0 = time.time()
    zeros = find_zeros_gl2(gl2, t_min=2.0, t_max=100.0, n_scan=3000)
    print(f"  영점 {len(zeros)}개 발견 ({time.time()-t0:.1f}초)", flush=True)
    if len(zeros) > 0:
        print(f"  첫 5개: {[f'{z:.4f}' for z in zeros[:5]]}", flush=True)
        # LMFDB 첫 영점 6.3622와 비교
        print(f"  첫 영점: {zeros[0]:.6f} (LMFDB: 6.362190)", flush=True)

    if len(zeros) < 5:
        print("❌ 영점 부족. 중단.", flush=True)
        return

    # ── 4. 4성질 검증 ──
    n_test = min(20, len(zeros))
    test_zeros = zeros[:n_test]
    print(f"\n[4/5] 4성질 검증 ({n_test}개 영점)...", flush=True)

    # 4-1. κδ²
    print("\n  [4-1] κδ² 곡률 집중 (δ=0.03)...", flush=True)
    t0 = time.time()
    kd2 = verify_kd2(gl2, test_zeros, delta=0.03)
    kd2_med = float(np.median(kd2))
    print(f"    중앙값: {kd2_med:.4f} (기대 ~1.0)", flush=True)
    print(f"    범위: [{np.min(kd2):.4f}, {np.max(kd2):.4f}]", flush=True)
    print(f"    {time.time()-t0:.1f}초", flush=True)

    # 4-2. 모노드로미
    print("\n  [4-2] 모노드로미 (r=0.5)...", flush=True)
    t0 = time.time()
    mono = verify_mono(gl2, test_zeros, radius=0.5)
    mono_abs = np.abs(mono)
    mono_dev = np.abs(mono_abs / np.pi - 1.0)
    print(f"    |mono|/π 평균: {np.mean(mono_abs/np.pi):.6f}", flush=True)
    print(f"    편차 최대: {np.max(mono_dev):.6f}", flush=True)
    n_pass_mono = int(np.sum(mono_dev < 0.01))
    print(f"    PASS: {n_pass_mono}/{n_test}", flush=True)
    print(f"    {time.time()-t0:.1f}초", flush=True)

    # 4-3. 검출률 (t∈[5,50])
    print("\n  [4-3] 곡률 검출...", flush=True)
    t0 = time.time()
    recall, detected, eligible = verify_detect(gl2, zeros, t_min=5.0, t_max=50.0, n_profile=2000)
    print(f"    recall: {recall:.1%} ({detected}/{eligible})", flush=True)
    print(f"    {time.time()-t0:.1f}초", flush=True)

    # 4-4. E비
    print("\n  [4-4] E(1.0)/E(0.7)...", flush=True)
    t0 = time.time()
    er = verify_eratio(gl2, test_zeros, n_per_zero=min(10, n_test))
    print(f"    중앙값: {np.median(er):.1f}×", flush=True)
    print(f"    범위: [{np.min(er):.1f}, {np.max(er):.1f}]×", flush=True)
    print(f"    {time.time()-t0:.1f}초", flush=True)

    # ── 5. 결과 파일 ──
    elapsed = time.time() - T_START

    kd2_pass = 0.8 < kd2_med < 1.2
    mono_pass = n_pass_mono / n_test >= 0.9
    detect_pass = recall >= 0.8
    eratio_pass = float(np.median(er)) > 2.0

    all_pass = kd2_pass and mono_pass and detect_pass and eratio_pass
    verdict = "ALL-PASS" if all_pass else "PARTIAL"

    with open(RESULT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("[C-351] GL(2) 다발 성질 검증 — 11a1\n")
        f.write("=" * 70 + "\n\n")

        f.write("설정:\n")
        f.write(f"  대상: 11a1 (y² + y = x³ - x² - 10x - 20)\n")
        f.write(f"  N=11, ε=+1, rank=0\n")
        f.write(f"  정규화: 클래식 (임계선 σ=1)\n")
        f.write(f"  계수: {NUM_COEFFS}개, dps={mpmath.mp.dps}\n")
        f.write(f"  영점: {len(zeros)}개 (t∈[2,100])\n")
        f.write(f"  검증 영점: {n_test}개\n")
        f.write(f"  함수방정식 오차: {fe_err:.2e}\n\n")

        f.write("a_p 검증:\n")
        for p in sorted(ref.keys()):
            f.write(f"  a_{p} = {an[p]} ✅\n")
        f.write("\n")

        f.write(f"영점 (처음 {min(10, len(zeros))}개):\n")
        for i, z in enumerate(zeros[:10]):
            f.write(f"  {i+1:3d}. t = {z:.6f}\n")
        f.write(f"  ... (총 {len(zeros)}개)\n\n")

        f.write("-" * 70 + "\n")
        f.write("4성질 검증 결과\n")
        f.write("-" * 70 + "\n\n")

        f.write(f"[1] κδ² (δ=0.03):\n")
        f.write(f"    중앙값: {kd2_med:.4f}\n")
        f.write(f"    범위: [{np.min(kd2):.4f}, {np.max(kd2):.4f}]\n")
        f.write(f"    판정: {'PASS ✅' if kd2_pass else 'FAIL ❌'}\n\n")

        f.write(f"[2] 모노드로미 (r=0.5):\n")
        f.write(f"    |mono|/π 평균: {np.mean(mono_abs/np.pi):.6f}\n")
        f.write(f"    편차 최대: {np.max(mono_dev):.6f}\n")
        f.write(f"    PASS: {n_pass_mono}/{n_test}\n")
        f.write(f"    판정: {'PASS ✅' if mono_pass else 'FAIL ❌'}\n\n")

        f.write(f"[3] 곡률 검출 (t∈[5,50]):\n")
        f.write(f"    recall: {recall:.1%} ({detected}/{eligible})\n")
        f.write(f"    판정: {'PASS ✅' if detect_pass else 'FAIL ❌'}\n\n")

        f.write(f"[4] E(1.0)/E(0.7):\n")
        f.write(f"    중앙값: {np.median(er):.1f}×\n")
        f.write(f"    범위: [{np.min(er):.1f}, {np.max(er):.1f}]×\n")
        f.write(f"    판정: {'PASS ✅' if eratio_pass else 'FAIL ❌'}\n\n")

        f.write("=" * 70 + "\n")
        f.write(f"종합: {verdict}\n")
        f.write(f"시간: {elapsed:.1f}초 ({elapsed/60:.1f}분)\n")
        f.write("=" * 70 + "\n")

        # 상세 데이터
        f.write("\n\n상세 데이터\n" + "-" * 40 + "\n\n")
        f.write("κδ² per zero:\n")
        for i, (z, v) in enumerate(zip(test_zeros, kd2)):
            f.write(f"  {i+1:3d}. t={z:.4f}  κδ²={v:.4f}\n")

        f.write("\nmono per zero:\n")
        for i, (z, m) in enumerate(zip(test_zeros, mono)):
            f.write(f"  {i+1:3d}. t={z:.4f}  mono={m/np.pi:.6f}π  dev={abs(abs(m)/np.pi-1):.6f}\n")

        f.write("\nE비 per zero:\n")
        for i, (z, e) in enumerate(zip(test_zeros[:len(er)], er)):
            f.write(f"  {i+1:3d}. t={z:.4f}  E비={e:.1f}×\n")

    print(f"\n{'='*70}", flush=True)
    print(f"결과: {RESULT_FILE}", flush=True)
    print(f"판정: {verdict}", flush=True)
    print(f"  κδ²={kd2_med:.4f} {'✅' if kd2_pass else '❌'}", flush=True)
    print(f"  mono PASS={n_pass_mono}/{n_test} {'✅' if mono_pass else '❌'}", flush=True)
    print(f"  detect={recall:.0%} {'✅' if detect_pass else '❌'}", flush=True)
    print(f"  E비={np.median(er):.1f}× {'✅' if eratio_pass else '❌'}", flush=True)
    print(f"시간: {elapsed/60:.1f}분", flush=True)


if __name__ == '__main__':
    main()
