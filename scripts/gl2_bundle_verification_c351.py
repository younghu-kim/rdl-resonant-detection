"""
=============================================================================
[C-351] GL(2) 다발 성질 검증 — 타원곡선 11a1
=============================================================================
목적: GL(1) Dirichlet에서 확립된 4가지 다발 성질이 GL(2) (타원곡선)에도 성립하는지 검증.
대상: 11a1 (y² + y = x³ - x² - 10x - 20, conductor N=11, root number ε=+1)

4가지 검증 성질:
  1. κδ² ≈ 1 (곡률 집중)
  2. 모노드로미 = ±π (단순 영점에서)
  3. detect = 100% (곡률 피크에서 영점 검출)
  4. E(0.5)/E(0.3) >> 1 (에너지 비율)

비교 대상: ζ(s) 결과 (기존 확립)

참조:
  - bundle_utils.py (제타/디리클레 공용 함수)
  - dirichlet_chi11_c328.py (GL(1) 검증 패턴)
  - checklist: 80dps, t∈[10,200], δ=0.03, python -u
=============================================================================
"""

import numpy as np
import mpmath
import time
import sys

mpmath.mp.dps = 80  # 고정밀도 (t>100 필수)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 11a1 Dirichlet 계수 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_ap(p, N=11):
    """
    a_p 계산 for 11a1: y² + y = x³ - x² - 10x - 20.
    Legendre 기호 사용 O(p) 알고리즘.
    y² + y = rhs ⟺ (2y+1)² = 4*rhs + 1 mod p
    a_p = -Σ_{x=0}^{p-1} (4*rhs(x) + 1 | p)
    """
    if p == N:
        return 1  # split multiplicative reduction at p=11
    if p == 2:
        # 특수 처리: p=2에서 직접 점 계수
        count = 0
        for x in range(2):
            for y in range(2):
                if (y*y + y - x*x*x + x*x + 10*x + 20) % 2 == 0:
                    count += 1
        return 2 + 1 - (count + 1)

    # Legendre 기호 활용: O(p)
    ap = 0
    for x in range(p):
        rhs = (pow(x, 3, p) - pow(x, 2, p) + (-10 % p) * x + (-20 % p)) % p
        disc = (4 * rhs + 1) % p
        # Legendre symbol via Euler's criterion
        if disc == 0:
            pass  # contributes 0
        else:
            leg = pow(disc, (p - 1) // 2, p)
            if leg == 1:
                ap -= 1  # -(+1)
            else:  # leg == p-1
                ap += 1  # -(-1)
    return ap


def compute_an_table(max_n):
    """a_n 테이블 계산 (1-indexed). 반환: list[0..max_n]"""
    N = 11
    a = [0] * (max_n + 1)
    a[1] = 1

    # Sieve of Eratosthenes for primes
    is_prime = [True] * (max_n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(max_n**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, max_n+1, i):
                is_prime[j] = False

    primes = [p for p in range(2, max_n+1) if is_prime[p]]
    ap_cache = {}

    # Compute a_p
    print(f"a_p 계산 중 (소수 {len(primes)}개)...", flush=True)
    for p in primes:
        ap_cache[p] = compute_ap(p, N)
    print(f"  a_2={ap_cache[2]}, a_3={ap_cache[3]}, a_5={ap_cache[5]}, a_7={ap_cache[7]}", flush=True)

    # Fill a_n using multiplicativity
    for p in primes:
        a[p] = ap_cache[p]
        # Powers of p
        pk = p
        k = 1
        while pk * p <= max_n:
            pk *= p
            k += 1
            if p == N:
                # Bad prime: a_{p^k} = a_p^k
                a[pk] = ap_cache[p] ** k
            else:
                # Good prime: a_{p^k} = a_p * a_{p^{k-1}} - p * a_{p^{k-2}}
                a[pk] = ap_cache[p] * a[pk // p] - p * a[pk // (p*p)] if k >= 2 else ap_cache[p]

    # Multiplicativity for composite n
    for n in range(2, max_n + 1):
        if a[n] != 0 or is_prime[n]:
            continue
        # Factor n and use multiplicativity
        temp = n
        factors = {}
        for p in primes:
            if p * p > temp:
                break
            while temp % p == 0:
                factors[p] = factors.get(p, 0) + 1
                temp //= p
        if temp > 1:
            factors[temp] = factors.get(temp, 0) + 1

        result = 1
        for p, e in factors.items():
            pe = p ** e
            result *= a[pe]
        a[n] = result

    return a


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# lcalc 기반 GL(2) L-함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def setup_lcalc_gl2(an_table, N=11, epsilon=1, num_coeffs=10000):
    """lcalc L-function 객체 생성 (analytic normalization)"""
    from sage.libs.lcalc.lcalc_Lfunction import Lfunction_D

    Q = float(np.sqrt(N) / (2 * np.pi))  # sqrt(N)/(2π)
    OMEGA = float(epsilon)  # root number
    kappa = [0.5]  # Γ(s/2 + 1/2)
    lambd = [0.5]

    # Dirichlet coefficients: a_n / sqrt(n) (analytic normalization)
    dirichlet_coeffs = [float(an_table[n]) / np.sqrt(n) for n in range(1, num_coeffs + 1)]

    L = Lfunction_D("11a1", 0, dirichlet_coeffs, 0, Q, OMEGA, kappa, lambd, [], [])
    return L


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 완비 L-함수 (mpmath 고정밀도)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

class GL2LFunction:
    """GL(2) 완비 L-함수: mpmath 기반 고정밀도 계산"""

    def __init__(self, an_table, N=11, epsilon=1):
        self.an = an_table
        self.N = mpmath.mpf(N)
        self.epsilon = epsilon
        self.max_n = len(an_table) - 1
        # Precompute: Q = sqrt(N)/(2π)
        self.Q = mpmath.sqrt(self.N) / (2 * mpmath.pi)

    def completed_L(self, s):
        """
        Λ(s) = Q^s Γ(s/2 + 1/2) L*(s)
        where L*(s) = Σ (a_n/√n) n^{-s}
        Analytic normalization: critical line at Re(s) = 1/2
        """
        # L*(s) via Dirichlet series
        L_val = mpmath.mpf(0)
        for n in range(1, min(self.max_n + 1, 20001)):
            if self.an[n] == 0:
                continue
            term = mpmath.mpf(self.an[n]) / mpmath.power(n, s + mpmath.mpf('0.5'))
            L_val += term
            # 조기 종료: 항이 충분히 작으면
            if n > 100 and abs(term) < mpmath.mpf(10)**(-mpmath.mp.dps + 5):
                break

        gamma_val = mpmath.gamma(s / 2 + mpmath.mpf('0.5'))
        Q_s = mpmath.power(self.Q, s)
        return Q_s * gamma_val * L_val

    def connection(self, s):
        """접속 Λ'/Λ (수치 미분)"""
        h = mpmath.mpf(1) / mpmath.mpf(10**20)
        val = self.completed_L(s)
        if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
            return mpmath.mpc(1e10, 0)
        val_d = (self.completed_L(s + h) - self.completed_L(s - h)) / (2 * h)
        return val_d / val

    def curvature(self, s):
        """곡률 κ = |Λ'/Λ|²"""
        L = self.connection(s)
        return float(abs(L)**2)

    def monodromy(self, t, radius=0.5, n_steps=64):
        """영점 주위 폐곡선 적분으로 모노드로미 계산"""
        s_center = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        total_delta = 0.0
        prev_arg = None

        for k in range(n_steps + 1):
            theta = 2 * np.pi * k / n_steps
            s = s_center + radius * mpmath.exp(1j * theta)
            val = self.completed_L(s)

            if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
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

def find_zeros_gl2_lcalc(L_lcalc, n_zeros=200):
    """lcalc로 GL(2) 영점 탐색"""
    zeros = L_lcalc.find_zeros_via_N(n_zeros)
    return [float(z) for z in zeros]


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4-성질 검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def verify_curvature_concentration(gl2_func, zeros, delta=0.03):
    """성질 1: κδ² ≈ 1 at zeros"""
    kd2_values = []
    for t_zero in zeros:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero))
        # δ만큼 떨어진 점에서 곡률 측정
        s_near = s + mpmath.mpf(str(delta))
        kappa = gl2_func.curvature(s_near)
        kd2 = kappa * delta**2
        kd2_values.append(kd2)
    return np.array(kd2_values)


def verify_monodromy(gl2_func, zeros, radius=0.5):
    """성질 2: 모노드로미 ≈ ±π"""
    mono_values = []
    for t_zero in zeros:
        mono = gl2_func.monodromy(t_zero, radius=radius)
        mono_values.append(mono)
    return np.array(mono_values)


def verify_detection(gl2_func, zeros, t_min, t_max, n_profile=5000):
    """성질 3: 곡률 피크로 영점 검출"""
    ts = np.linspace(t_min, t_max, n_profile)
    kappas = []
    for t in ts:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        k = gl2_func.curvature(s)
        kappas.append(k)
    kappas = np.array(kappas)

    # 곡률 피크 탐색 (로컬 극대)
    peaks = []
    for i in range(1, len(kappas) - 1):
        if kappas[i] > kappas[i-1] and kappas[i] > kappas[i+1]:
            peaks.append(ts[i])
    peaks = np.array(peaks)

    # 영점과 매칭 (tolerance = 1.0)
    tol = 1.0
    detected = 0
    for z in zeros:
        if t_min + tol < z < t_max - tol:
            if len(peaks) > 0 and np.min(np.abs(peaks - z)) < tol:
                detected += 1
    eligible = sum(1 for z in zeros if t_min + tol < z < t_max - tol)
    recall = detected / eligible if eligible > 0 else 0
    return recall, detected, eligible


def verify_energy_ratio(gl2_func, zeros, t_range=2.0, n_sigma=101):
    """성질 4: E(0.5)/E(0.3) >> 1"""
    sigma_05 = np.linspace(0.5 - 0.3, 0.5 + 0.3, n_sigma)
    sigma_03 = np.linspace(0.3 - 0.3, 0.3 + 0.3, n_sigma)

    E_ratios = []
    for t_zero in zeros[:20]:  # 처음 20개만 (속도)
        # E(σ₀) = Σ κ(σ₀ + it) for t near t_zero
        t_range_pts = np.linspace(t_zero - t_range, t_zero + t_range, 21)

        E_05 = 0.0
        for t in t_range_pts:
            s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
            E_05 += gl2_func.curvature(s)

        E_03 = 0.0
        for t in t_range_pts:
            s = mpmath.mpf('0.3') + 1j * mpmath.mpf(str(t))
            E_03 += gl2_func.curvature(s)

        if E_03 > 0:
            E_ratios.append(E_05 / E_03)
    return np.array(E_ratios)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    T_START = time.time()
    RESULT_FILE = "results/gl2_bundle_verification_c351.txt"

    print("=" * 70, flush=True)
    print("[C-351] GL(2) 다발 성질 검증 — 11a1", flush=True)
    print("=" * 70, flush=True)

    # 1. a_n 계산
    print("\n[1/5] a_n 계수 계산...", flush=True)
    t0 = time.time()
    NUM_COEFFS = 15000
    an_table = compute_an_table(NUM_COEFFS)
    print(f"  완료: {time.time()-t0:.1f}초", flush=True)
    print(f"  검증: a_2={an_table[2]}, a_3={an_table[3]}, a_5={an_table[5]}, "
          f"a_7={an_table[7]}, a_11={an_table[11]}, a_13={an_table[13]}", flush=True)

    # 검증: LMFDB 참조값
    expected_ap = {2: -2, 3: -1, 5: 1, 7: -2, 11: 1, 13: 4, 17: -2, 19: 0, 23: -1}
    all_match = True
    for p, expected in expected_ap.items():
        if an_table[p] != expected:
            print(f"  ⚠️ a_{p} 불일치: {an_table[p]} ≠ {expected}", flush=True)
            all_match = False
    if all_match:
        print("  ✅ 모든 a_p 참조값과 일치", flush=True)
    else:
        print("  ❌ a_p 불일치 — 중단", flush=True)
        return

    # 2. lcalc L-function 설정 + 영점 탐색
    print("\n[2/5] lcalc 설정 + 영점 탐색...", flush=True)
    t0 = time.time()
    L_lcalc = setup_lcalc_gl2(an_table, N=11, epsilon=1, num_coeffs=NUM_COEFFS)

    # lcalc로 영점 탐색
    zeros_lcalc = find_zeros_gl2_lcalc(L_lcalc, n_zeros=200)
    # t∈[10,200] 필터
    zeros = [z for z in zeros_lcalc if 10.0 <= z <= 200.0]
    print(f"  총 영점: {len(zeros_lcalc)}, t∈[10,200]: {len(zeros)}", flush=True)
    if zeros:
        print(f"  첫 5개: {zeros[:5]}", flush=True)
    print(f"  완료: {time.time()-t0:.1f}초", flush=True)

    if len(zeros) < 5:
        print("영점 부족 — 계수 수 증가 필요. 중단.", flush=True)
        return

    # 3. mpmath GL(2) 함수 설정
    print("\n[3/5] mpmath GL(2) L-function 설정...", flush=True)
    gl2 = GL2LFunction(an_table, N=11, epsilon=1)

    # Λ(s) at first zero validation
    z0 = zeros[0]
    s0 = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(z0))
    val_at_zero = gl2.completed_L(s0)
    print(f"  Λ(0.5+i*{z0:.4f}) = {float(abs(val_at_zero)):.2e} (should be ~0)", flush=True)

    # 4. 4성질 검증
    print("\n[4/5] 4성질 검증...", flush=True)

    # 4-1. κδ²
    print("  [4-1] κδ² 곡률 집중...", flush=True)
    t0 = time.time()
    # 처음 30개 영점으로 시작
    test_zeros = zeros[:30]
    kd2_vals = verify_curvature_concentration(gl2, test_zeros, delta=0.03)
    kd2_med = np.median(kd2_vals)
    print(f"    κδ² 중앙값: {kd2_med:.4f} (기대: ~1.0)", flush=True)
    print(f"    범위: [{np.min(kd2_vals):.4f}, {np.max(kd2_vals):.4f}]", flush=True)
    print(f"    {time.time()-t0:.1f}초", flush=True)

    # 4-2. 모노드로미
    print("  [4-2] 모노드로미...", flush=True)
    t0 = time.time()
    mono_vals = verify_monodromy(gl2, test_zeros, radius=0.5)
    mono_abs = np.abs(mono_vals)
    mono_dev = np.abs(mono_abs / np.pi - 1.0)
    print(f"    |mono|/π 평균: {np.mean(mono_abs/np.pi):.6f}", flush=True)
    print(f"    편차 최대: {np.max(mono_dev):.6f}", flush=True)
    n_pass_mono = np.sum(mono_dev < 0.01)
    print(f"    PASS (dev<0.01): {n_pass_mono}/{len(test_zeros)}", flush=True)
    print(f"    {time.time()-t0:.1f}초", flush=True)

    # 4-3. 검출률
    print("  [4-3] 곡률 기반 검출...", flush=True)
    t0 = time.time()
    recall, detected, eligible = verify_detection(gl2, zeros, t_min=10.0, t_max=60.0, n_profile=3000)
    print(f"    recall: {recall:.1%} ({detected}/{eligible})", flush=True)
    print(f"    {time.time()-t0:.1f}초", flush=True)

    # 4-4. 에너지 비율
    print("  [4-4] E(0.5)/E(0.3)...", flush=True)
    t0 = time.time()
    e_ratios = verify_energy_ratio(gl2, test_zeros[:10])
    print(f"    E비 중앙값: {np.median(e_ratios):.1f}×", flush=True)
    print(f"    E비 범위: [{np.min(e_ratios):.1f}, {np.max(e_ratios):.1f}]×", flush=True)
    print(f"    {time.time()-t0:.1f}초", flush=True)

    # 5. 결과 파일
    print("\n[5/5] 결과 저장...", flush=True)
    elapsed = time.time() - T_START

    # 판정
    kd2_pass = 0.8 < kd2_med < 1.2
    mono_pass = n_pass_mono / len(test_zeros) >= 0.9
    detect_pass = recall >= 0.8
    eratio_pass = np.median(e_ratios) > 2.0

    all_pass = kd2_pass and mono_pass and detect_pass and eratio_pass
    verdict = "ALL-PASS" if all_pass else "PARTIAL"

    with open(RESULT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("[C-351] GL(2) 다발 성질 검증 — 11a1\n")
        f.write("=" * 70 + "\n\n")

        f.write("설정:\n")
        f.write(f"  대상: 11a1 (y² + y = x³ - x² - 10x - 20)\n")
        f.write(f"  conductor: N=11, root number: ε=+1\n")
        f.write(f"  계수: {NUM_COEFFS}개, dps: {mpmath.mp.dps}\n")
        f.write(f"  영점: {len(zeros)}개 (t∈[10,200])\n")
        f.write(f"  검증 영점: {len(test_zeros)}개 (처음 30개)\n\n")

        f.write("a_p 검증:\n")
        for p in sorted(expected_ap.keys()):
            f.write(f"  a_{p} = {an_table[p]} (LMFDB: {expected_ap[p]}) ✅\n")
        f.write("\n")

        f.write("영점 (처음 10개):\n")
        for i, z in enumerate(zeros[:10]):
            f.write(f"  {i+1:3d}. t = {z:.6f}\n")
        f.write(f"  ... (총 {len(zeros)}개)\n\n")

        f.write("-" * 70 + "\n")
        f.write("4성질 검증 결과\n")
        f.write("-" * 70 + "\n\n")

        f.write(f"[1] κδ² 곡률 집중 (δ=0.03):\n")
        f.write(f"    중앙값: {kd2_med:.4f}\n")
        f.write(f"    범위: [{np.min(kd2_vals):.4f}, {np.max(kd2_vals):.4f}]\n")
        f.write(f"    판정: {'PASS' if kd2_pass else 'FAIL'} ({'✅' if kd2_pass else '❌'})\n\n")

        f.write(f"[2] 모노드로미 (r=0.5):\n")
        f.write(f"    |mono|/π 평균: {np.mean(mono_abs/np.pi):.6f}\n")
        f.write(f"    편차 최대: {np.max(mono_dev):.6f}\n")
        f.write(f"    PASS: {n_pass_mono}/{len(test_zeros)}\n")
        f.write(f"    판정: {'PASS' if mono_pass else 'FAIL'} ({'✅' if mono_pass else '❌'})\n\n")

        f.write(f"[3] 곡률 기반 검출 (t∈[10,60]):\n")
        f.write(f"    recall: {recall:.1%} ({detected}/{eligible})\n")
        f.write(f"    판정: {'PASS' if detect_pass else 'FAIL'} ({'✅' if detect_pass else '❌'})\n\n")

        f.write(f"[4] E(0.5)/E(0.3) 에너지 비율:\n")
        f.write(f"    중앙값: {np.median(e_ratios):.1f}×\n")
        f.write(f"    범위: [{np.min(e_ratios):.1f}, {np.max(e_ratios):.1f}]×\n")
        f.write(f"    판정: {'PASS' if eratio_pass else 'FAIL'} ({'✅' if eratio_pass else '❌'})\n\n")

        f.write("=" * 70 + "\n")
        f.write(f"종합 판정: {verdict}\n")
        f.write(f"소요 시간: {elapsed:.1f}초 ({elapsed/60:.1f}분)\n")
        f.write("=" * 70 + "\n\n")

        # 상세 데이터
        f.write("\n상세 데이터:\n\n")
        f.write("κδ² per zero:\n")
        for i, (z, kd2) in enumerate(zip(test_zeros, kd2_vals)):
            f.write(f"  {i+1:3d}. t={z:.4f}  κδ²={kd2:.4f}\n")

        f.write("\n모노드로미 per zero:\n")
        for i, (z, m) in enumerate(zip(test_zeros, mono_vals)):
            f.write(f"  {i+1:3d}. t={z:.4f}  mono={m/np.pi:.6f}π  dev={abs(abs(m)/np.pi - 1):.6f}\n")

        f.write("\n에너지 비율 per zero:\n")
        for i, (z, er) in enumerate(zip(test_zeros[:10], e_ratios)):
            f.write(f"  {i+1:3d}. t={z:.4f}  E비={er:.1f}×\n")

    print(f"\n결과 저장: {RESULT_FILE}", flush=True)
    print(f"종합 판정: {verdict}", flush=True)
    print(f"소요 시간: {elapsed:.1f}초", flush=True)


if __name__ == '__main__':
    main()
