"""
=============================================================================
[Project RDL] 고 conductor 디리클레 L-함수 다발 검증
=============================================================================
목표: mod 3,4,5를 넘어 mod 7, 8, 11로 보편성 확장
  - χ mod 7 (복소, order 6): 첫 번째 "큰" 복소 지표
  - χ mod 8 (실수): χ₄와 비교하여 실수 지표 에너지 이상치 검증
  - χ mod 11 (복소, order 10): 소수 conductor

각 지표에 대해 5성질 검증:
  1. 영점 수
  2. κ 집중도
  3. 모노드로미 양자화
  4. 에너지 집중도
  5. σ=0.5 위상 점프 유일성
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
import cmath

mpmath.mp.dps = 80

T_MIN, T_MAX = 10.0, 40.0

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 고 conductor 지표 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# mod 7: primitive root = 3, order 6 character
w6 = cmath.exp(2j * cmath.pi / 6)
_chi7 = [0, 1, w6**2, w6, w6**4, w6**5, w6**3]  # 3^k mod 7: 1→0,2→2,3→1,4→4,5→5,6→3

# mod 8: real character (Kronecker (-1/n))
_chi8 = [0, 1, 0, -1, 0, -1, 0, 1]

# mod 11: primitive root = 2, order 10 character
w10 = cmath.exp(2j * cmath.pi / 10)
_inds11 = {1:0, 2:1, 4:2, 8:3, 5:4, 10:5, 9:6, 7:7, 3:8, 6:9}
_chi11 = [0] + [w10**_inds11[k] for k in range(1, 11)]

CHARACTERS = {
    'χ₇ (mod 7)': {
        'chi': _chi7, 'q': 7, 'a': 1, 'type': 'complex',  # χ(-1)=w6³=-1 → odd → a=1
    },
    'χ₈ (mod 8)': {
        'chi': _chi8, 'q': 8, 'a': 0, 'type': 'real',  # χ(-1)=χ(7)=+1 → even → a=0
    },
    'χ₁₁ (mod 11)': {
        'chi': _chi11, 'q': 11, 'a': 1, 'type': 'complex',  # χ(-1)=w10⁵=-1 → odd → a=1
    },
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 함수 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def completed_L(s, char_info):
    """완비 L-함수: Λ(s, χ) = (q/π)^{s/2} Γ((s+a)/2) L(s, χ)"""
    q = mpmath.mpf(char_info['q'])
    a = mpmath.mpf(char_info['a'])
    L_val = mpmath.dirichlet(s, char_info['chi'])
    gamma_val = mpmath.gamma((s + a) / 2)
    prefactor = mpmath.power(q / mpmath.pi, s / 2)
    return prefactor * gamma_val * L_val


def connection(s, char_info):
    """접속 L = Λ'/Λ (수치 미분)"""
    h = mpmath.mpf(1) / mpmath.mpf(10**20)
    val = completed_L(s, char_info)
    if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e10, 0)
    d = (completed_L(s + h, char_info) - completed_L(s - h, char_info)) / (2 * h)
    return d / val


def curvature(s, char_info):
    """곡률 κ = |L|²"""
    L = connection(s, char_info)
    return float(abs(L)**2)


def monodromy_contour(t, char_info, radius=0.5, n_steps=64):
    """폐곡선 적분으로 모노드로미 계산"""
    center = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
    phase_accum = mpmath.mpf(0)
    prev_val = completed_L(center + radius, char_info)

    for k in range(1, n_steps + 1):
        theta = 2 * mpmath.pi * k / n_steps
        point = center + radius * mpmath.exp(1j * theta)
        curr_val = completed_L(point, char_info)
        if abs(prev_val) > 0 and abs(curr_val) > 0:
            ratio = curr_val / prev_val
            phase_accum += mpmath.im(mpmath.log(ratio))
        prev_val = curr_val

    return float(phase_accum / mpmath.pi)


def find_zeros(char_info, t_min, t_max, n_scan=2000):
    """부호 변화 + findroot로 영점 탐색"""
    ts = np.linspace(t_min, t_max, n_scan)
    zeros = []
    fail_count = 0

    prev_re, prev_t = None, None
    for t in ts:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        curr_re = mpmath.re(completed_L(s, char_info))
        if prev_re is not None and prev_re * curr_re < 0:
            try:
                def f_real(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    return mpmath.re(completed_L(sv, char_info))
                mid = mpmath.mpf(str((prev_t + float(t)) / 2))
                tz = float(mpmath.findroot(f_real, mid))
                if not zeros or abs(tz - zeros[-1]) > 0.1:
                    zeros.append(tz)
            except Exception as e:
                fail_count += 1
                print(f"  WARNING: findroot failed at t≈{(prev_t+float(t))/2:.2f}: {e}",
                      flush=True)
        prev_re, prev_t = curr_re, float(t)

    if fail_count > 0:
        print(f"  ⚠️ findroot 실패 {fail_count}회", flush=True)
    if len(zeros) == 0:
        print(f"  ⚠️ 영점 0개 — 탐색 로직 점검 필요!", flush=True)

    zeros.sort()
    return np.array(zeros)


def count_phase_jumps(char_info, sigma, t_min, t_max, n_points=1500,
                      threshold=2.5):
    """σ 고정 직선에서 위상 점프 수"""
    ts = np.linspace(t_min, t_max, n_points)
    jumps = 0
    prev_phase = None
    for t in ts:
        s = mpmath.mpc(sigma, t)
        val = completed_L(s, char_info)
        if abs(val) > 0:
            phase = float(mpmath.arg(val))
            if prev_phase is not None:
                diff = abs(phase - prev_phase)
                if diff > threshold:
                    jumps += 1
            prev_phase = phase
    return jumps


def energy_at_sigma(char_info, sigma, t_min, t_max, n_points=500):
    """σ 고정선에서 에너지 E = ∑|κ|"""
    ts = np.linspace(t_min, t_max, n_points)
    total = 0.0
    for t in ts:
        s = mpmath.mpc(sigma, t)
        k = curvature(s, char_info)
        total += min(k, 1e10)
    return total


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

if __name__ == '__main__':
    results = {}

    for name, ci in CHARACTERS.items():
        print(f"\n{'─'*50}", flush=True)
        print(f"  {name} ({ci['type']}) 분석 중...", flush=True)
        print(f"{'─'*50}", flush=True)
        t0 = time.time()
        r = {}

        # 1. 영점 탐색
        print("  [1/5] 영점 탐색...", flush=True)
        zeros = find_zeros(ci, T_MIN, T_MAX)
        r['n_zeros'] = len(zeros)
        print(f"        영점 {len(zeros)}개 발견", flush=True)

        # 2. κ 집중도 (영점 근방 δ=0.05 오프셋에서 측정, 정확히 영점 위 발산 방지)
        print("  [2/5] 곡률 집중도 측정...", flush=True)
        if len(zeros) > 0:
            near_kappa = []
            delta_offset = 0.05
            for z in zeros:
                s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(z + delta_offset))
                near_kappa.append(curvature(s, ci))

            far_ts = np.linspace(T_MIN, T_MAX, 200)
            far_kappa = []
            for t in far_ts:
                if all(abs(t - z) > 1.0 for z in zeros):
                    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
                    far_kappa.append(curvature(s, ci))

            if far_kappa and np.median(far_kappa) > 0:
                r['kappa_ratio'] = np.median(near_kappa) / np.median(far_kappa)
            else:
                r['kappa_ratio'] = float('nan')
            print(f"        κ 비율: {r['kappa_ratio']:.1f}×", flush=True)
        else:
            r['kappa_ratio'] = float('nan')
            print("        κ 비율: N/A (영점 없음)", flush=True)

        # 3. 모노드로미
        print("  [3/5] 모노드로미 양자화 검증...", flush=True)
        if len(zeros) > 0:
            monos = [monodromy_contour(z, ci) for z in zeros]
            r['mono_mean'] = np.mean(np.abs(monos))
            r['mono_std'] = np.std([abs(abs(m) - 1.0) for m in monos])
            print(f"        평균 |mono|/π = {r['mono_mean']:.4f}, "
                  f"편차 = {r['mono_std']:.4f}", flush=True)
        else:
            r['mono_mean'] = float('nan')
            r['mono_std'] = float('nan')
            print("        모노드로미: N/A", flush=True)

        # 4. 에너지 집중도
        print("  [4/5] 에너지 집중도 측정...", flush=True)
        e_half = energy_at_sigma(ci, 0.5, T_MIN, T_MAX)
        e_off = energy_at_sigma(ci, 0.3, T_MIN, T_MAX)
        r['energy_ratio'] = e_half / e_off if e_off > 0 else float('nan')
        print(f"        E(0.5)/E(0.3) = {r['energy_ratio']:.1f}×", flush=True)

        # 5. 위상 점프 유일성
        print("  [5/5] 위상 점프 유일성 검사...", flush=True)
        j_half = count_phase_jumps(ci, 0.5, T_MIN, T_MAX)
        j_others = max(
            count_phase_jumps(ci, 0.3, T_MIN, T_MAX),
            count_phase_jumps(ci, 0.7, T_MIN, T_MAX),
            count_phase_jumps(ci, 0.25, T_MIN, T_MAX),
        )
        r['jumps_half'] = j_half
        r['jumps_other'] = j_others
        print(f"        σ=0.5: {j_half} jumps, 타 σ 최대: {j_others} jumps",
              flush=True)

        elapsed = time.time() - t0
        r['time'] = elapsed
        print(f"  [소요 시간: {elapsed:.0f}초]", flush=True)

        results[name] = r

    # ━━━ 결과 출력 + 저장 ━━━

    print("\n" + "=" * 72, flush=True)
    print("고 conductor 디리클레 다발 검증 결과", flush=True)
    print("=" * 72, flush=True)

    header = f"{'측정 항목':<30}" + "".join(f"{n:>14}" for n in results.keys())
    print(header, flush=True)
    print("-" * len(header), flush=True)

    def row(label, key, fmt=".1f"):
        vals = "".join(f"{results[n].get(key, 'N/A'):>14{fmt}}"
                       if isinstance(results[n].get(key), (int, float))
                       and not np.isnan(results[n].get(key, float('nan')))
                       else f"{'N/A':>14}"
                       for n in results.keys())
        print(f"{label:<30}{vals}", flush=True)

    for name, r in results.items():
        print(f"\n{name} ({CHARACTERS[name]['type']}):", flush=True)
        print(f"  영점: {r['n_zeros']}", flush=True)
        print(f"  κ 비율: {r['kappa_ratio']:.1f}×" if not np.isnan(r['kappa_ratio']) else "  κ 비율: N/A", flush=True)
        print(f"  |mono|/π: {r['mono_mean']:.4f} ± {r['mono_std']:.4f}" if not np.isnan(r['mono_mean']) else "  mono: N/A", flush=True)
        print(f"  E(0.5)/E(0.3): {r['energy_ratio']:.1f}×" if not np.isnan(r['energy_ratio']) else "  E ratio: N/A", flush=True)
        print(f"  σ=0.5 점프: {r['jumps_half']}, 타σ 최대: {r['jumps_other']}", flush=True)

    # χ₄ vs χ₈ 에너지 비교 (실수 지표 이상치 검증)
    if 'χ₈ (mod 8)' in results:
        e8 = results['χ₈ (mod 8)']['energy_ratio']
        print(f"\n실수 지표 에너지 비교:", flush=True)
        print(f"  χ₄ (mod 4): 132.5×  [이전 결과]", flush=True)
        print(f"  χ₈ (mod 8): {e8:.1f}×  [이번 결과]", flush=True)
        if not np.isnan(e8) and e8 > 50:
            print(f"  → 실수 지표 에너지 돌출 확인!", flush=True)
        elif not np.isnan(e8):
            print(f"  → conductor 4 특수성 가능 (χ₈은 돌출 아님)", flush=True)

    # 판정
    print(f"\n{'='*72}", flush=True)
    all_pass = True
    for name, r in results.items():
        # 폐곡선 적분 모노드로미: 단순영점 → 2π → |mono|/π = 2.0
        passed = (r['n_zeros'] > 0
                  and not np.isnan(r['kappa_ratio']) and r['kappa_ratio'] > 10
                  and not np.isnan(r['mono_mean']) and abs(r['mono_mean'] - 2.0) < 0.1
                  and r['jumps_half'] > 0 and r['jumps_other'] == 0)
        status = "PASS" if passed else "PARTIAL"
        if not passed:
            all_pass = False
        print(f"  {name}: {status}", flush=True)

    if all_pass:
        print(f"\n  ✓ 보편성 확장 성공: mod 7, 8, 11 모두 다발 구조 확인", flush=True)
    print("=" * 72, flush=True)

    # 파일 저장
    out_path = os.path.expanduser(
        '~/Desktop/gdl_unified/results/dirichlet_high_conductor.txt')
    with open(out_path, 'w') as f:
        from datetime import datetime
        f.write("=" * 72 + "\n")
        f.write(f"고 conductor 디리클레 다발 검증\n")
        f.write(f"구간: t∈[{T_MIN}, {T_MAX}]\n")
        f.write(f"날짜: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"정밀도: {mpmath.mp.dps} 자릿수\n")
        f.write("=" * 72 + "\n\n")

        for name, r in results.items():
            f.write(f"{name} ({CHARACTERS[name]['type']}):\n")
            f.write(f"  영점: {r['n_zeros']}\n")
            if not np.isnan(r['kappa_ratio']):
                f.write(f"  κ 비율: {r['kappa_ratio']:.1f}×\n")
            else:
                f.write(f"  κ 비율: N/A\n")
            if not np.isnan(r['mono_mean']):
                f.write(f"  |mono|/π: {r['mono_mean']:.4f} ± {r['mono_std']:.4f}\n")
            else:
                f.write(f"  mono: N/A\n")
            if not np.isnan(r['energy_ratio']):
                f.write(f"  E(0.5)/E(0.3): {r['energy_ratio']:.1f}×\n")
            else:
                f.write(f"  E ratio: N/A\n")
            f.write(f"  σ=0.5 점프: {r['jumps_half']}, 타σ 최대: {r['jumps_other']}\n")
            f.write(f"  소요: {r['time']:.0f}초\n\n")

        if 'χ₈ (mod 8)' in results:
            e8 = results['χ₈ (mod 8)']['energy_ratio']
            f.write(f"실수 지표 에너지 비교:\n")
            f.write(f"  χ₄ (mod 4): 132.5× [기존]\n")
            f.write(f"  χ₈ (mod 8): {e8:.1f}× [신규]\n\n")

    print(f"\n결과 저장: {out_path}", flush=True)
