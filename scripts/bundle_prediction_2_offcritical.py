"""
=============================================================================
[Project RDL] 다발 예측 실험 #2: 임계선 밖 확장 (σ ≠ 1/2)
=============================================================================
Conjecture 1 검증: 임계선 밖 천 수 = 0

방법:
  1. σ ∈ {0.1, 0.2, 0.3, 0.4, 0.49, 0.51, 0.6, 0.7, 0.8, 0.9}에서
     t ∈ [10, 50] 구간의 위상 변화 Δarg(ξ) 계산
  2. σ = 1/2에서만 ±π 점프가 발생하는지 확인
  3. 곡률 밀도 κ(σ, t)의 2D 맵 작성
  4. 에너지 E(σ) = ∫|L(σ+it)|² dt의 σ-프로파일 → σ=1/2에서 극대

예측:
  - σ ≠ 1/2: Δarg ≈ 0 (위상 점프 없음), κ 유한
  - σ = 1/2: Δarg = ±π 영점마다, κ → ∞ 영점에서
  - E(σ) 프로파일이 σ=1/2에서 뾰족한 피크
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))

mpmath.mp.dps = 30

# ─── 기본 함수 ───

def xi_func(s):
    half = mpmath.mpf('0.5')
    return half * s * (s - 1) * mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s)

def L_func(s):
    h = mpmath.mpf('1e-15')
    xi_val = xi_func(s)
    if abs(xi_val) < 1e-40:
        return mpmath.mpc('inf')
    xi_d = (xi_func(s + h) - xi_func(s - h)) / (2 * h)
    return xi_d / xi_val

def curvature(s):
    L = L_func(s)
    return float(abs(L)**2)


# ─── 실험 1: σ별 위상 점프 수 ───

def count_phase_jumps(sigma, t_min=10.0, t_max=50.0, n_points=2000):
    """σ 고정, t를 스캔하며 |Δarg| > π/2인 점프 수 세기"""
    ts = np.linspace(t_min, t_max, n_points)
    jumps = 0
    prev_arg = None

    for t in ts:
        s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
        xi_val = xi_func(s)
        if abs(xi_val) < 1e-40:
            prev_arg = None
            continue
        curr_arg = float(mpmath.arg(xi_val))

        if prev_arg is not None:
            delta = curr_arg - prev_arg
            # branch cut 보정
            while delta > np.pi:
                delta -= 2 * np.pi
            while delta < -np.pi:
                delta += 2 * np.pi

            if abs(delta) > np.pi / 2:
                jumps += 1

        prev_arg = curr_arg

    return jumps


# ─── 실험 2: 곡률 2D 맵 ───

def curvature_2d_map(sigmas, t_points):
    """σ × t 격자에서 곡률 밀도 맵 계산"""
    kappa_map = np.zeros((len(sigmas), len(t_points)))
    for i, sigma in enumerate(sigmas):
        for j, t in enumerate(t_points):
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            try:
                kappa_map[i, j] = curvature(s)
            except Exception as e:
                print(f"WARNING: curvature({sigma:.3f},{t:.3f}) 실패: {e}")
                kappa_map[i, j] = np.nan
    return kappa_map


# ─── 실험 3: 에너지 σ-프로파일 ───

def energy_profile(sigmas, t_min=13.0, t_max=34.0, n_t=200):
    """E(σ) = ∫|L(σ+it)|² dt (사다리꼴 적분)"""
    ts = np.linspace(t_min, t_max, n_t)
    dt = (t_max - t_min) / (n_t - 1)
    energies = []

    for sigma in sigmas:
        vals = []
        for t in ts:
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            try:
                kappa = curvature(s)
                if np.isfinite(kappa):
                    vals.append(kappa)
                else:
                    vals.append(1e6)  # cap
            except Exception as e:
                print(f"WARNING: energy({sigma:.3f},{t:.3f}) 실패: {e}")
                vals.append(0.0)
        E = np.trapezoid(vals, dx=dt)
        energies.append(E)

    return np.array(energies)


# ─── 메인 ───

def main():
    out_path = os.path.expanduser(
        '~/Desktop/gdl_unified/outputs/analysis/bundle_prediction_offcritical.txt'
    )

    print("=" * 70)
    print("다발 예측 실험 #2: 임계선 밖 확장")
    print("=" * 70)

    # ── 실험 1: 위상 점프 수 ──
    print("\n[1] σ별 위상 점프 수 (t∈[10,50], 2000점)")
    sigmas_jump = [0.1, 0.2, 0.3, 0.4, 0.49, 0.5, 0.51, 0.6, 0.7, 0.8, 0.9]
    jump_results = {}

    for sigma in sigmas_jump:
        n_jumps = count_phase_jumps(sigma, n_points=2000)
        jump_results[sigma] = n_jumps
        marker = "★" if sigma == 0.5 else " "
        print(f"  {marker} σ={sigma:.2f}: {n_jumps} jumps")

    # ── 실험 2: 곡률 2D 맵 (조밀) ──
    print("\n[2] 곡률 2D 맵 (σ×t 격자)")
    sigmas_map = np.linspace(0.1, 0.9, 17)
    ts_map = np.linspace(13.0, 16.0, 50)  # 영점 1개 포함 구간

    print(f"  격자: {len(sigmas_map)}σ × {len(ts_map)}t = {len(sigmas_map)*len(ts_map)} 점")
    kappa_map = curvature_2d_map(sigmas_map, ts_map)

    # σ=0.5 행 vs 나머지 비교
    idx_half = np.argmin(np.abs(sigmas_map - 0.5))
    kappa_on = kappa_map[idx_half, :]
    kappa_off = np.nanmean(kappa_map[np.arange(len(sigmas_map)) != idx_half, :], axis=0)

    ratio = np.nanmax(kappa_on) / np.nanmax(kappa_off)
    print(f"  max κ(σ=0.5) / max κ(σ≠0.5) = {ratio:.1f}×")

    # ── 실험 3: 에너지 σ-프로파일 ──
    print("\n[3] 에너지 σ-프로파일 (t∈[13,34], 200점)")
    sigmas_energy = np.linspace(0.1, 0.9, 33)
    energies = energy_profile(sigmas_energy, n_t=200)

    idx_peak = np.argmax(energies)
    print(f"  에너지 피크: σ={sigmas_energy[idx_peak]:.3f}, "
          f"E={energies[idx_peak]:.2e}")
    print(f"  E(σ=0.5)/E(σ=0.3) = {energies[16]/energies[8]:.1f}×")

    # 결과 저장
    print(f"\n결과 저장: {out_path}")
    with open(out_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("다발 예측 실험 #2: 임계선 밖 확장\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write("=" * 70 + "\n")

        f.write("\n[1] 위상 점프 수 (t∈[10,50], 2000점)\n")
        f.write(f"{'σ':>8} {'jumps':>8}\n")
        f.write("-" * 20 + "\n")
        for sigma in sigmas_jump:
            f.write(f"{sigma:>8.2f} {jump_results[sigma]:>8}\n")

        f.write(f"\n[2] 곡률 2D 맵\n")
        f.write(f"격자: {len(sigmas_map)}σ × {len(ts_map)}t\n")
        f.write(f"max κ(σ=0.5) / max κ(σ≠0.5) = {ratio:.1f}×\n")

        f.write(f"\n[3] 에너지 σ-프로파일 (t∈[13,34])\n")
        f.write(f"{'σ':>8} {'E':>14}\n")
        f.write("-" * 25 + "\n")
        for sigma, E in zip(sigmas_energy, energies):
            marker = "★" if abs(sigma - 0.5) < 0.02 else " "
            f.write(f"{marker}{sigma:>7.3f} {E:>14.4e}\n")

        f.write(f"\n에너지 피크: σ={sigmas_energy[idx_peak]:.3f}\n")
        f.write(f"E(0.5)/E(0.3) = {energies[16]/energies[8]:.1f}×\n")

    # 수학자 지시: results/offcritical_c297.txt에도 복사
    import shutil
    copy_path = os.path.expanduser(
        '~/Desktop/gdl_unified/results/offcritical_c297.txt'
    )
    shutil.copy2(out_path, copy_path)
    print(f"복사 완료: {copy_path}")

    print("완료.")


if __name__ == '__main__':
    main()
