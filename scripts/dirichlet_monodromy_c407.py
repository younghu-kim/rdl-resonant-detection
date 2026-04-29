"""
=============================================================================
[C-407] Dirichlet L-함수 모노드로미 검증 — χ mod 5 (짝수, conductor=5)
=============================================================================
목표: 모노드로미 TP/FP 분리가 ζ(s)만의 현상인지 보편적인지 판별.
  - ζ(s)에서 25 TP/60 FP 완벽 분리 확립 (C-403 + C-406).
  - 이제 degree-1 Dirichlet L-함수에서 동일 현상 확인.

방법:
  1. L(s, χ₅_even)의 영점 수집 (T∈[10,100]) → TP
  2. 영점 간 중점 + 랜덤 비영점 → FP
  3. TP/FP에서 monodromy_contour(func='dirichlet') 계산
  4. TP mono≈2π, FP mono≈0 분리 검증

χ mod 5 (짝수, 이차): χ(0)=0, χ(1)=1, χ(2)=-1, χ(3)=-1, χ(4)=1, a=0
  → Legendre symbol (·/5)

성공 기준:
  - TP |mono|/π > 1.5 비율 ≥ 70%
  - FP |mono|/π < 0.3 비율 ≥ 70%
  - KS p < 0.01
  - 3/4 충족 → 양성
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))
from bundle_utils import (
    completed_L, monodromy_contour, find_zeros_dirichlet, curvature_dirichlet,
)

# ─── 설정 ───
T_MIN, T_MAX = 10.0, 100.0
N_SCAN = 5000  # 영점 탐색 해상도
MONO_RADII = [0.1, 0.01, 0.001]
MONO_NSTEPS = 256
SEEDS = [42, 7, 123, 314, 2024]  # FP 랜덤 위치 시드

mpmath.mp.dps = 80

# 짝수 원시 χ mod 5 (이차, Legendre symbol)
# χ(n): 0→0, 1→1, 2→-1, 3→-1, 4→1
CHI5_EVEN = {
    'chi': [0, 1, -1, -1, 1],
    'q': 5,
    'a': 0,  # 짝수 (even)
    'label': 'χ₅_even (mod 5, 이차, a=0)',
}

OUT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/dirichlet_monodromy_c407.txt'
)


def compute_curvature_at(t, char_info):
    """임계선 위 점 t에서 Dirichlet 곡률 κ = |Λ'/Λ|²"""
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
    return curvature_dirichlet(s, char_info)


def generate_fp_points(true_zeros, t_min, t_max, seed, n_midpoints=None, n_random=10):
    """
    FP 후보 생성:
      1. 인접 영점 중점 (영점이 없는 위치)
      2. 랜덤 비영점 위치 (최소 거리 확인)
    """
    rng = np.random.RandomState(seed)
    fps = []

    # 1. 인접 영점 중점
    if len(true_zeros) >= 2:
        for i in range(len(true_zeros) - 1):
            mid = (true_zeros[i] + true_zeros[i + 1]) / 2.0
            fps.append(mid)

    # 2. 랜덤 위치 (영점에서 최소 1.0 이상 떨어진)
    min_dist_from_zero = 1.0
    attempts = 0
    while len(fps) < len(fps) + n_random and attempts < n_random * 10:
        t_rand = rng.uniform(t_min, t_max)
        if len(true_zeros) > 0:
            dists = np.abs(true_zeros - t_rand)
            if dists.min() < min_dist_from_zero:
                attempts += 1
                continue
        # 중복 체크
        if fps and min(abs(t_rand - f) for f in fps) < 0.5:
            attempts += 1
            continue
        fps.append(t_rand)
        n_random -= 1
        attempts += 1

    return np.array(sorted(fps))


def main():
    START = time.time()

    print("=" * 70)
    print("[C-407] Dirichlet L-함수 모노드로미 검증")
    print(f"  χ: {CHI5_EVEN['label']}")
    print(f"  T: [{T_MIN}, {T_MAX}], dps={mpmath.mp.dps}, n_steps={MONO_NSTEPS}")
    print("=" * 70)

    # ─── 1. 영점 수집 ───
    print("\n[1단계] 영점 수집 중...")
    true_zeros = find_zeros_dirichlet(CHI5_EVEN, t_min=T_MIN, t_max=T_MAX, n_scan=N_SCAN)
    print(f"  발견 영점 수: {len(true_zeros)}")

    if len(true_zeros) == 0:
        print("⚠️ 영점 0개 — 탐색 로직 점검 필요")
        # 대비: n_scan 올려서 재시도
        print("  n_scan=10000으로 재시도...")
        true_zeros = find_zeros_dirichlet(CHI5_EVEN, t_min=T_MIN, t_max=T_MAX, n_scan=10000)
        print(f"  재시도 결과: {len(true_zeros)}개")
        if len(true_zeros) == 0:
            print("❌ 영점 여전히 0개. 완비 L-함수 검증 필요. 중단.")
            return

    if len(true_zeros) < 5:
        print(f"⚠️ 영점 {len(true_zeros)}개 — 통계적으로 부족. 주의.")

    print(f"  영점 목록 (처음 10개):")
    for i, z in enumerate(true_zeros[:10]):
        print(f"    [{i+1}] t = {z:.6f}")
    if len(true_zeros) > 10:
        print(f"    ... ({len(true_zeros)}개 총)")

    # ─── 1.5. 완비 L-함수 건전성 검사 ───
    print("\n[1.5단계] 완비 L-함수 건전성 검사...")
    n_check = min(5, len(true_zeros))
    for i in range(n_check):
        t = true_zeros[i]
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        val = completed_L(s, CHI5_EVEN)
        print(f"  Λ(1/2+i·{t:.4f}) = {float(mpmath.re(val)):.2e} + {float(mpmath.im(val)):.2e}i  "
              f"|Λ| = {float(abs(val)):.2e}")

    # ─── 2. FP 생성 ───
    print("\n[2단계] FP 후보 생성...")
    # 5시드에서 각각 FP 생성 후 합산·중복제거
    all_fp = []
    for seed in SEEDS:
        fps = generate_fp_points(true_zeros, T_MIN, T_MAX, seed, n_random=10)
        all_fp.extend(fps.tolist())

    # 중복 제거 (±0.5)
    all_fp = sorted(set(all_fp))
    fp_dedup = [all_fp[0]] if all_fp else []
    for f in all_fp[1:]:
        if f - fp_dedup[-1] > 0.5:
            fp_dedup.append(f)
    fp_unique = np.array(fp_dedup)

    # 영점과 겹치는 FP 제거
    fp_clean = []
    for f in fp_unique:
        if len(true_zeros) > 0:
            if np.min(np.abs(true_zeros - f)) > 0.5:
                fp_clean.append(f)
        else:
            fp_clean.append(f)
    fp_unique = np.array(fp_clean)

    print(f"  FP 후보: {len(fp_unique)}개 (영점에서 >0.5 거리)")

    # ─── 3. 모노드로미 계산 ───
    print("\n[3단계] 모노드로미 계산...")
    results = []

    # TP: true zeros에서 모노드로미
    n_tp = min(20, len(true_zeros))
    print(f"\n  [True Positives] — {n_tp}개 영점")
    for i in range(n_tp):
        t = float(true_zeros[i])
        monos_by_r = {}
        for r in MONO_RADII:
            try:
                monos_by_r[r] = monodromy_contour(t, radius=r, n_steps=MONO_NSTEPS,
                                                   func='dirichlet', char_info=CHI5_EVEN)
            except Exception as e:
                print(f"    WARNING: monodromy r={r} t={t:.4f} 실패: {e}")
                monos_by_r[r] = 0.0

        try:
            kappa = compute_curvature_at(t, CHI5_EVEN)
        except Exception as e:
            print(f"    WARNING: curvature t={t:.4f} 실패: {e}")
            kappa = float('inf')

        mono_main = monos_by_r[0.1]
        results.append(('TP', t, mono_main, kappa, monos_by_r))
        print(f"    t={t:.4f}: "
              f"mono(r=0.1)={monos_by_r[0.1]/np.pi:.4f}π, "
              f"mono(r=0.01)={monos_by_r[0.01]/np.pi:.4f}π, "
              f"mono(r=0.001)={monos_by_r[0.001]/np.pi:.4f}π, κ={kappa:.2e}")

    # FP: 비영점에서 모노드로미
    n_fp = min(30, len(fp_unique))
    print(f"\n  [False Positives] — {n_fp}개 비영점")
    for i in range(n_fp):
        t = float(fp_unique[i])
        monos_by_r = {}
        for r in MONO_RADII:
            try:
                monos_by_r[r] = monodromy_contour(t, radius=r, n_steps=MONO_NSTEPS,
                                                   func='dirichlet', char_info=CHI5_EVEN)
            except Exception as e:
                print(f"    WARNING: monodromy r={r} t={t:.4f} 실패: {e}")
                monos_by_r[r] = 0.0

        try:
            kappa = compute_curvature_at(t, CHI5_EVEN)
        except Exception as e:
            print(f"    WARNING: curvature t={t:.4f} 실패: {e}")
            kappa = 0.0

        mono_main = monos_by_r[0.1]
        results.append(('FP', t, mono_main, kappa, monos_by_r))
        print(f"    t={t:.4f}: "
              f"mono(r=0.1)={monos_by_r[0.1]/np.pi:.4f}π, "
              f"mono(r=0.01)={monos_by_r[0.01]/np.pi:.4f}π, "
              f"mono(r=0.001)={monos_by_r[0.001]/np.pi:.4f}π, κ={kappa:.2e}")

    # ─── 4. 통계 ───
    print("\n" + "=" * 70)
    print("통계 요약")
    print("=" * 70)

    tp_results = [r for r in results if r[0] == 'TP']
    fp_results = [r for r in results if r[0] == 'FP']

    if not tp_results or not fp_results:
        print("❌ TP 또는 FP가 0개 — 통계 불가")
        return

    # 반지름별 통계
    for r in MONO_RADII:
        tp_m = [abs(res[4][r]) for res in tp_results]
        fp_m = [abs(res[4][r]) for res in fp_results]
        if tp_m and fp_m:
            ks_stat, ks_p = stats.ks_2samp(tp_m, fp_m)
            print(f"\n  radius={r}:")
            print(f"    TP |mono|/π: mean={np.mean(tp_m)/np.pi:.4f}, std={np.std(tp_m)/np.pi:.4f}")
            print(f"    FP |mono|/π: mean={np.mean(fp_m)/np.pi:.4f}, std={np.std(fp_m)/np.pi:.4f}")
            print(f"    KS test: stat={ks_stat:.4f}, p={ks_p:.6e}")

    # 주 통계 (r=0.1)
    tp_monos = [abs(r[2]) for r in tp_results]
    fp_monos = [abs(r[2]) for r in fp_results]

    print(f"\n주 통계 (r=0.1, n_steps={MONO_NSTEPS}):")
    print(f"  TP |mono|/π: mean={np.mean(tp_monos)/np.pi:.4f}, std={np.std(tp_monos)/np.pi:.4f}")
    print(f"  FP |mono|/π: mean={np.mean(fp_monos)/np.pi:.4f}, std={np.std(fp_monos)/np.pi:.4f}")

    ks_stat_main, ks_p_main = stats.ks_2samp(tp_monos, fp_monos)
    print(f"  KS test: stat={ks_stat_main:.4f}, p={ks_p_main:.6e}")

    # ─── 5. 성공 기준 ───
    print("\n" + "=" * 70)
    print("성공 기준 체크")
    print("=" * 70)

    # 기준 1: TP |mono|/π > 1.5 비율 ≥ 70%
    tp_high = sum(1 for m in tp_monos if m / np.pi > 1.5) / len(tp_monos)
    crit1 = tp_high >= 0.7
    print(f"  TP |mono|/π > 1.5 비율: {tp_high:.1%} ({'✅' if crit1 else '❌'} ≥ 70%)")

    # 기준 2: FP |mono|/π < 0.3 비율 ≥ 70%
    fp_low = sum(1 for m in fp_monos if m / np.pi < 0.3) / len(fp_monos)
    crit2 = fp_low >= 0.7
    print(f"  FP |mono|/π < 0.3 비율: {fp_low:.1%} ({'✅' if crit2 else '❌'} ≥ 70%)")

    # 기준 3: KS p < 0.01
    crit3 = ks_p_main < 0.01
    print(f"  KS p < 0.01: {'✅' if crit3 else '❌'} (p={ks_p_main:.6e})")

    # 기준 4: 이중 기준 정밀도
    mono_threshold = np.pi * 0.5
    tp_pass = sum(1 for m in tp_monos if m > mono_threshold)
    fp_pass = sum(1 for m in fp_monos if m > mono_threshold)
    if tp_pass + fp_pass > 0:
        dual_precision = tp_pass / (tp_pass + fp_pass)
    else:
        dual_precision = 0.0
    crit4 = dual_precision > 0.4 if (tp_pass + fp_pass) > 0 else False

    print(f"\n  이중 기준 (|mono| > π/2):")
    print(f"    TP 통과: {tp_pass}/{len(tp_monos)}")
    print(f"    FP 통과: {fp_pass}/{len(fp_monos)}")
    print(f"    정밀도: {dual_precision:.1%} ({'✅' if crit4 else '❌'} > 40%)")

    n_pass = sum([crit1, crit2, crit3, crit4])
    if n_pass >= 4:
        overall = "강한 양성"
    elif n_pass >= 3:
        overall = "양성"
    elif n_pass >= 2:
        overall = "약한 양성"
    else:
        overall = "음성"

    print(f"\n  *** 종합 판정: {overall} ({n_pass}/4 기준 충족) ***")

    # ─── 6. 결과 저장 ───
    elapsed = time.time() - START
    print(f"\n소요: {elapsed:.1f}초 ({elapsed/60:.1f}분)")

    with open(OUT_PATH, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("[C-407] Dirichlet L-함수 모노드로미 검증\n")
        f.write(f"χ: {CHI5_EVEN['label']}\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"설정:\n")
        f.write(f"  T: [{T_MIN}, {T_MAX}]\n")
        f.write(f"  dps: {mpmath.mp.dps}\n")
        f.write(f"  n_steps: {MONO_NSTEPS}\n")
        f.write(f"  반지름: {MONO_RADII}\n")
        f.write(f"  FP 시드: {SEEDS}\n")
        f.write(f"  χ: {CHI5_EVEN['chi']} (q={CHI5_EVEN['q']}, a={CHI5_EVEN['a']})\n\n")

        f.write(f"영점 수: {len(true_zeros)}\n")
        f.write(f"  목록: {[f'{z:.6f}' for z in true_zeros]}\n\n")

        # 결과 테이블
        header = f"{'Type':>4} {'t':>12}"
        for r in MONO_RADII:
            header += f" {'mono_r='+str(r):>14}"
        header += f" {'κ':>12}"
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")
        for res in results:
            line = f"{res[0]:>4} {res[1]:>12.6f}"
            for r in MONO_RADII:
                line += f" {res[4][r]/np.pi:>13.4f}π"
            line += f" {res[3]:>12.2e}"
            f.write(line + "\n")

        f.write(f"\n\n반지름별 통계:\n")
        for r in MONO_RADII:
            tp_m = [abs(res[4][r]) for res in tp_results]
            fp_m = [abs(res[4][r]) for res in fp_results]
            if tp_m and fp_m:
                ks_s, ks_p = stats.ks_2samp(tp_m, fp_m)
                f.write(f"\n  radius={r}:\n")
                f.write(f"    TP |mono|/π: {np.mean(tp_m)/np.pi:.4f} ± {np.std(tp_m)/np.pi:.4f}\n")
                f.write(f"    FP |mono|/π: {np.mean(fp_m)/np.pi:.4f} ± {np.std(fp_m)/np.pi:.4f}\n")
                f.write(f"    KS: stat={ks_s:.4f}, p={ks_p:.6e}\n")

        f.write(f"\n\n성공 기준:\n")
        f.write(f"  TP |mono|/π > 1.5 비율: {tp_high:.1%} ({'PASS' if crit1 else 'FAIL'})\n")
        f.write(f"  FP |mono|/π < 0.3 비율: {fp_low:.1%} ({'PASS' if crit2 else 'FAIL'})\n")
        f.write(f"  KS p-value: {ks_p_main:.6e} ({'PASS' if crit3 else 'FAIL'})\n")
        f.write(f"  이중기준 정밀도: {dual_precision:.1%} ({'PASS' if crit4 else 'FAIL'})\n")
        f.write(f"\n  종합 판정: {overall} ({n_pass}/4)\n")
        f.write(f"\n소요: {elapsed:.1f}초\n")

    print(f"\n결과 저장: {OUT_PATH}")
    print("완료.")


if __name__ == '__main__':
    main()
