"""
=============================================================================
[C-409] EC L-함수 모노드로미 TP/FP 분리 검증 — 37a1 (degree 2, rank 1)
=============================================================================
목표: 모노드로미 TP/FP 분리가 rank 1 EC에서도 성립하는지 확인.
  - C-408: 11a1 (rank 0, N=11) → 강한 양성 4/4 (20TP/30FP)
  - 이번: 37a1 (rank 1, N=37) → rank가 달라도 동일한가?

방법: C-408과 동일.
  1. 37a1 EC L-함수 영점 수집 (PARI lfunzeros)
  2. FP 후보 생성 (영점 간 중점 + 랜덤)
  3. TP/FP 모노드로미 계산 (폐곡선 적분, center=1.0)
  4. TP mono≈2π, FP mono≈0 분리 검증

성공 기준:
  - TP |mono|/π > 1.5 비율 ≥ 70%
  - FP |mono|/π < 0.3 비율 ≥ 70%
  - KS p < 0.01
  - 이중기준 정밀도 > 40%
  - 3/4 충족 → 양성

LMFDB 확인: 37a1 = [0, 0, 1, -1, 0], rank=1, conductor=37
=============================================================================
"""

import sys, os, time
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)
print(f"PARI 초기화: 2GB 메모리, realprecision=100")

# ─── 설정 ───
T_MAX = 50.0
MONO_RADII = [0.1, 0.05, 0.01]
MONO_NSTEPS = 128
SEEDS = [42, 7, 123, 314, 2024]

CURVE_INFO = {
    'label': '37a1',
    'coeffs': [0, 0, 1, -1, 0],   # LMFDB 확인: y² + y = x³ - x
    'rank': 1,
    'N': 37,
}

OUT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/ec_monodromy_c409.txt'
)


def init_ec_lfun(info):
    """PARI EC L-함수 초기화"""
    E = pari.ellinit(info['coeffs'])
    lf = pari.lfuncreate(E)
    eps = int(pari.ellrootno(E))
    # lf 직접 사용 (정확성 우선, C-408과 동일)
    linit = lf
    print(f"  lf 직접 사용 (정확성 우선)")
    return E, lf, linit, eps


def ec_lambda(linit_or_lf, sigma, t):
    """PARI로 Λ(s) 계산. s = sigma + i*t"""
    s = pari(f"{sigma} + {t}*I")
    try:
        val = complex(pari.lfunlambda(linit_or_lf, s))
        return val
    except Exception:
        return None


def monodromy_contour_ec(linit, t_zero, radius=0.1, n_steps=256, center=1.0):
    """
    EC L-함수의 점 s=center+i*t_zero 주위 반지름 radius 원에서 모노드로미 계산.
    EC 임계선 = σ=1 (weight 2, FE: Λ(s)=ε·Λ(2-s)).
    """
    total_delta = 0.0
    prev_arg = None

    for k in range(n_steps + 1):
        theta = 2 * np.pi * k / n_steps
        sigma = center + radius * np.cos(theta)
        t = t_zero + radius * np.sin(theta)

        val = ec_lambda(linit, sigma, t)
        if val is None or abs(val) < 1e-300:
            continue

        curr_arg = np.angle(val)

        if prev_arg is not None:
            delta = curr_arg - prev_arg
            while delta > np.pi:
                delta -= 2 * np.pi
            while delta < -np.pi:
                delta += 2 * np.pi
            total_delta += delta

        prev_arg = curr_arg

    return total_delta


def curvature_ec(linit, t, delta=0.01):
    """EC L-함수의 곡률 κ = |Λ'/Λ|² 근사 (차분). EC 임계선=σ=1"""
    s_center = ec_lambda(linit, 1.0 + delta, t)
    s_plus = ec_lambda(linit, 1.0 + delta, t + delta)
    s_minus = ec_lambda(linit, 1.0 + delta, t - delta)
    if s_center is None or s_plus is None or s_minus is None:
        return 0.0
    if abs(s_center) < 1e-300:
        return float('inf')
    deriv = (s_plus - s_minus) / (2 * delta * 1j)
    conn = deriv / s_center
    return abs(conn) ** 2


def generate_fp_points(true_zeros, t_min, t_max, seed, n_random=10):
    """FP 후보 생성: 인접 영점 중점 + 랜덤 비영점"""
    rng = np.random.RandomState(seed)
    fps = []

    # 1. 인접 영점 중점
    if len(true_zeros) >= 2:
        for i in range(len(true_zeros) - 1):
            mid = (true_zeros[i] + true_zeros[i + 1]) / 2.0
            fps.append(mid)

    # 2. 랜덤 (영점에서 >1.0 거리)
    min_dist = 1.0
    added = 0
    attempts = 0
    while added < n_random and attempts < n_random * 20:
        t_rand = rng.uniform(t_min, t_max)
        if len(true_zeros) > 0:
            dists = np.abs(np.array(true_zeros) - t_rand)
            if dists.min() < min_dist:
                attempts += 1
                continue
        if fps and min(abs(t_rand - f) for f in fps) < 0.5:
            attempts += 1
            continue
        fps.append(t_rand)
        added += 1
        attempts += 1

    return np.array(sorted(fps))


def main():
    START = time.time()

    print("=" * 70)
    print(f"[C-409] EC L-함수 모노드로미 TP/FP 검증")
    print(f"  곡선: {CURVE_INFO['label']} (N={CURVE_INFO['N']}, rank={CURVE_INFO['rank']})")
    print(f"  T_MAX: {T_MAX}, n_steps={MONO_NSTEPS}")
    print("=" * 70)

    # ─── 1. 초기화 ───
    print("\n[1단계] EC L-함수 초기화...")
    E, lf, linit, eps = init_ec_lfun(CURVE_INFO)
    print(f"  root number ε = {eps}")

    # ─── 2. 영점 수집 ───
    print("\n[2단계] 영점 수집 중...")
    zeros_raw = pari.lfunzeros(lf, T_MAX)
    true_zeros = [float(z) for z in zeros_raw if float(z) > 1.0]
    print(f"  발견 영점 수: {len(true_zeros)}")

    if len(true_zeros) < 3:
        print("  영점 부족. T_MAX 확장 필요.")
        return

    print(f"  영점 목록 (처음 10개):")
    for i, z in enumerate(true_zeros[:10]):
        print(f"    [{i+1}] t = {z:.6f}")
    if len(true_zeros) > 10:
        print(f"    ... ({len(true_zeros)}개 총)")

    # ─── 2.5. 건전성 검사 ───
    print("\n[2.5단계] Λ 건전성 검사...")
    print("  EC 임계선: σ=1 (FE Λ(s)=ε·Λ(2-s), weight 2)")
    for z in true_zeros[:5]:
        val = ec_lambda(lf, 1.0, z)
        if val is not None:
            print(f"  Λ(1+i·{z:.4f}) = {val.real:.2e} + {val.imag:.2e}i  |Λ| = {abs(val):.2e}")
    lfunc = lf

    # ─── 3. FP 생성 ───
    print("\n[3단계] FP 후보 생성...")
    all_fp = []
    for seed in SEEDS:
        fps = generate_fp_points(true_zeros, 2.0, T_MAX, seed, n_random=10)
        all_fp.extend(fps.tolist())

    # 중복 제거
    all_fp = sorted(set(all_fp))
    fp_dedup = [all_fp[0]] if all_fp else []
    for f in all_fp[1:]:
        if f - fp_dedup[-1] > 0.5:
            fp_dedup.append(f)
    fp_unique = np.array(fp_dedup)

    # 영점과 겹치는 FP 제거
    fp_clean = []
    for f in fp_unique:
        if np.min(np.abs(np.array(true_zeros) - f)) > 0.5:
            fp_clean.append(f)
    fp_unique = np.array(fp_clean)
    print(f"  FP 후보: {len(fp_unique)}개 (영점에서 >0.5 거리)")

    # ─── 4. 모노드로미 계산 ───
    print("\n[4단계] 모노드로미 계산...")
    results = []

    n_tp = min(20, len(true_zeros))
    print(f"\n  [True Positives] — {n_tp}개 영점")
    for i in range(n_tp):
        t = true_zeros[i]
        monos_by_r = {}
        for r in MONO_RADII:
            try:
                monos_by_r[r] = monodromy_contour_ec(lfunc, t, radius=r, n_steps=MONO_NSTEPS)
            except Exception as e:
                print(f"    WARNING: monodromy r={r} t={t:.4f} 실패: {e}")
                monos_by_r[r] = 0.0

        kappa = curvature_ec(lfunc, t)
        mono_main = monos_by_r[MONO_RADII[0]]
        results.append(('TP', t, mono_main, kappa, monos_by_r))
        parts = [f"mono(r={r})={monos_by_r[r]/np.pi:.4f}π" for r in MONO_RADII]
        print(f"    t={t:.4f}: {', '.join(parts)}, κ={kappa:.2e}")

    n_fp = min(30, len(fp_unique))
    print(f"\n  [False Positives] — {n_fp}개 비영점")
    for i in range(n_fp):
        t = float(fp_unique[i])
        monos_by_r = {}
        for r in MONO_RADII:
            try:
                monos_by_r[r] = monodromy_contour_ec(lfunc, t, radius=r, n_steps=MONO_NSTEPS)
            except Exception as e:
                print(f"    WARNING: monodromy r={r} t={t:.4f} 실패: {e}")
                monos_by_r[r] = 0.0

        kappa = curvature_ec(lfunc, t)
        mono_main = monos_by_r[MONO_RADII[0]]
        results.append(('FP', t, mono_main, kappa, monos_by_r))
        parts = [f"mono(r={r})={monos_by_r[r]/np.pi:.4f}π" for r in MONO_RADII]
        print(f"    t={t:.4f}: {', '.join(parts)}, κ={kappa:.2e}")

    # ─── 5. 통계 ───
    print("\n" + "=" * 70)
    print("통계 요약")
    print("=" * 70)

    tp_results = [r for r in results if r[0] == 'TP']
    fp_results = [r for r in results if r[0] == 'FP']

    if not tp_results or not fp_results:
        print("  TP 또는 FP가 0개 — 통계 불가")
        return

    for r in MONO_RADII:
        tp_m = [abs(res[4][r]) for res in tp_results]
        fp_m = [abs(res[4][r]) for res in fp_results]
        if tp_m and fp_m:
            ks_stat, ks_p = stats.ks_2samp(tp_m, fp_m)
            print(f"\n  radius={r}:")
            print(f"    TP |mono|/π: mean={np.mean(tp_m)/np.pi:.4f}, std={np.std(tp_m)/np.pi:.4f}")
            print(f"    FP |mono|/π: mean={np.mean(fp_m)/np.pi:.4f}, std={np.std(fp_m)/np.pi:.4f}")
            print(f"    KS test: stat={ks_stat:.4f}, p={ks_p:.6e}")

    tp_monos = [abs(r[2]) for r in tp_results]
    fp_monos = [abs(r[2]) for r in fp_results]

    print(f"\n주 통계 (r=0.1, n_steps={MONO_NSTEPS}):")
    print(f"  TP |mono|/π: mean={np.mean(tp_monos)/np.pi:.4f}, std={np.std(tp_monos)/np.pi:.4f}")
    print(f"  FP |mono|/π: mean={np.mean(fp_monos)/np.pi:.4f}, std={np.std(fp_monos)/np.pi:.4f}")

    ks_stat_main, ks_p_main = stats.ks_2samp(tp_monos, fp_monos)
    print(f"  KS test: stat={ks_stat_main:.4f}, p={ks_p_main:.6e}")

    # ─── 6. 성공 기준 ───
    print("\n" + "=" * 70)
    print("성공 기준 체크")
    print("=" * 70)

    tp_high = sum(1 for m in tp_monos if m / np.pi > 1.5) / len(tp_monos)
    crit1 = tp_high >= 0.7
    print(f"  TP |mono|/π > 1.5 비율: {tp_high:.1%} ({'PASS' if crit1 else 'FAIL'} >= 70%)")

    fp_low = sum(1 for m in fp_monos if m / np.pi < 0.3) / len(fp_monos)
    crit2 = fp_low >= 0.7
    print(f"  FP |mono|/π < 0.3 비율: {fp_low:.1%} ({'PASS' if crit2 else 'FAIL'} >= 70%)")

    crit3 = ks_p_main < 0.01
    print(f"  KS p < 0.01: {'PASS' if crit3 else 'FAIL'} (p={ks_p_main:.6e})")

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
    print(f"    정밀도: {dual_precision:.1%} ({'PASS' if crit4 else 'FAIL'} > 40%)")

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

    # ─── 7. 결과 저장 ───
    elapsed = time.time() - START
    print(f"\n소요: {elapsed:.1f}초 ({elapsed/60:.1f}분)")

    with open(OUT_PATH, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("[C-409] EC L-함수 모노드로미 TP/FP 검증\n")
        f.write(f"곡선: {CURVE_INFO['label']} (N={CURVE_INFO['N']}, rank={CURVE_INFO['rank']})\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"설정:\n")
        f.write(f"  T_MAX: {T_MAX}\n")
        f.write(f"  PARI realprecision: 100\n")
        f.write(f"  n_steps: {MONO_NSTEPS}\n")
        f.write(f"  반지름: {MONO_RADII}\n")
        f.write(f"  FP 시드: {SEEDS}\n")
        f.write(f"  root number ε: {eps}\n")
        f.write(f"  LMFDB 계수: {CURVE_INFO['coeffs']}\n\n")

        f.write(f"영점 수: {len(true_zeros)}\n")
        f.write(f"  목록: {[f'{z:.6f}' for z in true_zeros]}\n\n")

        header = f"{'Type':>4} {'t':>12}"
        for r in MONO_RADII:
            header += f" {'mono_r=' + str(r):>14}"
        header += f" {'kappa':>12}"
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")
        for res in results:
            line = f"{res[0]:>4} {res[1]:>12.6f}"
            for r in MONO_RADII:
                line += f" {res[4][r] / np.pi:>13.4f}π"
            line += f" {res[3]:>12.2e}"
            f.write(line + "\n")

        f.write(f"\n\n반지름별 통계:\n")
        for r in MONO_RADII:
            tp_m = [abs(res[4][r]) for res in tp_results]
            fp_m = [abs(res[4][r]) for res in fp_results]
            if tp_m and fp_m:
                ks_s, ks_p = stats.ks_2samp(tp_m, fp_m)
                f.write(f"\n  radius={r}:\n")
                f.write(f"    TP |mono|/π: {np.mean(tp_m) / np.pi:.4f} +/- {np.std(tp_m) / np.pi:.4f}\n")
                f.write(f"    FP |mono|/π: {np.mean(fp_m) / np.pi:.4f} +/- {np.std(fp_m) / np.pi:.4f}\n")
                f.write(f"    KS: stat={ks_s:.4f}, p={ks_p:.6e}\n")

        f.write(f"\n\n성공 기준:\n")
        f.write(f"  TP |mono|/pi > 1.5 비율: {tp_high:.1%} ({'PASS' if crit1 else 'FAIL'})\n")
        f.write(f"  FP |mono|/pi < 0.3 비율: {fp_low:.1%} ({'PASS' if crit2 else 'FAIL'})\n")
        f.write(f"  KS p-value: {ks_p_main:.6e} ({'PASS' if crit3 else 'FAIL'})\n")
        f.write(f"  이중기준 정밀도: {dual_precision:.1%} ({'PASS' if crit4 else 'FAIL'})\n")
        f.write(f"\n  종합 판정: {overall} ({n_pass}/4)\n")

        f.write(f"\n누적 (C-403+406+407+408+409):\n")
        f.write(f"  L-함수: zeta + L(s,chi5) + L(s,11a1) + L(s,37a1)\n")
        f.write(f"  degree: 1-2, rank 0-1\n")
        f.write(f"  이전: 65 TP / 120 FP\n")
        f.write(f"  이번: {n_tp} TP / {n_fp} FP\n")
        f.write(f"  합산: {65 + n_tp} TP / {120 + n_fp} FP\n")

        f.write(f"\n소요: {elapsed:.1f}초\n")

    print(f"\n결과 저장: {OUT_PATH}")
    print("완료.")


if __name__ == '__main__':
    main()
