"""
=============================================================================
[C-415] Dedekind ζ_K 모노드로미 TP/FP 분리 검증
=============================================================================
목표: Conjecture 1 (모노드로미 보편성)의 첫 번째 예측-검증 테스트.
  기존 9 L-함수(GL(1)-GL(5), sym^k 타워)와 완전히 다른 구성 경로.
  수체 Q(α), α³-α-1=0의 Dedekind zeta ζ_K(s).

수학적 배경:
  - K = Q(α), f(x) = x³-x-1, disc(K) = -23
  - Gal(K̃/Q) ≅ S₃ (비가환), [K:Q] = 3
  - ζ_K(s) = ζ(s) · L(s, ρ₂) — degree 3 L-함수
  - PARI: lfuncreate(nfinit(x³-x-1)) → Dedekind zeta 직접 구성
  - conductor: |disc(K)| = 23
  - gammaV: [0, 0, 1] (1 real + 1 complex place)
  - 임계선: σ = 1/2 (표준)

예측 (Conjecture 1):
  - TP (영점): mono/π = 2.000, κδ² → 1
  - FP (비영점): mono/π = 0.000, κδ² ≠ 1
  - 이중기준 100% 분리

구성 경로의 독립성:
  - 기존: GL(1) Dirichlet, GL(2) EC/Maass/Δ, GL(3-5) sym^k
  - 이번: Dedekind zeta — 수체의 산술에서 직접 유래, sym^k와 무관
  - ζ_K의 영점은 ζ(s)의 영점과 L(s,ρ₂)의 영점을 모두 포함
  - 논증 원리에 의해 모든 단순 영점에서 mono=2π 예측

성공 기준:
  1. TP ≥ 15개, mono/π = 2.000 ± 0.01 (r=0.01 기준)
  2. FP ≥ 20개, mono/π = 0.000 ± 0.01 (r=0.01 기준)
  3. KS p < 1e-6 (r=0.01 기준)
  4. 이중기준(κ + mono) 100% 분리 (r=0.01 기준)
=============================================================================
"""

import sys, os, time
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 80)
print(f"PARI 초기화: 2GB 메모리, realprecision=80")

# ─── 설정 ───
T_MAX = 50.0
T_ZEROS_MAX = 70.0
MONO_RADII = [0.1, 0.01, 0.001]
MONO_NSTEPS = 128
PRIMARY_RADIUS_IDX = 1  # r=0.01

INFO = {
    'label': 'dedekind_s3',
    'source': 'ζ_K, K=Q(α), α³-α-1=0, Gal≅S₃',
    'polynomial': 'x^3-x-1',
    'discriminant': -23,
    'degree': 3,
    'sigma_c': 0.5,  # 표준 Dedekind zeta 임계선
}

OUT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/dedekind_s3_monodromy_c415.txt'
)


def init_dedekind_lfun():
    """PARI Dedekind zeta 초기화"""
    nf = pari(f'nfinit({INFO["polynomial"]})')
    lf = pari(f'lfuncreate({nf})')

    params = pari.lfunparams(lf)
    print(f"  lfunparams: {params}")
    print(f"  conductor: {INFO['discriminant']}")
    print(f"  gammaV: [0, 0, 1] (1 real + 1 complex embedding)")

    return lf, nf


def dedekind_lambda(lf, sigma, t):
    """PARI로 Λ(s) 계산. s = sigma + i*t"""
    s = pari(f"{sigma} + {t}*I")
    try:
        val = complex(pari.lfunlambda(lf, s))
        return val
    except Exception:
        return None


def monodromy_contour(lf, t_zero, radius=0.1, n_steps=128, center=0.5):
    """
    Dedekind zeta의 점 s=center+i*t_zero 주위 폐곡선 모노드로미.
    arg(Λ(s)) 누적 변화 계산 (인수 원리).
    """
    total_delta = 0.0
    prev_arg = None

    for k in range(n_steps + 1):
        theta = 2 * np.pi * k / n_steps
        sigma = center + radius * np.cos(theta)
        t = t_zero + radius * np.sin(theta)

        val = dedekind_lambda(lf, sigma, t)
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


def curvature(lf, t, center=0.5, delta_h=0.01):
    """κ = |Λ'/Λ|² 근사 (차분)"""
    s_center = dedekind_lambda(lf, center + delta_h, t)
    s_plus = dedekind_lambda(lf, center + delta_h, t + delta_h)
    s_minus = dedekind_lambda(lf, center + delta_h, t - delta_h)
    if s_center is None or s_plus is None or s_minus is None:
        return 0.0
    if abs(s_center) < 1e-300:
        return float('inf')
    deriv = (s_plus - s_minus) / (2 * delta_h * 1j)
    conn = deriv / s_center
    return abs(conn) ** 2


def generate_fp_points(all_zeros, t_min, t_max, n_random=30):
    """
    FP 후보 생성: 인접 영점 중점 + 저t 영역 + 랜덤 비영점.
    all_zeros는 T_ZEROS_MAX까지의 전체 영점 (FP가 영점 근방에 오지 않도록).
    min_dist=0.2: r=0.1 contour 안전 거리.
    """
    fps = []
    zeros_arr = np.array(all_zeros)

    # T_MAX 이내 영점만 중점 계산에 사용
    tp_zeros = [z for z in all_zeros if z <= t_max]

    # 1. 인접 영점 중점 (갭 > 0.4이면 포함)
    if len(tp_zeros) >= 2:
        for i in range(len(tp_zeros) - 1):
            mid = (tp_zeros[i] + tp_zeros[i + 1]) / 2.0
            gap = tp_zeros[i + 1] - tp_zeros[i]
            if gap > 0.4:
                fps.append(mid)

    # 2. 저t 영역 (첫 영점 이전)
    if len(all_zeros) > 0:
        first_zero = all_zeros[0]
        if first_zero > 1.0:
            low_t_max = first_zero - 0.2
            if low_t_max > 1.0:
                n_low = min(5, int((low_t_max - 0.5) / 0.3))
                for i in range(n_low):
                    t = 0.5 + (low_t_max - 0.5) * (i + 0.5) / n_low
                    fps.append(t)

    # 3. 랜덤: min_dist=0.2 (전체 all_zeros 기준)
    min_dist = 0.2
    for seed in [42, 7, 123, 314, 2024, 999, 55, 77, 1234, 5678]:
        rng = np.random.RandomState(seed)
        added = 0
        attempts = 0
        target = n_random // 10 + 2
        while added < target and attempts < 500:
            t_rand = rng.uniform(t_min, t_max)
            if len(all_zeros) > 0:
                dists = np.abs(zeros_arr - t_rand)
                if dists.min() < min_dist:
                    attempts += 1
                    continue
            if fps and min(abs(t_rand - f) for f in fps) < 0.15:
                attempts += 1
                continue
            fps.append(t_rand)
            added += 1
            attempts += 1

    return np.array(sorted(fps))


def main():
    START = time.time()
    primary_r = MONO_RADII[PRIMARY_RADIUS_IDX]

    print("=" * 70)
    print("[C-415] Dedekind ζ_K (S₃, x³-x-1) 모노드로미 TP/FP 검증")
    print(f"  source: {INFO['source']}")
    print(f"  polynomial: {INFO['polynomial']}")
    print(f"  discriminant: {INFO['discriminant']}")
    print(f"  degree: {INFO['degree']}")
    print(f"  σ_c: {INFO['sigma_c']}")
    print(f"  T_MAX: {T_MAX}, T_ZEROS_MAX: {T_ZEROS_MAX}")
    print(f"  n_steps: {MONO_NSTEPS}, radii: {MONO_RADII}")
    print(f"  주 판정 반지름: r={primary_r}")
    print(f"  핵심: Conjecture 1 예측-검증. sym^k와 독립된 구성 경로.")
    print(f"        ζ_K = ζ·L(s,ρ₂). 비가환 S₃ Galois 군 기원.")
    print("=" * 70)

    # ─── 1. L-함수 초기화 ───
    print("\n[1단계] Dedekind zeta 초기화...")
    lf, nf = init_dedekind_lfun()

    # FE 확인
    feq = pari.lfuncheckfeq(lf)
    print(f"  lfuncheckfeq: {feq} (음수 = 일치하는 자릿수)")
    feq_val = float(feq)
    if feq_val > -40:
        print(f"  ⚠️ lfuncheckfeq = {feq_val} (>-40). 정밀도 부족 가능성.")

    # ─── 2. 영점 수집 ───
    print(f"\n[2단계] 영점 수집 (T ≤ {T_ZEROS_MAX}, TP용 T ≤ {T_MAX})...")
    zeros_raw = pari.lfunzeros(lf, T_ZEROS_MAX)
    all_zeros = [float(z) for z in zeros_raw if float(z) > 1.0]
    true_zeros = [z for z in all_zeros if z <= T_MAX]
    print(f"  전체 영점: {len(all_zeros)} (T≤{T_ZEROS_MAX})")
    print(f"  TP용 영점: {len(true_zeros)} (T≤{T_MAX})")

    if len(true_zeros) < 15:
        print(f"  ⚠️ TP 영점 부족 ({len(true_zeros)} < 15). T_MAX 확장...")
        T_MAX_NEW = 100.0
        zeros_raw = pari.lfunzeros(lf, T_MAX_NEW + 20)
        all_zeros = [float(z) for z in zeros_raw if float(z) > 1.0]
        true_zeros = [z for z in all_zeros if z <= T_MAX_NEW]
        print(f"  확장 후: 전체 {len(all_zeros)}, TP용 {len(true_zeros)}")

    # TP를 최대 20개로 제한
    tp_zeros = true_zeros[:20]

    print(f"  사용할 TP 영점: {len(tp_zeros)}개")
    print(f"  영점 목록 (처음 10개):")
    for i, z in enumerate(tp_zeros[:10]):
        print(f"    [{i+1}] t = {z:.6f}")
    if len(tp_zeros) > 10:
        print(f"    ... ({len(tp_zeros)}개 총)")

    # ─── 2.5. 영점 기원 식별 ───
    # ζ_K = ζ·L(s,ρ₂) 이므로 일부 영점은 ζ(s)에서, 일부는 L(s,ρ₂)에서 유래
    print(f"\n[2.5단계] 영점 기원 식별 (ζ vs L(s,ρ₂))...")
    zeta_zeros_raw = pari('lfunzeros(1, %s)' % T_ZEROS_MAX)
    zeta_zeros = [float(z) for z in zeta_zeros_raw if float(z) > 1.0]
    print(f"  ζ(s) 영점 (T≤{T_ZEROS_MAX}): {len(zeta_zeros)}개")

    # 각 ζ_K 영점이 ζ 영점과 매칭되는지 확인
    origin_map = {}
    for z in all_zeros:
        dists = [abs(z - zz) for zz in zeta_zeros] if zeta_zeros else [999]
        min_dist = min(dists)
        if min_dist < 0.01:
            origin_map[z] = 'ζ'
        else:
            origin_map[z] = 'L(ρ₂)'

    n_from_zeta = sum(1 for v in origin_map.values() if v == 'ζ')
    n_from_artin = sum(1 for v in origin_map.values() if v == 'L(ρ₂)')
    print(f"  ζ(s) 기원: {n_from_zeta}개, L(s,ρ₂) 기원: {n_from_artin}개")

    for z in tp_zeros[:10]:
        print(f"    t={z:.4f}: {origin_map.get(z, '?')}")

    # ─── 3. 건전성 검사 ───
    sigma_c = INFO['sigma_c']
    print(f"\n[3단계] Λ 건전성 검사 (σ_c = {sigma_c})...")

    print("  영점에서 |Λ|:")
    for z in tp_zeros[:5]:
        val = dedekind_lambda(lf, sigma_c, z)
        if val is not None:
            print(f"    Λ({sigma_c}+i·{z:.4f}) = {val.real:.2e} + {val.imag:.2e}i  |Λ| = {abs(val):.2e}")

    # FE 대칭: Λ(s) = ε·Λ(1-s)
    print(f"\n  FE 대칭 검증:")
    for t_test in [15.0, 25.0, 35.0]:
        s_val = dedekind_lambda(lf, sigma_c, t_test)
        s_conj = dedekind_lambda(lf, 1.0 - sigma_c, t_test)
        if s_val is not None and s_conj is not None:
            ratio = abs(s_val) / abs(s_conj) if abs(s_conj) > 1e-300 else float('inf')
            print(f"    |Λ({sigma_c}+{t_test}i)| = {abs(s_val):.6e}, "
                  f"|Λ({1-sigma_c}+{t_test}i)| = {abs(s_conj):.6e}, ratio = {ratio:.6f}")

    # ─── 4. FP 생성 ───
    print(f"\n[4단계] FP 후보 생성...")
    fp_points = generate_fp_points(all_zeros, 2.0, T_MAX)
    print(f"  FP 후보: {len(fp_points)}개")
    for i, fp in enumerate(fp_points[:5]):
        print(f"    FP[{i+1}]: t = {fp:.4f}")
    if len(fp_points) > 5:
        print(f"    ... ({len(fp_points)}개 총)")

    # ─── 5. 모노드로미 + 곡률 측정 ───
    print(f"\n[5단계] 모노드로미 + 곡률 측정...")

    tp_results = []
    print(f"\n  --- TP ({len(tp_zeros)}개) ---")
    for i, tz in enumerate(tp_zeros):
        monos = {}
        for ri, r in enumerate(MONO_RADII):
            mono = monodromy_contour(lf, tz, radius=r, n_steps=MONO_NSTEPS, center=sigma_c)
            monos[r] = mono
        kappa = curvature(lf, tz, center=sigma_c, delta_h=0.01)
        origin = origin_map.get(tz, '?')
        tp_results.append({'t': tz, 'monos': monos, 'kappa': kappa, 'origin': origin})

        mono_primary = monos[primary_r] / np.pi
        mono_strs = ', '.join(f"r={r}:{monos[r]/np.pi:.4f}π" for r in MONO_RADII)
        print(f"  TP[{i+1:2d}] t={tz:8.4f} ({origin:5s}): mono/π = {mono_primary:.4f}  κ = {kappa:.2e}  [{mono_strs}]")

    fp_results = []
    print(f"\n  --- FP ({len(fp_points)}개) ---")
    for i, fp in enumerate(fp_points):
        monos = {}
        for ri, r in enumerate(MONO_RADII):
            mono = monodromy_contour(lf, fp, radius=r, n_steps=MONO_NSTEPS, center=sigma_c)
            monos[r] = mono
        kappa = curvature(lf, fp, center=sigma_c, delta_h=0.01)
        fp_results.append({'t': fp, 'monos': monos, 'kappa': kappa})

        mono_primary = monos[primary_r] / np.pi
        mono_strs = ', '.join(f"r={r}:{monos[r]/np.pi:.4f}π" for r in MONO_RADII)
        print(f"  FP[{i+1:2d}] t={fp:8.4f}: mono/π = {mono_primary:.4f}  κ = {kappa:.2e}  [{mono_strs}]")

    # ─── 6. 통계 분석 ───
    print(f"\n{'=' * 70}")
    print(f"[6단계] 통계 분석")
    print(f"{'=' * 70}")

    # 결과 파일 열기
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, 'w') as f:
        def out(msg=''):
            print(msg)
            f.write(msg + '\n')

        out("=" * 72)
        out(f"[C-415] Dedekind ζ_K (S₃, x³-x-1) 모노드로미 TP/FP 검증")
        out("=" * 72)
        out(f"수체: K = Q(α), α³-α-1=0")
        out(f"Galois 군: S₃ (비가환)")
        out(f"판별식: {INFO['discriminant']}")
        out(f"conductor: |disc| = 23")
        out(f"degree: {INFO['degree']}")
        out(f"gammaV: [0, 0, 1]")
        out(f"임계선: σ = {sigma_c}")
        out(f"lfuncheckfeq: {feq}")
        out(f"T_MAX: {T_MAX}, T_ZEROS_MAX: {T_ZEROS_MAX}")
        out(f"n_steps: {MONO_NSTEPS}, radii: {MONO_RADII}")
        out(f"주 판정 반지름: r={primary_r}")
        out(f"ζ_K = ζ(s)·L(s,ρ₂) 분해: {n_from_zeta} from ζ, {n_from_artin} from L(ρ₂)")
        out()

        # 반지름별 분석
        for ri, r in enumerate(MONO_RADII):
            out(f"\n{'─' * 60}")
            out(f"반지름 r = {r}")
            out(f"{'─' * 60}")

            tp_monos = np.array([res['monos'][r] / np.pi for res in tp_results])
            fp_monos = np.array([res['monos'][r] / np.pi for res in fp_results])
            tp_kappas = np.array([res['kappa'] for res in tp_results])
            fp_kappas = np.array([res['kappa'] for res in fp_results])

            out(f"\n  TP (n={len(tp_monos)}):")
            out(f"    mono/π: mean = {np.mean(tp_monos):.6f} ± {np.std(tp_monos):.6f}")
            out(f"    mono/π: min = {np.min(tp_monos):.6f}, max = {np.max(tp_monos):.6f}")
            out(f"    κ: mean = {np.mean(tp_kappas):.2e} ± {np.std(tp_kappas):.2e}")
            out(f"    κ: min = {np.min(tp_kappas):.2e}, max = {np.max(tp_kappas):.2e}")

            out(f"\n  FP (n={len(fp_monos)}):")
            out(f"    mono/π: mean = {np.mean(fp_monos):.6f} ± {np.std(fp_monos):.6f}")
            out(f"    mono/π: min = {np.min(fp_monos):.6f}, max = {np.max(fp_monos):.6f}")
            out(f"    κ: mean = {np.mean(fp_kappas):.2e} ± {np.std(fp_kappas):.2e}")
            out(f"    κ: min = {np.min(fp_kappas):.2e}, max = {np.max(fp_kappas):.2e}")

            # KS 검정
            if len(tp_monos) > 0 and len(fp_monos) > 0:
                ks_stat, ks_p = stats.ks_2samp(tp_monos, fp_monos)
                out(f"\n  KS 검정: D = {ks_stat:.6f}, p = {ks_p:.2e}")
            else:
                ks_p = 1.0
                out(f"\n  KS 검정: 데이터 부족")

            # 이중기준: mono > 1.0 AND κ > median(κ_all) → TP 예측
            all_kappas = np.concatenate([tp_kappas, fp_kappas])
            kappa_threshold = np.median(all_kappas) * 2
            tp_dual = sum(1 for res in tp_results if abs(res['monos'][r]/np.pi - 2.0) < 0.5 and res['kappa'] > kappa_threshold)
            fp_dual = sum(1 for res in fp_results if abs(res['monos'][r]/np.pi - 2.0) < 0.5 and res['kappa'] > kappa_threshold)
            dual_prec = tp_dual / (tp_dual + fp_dual) * 100 if (tp_dual + fp_dual) > 0 else 0
            dual_recall = tp_dual / len(tp_results) * 100 if len(tp_results) > 0 else 0
            out(f"\n  이중기준 (mono≈2π + κ>{kappa_threshold:.1f}):")
            out(f"    TP 포착: {tp_dual}/{len(tp_results)}")
            out(f"    FP 오판: {fp_dual}/{len(fp_results)}")
            out(f"    Precision: {dual_prec:.1f}%")
            out(f"    Recall: {dual_recall:.1f}%")

            # 주 판정 기준에서 최종 판정
            if ri == PRIMARY_RADIUS_IDX:
                pass_mono_tp = np.all(np.abs(tp_monos - 2.0) < 0.01)
                pass_mono_fp = np.all(np.abs(fp_monos) < 0.01)
                pass_ks = ks_p < 1e-6
                pass_dual = (dual_prec == 100.0) and (dual_recall == 100.0)

        # 기원별 분석
        out(f"\n{'─' * 60}")
        out(f"영점 기원별 분석 (주 반지름 r={primary_r})")
        out(f"{'─' * 60}")
        for origin_type in ['ζ', 'L(ρ₂)']:
            origin_results = [r for r in tp_results if r['origin'] == origin_type]
            if origin_results:
                monos = [r['monos'][primary_r] / np.pi for r in origin_results]
                kappas = [r['kappa'] for r in origin_results]
                out(f"\n  {origin_type} 기원 (n={len(origin_results)}):")
                out(f"    mono/π: {np.mean(monos):.6f} ± {np.std(monos):.6f}")
                out(f"    κ: mean={np.mean(kappas):.2e}, range=[{np.min(kappas):.2e}, {np.max(kappas):.2e}]")

        # 개별 TP 상세
        out(f"\n{'─' * 60}")
        out(f"TP 상세 (r={primary_r})")
        out(f"{'─' * 60}")
        for i, res in enumerate(tp_results):
            mono_p = res['monos'][primary_r] / np.pi
            out(f"  [{i+1:2d}] t={res['t']:8.4f} ({res['origin']:5s}) mono/π={mono_p:.6f}  κ={res['kappa']:.2e}")

        out(f"\n{'─' * 60}")
        out(f"FP 상세 (r={primary_r})")
        out(f"{'─' * 60}")
        for i, res in enumerate(fp_results):
            mono_p = res['monos'][primary_r] / np.pi
            out(f"  [{i+1:2d}] t={res['t']:8.4f} mono/π={mono_p:.6f}  κ={res['kappa']:.2e}")

        # 최종 판정
        out(f"\n{'=' * 72}")
        out(f"최종 판정 (r={primary_r} 기준)")
        out(f"{'=' * 72}")

        tp_monos_p = np.array([res['monos'][primary_r] / np.pi for res in tp_results])
        fp_monos_p = np.array([res['monos'][primary_r] / np.pi for res in fp_results])

        pass1 = len(tp_results) >= 15 and np.all(np.abs(tp_monos_p - 2.0) < 0.01)
        pass2 = len(fp_results) >= 20 and np.all(np.abs(fp_monos_p) < 0.01)
        ks_stat, ks_p = stats.ks_2samp(tp_monos_p, fp_monos_p)
        pass3 = ks_p < 1e-6

        tp_kappas = np.array([res['kappa'] for res in tp_results])
        fp_kappas = np.array([res['kappa'] for res in fp_results])
        all_k = np.concatenate([tp_kappas, fp_kappas])
        k_thresh = np.median(all_k) * 2
        tp_correct = sum(1 for r in tp_results if abs(r['monos'][primary_r]/np.pi - 2.0) < 0.5 and r['kappa'] > k_thresh)
        fp_correct = sum(1 for r in fp_results if abs(r['monos'][primary_r]/np.pi - 2.0) < 0.5 and r['kappa'] > k_thresh)
        pass4 = (tp_correct == len(tp_results)) and (fp_correct == 0)

        out(f"  SC1: TP ≥ 15, mono/π ∈ [1.99, 2.01]  → {'✅ PASS' if pass1 else '❌ FAIL'}")
        out(f"       TP = {len(tp_results)}, mono/π = {np.mean(tp_monos_p):.6f} ± {np.std(tp_monos_p):.6f}")
        out(f"  SC2: FP ≥ 20, mono/π ∈ [-0.01, 0.01] → {'✅ PASS' if pass2 else '❌ FAIL'}")
        out(f"       FP = {len(fp_results)}, mono/π = {np.mean(fp_monos_p):.6f} ± {np.std(fp_monos_p):.6f}")
        out(f"  SC3: KS p < 1e-6                      → {'✅ PASS' if pass3 else '❌ FAIL'}")
        out(f"       KS D = {ks_stat:.6f}, p = {ks_p:.2e}")
        out(f"  SC4: 이중기준 100% 분리              → {'✅ PASS' if pass4 else '❌ FAIL'}")
        out(f"       TP 포착: {tp_correct}/{len(tp_results)}, FP 오판: {fp_correct}/{len(fp_results)}")

        n_pass = sum([pass1, pass2, pass3, pass4])
        if n_pass == 4:
            verdict = "강한 양성 (4/4)"
        elif n_pass >= 3:
            verdict = f"양성 ({n_pass}/4)"
        elif n_pass >= 2:
            verdict = f"조건부 양성 ({n_pass}/4)"
        else:
            verdict = f"음성 ({n_pass}/4)"

        out(f"\n  *** 판정: {verdict} ***")
        out(f"\n  Conjecture 1 예측-검증: {'PASS' if n_pass == 4 else 'CONDITIONAL' if n_pass >= 2 else 'FAIL'}")
        out(f"  구성 경로: Dedekind zeta (sym^k 타워와 독립)")
        out(f"  Galois 군: S₃ (비가환 — 기존 검증은 모두 가환)")
        out(f"  ζ_K 분해: ζ(s)·L(s,ρ₂). 양쪽 기원 영점 모두 검증.")

        # 3반지름 비교 요약
        out(f"\n{'─' * 60}")
        out(f"3반지름 비교 (TP mono/π 평균 ± std)")
        out(f"{'─' * 60}")
        for r in MONO_RADII:
            tp_m = np.array([res['monos'][r] / np.pi for res in tp_results])
            fp_m = np.array([res['monos'][r] / np.pi for res in fp_results])
            out(f"  r={r}: TP {np.mean(tp_m):.6f}±{np.std(tp_m):.6f}, FP {np.mean(fp_m):.6f}±{np.std(fp_m):.6f}")

        # κ 분리 상세
        out(f"\n{'─' * 60}")
        out(f"κ 분리 상세")
        out(f"{'─' * 60}")
        if len(tp_kappas) > 0 and len(fp_kappas) > 0:
            out(f"  κ_TP: mean={np.mean(tp_kappas):.2e}, min={np.min(tp_kappas):.2e}")
            out(f"  κ_FP: mean={np.mean(fp_kappas):.2e}, max={np.max(fp_kappas):.2e}")
            if np.max(fp_kappas) > 0:
                ratio = np.min(tp_kappas) / np.max(fp_kappas)
                out(f"  min(κ_TP)/max(κ_FP) = {ratio:.1f}×")

        elapsed = time.time() - START
        out(f"\n  실행 시간: {elapsed:.1f}초 ({elapsed/60:.1f}분)")


if __name__ == '__main__':
    main()
