"""
=============================================================================
[Project RDL] False Positive 해부 실험
=============================================================================
수학자 질문: "FP(false positive)의 수학적 정체는 무엇인가?"

3가지 가설 검증:
  (A) 유사-영점: |F₂| 풍경에 영점이 아닌 극소점이 존재 → FP 위치에서 |F₂|'' > 0 (U자)
  (B) 해상도 한계: 인접 영점 쌍의 중간점이 FP → FP→true zero 거리 ~ 해상도
  (C) 무작위 잡음: 시드 간 FP 위치 비일관 → 낮은 겹침률

출력:
  1. 각 FP에서 |F₂|''(t) 부호 분포 (극소 vs 변곡점)
  2. FP→nearest true zero 거리 히스토그램
  3. GUE pair spacing과 FP-zero 간격 비교
  4. 5시드에서 FP 위치 겹침률

설정: t∈[14,50], K=128, 시드 5개, Baseline (L_tgt)
"""

import sys, os, time
import numpy as np
import torch
import torch.nn as nn

os.environ.setdefault("OMP_NUM_THREADS", "10")
os.environ.setdefault("MKL_NUM_THREADS", "10")
torch.set_num_threads(10)

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
from gdl.rdl.constants import PrecisionManager
from gdl.rdl.pipeline.xi_feature_dataset import (
    get_or_build_cache, compute_zeros_in_range, XiFeatureDataset
)
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss

PrecisionManager.setup_precision()

# ─── 설정 ───
T_MIN, T_MAX = 14.0, 50.0
N_POINTS = 1000
HIDDEN = 64
EPOCHS = 100
LR = 1e-3
BATCH = 32
IN_FEATURES = 128
SEEDS = [42, 7, 123, 314, 2024]


def train_baseline(dataset, seed):
    """표준 TotalResonanceLoss Baseline 훈련"""
    torch.manual_seed(seed); np.random.seed(seed)
    val_size = int(len(dataset) * 0.2)
    train_size = len(dataset) - val_size
    train_ds, val_ds = torch.utils.data.random_split(
        dataset, [train_size, val_size],
        generator=torch.Generator().manual_seed(seed)
    )
    train_loader = torch.utils.data.DataLoader(train_ds, batch_size=BATCH, shuffle=True, drop_last=True)
    val_loader = torch.utils.data.DataLoader(val_ds, batch_size=BATCH, shuffle=False)

    model = MasterResonantNetwork(
        in_features=dataset.in_features, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5,
        pqo_mode='cos2'
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)

    best_val = float('inf')
    t_start = time.time()
    for ep in range(EPOCHS):
        model.train()
        ep_loss = 0.0; ep_n = 0
        for X_batch, _ in train_loader:
            X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            optimizer.zero_grad(set_to_none=True)
            outputs = model(X_in)
            total_loss, _ = loss_fn(**outputs)
            if torch.isnan(total_loss) or torch.isinf(total_loss):
                continue
            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            ep_loss += total_loss.item(); ep_n += 1

        if (ep + 1) % 25 == 0 or ep == EPOCHS - 1:
            model.eval()
            va = 0.0; nv = 0
            with torch.enable_grad():
                for X_batch, _ in val_loader:
                    X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
                    X_in.requires_grad_(True)
                    o = model(X_in)
                    vl, _ = loss_fn(**o)
                    if not (torch.isnan(vl) or torch.isinf(vl)):
                        va += vl.item(); nv += 1
            vl_m = va / max(1, nv)
            best_val = min(best_val, vl_m)
            elapsed = time.time() - t_start
            print(f"    ep {ep+1}/{EPOCHS}: train={ep_loss/max(1,ep_n):.5f}, val={vl_m:.5f}, {elapsed:.0f}s", flush=True)

    return model, best_val


def compute_f2_landscape(model, dataset):
    """전체 데이터에서 |F₂| 풍경 계산 (phi/psi/L_G 기반)"""
    model.eval()
    loader = torch.utils.data.DataLoader(dataset, batch_size=64, shuffle=False)
    f2_vals = []
    with torch.enable_grad():
        for X_batch, _ in loader:
            X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            out = model(X_in)
            phi = out["phi"].detach()
            psi = out["psi"].detach()
            L_G = out["L_G"].detach()
            phi_real = phi.to(dtype=PrecisionManager.REAL_DTYPE)
            rot = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
            psi_c = psi.to(dtype=PrecisionManager.COMPLEX_DTYPE)
            f2_batch = (rot * (L_G - psi_c)).imag.mean(dim=-1).cpu().numpy()
            f2_vals.append(f2_batch)
    return np.abs(np.concatenate(f2_vals))


def find_fps_and_fns(f2_arr, is_zero, threshold):
    """FP/FN 분류"""
    predicted_zero = f2_arr < threshold
    is_z = is_zero.numpy() if hasattr(is_zero, 'numpy') else is_zero

    fp_mask = predicted_zero & (~is_z)  # 영점으로 예측했지만 실제로는 아닌 것
    fn_mask = (~predicted_zero) & is_z  # 영점이지만 놓친 것
    tp_mask = predicted_zero & is_z     # 올바르게 검출

    return fp_mask, fn_mask, tp_mask


def compute_f2_curvature_at_points(model, dataset, t_indices, delta=0.05, n_pts=21):
    """특정 t 위치 주변에서 |F₂|의 2차 미분(곡률) 계산

    각 FP 위치 주변 [-delta, +delta] 구간에서 n_pts개 점을 평가하고
    중심에서의 2차 미분 부호를 반환.
    양수 = U자 (극소점, 가설 A)
    음수 = V자 또는 변곡점
    """
    model.eval()
    curvatures = []
    f2_profiles = []

    for idx in t_indices:
        t_center = dataset.t[idx].item()
        t_local = torch.linspace(t_center - delta, t_center + delta, n_pts,
                                  dtype=PrecisionManager.REAL_DTYPE)

        # 각 t에 대해 특징 벡터 생성 + 모델 평가
        features_list = []
        for t_val in t_local:
            feat = dataset._build_features(t_val)
            features_list.append(feat)
        X = torch.stack(features_list)

        with torch.enable_grad():
            X_in = X.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            out = model(X_in)
            phi = out["phi"].detach()
            psi = out["psi"].detach()
            L_G = out["L_G"].detach()
            phi_real = phi.to(dtype=PrecisionManager.REAL_DTYPE)
            rot = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
            psi_c = psi.to(dtype=PrecisionManager.COMPLEX_DTYPE)
            f2_local = torch.abs((rot * (L_G - psi_c)).imag.mean(dim=-1)).cpu().numpy()

        f2_profiles.append((t_local.numpy(), f2_local))

        # 2차 미분: 중앙 차분 (center = n_pts//2)
        dt = (2 * delta) / (n_pts - 1)
        center = n_pts // 2
        if center > 0 and center < len(f2_local) - 1:
            d2f = (f2_local[center + 1] - 2 * f2_local[center] + f2_local[center - 1]) / (dt ** 2)
        else:
            d2f = np.nan
        curvatures.append(d2f)

    return np.array(curvatures), f2_profiles


def compute_gue_pair_spacing(zeros_list):
    """GUE 쌍 간격 분포 (정규화된 영점 간격)"""
    if len(zeros_list) < 2:
        return np.array([])
    spacings = np.diff(zeros_list)
    mean_spacing = np.mean(spacings)
    return spacings / mean_spacing  # 정규화


def main():
    out = []
    def log(msg):
        print(msg, flush=True); out.append(msg)

    log("=" * 72)
    log("  FP (False Positive) 해부 실험")
    log("=" * 72)
    log(f"  구간: t∈[{T_MIN},{T_MAX}], K={IN_FEATURES}, 시드={SEEDS}")
    log(f"  질문: FP의 수학적 정체는 무엇인가?")
    log(f"  가설 A: 유사-영점 (|F₂|'' > 0, U자 극소)")
    log(f"  가설 B: 해상도 한계 (FP→zero 거리 ~ resolution)")
    log(f"  가설 C: 무작위 잡음 (시드 간 FP 위치 비일관)")
    log("")

    start_time = time.time()

    # ─── 데이터 준비 ───
    zeros_list = compute_zeros_in_range(T_MIN, T_MAX, dps=25)
    n_zeros = len(zeros_list)
    zeros_arr = np.array([float(z) for z in zeros_list])
    log(f"  영점 수: {n_zeros}, 평균 간격: {(T_MAX-T_MIN)/n_zeros:.3f}")

    cache_path = os.path.expanduser(
        f'~/Desktop/gdl_unified/outputs/xi_cache_t{T_MIN}-{T_MAX}_n{N_POINTS}.pt'
    )
    cache_data = get_or_build_cache(cache_path, T_MIN, T_MAX, N_POINTS,
                                     mp_dps=50, zeros_list=zeros_list)
    ds = XiFeatureDataset(cache_data, in_features=IN_FEATURES)

    t_arr = ds.t.numpy()
    is_zero = ds.is_near_zero.numpy()
    n_true_zeros = int(is_zero.sum())
    log(f"  데이터 점: {len(ds)}, 영점 마스크 양성: {n_true_zeros}")

    # 정보론적 해상도 하한
    resolution = (T_MAX - T_MIN) / (2 * IN_FEATURES)
    log(f"  해상도 하한 (T_range/2K): {resolution:.4f}")

    # GUE pair spacing
    gue_spacings = compute_gue_pair_spacing(zeros_arr)
    if len(gue_spacings) > 0:
        log(f"  GUE 정규화 간격: mean={gue_spacings.mean():.4f}, std={gue_spacings.std():.4f}")

    # ─── 시드별 훈련 + FP 분석 ───
    all_fp_indices = []  # 시드별 FP 인덱스 집합
    all_fp_curvatures = []
    all_fp_distances = []
    all_fp_counts = []
    all_tp_counts = []
    per_seed_fp_curvature = []  # [{idx: curv}, ...]  구조적 FP 분리용
    per_seed_fp_distance = []   # [{idx: dist}, ...]

    for seed in SEEDS:
        log(f"\n  --- seed={seed} ---")

        # 훈련
        log(f"    [훈련] Baseline (L_tgt)...")
        model, val_loss = train_baseline(ds, seed)

        # |F₂| 풍경 계산
        f2_arr = compute_f2_landscape(model, ds)
        threshold = np.median(f2_arr) * 0.1

        # FP/TP 분류
        fp_mask, fn_mask, tp_mask = find_fps_and_fns(f2_arr, is_zero, threshold)
        n_fp = int(fp_mask.sum())
        n_tp = int(tp_mask.sum())
        n_fn = int(fn_mask.sum())
        n_pred = int((f2_arr < threshold).sum())
        precision = n_tp / max(n_pred, 1)
        recall = n_tp / max(n_true_zeros, 1)

        log(f"    val={val_loss:.5f}, threshold={threshold:.6f}")
        log(f"    TP={n_tp}, FP={n_fp}, FN={n_fn}, precision={precision:.3f}, recall={recall:.3f}")

        fp_indices = np.where(fp_mask)[0]
        all_fp_indices.append(set(fp_indices.tolist()))
        all_fp_counts.append(n_fp)
        all_tp_counts.append(n_tp)

        if n_fp == 0:
            log(f"    FP=0, 분석 생략")
            per_seed_fp_curvature.append({})
            per_seed_fp_distance.append({})
            continue

        # ─── 가설 A: 곡률 분석 ───
        log(f"    [가설 A] FP 위치 곡률 분석 ({n_fp}개)...")
        curvatures, profiles = compute_f2_curvature_at_points(model, ds, fp_indices)
        n_positive = int(np.sum(curvatures > 0))  # U자 (극소점)
        n_negative = int(np.sum(curvatures < 0))  # V자/변곡점
        n_nan = int(np.sum(np.isnan(curvatures)))
        all_fp_curvatures.extend(curvatures[~np.isnan(curvatures)].tolist())
        # 인덱스별 곡률 저장 (구조적 FP 분리용)
        curv_dict = {}
        for k, idx in enumerate(fp_indices):
            if not np.isnan(curvatures[k]):
                curv_dict[idx] = curvatures[k]
        per_seed_fp_curvature.append(curv_dict)
        log(f"    |F₂|'' > 0 (U자/극소): {n_positive}/{n_fp} ({100*n_positive/n_fp:.1f}%)")
        log(f"    |F₂|'' < 0 (V자/변곡): {n_negative}/{n_fp} ({100*n_negative/n_fp:.1f}%)")

        # ─── 가설 B: FP→true zero 거리 ───
        fp_t = t_arr[fp_indices]
        distances = []
        for t_fp in fp_t:
            dists = np.abs(zeros_arr - t_fp)
            distances.append(dists.min())
        distances = np.array(distances)
        all_fp_distances.extend(distances.tolist())
        # 인덱스별 거리 저장 (구조적 FP 분리용)
        dist_dict = {idx: distances[k] for k, idx in enumerate(fp_indices)}
        per_seed_fp_distance.append(dist_dict)

        near_res = np.sum(distances < resolution)
        near_2res = np.sum(distances < 2 * resolution)
        log(f"    [가설 B] FP→nearest zero 거리:")
        log(f"      mean={distances.mean():.4f}, median={np.median(distances):.4f}")
        log(f"      < resolution ({resolution:.4f}): {near_res}/{n_fp} ({100*near_res/n_fp:.1f}%)")
        log(f"      < 2×resolution ({2*resolution:.4f}): {near_2res}/{n_fp} ({100*near_2res/n_fp:.1f}%)")

    # ─── 가설 C: 시드 간 FP 겹침률 ───
    log(f"\n{'='*72}")
    log("  가설 C: 시드 간 FP 위치 일관성")
    log(f"{'='*72}")

    # 인접 인덱스도 겹침으로 간주 (±4 tolerance, 해상도 하한 ~0.14에 대응하는 인덱스 수)
    def expand_set(s, tol=4):
        expanded = set()
        for idx in s:
            for d in range(-tol, tol + 1):
                expanded.add(idx + d)
        return expanded

    if len(all_fp_indices) >= 2:
        # 모든 쌍에 대해 Jaccard 유사도
        pair_overlaps = []
        for i in range(len(all_fp_indices)):
            for j in range(i + 1, len(all_fp_indices)):
                s1 = expand_set(all_fp_indices[i])
                s2 = expand_set(all_fp_indices[j])
                if len(s1 | s2) > 0:
                    jaccard = len(s1 & s2) / len(s1 | s2)
                    pair_overlaps.append(jaccard)

        if pair_overlaps:
            log(f"  Jaccard 유사도 (±4 tolerance): mean={np.mean(pair_overlaps):.4f}, std={np.std(pair_overlaps):.4f}")
            log(f"  개별 쌍: {['%.3f' % j for j in pair_overlaps]}")

        # 3시드 이상 합의: 각 FP 위치에 대해 ±4 이내에 FP가 있는 시드 수
        all_fp_union = set()
        for s in all_fp_indices:
            all_fp_union.update(s)
        consensus_counts = {}  # idx -> 동의 시드 수
        for idx in all_fp_union:
            n_agree = sum(1 for s in all_fp_indices if any(abs(idx - j) <= 4 for j in s))
            consensus_counts[idx] = n_agree
        # 겹치는 인접 위치 제거 (±4 내 최고 합의만 유지)
        visited = set()
        consensus_unique = []
        for idx in sorted(consensus_counts, key=lambda x: -consensus_counts[x]):
            if idx in visited:
                continue
            consensus_unique.append((idx, consensus_counts[idx]))
            for d in range(-4, 5):
                visited.add(idx + d)
        consensus_3 = sum(1 for _, c in consensus_unique if c >= 3)
        consensus_4 = sum(1 for _, c in consensus_unique if c >= 4)
        consensus_5 = sum(1 for _, c in consensus_unique if c >= 5)
        log(f"  ≥3시드 합의 FP 위치: {consensus_3}")
        log(f"  ≥4시드 합의 FP 위치: {consensus_4}")
        log(f"  ≥5시드 합의 FP 위치 (전원): {consensus_5}")

        # ─── 구조적 FP 분리 분석 (수학자 권장) ───
        structural_fp_idx = [idx for idx, c in consensus_unique if c >= 3]
        if structural_fp_idx:
            log(f"\n  [구조적 FP 분리] ≥3시드 합의 위치의 곡률/거리:")
            struct_curvs = []
            struct_dists = []
            for idx in structural_fp_idx:
                for sd in range(len(per_seed_fp_curvature)):
                    # ±4 tolerance 내에서 해당 시드의 FP 매칭
                    for fp_idx in per_seed_fp_curvature[sd]:
                        if abs(fp_idx - idx) <= 4:
                            struct_curvs.append(per_seed_fp_curvature[sd][fp_idx])
                            if fp_idx in per_seed_fp_distance[sd]:
                                struct_dists.append(per_seed_fp_distance[sd][fp_idx])
                            break
            if struct_curvs:
                sc = np.array(struct_curvs)
                frac_u = np.mean(sc > 0)
                log(f"    곡률 U자(극소) 비율: {frac_u:.3f} ({len(sc)}개)")
            if struct_dists:
                sd_arr = np.array(struct_dists)
                log(f"    →zero 거리: mean={sd_arr.mean():.4f}, median={np.median(sd_arr):.4f}")
                log(f"    < resolution: {np.mean(sd_arr < resolution):.3f}")
                if struct_curvs:
                    if frac_u > 0.7 and np.mean(sd_arr < 2*resolution) > 0.5:
                        log(f"    → 구조적 FP = 영점 근처 유사-영점 (가설 A+B 결합)")
                    elif frac_u > 0.7:
                        log(f"    → 구조적 FP = 유사-영점 (가설 A 지배)")
                    elif np.mean(sd_arr < 2*resolution) > 0.5:
                        log(f"    → 구조적 FP = 해상도 한계 (가설 B 지배)")
                    else:
                        log(f"    → 구조적 FP 정체 불명확")

    # ─── 종합 분석 ───
    log(f"\n{'='*72}")
    log("  종합 분석")
    log(f"{'='*72}")

    log(f"\n  앙상블 통계 (5시드):")
    log(f"  FP 수: {np.mean(all_fp_counts):.1f}±{np.std(all_fp_counts):.1f}")
    log(f"  TP 수: {np.mean(all_tp_counts):.1f}±{np.std(all_tp_counts):.1f}")

    # 가설 A 종합
    if all_fp_curvatures:
        curv_arr = np.array(all_fp_curvatures)
        frac_positive = np.mean(curv_arr > 0)
        log(f"\n  [가설 A — 유사-영점] |F₂|'' > 0 비율: {frac_positive:.3f} ({len(curv_arr)}개)")
        if frac_positive > 0.7:
            log(f"    → 지지: 대부분의 FP가 |F₂| 풍경의 극소점 (유사-영점)")
        elif frac_positive < 0.3:
            log(f"    → 기각: FP는 극소점이 아님 (변곡점 또는 단조 구간)")
        else:
            log(f"    → 불분명: 혼합 패턴")

    # 가설 B 종합
    if all_fp_distances:
        dist_arr = np.array(all_fp_distances)
        frac_near = np.mean(dist_arr < resolution)
        frac_near_2x = np.mean(dist_arr < 2 * resolution)
        log(f"\n  [가설 B — 해상도 한계] FP→zero 거리:")
        log(f"    mean={dist_arr.mean():.4f}, median={np.median(dist_arr):.4f}")
        log(f"    < resolution: {frac_near:.3f}")
        log(f"    < 2×resolution: {frac_near_2x:.3f}")

        # GUE 비교
        if len(gue_spacings) > 0:
            mean_spacing = np.mean(np.diff(zeros_arr))
            normalized_fp_dist = dist_arr / mean_spacing
            log(f"    정규화 FP-zero 거리 (영점 간격 단위): mean={normalized_fp_dist.mean():.4f}")
            log(f"    GUE pair spacing mean: {gue_spacings.mean():.4f}")
            if normalized_fp_dist.mean() < 0.5:
                log(f"    → 지지: FP가 영점 근처에 밀집 (해상도 한계)")
            else:
                log(f"    → 기각: FP가 영점에서 먼 위치 (phantom minimum)")

    # 가설 C 종합
    if len(all_fp_indices) >= 2 and pair_overlaps:
        mean_jaccard = np.mean(pair_overlaps)
        log(f"\n  [가설 C — 무작위 잡음] 시드 간 Jaccard: {mean_jaccard:.4f}")
        if mean_jaccard > 0.5:
            log(f"    → 기각: FP 위치가 시드 간 일관적 (구조적 원인)")
        elif mean_jaccard < 0.1:
            log(f"    → 지지: FP 위치가 시드 간 비일관 (잡음)")
        else:
            log(f"    → 부분적: 일부 FP는 구조적, 일부는 잡음")

    # 최종 판정
    log(f"\n  --- 최종 판정 ---")

    verdicts = []
    if all_fp_curvatures:
        frac_pos = np.mean(np.array(all_fp_curvatures) > 0)
        if frac_pos > 0.5:
            verdicts.append("A(유사-영점)")
    if all_fp_distances:
        frac_nr = np.mean(np.array(all_fp_distances) < 2 * resolution)
        if frac_nr > 0.5:
            verdicts.append("B(해상도)")
    if len(all_fp_indices) >= 2 and pair_overlaps:
        if np.mean(pair_overlaps) < 0.2:
            verdicts.append("C(잡음)")

    if verdicts:
        log(f"  지배적 메커니즘: {' + '.join(verdicts)}")
    else:
        log(f"  지배적 메커니즘: 불명확 — 추가 분석 필요")

    elapsed = time.time() - start_time
    log(f"\n  총 실행 시간: {elapsed:.0f}s")

    # 저장
    os.makedirs(os.path.expanduser('~/Desktop/gdl_unified/results'), exist_ok=True)
    p = os.path.expanduser('~/Desktop/gdl_unified/results/fp_anatomy.txt')
    with open(p, 'w') as f:
        f.write('\n'.join(out))
    log(f"\n  저장: {p}")


if __name__ == '__main__':
    main()
