"""
=============================================================================
[Project RDL] 정밀도 개선 실험 — 후처리 필터링
=============================================================================
블라인드 영점 예측 recall 99% / precision 15% → 후처리로 정밀도 개선.

전략: 앙상블 합의, 깊이 필터, 간격 필터
"""

import sys, os, time, math
import numpy as np
import torch
from scipy.signal import argrelmin

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
from gdl.rdl.constants import PrecisionManager
from gdl.rdl.pipeline.xi_feature_dataset import (
    get_or_build_cache, compute_zeros_in_range, XiFeatureDataset
)
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss

PrecisionManager.setup_precision()

HIDDEN = 64
EPOCHS = 200
LR = 1e-3
BATCH = 32


def train_model(dataset, seed):
    torch.manual_seed(seed); np.random.seed(seed)

    val_size = int(len(dataset) * 0.2)
    train_size = len(dataset) - val_size
    train_ds, _ = torch.utils.data.random_split(
        dataset, [train_size, val_size],
        generator=torch.Generator().manual_seed(seed)
    )
    train_loader = torch.utils.data.DataLoader(train_ds, batch_size=BATCH, shuffle=True, drop_last=True)

    model = MasterResonantNetwork(
        in_features=dataset.in_features, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5,
        pqo_mode='cos2'
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)

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

        if (ep + 1) % 50 == 0 or ep == EPOCHS - 1:
            elapsed = time.time() - t_start
            print(f"      ep {ep+1}/{EPOCHS}: loss={ep_loss/max(1,ep_n):.5f}, {elapsed:.0f}s", flush=True)

    return model


def compute_dense_f2(model, dataset, t_min, t_max, n_grid=2000, batch_size=64):
    model.eval()
    t_grid = torch.linspace(t_min, t_max, n_grid, dtype=PrecisionManager.REAL_DTYPE)
    f2_grid = []

    # 배치 단위로 특징 벡터 구성 및 평가
    for start in range(0, n_grid, batch_size):
        end = min(start + batch_size, n_grid)
        X_batch = torch.stack([dataset._build_features(t_grid[i]) for i in range(start, end)])
        X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
        X_in.requires_grad_(True)
        with torch.enable_grad():
            out = model(X_in)
        f2_batch = out['Z_out'].abs().mean(dim=-1).detach().numpy()
        f2_grid.append(f2_batch)

    return t_grid.numpy(), np.concatenate(f2_grid)


def find_minima_raw(t_grid, f2_grid, order=5):
    idx = argrelmin(f2_grid, order=order)[0]
    return t_grid[idx], f2_grid[idx]


def filter_ensemble(all_minima_t, merge_radius=0.3, min_votes=2):
    if not all_minima_t:
        return np.array([])
    pool = np.sort(np.concatenate(all_minima_t))
    if len(pool) == 0:
        return np.array([])

    clusters = []
    current = [pool[0]]
    for t in pool[1:]:
        if t - current[-1] < merge_radius:
            current.append(t)
        else:
            clusters.append(current)
            current = [t]
    clusters.append(current)

    consensus = []
    for c in clusters:
        vote_count = 0
        for seed_minima in all_minima_t:
            for t in c:
                if np.min(np.abs(seed_minima - t)) < merge_radius:
                    vote_count += 1
                    break
        if vote_count >= min_votes:
            consensus.append(np.mean(c))

    return np.array(consensus)


def filter_spacing(t_arr, min_gap=1.5):
    if len(t_arr) <= 1:
        return t_arr
    sorted_t = np.sort(t_arr)
    keep = [sorted_t[0]]
    for t in sorted_t[1:]:
        if t - keep[-1] >= min_gap:
            keep.append(t)
    return np.array(keep)


def evaluate(pred_t, true_zeros, hit_threshold=0.3):
    if len(pred_t) == 0:
        return 0, 0, 0
    true_arr = np.array(true_zeros)
    hits = 0; matched = set()
    for t in pred_t:
        dists = np.abs(true_arr - t)
        idx = np.argmin(dists)
        if dists[idx] < hit_threshold and idx not in matched:
            hits += 1; matched.add(idx)
    precision = hits / len(pred_t) if len(pred_t) > 0 else 0
    recall = hits / len(true_arr) if len(true_arr) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall + 1e-12)
    return precision, recall, f1


def main():
    t_train_min, t_train_max = 100.0, 150.0
    t_pred_min, t_pred_max = 150.0, 200.0
    num_points = 1000
    in_features = 128
    seeds = [42, 7, 123, 2026, 314]
    n_grid = 2000

    out = []
    def log(msg):
        print(msg); out.append(msg)

    log("=" * 72)
    log("  정밀도 개선 실험 — 후처리 필터링")
    log("=" * 72)
    log(f"  학습: t∈[{t_train_min},{t_train_max}], 예측: t∈[{t_pred_min},{t_pred_max}]")
    log(f"  seeds={seeds}, ep={EPOCHS}, if={in_features}")
    log("")

    all_zeros = compute_zeros_in_range(100.0, 200.0, dps=25)
    pred_zeros = [z for z in all_zeros if t_pred_min <= z <= t_pred_max]

    cache_path = os.path.expanduser(
        f'~/Desktop/gdl_unified/outputs/xi_cache_t100.0-200.0_n{num_points}.pt'
    )
    cache_data = get_or_build_cache(cache_path, 100.0, 200.0, num_points, zeros_list=all_zeros)

    train_mask = (cache_data['t'] >= t_train_min) & (cache_data['t'] <= t_train_max)
    train_cache = {k: v[train_mask] for k, v in cache_data.items()}

    ds_full = XiFeatureDataset(cache_data, in_features=in_features)

    start = time.time()

    log("  [1단계] 모델 훈련 + 밀집 격자 F₂")
    all_minima_t = []
    all_f2_grids = []

    for s in seeds:
        log(f"    seed={s} 훈련 중...")
        ds_train = XiFeatureDataset(train_cache, in_features=in_features)
        model = train_model(ds_train, s)

        t_grid, f2_grid = compute_dense_f2(model, ds_full, t_pred_min, t_pred_max, n_grid)
        all_f2_grids.append(f2_grid)

        t_min_arr, f2_min_arr = find_minima_raw(t_grid, f2_grid)
        all_minima_t.append(t_min_arr)
        log(f"    → {len(t_min_arr)}개 극소")

    log(f"\n  [2단계] 필터 적용 (예측 영점: {len(pred_zeros)}개)")
    log(f"  {'필터':<30} {'정밀도':>8} {'재현율':>8} {'F1':>8} {'예측수':>6}")
    log(f"  {'-'*30} {'-'*8} {'-'*8} {'-'*8} {'-'*6}")

    results = {}

    # (a) 원시
    raw_p, raw_r, raw_f1 = [], [], []
    for mt in all_minima_t:
        p, r, f = evaluate(mt, pred_zeros)
        raw_p.append(p); raw_r.append(r); raw_f1.append(f)
    results['raw'] = (np.mean(raw_p), np.mean(raw_r), np.mean(raw_f1), np.mean([len(m) for m in all_minima_t]))
    log(f"  {'(a) 원시 (시드평균)':<30} {np.mean(raw_p):>8.3f} {np.mean(raw_r):>8.3f} {np.mean(raw_f1):>8.3f} {results['raw'][3]:>6.0f}")

    # (b) 앙상블 ≥2
    c2 = filter_ensemble(all_minima_t, merge_radius=0.3, min_votes=2)
    p, r, f = evaluate(c2, pred_zeros)
    results['ensemble_2'] = (p, r, f, len(c2))
    log(f"  {'(b) 앙상블 ≥2 시드':<30} {p:>8.3f} {r:>8.3f} {f:>8.3f} {len(c2):>6}")

    # (c) 앙상블 ≥3
    c3 = filter_ensemble(all_minima_t, merge_radius=0.3, min_votes=3)
    p, r, f = evaluate(c3, pred_zeros)
    results['ensemble_3'] = (p, r, f, len(c3))
    log(f"  {'(c) 앙상블 ≥3 시드':<30} {p:>8.3f} {r:>8.3f} {f:>8.3f} {len(c3):>6}")

    # (d) 앙상블 ≥3 + 깊이
    mean_f2 = np.mean(all_f2_grids, axis=0)
    if len(c3) > 0:
        c3_f2 = np.interp(c3, t_grid, mean_f2)
        depth_mask = c3_f2 < np.median(mean_f2) * 0.5
        c3_deep = c3[depth_mask]
    else:
        c3_deep = np.array([])
    p, r, f = evaluate(c3_deep, pred_zeros)
    results['ensemble_3_depth'] = (p, r, f, len(c3_deep))
    log(f"  {'(d) 앙상블 ≥3 + 깊이':<30} {p:>8.3f} {r:>8.3f} {f:>8.3f} {len(c3_deep):>6}")

    # (e) 앙상블 ≥3 + 간격
    c3_spaced = filter_spacing(c3, min_gap=1.5) if len(c3) > 0 else np.array([])
    p, r, f = evaluate(c3_spaced, pred_zeros)
    results['ensemble_3_spacing'] = (p, r, f, len(c3_spaced))
    log(f"  {'(e) 앙상블 ≥3 + 간격':<30} {p:>8.3f} {r:>8.3f} {f:>8.3f} {len(c3_spaced):>6}")

    # (f) 종합
    final = filter_spacing(c3_deep, min_gap=1.5) if len(c3_deep) > 0 else np.array([])
    p, r, f = evaluate(final, pred_zeros)
    results['final'] = (p, r, f, len(final))
    log(f"  {'(f) 종합 (≥3+깊이+간격)':<30} {p:>8.3f} {r:>8.3f} {f:>8.3f} {len(final):>6}")

    elapsed = time.time() - start
    log(f"\n  총 실행 시간: {elapsed:.0f}s")

    best_key = max(results, key=lambda k: results[k][2])
    best = results[best_key]
    log(f"\n  최고 F1: {best_key} → P={best[0]:.3f}, R={best[1]:.3f}, F1={best[2]:.3f}")

    if best[0] > results['raw'][0] * 1.5:
        log(f"  판정: 양성 — 정밀도 {best[0]/results['raw'][0]:.1f}배 개선")
    else:
        log(f"  판정: 중립 — 정밀도 유의미한 개선 없음")

    p_path = os.path.expanduser('~/Desktop/gdl_unified/results/precision_filter.txt')
    with open(p_path, 'w') as f:
        f.write('\n'.join(out))
    log(f"  저장: {p_path}")


if __name__ == '__main__':
    main()
