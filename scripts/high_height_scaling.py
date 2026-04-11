"""
=============================================================================
[Project RDL] 고높이 스케일링 실험
=============================================================================
t∈[100,200] 에서 검증된 F₂를 t∈[500,600], t∈[1000,1100] 에서도 테스트.
높이가 올라가면 영점 밀도 증가 → K가 부족해지거나 F₂ 성능 저하되는가?
"""

import sys, os, time, math
import numpy as np
import torch

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
from gdl.rdl.constants import PrecisionManager
from gdl.rdl.pipeline.xi_feature_dataset import (
    get_or_build_cache, compute_zeros_in_range, XiFeatureDataset
)
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss

PrecisionManager.setup_precision()

HIDDEN = 64
EPOCHS = 150
LR = 1e-3
BATCH = 32


def train_and_eval(dataset, seed):
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

        if (ep + 1) % 30 == 0 or ep == EPOCHS - 1:
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
            print(f"      ep {ep+1}/{EPOCHS}: train={ep_loss/max(1,ep_n):.5f}, val={vl_m:.5f}, {elapsed:.0f}s", flush=True)

    # F₂ 평가 (배치 처리)
    model.eval()
    eval_loader = torch.utils.data.DataLoader(dataset, batch_size=64, shuffle=False)
    f2_vals = []
    with torch.enable_grad():
        for X_batch, _ in eval_loader:
            X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            out = model(X_in)
            f2_batch = out['Z_out'].abs().mean(dim=-1).detach().numpy()
            f2_vals.append(f2_batch)

    f2_arr = np.concatenate(f2_vals)
    is_zero = dataset.is_near_zero.numpy()

    f2_zero = f2_arr[is_zero].mean() if is_zero.any() else 0
    f2_nonzero = f2_arr[~is_zero].mean() if (~is_zero).any() else 1
    ratio = f2_zero / (f2_nonzero + 1e-12)

    detected = int(np.sum(f2_arr[is_zero] < 0.01)) if is_zero.any() else 0
    total_zeros = int(is_zero.sum())

    return best_val, f2_zero, f2_nonzero, ratio, detected, total_zeros


def main():
    ranges = [
        (100.0, 200.0),
        (500.0, 600.0),
        (1000.0, 1100.0),
    ]
    num_points = 1000
    in_features_list = [128, 256]
    seeds = [42, 7, 123]

    out = []
    def log(msg):
        print(msg); out.append(msg)

    log("=" * 72)
    log("  고높이 스케일링 실험")
    log("=" * 72)
    log(f"  구간: {ranges}")
    log(f"  in_features: {in_features_list}, seeds={seeds}, ep={EPOCHS}")
    log("")

    start = time.time()
    all_results = {}

    for t_min, t_max in ranges:
        log(f"\n{'='*72}")
        log(f"  구간: t∈[{t_min}, {t_max}]")
        log(f"{'='*72}")

        zeros_list = compute_zeros_in_range(t_min, t_max, dps=25)
        n_zeros = len(zeros_list)
        density = n_zeros / (t_max - t_min)

        log(f"  영점 수: {n_zeros}, 밀도: {density:.3f}/unit")

        cache_path = os.path.expanduser(
            f'~/Desktop/gdl_unified/outputs/xi_cache_t{t_min}-{t_max}_n{num_points}.pt'
        )
        cache_data = get_or_build_cache(cache_path, t_min, t_max, num_points,
                                         mp_dps=50, zeros_list=zeros_list)

        for K in in_features_list:
            log(f"\n  --- K={K} ---")
            ds = XiFeatureDataset(cache_data, in_features=K)

            ratios = []; dets = []
            for s in seeds:
                t0 = time.time()
                val, f2z, f2nz, ratio, det, tot = train_and_eval(ds, s)
                dt = time.time() - t0
                log(f"    seed={s}: val={val:.5f}, ratio={ratio:.4f}, det={det}/{tot}, time={dt:.0f}s")
                ratios.append(ratio)
                dets.append(det)

            key = f"t[{int(t_min)},{int(t_max)}]_K{K}"
            all_results[key] = {
                'n_zeros': n_zeros, 'density': density,
                'ratio_mean': np.mean(ratios), 'ratio_std': np.std(ratios),
                'det_mean': np.mean(dets), 'det_total': tot
            }
            log(f"  → K={K}: ratio={np.mean(ratios):.4f}±{np.std(ratios):.4f}, det={np.mean(dets):.1f}/{tot}")

    elapsed = time.time() - start

    log(f"\n{'='*72}")
    log("  스케일링 요약")
    log(f"{'='*72}")
    log(f"  {'구간_K':<25} {'영점수':>6} {'밀도':>8} {'ratio':>15} {'검출':>10}")
    log(f"  {'-'*25} {'-'*6} {'-'*8} {'-'*15} {'-'*10}")

    for key, r in all_results.items():
        log(f"  {key:<25} {r['n_zeros']:>6} {r['density']:>8.3f} "
            f"{r['ratio_mean']:>7.4f}±{r['ratio_std']:.4f} "
            f"{r['det_mean']:>5.1f}/{r['det_total']}")

    log(f"\n  총 실행 시간: {elapsed:.0f}s")

    base_key = [k for k in all_results if 't[100,200]' in k and 'K128' in k][0]
    base_ratio = all_results[base_key]['ratio_mean']
    degraded = any(
        r['ratio_mean'] > base_ratio * 2.0
        for key, r in all_results.items() if 't[100,200]' not in key
    )

    if degraded:
        log("\n  판정: 음성 — 고높이에서 성능 저하 관찰")
    else:
        log("\n  판정: 양성 — 고높이에서도 성능 유지")

    p = os.path.expanduser('~/Desktop/gdl_unified/results/high_height_scaling.txt')
    with open(p, 'w') as f:
        f.write('\n'.join(out))
    log(f"  저장: {p}")


if __name__ == '__main__':
    main()
