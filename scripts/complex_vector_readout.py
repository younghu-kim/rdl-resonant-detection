"""
=============================================================================
[Project RDL] 복소 벡터 읽기 실험
=============================================================================
±π 불연속 병목의 아키텍처적 해결: arg ξ를 직접 예측하는 대신
ξ(1/2+it)를 복소 벡터 (Re, Im)로 예측하고 사후에 arg를 계산.

비교: (A) 기존 파이프라인 (TotalResonanceLoss)
      (B) 복소 벡터 읽기 (L_res+L_curv+L_pqo + MSE on Re/Im)
"""

import sys, os, time, math
import numpy as np
import torch
import torch.nn as nn

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


def eval_F2(model, dataset, threshold=0.01, batch_size=64):
    """F₂ 기반 영점 탐지 평가 (배치 처리)"""
    model.eval()
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)
    f2_vals = []
    with torch.enable_grad():
        for X_batch, _ in loader:
            X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            out = model(X_in)
            Z = out['Z_out']
            f2_batch = Z.abs().mean(dim=-1).detach().numpy()
            f2_vals.append(f2_batch)

    f2_arr = np.concatenate(f2_vals)
    is_zero = dataset.is_near_zero.numpy()

    f2_zero = f2_arr[is_zero].mean() if is_zero.any() else 0
    f2_nonzero = f2_arr[~is_zero].mean() if (~is_zero).any() else 1
    ratio = f2_zero / (f2_nonzero + 1e-12)

    detected = int(np.sum(f2_arr[is_zero] < threshold)) if is_zero.any() else 0
    total_zeros = int(is_zero.sum())

    return f2_zero, f2_nonzero, ratio, detected, total_zeros


def make_loaders(dataset, seed):
    val_size = int(len(dataset) * 0.2)
    train_size = len(dataset) - val_size
    train_ds, val_ds = torch.utils.data.random_split(
        dataset, [train_size, val_size],
        generator=torch.Generator().manual_seed(seed)
    )
    train_loader = torch.utils.data.DataLoader(train_ds, batch_size=BATCH, shuffle=True, drop_last=True)
    val_loader = torch.utils.data.DataLoader(val_ds, batch_size=BATCH, shuffle=False)
    return train_loader, val_loader


def train_baseline(dataset, seed):
    """기존 방식: TotalResonanceLoss 전체 파이프라인"""
    torch.manual_seed(seed); np.random.seed(seed)
    train_loader, val_loader = make_loaders(dataset, seed)

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
            print(f"      ep {ep+1}/{EPOCHS}: train={ep_loss/max(1,ep_n):.5f}, val={vl_m:.5f}, best={best_val:.5f}, {elapsed:.0f}s", flush=True)

    return model, best_val


def train_complex_vector(dataset, seed):
    """복소 벡터 읽기: Re/Im MSE + L_res/L_curv/L_pqo 보조 (L_tgt=0)
    프로젝션 전 hidden state (HIDDEN complex = HIDDEN*2 real)에서 (Re ξ, Im ξ) 예측.
    """
    torch.manual_seed(seed); np.random.seed(seed)
    train_loader, val_loader = make_loaders(dataset, seed)

    model = MasterResonantNetwork(
        in_features=dataset.in_features, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )

    # forward hook: out_projection 입력(프로젝션 전 Z)을 캡처
    _pre_proj = {}
    def _capture_pre_proj(module, inp, out):
        _pre_proj['Z'] = inp[0]
    hook_handle = model.out_projection.register_forward_hook(_capture_pre_proj)

    # 복소 벡터 읽기 헤드: HIDDEN complex → HIDDEN*2 real → 2 (Re ξ, Im ξ)
    readout = nn.Linear(HIDDEN * 2, 2, dtype=PrecisionManager.REAL_DTYPE)

    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=0.0, lambda_pqo=0.5,
        pqo_mode='cos2'
    )
    mse_fn = nn.MSELoss()

    params = list(model.parameters()) + list(readout.parameters())
    optimizer = torch.optim.Adam(params, lr=LR)
    lambda_mse = 2.0

    best_val = float('inf')
    t_start = time.time()
    for ep in range(EPOCHS):
        model.train(); readout.train()
        ep_loss = 0.0; ep_n = 0
        for X_batch, y_batch in train_loader:
            X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            optimizer.zero_grad(set_to_none=True)
            outputs = model(X_in)

            # 프로젝션 전 hidden state로 복소 벡터 읽기
            Z_hidden = _pre_proj['Z']  # [batch, HIDDEN] complex
            z_flat = torch.cat([Z_hidden.real, Z_hidden.imag], dim=-1).to(dtype=PrecisionManager.REAL_DTYPE)
            xi_pred = readout(z_flat)
            y_target = y_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            mse_loss = mse_fn(xi_pred, y_target)

            # 보조 공명 손실
            res_loss, _ = loss_fn(**outputs)
            if torch.isnan(res_loss) or torch.isinf(res_loss):
                res_loss = torch.tensor(0.0, dtype=PrecisionManager.REAL_DTYPE)

            total = res_loss + lambda_mse * mse_loss
            if torch.isnan(total) or torch.isinf(total):
                continue
            total.backward()
            torch.nn.utils.clip_grad_norm_(params, 5.0)
            optimizer.step()
            ep_loss += total.item(); ep_n += 1

        if (ep + 1) % 30 == 0 or ep == EPOCHS - 1:
            model.eval(); readout.eval()
            va = 0.0; nv = 0
            with torch.enable_grad():
                for X_batch, y_batch in val_loader:
                    X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
                    X_in.requires_grad_(True)
                    o = model(X_in)
                    Z_hidden = _pre_proj['Z']
                    z_flat = torch.cat([Z_hidden.real, Z_hidden.imag], dim=-1).to(dtype=PrecisionManager.REAL_DTYPE)
                    xi_pred = readout(z_flat)
                    y_target = y_batch.to(dtype=PrecisionManager.REAL_DTYPE)
                    mse_loss = mse_fn(xi_pred, y_target)
                    res_loss, _ = loss_fn(**o)
                    if torch.isnan(res_loss) or torch.isinf(res_loss):
                        res_loss = torch.tensor(0.0)
                    loss = res_loss + lambda_mse * mse_loss
                    if not (torch.isnan(loss) or torch.isinf(loss)):
                        va += loss.item(); nv += 1
            vl_m = va / max(1, nv)
            best_val = min(best_val, vl_m)
            elapsed = time.time() - t_start
            print(f"      ep {ep+1}/{EPOCHS}: train={ep_loss/max(1,ep_n):.5f}, val={vl_m:.5f}, best={best_val:.5f}, {elapsed:.0f}s", flush=True)

    hook_handle.remove()
    return model, best_val


def main():
    t_min, t_max = 100.0, 200.0
    num_points = 1000
    in_features = 128
    seeds = [42, 7, 123]

    out = []
    def log(msg):
        print(msg); out.append(msg)

    log("=" * 72)
    log("  복소 벡터 읽기 실험: 기존 vs 복소 벡터 readout")
    log("=" * 72)
    log(f"  t∈[{t_min},{t_max}], {num_points}점, if={in_features}, ep={EPOCHS}")
    log(f"  seeds={seeds}")
    log("")

    zeros_list = compute_zeros_in_range(t_min, t_max)
    cache_path = os.path.expanduser(
        f'~/Desktop/gdl_unified/outputs/xi_cache_t{t_min}-{t_max}_n{num_points}.pt'
    )
    cache_data = get_or_build_cache(cache_path, t_min, t_max, num_points, zeros_list=zeros_list)

    start = time.time()
    baseline_results = []
    complex_results = []

    for s in seeds:
        log(f"\n{'─'*72}")
        log(f"  seed={s}")
        log(f"{'─'*72}")

        ds = XiFeatureDataset(cache_data, in_features=in_features)

        log("  [Baseline] 훈련...")
        t0 = time.time()
        model_b, val_b = train_baseline(ds, s)
        f2z, f2nz, ratio_b, det_b, tot_z = eval_F2(model_b, ds)
        dt = time.time() - t0
        log(f"    val={val_b:.5f}, |F₂| ratio={ratio_b:.4f}, 검출={det_b}/{tot_z}, time={dt:.0f}s")
        baseline_results.append((val_b, ratio_b, det_b, tot_z))

        log("  [Complex Vector] 훈련...")
        t0 = time.time()
        model_c, val_c = train_complex_vector(ds, s)
        f2z, f2nz, ratio_c, det_c, tot_z = eval_F2(model_c, ds)
        dt = time.time() - t0
        log(f"    val={val_c:.5f}, |F₂| ratio={ratio_c:.4f}, 검출={det_c}/{tot_z}, time={dt:.0f}s")
        complex_results.append((val_c, ratio_c, det_c, tot_z))

    elapsed = time.time() - start

    log(f"\n{'='*72}")
    log("  요약")
    log(f"{'='*72}")

    br = [r[1] for r in baseline_results]
    bd = [r[2] for r in baseline_results]
    cr = [r[1] for r in complex_results]
    cd = [r[2] for r in complex_results]
    tot = baseline_results[0][3]

    log(f"  {'방식':<20} {'|F₂| ratio':>15} {'검출':>10}")
    log(f"  {'-'*20} {'-'*15} {'-'*10}")
    log(f"  {'Baseline':<20} {np.mean(br):>7.4f}±{np.std(br):.4f} {np.mean(bd):>6.1f}/{tot}")
    log(f"  {'Complex Vector':<20} {np.mean(cr):>7.4f}±{np.std(cr):.4f} {np.mean(cd):>6.1f}/{tot}")

    improvement = (np.mean(br) - np.mean(cr)) / (np.mean(br) + 1e-12) * 100
    det_improve = np.mean(cd) - np.mean(bd)

    log(f"\n  ratio 개선: {improvement:+.1f}%")
    log(f"  검출 개선: {det_improve:+.1f}개")
    log(f"  총 실행 시간: {elapsed:.0f}s")

    if np.mean(cr) < np.mean(br) * 0.8 or np.mean(cd) > np.mean(bd) * 1.2:
        verdict = "양성: 복소 벡터 읽기가 유의미하게 개선"
    elif np.mean(cr) > np.mean(br) * 1.2:
        verdict = "음성: 기존 방식이 우세"
    else:
        verdict = "중립: 유의미한 차이 없음"
    log(f"\n  판정: {verdict}")

    os.makedirs(os.path.expanduser('~/Desktop/gdl_unified/results'), exist_ok=True)
    p = os.path.expanduser('~/Desktop/gdl_unified/results/complex_vector_readout.txt')
    with open(p, 'w') as f:
        f.write('\n'.join(out))
    log(f"  저장: {p}")


if __name__ == '__main__':
    main()
