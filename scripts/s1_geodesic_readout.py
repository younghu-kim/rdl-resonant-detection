"""
=============================================================================
[Project RDL] S¹ 기하학적 각도 예측 실험
=============================================================================
±π 불연속의 기하학적 해결: arg ξ를 스칼라로 예측하는 대신
(cos θ, sin θ) ∈ S¹ 으로 예측하고 geodesic(대원) 손실 사용.

핵심 차이점 vs complex_vector_readout:
  - complex_vector: (Re ξ, Im ξ) → |ξ| 변동, 영점에서 0으로 수축
  - S¹: (cos θ, sin θ) → 단위원 위, ±π 분기점 자체가 없음

비교: (A) 기존 파이프라인 (TotalResonanceLoss)
      (B) S¹ geodesic readout (L_res+L_curv+L_pqo + geodesic on S¹)
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


def eval_F2(model, dataset, batch_size=64):
    """올바른 F₂ 잔차 기반 영점 탐지 평가 (배치 처리)
    F₂ = (e^{-iφ} · (L_G - ψ)).imag — blind_zero_prediction.py 와 동일
    """
    model.eval()
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)
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

    f2_arr = np.abs(np.concatenate(f2_vals))
    is_zero = dataset.is_near_zero.numpy()

    f2_zero = f2_arr[is_zero].mean() if is_zero.any() else 0
    f2_nonzero = f2_arr[~is_zero].mean() if (~is_zero).any() else 1
    ratio = f2_zero / (f2_nonzero + 1e-12)

    # 영점 검출: F₂ 잔차가 중앙값의 10% 미만
    threshold = np.median(f2_arr) * 0.1 if len(f2_arr) > 0 else 0.01
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


class GeodesicLoss(nn.Module):
    """S¹ 위의 geodesic 손실: L = 1 - cos(θ_pred - θ_true)
    (cos θ, sin θ) 표현에서 직접 계산: L = 1 - (c_p*c_t + s_p*s_t)
    """
    def forward(self, pred_cs, true_phase):
        # pred_cs: [batch, 2] = (cos θ_pred, sin θ_pred) — 정규화 전
        # true_phase: [batch, ...] = θ_true (라디안)

        # 예측을 단위원 위로 정규화
        norm = torch.sqrt(pred_cs[:, 0]**2 + pred_cs[:, 1]**2 + 1e-12)
        cos_pred = pred_cs[:, 0] / norm
        sin_pred = pred_cs[:, 1] / norm

        # 타겟의 cos/sin
        cos_true = torch.cos(true_phase[:, 0])
        sin_true = torch.sin(true_phase[:, 0])

        # Geodesic: 1 - cos(Δθ) = 1 - (cos_p*cos_t + sin_p*sin_t)
        cos_delta = cos_pred * cos_true + sin_pred * sin_true
        loss = 1.0 - cos_delta  # ∈ [0, 2]

        return loss.mean()


def train_s1_geodesic(dataset, seed):
    """S¹ geodesic readout: (cos θ, sin θ) 예측 + geodesic 손실"""
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

    # S¹ readout: HIDDEN complex → HIDDEN*2 real → 2 (cos θ, sin θ)
    readout = nn.Linear(HIDDEN * 2, 2, dtype=PrecisionManager.REAL_DTYPE)

    res_loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=0.0, lambda_pqo=0.5,
        pqo_mode='cos2'
    )
    geo_loss_fn = GeodesicLoss()

    params = list(model.parameters()) + list(readout.parameters())
    optimizer = torch.optim.Adam(params, lr=LR)
    lambda_geo = 2.0

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

            # S¹ readout: 프로젝션 전 hidden → (cos θ, sin θ)
            Z_hidden = _pre_proj['Z']
            z_flat = torch.cat([Z_hidden.real, Z_hidden.imag], dim=-1).to(dtype=PrecisionManager.REAL_DTYPE)
            cs_pred = readout(z_flat)  # [batch, 2]

            # 타겟: y_batch의 arg(ξ) → Ψ_target에서 위상각 추출
            # y_batch = (Re ξ, Im ξ) → θ = atan2(Im, Re)
            y_target = y_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            theta_true = torch.atan2(y_target[:, 1:2], y_target[:, 0:1])
            geo_loss = geo_loss_fn(cs_pred, theta_true)

            # 보조 공명 손실
            res_loss, _ = res_loss_fn(**outputs)
            if torch.isnan(res_loss) or torch.isinf(res_loss):
                res_loss = torch.tensor(0.0, dtype=PrecisionManager.REAL_DTYPE)

            total = res_loss + lambda_geo * geo_loss
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
                    cs_pred = readout(z_flat)
                    y_target = y_batch.to(dtype=PrecisionManager.REAL_DTYPE)
                    theta_true = torch.atan2(y_target[:, 1:2], y_target[:, 0:1])
                    geo_loss = geo_loss_fn(cs_pred, theta_true)
                    res_loss, _ = res_loss_fn(**o)
                    if torch.isnan(res_loss) or torch.isinf(res_loss):
                        res_loss = torch.tensor(0.0)
                    loss = res_loss + lambda_geo * geo_loss
                    if not (torch.isnan(loss) or torch.isinf(loss)):
                        va += loss.item(); nv += 1
            vl_m = va / max(1, nv)
            best_val = min(best_val, vl_m)
            elapsed = time.time() - t_start
            print(f"      ep {ep+1}/{EPOCHS}: train={ep_loss/max(1,ep_n):.5f}, val={vl_m:.5f}, best={best_val:.5f}, {elapsed:.0f}s", flush=True)

    hook_handle.remove()
    return model, readout, best_val


def eval_s1_phase_error(model, readout, dataset, _pre_proj_dict, batch_size=64):
    """S¹ readout의 위상 예측 정확도 평가"""
    model.eval(); readout.eval()
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)
    all_errors = []
    with torch.enable_grad():
        for X_batch, y_batch in loader:
            X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            model(X_in)
            Z_hidden = _pre_proj_dict['Z']
            z_flat = torch.cat([Z_hidden.real, Z_hidden.imag], dim=-1).to(dtype=PrecisionManager.REAL_DTYPE)
            cs_pred = readout(z_flat).detach()
            norm = torch.sqrt(cs_pred[:, 0]**2 + cs_pred[:, 1]**2 + 1e-12)
            cos_p = cs_pred[:, 0] / norm
            sin_p = cs_pred[:, 1] / norm
            theta_pred = torch.atan2(sin_p, cos_p)

            y_target = y_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            theta_true = torch.atan2(y_target[:, 1], y_target[:, 0])

            # 최단 각도 거리
            diff = theta_pred - theta_true
            angular_error = torch.abs(torch.atan2(torch.sin(diff), torch.cos(diff)))
            all_errors.append(angular_error.numpy())

    errors = np.concatenate(all_errors)
    is_zero = dataset.is_near_zero.numpy()

    err_zero = errors[is_zero].mean() if is_zero.any() else 0
    err_nonzero = errors[~is_zero].mean() if (~is_zero).any() else 0
    err_all = errors.mean()

    return err_all, err_zero, err_nonzero


def main():
    t_min, t_max = 100.0, 200.0
    num_points = 1000
    in_features = 128
    seeds = [42, 7, 123]

    out = []
    def log(msg):
        print(msg); out.append(msg)

    log("=" * 72)
    log("  S¹ 기하학적 각도 예측 실험: 기존 vs S¹ geodesic readout")
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
    s1_results = []

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

        log("  [S¹ Geodesic] 훈련...")
        t0 = time.time()
        # hook을 재설정해야 하므로 train_s1_geodesic 내부에서 관리
        model_s, readout_s, val_s = train_s1_geodesic(ds, s)
        f2z, f2nz, ratio_s, det_s, tot_z = eval_F2(model_s, ds)

        # S¹ 위상 오차 평가를 위한 hook 재설정
        _pre_proj_s = {}
        def _capture(module, inp, out):
            _pre_proj_s['Z'] = inp[0]
        h = model_s.out_projection.register_forward_hook(_capture)
        err_all, err_zero, err_nonzero = eval_s1_phase_error(model_s, readout_s, ds, _pre_proj_s)
        h.remove()

        dt = time.time() - t0
        log(f"    val={val_s:.5f}, |F₂| ratio={ratio_s:.4f}, 검출={det_s}/{tot_z}, time={dt:.0f}s")
        log(f"    S¹ 위상오차: 전체={err_all:.4f}rad, 영점={err_zero:.4f}rad, 비영점={err_nonzero:.4f}rad")
        s1_results.append((val_s, ratio_s, det_s, tot_z, err_all, err_zero, err_nonzero))

    elapsed = time.time() - start

    log(f"\n{'='*72}")
    log("  요약")
    log(f"{'='*72}")

    br = [r[1] for r in baseline_results]
    bd = [r[2] for r in baseline_results]
    sr = [r[1] for r in s1_results]
    sd = [r[2] for r in s1_results]
    tot = baseline_results[0][3]

    log(f"  {'방식':<20} {'|F₂| ratio':>15} {'검출':>10}")
    log(f"  {'-'*20} {'-'*15} {'-'*10}")
    log(f"  {'Baseline':<20} {np.mean(br):>7.4f}±{np.std(br):.4f} {np.mean(bd):>6.1f}/{tot}")
    log(f"  {'S¹ Geodesic':<20} {np.mean(sr):>7.4f}±{np.std(sr):.4f} {np.mean(sd):>6.1f}/{tot}")

    improvement = (np.mean(br) - np.mean(sr)) / (np.mean(br) + 1e-12) * 100
    det_improve = np.mean(sd) - np.mean(bd)

    log(f"\n  ratio 개선: {improvement:+.1f}%")
    log(f"  검출 개선: {det_improve:+.1f}개")

    # S¹ 위상 오차 요약
    err_z = [r[5] for r in s1_results]
    err_nz = [r[6] for r in s1_results]
    log(f"\n  S¹ 위상 오차 (rad):")
    log(f"    영점 근방:  {np.mean(err_z):.4f}±{np.std(err_z):.4f}")
    log(f"    비영점:    {np.mean(err_nz):.4f}±{np.std(err_nz):.4f}")
    log(f"    비율:     {np.mean(err_z)/(np.mean(err_nz)+1e-12):.2f}x")

    log(f"\n  총 실행 시간: {elapsed:.0f}s")

    if np.mean(sr) < np.mean(br) * 0.8 or np.mean(sd) > np.mean(bd) * 1.2:
        verdict = "양성: S¹ geodesic이 유의미하게 개선"
    elif np.mean(sr) > np.mean(br) * 1.2:
        verdict = "음성: 기존 방식이 우세"
    else:
        verdict = "중립: 유의미한 차이 없음"
    log(f"\n  판정: {verdict}")

    os.makedirs(os.path.expanduser('~/Desktop/gdl_unified/results'), exist_ok=True)
    p = os.path.expanduser('~/Desktop/gdl_unified/results/s1_geodesic_readout.txt')
    with open(p, 'w') as f:
        f.write('\n'.join(out))
    log(f"  저장: {p}")


if __name__ == '__main__':
    main()
