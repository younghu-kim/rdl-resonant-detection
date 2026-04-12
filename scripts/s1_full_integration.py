"""
=============================================================================
[Project RDL] S¹ Geodesic 전체 통합 실험
=============================================================================
L_tgt (target phase matching loss)를 S¹ geodesic 손실로 대체.

핵심 가설:
  기존 L_tgt = ||atan2(sin(Δφ), cos(Δφ))||² 는 ±π에서 atan2 기울기 불안정.
  S¹ geodesic L_geo = 1 - cos(arg(Z) - Ψ) 는 모든 곳에서 매끄러움.
  → 네트워크 내부 Z_out 구조가 ±π를 자연스럽게 우회하여
    F₂ ratio 개선 + 영점 검출 증가 기대.

비교:
  (A) Baseline: 표준 TotalResonanceLoss (L_res + L_curv + L_tgt + L_pqo)
  (B) S¹ Integration: L_tgt → L_geo 대체 (L_res + L_curv + L_geo + L_pqo)

이전 S¹ readout 결과와의 차이:
  - readout: 별도 head, Z_out 무영향 → F₂ ratio 무변화
  - integration: L_tgt 자체를 대체, Z_out 학습에 직접 영향 → F₂ ratio 변화 가능
"""

import sys, os, time, math
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

HIDDEN = 64
EPOCHS = 150
LR = 1e-3
BATCH = 32


def eval_F2(model, dataset, batch_size=64):
    """올바른 F₂ 잔차 기반 영점 탐지 평가 (blind_zero_prediction.py 방식)"""
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
    """기존 방식: TotalResonanceLoss 전체 파이프라인 (L_tgt 포함)"""
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


class GeodesicTargetLoss(nn.Module):
    """S¹ geodesic 손실을 Z_out의 위상에 직접 적용.

    L_geo = 1 - cos(arg(Z_out) - Psi_target)

    vs 기존 L_tgt = ||atan2(sin(Δφ), cos(Δφ))||²:
      - 둘 다 ±π 래핑 처리하지만, geodesic은 gradient가 sin(Δφ)로 모든 곳에서 연속
      - atan2 기반은 Δφ → ±π 근방에서 기울기가 0으로 collapse 가능
    """
    def __init__(self, reduction='mean'):
        super().__init__()
        self.reduction = reduction

    def forward(self, Z_out, Psi_target):
        if Z_out.dtype != PrecisionManager.COMPLEX_DTYPE:
            Z_out = Z_out.to(dtype=PrecisionManager.COMPLEX_DTYPE)
        if Psi_target.dtype != PrecisionManager.REAL_DTYPE:
            Psi_target = Psi_target.to(dtype=PrecisionManager.REAL_DTYPE)

        current_phase = torch.angle(Z_out)

        while Psi_target.dim() < current_phase.dim():
            Psi_target = Psi_target.unsqueeze(-1)

        # geodesic: 1 - cos(Δφ) ∈ [0, 2]
        delta = current_phase - Psi_target
        loss = 1.0 - torch.cos(delta)

        if self.reduction == 'mean':
            return loss.mean()
        elif self.reduction == 'sum':
            return loss.sum()
        return loss


def train_s1_integrated(dataset, seed):
    """S¹ 통합: L_tgt를 GeodesicTargetLoss로 대체

    전체 손실: L = λ_res·L_res + λ_curv·L_curv + λ_geo·L_geo + λ_pqo·L_pqo
    L_geo 가중치는 L_tgt와 동일하게 1.0 (공정 비교)
    """
    torch.manual_seed(seed); np.random.seed(seed)
    train_loader, val_loader = make_loaders(dataset, seed)

    model = MasterResonantNetwork(
        in_features=dataset.in_features, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    # L_tgt=0으로 비활성화하고, geodesic을 별도로 추가
    res_loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=0.0, lambda_pqo=0.5,
        pqo_mode='cos2'
    )
    geo_loss_fn = GeodesicTargetLoss(reduction='mean')
    lambda_geo = 1.0  # L_tgt와 동일 가중치

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

            # 공명 손실 (L_res + L_curv + L_pqo, L_tgt=0)
            res_loss, _ = res_loss_fn(**outputs)
            if torch.isnan(res_loss) or torch.isinf(res_loss):
                continue

            # S¹ geodesic 타겟 손실
            geo_loss = geo_loss_fn(outputs['Z_out'], outputs['Psi_target'])

            total_loss = res_loss + lambda_geo * geo_loss
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
                    vr, _ = res_loss_fn(**o)
                    vg = geo_loss_fn(o['Z_out'], o['Psi_target'])
                    vl = vr + lambda_geo * vg
                    if not (torch.isnan(vl) or torch.isinf(vl)):
                        va += vl.item(); nv += 1
            vl_m = va / max(1, nv)
            best_val = min(best_val, vl_m)
            elapsed = time.time() - t_start
            print(f"      ep {ep+1}/{EPOCHS}: train={ep_loss/max(1,ep_n):.5f}, val={vl_m:.5f}, best={best_val:.5f}, {elapsed:.0f}s", flush=True)

    return model, best_val


def main():
    t_min, t_max = 100.0, 200.0
    num_points = 1000
    in_features = 128
    seeds = [42, 7, 123]

    out = []
    def log(msg):
        print(msg, flush=True); out.append(msg)

    log("=" * 72)
    log("  S¹ Geodesic 전체 통합 실험: L_tgt vs L_geo (Z_out 직접 적용)")
    log("=" * 72)
    log(f"  t∈[{t_min},{t_max}], {num_points}점, if={in_features}, ep={EPOCHS}")
    log(f"  seeds={seeds}")
    log(f"  L_tgt = ||atan2(sin(Δφ), cos(Δφ))||²")
    log(f"  L_geo = 1 - cos(arg(Z) - Ψ)")
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

        log("  [A: Baseline L_tgt] 훈련...")
        t0 = time.time()
        model_b, val_b = train_baseline(ds, s)
        f2z_b, f2nz_b, ratio_b, det_b, tot_z = eval_F2(model_b, ds)
        dt = time.time() - t0
        log(f"    val={val_b:.5f}, |F₂| ratio={ratio_b:.4f}, 검출={det_b}/{tot_z}, time={dt:.0f}s")
        baseline_results.append((val_b, ratio_b, det_b, tot_z, f2z_b, f2nz_b))

        log("  [B: S¹ Geodesic L_geo] 훈련...")
        t0 = time.time()
        model_s, val_s = train_s1_integrated(ds, s)
        f2z_s, f2nz_s, ratio_s, det_s, tot_z = eval_F2(model_s, ds)
        dt = time.time() - t0
        log(f"    val={val_s:.5f}, |F₂| ratio={ratio_s:.4f}, 검출={det_s}/{tot_z}, time={dt:.0f}s")
        s1_results.append((val_s, ratio_s, det_s, tot_z, f2z_s, f2nz_s))

    elapsed = time.time() - start

    log(f"\n{'='*72}")
    log("  요약")
    log(f"{'='*72}")

    br = [r[1] for r in baseline_results]
    bd = [r[2] for r in baseline_results]
    sr = [r[1] for r in s1_results]
    sd = [r[2] for r in s1_results]
    tot = baseline_results[0][3]

    log(f"  {'방식':<25} {'|F₂| ratio':>15} {'검출':>10}")
    log(f"  {'-'*25} {'-'*15} {'-'*10}")
    log(f"  {'Baseline (L_tgt)':<25} {np.mean(br):>7.4f}±{np.std(br):.4f} {np.mean(bd):>6.1f}/{tot}")
    log(f"  {'S¹ Geodesic (L_geo)':<25} {np.mean(sr):>7.4f}±{np.std(sr):.4f} {np.mean(sd):>6.1f}/{tot}")

    # 비율이 낮을수록 좋음 (영점에서 F₂가 작으면 검출 용이)
    ratio_change = (np.mean(sr) - np.mean(br)) / (np.mean(br) + 1e-12) * 100
    det_improve = np.mean(sd) - np.mean(bd)

    log(f"\n  |F₂| ratio 변화: {ratio_change:+.1f}% ({'개선↓' if ratio_change < 0 else '악화↑' if ratio_change > 0 else '무변화'})")
    log(f"  검출 변화: {det_improve:+.1f}개")

    # F₂ 절대값 비교
    f2z_b_mean = np.mean([r[4] for r in baseline_results])
    f2nz_b_mean = np.mean([r[5] for r in baseline_results])
    f2z_s_mean = np.mean([r[4] for r in s1_results])
    f2nz_s_mean = np.mean([r[5] for r in s1_results])
    log(f"\n  |F₂| 절대값:")
    log(f"    Baseline:  영점={f2z_b_mean:.6f}, 비영점={f2nz_b_mean:.6f}")
    log(f"    S¹ Geo:    영점={f2z_s_mean:.6f}, 비영점={f2nz_s_mean:.6f}")

    log(f"\n  총 실행 시간: {elapsed:.0f}s")

    # 판정 기준:
    # - ratio < 0.8*baseline OR 검출 > 1.5*baseline → 양성
    # - ratio > 1.2*baseline → 음성
    if (np.mean(sr) < np.mean(br) * 0.8) or (np.mean(sd) > max(np.mean(bd) * 1.5, np.mean(bd) + 5)):
        verdict = "양성: S¹ geodesic 통합이 F₂ 구조를 유의미하게 개선"
    elif np.mean(sr) > np.mean(br) * 1.2 and np.mean(sd) <= np.mean(bd):
        verdict = "음성: 기존 L_tgt가 우세"
    else:
        verdict = "중립: 유의미한 차이 없음"
    log(f"\n  판정: {verdict}")

    # S¹ readout 결과와 비교
    log(f"\n  참고 — 이전 S¹ readout 결과:")
    log(f"    Readout: ratio 1.054±0.103, 검출 13.7/463 (별도 head, Z_out 무영향)")
    log(f"    Integration: ratio {np.mean(sr):.3f}±{np.std(sr):.3f}, 검출 {np.mean(sd):.1f}/{tot} (L_tgt 대체)")

    os.makedirs(os.path.expanduser('~/Desktop/gdl_unified/results'), exist_ok=True)
    p = os.path.expanduser('~/Desktop/gdl_unified/results/s1_full_integration.txt')
    with open(p, 'w') as f:
        f.write('\n'.join(out))
    log(f"\n  저장: {p}")


if __name__ == '__main__':
    main()
