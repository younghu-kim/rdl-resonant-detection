"""
=============================================================================
[Project RDL] S¹ Geodesic 고높이 검증 실험
=============================================================================
수학자 질문: "L_geo가 고높이(t>500)에서도 개선을 주는가?"

s1_full_integration에서 t∈[100,200] 양성 확인 (검출 4배, ratio -11%).
이를 t∈[500,600], t∈[1000,1100]에서 재검증.

비교:
  (A) Baseline: 표준 TotalResonanceLoss (L_tgt 포함)
  (B) S¹ Integration: L_tgt → L_geo 대체

핵심: 고높이에서 영점 밀도 증가 + xi 언더플로 → L_geo 우위가 유지되는가?
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

HIDDEN = 64
EPOCHS = 100
LR = 1e-3
BATCH = 32
IN_FEATURES = 128


class GeodesicTargetLoss(nn.Module):
    """S¹ geodesic 손실: L_geo = 1 - cos(arg(Z_out) - Psi_target)"""
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
        delta = current_phase - Psi_target
        loss = 1.0 - torch.cos(delta)
        if self.reduction == 'mean':
            return loss.mean()
        elif self.reduction == 'sum':
            return loss.sum()
        return loss


def eval_F2(model, dataset, batch_size=64):
    """올바른 F₂ 잔차 기반 영점 탐지 평가 (phi/psi/L_G 기반, Z_out.abs() 사용 금지)"""
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
    """Baseline: 표준 TotalResonanceLoss (L_tgt 포함)"""
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
            print(f"      ep {ep+1}/{EPOCHS}: train={ep_loss/max(1,ep_n):.5f}, val={vl_m:.5f}, {elapsed:.0f}s", flush=True)

    return model, best_val


def train_s1_integrated(dataset, seed):
    """S¹ 통합: L_tgt → GeodesicTargetLoss 대체"""
    torch.manual_seed(seed); np.random.seed(seed)
    train_loader, val_loader = make_loaders(dataset, seed)

    model = MasterResonantNetwork(
        in_features=dataset.in_features, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    res_loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=0.0, lambda_pqo=0.5,
        pqo_mode='cos2'
    )
    geo_loss_fn = GeodesicTargetLoss(reduction='mean')
    lambda_geo = 1.0

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
            res_loss, _ = res_loss_fn(**outputs)
            if torch.isnan(res_loss) or torch.isinf(res_loss):
                continue
            geo_loss = geo_loss_fn(outputs['Z_out'], outputs['Psi_target'])
            total_loss = res_loss + lambda_geo * geo_loss
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
                    vr, _ = res_loss_fn(**o)
                    vg = geo_loss_fn(o['Z_out'], o['Psi_target'])
                    vl = vr + lambda_geo * vg
                    if not (torch.isnan(vl) or torch.isinf(vl)):
                        va += vl.item(); nv += 1
            vl_m = va / max(1, nv)
            best_val = min(best_val, vl_m)
            elapsed = time.time() - t_start
            print(f"      ep {ep+1}/{EPOCHS}: train={ep_loss/max(1,ep_n):.5f}, val={vl_m:.5f}, {elapsed:.0f}s", flush=True)

    return model, best_val


def main():
    ranges = [
        (500.0, 600.0),
        (1000.0, 1100.0),
    ]
    seeds = [42, 7, 123]

    out = []
    def log(msg):
        print(msg, flush=True); out.append(msg)

    log("=" * 72)
    log("  S¹ Geodesic 고높이 검증: L_tgt vs L_geo at t>500")
    log("=" * 72)
    log(f"  구간: {ranges}")
    log(f"  K={IN_FEATURES}, seeds={seeds}, ep={EPOCHS}")
    log(f"  참조: t∈[100,200] 결과 — Baseline 검출 5.7/463, L_geo 검출 22.7/463")
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
            f'~/Desktop/gdl_unified/outputs/xi_cache_t{t_min}-{t_max}_n1000.pt'
        )
        cache_data = get_or_build_cache(cache_path, t_min, t_max, 1000,
                                         mp_dps=50, zeros_list=zeros_list)
        ds = XiFeatureDataset(cache_data, in_features=IN_FEATURES)

        baseline_results = []
        s1_results = []

        for s in seeds:
            log(f"\n  --- seed={s} ---")

            log("    [A: Baseline L_tgt] 훈련...")
            t0 = time.time()
            model_b, val_b = train_baseline(ds, s)
            f2z_b, f2nz_b, ratio_b, det_b, tot_z = eval_F2(model_b, ds)
            dt = time.time() - t0
            log(f"    val={val_b:.5f}, |F₂| ratio={ratio_b:.4f}, 검출={det_b}/{tot_z}, time={dt:.0f}s")
            baseline_results.append((val_b, ratio_b, det_b, tot_z, f2z_b, f2nz_b))

            log("    [B: S¹ Geodesic L_geo] 훈련...")
            t0 = time.time()
            model_s, val_s = train_s1_integrated(ds, s)
            f2z_s, f2nz_s, ratio_s, det_s, tot_z = eval_F2(model_s, ds)
            dt = time.time() - t0
            log(f"    val={val_s:.5f}, |F₂| ratio={ratio_s:.4f}, 검출={det_s}/{tot_z}, time={dt:.0f}s")
            s1_results.append((val_s, ratio_s, det_s, tot_z, f2z_s, f2nz_s))

        # 구간 요약
        br = [r[1] for r in baseline_results]
        bd = [r[2] for r in baseline_results]
        sr = [r[1] for r in s1_results]
        sd = [r[2] for r in s1_results]
        tot = baseline_results[0][3]

        log(f"\n  구간 요약: t∈[{t_min},{t_max}]")
        log(f"  {'방식':<25} {'|F₂| ratio':>15} {'검출':>10}")
        log(f"  {'-'*25} {'-'*15} {'-'*10}")
        log(f"  {'Baseline (L_tgt)':<25} {np.mean(br):>7.4f}±{np.std(br):.4f} {np.mean(bd):>6.1f}/{tot}")
        log(f"  {'S¹ Geodesic (L_geo)':<25} {np.mean(sr):>7.4f}±{np.std(sr):.4f} {np.mean(sd):>6.1f}/{tot}")

        ratio_change = (np.mean(sr) - np.mean(br)) / (np.mean(br) + 1e-12) * 100
        det_improve = np.mean(sd) - np.mean(bd)
        log(f"  |F₂| ratio 변화: {ratio_change:+.1f}%")
        log(f"  검출 변화: {det_improve:+.1f}개")

        all_results[(t_min, t_max)] = {
            'baseline': baseline_results,
            's1': s1_results,
            'n_zeros': n_zeros,
            'density': density,
        }

    elapsed = time.time() - start

    # 전체 요약
    log(f"\n{'='*72}")
    log("  전체 요약: L_geo 고높이 효과")
    log(f"{'='*72}")
    log(f"  {'구간':<20} {'Baseline 검출':>15} {'L_geo 검출':>15} {'검출 배율':>10} {'ratio 변화':>12}")
    log(f"  {'-'*20} {'-'*15} {'-'*15} {'-'*10} {'-'*12}")

    # 기존 t∈[100,200] 참조 데이터
    log(f"  {'t∈[100,200]*':<20} {'5.7/463':>15} {'22.7/463':>15} {'4.0x':>10} {'-11.0%':>12}")

    for (t_min, t_max), res in all_results.items():
        bd_m = np.mean([r[2] for r in res['baseline']])
        sd_m = np.mean([r[2] for r in res['s1']])
        tot = res['baseline'][0][3]
        ratio_b = np.mean([r[1] for r in res['baseline']])
        ratio_s = np.mean([r[1] for r in res['s1']])
        ratio_chg = (ratio_s - ratio_b) / (ratio_b + 1e-12) * 100
        det_mult = sd_m / max(bd_m, 0.1)
        log(f"  {'t∈['+str(int(t_min))+','+str(int(t_max))+']':<20} {bd_m:>6.1f}/{tot:<8} {sd_m:>6.1f}/{tot:<8} {det_mult:>6.1f}x {ratio_chg:>+10.1f}%")

    log(f"\n  * t∈[100,200]은 s1_full_integration.txt에서 인용 (동일 설계, ep=150)")

    # 종합 판정
    all_det_improve = []
    all_ratio_improve = []
    for res in all_results.values():
        bd_m = np.mean([r[2] for r in res['baseline']])
        sd_m = np.mean([r[2] for r in res['s1']])
        ratio_b = np.mean([r[1] for r in res['baseline']])
        ratio_s = np.mean([r[1] for r in res['s1']])
        all_det_improve.append(sd_m - bd_m)
        all_ratio_improve.append((ratio_s - ratio_b) / (ratio_b + 1e-12))

    avg_det_improve = np.mean(all_det_improve)
    avg_ratio_improve = np.mean(all_ratio_improve) * 100

    log(f"\n  고높이 평균 검출 개선: {avg_det_improve:+.1f}개")
    log(f"  고높이 평균 ratio 변화: {avg_ratio_improve:+.1f}%")

    if avg_det_improve > 3 or avg_ratio_improve < -5:
        verdict = "양성: L_geo가 고높이에서도 유의미한 개선을 제공"
    elif avg_det_improve < -3 or avg_ratio_improve > 10:
        verdict = "음성: L_geo 효과가 고높이에서 소실"
    else:
        verdict = "중립: 고높이에서 L_geo 효과가 불분명"

    log(f"\n  판정: {verdict}")
    log(f"  총 실행 시간: {elapsed:.0f}s")

    os.makedirs(os.path.expanduser('~/Desktop/gdl_unified/results'), exist_ok=True)
    p = os.path.expanduser('~/Desktop/gdl_unified/results/s1_geo_high_height.txt')
    with open(p, 'w') as f:
        f.write('\n'.join(out))
    log(f"\n  저장: {p}")


if __name__ == '__main__':
    main()
