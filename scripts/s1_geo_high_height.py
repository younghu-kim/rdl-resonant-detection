"""
=============================================================================
[Project RDL] S¹ Geodesic 고높이 검증 실험 v2
=============================================================================
수학자 질문: "L_geo가 고높이(t>500)에서도 개선을 주는가?"

s1_full_integration에서 t∈[100,200] 양성 확인 (검출 4배, ratio -11%).
이를 t∈[500,600], t∈[1000,1100]에서 재검증.

비교:
  (A) Baseline: 표준 TotalResonanceLoss (L_tgt 포함)
  (B) S¹ Integration: L_tgt → L_geo 대체

평가 방법: is_near_zero 대신 극소값 탐색 + 매칭 (blind_zero_prediction 방식).
고높이에서 xi 언더플로로 is_near_zero 마스크가 작동하지 않으므로,
블라인드 예측 방식(전반부 훈련 → 후반부 dense grid 극소값 탐색)을 사용한다.

v1 → v2 수정:
  - eval_F2: is_near_zero 마스크 → 극소값 탐색 + 매칭
  - 시간 분할: 전반부 훈련, 후반부 블라인드 예측
  - dense grid에서 F₂ 평가 → recall/precision/F1 보고
"""

import sys, os, time
import numpy as np
import torch
import torch.nn as nn
from scipy.signal import argrelmin

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
IN_FEATURES = 128
DENSE_GRID_N = 2000
HIT_THRESHOLD = 0.3   # 전형적 영점 간격의 절반


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


def eval_F2_at_points(model, features):
    """올바른 F₂ 잔차 계산 (phi/psi/L_G 기반, Z_out.abs() 사용 금지)"""
    model.eval()
    with torch.enable_grad():
        X = features.clone().requires_grad_(True)
        out = model(X)
    phi = out["phi"].detach()
    psi = out["psi"].detach()
    L_G = out["L_G"].detach()
    phi_real = phi.to(dtype=PrecisionManager.REAL_DTYPE)
    rot = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
    psi_c = psi.to(dtype=PrecisionManager.COMPLEX_DTYPE)
    return (rot * (L_G - psi_c)).imag.mean(dim=-1).cpu().numpy()


def find_local_minima(t_grid, abs_f2, order=5):
    """극소값 위치 탐색."""
    indices = argrelmin(abs_f2, order=order)[0]
    return t_grid[indices], abs_f2[indices]


def match_predictions(predicted_t, actual_zeros, threshold):
    """예측 영점과 실제 영점 매칭."""
    n_pred = len(predicted_t)
    if n_pred == 0:
        return 0, 0.0, 0.0, 0.0

    actual = np.array(actual_zeros)
    n_actual = len(actual)
    if n_actual == 0:
        return 0, 0.0, 0.0, 0.0

    hit_count = 0
    matched_pred = set()
    for i, z in enumerate(actual):
        dists = np.abs(predicted_t - z)
        j = np.argmin(dists)
        if dists[j] < threshold and j not in matched_pred:
            hit_count += 1
            matched_pred.add(j)

    precision = len(matched_pred) / n_pred if n_pred > 0 else 0.0
    recall = hit_count / n_actual if n_actual > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
    return hit_count, precision, recall, f1


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
    seeds = [42, 7, 123, 314, 2024]

    out = []
    def log(msg):
        print(msg, flush=True); out.append(msg)

    log("=" * 72)
    log("  S¹ Geodesic 고높이 검증 v2: L_tgt vs L_geo at t>500")
    log("=" * 72)
    log(f"  구간: {ranges}")
    log(f"  K={IN_FEATURES}, seeds={seeds}, ep={EPOCHS}")
    log(f"  평가: 전반부 훈련 → 후반부 블라인드 예측 (극소값 탐색 + 매칭)")
    log(f"  참조: t∈[100,200] 결과 — Baseline 검출 5.7/463, L_geo 검출 22.7/463")
    log(f"  v1 버그 수정: is_near_zero 마스크 → 극소값 탐색 (xi 언더플로 우회)")
    log("")

    start = time.time()
    all_results = {}

    for t_min, t_max in ranges:
        log(f"\n{'='*72}")
        log(f"  구간: t∈[{t_min}, {t_max}]")
        log(f"{'='*72}")

        # 영점 계산
        zeros_list = compute_zeros_in_range(t_min, t_max, dps=25)
        n_zeros = len(zeros_list)
        density = n_zeros / (t_max - t_min)
        log(f"  영점 수: {n_zeros}, 밀도: {density:.3f}/unit")

        # 시간 분할: 전반부 훈련, 후반부 블라인드 예측
        t_mid = (t_min + t_max) / 2
        train_zeros = [z for z in zeros_list if z <= t_mid]
        test_zeros = [z for z in zeros_list if z > t_mid]
        log(f"  학습: t∈[{t_min},{t_mid}] ({len(train_zeros)}개 영점)")
        log(f"  예측: t∈[{t_mid},{t_max}] ({len(test_zeros)}개 영점)")

        # 학습 범위 캐시
        train_cache_path = os.path.expanduser(
            f'~/Desktop/gdl_unified/outputs/xi_cache_t{t_min}-{t_mid}_n1000.pt'
        )
        train_cache = get_or_build_cache(train_cache_path, t_min, t_mid, 1000,
                                          mp_dps=50, zeros_list=train_zeros)
        train_ds = XiFeatureDataset(train_cache, in_features=IN_FEATURES)

        # 블라인드 예측용 밀집 격자
        dense_grid = np.linspace(t_mid, t_max, DENSE_GRID_N)
        dense_grid_list = dense_grid.tolist()

        baseline_results = []
        s1_results = []

        for s in seeds:
            log(f"\n  --- seed={s} ---")

            # (A) Baseline 훈련 + 평가
            log("    [A: Baseline L_tgt] 훈련...")
            t0 = time.time()
            model_b, val_b = train_baseline(train_ds, s)

            # 블라인드 예측 평가
            dense_features = train_ds.get_features_at_t(dense_grid_list)
            f2_b = eval_F2_at_points(model_b, dense_features)
            abs_f2_b = np.abs(f2_b)
            pred_t_b, _ = find_local_minima(dense_grid, abs_f2_b, order=5)
            hits_b, prec_b, rec_b, f1_b = match_predictions(
                pred_t_b, test_zeros, HIT_THRESHOLD)
            dt = time.time() - t0
            log(f"    val={val_b:.5f}, 예측={len(pred_t_b)}개, "
                f"적중={hits_b}/{len(test_zeros)}, "
                f"P={prec_b:.3f} R={rec_b:.3f} F1={f1_b:.3f}, {dt:.0f}s")
            baseline_results.append({
                'val': val_b, 'hits': hits_b, 'prec': prec_b,
                'rec': rec_b, 'f1': f1_b, 'n_pred': len(pred_t_b),
            })

            # (B) S¹ Geodesic 훈련 + 평가
            log("    [B: S¹ Geodesic L_geo] 훈련...")
            t0 = time.time()
            model_s, val_s = train_s1_integrated(train_ds, s)

            dense_features = train_ds.get_features_at_t(dense_grid_list)
            f2_s = eval_F2_at_points(model_s, dense_features)
            abs_f2_s = np.abs(f2_s)
            pred_t_s, _ = find_local_minima(dense_grid, abs_f2_s, order=5)
            hits_s, prec_s, rec_s, f1_s = match_predictions(
                pred_t_s, test_zeros, HIT_THRESHOLD)
            dt = time.time() - t0
            log(f"    val={val_s:.5f}, 예측={len(pred_t_s)}개, "
                f"적중={hits_s}/{len(test_zeros)}, "
                f"P={prec_s:.3f} R={rec_s:.3f} F1={f1_s:.3f}, {dt:.0f}s")
            s1_results.append({
                'val': val_s, 'hits': hits_s, 'prec': prec_s,
                'rec': rec_s, 'f1': f1_s, 'n_pred': len(pred_t_s),
            })

        # 구간 요약
        n_test = len(test_zeros)
        log(f"\n  구간 요약: t∈[{t_min},{t_max}], 테스트 영점 {n_test}개")
        log(f"  {'방식':<25} {'Recall':>10} {'Precision':>12} {'F1':>10} {'적중':>10}")
        log(f"  {'-'*25} {'-'*10} {'-'*12} {'-'*10} {'-'*10}")

        br = [r['rec'] for r in baseline_results]
        bp = [r['prec'] for r in baseline_results]
        bf = [r['f1'] for r in baseline_results]
        bh = [r['hits'] for r in baseline_results]
        sr = [r['rec'] for r in s1_results]
        sp = [r['prec'] for r in s1_results]
        sf = [r['f1'] for r in s1_results]
        sh = [r['hits'] for r in s1_results]

        log(f"  {'Baseline (L_tgt)':<25} {np.mean(br):>5.3f}±{np.std(br):.3f} "
            f"{np.mean(bp):>7.3f}±{np.std(bp):.3f} "
            f"{np.mean(bf):>5.3f}±{np.std(bf):.3f} "
            f"{np.mean(bh):>5.1f}/{n_test}")
        log(f"  {'S¹ Geodesic (L_geo)':<25} {np.mean(sr):>5.3f}±{np.std(sr):.3f} "
            f"{np.mean(sp):>7.3f}±{np.std(sp):.3f} "
            f"{np.mean(sf):>5.3f}±{np.std(sf):.3f} "
            f"{np.mean(sh):>5.1f}/{n_test}")

        rec_diff = np.mean(sr) - np.mean(br)
        f1_diff = np.mean(sf) - np.mean(bf)
        log(f"  Recall 차이: {rec_diff:+.3f}, F1 차이: {f1_diff:+.3f}")

        all_results[(t_min, t_max)] = {
            'baseline': baseline_results,
            's1': s1_results,
            'n_zeros': n_zeros,
            'n_test': n_test,
            'density': density,
        }

    elapsed = time.time() - start

    # 전체 요약
    log(f"\n{'='*72}")
    log("  전체 요약: L_geo 고높이 효과 (블라인드 예측)")
    log(f"{'='*72}")
    log(f"  {'구간':<20} {'B Recall':>10} {'S¹ Recall':>12} {'Δ Recall':>10} {'B F1':>8} {'S¹ F1':>8} {'Δ F1':>8}")
    log(f"  {'-'*20} {'-'*10} {'-'*12} {'-'*10} {'-'*8} {'-'*8} {'-'*8}")

    all_rec_diffs = []
    all_f1_diffs = []

    for (t_min, t_max), res in all_results.items():
        br_m = np.mean([r['rec'] for r in res['baseline']])
        sr_m = np.mean([r['rec'] for r in res['s1']])
        bf_m = np.mean([r['f1'] for r in res['baseline']])
        sf_m = np.mean([r['f1'] for r in res['s1']])
        label = f"t∈[{int(t_min)},{int(t_max)}]"
        log(f"  {label:<20} {br_m:>10.3f} {sr_m:>12.3f} {sr_m-br_m:>+10.3f} "
            f"{bf_m:>8.3f} {sf_m:>8.3f} {sf_m-bf_m:>+8.3f}")
        all_rec_diffs.append(sr_m - br_m)
        all_f1_diffs.append(sf_m - bf_m)

    avg_rec_diff = np.mean(all_rec_diffs)
    avg_f1_diff = np.mean(all_f1_diffs)

    log(f"\n  고높이 평균 Recall 차이: {avg_rec_diff:+.3f}")
    log(f"  고높이 평균 F1 차이: {avg_f1_diff:+.3f}")

    # 판정: recall/F1 기준
    if avg_rec_diff > 0.05 or avg_f1_diff > 0.03:
        verdict = "양성: L_geo가 고높이에서도 유의미한 개선을 제공"
    elif avg_rec_diff < -0.05 or avg_f1_diff < -0.03:
        verdict = "음성: L_geo 효과가 고높이에서 소실 또는 악화"
    else:
        verdict = "중립: 고높이에서 L_geo 효과가 불분명 (±5%p 이내)"

    log(f"\n  판정: {verdict}")
    log(f"  5시드, 2구간. 헌법 §3 준수 (≥5시드).")
    log(f"  총 실행 시간: {elapsed:.0f}s ({elapsed/3600:.1f}h)")

    os.makedirs(os.path.expanduser('~/Desktop/gdl_unified/results'), exist_ok=True)
    p = os.path.expanduser('~/Desktop/gdl_unified/results/s1_geo_high_height.txt')
    with open(p, 'w') as f:
        f.write('\n'.join(out))
    log(f"\n  저장: {p}")


if __name__ == '__main__':
    main()
