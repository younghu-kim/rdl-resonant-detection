"""
=============================================================================
[Project RDL] 고높이 스케일링 실험 v2
=============================================================================
t∈[100,200] 에서 검증된 F₂ 블라인드 영점 검출을
t∈[500,600], t∈[1000,1100] 에서도 테스트.

v1 대비 수정:
  1) F₂ 평가: Z_out.abs() (잘못됨) → (rot*(L_G-psi_c)).imag (올바른 F₂ 잔차)
  2) 영점 검출: is_near_zero 마스크 → 극소값 탐색 + 매칭 (blind_zero_prediction 방식)
  3) xi 정밀도 언더플로: 학습에서 xi 타겟을 사용하지 않으므로 무관.
     is_near_zero도 사용하지 않으므로 amp 언더플로도 무관.
"""

import sys, os, time, copy, math
import numpy as np
import torch
from scipy.signal import argrelmin

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
from gdl.rdl.constants import PrecisionManager
from gdl.rdl.pipeline.xi_feature_dataset import (
    get_or_build_cache, compute_zeros_in_range, XiFeatureDataset,
    get_xi_feature_dataloaders
)
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss

PrecisionManager.setup_precision()

HIDDEN = 64
EPOCHS = 150
LR = 1e-3
BATCH = 32
DENSE_GRID_N = 2000
HIT_THRESHOLD = 0.3   # 전형적 영점 간격의 절반


def eval_F2(model, features):
    """올바른 F₂ 잔차 계산 — blind_zero_prediction.py 와 동일."""
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


def train_one(train_loader, val_loader, in_features, seed, label):
    """단일 모델 학습 — early stopping 포함."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = MasterResonantNetwork(
        in_features=in_features, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5,
        pqo_mode='cos2'
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    best_val = float('inf')
    best_state = None
    best_ep = 0
    t_start = time.time()

    for ep in range(1, EPOCHS + 1):
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

        # 매 에폭 검증 (early stopping)
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
        if vl_m < best_val:
            best_val = vl_m
            best_state = copy.deepcopy(model.state_dict())
            best_ep = ep

        if ep % 30 == 0 or ep == EPOCHS:
            elapsed = time.time() - t_start
            print(f"      {label} ep {ep}/{EPOCHS}: "
                  f"train={ep_loss/max(1,ep_n):.5f}, val={vl_m:.5f}, "
                  f"best={best_val:.5f}@{best_ep}, {elapsed:.0f}s", flush=True)

    if best_state:
        model.load_state_dict(best_state)
    return model, best_val, best_ep, time.time() - t_start


def main():
    ranges = [
        (100.0, 200.0),
        (500.0, 600.0),
        (1000.0, 1100.0),
    ]
    num_points = 1000
    in_features_list = [128, 256]
    seeds = [42, 7, 123]

    # dps 스케일링: 높이에 비례하여 정밀도 증가
    def get_dps(t_max):
        # xi(1/2+it) ~ exp(-πt/4), 안전 마진 포함
        return max(50, int(t_max * 0.5) + 30)

    out = []
    def log(msg):
        print(msg, flush=True); out.append(msg)

    log("=" * 72)
    log("  고높이 스케일링 실험 v2 (올바른 F₂ 잔차 사용)")
    log("=" * 72)
    log(f"  구간: {ranges}")
    log(f"  in_features: {in_features_list}, seeds={seeds}, ep={EPOCHS}")
    log(f"  밀집 격자: {DENSE_GRID_N}점, 적중 임계: {HIT_THRESHOLD}")
    log("")

    start = time.time()
    all_results = {}

    for t_min, t_max in ranges:
        log(f"\n{'='*72}")
        log(f"  구간: t∈[{t_min}, {t_max}]")
        log(f"{'='*72}")

        # 영점 계산
        dps = get_dps(t_max)
        zeros_list = compute_zeros_in_range(t_min, t_max, dps=25)
        n_zeros = len(zeros_list)
        density = n_zeros / (t_max - t_min)
        log(f"  영점 수: {n_zeros}, 밀도: {density:.3f}/unit")

        # 학습/검증 분할: 전반부 학습, 후반부 블라인드 예측
        t_mid = (t_min + t_max) / 2
        train_zeros = [z for z in zeros_list if z <= t_mid]
        test_zeros = [z for z in zeros_list if z > t_mid]
        log(f"  학습: t∈[{t_min},{t_mid}] ({len(train_zeros)}개 영점)")
        log(f"  예측: t∈[{t_mid},{t_max}] ({len(test_zeros)}개 영점)")

        # 캐시 구축 (dps 스케일링)
        cache_path = os.path.expanduser(
            f'~/Desktop/gdl_unified/outputs/xi_cache_t{t_min}-{t_max}_n{num_points}_dps{dps}.pt'
        )
        cache_data = get_or_build_cache(cache_path, t_min, t_max, num_points,
                                         mp_dps=dps, zeros_list=zeros_list)

        # 밀집 격자 (블라인드 영역)
        dense_grid = np.linspace(t_mid, t_max, DENSE_GRID_N)
        dense_grid_list = dense_grid.tolist()

        for K in in_features_list:
            log(f"\n  --- K={K} ---")

            # 학습 범위 데이터셋
            train_cache_path = os.path.expanduser(
                f'~/Desktop/gdl_unified/outputs/xi_cache_t{t_min}-{t_mid}_n{num_points}_dps{dps}.pt'
            )
            train_cache = get_or_build_cache(train_cache_path, t_min, t_mid,
                                              num_points, mp_dps=dps,
                                              zeros_list=train_zeros)
            dataset = XiFeatureDataset(train_cache, in_features=K)

            val_size = int(len(dataset) * 0.2)
            train_size = len(dataset) - val_size
            train_ds, val_ds = torch.utils.data.random_split(
                dataset, [train_size, val_size],
                generator=torch.Generator().manual_seed(0)
            )
            train_loader = torch.utils.data.DataLoader(
                train_ds, batch_size=BATCH, shuffle=True, drop_last=True)
            val_loader = torch.utils.data.DataLoader(
                val_ds, batch_size=BATCH, shuffle=False)

            hits_list = []; prec_list = []; rec_list = []; f1_list = []

            for s in seeds:
                label = f"t[{int(t_min)},{int(t_max)}]/K{K}/s{s}"
                model, bv, be, dt = train_one(
                    train_loader, val_loader, K, s, label)

                # 밀집 격자에서 F₂ 평가
                dense_features = dataset.get_features_at_t(dense_grid_list)
                f2_dense = eval_F2(model, dense_features)
                abs_f2 = np.abs(f2_dense)

                # 극소값 탐색 + 매칭
                pred_t, pred_vals = find_local_minima(dense_grid, abs_f2, order=5)
                hits, prec, rec, f1 = match_predictions(
                    pred_t, test_zeros, HIT_THRESHOLD)

                log(f"    seed={s}: val={bv:.5f}@{be}, "
                    f"예측={len(pred_t)}개, 적중={hits}/{len(test_zeros)}, "
                    f"P={prec:.3f} R={rec:.3f} F1={f1:.3f}, {dt:.0f}s")

                hits_list.append(hits)
                prec_list.append(prec)
                rec_list.append(rec)
                f1_list.append(f1)

            key = f"t[{int(t_min)},{int(t_max)}]_K{K}"
            all_results[key] = {
                'n_zeros': n_zeros,
                'n_test': len(test_zeros),
                'density': density,
                'recall_mean': np.mean(rec_list),
                'recall_std': np.std(rec_list),
                'precision_mean': np.mean(prec_list),
                'f1_mean': np.mean(f1_list),
                'hits_mean': np.mean(hits_list),
            }
            log(f"  → K={K}: recall={np.mean(rec_list):.3f}±{np.std(rec_list):.3f}, "
                f"precision={np.mean(prec_list):.3f}, F1={np.mean(f1_list):.3f}")

    elapsed = time.time() - start

    log(f"\n{'='*72}")
    log("  스케일링 요약")
    log(f"{'='*72}")
    log(f"  {'구간_K':<25} {'영점':>5} {'밀도':>7} {'recall':>15} {'prec':>7} {'F1':>7}")
    log(f"  {'-'*25} {'-'*5} {'-'*7} {'-'*15} {'-'*7} {'-'*7}")

    for key, r in all_results.items():
        log(f"  {key:<25} {r['n_test']:>5} {r['density']:>7.3f} "
            f"{r['recall_mean']:>7.3f}±{r['recall_std']:.3f} "
            f"{r['precision_mean']:>7.3f} {r['f1_mean']:>7.3f}")

    log(f"\n  총 실행 시간: {elapsed:.0f}s")

    # 판정
    base_keys = [k for k in all_results if 't[100,200]' in k]
    if base_keys:
        base_recall = np.mean([all_results[k]['recall_mean'] for k in base_keys])
        high_keys = [k for k in all_results if 't[100,200]' not in k]
        high_recall = np.mean([all_results[k]['recall_mean'] for k in high_keys]) if high_keys else 0

        if high_recall >= base_recall * 0.8:
            log("\n  판정: 양성 — 고높이에서도 영점 검출 성능 유지")
        elif high_recall >= base_recall * 0.5:
            log("\n  판정: 중립 — 고높이에서 성능 저하 관찰, K 스케일링 필요 가능")
        else:
            log("\n  판정: 음성 — 고높이에서 성능 대폭 저하, K 스케일링 법칙 도출 필요")
    else:
        log("\n  판정: 데이터 부족")

    p = os.path.expanduser('~/Desktop/gdl_unified/results/high_height_scaling.txt')
    with open(p, 'w') as f:
        f.write('\n'.join(out))
    log(f"\n  저장: {p}")


if __name__ == '__main__':
    main()
