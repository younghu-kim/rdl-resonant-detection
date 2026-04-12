"""
=============================================================================
[Project RDL] Kuramoto + S¹ 통합 실험: L_geo가 r 궤적을 바꾸는가?
=============================================================================
수학자 질문:
  (1) L_geo 사용 시 r 궤적이 달라지는가?
  (2) cos 손실이 von Mises 사전분포를 암묵적으로 부과 → 더 빠른 동기화 예상
  (3) von Mises κ와 검출 성능 상관 (선행 데이터 수집)

설계:
  (A) Baseline: 표준 TotalResonanceLoss (L_tgt 포함) — kuramoto 실험 재현
  (B) S¹ Integration: L_tgt → L_geo 대체 — 매 에포크 r/std(φ) 추적
  두 방식의 r 궤적, std(φ) 궤적, κ 궤적을 직접 비교.

참조 데이터 (kuramoto_order_parameter.txt):
  Baseline r: 0.995 → 0.9994, std ratio 3.43x, t∈[100,200]
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
EPOCHS = 100
LR = 1e-3
BATCH = 32
IN_FEATURES = 128
SEEDS = [42, 7, 123]
T_RANGE = (100.0, 200.0)  # 기존 Kuramoto 결과와 직접 비교


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


def collect_gauge_phases(model, dataset, batch_size=64):
    """전체 데이터셋에 대해 게이지 위상 φ를 수집."""
    model.eval()
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)
    phi_list = []
    with torch.enable_grad():
        for X_batch, _ in loader:
            X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            out = model(X_in)
            phi_list.append(out["phi"].detach().cpu().numpy())
    return np.concatenate(phi_list, axis=0)  # [N, hidden_features]


def compute_kuramoto_metrics(phi_all):
    """Kuramoto 순서 매개변수, std, von Mises κ 계산."""
    exp_iphi = np.exp(1j * phi_all)  # [N, H]
    r_per_dim = np.abs(exp_iphi.mean(axis=0))  # [H]
    r_global = r_per_dim.mean()

    std_per_dim = phi_all.std(axis=0)  # [H]
    std_global = std_per_dim.mean()

    mean_phi = np.angle(exp_iphi.mean(axis=0)).mean()

    # von Mises κ 추정: κ ≈ 1/σ² (강결합 근사)
    kappa = 1.0 / (std_global**2 + 1e-12)

    # 정밀한 κ 추정: r ≈ I_1(κ)/I_0(κ) 역함수 근사
    # 강결합(r>0.9)에서 κ ≈ 1/(1-r²) 근사 사용
    kappa_from_r = 1.0 / (1.0 - r_global**2 + 1e-12) if r_global < 1.0 else float('inf')

    return {
        "r_global": r_global,
        "std_global": std_global,
        "mean_phi": mean_phi,
        "kappa_from_std": kappa,
        "kappa_from_r": kappa_from_r,
    }


def eval_F2(model, dataset, batch_size=64):
    """F₂ 잔차 기반 영점 탐지 평가 (phi/psi/L_G 기반)"""
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

    return ratio, detected, total_zeros


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


def train_with_kuramoto(dataset, seed, mode="baseline"):
    """학습 + 매 에포크 Kuramoto 추적.

    mode="baseline": 표준 TotalResonanceLoss
    mode="s1_geo":   L_tgt→L_geo 대체
    """
    torch.manual_seed(seed); np.random.seed(seed)
    train_loader, val_loader = make_loaders(dataset, seed)

    model = MasterResonantNetwork(
        in_features=dataset.in_features, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)

    if mode == "baseline":
        loss_fn = TotalResonanceLoss(
            lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5,
            pqo_mode='cos2'
        )
        geo_loss_fn = None
    else:
        loss_fn = TotalResonanceLoss(
            lambda_res=1.0, lambda_curv=0.1, lambda_tgt=0.0, lambda_pqo=0.5,
            pqo_mode='cos2'
        )
        geo_loss_fn = GeodesicTargetLoss(reduction='mean')

    # 에포크별 기록
    history = {"r": [], "std": [], "kappa_std": [], "kappa_r": [], "loss": []}

    t_start = time.time()

    # 초기 상태 (학습 전)
    m0 = compute_kuramoto_metrics(collect_gauge_phases(model, dataset))
    history["r"].append(m0["r_global"])
    history["std"].append(m0["std_global"])
    history["kappa_std"].append(m0["kappa_from_std"])
    history["kappa_r"].append(m0["kappa_from_r"])
    history["loss"].append(float('nan'))

    for ep in range(EPOCHS):
        model.train()
        ep_loss = 0.0; ep_n = 0
        for X_batch, _ in train_loader:
            X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            optimizer.zero_grad(set_to_none=True)
            outputs = model(X_in)

            if mode == "baseline":
                total_loss, _ = loss_fn(**outputs)
            else:
                res_loss, _ = loss_fn(**outputs)
                if torch.isnan(res_loss) or torch.isinf(res_loss):
                    continue
                geo_loss = geo_loss_fn(outputs['Z_out'], outputs['Psi_target'])
                total_loss = res_loss + 1.0 * geo_loss

            if torch.isnan(total_loss) or torch.isinf(total_loss):
                continue
            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            ep_loss += total_loss.item(); ep_n += 1

        avg_loss = ep_loss / max(1, ep_n)

        # 매 에포크 Kuramoto 측정
        metrics = compute_kuramoto_metrics(collect_gauge_phases(model, dataset))
        history["r"].append(metrics["r_global"])
        history["std"].append(metrics["std_global"])
        history["kappa_std"].append(metrics["kappa_from_std"])
        history["kappa_r"].append(metrics["kappa_from_r"])
        history["loss"].append(avg_loss)

        if (ep + 1) % 10 == 0:
            elapsed = time.time() - t_start
            print(f"    [{mode}] ep {ep+1}/{EPOCHS}: loss={avg_loss:.5f}, "
                  f"r={metrics['r_global']:.4f}, std={metrics['std_global']:.4f}, "
                  f"κ(std)={metrics['kappa_from_std']:.0f}, {elapsed:.0f}s", flush=True)

    elapsed = time.time() - t_start

    # 최종 F2 평가
    ratio, detected, total_zeros = eval_F2(model, dataset)

    return history, elapsed, ratio, detected, total_zeros


def main():
    results_path = os.path.expanduser("~/Desktop/gdl_unified/results/kuramoto_s1_integration.txt")
    out = []

    def log(msg):
        print(msg, flush=True)
        out.append(msg)

    log("=" * 72)
    log("  Kuramoto + S¹ 통합: L_geo가 위상 동기화 궤적을 바꾸는가?")
    log("=" * 72)
    log(f"  수학자 예측: L_geo = 1-cos(Δφ)는 von Mises 음의 로그우도")
    log(f"  → 동기화 가속 효과 예상 (더 빠른 r→1, 더 큰 κ)")
    log(f"  HIDDEN={HIDDEN}, EPOCHS={EPOCHS}, LR={LR}, BATCH={BATCH}")
    log(f"  t∈[{T_RANGE[0]},{T_RANGE[1]}], seeds={SEEDS}")
    log(f"  참조: Baseline r 0.995→0.9994, std ratio 3.43x (기존 실험)")
    log("")

    t_lo, t_hi = T_RANGE
    zeros = compute_zeros_in_range(t_lo, t_hi)
    n_zeros = len(zeros)
    log(f"  [Data] t∈[{t_lo},{t_hi}]: {n_zeros}개 영점")

    cache_path = os.path.expanduser(
        f"~/Desktop/gdl_unified/cache/xi_cache_t{t_lo}-{t_hi}_n1000.pt"
    )
    cache = get_or_build_cache(cache_path, t_lo, t_hi, 1000)
    dataset = XiFeatureDataset(cache, IN_FEATURES)
    log(f"  [Data] {len(dataset)}개 점, in_features={dataset.in_features}")

    start = time.time()
    all_baseline = []
    all_s1 = []

    for seed in SEEDS:
        log(f"\n{'─' * 72}")
        log(f"  seed={seed}")
        log(f"{'─' * 72}")

        # (A) Baseline
        log(f"\n  [A] Baseline (L_tgt) 학습...")
        h_b, t_b, ratio_b, det_b, tot_z = train_with_kuramoto(dataset, seed, mode="baseline")
        log(f"    완료: {t_b:.0f}s, F₂ ratio={ratio_b:.4f}, 검출={det_b}/{tot_z}")
        all_baseline.append((h_b, ratio_b, det_b, tot_z))

        # (B) S¹ Geodesic
        log(f"\n  [B] S¹ Geodesic (L_geo) 학습...")
        h_s, t_s, ratio_s, det_s, _ = train_with_kuramoto(dataset, seed, mode="s1_geo")
        log(f"    완료: {t_s:.0f}s, F₂ ratio={ratio_s:.4f}, 검출={det_s}/{tot_z}")
        all_s1.append((h_s, ratio_s, det_s, tot_z))

        # 시드별 궤적 비교
        r_b = np.array(h_b["r"])
        r_s = np.array(h_s["r"])
        std_b = np.array(h_b["std"])
        std_s = np.array(h_s["std"])
        k_b = np.array(h_b["kappa_std"])
        k_s = np.array(h_s["kappa_std"])

        log(f"\n  시드={seed} 궤적 비교:")
        log(f"  {'에포크':<8} {'r(Base)':>10} {'r(S¹)':>10} {'Δr':>10} "
            f"{'std(B)':>10} {'std(S¹)':>10} {'κ(B)':>10} {'κ(S¹)':>10}")
        log(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")

        for ep_i in [0, 1, 5, 10, 20, 50, EPOCHS]:
            if ep_i < len(r_b):
                dr = r_s[ep_i] - r_b[ep_i]
                log(f"  {ep_i:<8} {r_b[ep_i]:>10.5f} {r_s[ep_i]:>10.5f} {dr:>+10.5f} "
                    f"{std_b[ep_i]:>10.4f} {std_s[ep_i]:>10.4f} "
                    f"{k_b[ep_i]:>10.0f} {k_s[ep_i]:>10.0f}")

        # r 수렴 속도: r > 0.999 도달 에포크
        ep_999_b = next((i for i, r in enumerate(r_b) if r > 0.999), -1)
        ep_999_s = next((i for i, r in enumerate(r_s) if r > 0.999), -1)
        log(f"\n    r>0.999 도달: Baseline ep{ep_999_b}, S¹ ep{ep_999_s}")
        if ep_999_b > 0 and ep_999_s > 0:
            log(f"    가속 비율: {ep_999_b/ep_999_s:.2f}x")

    # ────────────────────────────────────────────────────────────────
    # 앙상블 요약
    # ────────────────────────────────────────────────────────────────
    log(f"\n{'=' * 72}")
    log(f"  앙상블 요약 ({len(SEEDS)} seeds)")
    log(f"{'=' * 72}")

    # r 궤적 평균
    r_b_all = np.array([h[0]["r"] for h in all_baseline])  # [n_seeds, EPOCHS+1]
    r_s_all = np.array([h[0]["r"] for h in all_s1])
    std_b_all = np.array([h[0]["std"] for h in all_baseline])
    std_s_all = np.array([h[0]["std"] for h in all_s1])
    k_b_all = np.array([h[0]["kappa_std"] for h in all_baseline])
    k_s_all = np.array([h[0]["kappa_std"] for h in all_s1])

    log(f"\n  에포크별 평균 r:")
    log(f"  {'에포크':<8} {'r(Base)':>15} {'r(S¹)':>15} {'Δr':>12}")
    log(f"  {'-'*8} {'-'*15} {'-'*15} {'-'*12}")
    for ep_i in [0, 1, 5, 10, 20, 50, EPOCHS]:
        if ep_i < r_b_all.shape[1]:
            rb_m = r_b_all[:, ep_i].mean()
            rb_s = r_b_all[:, ep_i].std()
            rs_m = r_s_all[:, ep_i].mean()
            rs_s = r_s_all[:, ep_i].std()
            dr = rs_m - rb_m
            log(f"  {ep_i:<8} {rb_m:.5f}±{rb_s:.5f} {rs_m:.5f}±{rs_s:.5f} {dr:>+.5f}")

    # 최종 지표
    log(f"\n  최종 (ep {EPOCHS}):")
    rb_final = r_b_all[:, -1]
    rs_final = r_s_all[:, -1]
    stdb_final = std_b_all[:, -1]
    stds_final = std_s_all[:, -1]
    kb_final = k_b_all[:, -1]
    ks_final = k_s_all[:, -1]

    stdb_init = std_b_all[:, 0]
    stds_init = std_s_all[:, 0]
    ratio_b_arr = stdb_init / (stdb_final + 1e-12)
    ratio_s_arr = stds_init / (stds_final + 1e-12)

    log(f"    Baseline: r={rb_final.mean():.5f}±{rb_final.std():.5f}, "
        f"std={stdb_final.mean():.4f}±{stdb_final.std():.4f}, "
        f"κ={kb_final.mean():.0f}±{kb_final.std():.0f}")
    log(f"    S¹ Geo:   r={rs_final.mean():.5f}±{rs_final.std():.5f}, "
        f"std={stds_final.mean():.4f}±{stds_final.std():.4f}, "
        f"κ={ks_final.mean():.0f}±{ks_final.std():.0f}")
    log(f"    std ratio: Baseline {ratio_b_arr.mean():.2f}±{ratio_b_arr.std():.2f}x, "
        f"S¹ {ratio_s_arr.mean():.2f}±{ratio_s_arr.std():.2f}x")

    # F2 비교
    log(f"\n  F₂ 검출 비교:")
    ratio_b_f2 = [h[1] for h in all_baseline]
    det_b_f2 = [h[2] for h in all_baseline]
    ratio_s_f2 = [h[1] for h in all_s1]
    det_s_f2 = [h[2] for h in all_s1]
    tot = all_baseline[0][3]

    log(f"    Baseline: |F₂| ratio={np.mean(ratio_b_f2):.4f}±{np.std(ratio_b_f2):.4f}, "
        f"검출={np.mean(det_b_f2):.1f}/{tot}")
    log(f"    S¹ Geo:   |F₂| ratio={np.mean(ratio_s_f2):.4f}±{np.std(ratio_s_f2):.4f}, "
        f"검출={np.mean(det_s_f2):.1f}/{tot}")

    # von Mises κ 궤적 (수학자 질문 #3 선행 데이터)
    log(f"\n  von Mises κ 궤적 (에포크별 평균):")
    log(f"  {'에포크':<8} {'κ(Base)':>12} {'κ(S¹)':>12} {'κ비율':>10}")
    log(f"  {'-'*8} {'-'*12} {'-'*12} {'-'*10}")
    for ep_i in [0, 1, 5, 10, 20, 50, EPOCHS]:
        if ep_i < k_b_all.shape[1]:
            kb_m = k_b_all[:, ep_i].mean()
            ks_m = k_s_all[:, ep_i].mean()
            ratio_k = ks_m / (kb_m + 1e-12)
            log(f"  {ep_i:<8} {kb_m:>12.0f} {ks_m:>12.0f} {ratio_k:>10.2f}x")

    # ────────────────────────────────────────────────────────────────
    # 판정
    # ────────────────────────────────────────────────────────────────
    dr_final = rs_final.mean() - rb_final.mean()
    dstd_ratio = ratio_s_arr.mean() - ratio_b_arr.mean()
    dk_final = ks_final.mean() - kb_final.mean()
    det_diff = np.mean(det_s_f2) - np.mean(det_b_f2)

    log(f"\n{'=' * 72}")
    log(f"  판정 기준:")
    log(f"    Δr (final): {dr_final:+.5f}")
    log(f"    Δ(std ratio): {dstd_ratio:+.2f}x")
    log(f"    Δκ (final): {dk_final:+.0f}")
    log(f"    검출 차이: {det_diff:+.1f}개")

    # 수렴 속도 비교: r > 0.999 도달 에포크
    ep999_b = []
    ep999_s = []
    for i in range(len(SEEDS)):
        rb = r_b_all[i]
        rs = r_s_all[i]
        e_b = next((j for j, r in enumerate(rb) if r > 0.999), EPOCHS)
        e_s = next((j for j, r in enumerate(rs) if r > 0.999), EPOCHS)
        ep999_b.append(e_b)
        ep999_s.append(e_s)
    log(f"    r>0.999 도달: Baseline 평균 ep{np.mean(ep999_b):.1f}, "
        f"S¹ 평균 ep{np.mean(ep999_s):.1f}")

    # 판정 로직
    faster = np.mean(ep999_s) < np.mean(ep999_b) * 0.8  # 20% 이상 빠르면
    higher_k = dk_final > 50  # κ 유의미 증가
    better_det = det_diff > 3  # 검출 3개 이상 증가

    if faster or (higher_k and better_det):
        verdict = "양성: L_geo가 동기화를 가속 (수학자 예측 확인)"
        if faster:
            accel = np.mean(ep999_b) / max(np.mean(ep999_s), 1)
            verdict += f" — 수렴 {accel:.1f}x 가속"
    elif dk_final < -50 or det_diff < -3:
        verdict = "음성: L_geo가 동기화에 부정적 영향"
    else:
        verdict = "중립: L_geo와 L_tgt의 동기화 궤적에 유의미한 차이 없음"

    log(f"\n  판정: {verdict}")

    elapsed = time.time() - start
    log(f"  총 실행 시간: {elapsed:.0f}s")

    # 수학자를 위한 핵심 답변
    log(f"\n{'=' * 72}")
    log(f"  수학자 질문에 대한 답변:")
    log(f"{'=' * 72}")
    log(f"  Q1: L_geo가 r 궤적을 바꾸는가?")
    log(f"    → r(final) 차이: {dr_final:+.5f}, std ratio 차이: {dstd_ratio:+.2f}x")
    log(f"  Q2: 더 빠른 동기화?")
    log(f"    → r>0.999 도달: Baseline ep{np.mean(ep999_b):.1f} vs S¹ ep{np.mean(ep999_s):.1f}")
    log(f"  Q3: κ-검출 상관 (선행 데이터):")
    log(f"    → Baseline: κ={kb_final.mean():.0f}, 검출={np.mean(det_b_f2):.1f}")
    log(f"    → S¹ Geo:   κ={ks_final.mean():.0f}, 검출={np.mean(det_s_f2):.1f}")

    # 저장
    with open(results_path, "w") as f:
        f.write("\n".join(out) + "\n")
    log(f"\n  저장: {results_path}")


if __name__ == "__main__":
    main()
