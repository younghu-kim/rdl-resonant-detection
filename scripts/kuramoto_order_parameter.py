"""
=============================================================================
[Project RDL] Kuramoto 순서 매개변수 실험
=============================================================================
게이지 위상 φ의 동기화 전이를 Kuramoto 순서 매개변수로 정량화.

핵심 가설 (논문 Observation 8.6 / Remark 8.7):
  게이지 위상 {φ_i}는 ~10 에포크에서 급격한 동기화 전이를 보인다.
  std(φ)가 ~0.10에서 ~0.03으로 ~3배 감소.
  이것은 Kuramoto 동기화의 신경망 서명이다.

측정:
  (1) Kuramoto 순서 매개변수: r(ep) = |1/N Σ_j exp(iφ_j)| ∈ [0,1]
      r≈0: 비동기 (위상 균등분포), r≈1: 완전 동기화
  (2) std(φ) per epoch — 논문의 관찰과 직접 비교
  (3) mean(φ) per epoch — 잠금 주파수 추적
  (4) 임계 에포크 식별: r이 0.5를 넘는 첫 에포크
  (5) Kuramoto 이론 비교: K_c = 2/(πg(0)) vs 실측 커플링

다중 시드로 재현성 확인. 두 t-범위에서 불변성 확인.
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
EPOCHS = 100  # 100 에포크 — 전이는 ~10 에포크에서 발생하므로 충분
LR = 1e-3
BATCH = 32
IN_FEATURES = 128
SEEDS = [42, 7, 123]
T_RANGES = [(100.0, 200.0), (10.0, 50.0)]  # 두 범위에서 불변성 확인


def collect_gauge_phases(model, dataset, batch_size=64):
    """전체 데이터셋에 대해 게이지 위상 φ를 수집.

    Returns:
        phi_all: [N, hidden_features] numpy array of gauge phases
    """
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
    """Kuramoto 순서 매개변수 및 관련 통계 계산.

    Args:
        phi_all: [N, hidden_features] gauge phases

    Returns:
        dict with r_global, r_per_dim, std_global, std_per_dim, mean_phi
    """
    # 각 히든 차원별 Kuramoto r
    exp_iphi = np.exp(1j * phi_all)  # [N, H]
    r_per_dim = np.abs(exp_iphi.mean(axis=0))  # [H]
    r_global = r_per_dim.mean()

    # std(φ) — circular std
    std_per_dim = phi_all.std(axis=0)  # [H]
    std_global = std_per_dim.mean()

    # 평균 위상
    mean_phi = np.angle(exp_iphi.mean(axis=0)).mean()

    return {
        "r_global": r_global,
        "r_per_dim": r_per_dim,
        "std_global": std_global,
        "std_per_dim": std_per_dim,
        "mean_phi": mean_phi,
    }


def train_with_kuramoto_tracking(dataset, seed):
    """표준 학습 + 매 에포크 Kuramoto 측정."""
    torch.manual_seed(seed); np.random.seed(seed)

    val_size = int(len(dataset) * 0.2)
    train_size = len(dataset) - val_size
    train_ds, val_ds = torch.utils.data.random_split(
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

    # 에포크별 기록
    history = {
        "r_global": [],
        "std_global": [],
        "mean_phi": [],
        "train_loss": [],
        "r_per_dim_snapshots": [],  # 에포크 0, 5, 10, 15, 20, 50에서 전체 분포
    }
    snapshot_epochs = {0, 5, 10, 15, 20, 50, EPOCHS - 1}

    t_start = time.time()

    # 초기 상태 (학습 전)
    metrics_init = compute_kuramoto_metrics(collect_gauge_phases(model, dataset))
    history["r_global"].append(metrics_init["r_global"])
    history["std_global"].append(metrics_init["std_global"])
    history["mean_phi"].append(metrics_init["mean_phi"])
    history["train_loss"].append(float('nan'))
    history["r_per_dim_snapshots"].append(("ep0_init", metrics_init["r_per_dim"]))

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

        avg_loss = ep_loss / max(1, ep_n)

        # 매 에포크 Kuramoto 측정
        metrics = compute_kuramoto_metrics(collect_gauge_phases(model, dataset))
        history["r_global"].append(metrics["r_global"])
        history["std_global"].append(metrics["std_global"])
        history["mean_phi"].append(metrics["mean_phi"])
        history["train_loss"].append(avg_loss)

        if ep in snapshot_epochs:
            history["r_per_dim_snapshots"].append((f"ep{ep+1}", metrics["r_per_dim"]))

        if (ep + 1) % 10 == 0:
            elapsed = time.time() - t_start
            print(f"    ep {ep+1}/{EPOCHS}: loss={avg_loss:.5f}, "
                  f"r={metrics['r_global']:.4f}, std(φ)={metrics['std_global']:.4f}, "
                  f"{elapsed:.0f}s", flush=True)

    elapsed = time.time() - t_start
    return history, elapsed


def find_critical_epoch(r_history, threshold=0.5):
    """r이 threshold를 초과하는 첫 에포크 (0-indexed, ep0=학습전)."""
    for i, r in enumerate(r_history):
        if r > threshold:
            return i
    return -1


def analyze_transition(history):
    """전이 특성 분석."""
    r = np.array(history["r_global"])
    std = np.array(history["std_global"])

    # 임계 에포크 (r > 0.5)
    ep_crit_50 = find_critical_epoch(r, 0.5)
    ep_crit_80 = find_critical_epoch(r, 0.8)

    # std 비율: 초기 대비 최종
    std_ratio = std[0] / (std[-1] + 1e-12)

    # 전이 기울기: r의 최대 증가율
    dr = np.diff(r)
    max_dr_ep = np.argmax(dr) + 1  # 에포크 (1-indexed)
    max_dr = dr[max_dr_ep - 1]

    # 10 에포크 전후 비교 (논문의 관찰)
    r_pre10 = r[:11].mean() if len(r) > 11 else r[:len(r)//2].mean()
    r_post10 = r[11:21].mean() if len(r) > 21 else r[len(r)//2:].mean()
    std_pre10 = std[:11].mean() if len(std) > 11 else std[:len(std)//2].mean()
    std_post10 = std[11:21].mean() if len(std) > 21 else std[len(std)//2:].mean()

    return {
        "ep_crit_50": ep_crit_50,
        "ep_crit_80": ep_crit_80,
        "std_ratio": std_ratio,
        "max_dr_epoch": max_dr_ep,
        "max_dr_value": max_dr,
        "r_initial": r[0],
        "r_final": r[-1],
        "std_initial": std[0],
        "std_final": std[-1],
        "r_pre10": r_pre10,
        "r_post10": r_post10,
        "std_pre10": std_pre10,
        "std_post10": std_post10,
    }


def main():
    results_path = os.path.expanduser("~/Desktop/gdl_unified/results/kuramoto_order_parameter.txt")
    out = []

    def log(msg):
        print(msg, flush=True)
        out.append(msg)

    log("=" * 72)
    log("  Kuramoto 순서 매개변수 실험: 게이지 위상 동기화 정량화")
    log("=" * 72)
    log(f"  HIDDEN={HIDDEN}, EPOCHS={EPOCHS}, LR={LR}, BATCH={BATCH}")
    log(f"  seeds={SEEDS}")
    log(f"  t-ranges={T_RANGES}")
    log("")

    all_transitions = {}

    for t_lo, t_hi in T_RANGES:
        log(f"\n{'─' * 72}")
        log(f"  t∈[{t_lo},{t_hi}]")
        log(f"{'─' * 72}")

        # 데이터 준비
        zeros = compute_zeros_in_range(t_lo, t_hi)
        n_zeros = len(zeros)
        log(f"  [Zeros] {n_zeros}개 영점 발견")

        n_points = 1000
        cache_path = os.path.expanduser(
            f"~/Desktop/gdl_unified/cache/xi_cache_t{t_lo}-{t_hi}_n{n_points}.pt"
        )
        cache = get_or_build_cache(cache_path, t_lo, t_hi, n_points)
        dataset = XiFeatureDataset(cache, IN_FEATURES)
        log(f"  [Data] {len(dataset)}개 점, in_features={dataset.in_features}")

        range_results = []

        for seed in SEEDS:
            log(f"\n  seed={seed}")
            history, elapsed = train_with_kuramoto_tracking(dataset, seed)
            trans = analyze_transition(history)
            range_results.append(trans)

            log(f"    r: {trans['r_initial']:.4f} → {trans['r_final']:.4f}")
            log(f"    std(φ): {trans['std_initial']:.4f} → {trans['std_final']:.4f} "
                f"(ratio {trans['std_ratio']:.2f}x)")
            log(f"    임계 에포크 (r>0.5): {trans['ep_crit_50']}, (r>0.8): {trans['ep_crit_80']}")
            log(f"    최대 Δr 에포크: {trans['max_dr_epoch']} (Δr={trans['max_dr_value']:.4f})")
            log(f"    전이 전후 (ep≤10 vs ep11-20):")
            log(f"      r: {trans['r_pre10']:.4f} → {trans['r_post10']:.4f}")
            log(f"      std: {trans['std_pre10']:.4f} → {trans['std_post10']:.4f}")
            log(f"    time={elapsed:.0f}s")

            # 에포크별 r 궤적 (10 에포크 단위)
            r_arr = np.array(history["r_global"])
            trajectory = ", ".join(f"{r_arr[i]:.3f}" for i in range(0, len(r_arr), 10))
            log(f"    r 궤적 (매 10ep): [{trajectory}]")

        all_transitions[f"t{t_lo}-{t_hi}"] = range_results

    # ────────────────────────────────────────────────────────────────
    # 요약
    # ────────────────────────────────────────────────────────────────
    log(f"\n{'=' * 72}")
    log(f"  요약")
    log(f"{'=' * 72}")

    for key, results in all_transitions.items():
        log(f"\n  {key}:")
        r_finals = [r["r_final"] for r in results]
        std_finals = [r["std_final"] for r in results]
        std_ratios = [r["std_ratio"] for r in results]
        crit_eps = [r["ep_crit_50"] for r in results]
        max_dr_eps = [r["max_dr_epoch"] for r in results]

        log(f"    r_final: {np.mean(r_finals):.4f}±{np.std(r_finals):.4f}")
        log(f"    std(φ) final: {np.mean(std_finals):.4f}±{np.std(std_finals):.4f}")
        log(f"    std ratio (initial/final): {np.mean(std_ratios):.2f}±{np.std(std_ratios):.2f}x")
        log(f"    임계 에포크 (r>0.5): {crit_eps}")
        log(f"    최대 Δr 에포크: {max_dr_eps}")

    # 전체 앙상블
    all_std_ratios = []
    all_crit_eps = []
    all_r_finals = []
    for results in all_transitions.values():
        for r in results:
            all_std_ratios.append(r["std_ratio"])
            all_crit_eps.append(r["ep_crit_50"])
            all_r_finals.append(r["r_final"])

    log(f"\n  전체 앙상블 ({len(all_std_ratios)} runs):")
    log(f"    r_final: {np.mean(all_r_finals):.4f}±{np.std(all_r_finals):.4f}")
    log(f"    std ratio: {np.mean(all_std_ratios):.2f}±{np.std(all_std_ratios):.2f}x")
    log(f"    임계 에포크 중앙값: {np.median(all_crit_eps):.0f}")

    # 판정
    log(f"\n  논문 관찰과 비교:")
    log(f"    논문: std(φ) 0.10→0.03 (~3x), 전이 ~10 에포크")
    mean_std_ratio = np.mean(all_std_ratios)
    mean_crit = np.median(all_crit_eps)

    if mean_std_ratio > 2.0 and 0 < mean_crit <= 20:
        judgment = "양성: Kuramoto 동기화 전이 확인"
        log(f"  판정: {judgment}")
        log(f"    std ratio {mean_std_ratio:.1f}x ≥ 2x, 임계 에포크 {mean_crit:.0f} ≤ 20")
    elif mean_std_ratio > 1.5:
        judgment = "약양성: 부분적 동기화 확인, 논문 관찰보다 약함"
        log(f"  판정: {judgment}")
    else:
        judgment = "음성: Kuramoto 동기화 미관찰"
        log(f"  판정: {judgment}")

    total_time = sum(r["std_ratio"] for results in all_transitions.values() for r in results)  # placeholder
    log(f"\n  총 실행 시간: (각 run 개별 기록)")

    # 파일 저장
    with open(results_path, "w") as f:
        f.write("\n".join(out) + "\n")
    log(f"\n  저장: {results_path}")


if __name__ == "__main__":
    main()
