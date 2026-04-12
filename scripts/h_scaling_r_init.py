"""
=============================================================================
[Project RDL] H-스케일링 실험: r_init의 기원
=============================================================================
수학자 질문: r_init ≈ 0.995는 어디서 오는가?

가설 A: Xavier+ReLU 구조적 → r_init ~ 1 - c/H  (H가 클수록 r_init → 1)
가설 B: 유한 크기 우연 → r_init 불규칙

실험 설계:
  H ∈ {16, 32, 64, 128, 256}에서 각 5 시드
  - r_init (학습 전) 측정
  - 30 에포크 학습 후 r_final, std(φ) ratio 측정
  - log(1 - r_init) vs log(H) 기울기로 스케일링 법칙 추정

t∈[100,200] 단일 범위 (t-범위 불변성 이미 확인됨).
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

# ── 설정 ──
H_VALUES = [16, 32, 64, 128, 256]
EPOCHS = 30       # r_init이 핵심이지만, 정련 속도도 비교
LR = 1e-3
BATCH = 32
IN_FEATURES = 128
SEEDS = [42, 7, 123, 314, 2024]
T_LO, T_HI = 100.0, 200.0


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


def compute_kuramoto_r(phi_all):
    """Kuramoto 순서 매개변수 r과 std(φ) 계산."""
    exp_iphi = np.exp(1j * phi_all)  # [N, H]
    r_per_dim = np.abs(exp_iphi.mean(axis=0))  # [H]
    r_global = r_per_dim.mean()
    std_global = phi_all.std(axis=0).mean()
    return r_global, std_global, r_per_dim


def run_single(dataset, hidden, seed):
    """단일 H, seed 조합에서 r_init 및 정련 측정."""
    torch.manual_seed(seed)
    np.random.seed(seed)

    val_size = int(len(dataset) * 0.2)
    train_size = len(dataset) - val_size
    train_ds, val_ds = torch.utils.data.random_split(
        dataset, [train_size, val_size],
        generator=torch.Generator().manual_seed(seed)
    )
    # H=256에서 메모리 절약: 배치 크기 조정
    bs = min(BATCH, 16) if hidden >= 256 else BATCH
    train_loader = torch.utils.data.DataLoader(
        train_ds, batch_size=bs, shuffle=True, drop_last=True
    )

    model = MasterResonantNetwork(
        in_features=dataset.in_features, hidden_features=hidden, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5,
        pqo_mode='cos2'
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)

    # ── r_init 측정 (학습 전) ──
    phi_init = collect_gauge_phases(model, dataset)
    r_init, std_init, r_per_dim_init = compute_kuramoto_r(phi_init)

    # ── 학습 + 매 에포크 r 추적 ──
    r_trajectory = [r_init]
    std_trajectory = [std_init]
    t0 = time.time()

    for ep in range(EPOCHS):
        model.train()
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

        # 매 에포크 r 측정
        phi_ep = collect_gauge_phases(model, dataset)
        r_ep, std_ep, _ = compute_kuramoto_r(phi_ep)
        r_trajectory.append(r_ep)
        std_trajectory.append(std_ep)

    elapsed = time.time() - t0
    r_final = r_trajectory[-1]
    std_final = std_trajectory[-1]
    std_ratio = std_init / (std_final + 1e-12)

    return {
        "r_init": r_init,
        "std_init": std_init,
        "r_final": r_final,
        "std_final": std_final,
        "std_ratio": std_ratio,
        "r_per_dim_init": r_per_dim_init,
        "r_trajectory": r_trajectory,
        "std_trajectory": std_trajectory,
        "elapsed": elapsed,
    }


def main():
    results_path = os.path.expanduser(
        "~/Desktop/gdl_unified/results/h_scaling_r_init.txt"
    )
    out = []

    def log(msg):
        print(msg, flush=True)
        out.append(msg)

    log("=" * 72)
    log("  H-스케일링 실험: r_init의 기원")
    log("=" * 72)
    log(f"  H_VALUES={H_VALUES}")
    log(f"  EPOCHS={EPOCHS}, LR={LR}, BATCH={BATCH}")
    log(f"  seeds={SEEDS}")
    log(f"  t∈[{T_LO},{T_HI}]")
    log("")

    # 데이터 준비 (한 번만)
    zeros = compute_zeros_in_range(T_LO, T_HI)
    log(f"  [Zeros] {len(zeros)}개 영점")
    cache_path = os.path.expanduser(
        f"~/Desktop/gdl_unified/cache/xi_cache_t{T_LO}-{T_HI}_n1000.pt"
    )
    cache = get_or_build_cache(cache_path, T_LO, T_HI, 1000)
    dataset = XiFeatureDataset(cache, IN_FEATURES)
    log(f"  [Data] {len(dataset)}개 점, in_features={dataset.in_features}")
    log("")

    # ── 메인 루프 ──
    all_results = {}  # H -> list of dicts

    for H in H_VALUES:
        log(f"{'─' * 72}")
        log(f"  H={H}")
        log(f"{'─' * 72}")

        h_results = []
        for seed in SEEDS:
            log(f"  seed={seed} ...", )
            try:
                res = run_single(dataset, H, seed)
                h_results.append(res)
                log(f"    r_init={res['r_init']:.6f}, r_final={res['r_final']:.6f}, "
                    f"std_ratio={res['std_ratio']:.2f}x, {res['elapsed']:.0f}s")
            except Exception as e:
                log(f"    ERROR: {e}")

        all_results[H] = h_results

        # H별 요약
        if h_results:
            r_inits = [r["r_init"] for r in h_results]
            r_finals = [r["r_final"] for r in h_results]
            std_ratios = [r["std_ratio"] for r in h_results]
            one_minus_r = [1.0 - r for r in r_inits]
            log(f"\n  H={H} 요약 ({len(h_results)} 시드):")
            log(f"    r_init:     {np.mean(r_inits):.6f} ± {np.std(r_inits):.6f}")
            log(f"    1-r_init:   {np.mean(one_minus_r):.6f} ± {np.std(one_minus_r):.6f}")
            log(f"    r_final:    {np.mean(r_finals):.6f} ± {np.std(r_finals):.6f}")
            log(f"    std_ratio:  {np.mean(std_ratios):.2f} ± {np.std(std_ratios):.2f}x")
        log("")

    # ── 스케일링 분석 ──
    log("=" * 72)
    log("  스케일링 분석: log(1-r_init) vs log(H)")
    log("=" * 72)

    H_arr = []
    mean_1mr = []
    std_1mr = []
    mean_r_init_arr = []

    for H in H_VALUES:
        if all_results[H]:
            r_inits = [r["r_init"] for r in all_results[H]]
            omr = [1.0 - r for r in r_inits]
            H_arr.append(H)
            mean_1mr.append(np.mean(omr))
            std_1mr.append(np.std(omr))
            mean_r_init_arr.append(np.mean(r_inits))

    H_arr = np.array(H_arr)
    mean_1mr = np.array(mean_1mr)

    log(f"\n  {'H':>6}  {'r_init':>10}  {'1-r_init':>12}  {'log(H)':>8}  {'log(1-r)':>10}")
    log(f"  {'─'*6}  {'─'*10}  {'─'*12}  {'─'*8}  {'─'*10}")
    for i, H in enumerate(H_arr):
        log(f"  {H:6d}  {mean_r_init_arr[i]:10.6f}  {mean_1mr[i]:12.6f}  "
            f"{np.log(H):8.3f}  {np.log(mean_1mr[i]):10.3f}")

    # 선형 회귀: log(1-r_init) = a * log(H) + b → 1-r_init ~ H^a
    if len(H_arr) >= 3:
        log_H = np.log(H_arr)
        log_1mr = np.log(mean_1mr)
        # 최소자승 피팅
        A = np.vstack([log_H, np.ones(len(log_H))]).T
        slope, intercept = np.linalg.lstsq(A, log_1mr, rcond=None)[0]
        c_fit = np.exp(intercept)
        # 잔차
        residuals = log_1mr - (slope * log_H + intercept)
        rmse = np.sqrt(np.mean(residuals**2))
        r_squared = 1 - np.sum(residuals**2) / np.sum((log_1mr - log_1mr.mean())**2)

        log(f"\n  파워 법칙 피팅: 1-r_init ~ c * H^α")
        log(f"    α (기울기) = {slope:.3f}")
        log(f"    c (절편)   = {c_fit:.4f}")
        log(f"    R²         = {r_squared:.4f}")
        log(f"    RMSE       = {rmse:.4f}")
        log(f"\n  해석:")
        if slope < -0.5 and r_squared > 0.8:
            log(f"    가설 A 지지: 1-r_init ~ H^({slope:.2f})")
            if abs(slope - (-1.0)) < 0.3:
                log(f"    기울기 ≈ -1 → 1-r_init ~ 1/H (수학자 예측 c/H와 일치)")
            elif abs(slope - (-2.0)) < 0.3:
                log(f"    기울기 ≈ -2 → 1-r_init ~ 1/H² (수학자 추측 O(1/H²)와 일치)")
            else:
                log(f"    기울기 = {slope:.2f} (정수 지수와 근접하지 않음)")
        elif r_squared < 0.5:
            log(f"    가설 B 지지: 스케일링 없음 (R²={r_squared:.2f} < 0.5)")
        else:
            log(f"    불확실: 약한 스케일링 (α={slope:.2f}, R²={r_squared:.2f})")

    # ── 정련 속도 비교 ──
    log(f"\n{'=' * 72}")
    log(f"  정련 속도: H별 (r_final - r_init) / epochs")
    log(f"{'=' * 72}")
    for H in H_VALUES:
        if all_results[H]:
            deltas = [r["r_final"] - r["r_init"] for r in all_results[H]]
            rates = [d / EPOCHS for d in deltas]
            log(f"  H={H:3d}: Δr={np.mean(deltas):.6f}±{np.std(deltas):.6f}, "
                f"rate={np.mean(rates)*1000:.4f}±{np.std(rates)*1000:.4f} ×10⁻³/ep")

    # ── r_per_dim 분포 분석 ──
    log(f"\n{'=' * 72}")
    log(f"  r_per_dim 초기 분포 (히든 차원별 Kuramoto r)")
    log(f"{'=' * 72}")
    for H in H_VALUES:
        if all_results[H]:
            # 첫 시드의 r_per_dim 분포 통계
            rpd = all_results[H][0]["r_per_dim_init"]
            log(f"  H={H:3d}: min={rpd.min():.4f}, median={np.median(rpd):.4f}, "
                f"max={rpd.max():.4f}, mean={rpd.mean():.4f}, std={rpd.std():.4f}")

    # ── 판정 ──
    log(f"\n{'=' * 72}")
    log(f"  판정")
    log(f"{'=' * 72}")

    if len(H_arr) >= 3:
        if slope < -0.5 and r_squared > 0.8:
            log(f"  양성: r_init의 H-의존성 확인")
            log(f"    1-r_init ~ H^({slope:.2f}), R²={r_squared:.3f}")
            log(f"    초기 동기화는 Xavier+ReLU 구조에서 유래 (우연이 아님)")
            log(f"    수학자의 '자발적 게이지 고정' 명제 실험적 지지")
        elif r_squared < 0.5:
            log(f"  음성: r_init의 H-의존성 미확인 (R²={r_squared:.2f})")
            log(f"    가설 B (유한 크기 우연) 지지")
        else:
            log(f"  중립: 약한 스케일링 관찰, 추가 H 값 필요")
            log(f"    α={slope:.2f}, R²={r_squared:.2f}")

    # 저장
    with open(results_path, "w") as f:
        f.write("\n".join(out) + "\n")
    log(f"\n  저장: {results_path}")


if __name__ == "__main__":
    main()
