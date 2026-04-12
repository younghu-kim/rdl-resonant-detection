"""
스킬: F₂ 기반 영점 검출 평가
=============================
모든 RDL 실험에서 공통으로 사용하는 |F₂| 잔차 기반 영점 탐지 평가 함수.

사용법:
    from eval_f2_detection import evaluate_f2_detection, compute_f2_landscape

    # 검출 통계 (recall, precision, F1 등)
    stats = evaluate_f2_detection(model, dataset, batch_size=64)

    # |F₂| 풍경만 필요할 때
    f2_arr = compute_f2_landscape(model, dataset, batch_size=64)

입력:
    model: MasterResonantNetwork (학습 완료 상태)
    dataset: XiFeatureDataset (.is_near_zero 속성 필요)

출력 (evaluate_f2_detection):
    dict with keys: f2_zero, f2_nonzero, ratio, detected, total_zeros,
                    precision, recall, f1, threshold, fp_indices, tp_indices

주의:
    - Z_out.abs() 사용 금지! phi/psi/L_G 기반 잔차만 사용
    - eval 시 torch.enable_grad() 필수 (MasterResonantNetwork의 일부 연산에 필요)
    - X_in.requires_grad_(True) 필수
"""

import numpy as np
import torch

import sys, os
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
from gdl.rdl.constants import PrecisionManager


def compute_f2_landscape(model, dataset, batch_size=64):
    """전체 데이터에서 |F₂| 값 배열 반환.

    Returns:
        np.ndarray: shape (N,), 각 데이터 점의 |F₂| 잔차 절대값
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
    return np.abs(np.concatenate(f2_vals))


def evaluate_f2_detection(model, dataset, batch_size=64):
    """F₂ 기반 영점 검출의 전체 통계를 반환.

    Returns:
        dict: f2_zero, f2_nonzero, ratio, detected, total_zeros,
              precision, recall, f1, threshold, fp_indices, tp_indices
    """
    f2_arr = compute_f2_landscape(model, dataset, batch_size)
    is_zero = dataset.is_near_zero.numpy()

    f2_zero = f2_arr[is_zero].mean() if is_zero.any() else 0
    f2_nonzero = f2_arr[~is_zero].mean() if (~is_zero).any() else 1
    ratio = f2_zero / (f2_nonzero + 1e-12)

    threshold = np.median(f2_arr) * 0.1 if len(f2_arr) > 0 else 0.01
    predicted_zero = f2_arr < threshold

    tp_mask = predicted_zero & is_zero
    fp_mask = predicted_zero & (~is_zero)

    n_tp = int(tp_mask.sum())
    n_fp = int(fp_mask.sum())
    total_zeros = int(is_zero.sum())
    n_pred = int(predicted_zero.sum())

    precision = n_tp / max(n_pred, 1)
    recall = n_tp / max(total_zeros, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-12)

    return {
        'f2_zero': f2_zero,
        'f2_nonzero': f2_nonzero,
        'ratio': ratio,
        'detected': n_tp,
        'total_zeros': total_zeros,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'threshold': threshold,
        'fp_indices': np.where(fp_mask)[0],
        'tp_indices': np.where(tp_mask)[0],
        'f2_arr': f2_arr,
    }
