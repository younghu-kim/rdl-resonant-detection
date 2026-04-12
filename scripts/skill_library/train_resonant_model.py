"""
스킬: 표준 MasterResonantNetwork 훈련 루프
============================================
Baseline (L_tgt) 및 S¹ Geodesic (L_geo) 훈련을 위한 표준 루프.

사용법:
    from train_resonant_model import train_baseline, train_s1_integrated

    model, val_loss = train_baseline(dataset, seed, epochs=100, lr=1e-3)
    model, val_loss = train_s1_integrated(dataset, seed, epochs=100, lr=1e-3)

입력:
    dataset: XiFeatureDataset
    seed: int (랜덤 시드)

출력:
    (model, best_val_loss)

주의:
    - hidden_features=64 (hidden_dim 아님!)
    - loss_fn(**outputs) 패턴 (개별 인자 전달 금지)
    - python -u 플래그 필수 (stdout 버퍼링 방지)
"""

import time
import numpy as np
import torch
import torch.nn as nn

import sys, os
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
from gdl.rdl.constants import PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss


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


def _make_loaders(dataset, seed, batch_size=32, val_frac=0.2):
    val_size = int(len(dataset) * val_frac)
    train_size = len(dataset) - val_size
    train_ds, val_ds = torch.utils.data.random_split(
        dataset, [train_size, val_size],
        generator=torch.Generator().manual_seed(seed)
    )
    train_loader = torch.utils.data.DataLoader(train_ds, batch_size=batch_size, shuffle=True, drop_last=True)
    val_loader = torch.utils.data.DataLoader(val_ds, batch_size=batch_size, shuffle=False)
    return train_loader, val_loader


def train_baseline(dataset, seed, epochs=100, lr=1e-3, hidden=64, batch_size=32, verbose=True):
    """표준 TotalResonanceLoss Baseline 훈련."""
    torch.manual_seed(seed); np.random.seed(seed)
    train_loader, val_loader = _make_loaders(dataset, seed, batch_size)

    model = MasterResonantNetwork(
        in_features=dataset.in_features, hidden_features=hidden, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5,
        pqo_mode='cos2'
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    best_val = float('inf')
    t_start = time.time()
    for ep in range(epochs):
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

        if (ep + 1) % 25 == 0 or ep == epochs - 1:
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
            if verbose:
                elapsed = time.time() - t_start
                print(f"    ep {ep+1}/{epochs}: train={ep_loss/max(1,ep_n):.5f}, val={vl_m:.5f}, {elapsed:.0f}s", flush=True)

    return model, best_val


def train_s1_integrated(dataset, seed, epochs=100, lr=1e-3, hidden=64,
                         batch_size=32, lambda_geo=1.0, verbose=True):
    """S¹ Geodesic 통합 훈련: L_tgt → L_geo 대체."""
    torch.manual_seed(seed); np.random.seed(seed)
    train_loader, val_loader = _make_loaders(dataset, seed, batch_size)

    model = MasterResonantNetwork(
        in_features=dataset.in_features, hidden_features=hidden, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    res_loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=0.0, lambda_pqo=0.5,
        pqo_mode='cos2'
    )
    geo_loss_fn = GeodesicTargetLoss(reduction='mean')
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    best_val = float('inf')
    t_start = time.time()
    for ep in range(epochs):
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

        if (ep + 1) % 25 == 0 or ep == epochs - 1:
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
            if verbose:
                elapsed = time.time() - t_start
                print(f"    ep {ep+1}/{epochs}: train={ep_loss/max(1,ep_n):.5f}, val={vl_m:.5f}, {elapsed:.0f}s", flush=True)

    return model, best_val
