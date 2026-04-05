"""
=============================================================================
[Project RDL] Observability - Curvature Gate B 3D Landscape Animator
=============================================================================
곡률 손실(Gate B: v_σ=0, A'=0, A''<0)이 파라미터 공간의
진폭 지형을 어떻게 형성하는지 3D 지형도로 시각화한다.

섭동 그리드 위에서 |Z|를 높이로 삼아 극대점 돔 형성을 관찰한다.
"""

import os
import numpy as np
import torch
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.animation as animation
from typing import Callable, Dict, List

from gdl.rdl.constants import PrecisionManager

class CurvatureGateBAnimator:
    """데이터 주변 섭동 공간의 3D 진폭/곡률 지형을 캡처하여 애니메이션으로 렌더링하는 엔진"""

    def __init__(self, X_base: torch.Tensor, radius: float = 2.0, resolution: int = 40, save_dir: str = "outputs/animations"):
        """
        Args:
            X_base (torch.Tensor): 탐색의 중심이 될 단일 원본 데이터 (Shape: [Features])
            radius (float): X_base를 중심으로 스캔할 섭동 평면의 반경
            resolution (int): 2D 그리드 해상도 (예: 40x40 = 1600개 지점 스캔)
            save_dir (str): 애니메이션 저장 경로
        """
        self.save_dir = save_dir
        os.makedirs(self.save_dir, exist_ok=True)
        self.frames: List[Dict] = []

        # 사이버네틱 우주 테마 (다크 모드)
        plt.style.use('dark_background')

        # -------------------------------------------------------------------
        # [Step 1] 원본 데이터 평탄화 및 직교 섭동 벡터(v1, v2) 생성
        # -------------------------------------------------------------------
        if X_base.dim() > 1:
            X_base = X_base[0] # 첫 번째 샘플 강제 추출

        X_flat = X_base.detach().clone().flatten()
        dim = X_flat.size(0)
        device = X_base.device
        dtype = PrecisionManager.REAL_DTYPE

        # 일관된 관측 캔버스를 위해 난수 시드 고정
        gen = torch.Generator(device=device)
        gen.manual_seed(42)
        v1 = torch.randn(dim, generator=gen, dtype=dtype, device=device)
        v2 = torch.randn(dim, generator=gen, dtype=dtype, device=device)

        # 그람-슈미트(Gram-Schmidt) 직교화를 통해 v1과 v2를 완벽한 수직 기저로 변환
        v1 = v1 / torch.norm(v1)
        v2 = v2 - torch.dot(v1, v2) * v1
        v2 = v2 / torch.norm(v2)

        # -------------------------------------------------------------------
        # [Step 2] 2D 탐색 공간 메쉬그리드(Grid) 생성 및 배치 텐서화
        # -------------------------------------------------------------------
        u_vals = np.linspace(-radius, radius, resolution)
        v_vals = np.linspace(-radius, radius, resolution)
        self.U, self.V = np.meshgrid(u_vals, v_vals)

        u_flat = torch.tensor(self.U.flatten(), dtype=dtype, device=device).unsqueeze(1)
        v_flat = torch.tensor(self.V.flatten(), dtype=dtype, device=device).unsqueeze(1)

        # X_grid = X_base + u*v1 + v*v2 (중심점 기준 섭동 공간 텐서 생성, Shape: [Resolution^2, Features])
        X_grid = X_flat.unsqueeze(0) + (u_flat @ v1.unsqueeze(0)) + (v_flat @ v2.unsqueeze(0))

        # 모델의 원래 입력 형태로 복원
        if X_base.dim() == 1:
            self.X_grid = X_grid
        else:
            self.X_grid = X_grid.view(-1, *X_base.shape)

        self.resolution = resolution

    def capture(self, epoch: int, model_forward: Callable):
        """
        현재 에포크의 모델 상태를 바탕으로 2D 섭동 공간의 진폭(A)과 곡률(A'') 지형을 캡처합니다.
        (Autograd를 통해 2차 미분을 역추적한 후, VRAM 보호를 위해 Numpy 배열로 하역합니다.)
        """
        # [방어 기제 1] 원본 X_grid가 오염되지 않도록 복제 후 미분 추적기(requires_grad) 부착
        X_eval = self.X_grid.clone().detach().requires_grad_(True)

        # 1. 모델 순방향 전파
        Z_out = model_forward(X_eval)
        if isinstance(Z_out, dict):
            Z_out = Z_out['Z_out'] # MasterNetwork 호환성

        # 2. 안전한 진폭 에너지 도출
        A = torch.sqrt(Z_out.real**2 + Z_out.imag**2 + 1e-8)

        # -------------------------------------------------------------------
        # [극강 최적화] O(1) 초거대 병렬 Hessian 프록시 추출
        # -------------------------------------------------------------------
        # 1600개의 점을 일일이 미분하지 않고, 전체를 합산하여 단 한 번의 Backward로 병렬 추출
        A_sum = A.sum()

        # 1차 미분 (극점 추적용 기울기)
        A_prime = torch.autograd.grad(
            outputs=A_sum,
            inputs=X_eval,
            create_graph=True,
            retain_graph=True
        )[0]

        A_prime_sum = A_prime.sum()

        # 2차 미분 곡률 (A'', Hessian)
        A_dp = torch.autograd.grad(
            outputs=A_prime_sum,
            inputs=X_eval,
            create_graph=False # 캡처 모듈이므로 여기서 미분 트리를 끊어 VRAM 누수 방지
        )[0]

        # 방향 독립적 곡률 스칼라화 (Laplacian Proxy)
        if A_dp.dim() > 1:
            curvature = A_dp.sum(dim=list(range(1, A_dp.dim())))
        else:
            curvature = A_dp

        # [방어 기제 2] Autograd 차단 및 Numpy 변환
        A_np = A.detach().cpu().numpy().reshape(self.resolution, self.resolution)
        Curv_np = curvature.detach().cpu().numpy().reshape(self.resolution, self.resolution)

        self.frames.append({
            'epoch': epoch,
            'A': A_np,
            'Curv': Curv_np
        })

    def render(self, filename: str = "curvature_collapse.gif", fps: int = 15):
        """캡처된 지형 프레임들을 360도 회전하는 3D 애니메이션으로 렌더링합니다."""
        if not self.frames:
            print("[WARN] 캡처된 프레임이 없습니다.")
            return

        print(f"\n🌋 총 {len(self.frames)} 프레임 3D 곡률 붕괴 애니메이션 렌더링 시작...")

        fig = plt.figure(figsize=(12, 8), dpi=120)
        ax = fig.add_subplot(111, projection='3d')

        # 축 배경 투명화 및 디자인
        ax.set_facecolor('black')
        fig.patch.set_facecolor('black')
        ax.xaxis.set_pane_color((0.1, 0.1, 0.1, 1.0))
        ax.yaxis.set_pane_color((0.1, 0.1, 0.1, 1.0))
        ax.zaxis.set_pane_color((0.1, 0.1, 0.1, 1.0))
        ax.grid(color='grey', linestyle=':', linewidth=0.5, alpha=0.3)

        # 카메라 흔들림을 막기 위한 전역 Z축 및 곡률 Min/Max 고정
        all_Curv = np.array([f['Curv'] for f in self.frames])
        all_A = np.array([f['A'] for f in self.frames])

        curv_max = np.max(np.abs(all_Curv)) + 1e-5
        norm = Normalize(vmin=-curv_max, vmax=curv_max)
        cmap = cm.coolwarm # 음수(오목/안정)는 파란색, 양수(볼록/불안정)는 빨간색

        z_min, z_max = all_A.min(), all_A.max()
        z_margin = max((z_max - z_min) * 0.1, 0.1)

        ax.set_zlim(z_min - z_margin, z_max + z_margin)
        ax.set_xlabel(r"Orthogonal Vector $v_1$", color='cyan', labelpad=15)
        ax.set_ylabel(r"Orthogonal Vector $v_2$", color='cyan', labelpad=15)
        ax.set_zlabel(r"Amplitude Energy ($A=|Z|$)", color='white', labelpad=15)

        title = ax.set_title("", fontsize=16, fontweight='bold', color='white', pad=20)

        # 컬러바 부착
        m = cm.ScalarMappable(cmap=cmap, norm=norm)
        m.set_array([])
        cbar = fig.colorbar(m, ax=ax, shrink=0.5, aspect=15, pad=0.1)
        cbar.set_label("Curvature A'' (Blue: Stable/Concave, Red: Unstable/Convex)", color='white')
        cbar.ax.yaxis.set_tick_params(color='white')
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

        surf = None

        # 프레임 업데이트 함수
        def update(frame_idx):
            nonlocal surf

            # 이전 프레임의 3D 표면 메모리 삭제 (RAM 폭발 방지)
            if surf is not None:
                surf.remove()

            data = self.frames[frame_idx]
            colors = cmap(norm(data['Curv']))

            # [표면 플롯]
            surf = ax.plot_surface(self.U, self.V, data['A'], facecolors=colors,
                                   rstride=1, cstride=1, linewidth=0.2, antialiased=True, alpha=0.9)

            # [카메라 360도 궤도 비행 (Orbital View)]
            # 에포크가 지남에 따라 카메라가 랜드스케이프 주변을 한 바퀴 회전합니다.
            ax.view_init(elev=30, azim=45 + (frame_idx * 360 / len(self.frames)))

            title.set_text(f"[RDL Curvature Gate B] Epoch: {data['epoch']:03d}\n"
                           r"Objective: $A' \to 0, A'' \le 0$ (Stable Blue Crater)")

            return fig,

        anim = animation.FuncAnimation(fig, update, frames=len(self.frames), interval=1000 // fps, blit=False)

        save_path = os.path.join(self.save_dir, filename)
        try:
            anim.save(save_path, writer='pillow', fps=fps)
            print(f"[OK] 3D 곡률 지형 애니메이션 렌더링 완료: {os.path.abspath(save_path)}")
        except Exception as e:
            print(f"[FAIL] 렌더링 실패: {e}")
        finally:
            plt.close(fig)


# =================================================================
# 직접 실행 시 극점 지형 렌더링 가상 애니메이션 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()
    print("\n--- [RDL Observability] Curvature Gate B 3D Animator Test ---")

    # 1. 16차원 가상 데이터 기준점 세팅
    X_base = torch.zeros(16, dtype=PrecisionManager.REAL_DTYPE)
    animator = CurvatureGateBAnimator(X_base=X_base, save_dir="outputs/tests", resolution=40)

    # 2. 에포크에 따라 변화하는 가상의 텐서 물리 엔진(Model) 추론 함수
    total_epochs = 30
    print(f"▶ 가상의 {total_epochs} Epoch 3D 진폭/곡률 지형 붕괴 과정 캡처 중...")

    for epoch in range(1, total_epochs + 1):
        progress = epoch / total_epochs

        def mock_model_forward(X_eval):
            # 중심에서 거리 계산 (16차원 벡터의 노름)
            dist_sq = (X_eval ** 2).sum(dim=-1)

            # [베이스 돔] 에포크가 지날수록 진폭 에너지가 중앙으로 솟아오름
            height = 1.0 + progress * 9.0
            width = 3.0 - 2.0 * progress # 점차 가파른 돔으로 변함
            dome = height * torch.exp(-dist_sq / width)

            # [노이즈] 초반엔 울퉁불퉁하지만 진행될수록 평탄해짐
            noise = torch.randn_like(dome) * (1.0 - progress) * 0.5

            A_val = dome + noise
            # 이 함수를 통과하는 순간 PyTorch Autograd가 자동으로 진짜 1, 2차 미분(A, A'') 곡률을 추출합니다!
            return torch.complex(A_val, torch.zeros_like(A_val))

        # 캡처 모듈이 mock_forward를 실행하며 1,600개 그리드의 2차 미분(Hessian)을 스스로 산출합니다.
        animator.capture(epoch, mock_model_forward)

    # 3. 전체 궤적 GIF 렌더링!
    animator.render(filename="mock_curvature_collapse.gif", fps=15)
