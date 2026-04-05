"""
=============================================================================
[Project RDL] Observability - Phase-Target Alignment Animator
=============================================================================
위상 정렬(Phase Alignment) 수렴 과정을 시각화.

에포크별 예측 위상각 arg(Z)과 타겟 Ψ(n)의 정합 정도를 렌더링하고
.gif 애니메이션으로 저장한다.
"""

import os
import numpy as np
import torch
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from typing import List, Dict


class PhaseAnimator:
    """훈련 과정 중 위상의 궤적을 캡처하고 애니메이션으로 렌더링하는 시각화 엔진"""

    def __init__(self, save_dir: str = "outputs/animations"):
        """
        Args:
            save_dir (str): 시각화 이미지가 저장될 디렉토리 경로
        """
        self.save_dir = save_dir
        os.makedirs(self.save_dir, exist_ok=True)
        self.frames: List[Dict] = []

        # 일관된 다크 테마(사이버네틱 우주 느낌) 적용
        plt.style.use('dark_background')

    def capture(self, epoch: int, Z_out: torch.Tensor, Psi_target: torch.Tensor, batch_idx: int = 0):
        """
        한 에포크가 끝날 때마다 위상 텐서를 numpy로 변환하여 저장합니다.

        Args:
            epoch (int): 현재 훈련 에포크
            Z_out (torch.Tensor): 네트워크 최종 복소 출력 (Shape: [Batch, Seq_len])
            Psi_target (torch.Tensor): 목표 위상 (Shape: [Batch, Seq_len])
            batch_idx (int): 시각화할 배치 내 샘플 인덱스 (통상 0번 샘플 추적)
        """
        # [방어 기제 1] Autograd 그래프 단절 및 순수 Numpy 배열화 (VRAM 누수 원천 차단)
        # 텐서의 미분 사슬을 끊어내고 GPU에서 CPU로 내립니다.
        Z_np = Z_out[batch_idx].detach().cpu().numpy()
        Psi_np = Psi_target[batch_idx].detach().cpu().numpy()

        # [방어 기제 2] 위상 래핑 해제 (Phase Unwrapping)
        # arg(Z) 연산은 태생적으로 -π ~ π 사이를 진동하며 톱니파 형태로 끊어집니다.
        # np.unwrap() 알고리즘을 적용하여 불연속점들을 매끄러운 1D 누적 단조 곡선으로 복원합니다.
        pred_phase = np.unwrap(np.angle(Z_np))
        target_phase = np.unwrap(Psi_np)

        self.frames.append({
            'epoch': epoch,
            'pred': pred_phase,
            'target': target_phase
        })

    def render(self, filename: str = "phase_alignment.gif", fps: int = 10):
        """
        저장된 프레임들을 모아 .gif 애니메이션으로 렌더링합니다.
        (파이썬 기본 내장 Pillow 라이터를 사용하여 별도 설치 의존성 제거)
        """
        if not self.frames:
            print("[WARN] 캡처된 프레임이 없습니다. capture()를 먼저 호출하세요.")
            return

        print(f"\n🎬 총 {len(self.frames)} 프레임 애니메이션 렌더링 시작...")

        # 1. 캔버스 세팅
        fig, ax = plt.subplots(figsize=(10, 5), dpi=120)
        seq_len = len(self.frames[0]['target'])
        x = np.arange(seq_len)

        # 2. 동적 Y축 스케일링을 위한 전역 Min/Max 탐색
        all_preds = np.array([f['pred'] for f in self.frames])
        all_targets = np.array([f['target'] for f in self.frames])

        y_min = min(all_preds.min(), all_targets.min())
        y_max = max(all_preds.max(), all_targets.max())
        margin = max(abs(y_max - y_min) * 0.1, 0.5)

        ax.set_xlim(0, seq_len - 1)
        ax.set_ylim(y_min - margin, y_max + margin)

        # 3. 그래프 플롯 선언 (목표는 시안색 파선, 예측은 마젠타 실선)
        line_target, = ax.plot([], [], color='cyan', linestyle='--', lw=2, alpha=0.8, label=r'Target Phase ($\Psi$)')
        line_pred, = ax.plot([], [], color='magenta', linestyle='-', lw=2.5, alpha=0.9, label=r'Predicted Phase ($\Phi$)')
        fill_col = None

        title = ax.set_title("", fontsize=14, fontweight='bold', color='white')
        ax.set_xlabel("Time Step (Sequence Dimension)", fontsize=12, color='lightgray')
        ax.set_ylabel("Accumulated Phase Angle (Unwrapped Radian)", fontsize=12, color='lightgray')
        ax.legend(loc='upper right', facecolor='black', edgecolor='white')
        ax.grid(True, linestyle=':', alpha=0.3, color='gray')

        # 4. 애니메이션 프레임 업데이트 함수
        def update(frame_idx):
            nonlocal fill_col
            data = self.frames[frame_idx]

            line_pred.set_data(x, data['pred'])
            line_target.set_data(x, data['target'])

            # 예측선과 타겟선 사이의 잔차 오차 영역을 하얀색 반투명으로 칠함
            # 훈련이 성공적이라면 에포크가 지날수록 이 하얀 안개가 걷혀야 함
            if fill_col is not None:
                fill_col.remove()
            fill_col = ax.fill_between(x, data['target'], data['pred'], color='white', alpha=0.15)

            # 오차(MAE) 계산을 통해 시각적 수렴도 실시간 타이틀 로깅
            mae = np.mean(np.abs(data['pred'] - data['target']))
            title.set_text(f"[RDL Phase Alignment] Epoch: {data['epoch']:03d} | MAE: {mae:.4f}")

            return line_target, line_pred, title, fill_col

        # 5. 애니메이션 객체 생성 (blit=False로 호환성 유지)
        anim = animation.FuncAnimation(
            fig, update, frames=len(self.frames),
            interval=1000 // fps, blit=False
        )

        # 6. 디스크 저장
        save_path = os.path.join(self.save_dir, filename)
        try:
            anim.save(save_path, writer='pillow', fps=fps)
            print(f"[OK] 애니메이션 렌더링 완료. 저장 위치: {os.path.abspath(save_path)}")
        except Exception as e:
            print(f"[FAIL] 렌더링 실패: {e}")
        finally:
            plt.close(fig) # VRAM 및 시스템 RAM 누수 방지를 위한 캔버스 소멸


# =================================================================
# 직접 실행 시 애니메이터 프레임 렌더링 가상 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    from gdl.rdl.constants import PrecisionManager
    PrecisionManager.setup_precision()

    print("\n--- [RDL Observability] Phase Animator Rendering Test ---")

    seq_len = 64
    animator = PhaseAnimator(save_dir="outputs/tests")

    # 1. 불변하는 가상의 타겟 위상 곡선 (Ground Truth)
    x = np.linspace(0, 4 * np.pi, seq_len)
    target_np = np.sin(x) + 0.5 * np.cos(2 * x)
    Psi_target = torch.tensor(target_np, dtype=PrecisionManager.REAL_DTYPE).unsqueeze(0)

    # 2. 에포크 진화 시뮬레이션 (점차 노이즈가 감소하며 타겟으로 수렴)
    total_epochs = 30
    print(f"▶ 가상의 {total_epochs} Epoch 훈련 궤적 캡처 중...")

    for epoch in range(1, total_epochs + 1):
        # 초반에는 노이즈 스케일과 위상 편이가 크지만 지수적으로 감쇠
        decay = np.exp(-epoch / 5.0)
        noise = np.random.randn(seq_len) * decay * 1.5

        current_angle = target_np + noise + (decay * np.pi)

        # 복소 공간 Z로 포장하여 텐서 생성
        Z_out = torch.complex(
            torch.tensor(np.cos(current_angle), dtype=PrecisionManager.REAL_DTYPE),
            torch.tensor(np.sin(current_angle), dtype=PrecisionManager.REAL_DTYPE)
        ).unsqueeze(0)

        # 모델 출력을 훈련 루프 중간에 캡처하는 상황 모사
        animator.capture(epoch, Z_out, Psi_target)

    # 3. 전체 궤적을 엮어 GIF 렌더링 가동!
    test_filename = "mock_phase_alignment.gif"
    animator.render(filename=test_filename, fps=10)
