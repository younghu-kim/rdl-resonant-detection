"""
=============================================================================
[Project RDL] Observability - Resonance Convergence Terminal Dashboard
=============================================================================
훈련 중 4항 공명 조건(판별식, 곡률, 타겟, PQO) 수렴 상태를
실시간 터미널 UI로 렌더링하는 대시보드.

모든 지표가 허용 오차(epsilon) 이내로 수렴하면
조기 종료(Early Stopping) 시그널을 발생시킨다.
"""

import time
import math
from typing import Dict, Optional

# 터미널 색상 제어를 위한 ANSI Escape Codes
class Colors:
    CYAN = '\033[96m'
    MAGENTA = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    WHITE = '\033[97m'
    RESET = '\033[0m'
    BOLD = '\033[1m'

class ResonanceDashboard:
    """3중 공명 붕괴 조건을 감시하고 터미널 UI를 렌더링하는 마스터 관측소"""

    def __init__(self, eps_res: float = 1e-2, eps_curv: float = 1e-2, eps_tgt: float = 1e-2):
        """
        Args:
            eps_res (float): 위상 정렬 허용 오차
            eps_curv (float): 곡률 안정화 허용 오차
            eps_tgt (float): 목표 매칭 허용 오차
        """
        self.eps = {
            'L_res': eps_res,
            'L_curv': eps_curv,
            'L_tgt': eps_tgt
        }
        self.is_collapsed = False
        self.start_time = time.time()

        # 시작 헤더 출력
        print(f"\n{Colors.BOLD}{Colors.CYAN}={'='*80}{Colors.RESET}")
        print(f"{Colors.BOLD}{Colors.MAGENTA}  [Project RDL] 🌌 공명 특이점 터미널 관측소 (Resonance Observatory){Colors.RESET}")
        print(f"{Colors.BOLD}{Colors.CYAN}={'='*80}{Colors.RESET}")
        print(f"{Colors.WHITE}▶ 특이점 붕괴 허용 오차 임계선 (ε){Colors.RESET}")
        print(f"  - Phase Align (L_res)  < {eps_res}")
        print(f"  - Curvature   (L_curv) < {eps_curv}")
        print(f"  - Target Match(L_tgt)  < {eps_tgt}\n")

    def _get_status_color(self, val: float, target_eps: float) -> str:
        """오차 수렴 정도에 따라 텍스트 색상을 동적으로 반환합니다."""
        if val <= target_eps:
            return Colors.GREEN  # 완벽한 공명 (특이점 돌파)
        elif val <= target_eps * 5:
            return Colors.YELLOW # 근접 중 (경고)
        else:
            return Colors.RED    # 붕괴 전 (혼돈)

    def log_epoch(self, epoch: int, metrics: Dict[str, float], val_loss: Optional[float] = None) -> bool:
        """
        매 에포크마다 훈련 상태를 사이버네틱 UI로 출력하고 붕괴 여부를 검사합니다.

        Returns:
            bool: 3중 공명 조건이 달성되어 특이점을 돌파했으면 True, 아니면 False
        """
        l_res = metrics.get('L_res', float('inf'))
        l_curv = metrics.get('L_curv', float('inf'))
        l_tgt = metrics.get('L_tgt', float('inf'))
        t_loss = metrics.get('Total', float('inf'))

        c_res = self._get_status_color(l_res, self.eps['L_res'])
        c_curv = self._get_status_color(l_curv, self.eps['L_curv'])
        c_tgt = self._get_status_color(l_tgt, self.eps['L_tgt'])

        elapsed = time.time() - self.start_time
        val_str = f" | {Colors.MAGENTA}Val: {val_loss:.4f}{Colors.RESET}" if val_loss is not None else ""

        # 상태 표시자 마킹 (해당 지표가 특이점 범위를 뚫었으면 초록색 별, 아니면 빈 공간)
        m_res = f"{Colors.GREEN}★{Colors.RESET}" if l_res <= self.eps['L_res'] else " "
        m_curv = f"{Colors.GREEN}★{Colors.RESET}" if l_curv <= self.eps['L_curv'] else " "
        m_tgt = f"{Colors.GREEN}★{Colors.RESET}" if l_tgt <= self.eps['L_tgt'] else " "

        log_msg = (
            f"{Colors.BOLD}[Epoch {epoch:04d}]{Colors.RESET} ⏱ {elapsed:05.1f}s | "
            f"Tot: {Colors.CYAN}{t_loss:8.4f}{Colors.RESET} || "
            f"Res: {m_res}{c_res}{l_res:8.4f}{Colors.RESET} | "
            f"Curv: {m_curv}{c_curv}{l_curv:8.4f}{Colors.RESET} | "
            f"Tgt: {m_tgt}{c_tgt}{l_tgt:8.4f}{Colors.RESET}"
            f"{val_str}"
        )
        print(log_msg)

        # 3중 특이점 도달 여부 판별
        if (l_res <= self.eps['L_res']) and (l_curv <= self.eps['L_curv']) and (l_tgt <= self.eps['L_tgt']):
            if not self.is_collapsed:
                self.is_collapsed = True
                self._announce_singularity(epoch, elapsed, l_res, l_curv, l_tgt)
            return True

        return False

    def _announce_singularity(self, epoch: int, elapsed: float, res: float, curv: float, tgt: float):
        """3중 공명 조건이 모두 허용 오차를 깨고 0으로 수렴했을 때 선포하는 장엄한 승리 로그"""
        border = f"{Colors.MAGENTA}{Colors.BOLD}================================================================================{Colors.RESET}"

        msg = f"""
{border}
{Colors.MAGENTA}{Colors.BOLD}
  ____  _____ ____   ___  _   _    _    _   _  ____ _____
 |  _ \| ____/ ___| / _ \| \ | |  / \  | \ | |/ ___| ____|
 | |_) |  _| \___ \| | | |  \| | / _ \ |  \| | |   |  _|
 |  _ <| |___ ___) | |_| | |\  |/ ___ \| |\  | |___| |___
 |_| \_\_____|____/ \___/|_| \_/_/   \_\_| \_|\____|_____|

  ____ ___  _     _        _    ____  ____  _____ ____
 / ___/ _ \| |   | |      / \  |  _ \/ ___|| ____|  _ \
| |  | | | | |   | |     / _ \ | |_) \___ \|  _| | | | |
| |__| |_| | |___| |___ / ___ \|  __/ ___) | |___| |_| |
 \____\___/|_____|_____/_/   \_\_|   |____/|_____|____/
{Colors.RESET}
{border}
{Colors.YELLOW}{Colors.BOLD}   [SYSTEM ALERT] 기하학적 특이점(Singularity) 돌파 감지!{Colors.RESET}
{Colors.WHITE}   경과 시간: {elapsed:.2f}초  |  돌파 에포크: Epoch {epoch}{Colors.RESET}
{border}
{Colors.CYAN}   >>> 위상(phi)이 공명 조건에 수렴했습니다. (Phase: {res:.6f}){Colors.RESET}
{Colors.BLUE}   >>> 곡률(A'')이 절대 안정을 향해 오목하게 붕괴했습니다. (Curvature: {curv:.6f}){Colors.RESET}
{Colors.GREEN}   >>> 목표(Ψ) 나침반과의 완벽한 동기화가 확인되었습니다. (Target: {tgt:.6f}){Colors.RESET}
{border}
   {Colors.GREEN}{Colors.BOLD}공명 조건 수렴 완료 (Resonance Convergence Achieved).{Colors.RESET}
{border}
"""
        print(msg)


# =================================================================
# 직접 실행 시 터미널 시각화 및 붕괴 선언 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    print("\n--- [RDL Observability] Resonance Dashboard Terminal Test ---")

    # 1. 엄격한 특이점 허용치(0.05) 설정
    dashboard = ResonanceDashboard(eps_res=0.05, eps_curv=0.05, eps_tgt=0.05)

    # 2. 에포크 진화 시뮬레이션
    total_epochs = 40
    for epoch in range(1, total_epochs + 1):
        # Loss가 10.0에서 0.01로 지수적 감소하는 상황 모사
        decay = math.exp(-epoch / 5.0)

        mock_metrics = {
            'L_res': 10.0 * decay + 0.005,  # 30에포크쯤엔 0.01 부근
            'L_curv': 5.0 * decay + 0.008,
            'L_tgt': 20.0 * decay + 0.002,
            'Total': 35.0 * decay + 0.015
        }

        # 화면 갱신 및 특이점 판별
        is_collapsed = dashboard.log_epoch(epoch, mock_metrics, val_loss=mock_metrics['Total']*1.1)

        # 실제 환경처럼 보이기 위한 짧은 딜레이
        time.sleep(0.15)

        # 3중 붕괴에 도달하면 조기 종료 (Early Stopping 모사)
        if is_collapsed:
            print(f"\n{Colors.CYAN}[OK] 공명 수렴 조건 달성. 조기 종료.{Colors.RESET}\n")
            break
