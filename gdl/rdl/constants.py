"""
=============================================================================
[Project RDL] Resonant Deep Learning - Global Constants & Precision Manager
=============================================================================
공명 수론 프레임워크의 핵심 상수를 불변(Immutable) 객체로 관리하며,
모든 텐서 연산을 Float64/Complex128 정밀도로 강제한다.

논문 참조:
  - main.tex: QROP-Net (격자 암호 + 프레넬 광학 + 딥러닝 위상 복원)
  - unified_en.tex: 통합 프레임워크 (위상 양자화 원리)
"""

import os
import torch
import dataclasses
from dataclasses import dataclass

# =====================================================================
# 1. 전역 하이퍼파라미터 객체 (Frozen Dataclass)
# =====================================================================
@dataclass(frozen=True)
class ResonantConstants:
    """
    공명 딥러닝 시스템의 핵심 상수. frozen=True로 런타임 변경을 차단한다.
    각 상수는 논문 방정식 번호 또는 [확장] 표기로 출처를 명시한다.
    """
    # ---------------------------------------------------------
    # [1] 이상 위상 함수 Psi(n) 상수 (unified_en.tex, eq:psi_n)
    #     Psi(n) = alpha * log(n) + beta * sin(gamma * log(n)) / n^delta
    # ---------------------------------------------------------
    ALPHA: float = 1.0          # alpha: 대수 상승 기울기 (alpha > 0)
    BETA: float  = 1.0          # beta: 진동 진폭 (beta >= 0)
    GAMMA: float = 2.5          # gamma: 위상 회전 주파수 (gamma > 0)
    DELTA: float = 1.0          # delta: 진동 감쇠율 (delta > 0)

    # ---------------------------------------------------------
    # [2] 분리된 스텝 상수 (용도별 독립)
    # ---------------------------------------------------------
    H_DIFF: float  = 1.5e-7     # 중심차분 미소단계 (eq:central_diff)
    H_GAUGE: float = 0.04       # 게이지 ODE 적분 단계 dt_j (eq:ode_alpha)
    H_LIE: float   = 1e-4      # [확장] 리 군 지수 사상 단계 Z * exp(h * L_G)

    # 하위 호환용 (레거시 코드에서 R_CONST.H 참조 시)
    H: float = 1.5e-7           # [deprecated] H_DIFF로 대체 예정

    # ---------------------------------------------------------
    # [3] 다원군(Multigroup) 채널 구성
    # ---------------------------------------------------------
    RHO: float = 3e-3           # rho: 릿지 정규화 계수 (과결정 채널용)

    # 논문 권장 3-채널 구성 (eq:channel_set)
    PAPER_CHANNEL_ANGLES: tuple = (-30.0, 15.0, 60.0)
    PAPER_CHANNEL_LAMBDAS: tuple = (1.10, 0.85, 1.00)    # 비등방 축척 lambda_k
    PAPER_CHANNEL_WEIGHTS: tuple = (0.40, 0.35, 0.25)    # 정규화 가중치 sum(w_k)=1

    # ---------------------------------------------------------
    # [4] 게이지 상태 동역학 (eq:ode_alpha ~ eq:ode_phi)
    # ---------------------------------------------------------
    TAU_G: float = 0.04         # tau_g: 게이지 감쇠 상수 (unified_en.tex, tau=0.04)
    KAPPA: float = 1.0          # [확장] kappa: 상태 의존 감쇠 계수 (논문 외 독자적 기여)

    # ---------------------------------------------------------
    # [5] 위상 양자화 연산자 P(phi;L) = sin^2(L*phi/2) (eq:pqo)
    # ---------------------------------------------------------
    PQO_L_DEFAULT: int = 2      # L=2 -> Gate A: sin^2(v)=0 at v=m*pi

    # ---------------------------------------------------------
    # [6] 수치 안정성
    # ---------------------------------------------------------
    EPSILON: float = 1e-16      # epsilon: 영분모/특이점 방지 상수


# 전역에서 호출할 단일(Singleton) 인스턴스 생성
# 파이썬의 모듈 캐싱 특성에 의해, 어디서 임포트하든 동일한 메모리 주소의 R_CONST를 참조
R_CONST = ResonantConstants()


# =====================================================================
# 2. 정밀도 및 결정론적 연산 매니저
# =====================================================================
class PrecisionManager:
    """
    딥러닝 프레임워크의 정밀도 희생을 방지하고,
    64-bit 부동소수점과 128-bit 복소수 연산을 강제한다.
    """
    REAL_DTYPE = torch.float64
    COMPLEX_DTYPE = torch.complex128

    @classmethod
    def setup_precision(cls) -> None:
        """프로그램 진입 시 호출하여 정밀도 환경을 고정한다."""

        # 1. 기본 텐서 타입을 Float64로 강제
        torch.set_default_dtype(cls.REAL_DTYPE)

        # 2. TF32 가속 차단 (10-bit 정밀도는 미세 위상차 연산에 부적합)
        if torch.cuda.is_available():
            torch.backends.cuda.matmul.allow_tf32 = False
            torch.backends.cudnn.allow_tf32 = False

        # 3. 결정론적 연산 강제
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

        # CUBLAS 환경 변수 설정
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
        torch.use_deterministic_algorithms(True, warn_only=True)

        print("[RDL Precision Manager] Float64/Complex128 고정 완료")


# 모듈 임포트 시 자동으로 정밀도 환경 세팅 실행
PrecisionManager.setup_precision()


# =================================================================
# 직접 실행 시 무결성 테스트
# =================================================================
if __name__ == "__main__":
    print("\n--- [RDL Constants & Precision Test] ---")

    # 1. 상수 불변성 테스트
    try:
        R_CONST.ALPHA = 2.0
        print("[FAIL] 상수가 변경되었습니다!")
    except getattr(dataclasses, 'FrozenInstanceError', Exception):
        print("[OK] 상수는 Frozen 상태로 보호됨")

    # 2. 텐서 정밀도 승격 테스트
    dummy_real = torch.tensor([1.0, 2.0])
    dummy_complex = torch.complex(dummy_real, dummy_real)

    if dummy_real.dtype == torch.float64 and dummy_complex.dtype == torch.complex128:
        print("[OK] Float64/Complex128 정밀도 통제 확인")
    else:
        print("[FAIL] 정밀도 통제 실패")

    # 3. 상수 분리 테스트
    assert R_CONST.H_DIFF != R_CONST.H_GAUGE, "H_DIFF와 H_GAUGE가 같으면 안됨"
    assert R_CONST.H_DIFF != R_CONST.H_LIE, "H_DIFF와 H_LIE가 같으면 안됨"
    assert len(R_CONST.PAPER_CHANNEL_ANGLES) == 3, "논문 채널은 3개"
    assert abs(sum(R_CONST.PAPER_CHANNEL_WEIGHTS) - 1.0) < 1e-10, "가중치 합은 1"
    print("[OK] 상수 분리 및 논문 채널 구성 검증 완료")
    print(f"  H_DIFF={R_CONST.H_DIFF}, H_GAUGE={R_CONST.H_GAUGE}, H_LIE={R_CONST.H_LIE}")
    print(f"  채널 각도: {R_CONST.PAPER_CHANNEL_ANGLES}")
    print(f"  채널 축척: {R_CONST.PAPER_CHANNEL_LAMBDAS}")
    print(f"  채널 가중치: {R_CONST.PAPER_CHANNEL_WEIGHTS}")
    print(f"  PQO L 기본값: {R_CONST.PQO_L_DEFAULT}")
