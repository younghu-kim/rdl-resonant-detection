# GDL Unified — Geometric Deep Learning (5G + Time) + Resonant Deep Learning

Bronstein et al.의 기하학적 딥러닝 논문을 실행 가능한 코드로 구현하고,
공명소수법칙(RDL)을 5G Gauge 도메인의 구체적 인스턴스로 통합한 프레임워크.

## 구조

```
gdl/
├── core/                  GDL Blueprint ABC + GDLFactory
│   ├── base_domain.py     BaseDomain
│   ├── base_signal.py     BaseSignal
│   ├── base_transform.py  BaseTransform
│   ├── base_blueprint.py  EquivariantLayer, NonLinearity, Pooling, Invariant, BaseGDLModel
│   └── factory.py         GDLFactory — 8개 도메인 자동 조립
│
├── domains/               5G + Time 도메인 구현
│   ├── grids/             1G 평행이동 등변 (CNN)
│   ├── groups/            2G 회전 등변 (Spherical CNN, Wigner-D)
│   ├── graphs/            3G 순열 등변 (DeepSets, GNN, Transformer)
│   ├── geodesics/         4G 등거리 등변 (Intrinsic Mesh CNN, Chebyshev)
│   ├── gauges/            5G 게이지 등변 (SO(2) Parallel Transport)
│   └── temporal/          Time 시간워핑 등변 (Difference-Gated RNN)
│
├── rdl/                   공명 딥러닝 (Resonant Deep Learning)
│   ├── constants.py       H_DIFF, H_GAUGE, H_LIE, PrecisionManager
│   ├── core/              복소 연산, 다원군 프레임, 릿지 솔버
│   ├── layers/            임베딩, PQO (sin²/cos²), 로그 미분, 공명 블록
│   ├── dynamics/          게이지 ODE (분리기, 감쇠, Heun 적분기, 상태 버퍼)
│   ├── losses/            F₂ 판별식, 곡률, 타겟, PQO, 4항 통합 손실
│   ├── optim/             PGGD (Phase-Gauge Gradient Descent)
│   ├── models/            MasterResonantNetwork
│   └── pipeline/          XiFeatureDataset, ResonanceTrainer
│
├── bridge/                RDL ↔ GDL 5G 통합
│   ├── gauge_rdl_adapter.py   MasterResonantNetwork → EquivariantLayer[GaugeSignal]
│   └── rdl_gauge_domain.py    U(1) 시간격자 → GaugeDomain 매핑
│
└── viz/                   시각화 (위상, 곡률, 대시보드)
```

## 수학적 배경

### GDL Common Blueprint (Bronstein et al. 2021)

모든 기하학적 도메인에서 동일한 아키텍처:

```
f = A ∘ σ_J ∘ B_J ∘ P_{J-1} ∘ ... ∘ P_1 ∘ σ_1 ∘ B_1
```

- **B**: EquivariantLayer — 대칭 보존 정보 혼합
- **σ**: GeometricNonLinearity — 기하학 보존 비선형 활성화
- **P**: LocalPooling — 스케일 분리 (조립화)
- **A**: GlobalInvariant — 도메인 파괴, 불변 리드아웃

### RDL ↔ GDL 5G 대응

| GDL 5G (Gauge) | RDL (Resonant) |
|---|---|
| SO(2) 구조군 | U(1) 위상 회전 |
| 평행 이동 R(ρ_ij) | 게이지 ODE 적분 |
| Schur 가중치 [[a,-b],[b,a]] | 복소 선형 연산 |
| 노름 게이팅 | 진폭-위상 분리 |

### RDL 핵심 수식

- **PQO**: P(φ;L) = sin²(Lφ/2), Maslov 보정: cos²(Lφ/2)
- **F₂ 판별식**: F₂(t) = Im{ e^{-iφ} (L̂ - φ̇) } = 0
- **게이지 ODE**: ψ_{j+1} = (1-α)ψ_j + α·v_t, φ_{j+1} = φ_j + (dt/2)(ψ_j + ψ_{j+1})

## 설치

```bash
# 가상환경 (torch + mpmath 필요)
pip install torch numpy scipy mpmath einops tqdm

# 개발 의존성
pip install pytest pytest-cov ruff mypy
```

## 실행

```bash
# 테스트 (115개)
python -m pytest tests/ -q

# GDL 8개 도메인 데모
python scripts/run_gdl_demo.py

# RDL ↔ GDL 5G 브릿지 데모
python scripts/run_bridge_demo.py

# ξ-함수 영점 탐지 훈련 (mpmath 필요)
python scripts/run_xi_pipeline.py --optimizer adam --pqo-mode cos2 --epochs 200

# 고위 영점 (t ∈ [100,200], 50개 영점)
python scripts/run_xi_pipeline.py --optimizer adam --pqo-mode cos2 --epochs 200 \
  --t-min 100.0 --t-max 200.0 --hidden 64
```

## 실험 결과

| 범위 | 영점 수 | 탐지율 | cos² PQO | 옵티마이저 |
|---|---|---|---|---|
| t ∈ [10,50] | 10 | **10/10** | < 0.0002 | Adam |
| t ∈ [100,200] | 50 | **50/50** | < 0.0011 | Adam |
| t ∈ [10,50] | 10 | 0/10 | — | PGGD |

PGGD 실패는 버그가 아니라 ξ-함수 손실 지형이 위상-자명(phase-trivial)이기 때문.

## 소스 출처

- **GDL 도메인**: `~/Desktop/Goun_Greed/` AI 설계서에서 코드 추출
- **RDL 모듈**: `~/Desktop/XI_function/` 동작 코드에서 이전 + import 정규화
- **Bridge**: 새로 작성

## 라이선스

Kim Young-Hu (2026). 연구용.
