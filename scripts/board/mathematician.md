# 수학자 보드
## 최종 업데이트: 2026-04-13

## 현재 수학적 판단

### 확립된 결과
1. **Kuramoto 동기화 확인 (양성)**: r_final=0.9994±0.0001, std(φ) ratio 3.43±0.09x
   - 논문 Obs 8.6 정량적 재현 완료
   - **핵심 발견**: 전이(transition)가 아닌 정련(refinement). r_init≈0.995 (에포크 0부터)
   - von Mises 자기일관성: σ=0.032 → r=exp(-σ²/2)=0.9995 ≈ 실측 0.9994
   - t-범위 불변 확인: [10,50]와 [100,200]에서 동일 패턴

2. **S¹ geodesic 통합 (양성)**: 검출 5.7→22.7/463 (4배), F₂ ratio 0.913
   - L_geo = 1 - cos(Δφ) = von Mises NLL (상수 제외)
   - L_geo 최소화 ↔ von Mises MLE ↔ κ 최대화 ↔ 동기화 강화
   - ±π 불연속 제거 + 통계적으로 자연스러운 손실 → 이중 이점

3. **고높이 스케일링 (양성)**: t=1100까지 recall 97.5%+, precision 17→26%
   - K=128 충분: Nyquist K≥W/(2Δt) 만족

### 새 명제 후보
**Prop: Spontaneous Gauge Fixing**
> H개 은닉 채널 위상 {φ_j(x)}는 학습 전 r_init ≈ 1 - O(1/H)로 시작.
> 학습은 κ를 O(H) → O(H²)로 정련. U(1)^H → U(1) 자발적 대칭 파괴.
> 검증 필요: H-스케일링 실험 (H=16,32,128,256에서 r_init 측정)

### 논문 수정 제안
- Obs 8.6의 "~10 에포크 전이"를 "supercritical 정련"으로 재해석
- 실측에서 r>0.99은 에포크 0부터. 전이는 초기화에서 이미 발생.

## 다음 수학적 질문 (우선순위순)

### 1순위: r_init ≈ 0.995의 기원 (H/N 스케일링)
- 가설 A: Xavier+ReLU 구조적 → r_init ~ 1 - c/H (실수 가중치가 위상 분산 억제)
- 가설 B: 유한 크기 우연 → r_init 불규칙
- 실험: H=16,32,64,128,256에서 r_init 측정. 로그-로그 플롯에서 기울기 확인
- Marchenko-Pastur 분포 연결 가능성 탐색

### 2순위: FP 해부 (precision ~25% 원인)
- |F₂|''(t) 부호: 극소점 vs 영점 구분 가능한가?
- FP→nearest true zero 거리 분포 vs GUE pair correlation
- 정보론적 하한 추정: 해상도 ~ (b-a)/(2K) ≈ 0.39

### 3순위: L_geo ↔ 동기화 인과성
- L_geo 학습 시 r(epoch) 궤적 vs L_tgt 궤적 비교
- 가설: L_geo=von Mises NLL → 더 빠른 동기화 → 더 나은 검출
- 빠른 실험, 1순위와 병렬 가능

## 이론적 메모

### von Mises–Kuramoto–L_geo 삼각 연결
```
von Mises 분포: p(φ) ∝ exp(κ cos(φ - μ))
  ↕ MLE
L_geo = 1 - cos(Δφ)  (음의 로그우도)
  ↕ 최소화
Kuramoto r → 1  (순서 매개변수 최대화)
```
세 개념이 수학적으로 동일한 최적화 문제. L_geo를 쓰면 네트워크는
자연스럽게 von Mises MLE를 수행하며 Kuramoto 동기화를 강화한다.

### 초기 동기화 기원 추측
Xavier 초기화에서 W_{ij} ~ N(0, 1/n_in). 3층 ReLU 네트워크의 출력은
z_j = Σ w_{jk} ReLU(Σ ...) 형태. ReLU가 음수를 제거하므로 Im(z_j)/Re(z_j)의
분산이 1/H에 비례하여 억제될 수 있음. 이는 arg(z_j)의 분산이 O(1/H)임을 의미.
→ r = 1 - O(1/H²) 예측. H=64이면 1 - O(1/4096) ≈ 0.9998 (실측 0.995와 order 일치).
정밀한 계산에는 Marchenko-Pastur 고유값 분포 필요.
