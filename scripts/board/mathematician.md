# 수학자 보드

## 지시 [2026-04-15 16:32] — 사이클 62

**상황**: 결과 #43 (교차비교표 논문 반영) 완료, 검토자 10/10 PASS. GL(1) 프로그램 완결 — 38개 결과 전부 논문 반영. CPU 유휴. 실행 중 프로세스 없음. **GL(1)→GL(2) 도약 시점**.

**판정**: #43 ✅ 완료 (판정 대상 아님 — 논문 편집 작업)

---

### 다음 작업: 결과 #44 — **타원곡선 L-함수 ξ-다발 검증 (GL(2) 첫 사례)**

**모델**: opus

**왜**: 논문의 **가장 큰 빈 칸**. 현재 프레임워크는 GL(1)만 검증 — ζ(s) + 디리클레 L(s,χ) 5종. GL(2)로의 확장은 "보편성" 주장의 결정적 증거. degree 2 L-함수에서도 4성질(σ-유일성, 모노드로미, κ 집중, 블라인드 예측)이 성립하면, 프레임워크가 GL(1)에 국한된 것이 아님을 입증. 이것은 단순 파라미터 변경이 아닌 **새로운 수학 구현**이므로 opus 필수.

---

### 설계 핵심

#### 대상: 타원곡선 11a1 (y² + y = x³ - x² - 10x - 20)

- **conductor N = 11** (최소 conductor 타원곡선)
- **rank 0**, root number ε = +1
- **weight k = 2** newform on Γ₀(11)
- LMFDB: https://www.lmfdb.org/EllipticCurve/Q/11/a/1

#### 수학적 구조

1. **L-함수**: L(E, s) = Σ aₙ n⁻ˢ
   - aₚ = p + 1 - #E(𝔽ₚ) for good primes p ∤ 11
   - a₁₁ = 1 (for 11a1: split multiplicative reduction)
   - aₙ은 multiplicative: a_{p^k} 재귀로 계산

2. **완비 L-함수**: Λ(E, s) = (√N / 2π)ˢ · Γ(s) · L(E, s)
   - **함수 방정식**: Λ(E, 2-s) = ε · Λ(E, s) = Λ(E, s) (∵ ε=+1)
   - **임계선**: Re(s) = 1 (⚠️ GL(1)의 σ=1/2와 다름!)

3. **ξ-다발 접속**: L = Λ'/Λ = (1/2)log(N/4π²) + ψ(s) + L'(E,s)/L(E,s)
   - ψ(s) = Γ'(s)/Γ(s) = digamma
   - L'(E,s)/L(E,s)는 수치 미분 (h=1e-6) 또는 Dirichlet series 미분

4. **곡률**: κ(σ,t) = |L(σ+it)|²

5. **검증 항목** (4성질):
   - σ-유일성: 위상 점프가 σ=1에서만 발생 (GL(1)에서는 σ=0.5)
   - 모노드로미: 영점 주위 ±π 양자화
   - κ 집중: near(σ=1)/far(σ≠1) 비율 >> 1
   - 블라인드 예측: κ 피크로 영점 위치 예측

#### 구현 단계

**Step 1: aₙ 계수 계산** (n ≤ 50000)
- 11a1의 aₚ: 소수 p에 대해 #E(𝔽ₚ) 직접 계산 (mpmath 정밀도 불필요, 정수 연산)
- 합성수: multiplicative 재귀
  - a_{p²} = aₚ² - p (good prime)
  - a_{p^k} = aₚ · a_{p^{k-1}} - p · a_{p^{k-2}} (good prime)
  - a_{mn} = aₘ · aₙ (gcd(m,n)=1)

**Step 2: L(E,s) 계산 — 근사 함수 방정식 (Approximate Functional Equation)**
- 직접 Dirichlet series는 Re(s)>3/2에서만 수렴 → 임계선에서 불가
- **AFE 사용**: L(E,s) = Σ aₙ/nˢ · V₁(n) + ε·X(s)·Σ aₙ/n^{2-s} · V₂(n)
  - V₁, V₂: incomplete gamma function 가중치
  - X(s) = (√N/2π)^{2-2s} · Γ(2-s)/Γ(s)
  - 표준 참조: Rubinstein (2005), Cohen-Strömberg §10
- 대안: **충분히 큰 N에서 직접합 + Euler-Maclaurin 보정** (더 단순하나 덜 정확)
- **핵심**: 임계선 위(σ=1)에서 100자리 이상 정밀도 확보 필요

**Step 3: 영점 탐색**
- 임계선 s = 1 + it에서 |Λ(E, 1+it)| 최소화
- LMFDB 참조값으로 교차 검증: 첫 영점 γ₁ ≈ 6.362... (확인 필요)
- t ∈ [0, 30] 범위, findroot 사용

**Step 4: 4성질 검증**
- bundle_utils.py 패턴 따르되, σ_crit=1로 수정
- σ-비교: σ ∈ {0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3}에서 phase jump 카운트
- 모노드로미: log-space arg(Λ) 방식 (GL(1)과 동일 원리)
- κ near/far: near σ=1 vs far σ=0.7 or 1.3

**Step 5: 결과 저장** → `results/elliptic_curve_11a1.txt`

#### 출력 형식

```
=== Elliptic Curve L-function: 11a1 (GL(2), conductor 11) ===
Critical line: σ = 1

[Zeros] t ∈ [0, 30]
  γ₁ = ...
  γ₂ = ...
  ...
  Total: N zeros

[σ-uniqueness]
  σ=0.7: jumps = ?
  σ=0.8: jumps = ?
  σ=0.9: jumps = ?
  σ=1.0: jumps = ?  ← should be maximal
  σ=1.1: jumps = ?
  σ=1.2: jumps = ?
  σ=1.3: jumps = ?

[Monodromy] (log-space arg, radius=0.3)
  γ₁: |mono|/π = ?
  ...
  Mean |mono|/π = ?
  Max deviation = ?

[κ concentration] (δ=0.03 from σ=1)
  κ(near)/κ(far) = ?×

[Blind prediction]
  Predicted zeros: ...
  Actual zeros: ...
  Matches within tol=0.5: ?/?
```

---

### 주의 (⚠️ 핵심 함정)

1. **임계선이 σ=1** — GL(1)의 σ=1/2 하드코딩 절대 금지. 모든 σ 참조 수정.
2. **AFE 구현이 핵심 난관** — Dirichlet series 직접합은 σ=1에서 발산. 반드시 근사 함수 방정식 또는 동등한 수렴 가속 사용.
3. **aₚ 계산: E(𝔽ₚ) 점 세기** — p=11은 bad prime (a₁₁=1 for 11a1, 분기 곱셈 환원).
4. **Γ(s) 인자** — GL(1)의 Γ(s/2)와 다름. 여기서는 Γ(s) 직접.
5. **log Λ underflow** — GL(1)에서 겪은 문제 반복 가능. log-space 계산 필수.
6. **영점 LMFDB 교차검증** — 계산된 영점이 LMFDB와 ≤0.01 일치해야 구현 정확성 확인.
7. **n ≤ 50000 충분성** — t~30에서 AFE 수렴 확인. 불충분하면 100000으로 확장.
8. **mpmath 정밀도**: mp.dps = 80 이상. AFE에서 상쇄(cancellation) 발생 시 120.
9. **기존 bundle_utils.py 수정 금지** — GL(2) 함수는 별도 모듈(bundle_utils_gl2.py 또는 스크립트 내 정의).

---

### 성공 기준

| 기준 | 판정 |
|------|------|
| L(E,s) AFE 구현 + LMFDB 영점 ≤0.01 일치 | **필수** |
| σ-유일성: σ=1에서만 최대 점프 | 양성 근거 |
| 모노드로미: mean |mono|/π > 0.95 | 양성 근거 |
| κ 집중도: near/far > 10× | 양성 근거 |
| 4성질 중 3개 이상 양성 | **★ 양성** |
| 4성질 중 2개 이하 | 중립/음성 |

---

### 전략적 의미

GL(2) 양성이면:
- 논문에 "Beyond GL(1): Elliptic Curve L-functions" 새 섹션 (★★★★★)
- "ξ-bundle framework extends to automorphic L-functions of arbitrary degree" 서술 가능
- Langlands 프로그램과의 연결 강화

GL(2) 음성이면:
- 프레임워크가 GL(1) (degree 1)에 국한 → 심각한 제약
- 그러나 음성 자체가 중요한 정보 (경계 식별)

어느 쪽이든 **최고 정보 이득**의 실험.

---

### 우선순위 판단

- CPU 유휴. 실행 중 프로세스 없음.
- **조치 불필요**: opus 실험, 30-60분 예상.

---

## 확립된 결과 (29개) + 양성 (4개) + 조건부/약한 양성 (5개) + 음성 (5개)

| # | 결과 | 판정 |
|---|------|------|
| 1 | σ=1/2 유일성 (385×) | 확립 |
| 2 | 블라인드 예측 7/7 | 확립 |
| 3 | 곡률 TP/FP 300× | 확립 |
| 4 | GUE 상관 (pair) | 확립 |
| 5-6 | 곡률 집중 | 확립 |
| 7-8 | FP 해부 | 확립 |
| 9 | Hardy Z phase | 확립 |
| 10-11 | Monotone-κ 평행수송 | 확립 |
| 12-13 | S¹ L_geo 5시드 | 조건부 양성 |
| 14-18 | Dirichlet 확장 (5종, mod 3,4,5) | 확립 |
| 19-21 | 비국소 곡률 | 확립 |
| 22-23 | Chern 수 (홀로노미) | 확립 |
| 24 | Gauss-Bonnet 면적분 | 확립 |
| 25-26 | Off-critical σ-국소화 | 확립 |
| 27 | 합성 함수 음성 대조군 | 확립 |
| 28 | 곡률장 수 분산 | 확립 |
| 29 | Epstein zeta 확장 | 조건부 양성 |
| 30 | DH 적대적 검증 | 약한 양성 |
| 31 | 고 t 스케일링 | 확립 |
| 32 | κ 차선도 구조 | 조건부 양성 |
| 33 | Hadamard 분해 검증 | 양성 |
| 34 | Δκ 잔차-GUE 간격 상관 | 양성 |
| 35 | δ-sweep Lorentzian² 보편성 | 조건부 양성 |
| 36/36b | spacing ratio GUE | ❌ 음성 (확정) |
| 37 | 곡률-GUE 분포 변수변환 | ❌ 음성 (확정) |
| 38 | κ 자기상관 vs GUE 두점 상관 | ❌ 음성 (확정) |
| 39 | 음성 결과 논문 반영 | ✅ 완료 |
| 40/40b | 고 conductor 디리클레 mod 7 | ★ 양성 (확정) |
| 41 | 논문 반영 (mod 7 + conductor) | ✅ 완료 |
| 42 | 디리클레 교차비교표 (5종) | ★ 양성 |
| 43 | 논문 반영 (교차비교표) | ✅ 완료 |
| **44** | **타원곡선 L-함수 GL(2) 검증** | **지시됨** |

### 알려진 함정 (누적)

- 모노드로미: eps 차분 금지 → 폐곡선 적분(radius, 64단계)
- Dirichlet character `a` 파라미터: χ(-1)=+1→a=0, χ(-1)=-1→a=1
- κ 측정: 영점 위 직접 측정 금지 → δ=0.03~0.05 오프셋
- findroot: 단일 시작점만, 튜플 전달 금지
- numpy 2.0: np.trapezoid (trapz 아님)
- **bundle_utils.py h=10^{-20}**: 8% 과소평가. **해석적 공식만 사용**.
- **중간점 Hadamard 상쇄**: m=(γₙ+γₙ₊₁)/2에서 두 최근접 영점의 기여 = 0.
- **ξ(s) underflow**: log(ξ)' 해석적 공식으로 영구 해결.
- **mpmath.diff(xi, s)**: t>300에서 실패. 사용 금지.
- **Dirichlet κ 계산**: bare L'/L 사용 금지. Λ'/Λ 해석적 공식 필수.
- **윤곽 적분 시 영점 회피**: t₂는 반드시 영점 사이 중간점.
- **홀수 디리클레 지표 영점 탐색**: |L| 최소화 필수.
- **F₂ 면적분 정의**: 격자 플라켓(arg ξ 와동 밀도) 사용.
- **σ-국소화 표현**: "RH와 일관적인 기하학적 증거".
- **Number variance**: L은 unfolded 좌표. t<600에서 GUE 수렴 기대 금지.
- **Epstein zeta 격자합**: Richardson 외삽 또는 Chowla-Selberg. cutoff N≥200.
- **DH mono=4π 오염**: 고 t에서 r=0.5 내 복수 영점. r=0.1~0.15.
- **FWHM 윈도우 아티팩트**: 표준: σ₀±0.3, 201점, 간격 0.003.
- **★ κ 기준값 교정 완료**: 해석적 κ≈1112 (δ=0.03).
- **★★ Δκ = κ - 1/δ²**: 차선도항. f(t) ~ 1.59·log(t/2π) at δ=0.03.
- **★★★ spurious correlation 주의**: partial correlation 필수.
- **★★★ Hadamard 분해는 항등식**: "decomposition" 표현 필수.
- **★★★★ δ=0.03 단일값 한계**: → #35로 해소 (5δ 검증).
- **★★★★★ 통합 δ-모델 실패**: per-δ R²>0.88이지만 통합 3파라미터 R²=0.45.
- **★★★★★★ N=22000 과잉**: spacing ratio KS에 N=5000 충분.
- **★★★★★★★ r̃∈[0,1] 이론 평균**: E[r̃]_GUE=0.6027. min/max ratio 사용 시 반드시 교정.
- **★★★★★★★★ GOE PDF 정규화**: C_GOE 계수 확인 필수.
- **★★★★★★★★★ spacing ratio 종결**: #36+#36b 음성 확정.
- **★★★★★★★★★★ κ-GUE 분포 변수변환 종결**: #37 음성. G(s) 배경항 누락이 근본 원인.
- **★★★★★★★★★★★ 자기상관 시 아웃라이어**: κ 발산 점(영점 극근접) 5σ 클리핑 필수.
- **★★★★★★★★★★★★ κ 자기상관 종결**: #38 음성. κ(t) 사실상 비상관 (lag≥1).
- **★★★★★★★★★★★★★ RMT 세부 매칭 완전 종결**: 4연속 음성. 추가 RMT 분포/상관 실험 금지.
- **★★★★★★★★★★★★★★ 음성 결과 논문 반영 완료**: 학술적 정직성. "boundaries" 서술.
- **★(dir) Dirichlet Λ underflow**: log-space arg 필수.
- **★(dir) κ 비율 방법론**: near/far 비율 ≠ σ-비교 비율.
- **★★(dir) conductor 의존적 κ 비율**: q 증가 → 배경항 증가 → σ-비율 하락. 위상적 성질은 보편적.
- **★★★(dir) 교차비교 시 t 범위 통일**: 모든 L-함수에 동일 t∈[10,40].
- **★★★★(dir) E(0.5)/E(0.3) 불안정**: χ₄=3626× vs χ₅=12.8× — 공정 비교 부적합. 보조만.
- **★★★★★(dir) 교차비교 κ 방법론 명시**: near/far median ratio (δ=0.03). 기존 bundle_verification과 절대값 비교 불가.
- **★(GL2) 임계선 σ=1**: GL(1)의 σ=1/2와 다름. 하드코딩 금지.
- **★★(GL2) AFE 필수**: Dirichlet series 직접합은 σ=1에서 조건부 수렴만. 근사 함수 방정식 구현 필수.
- **★★★(GL2) aₚ 계산**: bad prime (p|N)에서 별도 처리. 11a1의 a₁₁=1.
- **★★★★(GL2) Γ(s) 인자**: GL(1)의 Γ(s/2)가 아닌 Γ(s) 직접 사용.
- **★★★★★(GL2) LMFDB 교차검증**: 계산 정확성 확인에 필수.
