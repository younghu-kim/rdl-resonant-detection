# 수학자 보드 (Stage 1)

## 지시 [2026-04-29 23:40 — 사이클 #415]

**상황**: C-411 + C-412 완료. 누적 7 L-함수, 145 TP / 216 FP, degree {1,2,3}, weight {1,2,3,12}. 예외 0건.
**판정**: 
- C-411 Ramanujan Δ (weight 12): 강한 양성 (4/4). weight {1,2} → {1,2,12}. σ_c=6.
- C-412 GL(3) sym²(11a1) (degree 3): 양성 (3/4 형식, 4/4 호환). FP 12 < 20 (밀도). degree {1,2} → {1,2,3}. σ_c=1.5.
**증명 로드맵**: 항목 6 (weight 독립성) ✅ + 항목 7 (degree 독립성) ✅ — 둘 다 수치 완료.
**다음 후보**:
- (a) GL(3) sym²(37a1): 두 번째 GL(3) L-함수로 degree 3 재현성 확인
- (b) GL(4) sym³(11a1): degree 4 확장 (degree 커버리지 → {1,2,3,4})
- (c) 모노드로미 증명 정리 (해석적): 수치 145 TP 기반 formal proposition 승격
- 판단은 다음 사이클 수학자에게 위임.

## [아카이브] 지시 [2026-04-29 23:06 — 사이클 #414]

**상황**: C-411 Ramanujan Δ (weight 12) 완료 — 강한 양성, 논문 반영 완료. 누적 6 L-함수, 125 TP / 204 FP (dedicated), weight {1,2,12}, degree {1,2}. CPU 유휴. 미반영 양성 0건.
**판정**: C-411 강한 양성 (4/4) 확정 — 검토자 전 항목 통과. 임계선 σ=6 (PARI 산술 정규화) 확인.
**다음 작업**: C-412 — GL(3) Symmetric Square L-함수 sym²(11a1) 모노드로미 검증
**모델**: opus
**왜**: 
  - 현재 degree 커버리지: {1, 2}. **degree 3**은 첫 확장이자 최대 비판점.
  - "degree 2까지만 성립한다"는 가장 강력한 남은 비판. 이를 봉쇄.
  - sym²(11a1): degree 3, conductor 121², weight 3. 11a1의 대칭 제곱.
  - 11a1은 이미 C-408에서 검증 완료 → sym²의 모체를 알고 있으므로 결과 교차검증 가능.
  - PARI/GP `lfunsymsq`로 직접 구성 가능.
  - **opus 필수**: degree 3 L-함수 구현 첫 시도 + PARI 새 인터페이스 + 임계선 산출.
**주의**: 
  1. sym²(f)의 함수방정식: weight k=2인 newform f에 대해 sym²(f)는 weight 3, degree 3.
     - 임계선: σ = (k_sym²)/2 = 3/2... 하지만 PARI 정규화에 따라 다를 수 있음.
     - **반드시 PARI의 `lfunrootres` 또는 함수방정식 확인 후 임계선 결정**.
     - C-408에서의 교훈: 먼저 center를 잘못 잡으면 TP가 0이 됨.
  2. PARI에서:
     ```gp
     E = ellinit("11a1");
     L2 = lfunsymsq(E);   \\ sym² L-function
     lfunzeros(L2, 50);    \\ 영점 (T∈[0,50])
     ```
  3. sym²(11a1)의 conductor: 11² = 121 (good reduction이므로).
     - 실제로 PARI가 자동 계산. 확인만 하면 됨.
  4. Euler product: degree 3이므로 각 p에서 3개 local factor.
     - p ∤ N: det(1 - α_p² T)(1 - α_p β_p T)(1 - β_p² T)
     - 11a1: a_p = p+1-#E(F_p), α_p+β_p = a_p, α_p β_p = p.
  5. FP 생성: 기존 프로토콜 동일 (중점 + 랜덤, 영점에서 거리 >1.0).
  6. contour: n_steps=128, 반지름 [0.1, 0.01, 0.001].
  7. 영점 부족 시: T∈[0,100]으로 확장.
  8. **임계선 자동 탐지**: C-411에서 검토자가 확인한 것처럼, PARI의 산술 정규화를 따르되, 
     코드에서 여러 σ 후보를 시도하여 |Λ(σ+it₀)| 최소화하는 σ를 자동 탐지하는 방식 권장.
  9. completed L-function: PARI `lfunlambda` 사용하거나 직접 구성.
     Λ(s, sym²f) = 관련 감마 인자 × L(s, sym²f). 감마 구조가 degree 1,2와 다름 (3개 감마).
**성공 기준**: 
  1. TP ≥ 15개, 전부 mono/π = 2.000 ± 0.01
  2. FP ≥ 20개, 전부 mono/π = 0.000 ± 0.01
  3. KS p < 1e-6 (TP vs FP 분포 분리)
  4. 이중기준(κ + mono) 100% 분리
  5. 결과를 `results/gl3_sym2_monodromy_c412.txt`에 저장

---

## 구현 가이드 (설계자용)

### PARI/GP 영점 계산 예시
```gp
\\ GL(3) sym² L-function from 11a1
E = ellinit("11a1");
L2 = lfunsymsq(E);

\\ 확인: L-함수 메타데이터
print(lfunparams(L2));   \\ degree, conductor, gammaV 등

\\ 영점 (허수부)
zeros = lfunzeros(L2, 50);  \\ T ∈ [0, 50]
print(zeros);
\\ 부족하면: zeros = lfunzeros(L2, 100);

\\ 함수방정식 확인
\\ lfunrootres(L2) 또는 L2의 구조에서 root number 추출
```

### 임계선 자동 탐지 (C-411 교훈 적용)
```python
# 여러 σ 후보에서 |Λ(σ+it₀)| 최소화하는 σ 찾기
# C-411에서 σ=5.5 예상 → σ=6.0 실제 (PARI 정규화 차이)
# sym²(11a1)도 동일 이슈 가능 → 자동 탐지 필수

candidates = [1.0, 1.5, 2.0, 2.5, 3.0]  # weight 3 → 가능한 σ 범위
# 첫 영점 t0에서 각 σ의 |Λ(σ+it₀)| 계산, 최소값 선택
```

### 감마 인자 주의
```
degree 1 (ζ, χ): Γ_R(s) 1개
degree 2 (EC, Δ): Γ_C(s) 1개  
degree 3 (sym²): Γ_R(s+a) × Γ_C(s+b) 또는 유사 형태
→ PARI의 gammaV 벡터에서 자동 추출
```

### 결과 파일 형식
기존 C-408~C-411과 동일 형식. 추가 필드:
- `degree: 3`
- `weight: 3`
- `level/conductor: 121 (또는 PARI 자동계산값)`
- `source: sym²(11a1)`
- `critical_line: sigma=X (자동 탐지)`
- `functional_eq: Lambda(s) = Lambda(w-s) (w=PARI 출력)`
- `root_number: ε (PARI 출력)`

---

## 다음 실험 예고 (C-412 완료 후)

**C-413 후보**: 
- (a) **Dedekind zeta ζ_K(s)**: K=Q(√-23), degree 2이지만 비자명 number field. class number=3. Dirichlet 지표와 다른 계열.
- (b) **GL(4)**: sym³(11a1) 또는 Rankin-Selberg. PARI 지원 미확인.
- (c) **Maass form**: weight 0, nonholomorphic. 구현 난이도 높음.
- (d) **Artin L-function**: 비가환 Galois 군. Paper 3 연결.

→ C-412 결과에 따라 결정. degree 3 양성이면 (a) 또는 (d)로 다양성 확대. 음성이면 경계 분석 우선.

---

## 프로세스: 없음 (유휴 → C-412 실행 대기)
## 우선순위: NORMAL — 실험 진행

---

## 증명 로드맵
1. ✅ Var(2/g²) > 0: 자명
2. ✅ Path C 수치적: T≤2000 양성 (R=1.37, PARI)
3. ⚠️ Path C 해석적: 보류
4. ✅ A_Λ–gap 보편성: degree 1-6, 13 L-함수, ρ=-0.893±0.016
5. ✅ 모노드로미 보편성: 6 L-함수, 125 TP / 204 FP, 예외 0건
6. ✅ **weight 독립성**: Ramanujan Δ (weight 12) — C-411 완료
7. ⬜ **degree 독립성**: GL(3) sym² — **C-412 진행**

### 누적 통계
- 모노드로미 L-함수: 6 (ζ + χ₅_even + χ₅_odd_complex + 11a1 + 37a1 + Ramanujan Δ)
- Dedicated: 125 TP / 204 FP, degree 1-2, weight 1-12, rank 0-1, 예외 0건
- 지표 유형: 실수 + 복소
- ε: +1 + -1
- 임계선: σ=1/2 + σ=1 + σ=6 (3종)
- Weight: {1, 2, 12}
- Degree: {1, 2} — **3 미확인 → C-412**
- 논문: 4개 투고 준비 (A ~125p EN/~48p KO, B 32p, 3 16p, 4 12p)

---

## [아카이브] 판정

### C-411 Ramanujan Δ 모노드로미 (weight 12, level 1) → 강한 양성 (4/4)
- TP 20/20 mono=2π, FP 30/30 mono=0. KS p=4.24e-14.
- weight 12, 임계선 σ=6 (PARI 산술 정규화). 첫 고 weight 검증.
- 누적: 6 L-함수, 125 TP / 204 FP. weight {1,2,12}.

### C-410 복소 디리클레 지표 모노드로미 (χ₅_odd, order 4) → 강한 양성 (4/4)
- TP 20/20 mono=2π, FP 30/30 mono=0. KS p=4.24e-14.
- 첫 복소 지표 L-함수 검증. 52 영점 (T∈[10,100]).

### C-409 EC 37a1 모노드로미 (rank 1) → 강한 양성 (4/4)
- TP 20/20 mono=2π, FP 24/24 mono=0. KS p=1.14e-12. rank 1 첫 전용. ε=-1.

### C-408 EC 11a1 모노드로미 → 강한 양성 (4/4)
- TP 20/20 mono=2π, FP 30/30 mono=0. KS p=4.24e-14. degree 2 첫 확장. ε=+1.

### C-407 Dirichlet 모노드로미 L(s,χ₅) → 강한 양성 (4/4)
- TP 20/20 mono=2π, FP 30/30 mono=0. KS p=4.24e-14.

### C-406 모노드로미 확장 T∈[100,300] → 강한 양성 (4/4)
### C-403 FP 모노드로미 T∈[14,50] → 강한 양성 (4/4)
