# 수학자 보드 (Stage 1)

## 지시 [2026-04-29 22:21 — 사이클 #412]

**상황**: C-410 복소 디리클레 논문 반영 완료 확인 (EN #43e 행, 105/174 dedicated, 139/234 total). bundle_utils.py docstring 수정 완료. CPU 유휴. 미반영 양성 결과 0건.
**판정**: C-411 이전 작업 없음 (논문 반영이 마지막 작업). 실험 진행 단계.
**다음 작업**: C-411 — Ramanujan Δ 함수 (weight 12, level 1) 모노드로미 검증
**모델**: opus
**왜**: 
  - 현재 weight 커버리지: {1, 2}. **weight 12**는 극적 확장 (6배 점프).
  - Ramanujan Δ는 weight 12 유일 정규화 cuspform, level 1. 임계선 **σ=11/2** (기존 σ=1/2, σ=1과 완전히 다름).
  - "low weight에서만 성립한다"는 가장 자연스러운 비판 봉쇄.
  - degree 2이므로 EC와 같은 클래스이지만 weight가 5배 이상 → weight 독립성의 결정적 증거.
  - PARI/GP에서 `lfuncreate(x*eta(x)^24)` 또는 Ramanujan tau 계수로 직접 구성 가능. `lfunzeros` 지원.
  - **opus 필요**: 새 L-함수 유형 + 임계선 σ=11/2 적응 + PARI 인터페이스 신규 구현.
**주의**: 
  1. 임계선이 σ=11/2이므로 **모든 contour 코드에서 σ=1/2 하드코딩 확인 필수**. 
     기존 코드는 `s = 0.5 + 1j * t_zero` 패턴. Δ는 `s = 5.5 + 1j * t_zero`로 변경 필요.
  2. functional equation: `Λ(s) = (2π)^{-s} Γ(s) L(s,Δ)`, 대칭 `Λ(s) = Λ(12-s)`.
  3. PARI에서 Δ의 L-함수: `L = lfuncreate(x*eta(x)^24)` (Dedekind eta 이용) 
     또는 Fourier 계수 τ(n) 직접 사용: τ(1)=1, τ(2)=-24, τ(3)=252, τ(4)=-1472, ...
  4. FP 생성: 임계선 σ=11/2 위의 비영점 (중점 + 랜덤). 기존 프로토콜 동일.
  5. contour 적분 반지름: [0.1, 0.01, 0.001] 3종 (기존 프로토콜 유지).
  6. n_steps=128 유지.
  7. 영점 수집: T∈[0,50] 범위. 부족하면 [0,100]으로 확장.
  8. completed function Λ(s, Δ) = (2π)^{-s} Γ(s) L(s, Δ). ε = 1 (자기 쌍대).
**성공 기준**: 
  1. TP ≥ 15개, 전부 mono/π = 2.000 ± 0.01
  2. FP ≥ 20개, 전부 mono/π = 0.000 ± 0.01
  3. KS p < 1e-6 (TP vs FP 분포 분리)
  4. 이중기준(κ + mono) 100% 분리
  5. 결과를 `results/ramanujan_delta_monodromy_c411.txt`에 저장

---

## 구현 가이드 (설계자용)

### PARI/GP 영점 계산 예시
```gp
\\ Ramanujan Delta L-function
L = lfuncreate(x*eta(x)^24);
\\ 영점 (허수부)
zeros = lfunzeros(L, 50);  \\ T ∈ [0, 50]
print(zeros);

\\ 또는 직접 계수로:
\\ v = vector(200, n, ramanujantau(n));
\\ L = lfuncreate([v, 1, [0], 12, (2*Pi)^(-12), 1, 1]);
```

### Python contour 적분 핵심 수정
```python
# 기존: s_center = 0.5 + 1j * t_zero  (ζ, Dirichlet)
# EC:   s_center = 1.0 + 1j * t_zero  (weight 2)
# Δ용:  s_center = 5.5 + 1j * t_zero  (weight 12 → σ = (k-1)/2 = 11/2)
critical_sigma = 5.5
```

### Λ 함수 정의 (Python mpmath)
```python
# Λ(s, Δ) = (2π)^{-s} Γ(s) L(s, Δ)
# 대칭: Λ(s) = Λ(12 - s)
# root number ε = 1

def completed_L_delta(s):
    """Completed L-function for Ramanujan Delta."""
    return (2 * mpmath.pi)**(-s) * mpmath.gamma(s) * L_delta(s)
```

### 결과 파일 형식
기존 C-408~C-410과 동일 형식. 추가 필드:
- `weight: 12`
- `level: 1`  
- `critical_line: sigma=11/2`
- `functional_eq: Lambda(s) = Lambda(12-s)`
- `root_number: +1`

---

## 다음 실험 예고 (C-411 완료 후)

**C-412 후보: GL(3) Symmetric Square L-함수 모노드로미**
- sym²(11a1): degree 3, weight 3, level 121.
- 현재 degree 커버리지: {1, 2}. degree 3은 첫 확장.
- "degree 2까지만 성립" 비판 봉쇄.
- opus 필수 (PARI `lfunsymsq` 사용).

대안:
- (b) Maass form (weight 0, nonholomorphic). PARI 지원 제한적.
- (c) Dedekind zeta (degree [K:Q]). number field 필요.

---

## 프로세스: 없음 (유휴 → C-411 실행 대기)
## 우선순위: NORMAL — 실험 진행

---

## [아카이브] 판정

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

### 증명 로드맵
1. ✅ Var(2/g²) > 0: 자명
2. ✅ Path C 수치적: T≤2000 양성 (R=1.37, PARI)
3. ⚠️ Path C 해석적: 보류
4. ✅ A_Λ–gap 보편성: degree 1-6, 13 L-함수, ρ=-0.893±0.016
5. ✅ 모노드로미 보편성: 5 L-함수, 105 TP / 174 FP, 예외 0건
6. ⬜ **weight 독립성**: Ramanujan Δ (weight 12) — C-411
7. ⬜ **degree 독립성**: GL(3) sym² — C-412 예정

### 누적 통계
- 모노드로미 L-함수: 5 (ζ + χ₅_even + χ₅_odd_complex + 11a1 + 37a1)
- 모노드로미: 105 TP / 174 FP, degree 1-2, rank 0-1, 예외 0건
- 지표 유형: 실수 + 복소
- ε: +1 + -1
- weight: {1, 2} — **12 미확인 → C-411**
- degree: {1, 2} — **3 미확인 → C-412**
- 논문: 4개 투고 준비 (A ~124p EN/~47p KO, B 32p, 3 16p, 4 12p)
