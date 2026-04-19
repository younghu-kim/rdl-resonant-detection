# 제2논문 구조 설계 (Paper 2 Outline)

**작성**: 사이클 #185 (2026-04-20)
**데이터**: 실험 #107–#117 (11건)
**상태**: 초안 설계

---

## 가제목

**EN**: "ξ-Bundle Curvature as a Critical-Line Discriminator: Universality Across L-Function Families and an Analytic Detection Theorem"

**KO**: "ξ-다발 곡률의 임계선 감별 기능: L-함수 패밀리 보편성과 해석적 검출 정리"

---

## 논문 구조

### 1. Introduction (2–3p)
- Paper 1 요약: ξ-bundle framework (Thm 1–3, Gauss-Bonnet), ζ(s) + GL(1)–GL(5) 검증
- Paper 2 동기: (a) 자기쌍대 타원곡선, (b) Euler product 없는 L-함수, (c) **RH 위반 영점 감별**
- 핵심 기여 미리보기: Theorem 4 (κδ² discrimination) + 보편성

### 2. Framework Recap (1–2p)
- Definition 1–3, Theorem 1–3 요약 (Paper 1 참조)
- κ_near = |Λ'/Λ|²·δ² 정의 명확화
- **P3 tautology 한계 명시**: κδ² ≡ 1은 FE 직접 귀결. 비자명 정보는 c₁ 항에 있음.

### 3. Elliptic Curve L-Functions (3–4p)

#### 3.1. ε-보정과 rank-dependent 구조 (#107, #108)
- ε=−1 곡선: Im(Λ) 부호변환 방법 (B-27 해결)
- rank ≥ 2 center-zero 분리 (B-26 해소)
- **교훈**: P3 κδ²≡1은 FE tautology → 재정의 필요

#### 3.2. κ_near 보편성 (#109)
- 8곡선 (rank 0–3, conductor 11–5077): κδ² ≈ 1, CV = 0.47%
- rank/conductor/ε 무관 보편성

#### 3.3. 모노드로미 정밀도 (#110)
- 120영점, |Δarg/π| = 1.000000, mono/π = 2.000000
- 8/8 곡선 완벽 일치

#### 3.4. Conductor 스케일링 (#111)
- Δκ ≈ 2·log(N), R² = 0.73
- 유한-δ O(1) 항의 conductor 의존성

### 4. Beyond Euler Product (2–3p)

#### 4.1. Rankin-Selberg GL(4) (#112)
- 2쌍, 4/4 PASS. 다른 GL(4) 구성과의 비교 (sym³ vs R-S).

#### 4.2. Dedekind ζ (#113)
- 7수체 (D=−3~−163, h=1,3): 4/4 PASS
- 주의: Dedekind = ζ·L(χ_D), Euler product 있음

#### 4.3. Epstein ζ — FE-Only Sufficiency (#114)
- **핵심 결과**: 3 non-Euler Epstein ζ (h > 1), 47영점, 4/4 PASS
- Euler product 없이도 작동 → **ξ-bundle = FE 구조에만 의존**
- Euler/non-Euler 영점 비일치 (0/5) → 독립 확인

### 5. Critical-Line Discrimination: Theorem 4 (3–4p)

#### 5.1. Davenport-Heilbronn Function
- 정의, 함수방정식 + non-Euler
- 알려진 off-critical 영점: t<200에서 4개 (σ ∈ [0.57, 0.81])

#### 5.2. ξ-Bundle 4-Property on DH (#115)
- On-critical 10/10: κδ²=1.002±0.001, mono/π=2.000
- Off-critical 4/4: κδ²=1.210±0.099, mono=±π (3/4 정상)
- min(off) > max(on) → **완전 분리**

#### 5.3. Theorem 4: κδ² Critical-Line Discrimination
**Statement** + **Proof** (2줄: FE 미분 → cos(β−α)=0)

#### 5.4. δ-Sweep Validation (#116)
- 8 δ-값, log|κδ²−1| vs log(δ)
- On-critical: 기울기 2.00±0.003 (이론 2.0), R²=1.0000
- Off-critical: 기울기 0.96±0.04 (이론 1.0), R²≥0.999
- 모든 δ에서 완전 분리

#### 5.5. Corollary 4.1: c₁ Asymptotics (#117)
- c₁ = Re(Λ''/Λ') = 1/(σ₀−1/2) + O(1)
- 증명 (거울쌍 Hadamard 분해)
- 수치: c₁_fit/c₁_analytic = 0.998±0.002, N=4
- c₁·|σ−½| → 1 (mean=1.08)

### 6. Discussion (2p)
- 보편성 범위: GL(1)–GL(4), EC, Dedekind, Epstein, DH
- FE 의존성 확인 + Euler product 불필요
- Theorem 4: 최초의 ξ-bundle 기반 on/off-critical 감별기
- **한계**: c₁ 법칙 N=4, off-critical 영점 추가 확보 어려움
- **향후**: Artin L-함수, Selberg class 일반화, σ-국소화 증명

### 7. Conclusion (0.5p)

---

## 결과 매핑

| 실험 | 섹션 | 등급 | 논문 기여 |
|------|------|------|----------|
| #107 | 3.1 | ★★ | ε-보정 방법론 + P3 tautology 교훈 |
| #108 | 3.1 | 조건부 | ε-보정 검증 |
| #109 | 3.2 | ★★★ | κ_near 보편성 |
| #110 | 3.3 | ★★★ | 모노드로미 정밀도 |
| #111 | 3.4 | ★★ | conductor 스케일링 |
| #112 | 4.1 | ★★★ | R-S GL(4) |
| #113 | 4.2 | ★★★ | Dedekind ζ |
| #114 | 4.3 | ★★★ | Epstein non-Euler (FE-only) |
| #115 | 5.2 | ★★★ | DH 4-property + Theorem 4 발견 |
| #116 | 5.4 | ★★★ | δ-sweep 검증 |
| #117 | 5.5 | ★★ | c₁ 점근 법칙 |

## 예상 분량
- 본문: 15–20 페이지
- 부록 (수치표): 3–5 페이지
- 참고문헌: 2 페이지
- **총: ~25 페이지**

## 미결 사항
1. c₁ 법칙 N=4 → 추가 off-critical 영점? (★★ → ★★★ 승격 조건)
2. Artin L-함수 결과 포함 여부 (현재 데이터 없음)
3. Paper 1과의 중복 최소화 (framework recap 범위)
4. Theorem 4를 Paper 1에 포함할지 Paper 2에 둘지 → **Paper 2** (DH 맥락)
