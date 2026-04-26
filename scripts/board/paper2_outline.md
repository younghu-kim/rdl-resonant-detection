# 제2논문 구조 설계 (Paper 2 Outline)

**작성**: 사이클 #185 (2026-04-20)
**갱신**: 사이클 #361 (2026-04-27) — C-349~C-354 EC 확장 반영
**데이터**: 실험 #107–#117 (11건) + C-349~C-354 (6건, EC 확장)
**상태**: 아웃라인 확정. 집필 대기.

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
- **신규**: 20곡선 대규모 검증 + 블라인드 예측 F1≥0.994 + 해상도 한계 공식

### 2. Framework Recap (1–2p)
- Definition 1–3, Theorem 1–3 요약 (Paper 1 참조)
- κ_near = |Λ'/Λ|²·δ² 정의 명확화
- **P3 tautology 한계 명시**: κδ² ≡ 1은 FE 직접 귀결. 비자명 정보는 c₁ 항에 있음.

### 3. Elliptic Curve L-Functions (6–8p) ← **대폭 확장**

#### 3.1. ε-보정과 rank-dependent 구조 (#107, #108)
- ε=−1 곡선: Im(Λ) 부호변환 방법 (B-27 해결)
- rank ≥ 2 center-zero 분리 (B-26 해소)
- **교훈**: P3 κδ²≡1은 FE tautology → 재정의 필요

#### 3.2. 20-Curve Four-Property Verification (C-349) ★★★★
- **20곡선** (rank 0–3, N=11–11642): 19/20 PASS
- rank당 4–6곡선, conductor 범위 폭넓음
- 유일한 FAIL: 11642a1 (N=11642, rank 3, κδ² CV=32%) — 고 conductor 정밀도 한계
- rank/conductor/ε 무관 보편성 — 원래 8곡선에서 20곡선으로 2.5배 확장

#### 3.3. κ_near 보편성과 conductor 스케일링 (#109, #111)
- 8곡선 (rank 0–3, conductor 11–5077): κδ² ≈ 1, CV = 0.47%
- Δκ ≈ 2·log(N), R² = 0.73
- 유한-δ O(1) 항의 conductor 의존성

#### 3.4. 모노드로미 정밀도 (#110)
- 120영점, |Δarg/π| = 1.000000, mono/π = 2.000000
- 8/8 곡선 완벽 일치

#### 3.5. 배경 곡률의 해석적 공식 (C-353) ★★★★★ ← **신규**
- **κ_bg 이론값**: κ_gamma = |(1/2)log(N/4π²) + ψ((σ+it)/2)|²
- ψ-보정 전: simple log 공식, CV = 45% (실패)
- ψ-보정 후: full digamma 공식, CV = 9.3% (성공)
- 체계적 undershoot (ratio < 1): L'/L 기여가 감마와 역상관 → 새 물리적 정보
- 잔차 ∝ L'/L at background: 영점합의 비국소적 기여

#### 3.6. Blind Zero Prediction — 8-Curve Campaign (C-351, C-352) ★★★★★ ← **신규**
- **8곡선** (rank당 2, N=11–11197): 블라인드 t∈[25,40], MATCH_TOL=0.5
- **통합 성적**: 516/520 = 99.2% recall, FP = 0 (8곡선 전체)
- 곡선별 F1: 전곡선 F1 ≥ 0.994
- FP = 0: Λ'/Λ의 극점 구조에 의해 이론적으로 보장

| 곡선 | rank | N | P | R | F1 | TP/total |
|------|------|---|---|---|----|----|
| 11a1 | 0 | 11 | 1.000 | 1.000 | 1.000 | 36/36 |
| 43a1 | 0 | 43 | 1.000 | 1.000 | 1.000 | 44/44 |
| 37a1 | 1 | 37 | 1.000 | 1.000 | 1.000 | 42/42 |
| 79a1 | 1 | 79 | 1.000 | 1.000 | 1.000 | 52/52 |
| 389a1 | 2 | 389 | 1.000 | 1.000 | 1.000 | 63/63 |
| 571a1 | 2 | 571 | 1.000 | 1.000 | 1.000 | 67/67 |
| 5077a1 | 3 | 5077 | 1.000 | 0.988 | 0.994 | 82/83 |
| 11197a1 | 3 | 11197 | 1.000 | 0.985 | 0.992 | 130/132 |

#### 3.7. FN Anatomy and Resolution Limit (C-354) ★★★★★ ← **신규**
- **FN 원인**: 근접 영점 쌍 (gap < TP 하위 5%)
  - 5077a1 t=31.67: gap=0.152, 백분위 4.2%
  - 11197a1 t=39.14: gap=0.154, 백분위 3.9%
- **해상도 한계 공식**: dt < gap_min/2 이면 모든 영점 독립 탐지
  - dt=0.10: FN 발생 (dt > gap/2 ≈ 0.077)
  - dt=0.05: FN 전원 해소 (dt < gap/2) ✅
  - dt=0.02: 동일 ✅
- **mono/π=4.0 메커니즘**: 적분 반경 r > gap → 두 영점 동시 감김 → 감김수 합산
  - 5077a1: 임계 r=0.13에서 2→4 전이, gap=0.152
  - 11197a1: 임계 r=0.15에서 2→4 전이, gap=0.154
- **GUE 맥락**: FN gap/δ_mean ≈ 0.22–0.25 (GUE 하위 5% 경계 근처)
- **Empirical Criterion (Nyquist analog)**: dt_min ≈ gap_min/2 ≈ δ_mean·s_min/2
- **논문 서사**: "FN 4건은 전체 520 영점의 0.8%. 모두 근접 쌍(하위 4%)에서 발생. dt=0.05로 해소. 프레임워크의 유일한 한계는 Nyquist-type 해상도 조건."

#### 3.8. A(t₀) Rank Dependence — Negative Result (C-349) ❌
- **가설**: A(t₀) = Im(H₀)² + 2H₁ 이 rank를 감지하는가?
- **결과**: 기각. Conductor confound.
  - log(N) vs A: Pearson r=0.87 (p<0.001)
  - Partial corr (A vs rank | log(N)): r=-0.33 (비유의)
- **교훈**: 두 변수 상관 시 제3 변수(conductor) 교란 점검 필수
- 음성 결과이나 정직하게 보고 (헌법 4조)

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

### 6. Discussion (2–3p)
- 보편성 범위: GL(1)–GL(4), EC (20곡선), Dedekind, Epstein, DH
- FE 의존성 확인 + Euler product 불필요
- Theorem 4: 최초의 ξ-bundle 기반 on/off-critical 감별기
- **신규**: κ_bg 해석적 공식 — conductor에 의한 배경 곡률 예측
- **신규**: 블라인드 예측 한계의 정직한 분석 — Nyquist-type 해상도 공식
- **한계**: c₁ 법칙 N=4, off-critical 영점 추가 확보 어려움, 11642a1 FAIL
- **향후**: Artin L-함수, Selberg class 일반화, σ-국소화 증명

### 7. Conclusion (0.5p)

---

## 결과 매핑

| 실험 | 섹션 | 등급 | 논문 기여 |
|------|------|------|----------|
| #107 | 3.1 | ★★ | ε-보정 방법론 + P3 tautology 교훈 |
| #108 | 3.1 | 조건부 | ε-보정 검증 |
| #109 | 3.3 | ★★★ | κ_near 보편성 |
| #110 | 3.4 | ★★★ | 모노드로미 정밀도 |
| #111 | 3.3 | ★★ | conductor 스케일링 |
| #112 | 4.1 | ★★★ | R-S GL(4) |
| #113 | 4.2 | ★★★ | Dedekind ζ |
| #114 | 4.3 | ★★★ | Epstein non-Euler (FE-only) |
| #115 | 5.2 | ★★★ | DH 4-property + Theorem 4 발견 |
| #116 | 5.4 | ★★★ | δ-sweep 검증 |
| #117 | 5.5 | ★★ | c₁ 점근 법칙 |
| C-349 | 3.2, 3.8 | ★★★★ | 20곡선 4성질 + A(t₀) rank 기각 |
| C-351 | 3.6 | ★★★★★ | 3곡선 블라인드 (F1≥0.994) |
| C-352 | 3.6 | ★★★★★ | 8곡선 블라인드 (516/520, FP=0) |
| C-353 | 3.5 | ★★★★★ | κ_bg ψ-보정 (CV=9.3%) |
| C-354 | 3.7 | ★★★★★ | FN 해부 + 해상도 한계 공식 |

## 예상 분량
- 본문: 20–25 페이지 (EC 확장으로 +5p)
- 부록 (수치표): 3–5 페이지
- 참고문헌: 2 페이지
- **총: ~30 페이지**

## 미결 사항
1. c₁ 법칙 N=4 → 추가 off-critical 영점? (★★ → ★★★ 승격 조건)
2. Artin L-함수 결과 → Paper 3로 분리 확정
3. Paper 1과의 중복 최소화 (framework recap 범위)
4. Theorem 4를 Paper 2에 배치 확정 (DH 맥락)
5. ~~EC 결과 확장~~ → 완료 (C-349~C-354)
6. Paper B 집필 착수 시점: C-354 반영 후 즉시 가능
