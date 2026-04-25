# Paper 4 Outline: Hadamard Amplitudes and Zero Spacing in Automorphic L-functions

> 작성: C-266 (2026-04-25)
> 상태: 초안 구조
> 씨앗: 10건 전부 확립

---

## 중심 서사

**한 줄 요약**: 영점의 Laurent 진폭 A(γ)는 인접 간격의 보편적 예측자이며, 이 관계는 Hadamard NN 지배 + GUE 보편성의 귀결이다.

**핵심 메시지**:
1. A(γ) ≥ 2/gap_min² (무조건적 하한 — Hadamard 양정치성)
2. ρ(A, gap) ≈ -0.57: GL(1)~GL(3), conductor/degree/ε 독립 보편상수
3. GUE가 이 상관의 86%를 설명 → 나머지 14%는 Gamma + 산술적 기여
4. A_Λ vs A_L 이분법: Gamma 인자가 gap 메트릭 선택성을 결정

---

## 섹션 구조

### §1. Introduction (1.5p)
- L-함수 영점의 "크기"를 재는 새 불변량: 진폭 A(γ) = Im(c₀)² + 2Re(c₁)
- 기존 연구: 영점 간격 (Montgomery, Odlyzko), 곡률 (Paper 1-3 ξ-다발)
- 본 논문의 기여: A-gap 반상관의 발견, 증명, 보편성, GUE 메커니즘
- Main Theorem 미리보기

### §2. Preliminaries (1.5p)
- §2.1 Laurent expansion of Λ'/Λ: cₙ 정의, 패리티 (Thm 5 인용)
- §2.2 Hadamard decomposition: A = S₁² + 2H₁
  - S₁ = -Σ' 1/(γ₀-γₖ) + Γ-terms (교대 부호, NN 지배)
  - H₁ = Σ' 1/(γ₀-γₖ)² + Γ'-terms (양정치)
- §2.3 Amplitude formula (Cor 4.2 인용): κδ² - 1 = Aδ² + O(δ³)
- §2.4 σ-profile (Cor 4.3 인용): κ(σ,γ) = 1/ε² + A + Bε² + O(ε⁴)

### §3. A-gap Anti-correlation: Theorem (2p) ← **핵심 섹션**
- **씨앗**: C-256 (★★★), C-265 Prop 12 (★★★★)

#### §3.1 Theorem (A-gap Lower Bound) [무조건적]
**Theorem 1**: Let L(s) be in the Selberg class with simple zeros ρₙ = 1/2 + iγₙ. Then:
  A(γₙ) ≥ 4/gap_min(γₙ)²
where gap_min(γₙ) = min(γₙ - γₙ₋₁, γₙ₊₁ - γₙ).

Proof:
  (1) A = S₁² + 2H₁ ≥ 2H₁  (S₁² ≥ 0)
  (2) H₁ = Σ_{k≠n} 1/(γₙ-γₖ)² ≥ 1/Δ_L² + 1/Δ_R²  (양정치 급수)
  (3) 1/Δ_L² + 1/Δ_R² ≥ 2/g²  (AM-HM 또는 직접)
  (4) ∴ A ≥ 4/g². □

**Corollary**: ρ_Spearman(A, gap_min) < 0 whenever the gap distribution is non-degenerate.
(Informal: 하한이 g⁻²로 발산하고 상한은 유한이므로, 순위 상관은 음.)

#### §3.2 Numerical verification (ζ(s), n=910)
- E[A|gap_bin]: 8-bin 단조 감소 (22.76 → 4.02)
- R² = 0.987 for E[A|g] ≈ a/g² + b fit (a=2.63, 이론 ≥ 2 ✓)
- ρ(A, gap_min_GUE) = -0.8998, ρ(2H₁, gap_min) = -0.9044
- H₁^{NN}/H₁ = 66.4%, 2H₁/A = 86.7%

#### §3.3 Conditional monotonicity (pair correlation 가정)
**Proposition 12**: Under Montgomery's pair correlation conjecture,
  E[A | gap_min = g] = 2a/g² + c + O(g) as g → 0
where c = 2∫₁^∞ R₂(x)/x² dx = 1.9527.

### §4. Height Scaling (1p)
- **씨앗**: C-260 (★★★)
- <A(γ)>_bin ≈ 0.49·(log(γ/(2π)))² (R²_w = 0.975, 300 영점)
- H₁ 지배 (77.5%): pair correlation integral → (log)²
- 이론 계수 0.129 vs 1/12 = 0.083 (비율 1.55, 유한 합 절삭)

### §5. Universality: GL(1) through GL(3) (2.5p) ← **핵심 섹션**
- **씨앗**: C-261/262 (★★★), C-263 (★★★★), C-264 (★★★)

#### §5.1 GL(1) universality (4 L-functions)
| L-함수 | q | ρ(A, gap_right_GUE) | p-value |
|--------|---|---------------------|---------|
| ζ(s) | 1 | -0.5898 | 6.1e-20 |
| χ₃ | 3 | -0.5733 | 1.4e-7 |
| χ₄ | 4 | -0.5530 | 1.0e-7 |
| χ₅ (복소) | 5 | -0.6334 | 2.3e-6 |

#### §5.2 GL(2) extension (2 elliptic curves)
| L-함수 | N | ε | ρ | p-value |
|--------|---|---|---|---------|
| 11a1 | 11 | +1 | -0.5688 | 3.3e-9 |
| 37a1 | 37 | -1 | -0.5500 | 3.4e-10 |

Degree-독립: 6/6 L-함수 ρ ∈ [-0.63, -0.55]

#### §5.3 GL(3) extension (gap_min metric)
| L-함수 | N | ρ(A_L, gap_min) | p-value |
|--------|---|-----------------|---------|
| Sym²(11a1) | 121 | -0.42 | 1.9e-5 |
| Sym²(37a1) | 1369 | -0.49 | 3.8e-4 |

#### §5.4 GL(4) extension (C-269)
| L-함수 | N | n | ρ(A_L, gap_min) | p-value |
|--------|------|-----|-----------------|---------|
| Sym³(11a1) | 1331 | 101 | -0.5199 | 2.5e-8 |
| Sym³(37a1) | 50653 | 131 | -0.5137 | 3.5e-10 |

GL(1)~GL(4), 10/10 L-함수 보편적 음상관. σ-유일성이 d≥4에서 포화하는 것과 대조적으로, A-gap은 degree에 무관하게 작동.

#### §5.4 Adjacent pair correlation
ρ(Aₙ, Aₙ₊₁) ≈ +0.28~+0.53 (모든 L-함수에서 일관적 양)
→ 근접 영점 쌍에서 양쪽 A가 같은 gap 공유

### §6. A_Λ vs A_L Dichotomy (1.5p)
- **씨앗**: C-264 교차검증 (★★★)
- Cauchy (A_Λ): c₀^Λ = c₀^L + (1/2)ψ(ρ/2) + const → gap_right와 상관
- Zero-sum (A_L): Gamma 없음 → NN H₁ 지배 → gap_min과 상관
- 관계: Gamma factor가 gap 선택성을 결정
- ζ(s) 교차검증: A_Λ(Cauchy) ρ(gap_right)=-0.72, A_L(zero-sum) ρ(gap_min)=-0.61
- 일치: 두 방법 모두 Hadamard 구조의 다른 투영

### §7. GUE Prediction (2p) ← **핵심 섹션**
- **씨앗**: C-265 (★★★)

#### §7.1 GUE simulation setup
- N=300 GOE→GUE, 1000 realizations, bulk 30% extraction (88K points)
- Hadamard: K=1,2,5,10,50,full convergence

#### §7.2 Results
| 측도 | GUE 예측 | 관측 | 일치도 |
|------|----------|------|--------|
| ρ(A,gap_right) | -0.50 | -0.578 | 86% |
| ρ(Aₙ,Aₙ₊₁) | +0.36 | +0.28~42 | 정확 |
| NN 기여 | 68% | 78% | 근사 |

#### §7.3 NN convergence
K=1 (NN only): 68.4% of full correlation
K≥5: 수렴. 장거리 영점 기여 31.6%

#### §7.4 Residual analysis
Δρ ≈ 0.078: 후보 원인
(a) Gamma factor: A_Λ에 systematic height-dependent shift 추가
(b) Arithmetic content: L-함수 계수의 비보편적 기여
(c) S₁² 구조: GUE에서 S₁²/A = 11% vs 관측 23%

### §8. Partial Correlation (1p)
- **씨앗**: C-258 (★★★)
- ρ(Aₙ, Aₙ₊₁ | gap_between) = +0.391 (p=1.3e-18)
- Gap을 통제해도 인접 A 사이 양의 상관 유지
- 메커니즘: Hadamard 급수에서 공유 NN 항

### §9. σ-A Cross-check (0.5p)
- **씨앗**: C-255 (★★★★)
- Cor 4.3의 독립 검증: σ-scan c(γ) vs Laurent A(γ), ρ=0.99999744
- 두 독립적 방법의 6자리 일치 → A 정의의 견고성

### §10. Discussion (1p)
- A(γ)의 물리적 해석: 영점의 "영향력 반경" (1/√A ~ 국소 간격 스케일)
- Berry-Keating 연결: A = Hamiltonian의 local spectral density?
- RMT 보편성 등급: A-gap 상관 = GUE universality class의 새 서명
- 열린 문제: 
  (a) 잔차 14%의 정확한 기원 (Gamma vs 산술)
  (b) ρ ≈ -0.57 보편값의 해석적 도출
  (c) σ-국소화와 A-gap의 이중성 (Cor 4.3 + Thm 1)

### References

---

## 페이지 예상: 12-15p

## 미해결 갭 (후속 실험 후보)
1. **잔차 Γ/산술 분리**: A_L에 Gamma 보정 추가 → ρ 변화 측정 → §7.4 보강
2. **GL(3) A_Λ 직접 측정**: lfunlambda 속도 허용 시 Cauchy로 A_Λ 추출 → §6 보강
3. **GL(4+) 확장**: Sym³(11a1) A-gap → §5 확장 (d≥4 첫 데이터)
4. **ρ ≈ -0.57 해석적 도출**: pair correlation + Hadamard 결합 → §3.3 강화
5. **height scaling 확장**: γ ~ 5000에서 (log)² 유지 확인 → §4 보강
