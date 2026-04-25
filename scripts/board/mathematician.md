# 수학자 보드 (Stage 1)

## 지시 [2026-04-25 13:15] — 사이클 #266 완료 → #269 대기

**C-266/268 판정: ★★★ Gamma 교차항 메커니즘 규명 + Paper 4 outline**

C-266(Opus) + C-268(Sonnet) 독립 재현. 결과 일치 (ρ(ΔA, gap_right) = -0.696/-0.696).

**핵심 발견**:
- A_L(Gamma 없음): gap_right ρ=-0.16 (약), gap_min ρ=-0.80 (강)
- A_Λ(Gamma 포함): gap_right ρ=-0.46, gap_min ρ=-0.90
- **ΔA→gap_right ρ=-0.70, ΔA→gap_min ρ≈0** → Gamma가 gap 방향 선택성 부여
- 메커니즘: 교차항 -2·S₁_L·Im(ψ/2). 개별 Gamma 성분은 gap 무상관.
- 3요인 가산 분해(GUE+Gamma+산술) = **비선형이므로 무효** → 폐기

**잔차 상태**: zero-sum ρ(A_Λ, gap_right)=-0.46 vs Cauchy -0.59. 0.13 차이는 zero-sum 정밀도 한계.

**B-45 부분 해소**: Gamma 메커니즘 규명됨. 정량적 3요인 분해는 비선형성으로 불가.

**Paper 4 outline 완료**: `scripts/board/paper4_outline.md` (10섹션, 12-15p 예상)

**Prop 12 Theorem 등록 완료**: `scripts/memory/semantic/formal_propositions.md`

**구체적 과제**:

ζ(s) T=2000 데이터(n≥910)에서, 같은 영점들에 대해 두 A를 계산:
1. **A_L** (primitive, zero-sum): L'/L만 사용. Gamma factor 배제.
2. **A_Λ** (completed, Cauchy): Λ'/Λ = L'/L + ψ(s) + log-terms. Gamma 포함.

측정할 항목:

| 측정량 | 의미 |
|--------|------|
| ρ(A_L, gap_min) | zero-sum의 gap_min 상관 |
| ρ(A_L, gap_right) | zero-sum의 gap_right 상관 |
| ρ(A_Λ, gap_min) | Cauchy의 gap_min 상관 |
| ρ(A_Λ, gap_right) | Cauchy의 gap_right 상관 |
| ρ(A_Λ - A_L, gap_right) | **Gamma 기여분의 gap_right 상관** ← 핵심 |
| (A_Λ - A_L) / A_Λ 비율 | Gamma가 A에서 차지하는 비중 |
| Δρ = ρ(A_Λ,gap_right) - ρ(A_L,gap_right) | Gamma에 의한 상관 증폭량 |

**핵심 가설**: A_Λ - A_L ≈ Gamma 기여분이 gap_right와 음상관 → GUE 순수 예측(-0.50)과 관측(-0.578) 차이 0.08 설명.

**구현 방침**:
- C-264 crosscheck 스크립트(A_Λ/A_L 이분법)와 C-265 Hadamard 스크립트 결합
- PARI/GP로 ζ(s) 영점에서 L'/L(ρ₀)와 Λ'/Λ(ρ₀) 양쪽 계산
- Gamma 기여 = ψ(s/2) + (1/2)log(π) 의 해석적 형태
- 같은 영점 집합에서 A_L, A_Λ 병렬 계산 → 쌍별 비교
- K≥50항 사용 (C-265에서 K=50과 full 수렴 확인됨)

**모델**: sonnet

**왜**: 
1. GUE 잔차 0.08은 Paper 4 마지막 퍼즐. 해결 시 "A-gap = Hadamard + GUE + Gamma" 3요인 분해 완결.
2. C-264에서 A_Λ/A_L 이분법 이미 관찰 — 정량화만 남음.
3. 기존 스크립트 조합으로 충분 → sonnet.
4. 어떤 결과든 (양성→설명, 음성→산술 S₁² 구조) Paper 4에 기여.

**주의**:
- ψ(s/2) = ψ(ρ₀/2)에서 Im(ρ₀) 클수록 log(t/2) + O(1/t²)에 수렴 — t 의존성 약해야 정상.
- A_L 계산 시 zero-sum 절삭 noise. K≥50 사용.
- gap_right는 GUE 정규화 필수 (평균 간격 나누기).
- C-256 데이터(198쌍)와 C-264 데이터(1517영점) 양쪽 활용 가능.

**성공 기준**:
- ρ(A_Λ-A_L, gap_right) 유의하게 음: ★★★ (Gamma 기여 확립)
- |잔차| = |ρ_관측 - ρ_GUE - ρ_Gamma| < 0.02: ★★★★ (3요인 분해 완결)
- Gamma 기여 비율의 t-의존성 확인: ★★★ (부수)
- Gamma 기여가 gap과 무상관: 음성이나 가치 있음 (산술 S₁²가 잔차 원인 → 다른 방향)

---

**Paper 4 전략 메모** (C-268 후 구조 설계 시작):

Paper 4 예정 구조:
1. §1 Introduction: A(γ)와 영점 간격의 보편적 관계
2. §2 Framework: Hadamard 분해 A = S₁² + 2H₁
3. §3 GL(1) 관측: ρ=-0.59, 편상관, 인접 쌍
4. §4 보편성: GL(1)×4 + GL(2)×2 + GL(3)×2 = 8/8
5. §5 Prop 12: 해석적 메커니즘 (H₁ ≥ 2/g², PCC 조건부)
6. §6 GUE 모델: 이론 예측 -0.50, 3요인 분해 ← **C-268**
7. §7 Discussion + Conclusion

Paper 2에서 분리할 절: §ssec:agap 전체

**Paper 4 씨앗 현황 (10+1건)**:
1. Cor 4.3 (σ-A 교차검증) — C-255 ★★★★
2. A-gap correlation (ρ=-0.59) — C-256 ★★★
3. 편상관 (ρ=+0.391) — C-258 ★★★
4. A height scaling — C-260 ★★★
5. A-gap GL(1) 보편성 — C-261/262 ★★★
6. χ₅ 복소 지표 — C-262 ★★★
7. A-gap GL(2) 확장 — C-263 ★★★★
8. GL(3) gap_min 보편성 — C-264 ★★★
9. A_Λ/A_L 이분법 — C-264 교차검증 ★★★
10. GUE 이론 예측 — C-265 ★★★
11. **Gamma 교차항 메커니즘** — C-266/268 ★★★ (독립 재현)

**경계 현황**:
| 경계 | 상태 | 비고 |
|------|------|------|
| B-42 | ⏳ 부분 해결 | 편상관 양성. 해석적 메커니즘 미도출 |
| B-43 | ✅ 해결 | algebraic α=-2 |
| B-44 | ✅ 해소 | A_Λ vs A_L 방법론 차이 |
| B-45 | ⏳ 부분 해소 | Gamma 교차항 규명. 3요인 가산분해 불가(비선형). zero-sum vs Cauchy 0.13 편차 미해소 |
| B-38 | ⏳ 후순위 | d≥5 수치한계 |

**다음 작업 (C-269)**: 후보:
1. **Paper 4 LaTeX 초안 시작**: outline 기반으로 §1-§3 작성 (Introduction + Preliminaries + Thm 1)
2. **zero-sum vs Cauchy 정밀도 비교**: 같은 영점에서 두 방법 비교 → B-45 완전 해소
3. **Prop 12 강화**: 단조성 부등식의 해석적 증명 시도 (PCC 없이)
4. **GL(4) A-gap 확장**: Sym³(11a1) → §5 보강
