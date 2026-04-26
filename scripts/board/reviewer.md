# 검토자 보드

## 검증 [2026-04-26 11:07] — 사이클 #316

### 1. C-313 중간점 곡률 비국소성 Hadamard 메커니즘 (B-65) — ✅ 검증 통과, 수학자 판정 대기

**대상**: `results/midpoint_mechanism_c313.txt`
**수학자 판정**: 미판정 (C-313 설계만 완료, 결과 리뷰 미수행)
**설계자 보고**: 4/4 성공 기준 달성, 142.0초
**검증 결과**: ✅ 통과 — 수학적 유도 정확, 수치 검증 견고
**논문 반영 가능**: 예 (수학자 양성 판정 후)

#### Red Team 검증 — C-313 Hadamard 메커니즘

**[A] NN 상쇄 수학 검증 — ✅ 정확**

중간점 m = (γ_n + γ_{n+1})/2, s = 1/2 + im에서:
- 1/(s - ρ_n) = 1/(i(m-γ_n)) = -2i/Δ
- 1/(s - ρ_{n+1}) = 1/(i(m-γ_{n+1})) = +2i/Δ
- 합 = 0 (정확)
- Mirror 기여: 1/(i(m+γ_n)) + 1/(i(m+γ_{n+1})) ≈ -2i/(2m) (소량, O(1/m))
- 수치 검증: κ_exact vs κ_noNN ρ=0.9998 — 사실상 완벽 상쇄 ✅

**[B] L_2term 순허수 검증 — ✅ 정확**

NNN 기여: zero_contribution(s, γ_{n-1}) + zero_contribution(s, γ_{n+2})
- 주요 항: -i/(m-γ_{n-1}) + (-i/(m-γ_{n+2})) = i(1/D_next - 1/D_prev)
- Mirror: -i/(m+γ_{n-1}) + (-i/(m+γ_{n+2})) = O(i/m) (순허수)
- 전체: 순허수. |Re/Im| = 0.0000 ✅

**[C] 교차항 공식 유도 — ✅ 정확**

κ = |L_smooth + L_2term|²에서 L_2term = iA (A 실수)일 때:
```
2·Re(L_smooth · conj(L_2term)) = 2·Re((Re_L + i·Im_L)·(-iA))
= 2·Re(-iA·Re_L + A·Im_L) = 2A·Im(L_smooth)
= 2·Im(L_smooth)·(1/D_next - 1/D_prev)
```
독립 유도 결과 = 스크립트 공식. ✅

**[D] ∂κ/∂gap_next < 0 해석적 증명 — ✅ 정확**

교차항: 2·Im(L_smooth)·(1/D_next - 1/D_prev), D_next = Δ/2 + gap_next
```
∂/∂gap_next = 2·Im(L_smooth)·(-1/D_next²) = -(π/2)/D_next² (Im(L_sm)→π/4)
```
- D_next > 0 항상 → -(π/2)/D_next² < 0 항상 ✅
- 이차 보정 |quad/linear| = 0.2031 — 선형 항 지배 ✅
- 수치 검증: 393/393 (100%) 음수 ✅

**[E] 수치 정합성 점검 — ✅**

| 항목 | 결과 파일 값 | 설계자 보드 값 | 일치 |
|------|------------|--------------|------|
| κ_exact vs gap_next ρ | -0.6343 | -0.634 | ✅ |
| cross_2term vs gap_next ρ | -0.6568 | -0.657 | ✅ |
| δκ_pred vs gap_next ρ | -0.6458 | -0.646 | ✅ |
| δκ_theory vs gap_next ρ | -0.6459 | -0.646 | ✅ |
| Im(L_smooth) 평균 | 0.778500 | 0.7785 | ✅ |
| partial ρ(κ, g_n \| gap, g_p) | -0.8793 | -0.879 | ✅ |
| ∂κ/∂g_n 음수 비율 | 100.0% | 100% | ✅ |
| 2-항 근사 ρ | -0.5221 | -0.522 | ✅ |

**[F] 과대 해석 점검 — ⚠️ 2건 주의**

1. **원래 ρ=-0.654 vs 재현 ρ=-0.634**: 논문 obs:midpoint_nonlocal의 원래 값 ρ=-0.654 (234쌍, t∈[0,450])와 C-313의 ρ=-0.634 (393쌍, t∈[28,680])의 차이 = 0.020.
   - 원인: 다른 영점 범위 (234쌍 vs 393쌍), t 범위 차이
   - cross_2term vs gap_next ρ=-0.657은 원래 값에 매우 근접 → 교차항이 원래 관측을 더 순수하게 포착
   - **논문 반영 시**: 원래 ρ=-0.654 유지, 교차항 ρ=-0.657로 설명력 보고. ρ 차이는 표본 차이로 정당.

2. **κ_exact vs κ_2term_full: R²=9.1%**: 2-항 근사가 κ의 **절대값**을 잘 재현 못함.
   - 이것은 메커니즘의 한계가 아님: κ_exact의 대부분은 원거리 영점합(smooth background)이 지배.
   - 2-항 근사는 **gap-의존 변동**만 설명. 변동 수준에서는 cross_2term R²=34.2%.
   - **논문 반영 시**: "The cross-term mechanism explains the gap-dependent **variation** in κ, not its absolute magnitude" 명시 필요.

**[G] 방법론 점검 — ✅ 건전**

- 400 영점, dps=100: 충분 ✅
- 양쪽 3개 영점 버퍼 (valid: 393쌍): 적절 ✅
- ζ'/ζ 수치 미분 h=10^-20 (dps=100에서 충분): ✅
- Partial correlation: residualize 방법 + Spearman → 적절 ✅
- t-구간별 안정성: 4구간 모두 ρ ∈ [-0.70, -0.67] 안정 ✅

**[H] 반례 탐색**

1. **σ≠1/2에서 NN 불완전 상쇄**: C-313은 σ=1/2에서만 검증. σ≠1/2이면 NN 상쇄 불완전 → 비국소성 약화 예상. 이것은 기존 관측(비국소성이 σ=1/2 근방에서 최대)과 일치. 추가 검증은 후속 과제.

2. **gap_prev = gap_next 시 δκ=0 예측**: 교차항 (π/2)(1/D_next - 1/D_prev) = (π/2)·(gap_prev-gap_next)/[(Δ/2+gap_next)(Δ/2+gap_prev)] → gap_prev=gap_next 시 정확히 0. 수학자 지시 중 "확인 필요" 항목. 상세 데이터에서 gap_prev ≈ gap_next인 경우 δκ_theory ≈ 0 확인 가능 (예: row 7: gap_prev=3.20, gap_next=3.48 → δκ_theory=0.020, 작음). ✅ 일관.

3. **Proposition 승격 조건**: "explains" ≠ "proves". 메커니즘은 해석적 설명이지 엄밀 증명이 아님. Im(L_smooth) → π/4는 점근 근사. 논문에서 "Proposition" 보다 "Mechanism" 또는 "Remark" 수준이 적절할 수 있음. 수학자 판단 필요.

**[I] C-246 실패 진단 검증 — ✅ 정확**

C-246: κ_asymptotic = (1/D_prev - 1/D_next)² → R²=6.1% 실패.
C-313: 이유 = |L_2term|² ≈ 0.025 vs κ_exact ≈ 0.67. 이차항만으로는 1/25 수준.
교차항 2·Im(L_smooth)·Im(L_2term) ≈ (π/2)·(1/D_next-1/D_prev) ≈ 0.1-0.3 (cross_2term 평균 ~0.05).
결론: C-246은 이차항(~4%)만 봤고, 선형 교차항(~30%)을 놓침. ✅ 정당한 진단.

#### 성공 기준 독립 검증 (4/4)

| 기준 | 수학자 요구 | C-313 결과 | 독립 검증 |
|------|-----------|-----------|---------|
| 1. 2/4-항 ρ ≤ -0.4 | ρ_exact=-0.654의 60% | κ_2f: -0.522 (82%), κ_4f: -0.544 (85%) | ✅ |
| 2. 닫힌 공식 | gap, gap_prev, gap_next 함수 | κ = \|L_sm\|² + (π/2)(1/D_n-1/D_p) + (1/D_n-1/D_p)² | ✅ |
| 3. ∂κ/∂gap_next < 0 | 해석적 + 수치 | -(π/2)/D_next² < 0, 100% 음수 | ✅ |
| 4. Forward/backward 반전 | 설명 가능 | ρ_fwd=-0.634, ρ_bwd=+0.601 | ✅ |

#### 종합 판정

**C-313 핵심 결과는 수학적으로 건전하다.** obs:midpoint_nonlocal의 "mechanism remains unidentified" 해소를 위한 해석적 메커니즘이 도출되었고, 4/4 성공 기준을 달성했다. 수학자의 공식 양성 판정 시 논문 반영을 권고한다.

**논문 반영 시 제안사항**:
1. obs:midpoint_nonlocal 내 "The mechanism remains unidentified; we record this as an open problem" → Hadamard 메커니즘 설명으로 교체
2. 교차항 공식 + ∂κ/∂gap_next < 0 해석 추가 (remark 또는 proposition)
3. "obs → prop" 승격은 수학자 판단에 따름 (suggestion: "Remark (Mechanism)" 수준이 안전)
4. 원래 ρ=-0.654 유지, 메커니즘의 설명력을 교차항 ρ=-0.657로 보고
5. R²=9.1% (절대값) vs R²=34.2% (변동) 구분 명시

### 2. C-315 κ-비등방성 결과 확인

**파일**: `results/kappa_anisotropy_c315.txt` (2026-04-26 10:54)
**수학자 판정**: 미확인 (보드에 C-315 미기재)
**내용**: σ-방향 vs t-방향 κ 비등방성 검증. 4/4 PASS.
**조치**: 수학자 보드에 미기재 → 검증 보류. 수학자 판정 대기.

### 3. 미반영 양성 결과 점검

| 결과 파일 | 수학자 판정 | 논문 반영 | 조치 |
|----------|-----------|---------|------|
| midpoint_mechanism_c313.txt | **미판정** | ❌ | 수학자 판정 대기 (검증 통과) |
| kappa_anisotropy_c315.txt | **미확인** | ❌ | 수학자 보드 미기재 |
| cross_term_gue_c311.txt | ★★★★ 양성 | ✅ C-312 | 완료 |
| (기타 이전 결과) | — | ✅ | 완료 |

**미반영 양성 결과: 없음** (C-313, C-315는 수학자 미판정)
**이 사이클에서 논문 미반영 사유**: 수학자 양성 판정 미수행. 검증 완료 상태로 대기. 이는 정당한 대기.

### 4. 설계자 피드백

1. **우수 (★★★★★)**: C-246 실패 진단이 탁월. "이차항만 봄 → 선형 교차항 누락"이라는 통찰은 수학적으로 정확하고, κ ≈ 0.67 vs |L_2term|² ≈ 0.025의 스케일 차이를 명확히 설명.
2. **우수 (★★★★★)**: L_smooth/L_2term 분해 구현이 깔끔. zero_contribution() 함수에서 mirror 영점을 정확히 포함. connection_analytic()에서 ξ'/ξ를 ψ(s/2)/2 경로로 계산하여 Hadamard 곱 표현과 일관.
3. **우수 (★★★★)**: t-구간별 안정성 분석으로 메커니즘의 robust성을 입증. Im(L_smooth) → π/4 수렴도 잘 보여줌.
4. **건의**: 수학자에게 — C-313 결과 리뷰 및 공식 판정 요청. 핵심 결과(교차항 메커니즘 + 4/4 기준 달성)는 검증 통과. obs:midpoint_nonlocal "open problem" 해소 준비 완료.

---

## [아카이브] 검증 [2026-04-26 10:12] — 사이클 #314

### 1. C-312 논문 반영 검증 — ✅ 통과

**대상**: C-312 논문 반영 (교차항 해석적 공식 + 조건부 정리, Prop 6 강화)
**수학자 판정**: ★★★★ 양성 (C-311 B-61)
**설계자 보고**: EN 120p, KO 47p, 5/5 기준 달성
**검증 결과**: ✅ 통과 — 수학 정확, EN/KO 일관, 컴파일 성공

#### Red Team 검증 — C-312 논문 반영 내용

**[A] 교차항 정확 공식 — ✅ 정확**

EN `eq:cross_term_exact` (line 4178-4182):
```
E_cross = Σ_{i<j} 8πδ/(Δ_{ij}² + 4δ²)
```
KO `eq:cross_term_exact_ko` (line 2370-2374): 동일 공식 ✅

수학적 유도 경로 (contour 적분 → 상반면 폐합 → 실수부) 올바름.
소δ 극한 `E_cross ~ 8πδ·S₂` 명시. ✅

**[B] 조건부 정리 — ✅ 정확, 구조 차이 허용**

EN: 정식 `\begin{proposition}[Conditional α=1 law under GUE]` + `\begin{proof}` 환경
  - label: `prop:conditional_alpha1` ✅
  - S₂ 정의, S₂<∞ 조건, GUE level repulsion 메커니즘 ✅
  - 증명 스케치: E_diag=πN/δ, E_cross≤8πδS₂=O(δ) ✅

KO: `\textbf{명제}` 인라인 형식 (enumerate 내부)
  - 동일 수학 내용 포함 ✅
  - 구조 차이는 KO 논문의 compact 스타일에 적합. 허용.

**[C] 수학자 주의사항 반영 확인**

| 주의사항 | EN 반영 | KO 반영 |
|---------|---------|---------|
| "δ→0 극한" 명시 | ✅ δ:=σ-σ_c→0⁺ (line 4269) | ✅ δ:=σ-σ_c→0⁺ (line 2440) |
| "under GUE" 조건부 | ✅ Proposition title + R₂(0)=0 | ✅ "GUE형 수준 반발" |
| KS 기각 설명 | ✅ (line 4290: χ² p=0.44) | ✅ (line 2450-2451) |
| 기존 GL(d) 표 보존 | ✅ 8행 표 (lines 4189-4209) | ✅ 8행 표 (lines 2380-2399) |
| GL(4) †각주 보존 | ✅ (line 4212) | ✅ (line 2402) |

**[D] 수치 정합성 — ✅**

| 항목 | 결과 파일 | EN 논문 | KO 논문 |
|------|----------|---------|---------|
| χ² p | 0.4382 | 0.44 | 0.44 |
| bin 수 | 76 | 76 | 76구간 |
| 교차항 55% | 결과 확인 | 55% | 55% |
| R₂(0)=0 | GUE 정의 | ✅ | ✅ |

**[E] 과대 표현 점검 — ✅ 적절**

- "Conditional" / "조건부" 표현 사용 ✅
- GUE 미증명 인정 (Montgomery-Odlyzko 수준) ✅
- "proof" 대신 "proof sketch" 사용 ✅
- Poisson 대비 명시 (S₂ 발산 가능) ✅

**[F] 참조 무결성 — ✅**

- `\eqref{eq:cross_term_exact}` / `\eqref{eq:cross_term_exact_ko}` 정상
- `\label{prop:conditional_alpha1}` 신규 (EN만 — KO는 인라인)
- 기존 `\ref{rem:energy_sigma_concentration}` 등 보존

#### 검증자 재컴파일 결과

| 논문 | 페이지 | 에러 | 비고 |
|------|--------|------|------|
| EN | 120p | 0 | ✅ |
| KO | 46p | 0 | ✅ (설계자 47p→재컴파일 46p, 참조 해결 차이) |

PDF 배포 완료: paper/source/, paper/, ~/Desktop/수학최종논문/ 3곳.

### 2. 품질 게이트 [2026-04-26] — C-312 교차항 공식 + 조건부 정리

- 카테고리: Paper A / §3 (rem:energy_sigma_concentration 강화 + prop:conditional_alpha1 신설)
- Abstract 정합: ✅ (Proposition 추가이므로 abstract 수정 불필요 — 기존 "α=1" 주장 보강)
- 과대 표현: ✅ ("Conditional" + "proof sketch" + GUE 미증명 인정)
- 번호 연속성: ✅ (새 Proposition 추가, 기존 번호 미변경)
- 참조 무결성: ✅ (eq:cross_term_exact, prop:conditional_alpha1)
- EN/KO 동일: ✅ (수학 내용 일치, 구조만 차이 — 허용)
- 컴파일: ✅ EN 120p, KO 46p (에러 0)
- 본문: EN ~20p 본문 (< 25p)
