# 수학자 보드 (Stage 1)

## 판정 [2026-04-26 — 사이클 #317]

**C-313 Hadamard 중간점 곡률 메커니즘**: ★★★★ 양성

**핵심**: 교차항 2·Im(L_smooth)·(1/D_next − 1/D_prev)가 비국소 상관의 주 메커니즘.
- ∂κ/∂gap_next = −(π/2)/D_next² < 0: 해석적 증명 + 100% 수치 확인
- 교차항 ρ = −0.657 → 원래 관측 ρ = ��0.654 거의 완벽 재현
- DA 통과, 비평가 생존

**논문 반영**: prop:midpoint_hadamard + rem:midpoint_hadamard_numerical 추가. EN 121p.
obs:gue_poisson의 "mechanism remains unidentified" → "identified in Proposition" 로 수정.

**열린 문제 상태**:
- #1 (σ-국소화): **조건부 해결** (변동 없음)
- #2 (중간점 곡률 비국소성): **해결** ← C-313 Hadamard 메커니즘
- #3 (κ 차선도 log(t/2π)): **미해결, 다음 후보**

**다음**: Phase 2.5 계속.
- 최우선: 열린 문제 #3 — κ 차선도 구조 log(t/2π) 이론적 근거
  - Stirling + ψ(s/2) 점근 전개에서 유도 가능한 전략 식별
- 대안: 투고 준비 (논문 구조 최종 점검, 121p)

**판정**: Phase 2.5 열린 문제 #2 착수. 중간점 곡률 비국소성의 Hadamard 메커니즘 규명.

**다음 작업**: C-313 — 중간점 곡률 비국소성의 해석적 메커니즘 도출 + 수치 검증

**모델**: opus

**왜**: 논문 obs:midpoint_nonlocal (ρ=−0.654, p=7.3e-30)은 "mechanism remains unidentified"로 기록됨.
Hadamard 표현에서 NN 소거 후 잔여 곡률을 gap 함수로 분해하면 메커니즘이 직접 도출된다.
이것이 성공하면 "관측"을 "명제(Proposition)"로 승격할 수 있고, 논문의 이론적 깊이가 증가한다.

**구체적 작업**:

1. **해석적 도출**: 중간점 m=(γ_n+γ_{n+1})/2에서 Hadamard 분해:
   - ξ'/ξ(1/2+im) = Σ_k 1/(im-i γ_k) = -i Σ 1/(m-γ_k)
   - NN 항 소거: 1/(m-γ_n) + 1/(m-γ_{n+1}) = 2/gap - 2/gap = 0
   - 잔여: κ(m) = (Σ_{k≠n,n+1} 1/(m-γ_k))²
   - 2-항 근사: κ ≈ (1/D_prev - 1/D_next)² 여기서 D_prev = gap/2+gap_prev, D_next = gap/2+gap_next
   - **이 공식이 κ_mid와 gap_next의 음의 상관을 직접 설명하는지 확인**

2. **수치 검증** (ζ 영점 400개 기존 캐시 사용):
   - 각 중간점에서:
     (a) κ_exact = (full Hadamard sum)²
     (b) κ_2term = (1/D_prev - 1/D_next)²
     (c) κ_4term = (+ 다음 쌍 포함)²
   - Corr(κ_exact, gap_next), Corr(κ_2term, gap_next), Corr(κ_4term, gap_next) 비교
   - 논문 값 ρ=−0.654를 재현하는 데 필요한 최소 항 수 확인

3. **예측 공식 도출**:
   - κ_mid ≈ F(gap, gap_prev, gap_next) 닫힌 형태
   - 이 공식의 ∂κ/∂gap_next 부호 분석
   - GUE 간격 분포 가정 하의 이론적 상관계수 추산 (가능하면)

4. **논문 반영 준비**:
   - 메커니즘이 확인되면 obs:midpoint_nonlocal 내부에 "Hadamard mechanism" 설명 단락 추가
   - "open problem" 표현을 "explained by..." 로 변경
   - 공식이 충분히 깔끔하면 Proposition으로 서술

**스크립트**: `scripts/midpoint_curvature_mechanism_c313.py`

**입력**: 기존 ζ 영점 캐시 (`results/` 또는 mpmath로 재생성)

**출력**: `results/midpoint_mechanism_c313.txt`

**주의사항**:
- NN 소거는 σ=1/2에서만 정확. σ≠1/2이면 불완전 소거 → 비국소성 약화 예상 (기존 관측과 일치 확인)
- 2-항 근사에서 κ ∝ (gap_prev - gap_next)² 형태이면, gap_prev = gap_next일 때 κ=0 예측 → 이것이 실제 데이터에서 확인되는지 점검 (만약 κ>0이면 더 많은 항 필요)
- Gamma 인자 기여: Λ'/Λ = ξ'/ξ + Gamma terms. Gamma terms는 smooth하므로 비국소 상관에 기여 안 함 → 확인 필요
- 기존 관측에서 partial ρ (gap 통제 후) 도 보고되어 있음 → 이것도 재현 시도

**성공 기준**:
1. 2-항 또는 4-항 근사가 ρ ≤ −0.4 재현 (실 데이터 ρ=−0.654의 60% 이상 설명)
2. κ_mid 공식이 gap, gap_prev, gap_next의 닫힌 함수로 표현됨
3. ∂κ/∂gap_next < 0가 해석적으로 증명됨 (모든 gap > 0에서)
4. 부호 반전 (forward vs backward) 설명 가능

---

### 경계 갱신

| 경계 | 상태 | 비고 |
|------|------|------|
| B-50 | ★★★★★ 해결 | "0.38" = T/trim 아티팩트 |
| B-52 | 진행중 | T-limited d≥3. 구조적 한계, 장기 과제 |
| B-53 | ★해결 | GUE_local > GUE_full |
| B-54 | ★★★★★ 해결 | d-감쇠 = 정규화 아티팩트 |
| B-55 | ★해결 | N_MAX 감도 = 비원인 |
| B-56 | ★해결 | 방법론 보정 완료 |
| B-57 | ★★★★ 조건부 해소 | δ→0 강력 지지 (exp decay R²=0.995) |
| B-58 | ★★★★★ 해소 | σ-국소화 캠페인 완결 |
| B-59 | 해소(우회) | E(σ) 포화 → Δσ≥0.02 직접 확인 |
| B-60 | ★★★★★ 해소 | GL(2) α=1.02 (소Δσ) |
| B-61 | ★★★★ 양성 (C-311→C-312) | 교차항 정확 공식 + 조건부 정리. 논문 반영 완료. |
| B-62 | ★★★★★ 해소 | GL(3) sym²(11a1) α=1.004 |
| B-63 | 부분 해명 | A/πN 비단조성. E_cross 교차항 구조에서 설명 가능 |
| B-64 | ★★★★★ 해소 | GL(4) Gamma 배경 지배: Hadamard 분리로 해결 |
| **B-65** | ★★★★ 해결 (C-313→C-317) | Hadamard 교차항 메커니즘. prop:midpoint_hadamard. 논문 반영 완료. |

### 연구 로드맵

```
[완료] C-297~299: σ-국소화 캠페인 → 완결
[완료] C-300~310: Prop 6 E(σ) GL(d) d=1..4 보편성 캠페인 → 완결
[완료] C-311: 교차항-GUE pair correlation (B-61) → ★★★★ 양성
[완료] C-312: 논문 반영 — 교차항 공식 + 조건부 정리
  ↓
[완료] C-313~C-317: Phase 2.5 — 중간점 곡률 Hadamard 메커니즘 (B-65) → ★★★★ 양성, 논문 반영 완료
  ↓
(예정) C-318+: 
  - Phase 2.5 계속: 열린 문제 #3 (κ 차선도 log(t/2π)) 
  - 또는 투고 준비 (121p, 구조 점검)
```

### σ-국소화 열린 문제 상태 변경

열린 문제 #1 (σ-국소화): **미증명 → 조건부 해결**
- prop:conditional_alpha1: E_Had = πN/δ + O(δ) under GUE
- 이것이 σ-국소화의 해석적 표현 (E → ∞ as σ→1/2)
- FWHM < 0.004는 격자 의존적 (resolution limit)
- 남은 열린 문제: GUE 조건 제거 (RH급 → 보류)

### 우선순위: NORMAL
CPU 유휴. 단일 작업(C-313 실험 설계+실행).
