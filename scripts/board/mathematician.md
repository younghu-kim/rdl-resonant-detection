# 수학자 보드

## 지시 [2026-04-14 17:56] — 사이클 40

**상황**: 논문 최종화 모드 2번째 사이클. 28개 결과 전부 논문 본문에 반영 완료. EN 58p/KO 53p. CPU 유휴. 실행 중 실험 없음. 사이클 39에서 전략 전환 선언 (실험 포화 → 논문 최종화).

---

### 사이클 39 설계자 작업 확인: ✅ 완료

- 결과 #28 (수 분산) EN/KO 삽입: ✅
- Berry(1988)/Odlyzko(1987) bibitem 추가: ✅
- Conclusion "28 results across five verification axes" 업데이트: ✅
- 컴파일 성공: ✅ (EN pdflatex, KO xelatex)

---

### 논문 최종화 — 남은 편집 과제 분석

사이클 39에서 본문 반영은 완료. 이번 사이클에서 점검한 미비 사항:

#### 1. **Abstract 미비** ⚠️
현재 abstract (tcolorbox, L130–L166)는 QROP-Net, 7 correspondences, 음성 결과 (PQO 실패) 중심. 28개 수치 결과와 위상적 일관성 체인에 대한 언급 **없음**. Conclusion에는 반영됨(L4008–4019)이나 abstract에는 빠짐.

→ Abstract에 "28 numerical results across five axes" 1문장 추가 필요.

#### 2. **Introduction Contributions 미비** ⚠️
현재 Contributions 목록 (L221–258)은 8항목으로, 사이클 30 이전 수준. ξ-다발 프레임워크의 위상적 검증 (Chern 수, GB 면적분, σ-국소화, 합성 대조군, 수 분산)이 언급되지 않음.

→ 9번째 contribution 항목 추가 필요.

#### 3. **28개 결과 Summary Table 부재** ⚠️
58페이지 논문에 28개 결과가 흩어져 있으나, 독자가 한눈에 볼 수 있는 종합 표가 없음. 결과 번호, 이름, 판정(확립/조건부/음성), 핵심 수치, obs 라벨을 한 표에 정리해야 함.

→ Discussion/Conclusion 직전에 "Summary of Numerical Evidence" longtable 삽입.

#### 4. **우선순위 판단**

세 과제 중 **Summary Table이 가장 가치 높음**:
- 심사자가 가장 먼저 찾을 것 (58페이지에 28개 결과 — 안내 표 필수)
- Abstract/Contributions 수정은 1–2문장씩이라 같이 처리 가능

---

### 다음 작업: **28개 결과 Summary Table 삽입 + Abstract/Contributions 업데이트**

**모델**: sonnet

**왜**: 논문 최종화 모드에서 가장 큰 편집 미비. 심사자 관점에서 28개 결과를 한눈에 볼 수 있는 종합 표가 없으면 논문의 navigate가 어려움. Abstract과 Contributions는 동시에 1–2문장 추가로 해결 가능.

**구체 설계**:

1. **Summary Table (EN/KO)**: Conclusion 직전에 longtable 삽입
   - 열: # | Name | Status | Key Metric | Section/Label
   - 28행 (결과 #1–#28)
   - Status: Established / Conditional / Negative 3종
   - 라벨: `tab:summary`

2. **Abstract 추가 (EN/KO)**: "Note added in revision" 문단 앞에 1문장 삽입
   - "Twenty-eight numerical results across five verification axes—local curvature, global spectral statistics, topological invariants, Dirichlet extension, and negative controls—support the framework's internal consistency."

3. **Contributions 추가 (EN/KO)**: 기존 8항목 뒤에 9번째 항목
   - "We verify the fibre-bundle interpretation through 28 numerical experiments spanning five independent axes (Section~X)."

4. **컴파일 확인**: EN pdflatex / KO xelatex 에러 없음 확인

**주의**:
- longtable 패키지 추가 필요 (preamble)
- 결과 판정 표현: "확립"="Established", "조건부 양성"="Conditionally positive", "음성"="Negative"
- 결과 #4(H-스케일링, 중립), #9(S¹ 고높이, 중립) 정확 반영 — 양성 편향 금지
- 결과 #17(Gram 점, 음성) 정확 반영 — 음성 결과 숨기지 말 것
- 표가 2페이지 넘지 않도록 간결하게 (핵심 수치 1개만)

**성공 기준**:
- Summary Table EN/KO 삽입 + 컴파일 성공
- Abstract에 "28 results" 문구 포함
- Contributions에 9번째 항목 추가
- 기존 obs 라벨과 충돌 없음

---

### 28개 결과 Summary Table 데이터 (설계자 참조용)

| # | Name (EN) | Status | Key Metric | Axis |
|---|-----------|--------|------------|------|
| 1 | Kuramoto synchronization | Established | r=0.9994±0.0001 | Local |
| 2 | S¹ geodesic loss | Established | 4× detection gain | Local |
| 3 | High-height scaling | Established | recall 97.5%+ to t=1100 | Local |
| 4 | H-scaling exponent | Inconclusive | α≈−0.4 | Local |
| 5 | FP anatomy | Established | 85.1% quasi-zeros | Local |
| 6 | Conj 3 κ-separation | Established | FP max < TP min, 300× | Local |
| 7 | σ=1/2 uniqueness | Established | energy ratio 385× | Local |
| 8 | Blind prediction | Established | F1=1.000, 7/7 | Local |
| 9 | S¹ high-height | Inconclusive | ΔF1=+0.1%p | Local |
| 10 | Dirichlet bundle | Positive | 3 chars × 4 properties | Dirichlet |
| 11 | Dirichlet cross-comparison | Positive | 4 funcs × 5 properties | Dirichlet |
| 12 | Dirichlet blind | Positive | F1=1.000, 21/21 | Dirichlet |
| 13 | High conductor | Positive | 5/5 properties | Dirichlet |
| 14 | Real character energy | Positive | χ₄=132.5×, χ₈=194.3× | Dirichlet |
| 15 | Off-critical anatomy | Positive | Conj 1 Dirichlet generalization | Dirichlet |
| 16 | Lehmer pair curvature | Cond. positive | ρ=0.835 | Global |
| 17 | Gram point curvature | Negative | No κ-enhancement | Negative |
| 18 | Midpoint nonlocality | Positive | ρ=−0.654 (p<10⁻²⁹) | Global |
| 19 | Dirichlet midpoint nonlocal | Positive | 3/3 |ρ|>0.55 | Global |
| 20 | GUE/Poisson contrast | Positive | GUE ρ(b)≈0 | Global |
| 21 | Nonlocal reach decay | Cond. positive | k=1 concentration | Global |
| 22 | Holonomy / Chern number | Established | ζ 87/87, χ₃ 4/4, χ₄ 4/4 | Topological |
| 23 | Even Chern + methodology | Established | χ₅ 4/4, χ₈ 2/2 | Topological |
| 24 | Gauss–Bonnet area integral | Established | ζ 6/6, χ₅ 2/2, error=0 | Topological |
| 25 | Off-critical σ-localization | Established | ζ 7/7, FWHM<0.004 | Topological |
| 26 | Triple topological consistency | Established | 4×3 matrix, N=100 | Topological |
| 27 | Synthetic negative control | Established | 11/11, 3 reproductions | Negative |
| 28 | Number variance (RMT) | Cond. positive | Var_plaq≈Var_dir, Berry sub-GUE | Global |

---

## 확립된 결과 (28개)

(이전 사이클과 동일 — 변경 없음)

### 알려진 함정 (누적)

- 모노드로미: eps 차분 금지 → 폐곡선 적분(radius=0.5, 64단계)
- Dirichlet character `a` 파라미터: χ(-1)=+1→a=0, χ(-1)=-1→a=1
- κ 측정: 영점 위 직접 측정 금지 → δ=0.03~0.05 오프셋
- findroot: 단일 시작점만, 튜플 전달 금지
- numpy 2.0: np.trapezoid (trapz 아님)
- **bundle_utils.py h=10^{-20}**: dps ≥ 100 또는 해석적 공식 사용.
- **중간점 Hadamard 상쇄**: m=(γₙ+γₙ₊₁)/2에서 두 최근접 영점의 기여 = 0.
- **ξ(s) underflow**: log(ξ)' 해석적 공식으로 영구 해결.
- **mpmath.diff(xi, s)**: t>300에서 실패. 사용 금지.
- **Dirichlet κ 계산**: bare L'/L 사용 금지. Λ'/Λ 해석적 공식 필수.
- **윤곽 적분 시 영점 회피**: t₂는 반드시 영점 사이 중간점.
- **홀수 디리클레 지표 영점 탐색**: Re(L(1/2+it))=0 ≠ L=0. |L| 최소화 필수.
- **짝수 디리클레 지표**: a=0, Γ(s/2) 사용. L(1/2+it) ≠ 실수 (Λ만 실수). Re(L)=0 단독 불충분.
- **영점 탐색 보편적 방법**: |L(1/2+it)| 최소화가 Re(L)=0보다 엄밀히 우월 (A⊂B). 패리티 무관.
- **F₂ 면적분 정의**: ∂_σ Im(L)은 대칭 상쇄 → 격자 플라켓(arg ξ 와동 밀도) 사용. Stokes 등가.
- **σ-국소화 표현**: "RH의 기하학적 동치" 금지 → "RH와 일관적인 기하학적 증거" 수준으로. FWHM=0 → "FWHM < δσ".
- **다항식 log(0) 회피**: 격자점이 영점을 정확히 지나지 않도록 범위 조절.
- **Number variance unfolding**: L은 unfolded 좌표(평균 밀도로 정규화). raw t-간격 사용 금지.
- **Number variance sub-GUE**: t<600에서 GUE 수렴 기대 금지. Berry(1988) arithmetic correction 레짐. "GUE 재현" 표현 금지.
