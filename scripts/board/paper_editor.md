# 편집장 지시 [2026-04-17 02:36] — 논문 사이클 #1

## 판정: UPDATE_NEEDED

## 근거

GL(3) 확장 결과 #63~#67이 수학자 보드에서 모두 확립 판정을 받았으나, 논문에는 Discussion §Future Directions에 "GL(3) symmetric-square lifts is a natural question for a separate investigation" 한 줄만 존재. **5개 확립 결과가 전혀 반영되지 않음.** 이는 논문의 가장 큰 갭이다.

교정자 보고서: 초기 상태 (이전 사이클 없음). 교정자 지적사항 없음.

---

## 반영할 결과

| # | 결과명 | 카테고리 | 판정 | 반영 위치 (섹션) | 분량 | 우선순위 |
|---|--------|---------|------|----------------|------|---------|
| 63 | GL(3) sym²(11a1) 4성질 검증 | A (GL(3)) | ★★ 양성 | 부록: GL(3) Extension 신설 | 1.5p | 필수 |
| 64 | GL(3) sym²(37a1) 4성질 검증 | A (GL(3)) | ★★ 양성 | 부록: GL(3) Extension 신설 | 1p (63과 병합 표) | 필수 |
| 65 | GL(3) 블라인드 영점 예측 | A (GL(3)) | ★★★ 강양성 | 부록: GL(3) Extension | 0.5p | 필수 |
| 66b | GL(3) σ-유일성 5-Level 역보정 | A (GL(3)) | ★ 음성 (핵심) | 부록: GL(3) Extension + 한계 소절 | 0.5p | 필수 |
| 67 | GL(3) FP 21개 영점 검증 | A (GL(3)) | ★★★ 강양성 | 부록: GL(3) Extension (Zero-Finding) | 1p | 필수 |

**총 추가 분량**: ~4.5p (부록)

---

## 섹션 변경 지시

### 1. 신규 부록 섹션 생성: "Beyond GL(2): GL(3) Symmetric-Square $L$-functions"

**위치**: `\section{Beyond GL(1): GL(2)...}` (EN 5908행 / KO 4343행) 바로 뒤, `\section{Negative Result: Gram Point Curvature}` 앞.

**구조**:
```
\section{Beyond GL(2): GL(3) Symmetric-Square $L$-functions}
\label{app:gl3}

[도입부: GL(3) sym²(E) L-함수 설명, conductor N_E², degree 3, 감마 인자 2개, AFE 방법론]

\subsection{Four-Property Verification: sym²(11a1) and sym²(37a1)}
[#63, #64 결과 표]
- 함수방정식 잔차: 0.0 (두 곡선 모두)
- 영점 수: 63개 (11a1), 44개 (37a1)
- κ near/far: 283× (11a1), 222× (37a1) — 임계값 100× 초과
- 모노드로미/π: 12/12=2.000 (11a1), 15/15=2.000 (37a1)
- σ-유일성: FAIL (두 곡선 모두 — GL(3) 특성, B-01 참조)

\subsection{Blind Zero Prediction and the Zero-Finding Machine}
[#65, #67 결과]
- 블라인드: Recall=1.000 (LMFDB 15/15), 21개 "FP"
- #67 AFE 검증: 21/21 전부 실제 영점 (median |Λ|=2.5e-27)
- 삼중 확인: 부호변환(#63) + κ-피크(#65) + AFE(#67)
- **핵심 발견**: LMFDB 미수록 영점 독립 발견 → 교정 P=1.0, F1=1.0

\subsection{$\sigma$-Uniqueness Failure: Dirichlet Series Convergence Origin}
[#66b 결과]
- 5-Level 역보정 모두 FAIL
- 원인: 감마 인자가 아닌 Dirichlet 급수 수렴 특성
- GL(1)/GL(2)/GL(3) σ-유일성 패턴 비교
```

### 2. Summary Table에 GL(3) 결과 추가 (#50~#54)

**위치**: Summary Table 마지막 행 (#49) 뒤, `\end{longtable}` 앞에 추가.

**추가할 행**:
```latex
\midrule
50 & GL(3) sym$^2$(11a1) 4-property & \textbf{E} & $\kappa$ $283\times$; mono 12/12; FE${}=0.0$; $\sigma$-uniq.\ FAIL & GL(3) \\
51 & GL(3) sym$^2$(37a1) 4-property & \textbf{E} & $\kappa$ $222\times$; mono 15/15; FE${}=0.0$; $\sigma$-uniq.\ FAIL & GL(3) \\
52 & GL(3) blind zero prediction & \textbf{E} & $R{=}1.000$, $P{=}1.000$, $F_1{=}1.000$ (corrected); 21 new zeros & GL(3) \\
53 & GL(3) $\sigma$-uniqueness 5-level deconv. & \textbf{N} & 5/5 FAIL; origin${}={}$Dirichlet series convergence & GL(3) boundary \\
54 & GL(3) FP verification (AFE) & \textbf{E} & 21/21 zeros; median $|\Lambda|{=}2.5{\times}10^{-27}$ & GL(3) \\
```

또한 Table caption을 "Summary of 49 numerical results" → "Summary of 54 numerical results" 로 수정.

### 3. Discussion §Future Directions 수정

**EN 5434행**: "Extension to Maass forms and GL(3)" 문단을 업데이트.

**기존 내용** (삭제 금지, 확장):
> "Whether the framework extends to Maass forms ... and higher-rank L-functions such as GL(3) symmetric-square lifts is a natural question for a separate investigation."

**추가 내용**: GL(3) sym² 결과를 달성했음을 언급하고 app:gl3 참조. Maass form은 여전히 열린 문제로 남김.

변경 후:
> "For GL(3) symmetric-square lifts, we have now obtained affirmative results: see Appendix~\ref{app:gl3}. Extension to Maass forms (weight~0, non-holomorphic) remains an open problem."

---

## 구조 변경

- **신규 섹션 1개**: 부록에 GL(3) 섹션 추가 (GL(2) 섹션 뒤)
- **Part 재배치**: 불필요 (부록 내 추가)
- **분리 트리거 점검**: GL(3) 결과 5개 — paper_categories.md 기준으로는 Paper A에 유지 (GL(3) 독립 분리 조건인 "≥5개 확립 결과 또는 ≥15p 분량"에 도달했으나, #68 κ-conductor 스케일링 결과 대기 중이므로 이번 사이클에서는 분리하지 않고 Paper A 부록에 배치. 다음 사이클에서 재검토.)

---

## EN/KO 동기화 참고

- **GL(3) 내용**: EN/KO 양쪽 모두 없음 (동기 상태). 양쪽 모두에 동일하게 추가.
- **EN**: unified_master_en.tex 5908행 뒤 GL(3) 부록 추가
- **KO**: unified_master_ko.tex 4343행 뒤 GL(3) 부록 추가 (한국어 번역)
- **Summary Table**: EN/KO 양쪽 모두 #50~#54 동일하게 추가
- **Discussion**: EN 5434행 / KO 3937행 양쪽 수정

---

## 이전 교정자 지적사항 반영

해당 없음 (초기 상태, 이전 사이클 없음).

---

## 데이터 소스

집필자는 다음 결과 파일을 참조하여 정확한 수치를 인용할 것:

| 결과 파일 | 대상 |
|-----------|------|
| `results/gl3_sym2_11a1_63.txt` | #63 (11a1 4성질) |
| `results/gl3_sym2_37a1_64.txt` | #64 (37a1 4성질) |
| `results/gl3_blind_prediction_65.txt` | #65 (블라인드) |
| `results/gl3_gamma_deconv_66b.txt` | #66b (5-Level 역보정) |
| `results/gl3_fp_verification_67.txt` | #67 (FP 영점 검증) |

---

## 주의사항

1. **내용 삭제 절대 금지** — Discussion의 기존 GL(3) 언급은 삭제하지 말고 확장만 할 것
2. **수학 표기 일관성** — 기존 매크로 사용 (\xif, \Ftwo 등). GL(3) sym²는 `\mathrm{sym}^2` 사용
3. **GL(2) 섹션의 서술 패턴을 따를 것** — 동일한 구조 (도입 → Key differences → Validation → Four-property table → Discussion)
4. **κ 비율 해석 주의** — #63 far 샘플 3점뿐 (B-06 경계 존재). 이 한계를 각주로 언급
5. **#65 Precision 이중 표기** — 원래 P=0.417에서 교정 P=1.0 (#67에 의해). 두 값 모두 투명하게 보고
6. **#66b는 음성이지만 핵심 발견** — 부정적으로 쓰지 말고 "진단" 관점으로 서술 (GL(3) σ-uniqueness FAIL 원인 규명)
7. **#68 (κ-conductor 스케일링)은 아직 미완** — 이번 사이클에서는 "preliminary 2-point observation" 으로만 언급, 확정적 서술 금지

---
## 집필자 실행 보고 [2026-04-17 02:39]
- EN: GL(3) 부록 신설 (§Beyond GL(2): GL(3) Symmetric-Square), Summary Table #50–#54 추가, Discussion GL(3) 참조 추가, 6204행 → 6377행 (+173)
- KO: GL(3) 부록 신설 (§GL(2)를 넘어서: GL(3) 대칭 제곱), Summary Table #50–#54 추가, Discussion GL(3) 참조 추가, 4615행 → 4778행 (+163)
- 컴파일: EN PASS (85p, 에러 0), KO PASS (75p, 에러 0)
- 배포: PASS (paper/ + ~/Desktop/수학최종논문/)
- git: 3b02525
