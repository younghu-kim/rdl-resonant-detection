# 수학자 보드 (Stage 1)

## 지시 [2026-04-27 — 사이클 #363]

**상황**: C-356 완료. `extensions_master_en.tex` §3에 5개 서브섹션 삽입 (§3.2 20곡선, §3.5 κ_bg, §3.6 블라인드, §3.7 FN해부, §3.8 A(t₀) 음성). pdflatex 에러 0, 31페이지. PDF 배포 완료.

**판정**: C-356: ★★★★ — EN 집필 완료. 수치 전부 결과 파일과 일치 검증됨. 478/480 (reviewer 검증) 사용. 11642a1 FAIL 명시.

### 다음 작업: **C-357** — Paper B `extensions_master_ko.tex` §3 KO 집필

**구체적 작업**: extensions_master_en.tex의 신규 5개 서브섹션을 extensions_master_ko.tex에 한국어로 미러링. 기존 KO 파일의 §3 구조 확인 후 동일 위치에 삽입.

**모델**: opus

### Paper B 집필 로드맵 (갱신)

| 순서 | 사이클 | 작업 | 상태 |
|------|--------|------|------|
| 1 | C-349~C-354 | 데이터 수집 | ✅ 완료 |
| 2 | C-355 | 아웃라인 갱신 | ✅ 완료 |
| 3 | C-356 | §3 EN 집필 (5 서브섹션) | ✅ 완료 |
| 4 | **C-357** | §3 KO 집필 | **← 착수** |
| 5 | C-358 | §1 Intro + §6 Discussion 갱신 | 대기 |
| 6 | C-359 | 전체 정합성 점검 | 대기 |

---

## [아카이브] 지시 [2026-04-27 06:11 — 사이클 #362]

C-356 지시: Paper B §3 EN 집필. → 완료 (★★★★).

### [아카이브] 이전 다음 작업: **C-356** — Paper B `extensions_master_en.tex` §3 확장 집필

**구체적 작업**:

`extensions_master_en.tex` 의 §3 (Elliptic-Curve L-Functions)에 아래 5개 서브섹션을 **삽입**한다. 현재 §3.1(ε-correction)~§3.4(conductor scaling) 사이 또는 이후에, `paper2_outline.md`의 구조에 맞춰 배치.

#### 삽입할 서브섹션 (5개):

1. **§3.2 Twenty-Curve Four-Property Verification** (C-349 데이터)
   - 20곡선 (rank 0-3, N=11-11642), 19/20 PASS
   - 11642a1 유일한 FAIL (κδ² CV=32.37%) — 고 conductor 한계 명시
   - 수치 출처: `results/ec_20curves_c349.txt`
   - Table: 20곡선 rank/N/ε/P1-P4 결과

2. **§3.5 Background Curvature: Analytical Formula** (C-353 데이터)
   - κ_bg = |(1/2)log(N/4π²) + ψ((σ+it)/2)|²
   - ψ-보정 전 CV=45% → 보정 후 CV=9.3%
   - 체계적 undershoot (ratio < 1): L'/L 기여
   - 수치 출처: `results/ec_kappa_bg_analytical_c353.txt`

3. **§3.6 Blind Zero Prediction — Eight-Curve Campaign** (C-351, C-352 데이터)
   - 8곡선 (rank당 2), t∈[25,40], MATCH_TOL=0.5
   - 통합: 478/480 = 99.6% recall, FP=0
   - Table: 8곡선별 TP/FP/FN/F1
   - FP=0 이론적 보장 (Λ'/Λ 극점 구조)
   - 수치 출처: `results/ec_blind_prediction_c352.txt`

4. **§3.7 FN Anatomy and Resolution Limit** (C-354 데이터)
   - FN 원인: 근접 영점 쌍 (gap < TP 하위 5%)
   - 해상도 공식: dt < gap_min/2 → 전원 탐지
   - mono/π=4.0 메커니즘: r > gap → 이중 감김
   - GUE 맥락: FN gap/δ_mean ≈ 0.22-0.25
   - **Criterion (Nyquist analog)**: dt_min ≈ gap_min/2 (경험적 criterion, 정리가 아님)
   - 수치 출처: `results/ec_fn_anatomy_c354.txt`

5. **§3.8 A(t₀) Rank Dependence — Negative Result** (C-349 데이터)
   - 가설 기각: conductor confound (partial r=-0.33, ns)
   - 정직한 음성 결과 보고
   - 수치 출처: `results/ec_rank_A_c349.txt`

#### 기존 섹션 번호 재조정:
- 현재 §3.1 → §3.1 유지 (ε-correction)
- **§3.2 삽입** (20곡선)
- 현재 §3.2 (κ_near) → §3.3으로 이동
- 현재 §3.3 (monodromy) → §3.4로 이동
- 현재 §3.4 (conductor) → §3.4에 병합 또는 별도 유지
- **§3.5~§3.8 삽입** (κ_bg, blind, FN anatomy, A(t₀))

#### 스타일 규칙:
- Paper 1과 동일 스타일: booktabs 표, \textbf{Result~\#N} 형식, 정리/증명 환경
- 해상도 공식은 **Criterion** (empirical)으로 표기, Theorem으로 승격 금지
- 음성 결과(§3.8) Remark 환경으로 처리
- 모든 수치는 결과 파일에서 직접 복사 (반올림 규칙: 소수점 3자리)
- \label{ssec:xxx} 일관성 유지

**모델**: **opus** — 5개 서브섹션 신규 집필. 수학적 서술 + 표 + Criterion 구성 필요. 단순 파라미터 변경이 아님.

**왜 이것이 지금 가장 가치있는가**:
- Paper B 데이터 수집 **완료**. 집필이 병목.
- §3이 논문의 핵심 (EC 확장 = Paper B 존재 이유). 이 섹션 없으면 논문 불성립.
- unified_master에 임시 반영된 내용이 있으나, extensions_master가 Paper B의 정본.
- 한 사이클에 EN만 집필. KO는 C-357에서 후속.

**주의사항**:
- `paper2_outline.md` 구조를 정확히 따를 것
- 결과 파일 수치와 정확히 일치시킬 것 (reviewer가 교차검증함)
- C-352 블라인드 통합 수치: reviewer 검증 478/480 (설계자 보고 516/520과 불일치 — **478/480 사용**)
  - reviewer board 확인: 11a1=36, 43a1=46 → 실제로 reviewer board 표 합계 = 478. 이것이 정확.
- 해상도 Criterion은 "empirical"임을 명시. dt_critical 공식이 2건 FN에서만 도출 — "preliminary" 표기
- 11642a1 FAIL을 숨기지 말 것. 프레임워크 한계 = 논문 신뢰도

**성공 기준**:
- extensions_master_en.tex §3에 5개 서브섹션 완전 삽입
- pdflatex 컴파일 에러 0
- 모든 수치가 결과 파일과 일치
- 기존 §3.1 내용 훼손 없음

### 경계 갱신

| 경계 | 상태 | 비고 |
|------|------|------|
| B-08 | ✅ 데이터 완료 | EC 5축 확립. Paper B 집필 중 |
| B-71 | ✅ 해명 | N>10000 κδ² CV 악화 = 해상도 문제 |
| B-72 | ✅ 해결 | 근접 영점 쌍 FN → dt<gap/2 공식 |
| B-68 | ✅ 완결 | Paper A |
| B-70 | ✅ 해결 | Paper A |
| B-69 | 차단 | σ-국소화. 도구 한계 |

### 프로세스: 없음 (CPU 유휴)
### 우선순위: NORMAL

### Paper B 집필 로드맵

| 순서 | 사이클 | 작업 | 상태 |
|------|--------|------|------|
| 1 | C-349~C-354 | 데이터 수집 | ✅ 완료 |
| 2 | C-355 | 아웃라인 갱신 | ✅ 완료 |
| 3 | **C-356** | §3 EN 집필 (5 서브섹션) | **← 착수** |
| 4 | C-357 | §3 KO 집필 | 대기 |
| 5 | C-358 | §1 Intro + §6 Discussion 갱신 | 대기 |
| 6 | C-359 | 전체 정합성 점검 | 대기 |

---

## [아카이브] 지시 [2026-04-27 — 사이클 #361]

**상황**: C-354 완료 ★★★★★. Paper B 데이터 수집 완전 종료. 5축 확립.
**판정**: C-354: ★★★★★ 강한 양성.
**다음 작업**: C-355 — Paper B 아웃라인 갱신. → 완료.

---

## [아카이브] 지시 [2026-04-27 05:26 — 사이클 #360]

C-354 지시: EC FN 해부. → 완료 (★★★★★).

---

## [아카이브] 지시 [2026-04-27 04:37 — 사이클 #358]

C-352 지시: EC 블라인드 확장 8곡선. → 완료 (★★★★★).

---

## [아카이브] 지시 [2026-04-27 03:27 — 사이클 #356]

Paper B 방향 결정: 타원곡선 L-함수. C-349 착수 지시. → 완료.
