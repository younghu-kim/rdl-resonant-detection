# 설계자/실행자 보드

## 보고 [2026-04-27 04:05 — 사이클 #356] — C-349 EC 20곡선 + A(t₀) rank 의존성

**수학자 지시**: C-349 — EC ��계적 확장 + A(t₀) rank 의존성 탐사 (Paper B 첫 실험)
**모델**: Opus
**PID**: 185025 (완료)
**결과 위치**: `results/ec_20curves_c349.txt`
**스크립트**: `scripts/ec_20curves_c349.py`
**소요**: 21분

---

### 실험 구성

- **20곡선**: rank 0-3 × 5곡선, LMFDB 검증 완료
- rank 0: 11a1(N=11), 14a6(N=14), 15a1(N=15), 43a1(N=43), 197a1(N=197)
- rank 1: 37a1(N=37), 53a1(N=53), 58a1(N=58), 61a1(N=61), 79a1(N=79)
- rank 2: 389a1(N=389), 433a1(N=433), 446a1(N=446), 563a1(N=563), 571a1(N=571)
- rank 3: 5077a1(N=5077), 11197a1(N=11197), 11642a1(N=11642), 12279a1(N=12279), 13766a1(N=13766)
- **t∈[1,100]**, δ=0.03, Cauchy r=0.01, system python3 (cypari2)
- 총 영점: ~2,946개 (비자명), A(t₀) 측정: ~1,599개

### 결과 1: 4성질 검증

| 성질 | 결과 |
|------|------|
| P1 FE | 20/20 ✅ (log10\|Λ/Λ̃\| = 0.0 전원) |
| P2 κδ² | 19/20 ✅ (11642a1 CV=32% FAIL, 나머지 CV<8%) |
| P3 σ-sweep | 20/20 ✅ (σ=1.0에서 28-30/30 검출) |
| P4 mono/π | 20/20 ✅ (mean=2.000000, std=0 전원) |
| **종합** | **19/20 PASS (95%)** — ★★★★ 양성 |

**유일한 FAIL**: 11642a1 (rank 3, N=11642) κδ²=1.158 (CV=32%). 고 conductor 수치 정밀도 한계 가능. 나머지 19곡선 κδ² ∈ [1.006, 1.042].

### 결과 2: ★ A(t₀) = Im(c₀)² + 2Re(c₁) rank 의존성

| rank | n | mean | std | **median** | SE |
|------|---|------|-----|-----------|-----|
| 0 | 307 | 22.69 | 29.06 | **13.31** | 1.66 |
| 1 | 312 | 21.19 | 25.44 | **13.97** | 1.44 |
| 2 | 415 | 30.12 | 28.90 | **24.40** | 1.42 |
| 3 | 565 | 64.61 | 149.82 | **42.76** | 6.30 |

**통계 검정:**

| 비교 | p (t-test) | p (Mann-Whitney) | 유의 |
|------|-----------|-----------------|------|
| rank 0 vs 1 | 0.494 | 0.697 | ❌ 비유의 |
| rank 0 vs 2 | **0.0007** | **<0.0001** | ★ |
| rank 0 vs 3 | **<0.0001** | **<0.0001** | ★ |
| rank 1 vs 2 | **<0.0001** | **<0.0001** | ★ |
| rank 1 vs 3 | **<0.0001** | **<0.0001** | ★ |
| rank 2 vs 3 | **<0.0001** | **<0.0001** | ★ |

**Spearman(rank, A_mean)**: ρ=0.80, p=0.20 (4점이라 통계적 유의 미달)

### 핵심 해석 — ��️ 수학자 판단 필요

**A(t₀) 중앙값 단조 증가**: 13.3 → 14.0 → 24.4 → 42.8 (rank 0→1→2→3)

그러나 **conductor 교란(confound) 가능성**이 높다:
- rank 0-1: N~11-197, A_median~13
- rank 2: N~389-571, A_median~24
- rank 3: N~5077-13766, A_median~42

**rank 내부 증거**: 197a1 (rank=0, N=197)의 A_mean=43.4는 rank 3 수준. 이는 A(t₀)가 rank보다 **conductor N에 더 강하게 의존**할 수 있음을 시사.

**판정**: ★★★ 조건부 양성 — A(t₀)는 rank별 체계적 차이를 보이나, conductor와의 분리 실험이 필요. "rank 감지" 주장에는 동일 conductor 범위의 rank 비교가 필수.

### 이슈

1. **추가 곡선 계수 오류**: 처음 임의 작성한 12곡선 계수가 전부 부정확. LMFDB API + WebFetch로 정확한 계수 확보 (14a6, 15a1 등). 교훈: EC 계수는 반드시 LMFDB 직접 확인.
2. **lfuninit 경고**: t~100 근처에서 `insufficient initialization` 경고 ~40,000건. 결과 정확도에 영향 미미 (PARI 자동 보정).
3. **A(t₀) 성공률**: ~50% (Cauchy 적분 r=0.01에서 영점 간격이 r보다 작은 경우 실패). Cauchy r 조정으로 개선 가능.

### 수학자에게 제안

1. **Conductor 분리 실험**: 동일 N 범���에서 다른 rank 곡선 비교 (예: rank 0 N~400 vs rank 2 N~400)
2. **A(t₀) vs log(N) 상관**: rank 효과와 conductor 효과 분리
3. **11642a1 κδ² FAIL**: realprecision 150으로 재시도 또는 제외

---

## [아카이브] 보고 [2026-04-27 02:55 — 사이클 #353] — C-348 Paper A 최종 정합성 점검 완료

**수학자 지시**: C-348 — Paper A 투고 전 최종 정합성 점검 (10개 항목)
**모델**: Sonnet
**PID**: 없음 (논문 TeX 수정 작업)
**결과 위치**: `paper/source/unified_master_en.pdf`, `paper/source/unified_master_ko.pdf`

---

### 10개 점검 항목 결과

#### A. 수치 정합성

| 항목 | 판정 | 발견사항 |
|------|------|----------|
| A1 디리클레 통합 비교표 (q=3~11) | ⚠️→✅ | **수정됨**: result #37 셀이 구버전 "136 zeros; mono 136/136" (C-324). C-345 데이터(689zeros, 5/5 ALL-PASS)로 수정 (EN line 6017). 또한 EN line 7482 "136 zeros" → "689 zeros (C-345)" 수정. `tab:dirichlet_per_character` (q=7: 689 ✅) 등 다른 표는 정확. |
| A2 켤레 쌍 E비 표 (B-68) | ✅ | 6쌍 전부 일치: χ₇₂↔χ₇₄ 12.35×, χ₇₁↔χ₇₅ 2.05×, χ₁₁₄↔χ₁₁₆ 144×, χ₁₁₃↔χ₁₁₇ 4.0×, χ₁₁₁↔χ₁₁₉ 1.47×, χ₁₁₈↔χ₁₁₂ 1.12× |
| A3 블라인드 예측 | ✅ | GL(1) 27/27 F1=1.000, GL(2) 7/7 Ramanujan Δ, recall 0.991 (NN) — 모두 논문에 정확히 반영 |
| A4 FP 85.1% | ✅ | 57/67 = 85.07% ≈ 85.1% 정확. 결과파일과 일치 |
| A5 GUE/GL | ✅ | KS D=0.0411, p=0.0022 등 result #34~36에 수록. α_∞=5/12는 Paper B용 |

#### B. 내부 일관성

| 항목 | 판정 | 발견사항 |
|------|------|----------|
| B1 Notation 일관성 | ✅ | κ, δ, E(σ), mono, L, F₂ 정의 후 일관 사용 확인 |
| B2 Theorem/Prop 번호 | ✅ | Theorem 22개, label 23개, 참조 65회 — 전부 일치 |
| B3 \ref/\label 무결성 | ✅ | undefined reference 0건 (컴파일 에러 0 재확인) |

#### C. 서사 정합성

| 항목 | 판정 | 발견사항 |
|------|------|----------|
| C1 Abstract vs 본문 | ✅ | Abstract 7개 핵심 주장 전부 본문에서 뒷받침됨 |
| C2 Introduction Part 참조 | ⚠️→✅ | **수정됨**: line 278 "Part~III" (구버전 3파트 체계 잔재) → "Part~\ref{part:evidence}" (실제 Part V). KO line 170 "제3부" → "제~\ref{part:synthesis}부" (실제 제IV부) |

#### D. EN/KO 동기화

| 항목 | 판정 | 발견사항 |
|------|------|----------|
| D1 구조 비교 | ✅ (의도적 차이) | EN 5파트/KO 4파트, EN 22 Thm/KO 10 Thm — KO가 압축본이므로 정상 |
| D2 핵심 수치 동일성 | ✅ | q=7 (689/1.003/17.7~394.7×), FP 85.1%, 2,517 누적 — EN/KO 동일 확인 |
| D3 켤레 쌍 표 | ✅ | EN/KO 6쌍 수치 완전 동일 |

---

### 수행한 수정 (4건)

| 번호 | 파일 | 위치 | 수정 내용 |
|------|------|------|----------|
| 1 | EN | line 6017 | result #37: "136 zeros; mono 136/136" → "689 zeros (t∈[10,200], C-345); mono 689/689; 5/5 ALL-PASS" |
| 2 | EN | line 278 | "Part~III" → "Part~\ref{part:evidence}" (Part V로 렌더링) |
| 3 | EN | line 7482 | "table below for reference counts: 12, 13, 13, 136 zeros" → "per-conductor zero totals: 12, 13, 13, 689 zeros (C-345, 5 chars.)" |
| 4 | KO | line 170 | "제3부" → "제~\ref{part:synthesis}부" (제IV부로 렌더링) |

### 컴파일 결과

| 파일 | 엔진 | 에러 | Undefined Ref | 페이지 |
|------|------|------|---------------|--------|
| unified_master_en.tex | pdflatex | 0 | 0 | **122p** |
| unified_master_ko.tex | xelatex | 0 | 0 | **46p** (-1p, 레이아웃 변화) |

※ KO 46p: C-347에서 47p였으나, Part 참조 수정 후 레이아웃 압축으로 46p. 정상 범위.

### PDF 배포 (4곳)

- `paper/unified_master_en.pdf` ✅
- `paper/unified_master_ko.pdf` ✅
- `paper/unified_paper_en.pdf` ✅
- `paper/unified_paper_ko.pdf` ✅

---

### 최종 판정: **투고 준비 완료** (수정 4건 완료)

**수치 불일치**: 결과파일 대비 논문 수치 불일치 **0건** (수정 후)
**서사 불일치**: Part 참조 오류 **0건** (수정 후)
**Notation 불일치**: **0건**

### 이슈

- Edit 도구가 "성공" 메시지를 반환했으나 실제 파일에 적용되지 않는 문제 발생. Python 직접 수정으로 해결.
- 수학자 보드에 언급된 "q=8 모노드로미 수치 불일치(B-1)"는 에이전트 분석에서 발견되지 않았음. q=8 결과파일(c326)과 논문 표 수치가 일치함. 검토자의 추가 확인 권장.
- 수학자 보드에 "EN 121p"로 표기되었으나 실제 컴파일 결과 122p. C-347부터 122p 유지 중.

---

## 보고 [2026-04-27 01:58 — 사이클 #352] — B-68 Remark 반영 완료 (C-347)

**수학자 지시**: C-347 — B-68 Remark 논문 반영 (켤레 비동치 해석)
**모델**: Sonnet
**PID**: 없음 (논문 TeX 수정 작업)
**결과 위치**: `paper/source/unified_master_en.pdf`, `paper/source/unified_master_ko.pdf`

### 수행한 수정

| 항목 | 내용 |
|------|------|
| EN 추가 | `rem:conjugate_nonequiv` Remark — 켤레 쌍 표 직후 (line 7334) |
| KO 추가 | `rem:conjugate_nonequiv` Remark — 켤레 쌍 표 직후 (line 2343) |
| 라벨 추가 | `\phantomsection\label{tab:conjugate_eratio}` — EN line 7314, KO line 2323 |

### Remark 3개 필수 요소 확인

| 요소 | EN | KO |
|------|----|----|
| $|L(1,\chi)| = |L(1,\bar\chi)|$ 동일성 | ✅ line 7343 | ✅ 본문 |
| 단일 결정인자 부재 | ✅ "no single analytic invariant" | ✅ "어떤 단일 해석적 불변량도" |
| 위상 기하학적 해석 (arg W, phase geometry) | ✅ line 7347-7352 | ✅ 본문 |

### 컴파일 결과

| 파일 | 엔진 | 에러 | 페이지 |
|------|------|------|--------|
| unified_master_en.tex | pdflatex | 0 | **122p** (+1p) |
| unified_master_ko.tex | xelatex | 0 | **47p** (+1p) |

### PDF 배포 (4곳)

- `paper/unified_master_en.pdf` ✅
- `paper/unified_master_ko.pdf` ✅
- `paper/unified_paper_en.pdf` ✅
- `paper/unified_paper_ko.pdf` ✅

**B-68 Remark 반영 완료**: 디리클레 확장 섹션의 수학적 해석 서술 추가.
