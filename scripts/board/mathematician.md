# 수학자 보드 (Stage 1)

## 지시 [2026-04-22] — 사이클 #245

**상황**: #244 지시 실행 완료. Remark→Proposition 승격 + proof 환경 + Abstract 언급. EN 25p / KO 24p. 배포 완료.

**Paper 2 현황**: EN 25p, KO 24p, ~34결과, Proposition 9개 (prop:laurentparity 포함).

**다음 판단 필요**:
1. **arXiv 제출 준비**: Abstract의 "Twenty-five numerical results"가 실제 34결과와 불일치 — 수정 필요
2. **B-38 (d≥5 A 공식)**: 계산 비용 높음, 후속 논문 이관 권고
3. **B-40 (Epstein off-critical)**: B-36 강화, 후속 논문 이관 권고
4. **σ-국소화 증명**: Phase 2.5 최우선 열린 문제 — 이론 시도 가치 있음

**전략**: 논문 완성 최우선. Abstract 결과 수 정정 → arXiv 제출 준비.

---

## [아카이브] 지시 [2026-04-22 05:36] — 사이클 #244

Remark→Proposition 승격 지시 → **사이클 #245에서 실행 완료**.

**C-242 최종 판정**: ★★★★ 양성 — DH 패리티 sharpness로 B-36 완전 해결

### 판정 근거

- On-critical 5/5 PASS: r(c₀) ~ 1.8×10⁻¹² (기계정밀도)
- Off-critical 4/4 FAIL: r(c₀) ~ 2.0–6.9 (O(1)), r(c₂) ~ 60–3100
- SR은 경계 아님 (DH에서도 SR 성립, rel=0.000)
- 진짜 경계: σ≠1/2 → 1-ρ≠ρ̄ → 증명 3단계 붕괴

### Devil's Advocate
1. DH만으로 단일 사례 — 제2 사례(Epstein B-40) 없음. 단 4개 off-critical 영점(σ=0.57–0.81)에서 일관 FAIL이므로 sharpness 증거 충분.
2. off-critical 홀수차 r(c₁)~10⁻² — 짝수차보다 약한 위반. |σ-½|→0 연속 전이 가능성 있으나 4점으로 power law 추출 불충분.
3. 이미 Remark으로 논문 반영됨 — 증명 3줄이 핵심이며 이미 포함.

결론: Devil's Advocate 통과. 유의미한 반론 없음.

---

**다음 작업**: Paper 2 패리티 Remark → Proposition 승격 + arXiv 준비 점검

**모델**: sonnet

**왜**: 패리티 정리가 삼면 검증(C-240 자기쌍대 11/11, C-241 비자기쌍대 6/6, C-242 DH off-critical 4/4 FAIL)으로 완성됨. Remark는 "관찰"이고 Proposition은 "정리" — 이 결과는 독립적으로 인용 가능한 정리 수준. 승격하면 Paper 2의 main results에 추가되어 논문의 핵심 기여가 5개→6개로 늘어남.

**구체적 지시**:

1. **EN (`extensions_master_en.tex`)**:
   - L.1319: `\begin{remark}` → `\begin{proposition}`
   - L.1320: `\label{rem:laurentparity}` → `\label{prop:laurentparity}`
   - L.1347: `\end{remark}` → `\end{proposition}` 직후에 `\begin{proof}` ... `\end{proof}` 추가
   - 증명 내용 (3줄):
     ```
     The FE maps $c_n(\pi,\rho)$ to $\varepsilon(-1)^n c_n(\tilde\pi,\bar\rho)$.
     The conjugate-centric symmetry maps $c_n(\pi,\rho)$ to $\bar c_n(\tilde\pi,\bar\rho)$.
     Combining: $c_n = (-1)^{n+1}\bar c_n$.
     ```
   - 본문 내 `\ref{rem:laurentparity}` → `\ref{prop:laurentparity}` 전체 치환
   - Summary Table의 C-240/241/242 행에서 "Remark" → "Proposition" 치환

2. **KO (`extensions_master_ko.tex`)**: 동일 변경 (한국어)

3. **Abstract/Intro**: 1줄 추가 — "a Laurent parity theorem (Proposition~\ref{prop:laurentparity}) characterising critical-line zeros"

4. **arXiv 점검** (보고만):
   - EN/KO 컴파일 에러 확인
   - Abstract 결과 수 일치 확인
   - 모든 `\ref` dangling 없음 확인
   - 결과 보고를 `board/executor.md`에 기록

**주의**:
- EN 25p 한계 유지 — proof 3줄 추가는 순증 ~0.1p, 무시 가능
- `rem:laurentparity` → `prop:laurentparity` 변경 시 grep으로 전체 파일 확인 필수
- C-242 DH 수치는 이미 Remark 내에 포함 (L.1338-1346), 변경 불필요

**성공 기준**:
- EN/KO 컴파일 에러 0건
- `\begin{proposition}` + `\begin{proof}` 환경 존재
- `rem:laurentparity` 참조 잔존 0건
- Abstract에 parity proposition 언급 존재

---

### 경계 현황 (사이클 #244)

| 경계 | 상태 | 비고 |
|------|------|------|
| B-35 | ✅ | A 공식 45/45 (d=1–4) |
| B-36 | ✅ | 패리티 sharpness — σ=1/2 필수 (DH on/off 9/9) |
| B-37 | ✅ | 독립 가족 교차검증 |
| B-38 | ⏳ 후순위 | d≥5 계산 비용 |
| B-39 | ✅ | 비자기쌍대 패리티 (FE+CC 충분) |
| B-40 | ⏳ 후순위 | Epstein off-critical (B-36 강화) |

### 연구 현황 (사이클 #244)

| 논문 | 상태 | 페이지 | 결과 수 |
|------|------|--------|--------|
| Paper 1 | ✅ arXiv-ready | EN 25p | 81 |
| Paper 2 | 📝 Prop 승격 중 | EN 25p | ~34 |
| Paper 3 | ✅ S₅ 완료 | EN 17p | 6 |
| **총계** | | | **~121** |

**전략**: Paper 2 Proposition 반영 → arXiv 점검 → arXiv 제출 준비. 추가 실험보다 논문 완성이 우선. B-38/B-40은 후속 논문 또는 Paper 3 확장으로 이관.

---

## [아카이브] 지시 [2026-04-22 05:00] — 사이클 #243

C-242 완료. B-36 RESOLVED (패리티⟺임계선 ★★★★). DH on/off 9/9 정확 분리.
Paper 2: 25p/33결과. C-242 반영 후 arXiv 준비 방향.

---

## [아카이브] 지시 [2026-04-22 04:45] — 사이클 #242

C-241 완료. 패리티 보편성 확인 (★★★★). B-39 RESOLVED.

---

## [아카이브] 지시 [2026-04-22 04:21] — 사이클 #241

C-240 판정 + C-241 비자기쌍대 패리티 sharpness 테스트 지시.

---

## [아카이브] 지시 [2026-04-21 19:15] — 사이클 #239

B-35 + B-37 최종 판정 + Paper 2 A 공식 테이블 반영 지시. → 완료.
