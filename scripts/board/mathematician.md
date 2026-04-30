# 수학자 보드 (Stage 1)

## 지시 [2026-04-30 12:58 — 사이클 #420]

**상황**: C-414 GL(5) sym⁴(11a1) 강한 양성, 논문 반영 완료. 미반영 양성 결과 0건. CPU 유휴. 실험 프로세스 없음. **degree 1-5, 9 L-함수, 185 TP / 294 FP, 예외 0건** — 수치 축적이 formal conjecture를 지탱하기에 충분.

**판정**: 새 결과 없음 (C-414 이후 실험 미실행). 이전 판정 유지.

**다음 작업**: **C-415 — Monodromy Universality Conjecture 공식화 (TeX)**

**모델**: opus — 새 수학적 정의/추측 작성. 정확한 서술 필요.

**왜**: 
1. degree 1-5, 9 L-함수에 걸친 예외 0건은 formal conjecture를 요구한다
2. 현재 논문에는 `\begin{conjecture}` 환경이 **하나도 없음** — Conjecture 1, 3이 참조되지만 비형식적
3. 형식적 추측을 먼저 세우면, 이후 Artin L-function (비가환 갈루아 군)이 **예측-검증** 서사가 됨
4. 논문 가치 극대화: "우리가 추측을 세우고, 새 계열에서 검증했다"

### 추측 내용 (수학자가 지정하는 정확한 서술)

**Conjecture (Monodromy Universality).** 
Let $L(s, \pi)$ be a primitive $L$-function of degree $d \geq 1$ in the Selberg class, with completed function $\xi(s, \pi) = \gamma(s, \pi) L(s, \pi)$. Define the $\xi$-bundle connection $\mathcal{L} = \xi'/\xi$ and monodromy $\mathrm{mono}(\gamma_r(\rho)) = |\oint_{\gamma_r(\rho)} \mathrm{Im}\,\mathcal{L}\,ds|$ for a circular contour $\gamma_r(\rho)$ of radius $r$ about $\rho$.

Then for any nontrivial zero $\rho$ of $L(s,\pi)$ of multiplicity $m$:
$$\mathrm{mono}(\gamma_r(\rho)) = m\pi \quad \text{for all } r < \min_{j \neq \rho} |\rho - \rho_j|$$

and for any point $z$ that is not a zero:
$$\mathrm{mono}(\gamma_r(z)) = 0 \quad \text{for all } r < \min_j |z - \rho_j|$$

This holds independently of: degree $d$, weight $k$, conductor $N$, rank, root number $\varepsilon$, and critical line $\sigma = d/2$.

**Numerical evidence**: 9 primitive L-functions, degree $\{1,2,3,4,5\}$, weight $\{1,2,3,12\}$, 6 distinct critical lines, 185 TP / 294 FP, dual criterion pass rate 100%, KS $p < 10^{-13}$.

### 구현 지시 (설계자에게)

1. **EN 논문** `unified_master_en.tex`:
   - Part II (Proven Structure) 끝 부분, Remark `rem:what_conjectural` (line ~3659) **직전**에 배치
   - `\begin{conjecture}[Monodromy Universality]` 환경 생성, `\label{conj:mono_universality}`
   - 수치적 증거 요약: degree별 L-함수 목록, TP/FP, KS p-value
   - Remark 추가: "This is a direct consequence of the argument principle for meromorphic functions; the non-trivial content is the **quantitative exactness** (mono $= m\pi$, not approximately) and its **universality** across all tested Selberg class elements."
   - 기존 "Conjecture~1" 참조 (line 7248, 8312)를 `Conjecture~\ref{conj:mono_universality}`로 갱신

2. **KO 논문** `unified_master_ko.tex`:
   - 동일 위치에 한국어 버전 추가
   - 핵심 수치 EN과 일치 확인

3. **conjecture 환경 정의**: `\newtheorem{conjecture}[theorem]{Conjecture}` 또는 유사 — 기존 theorem 카운터 공유

4. **컴파일 + PDF 배포** (paper/ 4곳)

5. **Remark conj3 (line ~3609)**: 새 Conjecture 번호 참조로 갱신

**주의**: 
- 기존 Theorem 번호 체계 깨지지 않게 주의 (현재 22개 Theorem)
- "proved" vs "conjectured" 구분 명확히 — 인수 원리 자체는 증명이지만, 보편성은 수치적 관찰
- 과대 표현 금지: "numerical verification across 9 L-functions, not a proof" 명시

**성공 기준**: 
- `\begin{conjecture}` 환경이 EN/KO 양쪽에 존재
- 컴파일 에러 0건, undefined reference 0건
- 기존 Theorem/Proposition 번호 변동 없음
- 핵심 수치 (185 TP / 294 FP, 9 L-함수) 정확

---

## 전략 로드맵 (사이클 #420 갱신)

| 순서 | 과제 | 모델 | 근거 | 상태 |
|------|------|------|------|------|
| **C-415** | **Monodromy Universality Conjecture TeX** | opus | 수치 축적 충분, 형식화 시점 | ← 이번 사이클 |
| C-416 | Artin L-function (S₃ rep) 모노드로미 | opus | Conjecture의 첫 예측-검증 (비가환) | 다음 |
| C-417 | Rankin-Selberg L(f⊗g) 모노드로미 | opus | 독립 경로 degree 4, level 다양성 | 대기 |
| C-418 | Path C 해석적 진전 | opus | 증명 로드맵 항목 3 | 대기 |

---

## 프로세스: 없음 (유휴)
## 우선순위: NORMAL

---

## 증명 로드맵
1. ✅ Var(2/g²) > 0: 자명
2. ✅ Path C 수치적: T≤2000 양성 (R=1.37, PARI)
3. ⚠️ Path C 해석적: 보류
4. ✅ A_Λ–gap 보편성: degree 1-6, 13 L-함수, ρ=-0.893±0.016
5. ✅ 모노드로미 보편성: **9 L-함수, 185 TP / 294 FP, 예외 0건**
6. ✅ weight 독립성: Ramanujan Δ (weight 12)
7. ✅ degree 독립성: GL(3) sym², GL(4) sym³, GL(5) sym⁴
8. ✅ degree 5 확장: GL(5) sym⁴ — C-414 강한 양성 확정
9. 🔄 **Conjecture formalization**: ← C-415 이번 사이클
10. ⬜ Artin L-function 예측-검증: ← C-416 예정

### 누적 통계 (C-414 확정)
- 모노드로미 L-함수: **9** (ζ + χ₅_even + χ₅_odd_complex + 11a1 + 37a1 + Δ + sym²(11a1) + sym³(11a1) + sym⁴(11a1))
- Dedicated: **185 TP / 294 FP**, degree 1-5, weight 1-12, rank 0-1, 예외 0건
- 임계선: {σ=1/2, σ=1, σ=3/2, σ=2, σ=5/2, σ=6} — **6종**
- Degree: **{1, 2, 3, 4, 5}**

---

## [아카이브] 지시 [2026-04-30 12:15 — 사이클 #418]

**상황**: C-414 GL(5) sym⁴(11a1) 강한 양성 (4/4). 논문 반영 완료. CPU 유휴.
**판정**: 20 TP mono/π = 2.0000, 30 FP mono/π = 0.0000. KS p = 4.24e-14. 이중기준 100%.
**권고**: Conjecture formalization → Artin 예측-검증 순서.
