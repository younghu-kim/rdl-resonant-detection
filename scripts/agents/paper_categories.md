# 논문 카테고리 분류 체계

**작성**: 2026-04-15
**목적**: 결과를 연구 질문별로 분류하여 논문 라우팅 및 분리 판단에 사용

---

## 카테고리 정의

### Paper A: ξ-다발 기하학 (ξ-Bundle Geometry)

**연구 질문**: ξ-다발의 위상적 구조가 리만 영점을 어떻게 인코딩하는가?

**키워드**: Klein, monodromy, curvature, connection, σ-uniqueness, parallel transport, Gauss-Bonnet, Chern, holonomy, off-critical, blind prediction, FP anatomy, Hardy Z phase, gauge

**현재 결과**:
| # | 결과 | 판정 |
|---|------|------|
| 1 | σ=1/2 유일성 (385×) | 확립 |
| 2 | 블라인드 예측 7/7 | 확립 |
| 3 | 곡률 TP/FP 300× | 확립 |
| 5-6 | 곡률 집중 | 확립 |
| 7-8 | FP 해부 | 확립 |
| 9 | Hardy Z phase | 확립 |
| 10-11 | Monotone-κ 평행수송 | 확립 |
| 19-21 | 비국소 곡률 | 확립 |
| 22-23 | Chern 수 (홀로노미) | 확립 |
| 24 | Gauss-Bonnet 면적분 | 확립 |
| 25-26 | Off-critical σ-국소화 | 확립 |
| 27 | 합성 함수 음성 대조군 | 확립 |
| 30 | DH 적대적 검증 | 약한 양성 |
| 31 | 고 t 스케일링 | 확립 |
| 32 | κ 차선도 구조 | 조건부 양성 |
| 33 | Hadamard 분해 검증 | 양성 |
| 35 | δ-sweep Lorentzian² 보편성 | 조건부 양성 |
| 49 | FP 모노드로미 해부 (Conj 3) | 진행중 |

---

### Paper B: 스펙트럼 통계 (Spectral Statistics & RMT)

**연구 질문**: ξ-다발 곡률이 GUE 통계 및 영점 분포와 어떻게 관련되는가?

**키워드**: GUE, Poisson, spacing, Weyl, counting, pair correlation, number variance, RMT, Δκ, NNS, distribution

**현재 결과**:
| # | 결과 | 판정 |
|---|------|------|
| 4 | GUE 상관 (pair) | 확립 |
| 28 | 곡률장 수 분산 | 확립 |
| 34 | Δκ 잔차-GUE 간격 상관 | 양성 |
| 36/36b | spacing ratio GUE | 음성 |
| 37 | 곡률-GUE 분포 변수변환 | 음성 |
| 38 | κ 자기상관 vs GUE 두점 상관 | 음성 |

**참고**: 음성 결과 3개는 "경계 명확화" 역할 — GUE와의 직접 대응에는 한계가 있음을 보여준다. 이것도 논문에 포함 (negative results section).

---

### Paper B-ext: 디리클레 확장 (Dirichlet Extension)

**연구 질문**: ξ-다발 구조가 디리클레 L-함수에도 보편적으로 적용되는가?

**키워드**: Dirichlet, character, mod, conductor, L-function, GL(1), universality

**현재 결과**:
| # | 결과 | 판정 |
|---|------|------|
| 14-18 | Dirichlet 확장 (5종, mod 3,4,5) | 확립 |
| 29 | Epstein zeta 확장 | 조건부 양성 |
| 40/40b | 고 conductor 디리클레 mod 7 | 양성 |
| 42 | 디리클레 교차비교표 (5종) | 양성 |

**참고**: Paper A에 보편성 증거로 포함하거나, 분량이 크면 별도 섹션. 현재는 Paper A에 통합.

---

### Paper C: GL(2) 확장 (Elliptic Curve L-functions)

**연구 질문**: ξ-다발이 GL(2) (타원곡선 L-함수)까지 확장되는가? BSD 추측과의 연결은?

**키워드**: elliptic, GL(2), BSD, L-function, modular, 11a1, 37a1, Ramanujan, Δ, weight, conductor, AFE, σ-uniqueness

**현재 결과**:
| # | 결과 | 판정 |
|---|------|------|
| 44 | 타원곡선 11a1 GL(2) 검증 | 양성 (3/4) |
| 45 | 타원곡선 37a1 GL(2) 검증 | 양성 (3/4) |
| 46 | 라마누잔 Δ GL(2) 검증 | 양성 (3/4) |
| 47 | σ-유일성 분기 메커니즘 | 조건부 양성 |

**분리 상태**: 결과 4개 — 트리거 근접 (≥3). 현재는 Paper A 마지막 섹션에 포함, 축적되면 분리.

---

### Paper D: 양자 시뮬레이션 (Quantum Simulation)

**연구 질문**: ξ-다발의 위상적 구조를 양자 시뮬레이션으로 검증할 수 있는가?

**키워드**: Hamiltonian, DQPT, VQE, qubit, echo, Floquet, Bethe, Ansatz, topological code, stabilizer, Berry-Keating

**현재 결과**: 코드 완성 (quantum/ 디렉토리), 아직 논문 미반영

**분리 상태**: 독립 저장소 준비 단계. 코드 + 문서만 존재.

---

## 분류 규칙

### 키워드 매칭 우선순위
1. GL(2)/elliptic/BSD/modular → Paper C
2. Hamiltonian/DQPT/VQE/qubit/echo/Floquet → Paper D
3. GUE/Poisson/spacing/Weyl/RMT/distribution → Paper B
4. 나머지 → Paper A (기본)

### Bridge Result 처리 (병렬 기입)
복수 카테고리 해당 시 **양쪽 저장소에 동시 기입**:
- 주 카테고리: 상세 내용 (전체 증명/실험/데이터)
- 부 카테고리: 요약 1-2문단 + "see companion paper [ref] for details"

**병렬 기입 예시**:
- "GL(2)에서도 곡률 집중 성립" → 주: Paper C (GL(2)), 부: Paper A (요약만)
- "Bethe Ansatz ↔ Gauss-Bonnet 등가" → 주: Paper D (양자), 부: Paper A (요약만)
- "ξ-다발 프레임워크 정의 수정" → 모든 논문에 동기화

### S¹ 결과 (#12-13) 특수 처리
- S¹ L_geo 5시드: 조건부 양성, 재현성 경고
- Paper A에 포함하되 주의사항 명시

---

## 분리 트리거

| 조건 | 동작 |
|------|------|
| 단일 논문 본문 > 25p (부록 제외) | 분리 검토 |
| 한 카테고리 결과 ≥ 8개 | 분리 검토 |
| 새 카테고리(C/D) 결과 ≥ 3개 | 분리 검토 |

**현재 상태**:
- Paper A: ~20개 결과 → 분리 불필요 (메인 논문)
- Paper B: 6개 결과 → Paper A에 통합 유지
- Paper C: 4개 결과 → 트리거 충족, 분리 준비
- Paper D: 코드만 → 결과 생산 후 분리
