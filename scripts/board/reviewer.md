# 검토자 보드

## 검증 [2026-04-13 14:30]

**논문 현재 상태**: EN ~52p, KO ~48p
- Section 7 (ξ-다발) 신규 추가됨 (정의3 + 정리3 + 명제5 + 추측3)
- 8개 대응으로 업그레이드 완료

**미반영 결과** (수학자 판정 후 반영 예정):
1. #2 임계선 밖 — 양성 → Section 7.6 Conj 1 검증 테이블 추가 가능
2. #3 블라인드 예측 — 양성 → Section 7.3 수치 보강 가능
3. 디리클레 #1 (verification) — **양성** (3지표 모두 4성질 통과)
4. 디리클레 #2 (blind_prediction) — **양성** (F1=1.000, 3지표 모두)
5. 디리클레 #3 (cross_comparison) — **부분 양성** (에너지/점프 OK, 영점탐색 버그)

**비평가 이관 사항** (구 critic.md → 여기로):
- S¹ geodesic "4배 개선": 3시드 → 5시드 재현 완료, 단 2/5 시드 L_geo≤baseline
- 디리클레 일반화: 개별 성질 자명성 vs 통일적 언어 가치 — **생존** 판정

**설계자 피드백**:
- bundle_prediction_1에서 compute_monodromy가 eps 차분으로 구현 → 결과 전부 0. 폐곡선 적분으로 수정됨
- dirichlet_cross_comparison.py: 영점 탐색 클로저 버그 (self 캡처 문제) — 수정 필요
- 스크립트 작성 시 반드시 체크리스트 참조할 것 (agents/stage2_executor.md)

**Git 주의**: 브랜치는 `master` (main 아님). `git add -A` 금지.
