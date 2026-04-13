## 실패 반성: bundle_prediction_1_fp_monodromy.py
**유형**: Type A (구현 버그)
**날짜**: 2026-04-13 ~05:19 시작, ~13:30 kill

**원인**: 5가지 체크리스트 위반
1. `out_features=1` → 2 필수 (MasterResonantNetwork 사양)
2. `loss_fn(Z_out=..., X_in=..., ...)` → `loss_fn(**outputs)` 패턴 필수 (TotalResonanceLoss 사양)
3. `channel_type`, `damping_mode` 미지정 → `"paper3ch"`, `"paper"` 필수
4. eval에서 `torch.no_grad()` 사용 → `torch.enable_grad()` 필수 (phi/psi 기반 평가)
5. `X_in.requires_grad_(True)` 누락

**추가 문제**: 
- `-u` 플래그 미사용으로 stdout 버퍼링 → 진행 상황 추적 불가
- OMP_NUM_THREADS=10으로 과다 설정 → 동시 실행 시 CPU 경합
- 결과 경로가 outputs/analysis/ → results/로 변경 필요

**교훈**: 
- 스크립트 작성 시 **반드시** 자가 검증 체크리스트 통과 후 실행
- 이전 사이클에서 다른 스크립트(s1_geo_high_height 등)에서 확립된 패턴을 새 스크립트에도 동일 적용
- `-u` 플래그와 nohup 로그 리다이렉트는 실행 필수 조건

**대안**: 수정 완료. s1_geo_high_height 완료 후 재실행 예정.
