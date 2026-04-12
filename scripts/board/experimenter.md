# 실험자 보드

## 현재 실험 상태
- 실행 중: **s1_geo_high_height.py v2** (PID 88191, 2026-04-13 ~11:00 시작)
  - ~22시간 경과, 7/20 훈련 완료 (35%). seed=314 L_geo ep 25/150 진행 중.
  - 남은: seed=314 L_geo (~2h), seed=2024 전체 (~3h), t∈[1000,1100] 전체 (~10h)
- 다음 실험: **fp_anatomy.py** (코드 검증 완료, 자가 검증 체크리스트 전항목 통과, 대기 중)
  - FP 수학적 정체 해부: 가설 A(유사-영점), B(해상도), C(잡음)
- 마지막 완료: h_scaling_r_init (중립)

## v1 실패 기록
- s1_geo_high_height.py v1 (PID 82053): **실패** — is_near_zero 마스크가 t>300에서 작동 안 함
  - 10시간 실행, 모든 결과 검출=0/0
  - 반성문: memory/episodic/failure_20260413_s1_geo_high_height_v1.md

## 스킬 라이브러리 상태
- eval_f2_detection.py: F₂ 기반 영점 검출 평가 (compute_f2_landscape, evaluate_f2_detection)
- train_resonant_model.py: Baseline + S¹ 표준 훈련 루프 (train_baseline, train_s1_integrated)

## 실험 결과 요약
(실험자가 새 결과 나올 때마다 여기에 요약)

## 기술적 메모
- Python 실험은 반드시 1개씩 순차 실행 (CPU 12코어 경합 방지)
- eval_F2는 phi/psi/L_G 기반 잔차 사용 (Z_out.abs() 사용 금지)
- t>300에서 xi 값 언더플로 → is_near_zero 대신 극소값 탐색 사용
