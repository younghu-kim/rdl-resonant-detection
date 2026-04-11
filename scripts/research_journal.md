# RDL 자율 연구 일지

## 2026-04-11 23:48 사이클
**상황**: 3개 실험 모두 정상 실행 중, 새 결과 없음. complex_vector 벽시계 4h19m/CPU 17h13m (seed=42 baseline 아직 진행 중, 6 runs 중 1번째), high_height 벽시계 6h50m/CPU 27h07m, precision 벽시계 6h50m/CPU 26h59m. 모든 프로세스 Rl 상태, CPU 각 ~396%.
**판단**: 휴식. CPU 12코어 포화, 결과 없음, 모든 프로세스 건강. complex_vector의 run당 ~3.75h 추정 대비 4h19m 경과 — 첫 baseline run 완료 임박 가능. 로그 줄 수 변화 없음(17줄).
**실행**: 상황 파악 + 연구 일지 업데이트
**다음**: complex_vector 첫 run 완료 시 로그에 새 출력 예상. complex_vector_readout.txt 도착 시 즉시 Phase 1. 예상 잔여: complex_vector ~18h, precision ~12h, high_height ~20h+.

## 2026-04-11 23:15 사이클
**상황**: 3개 실험 모두 정상 실행 중, 새 결과 없음. complex_vector 벽시계 3h45m/CPU 15h (seed=42 baseline 아직 1/6), high_height 벽시계 6h17m/CPU 25h (t100-200 캐시 존재, t500-600 미생성), precision 벽시계 6h17m/CPU 25h. 시스템 load 28.88, 모든 프로세스 Rl/395%.
**판단**: 휴식. CPU 포화, 결과 없음, 프로세스 건강. high_height는 t100-200 학습 완료 후 t500-600 영점 계산 또는 아직 학습 중. complex_vector 총 예상 ~22h 벽시계 (3.75h × 6 runs).
**실행**: 상황 파악 + 연구 일지 업데이트
**다음**: complex_vector_readout.txt 도착 시 즉시 Phase 1. 예상 잔여: complex_vector ~18h, precision ~12-18h, high_height ~24h+.

## 2026-04-11 22:43 사이클
**상황**: 3개 실험 모두 정상 실행 중, 새 결과 없음. complex_vector 벽시계 3h14m/CPU 12h52m (seed=42 baseline 아직 진행 중, 6 runs 중 1번째), high_height 벽시계 5h45m/CPU 22h45m, precision 벽시계 5h45m/CPU 22h42m. 모든 프로세스 Rl 상태, CPU 각 ~395%. 로그 변화 없음.
**판단**: 휴식. CPU 포화, 결과 없음, 모든 프로세스 건강. complex_vector의 seed=42 baseline 1개에 3h14m+ 소요 중 → 6 runs 총 ~19h 벽시계 예상. precision/high_height는 CPU 22h+ 소비에도 결과 미생성, 아직 학습 중반~후반 추정.
**실행**: 상황 파악 + 연구 일지 업데이트
**다음**: complex_vector_readout.txt 도착 시 즉시 Phase 1. 예상 잔여: complex_vector ~16h, precision ~10-15h, high_height ~20h+.

## 2026-04-11 22:11 사이클
**상황**: 3개 실험 모두 정상 실행 중, 새 결과 없음. complex_vector 2h42m (seed=42 baseline 아직 진행 중, 6 runs 중 1번째), high_height 5h14m, precision 5h14m. 시스템 load 28.97, 12코어 포화. 로그 줄 수 변화 없음 (complex_vector 17줄, 나머지 2줄).
**판단**: 휴식. CPU 포화, 결과 없음, 모든 프로세스 정상 R 상태. 예상 완료: complex_vector ~3-4h 후 (총 ~6h), precision ~6-10h, high_height ~10-15h+.
**실행**: 상황 파악 + 연구 일지 업데이트
**다음**: complex_vector_readout.txt 도착 시 즉시 Phase 1. 양성이면 Tier1 항목1 [confirmed] + 전체 파이프라인 통합 설계, 음성이면 S¹ 기하학적 각도 예측 등 대안 검토.

## 2026-04-11 21:40 사이클
**상황**: 3개 실험 모두 정상 실행 중, 새 결과 없음. complex_vector 2h10m/CPU 8h39m (seed=42 baseline 아직 진행 중), high_height 4h42m/CPU 18h35m, precision 4h42m/CPU 18h29m. 시스템 load 29.89, 모두 R 상태. complex_vector 로그 17줄 변화 없음 — 훈련 내부 출력 없어 run 완료까지 대기.
**판단**: 휴식. CPU 포화, 결과 없음, 모든 프로세스 정상.
**실행**: 상황 파악 + 연구 일지 업데이트
**다음**: complex_vector_readout.txt 도착 시 즉시 Phase 1. 예상 잔여 4-8h.

## 2026-04-11 21:07 사이클
**상황**: 3개 실험 모두 정상 실행 중, 결과 파일 없음. complex_vector 1h38m (seed=42 baseline 진행 중, 6 runs 중 1번째), high_height 4h10m, precision 4h10m. 시스템 load 30.19, 12코어 포화. 로그 확인: complex_vector 17줄(초기 메시지만), 나머지 2줄(경고만).
**판단**: 휴식. CPU 완전 포화, 새 결과 없음. 훈련 함수에 중간 출력 없어 완료까지 로그 변화 없을 것. 예상 잔여: complex_vector ~7h, precision ~8-12h, high_height ~15h+.
**실행**: 상황 파악 + 연구 일지 업데이트
**다음**: complex_vector_readout.txt 도착 시 즉시 Phase 1. 로그 줄 수 변화(17→)로 진행 확인 가능.

## 2026-04-11 초기화
**상황**: 자율 연구 에이전트 시스템 구축. 6개 실험 완료, 3개 실행 중.
**확립된 결과**: Prop 8.5 필요조건, ±π 유일 병목, F₂ 일반화 99.1% recall, 매끄러운 타겟 전체 파이프라인 음성, 결정론적 푸리에 우위.
**실행 중**: complex_vector_readout, high_height_scaling, precision_filter (PID 40428, 40479, 40528)
**다음**: 3개 실험 완료 대기 → 결과 분석 → 논문 반영 → 후속 실험 판단

## 2026-04-11 18:48 에러
**에러**: exit code 127. 로그: /home/k0who029/Desktop/gdl_unified/outputs/auto_research_logs/run_20260411_184821.log

## 2026-04-11 18:50 사이클
**상황**: 3개 실험 모두 정상 실행 중 (PID 40429/40480/40529). 2시간 경과, 결과 파일 미생성. stdout 버퍼링으로 로그 출력 없음. embedding.py:67 경고는 초기화 시 1회 발생, 비위험. 12코어 전부 포화 상태.
**판단**: 휴식. 4번째 실험 시작은 CPU 경합으로 역효과. high_height_scaling이 가장 느릴 것 (t∈[500,600], [1000,1100] 영점 계산 + 캐시 빌드 필요). 논문 현재 상태 점검하여 결과 반영 준비.
**실행**: 상황 파악 + 연구 일지 업데이트. 논문 구조 확인 예정.
**다음**: 결과 파일 도착 대기. complex_vector_readout.txt가 가장 먼저 나올 가능성 높음 (기존 캐시 사용, 학습만 수행). 결과 도착 시 즉시 Phase 1 수행.

## 2026-04-12 15:10 사이클
**상황**: 3개 실험 모두 22시간째 실행 중 (PID 40429/40480/40529). CPU 각 390%, 메모리 377-411MB, R (running) 상태. stdout 버퍼링으로 /tmp/*.log에는 초기 경고만 기록. 결과 파일 미생성. t500-600/t1000-1100 캐시 파일도 미생성 → high_height_scaling은 아직 영점 계산 또는 초기 구간 학습 중.
**판단**: 휴식. 근거:
  - CPU 12코어 전부 포화 (3 프로세스 × ~4코어), 4번째 실험 시작 시 전체 성능 저하
  - 논문 구조 사전 검토 완료: complex_vector → Tier1 항목1 (L2820-2828), high_height → Tier2 항목1 (L2866-2871) + sec:f2_residual_open (L2320-2326), precision_filter → Tier2 항목2 근처 신규 추가 필요
  - 예상 완료 순서: complex_vector (가장 빠름, 기존 캐시+6 학습) > precision_filter (5시드×200ep+밀집격자) > high_height (영점계산+캐시빌드+18 학습, 가장 느림)
**실행**: 상황 파악 + 논문 삽입 지점 매핑 + 연구 일지 업데이트
**다음**: complex_vector_readout.txt 도착 시 즉시 Phase 1 수행. 양성이면 Tier1 항목1을 [confirmed]로 업데이트 + 전체 파이프라인 통합 실험 설계. 음성이면 대안 ±π 해결책 (예: 위상 해제 unwrapping, 코사인 기반 대리 등) 검토.

## 2026-04-11 25:30 사이클 (※ 이전 에이전트 날짜 오류 — 실제 프로세스는 4/11 16:57 시작)
**상황**: 이전 에이전트가 22시간/44시간 경과로 기록했으나 실제 `ps -o lstart` 확인 결과 프로세스 시작 시각은 2026-04-11 16:57. 날짜 혼동이 있었음.
**병목 분석**: CurvatureLoss가 매 배치마다 `torch.autograd.grad(create_graph=True)`로 2차 미분 계산 → 학습 속도 5-10배 감속. 3 프로세스 CPU 경합 추가 감속. t500-600/t1000-1100 캐시 미생성 → high_height는 영점 계산 또는 초기 학습 중.
**논문 삽입 지점**: complex_vector→L2820-2828, high_height→L2866-2871, precision→Tier1 근처 추가
**향후 실험 설계 시 고려**: `-u` 플래그(unbuffered stdout), 중간 체크포인트 저장, 단일 실험씩 순차 실행

## 2026-04-11 19:25 사이클
**상황**: 3개 실험 2시간 26분째 실행 중 (PID 40429/40480/40529). 모두 R 상태, CPU 각 390%, 메모리 377-421MB. 시스템 load 30.55 (12코어 완전 포화). 결과 파일 미생성. 로그에 초기 embedding.py 경고만 기록 (stdout 버퍼링). t100-200 캐시(27KB) 존재, t500-600/t1000-1100 캐시 미존재.
**판단**: 휴식. 근거:
  - CPU 12코어 포화 (load 30.55), 추가 실험 시작은 역효과
  - 이전 에이전트의 날짜 혼동(4/12, 4/13 항목) 정정: 실제 runtime은 2.5시간
  - 예상 완료 순서: complex_vector (기존 캐시+3시드×2방법×150ep) > precision_filter (5시드×200ep+밀집격자) > high_height (영점계산+캐시빌드+3범위×2K×3시드×150ep, 가장 느림)
  - 논문 삽입 지점 매핑 이미 완료
**실행**: 상황 파악 + 날짜 오류 정정 + 연구 일지 업데이트
**다음**: 결과 파일 도착 대기. complex_vector_readout.txt 가장 먼저 예상.
  - 양성: Tier1 항목1 [confirmed] + 전체 파이프라인 통합 실험 설계
  - 음성: 대안 ±π 해결책 탐색 (위상 해제, cos/sin 기반 대리 등)

## 2026-04-11 19:23 사이클
**상황**: 3개 실험 ~2.5시간 실행 중. 코드 리뷰 수행.
**발견**: complex_vector_readout.py에 **치명적 형상 불일치 버그** 발견!
  - `Z_out`은 `[batch, out_features=2]` complex → `z_flat=[batch,4]`
  - `readout = nn.Linear(128, 2)` → 128차원 입력 기대, 4차원 수신 → RuntimeError
  - 원인: `Z_out`은 out_projection 후(2차원)인데, readout은 pre-projection hidden(64차원) 기준
**수정**: `model.out_projection`에 forward hook → pre-projection Z 캡처 → z_flat=[batch,128] 형상 일치. `-u` 플래그 추가.
**실행**: PID 40429 종료 → 스크립트 수정 → 형상 검증 통과 → PID 46434로 재실행
**실행 중**: PID 46434 (complex_vector, 수정판), PID 40480 (high_height, 2.5h), PID 40529 (precision, 2.5h)
**다음**: complex_vector 로그 모니터링. 첫 baseline 완료 후 complex vector 단계 정상 진입 확인이 핵심.

## 2026-04-11 20:33 사이클
**상황**: 3개 실험 모두 정상 실행 중. 새 결과 파일 없음.
  - complex_vector (PID 46435): 1h 05min 경과, CPU 4h 19min. seed=42 baseline 훈련 중 (1/6 run). `-u` 플래그 작동, 로그 17줄 — 훈련 함수 내부에 진행 출력 없어 완료까지 출력 없을 것.
  - high_height (PID 40480): 3h 36min 벽시계, CPU 14h 14min. 18 runs × 150ep + t500/t1000 영점계산/캐시 빌드. stdout 버퍼링.
  - precision_filter (PID 40529): 3h 36min 벽시계, CPU 14h 07min. 5 runs × 200ep + 밀집 격자. stdout 버퍼링.
**추정 잔여**: complex_vector ~5-6h, precision_filter ~5-10h, high_height ~15-20h+
**판단**: 휴식. CPU 포화, 새 결과 없음, 논문 삽입 지점 매핑 완료. 대기가 최적.
**다음**: complex_vector_readout.txt 가장 먼저 도착 예상. 결과 도착 시 즉시 Phase 1.

## 2026-04-11 20:02 사이클
**상황**: 3개 실험 모두 정상 실행 중.
  - complex_vector (PID 46435): 32분 경과, seed=42 baseline 훈련 중 (150ep, 30ep마다 출력). `-u` 플래그 작동 확인.
  - high_height (PID 40480): 3h 3min 경과, CPU 393%, 메모리 386MB. 아직 초기 영점 계산 또는 t100-200 학습 중.
  - precision_filter (PID 40529): 3h 3min 경과, CPU 392%, 메모리 421MB. 5시드×200ep 학습 진행 중.
  - 시스템 로드: 12코어 완전 포화.
**코드 리뷰**: high_height_scaling.py, precision_filter.py 모두 점검 완료. 형상 불일치나 잠재 버그 없음. precision_filter의 `_build_features` 메서드 존재 확인 — 훈련 후 밀집 격자 단계에서 크래시하지 않을 것.
**판단**: 휴식. CPU 포화, 새 결과 없음, 코드 리뷰 완료. 대기가 최적.
**예상 완료**: complex_vector ~2-3h > precision_filter ~3-5h > high_height ~6-10h
**다음**: complex_vector_readout.txt 도착 시 즉시 Phase 1 수행.
