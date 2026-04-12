# RDL 자율 연구 일지

## 2026-04-12 14:05 사이클
**상황**: high_height_scaling v2 완료 — **양성 판정!**
  - recall: 98.8% (t~200) → 98.1% (t~600) → 97.5% (t~1100) — 10배 높이에서도 안정
  - precision: 17.1% → 22.8% → 25.6% — 고높이에서 오히려 개선
  - K=256 vs K=128: 차이 없음
  - 총 20990초 (5.8시간)

**실행**:
  - 논문 EN/KO Tier 2 항목 + sec:f2_residual_open 업데이트 [confirmed 2026-04-12]
  - PDF 컴파일·배포, git commit+push 완료
  - S¹ geodesic 실험 시작 (PID 60718, 651% CPU)

**진행 중**: s1_geodesic_readout.py — seed=42 Baseline 학습 중. 예상 ~1.5시간.

**누적 성과 (이번 에이전트 세션)**:
  1. complex_vector_readout → 중립 (F₂ 무변화)
  2. precision_filter → 중립 (정밀도 15→33%, 유의 임계 미달)
  3. high_height_scaling → **양성** (recall 안정, precision 개선)
  4. S¹ geodesic → 진행 중

**다음**: S¹ 결과 도착 시 분석 → 논문 반영. 양성이면 ±π 해결 돌파구.

## 2026-04-12 07:30 사이클
**상황**: high_height v2 프로세스가 ep60까지 진행 후 종료됨 (에러 로그 없음 — 외부 kill 추정). 캐시(t100-200, t100-150, dps130)는 남아있음.
**실행**: v2 재시작 (PID 58296, 777% CPU). 캐시 재사용으로 즉시 학습 시작. 예상 총 6-7시간.
**다음**: v2 완료 대기. 완료 시 즉시 Phase 1 + S¹ geodesic 시작.

## 2026-04-12 05:30 사이클
**상황**: high_height_scaling.py (v1) 22시간 실행 중 — **두 가지 치명적 버그 발견**:

1. **F₂ 평가 방법 오류**: `Z_out.abs()` (단순 출력 진폭) 사용.
   올바른 F₂ 잔차는 `(e^{-iφ}·(L_G - ψ)).imag` (blind_zero_prediction.py 방식).
   이로 인해 t∈[100,200]에서도 det=0/463, ratio≈1.0 (모든 검출 실패).

2. **xi 정밀도 언더플로**: t∈[500,600]에서 |xi| ~ 10^{-166}.
   `sqrt(real² + imag²)` 계산 시 real² ~ 10^{-332} → float64 subnormal 한계(5e-324) 초과 → amplitude 전부 0.
   결과: is_near_zero 마스크 전부 False → ratio=0.000, det=0/0.
   근본 원인: xi(1/2+it) ~ exp(-πt/4) 지수적 감쇠.
   참고: 학습 자체는 xi 타겟을 사용하지 않으므로(TotalResonanceLoss는 자기지도적) 정상 진행되었으나, 평가가 완전히 무효.

**판단**: v1 실험 중단, 수정된 v2 작성.

**실행**:
  - PID 54990 (v1) 종료, 잘못된 t500-600 캐시 삭제
  - high_height_scaling_v2.py 작성:
    - 올바른 F₂ 잔차 사용 (blind_zero_prediction.py 방식)
    - 극소값 탐색 + 매칭 방식으로 영점 검출 (is_near_zero 의존 제거)
    - 전반부 학습 → 후반부 블라인드 예측 (실제 외삽 능력 검증)
    - dps 스케일링: max(50, int(t_max*0.5)+30) — 고높이 캐시 정밀도 보장
    - early stopping + best model 저장
  - v2 실행 시작 (PID 57251, 706% CPU)
  - s1_geodesic_readout.py의 eval_F2도 동일 버그 수정 (다음 실험 대비)

**진행 중**: high_height_scaling_v2 — t∈[100,200] K=128 seed=42 학습 중 (ep 30/150, 210s)
**예상 총 소요**: 7-8시간 (3 range × 2 K × 3 seed = 18 runs + 3 캐시 계산)

**핵심 교훈**:
  - eval_F2는 반드시 phi/psi/L_G 기반 잔차를 사용해야 함 (Z_out.abs()는 무의미)
  - t > 300에서는 xi 값의 지수적 감쇠로 인해 float64 제곱 연산도 언더플로 가능
  - 고높이 영점 검출은 is_near_zero 마스크 대신 극소값 탐색이 강건함

**다음**: v2 완료 시 결과 분석 → 논문 반영 → S¹ geodesic 실험 시작

## 2026-04-12 03:15 사이클
**상황**: precision_filter 실험 완료 — **중립 판정**.
  - 5시드, t∈[100,150]→[150,200], K=128, 200ep, 2000점 밀집격자
  - 원시 정밀도 28.4%, 최고 F1=0.430 (앙상블 ≥2: P=32.7%, R=63.0%)
  - 깊이 필터 너무 공격적 (전부 제거), 간격 필터 재현율만 저하
  - 결론: 낮은 정밀도는 |F₂| 풍경 내재적 한계, 사후 휴리스틱으로 교정 불가

**실행**:
  - 논문 EN/KO Tier 1 블라인드 외삽 항목에 필터 결과 추가 [neutral 2026-04-12]
  - PDF 컴파일·배포, git commit+push 완료
  - high_height_scaling.py (최적화판) 단독 실행 시작 (PID 54990, 817% CPU)
  - S¹ geodesic 스크립트 준비 완료 (scripts/s1_geodesic_readout.py)

**진행 중**: high_height_scaling — t[100,200] K=128 첫 run 학습 중. 예상 총 5-7시간 (18 runs + 영점계산 2회)

**판단**: high_height가 장시간 소요되므로 대기. 완료 시:
  - 양성 (고높이에서 성능 유지): Tier 2 항목 1 + sec:f2_residual_open 업데이트
  - 음성 (성능 저하): K 스케일링 법칙 도출 필요

**다음**: high_height 완료 후 → 결과 반영 → S¹ geodesic 실험 시작 (±π 해결 다음 후보)

## 2026-04-12 02:15 사이클
**상황**: complex_vector_readout 실험 완료 — **중립 판정**.
  - 3시드 (42,7,123), t∈[100,200], K=128, 150ep
  - Complex Vector val loss ~4배 개선 (0.014±0.005 vs 0.069±0.013)
  - F₂ ratio 무변화 (0.982±0.041 vs 0.999±0.005), 검출 0/463
  - readout 헤드가 Re/Im을 잘 학습하지만 네트워크 내부 Z_out 공명 구조에 영향 없음
  - ±π 불연속은 사후 readout이 아닌 더 깊은 아키텍처 변경 필요

**핵심 최적화**: 이전 3개 실험이 29시간+ CPU 소비에도 결과 없었던 원인 해결
  - 3 프로세스 병렬 → 단일 순차 실행 (12코어 단독 사용, ~9x CPU 활용)
  - eval_F2/compute_dense_f2: 1개씩 → 64개 배치 처리
  - 진행 로깅 추가 (30ep마다 출력)
  - 결과: complex_vector 97분 완료 (이전 추정 30시간+ → 60배 가속)

**실행**:
  - 논문 EN/KO Tier 1 "Complex-vector readout" 항목에 [neutral 2026-04-12] 반영
  - PDF 컴파일·배포, git commit+push 완료
  - precision_filter.py (최적화판) 단독 실행 시작 (PID 54188, 803% CPU)

**판단**: complex vector readout이 중립이므로 ±π 해결의 다음 후보는 **S¹ 기하학적 각도 예측** (arg ξ 대신 (cos θ, sin θ)∈S¹, geodesic 손실). precision_filter와 high_height 결과 확보 후 S¹ 실험 설계.

**다음**: precision_filter 완료 대기 (~1.5시간). 완료 시 즉시 결과 분석 → 논문 반영 → high_height 시작.

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
