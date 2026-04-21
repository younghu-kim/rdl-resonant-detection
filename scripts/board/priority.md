# 수학자 우선순위 명령 보드

> 이 파일은 수학자(Stage 1)만 쓸 수 있습니다.
> run_cycle.sh가 매 사이클 시작 시 이 파일을 읽고 명령을 실행합니다.

## 사용 가능한 명령

```
KILL_PID:12345              — 특정 PID 즉시 종료
PAUSE_EXPLORATION:true      — 모든 탐사 프로세스 일시 중지 (SIGSTOP)
PAUSE_EXPLORATION:false     — 탐사 프로세스 재개 (SIGCONT)
PRIORITY:high               — 다음 실험 긴급 표시
```

## 판단 기준

| 상황 | 조치 |
|------|------|
| 탐사 3구간+ 안정, 수학 과제 대기 | KILL_PID |
| 긴급 수학 과제, 탐사 미완 | PAUSE_EXPLORATION:true |
| CPU 유휴, 수학 과제 없음 | 탐사 계속 허용 |
| 6시간+ 장기 실행, 새 정보 없음 | KILL_PID |

## 현재 명령

PRIORITY:high
# 사이클 #233 [2026-04-21 16:13]
# B-34 ★★★★ 돌파급: Slope Universality self-dual→일반 확장. FE+Schwarz 증명 발견.
# #216 Paper 2 Theorem 일반화 반영 지시 (opus). CPU 유휴. 즉시 실행.
