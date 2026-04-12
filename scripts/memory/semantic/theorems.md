# RDL 확립된 정리/결과 (의미 메모리)

## 확립된 명제

### Prop 8.5: Bandwidth-Density Bridge Inequality
- **내용**: K ≥ W/(2Δt) — 푸리에 기저 수의 정보론적 하한
- **증거**: 5-시드, ratio 1.42±0.73
- **확립일**: 2026-04-10

### ±π 유일 병목
- **내용**: ±π 불연속이 학습의 유일한 병목
- **증거**: isolation p=1.2e-15, normalized target으로 MSE 갭 99.8% 소멸
- **확립일**: 2026-04-11

### F₂ 블라인드 외삽
- **내용**: 블라인드 영점 예측 recall 99.1±1.7%
- **확립일**: 2026-04-10

### S¹ Geodesic 4배 개선 (3시드 양성 신호, 5시드 확인 대기)
- **내용**: L_geo가 L_tgt 대체 시 검출 5.7→22.7 (4배) — t∈[100,200], 3시드
- **메커니즘**: cos(Δφ) 기반 손실 = S¹ 위의 측지 거리, von Mises NLL과 동치
- **상태**: §3 준수를 위해 "확립"→"3시드 양성 신호"로 강등. s1_geo_high_height (5시드) 결과로 확인 예정
- **관찰일**: 2026-04-12

### 자발적 게이지 고정
- **내용**: r→1은 U(1)^H→U(1) 대칭 파괴, κ≈980 (von Mises)
- **확립일**: 2026-04-12

## 확립된 대응 관계
- |F₂| ↔ Hardy Z-function
- PQO ↔ Maslov index
- K ≥ W/(2Δt) ↔ Nyquist
- ±π 불연속 ↔ branch cut
- L_geo ↔ von Mises NLL
- r→1 ↔ Spontaneous Gauge Fixing

## 음성/중립 결과 (하지 말아야 할 것)
- complex_vector_readout: 중립 — val loss 개선이나 F₂/검출 변화 없음
- precision_filter: 중립 — 사후 필터로 정밀도 교정 불가 (F1=0.430 최대)
- smooth_target_full_loss: 음성 — 5.2% 악화, 시너지 없음
