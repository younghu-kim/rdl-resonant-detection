# 실패 반성문: s1_geo_high_height.py v1

**날짜**: 2026-04-13
**실험**: S¹ Geodesic 고높이 검증 (v1)
**실패 유형**: Type B — 방법론 오류 (코드 자체는 동작하나 평가가 무의미)
**낭비 시간**: ~10시간 (PID 82053, 00:21~10:50)

## 증상
- 모든 시드, 모든 구간에서 `검출=0/0`, `|F₂| ratio=0.0000`
- 모델은 정상 학습 (val loss 수렴), 하지만 평가 결과가 전부 0

## 근본 원인
- `eval_F2()` 함수가 `dataset.is_near_zero` 마스크 사용
- `is_near_zero = xi_amp < (xi_amp.median() * 0.01)`
- t>300에서 xi 진폭이 언더플로 → 모든 값이 ~0 → median도 ~0 → 임계값 이하인 점이 없음
- 결과: `total_zeros=0`, 어떤 검출도 불가능

## 교훈
1. **기존 성공 실험의 평가 방법론을 확인하라**: high_height_scaling_v2는 극소값 탐색+매칭을 사용했으나, 이 스크립트는 is_near_zero를 복사해 사용
2. **고높이(t>300) 실험에서는 항상 극소값 탐색 방식 사용**: is_near_zero는 t<300에서만 유효
3. **초반 출력에서 ratio=0.0000, 검출=0/0이 나오면 즉시 중단**: 정상이라면 최소한 일부 검출이 있어야 함
4. **자가 검증 체크리스트에 "고높이 eval 방식" 항목 추가 필요**

## 수정
- v2: is_near_zero → find_local_minima + match_predictions
- 블라인드 예측 방식 (전반부 훈련, 후반부 예측)
- high_height_scaling_v2와 동일한 검증된 방법론 적용

## 향후 방지
- t>300 실험 작성 시 반드시 `high_height_scaling_v2.py`의 평가 루틴 참조
- 스킬 라이브러리에 고높이 전용 평가 함수 추출 고려
