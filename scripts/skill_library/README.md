# RDL 스킬 라이브러리

성공한 실험에서 추출된 재사용 가능한 코드 루틴 모음.
Voyager(NVIDIA)의 스킬 라이브러리 패턴에서 영감.

## 사용법

새 실험 작성 전에 이 디렉토리를 확인하여 재사용 가능한 코드를 찾으세요.

```python
# 예시: 스킬 라이브러리에서 임포트
import sys; sys.path.insert(0, 'scripts/skill_library')
from eval_f2_detection import evaluate_f2_detection
```

## 스킬 추가 규칙

1. 실험이 **양성** 판정을 받은 후에만 스킬 추출
2. 각 스킬은 독립 실행 가능한 .py 파일
3. docstring에 사용법, 입출력, 주의사항 명시
4. 의존성은 표준 라이브러리 + torch + numpy + scipy만 허용

## 스킬 목록

| 파일 | 함수 | 설명 | 추출 원본 |
|------|------|------|-----------|
| eval_f2_detection.py | `compute_f2_landscape()` | 전체 데이터에서 \|F₂\| 배열 반환 | s1_geo_high_height, fp_anatomy |
| eval_f2_detection.py | `evaluate_f2_detection()` | 검출 통계 (TP/FP/precision/recall/F1) | 전 실험 공통 |
| train_resonant_model.py | `train_baseline()` | 표준 L_tgt Baseline 훈련 | 전 실험 공통 |
| train_resonant_model.py | `train_s1_integrated()` | S¹ L_geo 통합 훈련 | s1_full_integration, s1_geo_high_height |
| train_resonant_model.py | `GeodesicTargetLoss` | S¹ geodesic 손실 클래스 | s1_full_integration |
