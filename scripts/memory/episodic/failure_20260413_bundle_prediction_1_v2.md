## 실패 반성: bundle_prediction_1_fp_monodromy.py v2 (F₂ 계산 크래시)
**유형**: Type A (구현 버그)
**날짜**: 2026-04-13 ~10:02

**원인**: F₂ 계산에서 L_G 텐서 타입 불일치
- `L_G`는 ComplexDouble, `psi`는 Double
- `torch.complex(L_G - psi, zeros)` → ComplexDouble + Double 연산은 되나, 결과를 torch.complex()에 넘기면 "Expected Float/Double but got ComplexDouble" 에러

**수정**: fp_anatomy.py 패턴 적용
```python
# 잘못된 코드:
rot = torch.complex(torch.cos(phi.detach()), -torch.sin(phi.detach()))
diff = torch.complex(L_G.detach() - psi.detach(), torch.zeros_like(psi.detach()))
f2 = (rot * diff).imag

# 수정된 코드:
phi_real = phi.detach().to(dtype=PrecisionManager.REAL_DTYPE)
rot = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
psi_c = psi.detach().to(dtype=PrecisionManager.COMPLEX_DTYPE)
f2 = (rot * (L_G.detach() - psi_c)).imag.mean(dim=-1)
```

**교훈**:
1. 새 스크립트의 F₂ 계산은 반드시 fp_anatomy.py 패턴을 복사할 것
2. `torch.complex()`에 이미 complex인 텐서를 넘기지 말 것
3. `psi`를 ComplexDouble로 캐스팅 → `L_G - psi_c` 연산이 complex 공간에서 수행
4. `.mean(dim=-1)` 추가 — 채널 차원 축약 필요
