#!/usr/bin/env python3
"""
#106 arXiv 감사 FAIL 수정 스크립트
- FAIL 1: EN/KO Summary Table에 #67, #68 행 추가
- FAIL 2: EN Discussion 열거에 #67, #68, #76, #77 추가 + 번호 재조정
- FAIL 2: KO Discussion 열거에 #67, #68 추가 + 번호 재조정
- WARNING 1: B-12 상태 모순 제거 (EN/KO)
"""

import re

EN_PATH = "/home/k0who029/Desktop/gdl_unified/paper/source/unified_master_en.tex"
KO_PATH = "/home/k0who029/Desktop/gdl_unified/paper/source/unified_master_ko.tex"

def fix_en():
    with open(EN_PATH, 'r', encoding='utf-8') as f:
        content = f.read()

    original = content  # for diff reporting

    # =========================================================================
    # 1. EN Summary Table: Add #67 and #68 between #66 and #69
    # =========================================================================
    # Find the midrule after #66 row and before #69 row
    old_table = (
        "   & GL(4) $\\mathrm{sym}^3(11\\mathrm{a}1)$ \\\\\n"
        "\\midrule\n"
        "69 & GL(4) $\\xi$-bundle"
    )
    new_table = (
        "   & GL(4) $\\mathrm{sym}^3(11\\mathrm{a}1)$ \\\\\n"
        "\\midrule\n"
        "67 & GL(3) $\\mathrm{sym}^2(11\\mathrm{a}1)$ blind zero verification (21 FP $\\to$ TP) & \\textbf{E} &\n"
        "   21 extra peaks ($t{>}17.69$, beyond LMFDB) verified via AFE ($\\mathrm{dps}{=}80$,\n"
        "   $N_{\\mathrm{coeff}}{=}200$, 60-node GH quadrature);\n"
        "   median $|\\Lambda(1/2{+}it)|{=}2.5{\\times}10^{-27}$;\n"
        "   corrected: $P{=}R{=}F_1{=}1.000$; $\\xi$-bundle as zero-finding machine\n"
        "   (see \\S\\ref{subsec:gl3_blind})\n"
        "   & GL(3) sym$^2$(11a1) \\\\\n"
        "\\midrule\n"
        "68 & GL(3) $\\kappa$-conductor scaling (6 curves, sym$^2$) & \\textbf{E} &\n"
        "   log-log regression: $\\alpha{=}0.003$, $R^2{=}0.0002$ (no scaling);\n"
        "   $\\kappa_{\\mathrm{near}}{\\approx}1125$ (CV${=}0.15\\%$, 5 curves);\n"
        "   conductor-independent universal constant\n"
        "   (see \\S\\ref{subsec:gl3_kappa_scaling})\n"
        "   & GL(3) sym$^2$ \\\\\n"
        "\\midrule\n"
        "69 & GL(4) $\\xi$-bundle"
    )
    assert old_table in content, "ERROR: EN Summary Table target not found"
    content = content.replace(old_table, new_table, 1)
    print("[EN] Summary Table: #67 and #68 rows added")

    # =========================================================================
    # 2. EN Discussion: Add (xiii)/#67, (xiv)/#68, renumber, add #76/#77
    # =========================================================================
    # The discussion is on one giant line. Replace the key substring.
    # Pattern: from "result~\#66), and (xiii)~" to end of "(xxi)~...\#81...degree_univ})."

    old_disc = (
        "result~\\#66), and (xiii)~GL(4) $\\xi$-bundle $\\kappa_\\sigma$ ($\\sigma$-direction) "
        "$A(t_0)$ measurement with B-23 resolved ($\\kappa_\\sigma \\neq \\kappa_t$) and B-12 "
        "identified as normalization-dependent open question (result~\\#69, "
        "Remark~\\ref{rem:xi_bundle_b23}), (xiv)~GL(5) extension to $\\mathrm{sym}^4(11\\mathrm{a}1)$ "
        "($d{=}5$, $w{=}1$, $N{=}14641$) with PARI-verified functional equation to $331$-digit "
        "precision, $\\sigma$-uniqueness PASS ($\\sigma{=}2.5$, ${\\approx}10000\\times$ ratio), "
        "and $\\xi$-bundle $A = 14.88$ (per-zero CV $< 0.2\\%$, motivic normalization; "
        "results~\\#70, \\#71, Remark~\\ref{rem:gl5_sym4}), and (xv)~GL(3) $\\xi$-bundle $A(t_0)$ "
        "measurement at motivic centre $\\sigma = k/2 = 1.5$ with $\\sigma$-uniqueness 6/6 PASS "
        "($\\kappa \\approx 10^6$ at centre, ratio $> 8000\\times$) for both sym$^2$(11a1) "
        "($A{=}9.04$) and sym$^2$(37a1) ($A{=}14.51$), with $k$ auto-detected via $L[4]$ "
        "(B-26 resolved; results~\\#87, \\#88), completing GL(1)--GL(5) coverage, and (xvi)~Hadamard "
        "decomposition of $A(t_0)$ via global zero distribution: $A(t_0) = B(t_0)^2 + 2H_1(t_0)$ "
        "with $B$, $H_1$ expressed as explicit paired sums over all non-trivial zeros "
        "(Proposition~\\ref{prop:hadamard_A_decomp}); B-20 resolved: $A(t_0)$ variability "
        "analytically explained by near-zero $H_1$ contributions; 13/13 numerically verified "
        "($\\mathrm{err}_{\\mathrm{corr}}{<}0.003\\%$; result~\\#74, Remark~\\ref{rem:hadamard_A_numerical}), "
        "and (xvii)~degree universality of Proposition~\\ref{prop:hadamard_A_decomp} extended to "
        "GL(2) $L(s,\\Delta)$: 7/10 zeros below $1\\%$ error ($10/10$ below $2\\%$), $D{=}0$ "
        "confirmed ($|\\mathrm{Re}(c_0)|_{\\mathrm{max}} = 3.3\\times10^{-5}$), B-20 pattern "
        "reproduced; precision decrease consistent with Weyl asymptotics (result~\\#75, "
        "Remark~\\ref{rem:hadamard_gl2_numerical}), and (xviii)~quantification of the sign-change "
        "$\\sigma$-uniqueness boundary across degrees 1--5: at $d \\geq 4$, the sign-change count "
        "$S(\\sigma)$ saturates ($S(\\sigma) = \\text{const}$ for all tested $\\sigma$), rendering "
        "the PASS/FAIL distinction meaningless; GL(4) sym$^3$(11a1) achieves $S = 45$ with ratio "
        "$= 1.000$ (7/7~$\\sigma$), GL(5) sym$^4$(11a1) achieves $S \\in \\{59,60\\}$ with ratio "
        "$\\leq 1.017$ (7/7~$\\sigma$); saturation mechanism: uniform zero spacing $\\to$ sinusoidal "
        "Hardy $Z$ $\\to$ sign-change count $\\sigma$-invariant; boundary B-05 resolved: "
        "$\\sigma$-uniqueness valid for $d \\leq 3$ only (result~\\#78), (xix)~GL(3) $\\sigma$-sweep "
        "PARI re-measurement for sym$^2$(11a1): $\\Delta t{=}0.1$, 457 points, $S{\\in}\\{63,64\\}$ "
        "(7 offsets), $A{=}0.016$ (flat), confirming prior AFE asymmetry $A{=}0.222$ as "
        "coarse-sampling artifact (result~\\#79, Remark~\\ref{rem:sigma_saturation}), and (xx)~GL(3) "
        "Hadamard high-precision verification with $N{=}1189$ zeros ($T{\\leq}500$): $10/10{<}2\\%$, "
        "$6/10{<}1\\%$ relative error ($A_{\\text{true}}$ basis); Re$(c_0)$ Richardson correction "
        "method established ($A_{\\text{true}} = A_{\\text{dir}} - 2\\,\\text{Re}(c_0)/\\delta$, "
        "correcting $\\sim 2$--$5\\%$ numerical asymmetry); B-20 pattern reproduced "
        "($A_{\\max}{=}20.93$ at $t_8{=}12.07$ with 5 nearby zeros, $A_{\\min}{=}6.56$ at "
        "$t_2{=}4.73$); degree universality confirmed: $d{=}1$ (13/13, $<0.003\\%$), $d{=}2$ "
        "(7/10, $<1.3\\%$), $d{=}3$ (6/10, $<1\\%$, $10/10{<}2\\%$) (result~\\#80, "
        "Remark~\\ref{rem:hadamard_gl3_numerical}), and (xxi)~GL(4) $L(s,\\mathrm{sym}^3\\Delta)$ "
        "Hadamard $\\xi$-bundle $A{=}B^2{+}2H_1$: $7/10$ below $5\\%$, $8/10$ below $10\\%$, "
        "16~zeros ($t_0 \\in [0,25]$); degree universality $d{=}1$--$4$ confirmed "
        "(result~\\#81, Remark~\\ref{rem:hadamard_degree_univ})."
    )

    new_disc = (
        "result~\\#66), (xiii)~GL(3) sym$^2$(11a1) blind zero verification: 21 apparent FP "
        "reclassified as genuine zeros via AFE ($\\mathrm{dps}{=}80$, $N_{\\mathrm{coeff}}{=}200$, "
        "60-node GH quadrature); corrected $P{=}R{=}F_1{=}1.000$; $\\xi$-bundle independently "
        "discovers GL(3) zeros beyond LMFDB coverage (result~\\#67, see~\\S\\ref{subsec:gl3_blind}), "
        "(xiv)~GL(3) $\\kappa$-conductor scaling over six sym$^2$ $L$-functions: log-log regression "
        "$\\alpha{=}0.003$, $R^2{=}0.0002$ (no scaling); $\\kappa_{\\mathrm{near}}{\\approx}1125$ "
        "(CV${=}0.15\\%$, 5 curves) conductor-independent universal constant "
        "(result~\\#68, see~\\S\\ref{subsec:gl3_kappa_scaling}), and (xv)~GL(4) $\\xi$-bundle "
        "$\\kappa_\\sigma$ ($\\sigma$-direction) $A(t_0)$ measurement with B-23 resolved "
        "($\\kappa_\\sigma \\neq \\kappa_t$) and B-12 identified as normalization-dependent open "
        "question (result~\\#69, Remark~\\ref{rem:xi_bundle_b23}), (xvi)~GL(5) extension to "
        "$\\mathrm{sym}^4(11\\mathrm{a}1)$ ($d{=}5$, $w{=}1$, $N{=}14641$) with PARI-verified "
        "functional equation to $331$-digit precision, $\\sigma$-uniqueness PASS ($\\sigma{=}2.5$, "
        "${\\approx}10000\\times$ ratio), and $\\xi$-bundle $A = 14.88$ (per-zero CV $< 0.2\\%$, "
        "motivic normalization; results~\\#70, \\#71, Remark~\\ref{rem:gl5_sym4}), and (xvii)~GL(3) "
        "$\\xi$-bundle $A(t_0)$ measurement at motivic centre $\\sigma = k/2 = 1.5$ with "
        "$\\sigma$-uniqueness 6/6 PASS ($\\kappa \\approx 10^6$ at centre, ratio $> 8000\\times$) "
        "for both sym$^2$(11a1) ($A{=}9.04$) and sym$^2$(37a1) ($A{=}14.51$), with $k$ "
        "auto-detected via $L[4]$ (B-26 resolved; results~\\#87, \\#88), completing GL(1)--GL(5) "
        "coverage, and (xviii)~Hadamard decomposition of $A(t_0)$ via global zero distribution: "
        "$A(t_0) = B(t_0)^2 + 2H_1(t_0)$ with $B$, $H_1$ expressed as explicit paired sums over "
        "all non-trivial zeros (Proposition~\\ref{prop:hadamard_A_decomp}); B-20 resolved: $A(t_0)$ "
        "variability analytically explained by near-zero $H_1$ contributions; 13/13 numerically "
        "verified ($\\mathrm{err}_{\\mathrm{corr}}{<}0.003\\%$; result~\\#74, "
        "Remark~\\ref{rem:hadamard_A_numerical}), and (xix)~degree universality of "
        "Proposition~\\ref{prop:hadamard_A_decomp} extended to GL(2) $L(s,\\Delta)$: 7/10 zeros "
        "below $1\\%$ error ($10/10$ below $2\\%$), $D{=}0$ confirmed "
        "($|\\mathrm{Re}(c_0)|_{\\mathrm{max}} = 3.3\\times10^{-5}$), B-20 pattern reproduced; "
        "precision decrease consistent with Weyl asymptotics (result~\\#75, "
        "Remark~\\ref{rem:hadamard_gl2_numerical}), and (xx)~GL(3) Hadamard convergence limit: "
        "raw partial-sum mean error $22.6\\%$, $N{\\approx}15{,}000$ zeros required for ${<}2\\%$ "
        "(result~\\#76); convergence rate $\\alpha(d)$ measurement: GL(1) $\\alpha{=}0.558$, "
        "GL(2) $\\alpha{\\approx}0.36$, GL(3) $\\alpha{\\approx}0.5$; EM tail accuracy is "
        "differentiating factor (result~\\#77, Remark~\\ref{rem:hadamard_convergence_rate}), "
        "and (xxi)~quantification of the sign-change $\\sigma$-uniqueness boundary across "
        "degrees 1--5: at $d \\geq 4$, the sign-change count $S(\\sigma)$ saturates "
        "($S(\\sigma) = \\text{const}$ for all tested $\\sigma$), rendering the PASS/FAIL "
        "distinction meaningless; GL(4) sym$^3$(11a1) achieves $S = 45$ with ratio $= 1.000$ "
        "(7/7~$\\sigma$), GL(5) sym$^4$(11a1) achieves $S \\in \\{59,60\\}$ with ratio "
        "$\\leq 1.017$ (7/7~$\\sigma$); saturation mechanism: uniform zero spacing $\\to$ "
        "sinusoidal Hardy $Z$ $\\to$ sign-change count $\\sigma$-invariant; boundary B-05 "
        "resolved: $\\sigma$-uniqueness valid for $d \\leq 3$ only (result~\\#78), (xxii)~GL(3) "
        "$\\sigma$-sweep PARI re-measurement for sym$^2$(11a1): $\\Delta t{=}0.1$, 457 points, "
        "$S{\\in}\\{63,64\\}$ (7 offsets), $A{=}0.016$ (flat), confirming prior AFE asymmetry "
        "$A{=}0.222$ as coarse-sampling artifact (result~\\#79, Remark~\\ref{rem:sigma_saturation}), "
        "and (xxiii)~GL(3) Hadamard high-precision verification with $N{=}1189$ zeros "
        "($T{\\leq}500$): $10/10{<}2\\%$, $6/10{<}1\\%$ relative error ($A_{\\text{true}}$ basis); "
        "Re$(c_0)$ Richardson correction method established ($A_{\\text{true}} = A_{\\text{dir}} - "
        "2\\,\\text{Re}(c_0)/\\delta$, correcting $\\sim 2$--$5\\%$ numerical asymmetry); B-20 "
        "pattern reproduced ($A_{\\max}{=}20.93$ at $t_8{=}12.07$ with 5 nearby zeros, "
        "$A_{\\min}{=}6.56$ at $t_2{=}4.73$); degree universality confirmed: $d{=}1$ "
        "(13/13, $<0.003\\%$), $d{=}2$ (7/10, $<1.3\\%$), $d{=}3$ (6/10, $<1\\%$, "
        "$10/10{<}2\\%$) (result~\\#80, Remark~\\ref{rem:hadamard_gl3_numerical}), and "
        "(xxiv)~GL(4) $L(s,\\mathrm{sym}^3\\Delta)$ Hadamard $\\xi$-bundle $A{=}B^2{+}2H_1$: "
        "$7/10$ below $5\\%$, $8/10$ below $10\\%$, 16~zeros ($t_0 \\in [0,25]$); degree "
        "universality $d{=}1$--$4$ confirmed "
        "(result~\\#81, Remark~\\ref{rem:hadamard_degree_univ})."
    )

    assert old_disc in content, "ERROR: EN Discussion target not found"
    content = content.replace(old_disc, new_disc, 1)
    print("[EN] Discussion: (xiii)/#67, (xiv)/#68 added; #76/#77 added as (xx); renumbered to (xxiv)")

    # =========================================================================
    # 3. EN B-12: Fix "fully resolved" contradiction
    # =========================================================================
    # Line 7399: paragraph heading
    old_b12_heading = (
        "\\textbf{Theoretical basis for degree dependence (boundary B-12 resolved).}\\;"
    )
    new_b12_heading = (
        "\\textbf{Theoretical basis for degree dependence (boundary B-12 partially resolved).}\\;"
    )
    assert old_b12_heading in content, "ERROR: EN B-12 heading not found"
    content = content.replace(old_b12_heading, new_b12_heading, 1)
    print("[EN] B-12 heading: 'resolved' -> 'partially resolved'")

    # Line 7422: conclusion sentence
    old_b12_concl = (
        "\\emph{Boundary B-12 is therefore fully resolved}: $\\kappa_{\\mathrm{near}}(d)$ measures"
    )
    new_b12_concl = (
        "\\emph{Boundary B-12 is therefore partially resolved}: $\\kappa_{\\mathrm{near}}(d)$ "
        "degree-dependence is analytically explained (B-12a, resolved); cross-degree $A(t_0)$ "
        "normalization-dependence remains open (B-12b, open). Specifically, $\\kappa_{\\mathrm{near}}(d)$ measures"
    )
    assert old_b12_concl in content, "ERROR: EN B-12 conclusion not found"
    content = content.replace(old_b12_concl, new_b12_concl, 1)
    print("[EN] B-12 conclusion: 'fully resolved' -> 'partially resolved' + B-12a/B-12b note")

    with open(EN_PATH, 'w', encoding='utf-8') as f:
        f.write(content)
    print("[EN] File written.")


def fix_ko():
    with open(KO_PATH, 'r', encoding='utf-8') as f:
        content = f.read()

    # =========================================================================
    # 1. KO Summary Table: Add #67 and #68 between #66 and #69
    # =========================================================================
    old_table = (
        "   & GL(4) $\\mathrm{sym}^3(11\\mathrm{a}1)$ \\\\\n"
        "\\midrule\n"
        "69 & GL(4) $\\xi$-다발"
    )
    new_table = (
        "   & GL(4) $\\mathrm{sym}^3(11\\mathrm{a}1)$ \\\\\n"
        "\\midrule\n"
        "67 & GL(3) $\\mathrm{sym}^2(11\\mathrm{a}1)$ 블라인드 영점 검증 (21 FP $\\to$ TP) & \\textbf{E} &\n"
        "   LMFDB 범위 초과 ($t{>}17.69$) 추가 피크 21개를 AFE ($\\mathrm{dps}{=}80$,\n"
        "   $N_{\\mathrm{coeff}}{=}200$, 60-node GH 구적)로 검증;\n"
        "   중앙값 $|\\Lambda(1/2{+}it)|{=}2.5{\\times}10^{-27}$;\n"
        "   교정 후 $P{=}R{=}F_1{=}1.000$; $\\xi$-다발이 LMFDB 미등재 GL(3) 영점을 독립 발견\n"
        "   (\\S\\ref{subsec:gl3_blind} 참조)\n"
        "   & GL(3) sym$^2$(11a1) \\\\\n"
        "\\midrule\n"
        "68 & GL(3) $\\kappa$-도체 스케일링 (6곡선, sym$^2$) & \\textbf{E} &\n"
        "   로그-로그 회귀: $\\alpha{=}0.003$, $R^2{=}0.0002$ (스케일링 없음);\n"
        "   $\\kappa_{\\mathrm{near}}{\\approx}1125$ (CV${=}0.15\\%$, 5곡선);\n"
        "   도체 독립 보편상수\n"
        "   (\\S\\ref{subsec:gl3_kappa_scaling} 참조)\n"
        "   & GL(3) sym$^2$ \\\\\n"
        "\\midrule\n"
        "69 & GL(4) $\\xi$-다발"
    )
    assert old_table in content, "ERROR: KO Summary Table target not found"
    content = content.replace(old_table, new_table, 1)
    print("[KO] Summary Table: #67 and #68 rows added")

    # =========================================================================
    # 2. KO Discussion: Add (xiii)/#67 and (xiv)/#68 + renumber (xiii)→(xv)...(xxii)→(xxiv)
    # =========================================================================
    old_disc = (
        "결과~\\#66), (xiii)~GL(4) $\\xi$-다발 $\\kappa_\\sigma$ ($\\sigma$-방향) $A(t_0)$ 측정, "
        "B-23 해소 ($\\kappa_\\sigma \\neq \\kappa_t$), B-12는 정규화-의존 미결 문제로 확인 "
        "(결과~\\#69, 비고~\\ref{rem:xi_bundle_b23}), (xiv)~GL(5) 확장: "
        "$\\mathrm{sym}^4(11\\mathrm{a}1)$ ($d{=}5$, $w{=}1$, $N{=}14641$) PARI 함수방정식 "
        "$331$자리 검증, $\\sigma$-유일성 통과 ($\\sigma{=}2.5$, ${\\approx}10000\\times$ 비율), "
        "$\\xi$-다발 $A = 14.88$ (영점내 CV $< 0.2\\%$, motivic 정규화; 결과~\\#70, \\#71, "
        "비고~\\ref{rem:gl5_sym4}), (xv)~GL(3) $\\xi$-다발 $A(t_0)$ 측정 "
        "(motivic center $\\sigma = k/2 = 1.5$), $\\sigma$-유일성 6/6 PASS "
        "($\\kappa \\approx 10^6$, 비율 $> 8000\\times$), sym$^2$(11a1) $A{=}9.04$, "
        "sym$^2$(37a1) $A{=}14.51$, $k$ 자동 감지 (B-26 해결; 결과~\\#87, \\#88), "
        "GL(1)--GL(5) 완전 커버리지 달성, (xvi)~전역 영점 분포를 통한 $A(t_0)$의 Hadamard "
        "분해: $A(t_0) = B(t_0)^2 + 2H_1(t_0)$, $B$, $H_1$은 모든 비자명 영점에 대한 명시적 "
        "쌍합산 (명제~\\ref{prop:hadamard_A_decomp}); B-20 해결: 근접 영점 $H_1$ 기여로 "
        "$A(t_0)$ 변동 해석적 설명; 13/13 수치 검증 ($\\mathrm{err}_{\\mathrm{corr}}{<}0.003\\%$; "
        "결과~\\#74, 비고~\\ref{rem:hadamard_A_numerical}), (xvii)~GL(2) $L(s,\\Delta)$에서 "
        "명제~\\ref{prop:hadamard_A_decomp} 차수 보편성 확장: 7/10 영점 $1\\%$ 이하 오차 "
        "($10/10$은 $2\\%$ 이하), $D{=}0$ 확인 "
        "($|\\mathrm{Re}(c_0)|_{\\mathrm{max}} = 3.3\\times10^{-5}$), B-20 패턴 재현; "
        "정밀도 저하는 Weyl 점근 법칙과 일치 (결과~\\#75, 비고~\\ref{rem:hadamard_gl2_numerical}), "
        "(xviii)~GL(3) Hadamard 수렴 한계 규명: raw 부분합 수렴 속도 $\\alpha \\approx 0.5$ "
        "(degree 보편), Weyl 밀도 과보정으로 GL(3)에서 $N{\\approx}15{,}000$ 필요 (결과~\\#76), "
        "멱법칙 피팅으로 degree-보편 수렴 속도 확인: EM 꼬리 보정이 차별화 요인 "
        "(결과~\\#77, 비고~\\ref{rem:hadamard_convergence_rate}), (xix)~차수 1--5에 걸친 "
        "부호변환 $\\sigma$-유일성 경계 정량화: $d \\geq 4$에서 부호변환 수 $S(\\sigma)$가 "
        "포화 ($S(\\sigma) = \\text{const}$, 모든 시험 $\\sigma$에서), PASS/FAIL 구분 무의미; "
        "GL(4) sym$^3$(11a1)은 $S = 45$ (ratio $= 1.000$, 7/7~$\\sigma$), GL(5) "
        "sym$^4$(11a1)은 $S \\in \\{59,60\\}$ (ratio $\\leq 1.017$, 7/7~$\\sigma$); "
        "포화 메커니즘: 균일 영점 간격 $\\to$ 정현파적 Hardy $Z$ $\\to$ 부호변환 수 $\\sigma$-불변; "
        "경계 B-05 해결: $\\sigma$-유일성은 $d \\leq 3$에서만 유효 (결과~\\#78), (xx)~GL(3) "
        "sym$^2$(11a1) PARI $\\sigma$-스윕 재측정: $\\Delta t{=}0.1$, 457점, "
        "$S{\\in}\\{63,64\\}$ (7 오프셋), $A{=}0.016$ (평탄), 기존 AFE 비대칭 $A{=}0.222$가 "
        "조밀도 부족 아티팩트임을 확인 (결과~\\#79, 비고~\\ref{rem:sigma_saturation}), (xxi)~GL(3) "
        "Hadamard 고정밀 검증 ($N{=}1189$ 영점, $T{\\leq}500$): $10/10{<}2\\%$, $6/10{<}1\\%$ "
        "상대 오차 ($A_{\\text{true}}$ 기준); Re$(c_0)$ Richardson 보정 방법론 확립 "
        "($A_{\\text{true}} = A_{\\text{dir}} - 2\\,\\text{Re}(c_0)/\\delta$, ${\\sim}2$--$5\\%$ "
        "수치 비대칭 보정); B-20 패턴 재현 ($A_{\\max}{=}20.93$, $t_8{=}12.07$에서 근접 5영점; "
        "$A_{\\min}{=}6.56$, $t_2{=}4.73$); 차수 보편성 확립: $d{=}1$ (13/13, ${<}0.003\\%$), "
        "$d{=}2$ (7/10, ${<}1.3\\%$), $d{=}3$ (6/10, ${<}1\\%$, $10/10{<}2\\%$) "
        "(결과~\\#80, 비고~\\ref{rem:hadamard_gl3_numerical}), (xxii)~GL(4) "
        "$L(s,\\mathrm{sym}^3\\Delta)$ Hadamard $\\xi$-다발 $A{=}B^2{+}2H_1$: $7/10$은 "
        "$5\\%$ 이하, $8/10$은 $10\\%$ 이하, 16개 영점 ($t_0 \\in [0,25]$); "
        "차수 보편성 $d{=}1$--$4$ 확립 (결과~\\#81, 비고~\\ref{rem:hadamard_degree_univ})."
    )

    new_disc = (
        "결과~\\#66), (xiii)~GL(3) sym$^2$(11a1) 블라인드 영점 검증: 겉보기 FP 21개를 AFE "
        "($\\mathrm{dps}{=}80$, $N_{\\mathrm{coeff}}{=}200$)로 진짜 영점으로 재분류; "
        "교정 후 $P{=}R{=}F_1{=}1.000$; $\\xi$-다발이 LMFDB 범위 밖 GL(3) 영점을 독립 발견 "
        "(결과~\\#67, \\S\\ref{subsec:gl3_blind} 참조), "
        "(xiv)~GL(3) $\\kappa$-도체 스케일링 (6개 sym$^2$ $L$-함수): 로그-로그 회귀 "
        "$\\alpha{=}0.003$, $R^2{=}0.0002$ (스케일링 없음); $\\kappa_{\\mathrm{near}}{\\approx}1125$ "
        "(CV${=}0.15\\%$, 5곡선) 도체 독립 보편상수 "
        "(결과~\\#68, \\S\\ref{subsec:gl3_kappa_scaling} 참조), "
        "(xv)~GL(4) $\\xi$-다발 $\\kappa_\\sigma$ ($\\sigma$-방향) $A(t_0)$ 측정, "
        "B-23 해소 ($\\kappa_\\sigma \\neq \\kappa_t$), B-12는 정규화-의존 미결 문제로 확인 "
        "(결과~\\#69, 비고~\\ref{rem:xi_bundle_b23}), (xvi)~GL(5) 확장: "
        "$\\mathrm{sym}^4(11\\mathrm{a}1)$ ($d{=}5$, $w{=}1$, $N{=}14641$) PARI 함수방정식 "
        "$331$자리 검증, $\\sigma$-유일성 통과 ($\\sigma{=}2.5$, ${\\approx}10000\\times$ 비율), "
        "$\\xi$-다발 $A = 14.88$ (영점내 CV $< 0.2\\%$, motivic 정규화; 결과~\\#70, \\#71, "
        "비고~\\ref{rem:gl5_sym4}), (xvii)~GL(3) $\\xi$-다발 $A(t_0)$ 측정 "
        "(motivic center $\\sigma = k/2 = 1.5$), $\\sigma$-유일성 6/6 PASS "
        "($\\kappa \\approx 10^6$, 비율 $> 8000\\times$), sym$^2$(11a1) $A{=}9.04$, "
        "sym$^2$(37a1) $A{=}14.51$, $k$ 자동 감지 (B-26 해결; 결과~\\#87, \\#88), "
        "GL(1)--GL(5) 완전 커버리지 달성, (xviii)~전역 영점 분포를 통한 $A(t_0)$의 Hadamard "
        "분해: $A(t_0) = B(t_0)^2 + 2H_1(t_0)$, $B$, $H_1$은 모든 비자명 영점에 대한 명시적 "
        "쌍합산 (명제~\\ref{prop:hadamard_A_decomp}); B-20 해결: 근접 영점 $H_1$ 기여로 "
        "$A(t_0)$ 변동 해석적 설명; 13/13 수치 검증 ($\\mathrm{err}_{\\mathrm{corr}}{<}0.003\\%$; "
        "결과~\\#74, 비고~\\ref{rem:hadamard_A_numerical}), (xix)~GL(2) $L(s,\\Delta)$에서 "
        "명제~\\ref{prop:hadamard_A_decomp} 차수 보편성 확장: 7/10 영점 $1\\%$ 이하 오차 "
        "($10/10$은 $2\\%$ 이하), $D{=}0$ 확인 "
        "($|\\mathrm{Re}(c_0)|_{\\mathrm{max}} = 3.3\\times10^{-5}$), B-20 패턴 재현; "
        "정밀도 저하는 Weyl 점근 법칙과 일치 (결과~\\#75, 비고~\\ref{rem:hadamard_gl2_numerical}), "
        "(xx)~GL(3) Hadamard 수렴 한계 규명: raw 부분합 수렴 속도 $\\alpha \\approx 0.5$ "
        "(degree 보편), Weyl 밀도 과보정으로 GL(3)에서 $N{\\approx}15{,}000$ 필요 (결과~\\#76), "
        "멱법칙 피팅으로 degree-보편 수렴 속도 확인: EM 꼬리 보정이 차별화 요인 "
        "(결과~\\#77, 비고~\\ref{rem:hadamard_convergence_rate}), (xxi)~차수 1--5에 걸친 "
        "부호변환 $\\sigma$-유일성 경계 정량화: $d \\geq 4$에서 부호변환 수 $S(\\sigma)$가 "
        "포화 ($S(\\sigma) = \\text{const}$, 모든 시험 $\\sigma$에서), PASS/FAIL 구분 무의미; "
        "GL(4) sym$^3$(11a1)은 $S = 45$ (ratio $= 1.000$, 7/7~$\\sigma$), GL(5) "
        "sym$^4$(11a1)은 $S \\in \\{59,60\\}$ (ratio $\\leq 1.017$, 7/7~$\\sigma$); "
        "포화 메커니즘: 균일 영점 간격 $\\to$ 정현파적 Hardy $Z$ $\\to$ 부호변환 수 $\\sigma$-불변; "
        "경계 B-05 해결: $\\sigma$-유일성은 $d \\leq 3$에서만 유효 (결과~\\#78), (xxii)~GL(3) "
        "sym$^2$(11a1) PARI $\\sigma$-스윕 재측정: $\\Delta t{=}0.1$, 457점, "
        "$S{\\in}\\{63,64\\}$ (7 오프셋), $A{=}0.016$ (평탄), 기존 AFE 비대칭 $A{=}0.222$가 "
        "조밀도 부족 아티팩트임을 확인 (결과~\\#79, 비고~\\ref{rem:sigma_saturation}), (xxiii)~GL(3) "
        "Hadamard 고정밀 검증 ($N{=}1189$ 영점, $T{\\leq}500$): $10/10{<}2\\%$, $6/10{<}1\\%$ "
        "상대 오차 ($A_{\\text{true}}$ 기준); Re$(c_0)$ Richardson 보정 방법론 확립 "
        "($A_{\\text{true}} = A_{\\text{dir}} - 2\\,\\text{Re}(c_0)/\\delta$, ${\\sim}2$--$5\\%$ "
        "수치 비대칭 보정); B-20 패턴 재현 ($A_{\\max}{=}20.93$, $t_8{=}12.07$에서 근접 5영점; "
        "$A_{\\min}{=}6.56$, $t_2{=}4.73$); 차수 보편성 확립: $d{=}1$ (13/13, ${<}0.003\\%$), "
        "$d{=}2$ (7/10, ${<}1.3\\%$), $d{=}3$ (6/10, ${<}1\\%$, $10/10{<}2\\%$) "
        "(결과~\\#80, 비고~\\ref{rem:hadamard_gl3_numerical}), (xxiv)~GL(4) "
        "$L(s,\\mathrm{sym}^3\\Delta)$ Hadamard $\\xi$-다발 $A{=}B^2{+}2H_1$: $7/10$은 "
        "$5\\%$ 이하, $8/10$은 $10\\%$ 이하, 16개 영점 ($t_0 \\in [0,25]$); "
        "차수 보편성 $d{=}1$--$4$ 확립 (결과~\\#81, 비고~\\ref{rem:hadamard_degree_univ})."
    )

    assert old_disc in content, "ERROR: KO Discussion target not found"
    content = content.replace(old_disc, new_disc, 1)
    print("[KO] Discussion: (xiii)/#67, (xiv)/#68 added; renumbered to (xxiv)")

    # =========================================================================
    # 3. KO B-12: Fix "완전히 해소" contradiction (line 5697, 5717)
    # =========================================================================
    old_b12_heading = (
        "\\textbf{차수 의존성의 이론적 근거 (경계 B-12 해결).}\\;"
    )
    new_b12_heading = (
        "\\textbf{차수 의존성의 이론적 근거 (경계 B-12 부분 해결).}\\;"
    )
    assert old_b12_heading in content, "ERROR: KO B-12 heading not found"
    content = content.replace(old_b12_heading, new_b12_heading, 1)
    print("[KO] B-12 heading: '해결' -> '부분 해결'")

    old_b12_concl = (
        "\\emph{경계 B-12는 완전히 해소된다}: $\\kappa_{\\mathrm{near}}(d)$는 차수만이 아니라"
    )
    new_b12_concl = (
        "\\emph{경계 B-12는 부분적으로 해소된다}: $\\kappa_{\\mathrm{near}}(d)$ 차수 의존성은 "
        "해석적으로 설명됨 (B-12a, 해결); 교차-차수 $A(t_0)$ 정규화 의존성은 미결 (B-12b, 미결). "
        "구체적으로, $\\kappa_{\\mathrm{near}}(d)$는 차수만이 아니라"
    )
    assert old_b12_concl in content, "ERROR: KO B-12 conclusion not found"
    content = content.replace(old_b12_concl, new_b12_concl, 1)
    print("[KO] B-12 conclusion: '완전히 해소' -> '부분적으로 해소' + B-12a/B-12b 명시")

    with open(KO_PATH, 'w', encoding='utf-8') as f:
        f.write(content)
    print("[KO] File written.")


if __name__ == "__main__":
    print("=== #106 arXiv FAIL 수정 시작 ===")
    fix_en()
    fix_ko()
    print("=== 모든 수정 완료 ===")
