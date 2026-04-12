#!/usr/bin/env python3
"""메타 반성 — 연구 패턴 분석 및 고차원 통찰 추출

Stanford Generative Agents 영감. research_journal.md를 분석하여
반복 패턴, 생산성, 미해결 문제를 식별하고 권고사항을 생성한다.

실행: ~/qrop_env/bin/python scripts/meta_reflection.py --cycles 10
"""

import argparse
import re
from collections import Counter
from datetime import datetime
from pathlib import Path

# ─── 경로 설정 ──────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent
JOURNAL_FILE = BASE_DIR / "research_journal.md"
REFLECTION_DIR = BASE_DIR / "reflection_logs"
REFLECTION_DIR.mkdir(parents=True, exist_ok=True)

# ─── 판정 키워드 매핑 ───────────────────────────────────────────────
# 저널에서 양성/음성/중립 판정을 식별하는 패턴
POSITIVE_MARKERS = ["양성", "positive", "성공", "확인", "개선", "우수"]
NEGATIVE_MARKERS = ["음성", "negative", "실패", "기각", "악화", "열등"]
NEUTRAL_MARKERS = ["중립", "neutral", "보류", "미정"]

# 실패 유형 분류 키워드
FAILURE_TYPES = {
    "수렴 실패": ["수렴", "발산", "nan", "inf", "overflow"],
    "성능 미달": ["성능", "정확도", "MSE", "loss", "개선 없"],
    "실험 오류": ["오류", "에러", "버그", "크래시", "timeout"],
    "가설 기각": ["기각", "반증", "예측 불일치"],
    "데이터 부족": ["데이터", "샘플", "부족", "불충분"],
}


def parse_journal(text: str) -> list[dict]:
    """research_journal.md를 사이클 단위로 파싱.

    각 사이클은 '## ' 헤더로 시작하며, 본문에서 판정/산출물/실험 정보를 추출.
    """
    cycles = []

    # '## '로 시작하는 섹션 분리
    sections = re.split(r"(?=^## )", text, flags=re.MULTILINE)

    for section in sections:
        section = section.strip()
        if not section or not section.startswith("## "):
            continue

        # 헤더에서 날짜와 제목 추출
        header_line = section.split("\n", 1)[0]
        body = section.split("\n", 1)[1] if "\n" in section else ""

        # 날짜 추출 시도 (YYYY-MM-DD 패턴)
        date_match = re.search(r"(\d{4}-\d{2}-\d{2})", header_line)
        date_str = date_match.group(1) if date_match else "날짜 미상"

        # 판정 분류
        sentiment = classify_sentiment(body)

        # 실패 유형 탐지
        failure_types = detect_failure_types(body)

        # 산출물 카운트 (실행된 항목, 번호 리스트 등)
        outputs = count_outputs(body)

        # 미해결 질문 추출
        questions = extract_questions(body)

        # 연구 방향 추출
        direction = extract_direction(header_line, body)

        cycles.append({
            "header": header_line.replace("## ", ""),
            "date": date_str,
            "body": body,
            "sentiment": sentiment,
            "failure_types": failure_types,
            "outputs": outputs,
            "questions": questions,
            "direction": direction,
        })

    return cycles


def classify_sentiment(text: str) -> str:
    """텍스트의 양성/음성/중립 판정."""
    text_lower = text.lower()

    pos = sum(1 for m in POSITIVE_MARKERS if m in text_lower)
    neg = sum(1 for m in NEGATIVE_MARKERS if m in text_lower)
    neu = sum(1 for m in NEUTRAL_MARKERS if m in text_lower)

    if pos > neg and pos > neu:
        return "양성"
    elif neg > pos and neg > neu:
        return "음성"
    elif neu > 0:
        return "중립"
    else:
        # 기본값: 마커가 없으면 중립
        return "중립"


def detect_failure_types(text: str) -> list[str]:
    """실패 유형 탐지."""
    text_lower = text.lower()
    detected = []
    for ftype, keywords in FAILURE_TYPES.items():
        if any(kw in text_lower for kw in keywords):
            detected.append(ftype)
    return detected


def count_outputs(text: str) -> int:
    """사이클 내 산출물(실행/완료 항목) 수 추정.

    번호가 매겨진 항목(1. 2. 3.)이나 '실행', '완료', '저장' 등을 카운트.
    """
    # 번호 리스트 항목
    numbered = re.findall(r"^\d+\.", text, re.MULTILINE)
    # 실행/완료 키워드
    action_words = len(re.findall(r"(실행|완료|저장|생성|커밋|push|실험)", text))
    return max(len(numbered), action_words)


def extract_questions(text: str) -> list[str]:
    """미해결 질문 추출 — '?' 로 끝나는 문장."""
    questions = []
    for line in text.split("\n"):
        line = line.strip()
        if line.endswith("?") or line.endswith("?\""):
            # 인용부호 등 정리
            clean = re.sub(r"^[-*>]+\s*", "", line)
            if len(clean) > 5:
                questions.append(clean)
    return questions


def extract_direction(header: str, body: str) -> str:
    """연구 방향 키워드 추출."""
    # 헤더에서 괄호 안 내용이나 주요 키워드
    paren = re.search(r"\(([^)]+)\)", header)
    if paren:
        return paren.group(1)

    # 본문에서 '방향', '초점', '주제' 키워드 근처
    for line in body.split("\n")[:5]:
        if any(kw in line for kw in ["방향", "초점", "주제", "목표"]):
            clean = re.sub(r"^[-*>]+\s*", "", line).strip()
            return clean[:60]

    # 헤더 자체 반환 (60자 제한)
    return header[:60]


def analyze_cycles(cycles: list[dict], n: int) -> dict:
    """최근 N개 사이클 분석.

    Returns:
        분석 결과 딕셔너리
    """
    recent = cycles[-n:] if len(cycles) >= n else cycles
    total = len(recent)

    if total == 0:
        return {"error": "분석할 사이클이 없습니다."}

    # ── 1. 양성/음성/중립 비율 ──
    sentiments = Counter(c["sentiment"] for c in recent)
    pos_ratio = sentiments.get("양성", 0) / total * 100
    neg_ratio = sentiments.get("음성", 0) / total * 100
    neu_ratio = sentiments.get("중립", 0) / total * 100

    # ── 2. 반복 실패 패턴 (같은 유형 3회 이상) ──
    failure_counter = Counter()
    for c in recent:
        for ft in c["failure_types"]:
            failure_counter[ft] += 1
    repeated_failures = {k: v for k, v in failure_counter.items() if v >= 3}

    # ── 3. 가장 생산적인 방향 ──
    direction_output = {}
    direction_count = Counter()
    for c in recent:
        d = c["direction"]
        direction_count[d] += 1
        direction_output[d] = direction_output.get(d, 0) + c["outputs"]

    # 산출물 기준 상위 방향
    productive_dirs = sorted(
        direction_output.items(), key=lambda x: x[1], reverse=True
    )[:5]

    # ── 4. 오래된 미해결 질문 ──
    all_questions = []
    for c in recent:
        for q in c["questions"]:
            all_questions.append({"question": q, "date": c["date"]})

    # ── 5. 평균 산출물 수 ──
    total_outputs = sum(c["outputs"] for c in recent)
    avg_outputs = total_outputs / total if total > 0 else 0

    return {
        "total_cycles": total,
        "sentiments": dict(sentiments),
        "pos_ratio": pos_ratio,
        "neg_ratio": neg_ratio,
        "neu_ratio": neu_ratio,
        "repeated_failures": repeated_failures,
        "productive_dirs": productive_dirs,
        "questions": all_questions,
        "avg_outputs": avg_outputs,
        "total_outputs": total_outputs,
    }


def generate_report(analysis: dict) -> str:
    """분석 결과를 마크다운 보고서로 생성."""
    today = datetime.now().strftime("%Y-%m-%d")

    if "error" in analysis:
        return f"# 메타 반성 — {today}\n\n{analysis['error']}\n"

    lines = [
        f"# 메타 반성 — {today}",
        "",
        f"분석 대상: 최근 **{analysis['total_cycles']}개** 사이클",
        "",
    ]

    # ── 패턴 분석 ──
    lines.append("## 패턴 분석")
    lines.append("")

    if analysis["repeated_failures"]:
        lines.append("### 반복 실패 패턴 (3회 이상)")
        lines.append("")
        for ftype, count in sorted(
            analysis["repeated_failures"].items(), key=lambda x: -x[1]
        ):
            severity = "높음" if count >= 5 else "보통"
            lines.append(f"- **{ftype}**: {count}회 반복 (심각도: {severity})")
        lines.append("")
        lines.append("> 반복 실패는 접근 방식의 근본적 변경이 필요할 수 있음을 시사합니다.")
        lines.append("")
    else:
        lines.append("반복 실패 패턴이 탐지되지 않았습니다. (양호)")
        lines.append("")

    # ── 생산성 평가 ──
    lines.append("## 생산성 평가")
    lines.append("")
    lines.append(
        f"- 양성/음성/중립 비율: "
        f"**{analysis['pos_ratio']:.0f}%** / "
        f"**{analysis['neg_ratio']:.0f}%** / "
        f"**{analysis['neu_ratio']:.0f}%**"
    )
    lines.append(
        f"  - 양성: {analysis['sentiments'].get('양성', 0)}건, "
        f"음성: {analysis['sentiments'].get('음성', 0)}건, "
        f"중립: {analysis['sentiments'].get('중립', 0)}건"
    )
    lines.append(
        f"- 평균 사이클당 산출물: **{analysis['avg_outputs']:.1f}**개 "
        f"(총 {analysis['total_outputs']}개)"
    )
    lines.append("")

    if analysis["productive_dirs"]:
        lines.append("### 가장 생산적인 연구 방향 (산출물 기준)")
        lines.append("")
        for i, (direction, outputs) in enumerate(analysis["productive_dirs"], 1):
            lines.append(f"{i}. **{direction}** — 산출물 {outputs}개")
        lines.append("")

    # ── 방향 재검토 ──
    lines.append("## 방향 재검토")
    lines.append("")

    if analysis["questions"]:
        lines.append("### 미해결 질문")
        lines.append("")
        for q_info in analysis["questions"]:
            lines.append(f"- [{q_info['date']}] {q_info['question']}")
        lines.append("")

        # 오래된 질문 경고 (3일 이상)
        try:
            now = datetime.now()
            old_questions = []
            for q_info in analysis["questions"]:
                try:
                    q_date = datetime.strptime(q_info["date"], "%Y-%m-%d")
                    age = (now - q_date).days
                    if age >= 3:
                        old_questions.append((q_info, age))
                except ValueError:
                    continue

            if old_questions:
                lines.append("### 오래된 미해결 질문 경고")
                lines.append("")
                for q_info, age in old_questions:
                    lines.append(
                        f"- **{age}일 경과**: {q_info['question']}"
                    )
                lines.append("")
                lines.append(
                    "> 오래된 질문은 답을 구하거나, "
                    "현재 시점에서 더 이상 유효하지 않으면 폐기하세요."
                )
                lines.append("")
        except Exception:
            pass
    else:
        lines.append("미해결 질문이 발견되지 않았습니다.")
        lines.append("")

    # ── 권고사항 ──
    lines.append("## 권고사항")
    lines.append("")

    recommendations = []

    # 음성 비율이 높으면 경고
    if analysis["neg_ratio"] > 50:
        recommendations.append(
            "음성 결과 비율이 50%를 초과합니다. "
            "현재 접근 방식을 재검토하거나 새로운 방향을 탐색하세요."
        )

    # 반복 실패가 있으면 전환 권고
    if analysis["repeated_failures"]:
        top_failure = max(
            analysis["repeated_failures"].items(), key=lambda x: x[1]
        )
        recommendations.append(
            f"'{top_failure[0]}' 패턴이 {top_failure[1]}회 반복되고 있습니다. "
            f"근본 원인을 분석하고 다른 전략을 시도하세요."
        )

    # 산출물이 적으면 효율성 경고
    if analysis["avg_outputs"] < 2:
        recommendations.append(
            "사이클당 산출물이 적습니다. "
            "실험 설계를 단순화하거나 병렬 실행을 고려하세요."
        )

    # 양성 비율이 높으면 격려
    if analysis["pos_ratio"] > 60:
        recommendations.append(
            "양성 결과 비율이 높습니다. 현재 방향을 유지하되, "
            "확증 편향에 주의하세요."
        )

    # 중립이 많으면 결정력 부족 경고
    if analysis["neu_ratio"] > 60:
        recommendations.append(
            "중립 결과가 과반입니다. 실험 설계를 더 명확한 "
            "양성/음성 판정이 나오도록 개선하세요."
        )

    if not recommendations:
        recommendations.append("특별한 이상 패턴이 발견되지 않았습니다. 현재 방향을 유지하세요.")

    for i, rec in enumerate(recommendations, 1):
        lines.append(f"{i}. {rec}")
    lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="RDL 메타 반성 — 연구 패턴 분석 및 고차원 통찰 추출"
    )
    parser.add_argument(
        "--cycles", type=int, default=10,
        help="분석할 최근 사이클 수 (기본값: 10)"
    )
    args = parser.parse_args()

    # 저널 파일 읽기
    if not JOURNAL_FILE.exists():
        print(f"[오류] 연구 일지를 찾을 수 없습니다: {JOURNAL_FILE}")
        print("  research_journal.md가 scripts/ 디렉토리에 있어야 합니다.")
        return

    print(f"=== RDL 메타 반성 (최근 {args.cycles}개 사이클) ===\n")

    text = JOURNAL_FILE.read_text(encoding="utf-8")
    cycles = parse_journal(text)
    print(f"저널에서 {len(cycles)}개 사이클 파싱됨")

    # 분석 실행
    analysis = analyze_cycles(cycles, args.cycles)

    # 보고서 생성
    report = generate_report(analysis)

    # 파일 저장
    today = datetime.now().strftime("%Y-%m-%d")
    output_path = REFLECTION_DIR / f"reflection_{today}.md"
    output_path.write_text(report, encoding="utf-8")

    # 콘솔 출력
    print(report)
    print(f"\n저장 완료: {output_path}")


if __name__ == "__main__":
    main()
