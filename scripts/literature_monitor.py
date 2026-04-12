#!/usr/bin/env python3
"""RDL 문헌 자동 모니터링 — Semantic Scholar API 활용

Semantic Scholar API를 통해 RDL 관련 키워드로 최근 논문을 검색하고,
이전 스캔 결과와 비교하여 신규 논문만 리포트한다.

실행: ~/qrop_env/bin/python scripts/literature_monitor.py
"""

import json
import os
import sys
import time
import urllib.request
import urllib.parse
import urllib.error
from datetime import datetime, timedelta
from pathlib import Path

# ─── 설정 ───────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent
MEMORY_DIR = BASE_DIR / "memory" / "semantic"
MEMORY_DIR.mkdir(parents=True, exist_ok=True)

# 검색 키워드 목록
KEYWORDS = [
    "Riemann zeta zeros neural network",
    "phase quantization deep learning",
    "resonant detection zeros",
    "Hardy Z function machine learning",
    "Kuramoto synchronization neural",
    "gauge theory neural network",
    "Berry-Keating conjecture",
]

# 우리 연구와 겹치는 핵심 키워드 (유사 연구 감지용)
OUR_KEYWORDS = {
    "resonant", "phase quantization", "gauge", "zeta zeros",
    "hardy z", "kuramoto", "deep learning", "neural network",
    "zero detection", "riemann", "order parameter",
}

# API 설정
API_URL = "https://api.semanticscholar.org/graph/v1/paper/search"
FIELDS = "title,year,authors,abstract,citationCount"
RATE_LIMIT_SEC = 1.0  # 요청당 대기 시간 (초)


def fetch_papers(query: str, min_year: int) -> list[dict]:
    """Semantic Scholar API로 논문 검색.

    Args:
        query: 검색 키워드
        min_year: 이 연도 이상의 논문만 반환

    Returns:
        논문 딕셔너리 리스트
    """
    params = urllib.parse.urlencode({
        "query": query,
        "fields": FIELDS,
        "limit": 20,
        "year": f"{min_year}-",
    })
    url = f"{API_URL}?{params}"

    req = urllib.request.Request(url)
    req.add_header("User-Agent", "RDL-LiteratureMonitor/1.0")

    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            return data.get("data", [])
    except urllib.error.HTTPError as e:
        print(f"  [오류] HTTP {e.code}: {query[:40]}...")
        return []
    except urllib.error.URLError as e:
        print(f"  [오류] 연결 실패: {e.reason}")
        return []
    except Exception as e:
        print(f"  [오류] 예외: {e}")
        return []


def load_previous_ids() -> set[str]:
    """이전 스캔 파일들에서 이미 수집된 논문 ID 집합을 로드."""
    seen = set()
    for md_file in sorted(MEMORY_DIR.glob("literature_scan_*.md")):
        try:
            text = md_file.read_text(encoding="utf-8")
            # paperId 라인에서 ID 추출
            for line in text.splitlines():
                if line.startswith("- paperId: "):
                    pid = line.split("- paperId: ", 1)[1].strip()
                    if pid:
                        seen.add(pid)
        except Exception:
            continue
    return seen


def check_similarity(title: str, abstract: str) -> bool:
    """우리 연구와 유사한지 키워드 매칭으로 판별.

    제목+초록에 우리 핵심 키워드가 3개 이상 포함되면 유사 연구로 판단.
    """
    text = (title + " " + (abstract or "")).lower()
    matches = sum(1 for kw in OUR_KEYWORDS if kw in text)
    return matches >= 3


def format_authors(authors: list[dict]) -> str:
    """저자 리스트를 문자열로 포맷."""
    if not authors:
        return "N/A"
    names = [a.get("name", "?") for a in authors[:5]]
    suffix = f" 외 {len(authors) - 5}명" if len(authors) > 5 else ""
    return ", ".join(names) + suffix


def run_scan():
    """전체 스캔 실행."""
    today = datetime.now().strftime("%Y-%m-%d")
    print(f"=== RDL 문헌 모니터링 시작 ({today}) ===\n")

    # 6개월 전 기준 연도 계산
    six_months_ago = datetime.now() - timedelta(days=180)
    min_year = six_months_ago.year

    # 이전 스캔에서 이미 본 논문 ID 로드
    previous_ids = load_previous_ids()
    print(f"이전 스캔에서 {len(previous_ids)}개 논문 ID 로드됨\n")

    # 키워드별 검색
    all_papers = {}  # paperId -> paper dict
    new_papers = {}  # 신규 논문만

    for i, keyword in enumerate(KEYWORDS):
        print(f"[{i+1}/{len(KEYWORDS)}] 검색: {keyword}")
        papers = fetch_papers(keyword, min_year)
        print(f"  → {len(papers)}개 결과")

        for p in papers:
            pid = p.get("paperId")
            if not pid:
                continue

            # 6개월 필터 (year 기준 추가 확인)
            year = p.get("year")
            if year and year < min_year:
                continue

            p["_keyword"] = keyword
            all_papers[pid] = p

            if pid not in previous_ids:
                new_papers[pid] = p

        # API rate limit 준수
        if i < len(KEYWORDS) - 1:
            time.sleep(RATE_LIMIT_SEC)

    print(f"\n총 {len(all_papers)}개 논문 수집, 신규 {len(new_papers)}개\n")

    # 유사 연구 체크
    for pid, p in new_papers.items():
        title = p.get("title", "")
        abstract = p.get("abstract", "")
        if check_similarity(title, abstract):
            p["_similar"] = True

    # 결과 마크다운 생성
    lines = [
        f"# RDL 문헌 스캔 — {today}",
        "",
        f"검색 키워드 {len(KEYWORDS)}개, 총 {len(all_papers)}개 수집, "
        f"**신규 {len(new_papers)}개**",
        "",
    ]

    if new_papers:
        lines.append("## 신규 논문")
        lines.append("")

        # 유사 연구 먼저, 그 다음 인용수 내림차순
        sorted_papers = sorted(
            new_papers.values(),
            key=lambda p: (not p.get("_similar", False), -(p.get("citationCount") or 0)),
        )

        for p in sorted_papers:
            tag = "**[주의: 유사 연구]** " if p.get("_similar") else ""
            title = p.get("title", "제목 없음")
            year = p.get("year", "?")
            authors = format_authors(p.get("authors", []))
            citations = p.get("citationCount", 0)
            abstract = p.get("abstract", "")
            keyword = p.get("_keyword", "")
            pid = p.get("paperId", "")

            lines.append(f"### {tag}{title}")
            lines.append(f"- paperId: {pid}")
            lines.append(f"- 연도: {year}")
            lines.append(f"- 저자: {authors}")
            lines.append(f"- 인용수: {citations}")
            lines.append(f"- 검색어: {keyword}")
            if abstract:
                # 초록 첫 200자만 표시
                short = abstract[:200] + ("..." if len(abstract) > 200 else "")
                lines.append(f"- 초록: {short}")
            lines.append("")
    else:
        lines.append("## 신규 논문 없음")
        lines.append("")
        lines.append("이전 스캔 이후 새로운 논문이 발견되지 않았습니다.")
        lines.append("")

    # 통계 요약
    similar_count = sum(1 for p in new_papers.values() if p.get("_similar"))
    lines.append("## 요약 통계")
    lines.append(f"- 검색 키워드 수: {len(KEYWORDS)}")
    lines.append(f"- 총 수집 논문: {len(all_papers)}")
    lines.append(f"- 신규 논문: {len(new_papers)}")
    lines.append(f"- 유사 연구 경고: {similar_count}")
    lines.append(f"- 기준 연도: {min_year}~")
    lines.append("")

    # 파일 저장
    output_path = MEMORY_DIR / f"literature_scan_{today}.md"
    output_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"결과 저장: {output_path}")
    print(f"  신규 논문 {len(new_papers)}개 (유사 연구 경고 {similar_count}개)")


if __name__ == "__main__":
    run_scan()
