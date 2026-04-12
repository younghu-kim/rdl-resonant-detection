#!/usr/bin/env python3
"""가설 Elo 토너먼트 — Google AI Co-Scientist 영감

가설/전략을 Elo 레이팅으로 관리하고, 수동 대결 결과를 입력받아
순위를 업데이트하는 토너먼트 시스템.

사용법:
  python elo_tournament.py add "가설명" "설명"
  python elo_tournament.py match H1 H2 winner   (winner: 1, 2, draw)
  python elo_tournament.py rank
  python elo_tournament.py export

실행: ~/qrop_env/bin/python scripts/elo_tournament.py
"""

import json
import sys
import uuid
from datetime import datetime
from pathlib import Path

# ─── 경로 설정 ──────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent
HYPOTHESES_FILE = BASE_DIR / "memory" / "working" / "hypotheses.json"
EXPORT_DIR = BASE_DIR / "board"
EXPORT_FILE = EXPORT_DIR / "elo_ranking.md"

# 디렉토리 생성
HYPOTHESES_FILE.parent.mkdir(parents=True, exist_ok=True)
EXPORT_DIR.mkdir(parents=True, exist_ok=True)

# ─── Elo 설정 ───────────────────────────────────────────────────────
K_FACTOR = 32  # Elo K 계수
DEFAULT_ELO = 1500  # 초기 Elo


# ─── 데이터 입출력 ──────────────────────────────────────────────────

def load_data() -> dict:
    """가설 데이터 로드. 파일이 없으면 빈 구조 반환."""
    if HYPOTHESES_FILE.exists():
        try:
            text = HYPOTHESES_FILE.read_text(encoding="utf-8")
            return json.loads(text)
        except (json.JSONDecodeError, IOError):
            print("[경고] hypotheses.json 파싱 실패, 빈 데이터로 시작합니다.")
    return {"hypotheses": [], "matches": []}


def save_data(data: dict):
    """가설 데이터를 JSON으로 저장."""
    HYPOTHESES_FILE.write_text(
        json.dumps(data, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


# ─── Elo 계산 ───────────────────────────────────────────────────────

def expected_score(rating_a: float, rating_b: float) -> float:
    """A의 기대 승률 계산 (표준 Elo 공식)."""
    return 1.0 / (1.0 + 10.0 ** ((rating_b - rating_a) / 400.0))


def update_elo(
    elo_a: float, elo_b: float, result: str
) -> tuple[float, float]:
    """Elo 업데이트.

    Args:
        elo_a: 가설 A의 현재 Elo
        elo_b: 가설 B의 현재 Elo
        result: "1" (A 승), "2" (B 승), "draw" (무승부)

    Returns:
        (새 elo_a, 새 elo_b)
    """
    ea = expected_score(elo_a, elo_b)
    eb = 1.0 - ea

    if result == "1":
        sa, sb = 1.0, 0.0
    elif result == "2":
        sa, sb = 0.0, 1.0
    else:  # draw
        sa, sb = 0.5, 0.5

    new_a = elo_a + K_FACTOR * (sa - ea)
    new_b = elo_b + K_FACTOR * (sb - eb)
    return round(new_a, 1), round(new_b, 1)


# ─── 가설 검색 헬퍼 ─────────────────────────────────────────────────

def find_hypothesis(data: dict, identifier: str) -> dict | None:
    """ID 또는 이름(부분 매칭)으로 가설 검색."""
    # 정확한 ID 매칭
    for h in data["hypotheses"]:
        if h["id"] == identifier:
            return h

    # 이름 부분 매칭 (대소문자 무시)
    identifier_lower = identifier.lower()
    candidates = [
        h for h in data["hypotheses"]
        if identifier_lower in h["name"].lower()
    ]
    if len(candidates) == 1:
        return candidates[0]
    elif len(candidates) > 1:
        print(f"[오류] '{identifier}'에 해당하는 가설이 여러 개입니다:")
        for c in candidates:
            print(f"  - {c['id']}: {c['name']}")
        return None

    return None


# ─── 명령어 구현 ────────────────────────────────────────────────────

def cmd_add(name: str, description: str):
    """새 가설 추가."""
    data = load_data()

    # 중복 이름 체크
    for h in data["hypotheses"]:
        if h["name"].lower() == name.lower():
            print(f"[오류] 동일한 이름의 가설이 이미 존재합니다: {h['id']}")
            return

    # 짧은 ID 생성 (H + 6자리)
    short_id = "H" + uuid.uuid4().hex[:6].upper()

    hypothesis = {
        "id": short_id,
        "name": name,
        "description": description,
        "elo": DEFAULT_ELO,
        "matches": 0,
        "wins": 0,
        "losses": 0,
        "draws": 0,
        "created": datetime.now().isoformat(),
    }

    data["hypotheses"].append(hypothesis)
    save_data(data)
    print(f"[추가 완료] {short_id}: {name} (Elo: {DEFAULT_ELO})")


def cmd_match(id_a: str, id_b: str, winner: str):
    """두 가설의 대결 결과를 기록하고 Elo 업데이트."""
    if winner not in ("1", "2", "draw"):
        print("[오류] winner는 1, 2, draw 중 하나여야 합니다.")
        return

    data = load_data()

    hyp_a = find_hypothesis(data, id_a)
    hyp_b = find_hypothesis(data, id_b)

    if hyp_a is None:
        print(f"[오류] 가설을 찾을 수 없습니다: {id_a}")
        return
    if hyp_b is None:
        print(f"[오류] 가설을 찾을 수 없습니다: {id_b}")
        return
    if hyp_a["id"] == hyp_b["id"]:
        print("[오류] 같은 가설끼리는 대결할 수 없습니다.")
        return

    # Elo 업데이트
    old_a, old_b = hyp_a["elo"], hyp_b["elo"]
    new_a, new_b = update_elo(old_a, old_b, winner)

    hyp_a["elo"] = new_a
    hyp_b["elo"] = new_b
    hyp_a["matches"] += 1
    hyp_b["matches"] += 1

    # 승패 기록
    if winner == "1":
        hyp_a["wins"] += 1
        hyp_b["losses"] += 1
        result_str = f"{hyp_a['name']} 승리"
    elif winner == "2":
        hyp_a["losses"] += 1
        hyp_b["wins"] += 1
        result_str = f"{hyp_b['name']} 승리"
    else:
        hyp_a["draws"] += 1
        hyp_b["draws"] += 1
        result_str = "무승부"

    # 대결 기록 저장
    match_record = {
        "date": datetime.now().isoformat(),
        "hyp_a": hyp_a["id"],
        "hyp_b": hyp_b["id"],
        "winner": winner,
        "elo_before": [old_a, old_b],
        "elo_after": [new_a, new_b],
    }
    data.setdefault("matches", []).append(match_record)

    save_data(data)

    print(f"\n  대결: {hyp_a['name']} vs {hyp_b['name']}")
    print(f"  결과: {result_str}")
    print(f"  {hyp_a['name']}: {old_a} → {new_a} ({new_a - old_a:+.1f})")
    print(f"  {hyp_b['name']}: {old_b} → {new_b} ({new_b - old_b:+.1f})")


def cmd_rank():
    """현재 Elo 순위표 출력."""
    data = load_data()

    if not data["hypotheses"]:
        print("[정보] 등록된 가설이 없습니다. 'add' 명령으로 추가하세요.")
        return

    # Elo 내림차순 정렬
    ranked = sorted(data["hypotheses"], key=lambda h: h["elo"], reverse=True)

    print("\n╔════╦════════════════════════════════════╦══════╦═══════╦════════════╗")
    print("║ #  ║ 가설명                             ║  Elo ║ 경기수║ 승/패/무   ║")
    print("╠════╬════════════════════════════════════╬══════╬═══════╬════════════╣")

    for i, h in enumerate(ranked, 1):
        name = h["name"][:32].ljust(32)
        elo = str(int(h["elo"])).rjust(4)
        matches = str(h["matches"]).rjust(4)
        record = f"{h.get('wins',0)}/{h.get('losses',0)}/{h.get('draws',0)}"
        record = record.center(10)
        print(f"║ {i:<2} ║ {name}   ║ {elo} ║ {matches}  ║ {record} ║")

    print("╚════╩════════════════════════════════════╩══════╩═══════╩════════════╝\n")


def cmd_export():
    """순위표를 마크다운으로 내보내기."""
    data = load_data()

    if not data["hypotheses"]:
        print("[정보] 등록된 가설이 없습니다.")
        return

    ranked = sorted(data["hypotheses"], key=lambda h: h["elo"], reverse=True)
    today = datetime.now().strftime("%Y-%m-%d %H:%M")

    lines = [
        f"# Elo 가설 순위표",
        f"",
        f"최종 업데이트: {today}",
        f"",
        f"| 순위 | ID | 가설명 | Elo | 경기수 | 승/패/무 | 설명 |",
        f"|:---:|:---:|:---|:---:|:---:|:---:|:---|",
    ]

    for i, h in enumerate(ranked, 1):
        record = f"{h.get('wins',0)}/{h.get('losses',0)}/{h.get('draws',0)}"
        desc = h["description"][:50] + ("..." if len(h["description"]) > 50 else "")
        lines.append(
            f"| {i} | {h['id']} | {h['name']} | {int(h['elo'])} "
            f"| {h['matches']} | {record} | {desc} |"
        )

    lines.append("")

    # 최근 대결 기록 (최대 10개)
    recent = data.get("matches", [])[-10:]
    if recent:
        lines.append("## 최근 대결 기록")
        lines.append("")
        lines.append("| 날짜 | 가설 A | 가설 B | 결과 | Elo 변동 |")
        lines.append("|:---:|:---:|:---:|:---:|:---:|")

        # ID → 이름 매핑
        id_map = {h["id"]: h["name"] for h in data["hypotheses"]}

        for m in reversed(recent):
            date = m["date"][:10]
            name_a = id_map.get(m["hyp_a"], m["hyp_a"])
            name_b = id_map.get(m["hyp_b"], m["hyp_b"])
            if m["winner"] == "1":
                result = f"**{name_a}** 승"
            elif m["winner"] == "2":
                result = f"**{name_b}** 승"
            else:
                result = "무승부"
            delta_a = m["elo_after"][0] - m["elo_before"][0]
            delta_b = m["elo_after"][1] - m["elo_before"][1]
            elo_str = f"{delta_a:+.0f} / {delta_b:+.0f}"
            lines.append(f"| {date} | {name_a} | {name_b} | {result} | {elo_str} |")

        lines.append("")

    EXPORT_FILE.write_text("\n".join(lines), encoding="utf-8")
    print(f"[내보내기 완료] {EXPORT_FILE}")


def print_usage():
    """사용법 출력."""
    print("""
사용법:
  python elo_tournament.py add "가설명" "설명"
      새 가설을 추가합니다 (초기 Elo: 1500).

  python elo_tournament.py match <ID_A> <ID_B> <winner>
      두 가설의 대결 결과를 기록합니다.
      winner: 1 (A 승), 2 (B 승), draw (무승부)

  python elo_tournament.py rank
      현재 Elo 순위표를 출력합니다.

  python elo_tournament.py export
      순위표를 board/elo_ranking.md에 저장합니다.

예시:
  python elo_tournament.py add "S1 geodesic loss" "L_geo가 L_tgt보다 우수하다"
  python elo_tournament.py add "PGGD optimizer" "PGGD가 Adam보다 수렴이 빠르다"
  python elo_tournament.py match H1A2B3 H4C5D6 1
  python elo_tournament.py rank
""")


# ─── 메인 ───────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print_usage()
        sys.exit(1)

    command = sys.argv[1].lower()

    if command == "add":
        if len(sys.argv) < 4:
            print('[오류] 사용법: python elo_tournament.py add "가설명" "설명"')
            sys.exit(1)
        cmd_add(sys.argv[2], sys.argv[3])

    elif command == "match":
        if len(sys.argv) < 5:
            print("[오류] 사용법: python elo_tournament.py match <ID_A> <ID_B> <winner>")
            sys.exit(1)
        cmd_match(sys.argv[2], sys.argv[3], sys.argv[4])

    elif command == "rank":
        cmd_rank()

    elif command == "export":
        cmd_export()

    else:
        print(f"[오류] 알 수 없는 명령: {command}")
        print_usage()
        sys.exit(1)


if __name__ == "__main__":
    main()
