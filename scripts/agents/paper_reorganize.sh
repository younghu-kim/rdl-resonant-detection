#!/bin/bash
# =============================================================================
# 논문 구조 점검 및 자동 재배치 (24시간 단위)
# =============================================================================
# 사용법:
#   ./paper_reorganize.sh          # 1회 실행
#   ./paper_reorganize.sh daemon   # 24시간마다 자동 실행 (systemd timer)
#   ./paper_reorganize.sh status   # 현재 상태 확인
#
# 기능:
#   1. 카테고리별 결과 분류 점검 (paper_categories.md 기준)
#   2. 분리 트리거 확인 (25p / 8결과 / 3결과)
#   3. Summary Table 정합성 확인
#   4. EN/KO 내용 동기화 확인
#   5. 필요 시 Claude Code로 자동 재배치 실행
# =============================================================================

set -euo pipefail

PROJECT_DIR="$HOME/Desktop/gdl_unified"
PAPER_DIR="$PROJECT_DIR/paper/source"
AGENTS_DIR="$PROJECT_DIR/scripts/agents"
LOG_DIR="$PROJECT_DIR/outputs/paper_reorg_logs"
LOCK_FILE="/tmp/rdl_paper_reorg"

mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log() { echo -e "${BLUE}[$(date +%H:%M:%S)]${NC} $1"; }
ok()  { echo -e "${GREEN}[$(date +%H:%M:%S)] ✓${NC} $1"; }
warn(){ echo -e "${YELLOW}[$(date +%H:%M:%S)] ⚠${NC} $1"; }
err() { echo -e "${RED}[$(date +%H:%M:%S)] ✗${NC} $1"; }

# =============================================================================
# 상태 확인
# =============================================================================
check_status() {
    log "논문 구조 상태 확인..."

    EN_FILE="$PAPER_DIR/unified_master_en.tex"
    KO_FILE="$PAPER_DIR/unified_master_ko.tex"

    # 페이지 수 (마지막 컴파일 기준)
    EN_PAGES=$(grep -o 'Output written on.*(\([0-9]*\) pages' "$PAPER_DIR/unified_master_en.log" 2>/dev/null | grep -o '[0-9]* pages' | head -1 || echo "?")
    KO_PAGES=$(grep -o 'Output written on.*([0-9]* pages' "$PAPER_DIR/unified_master_ko.log" 2>/dev/null | grep -o '[0-9]* pages' | head -1 || echo "?")

    # 섹션 수
    EN_SECTIONS=$(grep -c '\\section' "$EN_FILE" 2>/dev/null || echo 0)
    KO_SECTIONS=$(grep -c '\\section' "$KO_FILE" 2>/dev/null || echo 0)

    # Part 수
    EN_PARTS=$(grep -c '\\part{' "$EN_FILE" 2>/dev/null || echo 0)
    KO_PARTS=$(grep -c '\\part{' "$KO_FILE" 2>/dev/null || echo 0)

    # Summary table 결과 수
    EN_RESULTS=$(grep -c '^[0-9]' "$EN_FILE" 2>/dev/null || echo 0)

    echo ""
    echo "  ┌─────────────────────────────────────┐"
    echo "  │       논문 구조 현황 보고서          │"
    echo "  ├─────────────────────────────────────┤"
    echo "  │ EN: ${EN_PAGES}, ${EN_SECTIONS} sections, ${EN_PARTS} parts"
    echo "  │ KO: ${KO_PAGES}, ${KO_SECTIONS} sections, ${KO_PARTS} parts"
    echo "  │ EN/KO 섹션 차이: $((EN_SECTIONS - KO_SECTIONS))"
    echo "  └─────────────────────────────────────┘"
    echo ""

    # 분리 트리거 확인
    log "분리 트리거 확인..."

    # GL(2) 결과 수
    GL2_COUNT=$(grep -c 'GL(2)\|elliptic\|11a1\|37a1\|Ramanujan.*Delta' "$PAPER_DIR/unified_master_en.tex" 2>/dev/null || echo 0)
    if [ "$GL2_COUNT" -gt 50 ]; then
        warn "GL(2) 관련 언급 ${GL2_COUNT}회 — 분리 검토 필요"
    else
        ok "GL(2) 분리 트리거 미충족"
    fi

    # 페이지 수 기반 트리거
    EN_PAGE_NUM=$(echo "$EN_PAGES" | grep -o '[0-9]*' || echo 0)
    if [ "$EN_PAGE_NUM" -gt 25 ]; then
        warn "EN 본문 ${EN_PAGE_NUM}p > 25p — 분리 검토 (부록 포함 수치)"
    fi
}

# =============================================================================
# 카테고리 정합성 검증
# =============================================================================
check_categories() {
    log "카테고리 정합성 검증..."

    EN_FILE="$PAPER_DIR/unified_master_en.tex"

    # Part IV 내 섹션 순서 확인
    MIDPOINT_LINE=$(grep -n 'Midpoint Curvature Nonlocality' "$EN_FILE" | head -1 | cut -d: -f1)
    CORR_LINE=$(grep -n 'Structural Correspondence Table' "$EN_FILE" | head -1 | cut -d: -f1)
    F2_LINE=$(grep -n 'The \$|\\Ftwo|\$' "$EN_FILE" | head -1 | cut -d: -f1)

    if [ -n "$MIDPOINT_LINE" ] && [ -n "$CORR_LINE" ] && [ -n "$F2_LINE" ]; then
        if [ "$MIDPOINT_LINE" -lt "$CORR_LINE" ] && [ "$CORR_LINE" -lt "$F2_LINE" ]; then
            ok "Part IV 순서 정상: Midpoint($MIDPOINT_LINE) → Corr($CORR_LINE) → F₂($F2_LINE)"
        else
            err "Part IV 순서 이상! Midpoint=$MIDPOINT_LINE, Corr=$CORR_LINE, F₂=$F2_LINE"
            return 1
        fi
    else
        warn "Part IV 섹션 위치 탐지 실패"
    fi

    # 카테고리 주석 존재 확인
    CAT_COMMENTS=$(grep -c '\[Paper [A-D]' "$EN_FILE" 2>/dev/null || echo 0)
    if [ "$CAT_COMMENTS" -ge 4 ]; then
        ok "카테고리 주석 ${CAT_COMMENTS}개 확인"
    else
        warn "카테고리 주석 부족: ${CAT_COMMENTS}개 (최소 4개 필요)"
    fi

    return 0
}

# =============================================================================
# EN/KO 동기화 확인
# =============================================================================
check_sync() {
    log "EN/KO 동기화 확인..."

    EN_FILE="$PAPER_DIR/unified_master_en.tex"
    KO_FILE="$PAPER_DIR/unified_master_ko.tex"

    # Part 이름 비교
    EN_PART4=$(grep '\\part{' "$EN_FILE" | sed -n '4p')
    KO_PART4=$(grep '\\part{' "$KO_FILE" | sed -n '4p')

    log "  EN Part IV: $EN_PART4"
    log "  KO Part IV: $KO_PART4"

    # \section 수 비교
    EN_SEC=$(grep -c '\\section' "$EN_FILE")
    KO_SEC=$(grep -c '\\section' "$KO_FILE")

    if [ "$EN_SEC" -eq "$KO_SEC" ]; then
        ok "섹션 수 일치: EN=$EN_SEC, KO=$KO_SEC"
    else
        warn "섹션 수 불일치: EN=$EN_SEC, KO=$KO_SEC (차이: $((EN_SEC - KO_SEC)))"
    fi

    # Summary table 결과 수 비교
    EN_TAB=$(grep -c '& \\textbf{[ECN]}' "$EN_FILE" 2>/dev/null || echo 0)
    KO_TAB=$(grep -c '& \\textbf{[ECN]}' "$KO_FILE" 2>/dev/null || echo 0)

    if [ "$EN_TAB" -eq "$KO_TAB" ]; then
        ok "Summary table 결과 수 일치: $EN_TAB"
    else
        warn "Summary table 불일치: EN=$EN_TAB, KO=$KO_TAB"
    fi
}

# =============================================================================
# 자동 재배치 실행 (Claude Code 호출)
# =============================================================================
run_reorganize() {
    log "논문 자동 재배치 시작..."

    LOGFILE="$LOG_DIR/reorg_${TIMESTAMP}.log"

    # 잠금 확인
    if [ -f "$LOCK_FILE" ]; then
        warn "이미 재배치 진행 중 (${LOCK_FILE}). 건너뜀."
        return 0
    fi
    touch "$LOCK_FILE"
    trap "rm -f $LOCK_FILE" EXIT

    # 1. 상태 확인
    check_status 2>&1 | tee -a "$LOGFILE"

    # 2. 카테고리 정합성
    if check_categories 2>&1 | tee -a "$LOGFILE"; then
        ok "카테고리 정합성 통과" | tee -a "$LOGFILE"
    else
        err "카테고리 정합성 실패 — Claude Code로 재배치 실행" | tee -a "$LOGFILE"

        # Claude Code 호출하여 재배치
        claude -p "
논문 TeX 파일의 구조가 올바르지 않다. paper_categories.md 기준으로 재배치하라.

1. unified_master_en.tex와 unified_master_ko.tex의 Part IV 내 섹션 순서를 확인
2. 올바른 순서: Midpoint Nonlocality → Structural Correspondence → |F₂| Residual
3. 잘못된 순서면 재배치
4. 재배치 후 pdflatex/xelatex 컴파일
5. PDF 배포 (~/Desktop/수학최종논문/ + paper/)
6. git commit + push

내용 누락 절대 금지.
" --model sonnet 2>&1 | tee -a "$LOGFILE"
    fi

    # 3. EN/KO 동기화
    check_sync 2>&1 | tee -a "$LOGFILE"

    # 4. 분리 트리거 감지 시 알림
    GL2_RESULTS=$(grep -c '\\\\textbf{E}.*GL(2)' "$PAPER_DIR/unified_master_en.tex" 2>/dev/null || echo 0)
    if [ "$GL2_RESULTS" -ge 5 ]; then
        warn "GL(2) 확립 결과 ${GL2_RESULTS}개 — 저장소 분리 시점" | tee -a "$LOGFILE"
    fi

    # 5. 특화 레포 동기화
    log "특화 레포 동기화 실행..." | tee -a "$LOGFILE"
    "$AGENTS_DIR/sync_repos.sh" 2>&1 | tee -a "$LOGFILE"

    ok "논문 재배치 + 동기화 점검 완료. 로그: $LOGFILE"
    rm -f "$LOCK_FILE"
}

# =============================================================================
# Daemon 모드 (systemd timer 설정)
# =============================================================================
setup_daemon() {
    log "24시간 논문 재배치 타이머 설정..."

    SERVICE_NAME="rdl-paper-reorg"
    SERVICE_FILE="$HOME/.config/systemd/user/${SERVICE_NAME}.service"
    TIMER_FILE="$HOME/.config/systemd/user/${SERVICE_NAME}.timer"

    mkdir -p "$HOME/.config/systemd/user"

    # Service 파일
    cat > "$SERVICE_FILE" << SEOF
[Unit]
Description=RDL Paper Reorganization Check

[Service]
Type=oneshot
WorkingDirectory=$PROJECT_DIR
ExecStart=$AGENTS_DIR/paper_reorganize.sh
Environment="PATH=$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin"
StandardOutput=append:$LOG_DIR/daemon.log
StandardError=append:$LOG_DIR/daemon_err.log
SEOF

    # Timer 파일 (24시간마다)
    cat > "$TIMER_FILE" << TEOF
[Unit]
Description=RDL Paper Reorganization Timer (24h)

[Timer]
OnCalendar=*-*-* 03:00:00
Persistent=true

[Install]
WantedBy=timers.target
TEOF

    systemctl --user daemon-reload
    systemctl --user enable --now "${SERVICE_NAME}.timer"

    ok "타이머 설정 완료: 매일 03:00에 실행"
    systemctl --user list-timers | grep "$SERVICE_NAME"
}

# =============================================================================
# 메인
# =============================================================================
case "${1:-}" in
    daemon)
        setup_daemon
        ;;
    status)
        check_status
        check_categories
        check_sync
        ;;
    *)
        run_reorganize
        ;;
esac
