#!/bin/bash
# =============================================================================
# RDL 자율 연구 코워커 데몬
# =============================================================================
# 사용법:
#   ./scripts/auto_research.sh          # 1회 실행
#   ./scripts/auto_research.sh --loop   # 30분마다 반복
#   ./scripts/auto_research.sh --loop 3600  # 1시간마다
#   ./scripts/auto_research.sh --watch  # 파일 변경 감지 시 실행
# =============================================================================

PROJ_DIR="$HOME/Desktop/gdl_unified"
PROMPT_FILE="$PROJ_DIR/scripts/auto_research_prompt.md"
JOURNAL="$PROJ_DIR/scripts/research_journal.md"
LOG_DIR="$PROJ_DIR/outputs/auto_research_logs"
mkdir -p "$LOG_DIR"

run_once() {
    local timestamp=$(date +%Y%m%d_%H%M%S)
    local logfile="$LOG_DIR/run_${timestamp}.log"

    echo "[$(date)] 자율 연구 코워커 시작" | tee "$logfile"

    cd "$PROJ_DIR" && claude -p "$(cat "$PROMPT_FILE")" \
        --model opus \
        --permission-mode bypassPermissions \
        --max-turns 50 \
        >> "$logfile" 2>&1

    local exit_code=$?
    echo "[$(date)] 코워커 완료 (exit: $exit_code)" | tee -a "$logfile"

    # 에러 시 일지에 기록
    if [ $exit_code -ne 0 ]; then
        echo "" >> "$JOURNAL"
        echo "## $(date +%Y-%m-%d\ %H:%M) 에러" >> "$JOURNAL"
        echo "**에러**: exit code $exit_code. 로그: $logfile" >> "$JOURNAL"
    fi
}

case "${1:-}" in
    --loop)
        INTERVAL=${2:-1800}  # 기본 30분
        echo "=== RDL 코워커 루프 모드 ==="
        echo "간격: ${INTERVAL}초, 프로젝트: $PROJ_DIR"
        echo "일지: $JOURNAL"
        echo "로그: $LOG_DIR"
        echo "중지: Ctrl+C 또는 kill $$"
        echo ""
        while true; do
            run_once
            echo "[$(date)] 다음 실행까지 ${INTERVAL}초 대기..."
            sleep "$INTERVAL"
        done
        ;;
    --watch)
        echo "=== RDL 코워커 감시 모드 ==="
        if ! command -v inotifywait &>/dev/null; then
            echo "inotify-tools 필요: sudo apt install inotify-tools"
            exit 1
        fi
        echo "results/ 감시 중... 새 .txt 파일 생성 시 자동 실행"
        while true; do
            inotifywait -q -e create -e modify "$PROJ_DIR/results/" \
                --include '\.txt$'
            echo "[$(date)] 새 결과 감지!"
            sleep 30  # 파일 쓰기 완료 대기
            run_once
        done
        ;;
    *)
        run_once
        ;;
esac
