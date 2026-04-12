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
CYCLE_COUNT_FILE="$LOG_DIR/.cycle_count"
MEMORY_DIR="$PROJ_DIR/scripts/memory"
REFLECTION_DIR="$PROJ_DIR/scripts/reflection_logs"
CLAUDE_MEM="$HOME/.claude/projects/-home-k0who029/memory"

mkdir -p "$LOG_DIR" "$MEMORY_DIR/working" "$MEMORY_DIR/episodic" "$MEMORY_DIR/semantic" "$REFLECTION_DIR"

# 사이클 카운터 관리
get_cycle_count() {
    if [ -f "$CYCLE_COUNT_FILE" ]; then
        cat "$CYCLE_COUNT_FILE"
    else
        echo 0
    fi
}

increment_cycle() {
    local count=$(get_cycle_count)
    echo $((count + 1)) > "$CYCLE_COUNT_FILE"
    echo $((count + 1))
}

# 메모리 정리: 3일 이상 된 working 메모리를 episodic으로 이관
cleanup_working_memory() {
    if [ -d "$MEMORY_DIR/working" ]; then
        find "$MEMORY_DIR/working" -name "*.md" -mtime +3 -exec mv {} "$MEMORY_DIR/episodic/" \; 2>/dev/null
    fi
}

run_once() {
    local timestamp=$(date +%Y%m%d_%H%M%S)
    local logfile="$LOG_DIR/run_${timestamp}.log"
    local cycle=$(increment_cycle)

    echo "[$(date)] 자율 연구 코워커 시작 (사이클 #$cycle)" | tee "$logfile"

    # 메모리 정리
    cleanup_working_memory

    # 클로드 코드 메모리 동기화 추적: 실행 전 타임스탬프 기록
    local mem_mtime_before=$(stat -c %Y "$CLAUDE_MEM/MEMORY.md" 2>/dev/null || echo 0)

    # 메타 반성 체크 (10사이클마다)
    if [ $((cycle % 10)) -eq 0 ] && [ -f "$PROJ_DIR/scripts/meta_reflection.py" ]; then
        echo "[$(date)] 메타 반성 실행 (사이클 #$cycle)" | tee -a "$logfile"
        ~/qrop_env/bin/python "$PROJ_DIR/scripts/meta_reflection.py" --cycles 10 >> "$logfile" 2>&1
    fi

    # 문헌 모니터링 체크 (주 1회, 약 336사이클 = 7일 * 48사이클/일)
    if [ $((cycle % 336)) -eq 0 ] && [ -f "$PROJ_DIR/scripts/literature_monitor.py" ]; then
        echo "[$(date)] 문헌 모니터링 실행" | tee -a "$logfile"
        ~/qrop_env/bin/python "$PROJ_DIR/scripts/literature_monitor.py" >> "$logfile" 2>&1
    fi

    cd "$PROJ_DIR" && claude -p "$(cat "$PROMPT_FILE")" \
        --model opus \
        --permission-mode bypassPermissions \
        --max-turns 50 \
        >> "$logfile" 2>&1

    local exit_code=$?
    echo "[$(date)] 코워커 완료 (exit: $exit_code, 사이클 #$cycle)" | tee -a "$logfile"

    # 클로드 코드 메모리 동기화 확인 (로그에 기록)
    if [ -f "$CLAUDE_MEM/MEMORY.md" ]; then
        local mem_mtime_after=$(stat -c %Y "$CLAUDE_MEM/MEMORY.md" 2>/dev/null || echo 0)
        if [ "$mem_mtime_after" -gt "$mem_mtime_before" ]; then
            echo "[$(date)] 클로드 코드 메모리 동기화됨 (사이클 #$cycle)" | tee -a "$logfile"
        else
            echo "[$(date)] 클로드 코드 메모리 변경 없음 (사이클 #$cycle)" | tee -a "$logfile"
        fi
    fi

    # 에러 시 일지에 기록
    if [ $exit_code -ne 0 ]; then
        echo "" >> "$JOURNAL"
        echo "## $(date +%Y-%m-%d\ %H:%M) 에러 (사이클 #$cycle)" >> "$JOURNAL"
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
