#!/bin/bash
# RDL 연구팀 3-에이전트 시스템 v2
# 보완: 파일 충돌 방지, 순차 실행, lock, 일일 상한, journal 정리
# Usage: ./team_run.sh [--loop SECONDS]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CLAUDE="/home/k0who029/.npm-global/bin/claude"
LOG_DIR="$PROJECT_DIR/outputs/team_logs"
STATE_FILE="$PROJECT_DIR/outputs/.team_state"
LOCK_FILE="$PROJECT_DIR/outputs/.team_lock"
CYCLE_COUNT_FILE="$PROJECT_DIR/outputs/.team_cycle_count"
MAX_CYCLES_PER_DAY=12

mkdir -p "$LOG_DIR"

# ─── Lock 관리 ───
acquire_lock() {
    if [ -f "$LOCK_FILE" ]; then
        local lock_age=$(( $(date +%s) - $(stat -c %Y "$LOCK_FILE" 2>/dev/null || echo 0) ))
        if [ "$lock_age" -gt 3600 ]; then
            echo "[$(date)] 오래된 lock 제거 (${lock_age}초)"
            rm -f "$LOCK_FILE"
        else
            echo "[$(date)] 다른 사이클 실행 중 — 스킵"
            return 1
        fi
    fi
    echo $$ > "$LOCK_FILE"
    return 0
}

release_lock() {
    rm -f "$LOCK_FILE"
}

# ─── 일일 사이클 상한 ───
check_daily_limit() {
    local today=$(date +%Y%m%d)
    local count=0
    if [ -f "$CYCLE_COUNT_FILE" ]; then
        local saved_date=$(head -1 "$CYCLE_COUNT_FILE" | cut -d: -f1)
        if [ "$saved_date" = "$today" ]; then
            count=$(head -1 "$CYCLE_COUNT_FILE" | cut -d: -f2)
        fi
    fi

    if [ "$count" -ge "$MAX_CYCLES_PER_DAY" ]; then
        echo "[$(date)] 일일 상한 도달 ($MAX_CYCLES_PER_DAY회) — 스킵"
        return 1
    fi

    # 카운트 증가
    echo "${today}:$((count + 1))" > "$CYCLE_COUNT_FILE"
    echo "[$(date)] 오늘 사이클: $((count + 1))/$MAX_CYCLES_PER_DAY"
    return 0
}

# ─── 상태 스냅샷 ───
snapshot() {
    local results=$(ls -l "$PROJECT_DIR/results/"*.txt 2>/dev/null | md5sum | cut -d' ' -f1)
    local running=$(ps aux | grep "qrop_env.*python.*scripts/" | grep -v grep | wc -l 2>/dev/null || echo 0)
    local boards=$(cat "$SCRIPT_DIR/board/"*.md 2>/dev/null | md5sum | cut -d' ' -f1)
    echo "${results}_${running}_${boards}"
}

has_changes() {
    local current=$(snapshot)
    local previous=""
    [ -f "$STATE_FILE" ] && previous=$(cat "$STATE_FILE")

    if [ "$current" = "$previous" ]; then
        return 1
    else
        echo "$current" > "$STATE_FILE"
        return 0
    fi
}

has_unreflected() {
    local reflected="$PROJECT_DIR/results/.reflected"
    local all_results=$(ls "$PROJECT_DIR/results/"*.txt 2>/dev/null | sort)
    local reflected_list=""
    [ -f "$reflected" ] && reflected_list=$(sort "$reflected" 2>/dev/null)

    if [ "$all_results" != "$reflected_list" ]; then
        return 0
    fi
    return 1
}

experiment_just_finished() {
    local running=$(ps aux | grep "qrop_env.*python.*scripts/" | grep -v grep | wc -l)
    local prev_running="$PROJECT_DIR/outputs/.prev_running"
    local was_running=0
    [ -f "$prev_running" ] && was_running=$(cat "$prev_running")
    echo "$running" > "$prev_running"

    if [ "$was_running" -gt 0 ] && [ "$running" -eq 0 ]; then
        return 0
    fi
    return 1
}

# ─── 좀비 검사 (토큰 0, bash만) ───
check_zombie() {
    local pids=$(ps aux | grep "qrop_env.*python.*scripts/" | grep -v grep | awk '{print $2}')
    for pid in $pids; do
        local cmdline=$(ps -p "$pid" -o args= 2>/dev/null || true)
        local script=$(echo "$cmdline" | grep -oP 'scripts/\K[^ ]+' || true)
        local logname=$(echo "$script" | sed 's/\.py$//')
        local logfile="/tmp/${logname}.log"

        [ ! -f "$logfile" ] && continue

        if tail -20 "$logfile" | grep -q "Traceback\|TypeError\|RuntimeError\|ValueError\|MemoryError"; then
            echo "[$(date)] 좀비 발견: PID=$pid ($script)"
            tail -3 "$logfile"
            kill "$pid" 2>/dev/null || true
            echo "[$(date)] PID=$pid kill 완료"

            # journal에 기록
            echo -e "\n## $(date '+%Y-%m-%d %H:%M') 좀비 자동 kill\n**PID**: $pid ($script)\n**원인**: $(tail -1 "$logfile")\n" \
                >> "$SCRIPT_DIR/research_journal.md"
        fi
    done
}

# ─── journal/board 크기 관리 ───
trim_files() {
    # research_journal.md: 최근 500줄만 유지
    local journal="$SCRIPT_DIR/research_journal.md"
    if [ -f "$journal" ]; then
        local lines=$(wc -l < "$journal")
        if [ "$lines" -gt 500 ]; then
            echo "[$(date)] journal 정리: ${lines}줄 → 500줄"
            head -500 "$journal" > "${journal}.tmp" && mv "${journal}.tmp" "$journal"
        fi
    fi

    # board 파일: 각 200줄 제한
    for board_file in "$SCRIPT_DIR/board/"*.md; do
        [ ! -f "$board_file" ] && continue
        local lines=$(wc -l < "$board_file")
        if [ "$lines" -gt 200 ]; then
            echo "[$(date)] board 정리: $(basename $board_file) ${lines}줄 → 200줄"
            head -200 "$board_file" > "${board_file}.tmp" && mv "${board_file}.tmp" "$board_file"
        fi
    done

    # 30일 이상 된 로그 삭제
    find "$LOG_DIR" -name "*.log" -mtime +30 -delete 2>/dev/null || true
}

# ─── 트리거 판단 ───
should_run() {
    if has_unreflected; then
        echo "[$(date)] 트리거: 미반영 결과 존재"
        return 0
    fi

    if experiment_just_finished; then
        echo "[$(date)] 트리거: 실험 완료 감지"
        return 0
    fi

    if has_changes; then
        echo "[$(date)] 트리거: 상태 변화 감지"
        return 0
    fi

    local running=$(ps aux | grep "qrop_env.*python.*scripts/" | grep -v grep | wc -l)
    if [ "$running" -eq 0 ]; then
        echo "[$(date)] 트리거: 실행 중 실험 없음 (새 실험 필요)"
        return 0
    fi

    return 1
}

# ─── 팀 사이클 (순차: 수학자 → 실험자+저술가) ───
run_team() {
    local TS=$(date +%Y%m%d_%H%M%S)

    if ! acquire_lock; then
        return 1
    fi

    if ! check_daily_limit; then
        release_lock
        return 1
    fi

    trap 'release_lock' EXIT

    echo "[$(date)] === 팀 사이클 시작 ==="

    # Phase 1: 수학자 (먼저 완료까지 대기)
    echo "[$(date)] [Phase 1] 수학자 시작..."
    $CLAUDE -p "$(cat $SCRIPT_DIR/agent_mathematician.md)" \
        --model opus \
        --permission-mode bypassPermissions \
        --max-turns 15 \
        > "$LOG_DIR/mathematician_${TS}.log" 2>&1
    local MATH_EXIT=$?
    echo "[$(date)] 수학자 완료 (exit: $MATH_EXIT)"

    if [ $MATH_EXIT -ne 0 ]; then
        echo "[$(date)] 수학자 실패 — 실험자/저술가도 중단"
        release_lock
        trap - EXIT
        return 1
    fi

    # Phase 2: 실험자 + 저술가 (수학자 완료 후 병렬)
    echo "[$(date)] [Phase 2] 실험자 + 저술가 시작..."

    $CLAUDE -p "$(cat $SCRIPT_DIR/agent_experimenter.md)" \
        --model opus \
        --permission-mode bypassPermissions \
        --max-turns 20 \
        > "$LOG_DIR/experimenter_${TS}.log" 2>&1 &
    PID_EXP=$!

    $CLAUDE -p "$(cat $SCRIPT_DIR/agent_writer.md)" \
        --model opus \
        --permission-mode bypassPermissions \
        --max-turns 15 \
        > "$LOG_DIR/writer_${TS}.log" 2>&1 &
    PID_WRITE=$!

    echo "[$(date)] 실험자=$PID_EXP, 저술가=$PID_WRITE"

    wait $PID_EXP;   echo "[$(date)] 실험자 완료 (exit: $?)"
    wait $PID_WRITE; echo "[$(date)] 저술가 완료 (exit: $?)"

    # 파일 크기 관리
    trim_files

    echo "[$(date)] === 팀 사이클 완료 ==="

    release_lock
    trap - EXIT
}

# ─── 메인 ───
cd "$PROJECT_DIR"

if [ "$1" = "--loop" ]; then
    INTERVAL=${2:-600}
    echo "[$(date)] 팀 루프 모드: ${INTERVAL}초 간격 체크, 일일 최대 ${MAX_CYCLES_PER_DAY}회"
    while true; do
        # 매번: 좀비 체크 (토큰 0)
        check_zombie

        # 조건 충족 시만: 팀 사이클 (토큰 사용)
        if should_run; then
            run_team || true
        else
            echo "[$(date)] 변화 없음 — 스킵"
        fi

        sleep "$INTERVAL"
    done
else
    check_zombie
    run_team
fi
