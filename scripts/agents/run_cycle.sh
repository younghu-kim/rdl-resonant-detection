#!/bin/bash
# =============================================================================
# RDL 3단 연구 사이클 오케스트레이터 v3
# =============================================================================
# 사용법:
#   ./run_cycle.sh          # 전체 3단계 실행
#   ./run_cycle.sh 1        # 1단계(수학자)만
#   ./run_cycle.sh 2        # 2단계(설계자)만
#   ./run_cycle.sh 3        # 3단계(검토자)만
#   ./run_cycle.sh daemon   # 30분마다 자동 실행
#   ./run_cycle.sh status   # 현재 상태 확인
#
# v3 변경사항 (구 auto_research.sh + team_run.sh 통합):
#   - 사용자 우선순위 잠금 (/tmp/rdl_user_active)
#   - 사이클 카운터 누적
#   - working 메모리 3일 자동 이관
#   - 메타 반성 (10사이클마다)
#   - 모델 동적 선택 (수학자 지정, Sonnet 실패→Opus 자동 승격)
#   - 검증 게이트 + 반성문 자동 기록
# =============================================================================

set -euo pipefail

PROJECT_DIR="$HOME/Desktop/gdl_unified"
AGENTS_DIR="$PROJECT_DIR/scripts/agents"
LOG_DIR="$PROJECT_DIR/outputs/auto_research_logs"
BOARD_DIR="$PROJECT_DIR/scripts/board"
RESULTS_DIR="$PROJECT_DIR/results"
ANALYSIS_DIR="$PROJECT_DIR/outputs/analysis"
MEMORY_DIR="$PROJECT_DIR/scripts/memory"
EPISODIC_DIR="$MEMORY_DIR/episodic"
JOURNAL="$PROJECT_DIR/scripts/research_journal.md"
CYCLE_COUNT_FILE="$LOG_DIR/.cycle_count"
LOCK_FILE="/tmp/rdl_user_active"

mkdir -p "$LOG_DIR" "$BOARD_DIR" "$RESULTS_DIR" "$EPISODIC_DIR" \
         "$MEMORY_DIR/working" "$MEMORY_DIR/semantic"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# 색상
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

log() { echo -e "${BLUE}[$(date +%H:%M:%S)]${NC} $1"; }
ok()  { echo -e "${GREEN}[$(date +%H:%M:%S)] ✓${NC} $1"; }
err() { echo -e "${RED}[$(date +%H:%M:%S)] ✗${NC} $1"; }
warn(){ echo -e "${YELLOW}[$(date +%H:%M:%S)] !${NC} $1"; }

# ─── 사이클 카운터 ───

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

# ─── 사용자 우선순위 잠금 ───

check_user_lock() {
    if [ -f "$LOCK_FILE" ]; then
        local lock_age
        lock_age=$(( $(date +%s) - $(stat -c %Y "$LOCK_FILE" 2>/dev/null || echo 0) ))
        # 잠금이 2시간 이상 오래되면 stale로 간주
        if [ "$lock_age" -gt 7200 ]; then
            warn "잠금 파일 오래됨 (${lock_age}초) — stale로 간주, 제거"
            rm -f "$LOCK_FILE"
            return 0
        fi
        warn "사용자 활성 세션 감지 (잠금: $LOCK_FILE, ${lock_age}초 전)"
        warn "사이클 건너뜀 — 사용자 오더 우선"
        return 1
    fi
    return 0
}

# 사용자 세션에서 잠금 설정/해제용 함수
lock()   { touch "$LOCK_FILE"; ok "사용자 잠금 설정"; }
unlock() { rm -f "$LOCK_FILE"; ok "사용자 잠금 해제"; }

# ─── 메모리 관리 ───

cleanup_working_memory() {
    if [ -d "$MEMORY_DIR/working" ]; then
        local moved=0
        while IFS= read -r f; do
            mv "$f" "$EPISODIC_DIR/" 2>/dev/null && ((moved++)) || true
        done < <(find "$MEMORY_DIR/working" -name "*.md" -mtime +3 2>/dev/null)
        if [ "$moved" -gt 0 ]; then
            log "working 메모리 ${moved}개 → episodic 이관"
        fi
    fi
}

# ─── 단계 실행 함수 ───

run_stage() {
    local stage=$1
    local name=$2
    local model="${3:-opus}"
    local prompt_file="$AGENTS_DIR/stage${stage}_${name}.md"
    local log_file="$LOG_DIR/${TIMESTAMP}_stage${stage}_${name}.log"

    if [ ! -f "$prompt_file" ]; then
        err "프롬프트 파일 없음: $prompt_file"
        return 1
    fi

    log "Stage $stage: $name 시작 (model: $model)..."

    local user_prompt="프로젝트 디렉토리: $PROJECT_DIR
현재 시각: $(date '+%Y-%m-%d %H:%M')
사이클: #$(get_cycle_count)

위 지침에 따라 이번 사이클을 수행하세요.
보드 파일을 읽고, 상황을 파악하고, 당신의 역할에 맞는 작업을 수행하세요.
작업이 끝나면 보드 파일에 결과를 기록하세요."

    cd "$PROJECT_DIR"
    claude -p "$user_prompt" \
        --system-prompt-file "$prompt_file" \
        --output-format text \
        --model "$model" \
        --allowedTools "Read,Write,Edit,Bash,Glob,Grep" \
        > "$log_file" 2>&1

    local exit_code=$?

    if [ $exit_code -eq 0 ]; then
        ok "Stage $stage: $name 완료 (로그: $log_file)"
    else
        err "Stage $stage: $name 실패 (exit=$exit_code, 로그: $log_file)"
        tail -10 "$log_file" 2>/dev/null || true
    fi

    return $exit_code
}

# ─── 검증 게이트 ───

gate_after_stage1() {
    log "검증 게이트 1: 수학자 지시 확인..."

    if ! grep -q "다음 작업" "$BOARD_DIR/mathematician.md" 2>/dev/null; then
        warn "수학자 지시 없음 — 대기 상태. Stage 2 건너뜀."
        return 1
    fi

    local running
    running=$(ps aux | grep qrop_env | grep -v grep | wc -l)
    if [ "$running" -gt 1 ]; then
        warn "실험 ${running}개 실행 중 — Stage 2 건너뜀."
        return 1
    fi

    ok "게이트 1 통과"
    return 0
}

gate_after_stage2() {
    log "검증 게이트 2: 결과 검증..."
    local stage2_log="$LOG_DIR/${TIMESTAMP}_stage2_executor.log"
    local issues=0

    if grep -qi "Traceback\|SyntaxError\|ImportError" "$stage2_log" 2>/dev/null; then
        warn "Stage 2 로그에 에러 발견"
        ((issues++))
    fi

    local new_results
    new_results=$(find "$RESULTS_DIR" "$ANALYSIS_DIR" -name "*.txt" -newer "$BOARD_DIR/executor.md" 2>/dev/null | wc -l)
    if [ "$new_results" -eq 0 ]; then
        warn "새 결과 파일 없음 (실험 실행 중일 수 있음)"
    fi

    for f in $(find "$RESULTS_DIR" "$ANALYSIS_DIR" -name "*.txt" -newer "$BOARD_DIR/executor.md" 2>/dev/null); do
        if grep -q "영점 0개\|0개 발견\|nan\|NaN\|P=0.0%.*R=0.0%" "$f" 2>/dev/null; then
            err "비정상 결과 감지: $f"
            ((issues++))
        fi
    done

    if [ "$issues" -gt 0 ]; then
        warn "게이트 2: ${issues}개 이슈 발견 — Stage 3에 경고 전달"
        echo "⚠️ 게이트 2 경고: ${issues}개 이슈 ($TIMESTAMP)" >> "$BOARD_DIR/reviewer.md"
    else
        ok "게이트 2 통과"
    fi

    return 0
}

# ─── 실패 시 반성문 ───

write_failure_report() {
    local stage=$1
    local log_file=$2
    local report="$EPISODIC_DIR/failure_$(date +%Y%m%d)_stage${stage}.md"

    cat > "$report" << EOF
# 실패 반성문: Stage $stage ($TIMESTAMP)

## 에러 로그 (마지막 20줄)
\`\`\`
$(tail -20 "$log_file" 2>/dev/null)
\`\`\`

## 원인 분석
(자동 생성 — 다음 사이클에서 수학자가 분석)

## 교훈
- Stage $stage에서 $(date '+%Y-%m-%d %H:%M')에 실패
- 로그: $log_file
EOF

    warn "반성문 작성: $report"
}

# ─── 전체 사이클 ───

run_full_cycle() {
    local cycle=$(increment_cycle)

    log "=========================================="
    log "RDL 연구 사이클 #$cycle 시작 ($TIMESTAMP)"
    log "=========================================="

    # 메모리 정리
    cleanup_working_memory

    # 실행 중인 실험 확인
    local running
    running=$(ps aux | grep qrop_env | grep -v grep | wc -l)
    if [ "$running" -gt 0 ]; then
        warn "실행 중인 실험 ${running}개 감지"
        ps aux | grep qrop_env | grep -v grep | awk '{print "  PID "$2": "$11" "$12}' || true
    fi

    # ── Stage 1: 수학자 (Opus) ──
    local s1_ok=0
    run_stage 1 mathematician opus && s1_ok=1

    if [ "$s1_ok" -eq 0 ]; then
        err "수학자 단계 실패"
        write_failure_report 1 "$LOG_DIR/${TIMESTAMP}_stage1_mathematician.log"
        return 1
    fi

    # 게이트 1
    gate_after_stage1 || {
        log "Stage 2,3 건너뜀 (수학자 대기 지시)"
        return 0
    }

    # ── Stage 2: 설계자 (모델은 수학자 지시 or 기본 sonnet) ──
    local s2_model
    s2_model=$(grep -i '^\*\*모델\*\*:' "$BOARD_DIR/mathematician.md" 2>/dev/null \
        | tail -1 | sed 's/.*: *//' | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')
    if [[ "$s2_model" != "opus" && "$s2_model" != "sonnet" ]]; then
        s2_model="sonnet"
    fi
    log "Stage 2 모델: $s2_model (수학자 지정)"

    local s2_ok=0
    local retry=0
    while [ $retry -lt 2 ] && [ $s2_ok -eq 0 ]; do
        if run_stage 2 executor "$s2_model"; then
            s2_ok=1
        else
            ((retry++))
            if [ $retry -eq 1 ] && [ "$s2_model" = "sonnet" ]; then
                warn "Sonnet 실패 → Opus로 업그레이드 재시도..."
                s2_model="opus"
            else
                warn "Stage 2 재시도 ($retry/2)..."
            fi
        fi
    done

    if [ "$s2_ok" -eq 0 ]; then
        err "설계자 2회 실패"
        write_failure_report 2 "$LOG_DIR/${TIMESTAMP}_stage2_executor.log"
    fi

    # 게이트 2
    gate_after_stage2

    # ── Stage 3: 검토자 (Opus) ──
    run_stage 3 reviewer opus || {
        warn "검토자 단계 실패"
        write_failure_report 3 "$LOG_DIR/${TIMESTAMP}_stage3_reviewer.log"
    }

    # ── 사이클 후 정리 ──

    # PDF 배포
    if [ -f "$PROJECT_DIR/paper/unified_master_en.pdf" ]; then
        mkdir -p "$HOME/Desktop/수학최종논문"
        cp "$PROJECT_DIR/paper/unified_master_en.pdf" "$HOME/Desktop/수학최종논문/" 2>/dev/null && \
            ok "EN PDF 배포 완료" || true
        cp "$PROJECT_DIR/paper/unified_master_ko.pdf" "$HOME/Desktop/수학최종논문/" 2>/dev/null && \
            ok "KO PDF 배포 완료" || true
    fi

    # 메모리 동기화 상태
    "$AGENTS_DIR/sync_memory.sh" status 2>/dev/null || true

    # 일지 기록
    echo "" >> "$JOURNAL"
    echo "## $(date +%Y-%m-%d\ %H:%M) 사이클 #$cycle (자동)" >> "$JOURNAL"
    echo "- Stage 1: $([ $s1_ok -eq 1 ] && echo '완료' || echo '실패')" >> "$JOURNAL"
    echo "- Stage 2: $([ $s2_ok -eq 1 ] && echo "완료 ($s2_model)" || echo '실패')" >> "$JOURNAL"
    echo "- 로그: $LOG_DIR/${TIMESTAMP}_*" >> "$JOURNAL"

    log "=========================================="
    ok "사이클 #$cycle 완료 ($TIMESTAMP)"
    log "=========================================="
}

# ─── 데몬 모드 ───

run_daemon() {
    local interval_min=${1:-30}
    local interval_sec=$((interval_min * 60))

    log "=========================================="
    log "RDL 3단 연구 데몬 v3 시작"
    log "간격: ${interval_min}분, 프로젝트: $PROJECT_DIR"
    log "잠금: $LOCK_FILE (사용자 우선)"
    log "사이클: #$(get_cycle_count)부터 계속"
    log "=========================================="

    while true; do
        if check_user_lock; then
            run_full_cycle || true
        fi
        log "다음 사이클: ${interval_min}분 후"
        sleep "$interval_sec"
    done
}

# ─── 상태 확인 ───

show_status() {
    echo -e "${CYAN}=== RDL 3단 연구 시스템 상태 ===${NC}"
    echo ""
    echo "사이클: #$(get_cycle_count)"
    echo "잠금: $([ -f "$LOCK_FILE" ] && echo "활성 (사용자 우선)" || echo "없음 (데몬 실행 가능)")"
    echo ""

    echo "--- 보드 ---"
    for f in "$BOARD_DIR"/*.md; do
        [ -f "$f" ] && echo "  $(basename $f): $(head -1 "$f" 2>/dev/null) ($(stat -c '%y' "$f" 2>/dev/null | cut -d. -f1))"
    done
    echo ""

    echo "--- 실행 중 실험 ---"
    local running
    running=$(ps aux | grep qrop_env | grep -v grep 2>/dev/null)
    if [ -n "$running" ]; then
        echo "$running" | awk '{print "  PID "$2": "$11" "$12}'
    else
        echo "  없음"
    fi
    echo ""

    echo "--- 최근 로그 ---"
    ls -lt "$LOG_DIR"/*.log 2>/dev/null | head -5 | awk '{print "  "$NF" ("$6" "$7" "$8")"}'
    echo ""

    echo "--- 서비스 ---"
    systemctl --user is-active rdl-cycle.service 2>/dev/null && \
        echo "  rdl-cycle: 활성" || echo "  rdl-cycle: 비활성"
}

# ─── 메인 ───

case "${1:-all}" in
    1)      run_stage 1 mathematician opus ;;
    2)
        s2m=$(grep -i '^\*\*모델\*\*:' "$BOARD_DIR/mathematician.md" 2>/dev/null \
            | tail -1 | sed 's/.*: *//' | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')
        [[ "$s2m" == "opus" || "$s2m" == "sonnet" ]] || s2m="sonnet"
        log "Stage 2 모델: $s2m (수학자 지정)"
        run_stage 2 executor "$s2m"
        ;;
    3)      run_stage 3 reviewer opus ;;
    all)    run_full_cycle ;;
    daemon) run_daemon "${2:-30}" ;;
    status) show_status ;;
    lock)   lock ;;
    unlock) unlock ;;
    *)
        echo "사용법: $0 [1|2|3|all|daemon|status|lock|unlock]"
        echo "  1       — 수학자만 (Opus)"
        echo "  2       — 설계자 (수학자 지정 모델)"
        echo "  3       — 검토자만 (Opus)"
        echo "  all     — 전체 사이클 (기본값)"
        echo "  daemon  — 자동 반복 (기본 30분)"
        echo "  status  — 현재 상태 확인"
        echo "  lock    — 사용자 잠금 (데몬 일시정지)"
        echo "  unlock  — 사용자 잠금 해제"
        exit 1
        ;;
esac
