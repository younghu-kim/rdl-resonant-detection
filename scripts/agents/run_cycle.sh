#!/bin/bash
# =============================================================================
# RDL 3단 연구 사이클 오케스트레이터 v3
# =============================================================================
# 사용법:
#   ./run_cycle.sh          # 전체 3단계 실행
#   ./run_cycle.sh 1        # 1단계(수학자)만
#   ./run_cycle.sh 2        # 2단계(설계자)만
#   ./run_cycle.sh 3        # 3단계(검토자)만
#   ./run_cycle.sh daemon   # 60분마다 자동 실행
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

# ─── 백그라운드 프로세스 감시 ───

# 수학자가 board에 우선순위 명령을 남기면 실행
PRIORITY_FILE="$BOARD_DIR/priority.md"

check_and_enforce_priority() {
    # 1) 수학자 우선순위 파일에 kill 지시가 있으면 실행
    if [ -f "$PRIORITY_FILE" ]; then
        local kill_pids
        kill_pids=$(grep -oP '(?<=KILL_PID:)\d+' "$PRIORITY_FILE" 2>/dev/null)
        if [ -n "$kill_pids" ]; then
            for pid in $kill_pids; do
                if ps -p "$pid" > /dev/null 2>&1; then
                    warn "수학자 우선순위 명령: PID $pid kill"
                    kill "$pid" 2>/dev/null && ok "PID $pid 종료" || warn "PID $pid kill 실패"
                fi
            done
        fi
    fi

    # 2) 장시간 brute-force 탐사 자동 감지
    #    ps -eo로 경과시간(etimes, 초 단위) 사용 — locale 무관, 정확
    local long_running=""
    while IFS= read -r line; do
        local pid etime cmd
        pid=$(echo "$line" | awk '{print $1}')
        etime=$(echo "$line" | awk '{print $2}')
        cmd=$(echo "$line" | awk '{$1=$2=""; print substr($0,3)}')
        # 경과 시간 > 6시간 (21600초)
        if [ "$etime" -gt 21600 ] 2>/dev/null; then
            long_running="$long_running $pid"
            warn "6시간+ 장기 탐사: PID $pid (${etime}초): $cmd"
        fi
    done < <(ps -eo pid,etimes,args 2>/dev/null | grep -E "python.*(overnight|exploration|continuous)" | grep -v grep)

    if [ -n "$long_running" ]; then
        # 보드에 1회만 기록 (중복 방지: 같은 PID가 이미 있으면 건너뜀)
        for pid in $long_running; do
            if ! grep -q "PID.*$pid" "$BOARD_DIR/mathematician.md" 2>/dev/null; then
                warn "수학자에게 보고: PID $pid"
            fi
        done
    fi
}

# 수학자가 우선순위 파일에 작성할 수 있는 명령 형식:
# KILL_PID:12345          — 특정 PID 즉시 종료
# PRIORITY:high           — 다음 실험 긴급 (daemon 간격 단축)
# PAUSE_EXPLORATION:true  — 모든 탐사 프로세스 일시 중지

process_priority_commands() {
    [ ! -f "$PRIORITY_FILE" ] && return 0

    # 우선순위 레벨 확인
    local priority_level
    priority_level=$(grep -oP '(?<=PRIORITY:)\w+' "$PRIORITY_FILE" 2>/dev/null | tail -1)

    # 탐사 일시 중지 명령
    if grep -q 'PAUSE_EXPLORATION:true' "$PRIORITY_FILE" 2>/dev/null; then
        log "수학자 명령: 모든 탐사 프로세스 일시 중지"
        local exploration_pids
        exploration_pids=$(ps aux | grep -E "python.*(overnight|exploration|continuous|sweep)" \
            | grep -v grep | awk '{print $2}')
        for pid in $exploration_pids; do
            kill -STOP "$pid" 2>/dev/null && ok "PID $pid 일시 중지 (SIGSTOP)" || true
        done
    fi

    # 탐사 재개 명령
    if grep -q 'PAUSE_EXPLORATION:false' "$PRIORITY_FILE" 2>/dev/null; then
        log "수학자 명령: 탐사 프로세스 재개"
        local stopped_pids
        stopped_pids=$(ps aux | grep -E "python.*(overnight|exploration|continuous|sweep)" \
            | grep -v grep | awk '$8 ~ /T/ {print $2}')
        for pid in $stopped_pids; do
            kill -CONT "$pid" 2>/dev/null && ok "PID $pid 재개 (SIGCONT)" || true
        done
    fi

    # 처리 완료 표시
    if [ -f "$PRIORITY_FILE" ]; then
        local processed_time
        processed_time=$(date '+%Y-%m-%d %H:%M')
        sed -i "1i# 처리됨: $processed_time" "$PRIORITY_FILE"
    fi
}

# 수렴 판단: 탐사 프로세스가 실행 중일 때만 판단
check_exploration_convergence() {
    # 탐사 프로세스가 없으면 판단할 필요 없음
    local exploration_running
    exploration_running=$(ps aux | grep -E "python.*(overnight|exploration|continuous)" \
        | grep -v grep | wc -l)
    if [ "$exploration_running" -eq 0 ]; then
        return 0
    fi

    local summary_file="$PROJECT_DIR/outputs/overnight/exploration_summary.json"
    [ ! -f "$summary_file" ] && return 0

    # JSON 파일이 1시간 이상 오래되면 stale — 건너뜀
    local file_age
    file_age=$(( $(date +%s) - $(stat -c %Y "$summary_file" 2>/dev/null || echo 0) ))
    if [ "$file_age" -gt 3600 ]; then
        return 0
    fi

    local recent_f2
    recent_f2=$(python3 -c "
import json, sys, statistics
try:
    with open('$summary_file') as f:
        data = json.load(f)
    if isinstance(data, list) and len(data) >= 3:
        last3 = [d.get('f2_mean', 999) for d in data[-3:]]
        if any(v >= 999 for v in last3):
            print('SKIP:invalid_data')
            sys.exit(0)
        mean_v = statistics.mean(last3)
        std_v = statistics.stdev(last3) if len(last3) > 1 else 0
        cv = std_v / mean_v if mean_v > 0 else 999
        if cv < 0.15:
            print(f'CONVERGED:cv={cv:.3f},f2={mean_v:.4f}')
        else:
            print(f'VARYING:cv={cv:.3f}')
    else:
        print('SKIP:insufficient_data')
except Exception as e:
    print(f'ERROR:{e}')
" 2>/dev/null)

    if echo "$recent_f2" | grep -q "CONVERGED"; then
        # 이미 보고했으면 중복 방지
        if ! grep -q "탐사 수렴.*$(date '+%Y-%m-%d')" "$BOARD_DIR/mathematician.md" 2>/dev/null; then
            warn "탐사 수렴 감지: $recent_f2"
            warn "→ 수학적 과제에 CPU 재할당 권고"
        fi
    fi
}

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
    /home/k0who029/.npm-global/bin/claude -p "$user_prompt" \
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

# ─── 변경 감지 (IDLE 판단) ───

has_board_changed() {
    # 보드 파일이 마지막 사이클 이후 변경되었는가?
    local marker="$LOG_DIR/.last_cycle_time"
    if [ ! -f "$marker" ]; then
        return 0  # 마커 없으면 실행
    fi
    # 보드, 결과, 스크립트 중 하나라도 변경되면 true
    local changed
    changed=$(find "$BOARD_DIR" "$RESULTS_DIR" -newer "$marker" -name "*.md" -o -name "*.txt" 2>/dev/null | head -1)
    if [ -n "$changed" ]; then
        return 0  # 변경 있음 → 실행
    fi
    # 실행 중 실험이 방금 끝났으면 true
    local finished_experiment
    finished_experiment=$(find "$RESULTS_DIR" -name "*.txt" -newer "$marker" 2>/dev/null | head -1)
    if [ -n "$finished_experiment" ]; then
        return 0
    fi
    return 1  # 변경 없음 → 건너뜀
}

has_new_results() {
    # 마지막 사이클 이후 새 결과 파일이 생겼는가?
    local marker="$LOG_DIR/.last_cycle_time"
    if [ ! -f "$marker" ]; then
        return 0
    fi
    local new
    new=$(find "$RESULTS_DIR" "$ANALYSIS_DIR" -name "*.txt" -newer "$marker" 2>/dev/null | wc -l)
    [ "$new" -gt 0 ]
}

is_experiment_running() {
    # 실험 프로세스가 실행 중인가?
    local running
    running=$(ps aux | grep -E "qrop_env|python.*gl[0-9]|python.*blind|python.*maass" \
        | grep -v grep | wc -l)
    [ "$running" -gt 0 ]
}

# ─── git 자동 동기화 (연구 사이클용) ───

git_sync_research() {
    local msg_prefix="${1:-auto-sync}"

    cd "$PROJECT_DIR"

    # 변경 사항 있는지 확인
    if git diff --quiet HEAD 2>/dev/null && [ -z "$(git ls-files --others --exclude-standard 2>/dev/null)" ]; then
        return 0
    fi

    log "git 동기화..."

    # 결과 파일 커밋
    local new_results
    new_results=$(git ls-files --others --exclude-standard -- results/ 2>/dev/null | wc -l)
    if [ "$new_results" -gt 0 ]; then
        git add results/*.txt results/*.py 2>/dev/null || true
        git commit -m "$(cat <<EOF
results: ${msg_prefix} — new experiment data

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
EOF
)" 2>/dev/null && ok "git: 결과 커밋" || true
    fi

    # 스크립트/보드 커밋
    local script_changes
    script_changes=$(git diff --name-only HEAD -- scripts/ 2>/dev/null | wc -l)
    script_changes=$((script_changes + $(git ls-files --others --exclude-standard -- scripts/ 2>/dev/null | wc -l)))
    if [ "$script_changes" -gt 0 ]; then
        git add scripts/ 2>/dev/null || true
        git commit -m "$(cat <<EOF
scripts: ${msg_prefix} — board/code updates

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
EOF
)" 2>/dev/null && ok "git: 스크립트 커밋" || true
    fi

    # push (unpushed가 있으면)
    if git log --oneline origin/main..HEAD 2>/dev/null | head -1 | grep -q .; then
        git push origin main 2>&1 && ok "git push 완료" || warn "git push 실패"
    fi
}

# ─── 논문 루프 인터페이스 (1h→3d 일방향 소통) ───

PAPER_PENDING="$BOARD_DIR/paper_pending.md"

update_paper_pending() {
    # 수학자 보드에서 새 양성 판정을 paper_pending.md에 축적
    # 3일 논문 루프(run_paper_cycle.sh)가 이 파일을 읽고 논문에 반영

    # pending 파일 초기화 (없으면)
    if [ ! -f "$PAPER_PENDING" ]; then
        cat > "$PAPER_PENDING" << 'PEOF'
# 논문 반영 대기 결과 (Paper Pending Queue)

> 1시간 연구 루프(run_cycle.sh)가 자동으로 결과를 추가합니다.
> 3일 논문 루프(run_paper_cycle.sh)가 반영 후 제거합니다.

| # | 결과명 | 판정 | 추가일 | 카테고리 |
|---|--------|------|--------|---------|
PEOF
    fi

    # 수학자 보드에서 양성 판정 추출 (★ 포함 행)
    local today
    today=$(date '+%Y-%m-%d')

    # 로드맵 테이블에서 완료+양성인 항목 찾기
    while IFS='|' read -r _ num name status verdict _; do
        num=$(echo "$num" | tr -d '[:space:]')
        verdict=$(echo "$verdict" | tr -d '[:space:]')

        # 숫자가 아니면 건너뜀
        [[ "$num" =~ ^[0-9]+$ ]] || continue

        # 양성 판정만 (★ 포함)
        echo "$verdict" | grep -q '★' || continue

        # 이미 pending에 있으면 건너뜀
        grep -q "^| $num " "$PAPER_PENDING" 2>/dev/null && continue

        # 카테고리 판정 (키워드 기반)
        local category="A"
        if echo "$name" | grep -qi "GL(2)\|elliptic\|11a1\|37a1\|Ramanujan"; then
            category="C (GL(2))"
        elif echo "$name" | grep -qi "GL(3)\|sym²\|conductor"; then
            category="A (GL(3))"
        elif echo "$name" | grep -qi "quantum\|VQE\|DQPT"; then
            category="D (양자)"
        elif echo "$name" | grep -qi "GUE\|spacing\|Weyl"; then
            category="B (스펙트럼)"
        fi

        # pending에 추가
        echo "| $num | $(echo "$name" | xargs) | $verdict | $today | $category |" >> "$PAPER_PENDING"
        log "논문 대기큐 추가: #$num ($verdict)"

    done < <(grep '|' "$BOARD_DIR/mathematician.md" 2>/dev/null | grep -i '완료\|✅')
}

# ─── 전체 사이클 ───

run_full_cycle() {
    local cycle=$(increment_cycle)

    log "=========================================="
    log "RDL 연구 사이클 #$cycle 시작 ($TIMESTAMP)"
    log "=========================================="

    # 메모리 정리
    cleanup_working_memory

    # ── 수학적 우선순위 시스템 ──
    check_and_enforce_priority
    process_priority_commands
    check_exploration_convergence

    # ── IDLE 감지: 변경 없으면 경량 사이클 ──
    if ! has_board_changed; then
        if is_experiment_running; then
            log "실험 실행 중 + 보드 변경 없음 → 대기 (API 절약)"
            touch "$LOG_DIR/.last_cycle_time"
            return 0
        fi
        log "보드/결과 변경 없음 → IDLE 사이클 건너뜀 (API 절약)"
        touch "$LOG_DIR/.last_cycle_time"
        return 0
    fi

    # 실행 중인 실험 확인
    local experiment_running=false
    if is_experiment_running; then
        experiment_running=true
        warn "실행 중인 실험 감지"
        ps aux | grep -E "qrop_env|python.*gl[0-9]|python.*blind|python.*maass" \
            | grep -v grep | awk '{print "  PID "$2": "$11" "$12}' || true
    fi

    # ── Stage 1: 수학자 (Opus) ──
    local s1_ok=0
    run_stage 1 mathematician opus && s1_ok=1

    if [ "$s1_ok" -eq 0 ]; then
        err "수학자 단계 실패"
        write_failure_report 1 "$LOG_DIR/${TIMESTAMP}_stage1_mathematician.log"
        touch "$LOG_DIR/.last_cycle_time"
        return 1
    fi

    # 게이트 1: 실험 실행 중이면 Stage 2 건너뜀
    if [ "$experiment_running" = true ]; then
        log "실험 실행 중 → Stage 2,3 건너뜀 (수학자 판정만 기록)"
        touch "$LOG_DIR/.last_cycle_time"
        return 0
    fi

    gate_after_stage1 || {
        log "Stage 2,3 건너뜀 (수학자 대기 지시)"
        touch "$LOG_DIR/.last_cycle_time"
        return 0
    }

    # ── Stage 2: 설계자 (모델은 수학자 지시 or 기본 opus) ──
    local s2_model
    s2_model=$(grep -i '\*\*모델\*\*' "$BOARD_DIR/mathematician.md" 2>/dev/null \
        | tail -1 | sed 's/.*[：:] *//' | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')
    if [[ "$s2_model" != "opus" && "$s2_model" != "sonnet" && "$s2_model" != "haiku" ]]; then
        s2_model="opus"
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

    # ── Stage 3: 검토자 — 새 결과 있을 때만 실행 (API 절약) ──
    if has_new_results || [ "$s2_ok" -eq 1 ]; then
        local tex_file="$PROJECT_DIR/paper/source/unified_master_en.tex"
        local tex_mtime_before
        tex_mtime_before=$(stat -c %Y "$tex_file" 2>/dev/null || echo 0)

        run_stage 3 reviewer opus || {
            warn "검토자 단계 실패"
            write_failure_report 3 "$LOG_DIR/${TIMESTAMP}_stage3_reviewer.log"
        }

        # 논문 반영 검증 게이트
        local tex_mtime_after
        tex_mtime_after=$(stat -c %Y "$tex_file" 2>/dev/null || echo 0)

        if [ "$tex_mtime_after" -gt "$tex_mtime_before" ]; then
            ok "논문 TeX 업데이트 감지 — PDF 재컴파일 확인"
        else
            local unreflected=0
            if [ -f "$RESULTS_DIR/.reflected" ]; then
                unreflected=$(diff <(ls "$RESULTS_DIR"/*.txt 2>/dev/null | sort) \
                    <(sort -u "$RESULTS_DIR/.reflected") 2>/dev/null | grep -c "^<" || echo 0)
            else
                unreflected=$(ls "$RESULTS_DIR"/*.txt 2>/dev/null | wc -l)
            fi
            if [ "$unreflected" -gt 0 ]; then
                warn "논문 미반영 결과 ${unreflected}개 감지"
                if ! grep -q "논문 미반영.*$(date '+%Y-%m-%d')" "$BOARD_DIR/reviewer.md" 2>/dev/null; then
                    echo "⚠️ 논문 미반영 결과 ${unreflected}개 ($(date '+%Y-%m-%d %H:%M'))" >> "$BOARD_DIR/reviewer.md"
                fi
            fi
        fi
    else
        log "새 결과 없음 → Stage 3 (검토자) 건너뜀 (API 절약)"
    fi

    # ── 논문 루프 인터페이스: paper_pending.md에 결과 축적 ──
    update_paper_pending

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

    # git 자동 동기화 (코드+결과+보드 분리 커밋 + push)
    git_sync_research "research-cycle #$cycle"

    # 사이클 타임스탬프 마커 갱신
    touch "$LOG_DIR/.last_cycle_time"

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
        s2m=$(grep -i '\*\*모델\*\*' "$BOARD_DIR/mathematician.md" 2>/dev/null \
            | tail -1 | sed 's/.*[：:] *//' | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')
        [[ "$s2m" == "opus" || "$s2m" == "sonnet" || "$s2m" == "haiku" ]] || s2m="opus"
        log "Stage 2 모델: $s2m (수학자 지정)"
        run_stage 2 executor "$s2m"
        ;;
    3)      run_stage 3 reviewer opus ;;
    all)    run_full_cycle ;;
    daemon) run_daemon "${2:-60}" ;;
    status) show_status ;;
    lock)   lock ;;
    unlock) unlock ;;
    *)
        echo "사용법: $0 [1|2|3|all|daemon|status|lock|unlock]"
        echo "  1       — 수학자만 (Opus)"
        echo "  2       — 설계자 (수학자 지정 모델)"
        echo "  3       — 검토자만 (Opus)"
        echo "  all     — 전체 사이클 (기본값)"
        echo "  daemon  — 자동 반복 (기본 60분)"
        echo "  status  — 현재 상태 확인"
        echo "  lock    — 사용자 잠금 (데몬 일시정지)"
        echo "  unlock  — 사용자 잠금 해제"
        exit 1
        ;;
esac
