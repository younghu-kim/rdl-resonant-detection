#!/bin/bash
# =============================================================================
# RDL 3단 논문 관리 사이클 오케스트레이터 v1
# =============================================================================
# 사용법:
#   ./run_paper_cycle.sh          # 전체 3단계 1회 실행
#   ./run_paper_cycle.sh 1        # Stage 1 (편집장)만
#   ./run_paper_cycle.sh 2        # Stage 2 (집필자)만
#   ./run_paper_cycle.sh 3        # Stage 3 (교정자)만
#   ./run_paper_cycle.sh daemon   # 3시간마다 자동 실행 (IDLE 시 API 0)
#   ./run_paper_cycle.sh status   # 현재 상태 확인
#
# 설계:
#   - 1시간 연구 루프(run_cycle.sh)가 board/paper_pending.md에 결과를 축적
#   - 3일 논문 루프(이 스크립트)가 축적된 결과를 논문에 반영
#   - 편집장 → 집필자 → 교정자 3단 파이프라인
#   - 교정자가 NEEDS_FIX 판정 시 → 집필자 재실행 (문답 루프, 최대 3회)
#   - IDLE 감지: 새 결과 없으면 API 호출 0
#
# 소통 구조:
#   [1h 루프] ──일방향──→ paper_pending.md (결과 축적)
#                              ↓
#   [3d 루프] 편집장 ← 읽음 + mathematician.md + boundaries.md
#                 ↓ 지시 (paper_editor.md)
#             집필자 → TeX 수정 + 컴파일 + 배포
#                 ↓ 결과
#             교정자 → 검증 (paper_proofreader.md)
#                 ↓ NEEDS_FIX?
#             집필자 재실행 (문답 루프, 최대 3회)
# =============================================================================

set -euo pipefail

PROJECT_DIR="$HOME/Desktop/gdl_unified"
AGENTS_DIR="$PROJECT_DIR/scripts/agents"
LOG_DIR="$PROJECT_DIR/outputs/paper_cycle_logs"
BOARD_DIR="$PROJECT_DIR/scripts/board"
RESULTS_DIR="$PROJECT_DIR/results"
PAPER_DIR="$PROJECT_DIR/paper/source"
DEPLOY_DIR="$HOME/Desktop/수학최종논문"
JOURNAL="$PROJECT_DIR/scripts/research_journal.md"
CYCLE_COUNT_FILE="$LOG_DIR/.paper_cycle_count"
LOCK_FILE="/tmp/rdl_paper_cycle"
USER_LOCK="/tmp/rdl_user_active"
LAST_UPDATE_MARKER="$LOG_DIR/.last_paper_update_time"
PENDING_FILE="$BOARD_DIR/paper_pending.md"

mkdir -p "$LOG_DIR" "$BOARD_DIR" "$DEPLOY_DIR"

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

get_paper_cycle() {
    if [ -f "$CYCLE_COUNT_FILE" ]; then
        cat "$CYCLE_COUNT_FILE"
    else
        echo 0
    fi
}

increment_paper_cycle() {
    local count=$(get_paper_cycle)
    echo $((count + 1)) > "$CYCLE_COUNT_FILE"
    echo $((count + 1))
}

# ─── 사용자 잠금 ───

check_user_lock() {
    if [ -f "$USER_LOCK" ]; then
        local lock_age
        lock_age=$(( $(date +%s) - $(stat -c %Y "$USER_LOCK" 2>/dev/null || echo 0) ))
        if [ "$lock_age" -gt 7200 ]; then
            warn "사용자 잠금 stale (${lock_age}초) — 제거"
            rm -f "$USER_LOCK"
            return 0
        fi
        warn "사용자 활성 세션 — 논문 사이클 건너뜀"
        return 1
    fi
    return 0
}

# ─── 연구 루프 충돌 방지 ───

is_research_cycle_running() {
    # run_cycle.sh가 claude를 호출 중인지 확인
    if pgrep -f "run_cycle.sh" > /dev/null 2>&1; then
        # claude CLI가 실제로 실행 중인지도 확인
        if pgrep -f "claude.*stage[0-9]" > /dev/null 2>&1; then
            return 0  # 실행 중
        fi
    fi
    return 1  # 미실행
}

# ─── IDLE 감지 (핵심 최적화) ───

has_new_material() {
    # paper_pending.md에 반영 대기 결과가 있는가?
    if [ -f "$PENDING_FILE" ]; then
        local pending_count
        pending_count=$(grep -c '^| [0-9]' "$PENDING_FILE" 2>/dev/null || echo 0)
        if [ "$pending_count" -gt 0 ]; then
            log "반영 대기 결과 ${pending_count}건 감지"
            return 0
        fi
    fi

    # 교정자가 NEEDS_FIX를 남겼는가?
    if [ -f "$BOARD_DIR/paper_proofreader.md" ]; then
        if grep -q "NEEDS_FIX\|동기화 필요\|컴파일 실패" "$BOARD_DIR/paper_proofreader.md" 2>/dev/null; then
            log "교정자 NEEDS_FIX 감지"
            return 0
        fi
    fi

    # 마커 이후 수학자 보드에 새 양성 판정이 있는가?
    if [ -f "$LAST_UPDATE_MARKER" ]; then
        if [ "$BOARD_DIR/mathematician.md" -nt "$LAST_UPDATE_MARKER" ]; then
            local new_stars
            new_stars=$(grep -c '★' "$BOARD_DIR/mathematician.md" 2>/dev/null || echo 0)
            local prev_stars
            prev_stars=$(cat "$LOG_DIR/.prev_star_count" 2>/dev/null || echo 0)
            if [ "$new_stars" -gt "$prev_stars" ]; then
                log "새 양성 판정 감지 (★ ${prev_stars}→${new_stars})"
                return 0
            fi
        fi
    else
        return 0  # 마커 없음 = 첫 실행
    fi

    return 1  # 새 자료 없음
}

# ─── 단계 실행 함수 ───

run_paper_stage() {
    local stage=$1
    local name=$2
    local model="${3:-opus}"
    local extra_context="${4:-}"
    local prompt_file="$AGENTS_DIR/paper_stage${stage}_${name}.md"
    local log_file="$LOG_DIR/${TIMESTAMP}_paper_stage${stage}_${name}.log"

    if [ ! -f "$prompt_file" ]; then
        err "프롬프트 파일 없음: $prompt_file"
        return 1
    fi

    log "Paper Stage $stage: $name 시작 (model: $model)..."

    local user_prompt="프로젝트 디렉토리: $PROJECT_DIR
현재 시각: $(date '+%Y-%m-%d %H:%M')
논문 사이클: #$(get_paper_cycle)

위 지침에 따라 이번 사이클을 수행하세요.
보드 파일을 읽고, 상황을 파악하고, 당신의 역할에 맞는 작업을 수행하세요.
작업이 끝나면 보드 파일에 결과를 기록하세요."

    # 문답 루프에서 전달된 추가 컨텍스트
    if [ -n "$extra_context" ]; then
        user_prompt="${user_prompt}

## 교정자 피드백 (이전 라운드에서 발견된 문제)
${extra_context}"
    fi

    cd "$PROJECT_DIR"
    /home/k0who029/.npm-global/bin/claude -p "$user_prompt" \
        --system-prompt-file "$prompt_file" \
        --output-format text \
        --model "$model" \
        --allowedTools "Read,Write,Edit,Bash,Glob,Grep" \
        > "$log_file" 2>&1

    local exit_code=$?

    if [ $exit_code -eq 0 ]; then
        ok "Paper Stage $stage: $name 완료 (로그: $log_file)"
    else
        err "Paper Stage $stage: $name 실패 (exit=$exit_code)"
        tail -10 "$log_file" 2>/dev/null || true
    fi

    return $exit_code
}

# ─── 검증 게이트 ───

gate_after_editor() {
    log "게이트 1: 편집장 지시 확인..."

    if ! [ -f "$BOARD_DIR/paper_editor.md" ]; then
        warn "편집장 보드 파일 없음"
        return 1
    fi

    if grep -q "NO_CHANGES_NEEDED\|변경 불필요" "$BOARD_DIR/paper_editor.md" 2>/dev/null; then
        log "편집장: 변경 불필요 → Stage 2,3 건너뜀"
        return 1
    fi

    # 오늘 날짜의 지시인지 확인
    if ! grep -q "$(date '+%Y-%m-%d')" "$BOARD_DIR/paper_editor.md" 2>/dev/null; then
        warn "편집장 보드가 오늘 날짜가 아님 — 오래된 지시 무시"
        return 1
    fi

    ok "게이트 1 (편집장) 통과"
    return 0
}

gate_after_writer() {
    log "게이트 2: 컴파일 확인..."

    local en_pdf="$PAPER_DIR/unified_master_en.pdf"
    local ko_pdf="$PAPER_DIR/unified_master_ko.pdf"

    # PDF 존재 확인
    if [ ! -f "$en_pdf" ] || [ ! -f "$ko_pdf" ]; then
        err "PDF 파일 누락"
        return 1
    fi

    # PDF가 최근 30분 이내인지
    local now
    now=$(date +%s)
    local en_age=$(( now - $(stat -c %Y "$en_pdf" 2>/dev/null || echo 0) ))
    if [ "$en_age" -gt 1800 ]; then
        warn "PDF가 30분 이상 오래됨 — 컴파일 실패 가능"
        return 1
    fi

    # LaTeX 치명적 에러 확인
    if grep -qi "Fatal error\|Emergency stop" "$PAPER_DIR/unified_master_en.log" 2>/dev/null; then
        err "EN TeX 컴파일 치명적 에러"
        return 1
    fi

    ok "게이트 2 (집필자) 통과"
    return 0
}

# ─── 실패 시 반성문 ───

write_paper_failure() {
    local stage=$1
    local log_file=$2
    local report="$LOG_DIR/failure_$(date +%Y%m%d)_paper_stage${stage}.md"

    cat > "$report" << EOF
# 논문 사이클 실패: Paper Stage $stage ($TIMESTAMP)

## 에러 로그 (마지막 20줄)
\`\`\`
$(tail -20 "$log_file" 2>/dev/null)
\`\`\`

## 원인
(다음 사이클에서 편집장이 분석)
EOF

    warn "반성문 작성: $report"
}

# ─── PDF 배포 ───

deploy_pdfs() {
    log "PDF 배포..."
    mkdir -p "$DEPLOY_DIR"

    for lang in en ko; do
        local src="$PAPER_DIR/unified_master_${lang}.pdf"
        if [ -f "$src" ]; then
            cp "$src" "$DEPLOY_DIR/" 2>/dev/null && ok "${lang} PDF → 수학최종논문/" || true
            cp "$src" "$PROJECT_DIR/paper/" 2>/dev/null || true
        fi
    done
}

# ─── pending 파일 정리 ───

clear_pending() {
    # 반영 완료된 결과를 pending에서 제거
    if [ -f "$PENDING_FILE" ]; then
        local backup="$LOG_DIR/paper_pending_${TIMESTAMP}.bak"
        cp "$PENDING_FILE" "$backup"

        # pending 파일을 빈 헤더만 남김
        cat > "$PENDING_FILE" << 'EOF'
# 논문 반영 대기 결과 (Paper Pending Queue)

> 1시간 연구 루프(run_cycle.sh)가 자동으로 결과를 추가합니다.
> 3일 논문 루프(run_paper_cycle.sh)가 반영 후 제거합니다.

| # | 결과명 | 판정 | 추가일 | 카테고리 |
|---|--------|------|--------|---------|
EOF
        ok "pending 큐 정리 완료 (백업: $backup)"
    fi
}

# ─── git 동기화 (논문+코드+결과 분리 커밋 + push) ───

git_sync_all() {
    local msg_prefix="${1:-auto-sync}"

    cd "$PROJECT_DIR"

    # 변경 사항 있는지 확인
    if git diff --quiet HEAD && [ -z "$(git ls-files --others --exclude-standard)" ]; then
        log "git: 변경 사항 없음"
        return 0
    fi

    log "git 동기화 시작..."

    # 1. 논문 커밋 (paper/)
    local paper_changes
    paper_changes=$(git diff --name-only HEAD -- paper/ 2>/dev/null | wc -l)
    paper_changes=$((paper_changes + $(git ls-files --others --exclude-standard -- paper/ 2>/dev/null | wc -l)))
    if [ "$paper_changes" -gt 0 ]; then
        git add paper/source/unified_master_en.tex paper/source/unified_master_ko.tex 2>/dev/null || true
        git add paper/source/unified_master_en.pdf paper/source/unified_master_ko.pdf 2>/dev/null || true
        git add paper/unified_master_en.pdf paper/unified_master_ko.pdf 2>/dev/null || true
        git commit -m "$(cat <<EOF
paper: ${msg_prefix} — TeX + PDF update

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
EOF
)" 2>/dev/null && ok "git: 논문 커밋 완료" || true
    fi

    # 2. 결과 커밋 (results/)
    local result_changes
    result_changes=$(git ls-files --others --exclude-standard -- results/ 2>/dev/null | wc -l)
    result_changes=$((result_changes + $(git diff --name-only HEAD -- results/ 2>/dev/null | wc -l)))
    if [ "$result_changes" -gt 0 ]; then
        git add results/*.txt 2>/dev/null || true
        git add results/*.py 2>/dev/null || true
        git commit -m "$(cat <<EOF
results: ${msg_prefix} — new experiment data

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
EOF
)" 2>/dev/null && ok "git: 결과 커밋 완료" || true
    fi

    # 3. 코드/스크립트 커밋 (scripts/)
    local script_changes
    script_changes=$(git diff --name-only HEAD -- scripts/ 2>/dev/null | wc -l)
    script_changes=$((script_changes + $(git ls-files --others --exclude-standard -- scripts/ 2>/dev/null | wc -l)))
    if [ "$script_changes" -gt 0 ]; then
        git add scripts/ 2>/dev/null || true
        git commit -m "$(cat <<EOF
scripts: ${msg_prefix} — board/agent/code updates

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
EOF
)" 2>/dev/null && ok "git: 스크립트 커밋 완료" || true
    fi

    # 4. 나머지 (board/ 등 루트 레벨)
    local other_changes
    other_changes=$(git diff --name-only HEAD 2>/dev/null | wc -l)
    other_changes=$((other_changes + $(git ls-files --others --exclude-standard 2>/dev/null | grep -v -E '^(paper|results|scripts)/' | wc -l)))
    if [ "$other_changes" -gt 0 ]; then
        git add -A 2>/dev/null || true
        git commit -m "$(cat <<EOF
chore: ${msg_prefix} — misc updates

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
EOF
)" 2>/dev/null && ok "git: 기타 커밋 완료" || true
    fi

    # 5. push
    if git log --oneline origin/master..HEAD 2>/dev/null | head -1 | grep -q .; then
        git push origin master 2>&1 && ok "git push 완료" || warn "git push 실패"
    fi
}

# ─── 문답 루프 (집필자 ↔ 교정자) ───

run_dialogue_loop() {
    local max_rounds=3
    local round=0

    while [ "$round" -lt "$max_rounds" ]; do
        ((round++))
        log "══════ 문답 라운드 $round/$max_rounds ══════"

        # ── 집필자 실행 ──
        local writer_context=""
        if [ "$round" -gt 1 ]; then
            # 이전 교정자 피드백을 집필자에게 전달
            writer_context=$(grep -A 50 "수정 필요 사항" "$BOARD_DIR/paper_proofreader.md" 2>/dev/null \
                | head -30 || echo "")
            log "교정자 피드백을 집필자에게 전달 (라운드 $round)"
        fi

        local s2_ok=0
        run_paper_stage 2 writer sonnet "$writer_context" && s2_ok=1

        if [ "$s2_ok" -eq 0 ]; then
            err "집필자 실패 (라운드 $round)"
            write_paper_failure 2 "$LOG_DIR/${TIMESTAMP}_paper_stage2_writer.log"
            return 1
        fi

        # 게이트 2: 컴파일 확인
        if ! gate_after_writer; then
            warn "컴파일 실패 — 교정자 건너뛰고 집필자 재시도"
            continue
        fi

        # ── 교정자 실행 ──
        run_paper_stage 3 proofreader sonnet || {
            warn "교정자 실패 (라운드 $round)"
            write_paper_failure 3 "$LOG_DIR/${TIMESTAMP}_paper_stage3_proofreader.log"
            # 교정자 실패해도 집필자 결과는 유지 — 루프 종료
            break
        }

        # ── 교정자 판정 확인 ──
        if grep -q "종합 판정: PASS" "$BOARD_DIR/paper_proofreader.md" 2>/dev/null; then
            ok "교정자 PASS — 문답 루프 종료 (라운드 $round)"
            return 0
        fi

        if grep -q "NEEDS_FIX" "$BOARD_DIR/paper_proofreader.md" 2>/dev/null; then
            if [ "$round" -lt "$max_rounds" ]; then
                warn "교정자 NEEDS_FIX — 집필자 재실행 (라운드 $((round+1)))"
                # 타임스탬프 갱신 (로그 파일명 중복 방지)
                TIMESTAMP=$(date +%Y%m%d_%H%M%S)
            else
                warn "문답 라운드 $max_rounds 소진 — NEEDS_FIX 잔존, 다음 사이클로 이월"
            fi
        else
            # PASS도 NEEDS_FIX도 아니면 (교정자가 이상한 출력) 루프 종료
            warn "교정자 판정 불명확 — 루프 종료"
            break
        fi
    done

    return 0
}

# ─── 전체 사이클 ───

run_full_paper_cycle() {
    local cycle=$(increment_paper_cycle)

    log "══════════════════════════════════════════"
    log "RDL 논문 관리 사이클 #$cycle 시작 ($TIMESTAMP)"
    log "══════════════════════════════════════════"

    # ── IDLE 감지 ──
    if ! has_new_material; then
        log "새 자료 없음 → IDLE (API 호출 0)"
        touch "$LAST_UPDATE_MARKER"
        return 0
    fi

    # ── 사용자 잠금 ──
    if ! check_user_lock; then
        return 0
    fi

    # ── 연구 루프 충돌 방지 ──
    local wait_count=0
    while is_research_cycle_running; do
        if [ "$wait_count" -ge 6 ]; then
            warn "연구 루프 30분 이상 실행 중 — 논문 사이클 건너뜀"
            return 0
        fi
        log "연구 루프 실행 중 — 5분 대기 ($((wait_count+1))/6)"
        sleep 300
        ((wait_count++))
    done

    # ── Stage 1: 편집장 (Editor-in-Chief) ──
    log "── Stage 1: 편집장 ──"
    local s1_ok=0
    run_paper_stage 1 editor opus && s1_ok=1

    if [ "$s1_ok" -eq 0 ]; then
        err "편집장 실패"
        write_paper_failure 1 "$LOG_DIR/${TIMESTAMP}_paper_stage1_editor.log"
        touch "$LAST_UPDATE_MARKER"
        return 1
    fi

    # 게이트 1: 변경 필요 여부
    if ! gate_after_editor; then
        log "변경 불필요 — 사이클 종료"
        touch "$LAST_UPDATE_MARKER"
        # ★ 카운터 갱신 (다음 IDLE 판단용)
        grep -c '★' "$BOARD_DIR/mathematician.md" > "$LOG_DIR/.prev_star_count" 2>/dev/null || true
        return 0
    fi

    # ── Stage 2+3: 문답 루프 (집필자 ↔ 교정자, 최대 3회) ──
    log "── Stage 2+3: 문답 루프 ──"

    # 변경 전 줄 수 기록 (안전장치)
    local en_lines_before ko_lines_before
    en_lines_before=$(wc -l < "$PAPER_DIR/unified_master_en.tex" 2>/dev/null || echo 0)
    ko_lines_before=$(wc -l < "$PAPER_DIR/unified_master_ko.tex" 2>/dev/null || echo 0)

    run_dialogue_loop

    # 변경 후 줄 수 확인 (무손실 안전장치)
    local en_lines_after ko_lines_after
    en_lines_after=$(wc -l < "$PAPER_DIR/unified_master_en.tex" 2>/dev/null || echo 0)
    ko_lines_after=$(wc -l < "$PAPER_DIR/unified_master_ko.tex" 2>/dev/null || echo 0)

    if [ "$en_lines_after" -lt "$en_lines_before" ]; then
        err "⚠ EN 줄 수 감소! ${en_lines_before} → ${en_lines_after} — 내용 손실 경고"
    fi
    if [ "$ko_lines_after" -lt "$ko_lines_before" ]; then
        err "⚠ KO 줄 수 감소! ${ko_lines_before} → ${ko_lines_after} — 내용 손실 경고"
    fi

    # ── 후처리 ──

    # PDF 배포
    deploy_pdfs

    # pending 큐 정리
    clear_pending

    # git commit + push (논문 + 코드 + 결과 + 보드 전부)
    git_sync_all "paper-cycle #$cycle"

    # 특화 레포 동기화
    if [ -x "$AGENTS_DIR/sync_repos.sh" ]; then
        log "특화 레포 동기화..."
        "$AGENTS_DIR/sync_repos.sh" 2>&1 | tee -a "$LOG_DIR/${TIMESTAMP}_sync.log" || true
    fi

    # 마커 갱신
    touch "$LAST_UPDATE_MARKER"
    grep -c '★' "$BOARD_DIR/mathematician.md" > "$LOG_DIR/.prev_star_count" 2>/dev/null || true

    # 일지 기록
    echo "" >> "$JOURNAL"
    echo "## $(date '+%Y-%m-%d %H:%M') 논문 사이클 #$cycle" >> "$JOURNAL"
    echo "- 편집장: $([ $s1_ok -eq 1 ] && echo '완료' || echo '실패')" >> "$JOURNAL"
    echo "- EN: ${en_lines_before}→${en_lines_after} lines" >> "$JOURNAL"
    echo "- KO: ${ko_lines_before}→${ko_lines_after} lines" >> "$JOURNAL"
    echo "- 로그: $LOG_DIR/${TIMESTAMP}_*" >> "$JOURNAL"

    log "══════════════════════════════════════════"
    ok "논문 사이클 #$cycle 완료"
    log "══════════════════════════════════════════"
}

# ─── 데몬 모드 ───

run_daemon() {
    local interval_min=${1:-180}
    local interval_sec=$((interval_min * 60))

    log "══════════════════════════════════════════"
    log "RDL 3단 논문 관리 데몬 v1 시작"
    log "간격: ${interval_min}분, 프로젝트: $PROJECT_DIR"
    log "잠금: $USER_LOCK (사용자 우선)"
    log "논문 사이클: #$(get_paper_cycle)부터"
    log "══════════════════════════════════════════"

    while true; do
        if check_user_lock; then
            run_full_paper_cycle || true
        fi
        log "다음 논문 사이클: ${interval_min}분 후"
        sleep "$interval_sec"
    done
}

# ─── 상태 확인 ───

show_paper_status() {
    echo -e "${CYAN}=== RDL 3단 논문 관리 시스템 상태 ===${NC}"
    echo ""
    echo "논문 사이클: #$(get_paper_cycle)"
    echo "사용자 잠금: $([ -f "$USER_LOCK" ] && echo "활성" || echo "없음")"
    echo ""

    echo "--- 논문 현황 ---"
    if [ -f "$PAPER_DIR/unified_master_en.tex" ]; then
        echo "  EN: $(wc -l < "$PAPER_DIR/unified_master_en.tex") lines, $(stat -c '%y' "$PAPER_DIR/unified_master_en.tex" 2>/dev/null | cut -d. -f1)"
        echo "  KO: $(wc -l < "$PAPER_DIR/unified_master_ko.tex") lines, $(stat -c '%y' "$PAPER_DIR/unified_master_ko.tex" 2>/dev/null | cut -d. -f1)"
    fi
    echo ""

    echo "--- 반영 대기 큐 ---"
    if [ -f "$PENDING_FILE" ]; then
        local pending
        pending=$(grep -c '^| [0-9]' "$PENDING_FILE" 2>/dev/null || echo 0)
        echo "  대기 결과: ${pending}건"
        grep '^| [0-9]' "$PENDING_FILE" 2>/dev/null | while read -r line; do
            echo "  $line"
        done
    else
        echo "  대기 큐 없음"
    fi
    echo ""

    echo "--- 보드 ---"
    for f in paper_editor.md paper_proofreader.md paper_pending.md; do
        if [ -f "$BOARD_DIR/$f" ]; then
            echo "  $f: $(head -1 "$BOARD_DIR/$f" 2>/dev/null) ($(stat -c '%y' "$BOARD_DIR/$f" 2>/dev/null | cut -d. -f1))"
        else
            echo "  $f: 없음"
        fi
    done
    echo ""

    echo "--- IDLE 판단 ---"
    if has_new_material; then
        echo "  상태: 실행 필요 (새 자료 있음)"
    else
        echo "  상태: IDLE (새 자료 없음)"
    fi
    echo ""

    echo "--- 최근 로그 ---"
    ls -lt "$LOG_DIR"/*.log 2>/dev/null | head -5 | awk '{print "  "$NF" ("$6" "$7" "$8")"}'
    echo ""

    echo "--- 서비스 ---"
    systemctl --user is-active rdl-paper-cycle.service 2>/dev/null && \
        echo "  rdl-paper-cycle: 활성" || echo "  rdl-paper-cycle: 비활성"
}

# ─── 메인 ───

case "${1:-all}" in
    1)      run_paper_stage 1 editor opus ;;
    2)      run_paper_stage 2 writer sonnet ;;
    3)      run_paper_stage 3 proofreader sonnet ;;
    all)    run_full_paper_cycle ;;
    daemon) run_daemon "${2:-180}" ;;
    status) show_paper_status ;;
    *)
        echo "사용법: $0 [1|2|3|all|daemon|status]"
        echo "  1       — 편집장만 (Opus)"
        echo "  2       — 집필자만 (Opus)"
        echo "  3       — 교정자만 (Opus)"
        echo "  all     — 전체 사이클 (기본값)"
        echo "  daemon  — 자동 반복 (기본 180분)"
        echo "  status  — 현재 상태 확인"
        exit 1
        ;;
esac
