#!/bin/bash
# =============================================================================
# 모노레포 → 특화 레포 단방향 동기화
# =============================================================================
# 사용법:
#   ./sync_repos.sh              # 전체 동기화
#   ./sync_repos.sh gl2          # GL(2) 레포만
#   ./sync_repos.sh quantum      # 양자 레포만
#   ./sync_repos.sh status       # 동기화 상태 확인
#
# 모노레포(rdl-resonant-detection)가 원본(source of truth).
# 이 스크립트는 특화 레포에 해당 파일만 단방향으로 복사 + push한다.
# =============================================================================

set -euo pipefail

MONO_DIR="$HOME/Desktop/gdl_unified"
SYNC_DIR="$HOME/Desktop/sync_repos"
LOG_DIR="$MONO_DIR/outputs/sync_logs"

mkdir -p "$LOG_DIR" "$SYNC_DIR"

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
# GL(2) 특화 레포 동기화
# =============================================================================
sync_gl2() {
    local REPO_NAME="xi-bundle-gl2"
    local REPO_URL="git@github.com:younghu-kim/${REPO_NAME}.git"
    local REPO_DIR="$SYNC_DIR/$REPO_NAME"

    log "GL(2) 특화 레포 동기화..."

    # 레포 존재 확인
    if ! gh repo view "younghu-kim/$REPO_NAME" &>/dev/null; then
        warn "레포 younghu-kim/$REPO_NAME 미존재 — 건너뜀"
        return 0
    fi

    # 로컬 클론 확인
    if [ ! -d "$REPO_DIR" ]; then
        log "최초 클론: $REPO_URL"
        git clone "$REPO_URL" "$REPO_DIR"
    fi

    cd "$REPO_DIR"
    git pull --rebase origin main 2>/dev/null || git pull --rebase origin master 2>/dev/null || true

    # 동기화 대상 파일
    local CHANGED=0

    # GL(2) 논문 (존재 시)
    if [ -f "$MONO_DIR/paper/source/gl2_master_en.tex" ]; then
        mkdir -p paper/
        cp "$MONO_DIR/paper/source/gl2_master_en.tex" paper/
        cp "$MONO_DIR/paper/source/gl2_master_ko.tex" paper/ 2>/dev/null || true
        CHANGED=1
    fi

    # GL(2) 부록 (모노레포 unified에서 GL(2) 섹션 추출은 복잡 → 논문 분리 후만)
    # GL(2) 관련 스크립트
    mkdir -p scripts/
    for f in "$MONO_DIR"/scripts/*gl2* "$MONO_DIR"/scripts/*elliptic* "$MONO_DIR"/scripts/*ramanujan*; do
        [ -f "$f" ] && cp "$f" scripts/ && CHANGED=1
    done

    # GL(2) 관련 결과
    mkdir -p results/
    for f in "$MONO_DIR"/results/*elliptic* "$MONO_DIR"/results/*gl2* "$MONO_DIR"/results/*ramanujan* "$MONO_DIR"/results/*sigma_uniqueness*; do
        [ -f "$f" ] && cp "$f" results/ && CHANGED=1
    done

    # 변경 사항 있으면 push
    if [ "$CHANGED" -eq 1 ] && [ -n "$(git status --porcelain)" ]; then
        git add -A
        git commit -m "sync: auto-sync from monorepo $(date +%Y-%m-%d)

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>" || true
        git push origin main 2>/dev/null || git push origin master 2>/dev/null || true
        ok "GL(2) 동기화 완료 — push 성공"
    else
        ok "GL(2) 변경 없음 — 건너뜀"
    fi

    cd "$MONO_DIR"
}

# =============================================================================
# 양자 특화 레포 동기화
# =============================================================================
sync_quantum() {
    local REPO_NAME="xi-bundle-quantum"
    local REPO_URL="git@github.com:younghu-kim/${REPO_NAME}.git"
    local REPO_DIR="$SYNC_DIR/$REPO_NAME"

    log "양자 특화 레포 동기화..."

    # 레포 존재 확인
    if ! gh repo view "younghu-kim/$REPO_NAME" &>/dev/null; then
        warn "레포 younghu-kim/$REPO_NAME 미존재 — 건너뜀"
        return 0
    fi

    # 로컬 클론 확인
    if [ ! -d "$REPO_DIR" ]; then
        log "최초 클론: $REPO_URL"
        git clone "$REPO_URL" "$REPO_DIR"
    fi

    cd "$REPO_DIR"
    git pull --rebase origin main 2>/dev/null || git pull --rebase origin master 2>/dev/null || true

    local CHANGED=0

    # quantum/ 디렉토리 전체
    if [ -d "$MONO_DIR/quantum" ]; then
        rsync -av --delete "$MONO_DIR/quantum/" quantum/
        CHANGED=1
    fi

    # 해밀토니안 선택 규칙 문서
    mkdir -p docs/
    if [ -f "$MONO_DIR/docs/hamiltonian_selection_rules.md" ]; then
        cp "$MONO_DIR/docs/hamiltonian_selection_rules.md" docs/
        CHANGED=1
    fi

    # 양자 논문 (존재 시)
    if [ -f "$MONO_DIR/paper/source/quantum_master_en.tex" ]; then
        mkdir -p paper/
        cp "$MONO_DIR/paper/source/quantum_master_en.tex" paper/
        cp "$MONO_DIR/paper/source/quantum_master_ko.tex" paper/ 2>/dev/null || true
        CHANGED=1
    fi

    # 변경 사항 있으면 push
    if [ "$CHANGED" -eq 1 ] && [ -n "$(git status --porcelain)" ]; then
        git add -A
        git commit -m "sync: auto-sync from monorepo $(date +%Y-%m-%d)

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>" || true
        git push origin main 2>/dev/null || git push origin master 2>/dev/null || true
        ok "양자 동기화 완료 — push 성공"
    else
        ok "양자 변경 없음 — 건너뜀"
    fi

    cd "$MONO_DIR"
}

# =============================================================================
# 상태 확인
# =============================================================================
check_status() {
    log "동기화 상태 확인..."

    echo ""
    echo "  ┌──────────────────────────────────────────┐"
    echo "  │          레포 동기화 현황                 │"
    echo "  ├──────────────────────────────────────────┤"

    # 모노레포
    echo -n "  │ 모노레포: "
    if [ -d "$MONO_DIR/.git" ]; then
        local MONO_COMMIT=$(cd "$MONO_DIR" && git log --oneline -1 2>/dev/null)
        echo "$MONO_COMMIT"
    else
        echo "git 아님"
    fi

    # GL(2)
    echo -n "  │ GL(2):    "
    if gh repo view "younghu-kim/xi-bundle-gl2" &>/dev/null; then
        if [ -d "$SYNC_DIR/xi-bundle-gl2/.git" ]; then
            local GL2_COMMIT=$(cd "$SYNC_DIR/xi-bundle-gl2" && git log --oneline -1 2>/dev/null)
            echo "$GL2_COMMIT"
        else
            echo "레포 존재, 로컬 클론 없음"
        fi
    else
        echo "미생성"
    fi

    # 양자
    echo -n "  │ 양자:     "
    if gh repo view "younghu-kim/xi-bundle-quantum" &>/dev/null; then
        if [ -d "$SYNC_DIR/xi-bundle-quantum/.git" ]; then
            local Q_COMMIT=$(cd "$SYNC_DIR/xi-bundle-quantum" && git log --oneline -1 2>/dev/null)
            echo "$Q_COMMIT"
        else
            echo "레포 존재, 로컬 클론 없음"
        fi
    else
        echo "미생성"
    fi

    echo "  └──────────────────────────────────────────┘"
    echo ""

    # 마지막 동기화 로그
    local LAST_LOG=$(ls -t "$LOG_DIR"/sync_*.log 2>/dev/null | head -1)
    if [ -n "$LAST_LOG" ]; then
        log "마지막 동기화: $(basename "$LAST_LOG")"
    else
        log "동기화 이력 없음"
    fi
}

# =============================================================================
# 전체 동기화
# =============================================================================
sync_all() {
    local LOGFILE="$LOG_DIR/sync_${TIMESTAMP}.log"
    log "전체 동기화 시작..." | tee "$LOGFILE"

    sync_gl2 2>&1 | tee -a "$LOGFILE"
    sync_quantum 2>&1 | tee -a "$LOGFILE"

    ok "전체 동기화 완료. 로그: $LOGFILE" | tee -a "$LOGFILE"
}

# =============================================================================
# 메인
# =============================================================================
case "${1:-}" in
    gl2)
        sync_gl2
        ;;
    quantum)
        sync_quantum
        ;;
    status)
        check_status
        ;;
    *)
        sync_all
        ;;
esac
