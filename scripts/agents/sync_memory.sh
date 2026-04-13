#!/bin/bash
# =============================================================================
# 에이전트 메모리 → Claude Code 메모리 동기화
# =============================================================================
# 에이전트 시스템(scripts/memory/)의 확립된 발견을
# Claude Code 세션 메모리(~/.claude/projects/.../memory/)에 요약 반영
#
# 사용법: ./sync_memory.sh [방향]
#   agent2claude  — 에이전트 → Claude Code (기본)
#   claude2agent  — Claude Code → 에이전트
#   status        — 양쪽 상태 확인
# =============================================================================

set -euo pipefail

PROJECT_DIR="$HOME/Desktop/gdl_unified"
AGENT_MEM="$PROJECT_DIR/scripts/memory"
CLAUDE_MEM="$HOME/.claude/projects/-home-k0who029/memory"

BLUE='\033[0;34m'
GREEN='\033[0;32m'
NC='\033[0m'

log() { echo -e "${BLUE}[sync]${NC} $1"; }
ok()  { echo -e "${GREEN}[sync] ✓${NC} $1"; }

case "${1:-status}" in
    status)
        echo "=== 에이전트 메모리 ==="
        echo "  semantic:"
        ls -1 "$AGENT_MEM/semantic/" 2>/dev/null | sed 's/^/    /'
        echo "  episodic: $(ls "$AGENT_MEM/episodic/" 2>/dev/null | wc -l)개"
        echo "  working:"
        ls -1 "$AGENT_MEM/working/" 2>/dev/null | sed 's/^/    /'
        echo ""
        echo "=== Claude Code 메모리 ==="
        echo "  총 파일: $(ls "$CLAUDE_MEM/"*.md 2>/dev/null | wc -l)개"
        echo "  최근 수정:"
        ls -lt "$CLAUDE_MEM/"*.md 2>/dev/null | head -5 | awk '{print "    "$NF" ("$6" "$7" "$8")"}'
        ;;

    agent2claude)
        log "에이전트 semantic → Claude Code 동기화..."

        # formal_propositions.md 요약을 Claude 메모리에 반영
        if [ -f "$AGENT_MEM/semantic/formal_propositions.md" ]; then
            local_ts=$(stat -c %Y "$AGENT_MEM/semantic/formal_propositions.md" 2>/dev/null || echo 0)
            claude_ts=0
            if [ -f "$CLAUDE_MEM/project_bundle_geometry.md" ]; then
                claude_ts=$(stat -c %Y "$CLAUDE_MEM/project_bundle_geometry.md" 2>/dev/null || echo 0)
            fi
            if [ "$local_ts" -gt "$claude_ts" ]; then
                log "formal_propositions가 더 최신 — Claude 메모리 갱신 필요"
                echo "  → 다음 Claude Code 세션에서 자동 반영됨"
            else
                ok "formal_propositions 동기화 상태 OK"
            fi
        fi

        # anti_patterns.md 확인
        if [ -f "$AGENT_MEM/semantic/anti_patterns.md" ]; then
            ok "anti_patterns.md 존재 ($(wc -l < "$AGENT_MEM/semantic/anti_patterns.md")줄)"
        fi

        # 최근 실패 반성문 요약
        recent_failures=$(find "$AGENT_MEM/episodic/" -name "failure_*.md" -newer "$AGENT_MEM/episodic/" -mmin -1440 2>/dev/null | wc -l)
        if [ "$recent_failures" -gt 0 ]; then
            log "최근 24시간 실패 반성문: ${recent_failures}개"
            find "$AGENT_MEM/episodic/" -name "failure_*.md" -mmin -1440 -exec basename {} \; 2>/dev/null | sed 's/^/  /'
        fi
        ;;

    claude2agent)
        log "Claude Code → 에이전트 동기화..."
        log "(현재 수동 — Claude Code 세션에서 에이전트 보드를 읽어 자동 반영)"
        ;;

    *)
        echo "사용법: $0 [agent2claude|claude2agent|status]"
        exit 1
        ;;
esac
