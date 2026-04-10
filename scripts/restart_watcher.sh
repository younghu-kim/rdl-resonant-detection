#!/bin/bash
# [5000,10000] 결과 파일 나타나면 프로세스 죽이고 재시작
RESULT_FILE=~/Desktop/gdl_unified/outputs/overnight/result_t5000-10000.json
LOG_FILE=~/Desktop/gdl_unified/outputs/overnight/run.log

echo "[watcher] result_t5000-10000.json 대기 중..."

while [ ! -f "$RESULT_FILE" ]; do
    sleep 30
done

echo "[watcher] [5000,10000] 완료 감지! 재시작합니다."

# 기존 프로세스 종료
pkill -f "overnight_exploration.py" 2>/dev/null
sleep 3

# 재시작 (수정된 코드 = 100 단위)
cd ~/Desktop/gdl_unified
nohup ~/qrop_env/bin/python3 -u scripts/overnight_exploration.py >> outputs/overnight/run.log 2>&1 &
echo "[watcher] 새 PID: $! — 100 단위 탐사 시작"
