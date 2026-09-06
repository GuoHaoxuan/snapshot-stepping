#!/bin/bash
# 第 $1 个 worker：tasks.txt 里 (行号-1) % 100 == $1 的那些 (星, 天)。已产出的天跳过，重跑安全。
cd /scratchfs2/gecam/guohx/gridrun
WORKERS=100; idx=$1; n=0
while read -r sat day; do
  if [ $((n % WORKERS)) -eq "$idx" ]; then
    y=${day:0:4}; m=${day:5:2}; dd=${day:8:2}
    case $sat in grid02) D=GRID-02;; grid03b) D=GRID-03B;; grid04) D=GRID-04;; grid07) D=GRID-07;; esac
    sig="data/$D/$y/$m/$y$m${dd}_signals.json"; hrs="data/$D/$y/$m/$y$m${dd}_hours.json"
    if [ ! -f "$sig" ] || [ ! -f "$hrs" ]; then ./blink search "$day" "$day" --instrument $sat 2>> farm_logs/worker_$idx.err; fi
  fi
  n=$((n + 1))
done < tasks.txt
