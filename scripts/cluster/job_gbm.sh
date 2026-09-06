#!/bin/bash
# 第 $1 个 worker：处理 days.txt 中 (行号-1) % WORKERS == $1 的那些天。
cd /scratchfs2/gecam/guohx/gbmrun
WORKERS=100
idx=$1
n=0
while read -r day; do
    if [ $((n % WORKERS)) -eq "$idx" ]; then
        y=${day:0:4}; m=${day:5:2}; dd=${day:8:2}
        sig="data/Fermi_GBM/$y/$m/$y$m${dd}_signals.json"
        hrs="data/Fermi_GBM/$y/$m/$y$m${dd}_hours.json"
        if [ ! -f "$sig" ] || [ ! -f "$hrs" ]; then
            ./blink search "$day" "$day" --instrument fermi-gbm
        fi
    fi
    n=$((n + 1))
done < days.txt
