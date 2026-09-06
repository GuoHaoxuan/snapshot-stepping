#!/bin/bash
cd /scratchfs2/gecam/guohx/v6run
n=$(($1+1))
day=$(sed -n "${n}p" days.txt)
[ -z "$day" ] && exit 0
y=${day:0:4}; m=${day:5:2}; dd=${day:8:2}
sig="data/Insight-HXMT_HE/$y/$m/$y$m${dd}_signals.json"
hrs="data/Insight-HXMT_HE/$y/$m/$y$m${dd}_hours.json"
[ -f "$sig" ] && [ -f "$hrs" ] && exit 0
# 计算节点上 /workfs2 只读、scratchfs2 inode 紧:错误日志成功即删,只有失败的天留档
log="farm_logs/e_${day}.log"
./blink search "$day" "$day" --workers 1 --worker 0 2> "$log"
code=$?
if [ "$code" -eq 0 ] && [ -f "$sig" ] && [ -f "$hrs" ]; then rm -f "$log"; fi
exit $code
