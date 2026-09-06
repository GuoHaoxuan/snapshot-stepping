#!/bin/bash
# 列出可搜索的天：BGO 逐小时 TTE 齐全（48 个文件 = 2 探头 x 24 小时）即可，
# NaI 有没有由 Chunk 自己决定（2017-10 之前没有，那时按单组搜）。
D=/hxmtfs/data/Fermi_GBM
out=/scratchfs2/gecam/guohx/gbmrun/days.txt
: > "$out"
for y in $(seq 2012 2020); do
    for d in $(ls $D/BGO/$y 2>/dev/null); do
        n=$(ls $D/BGO/$y/$d 2>/dev/null | grep -cE 'glg_tte_b[01]_[0-9]{6}_[0-9]{2}z')
        [ "$n" -eq 48 ] || continue
        # YYMMDD -> YYYY-MM-DD
        echo "20${d:0:2}-${d:2:2}-${d:4:2}" >> "$out"
    done
done
sort -u "$out" -o "$out"
echo "可搜天数: $(wc -l < "$out")"
echo "首: $(head -1 "$out")   末: $(tail -1 "$out")"
