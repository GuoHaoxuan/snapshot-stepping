#!/bin/bash
# 多段活时间的逐字节回归：GRID 4 天 vs v8、SVOM 1 天 vs v5、HXMT 1 天 vs attcheck（v8 二进制产物）
export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH
R=/scratchfs2/gecam/guohx; log=$R/gridrun/farm_logs/chain_livetime.log
say() { echo "$(date +%H:%M) $*" >> $log; }
say "chain start"
cd $R/gridrun; rm -f farm_logs/build_v10.out
hep_sub -g hxmt -mem 12000 -cpu 8 job_build_v10.sh >> $log 2>&1
while ! grep -q BUILD_EXIT farm_logs/build_v10.out 2>/dev/null; do sleep 30; done
grep -q "BUILD_EXIT=0" farm_logs/build_v10.out || { say "build failed, abort"; exit 1; }
say "built blink_v10 md5 $(md5sum blink_v10 | cut -c1-8)"
for d in $R/gridrun/v10test $R/svomrun5/v10test $R/v6run/v10test; do mkdir -p $d/farm_logs; cp $R/gridrun/blink_v10 $d/blink.new && mv -f $d/blink.new $d/blink; done
mk() { # dir day args name
  cat > $1/job_$4.sh <<EOS
#!/bin/bash
cd $1
./blink search $2 $2 $3 > farm_logs/$4.out 2> farm_logs/$4.err; echo "EXIT=\$?" >> farm_logs/$4.out
EOS
  chmod +x $1/job_$4.sh; (cd $1 && hep_sub -g hxmt -mem 16000 job_$4.sh >> $log 2>&1)
}
mk $R/gridrun/v10test 2023-08-12 "--instrument grid03b" g03b
mk $R/gridrun/v10test 2024-03-15 "--instrument grid07" g07
mk $R/gridrun/v10test 2022-03-11 "--instrument grid04" g04
mk $R/gridrun/v10test 2020-12-08 "--instrument grid02" g02
mk $R/svomrun5/v10test 2024-07-03 "--instrument svom-grm" svom
mk $R/v6run/v10test 2017-10-20 "--workers 1 --worker 0" hxmt
say "6 regression jobs submitted"
while [ $(cat $R/gridrun/v10test/farm_logs/*.out $R/svomrun5/v10test/farm_logs/*.out $R/v6run/v10test/farm_logs/*.out 2>/dev/null | grep -c "EXIT=") -lt 6 ]; do sleep 60; done
say "regression jobs done; byte comparison:"
cmp_one() { # new old label
  if cmp -s "$1" "$2"; then say "  $3: IDENTICAL ($(wc -c < "$1") bytes)"; else say "  $3: DIFFERENT"; python3 - "$1" "$2" >> $log 2>&1 <<PY
import json, sys
a, b = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
print("    candidates new/old:", len(a), len(b))
for x, y in zip(a, b):
    if x != y:
        print("    first diff:", {k: (x.get(k), y.get(k)) for k in set(x) | set(y) if x.get(k) != y.get(k)}); break
PY
  fi
}
cmp_one $R/gridrun/v10test/data/GRID-03B/2023/08/20230812_signals.json $R/gridrun/data/GRID-03B/2023/08/20230812_signals.json "GRID-03B 2023-08-12"
cmp_one $R/gridrun/v10test/data/GRID-07/2024/03/20240315_signals.json $R/gridrun/data/GRID-07/2024/03/20240315_signals.json "GRID-07 2024-03-15"
cmp_one $R/gridrun/v10test/data/GRID-04/2022/03/20220311_signals.json $R/gridrun/data/GRID-04/2022/03/20220311_signals.json "GRID-04 2022-03-11"
cmp_one $R/gridrun/v10test/data/GRID-02/2020/12/20201208_signals.json $R/gridrun/data/GRID-02/2020/12/20201208_signals.json "GRID-02 2020-12-08"
cmp_one $R/svomrun5/v10test/data/SVOM_GRM/2024/07/20240703_signals.json $R/svomrun5/data/SVOM_GRM/2024/07/20240703_signals.json "SVOM 2024-07-03"
cmp_one $R/v6run/v10test/data/Insight-HXMT_HE/2017/10/20171020_signals.json $R/v6run/attcheck/data/Insight-HXMT_HE/2017/10/20171020_signals.json "HXMT 2017-10-20 (vs v8 attcheck)"
for f in $R/gridrun/v10test/farm_logs/*.err $R/svomrun5/v10test/farm_logs/*.err $R/v6run/v10test/farm_logs/*.err; do grep -l "panick" $f >> $log 2>/dev/null; done
say "DONE"
