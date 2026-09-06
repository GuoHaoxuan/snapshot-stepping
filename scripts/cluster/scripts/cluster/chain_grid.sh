#!/bin/bash
# v6 链：构建 → 核对二进制内容 → 部署（临时名 + mv + md5 核对）→ 清掉 v4（与 v3 同一二进制的产物）→ 全量 → 汇总 → 特征 + 闪电关联
export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH
cd /scratchfs2/gecam/guohx/gridrun
log=farm_logs/chain_v8.log
say() { echo "$(date +%H:%M) $*" >> $log; }
say "chain start on $(hostname)"
rm -f farm_logs/build_v8.out
hep_sub -g hxmt -mem 12000 -cpu 8 job_build_v8.sh >> $log 2>&1
while ! grep -q BUILD_EXIT farm_logs/build_v8.out 2>/dev/null; do sleep 60; done
say "build: $(grep -h 'BUILD_EXIT\|binary copied' farm_logs/build_v8.out | tr '\n' ' ')"
grep -q "BUILD_EXIT=0" farm_logs/build_v8.out || { say "build failed, abort"; exit 1; }
{ [ -x blink_v8 ] && strings blink_v8 | grep -q dropped_single_detector && strings blink_v8 | grep -q without_attitude; } || { say "blink_v8 missing or stale, abort"; exit 1; }
cp blink_v8 blink.new 2>>$log && mv -f blink.new blink 2>>$log
[ "$(md5sum blink | cut -d' ' -f1)" = "$(md5sum blink_v8 | cut -d' ' -f1)" ] || { say "deploy failed: gridrun/blink md5 != blink_v8, abort"; exit 1; }
say "deployed gridrun/blink = blink_v8 (md5 $(md5sum blink | cut -c1-8))"
n=$(find data -type f | wc -l)
tar cf data_v7.tar data && [ $(tar tf data_v7.tar | grep -c -v '/$') -eq $n ] && rm -rf data && say "v7 archived: $n files" || { say "v5 archive mismatch, abort"; exit 1; }
mv -f tgfs_grid02.json tgfs_grid02_v7.json; mv -f tgfs_grid03b.json tgfs_grid03b_v7.json; mv -f tgfs_grid04.json tgfs_grid04_v7.json; mv -f tgfs_grid07.json tgfs_grid07_v7.json; cp features_sig.csv features_sig_v7.csv
rm -f farm_logs/worker_*.err job_grid.sh.out.* job_grid.sh.err.*
hep_sub -g hxmt -mem 8000 job_grid.sh -argu "%{ProcId}" -n 100 >> $log 2>&1
say "v8 submitted"
while :; do
  d=$(ls data/*/*/*/*_signals.json 2>/dev/null | wc -l); q=$(hep_q -u guohx 2>/dev/null | grep "job_grid.sh" | grep -vc " X ")
  say "v8 days=$d/1586 running=$q"
  [ "$d" -ge 1586 ] && break
  [ "$q" -eq 0 ] && { say "workers gone but days=$d"; break; }
  sleep 300
done
say "v8 run finished; panic/quota lines: $(cat farm_logs/worker_*.err 2>/dev/null | grep -c 'panick\|quota')"
say "null attitude files: $(grep -l ': null' data/GRID-*/*/*/*_signals.json 2>/dev/null | wc -l)"
python3 grid_summary.py >> $log 2>&1
rm -f farm_logs/feat.out wwlln_grid.log
hep_sub -g hxmt -mem 8000 -cpu 8 job_wwlln_v3.sh >> $log 2>&1
hep_sub -g hxmt -mem 8000 job_feat.sh >> $log 2>&1
while :; do a=$(grep -c "EXIT=" farm_logs/feat.out 2>/dev/null); b=$(grep -c "^end" wwlln_grid.log 2>/dev/null); [ "${a:-0}" -ge 1 ] && [ "${b:-0}" -ge 1 ] && break; sleep 60; done
say "feat: $(cat farm_logs/feat.out | tr '\n' ' ')"
say "wwlln:"; grep -v "^filter: [0-9]*/[0-9]* *$" wwlln_grid.log >> $log
ls -la tgfs_grid*.json >> $log 2>&1
say "DONE"
