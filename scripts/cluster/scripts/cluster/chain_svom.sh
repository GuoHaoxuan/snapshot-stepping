#!/bin/bash
# SVOM v4 = v3 + 单路 GRD 占比 ≤ 0.8。等 gridrun 的 v5 二进制部署好、v5 全量跑完再占农场。
export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH
cd /scratchfs2/gecam/guohx/svomrun5
log=farm_logs/chain_svom5.log
say() { echo "$(date +%H:%M) $*" >> $log; }
say "chain start on $(hostname)"
g=/scratchfs2/gecam/guohx/gridrun/farm_logs/chain_v7.log
while ! grep -q "deployed gridrun/blink\|abort" $g 2>/dev/null; do sleep 60; done
grep -q abort $g && { say "grid v5 chain aborted, abort"; exit 1; }
cp /scratchfs2/gecam/guohx/gridrun/blink_v7 blink.new && mv -f blink.new blink
[ "$(md5sum blink | cut -d' ' -f1)" = "$(md5sum /scratchfs2/gecam/guohx/gridrun/blink_v7 | cut -d' ' -f1)" ] || { say "deploy failed, abort"; exit 1; }
strings blink | grep -q dropped_single_detector || { say "binary lacks the rule, abort"; exit 1; }
say "deployed blink = gridrun/blink_v7 (md5 $(md5sum blink | cut -c1-8))"
while ! grep -q "v7 run finished\|abort" $g 2>/dev/null; do sleep 300; done
hep_sub -g hxmt -mem 8000 job_svom.sh -argu "%{ProcId}" -n 100 >> $log 2>&1
say "svom v5 submitted"
while :; do
  d=$(ls data/SVOM_GRM/*/*/*_signals.json 2>/dev/null | wc -l); q=$(hep_q -u guohx 2>/dev/null | grep "job_svom.sh" | grep -vc " X ")
  say "svom v5 days=$d/792 running=$q"
  [ "$d" -ge 792 ] && break
  [ "$q" -eq 0 ] && { say "workers gone but days=$d"; break; }
  sleep 300
done
say "svom v5 run finished; panic/quota lines: $(cat job_svom.sh.err.* 2>/dev/null | grep -c 'panick\|quota')"
python3 - >> $log 2>&1 <<PY
import json, glob, collections
for tag, root in (("v3", "/scratchfs2/gecam/guohx/svomrun3/data"), ("v4", "/scratchfs2/gecam/guohx/svomrun5/data")):
    sig = []; tot = collections.Counter()
    for f in sorted(glob.glob(root + "/SVOM_GRM/*/*/*_signals.json")): sig += json.load(open(f))
    for f in sorted(glob.glob(root + "/SVOM_GRM/*/*/*_hours.json")):
        for x in json.load(open(f))["hours"]:
            for k, v in (x.get("metrics") or {}).items():
                if k.startswith("dropped"): tot[k] += int(v)
    print(tag, "候选", len(sig), "显著(fa<=1e-5)", sum(c["false_positive_per_year"] <= 1e-5 for c in sig), "丢弃", dict(tot))
PY
rm -f wwlln_svom.log; hep_sub -g hxmt -mem 8000 -cpu 8 job_wwlln.sh >> $log 2>&1
while ! grep -q "^EXIT" wwlln_svom.log 2>/dev/null; do sleep 60; done
grep -v "^filter: [0-9]*/[0-9]* *$" wwlln_svom.log >> $log
say "DONE"
