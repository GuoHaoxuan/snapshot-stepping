#!/bin/bash
# HXMT 单探头占比审计：构建 → 核对 → 在 v6 显著候选上跑 acd-audit
export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH
cd /scratchfs2/gecam/guohx/gridrun
log=/scratchfs2/gecam/guohx/gridrun/farm_logs/chain_audit.log   # 绝对路径：脚本后面 cd 过，相对路径会写空
say() { echo "$(date +%H:%M) $*" >> $log; }
say "chain start on $(hostname)"
rm -f farm_logs/build_audit.out
hep_sub -g hxmt -mem 12000 -cpu 8 job_build_audit.sh >> $log 2>&1
while ! grep -q BUILD_EXIT farm_logs/build_audit.out 2>/dev/null; do sleep 60; done
say "build: $(grep -h 'BUILD_EXIT\|binary copied' farm_logs/build_audit.out | tr '\n' ' ')"
grep -q "BUILD_EXIT=0" farm_logs/build_audit.out || { say "build failed, abort"; exit 1; }
strings blink_audit | grep -q det_share_max || { say "blink_audit lacks det_share_max, abort"; exit 1; }
cd /scratchfs2/gecam/guohx/v6run/audit_det
cp /scratchfs2/gecam/guohx/gridrun/blink_audit blink.new && mv -f blink.new blink
rm -f audit.out
hep_sub -g hxmt -mem 16000 job_audit.sh >> $log 2>&1
say "audit submitted"
while ! grep -q "EXIT=" audit.out 2>/dev/null; do sleep 300; say "audit: $(tail -1 audit.err 2>/dev/null | cut -c1-80)"; done
say "audit: $(cat audit.out | tr '\n' ' ') $(tail -1 audit.err)"
say "DONE"
