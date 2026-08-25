#!/bin/bash
# 提交 SVOM/GRM 全量 TGF 搜索：100 个 worker，%{ProcId}(0..99) 取 days.txt 的行。
# 重跑安全：已产出且不比源数据旧的天会被跳过，中途失败直接重跑本脚本补齐。
# job_svom.sh 与 days.txt 都在运行目录 /scratchfs2/gecam/guohx/svomrun 下。

export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH

# 792 天 / 100 worker ≈ 8 天每作业，实测 ~4 分钟每天 ≈ 35 分钟，远低于墙钟上限。
hep_sub -g hxmt -mem 8192 job_svom.sh -argu "%{ProcId}" -n 100
