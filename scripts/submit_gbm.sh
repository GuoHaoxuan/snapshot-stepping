#!/bin/bash
# 提交 Fermi/GBM 全量 TGF 搜索：100 个 worker，%{ProcId}(0..99) 取 days.txt 的行。
# 重跑安全：已产出且不比源数据旧的天会被跳过，中途失败直接重跑本脚本补齐。
# job_gbm.sh 与 days.txt 都在运行目录 /scratchfs2/gecam/guohx/gbmrun 下。
#
# days.txt 由 gbm_days.sh 生成：BGO 逐小时 TTE 齐全的天即可入选，NaI 有没有
# 由 Chunk 自行判断（2017-10 之前没有，那些天按单组搜）。

export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH

# 2863 天 / 100 worker ≈ 29 天每作业。GBM 一天要读 336 个压缩 TTE（约 7.2 GB），
# 比 SVOM 重一个量级，内存也给足：一小时十几路事例流合起来有几千万条。
hep_sub -g hxmt -mem 16384 job_gbm.sh -argu "%{ProcId}" -n 100
