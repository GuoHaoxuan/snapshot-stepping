#!/bin/bash
cd /scratchfs2/gecam/guohx/svomrun5
export WWLLN_DB_PATH=/gecamfs/Exchange/GSDC/missions/AEfiles/WWLLN.db
echo "start $(date) node $(hostname) nproc $(nproc)" > wwlln_svom.log
./blink wwlln --instrument svom-grm >> wwlln_svom.log 2>&1
echo "EXIT=$? end $(date)" >> wwlln_svom.log
