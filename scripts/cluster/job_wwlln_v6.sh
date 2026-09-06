#!/bin/bash
cd /scratchfs2/gecam/guohx/v6run
export WWLLN_DB_PATH=/gecamfs/Exchange/GSDC/missions/AEfiles/WWLLN.db
echo "start $(date) node $(hostname) nproc $(nproc)" > wwlln_farm.log
./blink wwlln >> wwlln_farm.log 2>&1
echo "WWLLN_EXIT=$? end $(date)" >> wwlln_farm.log
