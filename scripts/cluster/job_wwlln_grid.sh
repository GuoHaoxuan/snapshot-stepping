#!/bin/bash
# 四颗星逐个关联；`blink wwlln` 把 tgfs.json 写在当前目录，跑完一颗改名留档。
cd /scratchfs2/gecam/guohx/gridrun
export WWLLN_DB_PATH=/gecamfs/Exchange/GSDC/missions/AEfiles/WWLLN.db
echo "start $(date) node $(hostname) nproc $(nproc)" > wwlln_grid.log
for sat in grid02 grid03b grid04 grid07; do
  ./blink wwlln --instrument $sat >> wwlln_grid.log 2>&1
  echo "$sat EXIT=$?" >> wwlln_grid.log
  mv tgfs.json tgfs_$sat.json
done
echo "end $(date)" >> wwlln_grid.log
