#!/bin/bash
cd /scratchfs2/gecam/guohx/v6run
log=catalog.log
: > "$log"
./blink catalog tgfs.json -o catalog_v6.csv >> "$log" 2>&1
echo "CATALOG_EXIT=$? end $(date)" >> "$log"
