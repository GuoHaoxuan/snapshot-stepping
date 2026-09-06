#!/bin/bash
cd /scratchfs2/gecam/guohx/gridrun
python3 grid_features.py > farm_logs/feat.out 2> farm_logs/feat.err; echo "EXIT=$?" >> farm_logs/feat.out
