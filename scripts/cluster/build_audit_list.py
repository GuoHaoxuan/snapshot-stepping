#!/usr/bin/env python3
"""Stream-scan v5 tgfs.json (pretty-printed) -> certified-sample audit list.

Keeps rows with associated==true (certified TGF pool) or fpy < 1e-5
(paper-significant tier: REP storm-day + 2025-09-30 check populations).
Line scanner, O(1) memory -- the 1.4 GB json never fully loads.
"""
SRC = "/scratchfs2/gecam/guohx/v5run/tgfs.json"
OUT = "certified_v5.csv"
n_kept = n_seen = 0
start = stop = fpy = None
with open(SRC) as f, open(OUT, "w") as out:
    out.write("date,start,stop,false_positive_per_year,assoc\n")
    for line in f:
        line = line.strip()
        if line.startswith('"start"'):
            start = line.split('"')[3]
        elif line.startswith('"stop"'):
            stop = line.split('"')[3]
        elif line.startswith('"false_positive_per_year"'):
            fpy = float(line.split(":")[1].rstrip(","))
        elif line.startswith('"associated"'):
            n_seen += 1
            assoc = "true" in line
            if assoc or fpy < 1e-5:
                out.write(f"{start[:10].replace('-', '')},{start},{stop},{fpy},{int(assoc)}\n")
                n_kept += 1
            start = stop = fpy = None
print(f"scanned {n_seen} candidates, kept {n_kept}")
