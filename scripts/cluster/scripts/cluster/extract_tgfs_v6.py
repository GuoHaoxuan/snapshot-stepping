#!/usr/bin/env python3
"""Stream-scan v6 tgfs.json -> full-pool candidate table with train fields.

Same line-scanner contract as sig_train_v6.py; one row per candidate:
day,fpy,assoc,coinc,neighbors_10min,is_train (matches the load_catalog
format of the fp-distribution scripts, plus the two train columns).
"""
SRC = "/scratchfs2/gecam/guohx/v6run/tgfs.json"
OUT = "tgfs_v6.csv"
n = 0
start = fpy = assoc = coinc = nb = None
with open(SRC) as f, open(OUT, "w") as out:
    out.write("day,fpy,assoc,coinc,neighbors_10min,is_train\n")
    for line in f:
        line = line.strip()
        if line.startswith('"start"'):
            start = line.split('"')[3]
        elif line.startswith('"false_positive_per_year"'):
            fpy = line.split(":")[1].strip().rstrip(",")
        elif line.startswith('"associated"'):
            assoc = "true" in line
        elif line.startswith('"coincidence_probability"'):
            coinc = line.split(":")[1].strip().rstrip(",")
        elif line.startswith('"neighbors_10min"'):
            nb = line.split(":")[1].strip().rstrip(",")
        elif line.startswith('"is_train"'):
            n += 1
            day = start[:10].replace("-", "")
            out.write(f"{day},{fpy},{int(assoc)},{coinc},{nb},{int('true' in line)}\n")
            start = fpy = assoc = coinc = nb = None
print(f"wrote {n} candidates -> {OUT}")
