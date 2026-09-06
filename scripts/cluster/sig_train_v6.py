#!/usr/bin/env python3
"""Stream-scan v6 tgfs.json -> train fields for the certified/significant subset.

Keeps rows with associated==true or fpy < 1e-5 (superset of the fpy<1e-8
significant tier). Same line-scanner contract as build_audit_list.py; the
train block serializes after lightning, so emit on the record-final is_train.
"""
SRC = "/scratchfs2/gecam/guohx/v6run/tgfs.json"
OUT = "certified_v6_train.csv"
n_seen = n_kept = 0
n_sig = n_sig_train = n_sig_assoc = n_sig_assoc_train = 0
start = stop = fpy = assoc = nb = None
with open(SRC) as f, open(OUT, "w") as out:
    out.write("date,start,stop,false_positive_per_year,assoc,neighbors_10min,is_train\n")
    for line in f:
        line = line.strip()
        if line.startswith('"start"'):
            start = line.split('"')[3]
        elif line.startswith('"stop"'):
            stop = line.split('"')[3]
        elif line.startswith('"false_positive_per_year"'):
            fpy = float(line.split(":")[1].rstrip(","))
        elif line.startswith('"associated"'):
            assoc = "true" in line
        elif line.startswith('"neighbors_10min"'):
            nb = int(line.split(":")[1].rstrip(","))
        elif line.startswith('"is_train"'):
            n_seen += 1
            is_train = "true" in line
            if fpy < 1e-8:
                n_sig += 1
                n_sig_train += is_train
                if assoc:
                    n_sig_assoc += 1
                    n_sig_assoc_train += is_train
            if assoc or fpy < 1e-5:
                day = start[:10].replace("-", "")
                out.write(f"{day},{start},{stop},{fpy},{int(assoc)},{nb},{int(is_train)}\n")
                n_kept += 1
            start = stop = fpy = assoc = nb = None
print(f"scanned {n_seen} candidates, kept {n_kept}")
print(f"significant(<1e-8): {n_sig}  train-flagged: {n_sig_train}  "
      f"kept after cut: {n_sig - n_sig_train}")
print(f"significant assoc: {n_sig_assoc}  of which train-flagged: {n_sig_assoc_train}")
