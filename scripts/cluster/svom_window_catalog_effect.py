"""窗长上限对目录判选（fa≤1e-5 ∪ (fa≤1 且关联)）的影响。"""
import json
from datetime import datetime, timezone


def load(path, until="2025-01-01"):
    out = []
    for r in json.load(open(path)):
        s, li = r["signal"], r["lightning"]
        if s["start"] >= until:
            continue
        out.append((s["false_positive_per_year"], bool(li.get("associated")),
                    bool(li.get("in_coverage", True)), li.get("coincidence_probability") or 0.0))
    return out


for lab, path in (("v6 (1 ms)", "scratch/tgfs_svom_v6.json"), ("v7 (5 ms)", "scratch/tgfs_svom_v7.json")):
    rows = load(path)
    direct = [r for r in rows if r[0] <= 1e-5]
    rescue = [r for r in rows if 1e-5 < r[0] <= 1.0 and r[1]]
    false_rescue = sum(r[3] for r in rows if 1e-5 < r[0] <= 1.0 and r[2])
    print("%-10s 池 %5d  目录 %4d = 直接 %4d + 仅关联 %3d   误救期望 %.2f"
          % (lab, len(rows), len(direct) + len(rescue), len(direct), len(rescue), false_rescue))
