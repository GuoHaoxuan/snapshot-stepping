"""针对性重跑 vs 权威目录：逐天按起始时间（±1 ms）配对，找只在新结果里出现且没有 attitude 字段的候选。"""
import json, glob, os, sys
from datetime import datetime, timezone
def t(iso):
    b = iso.rstrip("Z"); h, f = b.split("."); return datetime.strptime(h + "." + (f + "000000")[:6], "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc).timestamp()
for tag, new_root, old_root, sub in (("HXMT", "/scratchfs2/gecam/guohx/v6run/attcheck/data", "/scratchfs2/gecam/guohx/v6run/data", "Insight-HXMT_HE"), ("SVOM", "/scratchfs2/gecam/guohx/svomrun5/attcheck/data", "/scratchfs2/gecam/guohx/svomrun5/data", "SVOM_GRM")):
    days = [l.strip() for l in open(os.path.dirname(new_root) + "/days.txt") if l.strip()]
    done = 0; n_new = n_old = 0; only_new = []; only_old = []; no_att = 0; without = 0; noeph_new = 0; single = 0
    for d in days:
        y, m, dd = d.split("-"); rel = f"{sub}/{y}/{m}/{y}{m}{dd}"
        fn, fo = f"{new_root}/{rel}_signals.json", f"{old_root}/{rel}_signals.json"
        if not os.path.exists(fn): continue
        done += 1
        new = json.load(open(fn)); old = json.load(open(fo)) if os.path.exists(fo) else []
        for x in json.load(open(f"{new_root}/{rel}_hours.json"))["hours"]:
            mt = x.get("metrics") or {}; without += int(mt.get("without_attitude", 0)); noeph_new += int(mt.get("dropped_no_ephemeris", 0)); single += int(mt.get("dropped_single_detector", 0))
        n_new += len(new); n_old += len(old)
        to = [t(c["start"]) for c in old]; tn = [t(c["start"]) for c in new]
        for c, tc in zip(new, tn):
            if not any(abs(tc - x) < 1e-3 for x in to):
                only_new.append((d, c["start"][:23], c["false_positive_per_year"], "attitude" in c))
                if "attitude" not in c: no_att += 1
        for c, tc in zip(old, to):
            if not any(abs(tc - x) < 1e-3 for x in tn): only_old.append((d, c["start"][:23], c["false_positive_per_year"]))
    print(f"=== {tag}: 重跑 {done}/{len(days)} 天；旧 {n_old} 个候选，新 {n_new} 个；新结果 dropped_no_ephemeris {noeph_new}、without_attitude {without}、dropped_single_detector {single}")
    print(f"  只在新结果里: {len(only_new)}（其中无 attitude 字段 {no_att}）；只在旧结果里: {len(only_old)}")
    for r in only_new[:12]: print("   新:", r)
    for r in only_old[:6]: print("   旧:", r)
