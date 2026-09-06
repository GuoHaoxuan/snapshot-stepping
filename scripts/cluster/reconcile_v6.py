import json, glob, os, collections
v5 = "/scratchfs2/gecam/guohx/v5run/data/Insight-HXMT_HE"
v6 = "/scratchfs2/gecam/guohx/v6run/data/Insight-HXMT_HE"
d5 = {os.path.basename(p)[:8]: p for p in glob.glob(v5 + "/*/*/*_signals.json")}
d6 = {os.path.basename(p)[:8]: p for p in glob.glob(v6 + "/*/*/*_signals.json")}
only5, only6 = sorted(set(d5) - set(d6)), sorted(set(d6) - set(d5))
print("days: v5", len(d5), "v6", len(d6))
print("only-v5:", only5)
print("only-v6:", len(only6), only6)
n_same = 0
diff_days = []
n_cand = n_sig = n_acd_missing = 0
for day in sorted(d6):
    s6 = json.load(open(d6[day]))
    n_cand += len(s6)
    n_acd_missing += sum(1 for x in s6 if "acd" not in x)
    n_sig += sum(1 for x in s6 if x["false_positive_per_year"] < 1e-8)
    if day in d5:
        s5 = json.load(open(d5[day]))
        if [x["start"] for x in s5] == [x["start"] for x in s6]:
            n_same += 1
        else:
            diff_days.append((day, len(s5), len(s6)))
print("common days:", len(set(d5) & set(d6)), "identical:", n_same, "differing:", len(diff_days))
print("diff days:", diff_days[:30])
print("v6 candidates:", n_cand, "significant(<1e-8):", n_sig, "missing-acd:", n_acd_missing)
# 逐小时账本汇总
searched_h = excluded_h = 0
searched_s = 0.0
reasons = collections.Counter()
for hp in glob.glob(v6 + "/*/*/*_hours.json"):
    h = json.load(open(hp))
    searched_h += h["searched_hours"]
    excluded_h += h["excluded_hours"]
    searched_s += h["searched_seconds"]
    for k, v in h.get("excluded_by_reason", {}).items():
        reasons[k] += v
print("hours: searched", searched_h, "excluded", excluded_h,
      "exposure_days %.1f" % (searched_s / 86400.0))
print("excluded_by_reason:", dict(reasons))
