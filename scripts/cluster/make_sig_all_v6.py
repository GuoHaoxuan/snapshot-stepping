import json, glob, csv, os
root = "/scratchfs2/gecam/guohx/v6run/data/Insight-HXMT_HE"
out = "/scratchfs2/gecam/guohx/v6run/sig_all_v6.csv"
rows = []
for sp in sorted(glob.glob(root + "/*/*/*_signals.json")):
    date = os.path.basename(sp)[:8]
    for s in json.load(open(sp)):
        if s["false_positive_per_year"] < 1e-8:
            a = s.get("acd") or {}
            rows.append([
                date, s["start"], s["stop"], s["count"], s["mean"], s["sf"],
                s["false_positive_per_year"],
                s["position"]["latitude"], s["position"]["longitude"],
                s["position"]["altitude"],
                a.get("n", ""), a.get("n_acd", ""), a.get("n_acd_multi", ""),
                a.get("n_bg", ""), a.get("n_acd_bg", ""),
            ])
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["date","start","stop","count","mean","sf",
                "false_positive_per_year","lat","lon","alt",
                "n","n_acd","n_acd_multi","n_bg","n_acd_bg"])
    w.writerows(rows)
print("rows:", len(rows))
