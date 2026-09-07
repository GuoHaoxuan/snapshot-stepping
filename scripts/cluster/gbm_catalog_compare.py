"""把我们在 GBM 上搜出的候选与已发表的第二版 GBM TGF 目录逐个匹配，量完备性与纯度。

这是 HXMT 与 SVOM 都没有的外部真值。目录覆盖 2008-07-11 至 2016-07-31，所以只能对
那个区间内的年份做检验（2019 年不在覆盖内）。

MET 都是 Fermi MET（2001-01-01 起的连续 TT 秒），两边同一基准，直接按 MET 匹配。

用法: python3 gbm_catalog_compare.py <year> <run_data_dir> <catalog_dir> <out.csv>
"""
import csv, glob, json, os, sys
from datetime import datetime, timezone
import numpy as np

# 2001 年之后的闰秒（UTC 日期）；MET = (UTC − 2001-01-01) + 已过闰秒数
LEAPS = ("2006-01-01", "2009-01-01", "2012-07-01", "2015-07-01", "2017-01-01")
REF = datetime(2001, 1, 1, tzinfo=timezone.utc)
TOL_S = 0.01          # 目录给的是 TGF 的时刻，我们给的是候选窗起点，10 ms 足够宽


def met_of(iso):
    b = iso.rstrip("Z"); h, _, f = b.partition(".")
    t = datetime.strptime(h + "." + (f + "000000")[:6], "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)
    n = sum(1 for l in LEAPS if t >= datetime.fromisoformat(l).replace(tzinfo=timezone.utc))
    return (t - REF).total_seconds() + n


def main(year, run_dir, cat_dir, out):
    cat = []
    for r in csv.reader(open(os.path.join(cat_dir, "gbm_tgf_catalog_offline.csv"))):
        if not r or r[0].startswith("#"): continue
        if r[6].strip()[:4] != year: continue
        cat.append({"id": r[0].strip(), "met": float(r[1]), "nai": float(r[5]), "b0": float(r[3]), "b1": float(r[4]),
                    "width_ms": float(r[8]), "lon": float(r[10]), "lat": float(r[11]), "trig": r[14].strip()})
    wwlln = {r[0].strip() for r in csv.reader(open(os.path.join(cat_dir, "gbm_tgf_catalog_wwlln.csv"))) if r and not r[0].startswith("#")}
    print(f"目录 {year} 年：{len(cat)} 个 TGF，其中带 WWLLN 关联 {sum(1 for c in cat if c['id'] in wwlln)} 个")

    ours = []
    for f in sorted(glob.glob(os.path.join(run_dir, "Fermi_GBM", year, "*", "*_signals.json"))):
        for s in json.load(open(f)):
            ours.append({"met": met_of(s["start"]), "fa": s["false_positive_per_year"], "count": s["count"],
                         "dur": s["bin_size_best"], "lon": s["position"]["longitude"], "lat": s["position"]["latitude"]})
    days = len(glob.glob(os.path.join(run_dir, "Fermi_GBM", year, "*", "*_hours.json")))
    print(f"我们的搜索：{days} 天，{len(ours)} 个候选，显著（fa≤1e-5）{sum(1 for o in ours if o['fa'] <= 1e-5)} 个")
    if not ours: print("没有候选，先跑搜索"); return

    om = np.array([o["met"] for o in ours]); order = np.argsort(om); om = om[order]; ours = [ours[i] for i in order]
    # 目录里的每个 TGF：我们找到没有、fa 多少
    rows = []
    for c in cat:
        i = np.searchsorted(om, c["met"]); best = None
        for j in (i - 1, i, i + 1):
            if 0 <= j < len(om) and abs(om[j] - c["met"]) <= TOL_S:
                if best is None or ours[j]["fa"] < ours[best]["fa"]: best = j
        rows.append({"id": c["id"], "met": c["met"], "nai": c["nai"], "bgo": c["b0"] + c["b1"], "width_ms": c["width_ms"],
                     "wwlln": int(c["id"] in wwlln), "triggered": int(c["trig"] not in ("", "NULL")),
                     "found": int(best is not None), "our_fa": ours[best]["fa"] if best is not None else "",
                     "our_count": ours[best]["count"] if best is not None else "", "dt_ms": (om[best] - c["met"]) * 1e3 if best is not None else ""})
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    found = np.array([r["found"] for r in rows], bool)
    fa = np.array([r["our_fa"] if r["our_fa"] != "" else np.inf for r in rows], float)
    tot = np.array([r["nai"] + r["bgo"] for r in rows])
    print(f"完备性：目录 {len(rows)} 个里我们找到 {found.sum()}（{100*found.mean():.1f}%），其中达显著（fa≤1e-5）{int((fa <= 1e-5).sum())}")
    for lo, hi in ((0, 50), (50, 100), (100, 200), (200, 500), (500, 1e9)):
        m = (tot >= lo) & (tot < hi)
        if m.sum(): print(f"  目录计数 {lo:4.0f}–{hi:<6.0f} n={m.sum():4d} 找到 {100*found[m].mean():5.1f}% 达显著 {100*(fa[m] <= 1e-5).mean():5.1f}%")
    trig = np.array([r["triggered"] for r in rows], bool); ww = np.array([r["wwlln"] for r in rows], bool)
    print(f"  星上触发的 {trig.sum()} 个里找到 {100*found[trig].mean():.1f}%；带 WWLLN 的 {ww.sum()} 个里找到 {100*found[ww].mean():.1f}%")
    # 反向：我们的显著候选有多少不在目录里
    cm = np.array(sorted(c["met"] for c in cat))
    extra = 0; extra_sig = 0
    for o in ours:
        i = np.searchsorted(cm, o["met"])
        hit = any(0 <= j < len(cm) and abs(cm[j] - o["met"]) <= TOL_S for j in (i - 1, i, i + 1))
        if not hit:
            extra += 1
            if o["fa"] <= 1e-5: extra_sig += 1
    print(f"  我们有而目录没有：{extra} 个候选，其中显著 {extra_sig} 个（目录本身不完备，这些要另判）")
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
