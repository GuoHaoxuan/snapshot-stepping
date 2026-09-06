import json, glob, collections, numpy as np
R = "/scratchfs2/gecam/guohx/gridrun/data"
for sat in ("GRID-02", "GRID-03B", "GRID-04", "GRID-07"):
    sig = []; tot = collections.Counter(); corrupt = []; per_day = collections.Counter()
    for f in sorted(glob.glob(f"{R}/{sat}/*/*/*_signals.json")):
        s = json.load(open(f)); sig += s; per_day[f[-21:-13]] = len(s)
    for f in sorted(glob.glob(f"{R}/{sat}/*/*/*_hours.json")):
        h = json.load(open(f))
        for x in h["hours"]:
            for k, v in (x.get("metrics") or {}).items():
                if k.startswith("dropped"): tot[k] += int(v)
            if x.get("reason") == "corrupt_data": corrupt.append((f[-21:-13], x.get("hour"), (x.get("detail") or "")[:90]))
    fa = np.array([c["false_positive_per_year"] for c in sig]); lat = np.array([c["position"]["latitude"] for c in sig])
    cnt = np.array([c["count"] for c in sig]); mean = np.array([c["mean"] for c in sig]); dur = np.array([c["bin_size_best"] for c in sig]) * 1e6
    rate = mean / (dur * 1e-6)
    print(f"=== {sat}: {len(sig)} 候选, 丢弃 {dict(tot)}, corrupt_data {len(corrupt)} 小时 ===")
    for c in corrupt[:3]: print("   corrupt:", c)
    if not len(sig): continue
    print(f"  |lat| 直方(0..60,10°): 全部 {np.histogram(np.abs(lat), bins=range(0,70,10))[0].tolist()}  fa<=1e-5 {np.histogram(np.abs(lat[fa<=1e-5]), bins=range(0,70,10))[0].tolist()}")
    print(f"  候选所在处本底率 mean/dur (c/s): 中位 {np.median(rate):.0f}, p90 {np.percentile(rate,90):.0f}; fa<=1e-5 的: 中位 {np.median(rate[fa<=1e-5]) if (fa<=1e-5).any() else 0:.0f}")
    print(f"  count 中位 {np.median(cnt):.0f}, dur 中位 {np.median(dur):.0f} µs; fa<=1e-5 的 count 中位 {np.median(cnt[fa<=1e-5]) if (fa<=1e-5).any() else 0:.0f}")
    top = per_day.most_common(5); print(f"  候选最多的 5 天: {top}  （{len(per_day)} 天有候选, 中位 {np.median(list(per_day.values())):.0f}/天）")
    order = np.argsort(fa)[:6]
    print("  最显著 6 个:")
    for i in order:
        c = sig[i]; print(f"    {c['start'][:23]} count={c['count']:3d} mean={c['mean']:7.2f} rate={rate[i]:7.0f} fa={fa[i]:.1e} lon={c['position']['longitude']:7.1f} lat={c['position']['latitude']:6.1f} dur={dur[i]:.0f}us")
