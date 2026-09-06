"""检查天格试跑：账本里的小时状态、曝光是否等于归档 GTI 合计、候选、丢弃计数。"""
import json, glob, os, sys
from astropy.io import fits
R = "/scratchfs2/gecam/guohx/gridrun/data"; G = "/gecamfs/Exchange/GSDC/missions/GRID"
for sat, day in (("GRID-03B", "2023-08-12"), ("GRID-07", "2024-03-15"), ("GRID-04", "2022-03-11"), ("GRID-02", "2020-12-08")):
    y, m, d = day.split("-"); base = f"{R}/{sat}/{y}/{m}/{y}{m}{d}"
    if not os.path.exists(base + "_hours.json"):
        print(f"=== {sat} {day}: 还没有输出"); continue
    h = json.load(open(base + "_hours.json")); s = json.load(open(base + "_signals.json"))
    # 归档里这一天最高版本目录的 GTI 合计（只算落在当天 0–24h 内的部分）
    dd = f"{G}/{sat}/fits7/{y}/{m}/{d}"; v = sorted(os.listdir(dd))[-1]; gti = 0.0; npass = 0
    import datetime as dt
    ref = dt.datetime(2018, 1, 1, tzinfo=dt.timezone.utc); d0 = (dt.datetime(int(y), int(m), int(d), tzinfo=dt.timezone.utc) - ref).total_seconds(); d1 = d0 + 86400
    for f in glob.glob(f"{dd}/{v}/*.fits"):
        with fits.open(f) as hh:
            g = hh["GTI"].data; a, b = float(g["START"][0]), float(g["STOP"][0]); gti += max(0.0, min(b, d1) - max(a, d0)); npass += 1
    print(f"=== {sat} {day}  版本 {v}，{npass} 次过境，归档 GTI 合计 {gti/3600:.3f} h ===")
    print(f"  账本: searched={h['searched_hours']} excluded={h['excluded_hours']} {h['excluded_by_reason']} searched_seconds={h['searched_seconds']:.0f} ({h['searched_seconds']/3600:.3f} h)  候选 {h['n_signals']}")
    print(f"  曝光核对: 账本 {h['searched_seconds']:.1f} s vs GTI {gti:.1f} s → 差 {h['searched_seconds']-gti:+.1f} s")
    tot = {}
    for x in h["hours"]:
        if x.get("status") == "searched":
            for k, val in (x.get("metrics") or {}).items(): tot[k] = tot.get(k, 0) + val
    print("  诊断合计:", {k: int(val) for k, val in tot.items()})
    sig = sorted(s, key=lambda c: c["false_positive_per_year"])
    print(f"  候选 {len(s)}，fa<=1e-5: {sum(c['false_positive_per_year']<=1e-5 for c in s)}，最显著 3 个:")
    for c in sig[:3]:
        print(f"    {c['start'][11:23]} count={c['count']:3d} mean={c['mean']:6.2f} fa={c['false_positive_per_year']:.2e} lon={c['position']['longitude']:7.1f} lat={c['position']['latitude']:6.1f} dur={c['bin_size_best']*1e6:.0f}us")
    err = f"/scratchfs2/gecam/guohx/gridrun/farm_logs/test_{sat.lower().replace('-','')}.err"
    for line in open(err):
        if "Elapsed" in line or "Maximum resident" in line: print("  ", line.strip())
