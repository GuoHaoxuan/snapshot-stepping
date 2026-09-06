"""空洞占比 vs 计数率 / 纬度。向量化：F(x)=累计空洞长度，秒内空洞 = F(s+1)-F(s)。"""
from astropy.io import fits
import glob, os, json, sys, numpy as np, datetime as dt
G = "/gecamfs/Exchange/GSDC/missions/GRID"; R = "/scratchfs2/gecam/guohx/gridrun/data"
REF = dt.datetime(2018, 1, 1, tzinfo=dt.timezone.utc)
def met(iso):
    b = iso.rstrip("Z"); h, f = b.split("."); b = h + "." + (f + "000000")[:6]
    return (dt.datetime.strptime(b, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=dt.timezone.utc) - REF).total_seconds()
GAP = 0.003
for sat, day, tag in (("GRID-03B", "2023-08-12", "20230812"), ("GRID-02", "2020-12-08", "20201208"), ("GRID-04", "2022-03-11", "20220311")):
    y, m, d0 = day.split("-"); vdir = sorted(os.listdir(f"{G}/{sat}/fits7/{y}/{m}/{d0}"))[-1]
    pd = sorted(glob.glob(f"{G}/{sat}/fits8/{y}/{m}/{d0}/posatt*"))[-1]
    S = []
    for f in sorted(glob.glob(f"{G}/{sat}/fits7/{y}/{m}/{d0}/{vdir}/*.fits")):
        with fits.open(f) as h:
            t = np.sort(np.concatenate([np.asarray(h[f"EVENTS{k}"].data["TIME"], float)[np.asarray(h[f"EVENTS{k}"].data["EVT_TYPE"]) == 1] for k in range(4)]))
            g = h["GTI"].data; s0, s1 = float(g["START"][0]), float(g["STOP"][0])
        if len(t) < 2: continue
        pf = glob.glob(f"{pd}/*{os.path.basename(f)[8:29]}*.fits"); lat_t = lat_v = None
        if pf:
            with fits.open(pf[0]) as h: lat_t = np.asarray(h[1].data["TIME"], float); lat_v = np.asarray(h[1].data["Latitude"], float)
        d = np.diff(t); gs = t[:-1][d > GAP]; ge = gs + d[d > GAP]
        # 断点：每个空洞的起止；F 在空洞内线性上升
        bp = np.empty(2 * len(gs)); bp[0::2] = gs; bp[1::2] = ge
        fv = np.empty_like(bp); cum = np.concatenate([[0.0], np.cumsum(ge - gs)]); fv[0::2] = cum[:-1]; fv[1::2] = cum[1:]
        secs = np.arange(np.floor(s0), np.ceil(s1))
        F = lambda x: np.interp(x, bp, fv, left=0.0, right=cum[-1]) if len(bp) else np.zeros_like(x)
        gf = F(secs + 1) - F(secs)
        n = np.histogram(t, bins=np.append(secs, secs[-1] + 1))[0]
        lat = np.interp(secs, lat_t, lat_v) if lat_t is not None and np.isfinite(lat_v).any() else np.full(len(secs), np.nan)
        S.append(np.column_stack([secs, n, gf, lat]))
    S = np.vstack(S); rate, gf, lat = S[:, 1], S[:, 2], S[:, 3]
    print(f"=== {sat} {day}: {len(S)} 秒 ===")
    print("  按计数率分档：秒数 / 空洞占比中位 / 空洞占比>10% 的秒占比 / |lat| 中位")
    for lo, hi in ((0, 1000), (1000, 2000), (2000, 3000), (3000, 5000), (5000, 8000), (8000, 15000), (15000, 1e9)):
        mk = (rate >= lo) & (rate < hi)
        if mk.sum(): print(f"    {lo:5.0f}-{hi:5.0f} c/s: {mk.sum():6d} 秒  空洞中位 {100*np.median(gf[mk]):5.1f}%  >10% 占 {100*(gf[mk]>0.1).mean():5.1f}%  |lat| 中位 {np.nanmedian(np.abs(lat[mk])):5.1f}")
    print("  按 |lat| 分档：秒数 / 率中位 / 空洞占比中位")
    for lo, hi in ((0, 20), (20, 30), (30, 40), (40, 50), (50, 90)):
        mk = (np.abs(lat) >= lo) & (np.abs(lat) < hi)
        if mk.sum(): print(f"    |lat| {lo:2d}-{hi:2d}: {mk.sum():6d} 秒  率中位 {np.median(rate[mk]):6.0f} c/s  空洞中位 {100*np.median(gf[mk]):5.1f}%  >10% 占 {100*(gf[mk]>0.1).mean():5.1f}%")
    s = json.load(open(f"{R}/{sat}/{y}/{m}/{tag}_signals.json"))
    if s:
        cs = np.array([met(c["start"]) for c in s]); idx = np.clip(np.searchsorted(S[:, 0], np.floor(cs)), 0, len(S) - 1)
        cr, cg = rate[idx], gf[idx]
        print(f"  候选 {len(s)} 个所在秒：率中位 {np.median(cr):.0f} c/s（p10 {np.percentile(cr,10):.0f}），空洞占比中位 {100*np.median(cg):.0f}%，空洞==0 的候选 {int((cg==0).sum())} 个")
        quiet = [c for c, r_, g_ in zip(s, cr, cg) if r_ < 3000 and g_ == 0]
        print(f"  率<3000 且无空洞的候选: {len(quiet)} 个；其中 fa<=1e-5: {sum(c['false_positive_per_year']<=1e-5 for c in quiet)}")
        for c in sorted(quiet, key=lambda c: c['false_positive_per_year'])[:3]:
            print(f"     {c['start'][11:23]} count={c['count']} mean={c['mean']:.2f} fa={c['false_positive_per_year']:.1e} lat={c['position']['latitude']:.1f} lon={c['position']['longitude']:.1f}")
    sys.stdout.flush()
