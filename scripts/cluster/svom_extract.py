"""提取指定 UTC 时刻附近的 GRM 事例，输出 CSV 供本地画光变曲线。"""
from astropy.io import fits
import numpy as np, glob, sys, os
from datetime import datetime, timezone

REF = datetime(2017, 1, 1, tzinfo=timezone.utc)   # MET 零点（= UTC，2017 后无闰秒）
HALF = 0.030                                       # 取 ±30 ms

targets = [
    ("A_phase_2025-11-12", "2025-11-12T03:28:41.507010"),
    ("B_phase_2026-01-20", "2026-01-20T17:44:15.506838"),
    ("C_free_2025-06-08",  "2025-06-08T10:43:55.344500"),
    ("D_free_2025-03-31",  "2025-03-31T12:09:17.077212"),
]
out = open("/scratchfs2/gecam/guohx/svomrun/lightcurves.csv", "w")
out.write("tag,det,dt_s,pi\n")
for tag, iso in targets:
    t = datetime.strptime(iso, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)
    met = (t - REF).total_seconds()
    d = "/gecamfs/SVOM/Archived-DATA/GRM-DATA/L1B/daily/%s/grm_evt" % t.strftime("%Y/%m/%d")
    pat = "%s/svom_grm_evt_%s_%02d_v*.fits" % (d, t.strftime("%y%m%d"), t.hour)
    g = sorted(glob.glob(pat))
    if not g:
        print("  %s: 没有文件 %s" % (tag, pat)); continue
    n = 0
    with fits.open(g[-1]) as h:
        for i, det in zip((3, 4, 5), (1, 2, 3)):
            tt = np.asarray(h[i].data["TIME"], dtype=np.float64)
            pi = np.asarray(h[i].data["PI"])
            et = np.asarray(h[i].data["EVT_TYPE"])
            m = (tt >= met - HALF) & (tt <= met + HALF) & (et == 0) & (pi >= 15) & (pi < 256)
            for x, p in zip(tt[m] - met, pi[m]):
                out.write("%s,%d,%.9f,%d\n" % (tag, det, x, p))
            n += m.sum()
    print("  %s  MET=%.6f  ±30ms 内 keep 事例 %d" % (tag, met, n))
out.close()
print("written /scratchfs2/gecam/guohx/svomrun/lightcurves.csv")
