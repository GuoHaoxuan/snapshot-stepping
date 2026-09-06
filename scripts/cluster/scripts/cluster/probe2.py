from astropy.io import fits
import numpy as np, glob
D = "/hxmtfs/data/Fermi_GBM/2019/01/01/current"
DETS = ["n0","n1","n2","n3","n4","n5","n6","n7","n8","n9","na","nb","b0","b1"]
# 候选 2019-01-01T00:00:01.306729912Z, count=10, dur=6.08us
# 直接用 MET 找：UTC 秒 = 1.306729912, MET = (2019-01-01 - 2001-01-01) + 5 + 1.306729912
from datetime import datetime, timezone
REF = datetime(2001,1,1,tzinfo=timezone.utc)
base = (datetime(2019,1,1,tzinfo=timezone.utc)-REF).total_seconds() + 5
m0 = base + 1.306729912
print("base MET = %.6f, m0 = %.9f" % (base, m0))
rows=[]
for det in DETS:
    g = sorted(glob.glob("%s/glg_tte_%s_190101_00z_v*.fit.gz" % (D, det)))
    if not g: continue
    with fits.open(g[0]) as h:
        tt = np.asarray(h["EVENTS"].data["TIME"], dtype=np.float64)
        pha = np.asarray(h["EVENTS"].data["PHA"])
    k=(pha>0)&(pha<127); tt=tt[k]; pha=pha[k]
    lo=np.searchsorted(tt,m0-30e-6); hi=np.searchsorted(tt,m0+30e-6)
    for t,p in zip(tt[lo:hi],pha[lo:hi]):
        rows.append((t,det,int(p)))
rows.sort()
print("±30us 内共 %d 个事例：" % len(rows))
for t,det,p in rows:
    print("   %.9f  (%+9.3f us)  %s pha=%d" % (t, (t-m0)*1e6, det, p))
