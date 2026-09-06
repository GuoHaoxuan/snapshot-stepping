from astropy.io import fits
import numpy as np, glob, json
from datetime import datetime, timezone

REF = datetime(2001,1,1,tzinfo=timezone.utc); LEAPS = 5
D = "/hxmtfs/data/Fermi_GBM/2019/01/01/current"
DETS = ["n0","n1","n2","n3","n4","n5","n6","n7","n8","n9","na","nb","b0","b1"]

def met(iso):
    body = iso.rstrip("Z")
    head, frac = body.split(".")
    body = head + "." + (frac + "000000")[:6]
    t = datetime.strptime(body, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)
    return (t - REF).total_seconds() + LEAPS, t

s = json.load(open("/scratchfs2/gecam/guohx/gbmrun/data/Fermi_GBM/2019/01/20190101_signals.json"))
s.sort(key=lambda c: c["false_positive_per_year"])
probe = [s[0], s[1]] + [c for c in s if c["start"].startswith("2019-01-01T00:00:01")][:1]

for c in probe:
    m0, t = met(c["start"]); m1, _ = met(c["stop"])
    print("=== %s  count=%d mean=%.3f fa=%.2e dur=%.1fus" % (
        c["start"], c["count"], c["mean"], c["false_positive_per_year"], (m1-m0)*1e6))
    hh = t.hour
    tot_in, tot_wide = 0, 0
    for det in DETS:
        g = sorted(glob.glob("%s/glg_tte_%s_190101_%02dz_v*.fit.gz" % (D, det, hh)))
        if not g:
            print("   %s: NO FILE" % det); continue
        with fits.open(g[0]) as h:
            tt = np.asarray(h["EVENTS"].data["TIME"], dtype=np.float64)
            pha = np.asarray(h["EVENTS"].data["PHA"])
        k = (pha > 0) & (pha < 127)
        tt_k = tt[k]
        n_in = np.searchsorted(tt_k, m1, "right") - np.searchsorted(tt_k, m0, "left")
        # 宽窗 ±1ms，看看事例到底落在哪
        lo = np.searchsorted(tt_k, m0-1e-3, "left"); hi = np.searchsorted(tt_k, m1+1e-3, "right")
        tot_in += n_in; tot_wide += hi-lo
        if n_in or hi-lo:
            off = (tt_k[lo:hi] - m0)*1e6
            print("   %s: in=%2d wide=%3d  offsets_us=%s" % (
                det, n_in, hi-lo, np.array2string(off[:12], precision=1, max_line_width=200)))
    print("   >>> 窗内合计 %d (搜索说 %d)，±1ms 内合计 %d" % (tot_in, c["count"], tot_wide))
