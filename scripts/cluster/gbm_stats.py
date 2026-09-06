"""统计 GBM 候选窗内的逐探头计数分布。

一天的候选全在同一天的数据里，所以按小时载入一次事例、处理该小时的全部
候选。输出每个候选在 14 个探头上的计数，据此可以同时评估几种分组方案：
  3 组  = 12 个 NaI 合一组 + b0 + b1（GBM 团队的做法）
  14 组 = 每个探头一组
"""
from astropy.io import fits
import numpy as np, glob, json, csv
from collections import defaultdict
from datetime import datetime, timezone

REF = datetime(2001, 1, 1, tzinfo=timezone.utc)
LEAPS = 5
D = "/hxmtfs/data/Fermi_GBM/2019/01/01/current"
DETS = ["n0","n1","n2","n3","n4","n5","n6","n7","n8","n9","na","nb","b0","b1"]

def met(iso):
    # serde 会截掉小数末尾的零，所以秒的小数位数不固定，不能按定长切片。
    body = iso.rstrip("Z")
    if "." in body:
        head, frac = body.split(".")
        body = head + "." + (frac + "000000")[:6]
    else:
        body += ".000000"
    t = datetime.strptime(body, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)
    return (t - REF).total_seconds() + LEAPS, t

s = json.load(open("/scratchfs2/gecam/guohx/gbmrun/data/Fermi_GBM/2019/01/20190101_signals.json"))
by_hour = defaultdict(list)
for c in s:
    m0, t = met(c["start"]); m1, _ = met(c["stop"])
    by_hour[t.hour].append((c, m0, m1))
print("候选 %d，分布在 %d 小时" % (len(s), len(by_hour)), flush=True)

out = open("/scratchfs2/gecam/guohx/gbmrun/cand_detstats.csv", "w", newline="")
w = csv.writer(out)
w.writerow(["start", "fa", "count", "mean", "dur_us", "lon", "lat", "n_win",
            "pha_med"] + DETS)
for hh in sorted(by_hour):
    per_det = {}
    for det in DETS:
        g = sorted(glob.glob("%s/glg_tte_%s_190101_%02dz_v*.fit.gz" % (D, det, hh)))
        if not g:
            per_det[det] = (np.empty(0), np.empty(0, dtype=np.int16)); continue
        with fits.open(g[0]) as h:
            tt = np.asarray(h["EVENTS"].data["TIME"], dtype=np.float64)
            pha = np.asarray(h["EVENTS"].data["PHA"])
        k = (pha > 0) & (pha < 127)
        per_det[det] = (tt[k], pha[k])
    for c, m0, m1 in by_hour[hh]:
        counts = []; phas = []
        for det in DETS:
            tt, pha = per_det[det]
            i0, i1 = np.searchsorted(tt, [m0, m1], side="left"), None
            lo = np.searchsorted(tt, m0, side="left")
            hi = np.searchsorted(tt, m1, side="right")
            counts.append(hi - lo)
            if hi > lo: phas.append(pha[lo:hi])
        n_win = int(sum(counts))
        pm = int(np.median(np.concatenate(phas))) if phas else -1
        w.writerow([c["start"][11:26], "%.3e" % c["false_positive_per_year"],
                    c["count"], "%.3f" % c["mean"], "%.1f" % ((m1 - m0) * 1e6),
                    "%.2f" % c["position"]["longitude"], "%.2f" % c["position"]["latitude"],
                    n_win, pm] + counts)
    print("  T%02d 完成 %d 个候选" % (hh, len(by_hour[hh])), flush=True)
out.close()
print("done -> cand_detstats.csv")
