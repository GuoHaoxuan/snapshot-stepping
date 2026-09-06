"""统计 GBM 候选窗内的逐探头计数分布（修正版）。

对前一版的两处修正：
  1. 时间用全纳秒精度。候选窗只有几微秒宽，把 ISO 截到微秒会把边界上的
     整簇事例挡在窗外——实测有候选因此从 10 个数成 1 个。
  2. 窗口取最显著的那一格 [start+delay, start+delay+bin_size_best]，而不是
     合并后的包络 [start, stop]。包络比最佳格宽时会多数进邻近事例。
比较用半个 ulp 的容差：MET 量级 5.7e8，f64 的间隔是 119 ns，事例时间本身
就落在这个格子上，而 start 由 MET→UTC→ISO 转回来最多差不到 1 ns。
"""
from astropy.io import fits
import numpy as np, glob, json, csv, re
from collections import defaultdict
from datetime import datetime, timezone

REF = datetime(2001, 1, 1, tzinfo=timezone.utc)
LEAPS = 5
DAY = "2019-01-01"
D = "/hxmtfs/data/Fermi_GBM/2019/01/01/current"
STAMP = "190101"
DETS = ["n0","n1","n2","n3","n4","n5","n6","n7","n8","n9","na","nb","b0","b1"]
TOL = 6e-8  # 半个 ulp

BASE = (datetime(2019, 1, 1, tzinfo=timezone.utc) - REF).total_seconds() + LEAPS

ISO = re.compile(r"T(\d\d):(\d\d):(\d\d)\.(\d+)Z?$")

def met_of(iso):
    """ISO → MET，保住纳秒。"""
    hh, mm, ss, frac = ISO.search(iso).groups()
    ns = int((frac + "000000000")[:9])
    return BASE + int(hh) * 3600 + int(mm) * 60 + int(ss) + ns * 1e-9, int(hh)

s = json.load(open("/scratchfs2/gecam/guohx/gbmrun/data/Fermi_GBM/2019/01/20190101_signals.json"))
by_hour = defaultdict(list)
for c in s:
    m, hh = met_of(c["start"])
    m0 = m + c["delay"]
    m1 = m0 + c["bin_size_best"]
    by_hour[hh].append((c, m0, m1))
print("候选 %d，分布在 %d 小时" % (len(s), len(by_hour)), flush=True)

out = open("/scratchfs2/gecam/guohx/gbmrun/cand_detstats2.csv", "w", newline="")
w = csv.writer(out)
w.writerow(["start", "fa", "count", "mean", "dur_us", "lon", "lat",
            "n_win", "n_nai", "n_b0", "n_b1", "n_det", "n_times", "max_mult",
            "pha_med"] + DETS)
for hh in sorted(by_hour):
    per_det = {}
    for det in DETS:
        g = sorted(glob.glob("%s/glg_tte_%s_%s_%02dz_v*.fit.gz" % (D, det, STAMP, hh)))
        if not g:
            per_det[det] = (np.empty(0), np.empty(0, dtype=np.int16)); continue
        with fits.open(g[-1]) as h:   # 版本号最大的，与 Rust 侧一致
            tt = np.asarray(h["EVENTS"].data["TIME"], dtype=np.float64)
            pha = np.asarray(h["EVENTS"].data["PHA"])
        k = (pha > 0) & (pha < 127)
        per_det[det] = (tt[k], pha[k])
    for c, m0, m1 in by_hour[hh]:
        counts, phas, times = [], [], []
        for det in DETS:
            tt, pha = per_det[det]
            lo = np.searchsorted(tt, m0 - TOL, side="left")
            hi = np.searchsorted(tt, m1 + TOL, side="right")
            counts.append(hi - lo)
            if hi > lo:
                phas.append(pha[lo:hi]); times.append(tt[lo:hi])
        n_nai = int(sum(counts[:12])); n_b0, n_b1 = counts[12], counts[13]
        n_win = n_nai + n_b0 + n_b1
        n_det = int(sum(1 for x in counts if x > 0))
        if times:
            allt = np.concatenate(times)
            uniq, mult = np.unique(allt, return_counts=True)
            n_times, max_mult = len(uniq), int(mult.max())
        else:
            n_times, max_mult = 0, 0
        pm = int(np.median(np.concatenate(phas))) if phas else -1
        w.writerow([c["start"][11:30], "%.3e" % c["false_positive_per_year"],
                    c["count"], "%.3f" % c["mean"], "%.2f" % ((m1 - m0) * 1e6),
                    "%.2f" % c["position"]["longitude"], "%.2f" % c["position"]["latitude"],
                    n_win, n_nai, n_b0, n_b1, n_det, n_times, max_mult, pm] + counts)
    print("  T%02d 完成 %d 个候选" % (hh, len(by_hour[hh])), flush=True)
out.close()
print("done -> cand_detstats2.csv")
