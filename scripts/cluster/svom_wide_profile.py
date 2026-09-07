"""SVOM 证实 TGF 的真实时长：在 ±30 ms 宽窗里量，判断 1 ms 的搜索上限有没有截断它们。

搜索只测 ≤1 ms 的窗，一个真实时长 3 ms 的暴发会被最亮的 1 ms 子窗抓到、时长被记短。
这里绕开搜索，直接在事例流上量：100 µs 分格，找连续显著格构成的区间，再算背景扣除后的 T50/T90。

用法: python3 svom_wide_profile.py <sample.csv> <out.csv>
"""
import csv, glob, os, sys
from datetime import datetime, timedelta, timezone

import numpy as np
from astropy.io import fits
from scipy.stats import poisson

D = "/gecamfs/SVOM/Archived-DATA/GRM-DATA/L1B/daily"
REF = datetime(2017, 1, 1, tzinfo=timezone.utc)
HALF = 0.030          # ±30 ms 宽窗
BKG_IN, BKG_OUT = 0.005, 0.030   # 本底取 5–30 ms 环
BIN = 100e-6


def met(iso):
    b = iso.rstrip("Z"); h, _, f = b.partition(".")
    return (datetime.strptime(h + "." + (f + "000000")[:6], "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc) - REF).total_seconds()


def load_hour(m0):
    utc = REF + timedelta(seconds=m0)
    pat = "%s/%s/grm_evt/svom_grm_evt_%s_%s_v*.fits" % (D, utc.strftime("%Y/%m/%d"), utc.strftime("%y%m%d"), utc.strftime("%H"))
    g = sorted(glob.glob(pat))
    if not g:
        return None
    with fits.open(g[-1]) as h:
        T = []
        for i in (3, 4, 5):
            d = h[i].data
            t = np.asarray(d["TIME"], float); pi = np.asarray(d["PI"]); et = np.asarray(d["EVT_TYPE"])
            k = (et == 0) & (pi >= 25) & (pi < 256)      # 与搜索同款事例准入（能阈 ch25）
            T.append(t[k])
    T = np.concatenate(T); T.sort()
    return T


def main(sample, out):
    rows = list(csv.DictReader(open(sample)))
    w = csv.writer(open(out, "w", newline=""))
    w.writerow(["start", "span_ms", "best_ms", "n_core", "bkg_per_ms", "n_sig_bins", "run_ms",
                "excess", "t50_ms", "t90_ms", "excess_in_1ms", "frac_in_best_1ms"])
    for r in rows:
        m0 = met(r["start"])
        T = load_hour(m0)
        if T is None:
            print("缺文件", r["start"]); continue
        rel = T - m0
        near = np.abs(rel) <= BKG_OUT
        t = rel[near]
        ring = (np.abs(t) >= BKG_IN) & (np.abs(t) <= BKG_OUT)
        live_ring = 2 * (BKG_OUT - BKG_IN)
        bpm = ring.sum() / (live_ring * 1e3)                       # 每 ms 的本底计数
        edges = np.arange(-0.005, 0.005 + BIN / 2, BIN)
        cnt, _ = np.histogram(t, bins=edges)
        mu = bpm * (BIN * 1e3)
        sig = poisson.sf(cnt - 1, mu) < 1e-3
        if not sig.any():
            w.writerow([r["start"][:23], r["span"], r["best"], r["n"], "%.3f" % bpm, 0, 0, 0, "", "", "", ""])
            continue
        peak = int(np.argmax(cnt))
        lo = peak
        while lo - 1 >= 0 and sig[lo - 1]:
            lo -= 1
        hi = peak
        while hi + 1 < len(cnt) and sig[hi + 1]:
            hi += 1
        run_ms = (hi - lo + 1) * BIN * 1e3
        seg = (t >= edges[lo]) & (t < edges[hi + 1])
        excess = seg.sum() - bpm * run_ms
        # T50/T90：区间内背景扣除后的累积计数
        ts = np.sort(t[seg])
        if len(ts) > 2 and excess > 3:
            frac_bkg = bpm * run_ms / max(len(ts), 1)
            cum = np.arange(1, len(ts) + 1) - frac_bkg * len(ts) * (ts - ts[0]) / max(ts[-1] - ts[0], 1e-9)
            cum = np.maximum.accumulate(cum); cum /= cum[-1]
            t05, t25, t75, t95 = [float(np.interp(q, cum, ts)) for q in (0.05, 0.25, 0.75, 0.95)]
            t50 = (t75 - t25) * 1e3; t90 = (t95 - t05) * 1e3
        else:
            t50 = t90 = float("nan")
        # 最亮的 1 ms 子窗能收多少：上限截断损失的直接度量
        best1 = 0
        for s in np.arange(edges[lo] - 0.001, edges[hi + 1] + 1e-9, BIN):
            c1 = ((t >= s) & (t < s + 0.001)).sum()
            best1 = max(best1, c1 - bpm * 1.0)
        w.writerow([r["start"][:23], r["span"], r["best"], r["n"], "%.3f" % bpm, int(sig.sum()), "%.1f" % run_ms,
                    "%.1f" % excess, "%.3f" % t50, "%.3f" % t90, "%.1f" % best1,
                    "%.3f" % (best1 / excess if excess > 0 else float("nan"))])
        print("%s 跨度%.2f 最佳窗%.2f | 本底%.2f/ms 显著格%d 连续区间%.1f ms 超出%.0f T50=%.2f T90=%.2f 最亮1ms收%.0f(%.0f%%)"
              % (r["start"][:23], float(r["span"]), float(r["best"]), bpm, sig.sum(), run_ms, excess, t50, t90, best1,
                 100 * best1 / excess if excess > 0 else float("nan")), flush=True)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
