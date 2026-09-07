"""SVOM/GRM 能阈核对：证实 TGF、未证实显著候选、本底的 PI 谱，以及能阈扫描（不同能阈下证实样本的信噪）。
输入: svom_top795.csv（显著候选 start,stop,...）、assoc_svom_v5.csv（associated,in_coverage）
输出: spectra.csv（逐道计数）、thr_scan.csv（逐候选逐能阈的 S、B）、burst.csv（低道事例的成簇性）、ebounds.csv"""
from astropy.io import fits
import numpy as np, glob, csv, sys, os
from datetime import datetime, timezone
from scipy.stats import poisson
D = "/gecamfs/SVOM/Archived-DATA/GRM-DATA/L1B/daily"; REF = datetime(2017, 1, 1, tzinfo=timezone.utc)
R = "/scratchfs2/gecam/guohx/svomrun5"
def met(iso):
    b = iso.rstrip("Z"); h, _, f = b.partition("."); return (datetime.strptime(h + "." + (f + "000000")[:6], "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc) - REF).total_seconds()
cands = list(csv.DictReader(open(f"{R}/svom_top795.csv")))
assoc = {r["start"][:26]: (r["associated"] == "1", r["in_coverage"] == "1") for r in csv.DictReader(open(f"{R}/assoc_svom_v5.csv"))}
by_hour = {}
for c in cands:
    m0, m1 = met(c["start"]), met(c["stop"]); utc = REF + __import__("datetime").timedelta(seconds=m0)
    by_hour.setdefault((utc.strftime("%Y/%m/%d"), utc.strftime("%H")), []).append((c, m0, m1))
THR = [10, 12, 15, 18, 20, 25, 30, 40, 50]
spec = {"assoc": np.zeros(260), "nonassoc": np.zeros(260), "bkg": np.zeros(260), "bkg_seconds": 0.0, "hour_all": np.zeros(260), "hour_seconds": 0.0}
scan = open(f"{R}/thr_scan.csv", "w"); scan.write("start,fa,assoc,cov,dur_s," + ",".join(f"S{t}" for t in THR) + "," + ",".join(f"B{t}" for t in THR) + "\n")
burst = open(f"{R}/burst.csv", "w"); burst.write("hour,band,n_events,ms_bins_ge5,frac_events_in_bursty_bins\n")
done_hours = 0
for (day, hh), items in sorted(by_hour.items()):
    g = sorted(glob.glob(f"{D}/{day}/grm_evt/svom_grm_evt_{datetime.strptime(day, '%Y/%m/%d').strftime('%y%m%d')}_{hh}_v*.fits"))
    if not g: continue
    with fits.open(g[-1]) as h:
        if done_hours == 0:
            with open(f"{R}/ebounds.csv", "w") as f:
                f.write("channel,e_min,e_max\n")
                for r in h["EBOUNDS"].data: f.write(f"{int(r[0])},{float(r[1]):.4f},{float(r[2]):.4f}\n")
        T = []; P = []
        for i in (3, 4, 5):
            d = h[i].data; t = np.asarray(d["TIME"], float); pi = np.asarray(d["PI"]); et = np.asarray(d["EVT_TYPE"])
            k = (et == 0) & (pi < 256); T.append(t[k]); P.append(pi[k])
        T = np.concatenate(T); P = np.concatenate(P); o = np.argsort(T); T, P = T[o], P[o]
        gti = [(float(a), float(b)) for a, b in zip(h["GTI"].data["START"], h["GTI"].data["STOP"])]
    ingti = np.zeros(len(T), bool)
    for a, b in gti: ingti |= (T >= a) & (T <= b)
    T, P = T[ingti], P[ingti]
    # 整小时谱（前 20 个小时够了）
    if done_hours < 20:
        spec["hour_all"] += np.bincount(P, minlength=260)[:260]; spec["hour_seconds"] += sum(b - a for a, b in gti)
        # 低道事例成簇性：1 ms 格里 ≥5 个的格子占比
        for band, lo, hi in (("ch10-14", 10, 15), ("ch15-19", 15, 20), ("ch20-29", 20, 30), ("ch30-49", 30, 50), ("ch50-255", 50, 256)):
            m = (P >= lo) & (P < hi); tt = T[m]
            if len(tt) < 100: continue
            bins = np.floor((tt - tt[0]) * 1000).astype(np.int64); cnt = np.bincount(bins)
            big = cnt >= 5; n_in = int(cnt[big].sum())
            burst.write(f"{day} {hh},{band},{len(tt)},{int(big.sum())},{n_in/len(tt):.5f}\n")
    for c, m0, m1 in items:
        a_, cov = assoc.get(c["start"][:26], (False, True))
        core = (T >= m0) & (T <= m1); near = ((T >= m0 - 1) & (T < m0 - 0.01)) | ((T > m1 + 0.01) & (T <= m1 + 1))
        key = "assoc" if a_ else "nonassoc"
        spec[key] += np.bincount(P[core], minlength=260)[:260]; spec["bkg"] += np.bincount(P[near], minlength=260)[:260]; spec["bkg_seconds"] += 1.98
        dur = m1 - m0
        S = [int((core & (P >= t)).sum()) for t in THR]; B = [float((near & (P >= t)).sum()) / 1.98 * dur for t in THR]
        scan.write(f"{c['start'][:26]},{c['false_positive_per_year']},{int(a_)},{int(cov)},{dur:.6f}," + ",".join(map(str, S)) + "," + ",".join(f"{b:.4f}" for b in B) + "\n")
    done_hours += 1
    if done_hours % 50 == 0: print(done_hours, "hours", flush=True)
with open(f"{R}/spectra.csv", "w") as f:
    f.write("channel,assoc_core,nonassoc_core,bkg,hour_all\n")
    for ch in range(260): f.write(f"{ch},{int(spec['assoc'][ch])},{int(spec['nonassoc'][ch])},{int(spec['bkg'][ch])},{int(spec['hour_all'][ch])}\n")
    f.write(f"# bkg_seconds={spec['bkg_seconds']:.1f} hour_seconds={spec['hour_seconds']:.1f}\n")
print("done", done_hours, "hours")
