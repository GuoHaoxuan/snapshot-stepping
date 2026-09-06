"""GRID 显著候选（fa<=1e-5, 本底率<=5 kc/s）的事例级特征：几路探测器、时间戳簇集、能谱、本底。"""
from astropy.io import fits
import glob, os, json, csv, sys, numpy as np, datetime as dt
G = "/gecamfs/Exchange/GSDC/missions/GRID"; R = "/scratchfs2/gecam/guohx/gridrun/data"
REF = dt.datetime(2018, 1, 1, tzinfo=dt.timezone.utc)
def met(iso):
    b = iso.rstrip("Z"); h, f = b.split("."); b = h + "." + (f + "000000")[:6]
    return (dt.datetime.strptime(b, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=dt.timezone.utc) - REF).total_seconds()
w = csv.writer(open("/scratchfs2/gecam/guohx/gridrun/features_sig.csv", "w", newline=""))
w.writerow(["sat", "start", "fa", "count", "mean", "dur_us", "lon", "lat", "rate_win", "n_core", "n_det_hit", "det_frac_max", "max_mult", "n_times", "min_dt_us", "pi_core_med", "pi_bkg_med", "pi_ratio", "n_pm5ms", "pass_rate"])
n = 0
for sat in ("GRID-02", "GRID-03B", "GRID-04", "GRID-07"):
    sig = []
    for f in sorted(glob.glob(f"{R}/{sat}/*/*/*_signals.json")): sig += json.load(open(f))
    sig = [c for c in sig if c["false_positive_per_year"] <= 1e-5 and c["mean"] / c["bin_size_best"] <= 5000]
    cache = {}
    for c in sig:
        t0 = met(c["start"]) + c["delay"]; t1 = t0 + c["bin_size_best"]; day = c["start"][:10].replace("-", "/")
        key = day
        if key not in cache:
            dd = f"{G}/{sat}/fits7/{day}"; v = sorted(os.listdir(dd))[-1]; cache[key] = sorted(glob.glob(f"{dd}/{v}/*.fits"))
        for f in cache[key]:
            with fits.open(f) as h:
                g = h["GTI"].data; s0, s1 = float(g["START"][0]), float(g["STOP"][0])
                if not (s0 <= t0 <= s1): continue
                T = []; P = []; D = []
                for k in range(4):
                    x = h[f"EVENTS{k}"].data; t = np.asarray(x["TIME"], float); et = np.asarray(x["EVT_TYPE"]); pi = np.asarray(x["PI"])
                    m = et == 1; T.append(t[m]); P.append(pi[m]); D.append(np.full(m.sum(), k))
                T = np.concatenate(T); P = np.concatenate(P); D = np.concatenate(D); o = np.argsort(T); T, P, D = T[o], P[o], D[o]
                core = (T >= t0) & (T <= t1); near = (T >= t0 - 1) & (T <= t1 + 1) & ~((T >= t0 - 0.01) & (T <= t1 + 0.01))
                nc = int(core.sum())
                if nc == 0: break
                dets, cnts = np.unique(D[core], return_counts=True); _, mult = np.unique(T[core], return_counts=True)
                dts = np.diff(T[core]); mindt = dts[dts > 0].min() * 1e6 if (dts > 0).any() else 0.0
                win = (T >= t0 - 0.5) & (T <= t1 + 0.5)
                pm5 = int(((T >= t0 - 0.005) & (T <= t1 + 0.005)).sum()) - nc
                pb = np.median(P[near]) if near.any() else np.nan
                w.writerow([sat, c["start"][:23], f"{c['false_positive_per_year']:.3e}", c["count"], f"{c['mean']:.4f}", f"{c['bin_size_best']*1e6:.1f}",
                            f"{c['position']['longitude']:.2f}", f"{c['position']['latitude']:.2f}", f"{win.sum()/(1.0+c['bin_size_best']):.0f}", nc, len(dets), f"{cnts.max()/nc:.3f}",
                            int(mult.max()), len(mult), f"{mindt:.2f}", int(np.median(P[core])), pb if np.isnan(pb) else int(pb), f"{np.median(P[core])/max(pb,1):.3f}" if not np.isnan(pb) else "", pm5, f"{((T>=s0)&(T<=s1)).sum()/(s1-s0):.0f}"])
                n += 1
                break
print("rows:", n)
