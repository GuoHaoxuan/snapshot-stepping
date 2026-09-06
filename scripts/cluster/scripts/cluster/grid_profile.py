"""显著候选的宽窗时间轮廓：搜索把时长截在 1 ms 上限，这里在 ±30 ms 内用 1 ms 格、±3 ms 内用 100 µs 格量真实时长。
用法: python3 grid_profile.py <候选 CSV: sat,start,...> <输出 CSV>"""
from astropy.io import fits
from scipy.stats import poisson
import glob, os, csv, sys, numpy as np, datetime as dt
G = "/gecamfs/Exchange/GSDC/missions/GRID"
REF = dt.datetime(2018, 1, 1, tzinfo=dt.timezone.utc)
def met(iso):
    b = iso.rstrip("Z"); h, f = b.split("."); b = h + "." + (f + "000000")[:6]
    return (dt.datetime.strptime(b, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=dt.timezone.utc) - REF).total_seconds()
def pass_files(sat, t0):
    out = []
    for back in (0, 1):
        day = (REF + dt.timedelta(seconds=t0) - dt.timedelta(days=back)).strftime("%Y/%m/%d")
        dd = f"{G}/{sat}/fits7/{day}"
        if not os.path.isdir(dd): continue
        v = sorted(os.listdir(dd))[-1]; out += sorted(glob.glob(f"{dd}/{v}/*.fits"))
    return out
def run_len(sig, peak):
    # 从峰所在格向两边延伸的连续显著格数
    lo = peak
    while lo - 1 >= 0 and sig[lo - 1]: lo -= 1
    hi = peak
    while hi + 1 < len(sig) and sig[hi + 1]: hi += 1
    return lo, hi
rows = list(csv.DictReader(open(sys.argv[1])))
w = csv.writer(open(sys.argv[2], "w", newline=""))
w.writerow(["sat", "start", "dur_search_us", "lat", "bkg_per_ms", "dur_1ms_bins", "excess_1ms", "n_sig_bins_pm30ms", "dur_100us_bins", "excess_100us", "det_frac_max_wide", "n_det_wide", "pi_med_wide", "pi_med_bkg"])
n = 0
for r in rows:
    sat = r["sat"]; t0 = met(r["start"])
    for f in pass_files(sat, t0):
        with fits.open(f) as h:
            g = h["GTI"].data; s0, s1 = float(g["START"][0]), float(g["STOP"][0])
            if not (s0 <= t0 <= s1): continue
            T = []; P = []; D = []
            for k in range(4):
                x = h[f"EVENTS{k}"].data; t = np.asarray(x["TIME"], float); et = np.asarray(x["EVT_TYPE"]); pi = np.asarray(x["PI"])
                m = et == 1; T.append(t[m]); P.append(pi[m]); D.append(np.full(m.sum(), k))
            T = np.concatenate(T); P = np.concatenate(P); D = np.concatenate(D); o = np.argsort(T); T, P, D = T[o], P[o], D[o]
            bk = ((T >= t0 - 1) & (T <= t0 - 0.05)) | ((T >= t0 + 0.05) & (T <= t0 + 1))
            # 本底窗也夹到 GTI 内
            live = (min(t0 - 0.05, s1) - max(t0 - 1, s0)) + (min(t0 + 1, s1) - max(t0 + 0.05, s0))
            bpm = bk.sum() / (live * 1000.0)
            def profile(width_ms, half_ms):
                edges = t0 + np.arange(-half_ms, half_ms + width_ms / 2, width_ms) * 1e-3
                cnt, _ = np.histogram(T, bins=edges)
                mu = bpm * width_ms
                sig = poisson.sf(cnt - 1, mu) < 1e-3
                peak = int(np.argmax(cnt))
                if not sig[peak]: return 0, 0.0, int(sig.sum()), None
                lo, hi = run_len(sig, peak)
                excess = cnt[lo:hi + 1].sum() - mu * (hi - lo + 1)
                return hi - lo + 1, excess, int(sig.sum()), (edges[lo], edges[hi + 1])
            d1, e1, ns1, span1 = profile(1.0, 30.0)
            d01, e01, ns01, span01 = profile(0.1, 3.0)
            span = span1 or span01
            if span is None:
                w.writerow([sat, r["start"], r["dur_us"], r["lat"], f"{bpm:.3f}", 0, 0, ns1, 0, 0, "", "", "", ""]); n += 1; break
            wide = (T >= span[0]) & (T < span[1])
            dets, cnts = np.unique(D[wide], return_counts=True)
            w.writerow([sat, r["start"], r["dur_us"], r["lat"], f"{bpm:.3f}", d1, f"{e1:.1f}", ns1, d01, f"{e01:.1f}",
                        f"{cnts.max()/wide.sum():.2f}", len(dets), int(np.median(P[wide])), int(np.median(P[bk])) if bk.any() else ""])
            n += 1; break
print("rows:", n)
