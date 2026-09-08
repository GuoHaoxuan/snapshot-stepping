"""导出候选前后 ±50 ms 的逐事例（时间/能道/探头），供本地画光变。
用法: python3 grid_lightcurve.py <候选 CSV: sat,start,...> <输出目录>"""
from astropy.io import fits
import glob, os, csv, sys, numpy as np, datetime as dt

G = "/gecamfs/Exchange/GSDC/missions/GRID"
REF = dt.datetime(2018, 1, 1, tzinfo=dt.timezone.utc)
HALF = 0.05  # 导出半宽 50 ms
BKG = 1.0    # 本底窗半宽 1 s


def met(iso):
    b = iso.rstrip("Z")
    h, f = b.split(".")
    b = h + "." + (f + "000000")[:6]
    return (dt.datetime.strptime(b, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=dt.timezone.utc) - REF).total_seconds()


def pass_files(sat, t0):
    out = []
    for back in (0, 1):
        day = (REF + dt.timedelta(seconds=t0) - dt.timedelta(days=back)).strftime("%Y/%m/%d")
        dd = f"{G}/{sat}/fits7/{day}"
        if not os.path.isdir(dd):
            continue
        v = sorted(os.listdir(dd))[-1]
        out += sorted(glob.glob(f"{dd}/{v}/*.fits"))
    return out


rows = list(csv.DictReader(open(sys.argv[1])))
outdir = sys.argv[2]
os.makedirs(outdir, exist_ok=True)

for r in rows:
    sat, t0 = r["sat"], met(r["start"])
    tag = f"{sat}_{r['start'][:23].replace(':', '').replace('-', '').replace('.', '')}"
    done = False
    for f in pass_files(sat, t0):
        with fits.open(f) as h:
            g = h["GTI"].data
            s0, s1 = float(g["START"][0]), float(g["STOP"][0])
            if not (s0 <= t0 <= s1):
                continue
            T, P, D = [], [], []
            for k in range(4):
                x = h[f"EVENTS{k}"].data
                t = np.asarray(x["TIME"], float)
                m = np.asarray(x["EVT_TYPE"]) == 1
                keep = m & (t >= t0 - BKG) & (t <= t0 + BKG)
                T.append(t[keep])
                P.append(np.asarray(x["PI"])[keep])
                D.append(np.full(int(keep.sum()), k))
            T = np.concatenate(T); P = np.concatenate(P); D = np.concatenate(D)
            o = np.argsort(T); T, P, D = T[o], P[o], D[o]
            # 本底谱/本底率用整个 ±1 s（GTI 内），光变只导出 ±50 ms
            live = min(t0 + BKG, s1) - max(t0 - BKG, s0)
            near = np.abs(T - t0) <= HALF
            with open(f"{outdir}/{tag}.csv", "w", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["dt_ms", "pi", "det"])
                for tt, pp, dd_ in zip(T[near], P[near], D[near]):
                    wr.writerow([f"{(tt - t0) * 1e3:.6f}", int(pp), int(dd_)])
            far = np.abs(T - t0) > 0.005
            with open(f"{outdir}/{tag}_bkg.csv", "w", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["live_s", "n_bkg", "pi"])
                wr.writerow([f"{live:.3f}", int(far.sum()), ""])
                for pp in P[far]:
                    wr.writerow(["", "", int(pp)])
            print(tag, "events", int(near.sum()), "bkg_rate", f"{far.sum() / live:.1f}")
            done = True
            break
    if not done:
        print(tag, "NOT FOUND")
