"""计数守恒检验：真暴发在 ±100 ms 里净多出 excess 个计数；把积压事例重新盖时间戳的读出伪象总数不变（暴发的超出 = 周围的亏缺）。
用法: python3 grid_conserve.py <候选 CSV> <输出 CSV>"""
from astropy.io import fits
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
rows = list(csv.DictReader(open(sys.argv[1])))
w = csv.writer(open(sys.argv[2], "w", newline=""))
w.writerow(["sat", "start", "dur_search_us", "lat", "bkg_per_ms", "net_pm10ms", "net_pm100ms", "sigma_100", "net_pm500ms", "sigma_500", "gap_max_ms_pm100", "gap_rl", "gap_offset_ms"])
n = 0
for r in rows:
    sat = r["sat"]; t0 = met(r["start"])
    for f in pass_files(sat, t0):
        with fits.open(f) as h:
            g = h["GTI"].data; s0, s1 = float(g["START"][0]), float(g["STOP"][0])
            if not (s0 <= t0 <= s1): continue
            T = np.concatenate([np.asarray(h[f"EVENTS{k}"].data["TIME"], float)[np.asarray(h[f"EVENTS{k}"].data["EVT_TYPE"]) == 1] for k in range(4)]); T.sort()
            def cnt(a, b):
                a, b = max(a, s0), min(b, s1)
                return int(((T >= a) & (T <= b)).sum()), (b - a) * 1000.0
            # 本底：±1..±3 s（避开 ±1 s 内可能的结构），夹到 GTI
            nb, lb = cnt(t0 - 3, t0 - 1); nb2, lb2 = cnt(t0 + 1, t0 + 3); bpm = (nb + nb2) / max(lb + lb2, 1e-9)
            out = [sat, r["start"], r["dur_us"], r["lat"], f"{bpm:.3f}"]
            for half in (0.01, 0.1, 0.5):
                nn, ll = cnt(t0 - half, t0 + half); exp = bpm * ll
                out.append(f"{nn - exp:.1f}")
                if half > 0.01: out.append(f"{np.sqrt(exp):.1f}")
            m = (T >= t0 - 0.1) & (T <= t0 + 0.1); tt = T[m]
            if len(tt) > 1:
                d = np.diff(tt); i = int(np.argmax(d)); gap = d[i] * 1000
                out += [f"{gap:.2f}", f"{bpm*gap:.1f}", f"{(tt[i] - t0)*1000:.1f}"]
            else: out += ["", "", ""]
            w.writerow(out); n += 1; break
print("rows:", n)
