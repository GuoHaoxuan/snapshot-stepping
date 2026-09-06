import numpy as np, csv
rows = list(csv.DictReader(open("/scratchfs2/gecam/guohx/gbmrun/cand_detstats2.csv")))
g = lambda k, t=int: np.array([t(r[k]) for r in rows])
cnt, nw = g("count"), g("n_win")
nai, b0, b1 = g("n_nai"), g("n_b0"), g("n_b1")
ndet, ntimes, mult = g("n_det"), g("n_times"), g("max_mult")
dur, fa, pha = g("dur_us", float), g("fa", float), g("pha_med")
lon, lat = g("lon", float), g("lat", float)
frac = mult / np.maximum(nw, 1)

sig = np.where(fa <= 1e-5)[0]
sig = sig[np.argsort(fa[sig])]
print("fa<=1e-5 的 %d 个候选：" % len(sig))
print("  %-19s %8s %5s %4s %3s %3s %4s %4s %4s %5s %7s %4s %8s %7s" % (
    "start", "fa", "count", "NaI", "b0", "b1", "ndet", "ntim", "mult", "mult/n", "dur_us", "pha", "lon", "lat"))
for i in sig:
    print("  %-19s %8.1e %5d %4d %3d %3d %4d %4d %4d %5.2f %7.1f %4d %8.1f %7.1f" % (
        rows[i]["start"][:19], fa[i], cnt[i], nai[i], b0[i], b1[i], ndet[i], ntimes[i],
        mult[i], frac[i], dur[i], pha[i], lon[i], lat[i]))

print()
print("--- 组合判据（全体 47904 / 显著 41）---")
grp2 = np.sort(np.stack([nai, b0, b1], 1), 1)[:, -2]
tests = [
    ("mult/n_win <= 0.35", frac <= 0.35),
    ("mult/n_win <= 0.5 ", frac <= 0.5),
    ("dur>=50us          ", dur >= 50),
    ("dur>=50us & 组内>=3  ", (dur >= 50) & (grp2 >= 3)),
    ("dur>=50us & frac<=.35", (dur >= 50) & (frac <= 0.35)),
    ("组内>=3 & frac<=0.35 ", (grp2 >= 3) & (frac <= 0.35)),
    ("ntimes>=8          ", ntimes >= 8),
    ("ntimes>=8 & 组内>=3  ", (ntimes >= 8) & (grp2 >= 3)),
]
s = fa <= 1e-5
for name, m in tests:
    print("  %-22s : %6d (%5.2f%%) | %3d / 41" % (name, m.sum(), 100*m.mean(), (m & s).sum()))
