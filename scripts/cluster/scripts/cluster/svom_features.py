"""逐候选提取事例级特征，用来搞清楚这些候选到底是什么。

对每个候选取三个时间窗：
  core     [start, stop]              —— 候选本身
  near     core 两侧各 1 s，扣掉 ±10 ms —— 本地本底
输出能谱、探头分布、反符合占比、持续时间等，供离线分类。
"""
from astropy.io import fits
import numpy as np, glob, csv, sys, os
from collections import defaultdict
from datetime import datetime, timezone

REF = datetime(2017, 1, 1, tzinfo=timezone.utc)
D = "/gecamfs/SVOM/Archived-DATA/GRM-DATA/L1B/daily"


def met(iso):
    # serde 会截掉小数末尾的零，秒的小数位数不固定，不能按定长切片
    body = iso.rstrip("Z")
    if "." in body:
        head, frac = body.split(".")
        body = head + "." + (frac + "000000")[:6]
    else:
        body += ".000000"
    t = datetime.strptime(body, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)
    return (t - REF).total_seconds(), t


rows = list(csv.DictReader(open(sys.argv[1])))
by_hour = defaultdict(list)
for r in rows:
    m, t = met(r["start"])
    by_hour[t.strftime("%Y/%m/%d %H")].append((r, m, met(r["stop"])[0], t))

out = open(sys.argv[2], "w", newline="")
w = csv.writer(out)
w.writerow(["start", "fa", "lon", "lat", "phase", "n_core", "dur_ms",
            "pi_core_med", "pi_bkg_med", "pi_med_ratio", "acd_core", "acd_bkg",
            "det_frac_max", "rate_bkg", "n_det_hit"])
done = 0
WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 1
IDX = int(sys.argv[4]) if len(sys.argv) > 4 else 0
for n_hour, (hour, items) in enumerate(sorted(by_hour.items())):
    if n_hour % WORKERS != IDX:
        continue
    day, hh = hour.split()
    pat = "%s/%s/grm_evt/svom_grm_evt_%s_%s_v*.fits" % (
        D, day, datetime.strptime(day, "%Y/%m/%d").strftime("%y%m%d"), hh)
    g = sorted(glob.glob(pat))
    if not g:
        print("  缺文件", pat, flush=True); continue
    with fits.open(g[-1]) as h:
        T, P, A, Dt = [], [], [], []
        for i, det in zip((3, 4, 5), (1, 2, 3)):
            d = h[i].data
            t = np.asarray(d["TIME"], dtype=np.float64)
            pi = np.asarray(d["PI"]); et = np.asarray(d["EVT_TYPE"])
            ac = np.asarray(d["ANTI_COIN"])
            k = (et == 0) & (pi >= 15) & (pi < 256)          # 与搜索同款 keep
            T.append(t[k]); P.append(pi[k]); A.append(ac[k])
            Dt.append(np.full(k.sum(), det, dtype=np.int8))
        T = np.concatenate(T); P = np.concatenate(P)
        A = np.concatenate(A); Dt = np.concatenate(Dt)
        o = np.argsort(T); T, P, A, Dt = T[o], P[o], A[o], Dt[o]
    for r, m0, m1, _ in items:
        core = (T >= m0) & (T <= m1)
        near = ((T >= m0 - 1.0) & (T <= m1 + 1.0)) & ~((T >= m0 - 0.01) & (T <= m1 + 0.01))
        nc = int(core.sum()); nb = int(near.sum())
        if nc == 0:
            continue
        dets, cnts = np.unique(Dt[core], return_counts=True)
        w.writerow([
            r["start"][:26], r["false_positive_per_year"], r["lon"], r["lat"],
            "%.6f" % (float("0" + r["start"][19:].rstrip("Z"))),
            nc, "%.4f" % ((m1 - m0) * 1e3),
            int(np.median(P[core])), int(np.median(P[near])) if nb else -1,
            "%.3f" % (np.median(P[core]) / max(np.median(P[near]), 1) if nb else -1),
            "%.4f" % A[core].mean(), "%.4f" % (A[near].mean() if nb else -1),
            "%.3f" % (cnts.max() / nc), "%.0f" % (nb / 2.0), len(dets),
        ])
    done += len(items)
    print("  %s  累计 %d/%d" % (hour, done, len(rows)), flush=True)
out.close()
print("done ->", sys.argv[2])
