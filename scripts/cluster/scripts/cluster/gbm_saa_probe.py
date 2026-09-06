"""GBM 试跑日（2019-01-01）的 SAA 边缘：TTE 的 GTI 是什么、事例流里的缺口在哪、缺口两端有没有 GTI 外的尾巴、候选离缺口多远。"""
from astropy.io import fits
import numpy as np, glob, json, sys
from datetime import datetime, timezone
D = "/hxmtfs/data/Fermi_GBM/2019/01/01/current"
DETS = ["n0","n1","n2","n3","n4","n5","n6","n7","n8","n9","na","nb","b0","b1"]
REF = datetime(2001, 1, 1, tzinfo=timezone.utc)
def met(iso):
    b = iso.rstrip("Z"); h, f = b.split("."); return (datetime.strptime(h + "." + (f + "000000")[:6], "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc) - REF).total_seconds() + 5.0  # 2019 年 GBM MET 比 UTC 多 5 个闰秒
# poshist SAA 段
ph = sorted(glob.glob("%s/glg_poshist_all_190101_v*.fit" % D))[-1]
with fits.open(ph) as h:
    pt = np.asarray(h[1].data["SCLK_UTC"], float); fl = np.asarray(h[1].data["FLAGS"])
saa = (fl >> 1) & 1 == 1
segs = []; i = 0
while i < len(saa):
    if saa[i]:
        j = i
        while j < len(saa) and saa[j]: j += 1
        segs.append((pt[i], pt[j-1] + 1.0)); i = j
    else: i += 1
print("poshist SAA 段 %d 个（总 %.0f s）:" % (len(segs), sum(b - a for a, b in segs)))
for a, b in segs: print("   %.0f .. %.0f  (%.0f s)" % (a, b, b - a))
# 候选
try:
    cands = json.load(open("/scratchfs2/gecam/guohx/gbmrun/data/Fermi_GBM/2019/01/20190101_signals.json"))
except Exception as e:
    cands = []; print("no candidates file:", e)
ct = np.array([met(c["start"]) for c in cands]); cfa = np.array([c["false_positive_per_year"] for c in cands])
print("候选 %d 个" % len(cands))
# 逐小时逐探头
edge_rows = []
for hh in range(24):
    for det in DETS:
        g = sorted(glob.glob("%s/glg_tte_%s_190101_%02dz_v*.fit.gz" % (D, det, hh)))
        if not g: continue
        with fits.open(g[-1]) as h:
            gti = [(float(a), float(b)) for a, b in zip(h["GTI"].data["START"], h["GTI"].data["STOP"])]
            t = np.asarray(h["EVENTS"].data["TIME"], float)
        if hh == 0 and det in ("n0", "b0"): print("  %s %02dz GTI 行数 %d: %s；事例 %d，首末 %.1f..%.1f" % (det, hh, len(gti), [(round(a, 1), round(b, 1)) for a, b in gti][:3], len(t), t[0], t[-1]))
        d = np.diff(t); gaps = np.where(d > 5.0)[0]
        for k in gaps:
            a, b = t[k], t[k+1]   # 缺口 [a, b]
            # 缺口是否在 GTI 里被排除：GTI 是否有一段结束在 a 附近、下一段起于 b 附近
            in_gti = any(abs(gs - b) < 2 or abs(ge - a) < 2 for gs, ge in gti)
            # 缺口两侧 2 s 与 20–22 s 的速率
            r_before_edge = ((t >= a - 2) & (t < a)).sum() / 2.0; r_before_far = ((t >= a - 22) & (t < a - 20)).sum() / 2.0
            r_after_edge = ((t > b) & (t <= b + 2)).sum() / 2.0; r_after_far = ((t > b + 20) & (t <= b + 22)).sum() / 2.0
            # 是否落在 poshist SAA 段
            in_saa = any(sa - 60 < a < sb + 60 or sa - 60 < b < sb + 60 for sa, sb in segs)
            edge_rows.append((hh, det, a, b, b - a, in_gti, in_saa, r_before_far, r_before_edge, r_after_edge, r_after_far))
print("事例缺口（>5 s）共 %d 个（14 路合计）" % len(edge_rows))
print("%3s %3s %12s %8s %6s %6s | %8s %8s %8s %8s" % ("h", "det", "gap_start", "len_s", "inGTI", "inSAA", "r-far", "r-edge", "r+edge", "r+far"))
for r in edge_rows[:60]:
    print("%3d %3s %12.1f %8.1f %6s %6s | %8.0f %8.0f %8.0f %8.0f" % r)
# 候选离最近缺口边的距离
if len(ct):
    edges = np.array(sorted(set([r[2] for r in edge_rows] + [r[3] for r in edge_rows])))
    dist = np.array([np.min(np.abs(edges - x)) if len(edges) else np.inf for x in ct])
    print("候选离最近缺口边的距离: <1 s %d, <2 s %d, <10 s %d, 总 %d；fa<=1e-5 的: <2 s %d / %d" % ((dist < 1).sum(), (dist < 2).sum(), (dist < 10).sum(), len(ct), ((dist < 2) & (cfa <= 1e-5)).sum(), (cfa <= 1e-5).sum()))
    for x, f_, dd in sorted(zip(ct, cfa, dist), key=lambda z: z[2])[:15]: print("   MET %.3f fa=%.1e 距缺口 %.2f s" % (x, f_, dd))
