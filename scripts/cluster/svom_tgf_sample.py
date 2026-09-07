"""导出 SVOM 闪电证实 TGF 样本的事例级数据，供离线做时长、能谱、死时间、地理分布分析。

每个 TGF 导出三样：
  events   最显著格中心 ±50 ms 内的逐事例（时间、能道、探头、增益档、死时间、反符合）
  bkg      本底窗（候选两侧各 1 s、扣掉 ±10 ms，并夹到 GTI 内）的能道直方与死时间合计
  meta     位置、本底计数、本底窗活时间、逐探头计数与死时间

口径与搜索侧一致：事例准入 EVT_TYPE==0（能道不在这里切，留给离线分析）；
本底窗与 `svom_features.py` 同款（±1 s 挖 ±10 ms）。

**ANTI_COIN==1 是星上标定源事例**（49–57 keV 的线，三路合计约 36 c/s），不是反符合
标志；本底谱按 AC==0 统计，AC==1 的谱单独导出一份用来量化它的影响。不扣的话，扣本底
时会在 50 keV 附近过扣，低能段的谱型会被带偏。

用法: python3 svom_tgf_sample.py <confirmed.csv> <outdir>
"""
from astropy.io import fits
import numpy as np, glob, csv, sys, os
from collections import defaultdict
from datetime import datetime, timezone

REF = datetime(2017, 1, 1, tzinfo=timezone.utc)
D = "/gecamfs/SVOM/Archived-DATA/GRM-DATA/L1B/daily"
HALF_EVENTS = 0.050      # 逐事例导出的半窗
BKG_HALF = 1.0           # 本底窗半宽
BKG_HOLLOW = 0.010       # 本底窗中间挖掉的半宽
NCHAN = 259


def met(iso):
    body = iso.rstrip("Z")
    if "." in body:
        head, frac = body.split("."); body = head + "." + (frac + "000000")[:6]
    else:
        body += ".000000"
    t = datetime.strptime(body, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)
    return (t - REF).total_seconds(), t


def load_hour(day, hh):
    """一小时三路事例合并（按时间排序）+ GTI + EBOUNDS。找不到文件返回 None。"""
    pat = "%s/%s/grm_evt/svom_grm_evt_%s_%s_v*.fits" % (
        D, day, datetime.strptime(day, "%Y/%m/%d").strftime("%y%m%d"), hh)
    g = sorted(glob.glob(pat))
    if not g:
        return None
    cols = {k: [] for k in ("t", "pi", "gain", "dead", "evt", "anti", "det")}
    with fits.open(g[-1]) as h:
        gti = [(float(a), float(b)) for a, b in zip(h["GTI"].data["START"], h["GTI"].data["STOP"])]
        ebounds = np.array([(float(r[1]), float(r[2])) for r in h["EBOUNDS"].data])
        for i, det in zip((3, 4, 5), (1, 2, 3)):
            d = h[i].data
            n = len(d["TIME"])
            cols["t"].append(np.asarray(d["TIME"], float))
            cols["pi"].append(np.asarray(d["PI"], np.int16))
            cols["gain"].append(np.asarray(d["GAIN_TYPE"], np.int8))
            cols["dead"].append(np.asarray(d["DEAD_TIME"], np.float32))
            cols["evt"].append(np.asarray(d["EVT_TYPE"], np.int8))
            cols["anti"].append(np.asarray(d["ANTI_COIN"], np.int8))
            cols["det"].append(np.full(n, det, np.int8))
    arr = {k: np.concatenate(v) for k, v in cols.items()}
    o = np.argsort(arr["t"], kind="stable")
    return {k: v[o] for k, v in arr.items()}, gti, ebounds


def live_length(gti, a, b):
    return sum(max(0.0, min(y, b) - max(x, a)) for x, y in gti)


def main(inp, outdir):
    os.makedirs(outdir, exist_ok=True)
    rows = list(csv.DictReader(open(inp)))
    by_hour = defaultdict(list)
    for r in rows:
        m0, t0 = met(r["start"])
        # 最显著格中心：start + delay + bin_size_best/2
        peak = m0 + float(r["delay"]) + float(r["bin_size_best"]) / 2.0
        by_hour[t0.strftime("%Y/%m/%d %H")].append((r, m0, met(r["stop"])[0], peak))

    ev = open(os.path.join(outdir, "sample_events.csv"), "w", newline="")
    we = csv.writer(ev)
    we.writerow(["idx", "dt_s", "pi", "det", "gain", "dead_us", "anti", "evt_type"])
    bk = open(os.path.join(outdir, "sample_bkg_spec.csv"), "w", newline="")
    wb = csv.writer(bk)
    wb.writerow(["idx", "live_s"] + ["c%d" % i for i in range(NCHAN)])
    bk1 = open(os.path.join(outdir, "sample_bkg_spec_ac1.csv"), "w", newline="")
    wb1 = csv.writer(bk1)
    wb1.writerow(["idx", "live_s"] + ["c%d" % i for i in range(NCHAN)])
    me = open(os.path.join(outdir, "sample_meta.csv"), "w", newline="")
    wm = csv.writer(me)
    wm.writerow(["idx", "start", "fa", "lon", "lat", "alt_m", "peak_met", "bkg_live_s", "bkg_counts", "bkg_dead_us",
                 "bkg_counts_d1", "bkg_counts_d2", "bkg_counts_d3",
                 "bkg_dead_us_d1", "bkg_dead_us_d2", "bkg_dead_us_d3", "n_gti_seg"])

    eb_written = False
    done = 0
    for hour, items in sorted(by_hour.items()):
        day, hh = hour.split()
        loaded = load_hour(day, hh)
        if loaded is None:
            print("missing", hour, flush=True)
            continue
        arr, gti, ebounds = loaded
        if not eb_written:
            with open(os.path.join(outdir, "ebounds.csv"), "w") as f:
                f.write("channel,e_min,e_max\n")
                for i, (lo, hi) in enumerate(ebounds):
                    f.write("%d,%.6f,%.6f\n" % (i, lo, hi))
            eb_written = True
        t = arr["t"]
        keep_evt = arr["evt"] == 0
        ingti = np.zeros(len(t), bool)
        for x, y in gti:
            ingti |= (t >= x) & (t <= y)
        for r, m0, m1, peak in items:
            idx = int(r["idx"])
            lo, hi = np.searchsorted(t, [peak - HALF_EVENTS, peak + HALF_EVENTS])
            for j in range(lo, hi):
                we.writerow([idx, "%.9f" % (t[j] - peak), int(arr["pi"][j]), int(arr["det"][j]),
                             int(arr["gain"][j]), "%.3f" % arr["dead"][j], int(arr["anti"][j]), int(arr["evt"][j])])
            a, b = m0 - BKG_HALF, m1 + BKG_HALF
            ha, hb = m0 - BKG_HOLLOW, m1 + BKG_HOLLOW
            live = live_length(gti, a, ha) + live_length(gti, hb, b)
            win = (((t >= a) & (t < ha)) | ((t > hb) & (t <= b))) & ingti & keep_evt
            sel = win & (arr["anti"] == 0)      # 标定源事例不计入本底
            spec = np.bincount(np.clip(arr["pi"][sel], 0, NCHAN - 1), minlength=NCHAN)
            wb.writerow([idx, "%.6f" % live] + [int(x) for x in spec])
            sel1 = win & (arr["anti"] == 1)
            spec1 = np.bincount(np.clip(arr["pi"][sel1], 0, NCHAN - 1), minlength=NCHAN)
            wb1.writerow([idx, "%.6f" % live] + [int(x) for x in spec1])
            dsel = arr["det"][sel]
            dead = arr["dead"][sel]
            wm.writerow([idx, r["start"], r["fa"], r["lon"], r["lat"], r["alt_m"], "%.6f" % peak,
                         "%.6f" % live, int(sel.sum()), "%.1f" % dead.sum()]
                        + [int((dsel == d).sum()) for d in (1, 2, 3)]
                        + ["%.1f" % dead[dsel == d].sum() for d in (1, 2, 3)]
                        + [len(gti)])
            done += 1
        print("hour", hour, "done", done, "/", len(rows), flush=True)
    for f in (ev, bk, bk1, me):
        f.close()
    print("wrote", done, "TGFs to", outdir)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
