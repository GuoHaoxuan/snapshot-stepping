"""从 L1B 的 GTI 缺口反推 SVOM 实际使用的 SAA 区域。

SVOM 的 L1B 把 SAA 期间排除在 GTI 之外（这正是早先两类假信号的根因），所以
「在轨但不在 GTI 内」的时刻对应的星下点，就是这套数据里 SAA 的操作定义。
不依赖任何外部模型或文档。

用法: python3 svom_saa_region.py <n_days> <out.csv>
"""
from astropy.io import fits
import numpy as np, glob, os, sys, csv
from datetime import datetime, timedelta, timezone

D = "/gecamfs/SVOM/Archived-DATA/GRM-DATA/L1B/daily"
STEP = 10.0     # 每 10 s 采一个点


def hours_of(day):
    """每个整点只取版本号最高的那个事例文件（同一小时会有 v00/v01/v02）。"""
    best = {}
    for f in glob.glob(f"{D}/{day}/grm_evt/svom_grm_evt_*_v*.fits"):
        parts = os.path.basename(f).split("_")
        best[parts[4]] = max(best.get(parts[4], ""), f)
    return [best[k] for k in sorted(best)]


def main(n_days, out):
    # 只取真有事例文件的天：早期有些日目录只有 att/orb/hk
    days = [d for d in sorted(glob.glob(f"{D}/2024/*/*")) if glob.glob(f"{d}/grm_evt/*.fits") and glob.glob(f"{d}/orb/*.fits")][:n_days]
    in_saa, in_gti = [], []
    for dpath in days:
        day = "/".join(dpath.split("/")[-3:])
        for f in hours_of(day):
            base = os.path.basename(f)
            hh = base.split("_")[4]
            orbs = sorted(glob.glob(f"{D}/{day}/orb/svom_orb_{base.split('_')[3]}_{hh}_v*.fits"))
            if not orbs: continue
            try:
                with fits.open(f) as h:
                    gti = [(float(a), float(b)) for a, b in zip(h["GTI"].data["START"], h["GTI"].data["STOP"])]
                with fits.open(orbs[-1]) as h:
                    d = h[1].data
                    cols = {c.upper(): c for c in d.columns.names}
                    t = np.asarray(d[cols["TIME"]], float)
                    lon = np.asarray(d[cols.get("LONGITUDE", cols.get("LON", "LONGITUDE"))], float)
                    lat = np.asarray(d[cols.get("LATITUDE", cols.get("LAT", "LATITUDE"))], float)
            except Exception as e:
                print("skip", base, e, file=sys.stderr); continue
            if len(t) < 10 or not gti: continue
            grid = np.arange(t.min(), t.max(), STEP)
            inside = np.zeros(len(grid), bool)
            for a, b in gti: inside |= (grid >= a) & (grid <= b)
            glon = np.interp(grid, t, np.unwrap(np.radians(lon))); glon = ((np.degrees(glon) + 180) % 360) - 180
            glat = np.interp(grid, t, lat)
            in_saa.extend(zip(glon[~inside], glat[~inside]))
            in_gti.extend(zip(glon[inside][::20], glat[inside][::20]))
    print(f"{len(days)} 天：GTI 外采样点 {len(in_saa)}，GTI 内（抽样）{len(in_gti)}")
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["lon", "lat", "in_gti"])
        for lo, la in in_saa: w.writerow(["%.3f" % lo, "%.3f" % la, 0])
        for lo, la in in_gti: w.writerow(["%.3f" % lo, "%.3f" % la, 1])
    if in_saa:
        a = np.array(in_saa)
        print("GTI 外点的经度 %.1f..%.1f，纬度 %.1f..%.1f" % (a[:, 0].min(), a[:, 0].max(), a[:, 1].min(), a[:, 1].max()))
    print("wrote", out)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20, sys.argv[2] if len(sys.argv) > 2 else "saa_region.csv")
