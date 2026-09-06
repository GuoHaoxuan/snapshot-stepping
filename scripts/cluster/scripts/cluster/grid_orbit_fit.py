"""天格轨道自定位：太阳同步圆轨道的轨道面从有位置解的数据定标，缺位置解时用本底计数率随磁纬的变化拟合沿轨相位。

子命令:
  calib    <sat> <day YYYY/MM/DD>... -o params.json   # 轨道面、周期、高度 + 本底率-偶极磁纬模板
  validate <sat> <day> params.json                    # 盲拟合有真值的一天，报告位置误差
  fit      <sat> <day> params.json -o orbit.csv       # 无位置解的一天：输出 10 s 一行的 time,lon,lat,alt_m
"""
import argparse, glob, json, os, sys
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import get_sun
from datetime import datetime, timezone, timedelta

G = "/gecamfs/Exchange/GSDC/missions/GRID"; REF = datetime(2018, 1, 1, tzinfo=timezone.utc)
R_E = 6371.0
BIN = 10.0
# 地磁偶极（IGRF-13 2020 近似）：北磁极 80.7°N, 72.7°W
POLE_LAT, POLE_LON = np.radians(80.7), np.radians(-72.7)


def files(sat, prod, day):
    d = f"{G}/{sat}/{prod}/{day}"
    if not os.path.isdir(d): return []
    v = sorted(os.listdir(d))[-1]; return sorted(glob.glob(f"{d}/{v}/*.fits"))


def dipole_lat(lat_deg, lon_deg):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    s = np.sin(lat) * np.sin(POLE_LAT) + np.cos(lat) * np.cos(POLE_LAT) * np.cos(lon - POLE_LON)
    return np.degrees(np.arcsin(np.clip(s, -1, 1)))


def in_saa(lat, lon):
    return (lon > -90) & (lon < 20) & (lat > -45) & (lat < 0)


def load_rates(sat, day):
    """每次过境的 10 s 计数率序列（四路 EVT_TYPE=1 合计）。"""
    out = []
    for f in files(sat, "fits7", day):
        with fits.open(f) as h:
            g = h["GTI"].data; a, b = float(g["START"][0]), float(g["STOP"][0])
            T = np.concatenate([np.asarray(h[f"EVENTS{k}"].data["TIME"], float)[np.asarray(h[f"EVENTS{k}"].data["EVT_TYPE"]) == 1] for k in range(4)])
        edges = np.arange(a, b + BIN, BIN); cnt, _ = np.histogram(T, bins=edges)
        out.append((edges[:-1] + BIN / 2, cnt / BIN))
    return out


def load_positions(sat, day):
    T = []; LAT = []; LON = []; ALT = []; XYZ = []
    for f in files(sat, "fits8", day):
        with fits.open(f) as h:
            d = h[1].data
            T.append(np.asarray(d["TIME"], float)); LAT.append(np.asarray(d["Latitude"], float)); LON.append(np.asarray(d["Longitude"], float)); ALT.append(np.asarray(d["Altitude"], float))
            XYZ.append(np.stack([np.asarray(d["X_J2000"], float), np.asarray(d["Y_J2000"], float), np.asarray(d["Z_J2000"], float)], axis=1))
    if not T: return None
    T = np.concatenate(T); o = np.argsort(T)
    return T[o], np.concatenate(LAT)[o], np.concatenate(LON)[o], np.concatenate(ALT)[o], np.concatenate(XYZ)[o]


def gmst_rad(met):
    t = Time(REF.timestamp() + met, format="unix", scale="utc")
    return t.sidereal_time("mean", "greenwich").rad


def sun_ra_rad(met):
    t = Time(REF.timestamp() + met, format="unix", scale="utc")
    return get_sun(t).ra.rad


def model_latlon(met, p, t0, period):
    """圆轨道：升交点赤经跟着太阳走（太阳同步），沿轨相位 u = 2π (t − t0)/P。"""
    inc = np.radians(p["inc_deg"]); raan = sun_ra_rad(met) + np.radians(p["ltan_h"] * 15.0 - 180.0)
    u = 2 * np.pi * (met - t0) / period
    x = np.cos(u) * np.cos(raan) - np.sin(u) * np.sin(raan) * np.cos(inc)
    y = np.cos(u) * np.sin(raan) + np.sin(u) * np.cos(raan) * np.cos(inc)
    z = np.sin(u) * np.sin(inc)
    lat = np.degrees(np.arcsin(z)); lon = np.degrees(np.arctan2(y, x) - gmst_rad(met))
    return lat, ((lon + 180) % 360) - 180


def calib(sat, days, out):
    incs = []; ltans = []; periods = []; alts = []; tmpl_x = []; tmpl_y = []; saa_r = []
    for day in days:
        pos = load_positions(sat, day)
        if pos is None: continue
        T, LAT, LON, ALT, XYZ = pos; ok = np.isfinite(LAT) & np.isfinite(XYZ[:, 0])
        T, LAT, LON, ALT, XYZ = T[ok], LAT[ok], LON[ok], ALT[ok], XYZ[ok]
        if len(T) < 100: continue
        # 角动量方向 → 倾角、升交点赤经；只用相邻 1 s 的样本
        d = np.diff(T); m = np.where((d > 0.5) & (d < 1.5))[0]
        r = XYZ[m]; v = (XYZ[m + 1] - XYZ[m]) / d[m][:, None]
        hvec = np.cross(r, v); hn = hvec / np.linalg.norm(hvec, axis=1, keepdims=True)
        inc = np.degrees(np.arccos(hn[:, 2])); raan = np.degrees(np.arctan2(hn[:, 0], -hn[:, 1]))
        ra_sun = np.degrees(sun_ra_rad(T[m]))
        ltan = ((raan - ra_sun + 180.0) / 15.0) % 24
        incs.append(np.median(inc)); ltans.append(np.median(ltan)); alts.append(np.median(ALT) / 1e3)
        asc = np.where((LAT[:-1] < 0) & (LAT[1:] >= 0) & (np.diff(T) < 20))[0]
        if len(asc) > 1:
            dt = np.diff(T[asc]); dt = dt[(dt > 5000) & (dt < 6500)]
            if len(dt): periods.append(np.median(dt))
        # 模板：本底率 vs 偶极磁纬（SAA 单独）
        for t, r_ in load_rates(sat, day):
            lat = np.interp(t, T, LAT); lon = np.interp(t, T, LON)
            ml = dipole_lat(lat, lon); s = in_saa(lat, lon)
            tmpl_x.extend(ml[~s]); tmpl_y.extend(r_[~s]); saa_r.extend(r_[s])
        print(f"  {day}: inc {np.median(inc):.2f}° ltan {np.median(ltan):.2f} h alt {np.median(ALT)/1e3:.1f} km 升交点数 {len(asc)}")
    tmpl_x = np.array(tmpl_x); tmpl_y = np.array(tmpl_y)
    edges = np.arange(-90, 91, 3.0); centers = edges[:-1] + 1.5; med = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mm = (tmpl_x >= lo) & (tmpl_x < hi); med.append(float(np.median(tmpl_y[mm])) if mm.sum() >= 5 else np.nan)
    params = {"sat": sat, "inc_deg": float(np.median(incs)), "ltan_h": float(np.median(ltans)), "period_s": float(np.median(periods)), "alt_km": float(np.median(alts)),
              "template_maglat": centers.tolist(), "template_rate": med, "saa_rate_median": float(np.median(saa_r)) if saa_r else None, "days": days}
    json.dump(params, open(out, "w"), indent=1)
    print("params:", {k: v for k, v in params.items() if not k.startswith("template")})
    print("template (maglat: rate):", " ".join(f"{c:.0f}:{v:.0f}" for c, v in zip(centers, med) if np.isfinite(v)))


def predict_rate(p, lat, lon):
    x = np.array(p["template_maglat"]); y = np.array(p["template_rate"], float); ok = np.isfinite(y)
    r = np.interp(dipole_lat(lat, lon), x[ok], y[ok])
    if p.get("saa_rate_median"): r = np.where(in_saa(lat, lon), p["saa_rate_median"], r)
    return r


def fit_day(sat, day, p, verbose=True):
    """全天所有过境共用 (t0, period)：先粗扫相位，再细化。返回 (t0, period, 损失)。"""
    passes = load_rates(sat, day)
    if not passes: return None
    T = np.concatenate([t for t, _ in passes]); R = np.concatenate([r for _, r in passes])
    keep = R > 0; T, R = T[keep], R[keep]; logR = np.log(R)
    ref = float(np.floor(T.min() / p["period_s"]) * p["period_s"])
    # 预先算太阳赤经与 GMST（与相位无关），避免每个试探重复算
    def loss(t0, period):
        lat, lon = model_latlon(T, p, t0, period)
        pr = np.log(np.maximum(predict_rate(p, lat, lon), 1.0))
        # 每次过境允许一个增益常数（探测器状态可能不同）：按过境去均值后比对
        resid = []
        i = 0
        for t, r in passes:
            n = int((r > 0).sum()); seg = slice(i, i + n); i += n
            d = logR[seg] - pr[seg]; resid.append(d - np.median(d))
        return float(np.sqrt(np.mean(np.concatenate(resid) ** 2)))
    best = (None, None, np.inf)
    for period in (p["period_s"],):
        for t0 in np.arange(ref, ref + period, 20.0):
            l = loss(t0, period)
            if l < best[2]: best = (t0, period, l)
    t0, period, l = best
    for dt0, dp in ((5.0, 2.0), (1.0, 0.5)):
        for pp in np.arange(period - 5 * dp, period + 5 * dp + 1e-9, dp):
            for tt in np.arange(t0 - 5 * dt0, t0 + 5 * dt0 + 1e-9, dt0):
                l = loss(tt, pp)
                if l < best[2]: best = (tt, pp, l)
        t0, period, l = best
    if verbose: print(f"  {day}: t0 = ref + {t0 - ref:.1f} s, period {period:.1f} s, loss {l:.3f}, passes {len(passes)}, bins {len(T)}")
    return best


def validate(sat, day, params):
    p = json.load(open(params)); pos = load_positions(sat, day)
    if pos is None: print("no positions"); return
    T, LAT, LON, ALT, _ = pos; ok = np.isfinite(LAT)
    best = fit_day(sat, day, p)
    if best is None: return
    t0, period, _ = best
    lat, lon = model_latlon(T[ok], p, t0, period)
    dlat = lat - LAT[ok]; dlon = ((lon - LON[ok] + 180) % 360) - 180
    dist = 2 * (R_E + p["alt_km"]) * np.arcsin(np.sqrt(np.sin(np.radians(dlat) / 2) ** 2 + np.cos(np.radians(LAT[ok])) * np.cos(np.radians(lat)) * np.sin(np.radians(dlon) / 2) ** 2))
    print(f"  验证 {sat} {day}: 位置误差 km 中位 {np.median(dist):.0f} p90 {np.percentile(dist, 90):.0f} max {dist.max():.0f}；样本 {ok.sum()}")
    # 也报一下若只用相位 0（不拟合）的误差，作对照
    return np.median(dist)


def fit(sat, day, params, out):
    p = json.load(open(params)); best = fit_day(sat, day, p)
    if best is None: print("no passes"); return
    t0, period, _ = best
    with open(out, "w") as f:
        f.write("time,lon,lat,alt_m\n")
        for t, _ in load_rates(sat, day):
            tt = np.arange(t[0] - BIN / 2, t[-1] + BIN / 2 + 1e-9, BIN)
            lat, lon = model_latlon(tt, p, t0, period)
            for a, b, c in zip(tt, lon, lat): f.write(f"{a:.1f},{b:.4f},{c:.4f},{p['alt_km']*1e3:.0f}\n")
    print("wrote", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("calib"); a.add_argument("sat"); a.add_argument("days", nargs="+"); a.add_argument("-o", required=True)
    b = sub.add_parser("validate"); b.add_argument("sat"); b.add_argument("day"); b.add_argument("params")
    c = sub.add_parser("fit"); c.add_argument("sat"); c.add_argument("day"); c.add_argument("params"); c.add_argument("-o", required=True)
    args = ap.parse_args()
    if args.cmd == "calib": calib(args.sat, args.days, args.o)
    elif args.cmd == "validate": validate(args.sat, args.day, args.params)
    else: fit(args.sat, args.day, args.params, args.o)
