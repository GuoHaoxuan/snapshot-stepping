"""天格轨道自定位：太阳同步圆轨道的轨道面从有位置解的数据定标，缺位置解时用本底计数率随磁纬的变化拟合沿轨相位。

子命令:
  calib    <sat> <day YYYY/MM/DD>... -o params.json   # 轨道面、周期、高度 + 本底率-偶极磁纬模板
  validate <sat> <day> params.json                    # 盲拟合有真值的一天，报告位置误差
  fit      <sat> <day> params.json -o orbit.csv       # 无位置解的一天：输出 10 s 一行的 time,lon,lat,alt_m
"""
import argparse, csv, glob, json, os, sys
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import get_sun
from astropy.utils import iers
iers.conf.auto_download = False   # 农场节点没有外网；用内置的 IERS-B，精度对本用途足够
from datetime import datetime, timezone, timedelta

G = "/gecamfs/Exchange/GSDC/missions/GRID"; REF = datetime(2018, 1, 1, tzinfo=timezone.utc)
R_E = 6371.0
BIN = 10.0
# 损失函数的口味：极区/极光带的率逐日起伏大（2024 年太阳活动高），模板在那里不稳，
# 只用 |磁纬| 不超过 MAGLAT_MAX 的格子；残差用 "rms" 或 "mad"（中位绝对偏差，抗离群）
LOSS = "rms"
MAGLAT_MAX = None
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
    """每次过境的 10 s 计数率序列（四路 EVT_TYPE=1 合计）。坏文件（如 0 字节）跳过。"""
    out = []
    for f in files(sat, "fits7", day):
        try:
            with fits.open(f) as h:
                g = h["GTI"].data; a, b = float(g["START"][0]), float(g["STOP"][0])
                T = np.concatenate([np.asarray(h[f"EVENTS{k}"].data["TIME"], float)[np.asarray(h[f"EVENTS{k}"].data["EVT_TYPE"]) == 1] for k in range(4)])
        except Exception as e:
            print("  skip", os.path.basename(f), e, file=sys.stderr); continue
        edges = np.arange(a, b + BIN, BIN); cnt, _ = np.histogram(T, bins=edges)
        out.append((edges[:-1] + BIN / 2, cnt / BIN))
    return out


def load_positions(sat, day):
    T = []; LAT = []; LON = []; ALT = []; XYZ = []
    for f in files(sat, "fits8", day):
        try:
            with fits.open(f) as h:
                d = h[1].data
                T.append(np.asarray(d["TIME"], float)); LAT.append(np.asarray(d["Latitude"], float)); LON.append(np.asarray(d["Longitude"], float)); ALT.append(np.asarray(d["Altitude"], float))
                XYZ.append(np.stack([np.asarray(d["X_J2000"], float), np.asarray(d["Y_J2000"], float), np.asarray(d["Z_J2000"], float)], axis=1))
        except Exception as e:
            print("  skip", os.path.basename(f), e, file=sys.stderr); continue
    if not T: return None
    T = np.concatenate(T); o = np.argsort(T)
    return T[o], np.concatenate(LAT)[o], np.concatenate(LON)[o], np.concatenate(ALT)[o], np.concatenate(XYZ)[o]


def gmst_rad(met):
    t = Time(REF.timestamp() + met, format="unix", scale="utc")
    return t.sidereal_time("mean", "greenwich").rad


def sun_ra_rad(met):
    t = Time(REF.timestamp() + met, format="unix", scale="utc")
    return get_sun(t).ra.rad


def model_latlon(met, p, t0, period, sky=None):
    """圆轨道：升交点赤经跟着太阳走（太阳同步），沿轨相位 u = 2π (t − t0)/P。

    `sky` = (sun_ra, gmst) 是与相位无关的量，拟合时对同一批时刻预先算一次
    （astropy 每次算都会试着下载 IERS 表，农场节点没有外网，白等超时）。"""
    sun_ra, gmst = sky if sky is not None else (sun_ra_rad(met), gmst_rad(met))
    # 升交点地方时不是严格常数（太阳同步只是近似），按定标窗测出的线性漂移外推
    ltan = p["ltan_h"] + p.get("ltan_rate_h_per_day", 0.0) * (met - p.get("ltan_epoch_met", met)) / 86400.0
    inc = np.radians(p["inc_deg"]); raan = sun_ra + np.radians(ltan * 15.0 - 180.0)
    u = 2 * np.pi * (met - t0) / period
    x = np.cos(u) * np.cos(raan) - np.sin(u) * np.sin(raan) * np.cos(inc)
    y = np.cos(u) * np.sin(raan) + np.sin(u) * np.cos(raan) * np.cos(inc)
    z = np.sin(u) * np.sin(inc)
    lat = np.degrees(np.arcsin(z)); lon = np.degrees(np.arctan2(y, x) - gmst)
    return lat, ((lon + 180) % 360) - 180


def stk_plane(path):
    """任务组的 STK 星历（source_orbit/orb_*.txt，10 s 一行）：第一组 xyz 是地固系、第二组是惯性系。
    返回 (历元 MET, 倾角, 升交点地方时, 周期, 高度)；没有事例，只能定轨道面，不能定模板。"""
    from datetime import datetime as DT
    T = []; XYZ = []; LAT = []; ALT = []
    for l in open(path, errors="ignore"):
        parts = l.split()
        if len(parts) < 15 or not parts[0].isdigit(): continue
        try:
            t = DT.strptime(" ".join(parts[:4]), "%d %b %Y %H:%M:%S.%f").replace(tzinfo=timezone.utc)
        except ValueError: continue
        T.append((t - REF).total_seconds()); LAT.append(float(parts[4])); ALT.append(float(parts[6]))
        XYZ.append([float(parts[12]), float(parts[13]), float(parts[14])])
    T = np.array(T); XYZ = np.array(XYZ); LAT = np.array(LAT); ALT = np.array(ALT)
    if len(T) < 100: return None
    d = np.diff(T); m = np.where((d > 5) & (d < 15))[0]
    r = XYZ[m]; v = (XYZ[m + 1] - XYZ[m]) / d[m][:, None]
    hvec = np.cross(r, v); hn = hvec / np.linalg.norm(hvec, axis=1, keepdims=True)
    inc = np.degrees(np.arccos(hn[:, 2])); raan = np.degrees(np.arctan2(hn[:, 0], -hn[:, 1]))
    ltan = ((raan - np.degrees(sun_ra_rad(T[m])) + 180.0) / 15.0) % 24
    asc = np.where((LAT[:-1] < 0) & (LAT[1:] >= 0))[0]; dt = np.diff(T[asc]); dt = dt[(dt > 5000) & (dt < 6500)]
    return float(np.median(T[m])), float(np.nanmedian(inc)), float(np.nanmedian(ltan)), float(np.median(dt)) if len(dt) else None, float(np.median(ALT))


def calib(sat, days, out, stk_files=()):
    incs = []; ltans = []; periods = []; alts = []; tmpl_x = []; tmpl_y = []; saa_r = []; anchors = []; ltan_epochs = []
    for path in stk_files:
        r = stk_plane(path)
        if r is None: print("  stk skip", os.path.basename(path)); continue
        epoch, inc, ltan, period, alt = r
        incs.append(inc); ltans.append(ltan); ltan_epochs.append(epoch); alts.append(alt)
        if period: periods.append((epoch, period))
        print(f"  {os.path.basename(path)}: inc {inc:.2f}° ltan {ltan:.2f} h alt {alt:.1f} km period {period}")
    # 模板只用最近的三个有事例的日子：本底率随太阳活动逐年变，早的日子形状也会走样
    template_days = set(sorted(days)[-3:])
    for day in days:
        pos = load_positions(sat, day)
        if pos is None: continue
        T, LAT, LON, ALT, XYZ = pos; ok = np.isfinite(LAT) & np.isfinite(XYZ[:, 0]) & np.isfinite(XYZ[:, 1]) & np.isfinite(XYZ[:, 2])
        T, LAT, LON, ALT, XYZ = T[ok], LAT[ok], LON[ok], ALT[ok], XYZ[ok]
        if len(T) < 100:
            print(f"  {day}: 有效 J2000 行 {len(T)}，跳过"); continue
        # 角动量方向 → 倾角、升交点赤经；只用相邻 1 s 的样本
        d = np.diff(T); m = np.where((d > 0.5) & (d < 1.5))[0]
        r = XYZ[m]; v = (XYZ[m + 1] - XYZ[m]) / d[m][:, None]
        hvec = np.cross(r, v); hn = hvec / np.linalg.norm(hvec, axis=1, keepdims=True)
        inc = np.degrees(np.arccos(hn[:, 2])); raan = np.degrees(np.arctan2(hn[:, 0], -hn[:, 1]))
        ra_sun = np.degrees(sun_ra_rad(T[m]))
        ltan = ((raan - ra_sun + 180.0) / 15.0) % 24
        if not (np.isfinite(np.nanmedian(inc)) and np.isfinite(np.nanmedian(ltan))):
            print(f"  {day}: 角动量算不出（速度差分含 NaN），跳过"); continue
        incs.append(np.nanmedian(inc)); ltans.append(np.nanmedian(ltan)); alts.append(np.nanmedian(ALT) / 1e3); ltan_epochs.append(float(np.median(T[m])))
        asc = np.where((LAT[:-1] < 0) & (LAT[1:] >= 0) & (np.diff(T) < 20))[0]
        if len(asc) > 1:
            dt = np.diff(T[asc]); dt = dt[(dt > 5000) & (dt < 6500)]
            if len(dt): periods.append((float(np.median(T[asc])), float(np.median(dt))))
        if len(asc):
            # 绝对相位锚：最后一个升交点的时刻（线性插到 lat=0）
            i = asc[-1]; frac = -LAT[i] / (LAT[i + 1] - LAT[i]); anchors.append(float(T[i] + frac * (T[i + 1] - T[i])))
        # 模板：本底率 vs 偶极磁纬（SAA 单独）
        if day not in template_days: continue
        for t, r_ in load_rates(sat, day):
            lat = np.interp(t, T, LAT); lon = np.interp(t, T, LON)
            ml = dipole_lat(lat, lon); s = in_saa(lat, lon)
            tmpl_x.extend(ml[~s]); tmpl_y.extend(r_[~s]); saa_r.extend(r_[s])
        print(f"  {day}: inc {np.median(inc):.2f}° ltan {np.median(ltan):.2f} h alt {np.median(ALT)/1e3:.1f} km 升交点数 {len(asc)}")
    if not incs:
        sys.exit(f"{sat}: 定标日里没有一天算得出轨道面，换几天")
    tmpl_x = np.array(tmpl_x); tmpl_y = np.array(tmpl_y)
    edges = np.arange(-90, 91, 3.0); centers = edges[:-1] + 1.5; med = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mm = (tmpl_x >= lo) & (tmpl_x < hi); med.append(float(np.median(tmpl_y[mm])) if mm.sum() >= 5 else np.nan)
    # 升交点地方时的线性漂移：定标日跨度够大（两段窗）才定斜率，否则当常数
    ltan_epochs = np.array(ltan_epochs); ltans_arr = np.array(ltans)
    if len(ltans_arr) >= 4 and (ltan_epochs.max() - ltan_epochs.min()) > 10 * 86400:
        x = (ltan_epochs - ltan_epochs.mean()) / 86400.0
        slope, intercept = np.polyfit(x, ltans_arr, 1); res1 = ltans_arr - (intercept + slope * x)
        ltan_h, ltan_rate, ltan_epoch = float(intercept), float(slope), float(ltan_epochs.mean())
        print(f"  升交点地方时漂移 {ltan_rate*60:+.2f} min/day（跨 {(x.max()-x.min()):.0f} 天，线性残差 rms {np.sqrt(np.mean(res1**2))*60:.2f} min）")
        if len(ltans_arr) >= 6 and (x.max() - x.min()) > 120:
            c2 = np.polyfit(x, ltans_arr, 2); res2 = ltans_arr - np.polyval(c2, x)
            print(f"  二次拟合残差 rms {np.sqrt(np.mean(res2**2))*60:.2f} min（二次项 {c2[0]*60:+.4f} min/day²）——只作参考，模型仍用线性")
    else:
        ltan_h, ltan_rate, ltan_epoch = float(np.median(ltans_arr)), 0.0, float(np.median(ltan_epochs)) if len(ltan_epochs) else 0.0
    params = {"sat": sat, "inc_deg": float(np.median(incs)), "ltan_h": ltan_h, "ltan_rate_h_per_day": ltan_rate, "ltan_epoch_met": ltan_epoch,
              # 周期取最晚三天的中位（阻力让它随时间变，早的窗不代表现在）
              "period_s": float(np.median([pp for _, pp in sorted(periods)[-3:]])), "alt_km": float(np.median(alts)),
              "template_maglat": centers.tolist(), "template_rate": med, "saa_rate_median": float(np.median(saa_r)) if saa_r else None, "days": days,
              # 相位锚：定标日里最晚的升交点；之后每天的 t0 以前一天为先验、只在 ±600 s 内找，
              # 防止半圈翻转（南北对称的模板下另一解的 loss 只差一点，GRID-07 盲拟合就翻到了对面）
              "t0_anchor": float(max(anchors)) if anchors else None}
    json.dump(params, open(out, "w"), indent=1)
    print("params:", {k: v for k, v in params.items() if not k.startswith("template")})
    print("template (maglat: rate):", " ".join(f"{c:.0f}:{v:.0f}" for c, v in zip(centers, med) if np.isfinite(v)))


def predict_rate(p, lat, lon):
    x = np.array(p["template_maglat"]); y = np.array(p["template_rate"], float); ok = np.isfinite(y)
    r = np.interp(dipole_lat(lat, lon), x[ok], y[ok])
    if p.get("saa_rate_median"): r = np.where(in_saa(lat, lon), p["saa_rate_median"], r)
    return r


def fit_day(sat, day, p, verbose=True, period_prior=None, t0_prior=None):
    """全天所有过境共用 (t0, period)：先粗扫相位，再细化。返回 (t0, period, 损失)。

    周期随大气阻力慢慢变短，`period_prior`（前一天拟出的周期）给出搜索中心，±20 s 内找。
    `t0_prior`（前一天的 t0，或定标的升交点锚）给出相位中心，只在 ±600 s 内找：
    模板南北近似对称，另一解只差半圈、loss 只差一点，没有先验时可能翻到对面。"""
    passes = load_rates(sat, day)
    if not passes: return None
    T = np.concatenate([t for t, _ in passes]); R = np.concatenate([r for _, r in passes])
    keep = R > 0; T, R = T[keep], R[keep]; logR = np.log(R)
    ref = float(np.floor(T.min() / p["period_s"]) * p["period_s"])
    sky = (sun_ra_rad(T), gmst_rad(T))
    def loss(t0, period):
        lat, lon = model_latlon(T, p, t0, period, sky)
        pr = np.log(np.maximum(predict_rate(p, lat, lon), 1.0))
        use = np.ones(len(T), bool) if MAGLAT_MAX is None else (np.abs(dipole_lat(lat, lon)) <= MAGLAT_MAX)
        # 每次过境允许一个增益常数（探测器状态可能不同）：按过境去均值后比对
        resid = []
        i = 0
        for t, r in passes:
            n = int((r > 0).sum()); seg = slice(i, i + n); i += n
            m = use[seg]
            if m.sum() < 3: continue
            d = (logR[seg] - pr[seg])[m]; resid.append(d - np.median(d))
        if not resid: return np.inf
        res = np.concatenate(resid)
        return float(np.median(np.abs(res))) if LOSS == "mad" else float(np.sqrt(np.mean(res ** 2)))
    best = (None, None, np.inf)
    center = period_prior if period_prior else p["period_s"]
    if t0_prior is not None:
        t0_grid = np.arange(t0_prior - 600.0, t0_prior + 600.0 + 1e-9, 10.0)
    else:
        t0_grid = None
    # 周期一天能变 1 s（2024 年太阳活动高、440 km 的星阻力大），搜索范围要够
    for period in np.arange(center - 24.0, center + 24.0 + 1e-9, 3.0):
        grid = t0_grid if t0_grid is not None else np.arange(ref, ref + period, 20.0)
        for t0 in grid:
            l = loss(t0, period)
            if l < best[2]: best = (t0, period, l)
    t0, period, l = best
    if t0 is None:
        print(f"  {day}: 损失全为无穷（模板/参数有 NaN？），放弃", file=sys.stderr); return None
    for dt0, dp in ((5.0, 1.0), (1.0, 0.25)):
        for pp in np.arange(period - 5 * dp, period + 5 * dp + 1e-9, dp):
            for tt in np.arange(t0 - 5 * dt0, t0 + 5 * dt0 + 1e-9, dt0):
                l = loss(tt, pp)
                if l < best[2]: best = (tt, pp, l)
        t0, period, l = best
    # 把相位历元挪到这一天的末尾：u = 2π(t − t0)/P 里 t0 离数据越远，周期误差的
    # 杠杆越长（6 天 92 圈，周期差 8 s 就是 740 s 的相位差，超出 ±600 s 的搜索）
    t0 = t0 + period * np.floor((T.max() - t0) / period)
    if verbose: print(f"  {day}: t0 = ref + {t0 - ref:.1f} s, period {period:.1f} s, loss {l:.3f}, passes {len(passes)}, bins {len(T)}")
    return (t0, period, l)


def validate(sat, day, params):
    """像生产那样从定标锚的次日起逐天跟踪相位到 `day`，再拿 `day` 的真实位置量误差。"""
    p = json.load(open(params)); pos = load_positions(sat, day)
    if pos is None: print("no positions"); return
    T, LAT, LON, ALT, _ = pos; ok = np.isfinite(LAT)
    from datetime import date
    t0_prior = p.get("t0_anchor"); period_prior = None
    if t0_prior is not None:
        d = (REF + timedelta(seconds=t0_prior)).date() + timedelta(days=1); target = date(*map(int, day.split("/")))
        while d < target:
            b = fit_day(sat, d.strftime("%Y/%m/%d"), p, verbose=False, period_prior=period_prior, t0_prior=t0_prior)
            if b is not None: t0_prior, period_prior, _ = b; print(f"    跟踪 {d}: period {period_prior:.1f} loss {b[2]:.3f}")
            d += timedelta(days=1)
    best = fit_day(sat, day, p, period_prior=period_prior, t0_prior=t0_prior)
    if best is None: return
    t0, period, _ = best
    lat, lon = model_latlon(T[ok], p, t0, period)
    dlat = lat - LAT[ok]; dlon = ((lon - LON[ok] + 180) % 360) - 180
    dist = 2 * (R_E + p["alt_km"]) * np.arcsin(np.sqrt(np.sin(np.radians(dlat) / 2) ** 2 + np.cos(np.radians(LAT[ok])) * np.cos(np.radians(lat)) * np.sin(np.radians(dlon) / 2) ** 2))
    print(f"  验证 {sat} {day}: 位置误差 km 中位 {np.median(dist):.0f} p90 {np.percentile(dist, 90):.0f} max {dist.max():.0f}；样本 {ok.sum()}")
    # 也报一下若只用相位 0（不拟合）的误差，作对照
    return np.median(dist)


def write_orbit(sat, day, p, t0, period, out):
    with open(out, "w") as f:
        f.write("time,lon,lat,alt_m\n")
        for t, _ in load_rates(sat, day):
            tt = np.arange(t[0] - BIN / 2, t[-1] + BIN / 2 + 1e-9, BIN)
            lat, lon = model_latlon(tt, p, t0, period)
            for a, b, c in zip(tt, lon, lat): f.write(f"{a:.1f},{b:.4f},{c:.4f},{p['alt_km']*1e3:.0f}\n")


def fit(sat, day, params, out):
    p = json.load(open(params)); best = fit_day(sat, day, p)
    if best is None: print("no passes"); return
    t0, period, _ = best
    write_orbit(sat, day, p, t0, period, out); print("wrote", out)


def fit_range(sat, start, end, params, outdir):
    """逐天拟合并写 <outdir>/<sat>/YYYYMMDD.csv；周期以前一天的拟合值为先验；日志一行一天（loss 高的天要人看）。"""
    from datetime import date
    p = json.load(open(params)); os.makedirs(os.path.join(outdir, sat), exist_ok=True)
    d0 = date(*map(int, start.split("-"))); d1 = date(*map(int, end.split("-"))); prior = None; t0_prior = p.get("t0_anchor")
    logpath = os.path.join(outdir, sat, "fit_log.csv"); log = open(logpath, "a")
    if os.path.getsize(logpath) == 0: log.write("day,t0_offset,period,loss,passes,t0_abs\n")
    d = d0
    while d <= d1:
        day = d.strftime("%Y/%m/%d"); out = os.path.join(outdir, sat, d.strftime("%Y%m%d") + ".csv")
        if not os.path.exists(out):
            best = fit_day(sat, day, p, verbose=False, period_prior=prior, t0_prior=t0_prior)
            if best is not None:
                t0, period, l = best; prior = period; t0_prior = t0
                write_orbit(sat, day, p, t0, period, out)
                ref = float(np.floor(t0 / period) * period)
                log.write(f"{d.strftime('%Y-%m-%d')},{t0 - ref:.1f},{period:.2f},{l:.3f},{len(load_rates(sat, day))},{t0:.2f}\n"); log.flush()
                print(f"{day}: period {period:.1f} loss {l:.3f}", flush=True)
        d += timedelta(days=1)


def smooth(sat, params, outdir, max_loss=None, max_resid=120.0, degree=2):
    """把逐日拟合的相位历元串成一条平滑的轨道：升交点时刻 T 对圈数 N 做多项式（周期 = dT/dN 单调缓变），
    剔除 loss 高或偏离曲线的天，再用平滑模型重写每天的轨道表；被剔除的天删掉表（候选按无星历丢，不造假位置）。"""
    p = json.load(open(params)); d = os.path.join(outdir, sat)
    rows = [r for r in csv.DictReader(open(os.path.join(d, "fit_log.csv"))) if r["day"] != "day" and r.get("t0_abs")]
    rows.sort(key=lambda r: r["day"])
    T = np.array([float(r["t0_abs"]) for r in rows]); L = np.array([float(r["loss"]) for r in rows]); P = np.array([float(r["period"]) for r in rows])
    # 圈数：相邻历元之差除以当时的周期取整，累加
    N = np.zeros(len(T)); 
    for i in range(1, len(T)): N[i] = N[i - 1] + np.round((T[i] - T[i - 1]) / P[i - 1])
    # loss 门槛自适应：中位 + 3 倍中位绝对偏差（不同星、不同季节的 loss 水平不同），上限 0.9
    if max_loss is None:
        max_loss = min(0.9, float(np.median(L) + 3 * np.median(np.abs(L - np.median(L)))))
    good = L <= max_loss
    print(f"  loss 门槛 {max_loss:.3f}，超出的天 {int((~good).sum())}")
    for _ in range(6):
        if good.sum() < degree + 2: break
        coef = np.polyfit(N[good], T[good], degree); resid = T - np.polyval(coef, N)
        new_good = good & (np.abs(resid) <= max_resid)
        if new_good.sum() == good.sum(): break
        good = new_good
    coef = np.polyfit(N[good], T[good], degree); resid = T - np.polyval(coef, N)
    dcoef = np.polyder(coef)
    print(f"{sat}: {len(rows)} 天，保留 {good.sum()}；残差 |中位| {np.median(np.abs(resid[good])):.1f} s，周期 {np.polyval(dcoef, N[0]):.1f} → {np.polyval(dcoef, N[-1]):.1f} s")
    with open(os.path.join(d, "smooth_log.csv"), "w") as f:
        f.write("day,loss,period_fit,period_smooth,resid_s,kept\n")
        for r, l, pp, n, rs, g in zip(rows, L, P, N, resid, good):
            f.write(f"{r['day']},{l:.3f},{pp:.1f},{np.polyval(dcoef, n):.2f},{rs:.1f},{int(g)}\n")
    for r, n, g in zip(rows, N, good):
        y, m, dd = r["day"].split("-"); out = os.path.join(d, f"{y}{m}{dd}.csv")
        if not g:
            if os.path.exists(out): os.remove(out)
            continue
        t0 = float(np.polyval(coef, n)); period = float(np.polyval(dcoef, n))
        write_orbit(sat, f"{y}/{m}/{dd}", p, t0, period, out)
    print("  剔除的天:", [r["day"] for r, g in zip(rows, good) if not g][:20])


def validate_range(sat, start, end, params, outdir):
    """在有真值的一段上走完整的生产流程（逐日拟合 → 平滑），再逐天量误差。"""
    from datetime import date
    fit_range(sat, start, end, params, outdir)
    smooth(sat, params, outdir)
    p = json.load(open(params)); d0 = date(*map(int, start.split("-"))); d1 = date(*map(int, end.split("-"))); d = d0
    while d <= d1:
        path = os.path.join(outdir, sat, d.strftime("%Y%m%d") + ".csv")
        pos = load_positions(sat, d.strftime("%Y/%m/%d"))
        if os.path.exists(path) and pos is not None:
            T, LAT, LON, ALT, _ = pos; ok = np.isfinite(LAT)
            tab = np.loadtxt(path, delimiter=",", skiprows=1)
            if ok.sum() > 10 and len(tab) > 2:
                lon = np.interp(T[ok], tab[:, 0], np.unwrap(np.radians(tab[:, 1]))); lon = np.degrees(lon); lat = np.interp(T[ok], tab[:, 0], tab[:, 2])
                inside = (T[ok] >= tab[0, 0]) & (T[ok] <= tab[-1, 0])
                def gc(lat1, lon1, lat2, lon2):
                    return 2 * (R_E + p["alt_km"]) * np.arcsin(np.sqrt(np.sin(np.radians(lat1 - lat2) / 2) ** 2 + np.cos(np.radians(lat2)) * np.cos(np.radians(lat1)) * np.sin(np.radians(((lon1 - lon2 + 180) % 360) - 180) / 2) ** 2))
                dist = gc(lat, lon, LAT[ok], LON[ok])[inside]
                if len(dist):
                    # 沿轨/横轨分解：把表的时间平移 δ 再比，取误差最小的 δ 就是沿轨（相位）误差，剩下的是横轨（轨道面）误差
                    best = (np.median(dist), 0.0)
                    for delta in np.arange(-60, 60.1, 2.0):
                        lo2 = np.degrees(np.interp(T[ok] + delta, tab[:, 0], np.unwrap(np.radians(tab[:, 1])))); la2 = np.interp(T[ok] + delta, tab[:, 0], tab[:, 2])
                        m2 = np.median(gc(la2, lo2, LAT[ok], LON[ok])[inside])
                        if m2 < best[0]: best = (m2, delta)
                    print(f"  {d}: 误差中位 {np.median(dist):.0f} km p90 {np.percentile(dist, 90):.0f} max {dist.max():.0f}；沿轨相位差 {best[1]:+.0f} s（{best[1]*7.6:+.0f} km），去掉后横轨中位 {best[0]:.0f} km")
        elif pos is not None: print(f"  {d}: 无轨道表（被剔除或没数据）")
        d += timedelta(days=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("calib"); a.add_argument("sat"); a.add_argument("days", nargs="+"); a.add_argument("-o", required=True); a.add_argument("--stk", nargs="*", default=[], help="任务组 STK 星历文件，只定轨道面")
    b = sub.add_parser("validate"); b.add_argument("sat"); b.add_argument("day"); b.add_argument("params")
    c = sub.add_parser("fit"); c.add_argument("sat"); c.add_argument("day"); c.add_argument("params"); c.add_argument("-o", required=True)
    r = sub.add_parser("fit-range"); r.add_argument("sat"); r.add_argument("start"); r.add_argument("end"); r.add_argument("params"); r.add_argument("-o", required=True)
    sm = sub.add_parser("smooth"); sm.add_argument("sat"); sm.add_argument("params"); sm.add_argument("-o", required=True)
    vr = sub.add_parser("validate-range"); vr.add_argument("sat"); vr.add_argument("start"); vr.add_argument("end"); vr.add_argument("params"); vr.add_argument("-o", required=True)
    args = ap.parse_args()
    if args.cmd == "smooth": smooth(args.sat, args.params, args.o); sys.exit(0)
    if args.cmd == "validate-range": validate_range(args.sat, args.start, args.end, args.params, args.o); sys.exit(0)
    if args.cmd == "calib": calib(args.sat, args.days, args.o, args.stk)
    elif args.cmd == "validate": validate(args.sat, args.day, args.params)
    elif args.cmd == "fit": fit(args.sat, args.day, args.params, args.o)
    else: fit_range(args.sat, args.start, args.end, args.params, args.o)
