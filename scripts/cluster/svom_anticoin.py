"""SVOM/GRM 的 ANTI_COIN 列是什么：取值、与其它列的联合分布、空间依赖、三路时间相关。

判据能不能用，取决于它是不是真的在标带电粒子。三条线索：
  1. 取值与联合分布——是位标志还是计数，跟 PI/GAIN_TYPE/EVT_TYPE/FLAG/DEAD_TIME 怎么关联；
  2. 空间依赖——按偶极磁纬与 SAA 分箱的触发率，带电粒子应在 SAA 与高磁纬显著升高；
  3. 时间相关——AC=1 的事例在另外两路里有没有同时刻的伙伴（穿整星的粒子会让多路同时响应），
     以 AC=0 的事例做对照。

用法: python3 svom_anticoin.py <YYYY/MM/DD> <hour_step> <out_prefix>
"""
import glob
import sys
from datetime import datetime, timedelta, timezone

import numpy as np
from astropy.io import fits

ARCHIVE = "/gecamfs/SVOM/Archived-DATA/GRM-DATA/L1B/daily"
REF = datetime(2017, 1, 1, tzinfo=timezone.utc)          # SVOM MET 零点
POLE_LAT, POLE_LON = np.radians(80.7), np.radians(-72.7)  # 偶极北磁极
COINC_US = 20.0                                           # 三路符合窗（µs）
MAX_PAIRS = 40000                                         # 相关性检验的抽样上限


def dipole_lat(lat_deg, lon_deg):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    s = np.sin(lat) * np.sin(POLE_LAT) + np.cos(lat) * np.cos(POLE_LAT) * np.cos(lon - POLE_LON)
    return np.degrees(np.arcsin(np.clip(s, -1, 1)))


def in_saa(lat, lon):
    return (lon > -90) & (lon < 40) & (lat > -40) & (lat < 5)


def latest(pattern):
    g = sorted(glob.glob(pattern))
    return g[-1] if g else None


def read_hour(day, hh):
    ymd = datetime.strptime(day, "%Y/%m/%d").strftime("%y%m%d")
    evt = latest(f"{ARCHIVE}/{day}/grm_evt/svom_grm_evt_{ymd}_{hh:02d}_v*.fits")
    orb = latest(f"{ARCHIVE}/{day}/orb/svom_orb_{ymd}_{hh:02d}_v*.fits")
    if evt is None:
        return None
    dets = []
    with fits.open(evt) as h:
        gti = [(float(a), float(b)) for a, b in zip(h["GTI"].data["START"], h["GTI"].data["STOP"])]
        for i in (3, 4, 5):
            d = h[i].data
            dets.append(dict(
                t=np.asarray(d["TIME"], float), pi=np.asarray(d["PI"]),
                et=np.asarray(d["EVT_TYPE"]), gt=np.asarray(d["GAIN_TYPE"]),
                ac=np.asarray(d["ANTI_COIN"]), fl=np.asarray(d["FLAG"]),
                dt=np.asarray(d["DEAD_TIME"], float)))
    pos = None
    if orb is not None:
        with fits.open(orb) as h:
            o = h["ORB"].data
            pos = (np.asarray(o["TIME"], float), np.asarray(o["LAT"], float), np.asarray(o["LON"], float))
    return dets, gti, pos


def main():
    day, step, prefix = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    hours = list(range(0, 24, step))

    joint = {}                      # (ac, et, gt, fl) -> count
    pi_by_ac = {0: np.zeros(260), 1: np.zeros(260)}
    dt_by_ac = {0: [], 1: []}
    mag_edges = np.arange(-60, 61, 5.0)
    mag_n = np.zeros(len(mag_edges) - 1); mag_ac = np.zeros_like(mag_n)
    saa_n = np.zeros(2); saa_ac = np.zeros(2)
    rate_bins = np.zeros(0)
    corr = {0: [], 1: []}           # AC 值 -> 到另一路最近事例的 |dt|（µs）
    partner_ac = {0: [0, 0], 1: [0, 0]}   # [伙伴 AC=0 数, 伙伴 AC=1 数]，按本事例 AC 分
    per_det_n = np.zeros(3); per_det_ac = np.zeros(3)
    rate_vs_ac = []                 # (每秒计数, 该秒 AC 率)
    n_hours = 0

    for hh in hours:
        got = read_hour(day, hh)
        if got is None:
            print("  缺文件", day, hh, flush=True); continue
        dets, gti, pos = got
        n_hours += 1
        # 只用 GTI 内、EVT_TYPE=0（正常事例）；PI 阈值留到分箱时再看，先要全谱
        keep = []
        for k, d in enumerate(dets):
            m = np.zeros(len(d["t"]), bool)
            for a, b in gti:
                m |= (d["t"] >= a) & (d["t"] <= b)
            keep.append(m)
            # 联合分布用全部事例（含 EVT_TYPE=1）
            for ac, et, gt, fl in zip(d["ac"][m], d["et"][m], d["gt"][m], d["fl"][m]):
                joint[(int(ac), int(et), int(gt), int(fl))] = joint.get((int(ac), int(et), int(gt), int(fl)), 0) + 1
            good = m & (d["et"] == 0)
            per_det_n[k] += good.sum(); per_det_ac[k] += d["ac"][good].sum()
            for v in (0, 1):
                sel = good & (d["ac"] == v)
                pi_by_ac[v] += np.bincount(np.clip(d["pi"][sel], 0, 259), minlength=260)[:260]
                dt_by_ac[v].append(d["dt"][sel][::37])       # 抽样，够看分布
        # 空间依赖：用第一路（统计足够）按秒插值位置
        if pos is not None:
            ot, olat, olon = pos
            d = dets[0]; good = keep[0] & (d["et"] == 0)
            t = d["t"][good]; ac = d["ac"][good]
            if len(t) and len(ot) > 2:
                lat = np.interp(t, ot, olat); lon = np.interp(t, ot, olon)
                ml = dipole_lat(lat, lon); saa = in_saa(lat, lon)
                idx = np.clip(np.digitize(ml, mag_edges) - 1, 0, len(mag_n) - 1)
                np.add.at(mag_n, idx, 1); np.add.at(mag_ac, idx, ac)
                saa_n[0] += (~saa).sum(); saa_ac[0] += ac[~saa].sum()
                saa_n[1] += saa.sum(); saa_ac[1] += ac[saa].sum()
                # 每秒计数率与该秒 AC 率
                sec = np.floor(t).astype(np.int64); u, inv = np.unique(sec, return_inverse=True)
                cnt = np.bincount(inv); acs = np.bincount(inv, weights=ac)
                rate_vs_ac.extend(zip(cnt.tolist(), (acs / np.maximum(cnt, 1)).tolist()))
        # 时间相关：第一路的事例，到第二/三路最近事例的时间差
        d0 = dets[0]; g0 = keep[0] & (d0["et"] == 0)
        t0 = d0["t"][g0]; a0 = d0["ac"][g0]
        others = []
        for k in (1, 2):
            dk = dets[k]; gk = keep[k] & (dk["et"] == 0)
            others.append((np.sort(dk["t"][gk]), dk["ac"][gk][np.argsort(dk["t"][gk])]))
        for v in (0, 1):
            sel = np.where(a0 == v)[0]
            if len(sel) > MAX_PAIRS:
                sel = np.random.default_rng(0).choice(sel, MAX_PAIRS, replace=False)
            ts = t0[sel]
            best = np.full(len(ts), np.inf); best_ac = np.zeros(len(ts), int)
            for ot_, oac in others:
                if len(ot_) < 2:
                    continue
                j = np.searchsorted(ot_, ts)
                for off in (-1, 0):
                    jj = np.clip(j + off, 0, len(ot_) - 1)
                    dd = np.abs(ot_[jj] - ts)
                    upd = dd < best
                    best[upd] = dd[upd]; best_ac[upd] = oac[jj][upd]
            corr[v].append(best * 1e6)
            close = best * 1e6 <= COINC_US
            partner_ac[v][0] += int((close & (best_ac == 0)).sum())
            partner_ac[v][1] += int((close & (best_ac == 1)).sum())
        print("  %s %02dz 完成" % (day, hh), flush=True)

    # 输出
    with open(prefix + "_summary.txt", "w") as f:
        def out(s):
            print(s, flush=True); f.write(s + "\n")
        out("小时数 %d（%s，每 %d 小时取一个）" % (n_hours, day, step))
        tot = sum(joint.values())
        out("\n[1] 联合分布 (ANTI_COIN, EVT_TYPE, GAIN_TYPE, FLAG) -> 计数 / 占比")
        for k in sorted(joint, key=lambda k: -joint[k]):
            out("   AC=%d ET=%d GAIN=%d FLAG=%d : %10d  %6.3f%%" % (k + (joint[k], 100 * joint[k] / tot)))
        out("\n   逐路触发率: " + " ".join("GRD%d %.4f%%" % (i + 1, 100 * per_det_ac[i] / max(per_det_n[i], 1)) for i in range(3)))
        for v in (0, 1):
            p = pi_by_ac[v]; c = np.arange(260)
            m = p.sum()
            out("   AC=%d 的 PI 谱: n=%.0f 中位 %.0f, PI<25 占 %.1f%%, PI>=100 占 %.1f%%" % (
                v, m, c[np.searchsorted(np.cumsum(p), m / 2)], 100 * p[:25].sum() / m, 100 * p[100:].sum() / m))
            dd = np.concatenate(dt_by_ac[v]) if dt_by_ac[v] else np.array([0.0])
            out("   AC=%d 的 DEAD_TIME(µs): 中位 %.3f p10 %.3f p90 %.3f" % (v, np.median(dd), np.percentile(dd, 10), np.percentile(dd, 90)))

        out("\n[2] 空间依赖（GRD1，EVT_TYPE=0）")
        out("   磁纬     事例数      AC 率")
        for i in range(len(mag_n)):
            if mag_n[i] < 1000:
                continue
            out("   %+4.0f..%+4.0f %9.0f %8.4f%%" % (mag_edges[i], mag_edges[i + 1], mag_n[i], 100 * mag_ac[i] / mag_n[i]))
        out("   SAA 外 %9.0f %8.4f%% ；SAA 内 %9.0f %8.4f%%" % (saa_n[0], 100 * saa_ac[0] / max(saa_n[0], 1), saa_n[1], 100 * saa_ac[1] / max(saa_n[1], 1)))
        if rate_vs_ac:
            rv = np.array(rate_vs_ac)
            out("   每秒计数率 vs 该秒 AC 率（分位）:")
            q = np.quantile(rv[:, 0], [0.1, 0.3, 0.5, 0.7, 0.9, 0.99])
            for lo, hi in zip([0] + list(q), list(q) + [1e9]):
                m = (rv[:, 0] >= lo) & (rv[:, 0] < hi)
                if m.sum() < 20:
                    continue
                out("     率 %6.0f–%6.0f c/s: n=%6d AC 率 %.4f%%" % (lo, hi, m.sum(), 100 * np.average(rv[m, 1], weights=rv[m, 0])))

        out("\n[3] 三路时间相关（GRD1 的事例，到 GRD2/3 最近事例的 |dt|）")
        for v in (0, 1):
            c = np.concatenate(corr[v]) if corr[v] else np.array([1e9])
            out("   本事例 AC=%d: n=%d 最近 |dt| 中位 %.1f µs; ≤%.0f µs 的占 %.2f%%; ≤1 µs 的占 %.3f%%" % (
                v, len(c), np.median(c), COINC_US, 100 * (c <= COINC_US).mean(), 100 * (c <= 1.0).mean()))
            tot_close = sum(partner_ac[v])
            if tot_close:
                out("     其中伙伴事例 AC=1 的占 %.2f%%（%d/%d）" % (100 * partner_ac[v][1] / tot_close, partner_ac[v][1], tot_close))
    np.savez(prefix + "_arrays.npz", mag_edges=mag_edges, mag_n=mag_n, mag_ac=mag_ac,
             saa_n=saa_n, saa_ac=saa_ac, pi_ac0=pi_by_ac[0], pi_ac1=pi_by_ac[1],
             corr0=np.concatenate(corr[0]) if corr[0] else np.zeros(0),
             corr1=np.concatenate(corr[1]) if corr[1] else np.zeros(0),
             rate_vs_ac=np.array(rate_vs_ac) if rate_vs_ac else np.zeros((0, 2)))
    print("wrote", prefix + "_summary.txt", prefix + "_arrays.npz")


if __name__ == "__main__":
    main()
