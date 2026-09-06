"""把候选所在处的地磁场线追到南北两个 100 km 高度的足点（IGRF，RK4），并估算电子沿场线的飞行时间。

用途：天格的毫秒级软暴若是 TGF 的电子束（TEB），闪电应在场线足点附近、领先几十毫秒，
而不在星下点。输出 CSV 供造"足点版"候选表去跑 `blink wwlln --window-ms`。

用法:
    python3 scripts/grid_footpoints.py <features.csv> <tgfs_*.json ...> -o <CSV>
"""
import argparse, csv, json
import numpy as np
import ppigrf
from datetime import datetime, timezone

R_E = 6371.2
FOOT_KM = 100.0
STEP_KM = 5.0
MAX_STEPS = 12000   # 60000 km：再长就是开放场线
C_KM_S = 299792.458


def unit_b(lon, lat, h, date):
    """一批点上的单位磁场向量 (E, N, U)。IGRF 逐点算勒让德函数很慢，所以整批一起算。"""
    Be, Bn, Bu = ppigrf.igrf(lon, lat, h, date)
    v = np.stack([np.ravel(Be), np.ravel(Bn), np.ravel(Bu)], axis=1).astype(float)
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def deriv(state, date, sign):
    lon, lat, h = state[:, 0], state[:, 1], state[:, 2]
    b = sign * unit_b(lon, lat, h, date)
    r = R_E + h
    return np.stack([np.degrees(b[:, 0] / (r * np.cos(np.radians(lat)))), np.degrees(b[:, 1] / r), b[:, 2]], axis=1)


def trace(lon, lat, h, date, sign):
    """整批候选沿 sign*B 同步走到 100 km；返回 (lon, lat, 路径长度 km, 是否到达) 的数组。"""
    x = np.stack([lon, lat, h], axis=1).astype(float); n = len(x)
    path = np.zeros(n); done = np.zeros(n, bool); reached = np.zeros(n, bool)
    for _ in range(MAX_STEPS):
        act = ~done
        if not act.any():
            break
        xa = x[act]
        k1 = deriv(xa, date, sign)
        k2 = deriv(xa + 0.5 * STEP_KM * k1, date, sign)
        k3 = deriv(xa + 0.5 * STEP_KM * k2, date, sign)
        k4 = deriv(xa + STEP_KM * k3, date, sign)
        step = STEP_KM * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        idx = np.where(act)[0]
        for j, i in enumerate(idx):
            if xa[j, 2] + step[j, 2] <= FOOT_KM:
                frac = (xa[j, 2] - FOOT_KM) / max(-step[j, 2], 1e-9)
                x[i] = xa[j] + frac * step[j]; path[i] += frac * STEP_KM; done[i] = True; reached[i] = True
            else:
                x[i] = xa[j] + step[j]; path[i] += STEP_KM
                if abs(x[i, 1]) > 89.5 or x[i, 2] > 60000:
                    done[i] = True
    return ((x[:, 0] + 180) % 360) - 180, x[:, 1], path, reached


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv"); ap.add_argument("tgfs", nargs="+"); ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    alt = {}
    for path in args.tgfs:
        for rec in json.load(open(path)):
            s = rec["signal"]; alt[(s["instrument"], s["start"][:23])] = s["position"]["altitude"] / 1000.0
    rows = [r for r in csv.DictReader(open(args.csv)) if (r["sat"], r["start"][:23]) in alt]
    lon = np.array([float(r["lon"]) for r in rows]); lat = np.array([float(r["lat"]) for r in rows]); h = np.array([alt[(r["sat"], r["start"][:23])] for r in rows])
    # IGRF 的历元按候选的中位年份取一次：一两年内场线足点的变化远小于 800 km 的关联半径
    years = sorted(int(r["start"][:4]) for r in rows); date = datetime(years[len(years) // 2], 7, 1)
    fn_lon, fn_lat, pn, okn = trace(lon, lat, h, date, +1.0)
    fs_lon, fs_lat, ps, oks = trace(lon, lat, h, date, -1.0)
    swap = fn_lat < fs_lat   # 保证 N 是北足点
    for a, b in ((fn_lon, fs_lon), (fn_lat, fs_lat), (pn, ps), (okn, oks)):
        a[swap], b[swap] = b[swap].copy(), a[swap].copy()
    w = csv.writer(open(args.output, "w", newline=""))
    w.writerow(["sat", "start", "fa", "dur_us", "lon", "lat", "alt_km", "cls", "footN_lon", "footN_lat", "pathN_km", "reachedN", "footS_lon", "footS_lat", "pathS_km", "reachedS", "travelN_ms", "travelS_ms"])
    for i, r in enumerate(rows):
        cls = "short" if float(r["dur_us"]) < 500 else "long"
        w.writerow([r["sat"], r["start"][:23], r["fa"], r["dur_us"], f"{lon[i]:.3f}", f"{lat[i]:.3f}", f"{h[i]:.1f}", cls,
                    f"{fn_lon[i]:.3f}", f"{fn_lat[i]:.3f}", f"{pn[i]:.0f}", int(okn[i]), f"{fs_lon[i]:.3f}", f"{fs_lat[i]:.3f}", f"{ps[i]:.0f}", int(oks[i]),
                    f"{pn[i] / (0.95 * C_KM_S) * 1e3:.1f}", f"{ps[i] / (0.95 * C_KM_S) * 1e3:.1f}"])
        print("%s %s %-5s (%7.2f,%6.2f) -> N (%7.2f,%6.2f) %5.0f km %s | S (%7.2f,%6.2f) %5.0f km %s" % (
            r["sat"], r["start"][:19], cls, lon[i], lat[i], fn_lon[i], fn_lat[i], pn[i], "ok" if okn[i] else "OPEN", fs_lon[i], fs_lat[i], ps[i], "ok" if oks[i] else "OPEN"))


if __name__ == "__main__":
    main()
