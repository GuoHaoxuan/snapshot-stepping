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
    Be, Bn, Bu = ppigrf.igrf(lon, lat, h, date)
    v = np.array([float(np.squeeze(Be)), float(np.squeeze(Bn)), float(np.squeeze(Bu))])
    return v / np.linalg.norm(v)


def deriv(state, date, sign):
    lon, lat, h = state
    e, n, u = sign * unit_b(lon, lat, h, date)
    r = R_E + h
    return np.array([np.degrees(e / (r * np.cos(np.radians(lat)))), np.degrees(n / r), u])


def trace(lon, lat, h, date, sign):
    """沿 sign*B 走到 100 km；返回 (lon, lat, 路径长度 km, 是否到达)。"""
    x = np.array([lon, lat, h], float); path = 0.0
    for _ in range(MAX_STEPS):
        k1 = deriv(x, date, sign)
        k2 = deriv(x + 0.5 * STEP_KM * k1, date, sign)
        k3 = deriv(x + 0.5 * STEP_KM * k2, date, sign)
        k4 = deriv(x + STEP_KM * k3, date, sign)
        step = STEP_KM * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        if x[2] + step[2] <= FOOT_KM:
            frac = (x[2] - FOOT_KM) / max(-step[2], 1e-9)
            x = x + frac * step; path += frac * STEP_KM
            return ((x[0] + 180) % 360) - 180, x[1], path, True
        x = x + step; path += STEP_KM
        if abs(x[1]) > 89.5 or x[2] > 60000:
            break
    return ((x[0] + 180) % 360) - 180, x[1], path, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv"); ap.add_argument("tgfs", nargs="+"); ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    alt = {}
    for path in args.tgfs:
        for rec in json.load(open(path)):
            s = rec["signal"]; alt[(s["instrument"], s["start"][:23])] = s["position"]["altitude"] / 1000.0
    rows = list(csv.DictReader(open(args.csv)))
    w = csv.writer(open(args.output, "w", newline=""))
    w.writerow(["sat", "start", "fa", "dur_us", "lon", "lat", "alt_km", "cls", "footN_lon", "footN_lat", "pathN_km", "reachedN", "footS_lon", "footS_lat", "pathS_km", "reachedS", "travelN_ms", "travelS_ms"])
    for r in rows:
        key = (r["sat"], r["start"][:23])
        if key not in alt:
            print("no altitude for", key); continue
        lon, lat, h = float(r["lon"]), float(r["lat"]), alt[key]
        date = datetime.strptime(r["start"][:10], "%Y-%m-%d")
        cls = "short" if float(r["dur_us"]) < 500 else "long"
        # 北半球里 B 向下：沿 +B 到北足点，沿 -B 到南足点（IGRF 的 Bu 在北半球为负）
        fn_lon, fn_lat, pn, okn = trace(lon, lat, h, date, +1.0)
        fs_lon, fs_lat, ps, oks = trace(lon, lat, h, date, -1.0)
        if fn_lat < fs_lat:   # 保证 N 是北足点
            fn_lon, fn_lat, pn, okn, fs_lon, fs_lat, ps, oks = fs_lon, fs_lat, ps, oks, fn_lon, fn_lat, pn, okn
        w.writerow([r["sat"], r["start"][:23], r["fa"], r["dur_us"], f"{lon:.3f}", f"{lat:.3f}", f"{h:.1f}", cls,
                    f"{fn_lon:.3f}", f"{fn_lat:.3f}", f"{pn:.0f}", int(okn), f"{fs_lon:.3f}", f"{fs_lat:.3f}", f"{ps:.0f}", int(oks),
                    f"{pn / (0.95 * C_KM_S) * 1e3:.1f}", f"{ps / (0.95 * C_KM_S) * 1e3:.1f}"])
        print("%s %s %-5s (%7.2f,%6.2f) -> N (%7.2f,%6.2f) %5.0f km %s | S (%7.2f,%6.2f) %5.0f km %s" % (r["sat"], r["start"][:19], cls, lon, lat, fn_lon, fn_lat, pn, "ok" if okn else "OPEN", fs_lon, fs_lat, ps, "ok" if oks else "OPEN"))


if __name__ == "__main__":
    main()
