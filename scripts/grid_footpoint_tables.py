"""把候选搬到场线足点、把时刻提前电子飞行时间，做成 `blink wwlln` 能读的候选表。

对每个候选出两份：近足点（几百 km，飞行 1–3 ms）与远足点（上万 km，飞行几十到两百 ms）。
时刻前移后，关联窗 ±W 就等价于"闪电领先探测 travel ± W"。两份表都挂在 GRID-02 的目录名下
（它的首日最早，能覆盖四颗星的全部日期），候选的 instrument 字段保留原星名。

用法:
    python3 scripts/grid_footpoint_tables.py <footpoints.csv> <tgfs_*.json ...> -o <DIR>
输出 <DIR>/near/data/GRID-02/... 与 <DIR>/far/data/GRID-02/...，以及 <DIR>/map.csv。
"""
import argparse, csv, json, os
from collections import defaultdict
from datetime import datetime, timedelta, timezone


def parse(iso):
    b = iso.rstrip("Z"); h, f = b.split("."); f = (f + "000000000")[:9]
    return datetime.strptime(h, "%Y-%m-%dT%H:%M:%S").replace(tzinfo=timezone.utc), int(f)


def fmt(dt, ns):
    return dt.strftime("%Y-%m-%dT%H:%M:%S") + ".%09dZ" % ns


def shift(iso, ms):
    dt, ns = parse(iso)
    total = dt.timestamp() * 1e9 + ns - ms * 1e6
    sec = int(total // 1e9); ns2 = int(total - sec * 1e9)
    return fmt(datetime.fromtimestamp(sec, tz=timezone.utc), ns2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("footpoints"); ap.add_argument("tgfs", nargs="+"); ap.add_argument("-o", "--outdir", required=True)
    args = ap.parse_args()
    sig = {}
    for path in args.tgfs:
        for rec in json.load(open(path)):
            s = rec["signal"]; sig[(s["instrument"], s["start"][:23])] = s
    rows = list(csv.DictReader(open(args.footpoints)))
    tables = {"near": defaultdict(list), "far": defaultdict(list)}
    mapping = []
    for r in rows:
        s = sig[(r["sat"], r["start"][:23])]
        pn, ps = float(r["pathN_km"]), float(r["pathS_km"])
        near, far = (("N", r["footN_lon"], r["footN_lat"], float(r["travelN_ms"])), ("S", r["footS_lon"], r["footS_lat"], float(r["travelS_ms"]))) if pn < ps else (("S", r["footS_lon"], r["footS_lat"], float(r["travelS_ms"])), ("N", r["footN_lon"], r["footN_lat"], float(r["travelN_ms"])))
        for tag, (hemi, lon, lat, travel) in (("near", near), ("far", far)):
            c = dict(s)
            c["start"] = shift(s["start"], travel); c["stop"] = shift(s["stop"], travel)
            c["position"] = {"longitude": float(lon), "latitude": float(lat), "altitude": 100000.0}
            day = c["start"][:10]
            tables[tag][day].append(c)
            mapping.append([r["sat"], r["start"][:23], r["cls"], tag, hemi, lon, lat, f"{travel:.1f}", c["start"][:23]])
    for tag, days in tables.items():
        for day, cands in days.items():
            y, m, d = day.split("-")
            outdir = os.path.join(args.outdir, tag, "data", "GRID-02", y, m); os.makedirs(outdir, exist_ok=True)
            json.dump(sorted(cands, key=lambda c: c["start"]), open(os.path.join(outdir, f"{y}{m}{d}_signals.json"), "w"), indent=1)
    with open(os.path.join(args.outdir, "map.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["sat", "start", "cls", "table", "hemi", "lon", "lat", "travel_ms", "shifted_start"]); w.writerows(mapping)
    print("candidates", len(rows), "-> near/far tables in", args.outdir)


if __name__ == "__main__":
    main()
