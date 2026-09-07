"""SVOM 窗长上限 1 ms（v6）与 5 ms（v7）的对照，覆盖期同口径。

判据：新增候选的闪电关联率是否与现有一致。一致说明新增的是真 TGF；掉到偶然水平
说明放宽窗长只是把本底涨落收了进来。

配对按区间重叠（放宽 10 ms）：窗长放宽后相邻触发窗合并得更多，候选的起始时刻会平移几毫秒，
按起点严格配对会把同一个暴发误判成"v6 丢了、v7 新增"。
"""
import json
import sys
from datetime import datetime, timezone

import numpy as np

TOL = 10e-3


def met(iso):
    b = iso.rstrip("Z"); h, _, f = b.partition(".")
    return datetime.strptime(h + "." + (f + "000000")[:6], "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc).timestamp()


def load(path, until="2025-01-01"):
    """只取 WWLLN 覆盖期（v7 只跑了这一段，两边必须同口径）。"""
    out = []
    for r in json.load(open(path)):
        s, li = r["signal"], r["lightning"]
        if s["start"] >= until:
            continue
        t0, t1 = met(s["start"]), met(s["stop"])
        out.append(dict(t0=t0, t1=t1, start=s["start"][:23], fa=s["false_positive_per_year"],
                        best=s["bin_size_best"] * 1e3, span=(t1 - t0) * 1e3, count=s["count"],
                        assoc=bool(li.get("associated")), cov=bool(li.get("in_coverage", True)),
                        prob=li.get("coincidence_probability") or 0.0))
    out.sort(key=lambda x: x["t0"])
    return out


def match(a, b):
    """a 的每条在 b 中是否有区间重叠（各放宽 TOL）的对应。"""
    tb0 = np.array([x["t0"] for x in b]); tb1 = np.array([x["t1"] for x in b])
    order = np.argsort(tb0); tb0, tb1 = tb0[order], tb1[order]
    bs = [b[i] for i in order]
    out = []
    for x in a:
        lo = np.searchsorted(tb0, x["t0"] - TOL - 0.05)
        hi = np.searchsorted(tb0, x["t1"] + TOL)
        hit = None
        for i in range(lo, hi):
            if tb1[i] + TOL >= x["t0"] - TOL and tb0[i] - TOL <= x["t1"] + TOL:
                hit = bs[i]; break
        out.append(hit)
    return out


def line(label, rows):
    cov = [r for r in rows if r["cov"]]
    a = [r for r in cov if r["assoc"]]
    p = sum(r["prob"] for r in cov)
    print("  %-26s 覆盖内 %5d  关联 %4d (%5.1f%%)  偶然期望 %6.2f"
          % (label, len(cov), len(a), 100 * len(a) / max(len(cov), 1), p))


def main(p6, p7):
    v6 = load(p6); v7 = load(p7)
    print("=== 总体（WWLLN 覆盖期 191 天）")
    for lab, v in (("v6 (1 ms)", v6), ("v7 (5 ms)", v7)):
        sig = [r for r in v if r["fa"] <= 1e-5]
        siga = [r for r in sig if r["cov"] and r["assoc"]]
        alla = [r for r in v if r["cov"] and r["assoc"]]
        print("  %-10s 候选 %5d  显著 %4d  显著且关联 %3d (%.1f%%)  全部关联 %3d"
              % (lab, len(v), len(sig), len(siga), 100 * len(siga) / max(len(sig), 1), len(alla)))

    m76 = match(v7, v6); m67 = match(v6, v7)
    new = [x for x, h in zip(v7, m76) if h is None]
    both7 = [x for x, h in zip(v7, m76) if h is not None]
    lost = [x for x, h in zip(v6, m67) if h is None]
    print("\n=== 逐条配对（区间重叠 ±10 ms）")
    print("  v7 新增 %d，两版都有 %d，v6 独有 %d" % (len(new), len(both7), len(lost)))

    print("\n=== 关键判据：新增候选的关联率")
    line("v7 新增（全部）", new)
    line("v7 新增（显著 fa≤1e-5）", [r for r in new if r["fa"] <= 1e-5])
    line("两版都有（显著）", [r for r in both7 if r["fa"] <= 1e-5])
    line("v6 显著（基准）", [r for r in v6 if r["fa"] <= 1e-5])
    line("v6 独有（丢掉的）", lost)

    print("\n=== v7 显著候选用了多长的窗")
    sig7 = [r for r in v7 if r["fa"] <= 1e-5]
    b = np.array([r["best"] for r in sig7])
    print("  bin_size_best ms: 中位 %.3f p90 %.3f max %.3f；>1 ms 的 %d 个 (%.1f%%)"
          % (np.median(b), np.percentile(b, 90), b.max(), (b > 1.0).sum(), 100 * (b > 1.0).mean()))
    ba = np.array([r["best"] for r in sig7 if r["cov"] and r["assoc"]])
    if len(ba):
        print("  其中已证实 %d 个: 中位 %.3f max %.3f；>1 ms 的 %d 个"
              % (len(ba), np.median(ba), ba.max(), (ba > 1.0).sum()))

    print("\n=== 已证实 TGF 在两版之间的对应")
    a6 = [r for r in v6 if r["fa"] <= 1e-5 and r["cov"] and r["assoc"]]
    hits = match(a6, v7)
    miss = [x for x, h in zip(a6, hits) if h is None]
    grew = [(x, h) for x, h in zip(a6, hits) if h is not None and h["best"] > x["best"] * 1.5]
    print("  v6 已证实 %d 个：v7 里对应不上的 %d，窗长增大 50%% 以上的 %d" % (len(a6), len(miss), len(grew)))
    for x, h in sorted(grew, key=lambda z: -z[1]["best"]):
        print("     %s  窗 %.3f → %.3f ms   fa %.1e → %.1e   计数 %d → %d"
              % (x["start"], x["best"], h["best"], x["fa"], h["fa"], x["count"], h["count"]))
    a7 = [r for r in v7 if r["fa"] <= 1e-5 and r["cov"] and r["assoc"]]
    back = match(a7, v6)
    only7 = [x for x, h in zip(a7, back) if h is None]
    print("  v7 已证实 %d 个：v6 里没有的 %d" % (len(a7), len(only7)))
    for x in only7:
        print("     %s 窗 %.3f ms fa %.1e 计数 %d" % (x["start"], x["best"], x["fa"], x["count"]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
