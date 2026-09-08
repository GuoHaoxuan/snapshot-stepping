"""天格 GRID TGF 搜索的讲图：给报告用，一张图讲一件事，字号按投影调大。

图 1：显著候选分成两群，短硬暴关联闪电、毫秒级软暴不关联。
图 2：地理分布 + 磁力线足点检验（软暴不是 TGF 的电子束）。

用法:
    python3 scripts/plot_grid_talk.py <features_sig.csv> <tgfs_*.json ...> -o <目录> [--footpoint <dir>]
"""
import argparse, csv, glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

plt.rcParams.update({
    "font.sans-serif": ["PingFang SC", "Arial Unicode MS"], "font.family": "sans-serif",
    "axes.unicode_minus": False, "font.size": 15, "axes.titlesize": 17, "axes.labelsize": 15,
    "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13, "lines.linewidth": 2,
})
SHORT_US = 500.0
RED, BLUE, GREY = "#c53030", "#2b6cb0", "0.55"
POLE_LAT, POLE_LON = np.radians(80.7), np.radians(-72.7)


def dipole_lat(lat_deg, lon_deg):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    s = np.sin(lat) * np.sin(POLE_LAT) + np.cos(lat) * np.cos(POLE_LAT) * np.cos(lon - POLE_LON)
    return np.degrees(np.arcsin(np.clip(s, -1, 1)))


def load(features, tgfs):
    assoc = {}
    for path in tgfs:
        for rec in json.load(open(path)):
            s = rec["signal"]
            assoc[(s["instrument"], s["start"][:23])] = bool(rec["lightning"].get("associated"))
    rows = list(csv.DictReader(open(features)))
    d = dict(
        sat=np.array([r["sat"] for r in rows]),
        start=np.array([r["start"][:23] for r in rows]),
        fa=np.array([float(r["fa"]) for r in rows]),
        dur=np.array([float(r["dur_us"]) for r in rows]),
        lon=np.array([float(r["lon"]) for r in rows]),
        lat=np.array([float(r["lat"]) for r in rows]),
        pir=np.array([float(r["pi_ratio"]) if r["pi_ratio"] not in ("", "nan") else np.nan for r in rows]),
        rate=np.array([float(r["rate_win"]) for r in rows]),
        n=np.array([float(r["n_core"]) for r in rows]),
    )
    d["assoc"] = np.array([assoc.get((s, t), False) for s, t in zip(d["sat"], d["start"])])
    d["mlat"] = dipole_lat(d["lat"], d["lon"])
    d["short"] = d["dur"] < SHORT_US
    return d


def fig_two_populations(d, out):
    """图 1：三个维度上两群分开，闪电只认其中一群。"""
    s, l = d["short"], ~d["short"]
    a = d["assoc"]
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6))
    ax = axes[0]
    bins = np.logspace(np.log10(20), np.log10(1100), 22)
    ax.hist(d["dur"][l], bins=bins, color=BLUE, alpha=0.55, label="毫秒级软暴 (%d)" % l.sum())
    ax.hist(d["dur"][s], bins=bins, color=RED, alpha=0.65, label="短硬暴 (%d)" % s.sum())
    ax.hist(d["dur"][a], bins=bins, histtype="step", edgecolor="k", lw=2.4, label="闪电证实 (%d)" % a.sum())
    ax.axvline(SHORT_US, color="k", ls="--", lw=1.5)
    ax.set_xscale("log"); ax.set_xlabel("搜索窗长 (µs)"); ax.set_ylabel("候选数")
    # 软暴堆在 1 ms 是搜索窗上限，不是真实时长：宽窗量出来是 3–4 ms
    ax.annotate("顶在 1 ms 搜索上限\n宽窗量出真实时长 3–4 ms", xy=(950, 43), xytext=(150, 34),
                fontsize=12, color="#2c5282", ha="center",
                arrowprops=dict(arrowstyle="->", color="#2c5282", lw=1.4))
    ax.set_title("(a) 时长", pad=10); ax.legend(loc="upper left", fontsize=12)
    ax = axes[1]
    bins = np.linspace(0, 8, 21)
    ax.hist(np.clip(d["pir"][l], 0, 8), bins=bins, color=BLUE, alpha=0.55)
    ax.hist(np.clip(d["pir"][s], 0, 8), bins=bins, color=RED, alpha=0.65)
    ax.hist(np.clip(d["pir"][a], 0, 8), bins=bins, histtype="step", edgecolor="k", lw=2.4)
    ax.axvline(1, color="k", ls=":", lw=1.5)
    ax.set_xlabel("窗内能道中位数 ÷ 本底"); ax.set_ylabel("候选数")
    ax.set_title("(b) 谱硬度：短暴比本底硬得多", pad=10)
    ax = axes[2]
    bins = np.linspace(0, 70, 15)
    ax.hist(np.abs(d["mlat"][l]), bins=bins, color=BLUE, alpha=0.55)
    ax.hist(np.abs(d["mlat"][s]), bins=bins, color=RED, alpha=0.65)
    ax.hist(np.abs(d["mlat"][a]), bins=bins, histtype="step", edgecolor="k", lw=2.4)
    ax.set_xlabel("|偶极磁纬| (°)"); ax.set_ylabel("候选数")
    ax.set_title("(c) 磁纬：软暴在高磁纬", pad=10)
    fig.suptitle("天格 %d 个显著候选分成两群：短硬暴关联闪电（%d/%d），毫秒级软暴一个都不关联（0/%d）"
                 % (len(d["fa"]), int((a & s).sum()), int(s.sum()), int(l.sum())), fontsize=17)
    fig.tight_layout(rect=(0, 0, 1, 0.92)); fig.savefig(out, dpi=160); print("wrote", out)


def fig_map_footpoint(d, out, fp_dir=None):
    """图 2：地理分布，以及足点检验（软暴不是 TGF 的电子束）。"""
    s, l, a = d["short"], ~d["short"], d["assoc"]
    has_fp = fp_dir and os.path.exists(os.path.join(fp_dir, "map.csv"))
    fig = plt.figure(figsize=(14, 8.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1], hspace=0.42)
    ax = fig.add_subplot(gs[0], projection=ccrs.PlateCarree())
    span = min(90, np.ceil(np.abs(d["lat"]).max()) + 4)
    ax.set_extent([-180, 180, -span, span], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="0.94")
    ax.add_feature(cfeature.COASTLINE, lw=0.5, edgecolor="0.45")
    ax.gridlines(draw_labels=False, lw=0.3, color="0.9")
    ax.scatter(d["lon"][l], d["lat"][l], s=52, c=BLUE, marker="s", lw=0.4, edgecolor="k", alpha=0.8,
               transform=ccrs.PlateCarree(), zorder=4, label="毫秒级软暴 (%d)" % l.sum())
    ax.scatter(d["lon"][s & ~a], d["lat"][s & ~a], s=48, c="0.75", lw=0.4, edgecolor="0.4",
               transform=ccrs.PlateCarree(), zorder=5, label="短硬暴，未关联 (%d)" % int((s & ~a).sum()))
    ax.scatter(d["lon"][a], d["lat"][a], s=150, c=RED, marker="*", lw=0.6, edgecolor="k",
               transform=ccrs.PlateCarree(), zorder=6, label="闪电证实的 TGF (%d)" % a.sum())
    ax.set_title("短硬暴在低纬雷暴区，软暴在高纬", pad=10)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.30), ncol=3, frameon=False, fontsize=12)

    ax = fig.add_subplot(gs[1])
    if has_fp:
        m = {(r["table"], r["shifted_start"]): r for r in csv.DictReader(open(os.path.join(fp_dir, "map.csv")))}
        res = {}
        for tag, w in (("near", 10), ("far", 25)):
            path = os.path.join(fp_dir, "tgfs_%s.json" % tag)
            if not os.path.exists(path): continue
            for rec in json.load(open(path)):
                sig, li = rec["signal"], rec["lightning"]
                r = m.get((tag, sig["start"][:23]))
                if r is None: continue
                key = (tag, r["cls"])
                n, k = res.get(key, (0, 0))
                res[key] = (n + 1, k + int(bool(li.get("associated"))))
        labels, vals, colors = [], [], []
        for cls, cname, color in (("short", "短硬暴", RED), ("long", "毫秒级软暴", BLUE)):
            for tag, tname in (("near", "近足点"), ("far", "远足点")):
                n, k = res.get((tag, cls), (0, 0))
                labels.append("%s\n%s" % (cname, tname)); vals.append(k); colors.append(color)
        x = np.arange(len(labels))
        bars = ax.bar(x, vals, 0.55, color=colors, alpha=0.85)
        for xi, v, (cls, tag) in zip(x, vals, [(c, t) for c in ("short", "long") for t in ("near", "far")]):
            n = res.get((tag, cls), (0, 0))[0]
            ax.text(xi, v + 0.15, "%d / %d" % (v, n), ha="center", fontsize=14)
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=13)
        ax.set_ylabel("关联到闪电的个数"); ax.set_ylim(0, max(max(vals) + 1.5, 3))
        ax.set_title("足点检验：把候选沿磁力线追到 100 km 足点、时刻前移电子飞行时间后再查闪电", pad=10)
        ax.text(0.99, 0.93, "软暴若是 TGF 的电子束，按短暴的关联率应有 5–6 个\n两个足点都是 0（P ≈ 0.3%）→ 不是 TEB",
                transform=ax.transAxes, ha="right", va="top", fontsize=13,
                bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    else:
        ax.axis("off"); ax.text(0.5, 0.5, "足点检验结果未就绪", ha="center", va="center")
    fig.suptitle("天格：短硬暴是 TGF，毫秒级软暴是磁层电子沉降", fontsize=18)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("features"); ap.add_argument("tgfs", nargs="+")
    ap.add_argument("-o", "--outdir", required=True); ap.add_argument("--footpoint")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    d = load(args.features, args.tgfs)
    fig_two_populations(d, os.path.join(args.outdir, "grid_talk_1_populations.png"))
    fig_map_footpoint(d, os.path.join(args.outdir, "grid_talk_2_map.png"), args.footpoint)


if __name__ == "__main__":
    main()
