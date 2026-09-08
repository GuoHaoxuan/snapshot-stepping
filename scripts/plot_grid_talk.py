"""天格 GRID TGF 搜索的讲图：给报告用，一张图讲一件事，字号按投影调大。

图 1：显著候选分成两群，短硬暴关联闪电、毫秒级软暴不关联。
图 2：地理分布 + 磁力线足点检验（软暴不是 TGF 的电子束）。
图 3：两群各举一例的光变与逐事例沉积能量。

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
    ax.annotate("顶在 1 ms 搜索上限\n宽窗量出真实时长 3–4 ms", xy=(950, 42), xytext=(120, 24),
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
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1], hspace=0.62)
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
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.24), ncol=3, frameon=False, fontsize=12)

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
        ax.set_ylabel("关联到闪电的个数"); ax.set_ylim(0, max(max(vals) * 2.1, 4))
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        ax.set_title("足点检验：沿磁力线追到 100 km 足点、时刻前移电子飞行时间后再查闪电", pad=10, fontsize=15)
        # 零结果的说服力按图上实测值现算：拿短硬暴近足点的关联率当先验，看软暴全零有多不像话
        k_s, n_s = res.get(("near", "short"), (0, 0))[1], res.get(("near", "short"), (0, 0))[0]
        n_l = res.get(("near", "long"), (0, 0))[0]
        if n_s and n_l:
            rho = k_s / n_s
            # 只写期望个数，不写 P：先验取自低纬雷暴区，搬到高纬远洋会高估 WWLLN 的效率
            note = ("软暴若是 TGF 的电子束，按短硬暴近足点的关联率 %.0f%% 算，%d 个应有约 %.0f 个\n"
                    "两个足点都是 0；软暴的谱又与本底重合 → 不是 TEB" % (rho * 100, n_l, rho * n_l))
            ax.text(0.99, 0.93, note, transform=ax.transAxes, ha="right", va="top", fontsize=13,
                    bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    else:
        ax.axis("off"); ax.text(0.5, 0.5, "足点检验结果未就绪", ha="center", va="center")
    fig.suptitle("天格：短硬暴是 TGF，毫秒级软暴是磁层电子沉降", fontsize=18)
    # 地图受投影比例限制、比默认栅格窄，把下面的柱状图收到同样宽度免得头重脚轻
    fig.canvas.draw()
    box_map = fig.axes[0].get_position()
    box_bar = ax.get_position()
    ax.set_position([box_map.x0, box_bar.y0, box_map.width, box_bar.height])
    fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)


# 图 3 的两个例子：短硬暴取闪电证实里最显著的一个，毫秒级软暴取参数最贴近该
# 群中位（硬度比 1.2、|磁纬| 44°、本底 1.6 kc/s）的一个。
EXAMPLES = [
    dict(tag="GRID-03B_20221004T000956543", sat="GRID-03B", color=RED,
         title="短硬暴：GRID-03B 2022-10-04 00:09:56.5 UTC",
         note="152.9°E 11.96°S｜偶极磁纬 −20°｜同时刻有闪电，巧合概率 3×10⁻⁴"),
    dict(tag="GRID-02_20210124T195431291", sat="GRID-02", color=BLUE,
         title="毫秒级软暴：GRID-02 2021-01-24 19:54:31.3 UTC",
         note="149.3°W 47.2°S｜偶极磁纬 −44°｜WWLLN 有覆盖但无闪电"),
]


def _load_example(lc_dir, ex):
    """读逐事例 CSV 与 EBOUNDS，返回相对暴中心的时间、沉积能量、本底率、T90。"""
    from scipy.stats import poisson

    d = np.genfromtxt(os.path.join(lc_dir, ex["tag"] + ".csv"), delimiter=",", names=True)
    t = np.atleast_1d(d["dt_ms"]).astype(float)
    pi = np.atleast_1d(d["pi"]).astype(int)
    b = np.genfromtxt(os.path.join(lc_dir, ex["tag"] + "_bkg.csv"), delimiter=",", names=True)
    rate = float(np.atleast_1d(b["n_bkg"])[0]) / float(np.atleast_1d(b["live_s"])[0])

    eb = np.genfromtxt(os.path.join(lc_dir, "ebounds_" + ex["sat"] + ".csv"), delimiter=",", names=True)
    table = dict(zip(np.atleast_1d(eb["ch"]).astype(int),
                     np.sqrt(np.atleast_1d(eb["e_min"]) * np.atleast_1d(eb["e_max"]))))
    energy = np.array([table.get(int(c), np.nan) for c in pi])

    # 先用 1 ms 格粗框出暴的范围：含最大格、逐格泊松 p < 1e-4 的连续段，两边各放宽一格。
    # （更细的格里单个低计数格会把连续段打断，不能直接拿来量时长。）
    edges = np.arange(-15.0, 15.01, 1.0)
    cnt, _ = np.histogram(t, bins=edges)
    sig = poisson.sf(cnt - 1, rate / 1000.0) < 1e-4
    peak = int(np.argmax(cnt))
    lo = hi = peak
    while lo - 1 >= 0 and sig[lo - 1]:
        lo -= 1
    while hi + 1 < len(sig) and sig[hi + 1]:
        hi += 1
    b0, b1 = edges[max(lo - 1, 0)], edges[min(hi + 2, len(edges) - 1)]

    # 框内扣掉本底的累积计数，取 5%–95% 作为 T90；暴中心取 T90 区间中点。
    ts = np.sort(t[(t >= b0) & (t < b1)])
    cum = np.arange(1, len(ts) + 1) - rate / 1000.0 * (ts - b0)
    total = cum[-1]
    t05 = float(np.interp(0.05 * total, cum, ts))
    t95 = float(np.interp(0.95 * total, cum, ts))
    centre = 0.5 * (t05 + t95)
    inside = (t >= t05) & (t <= t95)
    excess = float(inside.sum() - rate / 1000.0 * (t95 - t05))
    return t - centre, energy, rate, inside, (t95 - t05), excess


def fig_lightcurves(lc_dir, out):
    """图 3：两群各举一例，上排光变、下排逐事例沉积能量。"""
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 9.4),
                             gridspec_kw=dict(height_ratios=[1.15, 1.0], hspace=0.34, wspace=0.22))
    half, w = 8.0, 0.2   # 上排统一 ±8 ms、200 µs 格，两群的时长直接可比

    for col, ex in enumerate(EXAMPLES):
        t, energy, rate, inside, span, excess = _load_example(lc_dir, ex)
        mu = rate * w / 1000.0

        ax = axes[0, col]
        edges = np.arange(-half, half + w / 2, w)
        ax.hist(t, bins=edges, color=ex["color"], alpha=0.85, edgecolor=ex["color"])
        ax.axhline(mu, color="0.35", ls="--", lw=1.5)
        ax.set_xlim(-half, half)
        ax.set_ylim(0, max(np.histogram(t, bins=edges)[0].max() * 1.28, 3))
        ax.set_xlabel("相对暴中心时间 (ms)")
        ax.set_ylabel(f"计数 / {w * 1000:.0f} µs")
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        ax.set_title(ex["title"], pad=10)
        ax.text(0.02, 0.95, f"({'ab'[col]})", transform=ax.transAxes, va="top", fontsize=15, weight="bold")
        ax.text(0.98, 0.95,
                f"超出本底 {excess:.0f} 计数    T90 = {span * 1000:.0f} µs\n本底 {rate / 1000:.2f} 计数/ms（虚线）",
                transform=ax.transAxes, ha="right", va="top", fontsize=13, color=ex["color"],
                linespacing=1.5)

        # 短硬暴在 200 µs 格里只剩一根针，插一个 20 µs 格的放大图
        if span < 0.5:   # 短硬暴
            ins = ax.inset_axes([0.53, 0.24, 0.43, 0.44])
            fine = 0.02
            ins.hist(t, bins=np.arange(-0.3, 0.3 + fine / 2, fine), color=ex["color"], edgecolor=ex["color"])
            ins.set_xlim(-0.3, 0.3)
            ins.set_title("放大：20 µs 格", fontsize=12, pad=3)
            ins.tick_params(labelsize=11)
            ins.set_xlabel("ms", fontsize=11, labelpad=1)

        ax = axes[1, col]
        far = np.abs(t) > 1.5 * max(span, 0.3)
        ax.scatter(t[far], energy[far], s=22, color=GREY, alpha=0.55, label="本底事例", zorder=2)
        ax.scatter(t[inside], energy[inside], s=52, color=ex["color"], edgecolor="white",
                   linewidth=0.6, label="暴内事例", zorder=3)
        med_b, med_i = np.nanmedian(energy[far]), np.nanmedian(energy[inside])
        ax.axhline(med_b, color="0.35", ls="--", lw=1.5)
        ax.text(-half * 0.98, med_b, f" 本底中位 {med_b:.0f} keV", ha="left", va="bottom",
                fontsize=12, color="0.35")
        ax.axhline(med_i, color=ex["color"], ls="-", lw=1.5, alpha=0.7)
        ax.text(half * 0.98, med_i, f"暴内中位 {med_i:.0f} keV ", ha="right",
                va="bottom" if med_i > med_b else "top", fontsize=13, color=ex["color"])
        ax.set_yscale("log")
        ax.set_ylim(8, 2500)
        ax.set_xlim(-half, half)
        ax.set_xlabel("相对暴中心时间 (ms)")
        ax.set_ylabel("沉积能量 (keV)")
        ax.text(0.02, 0.95, f"({'cd'[col]})", transform=ax.transAxes, va="top", fontsize=15, weight="bold")
        ax.legend(loc="upper right", framealpha=0.9)
        ax.text(0.5, -0.30, ex["note"], transform=ax.transAxes, ha="center", fontsize=13, color="0.3")

    fig.suptitle("两群各一例：短硬暴是 70 µs 的硬脉冲，毫秒级软暴是 3.0 ms 的平顶、能谱与本底同", y=0.985)
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("features"); ap.add_argument("tgfs", nargs="+")
    ap.add_argument("-o", "--outdir", required=True); ap.add_argument("--footpoint")
    ap.add_argument("--lightcurves", help="逐事例导出目录，见 scripts/cluster/grid_lightcurve.py")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    d = load(args.features, args.tgfs)
    fig_two_populations(d, os.path.join(args.outdir, "grid_talk_1_populations.png"))
    fig_map_footpoint(d, os.path.join(args.outdir, "grid_talk_2_map.png"), args.footpoint)
    if args.lightcurves:
        fig_lightcurves(args.lightcurves, os.path.join(args.outdir, "grid_talk_3_lightcurves.png"))


if __name__ == "__main__":
    main()
