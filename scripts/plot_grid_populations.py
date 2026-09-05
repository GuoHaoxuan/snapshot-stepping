"""天格 GRID 全量搜索结果里两个群体的分离：短硬暴（TGF 候选）与毫秒级软暴。

给了 `blink wwlln` 的 tgfs JSON 就把经 WWLLN 闪电证实的那一批单独标出来——它是这张
图里唯一有外部真值的子集。两群的分界（500 µs）来自宽窗轮廓：短群在 ±30 ms 里是单个
亚毫秒脉冲、谱比本底硬；长群真实时长 3–4 ms、谱与本底一样、都在高磁纬。

用法:
    python3 scripts/plot_grid_populations.py <features.csv> [tgfs_*.json ...] -o <PNG>
"""
import argparse, csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

plt.rcParams.update({
    "font.sans-serif": ["PingFang SC", "Arial Unicode MS"],
    "font.family": "sans-serif",
    "axes.unicode_minus": False,
})

SHORT_US = 500.0   # 搜索窗长的分界；长群大多顶在 1 ms 的搜索上限


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("tgfs", nargs="*", help="blink wwlln 写出的 tgfs JSON（每颗星一个）")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.csv)))
    assoc = {}
    for path in args.tgfs:
        for rec in json.load(open(path)):
            s, li = rec["signal"], rec["lightning"]
            assoc[(s["instrument"], s["start"][:23])] = bool(li.get("associated")) and bool(li.get("in_coverage", True))
    for r in rows:
        r["confirmed"] = int(assoc.get((r["sat"], r["start"][:23]), False))
    col = lambda k, t=float: np.array([t(r[k]) for r in rows])
    confirmed = col("confirmed", int).astype(bool)
    sat = np.array([r["sat"] for r in rows])
    lon, lat = col("lon"), col("lat")
    dur, rate, ratio = col("dur_us"), col("rate_win"), col("pi_ratio")
    fa = np.log10(col("fa"))

    short = dur < SHORT_US
    long_ = ~short
    classes = [
        (short, "短硬暴 (TGF 候选)", "tab:blue", "o"),
        (long_, "毫秒级软暴", "tab:red", "s"),
    ]

    fig = plt.figure(figsize=(15, 7.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.85, 1], hspace=0.72, wspace=0.26)

    ax = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
    # 天格是太阳同步轨道，什么纬度都到得了；按候选实际落点裁一下
    span = min(90, np.ceil(np.abs(lat).max()) + 4)
    ax.set_extent([-180, 180, -span, span], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="0.95")
    ax.add_feature(cfeature.COASTLINE, lw=0.4, edgecolor="0.4")
    ax.gridlines(draw_labels=False, lw=0.3, color="0.88")
    pcm = ax.scatter(lon[short], lat[short], c=np.clip(fa[short], -45, None), s=22,
                     cmap="viridis_r", vmin=-45, vmax=-5, lw=0.2, edgecolor="0.3",
                     transform=ccrs.PlateCarree(),
                     label="短硬暴 (%d)" % short.sum(), zorder=4)
    ax.scatter(lon[long_], lat[long_], s=26, c="tab:red", marker="s", lw=0.3,
               edgecolor="k", transform=ccrs.PlateCarree(),
               label="毫秒级软暴 (%d)" % long_.sum(), zorder=5)
    fig.canvas.draw()
    box = ax.get_position()
    cax = fig.add_axes([box.x0 + 0.15 * box.width, box.y0 - 0.085,
                        0.70 * box.width, 0.018])
    cb = fig.colorbar(pcm, cax=cax, orientation="horizontal")
    cb.set_label("短硬暴的 $\\log_{10}$ fa（越小越显著）", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    if confirmed.any():
        ax.scatter(lon[confirmed], lat[confirmed], s=110, facecolor="none",
                   edgecolor="crimson", lw=1.2, transform=ccrs.PlateCarree(),
                   label="其中经闪电证实 (%d)" % confirmed.sum(), zorder=6)
    ax.set_title("(a) fa$\\leq10^{-5}$ 的 %d 个候选的地理分布（%s）" %
                 (len(rows), "、".join(sorted(set(sat)))), fontsize=12)
    ax.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.04),
              ncol=3, frameon=False)

    ax = fig.add_subplot(gs[1, 0])
    bins = np.logspace(1, 3.1, 36)
    for mask, name, color, _ in classes:
        ax.hist(dur[mask], bins=bins, histtype="step", lw=1.6, color=color, label=name)
    if confirmed.any():
        ax.hist(dur[confirmed], bins=bins, color="crimson", alpha=0.55,
                label="经闪电证实 (%d)" % confirmed.sum())
    ax.axvline(1000, color="0.5", ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("搜索窗长 (µs)；1 ms 是搜索上限")
    ax.set_ylabel("候选数")
    ax.set_title("(b) 窗长：长群顶在上限，真实时长 3–4 ms", fontsize=10)
    ax.legend(fontsize=8)

    ax = fig.add_subplot(gs[1, 1])
    bins = np.linspace(0, 8, 33)
    for mask, name, color, _ in classes:
        ax.hist(np.clip(ratio[mask], 0, 8), bins=bins, histtype="step", lw=1.6, color=color, label=name)
    if confirmed.any():
        ax.hist(np.clip(ratio[confirmed], 0, 8), bins=bins, color="crimson", alpha=0.55,
                label="经闪电证实 (%d)" % confirmed.sum())
    ax.axvline(1, color="0.5", ls=":", lw=1)
    ax.set_xlabel("能道中位数比（窗内 / 本底）")
    ax.set_ylabel("候选数")
    ax.set_title("(c) 短群比本底硬，长群与本底一样", fontsize=10)
    ax.legend(fontsize=8)

    ax = fig.add_subplot(gs[1, 2])
    bins = np.linspace(0, 60, 25)
    for mask, name, color, _ in classes:
        ax.hist(np.abs(lat[mask]), bins=bins, histtype="step", lw=1.6, color=color, label=name)
    if confirmed.any():
        ax.hist(np.abs(lat[confirmed]), bins=bins, color="crimson", alpha=0.55,
                label="经闪电证实 (%d)" % confirmed.sum())
    ax.set_xlabel("|纬度| (°)")
    ax.set_ylabel("候选数")
    ax.set_title("(d) 短群在低纬雷暴区，长群在高磁纬", fontsize=10)
    ax.legend(fontsize=8)

    fig.suptitle("天格 GRID 全量搜索：显著候选里的两个群体\n"
                 "（分界来自 ±30 ms 宽窗轮廓与计数守恒检验；短群另经 WWLLN 闪电关联证实，长群 0 关联）",
                 fontsize=12)
    fig.savefig(args.output, dpi=140, bbox_inches="tight")
    print("wrote", args.output)
    for mask, name, _, _ in classes:
        print("  %-18s %3d  证实 %2d  能道中位数比 %.2f  本底率中位 %6.0f c/s  |lat| 中位 %4.1f" %
              (name, mask.sum(), (mask & confirmed).sum(), np.median(ratio[mask]),
               np.median(rate[mask]), np.median(np.abs(lat[mask]))))


if __name__ == "__main__":
    main()
