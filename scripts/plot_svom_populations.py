"""SVOM/GRM 全量搜索结果里三个群体的分离：1 Hz 假信号、高本底、TGF 候选。

给了关联结果就把经 WWLLN 闪电证实的那一批单独标出来——它是这张图里唯一有
外部真值的子集，其余的分类都是靠特征分布推出来的。

用法:
    python3 scripts/plot_svom_populations.py <features.csv> [assoc.csv] -o <PNG>
"""
import argparse, csv
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

RATE_HIGH = 2e4          # 本底率高于此即在 SAA 西缘，窗内能谱与本底无异
PHASE_LO, PHASE_HI = 0.500, 0.510  # 1 Hz 假信号聚集的整秒相位。窗宽 10 ms，
                                   # 实测 102 个而均匀期望 9.9 个，纯度约 90%；
                                   # 放宽到 ±30 ms 只多收 19 个而纯度掉到 75%


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("assoc", nargs="?", help="blink wwlln 的关联结果 CSV")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.csv)))
    if args.assoc:
        assoc = {r["start"][:26]: r for r in csv.DictReader(open(args.assoc))}
        for r in rows:
            a = assoc.get(r["start"][:26], {})
            r["confirmed"] = int(a.get("associated") == "1" and a.get("in_coverage") == "1")
    else:
        for r in rows:
            r["confirmed"] = 0
    col = lambda k, t=float: np.array([t(r[k]) for r in rows])
    associated = col("confirmed", int).astype(bool)
    lon, lat = col("lon"), col("lat")
    phase, rate, ratio = col("phase"), col("rate_bkg"), col("pi_med_ratio")
    dur = col("dur_ms")

    high = rate > RATE_HIGH
    onehz = (~high) & (phase >= PHASE_LO) & (phase <= PHASE_HI)
    # 这一类由排除法得到，但已有正面证据：WWLLN 覆盖内的 177 个里 67 个
    # 关联到闪电（±5 ms、800 km），偶然期望 0.74 个；同样处理的两个对照组
    # 落在期望上（高本底 0/18，1 Hz 2/23 可由相位窗的已知污染解释）。
    # 见 plot_svom_association.py。
    tgf = ~(high | onehz)
    # 只把归到 TGF 一类里的算作「已证实」。1 Hz 组里也有 2 个关联上了，但那
    # 两个是偶然落进相位窗的普通候选（见关联图），归在假信号名下会误导。
    confirmed = associated & tgf
    classes = [
        (tgf, "TGF 候选", "tab:blue"),
        (onehz, "1 Hz 假信号", "tab:orange"),
        (high, "高本底 (SAA 西缘)", "tab:red"),
    ]

    fig = plt.figure(figsize=(15, 7.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.85, 1], hspace=0.72, wspace=0.26)

    ax = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
    # 只画卫星到得了的纬度带。SVOM 轨道倾角 30°，候选实测落在 ±29.1°，
    # 留 4° 余量；高纬那一大片本来就不可能有候选，画出来只是压扁地图。
    span = np.ceil(np.abs(lat).max()) + 4
    ax.set_extent([-180, 180, -span, span], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="0.95")
    ax.add_feature(cfeature.COASTLINE, lw=0.4, edgecolor="0.4")
    ax.gridlines(draw_labels=False, lw=0.3, color="0.88")
    # 其余候选按显著性着色。fa 跨 260 个量级，取 log 后再把最显著的一端压到
    # −60，否则少数极端值会把整条色标拉平、看不出多数候选的差别。
    fa = np.log10(col("fa"))
    pcm = ax.scatter(lon[tgf], lat[tgf], c=np.clip(fa[tgf], -60, None), s=17,
                     cmap="viridis_r", vmin=-60, vmax=-5, lw=0.2, edgecolor="0.3",
                     transform=ccrs.PlateCarree(),
                     label="TGF 候选 (%d)" % tgf.sum(), zorder=4)
    # 色标放地图正下方。不能用 colorbar(ax=...)：地图有固定纵横比，让它
    # 从自己的格子里让出空间会把地图整个缩小。先画一次拿到地图落定后的
    # 位置，再按这个位置摆一根等宽的色标轴。
    fig.canvas.draw()
    box = ax.get_position()
    cax = fig.add_axes([box.x0 + 0.15 * box.width, box.y0 - 0.085,
                        0.70 * box.width, 0.018])
    cb = fig.colorbar(pcm, cax=cax, orientation="horizontal")
    cb.set_label("TGF 候选的 $\\log_{10}$ fa（越小越显著）", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    if confirmed.any():
        ax.scatter(lon[confirmed], lat[confirmed], s=95, facecolor="none",
                   edgecolor="crimson", lw=1.1, transform=ccrs.PlateCarree(),
                   label="其中经闪电证实 (%d)" % confirmed.sum(), zorder=6)
    for mask, name, color, marker in [
        (onehz, "1 Hz 假信号", "tab:orange", "^"),
        (high, "高本底 (SAA 西缘)", "tab:red", "s"),
    ]:
        ax.scatter(lon[mask], lat[mask], s=22, c=color, marker=marker, lw=0.3,
                   edgecolor="k", transform=ccrs.PlateCarree(),
                   label="%s (%d)" % (name, mask.sum()), zorder=5)
    ax.set_title("(a) fa$\\leq10^{-5}$ 的 %d 个候选的地理分布" % len(rows), fontsize=12)
    ax.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.04),
              ncol=4, frameon=False)

    ax = fig.add_subplot(gs[1, 0])
    bins = np.linspace(0, 1, 51)
    ax.hist(phase, bins=bins, color="0.75", label="全部")
    ax.hist(phase[onehz], bins=bins, color="tab:orange", label="1 Hz 假信号")
    ax.axhline(len(rows) / 50, color="k", ls="--", lw=1, label="均匀分布的期望")
    ax.set_xlabel("整秒内的相位")
    ax.set_ylabel("候选数")
    ax.set_title("(b) 相位聚在 0.505 s", fontsize=10)
    ax.legend(fontsize=8)

    ax = fig.add_subplot(gs[1, 1])
    bins = np.linspace(0, 4, 41)
    for mask, name, color in classes:
        ax.hist(ratio[mask], bins=bins, histtype="step", lw=1.6, color=color, label=name)
    if confirmed.any():
        ax.hist(ratio[confirmed], bins=bins, color="crimson", alpha=0.55,
                label="经闪电证实 (%d)" % confirmed.sum())
    ax.set_xlabel("能道中位数比（窗内 / 本底）")
    ax.set_ylabel("候选数")
    ax.set_title("(c) 窗内事例的能道中位数 ÷ 本底的\n"
                 "已证实的与整类分不开（KS $p=0.15$）", fontsize=10)
    ax.legend(fontsize=8)

    ax = fig.add_subplot(gs[1, 2])
    bins = np.logspace(-2, 2, 41)
    for mask, name, color in classes:
        ax.hist(dur[mask], bins=bins, histtype="step", lw=1.6, color=color, label=name)
    if confirmed.any():
        ax.hist(dur[confirmed], bins=bins, color="crimson", alpha=0.55,
                label="经闪电证实 (%d)" % confirmed.sum())
    ax.set_xscale("log")
    ax.set_xlabel("窗长 (ms)")
    ax.set_ylabel("候选数")
    ax.set_title("(d) 窗长分布\n已证实的偏短（1.31 vs 1.49 ms, $p=0.002$）",
                 fontsize=10)
    ax.legend(fontsize=8)

    fig.suptitle("SVOM/GRM 792 天全量搜索：显著候选里的三个群体\n"
             "（TGF 候选一类由排除法得到，另经 WWLLN 闪电关联证实，见关联图）",
             fontsize=12)
    fig.savefig(args.output, dpi=140, bbox_inches="tight")
    print("wrote", args.output)
    for mask, name, _ in classes:
        print("  %-16s %4d  能道中位数比 %.2f  本底率中位 %7.0f c/s" %
              (name, mask.sum(), np.median(ratio[mask]), np.median(rate[mask])))


if __name__ == "__main__":
    main()
