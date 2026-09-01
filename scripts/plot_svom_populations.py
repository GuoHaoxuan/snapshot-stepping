"""SVOM/GRM 全量搜索结果里三个群体的分离：1 Hz 假信号、高本底、TGF 候选。

用法:
    python3 scripts/plot_svom_populations.py <features.csv> -o <PNG>
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

RATE_HIGH = 2e4          # 本底率高于此即在 SAA 西缘，谱硬度恒为 1.00
PHASE_LO, PHASE_HI = 0.49, 0.52   # 1 Hz 假信号聚集的整秒相位


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.csv)))
    col = lambda k, t=float: np.array([t(r[k]) for r in rows])
    lon, lat = col("lon"), col("lat")
    phase, rate, hard = col("phase"), col("rate_bkg"), col("hardness")
    dur = col("dur_ms")

    high = rate > RATE_HIGH
    onehz = (~high) & (phase >= PHASE_LO) & (phase <= PHASE_HI)
    tgf = ~(high | onehz)
    classes = [
        (tgf, "TGF 候选", "tab:blue"),
        (onehz, "1 Hz 假信号", "tab:orange"),
        (high, "高本底 (SAA 西缘)", "tab:red"),
    ]

    fig = plt.figure(figsize=(15, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.28, wspace=0.26)

    ax = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor="0.95")
    ax.add_feature(cfeature.COASTLINE, lw=0.4, edgecolor="0.4")
    ax.gridlines(draw_labels=False, lw=0.3, color="0.88")
    for mask, name, color in classes:
        ax.scatter(lon[mask], lat[mask], s=16, c=color, lw=0.2, edgecolor="k",
                   alpha=0.85, transform=ccrs.PlateCarree(),
                   label="%s (%d)" % (name, mask.sum()), zorder=4)
    ax.set_title("(a) fa$\\leq10^{-5}$ 的 %d 个候选的地理分布" % len(rows), fontsize=12)
    ax.legend(fontsize=9, loc="lower left", ncol=3)

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
        ax.hist(hard[mask], bins=bins, histtype="step", lw=1.6, color=color, label=name)
    ax.set_xlabel("硬度比")
    ax.set_ylabel("候选数")
    ax.set_title("(c) 硬度比：TGF 候选偏硬", fontsize=10)
    ax.legend(fontsize=8)

    ax = fig.add_subplot(gs[1, 2])
    bins = np.logspace(-2, 2, 41)
    for mask, name, color in classes:
        ax.hist(dur[mask], bins=bins, histtype="step", lw=1.6, color=color, label=name)
    ax.set_xscale("log")
    ax.set_xlabel("窗长 (ms)")
    ax.set_ylabel("候选数")
    ax.set_title("(d) 窗长分布", fontsize=10)
    ax.legend(fontsize=8)

    fig.suptitle("SVOM/GRM 792 天全量搜索：显著候选里的三个群体", fontsize=13)
    fig.savefig(args.output, dpi=140, bbox_inches="tight")
    print("wrote", args.output)
    for mask, name, _ in classes:
        print("  %-16s %4d  硬度比中位 %.2f  本底率中位 %7.0f c/s" %
              (name, mask.sum(), np.median(hard[mask]), np.median(rate[mask])))


if __name__ == "__main__":
    main()
