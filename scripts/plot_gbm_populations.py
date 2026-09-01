"""GBM 一天的候选构成：带电粒子与 TGF 候选在同时性、时间结构与地理上的分离。

用法:
    python3 scripts/plot_gbm_populations.py <cand_detstats.csv> -o <PNG>
"""
import argparse, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.sans-serif": ["PingFang SC", "Arial Unicode MS"],
    "font.family": "sans-serif",
    "axes.unicode_minus": False,
})
import cartopy.crs as ccrs
import cartopy.feature as cfeature

CUT = 0.35

def load(path):
    rows = list(csv.DictReader(open(path)))
    col = lambda k, t=int: np.array([t(r[k]) for r in rows])
    d = {k: col(k) for k in ("count", "n_win", "n_nai", "n_b0", "n_b1",
                             "n_det", "n_times", "max_mult", "pha_med")}
    for k in ("dur_us", "fa", "lon", "lat"):
        d[k] = col(k, float)
    d["frac"] = d["max_mult"] / np.maximum(d["n_win"], 1)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    d = load(args.csv)

    keep = d["frac"] <= CUT
    sig = d["fa"] <= 1e-5
    fig = plt.figure(figsize=(17, 7.6))
    gs = fig.add_gridspec(2, 6, height_ratios=[1.5, 1], hspace=0.30, wspace=1.1)

    # (a) 同时性比例的分布
    ax = fig.add_subplot(gs[0, 0:2])
    bins = np.linspace(0, 1, 51)
    ax.hist(d["frac"], bins=bins, color="0.75", label="全部候选 (%d)" % len(d["frac"]))
    ax.hist(d["frac"][sig], bins=bins, color="crimson", alpha=0.85,
            label="fa$\\leq10^{-5}$ (%d)" % sig.sum())
    ax.axvline(CUT, color="k", ls="--", lw=1)
    ax.set_yscale("log")
    ax.set_xlabel("最大同时重数 / 窗内计数")
    ax.set_ylabel("候选数")
    ax.set_title("(a) 同时性比例：池子是连续的,\n显著候选才分成两群", fontsize=10)
    ax.legend(fontsize=8)

    # (b) 时长 vs 同时性
    ax = fig.add_subplot(gs[0, 2:4])
    ax.scatter(d["dur_us"][~sig], d["frac"][~sig], s=2, c="0.8", lw=0, rasterized=True)
    ax.scatter(d["dur_us"][sig & ~keep], d["frac"][sig & ~keep], s=42, c="tab:orange",
               marker="x", label="显著·判为粒子 (%d)" % (sig & ~keep).sum())
    ax.scatter(d["dur_us"][sig & keep], d["frac"][sig & keep], s=52, c="crimson",
               edgecolor="k", lw=0.5, label="显著·TGF 候选 (%d)" % (sig & keep).sum())
    ax.axhline(CUT, color="k", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("窗长 (µs)")
    ax.set_ylabel("最大同时重数 / 窗内计数")
    ax.set_title("(b) 两群在时长—同时性平面上分开", fontsize=10)
    ax.legend(fontsize=8, loc="upper right")

    # (c) 粒子占比随纬度——比值不需要曝光归一，宇宙线的地磁依赖应当在这里露出来
    ax = fig.add_subplot(gs[0, 4:6])
    edges = np.linspace(-27, 27, 28)
    mid = 0.5 * (edges[:-1] + edges[1:])
    total, _ = np.histogram(d["lat"], bins=edges)
    particle, _ = np.histogram(d["lat"][~keep], bins=edges)
    good = total > 50
    ratio = np.where(good, particle / np.maximum(total, 1), np.nan)
    err = np.where(good, np.sqrt(particle) / np.maximum(total, 1), np.nan)
    ax.errorbar(mid, ratio, yerr=err, fmt="o-", ms=3, lw=1, color="tab:orange")
    ax.set_xlabel("地理纬度 (deg)")
    ax.set_ylabel("判为粒子的占比")
    ax.set_title("(c) 粒子占比随纬度的变化", fontsize=10)
    ax.grid(alpha=0.3)

    # (d)(e) 地理密度
    lon_edges = np.linspace(-180, 180, 73)
    lat_edges = np.linspace(-30, 30, 25)
    for i, (mask, title, cmap) in enumerate([
        (~keep, "(d) 判为带电粒子的候选 (%d)" % (~keep).sum(), "Oranges"),
        (keep, "(e) 通过同时性判据的候选 (%d)" % keep.sum(), "Blues"),
    ]):
        ax = fig.add_subplot(gs[1, i * 3 : i * 3 + 3], projection=ccrs.PlateCarree())
        ax.set_extent([-180, 180, -35, 35], crs=ccrs.PlateCarree())
        hist, _, _ = np.histogram2d(d["lon"][mask], d["lat"][mask],
                                    bins=[lon_edges, lat_edges])
        pcm = ax.pcolormesh(lon_edges, lat_edges, hist.T, cmap=cmap,
                            transform=ccrs.PlateCarree(), shading="auto")
        ax.add_feature(cfeature.COASTLINE, lw=0.4, edgecolor="0.3")
        both = mask & sig
        ax.scatter(d["lon"][both], d["lat"][both], s=60, c="crimson", marker="*",
                   edgecolor="k", lw=0.4, transform=ccrs.PlateCarree(),
                   label="fa$\\leq10^{-5}$ (%d)" % both.sum(), zorder=5)
        fig.colorbar(pcm, ax=ax, orientation="horizontal", pad=0.04, fraction=0.05,
                     label="每 5°×2.5° 格子的候选数")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8, loc="lower left")

    fig.suptitle("Fermi/GBM 2019-01-01 单日候选：带电粒子与 TGF 候选的分离", fontsize=13)
    fig.savefig(args.output, dpi=140, bbox_inches="tight")
    print("wrote", args.output)


if __name__ == "__main__":
    main()
