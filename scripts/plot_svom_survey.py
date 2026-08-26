#!/usr/bin/env python3
"""SVOM/GRM 首轮全量搜索的概览图：fp 分布、时间分布、地理分布。

输入 sig_all_svom.csv（blink search 的逐日 signals.json 汇总）与
exposure_svom.csv（逐年曝光秒数，率的分母）。

usage: plot_svom_survey.py <sig_csv> <exposure_csv> -o <前缀>
"""
import argparse, csv, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).parent))
import pubstyle

SELECT = 1e-5          # 判选 A 层阈值（沿用论文口径）
C_ALL = "0.55"
C_SEL = pubstyle.C_EXT2   # SVOM/GRM 的配色角色


def load(path):
    rows = list(csv.DictReader(open(path)))
    t = np.array([datetime.strptime(r["start"][:26] + "Z", "%Y-%m-%dT%H:%M:%S.%fZ")
                  .replace(tzinfo=timezone.utc) for r in rows])
    return dict(
        t=t,
        fa=np.array([float(r["false_positive_per_year"]) for r in rows]),
        count=np.array([int(r["count"]) for r in rows]),
        bin_us=np.array([float(r["bin_size_best"]) for r in rows]) * 1e6,
        lon=np.array([float(r["lon"]) for r in rows]),
        lat=np.array([float(r["lat"]) for r in rows]),
    )


def fig_fp(d, out):
    # 绝大多数候选落在 1e-60 以右；更显著的那条长尾单独用注记交代，
    # 免得把横轴拉到 1e-160 之后所有结构都挤成一根线。
    floor = 1e-60
    fa = np.clip(d["fa"], floor, None)
    bins = np.logspace(np.log10(floor), np.log10(20.0), 120)
    fig, ax = plt.subplots(figsize=(pubstyle.FULL_W * 0.62, 2.6))
    n, _ = np.histogram(fa, bins=bins)
    ax.stairs(n, bins, color=C_SEL, lw=1.1, label="all candidates (%d)" % len(fa))
    ax.axvline(SELECT, color="crimson", lw=0.7, ls=":")
    ax.annotate("selection  $fa\\leq10^{-5}$  (%d)" % int((d["fa"] <= SELECT).sum()),
                xy=(SELECT, 0.97), xycoords=("data", "axes fraction"),
                xytext=(-4, 0), textcoords="offset points",
                color="crimson", fontsize=6.5, ha="right", va="top")
    n_tail = int((d["fa"] < floor).sum())
    ax.annotate("%d candidates below $10^{-60}$\n(most significant $4\\times10^{-266}$)" % n_tail,
                xy=(0.02, 0.62), xycoords="axes fraction", fontsize=6.5, color="0.35")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(floor / 3, 40); ax.set_ylim(0.7, None)
    ax.set_xlabel("false positives per year")
    ax.set_ylabel("candidates per bin")
    ax.set_title("SVOM/GRM first full search: fp-distribution", fontsize=8)
    ax.legend(loc="upper left")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print("saved", out)


def fig_time(d, exposure_csv, out):
    sel = d["fa"] <= SELECT
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(pubstyle.FULL_W * 0.62, 3.4),
                                  sharex=True, height_ratios=[2, 1])
    days_all = np.array([t.date().toordinal() for t in d["t"]])
    lo, hi = days_all.min(), days_all.max()
    edges = np.arange(lo, hi + 2)
    n_all, _ = np.histogram(days_all, bins=edges)
    n_sel, _ = np.histogram(days_all[sel], bins=edges)
    x = np.array([datetime.fromordinal(int(e)) for e in edges[:-1]])
    ax.plot(x, n_all, color=C_ALL, lw=0.5, label="all (%d)" % len(days_all))
    ax.set_ylabel("candidates / day")
    ax.set_title("SVOM/GRM candidates over the mission", fontsize=8)
    ax.legend(loc="upper left")
    top = np.argsort(n_all)[-3:][::-1]
    for i in top:
        ax.annotate(x[i].strftime("%Y-%m-%d"), xy=(x[i], n_all[i]),
                    xytext=(0, 3), textcoords="offset points",
                    fontsize=6, ha="center", color="crimson")
    ax2.plot(x, n_sel, color=C_SEL, lw=0.6,
             label="$fa\\leq10^{-5}$ (%d)" % int(sel.sum()))
    ax2.set_ylabel("selected / day")
    ax2.legend(loc="upper left")
    ax2.set_xlabel("date (UTC)")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print("saved", out)


def fig_geo(d, out):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    sel = d["fa"] <= SELECT
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(pubstyle.FULL_W * 0.8, 3.9))
    # 地图与直方图共用同一条经度轴，所以手工给定两块画布的位置，
    # 让它们左右边界严格对齐（GeoAxes 的 aspect 会自作主张，得关掉）。
    ax = fig.add_axes([0.08, 0.50, 0.90, 0.42], projection=proj)
    ax2 = fig.add_axes([0.08, 0.13, 0.90, 0.30])
    ax.add_feature(cfeature.LAND, facecolor="0.94", zorder=-2)
    ax.add_feature(cfeature.COASTLINE, edgecolor="0.55", linewidth=0.35, zorder=-1)
    ax.set_extent([-180, 180, -33, 33], crs=proj)
    ax.set_aspect("auto")
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="0.85", alpha=0.8)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": 6.5}
    ax.scatter(d["lon"][~sel], d["lat"][~sel], s=0.7, c=C_ALL, lw=0,
               alpha=0.30, rasterized=True, transform=proj,
               label="sub-threshold (%d)" % int((~sel).sum()))
    ax.scatter(d["lon"][sel], d["lat"][sel], s=4.0, c=C_SEL, lw=0,
               transform=proj, label="$fa\\leq10^{-5}$ (%d)" % int(sel.sum()))
    ax.set_title("SVOM/GRM candidate positions (sub-satellite point)", fontsize=8)
    ax.legend(loc="lower left", ncol=2, framealpha=0.9, fontsize=6)
    # 经度直方图：真 TGF 应当在非洲/南美/东南亚三个雷暴区堆出峰
    b = np.arange(-180, 181, 5)
    ax2.set_xlim(-180, 180)
    ax2.stairs(np.histogram(d["lon"][sel], bins=b)[0], b, color=C_SEL, lw=1.0,
               label="selected, per 5$^\\circ$ longitude")
    for name, x0 in (("America", -75), ("Africa", 22), ("Maritime Cont.", 122)):
        ax2.axvline(x0, color="0.7", lw=0.5, ls="--", zorder=-1)
        ax2.annotate(name, xy=(x0, 0.97), xycoords=("data", "axes fraction"),
                     fontsize=6, ha="center", va="top", color="0.35")
    ax2.set_xlabel("longitude (deg)"); ax2.set_ylabel("selected / 5$^\\circ$")
    ax2.set_ylim(0, None)
    ax2.legend(loc="upper right", fontsize=6, framealpha=0.9)
    fig.savefig(out, dpi=200)
    print("saved", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sig_csv")
    ap.add_argument("exposure_csv")
    ap.add_argument("-o", "--out", default="svom_survey")
    args = ap.parse_args()
    pubstyle.apply()
    d = load(args.sig_csv)
    fig_fp(d, args.out + "_fp.png")
    fig_time(d, args.exposure_csv, args.out + "_time.png")
    fig_geo(d, args.out + "_geo.png")

    sel = d["fa"] <= SELECT
    print("\n== 概览 ==")
    print("候选 %d，其中 fa<=1e-5 的 %d" % (len(d["fa"]), int(sel.sum())))
    print("选中候选 bin_size_best: 中位 %.0f us, 四分位 %.0f–%.0f us"
          % (np.median(d["bin_us"][sel]), *np.percentile(d["bin_us"][sel], [25, 75])))
    print("选中候选 |lat|<=30 占比 %.1f%%" % (100 * (np.abs(d["lat"][sel]) <= 30).mean()))


if __name__ == "__main__":
    main()
