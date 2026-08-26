#!/usr/bin/env python3
"""SVOM/GRM 候选的亚秒相位诊断。

真事件的发生时刻在一秒内均匀分布；固定相位的堆积只可能来自星上的周期性
动作。本图把相位分布与地理位置对照，用来定位伪信号的来源。

usage: plot_svom_phase.py <sig_csv> -o <png>
"""
import argparse, csv, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import pubstyle

SELECT = 1e-5
PHASE_LO, PHASE_HI = 0.500, 0.510


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sig_csv"); ap.add_argument("-o", "--out", default="svom_phase.png")
    args = ap.parse_args()
    pubstyle.apply()

    r = list(csv.DictReader(open(args.sig_csv)))
    fa = np.array([float(x["false_positive_per_year"]) for x in r])
    lon = np.array([float(x["lon"]) for x in r])
    lat = np.array([float(x["lat"]) for x in r])
    ph = np.array([float(x["start"][19:26]) for x in r])
    sel = fa <= SELECT
    onph = sel & (ph >= PHASE_LO) & (ph < PHASE_HI)
    offph = sel & ~((ph >= PHASE_LO) & (ph < PHASE_HI))

    import cartopy.crs as ccrs, cartopy.feature as cfeature
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(pubstyle.FULL_W * 0.8, 4.3))
    ax1 = fig.add_axes([0.09, 0.60, 0.88, 0.33])
    ax2 = fig.add_axes([0.09, 0.10, 0.88, 0.36], projection=proj)

    # 上：相位直方图
    b = np.arange(0, 1.001, 0.005)
    ax1.stairs(np.histogram(ph[sel], bins=b)[0], b, color=pubstyle.C_EXT2, lw=0.9,
               label="$fa\\leq10^{-5}$ (%d)" % int(sel.sum()))
    ax1.axhline(sel.sum() / 200, color="0.6", lw=0.6, ls="--",
                label="uniform expectation")
    ax1.set_xlim(0, 1); ax1.set_xlabel("sub-second phase of candidate start (s)")
    ax1.set_ylabel("candidates / 5 ms")
    ax1.annotate("%d candidates in [%.3f, %.3f]\n= %.0f$\\times$ uniform"
                 % (int(onph.sum()), PHASE_LO, PHASE_HI,
                    onph.sum() / (sel.sum() * (PHASE_HI - PHASE_LO))),
                 xy=(PHASE_HI, 0.92), xycoords=("data", "axes fraction"),
                 xytext=(8, 0), textcoords="offset points",
                 fontsize=6.5, color="crimson", va="top")
    ax1.legend(loc="upper left", fontsize=6)
    ax1.set_title("SVOM/GRM: a 1 Hz artefact hides among the significant candidates",
                  fontsize=8)

    # 下：按相位着色的地理分布
    ax2.add_feature(cfeature.LAND, facecolor="0.94", zorder=-2)
    ax2.add_feature(cfeature.COASTLINE, edgecolor="0.55", linewidth=0.35, zorder=-1)
    ax2.set_extent([-180, 180, -33, 33], crs=proj); ax2.set_aspect("auto")
    gl = ax2.gridlines(draw_labels=True, linewidth=0.3, color="0.85")
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": 6.5}
    ax2.scatter(lon[offph], lat[offph], s=3.5, c=pubstyle.C_EXT2, lw=0, transform=proj,
                label="phase elsewhere (%d)" % int(offph.sum()))
    ax2.scatter(lon[onph], lat[onph], s=9, c="crimson", lw=0, marker="X", transform=proj,
                label="phase in [%.2f, %.2f] (%d)" % (PHASE_LO, PHASE_HI, int(onph.sum())))
    ax2.legend(loc="lower left", fontsize=6, framealpha=0.9, ncol=2)
    fig.savefig(args.out, dpi=200)
    print("saved", args.out)

    print("\n相位内候选的经度分布（每 30 度）:")
    h, e = np.histogram(lon[onph], bins=np.arange(-180, 181, 30))
    for c, lo in zip(h, e[:-1]):
        if c: print("   %4d..%4d : %s (%d)" % (lo, lo + 30, "#" * c, c))


if __name__ == "__main__":
    main()
