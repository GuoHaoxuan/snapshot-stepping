#!/usr/bin/env python3
"""SVOM/GRM 候选的光变曲线与能谱。

上排是亚秒相位落在 0.5067 s 附近的候选（疑似 1 Hz 周期性伪信号），
下排是相位自由的候选。同一版式便于直接对比形态。

usage: plot_svom_lightcurves.py <lightcurves_csv> -o <png>
"""
import argparse, csv, sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import pubstyle

META = {
    "A_phase_2025-11-12": ("2025-11-12 03:28:41.507", "phase 0.5070", 164, 4.4e-266, -46.3, 4.9),
    "B_phase_2026-01-20": ("2026-01-20 17:44:15.507", "phase 0.5068", 91, 1.1e-115, -46.6, 5.8),
    "C_free_2025-06-08": ("2025-06-08 10:43:55.345", "phase 0.3445", 118, 2.9e-98, 137.4, -25.4),
    "D_free_2025-03-31": ("2025-03-31 12:09:17.077", "phase 0.0772", 76, 4.2e-90, 13.8, -3.6),
}
ORDER = ["A_phase_2025-11-12", "B_phase_2026-01-20", "C_free_2025-06-08", "D_free_2025-03-31"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("-o", "--out", default="svom_lightcurves.png")
    ap.add_argument("--bin-ms", type=float, default=0.2)
    args = ap.parse_args()
    pubstyle.apply()

    ev = defaultdict(lambda: ([], []))
    for row in csv.DictReader(open(args.csv)):
        t, p = ev[row["tag"]]
        t.append(float(row["dt_s"]) * 1e3)
        p.append(int(row["pi"]))

    fig, axes = plt.subplots(2, 2, figsize=(pubstyle.FULL_W, 3.6), sharex=True)
    bw = args.bin_ms
    bins = np.arange(-30, 30 + bw, bw)
    for ax, tag in zip(axes.ravel(), ORDER):
        t = np.array(ev[tag][0])
        label, phase, count, fa, lon, lat = META[tag]
        colour = pubstyle.C_EXT2 if "phase 0.50" in phase else pubstyle.C_OBS
        ax.stairs(np.histogram(t, bins=bins)[0], bins, color=colour, lw=0.8)
        rate = len(t) / 0.060
        ax.set_title("%s   %s\n$fa$=%.0e  count=%d  (%.0f, %.0f)$^\\circ$  bkg %.0f c/s"
                     % (label, phase, fa, count, lon, lat, rate), fontsize=6.5)
        ax.set_xlim(-30, 30)
    for ax in axes[1]:
        ax.set_xlabel("time since candidate start (ms)")
    for ax in axes[:, 0]:
        ax.set_ylabel("counts / %.1f ms" % bw)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print("saved", args.out)

    print("\n== 窗内特征 ==")
    for tag in ORDER:
        t = np.array(ev[tag][0]); p = np.array(ev[tag][1])
        core = np.abs(t) <= 1.0
        print("%-20s 本底 %5.0f c/s | 核心±1ms %3d 事例  PI中位 %3d | 窗外 PI中位 %3d"
              % (tag, len(t)/0.060, core.sum(),
                 int(np.median(p[core])) if core.sum() else -1,
                 int(np.median(p[~core]))))


if __name__ == "__main__":
    main()
