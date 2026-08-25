#!/usr/bin/env python3
"""TGF counts over time under the paper's selection, v5 catalog vs paper-era.

Selection (identical to the catalog paper): region 1 = fpy < 1e-5;
region 2 = lightning-associated with fpy in [1e-5, 1). Region 2 is empty
after 2024 (WWLLN library ends there).

Top panel: monthly selected-TGF counts, v5 (full mission) with the paper-era
catalog overlaid (its window ends 2024-12-31). Bottom panel: per-day counts,
log scale -- REP contamination shows up as isolated single-day spikes, real
TGFs as the smooth seasonal band.

usage: plot_tgf_time_distribution.py <paper.csv> <v5.csv> -o <png>
"""
import argparse
import datetime

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    day, fpy, assoc = [], [], []
    with open(path) as f:
        first = f.readline()
        if not first.startswith("day"):
            f.seek(0)
        for line in f:
            d, y, a, _ = line.rstrip("\n").split(",")
            day.append(int(d))
            fpy.append(float(y))
            assoc.append(a == "1")
    return np.array(day), np.array(fpy), np.array(assoc)


def select(day, fpy, assoc):
    r1 = fpy < 1e-5
    r2 = (fpy >= 1e-5) & (fpy < 1.0) & assoc
    return day[r1 | r2]


def to_dates(days):
    return np.array([datetime.date(d // 10000, d // 100 % 100, d % 100) for d in days])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paper_csv")
    ap.add_argument("v5_csv")
    ap.add_argument("-o", "--out", default="tgf_time_distribution.png")
    args = ap.parse_args()

    sel_paper = select(*load(args.paper_csv))
    sel_v5 = select(*load(args.v5_csv))
    print(f"paper selected: {len(sel_paper)}  v5 selected: {len(sel_v5)}")

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(12.5, 7.5), sharex=True,
        gridspec_kw=dict(height_ratios=[1.0, 1.0], hspace=0.08),
    )

    # 上图：月计数
    month_edges = [datetime.date(y, m, 1)
                   for y in range(2017, 2028) for m in range(1, 13)]
    month_nums = mdates.date2num(month_edges)
    for label, days, style in (
        ("v5, paper selection (%d)" % len(sel_v5), sel_v5,
         dict(color="C0", lw=1.3)),
        ("paper catalog (%d, ends 2024-12)" % len(sel_paper), sel_paper,
         dict(color="k", lw=1.0)),
    ):
        n, _ = np.histogram(mdates.date2num(to_dates(days)), bins=month_nums)
        ax1.stairs(n, month_nums, label=label, **style)
    ax1.set_ylabel("selected TGFs per month")
    ax1.legend(loc="upper left", fontsize=9)
    ax1.grid(alpha=0.25)

    # 下图：逐日计数（对数），风暴日以孤立尖峰现形
    days_u, counts = np.unique(sel_v5, return_counts=True)
    dates_u = to_dates(days_u)
    ax2.vlines(mdates.date2num(dates_u), 0.7, counts, color="C0", lw=0.8, alpha=0.8)
    ax2.set_yscale("log")
    ax2.set_ylim(0.7, counts.max() * 2)
    ax2.set_ylabel("selected TGFs per day (v5)")
    ax2.grid(alpha=0.25)

    order = np.argsort(counts)[::-1]
    top = order[:12]
    for i in top[:6]:
        ax2.annotate(str(days_u[i]), (mdates.date2num(dates_u[i]), counts[i]),
                     textcoords="offset points", xytext=(0, 4),
                     ha="center", fontsize=7, rotation=45)
    ax2.xaxis.set_major_locator(mdates.YearLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print("saved", args.out)

    print("\ntop days (v5, paper selection):")
    csum = 0
    for rank, i in enumerate(top, 1):
        csum += counts[i]
        print(f"  {rank:2d}. {days_u[i]}  {counts[i]:4d}")
    print(f"top 12 days: {counts[top].sum()} of {len(sel_v5)} "
          f"({100.0 * counts[top].sum() / len(sel_v5):.1f}%)")
    for k in (1, 5, 10, 20, 50):
        s = counts[order[:k]].sum()
        print(f"  top {k:3d} days -> {s:5d} ({100.0 * s / len(sel_v5):.1f}%)")


if __name__ == "__main__":
    main()
