#!/usr/bin/env python3
"""fp-distribution for v6 with the frozen train cut (threshold 34).

Same conventions as plot_fp_distribution_v5_vs_paper.py (log-bin histogram of
false_positive_per_year, power-law fit on the associated curve), full mission
window (the REP train story includes 2025+, unlike the paper-comparison plots
which stop at 2025-01-01).  Curves:

  v5 all           -- reference (time-base-fixed catalog, no train fields)
  v6 all           -- identical search + 2026-08 days, train/acd enrichment
  v6 train members -- neighbors_10min > 34, what the frozen cut removes
  v6 after cut     -- the catalog population
  v6 associated    -- WWLLN-associated, with power-law fit (paper b = 0.039)

Inputs: tgfs_v5.csv (day,fpy,assoc,coinc) and tgfs_v6.csv (same +
neighbors_10min,is_train), both extracted from the farm tgfs.json.

usage: plot_fp_distribution_v6.py <v5.csv> <v6.csv> -o <png>
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def power_law(x, a, b):
    return a * x**b


def load(path, with_train):
    cols = 6 if with_train else 4
    day, fpy, assoc, train = [], [], [], []
    with open(path) as f:
        first = f.readline()
        if not first.startswith("day"):
            f.seek(0)
        for line in f:
            parts = line.rstrip("\n").split(",")[:cols]
            day.append(int(parts[0]))
            fpy.append(float(parts[1]))
            assoc.append(parts[2] == "1")
            train.append(with_train and parts[5] == "1")
    return (np.array(day), np.array(fpy), np.array(assoc), np.array(train))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("v5_csv")
    ap.add_argument("v6_csv")
    ap.add_argument("-o", "--out", default="fp_distribution_v6.png")
    args = ap.parse_args()

    _, fpy5, assoc5, _ = load(args.v5_csv, with_train=False)
    _, fpy6, assoc6, train6 = load(args.v6_csv, with_train=True)

    floor = 1e-160
    clip = lambda a: np.clip(a, floor, None)
    bins = np.logspace(np.log10(floor), np.log10(20.0), 200)
    centers = np.sqrt(bins[:-1] * bins[1:])

    fig, ax = plt.subplots(figsize=(11.5, 6))

    curves = [
        ("v5 all (%d)" % len(fpy5), clip(fpy5), dict(color="0.65", lw=1.0)),
        ("v6 all (%d)" % len(fpy6), clip(fpy6), dict(color="C0", lw=1.4)),
        ("v6 train members (%d)" % int(train6.sum()), clip(fpy6[train6]),
         dict(color="C3", lw=1.2)),
        ("v6 after cut (%d)" % int((~train6).sum()), clip(fpy6[~train6]),
         dict(color="k", lw=1.2)),
        ("v6 associated (%d)" % int(assoc6.sum()), clip(fpy6[assoc6]),
         dict(color="C2", lw=1.2)),
    ]
    hists = {}
    for label, arr, style in curves:
        n, _ = np.histogram(arr, bins=bins)
        hists[label] = n
        ax.stairs(n, bins, label=label, **style)

    # 幂律拟合：关联曲线显著侧（同 v4/v5 脚本的 1e-50..1e-2 窗）
    label = "v6 associated (%d)" % int(assoc6.sum())
    n = hists[label]
    sel = (centers > 1e-50) & (centers < 1e-2)
    ok = n[sel] > 0
    pa, _ = curve_fit(power_law, centers[sel][ok], n[sel][ok], p0=(10.0, 0.05))
    x_fit = np.logspace(np.log10(floor), np.log10(20.0), 200)
    ax.plot(x_fit, power_law(x_fit, *pa), ls="--", color="C2", lw=0.8,
            alpha=0.7, zorder=-1,
            label=r"associated power law ($b=%.4f$)" % pa[1])

    ax.axvline(1e-8, color="crimson", lw=0.7, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(floor / 3, 40)
    ax.set_ylim(0.7, 4e6)
    ax.set_xlabel("false positives per year")
    ax.set_ylabel("candidates per bin")
    ax.set_title("fp-distribution, v6 full mission, train cut frozen at 34")
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print("saved", args.out)

    for label, arr, _ in curves:
        print("%-28s fpy<1e-8: %6d   fpy<1e-30: %5d"
              % (label, int((arr < 1e-8).sum()), int((arr < 1e-30).sum())))


if __name__ == "__main__":
    main()
