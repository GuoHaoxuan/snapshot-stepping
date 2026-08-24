#!/usr/bin/env python3
"""Estimated TGFs, v6 final: the paper extrapolation with the train cut applied.

Method is identical to plot_fp_distribution_v4_paperfmt.py (the run that gave
the REP-inflated 87,083): 200 log bins from the data minimum to fpy=20,
significant-side power-law fit (1e-50..1e-8) to the all-candidates curve,
then the fit summed over same-density log bins extended out to fpy=3.16e11
(the physical ceiling: one 100-us TGF per moment).  The only change is the
input population -- v6 candidates with train members (neighbors_10min > 34,
the frozen REP cut) removed.  The with-REP number is printed alongside on the
same binning to expose the inflation factor.

Figure follows the paper format (20x7 cm, reversed log x axis, shaded
extrapolation), drawn for the after-cut catalog.

usage: estimate_tgfs_v6.py <tgfs_v6.csv> [--full] [-o out.pdf]
       (default window < 2025-01-01 to compare with the paper and the v4
        87,083; --full uses the whole mission)
"""
import argparse
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def power_law(x, a, b):
    return a * x**b


def load(path):
    day, fpy, assoc, coinc, train = [], [], [], [], []
    with open(path) as f:
        f.readline()
        for line in f:
            d, y, a, c, _, t = line.rstrip("\n").split(",")
            day.append(int(d))
            fpy.append(float(y))
            assoc.append(a == "1")
            coinc.append(float(c))
            train.append(t == "1")
    return (np.array(day), np.array(fpy), np.array(assoc),
            np.array(coinc), np.array(train))


def estimate(n_all, centers, min_fp, max_fp, bins, p0):
    condition = (centers < 1e-8) & (centers > 1e-50)
    fit, _ = curve_fit(power_law, centers[condition], n_all[condition], p0=p0)
    n_ext = int(np.round((bins + 1) * (np.log10(max_fp) - np.log10(min_fp))
                         / (np.log10(20) - np.log10(min_fp))))
    fp_bins_ext = np.logspace(np.log10(min_fp), np.log10(max_fp), n_ext)
    total = float(power_law((fp_bins_ext[:-1] + fp_bins_ext[1:]) / 2, *fit).sum())
    return fit, total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("v6_csv")
    ap.add_argument("--full", action="store_true",
                    help="full mission window instead of < 2025-01-01")
    ap.add_argument("-o", "--out", default="fp-distribution-v6-final.pdf")
    args = ap.parse_args()

    day, fpy, assoc, coinc, train = load(args.v6_csv)
    window = np.ones_like(train) if args.full else day < 20250101
    tag = "full mission" if args.full else "< 2025-01-01"

    keep = window & ~train
    fpy_all = fpy[keep]
    fpy_rep = fpy[window]          # REP included, same window -- for comparison
    fpy_assoc = fpy[keep & assoc]
    assoc_pairs = list(zip(fpy[keep], coinc[keep]))

    n_zero = int((fpy_all == 0).sum())
    positive_min = min(fpy_all[fpy_all > 0].min(), fpy_assoc[fpy_assoc > 0].min())
    min_fp = float(positive_min)
    fpy_all = np.clip(fpy_all, min_fp, None)
    fpy_rep = np.clip(fpy_rep, min_fp, None)
    max_fp = 3.16e11
    bins = 200
    fp_bins = np.logspace(np.log10(min_fp), np.log10(20), bins + 1)
    centers = (fp_bins[:-1] + fp_bins[1:]) / 2

    plt.rcParams.update(
        {
            "text.usetex": shutil.which("latex") is not None,
            "font.family": "serif",
            "font.serif": ["Computer Modern"],
            "text.latex.preamble": "\\usepackage{amsmath}",
            "lines.linewidth": 1,
        }
    )
    cm = 1 / 2.54
    plt.figure(figsize=(20 * cm, 7 * cm), dpi=1200)

    n_rep, _, _ = plt.hist(fpy_rep, bins=fp_bins, histtype="step",
                           edgecolor="#BBBBBB", lw=0.8)
    n_all, _, _ = plt.hist(fpy_all, bins=fp_bins, histtype="step",
                           edgecolor="C0", lw=1.1)

    misassociated_count = np.zeros_like(fp_bins[:-1])
    for fp_year, prob in assoc_pairs:
        index = np.digitize(fp_year, fp_bins) - 1
        if 0 <= index < len(misassociated_count):
            misassociated_count[index] += prob
    plt.stairs(misassociated_count, fp_bins, fill=False, edgecolor="C1")

    n_associated, _, _ = plt.hist(
        fpy_assoc, bins=fp_bins, histtype="step", edgecolor="C2", alpha=0.5
    )
    plt.stairs(n_associated - misassociated_count, fp_bins, fill=False, edgecolor="C2")

    # 弱显著侧拟合（画参考线用）
    condition = (centers > 1e-4) & (centers < 20)
    fit_left, _ = curve_fit(power_law, centers[condition], n_all[condition])
    x_fit = np.logspace(np.log10(min_fp), np.log10(max_fp), 100)
    plt.plot(x_fit, power_law(x_fit, *fit_left), color="#CCCCCC", ls="--", zorder=-1)

    # 显著侧拟合 + 外推积分：目录版（train 摘除）与含 REP 版同 binning 对比
    fit_right, estimated = estimate(n_all, centers, min_fp, max_fp, bins, fit_left)
    fit_rep, estimated_rep = estimate(n_rep, centers, min_fp, max_fp, bins, fit_left)
    y_fit = power_law(x_fit, *fit_right)
    plt.plot(x_fit, y_fit, color="#CCCCCC", ls="--", zorder=-1)
    plt.fill_between(x_fit, y_fit, 1e-1, color="C0", alpha=0.1, zorder=-2)

    # 关联曲线拟合
    condition = (centers < 1e-2) & (centers > 1e-50)
    fit_assoc, _ = curve_fit(
        power_law, centers[condition], n_associated[condition], p0=fit_right
    )
    plt.plot(x_fit, power_law(x_fit, *fit_assoc), color="#CCCCCC", ls="--", zorder=-1)

    plt.ylim(0.5, 1e6)
    legend_handles = [
        mpatches.Patch(edgecolor="C0", facecolor="None",
                       label="TGF Candidates (v6, train cut)"),
        mpatches.Patch(edgecolor="#BBBBBB", facecolor="None",
                       label="$\\cdots$ before train cut"),
        mpatches.Patch(edgecolor="C1", facecolor="None",
                       label="Mis-associated TGF Candidates"),
        plt.Line2D([0], [0], color="#CCCCCC", linestyle="--", label="Power Law Fits"),
        mpatches.Patch(edgecolor="C2", facecolor="None", alpha=0.5,
                       label="TGF Candidates with Lightning"),
        mpatches.Patch(edgecolor="C2", facecolor="None",
                       label="$\\cdots$ But Mis-associated Excluded"),
        mpatches.Patch(facecolor="C0", edgecolor="None", alpha=0.1,
                       label="Estimated TGFs"),
    ]
    plt.legend(handles=legend_handles, ncols=2, loc="upper right", fontsize=7)
    plt.xlabel("Expected Annual False Positive Under Poisson Assumption")
    plt.ylabel("Number")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(max_fp, min_fp)
    plt.savefig(args.out, bbox_inches="tight")
    if args.out.endswith(".pdf"):
        plt.savefig(args.out[:-4] + ".png", dpi=200, bbox_inches="tight")
    print("saved", args.out)

    print(f"window: {tag}")
    print(f"candidates: {len(fpy_all)} after cut / {len(fpy_rep)} with REP "
          f"(zeros clipped: {n_zero})")
    print("fit right (a,b)  after cut =", fit_right)
    print("fit right (a,b)  with REP  =", fit_rep)
    print("fit assoc (a,b)            =", fit_assoc)
    print(f"Estimated TGFs  after cut : {estimated:.0f}")
    print(f"Estimated TGFs  with REP  : {estimated_rep:.0f}  "
          f"(inflation x{estimated_rep / estimated:.2f})")
    print(f"associated: {len(fpy_assoc)}  expected misassoc: "
          f"{float(np.sum(coinc[keep & assoc])):.1f}")


if __name__ == "__main__":
    main()
