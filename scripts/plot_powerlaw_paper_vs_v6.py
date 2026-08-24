#!/usr/bin/env python3
"""Paper-era power law against the v6 train-cut catalog, paper figure format.

Both catalogs restricted to < 2025-01-01 (paper window), same binning (200 log
bins from the joint positive minimum to fpy=20).  For each all-candidates
curve: the significant-side power-law fit (1e-50..1e-8) and its extrapolated
integral out to fpy=3.16e11 (the "Estimated TGFs" of the paper); the v6
associated curve and fit are drawn as the independent slope anchor.

Inputs are day,fpy,assoc,coinc CSVs:
  paper_catalog.csv  from /Volumes/Graphite/blink.db (paper-era catalog)
  tgfs_v6.csv        adds neighbors_10min,is_train (train cut applied here)

usage: plot_powerlaw_paper_vs_v6.py <paper.csv> <v6.csv> [-o out.pdf]
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


def load(path, with_train):
    day, fpy, assoc, train = [], [], [], []
    with open(path) as f:
        first = f.readline()
        if not first[:1].isdigit():
            pass  # header consumed
        else:
            f.seek(0)
        for line in f:
            parts = line.rstrip("\n").split(",")
            day.append(int(parts[0]))
            fpy.append(float(parts[1]))
            assoc.append(parts[2] == "1")
            train.append(with_train and parts[5] == "1")
    return np.array(day), np.array(fpy), np.array(assoc), np.array(train)


def fit_and_estimate(fpy_arr, centers, min_fp, max_fp, bins, fp_bins):
    n, _ = np.histogram(fpy_arr, bins=fp_bins)
    condition = (centers > 1e-4) & (centers < 20)
    p0, _ = curve_fit(power_law, centers[condition], n[condition])
    condition = (centers < 1e-8) & (centers > 1e-50)
    fit, _ = curve_fit(power_law, centers[condition], n[condition], p0=p0)
    n_ext = int(np.round((bins + 1) * (np.log10(max_fp) - np.log10(min_fp))
                         / (np.log10(20) - np.log10(min_fp))))
    ext = np.logspace(np.log10(min_fp), np.log10(max_fp), n_ext)
    total = float(power_law((ext[:-1] + ext[1:]) / 2, *fit).sum())
    return n, fit, total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paper_csv")
    ap.add_argument("v6_csv")
    ap.add_argument("-o", "--out", default="powerlaw-paper-vs-v6.pdf")
    args = ap.parse_args()

    day_p, fpy_p, _, _ = load(args.paper_csv, with_train=False)
    day_6, fpy_6, assoc_6, train_6 = load(args.v6_csv, with_train=True)

    fpy_paper = fpy_p[day_p < 20250101]
    keep = (day_6 < 20250101) & ~train_6
    fpy_v6 = fpy_6[keep]
    fpy_assoc = fpy_6[keep & assoc_6]

    min_fp = float(min(fpy_paper[fpy_paper > 0].min(), fpy_v6[fpy_v6 > 0].min()))
    fpy_paper = np.clip(fpy_paper, min_fp, None)
    fpy_v6 = np.clip(fpy_v6, min_fp, None)
    max_fp, bins = 3.16e11, 200
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

    n_paper, fit_paper, est_paper = fit_and_estimate(
        fpy_paper, centers, min_fp, max_fp, bins, fp_bins)
    n_v6, fit_v6, est_v6 = fit_and_estimate(
        fpy_v6, centers, min_fp, max_fp, bins, fp_bins)

    plt.stairs(n_paper, fp_bins, fill=False, edgecolor="k", lw=0.7, alpha=0.75)
    plt.stairs(n_v6, fp_bins, fill=False, edgecolor="C0", lw=1.1)
    n_assoc, _, _ = plt.hist(np.clip(fpy_assoc, min_fp, None), bins=fp_bins,
                             histtype="step", edgecolor="C2", alpha=0.6)

    x_fit = np.logspace(np.log10(min_fp), np.log10(max_fp), 100)
    plt.plot(x_fit, power_law(x_fit, *fit_paper), color="k", ls="--", lw=0.7,
             alpha=0.5, zorder=-1)
    plt.plot(x_fit, power_law(x_fit, *fit_v6), color="C0", ls="--", lw=0.8,
             zorder=-1)
    plt.fill_between(x_fit, power_law(x_fit, *fit_v6), 1e-1,
                     color="C0", alpha=0.1, zorder=-2)

    condition = (centers < 1e-2) & (centers > 1e-50)
    fit_assoc, _ = curve_fit(power_law, centers[condition], n_assoc[condition],
                             p0=fit_v6)
    plt.plot(x_fit, power_law(x_fit, *fit_assoc), color="C2", ls="--", lw=0.7,
             alpha=0.6, zorder=-1)

    legend_handles = [
        mpatches.Patch(edgecolor="k", facecolor="None", alpha=0.75,
                       label="Paper Search (%d)" % len(fpy_paper)),
        mpatches.Patch(edgecolor="C0", facecolor="None",
                       label="v6, train cut (%d)" % len(fpy_v6)),
        mpatches.Patch(edgecolor="C2", facecolor="None", alpha=0.6,
                       label="v6 with Lightning (%d)" % len(fpy_assoc)),
        plt.Line2D([0], [0], color="k", ls="--", lw=0.7, alpha=0.5,
                   label="Paper fit $b=%.4f$, Est.\\ %d" % (fit_paper[1], round(est_paper))),
        plt.Line2D([0], [0], color="C0", ls="--", lw=0.8,
                   label="v6 fit $b=%.4f$, Est.\\ %d" % (fit_v6[1], round(est_v6))),
        plt.Line2D([0], [0], color="C2", ls="--", lw=0.7, alpha=0.6,
                   label="Assoc.\\ fit $b=%.4f$" % fit_assoc[1]),
    ]
    plt.legend(handles=legend_handles, ncols=2, loc="upper right", fontsize=7)
    plt.xlabel("Expected Annual False Positive Under Poisson Assumption")
    plt.ylabel("Number")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(0.5, 1e6)
    plt.xlim(max_fp, min_fp)
    plt.savefig(args.out, bbox_inches="tight")
    if args.out.endswith(".pdf"):
        plt.savefig(args.out[:-4] + ".png", dpi=200, bbox_inches="tight")
    print("saved", args.out)
    print("paper: n=%d  fit(a,b)=%s  Estimated=%d" % (len(fpy_paper), fit_paper, round(est_paper)))
    print("v6   : n=%d  fit(a,b)=%s  Estimated=%d" % (len(fpy_v6), fit_v6, round(est_v6)))
    print("assoc: n=%d  fit(a,b)=%s" % (len(fpy_assoc), fit_assoc))


if __name__ == "__main__":
    main()
