#!/usr/bin/env python3
"""fp-distribution: the paper-era catalog against v4, same window (< 2025-01-01).

Reproduces the paper figure's conventions (log-bin histogram of
false_positive_per_year, power-law fits, lightning-associated curve) and lays
the old catalog over the rebuilt one. The three catalogs differ by pipeline
generation, not just by data quality cuts:

  paper  -- original search: old saturation heuristics (6.9 ms packet-gap +
            the 10-in-10s continuous() cluster veto, since removed) and the
            paper-era WWLLN association (2042 associated candidates);
  v3     -- rebuilt search (FIFO-reset saturation, no cluster veto), before
            the hour-exclusion guards: flood + duplicated hours included;
  v4     -- same search with bad hours excluded at the source.

usage: plot_fp_distribution_v4_vs_paper.py <fpy_paper.npz> <fpy_v3v4_pre25.npz> -o <png>
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def power_law(x, a, b):
    return a * x**b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paper_npz")
    ap.add_argument("new_npz")
    ap.add_argument("-o", "--out", default="fp_distribution_v4_vs_paper.png")
    args = ap.parse_args()

    dp = np.load(args.paper_npz)
    dn = np.load(args.new_npz)
    paper, assoc = dp["fpy"], dp["assoc"]
    v4, v4_storm, v3 = dn["v4"], dn["v4_storm"], dn["v3"]

    floor = 1e-160
    clip = lambda a: np.clip(a, floor, None)
    bins = np.logspace(np.log10(floor), np.log10(20.0), 200)
    centers = np.sqrt(bins[:-1] * bins[1:])

    fig, ax = plt.subplots(figsize=(11.5, 6))

    curves = [
        ("paper all (old search, %d)" % len(paper), clip(paper), dict(color="k", lw=1.0)),
        ("v3 all (%d)" % len(v3), clip(v3), dict(color="0.7", lw=1.2)),
        ("v4 all (%d)" % len(v4), clip(v4), dict(color="C0", lw=1.4)),
        (
            "v4 minus storm days (%d)" % int((~v4_storm).sum()),
            clip(v4[~v4_storm]),
            dict(color="C2", lw=1.1),
        ),
        (
            "paper lightning-associated (%d)" % int((assoc > 0).sum()),
            clip(paper[assoc > 0]),
            dict(color="C3", lw=1.2),
        ),
    ]
    hists = {}
    for label, arr, style in curves:
        n, _ = np.histogram(arr, bins=bins)
        hists[label] = n
        ax.stairs(n, bins, label=label, **style)

    # 幂律拟合：论文关联曲线（显著侧, 同原脚本 1e-50..1e-2 窗）
    n_assoc = hists["paper lightning-associated (%d)" % int((assoc > 0).sum())]
    sel = (centers > 1e-50) & (centers < 1e-2)
    ok = n_assoc[sel] > 0
    pa, _ = curve_fit(power_law, centers[sel][ok], n_assoc[sel][ok], p0=(10.0, 0.05))
    x_fit = np.logspace(np.log10(floor), np.log10(20.0), 200)
    ax.plot(
        x_fit,
        power_law(x_fit, *pa),
        ls="--",
        color="C3",
        lw=0.8,
        alpha=0.7,
        zorder=-1,
        label=r"associated power law ($b=%.3f$)" % pa[1],
    )

    ax.axvline(1e-8, color="crimson", lw=0.7, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(floor / 3, 40)
    ax.set_ylim(0.7, 4e6)
    ax.set_xlabel("false positives per year")
    ax.set_ylabel("candidates per bin")
    ax.set_title("fp-distribution, all catalogs restricted to < 2025-01-01")
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print("saved", args.out)

    print("associated power-law b =", pa[1])
    for label, arr, _ in curves:
        print(
            "%-38s fpy<1e-8: %6d   fpy<1e-30: %5d"
            % (label, int((arr < 1e-8).sum()), int((arr < 1e-30).sum()))
        )


if __name__ == "__main__":
    main()
