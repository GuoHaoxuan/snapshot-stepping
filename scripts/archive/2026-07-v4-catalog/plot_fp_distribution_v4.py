#!/usr/bin/env python3
"""fp-distribution of the v4 catalog, against v3.

Same conventions as the original fp-distribution3.py: histogram of
false_positive_per_year in log bins, log-log axes, power-law fit
N(bin) = a * fpy^b on the weak-significance side (1e-4 .. 20) extrapolated
across the plot. What changed since that plot was made:

  * v3's "all candidates" dragged a tail of fabricated ultra-significant
    candidates down to fpy = 0 (float underflow) -- the flood + duplicated
    hours. v4 excludes those hours at the source.
  * v4 still reaches fpy ~ 1e-197, but that is the storm-day REP
    microbursts -- real signals, split out as their own curve here.
  * v4-minus-storm-days is the population the TGF power law lives in.

usage: plot_fp_distribution_v4.py <fpy_v3v4.npz> -o <png>
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
    ap.add_argument("npz")
    ap.add_argument("-o", "--out", default="fp_distribution_v4.png")
    args = ap.parse_args()

    d = np.load(args.npz)
    v4, v4_storm, v3 = d["v4"], d["v4_storm"], d["v3"]

    # v3 里 flood 造出的 fpy 会下溢成字面的 0，压到最低 bin 边上并单独交代
    n_zero_v3 = int((v3 == 0).sum())
    floor = 1e-160
    v3 = np.clip(v3, floor, None)
    v4c = np.clip(v4, floor, None)

    bins = np.logspace(np.log10(floor), np.log10(20.0), 200)
    centers = np.sqrt(bins[:-1] * bins[1:])

    fig, ax = plt.subplots(figsize=(11, 5.5))
    n_v3, _ = np.histogram(v3, bins=bins)
    n_v4, _ = np.histogram(v4c, bins=bins)
    n_v4q, _ = np.histogram(v4c[~v4_storm], bins=bins)

    ax.stairs(n_v3, bins, color="0.65", lw=1.2, label="v3 all (flood + duplicated hours in)")
    ax.stairs(n_v4, bins, color="C0", lw=1.4, label="v4 all")
    ax.stairs(
        n_v4q, bins, color="C2", lw=1.2, label="v4 minus 8 storm days (TGF-side population)"
    )

    # 幂律拟合：弱显著侧 1e-4..20（与原脚本同窗），在无风暴曲线上拟
    sel = (centers > 1e-4) & (centers < 20)
    params, _ = curve_fit(power_law, centers[sel], n_v4q[sel])
    x_fit = np.logspace(np.log10(floor), np.log10(20.0), 200)
    ax.plot(
        x_fit,
        power_law(x_fit, *params),
        ls="--",
        color="0.3",
        lw=1,
        zorder=-1,
        label=r"power law fit ($b=%.3f$), extrapolated" % params[1],
    )

    ax.axvline(1e-8, color="crimson", lw=0.8, ls=":")
    ax.text(
        1e-8,
        1e6,
        "  census cut fpy=1e-8",
        color="crimson",
        fontsize=8,
        rotation=90,
        va="top",
    )
    ax.annotate(
        "v3 flood/dup tail\n(%d candidates at fpy=0, clipped here)" % n_zero_v3,
        xy=(floor * 3, 2e3),
        fontsize=8,
        color="0.4",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(floor / 3, 40)
    ax.set_ylim(0.7, 4e6)
    ax.set_xlabel("false positives per year")
    ax.set_ylabel("candidates per bin")
    ax.set_title(
        "fp-distribution: v3 (%.2fM) vs v4 (%.2fM, excluded hours removed at source)"
        % (len(v3) / 1e6, len(v4) / 1e6)
    )
    ax.legend(loc="upper left", fontsize=9)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print("saved", args.out)
    print("power-law b (v4 non-storm, 1e-4..20) =", params[1])
    for name, arr in (("v3", v3), ("v4", v4c), ("v4-nonstorm", v4c[~v4_storm])):
        print(
            "%-12s n=%8d   fpy<1e-8: %6d   fpy<1e-50: %6d"
            % (name, len(arr), int((arr < 1e-8).sum()), int((arr < 1e-50).sum()))
        )


if __name__ == "__main__":
    main()
