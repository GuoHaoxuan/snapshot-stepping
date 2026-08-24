#!/usr/bin/env python3
"""fp-distribution in the exact format of the paper figure (fp-distribution3.py),
with all three catalog generations overlaid: the paper-era search (old
saturation heuristics incl. the 10-in-10s cluster veto), v3 (rebuilt search,
before the hour-exclusion guards: flood and duplicated hours included) and v4
(bad hours excluded at the source). All restricted to < 2025-01-01.

Everything follows the paper's conventions: 20x7 cm, LaTeX Computer Modern,
log-log with the x axis reversed (weak significance on the left), 200 log bins
from the data minimum up to fpy=20, two power-law fits to the all-candidates
curve (weak side 1e-4..20, significant side 1e-50..1e-8), the associated fit
(1e-50..1e-2), and the shaded extrapolation ("Estimated TGFs") integrated out
to fpy=3.16e11.

The lightning curves come from the v4 WWLLN association (tgfs_v4.npz);
the paper-era and v3 all-candidate curves are overlaid for comparison.

usage: plot_fp_distribution_v4_paperfmt.py <fpy_paper.npz> <fpy_v3v4_pre25.npz> <tgfs_v4.npz> [-o out.pdf]
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def power_law(x, a, b):
    return a * x**b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paper_npz")
    ap.add_argument("new_npz")
    ap.add_argument("tgfs_npz")
    ap.add_argument("-o", "--out", default="fp-distribution-v4.pdf")
    args = ap.parse_args()

    dp = np.load(args.paper_npz)
    dn = np.load(args.new_npz)
    dt = np.load(args.tgfs_npz)
    pre25 = dt["day"] < 20250101
    fpy_all = dt["fpy"][pre25]
    fpy_v3 = dn["v3"]
    fpy_paper = dp["fpy"]
    fpy_assoc = dt["fpy"][pre25 & dt["assoc"]]
    assoc_pairs = list(zip(dt["fpy"][pre25], dt["coinc"][pre25]))

    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern"],
            "text.latex.preamble": "\\usepackage{amsmath}",
            "lines.linewidth": 1,
        }
    )
    cm = 1 / 2.54
    plt.figure(figsize=(20 * cm, 7 * cm), dpi=1200)

    positive_min = min(
        fpy_all[fpy_all > 0].min(),
        fpy_v3[fpy_v3 > 0].min(),
        fpy_paper[fpy_paper > 0].min(),
        fpy_assoc[fpy_assoc > 0].min(),
    )
    min_fp = float(positive_min)
    # v3 里 flood 造出的 fpy 下溢成 0，压到最低 bin
    n_zero_v3 = int((fpy_v3 == 0).sum())
    fpy_v3 = np.clip(fpy_v3, min_fp, None)
    max_fp = 3.16e11
    bins = 200
    fp_bins = np.logspace(np.log10(min_fp), np.log10(20), bins + 1)
    centers = (fp_bins[:-1] + fp_bins[1:]) / 2

    # 三代目录的全部候选
    plt.hist(fpy_v3, bins=fp_bins, histtype="step", edgecolor="#BBBBBB", lw=0.8)
    plt.hist(fpy_paper, bins=fp_bins, histtype="step", edgecolor="k", lw=0.7, alpha=0.75)
    n_all, _, _ = plt.hist(fpy_all, bins=fp_bins, histtype="step", edgecolor="C0", lw=1.1)

    # 误关联（v4 coincidence_probability 按 bin 求和）
    misassociated_count = np.zeros_like(fp_bins[:-1])
    for fp_year, prob in assoc_pairs:
        index = np.digitize(fp_year, fp_bins) - 1
        if 0 <= index < len(misassociated_count):
            misassociated_count[index] += prob
    plt.stairs(misassociated_count, fp_bins, fill=False, edgecolor="C1")

    # 闪电关联（v4）
    n_associated, _, _ = plt.hist(
        fpy_assoc, bins=fp_bins, histtype="step", edgecolor="C2", alpha=0.5
    )
    plt.stairs(n_associated - misassociated_count, fp_bins, fill=False, edgecolor="C2")

    # 拟合：弱显著侧
    condition = (centers > 1e-4) & (centers < 20)
    fit_left, _ = curve_fit(power_law, centers[condition], n_all[condition])
    x_fit = np.logspace(np.log10(min_fp), np.log10(max_fp), 100)
    plt.plot(x_fit, power_law(x_fit, *fit_left), color="#CCCCCC", ls="--", zorder=-1)

    # 拟合：显著侧 + Estimated TGFs 填充与外推积分
    condition = (centers < 1e-8) & (centers > 1e-50)
    fit_right, _ = curve_fit(power_law, centers[condition], n_all[condition], p0=fit_left)
    y_fit = power_law(x_fit, *fit_right)
    plt.plot(x_fit, y_fit, color="#CCCCCC", ls="--", zorder=-1)
    plt.fill_between(x_fit, y_fit, 1e-1, color="C0", alpha=0.1, zorder=-2)

    fp_bins_ext = np.logspace(
        np.log10(min_fp),
        np.log10(max_fp),
        int(
            np.round(
                (bins + 1)
                * (np.log10(max_fp) - np.log10(min_fp))
                / (np.log10(20) - np.log10(min_fp))
            )
        ),
    )
    estimated = float(
        power_law((fp_bins_ext[:-1] + fp_bins_ext[1:]) / 2, *fit_right).sum()
    )

    # 拟合：关联曲线
    condition = (centers < 1e-2) & (centers > 1e-50)
    ok = n_associated[condition] >= 0
    fit_assoc, _ = curve_fit(
        power_law, centers[condition][ok], n_associated[condition][ok], p0=fit_right
    )
    plt.plot(x_fit, power_law(x_fit, *fit_assoc), color="#CCCCCC", ls="--", zorder=-1)

    plt.ylim(0.5, 1e6)
    legend_handles = [
        mpatches.Patch(edgecolor="C0", facecolor="None", label="TGF Candidates (v4)"),
        mpatches.Patch(edgecolor="k", facecolor="None", alpha=0.75, label="$\\cdots$ Paper Search"),
        mpatches.Patch(edgecolor="#BBBBBB", facecolor="None", label="$\\cdots$ v3 (no hour guards)"),
        mpatches.Patch(
            edgecolor="C1",
            facecolor="None",
            label="Mis-associated TGF Candidates",
        ),
        plt.Line2D([0], [0], color="#CCCCCC", linestyle="--", label="Power Law Fits"),
        mpatches.Patch(
            edgecolor="C2",
            facecolor="None",
            alpha=0.5,
            label="TGF Candidates with Lightning",
        ),
        mpatches.Patch(
            edgecolor="C2", facecolor="None", label="$\\cdots$ But Mis-associated Excluded"
        ),
        mpatches.Patch(
            facecolor="C0", edgecolor="None", alpha=0.1, label="Estimated TGFs"
        ),
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
    print("fit left  (a,b) =", fit_left)
    print("fit right (a,b) =", fit_right)
    print("fit assoc (a,b) =", fit_assoc)
    print("Estimated TGFs:", estimated)
    print("v3 zeros clipped:", n_zero_v3)
    print("v4 associated (<2025):", len(fpy_assoc))
    print("expected misassociations (sum coinc):", float(np.sum(dt["coinc"][pre25])))


if __name__ == "__main__":
    main()
