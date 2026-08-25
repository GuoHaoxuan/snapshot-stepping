#!/usr/bin/env python3
"""fp-distribution: the paper-era catalog against v5, same window (< 2025-01-01).

Same conventions as archive/2026-07-v4-catalog/plot_fp_distribution_v4_vs_paper.py (log-bin histogram of
false_positive_per_year, power-law fits, lightning-associated curve), with the
v5 catalog (time-base fix: stime-offset mode + freeze guard + edge UTC gate)
laid over v4 and the paper-era search. The generations differ by pipeline:

  paper  -- original search: old saturation heuristics (6.9 ms packet-gap +
            the 10-in-10s continuous() cluster veto, since removed) and the
            paper-era WWLLN association (2042 associated candidates);
  v4     -- rebuilt search (FIFO-reset saturation masks, hour-exclusion
            ledger), but 2897 hours carry a poisoned stime offset that
            displaced the saturation masks;
  v5     -- same search with the time base fixed: fake significant candidates
            from mask-displaced hours are gone (6945 -> 4248 significant).

Inputs are day,fpy,assoc,coinc CSVs (one row per candidate; day = YYYYMMDD):
  paper_catalog.csv  from /Volumes/Graphite/blink.db (paper-era catalog)
  tgfs_v4.csv        from the v4 tgfs.json (server extraction)
  tgfs_v5.csv        from the v5 tgfs.json (server extraction)

usage: plot_fp_distribution_v5_vs_paper.py <paper.csv> <v4.csv> <v5.csv> -o <png>
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def power_law(x, a, b):
    return a * x**b


def load_catalog(path):
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paper_csv")
    ap.add_argument("v5_csv")
    ap.add_argument("--v4", dest="v4_csv", default=None,
                    help="optionally overlay the v4 all-candidates curve")
    ap.add_argument("-o", "--out", default="fp_distribution_v5_vs_paper.png")
    args = ap.parse_args()

    sources = [("paper", args.paper_csv), ("v5", args.v5_csv)]
    if args.v4_csv:
        sources.insert(1, ("v4", args.v4_csv))
    cats = {}
    for tag, path in sources:
        day, fpy, assoc = load_catalog(path)
        pre25 = day < 20250101
        cats[tag] = (fpy[pre25], assoc[pre25])
        print(f"{tag}: {pre25.sum()} candidates pre-2025, "
              f"{int(assoc[pre25].sum())} associated, "
              f"{int((fpy[pre25] < 1e-8).sum())} significant")

    floor = 1e-160
    clip = lambda a: np.clip(a, floor, None)
    bins = np.logspace(np.log10(floor), np.log10(20.0), 200)
    centers = np.sqrt(bins[:-1] * bins[1:])

    fig, ax = plt.subplots(figsize=(11.5, 6))

    curves = [
        ("paper all (%d)" % len(cats["paper"][0]), clip(cats["paper"][0]),
         dict(color="k", lw=1.0)),
        ("v5 all (%d)" % len(cats["v5"][0]), clip(cats["v5"][0]),
         dict(color="C0", lw=1.4)),
        ("paper associated (%d)" % int(cats["paper"][1].sum()),
         clip(cats["paper"][0][cats["paper"][1]]), dict(color="C3", lw=1.2)),
        ("v5 associated (%d)" % int(cats["v5"][1].sum()),
         clip(cats["v5"][0][cats["v5"][1]]), dict(color="C2", lw=1.2)),
    ]
    if "v4" in cats:
        curves.insert(1, ("v4 all (%d)" % len(cats["v4"][0]), clip(cats["v4"][0]),
                          dict(color="0.65", lw=1.1)))
    hists = {}
    for label, arr, style in curves:
        n, _ = np.histogram(arr, bins=bins)
        hists[label] = n
        ax.stairs(n, bins, label=label, **style)

    # 幂律拟合：关联曲线显著侧（同 v4 脚本的 1e-50..1e-2 窗）
    x_fit = np.logspace(np.log10(floor), np.log10(20.0), 200)
    fits = {}
    for label, color in ((("paper associated (%d)" % int(cats["paper"][1].sum())), "C3"),
                         (("v5 associated (%d)" % int(cats["v5"][1].sum())), "C2")):
        n = hists[label]
        sel = (centers > 1e-50) & (centers < 1e-2)
        ok = n[sel] > 0
        pa, _ = curve_fit(power_law, centers[sel][ok], n[sel][ok], p0=(10.0, 0.05))
        fits[label] = pa
        ax.plot(x_fit, power_law(x_fit, *pa), ls="--", color=color, lw=0.8,
                alpha=0.7, zorder=-1,
                label=r"%s power law ($b=%.4f$)" % (label.split(" (")[0], pa[1]))

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

    for label, arr, _ in curves:
        print("%-28s fpy<1e-8: %6d   fpy<1e-30: %5d"
              % (label, int((arr < 1e-8).sum()), int((arr < 1e-30).sum())))


if __name__ == "__main__":
    main()
