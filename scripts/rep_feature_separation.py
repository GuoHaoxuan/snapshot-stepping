#!/usr/bin/env python3
"""REP vs TGF feature separation on certified samples, v5 catalog.

Certified TGF = lightning-associated candidates (WWLLN, <=2024).
Certified REP = paper-selected candidates on the pre-2025 storm days
(>=15 selected/day, all known geomagnetic storms) that are NOT associated.
2025-09-30 (the largest storm cluster, outside WWLLN coverage) is overlaid
as a check population.

First-pass features, all from the candidate records themselves:
  bg_rate  = mean / dur   -- local background rate the search measured
                             around the candidate (1-s window), the cheap
                             proxy for the REP envelope;
  dur      = stop - start of the candidate interval;
  train    = number of other candidates (any fpy) within +-600 s;
  |lat|    = absolute geographic latitude (reference only).

usage: rep_feature_separation.py <tgfs_v5_feat.csv> -o <png>
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STORM_DAYS = {20230424, 20240428, 20240813, 20240918,
              20241009, 20241010, 20241023, 20241024}
CHECK_DAY = 20250930


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("feat_csv")
    ap.add_argument("-o", "--out", default="rep_feature_separation.png")
    args = ap.parse_args()

    day, t, dur, count, mean, fpy, assoc, lat = [], [], [], [], [], [], [], []
    with open(args.feat_csv) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(",")
            day.append(int(p[0])); t.append(float(p[1])); dur.append(float(p[2]))
            count.append(int(p[3])); mean.append(float(p[4])); fpy.append(float(p[5]))
            assoc.append(p[6] == "1"); lat.append(float(p[8]))
    day = np.array(day); t = np.array(t); dur = np.array(dur)
    count = np.array(count); mean = np.array(mean); fpy = np.array(fpy)
    assoc = np.array(assoc); lat = np.array(lat)

    sel = (fpy < 1e-5) | ((fpy >= 1e-5) & (fpy < 1.0) & assoc)
    is_storm = np.isin(day, list(STORM_DAYS))
    lab_tgf = assoc
    lab_rep = sel & is_storm & ~assoc
    lab_chk = sel & (day == CHECK_DAY)
    print(f"certified TGF: {lab_tgf.sum()}  certified REP: {lab_rep.sum()}  "
          f"check 20250930: {lab_chk.sum()}")

    # 列车密度：全候选（不限显著）在 ±600s 内的其他候选数
    order = np.argsort(t)
    ts = t[order]
    lo = np.searchsorted(ts, t - 600.0)
    hi = np.searchsorted(ts, t + 600.0)
    train = (hi - lo - 1).astype(float)

    bg_rate = mean / dur

    feats = [
        ("background rate around candidate [cts/s]", np.maximum(bg_rate, 1e-1), "log"),
        ("candidate duration [s]", np.maximum(dur, 1e-6), "log"),
        (r"candidates within $\pm$10 min", train + 1.0, "log"),
        ("|geographic latitude| [deg]", np.abs(lat), "linear"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (name, x, scale) in zip(axes.ravel(), feats):
        if scale == "log":
            bins = np.logspace(np.log10(x[lab_tgf | lab_rep | lab_chk].min()),
                               np.log10(x[lab_tgf | lab_rep | lab_chk].max()), 60)
            ax.set_xscale("log")
        else:
            bins = np.linspace(0, 90, 60)
        for m, label, style in ((lab_tgf, "TGF (assoc, %d)" % lab_tgf.sum(),
                                 dict(color="C2", lw=1.5)),
                                (lab_rep, "REP (storm, %d)" % lab_rep.sum(),
                                 dict(color="C3", lw=1.5)),
                                (lab_chk, "2025-09-30 (%d)" % lab_chk.sum(),
                                 dict(color="C1", lw=1.1, ls="--"))):
            n, _ = np.histogram(x[m], bins=bins)
            ax.stairs(n / max(m.sum(), 1), bins, label=label, **style)
        ax.set_xlabel(name)
        ax.set_ylabel("fraction per bin")
        ax.legend(fontsize=8)
        q = lambda m: np.percentile(x[m], [5, 50, 95])
        print(f"{name}")
        print(f"  TGF  p5/50/95: {q(lab_tgf)}")
        print(f"  REP  p5/50/95: {q(lab_rep)}")
        print(f"  chk  p5/50/95: {q(lab_chk)}")
        # 分离度：保 99% TGF 的单侧切割能杀掉多少 REP
        for side in ("above", "below"):
            if side == "above":
                cut = np.percentile(x[lab_tgf], 99)
                killed = (x[lab_rep] > cut).mean()
            else:
                cut = np.percentile(x[lab_tgf], 1)
                killed = (x[lab_rep] < cut).mean()
            print(f"  cut {side} TGF-99% ({cut:.3g}): kills {100 * killed:.1f}% of REP")

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print("saved", args.out)


if __name__ == "__main__":
    main()
