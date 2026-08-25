#!/usr/bin/env python3
"""Light curve of the Gannon-storm REP episode (2024-05-11, HXMT/HE).

Three zoom levels of the same keep()-filtered event stream the TGF search
sees (CsI, non-Am241, channel >= 38):

  A. the whole hour at 1 s -- the precipitation pass stands on the normal
     orbital background modulation, latitude overlaid;
  B. the episode at 50 ms -- a quasi-continuous train of microbursts;
  C. 60 s at 5 ms -- individual sub-second microbursts, the structures the
     search fragments into "candidates".

usage: plot_gannon_lightcurve.py <gannon_lc.npz> -o <png>
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def utc(sec):
    return "%02d:%02d:%02d" % (20 + sec // 3600, sec % 3600 // 60, sec % 60)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("-o", "--out", default="gannon_lightcurve.png")
    args = ap.parse_args()
    d = np.load(args.npz)

    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        3, 1, figsize=(11, 10), gridspec_kw=dict(hspace=0.35)
    )

    # ── A: 整小时, 1 s ──────────────────────────────────────────────────
    x = np.arange(0, 3600) + 0.5
    ax_a.plot(x, d["hour"], lw=0.5, color="steelblue")
    ax_a.axvspan(2450, 3050, color="crimson", alpha=0.10)
    ax_a.set_ylabel("counts / s (keep-filtered)")
    ax_a.set_title("2024-05-11 T20 (Gannon storm): full hour, 1 s bins")
    ax_a.set_xlim(0, 3600)
    ticks = np.arange(0, 3601, 600)
    ax_a.set_xticks(ticks)
    ax_a.set_xticklabels([utc(s) for s in ticks])
    ax_lat = ax_a.twinx()
    ax_lat.plot(d["olat_t"], d["olat"], color="gray", lw=0.8, alpha=0.7)
    ax_lat.set_ylabel("latitude (deg)", color="gray")
    ax_lat.tick_params(axis="y", colors="gray")

    # ── B: episode, 50 ms ───────────────────────────────────────────────
    xb = np.arange(2450, 3050, 0.05)[: len(d["epi"])] + 0.025
    ax_b.plot(xb, d["epi"] / 0.05, lw=0.4, color="crimson")
    ax_b.axvspan(2900, 2960, color="k", alpha=0.08)
    ax_b.set_ylabel("counts / s")
    ax_b.set_title("the precipitation pass, 50 ms bins -- a train of microbursts")
    ax_b.set_xlim(2450, 3050)
    ticks = np.arange(2500, 3050, 100)
    ax_b.set_xticks(ticks)
    ax_b.set_xticklabels([utc(s) for s in ticks])

    # ── C: 60 s, 5 ms ───────────────────────────────────────────────────
    xc = np.arange(2900, 2960, 0.005)[: len(d["zoom"])] + 0.0025
    ax_c.plot(xc, d["zoom"] / 0.005, lw=0.4, color="k")
    ax_c.set_ylabel("counts / s")
    ax_c.set_xlabel("UTC")
    ax_c.set_title("60 s zoom, 5 ms bins -- individual sub-second microbursts")
    ax_c.set_xlim(2900, 2960)
    ticks = np.arange(2900, 2961, 10)
    ax_c.set_xticks(ticks)
    ax_c.set_xticklabels([utc(s) for s in ticks])

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print("saved", args.out)


if __name__ == "__main__":
    main()
