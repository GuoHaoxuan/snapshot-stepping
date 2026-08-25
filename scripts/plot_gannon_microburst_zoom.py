#!/usr/bin/env python3
"""Short-timescale zoom into the Gannon REP microburst train (2024-05-11 T20).

Two panels off the raw keep()-filtered event times:
  A. 10 s at 2 ms bins -- the train structure, individually resolved bursts;
  B. 2 s at 1 ms bins  -- single microbursts with their rise and decay.

usage: plot_gannon_microburst_zoom.py <gannon_raw20s.npz> -o <png>
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def utc(sec):
    m, s = int(sec) % 3600 // 60, sec % 60
    return "20:%02d:%06.3f" % (m, s) if s % 1 else "20:%02d:%02d" % (m, int(s))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("-o", "--out", default="gannon_microburst_zoom.png")
    args = ap.parse_args()
    rel = np.load(args.npz)["rel"]

    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, figsize=(11, 7), gridspec_kw=dict(hspace=0.3)
    )

    # ── A: 10 s, 2 ms ───────────────────────────────────────────────────
    lo, hi, bw = 2937.0, 2947.0, 0.002
    bins = np.arange(lo, hi + bw / 2, bw)
    h, _ = np.histogram(rel, bins=bins)
    ax_a.plot(bins[:-1] + bw / 2, h / bw, lw=0.5, color="crimson")
    ax_a.axvspan(2940.9, 2942.9, color="k", alpha=0.08)
    ax_a.set_xlim(lo, hi)
    ax_a.set_ylabel("counts / s")
    ax_a.set_title("microburst train: 10 s at 2 ms bins (20:48:57 - 20:49:07 UTC)")
    ticks = np.arange(lo, hi + 0.5, 1.0)
    ax_a.set_xticks(ticks)
    ax_a.set_xticklabels([utc(t) for t in ticks], fontsize=8)

    # ── B: 2 s, 1 ms ────────────────────────────────────────────────────
    lo, hi, bw = 2940.9, 2942.9, 0.001
    bins = np.arange(lo, hi + bw / 2, bw)
    h, _ = np.histogram(rel, bins=bins)
    ax_b.plot(bins[:-1] + bw / 2, h / bw, lw=0.6, color="k")
    ax_b.set_xlim(lo, hi)
    ax_b.set_ylabel("counts / s")
    ax_b.set_xlabel("UTC")
    ax_b.set_title("single microbursts: 2 s at 1 ms bins")
    ticks = np.arange(2941.0, 2943.0, 0.25)
    ax_b.set_xticks(ticks)
    ax_b.set_xticklabels([utc(t) for t in ticks], fontsize=8)

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print("saved", args.out)


if __name__ == "__main__":
    main()
