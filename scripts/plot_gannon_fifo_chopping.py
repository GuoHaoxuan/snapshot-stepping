#!/usr/bin/env python3
"""Short-timescale look at the Gannon episode -- and what it actually shows.

At millisecond resolution the "microburst train" resolves into per-box FIFO
saturation cycling: each of the three HE boxes fills its FIFO in ~7 ms at
this flux, resets, sits dead for ~14 ms and repeats -- a strict ~21 ms comb
(measured per-box zero-gap spacing: median 21 ms, range 16-22 ms), with the
three boxes' phases nearly independent. The merged stream drops to exact
zero only where the three off-phases overlap (38% observed vs 33% expected
for independent phases).

  A. 10 s at 2 ms, merged -- the apparent "burst train";
  B. 2 s at 1 ms, per box -- three interleaved combs: the instrument is
     deeply saturated, and the envelope is a lower bound on the true flux
     (duty cycle ~1/3 per box).

usage: plot_gannon_fifo_chopping.py <gannon_raw20s_box.npz> -o <png>
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def utc(sec):
    m, s = int(sec) % 3600 // 60, sec % 60
    return "20:%02d:%02d" % (m, int(s))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("-o", "--out", default="gannon_fifo_chopping.png")
    args = ap.parse_args()
    d = np.load(args.npz)
    rel, box = d["rel"], d["box"]

    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, figsize=(11, 7.5), gridspec_kw=dict(hspace=0.35)
    )

    # ── A: 10 s, 2 ms, 合并流 ───────────────────────────────────────────
    lo, hi, bw = 2937.0, 2947.0, 0.002
    bins = np.arange(lo, hi + bw / 2, bw)
    h, _ = np.histogram(rel, bins=bins)
    ax_a.plot(bins[:-1] + bw / 2, h / bw, lw=0.4, color="crimson")
    ax_a.axvspan(2940.9, 2942.9, color="k", alpha=0.08)
    ax_a.set_xlim(lo, hi)
    ax_a.set_ylabel("counts / s (3 boxes merged)")
    ax_a.set_title(
        "10 s at 2 ms: the apparent burst train (20:48:57 - 20:49:07 UTC)"
    )
    ticks = np.arange(lo, hi + 0.5, 1.0)
    ax_a.set_xticks(ticks)
    ax_a.set_xticklabels([utc(t) for t in ticks], fontsize=8)

    # ── B: 2 s, 1 ms, 分机箱 ────────────────────────────────────────────
    lo, hi, bw = 2940.9, 2942.9, 0.001
    bins = np.arange(lo, hi + bw / 2, bw)
    x = bins[:-1] + bw / 2
    colors = ["C0", "C2", "C4"]
    for b, (name, color) in enumerate(zip("ABC", colors)):
        hb, _ = np.histogram(rel[box == b], bins=bins)
        ax_b.plot(
            x,
            hb / bw + b * 45_000,          # 逐箱竖直错开
            lw=0.5,
            color=color,
            label="box %s (+%dk offset)" % (name, b * 45),
        )
    ax_b.set_xlim(lo, hi)
    ax_b.set_yticks([0, 45_000, 90_000])
    ax_b.set_yticklabels(["0", "0", "0"])
    ax_b.set_ylabel("counts / s per box (offset)")
    ax_b.set_xlabel("UTC 20:49:00.9 + seconds")
    ax_b.set_title(
        "2 s at 1 ms, per box: three independent ~21 ms FIFO saturation combs "
        "(fill ~7 ms, dead ~14 ms)"
    )
    ax_b.legend(loc="upper right", fontsize=8, ncols=3)
    ticks = np.arange(2941.0, 2943.0, 0.25)
    ax_b.set_xticks(ticks)
    ax_b.set_xticklabels(["%.2f" % (t - 2940.9) for t in ticks], fontsize=8)

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print("saved", args.out)


if __name__ == "__main__":
    main()
