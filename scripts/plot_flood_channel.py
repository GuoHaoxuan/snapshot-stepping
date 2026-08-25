#!/usr/bin/env python3
"""Raw-channel spectrum: a flood episode against a normal hour and a real burst.

The three histograms are frozen counts measured with channel_hist.py (which
decodes byte0 straight out of the CCSDS packets, no time solve). They are the
evidence behind the flood criterion in `config_guard.rs`: channels 13-19 are
normally an almost-empty gap between the folded low-energy peak and the start
of the nominal spectrum at ch20, and a flood fills that gap by three orders of
magnitude. GRB 221009A is here to show a genuinely bright burst does NOT --
its real high-energy events fold into the ch0-12 tail and never touch the gap.

usage: plot_flood_channel.py -o <png>
"""
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ch = np.arange(46)
NORM = np.array([22394,22252,21496,21635,21512,21827,21398,21090,21249,20725,19816,19207,17778,13207,6101,3136,1721,215,74,293,100203,374034,236243,227254,221705,221333,216252,206289,189491,171548,151636,136868,120896,112461,105441,103013,100187,98426,99623,100491,102673,103982,105865,110477,114480,120107],float)
FLOOD= np.array([22727,22518,22342,22554,21946,22137,21735,21724,21730,21024,20294,19554,42275,1838706,2207028,1013928,745748,507259,351947,258486,200412,172705,161639,161654,164299,171263,173191,171169,160068,147076,130465,117402,104527,97426,91865,88276,87655,87401,88164,90555,92083,94810,97235,102308,107082,113362],float)
GRB  = np.array([11003,10737,10792,10647,10681,10628,10482,10432,10329,10228,9330,9226,8615,6408,2618,1453,816,113,71,44,187727,625775,362266,367194,351274,339110,324786,311653,294737,279111,264417,254798,241253,239531,234172,234924,234270,233448,237218,242545,246474,254400,260172,268826,267576,257159],float)
def nrm(a): return a/a.sum()


ap = argparse.ArgumentParser()
ap.add_argument("-o", "--out", required=True, help="output PNG path")
args = ap.parse_args()

fig, ax = plt.subplots(figsize=(11,6))
ax.axvspan(12.5,19.5, color="orange", alpha=0.13, zorder=0)
ax.text(16, 2.6e-1, "ch 13-19: normally an\nempty gap; flood floods it", ha="center", va="top", fontsize=10, color="#b35900")
ax.step(ch, nrm(FLOOD), where="mid", color="#d62728", lw=2.3, label="FLOOD  E3  (2019-12-06 08h)")
ax.step(ch, nrm(NORM),  where="mid", color="#1f77b4", lw=1.8, label="Normal      (2019-12-04 10h)")
ax.step(ch, nrm(GRB),   where="mid", color="#2ca02c", lw=1.8, ls="--", label="GRB 221009A real bright burst (2022-10-09 13h)")
ax.axvline(19.5, color="k", lw=1, ls=":", alpha=0.6)
ax.text(20.2, 1.6e-5, "nominal fold edge ch19/20\n(normal spectrum starts at ch20)", fontsize=8.5, color="k", alpha=0.75)
ax.set_yscale("log"); ax.set_xlim(-0.5,45.5); ax.set_ylim(1e-6,3e-1)
ax.set_xlabel("raw channel (byte0)"); ax.set_ylabel("normalized fraction (log)")
ax.set_title("HXMT/HE raw-channel spectrum: flood vs normal vs real bright burst")
ax.legend(loc="upper right", fontsize=9.5); ax.grid(alpha=0.25, which="both")
fig.tight_layout(); fig.savefig(args.out, dpi=140)
print("saved", args.out)
