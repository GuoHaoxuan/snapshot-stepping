#!/usr/bin/env python3
"""Overview of the v4 significant-candidate census (fpy < 1e-8).

Three views that together answer "why do storm days carry hundreds of
significant candidates":

  A. daily counts across the mission -- flat 2017-2023 baseline, storm-time
     spikes in the 2024-25 solar maximum;
  B. |latitude| distributions -- storm-day candidates sit at the outer
     radiation belt footpoints, the baseline does not;
  C. the Gannon day resolved in time -- all 1049 candidates sit inside ONE
     7-minute southern high-latitude pass: a quasi-continuous train of
     sub-second microbursts at ~2.5 candidates/s. Storm-day counts measure
     episode intensity, not a number of independent events; every top storm
     day collapses to one or two such minutes-long episodes.

usage: plot_sig_v4_overview.py <sig_all_v4.csv> -o <png>
"""
import argparse
import collections
import datetime

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

STORM_TOP = {
    "20240511": "Gannon",
    "20250613": "",
    "20240628": "",
    "20250930": "",
}


def load(path):
    rows = []
    with open(path) as f:
        header = f.readline().strip().split(",")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split(",")
            rows.append(
                dict(
                    date=p[idx["date"]],
                    start=p[idx["start"]],
                    lat=float(p[idx["lat"]]),
                    lon=float(p[idx["lon"]]),
                    fpy=float(p[idx["false_positive_per_year"]]),
                )
            )
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("-o", "--out", default="sig_v4_overview.png")
    args = ap.parse_args()

    rows = load(args.csv)
    per_day = collections.Counter(r["date"] for r in rows)

    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        3, 1, figsize=(11, 11), gridspec_kw=dict(hspace=0.32)
    )

    # ── A: daily counts over the mission ────────────────────────────────
    days = sorted(per_day)
    dts = [datetime.datetime.strptime(d, "%Y%m%d") for d in days]
    counts = [per_day[d] for d in days]
    ax_a.vlines(dts, 0.7, counts, color="steelblue", lw=0.7, alpha=0.8)
    ax_a.set_yscale("log")
    ax_a.set_ylim(0.7, 2000)
    ax_a.set_ylabel("significant candidates / day")
    ax_a.set_title(
        "v4 census (fpy < 1e-8): %d candidates on %d days" % (len(rows), len(days))
    )
    for d, label in STORM_TOP.items():
        dt = datetime.datetime.strptime(d, "%Y%m%d")
        ax_a.annotate(
            (label + " " if label else "") + "%s (%d)" % (dt.strftime("%y-%m-%d"), per_day[d]),
            xy=(dt, per_day[d]),
            xytext=(dt, per_day[d] * 1.8),
            ha="center",
            fontsize=8,
            arrowprops=dict(arrowstyle="-", lw=0.6),
        )
    ax_a.xaxis.set_major_locator(mdates.YearLocator())
    ax_a.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    # ── B: |latitude| distributions ─────────────────────────────────────
    storm_days = {d for d, c in per_day.items() if c >= 50}
    lat_storm = [abs(r["lat"]) for r in rows if r["date"] in storm_days]
    lat_base = [
        abs(r["lat"]) for r in rows if r["date"] not in storm_days and r["date"][:4] < "2024"
    ]
    bins = np.linspace(0, 45, 46)
    ax_b.hist(
        lat_base,
        bins=bins,
        density=True,
        histtype="stepfilled",
        alpha=0.45,
        color="gray",
        label="2017-2023 baseline days (n=%d)" % len(lat_base),
    )
    ax_b.hist(
        lat_storm,
        bins=bins,
        density=True,
        histtype="step",
        lw=2,
        color="crimson",
        label="storm days, >=50/day (n=%d)" % len(lat_storm),
    )
    frac_storm = np.mean(np.array(lat_storm) > 40)
    frac_base = np.mean(np.array(lat_base) > 40)
    ax_b.set_xlabel("|latitude| (deg)")
    ax_b.set_ylabel("normalized")
    ax_b.set_title(
        "|lat|>40 deg: storm days %.0f%%  vs  baseline %.0f%% "
        "(outer-belt footpoints, HXMT inclination 43 deg)"
        % (100 * frac_storm, 100 * frac_base)
    )
    ax_b.legend(loc="upper left", fontsize=9)

    # ── C: the Gannon episode -- one 7-minute pass ──────────────────────
    gannon = [r for r in rows if r["date"] == "20240511"]
    secs = np.array(
        [
            int(r["start"][11:13]) * 3600
            + int(r["start"][14:16]) * 60
            + float(r["start"][17:].rstrip("Z"))
            for r in gannon
        ]
    )
    lats = np.array([r["lat"] for r in gannon])
    t0 = secs.min()
    rel = secs - t0
    bins = np.arange(0, rel.max() + 2, 1.0)
    ax_c.hist(rel, bins=bins, color="crimson", alpha=0.85)
    ax_c.set_xlabel(
        "seconds after %02d:%02d:%02d UTC (2024-05-11)"
        % (t0 // 3600, t0 % 3600 // 60, t0 % 60)
    )
    ax_c.set_ylabel("significant candidates / s")
    ax_c.set_title(
        "Gannon: all %d candidates inside ONE 7-min southern pass "
        "(lat %.0f..%.0f) -- a train of sub-second microbursts, not %d events"
        % (len(gannon), lats.min(), lats.max(), len(gannon))
    )
    ax_lat = ax_c.twinx()
    order = np.argsort(rel)
    ax_lat.plot(rel[order], lats[order], color="gray", lw=1.0, alpha=0.8)
    ax_lat.set_ylabel("latitude (deg)", color="gray")
    ax_lat.tick_params(axis="y", colors="gray")

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print("saved", args.out)
    # 给正文引用的数字
    print("storm days (>=50/day): %d days, %d candidates (%.0f%% of census)"
          % (len(storm_days), len(lat_storm), 100 * len(lat_storm) / len(rows)))


if __name__ == "__main__":
    main()
