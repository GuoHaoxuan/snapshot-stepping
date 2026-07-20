#!/usr/bin/env python3
"""Band-resolved HXMT/HE (NaI) vs INTEGRAL/IBIS-ISGRI light curves for
SGR 1935+2154 / FRB 200428, in three DEPOSITED-energy bands.

HXMT/HE: NaI events (pulse_width in [54,70], per the HXMT handbook) taken from
the 1B FIFO-gap reconstruction cache that now carries per-event channel AND
pulse_width; observed events have a hole during the FIFO gap, the reconstruction
(fillers, with jointly-recovered channel+pulse_width) fills it.  Deposited
energy comes from the CALDB detector-wise 3-piece-quadratic E-C
(hxmt_he_gain_20171030_v1.fits, averaged over the 18 NaI units).

INTEGRAL/IBIS-ISGRI: photon events with ISGRI_ENERGY (deposited energy in CdTe).

CAVEAT: NaI and CdTe redistribute an incident spectrum into deposited energy
differently, so a "deposited keV band" does NOT select identical incident
photons across the two instruments (largest at the low edge and near escape
features).  This is the standard deposited/PI-energy comparison, not an
incident-energy unfold.

Reuses the light-travel / time-alignment machinery from plot_hxmt_vs_ibis.py.

Rendered in the paper's publication style (pubstyle: full text width,
STIX fonts, paper-wide colour roles). Four-panel layout matching the
cross-satellite figures f7/f15: an all-events overview panel (no energy
selection) carrying the engineering-counter channel, over three deposited
NaI sub-bands.
IBIS shares the C_EXT2 purple used for non-GBM externals (SVOM/GRM); the
engineering channel keeps its C_ENG green, so the CHIME radio-pulse
markers are drawn grey to avoid a colour clash.

Standard command (200428 companion seed figure):
  .venv/bin/python scripts/plot_hxmt_vs_ibis_bands.py \
      -o ../paper-hxmt-saturation/figures/companion/f13_xsat_200428_bands.pdf
"""
import argparse
import os
import sys
from pathlib import Path

# C25 baseline path defaults to /tmp (wiped on reboot); point it at the
# persistent copy BEFORE any import can freeze the default, so the
# engineering channel survives a fresh session with no env setup.
os.environ.setdefault("C25_JSON", "data/hxmt_aux/per_det_25param.json")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.lines import Line2D
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).parent))
from plot_hxmt_vs_ibis import (  # noqa: E402
    compute_hxmt_light_travel, compute_integral_light_travel, load_ibis_events,
    HXMT_TRIGGER_UTC_STR, IBIS_FILE, INTEGRAL_ORBIT_FILE, HXMT_ORBIT_FILE,
)
sys.path.insert(0, str(Path(__file__).parent / "he_nai_cal"))
from nai_pha2pi import HEGainCalibration  # noqa: E402
from plot_hxmt_csi_multi import hxmt_met  # noqa: E402
from engineering_prediction import load_engineering_prediction, T_REF  # noqa: E402
import pubstyle  # noqa: E402

CACHE_PW = "data/cache_frb200428_reconstruct_3box_pw.csv"
GAINFILE = "data/hxmt_aux/hxmt_he_gain_20171030_v1.fits"
NAI_PW = (54, 70)                      # NaI pulse-width window
BANDS = [(20, 50), (50, 100), (100, 200)]  # deposited keV

# Reference peak positions (s, on our geocentric axis, T0 = 14:34:24.011 UTC).
# Mereghetti+2020 three IBIS X-ray peaks: t1,t2,t3 = 0.434/0.462/0.493 s w.r.t.
# their T0 = 14:34:24 UTC (geocentric) → subtract 0.011 s for our T0.
XRAY_PEAKS = [(0.434 - 0.011, "t1"), (0.462 - 0.011, "t2"), (0.493 - 0.011, "t3")]
# CHIME/STARE2 radio pulses: P1 = 14:34:24.4265 UTC (→ +0.4155 s), P2 = P1+28.9 ms.
RADIO_PEAKS = [(0.4155, "P1"), (0.4155 + 0.0289, "P2")]

# Engineering-counter channel inputs (200428; 1B counters at 14:00 UTC).
ENG_ORBIT = "data/hxmt_aux/HXMT_20200428T14_Orbit_FFFFFF_V1_1K.FITS"
ENG_DATE, ENG_HOUR = "20200428", "140000"
ENG_LOAD = 60.0    # ± s loaded so the 1 Hz channel has clean background samples


def channel_to_kev_lut():
    """Average (over 18 NaI dets) Normal-mode channel -> deposited keV LUT."""
    cal = HEGainCalibration(GAINFILE)
    ch = np.arange(256, dtype=float)
    E = np.vstack([cal.channel_to_energy(d, ch, obs_mode="Normal") for d in range(18)])
    return np.nanmean(E, axis=0)


def load_hxmt_nai(cache, trig_met, E_lut, light_travel):
    """NaI-selected HXMT events -> (obs_t, obs_e, fill_t, fill_e) rel. trigger."""
    obs_t, obs_e, fill_t, fill_e = [], [], [], []
    with open(cache) as f:
        next(f)
        for line in f:
            p = line.split(",")
            typ, met, ch, pw = p[1], float(p[2]), int(p[3]), int(p[4])
            if not (NAI_PW[0] <= pw <= NAI_PW[1]):
                continue
            t = met - trig_met + light_travel
            e = float(E_lut[min(max(ch, 0), 255)])
            if typ == "EVT":
                obs_t.append(t); obs_e.append(e)
            elif typ == "FILL_GAP":
                fill_t.append(t); fill_e.append(e)
    return (np.asarray(obs_t), np.asarray(obs_e),
            np.asarray(fill_t), np.asarray(fill_e))


def load_hxmt_all(cache, trig_met, light_travel):
    """All reconstructed events (no NaI / energy selection) -> (obs_t, fill_t)."""
    obs_t, fill_t = [], []
    with open(cache) as f:
        next(f)
        for line in f:
            p = line.split(",")
            typ, met = p[1], float(p[2])
            t = met - trig_met + light_travel
            if typ == "EVT":
                obs_t.append(t)
            elif typ == "FILL_GAP":
                fill_t.append(t)
    return np.asarray(obs_t), np.asarray(fill_t)


def fit_background(x, rate, bkgm, deg=1):
    """Polynomial background fit (default linear) over the off-burst windows,
    evaluated across the whole time axis — interpolates the background trend
    through the burst instead of using a flat mean."""
    coef = np.polyfit(x[bkgm], rate[bkgm], deg)
    return np.polyval(coef, x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", type=float, default=0.005)
    ap.add_argument("--bkg-deg", type=int, default=1,
                    help="background polynomial degree (1 = linear)")
    ap.add_argument("--before", type=float, default=0.3)
    ap.add_argument("--after", type=float, default=0.9)
    ap.add_argument("--bkg", type=float, nargs=4, default=[-0.3, -0.05, 0.7, 0.9])
    ap.add_argument("--scale-range", type=float, nargs=2, default=[0.35, 0.65])
    ap.add_argument("--xlim", type=float, nargs=2, default=None,
                    metavar=("T1", "T2"),
                    help="display window (s); rates/background still use full range")
    ap.add_argument("-o", "--output", default="hxmt_vs_ibis_bands.png")
    args = ap.parse_args()

    pubstyle.apply()

    trig = Time(HXMT_TRIGGER_UTC_STR, scale="utc")
    trig_met = (trig - Time("2012-01-01T00:00:00", scale="utc")).sec
    trig_ijd = (trig.tt - Time("2000-01-01T00:00:00", scale="tt")).to("day").value
    hlt = compute_hxmt_light_travel(trig_met, trig, HXMT_ORBIT_FILE)
    ilt = compute_integral_light_travel(trig_ijd, IBIS_FILE, INTEGRAL_ORBIT_FILE)
    print(f"  HXMT light-travel {hlt*1e3:+.1f} ms, INTEGRAL {ilt:+.3f} s", file=sys.stderr)

    E_lut = channel_to_kev_lut()
    obs_t, obs_e, fill_t, fill_e = load_hxmt_nai(CACHE_PW, trig_met, E_lut, hlt)
    all_t = np.concatenate([obs_t, fill_t])
    all_e = np.concatenate([obs_e, fill_e])
    # all events (no NaI / energy selection) for the overview panel
    allev_obs_t, allev_fill_t = load_hxmt_all(CACHE_PW, trig_met, hlt)
    allev_all_t = np.concatenate([allev_obs_t, allev_fill_t])
    ibis_t, ibis_e = load_ibis_events(IBIS_FILE, trig_ijd, ilt)
    print(f"  HXMT all: {len(allev_obs_t):,} obs + {len(allev_fill_t):,} fill "
          f"(overview panel);  NaI sub-bands: {len(obs_t):,} obs + "
          f"{len(fill_t):,} fill;  IBIS: {len(ibis_t):,} events", file=sys.stderr)

    # Engineering-counter prediction (1 Hz, 18-detector sum). Its trigger MET
    # must be the TAI-based HXMT MET (same clock as the 1B counters), NOT the
    # UTC-epoch trig_met used for the event axis, else the 1 Hz samples land
    # several leap-seconds off. The counters register the SAME photons hitting
    # HE, so to share the geocentric axis of the events and IBIS they get the
    # same light travel (hlt) the events do. The integer GPS/UTC-second cycle
    # edges (T0 - 0.011 s on the HXMT clock, since T0 = ...:24.011) then map to
    # ~T0 once shifted by hlt (+11.4 ms) — a coincidence of this burst (hlt
    # happens to match T0's 0.011 s sub-second part), not an alignment to T0.
    eng_t = eng_rate = net_eng = None
    try:
        ty = (np.datetime64("2020-04-28") - T_REF).astype(
            "timedelta64[D]").astype(float) / 365.25
        eng_t, eng_rate = load_engineering_prediction(
            date_str=ENG_DATE, hour_str=ENG_HOUR,
            trigger_met=hxmt_met(HXMT_TRIGGER_UTC_STR),
            before=ENG_LOAD, after=ENG_LOAD, t_years_const=ty,
            orbit_path=ENG_ORBIT)
    except Exception as exc:  # noqa: BLE001
        print(f"  WARN: engineering channel unavailable: {exc}", file=sys.stderr)
    if eng_t is not None:
        eng_t = eng_t + hlt
        ebm = (np.abs(eng_t) > 5) & (np.abs(eng_t) < ENG_LOAD)
        net_eng = eng_rate - (np.median(eng_rate[ebm]) if ebm.any() else 0.0)

    edges = np.arange(-args.before, args.after + args.bin, args.bin)
    x = edges[:-1] + args.bin / 2
    t1, t2, t3, t4 = args.bkg
    bkgm = ((x >= t1) & (x < t2)) | ((x >= t3) & (x < t4))
    s1, s2 = args.scale_range
    sm = (x >= s1) & (x < s2)

    xl = tuple(args.xlim) if args.xlim else (-args.before, args.after)
    vis = (x >= xl[0]) & (x < xl[1])

    # Top overview panel (all events, no energy selection) + three NaI
    # sub-band panels, matching the paper's cross-satellite figures (f7/f15).
    panels = [(20, 200)] + BANDS  # top (20,200) is a placeholder; see is_top
    fig, axes = plt.subplots(
        len(panels), 1, figsize=(pubstyle.FULL_W, 6.0), sharex=True,
        gridspec_kw={"hspace": 0.0})
    for i, (ax, (elo, ehi)) in enumerate(zip(axes, panels)):
        is_top = i == 0
        if is_top:  # overview: all events, no NaI / energy selection
            r_obs = np.histogram(allev_obs_t, bins=edges)[0] / args.bin
            r_all = np.histogram(allev_all_t, bins=edges)[0] / args.bin
            r_ibis = np.histogram(ibis_t, bins=edges)[0] / args.bin
            n_ibis_raw = np.histogram(ibis_t, bins=edges)[0]
        else:
            def rate(t, e):
                m = (e >= elo) & (e < ehi)
                return np.histogram(t[m], bins=edges)[0] / args.bin
            r_obs = rate(obs_t, obs_e)
            r_all = rate(all_t, all_e)
            r_ibis = rate(ibis_t, ibis_e)
            # IBIS 每 bin 原始计数 → 泊松误差(率误差 = sqrt(N)/bin)
            ibis_m = (ibis_e >= elo) & (ibis_e < ehi)
            n_ibis_raw = np.histogram(ibis_t[ibis_m], bins=edges)[0]
        ibis_err = np.sqrt(n_ibis_raw) / args.bin
        n_obs = r_obs - fit_background(x, r_obs, bkgm, args.bkg_deg)
        n_all = r_all - fit_background(x, r_all, bkgm, args.bkg_deg)
        n_ibis = r_ibis - fit_background(x, r_ibis, bkgm, args.bkg_deg)
        scale = n_all[sm].sum() / n_ibis[sm].sum() if n_ibis[sm].sum() > 0 else 1.0

        # reference-time markers behind the curves
        for tp, _ in XRAY_PEAKS:
            ax.axvline(tp, color=pubstyle.C_SAT, ls="--", lw=0.7, alpha=0.55, zorder=0)
        for rp, _ in RADIO_PEAKS:
            ax.axvline(rp, color="0.35", ls=":", lw=0.9, alpha=0.8, zorder=0)
        det = "" if is_top else "NaI "
        # HXMT observed / reconstructed (dark + light blue, filled between)
        ax.fill_between(x, n_obs, n_all, step="mid", color=pubstyle.C_RECON,
                        alpha=0.30, zorder=4)
        ax.step(x, n_obs, where="mid", color=pubstyle.C_OBS, lw=1.0,
                label=f"HXMT/HE {det}observed", zorder=6)
        ax.step(x, n_all, where="mid", color=pubstyle.C_RECON, lw=1.2,
                label=f"HXMT/HE {det}reconstructed", zorder=7)
        # INTEGRAL/IBIS (purple: external reference other than Fermi/GBM)
        ax.fill_between(x, (n_ibis - ibis_err) * scale, (n_ibis + ibis_err) * scale,
                        step="mid", color=pubstyle.C_EXT2, alpha=0.12, lw=0, zorder=2)
        ax.step(x, n_ibis * scale, where="mid", color=pubstyle.C_EXT2, lw=1.0,
                label=rf"INTEGRAL/IBIS-ISGRI $\times${scale:.1f} ($\pm\sqrt{{N}}$)",
                zorder=3)
        # engineering-counter channel: overview panel only (1 Hz, whole-burst)
        if is_top and eng_t is not None:
            ax.step(eng_t, net_eng, where="post", color=pubstyle.C_ENG, lw=1.1,
                    label=r"engineering $\widehat{S}_{\rm rec}^{\rm eng}$"
                          " (1 Hz, 18-det sum)", zorder=5)
        ax.axhline(0, color="grey", lw=0.5)
        ax.margins(x=0)
        ax.set_ylabel("net rate (counts/s)")
        tag = ("all events" if is_top else f"{elo}–{ehi} keV (deposited)")
        ax.text(0.02, 0.92, f"{tag}, {args.bin*1e3:.0f} ms bins",
                transform=ax.transAxes, fontweight="bold", va="top", fontsize=8)
        cand = [n_all[vis].max(), (n_ibis * scale)[vis].max()]
        if is_top and eng_t is not None:
            ev = (eng_t >= xl[0]) & (eng_t < xl[1])
            if ev.any():
                cand.append(net_eng[ev].max())
        ytop = max(cand) * 1.20
        ax.set_ylim(min(0, n_all[vis].min() * 1.1), ytop)
        ax.yaxis.set_major_locator(
            matplotlib.ticker.MaxNLocator(nbins=5, prune="both"))
        h, la = ax.get_legend_handles_labels()
        h += [Line2D([0], [0], color=pubstyle.C_SAT, ls="--", lw=0.7),
              Line2D([0], [0], color="0.35", ls=":", lw=0.9)]
        la += ["IBIS X-ray peaks (Mereghetti+20)", "radio pulses (CHIME)"]
        ax.legend(h, la, loc="upper right")

    axes[-1].set_xlim(*xl)
    axes[-1].set_xlabel(
        f"time since HXMT $T_0$ (s)   [$T_0$ = {HXMT_TRIGGER_UTC_STR} UTC]")
    fig.subplots_adjust(left=0.077, right=0.988, top=0.99, bottom=0.062)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output)
    print(f"wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
