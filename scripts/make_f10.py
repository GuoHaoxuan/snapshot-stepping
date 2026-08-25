#!/usr/bin/env python3
"""GRB 221009A main-pulse / full-burst light curve: event-level + engineering
D/tau vs GECAM-C, with the 1-Hz dead-time (live-fraction) correction.

Reproduces the former (uncommitted) f10 figures. Six traces:
  1. HXMT/HE observed (event-level)              -- navy solid
  2. HXMT/HE + reconstructed (event-level)       -- sky-blue solid (+ fill)
  3. HXMT/HE event-level / f_live (DT corrected) -- navy dashed
  4. Eng recovered D/tau (post-DT, ceiling-free) -- green solid
  5. Eng recovered D/(tau*f_live) (DT corrected) -- green dashed
  6. GECAM-C GRD01 LG (x SCALE)                  -- orange solid

Engineering channel (traces 4/5): the per-detector physical-event counter
Cnt_PHODet saturates at a 14-bit ceiling (16384) below the BOAT main-pulse
rate, so the processed-event count is instead recovered from the dead-time
counter as N_proc = D_sec / tau_i, with tau_i ~ 18 us calibrated per detector
on the post-flare tail (where Cnt_PHODet is not at ceiling). Dividing by the
live fraction f_live = 1 - D/L_cyc lifts the processed rate to a source rate.

Usage:
  python3 scripts/make_f10.py --mode main -o /tmp           # main-pulse zoom
  python3 scripts/make_f10.py --mode full -o /tmp           # full burst, symlog
  python3 scripts/make_f10.py --mode both -o ../paper-hxmt-221009a/figures
"""
import argparse, csv, colorsys, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).parent))
from plot_hxmt_vs_gecam import compute_time_offset, load_gecam_btime
from engineering_prediction import BOX_CODE, BOX_OFFSET, MET_CORRECTION, L_CYC_TO_SEC

CACHE = "data/cache_221009a_reconstruct.csv"
TRIGGER_MET_ASTROPY = 339945423.0
TRIGGER_MET_PY = 339945422.0
TRIGGER_UTC = "2022-10-09T13:17:00"
FIFO_CEIL = 38000          # observed-rate FIFO clip, for annotation only
GECAM_SCALE = 200.0        # x200, from the 1/200 precursor calibration (see f9)
TAIL = (300.0, 700.0)      # tau_i calibration window (s, rel to trigger)
PHO_CEILING = 16384        # 14-bit Cnt_PHODet ceiling

# ── blue family (match plot_hxmt_vs_gecam) ──
_h, _l, _s = colorsys.rgb_to_hls(*mcolors.to_rgb("C0"))
NAVY = colorsys.hls_to_rgb(_h, 0.25, _s)
SKY = colorsys.hls_to_rgb(_h, 0.58, _s)
GREEN, ORANGE = "C2", "C1"


def load_hxmt(before, after, bin_w):
    """Event-level observed and observed+reconstructed net rate on a bin grid."""
    obs_t, fill_t = [], []
    with open(CACHE) as f:
        r = csv.reader(f); next(r)
        for row in r:
            t = float(row[2]) - TRIGGER_MET_ASTROPY
            (obs_t if row[1] == "EVT" else fill_t).append(t)
    obs_t = np.array(obs_t); fill_t = np.array(fill_t)
    edges = np.arange(-before, after + bin_w, bin_w)
    x = edges[:-1]
    r_obs = np.histogram(obs_t, bins=edges)[0] / bin_w
    r_all = np.histogram(np.concatenate([obs_t, fill_t]), bins=edges)[0] / bin_w
    bkg = r_all[x < -10].mean()  # gross quiescent background (f_live ~ 1 there)
    return x, edges, r_obs - bkg, r_all - bkg, len(fill_t), bkg


def load_eng_dtau(before, after, tail=TAIL):
    """Per-second engineering recovery, summed over 18 detectors.

    Returns (t_rel, net_dtau, flive, tau_us) where t_rel is at 1-Hz cadence
    relative to the astropy trigger, net_dtau is the background-subtracted
    D/tau processed-event rate (evt/s), flive is the 18-detector mean live
    fraction, and tau_us is the per-detector tau_i in microseconds.
    """
    t_lo = TRIGGER_MET_PY - before
    t_hi = TRIGGER_MET_PY + after
    sec_dtau, sec_flive = {}, {}
    tau_us = []
    for box, code in BOX_CODE.items():
        folder = Path(f"data/1B/2022/20221009/{code}")
        matches = sorted(folder.glob(f"HXMT_1B_{code}_20221009T130000*.fits"))
        if not matches:
            print(f"  WARN: no 1B HE_Eng for box {box}", file=sys.stderr)
            continue
        fe = fits.open(matches[0], memmap=True)
        d = fe["HE_Eng"].data
        offset = d["UTC_Last_Bdc"][0] - d["sTime_Last_Bdc"][0]
        met = d["Time"].astype(float) + offset + MET_CORRECTION
        lc = d["Length_Time_Cycle"].astype(float)
        mask = (met >= t_lo) & (met <= t_hi)
        met_m = met[mask]
        lc_m = np.where(lc[mask] > 0, lc[mask], np.nan)
        trel_m = met_m - TRIGGER_MET_PY
        tail_mask = (trel_m >= tail[0]) & (trel_m <= tail[1])
        for det_local in range(6):
            g = BOX_OFFSET[box] + det_local
            pho = d[f"Cnt_PHODet_{g}"].astype(float)[mask]
            dtk = d[f"DeadTime_PHODet_{g}"].astype(float)[mask]
            d_sec = dtk * L_CYC_TO_SEC
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(pho > 0, d_sec / pho, np.nan)
                # only calibrate where Cnt_PHODet is below its 14-bit ceiling
                ratio = np.where(pho < PHO_CEILING - 1, ratio, np.nan)
            tau_i = np.nanmedian(ratio[tail_mask])
            tau_us.append(tau_i * 1e6)
            n_proc = d_sec / tau_i
            lf = 1.0 - dtk / lc_m
            for i, t in enumerate(met_m):
                k = int(round(t))
                sec_dtau[k] = sec_dtau.get(k, 0.0) + n_proc[i]
                sec_flive.setdefault(k, []).append(lf[i])
        fe.close()
    if not sec_dtau:
        return None, None, None, None
    secs = np.array(sorted(sec_dtau.keys()), dtype=float)
    dtau = np.array([sec_dtau[int(s)] for s in secs])
    flive = np.array([np.nanmean(sec_flive[int(s)]) for s in secs])
    t_rel = secs - TRIGGER_MET_ASTROPY
    bkg = dtau[t_rel < -10].mean() if (t_rel < -10).any() else 0.0  # gross
    return t_rel, dtau - bkg, flive, np.array(tau_us), bkg


def to_grid(x, t1hz, y1hz):
    """Hold a 1-Hz series onto the (1-s) bin-left-edge grid x."""
    m = {int(round(t)): v for t, v in zip(t1hz, y1hz)}
    return np.array([m.get(int(round(xx)), np.nan) for xx in x])


def robust(r):
    rv = r[np.isfinite(r)]
    if not len(rv):
        return None
    q75, q25 = np.percentile(rv, [75, 25])
    return np.median(rv), (q75 - q25) / 1.349, len(rv)


def draw_traces(ax, x, obs, allr, all_dt, dtau, dtau_dt, gecam_s, n_fill,
                fill=True, lw=1.1, show_eng=True):
    if fill:
        ax.fill_between(x, 0, np.nan_to_num(allr), step="post", color="C0",
                        alpha=0.30, edgecolor="none", zorder=1)
        ax.fill_between(x, 0, np.nan_to_num(obs), step="post", color="C0",
                        alpha=0.30, edgecolor="none", zorder=1)
    ax.step(x, obs, where="post", color=NAVY, lw=lw,
            label="HXMT/HE observed (event-level)", zorder=4)
    ax.step(x, allr, where="post", color=SKY, lw=lw,
            label="HXMT/HE + reconstructed (event-level)", zorder=5)
    ax.step(x, all_dt, where="post", color=NAVY, lw=lw, ls="--",
            label=r"HXMT/HE event-level $\div$ live-fraction (DT corrected)",
            zorder=6)
    if show_eng:
        ax.step(x, dtau, where="post", color=GREEN, lw=lw,
                label=r"Eng recovered $D/\tau$ (post-DT, ceiling-free)", zorder=3)
        ax.step(x, dtau_dt, where="post", color=GREEN, lw=lw, ls="--",
                label=r"Eng recovered $D/(\tau\cdot\bar f_{\rm live})$ (DT corrected)",
                zorder=3)
    ax.step(x, np.nan_to_num(gecam_s), where="post", color=ORANGE, lw=lw,
            label=rf"GECAM-C GRD01 LG ($\times${GECAM_SCALE:.0f})", zorder=2)


def make_figure(mode, x, obs, allr, all_dt, dtau, dtau_dt, gecam_s, n_fill,
                out_dir, bin_w, show_eng=True, show_title=True):
    is_full = (mode == "full")
    is_flare = (mode == "flare")
    eng_str = " + engineering" if show_eng else ""
    if is_full:
        xlo, xhi = -50, 700
    elif is_flare:
        xlo, xhi = 440, 610
    else:
        xlo, xhi = 170, 310
    fig, (axc, axr) = plt.subplots(
        2, 1, figsize=(14.5, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})

    draw_traces(axc, x, obs, allr, all_dt, dtau, dtau_dt, gecam_s, n_fill,
                fill=not is_full, show_eng=show_eng)

    if is_full:
        axc.set_yscale("symlog", linthresh=1e3)
        axc.set_ylim(1e2, 1e7)
        axc.set_ylabel("Net count rate (evt/s)  [symlog, linthresh=1e3]")
        title = (f"GRB 221009A full burst: event-level{eng_str} vs "
                 f"GECAM-C  [{bin_w} s bins, GECAM ×{GECAM_SCALE:.0f}]")
        # GECAM-C crosses the South Atlantic Anomaly over T0-50..+172 s; its
        # low-gain rate there is particle background, not source flux, so the
        # early orange hump is a data-quality flag, not a real disagreement.
        axc.axvspan(xlo, 172.0, color="0.55", alpha=0.10, zorder=0)
        axc.annotate("GECAM-C in SAA\n(data unusable)", xy=(55, 2.3e6),
                     ha="center", va="center", fontsize=13, color="0.30",
                     style="italic", fontweight="bold")
    elif is_flare:
        w = (x >= xlo) & (x < xhi)
        top = np.nanmax([np.nanmax(np.nan_to_num(allr[w])),
                         np.nanmax(np.nan_to_num(all_dt[w])),
                         np.nanmax(np.nan_to_num(gecam_s[w]))])
        axc.set_ylim(-0.05 * top, 1.15 * top)
        axc.set_ylabel("Net count rate (evt/s)")
        title = (f"GRB 221009A post-main-pulse flare: event-level{eng_str} vs "
                 f"GECAM-C  [{bin_w} s bins, GECAM ×{GECAM_SCALE:.0f}]")
    else:
        axc.set_ylim(-0.2e6, 7.7e6)
        axc.set_ylabel("Net count rate (evt/s)")
        title = (f"GRB 221009A main pulse: event-level{eng_str} vs "
                 f"GECAM-C  [{bin_w} s bins, GECAM ×{GECAM_SCALE:.0f}]")
        # inset: early main emission, where the event level is not yet clipped
        axin = axc.inset_axes([0.06, 0.40, 0.34, 0.55])
        draw_traces(axin, x, obs, allr, all_dt, dtau, dtau_dt, gecam_s, n_fill,
                    fill=True, lw=0.9, show_eng=show_eng)
        axin.set_xlim(175, 212)
        axin.set_ylim(-3000, 88000)
        axin.set_title("early main emission (T+175..212)", fontsize=9)
        axin.tick_params(labelsize=8)
        axin.axhline(0, color="gray", lw=0.4, ls="--")

    axc.set_xlim(xlo, xhi)
    axc.axhline(0, color="gray", lw=0.5, ls="--")
    if show_title:
        axc.set_title(title, fontweight="bold", fontsize=12)
    axc.legend(loc="upper right", fontsize=9.5, framealpha=0.92)

    # ── ratio panel ──
    # significance: inside the window, GECAM net above a few-sigma floor and the
    # event-level reconstruction present (finite, positive). The floor is set on
    # the GECAM *net* rate (cts/s, unscaled) so it is instrument-physical.
    win = (x >= xlo) & (x < xhi)
    g_net = gecam_s / GECAM_SCALE
    G_FLOOR = 50.0  # GECAM net cts/s; yields 76 significant bins on the main pulse
    # HXMT floor: 0 on the main pulse (76 bins), but on the full burst it removes
    # the soft precursor/quiescent region where GECAM x200 is large yet HXMT has
    # no burst signal (ratio -> 0), which is why the full view is exploratory.
    H_FLOOR = 5e3 if is_full else 0.0
    sig = (win & np.isfinite(g_net) & (g_net > G_FLOOR)
           & np.isfinite(allr) & (allr > H_FLOOR))
    ratio_specs = [
        (allr,    SKY,   "-",  "evt_rec / GECAM"),
        (all_dt,  NAVY,  "--", "evt_rec_DT / GECAM"),
    ]
    if show_eng:
        ratio_specs += [
            (dtau,    GREEN, "-",  r"eng $D/\tau$ / GECAM"),
            (dtau_dt, GREEN, "--", r"eng $D/(\tau\bar f)$ / GECAM"),
        ]
    ratios = []
    for num, color, ls, lbl in ratio_specs:
        with np.errstate(divide="ignore", invalid="ignore"):
            rr = np.where(sig & np.isfinite(num), num / gecam_s, np.nan)
        axr.step(x, rr, where="post", color=color, lw=1.0, ls=ls, label=lbl)
        ratios.append((lbl, rr))

    axr.axhline(1.0, color="gray", lw=0.5, ls="--")
    axr.set_ylim(0, 2.5)
    axr.set_ylabel("HXMT / GECAM")
    axr.set_xlabel(f"Time since trigger (s)  [$T_0$ = {TRIGGER_UTC} UTC]")
    axr.legend(loc="upper right", ncol=2, fontsize=8.5, framealpha=0.9)

    names = ["evt_rec      ", "evt_rec_DT   ", "eng D/tau    ", "eng D/tau/f  "]
    lines = []
    for (lbl, rr), nm in zip(ratios, names):
        st = robust(rr)
        if st:
            med, sig_, n = st
            lines.append(f"{nm}/ GECAM = {med:.2f} ± {sig_:.2f} ({n} bins)")
    if lines:
        axr.text(0.012, 0.95, "\n".join(lines), transform=axr.transAxes,
                 ha="left", va="top", fontsize=8.5, family="monospace",
                 bbox=dict(facecolor="white", alpha=0.9, edgecolor="lightgray"))

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = {"full": "f10_221009_full_burst",
            "flare": "f10_221009_flare",
            "main": "f10_with_eng"}[mode]
    for ext in ("pdf", "png"):
        out = out_dir / f"{stem}.{ext}"
        fig.savefig(out, dpi=160 if ext == "png" else None, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)
    return [(lbl, robust(rr)) for lbl, rr in ratios]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["main", "full", "flare", "both"], default="both")
    ap.add_argument("--bin", type=float, default=1.0)
    ap.add_argument("--before", type=float, default=50.0)
    ap.add_argument("--after", type=float, default=700.0)
    ap.add_argument("-o", "--out-dir", default="/tmp")
    ap.add_argument("--no-eng", action="store_true",
                    help="drop the two engineering D/tau recovery traces (declutter)")
    ap.add_argument("--no-title", action="store_true",
                    help="suppress the in-figure title (paper output)")
    args = ap.parse_args()

    print("Loading HXMT event-level cache...", file=sys.stderr)
    x, edges, obs, allr, n_fill, B_evt = load_hxmt(args.before, args.after, args.bin)

    print("Loading engineering D/tau recovery...", file=sys.stderr)
    et, dtau, flive, tau_us, B_eng = load_eng_dtau(args.before, args.after)
    print(f"  tau_i = {np.nanmean(tau_us):.1f} us mean "
          f"(range {np.nanmin(tau_us):.1f}..{np.nanmax(tau_us):.1f})", file=sys.stderr)
    print(f"  gross background: event-level {B_evt:,.0f} evt/s, "
          f"eng D/tau {B_eng:,.0f} evt/s", file=sys.stderr)

    print("Loading GECAM-C GRD01 LG...", file=sys.stderr)
    g_met, _, _ = compute_time_offset()
    _, g_rate, _ = load_gecam_btime(g_met, args.before, args.after, args.bin,
                                    "lg", True)
    g_rate = np.nan_to_num(g_rate, nan=0.0)
    gecam_s = g_rate * GECAM_SCALE

    # map 1-Hz engineering series onto the bin grid
    flive_g = to_grid(x, et, flive)
    dtau_g = to_grid(x, et, dtau)
    # Dead-time correction. Live fraction gates the GROSS rate (source +
    # background), so the correct net source rate is gross/f_live - B, i.e.
    # (net + B)/f_live - B, NOT net/f_live. The difference is the missing term
    # B*(1-f_live)/f_live, which is rate-dependent (largest at the deepest
    # peaks) and so is not absorbed by the constant GECAM scale.
    with np.errstate(divide="ignore", invalid="ignore"):
        all_dt = np.where(flive_g > 0.05, (allr + B_evt) / flive_g - B_evt, np.nan)
        dtau_dt = np.where(flive_g > 0.05, (dtau_g + B_eng) / flive_g - B_eng, np.nan)

    modes = ["main", "full"] if args.mode == "both" else [args.mode]
    for m in modes:
        print(f"\n=== mode={m} ===", file=sys.stderr)
        stats = make_figure(m, x, obs, allr, all_dt, dtau_g, dtau_dt, gecam_s,
                            n_fill, args.out_dir, args.bin,
                            show_eng=not args.no_eng, show_title=not args.no_title)
        for lbl, st in stats:
            if st:
                print(f"  {lbl:32s} median={st[0]:.2f} IQRsig={st[1]:.2f} "
                      f"n={st[2]}")


if __name__ == "__main__":
    main()
