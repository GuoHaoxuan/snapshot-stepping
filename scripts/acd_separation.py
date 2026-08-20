#!/usr/bin/env python3
"""ACD-coincidence separation between certified TGF and REP populations.

Second, physics-based discriminant to cross-validate the train-density cut:
REP electrons must cross the plastic ACD to reach NaI/CsI (fires ~always),
gammas only Compton in it at the few-percent level.  Input is the output of

    blink acd-audit sig_all_v5.csv -o sig_all_v5_acd.csv

run on a machine with 1K archive access.  Certification mirrors
rep_feature_separation.py: TGF = WWLLN-associated; REP = paper-selected,
pre-2025 storm-day, not associated; 2025-09-30 overlaid as check population.

Association labels: either an `assoc` column already present in the audited
CSV (0/1), or --assoc <csv> with `start` and `assoc` columns joined on the
exact `start` string (e.g. exported from tgfs.json).

Features:
  f_sig  = n_acd / n          -- ACD-fired fraction inside the candidate window
  excess = f_sig - n_acd_bg/n_bg  -- baseline-subtracted fraction
  f_bg   = n_acd_bg / n_bg    -- local baseline (the REP envelope itself)
  multi  = n_acd_multi / n    -- >=2-paddle fraction (chance-coincidence guard)

usage: acd_separation.py <audited_csv> [--assoc <csv>] [-o <png>]
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STORM_DAYS = {20230424, 20240428, 20240813, 20240918,
              20241009, 20241010, 20241023, 20241024}
CHECK_DAY = 20250930


def read_csv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split(",")
        rows = [line.rstrip("\n").split(",") for line in f if line.strip()]
    return header, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("audited_csv")
    ap.add_argument("--assoc", help="csv with `start` and `assoc` columns")
    ap.add_argument("-o", "--out", default="acd_separation.png")
    args = ap.parse_args()

    header, rows = read_csv(args.audited_csv)
    col = {name: i for i, name in enumerate(header)}
    for name in ("date", "start", "false_positive_per_year",
                 "n", "n_acd", "n_acd_multi", "n_bg", "n_acd_bg"):
        if name not in col:
            raise SystemExit(f"missing column `{name}` in {args.audited_csv}")

    # 审计缺 1K 小时的行是空串,读成 nan 后统一丢弃
    def f(row, name):
        value = row[col[name]]
        return float(value) if value else np.nan

    day = np.array([int(r[col["date"]]) for r in rows])
    start = np.array([r[col["start"]] for r in rows])
    fpy = np.array([f(r, "false_positive_per_year") for r in rows])
    n = np.array([f(r, "n") for r in rows])
    n_acd = np.array([f(r, "n_acd") for r in rows])
    n_multi = np.array([f(r, "n_acd_multi") for r in rows])
    n_bg = np.array([f(r, "n_bg") for r in rows])
    n_acd_bg = np.array([f(r, "n_acd_bg") for r in rows])

    if "assoc" in col:
        assoc = np.array([r[col["assoc"]] == "1" for r in rows])
    elif args.assoc:
        a_header, a_rows = read_csv(args.assoc)
        a_col = {name: i for i, name in enumerate(a_header)}
        table = {r[a_col["start"]]: r[a_col["assoc"]] == "1" for r in a_rows}
        assoc = np.array([table.get(s, False) for s in start])
        n_hit = sum(s in table for s in start)
        print(f"assoc join: {n_hit}/{len(start)} rows matched")
    else:
        raise SystemExit("no `assoc` column and no --assoc file given")

    ok = np.isfinite(n) & (n > 0) & np.isfinite(n_bg) & (n_bg > 0)
    print(f"rows: {len(rows)}  usable: {ok.sum()} "
          f"(dropped {len(rows) - ok.sum()}: unaudited or empty windows)")

    sel = (fpy < 1e-5) | ((fpy >= 1e-5) & (fpy < 1.0) & assoc)
    is_storm = np.isin(day, list(STORM_DAYS))
    lab_tgf = assoc & ok
    lab_rep = sel & is_storm & ~assoc & ok
    lab_chk = sel & (day == CHECK_DAY) & ok
    print(f"certified TGF: {lab_tgf.sum()}  certified REP: {lab_rep.sum()}  "
          f"check {CHECK_DAY}: {lab_chk.sum()}")

    f_sig = np.where(ok, n_acd / np.maximum(n, 1), np.nan)
    f_bg = np.where(ok, n_acd_bg / np.maximum(n_bg, 1), np.nan)
    feats = [
        ("ACD-fired fraction in candidate window", f_sig, "linear"),
        ("baseline-subtracted ACD fraction", f_sig - f_bg, "linear"),
        ("baseline ACD fraction (local envelope)", f_bg, "linear"),
        (r"$\geq$2-paddle fraction in window", np.where(ok, n_multi / np.maximum(n, 1), np.nan),
         "linear"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (name, x, _) in zip(axes.ravel(), feats):
        pool = x[lab_tgf | lab_rep | lab_chk]
        bins = np.linspace(np.nanmin(pool), np.nanmax(pool), 60)
        for m, label, style in ((lab_tgf, "TGF (assoc, %d)" % lab_tgf.sum(),
                                 dict(color="C2", lw=1.5)),
                                (lab_rep, "REP (storm, %d)" % lab_rep.sum(),
                                 dict(color="C3", lw=1.5)),
                                (lab_chk, "%d (%d)" % (CHECK_DAY, lab_chk.sum()),
                                 dict(color="C1", lw=1.1, ls="--"))):
            counts, _ = np.histogram(x[m], bins=bins)
            ax.stairs(counts / max(m.sum(), 1), bins, label=label, **style)
        ax.set_xlabel(name)
        ax.set_ylabel("fraction per bin")
        ax.legend(fontsize=8)
        q = lambda m: np.nanpercentile(x[m], [5, 50, 95])
        print(f"{name}")
        print(f"  TGF  p5/50/95: {q(lab_tgf)}")
        print(f"  REP  p5/50/95: {q(lab_rep)}")
        if lab_chk.sum():
            print(f"  chk  p5/50/95: {q(lab_chk)}")
        # 分离度：保 99% TGF 的单侧切割能杀掉多少 REP（阈值仍只由 TGF 定标）
        for side in ("above", "below"):
            if side == "above":
                cut = np.nanpercentile(x[lab_tgf], 99)
                killed = np.nanmean(x[lab_rep] > cut)
            else:
                cut = np.nanpercentile(x[lab_tgf], 1)
                killed = np.nanmean(x[lab_rep] < cut)
            print(f"  cut {side} TGF-99% ({cut:.3g}): kills {100 * killed:.1f}% of REP")

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print("saved", args.out)


if __name__ == "__main__":
    main()
