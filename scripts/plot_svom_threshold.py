"""SVOM/GRM 能阈核对图：(a) 证实 TGF / 未证实显著候选 / 本底的 PI 谱（keV 轴），(b) 能阈扫描——
不同能阈下证实样本与未证实样本的信噪（S/√B 中位）和"仍达显著"的比例，(c) 低道事例的成簇性。

用法: python3 scripts/plot_svom_threshold.py <spectra.csv> <thr_scan.csv> <burst.csv> <ebounds.csv> -o <PNG>
"""
import argparse, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import poisson

plt.rcParams.update({"font.sans-serif": ["PingFang SC", "Arial Unicode MS"], "font.family": "sans-serif", "axes.unicode_minus": False})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spectra"); ap.add_argument("scan"); ap.add_argument("burst"); ap.add_argument("ebounds"); ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    eb = {int(r["channel"]): (float(r["e_min"]), float(r["e_max"])) for r in csv.DictReader(open(args.ebounds))}
    rows = [r for r in csv.DictReader(l for l in open(args.spectra) if not l.startswith("#"))]
    meta = [l for l in open(args.spectra) if l.startswith("#")][0]
    bkg_s = float(meta.split("bkg_seconds=")[1].split()[0]); hour_s = float(meta.split("hour_seconds=")[1].split()[0])
    ch = np.array([int(r["channel"]) for r in rows]); a = np.array([float(r["assoc_core"]) for r in rows]); n = np.array([float(r["nonassoc_core"]) for r in rows]); b = np.array([float(r["bkg"]) for r in rows]); hall = np.array([float(r["hour_all"]) for r in rows])
    last = max(eb); eb.setdefault(last + 1, (eb[last][1], eb[last][1] * 1.1))   # 溢出道 259 不在 EBOUNDS 里
    emin = np.array([eb[c][0] for c in ch]); emax = np.array([eb[c][1] for c in ch]); width = np.maximum(emax - emin, 1e-3)
    scan = list(csv.DictReader(open(args.scan)))
    thr = [int(k[1:]) for k in scan[0] if k.startswith("S") and k[1:].isdigit()]
    assoc = np.array([r["assoc"] == "1" for r in scan]); cov = np.array([r["cov"] == "1" for r in scan])
    S = {t: np.array([float(r[f"S{t}"]) for r in scan]) for t in thr}; B = {t: np.array([float(r[f"B{t}"]) for r in scan]) for t in thr}

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    ax = axes[0]
    sel = (ch >= 10) & (ch <= 255)
    # 横轴保持能道；本底按秒归一，候选窗内谱按峰归一到本底的峰（只比形状）
    bkg_per_s = hall / max(hour_s, 1)
    scale = bkg_per_s[sel].max()
    ax.step(ch[sel], bkg_per_s[sel], where="mid", color="0.4", lw=1.2, label="整小时本底（%.0f h）" % (hour_s / 3600))
    ax.step(ch[sel], (a / max(a[sel].max(), 1) * scale)[sel], where="mid", color="crimson", lw=1.6, label="证实 TGF 窗内（%d 个事例，形状）" % a.sum())
    ax.step(ch[sel], (n / max(n[sel].max(), 1) * scale)[sel], where="mid", color="tab:blue", lw=1.2, label="未证实显著候选窗内（%d 个，形状）" % n.sum())
    for t, c in ((15, "k"), (10, "0.6"), (25, "0.6")):
        ax.axvline(t, color=c, ls="--", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("PI 能道"); ax.set_ylabel("计数 / 道 / s（本底）；候选谱按峰归一")
    ax.set_xticks([10, 15, 25, 50, 100, 200]); ax.set_xticklabels(["10", "15", "25", "50", "100", "200"])
    ax.set_title("(a) PI 谱：竖线 = ch10（硬件下限 %.0f keV）、ch15（现用 %.1f keV）、ch25（%.0f keV）" % (eb[10][0], eb[15][0], eb[25][0]), fontsize=10); ax.legend(fontsize=8)

    ax = axes[1]
    xs = [eb[t][0] for t in thr]
    for mask, name, color in ((assoc & cov, "证实 TGF", "crimson"), (~assoc & cov, "覆盖内未证实", "tab:blue"), (~cov, "覆盖外", "0.5")):
        if mask.sum() == 0: continue
        snr = [np.median(S[t][mask] / np.sqrt(np.maximum(B[t][mask], 0.05))) for t in thr]
        keep = [np.mean(poisson.sf(S[t][mask] - 1, np.maximum(B[t][mask], 1e-6)) <= 1e-8) for t in thr]
        ax.plot(xs, snr, "o-", color=color, label="%s (n=%d) S/√B 中位" % (name, mask.sum()))
        ax.plot(xs, np.array(keep) * max(snr), "s--", color=color, alpha=0.6, label="%s 仍达 p≤1e-8 的比例 × 尺度" % name)
    ax.axvline(eb[15][0], color="k", ls="--", lw=1)
    ax.set_xscale("log"); ax.set_xlabel("能阈 (keV)"); ax.set_ylabel("S/√B 中位")
    ax.set_title("(b) 能阈扫描：窗内计数 S 与本底期望 B 都只算阈上事例", fontsize=10); ax.legend(fontsize=7)

    ax = axes[2]
    brows = list(csv.DictReader(open(args.burst)))
    bands = sorted(set(r["band"] for r in brows), key=lambda s: int(s[2:].split("-")[0]))
    vals = [[float(r["frac_events_in_bursty_bins"]) for r in brows if r["band"] == bnd] for bnd in bands]
    ax.boxplot(vals, tick_labels=bands, showfliers=False)
    ax.set_yscale("log"); ax.set_ylabel("落在 ≥5 个/ms 格子里的事例占比"); ax.set_xlabel("PI 道段")
    ax.set_title("(c) 低道事例是否成簇（毛刺会抬高低道）", fontsize=10)
    fig.suptitle("SVOM/GRM 能阈核对（现用 PI ≥ 15 ≈ 22.5 keV）", fontsize=12)
    fig.tight_layout(); fig.savefig(args.output, dpi=140, bbox_inches="tight"); print("wrote", args.output)
    print("能阈  keV   证实:S中位 B中位 S/√B  p≤1e-8比例 | 未证实(覆盖内):S中位 B中位 S/√B p≤1e-8比例")
    for t in thr:
        m1 = assoc & cov; m2 = ~assoc & cov
        f = lambda m: (np.median(S[t][m]), np.median(B[t][m]), np.median(S[t][m] / np.sqrt(np.maximum(B[t][m], 0.05))), np.mean(poisson.sf(S[t][m] - 1, np.maximum(B[t][m], 1e-6)) <= 1e-8))
        print("%3d %6.1f   %5.0f %6.2f %5.1f %5.2f | %5.0f %6.2f %5.1f %5.2f" % ((t, eb[t][0]) + f(m1) + f(m2)))


if __name__ == "__main__":
    main()
