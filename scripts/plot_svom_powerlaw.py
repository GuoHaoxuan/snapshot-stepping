"""SVOM/GRM 候选数随虚警率的分布（HXMT 论文里的 power-law 图，同一套画法），只用 WWLLN 覆盖内的时段。

四条阶梯：全部候选、误关联期望（各格偶然概率之和）、闪电关联的、关联减误关联。
三条幂律拟合：全部候选的本底段（fa > 1e-3）与 TGF 段（fa < 1e-8）、关联候选（fa < 1e-2）。
着色：fa ≤ 1e-5 直接接受；1e-5–1 仅关联者接受；1–20 拒绝。

用法: python3 scripts/plot_svom_powerlaw.py <tgfs.json> -o <PNG> [--until 2025-01-01]
"""
import argparse, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit

plt.rcParams.update({"font.sans-serif": ["PingFang SC", "Arial Unicode MS"], "font.family": "sans-serif", "axes.unicode_minus": False, "lines.linewidth": 1})


def power_law(x, a, b):
    return a * x ** b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tgfs"); ap.add_argument("-o", "--output", required=True); ap.add_argument("--until", default="2025-01-01", help="只用此日期之前（WWLLN 覆盖内）")
    args = ap.parse_args()
    data = [r for r in json.load(open(args.tgfs)) if r["signal"]["start"] < args.until and r["lightning"].get("in_coverage", True)]
    fa_all = np.array([r["signal"]["false_positive_per_year"] for r in data])
    assoc = np.array([bool(r["lightning"].get("associated")) for r in data])
    prob = np.array([r["lightning"].get("coincidence_probability") or 0.0 for r in data])
    min_fp, max_fp, bins = 1e-30, fa_all.max(), 100
    edges = np.logspace(np.log10(min_fp), np.log10(max_fp), bins + 1); centers = np.sqrt(edges[:-1] * edges[1:])
    fig = plt.figure(figsize=(20 / 2.54, 7.5 / 2.54), dpi=160)
    n_all, _, _ = plt.hist(fa_all, bins=edges, histtype="step", color="C0")
    mis = np.zeros(bins)
    idx = np.clip(np.digitize(fa_all, edges) - 1, 0, bins - 1)
    np.add.at(mis, idx, prob)
    plt.stairs(mis, edges, fill=False, edgecolor="C1")
    n_assoc, _, _ = plt.hist(fa_all[assoc], bins=edges, histtype="step", edgecolor="C2", alpha=0.5)
    plt.stairs(n_assoc - mis, edges, fill=False, edgecolor="C2")
    fits = {}
    for name, cond, y, p0 in (("all_bkg", centers > 1e-3, n_all, None), ("all_tgf", (centers < 1e-8) & (centers > 1e-30), n_all, None), ("assoc", (centers < 1e-2) & (centers > 1e-30), n_assoc, None)):
        m = cond & (y > 0)
        if m.sum() < 3: continue
        try:
            params, _ = curve_fit(power_law, centers[m], y[m], p0=p0 or (1.0, 0.1), maxfev=20000)
        except Exception as e:
            print("fit failed", name, e); continue
        fits[name] = params
        x = np.logspace(np.log10(min_fp), np.log10(max_fp), 200); plt.plot(x, power_law(x, *params), color="#AAAAAA", linestyle="--", zorder=-1)
    plt.axvspan(1e-30, 1e-5, facecolor="C2", edgecolor="None", alpha=0.1, zorder=-2)
    plt.axvspan(1e-5, 1, facecolor="C1", edgecolor="None", alpha=0.1, zorder=-2)
    plt.axvspan(1, 20, facecolor="C3", edgecolor="None", alpha=0.1, zorder=-2)
    handles = [mpatches.Patch(edgecolor="C0", facecolor="None", label="全部候选"),
               mpatches.Patch(edgecolor="C1", facecolor="None", label="误关联期望"),
               mpatches.Patch(edgecolor="C2", facecolor="None", alpha=0.5, label="闪电关联"),
               mpatches.Patch(edgecolor="C2", facecolor="None", label="关联 − 误关联"),
               plt.Line2D([0], [0], color="#AAAAAA", linestyle="--", label="幂律拟合"),
               mpatches.Patch(facecolor="C2", edgecolor="None", alpha=0.1, label="直接接受"),
               mpatches.Patch(facecolor="C1", edgecolor="None", alpha=0.1, label="仅关联者接受"),
               mpatches.Patch(facecolor="C3", edgecolor="None", alpha=0.1, label="拒绝")]
    plt.legend(handles=handles, ncols=2, loc="upper right", fontsize=7)
    plt.xlabel("泊松假设下的年虚警期望 fa"); plt.ylabel("候选数"); plt.xscale("log"); plt.yscale("log"); plt.xlim(max_fp, min_fp); plt.ylim(0.5, 1e5)
    plt.title("SVOM/GRM，WWLLN 覆盖内（%s 之前）：%d 个候选，%d 个关联" % (args.until, len(data), assoc.sum()), fontsize=9)
    plt.savefig(args.output, bbox_inches="tight"); print("wrote", args.output)
    for k, v in fits.items(): print("  fit %-8s a=%.3g b=%.4f" % (k, v[0], v[1]))
    accept = (fa_all < 1e-5) | ((fa_all < 1.0) & assoc)
    print("  接受的候选（论文判选）: %d = 直接 %d + 仅关联 %d；关联候选中 fa<=1e-5 的 %d / %d" % (accept.sum(), (fa_all < 1e-5).sum(), ((fa_all >= 1e-5) & (fa_all < 1) & assoc).sum(), (assoc & (fa_all <= 1e-5)).sum(), assoc.sum()))
    # 每个 fa 十倍程内的关联率与误关联期望
    print("  fa 区间           候选   关联  误关联期望  关联率")
    for lo, hi in ((1e-30, 1e-10), (1e-10, 1e-5), (1e-5, 1e-3), (1e-3, 1e-1), (1e-1, 1), (1, 20)):
        m = (fa_all >= lo) & (fa_all < hi)
        print("  [%7.0e, %7.0e) %6d %6d %8.2f   %5.1f%%" % (lo, hi, m.sum(), (assoc & m).sum(), prob[m].sum(), 100 * (assoc & m).sum() / max(m.sum(), 1)))


if __name__ == "__main__":
    main()
