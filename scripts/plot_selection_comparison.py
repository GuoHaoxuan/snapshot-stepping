"""HXMT 与 SVOM 的判选对比：两层判据（fa ≤ 1e-5 直接接受；1e-5 < fa ≤ 1 且闪电关联）在两台仪器上是否都成立。

(a) 候选数密度随 fa 的分布，按各自曝光归一——两台仪器是否有同样平坦的 TGF 段与陡峭的本底段。
(b) 分段闪电关联率与偶然期望——第一层的边界（1e-5）和第二层的上界（1）落在哪里才对。
(c) 两层各接受多少、误救期望多少。

输入: selection_bins.csv（scripts/cluster/selection_bins.py 在集群上生成）
用法: python3 scripts/plot_selection_comparison.py selection_bins.csv -o out.png
"""
import argparse, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({"font.sans-serif": ["PingFang SC", "Arial Unicode MS"], "font.family": "sans-serif", "axes.unicode_minus": False})
COLORS = {"HXMT/HE": "tab:purple", "SVOM/GRM": "tab:cyan"}
DIRECT, RESCUE = 1e-5, 1.0


def power_law(x, a, b):
    return a * x ** b


def load(path):
    out = {}
    for r in csv.DictReader(open(path)):
        d = out.setdefault(r["instrument"], {"lo": [], "hi": [], "n": [], "a": [], "p": [], "exp": float(r["exposure_s"])})
        d["lo"].append(float(r["fa_lo"])); d["hi"].append(float(r["fa_hi"]))
        d["n"].append(float(r["n_all"])); d["a"].append(float(r["n_assoc"])); d["p"].append(float(r["sum_prob"]))
    for d in out.values():
        for k in ("lo", "hi", "n", "a", "p"): d[k] = np.array(d[k])
    return out


def regroup(d, min_count=25):
    """从高 fa 往低合并细箱，直到每组至少 min_count 个候选（关联率才有统计意义）。"""
    groups = []
    i = len(d["n"]) - 1
    while i >= 0:
        j = i; n = d["n"][i]
        while n < min_count and j > 0:
            j -= 1; n += d["n"][j]
        if n > 0: groups.append((d["lo"][j], d["hi"][i], n, d["a"][j:i + 1].sum(), d["p"][j:i + 1].sum()))
        i = j - 1
    return groups[::-1]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("bins"); ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    data = load(args.bins)
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6))

    ax = axes[0]
    for name, d in data.items():
        days = d["exp"] / 86400.0
        width = np.log10(d["hi"]) - np.log10(d["lo"])          # 每箱的 decade 宽度
        y = d["n"] / width / max(days, 1)
        centers = np.sqrt(d["lo"] * d["hi"])
        m = d["n"] > 0
        ax.step(centers[m], y[m], where="mid", color=COLORS[name], lw=1.3, label="%s（%.0f 天曝光，%d 个候选）" % (name, days, d["n"].sum()))
        # TGF 段幂律（fa < 1e-8）
        f = m & (centers < 1e-8)
        if f.sum() > 5:
            par, _ = curve_fit(power_law, centers[f], y[f], p0=(1, 0.05), maxfev=40000)
            x = np.logspace(np.log10(centers[m].min()), 0, 100)
            ax.plot(x, power_law(x, *par), ls="--", lw=1, color=COLORS[name], alpha=0.6, label="  TGF 段幂律 指数 %.3f" % par[1])
    ax.axvspan(1e-62, DIRECT, facecolor="C2", alpha=0.08, zorder=-2)
    ax.axvspan(DIRECT, RESCUE, facecolor="C1", alpha=0.10, zorder=-2)
    ax.axvspan(RESCUE, 20, facecolor="C3", alpha=0.10, zorder=-2)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(20, 1e-62)
    ax.set_xlabel("泊松假设下的年虚警期望 fa"); ax.set_ylabel("候选数 / decade / 曝光天")
    ax.set_title("(a) 两台仪器的候选分布：同样平坦的 TGF 段", fontsize=10); ax.legend(fontsize=7, loc="lower left")

    ax = axes[1]
    for name, d in data.items():
        g = regroup(d)
        x = np.array([np.sqrt(a * b) for a, b, _, _, _ in g]); n = np.array([q[2] for q in g]); a_ = np.array([q[3] for q in g]); p = np.array([q[4] for q in g])
        rate = 100 * a_ / n; err = 100 * np.sqrt(np.maximum(a_, 1)) / n
        ax.errorbar(x, rate, yerr=err, fmt="o-", ms=3.5, lw=1.2, color=COLORS[name], label="%s 实测关联率" % name)
        ax.plot(x, 100 * p / n, ls=":", lw=1.4, color=COLORS[name], label="%s 偶然期望" % name)
    ax.axvline(DIRECT, color="k", ls="--", lw=1); ax.axvline(RESCUE, color="k", ls="--", lw=1)
    ax.axvspan(1e-62, DIRECT, facecolor="C2", alpha=0.08, zorder=-2)
    ax.axvspan(DIRECT, RESCUE, facecolor="C1", alpha=0.10, zorder=-2)
    ax.axvspan(RESCUE, 20, facecolor="C3", alpha=0.10, zorder=-2)
    ax.set_xscale("log"); ax.set_xlim(20, 1e-62); ax.set_ylim(0, None)
    ax.set_xlabel("泊松假设下的年虚警期望 fa"); ax.set_ylabel("闪电关联率 (%)")
    ax.set_title("(b) 关联率：第二层（1e-5–1）捞的是真 TGF，fa > 1 掉回偶然", fontsize=10); ax.legend(fontsize=7)

    ax = axes[2]
    names = list(data); xs = np.arange(len(names)); w = 0.34
    for k, (lab, lo, hi, need_assoc, color) in enumerate((("第一层 fa ≤ 1e-5", 0.0, DIRECT, False, "C2"), ("第二层 1e-5 < fa ≤ 1 且关联", DIRECT, RESCUE, True, "C1"))):
        vals = []; false_ = []
        for name in names:
            d = data[name]; m = (d["lo"] >= lo) & (d["hi"] <= hi) if lo else (d["hi"] <= hi)
            vals.append((d["a"] if need_assoc else d["n"])[m].sum())
            false_.append(d["p"][m].sum() if need_assoc else np.nan)
        b = ax.bar(xs + (k - 0.5) * w, vals, w, color=color, alpha=0.85, label=lab)
        for i, v in enumerate(vals):
            note = "%d" % v + ("\n(误救 %.1f)" % false_[i] if need_assoc else "")
            ax.text(xs[i] + (k - 0.5) * w, v * 1.05, note, ha="center", fontsize=8)
    ax.set_xticks(xs); ax.set_xticklabels(names); ax.set_yscale("log"); ax.set_ylabel("接受的候选数")
    ax.set_ylim(1, 2e5)
    ax.set_title("(c) 两层各接受多少（覆盖内、已去列车）", fontsize=10); ax.legend(fontsize=8, loc="upper right")

    fig.suptitle("目录判选在 HXMT 与 SVOM 上的对比：fa ≤ 1e-5 直接接受，1e-5 < fa ≤ 1 且闪电关联者接受", fontsize=12)
    fig.tight_layout(); fig.savefig(args.output, dpi=150, bbox_inches="tight"); print("wrote", args.output)
    for name, d in data.items():
        m1 = d["hi"] <= DIRECT; m2 = (d["lo"] >= DIRECT) & (d["hi"] <= RESCUE); m3 = d["lo"] >= RESCUE
        for lab, m in (("fa<=1e-5", m1), ("1e-5..1", m2), ("fa>1", m3)):
            print("  %-9s %-9s 候选 %6d 关联 %5d (%4.1f%%) 偶然期望 %6.2f" % (name, lab, d["n"][m].sum(), d["a"][m].sum(), 100 * d["a"][m].sum() / max(d["n"][m].sum(), 1), d["p"][m].sum()))


if __name__ == "__main__":
    main()
