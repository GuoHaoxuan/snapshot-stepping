"""SVOM/GRM 的 ANTI_COIN 列是什么：四条证据说明它标的是星上标定源事例，不是带电粒子。

(a) AC=1 的绝对速率与总计数率无关（恒定 ~15 c/s/路）——放射源的特征，不是随粒子通量变化的反符合。
(b) AC=1 的能谱是一条 49–57 keV 的线（外加 18–31 keV 次峰），AC=0 是连续本底。
(c) 触发率与偶极磁纬无关，且 SAA 内反而更低（恒定速率被高计数率稀释）——带电粒子应当相反。
(d) 候选窗内的 AC 计数与"恒定速率 × 窗长"相符，与"本底比例 × 窗内计数"差 4 倍。

用法: python3 scripts/plot_svom_anticoin.py <arrays.npz> <features.csv> <ebounds.csv> -o <PNG>
"""
import argparse
import csv

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.sans-serif": ["PingFang SC", "Arial Unicode MS"], "font.family": "sans-serif", "axes.unicode_minus": False})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz"); ap.add_argument("features"); ap.add_argument("ebounds")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    z = np.load(args.npz)
    eb = {int(r["channel"]): (float(r["e_min"]), float(r["e_max"])) for r in csv.DictReader(open(args.ebounds))}

    fig, axes = plt.subplots(1, 4, figsize=(19, 4.3))

    # (a) 绝对速率 vs 总计数率
    ax = axes[0]
    rv = z["rate_vs_ac"]
    if len(rv):
        q = np.quantile(rv[:, 0], [0.1, 0.3, 0.5, 0.7, 0.9, 0.99])
        xs, ys, fs = [], [], []
        for lo, hi in zip([0] + list(q), list(q) + [1e9]):
            m = (rv[:, 0] >= lo) & (rv[:, 0] < hi)
            if m.sum() < 20:
                continue
            r_mean = rv[m, 0].mean()
            frac = np.average(rv[m, 1], weights=rv[m, 0])
            xs.append(r_mean); ys.append(r_mean * frac); fs.append(100 * frac)
        ax.plot(xs, ys, "o-", color="tab:red", lw=1.6, label="AC=1 绝对速率 (c/s)")
        ax2 = ax.twinx()
        ax2.plot(xs, fs, "s--", color="0.5", lw=1.2, label="AC=1 占比 (%)")
        ax2.set_ylabel("AC=1 占比 (%)", color="0.4"); ax2.set_ylim(0, max(fs) * 1.2)
        ax.set_ylim(0, max(ys) * 1.6)
        ax.legend(fontsize=8, loc="upper left"); ax2.legend(fontsize=8, loc="upper right")
    ax.set_xlabel("该秒的总计数率 (c/s，GRD1)"); ax.set_ylabel("AC=1 速率 (c/s)", color="tab:red")
    ax.set_title("(a) 绝对速率恒定，占比只是被稀释", fontsize=10)

    # (b) 能谱
    ax = axes[1]
    p0, p1 = z["pi_ac0"], z["pi_ac1"]
    ch = np.arange(len(p0))
    lastc = max(eb)
    energy = np.array([0.5 * (eb[c][0] + eb[c][1]) if c in eb else eb[lastc][1] for c in ch])
    sel = (ch >= 5) & (ch <= 120)
    ax.step(energy[sel], (p1 / p1.sum())[sel], where="mid", color="tab:red", lw=1.5, label="ANTI_COIN=1（%.2g 个）" % p1.sum())
    ax.step(energy[sel], (p0 / p0.sum())[sel], where="mid", color="0.45", lw=1.2, label="ANTI_COIN=0（%.2g 个）" % p0.sum())
    ax.axvspan(eb[28][0], eb[30][1], color="tab:red", alpha=0.12)
    ax.annotate("%.0f–%.0f keV 线\n占 AC=1 的 %.0f%%" % (eb[28][0], eb[30][1], 100 * p1[28:31].sum() / p1.sum()),
                xy=(eb[29][0], (p1 / p1.sum())[29]), xytext=(90, 0.09), fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.axvline(eb[25][0], color="k", ls="--", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("能量 (keV)；竖线=搜索能阈 ch25"); ax.set_ylabel("归一化计数")
    ax.set_title("(b) AC=1 是一条线：标定源，不是粒子", fontsize=10); ax.legend(fontsize=7, loc="lower left")

    # (c) 空间依赖
    ax = axes[2]
    n, a, e = z["mag_n"], z["mag_ac"], z["mag_edges"]
    m = n > 1000
    c = 0.5 * (e[:-1] + e[1:])
    ax.step(c[m], 100 * a[m] / n[m], where="mid", color="tab:blue", lw=1.6, label="AC=1 占比 vs 磁纬")
    sn, sa = z["saa_n"], z["saa_ac"]
    ax.axhline(100 * sa[0] / max(sn[0], 1), color="0.5", ls="--", lw=1.2, label="SAA 外 %.2f%%" % (100 * sa[0] / max(sn[0], 1)))
    ax.axhline(100 * sa[1] / max(sn[1], 1), color="tab:orange", ls=":", lw=1.6, label="SAA 内 %.2f%%（更低）" % (100 * sa[1] / max(sn[1], 1)))
    ax.set_ylim(0, 1.6)
    ax.set_xlabel("偶极磁纬 (°)"); ax.set_ylabel("ANTI_COIN=1 占比 (%)")
    ax.set_title("(c) 与磁纬无关、SAA 内更低——排除带电粒子", fontsize=10); ax.legend(fontsize=7)

    # (d) 候选窗内的两种模型
    ax = axes[3]
    rows = list(csv.DictReader(open(args.features)))
    core = np.array([float(r["acd_core"]) for r in rows]); bkg = np.array([float(r["acd_bkg"]) for r in rows])
    nn = np.array([float(r["n_core"]) for r in rows]); dur = np.array([float(r["dur_ms"]) for r in rows]) * 1e-3
    brate = np.array([float(r["rate_bkg"]) for r in rows])
    obs = (core * nn).sum(); exp_cnt = (bkg * nn).sum(); exp_time = (bkg * brate * dur).sum()
    bars = ax.bar([0, 1, 2], [obs, exp_cnt, exp_time], color=["tab:red", "0.7", "tab:green"], width=0.6)
    for b, v in zip(bars, [obs, exp_cnt, exp_time]):
        ax.text(b.get_x() + b.get_width() / 2, v + 5, "%.0f" % v, ha="center", fontsize=9)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["实测\n（%d 个候选窗）" % len(rows), "模型 A\n本底比例 × 窗内计数", "模型 B\n恒定速率 × 窗长"], fontsize=8)
    ax.set_ylabel("候选窗内的 ANTI_COIN=1 计数")
    ax.set_title("(d) 窗内 AC 与窗长成正比，与信号无关", fontsize=10)

    fig.suptitle("SVOM/GRM 的 ANTI_COIN：星上标定源事例标记（恒定 ~15 c/s/路，49–57 keV 线），不是带电粒子反符合", fontsize=12)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150, bbox_inches="tight")
    print("wrote", args.output)


if __name__ == "__main__":
    main()
