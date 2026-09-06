"""SVOM/GRM 候选与 WWLLN 闪电的关联，以及两个对照组。

用法:
    python3 scripts/plot_svom_association.py <features.csv> <assoc.csv> -o <PNG>
"""
import argparse, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature

plt.rcParams.update({
    "font.sans-serif": ["PingFang SC", "Arial Unicode MS"],
    "font.family": "sans-serif",
    "axes.unicode_minus": False,
})

RATE_HIGH = 2e4
PHASE_LO, PHASE_HI = 0.500, 0.510


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("features")
    ap.add_argument("assoc")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    feat = {r["start"][:26]: r for r in csv.DictReader(open(args.features))}
    rows = [r for r in csv.DictReader(open(args.assoc)) if r["start"][:26] in feat]
    for r in rows:
        r.update(feat[r["start"][:26]])
    col = lambda k, t=float: np.array([t(r[k]) for r in rows])

    cov = col("in_coverage", int).astype(bool)
    asc = col("associated", int).astype(bool)
    prob = np.array([float(r["coincidence_probability"]) if r["coincidence_probability"]
                     else np.nan for r in rows])
    lon, lat = col("lon"), col("lat")
    phase, rate = col("phase"), col("rate_bkg")
    high = rate > RATE_HIGH
    onehz = (~high) & (phase >= PHASE_LO) & (phase < PHASE_HI)
    tgf = ~(high | onehz)

    groups = [(tgf, "TGF 候选"), (onehz, "1 Hz 假信号\n（对照）"), (high, "高本底\n（对照）")]

    # 地图裁到纬度带后是一条窄带，单独占一行才不会被挤小
    fig = plt.figure(figsize=(13, 8.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.45, wspace=0.24)

    # (a) 覆盖内候选的地理分布，关联上的单独标出
    ax = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
    span = np.ceil(np.abs(lat).max()) + 4
    ax.set_extent([-180, 180, -span, span], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="0.95")
    ax.add_feature(cfeature.COASTLINE, lw=0.4, edgecolor="0.4")
    m = tgf & cov
    ax.scatter(lon[m & ~asc], lat[m & ~asc], s=26, c="0.72", lw=0.3, edgecolor="0.4",
               transform=ccrs.PlateCarree(), zorder=4,
               label="未关联 (%d)" % (m & ~asc).sum())
    ax.scatter(lon[m & asc], lat[m & asc], s=52, c="crimson", marker="*",
               lw=0.4, edgecolor="k", transform=ccrs.PlateCarree(), zorder=5,
               label="关联到闪电 (%d)" % (m & asc).sum())
    ax.set_title("(a) WWLLN 覆盖内的 %d 个 TGF 候选\n（库止于 2024-12-31，"
                 "全部 %d 个里只有这些查得了）" % (m.sum(), tgf.sum()), fontsize=10)
    ax.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=2,
              frameon=False)
    ax.set_anchor("C")

    # (b) 实测关联数 vs 偶然期望
    ax = fig.add_subplot(gs[1, 0])
    x = np.arange(len(groups))
    obs = [int((g & cov & asc).sum()) for g, _ in groups]
    exp = [float(np.nansum(prob[g & cov])) for g, _ in groups]
    ax.bar(x - 0.19, obs, 0.38, color="crimson", label="实测关联数")
    ax.bar(x + 0.19, np.maximum(exp, 1e-3), 0.38, color="0.65", label="偶然期望")
    for i, (o, e) in enumerate(zip(obs, exp)):
        ax.text(i - 0.19, max(o, 1e-3) * 1.35, str(o), ha="center", fontsize=9)
        ax.text(i + 0.19, max(e, 1e-3) * 1.35, "%.2g" % e, ha="center", fontsize=8,
                color="0.35")
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 400)
    ax.set_xticks(x)
    ax.set_xticklabels([n for _, n in groups], fontsize=9)
    ax.set_ylabel("候选数")
    ax.set_title("(b) 关联数与偶然期望", fontsize=10)
    ax.legend(fontsize=8, loc="upper right")

    # (c) 关联率
    ax = fig.add_subplot(gs[1, 1])
    rate_obs, rate_err, rate_exp = [], [], []
    for g, _ in groups:
        n = int((g & cov).sum())
        k = int((g & cov & asc).sum())
        # 对照组在源头修掉以后可能一个都不剩，别除零
        rate_obs.append(100 * k / n if n else 0.0)
        rate_err.append(100 * np.sqrt(k) / n if n else 0.0)
        rate_exp.append(100 * np.nansum(prob[g & cov]) / n if n else 0.0)
    ax.bar(x, rate_obs, 0.5, yerr=rate_err, color=["crimson", "0.6", "0.6"], capsize=4)
    ax.plot(x, rate_exp, "k_", ms=22, label="偶然期望")
    for i, (r, e, n) in enumerate(zip(rate_obs, rate_err,
                                      [int((g & cov).sum()) for g, _ in groups])):
        ax.text(i, r + e + 2.0, "%.1f%%\n(n=%d)" % (r, n), ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([n for _, n in groups], fontsize=9)
    ax.set_ylabel("关联率 (%)")
    ax.set_ylim(0, 52)
    ax.set_title("(c) 关联率", fontsize=10)
    ax.legend(fontsize=8)

    n, mu = obs[0], exp[0]
    p = stats.poisson.sf(n - 1, mu)
    fig.suptitle("SVOM/GRM 候选与 WWLLN 闪电的关联（±5 ms，800 km）\n"
                 "TGF 候选 %d/%d 关联，偶然期望 %.2f 个，$P=%.0e$；两个对照组落在期望上"
                 % (n, int((tgf & cov).sum()), mu, p), fontsize=12)
    fig.savefig(args.output, dpi=140, bbox_inches="tight")
    print("wrote", args.output)


if __name__ == "__main__":
    main()
