"""SVOM/GRM TGF 搜索的讲图：给报告用，一张图讲一件事，字号按投影调大。

分析图（evidence/ 下那些）是给复核用的，面板多、字号小，投影看不清；这里只留主线。

用法:
    python3 scripts/plot_svom_talk.py <tgfs_v8.json> <sample_dir> <per_tgf.csv> -o <目录>
产出 talk_1_catalog.png / talk_2_criteria.png / talk_3_science.png / talk_4_map.png
"""
import argparse, csv, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.optimize import curve_fit

plt.rcParams.update({
    "font.sans-serif": ["PingFang SC", "Arial Unicode MS"], "font.family": "sans-serif",
    "axes.unicode_minus": False, "font.size": 15, "axes.titlesize": 17, "axes.labelsize": 15,
    "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13, "lines.linewidth": 2,
})
DIRECT, RESCUE = 1e-5, 1.0
BLUE, RED, GREY = "#2b6cb0", "#c53030", "0.55"


def power_law(x, a, b):
    return a * x ** b


def load_candidates(path):
    recs = json.load(open(path))
    fa = np.array([r["signal"]["false_positive_per_year"] for r in recs])
    assoc = np.array([bool(r["lightning"].get("associated")) for r in recs])
    cov = np.array([bool(r["lightning"].get("in_coverage", True)) for r in recs])
    prob = np.array([r["lightning"].get("coincidence_probability") or 0.0 for r in recs])
    lon = np.array([r["signal"]["position"]["longitude"] for r in recs])
    lat = np.array([r["signal"]["position"]["latitude"] for r in recs])
    train = np.array([bool(r["train"].get("is_train")) for r in recs])
    return fa, assoc, cov, prob, lon, lat, train


# 曝光（svomrun8 的 searched_seconds 汇总）：覆盖内 140.2 d / 191 天，覆盖外 469.0 d / 601 天
EXP_COV, EXP_OUT = 140.2, 469.0


def fig_catalog(d, out):
    """图 1：两个时期分开画——有闪电数据的能验证，没有的只能用第一层判据。

    这两段必须分开：第二层判选（fa ≤ 1 且关联）只在有闪电数据的时期能用，
    合在一张图上会让「关联」那条线看起来偏低，其实只是时间基线短了四分之三。
    """
    fa, assoc, cov, prob, _, _, train = d
    keep = ~train
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2), sharey=True)
    lo, hi = 1e-45, 20.0
    edges = np.logspace(np.log10(lo), np.log10(hi), 80)
    centers = np.sqrt(edges[:-1] * edges[1:])
    panels = (
        (axes[0], keep & cov, EXP_COV, "有闪电数据：2024-06 至 12（曝光 140 天）", True),
        (axes[1], keep & ~cov, EXP_OUT, "无闪电数据：2025-01 之后（曝光 469 天）", False),
    )
    for ax, m, exp, title, show_assoc in panels:
        n_all, _ = np.histogram(fa[m], edges)
        y = n_all / exp
        ax.step(centers, np.where(n_all > 0, y, np.nan), where="mid", color=BLUE, lw=2.2, label="全部候选")
        if show_assoc:
            n_as, _ = np.histogram(fa[m & assoc], edges)
            ax.step(centers, np.where(n_as > 0, n_as / exp, np.nan), where="mid", color=RED, lw=2.2, label="关联到闪电")
        f = (centers < 1e-8) & (n_all > 0)
        par, _ = curve_fit(power_law, centers[f], y[f], p0=(1, 0.05), maxfev=40000)
        x = np.logspace(np.log10(lo), 0, 200)
        ax.plot(x, power_law(x, *par), ls="--", color=GREY, lw=1.8, label="TGF 段幂律 指数 %.3f" % par[1])
        ax.axvspan(lo, DIRECT, facecolor="#2f855a", alpha=0.10, zorder=-2)
        ax.axvspan(DIRECT, RESCUE, facecolor="#dd6b20", alpha=0.13 if show_assoc else 0.04, zorder=-2)
        ax.axvspan(RESCUE, hi, facecolor="#c53030", alpha=0.10, zorder=-2)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(hi, lo); ax.set_ylim(2e-3, 3e2)
        ax.set_xlabel("虚警率 fa")
        ax.set_title(title, pad=10)
        ax.legend(loc="upper right", framealpha=0.9, fontsize=12)
    n_d_cov = int(((fa <= DIRECT) & keep & cov).sum()); n_r_cov = int(((fa > DIRECT) & (fa <= RESCUE) & assoc & keep & cov).sum())
    n_d_out = int(((fa <= DIRECT) & keep & ~cov).sum())
    axes[0].text(1e-24, 4e1, "直接接受 %d" % n_d_cov, ha="center", fontsize=16, color="#22543d")
    axes[0].text(2e-3, 4e1, "闪电救回 %d" % n_r_cov, ha="center", fontsize=16, color="#7b341e")
    axes[1].text(1e-24, 4e1, "直接接受 %d" % n_d_out, ha="center", fontsize=16, color="#22543d")
    axes[1].text(2e-3, 4e1, "没有闪电数据\n救不回来", ha="center", va="center", fontsize=15, color="#742a2a")
    axes[0].set_ylabel("候选数 / 曝光天")
    fig.suptitle("SVOM/GRM 目录 950 个 = 有闪电的 191 天里 %d + 没闪电的 601 天里 %d" % (n_d_cov + n_r_cov, n_d_out), fontsize=18)
    fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(out, dpi=160); print("wrote", out)


def fig_criteria(d, out):
    """图 2：分段闪电关联率——两条判选边界的依据。"""
    fa, assoc, cov, prob, _, _, train = d
    keep = ~train & cov
    edges = np.array([1e-30, 1e-15, 1e-10, 1e-7, 1e-5, 1e-3, 1e-1, 1.0, 20.0])
    x, rate, err, chance, ns = [], [], [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = keep & (fa >= a) & (fa < b)
        n = int(m.sum())
        if n < 5: continue
        k = int((assoc & m).sum())
        x.append(np.sqrt(a * b)); rate.append(100 * k / n); err.append(100 * np.sqrt(max(k, 1)) / n)
        chance.append(100 * prob[m].sum() / n); ns.append(n)
    x = np.array(x)
    fig, ax = plt.subplots(figsize=(11, 6.2))
    ax.axvspan(1e-62, DIRECT, facecolor="#2f855a", alpha=0.10, zorder=-2)
    ax.axvspan(DIRECT, RESCUE, facecolor="#dd6b20", alpha=0.13, zorder=-2)
    ax.axvspan(RESCUE, 20, facecolor="#c53030", alpha=0.10, zorder=-2)
    ax.errorbar(x, rate, yerr=err, fmt="o-", ms=9, lw=2.4, color=RED, capsize=4, label="实测闪电关联率")
    ax.plot(x, chance, ls=":", lw=2.4, color=GREY, label="偶然关联的期望")
    for xi, ri, ni in zip(x, rate, ns):
        ax.annotate("n=%d" % ni, (xi, ri), textcoords="offset points", xytext=(0, -26), ha="center", fontsize=12, color="0.35")
    ax.axvline(DIRECT, color="k", ls="--", lw=1.5); ax.axvline(RESCUE, color="k", ls="--", lw=1.5)
    ax.set_xscale("log"); ax.set_xlim(20, 1e-26); ax.set_ylim(-3, 74)
    ax.set_xlabel("虚警率 fa")
    ax.set_ylabel("闪电关联率 (%)")
    ax.set_title("判选的两条线来自数据：1e-5 处关联率陡升，fa > 1 掉回偶然", pad=12)
    ax.legend(loc="upper right")
    fig.tight_layout(); fig.savefig(out, dpi=160); print("wrote", out)


def fig_science(sample_dir, per_tgf, out):
    """图 3：证实 TGF 的时间结构与能谱。事例载入与脉冲拟合复用分析脚本，口径一致。"""
    import importlib.util
    spec = importlib.util.spec_from_file_location("tgfsample", os.path.join(os.path.dirname(__file__), "plot_svom_tgf_sample.py"))
    S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
    events, meta, bkg_spec, _bkg_ac1, eb = S.load(sample_dir)
    rows = list(csv.DictReader(open(per_tgf)))
    t50 = np.array([float(r["T50_us"]) for r in rows])
    t90 = np.array([float(r["T90_us"]) for r in rows])

    stack = []; core = np.zeros(259); bkg_scaled = np.zeros(259); total_width = 0.0
    for i, e in events.items():
        m = meta.get(i)
        if m is None: continue
        live = float(m["bkg_live_s"]); rate = float(m["bkg_counts"]) / live
        good = (e["evt"] == 0) & (e["anti"] == 0) & (e["pi"] >= S.PI_LO) & (e["pi"] < S.PI_HI)
        t = e["dt"][good]
        win = t[np.abs(t) <= S.ANALYSIS_HALF]
        fit = S.fit_pulse(win, rate)
        if fit is None: continue
        _, mu, sig = fit
        stack.append(win - mu)
        w = (e["dt"] >= mu - 2 * sig) & (e["dt"] <= mu + 2 * sig) & (e["evt"] == 0) & (e["anti"] == 0)
        core += np.bincount(np.clip(e["pi"][w], 0, 258), minlength=259)
        width = 4 * sig
        total_width += width
        b_live, b_spec = bkg_spec[i]
        bkg_scaled += b_spec * (width / b_live)
    stack = np.concatenate(stack) * 1e6      # µs

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.0))
    ax = axes[0]
    ax.hist(stack, bins=np.arange(-1000, 1001, 25), color=RED, alpha=0.85)
    ax.set_xlabel("相对脉冲中心的时间 (µs)"); ax.set_ylabel("计数 / 25 µs")
    ax.set_title("(a) %d 个证实 TGF 的叠加光变" % len(rows), pad=10)
    ax.text(0.97, 0.93, "T50 中位 %.0f µs\nT90 中位 %.0f µs" % (np.median(t50), np.median(t90)),
            transform=ax.transAxes, ha="right", va="top", fontsize=15,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))

    ax = axes[1]
    emin, emax = eb[:, 0], eb[:, 1]          # 分析脚本返回的是 (259, 2) 数组，不是字典
    net = core - bkg_scaled
    wk = emax - emin
    sel = (emin >= 42) & (emax <= 9000) & (net > 0) & (wk > 0)
    e_ = np.sqrt(emin * emax)[sel]
    y = (net / wk / max(total_width, 1e-9))[sel]
    yerr = (np.sqrt(np.maximum(core, 1)) / wk / max(total_width, 1e-9))[sel]
    ax.errorbar(e_, y, yerr=yerr, fmt="o", ms=5, color=RED, lw=1.2, capsize=0)
    # 斜率取分析脚本里逐道泊松（Cash）似然、模型在每道积分的拟合结果，这里只拟合归一化。
    # 逐点最小二乘加几何中心当能量会把指数带偏 0.1–0.2，讲图上不该显示那个数。
    for lo, hi, ls, idx, err, lab in ((42, 640, "-", -1.01, 0.04, "42–640 keV"),
                                      (640, 8041, "--", -0.90, 0.07, "640–8041 keV")):
        m = (e_ >= lo) & (e_ <= hi) & (y > 0)
        if m.sum() < 5: continue
        amp, _ = curve_fit(lambda x, a: a * x ** idx, e_[m], y[m], sigma=yerr[m], absolute_sigma=True, p0=(1e3,), maxfev=40000)
        xs = np.logspace(np.log10(lo), np.log10(hi), 50)
        ax.plot(xs, amp[0] * xs ** idx, ls=ls, color="k", lw=2.2, label="%s：计数谱指数 %.2f ± %.2f" % (lab, idx, err))
    ax.axvline(640, color=GREY, ls=":", lw=1.6)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("能道标称能量 (keV)"); ax.set_ylabel("计数 / keV / s")
    # 必须写明未解卷积：响应矩阵不可逆，这里画的是计数谱，横轴是 EBOUNDS 的标称能量，
    # 拟合出的指数不是光子谱指数。要给光子谱得用响应做正向折叠（卡在姿态四元数约定）。
    ax.set_title("(b) 叠加计数谱（未解卷积）：到 8 MeV 无截断", pad=10)
    ax.legend(loc="lower left")
    fig.suptitle("SVOM/GRM 闪电证实的 TGF：亚毫秒、硬谱（计数谱，未解卷积）", fontsize=18)
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(out, dpi=160); print("wrote", out)


def fig_map(d, per_tgf, out):
    """图 4：两个时期的地理分布分开画。

    上：有闪电数据的 191 天，红星是证实的 TGF；下：无闪电数据的 601 天，只能给出显著候选。
    两幅落在同样的三大雷暴区——这就是「覆盖外与覆盖内是同一人群」在地图上的样子。
    """
    fa, assoc, cov, prob, lon, lat, train = d
    sig = ~train & (fa <= DIRECT)
    conf = sig & cov & assoc
    unconf = sig & cov & ~assoc
    outside = sig & ~cov
    span = np.ceil(np.abs(lat[sig]).max()) + 3
    fig, axes = plt.subplots(2, 1, figsize=(13.5, 9.6), subplot_kw=dict(projection=ccrs.PlateCarree()))
    for ax in axes:
        ax.set_extent([-180, 180, -span, span], crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, facecolor="0.94")
        ax.add_feature(cfeature.COASTLINE, lw=0.5, edgecolor="0.45")
        ax.gridlines(draw_labels=False, lw=0.3, color="0.9")
    ax = axes[0]
    ax.scatter(lon[unconf], lat[unconf], s=34, c="0.72", lw=0.4, edgecolor="0.45",
               transform=ccrs.PlateCarree(), zorder=4, label="未关联 (%d)" % unconf.sum())
    ax.scatter(lon[conf], lat[conf], s=110, c=RED, marker="*", lw=0.5, edgecolor="k",
               transform=ccrs.PlateCarree(), zorder=5, label="闪电证实的 TGF (%d)" % conf.sum())
    ax.set_title("有闪电数据：2024-06 至 12（191 天）——%d 个显著候选，%d 个被闪电证实" % (int(conf.sum() + unconf.sum()), int(conf.sum())), pad=8)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.28), ncol=2, frameon=False)
    ax = axes[1]
    ax.scatter(lon[outside], lat[outside], s=34, c="#2b6cb0", alpha=0.65, lw=0.3, edgecolor="0.35",
               transform=ccrs.PlateCarree(), zorder=4, label="显著候选，无闪电数据可验证 (%d)" % outside.sum())
    ax.set_title("无闪电数据：2025-01 之后（601 天）——%d 个显著候选，落在同样的雷暴区" % int(outside.sum()), pad=8)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.28), ncol=1, frameon=False)
    fig.suptitle("闪电证实的 TGF 与无法验证的候选落在同样的三大雷暴区", fontsize=18)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tgfs"); ap.add_argument("sample_dir"); ap.add_argument("per_tgf")
    ap.add_argument("-o", "--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    d = load_candidates(args.tgfs)
    fig_catalog(d, os.path.join(args.outdir, "talk_1_catalog.png"))
    fig_criteria(d, os.path.join(args.outdir, "talk_2_criteria.png"))
    fig_science(args.sample_dir, args.per_tgf, os.path.join(args.outdir, "talk_3_science.png"))
    fig_map(d, args.per_tgf, os.path.join(args.outdir, "talk_4_map.png"))


if __name__ == "__main__":
    main()
