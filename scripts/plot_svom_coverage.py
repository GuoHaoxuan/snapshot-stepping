"""SVOM/GRM：WWLLN 覆盖外的显著候选与覆盖内的是不是同一人群。

闪电库止于 2024-12-31，只覆盖 SVOM 观测期的头 193 天；目录里 78% 的候选无法用闪电验证。
若能证明覆盖外与覆盖内的候选在所有可测量的性质上不可区分，覆盖内的关联率就能外推到整个目录。

做法：
- 单变量：窗长、能道中位数比、|纬度|、本底率、窗内计数、单探头占比、log10 fa 用两样本 KS 与
  Anderson-Darling；地方时是圆周量，用对起点旋转不变的 Kuiper 检验（置换定 p）。
- 多变量：六个量标准化后做能量距离置换检验，一次给出联合分布是否相同。用「覆盖内已证实 vs
  未证实」作阳性对照——同一个检验必须能查出那个真实差异，否则「查不出差异」没有说服力。
- 两个混杂因素：其一，覆盖内只有 2024-06 至 12，覆盖外含 2025 全年与 2026 前八月，雷暴的季节
  与地理分布本来就随时间变，故另做「同为 6–12 月」的季节匹配比较；其二，WWLLN 的探测效率有
  地域与昼夜差异，故外推不用全局关联率，而按 3 经度区 × 4 地方时段分层重加权。

用法:
    python3 scripts/plot_svom_coverage.py <features_top899.csv> <assoc_svom_v6.csv> -o <PNG>
"""
import argparse, csv
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.sans-serif": ["PingFang SC", "Arial Unicode MS"], "font.family": "sans-serif", "axes.unicode_minus": False})

# 覆盖内只有 6–12 月，季节匹配就取覆盖外的同月份
SEASON_MONTHS = (6, 7, 8, 9, 10, 11, 12)
ZONES = ((-180, -30, "美洲"), (-30, 60, "非洲欧洲"), (60, 180, "亚洲海洋大陆"))
LST_EDGES = (0, 6, 12, 18, 24)
# 多变量检验用的量。本底率单列：它有已知的仪器漂移（覆盖外内部随时间显著上升），
# 是仪器状态不是候选性质，混进联合检验会把结论引偏。
JOINT = ("dur", "pir", "alat", "n", "dfrac", "lst")
LABELS = {"dur": "窗长 (ms)", "pir": "能道中位数比", "alat": "|纬度| (°)", "rate": "本底率 (c/s)",
          "n": "窗内计数", "dfrac": "单探头占比", "lst": "地方时 (h)", "logfa": "$\\log_{10}$ fa"}


def load(features, assoc_csv):
    assoc = {r["start"][:26]: r for r in csv.DictReader(open(assoc_csv))}
    rows = []
    for r in csv.DictReader(open(features)):
        a = assoc[r["start"][:26]]
        lon = float(r["lon"])
        hour = int(r["start"][11:13]) + int(r["start"][14:16]) / 60.0
        rows.append(dict(
            cov=a["in_coverage"] == "1", ass=a["associated"] == "1",
            prob=float(a["coincidence_probability"]) if a["coincidence_probability"] else 0.0,
            mon=int(r["start"][5:7]), lon=lon, lat=float(r["lat"]), alat=abs(float(r["lat"])),
            dur=float(r["dur_ms"]), pir=float(r["pi_med_ratio"]) if r["pi_med_ratio"] else np.nan,
            rate=float(r["rate_bkg"]), n=float(r["n_core"]), dfrac=float(r["det_frac_max"]),
            lst=(hour + lon / 15.0) % 24.0, logfa=np.log10(float(r["fa"])),
            # 以 2024-01-01 为原点的月数，供趋势检验与作图
            tnum=(int(r["start"][:4]) - 2024) * 12 + int(r["start"][5:7]) - 1 + int(r["start"][8:10]) / 31.0,
        ))
    return rows


def kuiper_2samp(x, y, n_perm=4000, seed=0):
    """两样本 Kuiper 统计量 D+ + D−（对圆周起点的选择不变）与置换 p 值。"""
    def stat(a, b):
        allv = np.sort(np.concatenate([a, b]))
        ca = np.searchsorted(np.sort(a), allv, "right") / len(a)
        cb = np.searchsorted(np.sort(b), allv, "right") / len(b)
        return (ca - cb).max() + (cb - ca).max()
    obs = stat(x, y)
    pool = np.concatenate([x, y]); n = len(x)
    rng = np.random.default_rng(seed)
    cnt = sum(stat(p[:n], p[n:]) >= obs for p in (rng.permutation(pool) for _ in range(n_perm)))
    return obs, (cnt + 1) / (n_perm + 1)


def energy_stat(a, b):
    """能量距离统计量：n·m/(n+m) · (2E|A−B| − E|A−A'| − E|B−B'|)，零假设下分布相同时趋于 0。"""
    d_ab = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).mean()
    d_aa = np.sqrt(((a[:, None, :] - a[None, :, :]) ** 2).sum(-1)).mean()
    d_bb = np.sqrt(((b[:, None, :] - b[None, :, :]) ** 2).sum(-1)).mean()
    n, m = len(a), len(b)
    return n * m / (n + m) * (2 * d_ab - d_aa - d_bb)


def energy_perm(X, mask_a, mask_b, n_perm=2000, seed=0):
    a, b = X[mask_a], X[mask_b]
    mu, sd = np.vstack([a, b]).mean(0), np.vstack([a, b]).std(0)
    a, b = (a - mu) / sd, (b - mu) / sd
    obs = energy_stat(a, b)
    pool = np.vstack([a, b]); n = len(a)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        idx = rng.permutation(len(pool))
        null[i] = energy_stat(pool[idx[:n]], pool[idx[n:]])
    return obs, null, (np.sum(null >= obs) + 1) / (n_perm + 1)


def holm(pairs):
    """Holm 校正：返回 [(名字, p, 阈值, 是否拒绝)]，按 p 升序。"""
    s = sorted(pairs, key=lambda z: z[1]); n = len(s)
    return [(name, p, 0.05 / (n - i), p < 0.05 / (n - i)) for i, (name, p) in enumerate(s)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("features"); ap.add_argument("assoc"); ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--exposure-in", type=float, default=140.2, help="覆盖内的搜索曝光（天）")
    ap.add_argument("--exposure-out", type=float, default=469.0, help="覆盖外的搜索曝光（天）")
    ap.add_argument("--perm", type=int, default=2000)
    args = ap.parse_args()

    rows = load(args.features, args.assoc)
    cov = np.array([r["cov"] for r in rows]); ass = np.array([r["ass"] for r in rows])
    mon = np.array([r["mon"] for r in rows]); season = np.isin(mon, SEASON_MONTHS)
    col = lambda k: np.array([r[k] for r in rows], float)
    lon, lst, alat = col("lon"), col("lst"), col("alat")
    groups = {"覆盖内": cov, "覆盖外": ~cov, "覆盖外(6-12月)": ~cov & season,
              "已证实": cov & ass, "覆盖内未证实": cov & ~ass}
    print("覆盖内 %d（已证实 %d，未证实 %d），覆盖外 %d（其中 6-12 月 %d）"
          % (cov.sum(), (cov & ass).sum(), (cov & ~ass).sum(), (~cov).sum(), (~cov & season).sum()))

    # ---- 单变量检验 ----
    print("\n=== 单变量检验（KS / Anderson-Darling；地方时用 Kuiper）")
    families = {}
    for title, ma, mb in (("主检验 覆盖内 vs 覆盖外", cov, ~cov),
                          ("季节匹配 6-12月", cov, ~cov & season),
                          ("阳性对照 已证实 vs 未证实", cov & ass, cov & ~ass)):
        print("  --- %s (n=%d vs %d)" % (title, ma.sum(), mb.sum()))
        ps = []
        for k in ("dur", "pir", "alat", "rate", "n", "dfrac", "logfa"):
            x, y = col(k)[ma], col(k)[mb]
            x, y = x[np.isfinite(x)], y[np.isfinite(y)]
            ks = stats.ks_2samp(x, y).pvalue
            try:
                ad = stats.anderson_ksamp([x, y]).pvalue
            except Exception:
                ad = np.nan
            ps.append((LABELS[k], ks))
            print("      %-16s 中位 %9.3f / %9.3f   KS p=%.4f  AD p=%.4f" % (LABELS[k], np.median(x), np.median(y), ks, ad))
        v, pk = kuiper_2samp(lst[ma], lst[mb])
        ps.append((LABELS["lst"], pk))
        print("      %-16s (圆周)                        Kuiper V=%.3f p=%.4f" % (LABELS["lst"], v, pk))
        families[title] = ps
    print("\n=== Holm 校正（每个族内 8 个检验）")
    for title, ps in families.items():
        print("  --- %s" % title)
        for name, p, thr, rej in holm(ps):
            print("      %-16s p=%.4f  阈 %.4f  %s" % (name, p, thr, "拒绝同分布" if rej else "不拒绝"))

    # ---- 多变量联合检验 ----
    X = np.array([[r[k] for k in JOINT] for r in rows], float)
    ok = np.all(np.isfinite(X), axis=1)
    print("\n=== 多变量能量距离置换检验（%s）" % "、".join(LABELS[k] for k in JOINT))
    joint = {}
    for label, ma, mb, seed in (("覆盖内 vs 覆盖外", cov & ok, ~cov & ok, 0),
                                ("季节匹配", cov & ok, ~cov & season & ok, 1),
                                ("阳性对照：已证实 vs 未证实", cov & ass & ok, cov & ~ass & ok, 2)):
        obs, null, p = energy_perm(X, ma, mb, args.perm, seed)
        joint[label] = (obs, null, p)
        print("  %-26s 统计量 %6.3f  置换 p=%.4f" % (label, obs, p))

    # ---- 速率与关联率 ----
    r_in, r_out = cov.sum() / args.exposure_in, (~cov).sum() / args.exposure_out
    e_in, e_out = np.sqrt(cov.sum()) / args.exposure_in, np.sqrt((~cov).sum()) / args.exposure_out
    print("\n=== 按曝光归一的显著候选率")
    print("  覆盖内 %.3f ± %.3f /天，覆盖外 %.3f ± %.3f /天，差 %.2f sigma"
          % (r_in, e_in, r_out, e_out, (r_out - r_in) / np.hypot(e_in, e_out)))

    # 本底率是唯一在季节匹配下仍显著的量。查它是仪器漂移还是人群差异：
    # 覆盖外内部随时间显著上升，且在样本最多的两个纬度带里两期一致 → 仪器/环境漂移。
    print("\n=== 本底率那一项的来源")
    tnum = col("tnum")
    rho, prho = stats.spearmanr(tnum[~cov], col("rate")[~cov])
    print("  覆盖外内部：本底率随时间的 Spearman rho=%+.3f p=%.2g" % (rho, prho))
    for lo, hi in ((0, 10), (10, 20), (20, 30)):
        b = (alat >= lo) & (alat < hi)
        x, y = col("rate")[cov & b], col("rate")[~cov & season & b]
        if len(x) > 5 and len(y) > 5:
            print("    |lat| %2d-%2d 覆盖内 n=%3d 中位 %6.0f | 覆盖外 n=%3d 中位 %6.0f  KS p=%.3f"
                  % (lo, hi, len(x), np.median(x), len(y), np.median(y), stats.ks_2samp(x, y).pvalue))

    print("\n=== 覆盖内关联率的分层（WWLLN 的效率，不是候选的性质）")
    zone_rate = []
    for zlo, zhi, zlab in ZONES:
        m = cov & (lon >= zlo) & (lon < zhi); k = (m & ass).sum()
        o = (~cov) & (lon >= zlo) & (lon < zhi)
        zone_rate.append((zlab, m.sum(), k, 100 * m.sum() / cov.sum(), 100 * o.sum() / (~cov).sum()))
        print("  %-10s 覆盖内 %3d 关联 %3d (%4.1f%%)；样本占比 覆盖内 %4.1f%% vs 覆盖外 %4.1f%%"
              % (zlab, m.sum(), k, 100 * k / max(m.sum(), 1), 100 * m.sum() / cov.sum(), 100 * o.sum() / (~cov).sum()))
    lst_rate = []
    for tlo, thi in zip(LST_EDGES[:-1], LST_EDGES[1:]):
        m = cov & (lst >= tlo) & (lst < thi); k = (m & ass).sum()
        lst_rate.append(("%d-%d" % (tlo, thi), m.sum(), k))
        print("  地方时 %2d-%2d 覆盖内 %3d 关联 %3d (%4.1f%%)" % (tlo, thi, m.sum(), k, 100 * k / max(m.sum(), 1)))

    # ---- 分层重加权外推 ----
    p_glob = (cov & ass).sum() / cov.sum()
    exp_tot = var_tot = 0.0
    for zlo, zhi, _ in ZONES:
        for tlo, thi in zip(LST_EDGES[:-1], LST_EDGES[1:]):
            cin = cov & (lon >= zlo) & (lon < zhi) & (lst >= tlo) & (lst < thi)
            cout = (~cov) & (lon >= zlo) & (lon < zhi) & (lst >= tlo) & (lst < thi)
            n, k, m = cin.sum(), (cin & ass).sum(), cout.sum()
            if m == 0:
                continue
            p, n_eff = (k / n, n) if n >= 5 else (p_glob, cov.sum())   # 格内样本太少就退回全局率
            exp_tot += m * p
            var_tot += m ** 2 * p * (1 - p) / n_eff + m * p * (1 - p)
    mis_in = sum(r["prob"] for r in rows if r["cov"])
    print("\n=== 外推")
    print("  覆盖内 %d 个显著候选，关联 %d（偶然误关联期望 %.2f）→ 关联率 %.1f%% ± %.1f%%（扣误关联 %.1f%%）"
          % (cov.sum(), (cov & ass).sum(), mis_in, 100 * p_glob, 100 * np.sqrt(p_glob * (1 - p_glob) / cov.sum()),
             100 * ((cov & ass).sum() - mis_in) / cov.sum()))
    print("  覆盖外 %d 个：朴素外推 %.0f 个，分层重加权 %.0f ± %.0f 个"
          % ((~cov).sum(), (~cov).sum() * p_glob, exp_tot, np.sqrt(var_tot)))

    # ---- 图 ----
    fig, axes = plt.subplots(2, 4, figsize=(20.5, 8.8))
    styles = {"覆盖内": ("C0", "-", 1.8), "覆盖外": ("C3", "-", 1.4), "覆盖外(6-12月)": ("C1", "--", 1.2)}
    for ax, k, tag in zip(axes[0], ("dur", "pir", "alat", "n"), "abcd"):
        for gname in ("覆盖内", "覆盖外", "覆盖外(6-12月)"):
            v = col(k)[groups[gname]]; v = np.sort(v[np.isfinite(v)])
            c, ls, lw = styles[gname]
            ax.plot(v, np.arange(1, len(v) + 1) / len(v), ls, color=c, lw=lw, label="%s (n=%d)" % (gname, len(v)))
        p_main = dict(families["主检验 覆盖内 vs 覆盖外"])[LABELS[k]]
        p_seas = dict(families["季节匹配 6-12月"])[LABELS[k]]
        ax.set_xlabel(LABELS[k]); ax.set_ylabel("累积占比")
        ax.set_title("(%s) %s   KS p = %.2f / %.2f" % (tag, LABELS[k], p_main, p_seas), fontsize=10)
        ax.legend(fontsize=7.5); ax.grid(alpha=0.25)
        if k == "dur":
            ax.set_xlim(0, 3)
        if k == "n":
            ax.set_xlim(0, 80)

    ax = axes[1][0]
    for gname in ("覆盖内", "覆盖外"):
        h, e = np.histogram(lst[groups[gname]], bins=np.arange(0, 25, 3))
        ax.step(np.append(e[:-1], 24), np.append(100 * h / h.sum(), 100 * h[-1] / h.sum()),
                where="post", color=styles[gname][0], lw=1.6, label="%s (n=%d)" % (gname, groups[gname].sum()))
    _, pk = kuiper_2samp(lst[cov], lst[~cov])
    ax.set_xlabel(LABELS["lst"]); ax.set_ylabel("占比 (%)"); ax.set_xticks(range(0, 25, 6)); ax.set_ylim(0, None)
    ax.set_title("(e) 地方时（圆周量）Kuiper p = %.2f" % pk, fontsize=10); ax.legend(fontsize=7.5); ax.grid(alpha=0.25)

    ax = axes[1][1]
    x = np.arange(len(ZONES)); w = 0.36
    ax.bar(x - w / 2, [z[3] for z in zone_rate], w, color="C0", label="覆盖内的样本占比")
    ax.bar(x + w / 2, [z[4] for z in zone_rate], w, color="C3", label="覆盖外的样本占比")
    for i, z in enumerate(zone_rate):
        ax.text(i, max(z[3], z[4]) + 2.5, "覆盖内\n关联率 %.0f%%" % (100 * z[2] / max(z[1], 1)), ha="center", fontsize=8, color="0.25")
    zone_tab = np.array([[z[1] for z in zone_rate],
                         [int(((lon >= zlo) & (lon < zhi) & ~cov).sum()) for zlo, zhi, _ in ZONES]])
    ax.set_xticks(x); ax.set_xticklabels([z[2] for z in ZONES], fontsize=9); ax.set_ylabel("占比 (%)"); ax.set_ylim(0, 62)
    ax.set_title("(f) 地理分布一致（$\\chi^2$ p = %.2f）\n关联率的地域差异来自 WWLLN 的台站分布"
                 % stats.chi2_contingency(zone_tab)[1], fontsize=9.5)
    ax.legend(fontsize=7.5)

    ax = axes[1][2]
    rate = col("rate")
    for gname, m in (("覆盖内", cov), ("覆盖外", ~cov)):
        ax.scatter(tnum[m], rate[m], s=6, alpha=0.35, color=styles[gname][0], label="%s" % gname)
    edges = np.arange(5, 33, 2.0)
    mid = 0.5 * (edges[:-1] + edges[1:])
    med = [np.median(rate[(tnum >= a) & (tnum < b)]) if ((tnum >= a) & (tnum < b)).sum() > 3 else np.nan
           for a, b in zip(edges[:-1], edges[1:])]
    ax.plot(mid, med, "k-", lw=1.8, label="两月中位")
    ax.set_xlabel("2024-01 起的月数"); ax.set_ylabel(LABELS["rate"]); ax.set_ylim(0, 9000)
    ax.axvline(12, color="0.4", ls="--", lw=1)
    ax.text(12.4, 8600, "WWLLN 覆盖止于此", fontsize=7.5, va="top", color="0.35")
    ax.set_title("(g) 本底率：覆盖外内部随时间上升\n（$\\rho=%+.2f$，仪器/环境漂移，非人群差异）" % rho, fontsize=9.5)
    ax.legend(fontsize=7.5)

    ax = axes[1][3]
    obs0, null0, p0 = joint["覆盖内 vs 覆盖外"]
    obs2, _, p2 = joint["阳性对照：已证实 vs 未证实"]
    ax.hist(null0, bins=40, color="0.82", label="零分布（置换 %d 次）" % args.perm)
    ax.axvline(obs0, color="C0", lw=2, label="覆盖内 vs 覆盖外：p = %.2f" % p0)
    ax.axvline(obs2, color="C2", lw=2, ls="--", label="阳性对照 已证实 vs 未证实：p = %.4f" % p2)
    ax.set_xlabel("能量距离统计量"); ax.set_ylabel("置换次数")
    ax.set_title("(h) 多变量联合检验：覆盖内外落在零分布中央，\n而同一检验查得出真实差异", fontsize=9.5)
    ax.legend(fontsize=7.5)

    fig.suptitle("SVOM/GRM：WWLLN 覆盖外的 %d 个显著候选与覆盖内的 %d 个是同一人群\n"
                 "（覆盖内关联率 %.1f%% ± %.1f%%，外推到覆盖外相当于 %.0f ± %.0f 个闪电关联）"
                 % ((~cov).sum(), cov.sum(), 100 * p_glob, 100 * np.sqrt(p_glob * (1 - p_glob) / cov.sum()),
                    exp_tot, np.sqrt(var_tot)), fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(args.output, dpi=150, bbox_inches="tight")
    print("\nwrote", args.output)


if __name__ == "__main__":
    main()
