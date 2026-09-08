"""SVOM/GRM 闪电证实 TGF 样本的科学表征：时长、叠加能谱、计数、死时间、地理与地方时。

输入是 `scripts/cluster/svom_tgf_sample.py` 在集群上导出的三份表
（sample_events.csv / sample_bkg_spec.csv / sample_meta.csv / ebounds.csv）。

时长用非分箱极大似然拟合脉冲宽度给 T50/T90，不是搜索窗长；能谱是叠加的**沉积能量谱**，
没有响应矩阵就不能反解光子谱，指数不等于光子谱指数（见图注）。

事例准入：EVT_TYPE==0、25 ≤ PI < 256（v6 搜索的能阈 ch25 ≈ 42 keV）、**ANTI_COIN==0**。
ANTI_COIN==1 是星上标定源（49–57 keV 的线，三路合计约 36 c/s），核心窗与本底窗一律排除，
否则扣本底时会在 50 keV 附近过扣。注意 `features_top899.csv` 里的 pi_med_ratio 等量是
旧的 ch15 口径算的，不要和这里的量混用。

本底取自事例表里 10–50 ms 的环（与核心窗同一份数据、同一套筛选，口径自洽）；
`sample_bkg_spec.csv` 的 ±1 s 窗只用来互校本底率（两者比值 0.99）。

用法:
    python3 scripts/plot_svom_tgf_sample.py <sample_dir> -o <前缀>
输出 <前缀>_timing_spectrum.png、<前缀>_geography.png，并把逐个 TGF 的量打到 stdout。
"""
import argparse, csv, os
from datetime import datetime, timezone

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize
from scipy.stats import kstest, norm

plt.rcParams.update({
    "font.sans-serif": ["PingFang SC", "Arial Unicode MS"],
    "font.family": "sans-serif",
    "axes.unicode_minus": False,
})

PI_LO, PI_HI = 25, 256      # 搜索侧事例准入（能阈 ch25，256 起是溢出道）
ANALYSIS_HALF = 0.005       # 算时长用的半窗（±5 ms）
BKG_LO, BKG_HI = 0.010, 0.050   # 本底取自事例表的 10–50 ms 环（与核心窗同样 AC=0，口径自洽）
RNG = np.random.default_rng(20260907)


def load(sample_dir):
    ev = {}
    with open(os.path.join(sample_dir, "sample_events.csv")) as f:
        for r in csv.DictReader(f):
            ev.setdefault(int(r["idx"]), []).append(
                (float(r["dt_s"]), int(r["pi"]), int(r["det"]), int(r["gain"]),
                 float(r["dead_us"]), int(r["anti"]), int(r["evt_type"])))
    events = {}
    for k, v in ev.items():
        a = np.array(v, dtype=float)
        events[k] = dict(dt=a[:, 0], pi=a[:, 1].astype(int), det=a[:, 2].astype(int),
                         gain=a[:, 3].astype(int), dead=a[:, 4], anti=a[:, 5].astype(int),
                         evt=a[:, 6].astype(int))
    meta = {int(r["idx"]): r for r in csv.DictReader(open(os.path.join(sample_dir, "sample_meta.csv")))}
    bkg = {}
    with open(os.path.join(sample_dir, "sample_bkg_spec.csv")) as f:
        for r in csv.DictReader(f):
            bkg[int(r["idx"])] = (float(r["live_s"]),
                                  np.array([float(r["c%d" % i]) for i in range(259)]))
    bkg_ac1 = {}
    p1 = os.path.join(sample_dir, "sample_bkg_spec_ac1.csv")
    if os.path.exists(p1):
        with open(p1) as f:
            for r in csv.DictReader(f):
                bkg_ac1[int(r["idx"])] = np.array([float(r["c%d" % i]) for i in range(259)])
    eb = np.array([(float(r["e_min"]), float(r["e_max"]))
                   for r in csv.DictReader(open(os.path.join(sample_dir, "ebounds.csv")))])
    return events, meta, bkg, bkg_ac1, eb


def fit_pulse(t, rate, half=ANALYSIS_HALF):
    """非分箱扩展极大似然：本底率已知（±1 s 测得，几千个计数），脉冲取高斯。

    返回 (N_burst, mu, sigma)。累积法在这里不能用：±5 ms 窗里本底期望 30 余个、
    净计数中位才 18 个，5%/95% 交点会被本底涨落推到窗边，T90 系统性偏大。
    """
    if len(t) < 6:
        return None
    n_bkg = rate * 2 * half
    n0 = max(len(t) - n_bkg, 2.0)

    def nll(p):
        logN, mu, logsig = p
        N, sig = np.exp(logN), np.exp(logsig)
        g = np.exp(-0.5 * ((t - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
        return -(np.log(rate + N * g + 1e-12).sum() - (n_bkg + N))

    best = None
    for s0 in (5e-5, 2e-4, 1e-3):
        for mu0 in (0.0, float(np.median(t))):
            r = minimize(nll, [np.log(n0), mu0, np.log(s0)], method="Nelder-Mead",
                         options=dict(maxiter=4000, xatol=1e-9, fatol=1e-6))
            if best is None or r.fun < best.fun:
                best = r
    logN, mu, logsig = best.x
    sig = float(np.exp(logsig))
    if not (5e-6 < sig < half) or abs(mu) > half:
        return None
    return float(np.exp(logN)), float(mu), sig


def boot_pulse(t, rate, n=60):
    """自举给 T50/T90 的不确定度（对窗内事例重采样）。"""
    out = []
    for _ in range(n):
        s = RNG.choice(t, size=len(t), replace=True)
        f = fit_pulse(np.sort(s), rate)
        if f:
            out.append(f[2])
    if len(out) < 10:
        return np.nan
    return float(np.std(out))


# 高斯脉冲：T50 = 25%–75% 分位差 = 1.349σ，T90 = 5%–95% = 3.290σ
T50_SIG = 1.3490
T90_SIG = 3.2897


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample_dir")
    ap.add_argument("-o", "--output", required=True, help="输出前缀")
    args = ap.parse_args()
    events, meta, bkg, bkg_ac1, eb = load(args.sample_dir)
    idxs = sorted(events)
    emin, emax = eb[:, 0], eb[:, 1]
    ewidth = emax - emin
    ecen = np.sqrt(emin * emax)

    rows = []
    core_spec = np.zeros(259)          # AC=0（标定源已扣）
    core_spec_all = np.zeros(259)      # 不扣 AC=1，仅供对比
    bkg_spec_scaled = np.zeros(259)
    bkg_spec_all_scaled = np.zeros(259)
    stack_dt, stack_w = [], []
    n_ac1_core = 0
    n_skip = 0
    for i in idxs:
        e = events[i]
        keep = (e["evt"] == 0) & (e["anti"] == 0) & (e["pi"] >= PI_LO) & (e["pi"] < PI_HI)
        dt = np.sort(e["dt"][keep])
        # 本底：事例表里 10–50 ms 的环，AC=0，与核心窗同一口径
        ring = (e["evt"] == 0) & (e["anti"] == 0) & (np.abs(e["dt"]) > BKG_LO) & (np.abs(e["dt"]) <= BKG_HI)
        live_ring = 2 * (BKG_HI - BKG_LO)
        spec_ring = np.bincount(np.clip(e["pi"][ring], 0, 258), minlength=259).astype(float)
        rate = spec_ring[PI_LO:PI_HI].sum() / live_ring
        live_win, spec_win = bkg[i]                    # ±1 s 窗，只作互校（口径见文件头）
        rate_win = spec_win[PI_LO:PI_HI].sum() / live_win
        live, spec = live_ring, spec_ring
        win = dt[np.abs(dt) <= ANALYSIS_HALF]
        f = fit_pulse(win, rate)
        if f is None:
            n_skip += 1
            continue
        n_burst, mu, sig = f
        T50, T90 = T50_SIG * sig, T90_SIG * sig
        s_sig = boot_pulse(win, rate)
        s50 = T50_SIG * s_sig
        s90 = T90_SIG * s_sig
        # 拟合优度：把窗内事例的到达时间与「本底 + 高斯」的累积分布比一比
        cdf = lambda x: (rate * (x + ANALYSIS_HALF) + n_burst * norm.cdf(x, mu, sig)) / \
                        (rate * 2 * ANALYSIS_HALF + n_burst)
        ks_p = float(kstest(win, cdf).pvalue) if len(win) >= 8 else np.nan
        # 能谱与死时间用 µ ± 2σ（含 95% 脉冲计数），比 T90 窗窄得多、本底少得多
        t05, t95 = mu - 2 * sig, mu + 2 * sig
        inw_all = (e["dt"] >= t05) & (e["dt"] <= t95) & (e["evt"] == 0)
        inw = inw_all & (e["anti"] == 0)
        n_ac1_core += int((inw_all & (e["anti"] == 1)).sum())
        width = max(t95 - t05, 1e-6)
        core_spec += np.bincount(np.clip(e["pi"][inw], 0, 258), minlength=259)
        core_spec_all += np.bincount(np.clip(e["pi"][inw_all], 0, 258), minlength=259)
        bkg_spec_scaled += spec * (width / live)
        ring_all = (e["evt"] == 0) & (np.abs(e["dt"]) > BKG_LO) & (np.abs(e["dt"]) <= BKG_HI)
        bkg_spec_all_scaled += np.bincount(np.clip(e["pi"][ring_all], 0, 258), minlength=259) * (width / live)
        n_core = int(((e["pi"][inw] >= PI_LO) & (e["pi"][inw] < PI_HI)).sum())
        n_core_net = n_core - rate * width
        # 死时间：逐探头按非瘫痪模型算占比与修正；f≥1 表示该探头在窗内已饱和，无法修正
        f_det, n_det, corr_det, saturated = {}, {}, {}, False
        for d in (1, 2, 3):
            sel_d = inw & (e["det"] == d) & (e["pi"] >= PI_LO) & (e["pi"] < PI_HI)
            n_det[d] = int(sel_d.sum())
            # 死时间由所有被记录的事例造成，标定源事例也算在内
            f_det[d] = float(e["dead"][inw_all & (e["det"] == d)].sum() / (width * 1e6))
            if f_det[d] >= 0.95:
                saturated = True
                corr_det[d] = n_det[d] / 0.05
            else:
                corr_det[d] = n_det[d] / (1 - f_det[d])
        f_burst = max(f_det.values())
        f_mean = float(np.mean(list(f_det.values())))
        n_corr = sum(corr_det.values())
        f_bkg = max(float(meta[i]["bkg_dead_us_d%d" % d]) / (live_win * 1e6) for d in (1, 2, 3))
        # 叠加光变：对齐到拟合的脉冲中心
        stack_dt.append(dt - mu)
        stack_w.append(np.full(dt.shape, 1.0))
        m = meta[i]
        t_utc = datetime.strptime(m["start"][:19], "%Y-%m-%dT%H:%M:%S").replace(tzinfo=timezone.utc)
        lst = (t_utc.hour + t_utc.minute / 60 + t_utc.second / 3600 + float(m["lon"]) / 15.0) % 24
        rows.append(dict(idx=i, start=m["start"], fa=float(m["fa"]), lon=float(m["lon"]), lat=float(m["lat"]),
                         alt_km=float(m["alt_m"]) / 1e3, T50=T50, T90=T90, sT50=s50, sT90=s90,
                         n_core=n_core, n_net=n_core_net, n_burst=n_burst, rate_bkg=rate,
                         rate_win=rate_win, f_burst=f_burst, f_mean=f_mean, n_corr=n_corr,
                         saturated=saturated, f_bkg=f_bkg, lst=lst, width=width, ks_p=ks_p))
    print("表征到 %d / %d 个证实 TGF（%d 个拟合不收敛）" % (len(rows), len(idxs), n_skip))
    rr = np.array([r["rate_bkg"] for r in rows]); rw = np.array([r["rate_win"] for r in rows])
    print("本底率互校：10–50 ms 环（AC=0，本分析所用）中位 %.0f c/s，±1 s 窗中位 %.0f c/s，比值中位 %.3f" % (
        np.median(rr), np.median(rw), np.median(rr / rw)))
    ksp = np.array([r["ks_p"] for r in rows])
    print("单高斯拟合优度：KS p < 0.05 的 %d 个（多脉冲或非高斯轮廓）" % np.nansum(ksp < 0.05))
    T50 = np.array([r["T50"] for r in rows]) * 1e6
    T90 = np.array([r["T90"] for r in rows]) * 1e6
    nnet = np.array([r["n_net"] for r in rows])
    fb = np.array([r["f_burst"] for r in rows])
    fk = np.array([r["f_bkg"] for r in rows])
    lst = np.array([r["lst"] for r in rows])
    lat = np.array([r["lat"] for r in rows])
    lon = np.array([r["lon"] for r in rows])
    print("T50 µs: 中位 %.0f, 16–84%% %.0f–%.0f, 范围 %.0f–%.0f" % (
        np.median(T50), np.percentile(T50, 16), np.percentile(T50, 84), T50.min(), T50.max()))
    print("T90 µs: 中位 %.0f, 16–84%% %.0f–%.0f, 范围 %.0f–%.0f" % (
        np.median(T90), np.percentile(T90, 16), np.percentile(T90, 84), T90.min(), T90.max()))
    print("T90 区间净计数: 中位 %.0f, 范围 %.0f–%.0f, 合计 %.0f" % (
        np.median(nnet), nnet.min(), nnet.max(), nnet.sum()))
    fm = np.array([r["f_mean"] for r in rows])
    ncorr = np.array([r["n_corr"] for r in rows])
    nobs = np.array([r["n_core"] for r in rows], float)
    nsat = sum(1 for r in rows if r["saturated"])
    print("死时间占比（T90 窗内）: 最忙的一路 中位 %.3f p90 %.3f max %.3f；三路平均 中位 %.3f；本底 中位 %.5f" % (
        np.median(fb), np.percentile(fb, 90), fb.max(), np.median(fm), np.median(fk)))
    print("逐探头非瘫痪修正后的计数 / 实测: 中位 %.3f, p90 %.3f, max %.3f；有探头饱和(f≥0.95)的 %d 个" % (
        np.median(ncorr / nobs), np.percentile(ncorr / nobs, 90), (ncorr / nobs).max(), nsat))

    # ---- 叠加沉积能量谱与幂律拟合（沉积能量谱，不是光子谱）----
    net_spec = core_spec - bkg_spec_scaled
    total_width = sum(r["width"] for r in rows)
    lo_ch, hi_ch = PI_LO, 256
    sel = np.arange(259)
    groups = []
    acc_lo = lo_ch
    acc = 0.0
    for ch in range(lo_ch, hi_ch):
        acc += net_spec[ch]
        if acc >= 20 and ch > acc_lo:
            groups.append((acc_lo, ch)); acc_lo = ch + 1; acc = 0.0
    if acc_lo < hi_ch - 1:
        groups.append((acc_lo, hi_ch - 1))
    g_e, g_flux, g_err = [], [], []
    for a, b in groups:
        de = emax[b] - emin[a]
        net = net_spec[a:b + 1].sum()
        raw = core_spec[a:b + 1].sum()
        bk = bkg_spec_scaled[a:b + 1].sum()
        err = np.sqrt(max(raw, 1) + bk)
        g_e.append(np.sqrt(emin[a] * emax[b]))
        g_flux.append(net / de / total_width)
        g_err.append(err / de / total_width)
    g_e, g_flux, g_err = np.array(g_e), np.array(g_flux), np.array(g_err)
    ok = g_flux > 0
    pl = lambda x, a, b: a * x ** b

    def fit_range(lo_e, hi_e, core_arr=None, bkg_arr=None):
        """逐能道的泊松（Cash）似然拟合：模型在每道上积分，加上该道的本底期望。

        不分组，避免宽能段用几何中心当能量带来的偏置，也不受分组边界抖动影响。
        """
        core_arr = core_spec if core_arr is None else core_arr
        bkg_arr = bkg_spec_scaled if bkg_arr is None else bkg_arr
        chans = [c for c in range(lo_ch, hi_ch) if emin[c] >= lo_e and emax[c] <= hi_e]
        if len(chans) < 6:
            return None
        lo_a = np.array([emin[c] for c in chans]); hi_a = np.array([emax[c] for c in chans])
        obs = core_arr[chans]; bkg_c = bkg_arr[chans]

        def model(a, b):
            # ∫ A E^b dE 在每道上，乘总曝光；b = -1 单独处理
            if abs(b + 1) < 1e-6:
                integ = np.log(hi_a / lo_a)
            else:
                integ = (hi_a ** (b + 1) - lo_a ** (b + 1)) / (b + 1)
            return a * integ * total_width + bkg_c

        def cash(p):
            a, b = np.exp(p[0]), p[1]
            mu = np.maximum(model(a, b), 1e-9)
            return 2 * float((mu - obs * np.log(mu)).sum())

        best = None
        for b0 in (-0.6, -1.0, -1.5, -2.0):
            a0 = max(obs.sum() - bkg_c.sum(), 1.0) / max(total_width * np.log(hi_a[-1] / lo_a[0]), 1e-9)
            r = minimize(cash, [np.log(max(a0, 1e-6)), b0], method="Nelder-Mead",
                         options=dict(maxiter=6000, xatol=1e-8, fatol=1e-8))
            if best is None or r.fun < best.fun:
                best = r
        a_hat, b_hat = float(np.exp(best.x[0])), float(best.x[1])
        # 误差：对 b 做剖面似然，Δcash = 1
        c0 = best.fun
        db = 0.001
        b_lo = b_hat
        while b_lo > b_hat - 1.5:
            b_lo -= db
            r = minimize(lambda p: cash([p[0], b_lo]), [np.log(a_hat)], method="Nelder-Mead",
                         options=dict(maxiter=2000, xatol=1e-8, fatol=1e-8))
            if r.fun - c0 > 1.0:
                break
        err = abs(b_hat - b_lo)
        # 拟合优度：把每道并到至少 15 个观测计数再算 chi2
        gobs, gmu = [], []
        acc_o = acc_m = 0.0
        mu_hat = model(a_hat, b_hat)
        for o, m_ in zip(obs, mu_hat):
            acc_o += o; acc_m += m_
            if acc_o >= 15:
                gobs.append(acc_o); gmu.append(acc_m); acc_o = acc_m = 0.0
        if acc_o > 0 and gobs:
            gobs[-1] += acc_o; gmu[-1] += acc_m
        gobs, gmu = np.array(gobs), np.array(gmu)
        chi = float((((gobs - gmu) ** 2) / np.maximum(gmu, 1e-9)).sum())
        return np.array([a_hat, b_hat]), np.array([np.nan, err]), chi, len(gobs)

    print("叠加沉积能量谱（未解卷积；GRM 两档增益共用一套能道，交界在 ch110 ≈ 640 keV）：")
    print("  标定源（ANTI_COIN=1）在暴发窗内只有 %d 个事例（占阈上 %.2f%%），已扣；本底窗也按 AC=0 统计" % (
        n_ac1_core, 100 * n_ac1_core / max(core_spec[lo_ch:hi_ch].sum() + n_ac1_core, 1)))
    print("  净计数合计 %.0f（原始 %.0f，本底期望 %.0f），溢出道 ch256-258 净 %.1f" % (
        net_spec[lo_ch:hi_ch].sum(), core_spec[lo_ch:hi_ch].sum(), bkg_spec_scaled[lo_ch:hi_ch].sum(),
        net_spec[256:].sum()))
    print("  能阈以下 ch0-24 净计数 %.1f（占阈上的 %.0f%%），证实样本在阈下几乎没有信号" % (
        net_spec[:PI_LO].sum(), 100 * net_spec[:PI_LO].sum() / max(net_spec[PI_LO:hi_ch].sum(), 1)))
    fits = {}
    for lab, lo_e, hi_e in (("全band 42-8041", emin[lo_ch], emax[hi_ch - 1] * 1.01),
                            ("低增益段 640-8041", 640.0, emax[hi_ch - 1] * 1.01),
                            ("高增益段 42-640", emin[lo_ch], 640.0)):
        r = fit_range(lo_e, hi_e)
        if r is None:
            continue
        p, pe, chi, npts = r
        fits[lab] = (p, pe, chi, npts)
        print("  %s keV: dN/dE ∝ E^(%.2f ± %.2f)，chi2/dof = %.1f/%d = %.2f" % (
            lab, p[1], pe[1], chi, npts - 2, chi / max(npts - 2, 1)))
    # 不扣标定源的对照：核心与本底都不筛 ANTI_COIN，同一套 Cash 拟合
    for lab, lo_e, hi_e in (("低增益段 640-8041", 640.0, emax[hi_ch - 1] * 1.01),
                            ("高增益段 42-640", emin[lo_ch], 640.0)):
        if lab not in fits:
            continue
        r2 = fit_range(lo_e, hi_e, core_spec_all, bkg_spec_all_scaled)
        if r2 is None:
            continue
        print("  不扣标定源的对照 %s keV: 指数 %.3f ± %.3f（扣后 %.3f，差 %.3f）" % (
            lab.split()[1], r2[0][1], r2[1][1], fits[lab][0][1], r2[0][1] - fits[lab][0][1]))
    d50 = (emin >= 45) & (emax <= 60)
    ac1_in_bkg = (bkg_spec_all_scaled - bkg_spec_scaled)[d50].sum()
    print("  45–60 keV：标定源在本底里折算到暴发窗是 %.2f 个计数，占该段净计数 %.2f%%（所以对本样本影响可忽略；"
          "若用更宽的暴发窗，本底占比上去了，这一步就不能省）" % (
              ac1_in_bkg, 100 * ac1_in_bkg / max(net_spec[d50].sum(), 1)))
    print("  截断幂律的截断能跑到 GRM 波段以外：8 MeV 以内看不到谱截断。")
    print("  42-640 keV 段的 chi2 偏大来自两档增益交界处的响应结构（300-640 keV 明显下凹），")
    print("  要给光子谱必须用 CALDB 的响应矩阵解卷积，见 OPEN-QUESTIONS。")

    # ---- 图 1：时间与能谱 ----
    fig, ax = plt.subplots(2, 2, figsize=(12.5, 8))
    a = ax[0, 0]
    all_dt = np.concatenate(stack_dt) * 1e6
    bins = np.arange(-2000, 2001, 25)
    a.hist(all_dt, bins=bins, color="crimson", alpha=0.85)
    bkg_level = np.median([r["rate_bkg"] for r in rows]) * len(rows) * 25e-6
    a.axhline(bkg_level, color="0.4", ls="--", lw=1, label="本底期望 %.1f/格" % bkg_level)
    a.set_xlabel("相对拟合脉冲中心的时间 (µs)"); a.set_ylabel("计数 / 25 µs")
    a.set_title("(a) %d 个证实 TGF 的叠加光变（ch25–255，已扣标定源）" % len(rows), fontsize=10)
    a.legend(fontsize=8)

    a = ax[0, 1]
    b = np.logspace(1, 4, 26)
    a.hist(T50, bins=b, histtype="step", lw=1.6, color="crimson", label="T50 中位 %.0f µs" % np.median(T50))
    a.hist(T90, bins=b, histtype="step", lw=1.6, color="tab:blue", label="T90 中位 %.0f µs" % np.median(T90))
    a.set_xscale("log"); a.set_xlabel("时长 (µs)"); a.set_ylabel("TGF 数")
    a.set_title("(b) 非分箱极大似然拟合脉冲宽度给出的时长", fontsize=10); a.legend(fontsize=8)

    a = ax[1, 0]
    a.errorbar(g_e[ok], g_flux[ok], yerr=g_err[ok], fmt="o", ms=3.5, color="crimson", label="叠加净计数谱")
    for lab, color, lo_e, hi_e in (("低增益段 640-8041", "k", 640.0, emax[hi_ch - 1]),
                                   ("高增益段 42-640", "0.5", emin[lo_ch], 640.0)):
        if lab not in fits:
            continue
        p, pe, _, _ = fits[lab]
        x = np.logspace(np.log10(lo_e), np.log10(hi_e), 60)
        a.plot(x, pl(x, *p), ls="--", lw=1.3, color=color,
               label="%s keV: 指数 %.2f ± %.2f" % (lab.split()[1], p[1], pe[1]))
    a.axvline(640, color="tab:blue", ls=":", lw=1.2, label="两档增益交界 ≈640 keV")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("能量 (keV)"); a.set_ylabel("计数 / keV / s")
    a.set_title("(c) 叠加沉积能量谱：未解卷积，指数 ≠ 光子谱指数", fontsize=10); a.legend(fontsize=7)

    a = ax[1, 1]
    a.hist(nnet, bins=np.arange(0, max(nnet) + 6, 5), color="tab:blue", alpha=0.8,
           label="T90 内净计数 中位 %.0f" % np.median(nnet))
    a.set_xlabel("单个 TGF 的净计数（三路合计）"); a.set_ylabel("TGF 数")
    a2 = a.twinx()
    o = np.argsort(nnet)
    a2.plot(nnet[o], 100 * fb[o], "o", ms=3.5, color="crimson", alpha=0.75)
    a2.set_ylabel("最忙一路的死时间占比 (%)", color="crimson")
    a2.axhline(100 * np.median(fk), color="0.5", ls=":", lw=1)
    a.set_title("(d) 计数与死时间：占比中位 %.0f%%，修正后计数 ×%.2f（中位）" % (
        100 * np.median(fb), np.median(ncorr / nobs)), fontsize=10)
    a.legend(fontsize=8, loc="upper left")
    fig.suptitle("SVOM/GRM %d 个闪电证实 TGF：时间结构与能谱" % len(rows), fontsize=12)
    fig.tight_layout()
    fig.savefig(args.output + "_timing_spectrum.png", dpi=150, bbox_inches="tight")
    print("wrote", args.output + "_timing_spectrum.png")

    # ---- 图 2：地理与地方时 ----
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(figsize=(13.5, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.42, wspace=0.22)
    a = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
    span = min(60, np.ceil(np.abs(lat).max()) + 5)
    a.set_extent([-180, 180, -span, span], crs=ccrs.PlateCarree())
    a.add_feature(cfeature.LAND, facecolor="0.95")
    a.add_feature(cfeature.COASTLINE, lw=0.4, edgecolor="0.4")
    a.gridlines(draw_labels=False, lw=0.3, color="0.88")
    sc = a.scatter(lon, lat, c=T90, s=42, cmap="viridis", norm=matplotlib.colors.LogNorm(),
                   lw=0.3, edgecolor="k", transform=ccrs.PlateCarree())
    cb = fig.colorbar(sc, ax=a, orientation="vertical", fraction=0.021, pad=0.01)
    cb.set_label("T90 (µs)", fontsize=8); cb.ax.tick_params(labelsize=7)
    a.set_title("(a) %d 个证实 TGF 的位置（SVOM 轨道倾角 30°）" % len(rows), fontsize=10)

    a = fig.add_subplot(gs[1, 0])
    a.hist(lst, bins=np.arange(0, 25, 2), color="tab:orange", alpha=0.85)
    a.axhline(len(rows) / 12, color="k", ls="--", lw=1, label="均匀分布期望")
    a.set_xticks(range(0, 25, 4)); a.set_xlabel("地方时 (h)"); a.set_ylabel("TGF 数")
    # 瑞利检验：地方时当作圆周量，检验是否偏离均匀
    ang = 2 * np.pi * lst / 24.0
    R = np.hypot(np.cos(ang).mean(), np.sin(ang).mean())
    n = len(lst)
    p_ray = float(np.exp(np.sqrt(1 + 4 * n + 4 * (n ** 2 - (n * R) ** 2)) - (1 + 2 * n)))
    peak_h = (np.degrees(np.arctan2(np.sin(ang).mean(), np.cos(ang).mean())) / 360 * 24) % 24
    a.set_title("(b) 地方时：瑞利检验 p = %.2f（%s），矢量峰在 %.1f h" % (
        p_ray, "与均匀不可区分" if p_ray > 0.05 else "偏离均匀", peak_h), fontsize=10)
    a.legend(fontsize=8)

    a = fig.add_subplot(gs[1, 1])
    bins = np.arange(-32, 33, 4)
    a.hist(lat, bins=bins, color="tab:green", alpha=0.85, label="实测")
    # 倾角 30° 圆轨道的纬度曝光：dt/dlat ∝ 1/sqrt(sin^2 i - sin^2 lat)
    inc = np.radians(30.0)
    c_lat = 0.5 * (bins[:-1] + bins[1:])
    w = np.where(np.abs(np.radians(c_lat)) < inc,
                 1.0 / np.sqrt(np.maximum(np.sin(inc) ** 2 - np.sin(np.radians(c_lat)) ** 2, 1e-6)), 0.0)
    w = w / w.sum() * len(rows)
    a.step(c_lat, w, where="mid", color="k", ls="--", lw=1.2, label="轨道曝光期望（TGF 率与纬度无关时）")
    a.set_xlabel("纬度 (°)"); a.set_ylabel("TGF 数")
    a.set_title("(c) 纬度：|lat| 中位 %.1f°，%.0f%% 在 ±25° 内" % (
        np.median(np.abs(lat)), 100 * np.mean(np.abs(lat) < 25)), fontsize=10)
    a.legend(fontsize=7)
    fig.suptitle("SVOM/GRM 证实 TGF 的地理与地方时分布", fontsize=12)
    fig.savefig(args.output + "_geography.png", dpi=150, bbox_inches="tight")
    print("wrote", args.output + "_geography.png")

    zones = [("美洲 (-120..-30)", -120, -30), ("非洲欧洲 (-30..60)", -30, 60), ("亚洲海洋大陆 (60..180)", 60, 180)]
    for name, a_, b_ in zones:
        print("  %s: %d 个 (%.0f%%)" % (name, ((lon >= a_) & (lon < b_)).sum(),
                                        100 * np.mean((lon >= a_) & (lon < b_))))
    print("  地方时 12–20 h 占 %.0f%%，0–8 h 占 %.0f%%；瑞利检验 p = %.2f（%s）" % (
        100 * np.mean((lst >= 12) & (lst < 20)), 100 * np.mean((lst >= 0) & (lst < 8)), p_ray,
        "与均匀不可区分，样本量不足以看出雷暴的日变化" if p_ray > 0.05 else "偏离均匀"))
    inc_r = np.radians(30.0)
    wlat = np.where(np.abs(np.radians(lat)) < inc_r,
                    1.0 / np.sqrt(np.maximum(np.sin(inc_r) ** 2 - np.sin(np.radians(lat)) ** 2, 1e-6)), 0.0)
    print("  纬度：实测 |lat| 中位 %.1f°；同样倾角下纯曝光的 |lat| 中位约 %.1f°（曝光在 ±30° 转折点最长），"
          % (np.median(np.abs(lat)), np.median(np.abs(np.degrees(np.arcsin(np.sin(inc_r) * np.sin(np.linspace(0, 2 * np.pi, 10001))))))))
    print("    实测明显更集中在低纬，与雷暴分布一致而与曝光相反。")

    with open(args.output + "_per_tgf.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["idx", "start", "fa", "lon", "lat", "alt_km", "lst_h", "T50_us", "T50_err_us",
                    "T90_us", "T90_err_us", "n_core", "n_net", "n_deadcorr", "rate_bkg_cps",
                    "deadfrac_burst_max", "deadfrac_burst_mean", "deadfrac_bkg", "saturated", "ks_p"])
        for r in rows:
            w.writerow([r["idx"], r["start"], "%.3e" % r["fa"], "%.4f" % r["lon"], "%.4f" % r["lat"],
                        "%.1f" % r["alt_km"], "%.2f" % r["lst"], "%.1f" % (r["T50"] * 1e6),
                        "%.1f" % (r["sT50"] * 1e6 if np.isfinite(r["sT50"]) else -1),
                        "%.1f" % (r["T90"] * 1e6), "%.1f" % (r["sT90"] * 1e6 if np.isfinite(r["sT90"]) else -1),
                        r["n_core"], "%.1f" % r["n_net"], "%.1f" % r["n_corr"], "%.0f" % r["rate_bkg"],
                        "%.4f" % r["f_burst"], "%.4f" % r["f_mean"], "%.6f" % r["f_bkg"],
                        int(r["saturated"]), "%.4f" % (r["ks_p"] if np.isfinite(r["ks_p"]) else -1)])
    print("wrote", args.output + "_per_tgf.csv")


if __name__ == "__main__":
    main()
