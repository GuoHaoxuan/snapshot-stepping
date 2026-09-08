"""核验 SVOM/GRM 的 DEAD_TIME 列到底是什么，以及死时间该按哪种模型修正。

决定性检验：若 DEAD_TIME 是该事例造成的死时间，则同一路探测器里相邻事例的间隔不可能
短于它——间隔分布应在死时间处有硬截断，且截断位置随 DEAD_TIME 变化。
"""
from astropy.io import fits
import numpy as np, glob, sys

f = sorted(glob.glob("/gecamfs/SVOM/Archived-DATA/GRM-DATA/L1B/daily/2024/07/03/grm_evt/svom_grm_evt_240703_0*_v*.fits"))[0]
print("文件:", f.split("/")[-1])
with fits.open(f) as h:
    hdr = h[3].header
    for k in sorted(hdr):
        if "DEAD" in k or "EXPO" in k or "LIVE" in k or "TIME" in k.upper()[:4]:
            print("  头:", k, "=", hdr[k])
    for det, i in (("GRD1", 3), ("GRD2", 4), ("GRD3", 5)):
        d = h[i].data
        t = np.asarray(d["TIME"], float); dt_col = np.asarray(d["DEAD_TIME"], float)
        pi = np.asarray(d["PI"]); et = np.asarray(d["EVT_TYPE"]); gain = np.asarray(d["GAIN_TYPE"])
        o = np.argsort(t); t, dt_col, pi, et, gain = t[o], dt_col[o], pi[o], et[o], gain[o]
        gap = np.diff(t) * 1e6      # µs
        print(f"\n=== {det}: {len(t)} 事例，{t[-1]-t[0]:.0f} s")
        print("  DEAD_TIME 取值: 唯一值 %d 个，中位 %.3f，均值 %.3f，min %.3f max %.3f µs" % (len(np.unique(dt_col)), np.median(dt_col), dt_col.mean(), dt_col.min(), dt_col.max()))
        u, c = np.unique(np.round(dt_col, 3), return_counts=True)
        top = np.argsort(-c)[:6]; print("  最常见取值:", [(float(u[k]), int(c[k])) for k in top])
        print("  相邻事例间隔 µs: min %.3f  p0.1 %.3f  p1 %.3f  p5 %.3f  中位 %.1f" % (gap.min(), *np.percentile(gap, [0.1, 1, 5]), np.median(gap)))
        print("  间隔 < 该事例 DEAD_TIME 的比例: %.4f%%（若这一列是死时间，应≈0）" % (100 * np.mean(gap < dt_col[1:])))
        print("  间隔 < 4 µs 的比例: %.4f%%；< 2 µs: %.4f%%" % (100 * np.mean(gap < 4), 100 * np.mean(gap < 2)))
        # 间隔分布的低端形状
        hist, edges = np.histogram(gap[gap < 30], bins=np.arange(0, 30.5, 0.5))
        print("  间隔直方(0–15 µs, 0.5 µs 一格):", hist[:30].tolist())
        # DEAD_TIME 是否随脉冲高度变化
        for lo, hi in ((25, 60), (60, 110), (110, 160), (160, 256)):
            m = (pi >= lo) & (pi < hi) & (et == 0)
            if m.sum() > 1000: print("    PI %3d–%3d: DEAD_TIME 中位 %.3f µs (n=%d)" % (lo, hi, np.median(dt_col[m]), m.sum()))
        for g in np.unique(gain):
            m = gain == g
            if m.sum() > 1000: print("    GAIN_TYPE %d: DEAD_TIME 中位 %.3f µs (n=%d)" % (g, np.median(dt_col[m]), m.sum()))
        # 一秒内 DEAD_TIME 之和 = 死时间占比？
        sec = np.floor(t - t[0]).astype(int); n_sec = sec.max() + 1
        s_dt = np.bincount(sec, weights=dt_col, minlength=n_sec) * 1e-6
        s_n = np.bincount(sec, minlength=n_sec)
        ok = s_n > 100
        print("  逐秒 Σ DEAD_TIME: 中位 %.4f s（即死时间占比 %.2f%%），与该秒计数的相关 r=%.3f" % (np.median(s_dt[ok]), 100*np.median(s_dt[ok]), np.corrcoef(s_dt[ok], s_n[ok])[0,1]))
        # 非瘫痪模型自洽性：观测率 m 与 Σdt 的关系应满足 死时间占比 = m·τ
        tau = np.median(dt_col) * 1e-6
        print("  非瘫痪预期占比 = 观测率×τ = %.4f（实测 %.4f）" % (np.median(s_n[ok]) * tau, np.median(s_dt[ok])))
