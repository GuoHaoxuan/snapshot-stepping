"""看 47904 个 GBM 候选的分组计数构成，以及几种候选判据各能压掉多少。"""
import numpy as np, csv

rows = list(csv.DictReader(open("/scratchfs2/gecam/guohx/gbmrun/cand_detstats2.csv")))
g = lambda k, t=int: np.array([t(r[k]) for r in rows])
cnt, nw = g("count"), g("n_win")
nai, b0, b1 = g("n_nai"), g("n_b0"), g("n_b1")
ndet, ntimes, mult = g("n_det"), g("n_times"), g("max_mult")
dur, fa = g("dur_us", float), g("fa", float)
pha = g("pha_med")
n = len(rows)
print("候选 %d" % n)
print("提取与搜索计数一致: %.1f%%   (n_win-count 的 p10/p50/p90 = %d/%d/%d)"
      % (100 * (nw == cnt).mean(), *np.percentile(nw - cnt, [10, 50, 90]).astype(int)))
print()
print("--- 分组计数 ---")
for name, v in [("NaI 合计", nai), ("b0", b0), ("b1", b1)]:
    print("  %-8s p10=%2d p50=%2d p90=%2d  为0的占 %.1f%%" %
          (name, *np.percentile(v, [10, 50, 90]).astype(int), 100 * (v == 0).mean()))
print("  三组里非空的组数: " + "  ".join(
    "%d组=%.1f%%" % (k, 100 * (((nai > 0).astype(int) + (b0 > 0) + (b1 > 0)) == k).mean())
    for k in (1, 2, 3)))
print()
print("--- 同时性（宇宙线穿整星的直接指标）---")
print("  不同时间戳数 n_times p10=%d p50=%d p90=%d" % tuple(np.percentile(ntimes, [10, 50, 90])))
print("  单一时间戳最大重数 max_mult p50=%d p90=%d max=%d" % (
    np.percentile(mult, 50), np.percentile(mult, 90), mult.max()))
print("  count 全落在一个时间戳上(n_times==1): %.1f%%" % (100 * (ntimes == 1).mean()))
print("  一半以上计数挤在一个时间戳: %.1f%%" % (100 * (mult >= nw / 2).mean()))
print("  点亮>=6 个探头且每个<=2 计数: %.1f%%" % (100 * ((ndet >= 6) & (nw <= 2 * ndet)).mean()))
print()
print("--- 若干判据的留存率（全体 / fa<=1e-5 那一档）---")
sig = fa <= 1e-5
print("  基线                     : %6d  |  %4d" % (n, sig.sum()))
cuts = [
    ("各越线组各自 >=3 计数", (np.sort(np.stack([nai, b0, b1], 1), 1)[:, -2] >= 3)),
    ("各越线组各自 >=5 计数", (np.sort(np.stack([nai, b0, b1], 1), 1)[:, -2] >= 5)),
    ("各越线组各自 >=8 计数", (np.sort(np.stack([nai, b0, b1], 1), 1)[:, -2] >= 8)),
    ("max_mult <= 2         ", mult <= 2),
    ("max_mult <= 3         ", mult <= 3),
    ("n_times >= 4          ", ntimes >= 4),
    ("时长 >= 20us          ", dur >= 20),
    ("时长 >= 50us          ", dur >= 50),
]
for name, m in cuts:
    print("  %-24s : %6d (%5.1f%%) |  %4d (%5.1f%%)" %
          (name, m.sum(), 100 * m.mean(), (m & sig).sum(), 100 * (m & sig).sum() / max(sig.sum(), 1)))
print()
comb = (np.sort(np.stack([nai, b0, b1], 1), 1)[:, -2] >= 3) & (mult <= 3)
print("  组内>=3 且 max_mult<=3   : %6d (%5.1f%%) |  %4d (%5.1f%%)" %
      (comb.sum(), 100 * comb.mean(), (comb & sig).sum(), 100 * (comb & sig).sum() / max(sig.sum(), 1)))
print()
print("--- 最显著的 12 个候选在各判据下的表现 ---")
order = np.argsort(fa)[:12]
print("  %-19s %6s %5s %4s %4s %4s %5s %5s %5s %7s" %
      ("start", "fa", "count", "NaI", "b0", "b1", "ndet", "ntim", "mult", "dur_us"))
for i in order:
    print("  %-19s %6.0e %5d %4d %4d %4d %5d %5d %5d %7.1f" %
          (rows[i]["start"], fa[i], cnt[i], nai[i], b0[i], b1[i], ndet[i], ntimes[i], mult[i], dur[i]))
