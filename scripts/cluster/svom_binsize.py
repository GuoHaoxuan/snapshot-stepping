"""SVOM v6：搜索窗 bin_size_best 是否顶在 1 ms 的上限上（上限截断的直接证据）。"""
import json
import numpy as np

recs = json.load(open("scratch/tgfs_svom_v6.json"))
fa = np.array([r["signal"]["false_positive_per_year"] for r in recs])
best = np.array([r["signal"]["bin_size_best"] for r in recs]) * 1e3      # ms
bmax = np.array([r["signal"]["bin_size_max"] for r in recs]) * 1e3
span = np.array([(r["signal"]["bin_size_min"], r["signal"]["bin_size_max"]) for r in recs])
cov = np.array([bool(r["lightning"].get("in_coverage", True)) for r in recs])
asc = np.array([bool(r["lightning"].get("associated")) for r in recs])
sig = fa <= 1e-5
CAP = 1.0

print("候选 %d，显著 %d，覆盖内显著 %d，证实 %d" % (len(recs), sig.sum(), (sig & cov).sum(), (sig & cov & asc).sum()))
for name, m in (("证实 TGF", sig & cov & asc), ("覆盖内未证实", sig & cov & ~asc),
                ("覆盖外显著", sig & ~cov), ("全部显著", sig), ("全部候选", np.ones(len(recs), bool))):
    b = best[m]
    print("\n%s n=%d" % (name, m.sum()))
    print("   bin_size_best ms: min %.4f p25 %.4f 中位 %.4f p75 %.4f p90 %.4f max %.4f"
          % (b.min(), np.percentile(b, 25), np.median(b), np.percentile(b, 75), np.percentile(b, 90), b.max()))
    for thr in (0.90, 0.95, 0.99, 0.999):
        print("      >=%.3f ms（上限的 %.0f%%）: %d 个 (%.1f%%)" % (thr, thr * 100, (b >= thr).sum(), 100 * (b >= thr).mean()))
    bm = bmax[m]
    print("   bin_size_max（合并里最长的触发窗）>=0.99 ms: %d 个 (%.1f%%)" % ((bm >= 0.99).sum(), 100 * (bm >= 0.99).mean()))

# 均匀分布假设下的对照：若窗长不受上限影响，0.9-1.0 区间应只占 ~10%
b = best[sig & cov & asc]
h, e = np.histogram(b, bins=np.arange(0, 1.05, 0.1))
print("\n证实 TGF 的 bin_size_best 直方（0.1 ms 一格，0→1.0）:", h.tolist())
b2 = best[sig]
h2, _ = np.histogram(b2, bins=np.arange(0, 1.05, 0.1))
print("全部显著候选的 bin_size_best 直方:", h2.tolist())
