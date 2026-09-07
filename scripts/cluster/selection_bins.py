"""把 HXMT 与 SVOM 的候选按 fa 分箱，输出逐箱的候选数、闪电关联数、误关联期望，以及各自的曝光。

判选的两层（fa ≤ 1e-5 直接接受；1e-5 < fa ≤ 1 且关联）是否该在两台仪器上通用，要看的就是
这两条曲线：候选数随 fa 的分布（TGF 段是否同样平坦），以及分段关联率相对偶然期望的高低。
只统计 WWLLN 覆盖内的候选，池级去列车与判选一致。

用法: python3 selection_bins.py <out.csv>
"""
import json, glob, sys
import numpy as np

RUNS = (
    ("HXMT/HE", "/scratchfs2/gecam/guohx/v6run/tgfs.json", "/scratchfs2/gecam/guohx/v6run/data/Insight-HXMT_HE/*/*/*_hours.json"),
    ("SVOM/GRM", "/scratchfs2/gecam/guohx/svomrun6/tgfs.json", "/scratchfs2/gecam/guohx/svomrun6/data/SVOM_GRM/*/*/*_hours.json"),
)
EDGES = np.logspace(-60, np.log10(20.0), 121)


def main(out):
    rows = []
    for name, tgfs, hours_glob in RUNS:
        n_all = np.zeros(len(EDGES) - 1); n_assoc = np.zeros_like(n_all); s_prob = np.zeros_like(n_all)
        n_train = 0; n_total = 0
        with open(tgfs) as f:
            data = json.load(f)
        for r in data:
            n_total += 1
            if r["train"].get("is_train"):
                n_train += 1; continue          # 池级去列车，与判选同一口径
            li = r["lightning"]
            if not li.get("in_coverage", True):
                continue
            fa = r["signal"]["false_positive_per_year"]
            i = int(np.clip(np.digitize(fa, EDGES) - 1, 0, len(n_all) - 1))
            n_all[i] += 1
            if li.get("associated"): n_assoc[i] += 1
            s_prob[i] += li.get("coincidence_probability") or 0.0
        # 覆盖内曝光：只数 2025-01-01 之前的天
        secs = 0.0
        for f in glob.glob(hours_glob):
            if f.split("/")[-1][:8] >= "20250101": continue
            secs += json.load(open(f))["searched_seconds"]
        print(f"{name}: 候选 {n_total}，列车 {n_train}，覆盖内 {int(n_all.sum())}，关联 {int(n_assoc.sum())}，曝光 {secs/86400:.1f} d", flush=True)
        for i in range(len(n_all)):
            rows.append((name, EDGES[i], EDGES[i + 1], n_all[i], n_assoc[i], s_prob[i], secs))
    with open(out, "w") as f:
        f.write("instrument,fa_lo,fa_hi,n_all,n_assoc,sum_prob,exposure_s\n")
        for r in rows:
            f.write("%s,%.6e,%.6e,%d,%d,%.6f,%.1f\n" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6]))
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "selection_bins.csv")
