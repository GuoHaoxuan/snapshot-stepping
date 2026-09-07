"""决定性检验：事例表的 ANTI_COIN=1 是否就是与带电粒子探测器 GCD 的符合。

HK 文件（1 Hz）里有 GCD1CNT/GCD2CNT/GCD3CNT——GRM 的带电粒子探测器逐秒计数。若 ANTI_COIN
标的是"该事例与 GCD 符合"，那么逐秒的 AC=1 计数应与同一路的 GCDnCNT 强相关，且与本路的
GRDnCNT（总计数）弱相关。逐路配对（GRD1↔GCD1 …）还能验证是不是一一对应。

用法: python3 svom_anticoin_gcd.py <YYYY/MM/DD> <hour_step> <out_prefix>
"""
import glob
import sys
from datetime import datetime

import numpy as np
from astropy.io import fits

ARCHIVE = "/gecamfs/SVOM/Archived-DATA/GRM-DATA/L1B/daily"


def latest(pattern):
    g = sorted(glob.glob(pattern))
    return g[-1] if g else None


def main():
    day, step, prefix = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    ymd = datetime.strptime(day, "%Y/%m/%d").strftime("%y%m%d")
    rows = []          # 每秒一行: [ac1, ac2, ac3, n1, n2, n3, gcd1, gcd2, gcd3]
    for hh in range(0, 24, step):
        evt = latest(f"{ARCHIVE}/{day}/grm_evt/svom_grm_evt_{ymd}_{hh:02d}_v*.fits")
        hk = latest(f"{ARCHIVE}/{day}/hk/svom_hk_{ymd}_{hh:02d}_v*.fits")
        if evt is None or hk is None:
            print("  缺文件", day, hh, flush=True); continue
        with fits.open(evt) as h:
            gti = [(float(a), float(b)) for a, b in zip(h["GTI"].data["START"], h["GTI"].data["STOP"])]
            per_det = []
            for i in (3, 4, 5):
                d = h[i].data
                t = np.asarray(d["TIME"], float); ac = np.asarray(d["ANTI_COIN"]); et = np.asarray(d["EVT_TYPE"])
                m = et == 0
                keep = np.zeros(len(t), bool)
                for a, b in gti:
                    keep |= (t >= a) & (t <= b)
                m &= keep
                per_det.append((t[m], ac[m]))
        with fits.open(hk) as h:
            d = h["HK"].data
            hkt = np.asarray(d["TIME"], float)
            gcd = np.stack([np.asarray(d[f"GCD{k}CNT"], float) for k in (1, 2, 3)], axis=1)
        sec0 = int(np.floor(hkt.min())); sec1 = int(np.ceil(hkt.max()))
        edges = np.arange(sec0, sec1 + 1)
        n_sec = len(edges) - 1
        ac_sec = np.zeros((n_sec, 3)); n_sec_cnt = np.zeros((n_sec, 3))
        for k, (t, ac) in enumerate(per_det):
            if len(t) == 0:
                continue
            n_sec_cnt[:, k], _ = np.histogram(t, bins=edges)
            ac_sec[:, k], _ = np.histogram(t[ac == 1], bins=edges)
        # HK 对齐到秒
        gcd_sec = np.full((n_sec, 3), np.nan)
        idx = np.clip(np.searchsorted(edges[:-1], np.floor(hkt)) , 0, n_sec - 1)
        for k in range(3):
            gcd_sec[idx, k] = gcd[:, k]
        ok = np.isfinite(gcd_sec).all(axis=1) & (n_sec_cnt.sum(axis=1) > 100)
        rows.append(np.hstack([ac_sec[ok], n_sec_cnt[ok], gcd_sec[ok]]))
        print("  %s %02dz: %d 秒" % (day, hh, ok.sum()), flush=True)
    if not rows:
        print("no data"); return
    R = np.vstack(rows)
    ac, cnt, gcd = R[:, 0:3], R[:, 3:6], R[:, 6:9]
    np.savez(prefix + "_gcd.npz", ac=ac, cnt=cnt, gcd=gcd)
    with open(prefix + "_gcd.txt", "w") as f:
        def out(s):
            print(s, flush=True); f.write(s + "\n")
        out("秒数 %d（%s）" % (len(R), day))
        out("逐秒平均: AC=1 %.2f/%.2f/%.2f c/s；总计数 %.0f/%.0f/%.0f c/s；GCD %.1f/%.1f/%.1f c/s" % (
            *ac.mean(axis=0), *cnt.mean(axis=0), *gcd.mean(axis=0)))
        out("\n相关系数矩阵（逐秒，Pearson）：行=各路 AC=1 计数，列=各路 GCD 计数")
        out("            " + "".join("   GCD%d" % (k + 1) for k in range(3)) + "     GRD总计数(同路)")
        for i in range(3):
            cs = [np.corrcoef(ac[:, i], gcd[:, k])[0, 1] for k in range(3)]
            own = np.corrcoef(ac[:, i], cnt[:, i])[0, 1]
            out("   AC(GRD%d)  " % (i + 1) + "".join(" %6.3f" % c for c in cs) + "        %6.3f" % own)
        out("\n对照：各路总计数与 GCD 的相关（若 AC 只是随计数率涨落，这一行也会高）")
        for i in range(3):
            cs = [np.corrcoef(cnt[:, i], gcd[:, k])[0, 1] for k in range(3)]
            out("   N(GRD%d)   " % (i + 1) + "".join(" %6.3f" % c for c in cs))
        out("\n比值 AC=1 / GCD（同路，逐秒中位）：" + " ".join("%.4f" % np.median(ac[:, k] / np.maximum(gcd[:, k], 1)) for k in range(3)))
        # 分箱：GCD 高低两端的 AC 率
        out("\n按同路 GCD 计数分箱（GRD1）：")
        q = np.quantile(gcd[:, 0], [0.1, 0.3, 0.5, 0.7, 0.9])
        for lo, hi in zip([0] + list(q), list(q) + [1e9]):
            m = (gcd[:, 0] >= lo) & (gcd[:, 0] < hi)
            if m.sum() < 20:
                continue
            out("   GCD1 %6.0f-%6.0f c/s: n=%5d  AC=1 %6.2f c/s  总计数 %7.0f c/s  AC 率 %.4f%%" % (
                lo, hi, m.sum(), ac[m, 0].mean(), cnt[m, 0].mean(), 100 * ac[m, 0].sum() / max(cnt[m, 0].sum(), 1)))
    print("wrote", prefix + "_gcd.txt")


if __name__ == "__main__":
    main()
