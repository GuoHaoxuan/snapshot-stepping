"""放宽窗长上限对灵敏度的理论收益：对不同真实时长 T 的暴发，进目录（fa ≤ 1e-5）所需的最少超出计数。

搜索会试所有 ≤ 上限的窗，所以对某个 T，取各可用窗里门槛最低的那个。
本底 r = 3576 c/s（SVOM v6 显著候选处的中位本底率），min_number = 8 也要满足。
"""
import numpy as np
from scipy.stats import poisson

YEAR = 3600 * 24 * 365.25
R = 3576.0          # c/s
FA = 1e-5           # 目录的直接接受线
MIN_N = 8


def need_counts(W):
    """窗长 W 秒时，窗内需要多少总计数才能达到 fa ≤ FA。"""
    mu = R * W
    p_thr = FA * W / YEAR
    n = max(MIN_N, int(mu))
    while n < 100000:
        if poisson.sf(n - 1, mu) < p_thr:
            return n, mu
        n += 1
    return np.inf, mu


def need_excess(W, T):
    n, mu = need_counts(W)
    frac = min(W, T) / T
    return (n - mu) / frac


print("本底 %.0f c/s，目录线 fa ≤ %.0e，min_number = %d" % (R, FA, MIN_N))
print(" 真实时长 |    上限 1 ms：最优窗 / 所需超出 |    上限 5 ms：最优窗 / 所需超出 | 阈值降低")
grid = np.arange(0.02e-3, 5.001e-3, 0.02e-3)
for T_ms in (0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0):
    T = T_ms * 1e-3
    res = {}
    for cap, tag in ((1e-3, "1ms"), (5e-3, "5ms")):
        ws = grid[grid <= cap]
        e = np.array([need_excess(w, T) for w in ws])
        i = int(np.argmin(e))
        res[tag] = (ws[i] * 1e3, e[i])
    w1, e1 = res["1ms"]; w5, e5 = res["5ms"]
    print("  %5.1f ms |  %.2f ms / %6.1f 计数            |  %.2f ms / %6.1f 计数            | %5.0f%%"
          % (T_ms, w1, e1, w5, e5, 100 * (1 - e5 / e1)))
