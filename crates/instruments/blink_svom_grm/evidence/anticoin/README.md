# `ANTI_COIN` 列的语义（2026-09-07）

判定它是星上标定源事例标记、不是带电粒子反符合的依据。结论见上两级的
`OPEN-QUESTIONS.md` 第 3 条。

- `svom_anticoin.png`：四条证据——绝对速率恒定、49–57 keV 线状能谱、与磁纬/SAA 无关、
  候选窗内计数与窗长成正比。
- `ac_0703_summary.txt`：联合分布、逐路触发率、磁纬与 SAA 分箱、三路时间相关。
- `ac_0703_gcd.txt`：与 HK 里带电粒子/标定探测器计数 `GCD1/2/3CNT` 的逐秒相关（|r| ≤ 0.02）。
- `ac_0703_arrays.npz`：作图用数组（磁纬分箱、AC=0/1 的 PI 谱、逐秒速率、时间相关直方图）。

数据取 2024-07-03 的 12 个整点小时（每 2 小时一个），1.7 亿个事例。
生成：`scripts/cluster/svom_anticoin.py`、`scripts/cluster/svom_anticoin_gcd.py`
→ `scripts/plot_svom_anticoin.py`。
