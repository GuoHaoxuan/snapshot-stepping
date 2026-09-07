# Fermi/GBM 第二版 TGF 目录（外部真值）

来源：https://fermi.gsfc.nasa.gov/ssc/data/access/gbm/tgf/ （2026-09-07 下载）。
覆盖 2008-07-11 至 **2016-07-31**，之后没有发表的续编。

| 文件 | 内容 |
|---|---|
| `gbm_tgf_catalog_offline.csv` | 离线搜索表 4135 个（`OS_ID, MET, file, BGO_0_N, BGO_1_N, NAI_N, Date, UTC, Width_ms, P2, Lon, Lat, Alt, LST, TRIG_ID`） |
| `gbm_tgf_catalog_trig.csv` | 星上触发表 686 个 |
| `gbm_tgf_catalog_wwlln.csv` | WWLLN 关联 1544 条（含分离时间与距离） |
| `gbm_tgf_catalog_teb.csv` | 地球电子束 30 个 |

逐年（离线表）：2010 年 241、2011 年 492、2012 年 651、2013 年 779、2014 年 811、2015 年 725、
2016 年（到 7 月）405。

用途：本流水线在 GBM 上的**完备性与纯度检验**——这是 HXMT 与 SVOM 都给不了的外部真值。
MET 与 TTE 文件同一时间基准，可直接按 MET 匹配。**必须跑 2016-07 之前的年份**；
我们先跑的 2019 年不在目录覆盖内，只能作规模演示。
