# SVOM/GRM 能阈核对（2026-09-07）

能阈从 ch15（22.5 keV）改为 ch25（42 keV）的依据。图 `svom_threshold.png`，数据：

- `spectra.csv`：逐能道计数——72 个闪电证实 TGF 的候选窗内、未证实显著候选窗内、候选旁 ±1 s 本底、整小时本底（18 h）。
- `thr_scan.csv`：795 个显著候选逐个在 9 档能阈下的窗内计数 S 与本底期望 B。
- `burst.csv`：各道段事例落在"≥5 个/ms"格子里的占比（低道是否成簇）。
- `ebounds.csv`：EBOUNDS 道–keV 对照。

生成：`scripts/cluster/svom_spectra.py`（集群，读 L1B）→ `scripts/plot_svom_threshold.py`。
结论见上一级的 `OPEN-QUESTIONS.md` 第 14 条。
