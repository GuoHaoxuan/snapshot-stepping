# WWLLN 覆盖外的候选与覆盖内是同一人群（2026-09-07）

闪电库止于 2024-12-31，只覆盖 SVOM 观测期的头 193 天（曝光 140 天）；v6 目录里 78% 的候选
无法用闪电验证。这里的检验支撑「把覆盖内的关联率外推到整个目录」这一步。

- `svom_coverage.png`：八个面板，见下面的结论与上一级 `OPEN-QUESTIONS.md` 第 17 条。
- `features_top899.csv`：899 个显著候选（fa ≤ 1e-5）的特征，svomrun6 的 `svom_features.py` 产出。
- `assoc_svom_v6.csv`：同批候选的闪电关联结果（`blink wwlln` 写出的 tgfs.json 提取）。

复现：

    python3 scripts/plot_svom_coverage.py \
        crates/instruments/blink_svom_grm/evidence/coverage/features_top899.csv \
        crates/instruments/blink_svom_grm/evidence/coverage/assoc_svom_v6.csv \
        -o svom_coverage.png

曝光（覆盖内 140.2 天、覆盖外 469.0 天）由脚本的 `--exposure-in/--exposure-out` 传入，
数字取自 svomrun6 逐天账本 `searched_seconds` 的分段求和。
