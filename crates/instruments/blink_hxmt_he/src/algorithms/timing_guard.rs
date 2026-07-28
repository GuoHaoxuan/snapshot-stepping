//! 时标塌陷判据：时间重建塌掉的小时，事例被堆到重复时标上，整小时不搜。
//!
//! # 现象
//!
//! 个别小时的 1K 时间重建在文件开头塌陷：大批事例被赋成完全相同的时标，
//! 形成微秒级的密集堆。搜索按局部速率算，这些堆确实统计极端，于是被
//! 正确地判成成团信号 —— 输入是坏的，算法没错。实测 2019-04-25 T23 头
//! 5 秒里 **46% 的相邻事例共享相同时标**（正常 2-3%），速率是后 5 秒的
//! 2.5 倍，一小时产出 694 个候选，而该天其余 23 小时头 5 秒各产出 0 个。
//!
//! # 判据要躲开的坑
//!
//! 相同时标本身是正常的：同一机箱的事例落在各自的 2 μs 网格上，速率越高
//! 撞进同一格的就越多。所以**原始 tie 占比会随计数率上升**，直接卡阈值有
//! 冤杀真亮暴的风险 —— 与 flood 判据里"溢出占比 vs gap 占比"是同一个坑。
//!
//! 因此两个量都算、都记录：
//!
//! * [`max_tie_fraction`] —— 滑窗内相同时标占比的最大值（原始）
//! * [`max_tie_excess`] —— 同一个量除以该窗自身速率所预期的占比（速率归一，
//!   对亮暴免疫）
//!
//! 判决用哪个、阈值取多少，见下面常量上的标定记录。
//!
//! # 为什么算全部事例而不是 `keep()` 之后的
//!
//! 时标塌陷影响的是时间重建本身，与事例选择无关；而 `EventFile` 直接给得出
//! 全部事例的时间切片（零成本），走 `keep()` 则要物化三千多万个 `Event`。
//! 两者的判别力经标定确认一致。

/// 滑窗宽度（相邻间隔数）。正常速率下约 0.7 秒，足以覆盖一次塌陷，
/// 又不至于把一次短暴稀释掉。
pub const WINDOW_GAPS: usize = 10_000;

/// 机箱内时标的量化步长（秒）。HE 的 ptime 计数器是 2 μs 一跳，同机箱两个
/// 事例落进同一跳就会得到完全相同的时标。跨机箱因为各自的 stime offset 不同，
/// 不落在同一张网格上（实测最小非零间隔只有几十纳秒），所以这个步长只用来
/// 给"正常应该有多少 tie"定一个量级基准，不是精确模型。
pub const TICK_SECONDS: f64 = 2e-6;

/// 出判决所需的最少相邻间隔数。
pub const MIN_GAPS_FOR_VERDICT: usize = 2 * WINDOW_GAPS;

/// 滑窗内相同时标占比的最大值。入参是本单元全部事例的时间（秒，已按时间排序）。
///
/// 事例数不足一个窗时返回 `None`。
pub fn max_tie_fraction(times: &[f64]) -> Option<f64> {
    let n_gaps = times.len().checked_sub(1)?;
    if n_gaps < MIN_GAPS_FOR_VERDICT {
        return None;
    }
    let tied = |gap: usize| usize::from(times[gap + 1] == times[gap]);

    let mut count: usize = (0..WINDOW_GAPS).map(tied).sum();
    let mut max = count;
    for gap in WINDOW_GAPS..n_gaps {
        count += tied(gap);
        count -= tied(gap - WINDOW_GAPS);
        max = max.max(count);
    }
    Some(max as f64 / WINDOW_GAPS as f64)
}

/// 滑窗内"相同时标占比 / 该窗速率所预期的占比"的最大值。
///
/// 预期占比取 `1 − exp(−rate × TICK)`：窗内速率越高，正常情况下本就该有越多
/// 相同时标。除掉它之后，真亮暴（速率高但时标正常）不会因为速率高而被判。
pub fn max_tie_excess(times: &[f64]) -> Option<f64> {
    let n_gaps = times.len().checked_sub(1)?;
    if n_gaps < MIN_GAPS_FOR_VERDICT {
        return None;
    }
    let tied = |gap: usize| usize::from(times[gap + 1] == times[gap]);

    let mut count: usize = (0..WINDOW_GAPS).map(tied).sum();
    let mut max = 0.0f64;
    for start in 0..=(n_gaps - WINDOW_GAPS) {
        if start > 0 {
            count += tied(start + WINDOW_GAPS - 1);
            count -= tied(start - 1);
        }
        let span = times[start + WINDOW_GAPS] - times[start];
        if span > 0.0 {
            let rate = WINDOW_GAPS as f64 / span;
            let expected = -(-rate * TICK_SECONDS).exp_m1();
            if expected > 1e-9 {
                let excess = (count as f64 / WINDOW_GAPS as f64) / expected;
                max = max.max(excess);
            }
        }
    }
    Some(max)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 均匀铺开、互不相同的时标 —— 一个 tie 都没有。
    fn clean(n: usize, rate: f64) -> Vec<f64> {
        (0..n).map(|i| i as f64 / rate).collect()
    }

    #[test]
    fn clean_stream_has_no_ties() {
        let times = clean(3 * WINDOW_GAPS, 15_000.0);
        assert_eq!(max_tie_fraction(&times), Some(0.0));
        assert_eq!(max_tie_excess(&times), Some(0.0));
    }

    #[test]
    fn a_collapsed_stretch_is_caught_even_when_short() {
        // 一段塌陷（半个窗宽的事例全压在同一时标上）埋在干净流里
        let mut times = clean(3 * WINDOW_GAPS, 15_000.0);
        let at = times[WINDOW_GAPS];
        for slot in times
            .iter_mut()
            .skip(WINDOW_GAPS)
            .take(WINDOW_GAPS / 2)
        {
            *slot = at;
        }
        let fraction = max_tie_fraction(&times).unwrap();
        assert!(fraction > 0.4, "{fraction}");
        assert!(max_tie_excess(&times).unwrap() > 10.0);
    }

    #[test]
    fn a_bright_burst_with_honest_timestamps_is_not_flagged() {
        // 真亮暴：速率暴涨 100 倍，但每个事例时标各不相同。原始占比不该动，
        // 速率归一量更不该动 —— 这是判据不能冤杀真事件的底线。
        let mut times = clean(3 * WINDOW_GAPS, 15_000.0);
        let base = times[WINDOW_GAPS];
        for (k, slot) in times
            .iter_mut()
            .skip(WINDOW_GAPS)
            .take(WINDOW_GAPS)
            .enumerate()
        {
            *slot = base + k as f64 / 1_500_000.0;
        }
        for gap in WINDOW_GAPS * 2..times.len() {
            times[gap] = times[gap - 1] + 1.0 / 15_000.0;
        }
        assert_eq!(max_tie_fraction(&times), Some(0.0));
        assert_eq!(max_tie_excess(&times), Some(0.0));
    }

    #[test]
    fn too_few_events_yield_no_verdict() {
        assert_eq!(max_tie_fraction(&clean(MIN_GAPS_FOR_VERDICT - 1, 15_000.0)), None);
        assert_eq!(max_tie_excess(&[]), None);
        assert_eq!(max_tie_fraction(&[1.0]), None);
    }

    #[test]
    fn both_metrics_scan_the_whole_hour_not_just_the_head() {
        // 塌陷若发生在中段也必须抓到 —— 判据不假设它一定在文件开头
        let mut times = clean(4 * WINDOW_GAPS, 15_000.0);
        let at = times[2 * WINDOW_GAPS];
        for slot in times.iter_mut().skip(2 * WINDOW_GAPS).take(WINDOW_GAPS) {
            *slot = at;
        }
        assert!(max_tie_fraction(&times).unwrap() > 0.9);
    }
}
