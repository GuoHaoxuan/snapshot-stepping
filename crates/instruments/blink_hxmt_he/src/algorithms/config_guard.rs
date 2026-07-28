//! 非标准配置判据：星上能量阈值被调过的时段，能道解码是错的，整小时不搜。
//!
//! # 背景
//!
//! 星上能量阈值在任务早期被调过三段（2017-09-29 起约 66h、2018-03-18 起约
//! 24h、2019-12-06 起约 20h；2020 年后未再出现）。阈值一动，ADC 的折返点就
//! 不再落在 channel 19 上，而解码规则 [`Event::channel()`](crate::types::Event)
//! 里的 `raw < 20 → +256` 是写死的 —— 这些时段的能道被系统性解错，凭空造出
//! 大量"高能溢出暴"。这些候选全是解码假象，不是天体信号，也不是探测器饱和。
//!
//! # 为什么用直方而不用遥测寄存器 / 静态时间窗表
//!
//! 找具体寄存器失败过一次（曾锁定 `HE_ANBL` 的 `TMY005`，被第三段 episode 的
//! 交叉验证证伪：那只是它自己在 75/76 之间抖）。标准 1K 产品里也没有显式的
//! 阈值字段。而写死一张时间窗表挡不住"以后再调一次配置"。
//!
//! 所以判据直接量数据自己：正常配置下 raw channel 的分布是三段式 ——
//!
//! * `ch 0-12`：真高能事件折返后落下的平坦低尾
//! * `ch 13-19`：**近乎空的 gap**
//! * `ch 20` 起：正常能谱陡然开始（折返边界就在 19/20）
//!
//! 非标准配置把这段 gap 填满（峰堆在 ch 13-14）。所以"gap 里坐了多少事例"
//! 就是判据本身，且随数据自适应。
//!
//! # 实测标定（1K HE-Evt，即搜索实际吃的那份数据）
//!
//! | 小时 | frac(ch15-19) | frac(ch17-19) |
//! |---|---|---|
//! | 正常 2019-12-04T10 | 1.95e-4 | 3.28e-5 |
//! | 正常 2019-12-08T10 | 2.30e-4 | 3.35e-5 |
//! | GRB 221009A 2022-10-09T13 | 8.57e-5 | 1.16e-5 |
//! | **非标准 2019-12-06T08** | **1.28e-1** | **4.91e-2** |
//!
//! 取 ch15-19、阈值 1e-3：比最坏的正常小时高 4.3×，比坏配置低 128×，中间空
//! 两个量级。（ch17-19 分离更大但计数少 6 倍，低统计小时更抖，故取前者。）
//!
//! # 对真亮暴免疫
//!
//! 史上最亮的 GRB 221009A，gap 占比比普通小时还低 —— 真高能光子折返落在
//! ch0-12 的平尾里，根本碰不到 ch15-19 这段 gap。这正是判据用 gap 占比、
//! 而不用"溢出占比 frac(ch0-19)"的原因：后者对 221009A 只剩 17× 余量，
//! 而且余量会随暴的亮度继续收窄，早晚冤杀真事件。
//!
//! # 对部分污染小时的灵敏度
//!
//! 坏配置小时的事例率本身约为正常的 2 倍，所以一小时里只要有约 0.4%
//! （约 15 秒）落在坏配置里，占比就已越过 1e-3。episode 边界那两个小时
//! 不会漏。

/// gap 道段下界（含）。
pub const GAP_CHANNEL_LO: u8 = 15;

/// gap 道段上界（含）。
pub const GAP_CHANNEL_HI: u8 = 19;

/// 判为非标准配置的 gap 占比阈值。
pub const GAP_FRACTION_THRESHOLD: f64 = 1e-3;

/// 出判决所需的最少事例数。
///
/// 正常小时的 gap 占比约 2e-4，要靠泊松涨落假冒到 1e-3 需要向上涨 4.3 倍；
/// 在 1e5 事例上这种涨落的概率已经小到不可能（期望 23 个、需要 100 个）。
/// 低于这个数不出"非标准"判决 —— 而且这一侧是安全的：坏配置小时的事例数
/// 是正常的两倍（3.4e7 → 7.0e7），永远不可能落到低统计这一档里来。
pub const MIN_EVENTS_FOR_VERDICT: usize = 100_000;

/// gap 道段占比。入参是单元内**全部**事例的 raw channel（`Channel` 列原值，
/// 不做 `keep()` 过滤、不做 `<20 → +256` 折返），与标定时的量法一致。
///
/// 事例表为空时返回 `None`。
pub fn gap_fraction(channels: &[u8]) -> Option<f64> {
    if channels.is_empty() {
        return None;
    }
    let n_gap = channels
        .iter()
        .filter(|&&channel| (GAP_CHANNEL_LO..=GAP_CHANNEL_HI).contains(&channel))
        .count();
    Some(n_gap as f64 / channels.len() as f64)
}

/// 该单元是否处于非标准配置。事例数不足以出判决时返回 `false`（见
/// [`MIN_EVENTS_FOR_VERDICT`]）。
pub fn is_nonstandard(channels: &[u8]) -> bool {
    if channels.len() < MIN_EVENTS_FOR_VERDICT {
        return false;
    }
    gap_fraction(channels).is_some_and(|fraction| fraction > GAP_FRACTION_THRESHOLD)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 用逐道计数造一份 channel 序列。下标 0-45 是真实道号，
    /// 下标 46 是"ch46 以上整条谱"的合计（都在 gap 之外，不影响判据）。
    fn from_hist(counts: &[usize; 47]) -> Vec<u8> {
        let mut channels = Vec::new();
        for (channel, &count) in counts.iter().enumerate() {
            channels.extend(std::iter::repeat_n(channel as u8, count));
        }
        channels
    }

    // 实测直方（1K HE-Evt）按 1/100 等比缩小：占比保到 3 位有效数字，
    // 事例数仍高于判决门槛。
    const NORMAL_20191204T10: [usize; 47] = [
        453, 449, 442, 441, 444, 438, 429, 427, 418, 407, 398, 367, 313, 214, 123, 44, 12, 3, 1, 8,
        2214, 7460, 4759, 4620, 4589, 4609, 4543, 4263, 3892, 3482, 3050, 2763, 2448, 2269, 2125,
        2090, 2039, 2027, 2028, 2056, 2097, 2136, 2182, 2279, 2365, 2466, 256754,
    ];
    const FLOOD_20191206T08: [usize; 47] = [
        674, 667, 664, 656, 654, 649, 643, 647, 635, 613, 589, 558, 1429, 63844, 65397, 31776,
        23018, 15642, 10806, 7834, 5989, 5072, 4678, 4686, 4867, 5158, 5299, 5201, 4858, 4379,
        3862, 3490, 3095, 2851, 2686, 2614, 2567, 2558, 2575, 2623, 2685, 2768, 2843, 2998, 3138,
        3334, 373619,
    ];
    const GRB221009A_20221009T13: [usize; 47] = [
        325, 317, 319, 315, 313, 309, 307, 305, 300, 291, 275, 262, 225, 160, 79, 28, 11, 3, 2, 1,
        5936, 17860, 10669, 10739, 10347, 9924, 9537, 9099, 8676, 8210, 7743, 7518, 7141, 7003,
        6866, 6924, 6887, 6919, 7037, 7137, 7305, 7495, 7653, 7843, 7768, 7483, 304708,
    ];

    #[test]
    fn normal_hour_is_standard() {
        let channels = from_hist(&NORMAL_20191204T10);
        let fraction = gap_fraction(&channels).unwrap();
        assert!(!is_nonstandard(&channels));
        // 实测 1.95e-4，离阈值还有 4 倍以上余量
        assert!(fraction < GAP_FRACTION_THRESHOLD / 4.0, "{fraction}");
    }

    #[test]
    fn brightest_grb_is_standard() {
        // 真亮暴不能被冤杀：221009A 的 gap 占比比普通小时还低 —— 真高能光子
        // 折返落在 ch0-12 平尾，碰不到 gap。
        let channels = from_hist(&GRB221009A_20221009T13);
        let fraction = gap_fraction(&channels).unwrap();
        assert!(!is_nonstandard(&channels));
        assert!(fraction < gap_fraction(&from_hist(&NORMAL_20191204T10)).unwrap());
    }

    #[test]
    fn flood_hour_is_nonstandard() {
        let channels = from_hist(&FLOOD_20191206T08);
        let fraction = gap_fraction(&channels).unwrap();
        assert!(is_nonstandard(&channels));
        // 实测 1.28e-1：判据两侧差两个量级以上
        assert!(fraction > 100.0 * GAP_FRACTION_THRESHOLD, "{fraction}");
    }

    #[test]
    fn partially_contaminated_hour_is_caught() {
        // episode 边界那一小时只有一部分落在坏配置里。坏配置的事例率约为正常的
        // 两倍，所以按时间只占 1% 时，事例数约占 2%，照样远在阈值之上。
        let mut channels = from_hist(&NORMAL_20191204T10);
        let flood = from_hist(&FLOOD_20191206T08);
        // step_by 而不是 take：from_hist 按道号升序铺开，取前缀只会拿到 ch0
        channels.extend(flood.iter().step_by(50).copied());
        assert!(is_nonstandard(&channels), "{:?}", gap_fraction(&channels));
    }

    #[test]
    fn low_statistics_never_yields_a_verdict() {
        // 事例数不够就不出判决（这一侧是安全的：坏配置小时事例数只多不少）。
        let channels = vec![GAP_CHANNEL_LO; MIN_EVENTS_FOR_VERDICT - 1];
        assert!(!is_nonstandard(&channels));
        assert_eq!(gap_fraction(&channels), Some(1.0));
    }

    #[test]
    fn empty_hour_has_no_fraction() {
        assert_eq!(gap_fraction(&[]), None);
        assert!(!is_nonstandard(&[]));
    }
}
