//! 重复事例判据：1K 事例表里整行重复写入的小时，整小时不搜。
//!
//! # 现象
//!
//! 个别小时的 1K 事例表把事例整行写了不止一份 —— 时间、探头、能道、脉宽、
//! 类型、ACD 全都逐字段相同。实测：
//!
//! * 2019-04-25 T23 头 4 秒，一部分事例各写了 **4 份**；
//! * 2022-01-07 T02 全小时，**每个事例都写了 2 份**（重复行占比精确的 0.5）；
//! * 2020-10-18 T13 部分时段写到 3 份。
//!
//! # 为什么必须整小时排除
//!
//! 重复把暴发和本底一起放大，而泊松显著性不是尺度不变的：`count=8` 压在
//! `本底=1` 上的 sf 是 1.1e-6，两边都翻倍成 `16 / 2` 就成了 1.3e-11 ——
//! 凭空多出五个量级。所以整段重复会系统性地抬高该时段所有候选的显著性，
//! 造出一批看起来极显著、其实不存在的事件。2019-04-25 T23 那 4 秒就产出了
//! 694 个候选，而当天其余 23 小时在同样位置各产出 0 个。
//!
//! # 判据
//!
//! 判据直接量病因：**逐字段完全相同的重复行**。正常小时是精确的 0，不是
//! "很小的数" —— 59 个随机抽样小时里 58 个一行重复都没有。
//!
//! 早前试过用"相同时标占比"来判，失败了：相同时标是**症状**，真实高速率
//! 时段本来就会有大量事例落进同一个 2 μs 计数格，两个真实 REP 小时的相同
//! 时标占比比坏小时还高。重复行是**病因**，真事件产生不了逐字段相同的副本。
//!
//! # 免疫性（实测，全部为 0 重复行）
//!
//! | 小时 | 重复行占比 |
//! |---|---|
//! | GRB 221009A（史上最亮）| 0 |
//! | 真暴 2024-08-13 T13 / 2025-11-16 T07 | 0 |
//! | flood E1 / E3 | 0 |
//! | REP 风暴日 2024-05-11 T10（Gannon）| 0 |
//! | 随机抽样 59 个小时中的 58 个 | 0 |
//!
//! 判据不会因为"亮"或"事例多"而误判 —— 它问的是"同一份数据是不是被写了
//! 两遍"，与物理无关。

/// 判为重复写入的重复行占比阈值。
///
/// 正常小时是精确的 0；实测到的坏小时最低是 0.0024。取 1e-4 —— 比正常态高
/// 出一个不可能由涨落跨过的距离（正常态没有涨落，就是 0），又比最轻的坏小时
/// 低 24 倍。取一个非零值而不是"只要有一行重复就排除"，是留一点余量给
/// 未来可能出现的、真正偶发的同格碰撞。
pub const DUPLICATE_FRACTION_THRESHOLD: f64 = 1e-4;

/// 出判决所需的最少事例数。
pub const MIN_EVENTS_FOR_VERDICT: usize = 100_000;

/// 把一个事例的身份打包成一个整数，便于在 run 内排序查重。ACD 不参与：
/// 它是 18 位数组，取它要额外一份拷贝，而实测重复行连 ACD 都相同，
/// 这四个字段已经足够定身份。
fn identity(detector: u8, channel: u8, pulse_width: u8, event_type: u8) -> u32 {
    (detector as u32) << 24 | (channel as u32) << 16 | (pulse_width as u32) << 8 | event_type as u32
}

/// 逐字段完全相同的重复行占全部事例的比例。
///
/// 入参按时间升序。重复行必然共享时标，所以只需在每段等时标的 run 内部查重；
/// run 通常只有几个成员，整体是一遍线性扫描。
///
/// 事例数不足以出判决时返回 `None`。
pub fn duplicate_fraction(
    times: &[f64],
    detectors: &[u8],
    channels: &[u8],
    pulse_widths: &[u8],
    event_types: &[u8],
) -> Option<f64> {
    let n = times.len();
    if n < MIN_EVENTS_FOR_VERDICT
        || detectors.len() != n
        || channels.len() != n
        || pulse_widths.len() != n
        || event_types.len() != n
    {
        return None;
    }

    let mut duplicates = 0usize;
    let mut run_start = 0usize;
    // run 内排序查重。用排序而不是逐个回看，是为了让代价随 run 长度是
    // O(k log k) 而不是 O(k²) —— 实测最长的 run 只有几个成员，但那是数据
    // 碰巧，不该拿它当边界保证。
    let mut run: Vec<u32> = Vec::new();
    for i in 1..=n {
        if i < n && times[i] == times[run_start] {
            continue;
        }
        // [run_start, i) 是一段时标相同的事例
        if i - run_start > 1 {
            run.clear();
            run.extend((run_start..i).map(|slot| {
                identity(
                    detectors[slot],
                    channels[slot],
                    pulse_widths[slot],
                    event_types[slot],
                )
            }));
            run.sort_unstable();
            duplicates += run.windows(2).filter(|pair| pair[0] == pair[1]).count();
        }
        run_start = i;
    }

    Some(duplicates as f64 / n as f64)
}

/// 该单元的事例表是否含有重复写入。
pub fn has_duplicated_events(
    times: &[f64],
    detectors: &[u8],
    channels: &[u8],
    pulse_widths: &[u8],
    event_types: &[u8],
) -> bool {
    duplicate_fraction(times, detectors, channels, pulse_widths, event_types)
        .is_some_and(|fraction| fraction > DUPLICATE_FRACTION_THRESHOLD)
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Table {
        times: Vec<f64>,
        detectors: Vec<u8>,
        channels: Vec<u8>,
        pulse_widths: Vec<u8>,
        event_types: Vec<u8>,
    }

    impl Table {
        fn fraction(&self) -> Option<f64> {
            duplicate_fraction(
                &self.times,
                &self.detectors,
                &self.channels,
                &self.pulse_widths,
                &self.event_types,
            )
        }

        fn flagged(&self) -> bool {
            has_duplicated_events(
                &self.times,
                &self.detectors,
                &self.channels,
                &self.pulse_widths,
                &self.event_types,
            )
        }
    }

    /// n 个互不相同的事例，时标每 2 μs 一个（同格里放两个不同探头的事例，
    /// 模拟真实数据里大量共享时标但并不重复的情形）。
    fn honest(n: usize) -> Table {
        let mut t = Table {
            times: Vec::new(),
            detectors: Vec::new(),
            channels: Vec::new(),
            pulse_widths: Vec::new(),
            event_types: Vec::new(),
        };
        for i in 0..n {
            t.times.push((i / 2) as f64 * 2e-6);
            t.detectors.push((i % 18) as u8);
            t.channels.push((i % 251) as u8);
            t.pulse_widths.push(90 + (i % 7) as u8);
            t.event_types.push(0);
        }
        t
    }

    /// 把 [from, to) 这一段每个事例再写 copies-1 遍。
    fn duplicate_stretch(table: &Table, from: usize, to: usize, copies: usize) -> Table {
        let mut out = Table {
            times: Vec::new(),
            detectors: Vec::new(),
            channels: Vec::new(),
            pulse_widths: Vec::new(),
            event_types: Vec::new(),
        };
        for i in 0..table.times.len() {
            let n = if i >= from && i < to { copies } else { 1 };
            for _ in 0..n {
                out.times.push(table.times[i]);
                out.detectors.push(table.detectors[i]);
                out.channels.push(table.channels[i]);
                out.pulse_widths.push(table.pulse_widths[i]);
                out.event_types.push(table.event_types[i]);
            }
        }
        out
    }

    #[test]
    fn an_honest_table_has_no_duplicates_even_when_timestamps_collide() {
        // 共享时标本身不是罪 —— 真实高速率时段大量事例落在同一个计数格里
        let table = honest(200_000);
        assert_eq!(table.fraction(), Some(0.0));
        assert!(!table.flagged());
    }

    #[test]
    fn a_wholly_doubled_table_is_flagged() {
        // 2022-01-07T02 的形状：每个事例整整两份
        let table = duplicate_stretch(&honest(150_000), 0, 150_000, 2);
        let fraction = table.fraction().unwrap();
        assert!((fraction - 0.5).abs() < 1e-9, "{fraction}");
        assert!(table.flagged());
    }

    #[test]
    fn a_short_quadrupled_stretch_is_flagged() {
        // 2019-04-25T23 的形状：只有开头一小段各写 4 份，占全小时不到 0.3%
        let table = duplicate_stretch(&honest(200_000), 0, 3_000, 4);
        let fraction = table.fraction().unwrap();
        assert!(fraction > DUPLICATE_FRACTION_THRESHOLD, "{fraction}");
        assert!(fraction < 0.05, "{fraction}");
        assert!(table.flagged());
    }

    #[test]
    fn a_stretch_in_the_middle_is_flagged_too() {
        // 判据不假设重复一定在文件开头
        let table = duplicate_stretch(&honest(200_000), 120_000, 123_000, 3);
        assert!(table.flagged());
    }

    #[test]
    fn a_handful_of_duplicates_stays_below_the_threshold() {
        // 阈值留的余量：偶发几行重复不足以判掉一整个小时
        let table = duplicate_stretch(&honest(200_000), 0, 5, 2);
        let fraction = table.fraction().unwrap();
        assert!(fraction > 0.0);
        assert!(!table.flagged(), "{fraction}");
    }

    #[test]
    fn too_few_events_yield_no_verdict() {
        let table = honest(MIN_EVENTS_FOR_VERDICT - 1);
        assert_eq!(table.fraction(), None);
        assert!(!table.flagged());
    }

    #[test]
    fn mismatched_column_lengths_yield_no_verdict() {
        let mut table = honest(200_000);
        table.channels.pop();
        assert_eq!(table.fraction(), None);
    }
}
