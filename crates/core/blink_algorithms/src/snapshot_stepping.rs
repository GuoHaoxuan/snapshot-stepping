use crate::{constants::DAYS_PER_YEAR, types::candidate::Candidate};
use blink_core::{traits::Event, types::MissionElapsedTime};
use statrs::distribution::{DiscreteCDF, Poisson};
use std::cmp::Ordering;
use statrs::function::gamma::ln_gamma;
use uom::si::f64::*;

pub struct SearchConfig {
    pub min_duration: Time,
    pub max_duration: Time,
    pub neighbor: Time,
    pub hollow: Time,
    pub false_positive_per_year: f64,
    pub min_number: u32,
    /// 判定一次触发需要多少个组各自越线（符合数）。
    ///
    /// 单组仪器填 1，此时整套判据精确退化回"这一组自己显著就算"。多组时
    /// 要求 N 组同时越线，是为了挡掉只点亮一个组的东西——GBM 上实测一天
    /// 12.6 万个候选里有 72.6% 是宇宙线穿过整星，每个探头各留一个高能沉积，
    /// 12 个 NaI 合起来轻易凑够 min_number，但两个 BGO 各自只有一个，
    /// 要求 2 组符合就能把它们挡在外面。
    pub coincidence: usize,
}

impl Default for SearchConfig {
    fn default() -> Self {
        Self {
            min_duration: Time::new::<uom::si::time::microsecond>(10.0),
            max_duration: Time::new::<uom::si::time::millisecond>(1.0),
            neighbor: Time::new::<uom::si::time::second>(1.0),
            hollow: Time::new::<uom::si::time::millisecond>(10.0),
            false_positive_per_year: 20.0,
            min_number: 8,
            coincidence: 1,
        }
    }
}

/// 从 n 个里取 k 个的组合数。组数最多十几，不会溢出。
fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut result = 1.0_f64;
    for i in 0..k {
        result = result * (n - i) as f64 / (i + 1) as f64;
    }
    result
}

// fn coincidence_prob(probs: &[f64], n: usize) -> f64 {
//     let mut cache = vec![0.0; n + 1];
//     cache[0] = 1.0;

//     for m_i in 1..=probs.len() {
//         for n_i in (0..=n).rev() {
//             cache[n_i] = match (m_i, n_i) {
//                 (_, 0) => 1.0,
//                 // The following line can be removed, because cache[n_i] is 0.0 initially,
//                 // although it is meaningless mathematically
//                 // (m_i, n_i) if m_i == n_i => probs[m_i - 1] * cache[n_i - 1],
//                 (m_i, n_i) => probs[m_i - 1] * cache[n_i - 1] + (1.0 - probs[m_i - 1]) * cache[n_i],
//             }
//         }
//     }

//     cache[n]
// }

pub fn poisson_isf(p: f64, lambda: f64) -> u32 {
    let mut k = 0;
    let mut cumulative_prob = (-lambda).exp();
    let mut part = 0.0;

    while cumulative_prob < 1.0 - p {
        k += 1;
        part += (lambda / k as f64).ln();
        cumulative_prob += (-lambda + part).exp();
    }

    k
}

/// `gti`：活时间段，按时间排好、互不重叠。本底窗（候选两侧各 `neighbor/2`）
/// 与空窗（各 `hollow/2`）都夹到候选所在的那一段里再算时长，所以本底率的分母
/// 只数活时间。传 `[[start, stop]]` 一段即等价于原先按 chunk 边界夹取的行为，
/// 算术逐位相同。
///
/// 不这样做的后果实测过：SVOM 在 SAA 停机前是 130–210 kc/s 的平台后硬切断，
/// 停机前半秒内候选的本底窗一半伸进空区，分子少一半、分母按墙钟不变，期望从
/// 约 100 压到 48，普通计数就成了 fa=1e-10——一天冒出 89 个这样的假候选。
/// 分子不必动：事例流本来就该按 GTI 过滤，计数只来自活时间。
pub fn search_new<E: Event>(
    data: &[E],
    group_number: usize,
    start: MissionElapsedTime<E::Instrument>,
    stop: MissionElapsedTime<E::Instrument>,
    gti: &[[MissionElapsedTime<E::Instrument>; 2]],
    config: SearchConfig,
) -> Vec<Candidate<E::Instrument>> {
    let mut result: Vec<Candidate<E::Instrument>> = Vec::new();
    // let mut cache = vec![
    //     vec![None; CACHE_COUNT_MAX as usize];
    //     (CACHE_MEAN_MAX * CACHE_MEAN_HASH_FACTOR).ceil() as usize
    // ];
    // let mut cache = vec![0; 100_000];

    let mut cursor = data
        .binary_search_by(|event| event.time().cmp(&start))
        .unwrap_or_else(|index| index);
    if cursor == data.len() {
        return result;
    }

    // 两个快照窗的初始化与下面主循环里的维护（见 loop 尾部）用同一种写法：
    // 判 `[idx + 1]` 是否还在窗内，再前进并计数。早先这里判的是 `[idx]`、
    // 前进后才计数，于是（a）会把窗外的那一个事例也计进本底，（b）末元素仍
    // 满足时间条件时会索引 `data[len]` 越界 —— 只有当整段数据比半个窗还短
    // 才会踩到，一小时的量踩不到，但那是靠数据量挡着、不是靠边界检查挡着。
    let mut mean_start_snapshot = cursor;
    let mut mean_stop_snapshot = cursor;
    let mut mean_numbers_snapshot: Vec<u32> = vec![0; group_number];
    mean_numbers_snapshot[data[cursor].group() as usize] = 1;
    while mean_stop_snapshot + 1 < data.len()
        && data[mean_stop_snapshot + 1].time() - data[cursor].time() < config.neighbor / 2.0
    {
        mean_stop_snapshot += 1;
        mean_numbers_snapshot[data[mean_stop_snapshot].group() as usize] += 1;
    }

    let mut hollow_start_snapshot = cursor;
    let mut hollow_stop_snapshot = cursor;
    let mut hollow_numbers_snapshot: Vec<u32> = vec![0; group_number];
    hollow_numbers_snapshot[data[cursor].group() as usize] = 1;
    while hollow_stop_snapshot + 1 < data.len()
        && data[hollow_stop_snapshot + 1].time() - data[cursor].time() < config.hollow / 2.0
    {
        hollow_stop_snapshot += 1;
        hollow_numbers_snapshot[data[hollow_stop_snapshot].group() as usize] += 1;
    }

    // 三个逐组计数器在循环外分配一次、每轮重置。它们原先是每个 cursor 位置
    // 现做一个 vec![] 加两次 clone —— 外层循环有多少个事例就分配多少次三连，
    // GBM 一小时是三千万个事例，即九千万次 malloc/free，全是白付的。
    let mut numbers: Vec<u32> = vec![0; group_number];
    let mut mean_numbers: Vec<u32> = vec![0; group_number];
    let mut hollow_numbers: Vec<u32> = vec![0; group_number];
    // 符合判据要挑第 N 小的组 p 值，这块缓冲同样只分配一次。
    let mut group_fps: Vec<f64> = vec![0.0; group_number];

    // cursor 所在的活时间段。段按时间排好、cursor 单调前进，下标只增不减；
    // 每个 cursor 位置只算一次，内层循环里的四次夹取仍是原来那四次比较。
    // 候选窗最长 max_duration，远短于任何缺口，不会跨段。cursor 落在段外
    // （事例没按 GTI 过滤时会发生）就退回 chunk 边界，即原先的行为。
    let mut segment = 0usize;

    loop {
        let cursor_time = data[cursor].time();
        while segment + 1 < gti.len() && gti[segment + 1][0] <= cursor_time {
            segment += 1;
        }
        let (live_start, live_stop) = match gti.get(segment) {
            Some(seg) if seg[0] <= cursor_time && cursor_time <= seg[1] => {
                (seg[0].max(start), seg[1].min(stop))
            }
            _ => (start, stop),
        };

        let mut step = 0;
        numbers.fill(0);
        numbers[data[cursor].group() as usize] = 1;
        let mut mean_stop = mean_stop_snapshot;
        mean_numbers.copy_from_slice(&mean_numbers_snapshot);
        let mut hollow_stop = hollow_stop_snapshot;
        hollow_numbers.copy_from_slice(&hollow_numbers_snapshot);

        loop {
            // 窗内事例数恒等于 step + 1：进循环时记了 1 个，之后每前进一步
            // 正好收一个事例进某个组。原先每轮都把逐组计数器求和一遍，而内层
            // 循环 GBM 一小时要跑六亿次。
            let total_number = step as u32 + 1;
            debug_assert_eq!(total_number, numbers.iter().sum::<u32>());
            let duration = data[cursor + step].time() - data[cursor].time();
            if total_number >= config.min_number && duration >= config.min_duration {
                let mean_start_time =
                    (data[cursor].time() - config.neighbor / 2.0).max(live_start);
                let mean_stop_time =
                    (data[cursor + step].time() + config.neighbor / 2.0).min(live_stop);
                let hollow_start_time =
                    (data[cursor].time() - config.hollow / 2.0).max(live_start);
                let hollow_stop_time =
                    (data[cursor + step].time() + config.hollow / 2.0).min(live_stop);
                let pure_mean_duration =
                    (mean_stop_time - mean_start_time) - (hollow_stop_time - hollow_start_time);
                let pure_mean_percent =
                    (duration / pure_mean_duration).get::<uom::si::ratio::ratio>();
                let threshold = config.false_positive_per_year
                    / (uom::si::f64::Time::new::<uom::si::time::second>(3600.0)
                        * 24.0
                        * DAYS_PER_YEAR
                        / duration)
                        .get::<uom::si::ratio::ratio>();
                // 要求 n 组里有 N 组各自越线，零假设下的虚警率是 C(n,N)·alpha^N，
                // 令它等于 threshold 便反解出单组该用的门槛。N=1 时就是
                // threshold/n，即 Bonferroni；剪枝也按这个门槛来。
                let n_required = config.coincidence.max(1).min(group_number);
                let combinations = binomial(group_number, n_required);
                let group_threshold =
                    (threshold / combinations).powf(1.0 / n_required as f64);
                let ln_group_threshold = group_threshold.ln();
                // 逐组算尾概率并直接取最小值。这里曾经先 collect 成 Vec 再取 min，
                // 而这段在内层循环里 —— 每个 (cursor, step) 组合都要在堆上分配再
                // 释放一次，GBM 一天这种组合是几十亿次量级。fold 的次序与原先
                // 逐个取 min 完全一致，数值不变。
                let group_fp = |group: usize| -> f64 {
                    let pure_mean_number = mean_numbers[group] - hollow_numbers[group];
                    let equivalent_background_number =
                        pure_mean_number as f64 * pure_mean_percent;
                        match (equivalent_background_number, numbers[group]) {
                            (0.0, 0) => 1.0,
                            (0.0, _) => 1.0,
                            _ => {
                                let lambda = equivalent_background_number;
                                let count = numbers[group] as f64;
                                // 剪枝，严格等价，不是近似：
                                //   sf(count) = P(X > count) >= P(X = count+1)
                                // 所以只要单项 PMF 已经不小于本组门槛，这一组的
                                // 尾概率也不小于门槛，越不了线。此时返回 +inf 让它
                                // 退出 min 的竞争 —— 真正触发时取到的最小值必然来自
                                // 没被剪的组，那个值是精确算出来的，候选记录的显著性
                                // 因此不受影响。省掉的全是注定不触发的不完全伽马。
                                let ln_pmf = -lambda + (count + 1.0) * lambda.ln()
                                    - ln_gamma(count + 2.0);
                                if ln_pmf >= ln_group_threshold {
                                    f64::INFINITY
                                } else {
                                    Poisson::new(lambda).unwrap().sf(numbers[group] as u64)
                                }
                            }
                        }
                };
                // 分组判据的 Bonferroni 校正：每个窗口按组各做一次检验，就
                // 有几次中奖机会，纯本底下"至少一组超阈"的概率约为单组的
                // group_number 倍。取最显著的那组、乘以组数，误报率就回到
                // false_positive_per_year 标称的水平，`fa` 的语义不变。
                //
                // 单组时是 fps[0] * 1.0，浮点上精确等于原值 —— HXMT 与 SVOM
                // 逐位不受影响。
                // N=1 走原来的路：只要最小的那个，不必排序，数值与从前逐位相同。
                let fp = if n_required == 1 {
                    (0..group_number).map(group_fp).fold(f64::INFINITY, f64::min)
                        * group_number as f64
                } else {
                    // 需要第 N 小的组 p 值。组数最多十几个，就地插入排序最省。
                    for group in 0..group_number {
                        group_fps[group] = group_fp(group);
                    }
                    let slice = &mut group_fps[..group_number];
                    slice.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
                    combinations * slice[n_required - 1].powi(n_required as i32)
                };

                if fp < threshold {
                    let total_equivalent_background_number = (0..group_number)
                        .map(|group| mean_numbers[group] - hollow_numbers[group])
                        .sum::<u32>()
                        as f64
                        * pure_mean_percent;
                    // println!(
                    //     "Found trigger: total_number: {}, equivalent_background_number: {}, fp: {}, threshold: {}, duration: {}",
                    //     total_number,
                    //     total_equivalent_background_number,
                    //     fp,
                    //     threshold,
                    //     duration.to_seconds() * 1e6
                    // );
                    let current = Candidate::new(
                        data[cursor].time(),
                        data[cursor + step].time(),
                        total_number,
                        total_equivalent_background_number,
                        fp,
                    );
                    if let Some(last) = result.last_mut() {
                        if last.mergeable(&current, 0.0) {
                            *last = last.merge(&current);
                        } else {
                            result.push(current);
                        }
                    } else {
                        result.push(current);
                    }
                }
            }

            step += 1;
            if cursor + step >= data.len()
                || data[cursor + step].time() - data[cursor].time() >= config.max_duration
                || data[cursor + step].time() >= stop
            {
                break;
            }
            numbers[data[cursor + step].group() as usize] += 1;
            while mean_stop + 1 < data.len()
                && data[mean_stop + 1].time() - data[cursor + step].time() < config.neighbor / 2.0
            {
                mean_stop += 1;
                mean_numbers[data[mean_stop].group() as usize] += 1;
            }
            while hollow_stop + 1 < data.len()
                && data[hollow_stop + 1].time() - data[cursor + step].time() < config.hollow / 2.0
            {
                hollow_stop += 1;
                hollow_numbers[data[hollow_stop].group() as usize] += 1;
            }
        }

        cursor += 1;
        if cursor >= data.len() || data[cursor].time() >= stop {
            break;
        }
        while mean_start_snapshot + 1 < data.len()
            && data[cursor].time() - data[mean_start_snapshot + 1].time() > config.neighbor / 2.0
        {
            mean_numbers_snapshot[data[mean_start_snapshot].group() as usize] -= 1;
            mean_start_snapshot += 1;
        }
        while mean_stop_snapshot + 1 < data.len()
            && data[mean_stop_snapshot + 1].time() - data[cursor].time() < config.neighbor / 2.0
        {
            mean_stop_snapshot += 1;
            mean_numbers_snapshot[data[mean_stop_snapshot].group() as usize] += 1;
        }
        while hollow_start_snapshot + 1 < data.len()
            && data[cursor].time() - data[hollow_start_snapshot + 1].time() > config.hollow / 2.0
        {
            hollow_numbers_snapshot[data[hollow_start_snapshot].group() as usize] -= 1;
            hollow_start_snapshot += 1;
        }
        while hollow_stop_snapshot + 1 < data.len()
            && data[hollow_stop_snapshot + 1].time() - data[cursor].time() < config.hollow / 2.0
        {
            hollow_stop_snapshot += 1;
            hollow_numbers_snapshot[data[hollow_stop_snapshot].group() as usize] += 1;
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use blink_core::{
        error::Error,
        traits::{Chunk, Instrument},
        types::{Coverage, Signal},
    };
    use chrono::{DateTime, NaiveDate, TimeZone, Utc};
    use serde::Serialize;
    use std::sync::OnceLock;

    #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
    struct TestInstrument;

    struct TestChunk;

    impl Chunk for TestChunk {
        type Event = TestEvent;

        fn from_epoch(_: &DateTime<Utc>) -> Result<Self, Error> {
            Err(Error::Unknown)
        }
        fn search(&self) -> Vec<Signal<Self::Event>> {
            Vec::new()
        }
        fn last_modified(_: &DateTime<Utc>) -> Result<DateTime<Utc>, Error> {
            Err(Error::Unknown)
        }
        fn coverage(&self) -> Coverage {
            Coverage {
                span_seconds: 0.0,
                masked_seconds: 0.0,
            }
        }
    }

    impl Instrument for TestInstrument {
        type Chunk = TestChunk;

        fn ref_time() -> &'static DateTime<Utc> {
            static REF_TIME: OnceLock<DateTime<Utc>> = OnceLock::new();
            REF_TIME.get_or_init(|| Utc.with_ymd_and_hms(2012, 1, 1, 0, 0, 0).unwrap())
        }
        fn launch_day() -> NaiveDate {
            NaiveDate::from_ymd_opt(2017, 6, 15).unwrap()
        }
        fn name() -> &'static str {
            "test"
        }
    }

    #[derive(Serialize, Debug, Clone)]
    struct TestEvent {
        seconds: f64,
        group: u8,
    }

    impl Event for TestEvent {
        type Instrument = TestInstrument;
        type ChannelType = u16;

        fn time(&self) -> MissionElapsedTime<TestInstrument> {
            MissionElapsedTime::new(self.seconds)
        }
        fn channel(&self) -> u16 {
            100
        }
        fn group(&self) -> u8 {
            self.group
        }
        fn keep(&self) -> bool {
            true
        }
    }

    fn at(seconds: &[f64]) -> Vec<TestEvent> {
        seconds
            .iter()
            .map(|&s| TestEvent {
                seconds: s,
                group: 0,
            })
            .collect()
    }

    /// 本底 200 个事例铺在 20 s 上，再在 10.0 s 处塞 30 个事例进 100 us。
    fn burst_on_quiet_background(group: u8) -> Vec<TestEvent> {
        let mut seconds: Vec<f64> = (0..200).map(|i| i as f64 * 0.1).collect();
        seconds.extend((0..30).map(|i| 10.0 + i as f64 * 3e-6));
        seconds.sort_by(|a, b| a.partial_cmp(b).unwrap());
        seconds
            .iter()
            .map(|&s| TestEvent { seconds: s, group })
            .collect()
    }

    fn run_grouped(data: &[TestEvent], group_number: usize) -> Vec<crate::types::candidate::Candidate<TestInstrument>> {
        run_coincident(data, group_number, 1)
    }

    fn run_coincident(
        data: &[TestEvent],
        group_number: usize,
        coincidence: usize,
    ) -> Vec<crate::types::candidate::Candidate<TestInstrument>> {
        search_new(
            data,
            group_number,
            MissionElapsedTime::<TestInstrument>::new(0.0),
            MissionElapsedTime::<TestInstrument>::new(3600.0),
            &whole_span(),
            SearchConfig {
                coincidence,
                ..SearchConfig::default()
            },
        )
    }

    /// 本底铺在两组上，再在 10.0 s 处让指定的那些组各自出现一次暴发。
    fn burst_in_groups(groups: &[u8]) -> Vec<TestEvent> {
        let mut events: Vec<TestEvent> = Vec::new();
        for group in 0..2u8 {
            // 两组本底错开半个间隔。若只错开很小一点，第二组的某个本底事例会
            // 正好压在暴发窗边界上、凑成一次真的符合——那是数据的巧合，不是
            // 判据的毛病，但会让这个用例测不到想测的东西。
            events.extend((0..200).map(|i| TestEvent {
                seconds: i as f64 * 0.1 + group as f64 * 0.05,
                group,
            }));
        }
        for &group in groups {
            events.extend((0..30).map(|i| TestEvent {
                seconds: 10.0 + i as f64 * 3e-6,
                group,
            }));
        }
        events.sort_by(|a, b| a.seconds.partial_cmp(&b.seconds).unwrap());
        events
    }

    fn run(data: &[TestEvent]) -> Vec<crate::types::candidate::Candidate<TestInstrument>> {
        search_new(
            data,
            1,
            MissionElapsedTime::<TestInstrument>::new(0.0),
            MissionElapsedTime::<TestInstrument>::new(3600.0),
            &whole_span(),
            SearchConfig::default(),
        )
    }

    fn met(seconds: f64) -> MissionElapsedTime<TestInstrument> {
        MissionElapsedTime::new(seconds)
    }

    /// 整个 [0, 3600] 当活时间——即原先按 chunk 边界夹取的行为。
    fn whole_span() -> Vec<[MissionElapsedTime<TestInstrument>; 2]> {
        vec![[met(0.0), met(3600.0)]]
    }

    /// 每 0.5 µs 一个事例、只铺在 [10.000, 10.003] s 上：前面是空的，然后是
    /// 一块平台，正是 SAA 停机前那种形态（缩小版）。
    fn plateau_after_a_gap() -> Vec<TestEvent> {
        at(&(0..6000).map(|i| 10.0 + i as f64 * 5e-7).collect::<Vec<_>>())
    }

    fn run_with_live(
        data: &[TestEvent],
        gti: &[[f64; 2]],
    ) -> Vec<crate::types::candidate::Candidate<TestInstrument>> {
        let gti: Vec<[MissionElapsedTime<TestInstrument>; 2]> =
            gti.iter().map(|s| [met(s[0]), met(s[1])]).collect();
        search_new(
            data,
            1,
            met(0.0),
            met(3600.0),
            &gti,
            SearchConfig {
                min_duration: Time::new::<uom::si::time::microsecond>(0.0),
                max_duration: Time::new::<uom::si::time::microsecond>(100.0),
                neighbor: Time::new::<uom::si::time::millisecond>(4.0),
                hollow: Time::new::<uom::si::time::microsecond>(100.0),
                false_positive_per_year: 20.0,
                min_number: 8,
                coincidence: 1,
            },
        )
    }

    #[test]
    fn a_flat_plateau_next_to_a_gap_only_triggers_when_the_gap_is_counted_as_live() {
        let data = plateau_after_a_gap();
        // 把空区也当活时间（原先的行为）：平台头 2 ms 内的候选本底窗一半在空区，
        // 均值减半，200 个计数对期望 100，平台自己就触发。
        assert!(!run_with_live(&data, &[[0.0, 3600.0]]).is_empty());
        // 按活时间归一：期望就是 200，平台只是平台。
        assert!(run_with_live(&data, &[[10.0, 10.003]]).is_empty());
    }

    #[test]
    fn a_burst_right_after_a_gap_is_still_found() {
        // 平台上在恢复后 0.5 ms 处再塞 60 个事例进 30 µs：期望 60，实测 120。
        let mut data = plateau_after_a_gap();
        data.extend(at(&(0..60).map(|i| 10.0005 + i as f64 * 5e-7).collect::<Vec<_>>()));
        data.sort_by(|a, b| a.seconds.partial_cmp(&b.seconds).unwrap());
        let found = run_with_live(&data, &[[10.0, 10.003]]);
        assert_eq!(found.len(), 1);
        assert!((found[0].start.met() - 10.0005).abs() < 1e-4);
    }

    #[test]
    fn data_shorter_than_the_background_window_does_not_panic() {
        // 整段数据（49 ms）比半个 neighbor 窗（500 ms）还短。快照初始化会一路
        // 走到数组末尾 —— 早先的写法在这里索引 data[len] 越界 panic。
        let data = at(&(0..50).map(|i| i as f64 * 1e-3).collect::<Vec<_>>());
        // 等间隔 1 ms、max_duration 也是 1 ms，一个候选最多凑到 1 个事例，
        // 远不够 min_number=8，所以正确结果是没有候选。
        assert!(run(&data).is_empty());
    }

    #[test]
    fn a_single_event_does_not_panic() {
        assert!(run(&at(&[0.0])).is_empty());
    }

    #[test]
    fn events_before_start_are_skipped() {
        // start=0，负时刻的事例不该被当作起点
        let data = at(&[-2.0, -1.0, 10.0, 20.0]);
        assert!(run(&data).is_empty());
    }

    #[test]
    fn a_dense_burst_on_a_quiet_background_is_found() {
        let found = run(&burst_on_quiet_background(0));
        assert!(!found.is_empty(), "本底之上 30 个事例挤在 100us 里必须被找到");
    }

    #[test]
    fn a_single_group_degenerates_to_the_ungrouped_formula() {
        // 分组判据必须在 n=1 时精确退化：候选带的显著性要跟拿合并计数重算
        // 出来的一模一样。HXMT 与 SVOM 都跑在 n=1 上，这条守住它们逐位不变。
        let found = run_grouped(&burst_on_quiet_background(0), 1);
        assert!(!found.is_empty());
        for candidate in &found {
            assert_eq!(
                candidate.sf(),
                crate::poisson::sf(candidate.mean, candidate.count),
                "n=1 时存下来的 sf 必须与重算值逐位相同"
            );
        }
    }

    #[test]
    fn splitting_into_groups_costs_a_trials_factor() {
        // 同一份数据，事例全在第 0 组。分成两组等于每个窗口多了一次中奖机会，
        // Bonferroni 把显著性乘上组数，正好抵消翻倍的误报率。
        let one = run_grouped(&burst_on_quiet_background(0), 1);
        let two = run_grouped(&burst_on_quiet_background(0), 2);
        assert_eq!(one.len(), two.len());
        for (single, pair) in one.iter().zip(two.iter()) {
            assert_eq!(pair.sf(), single.sf() * 2.0, "两组应付 2 倍试验次数惩罚");
        }
    }

    #[test]
    fn requiring_two_groups_rejects_a_burst_that_lights_only_one() {
        // 只点亮一个组 —— GBM 上单探头爆发就是这个样子。要求两组符合时，
        // 第二小的组 p 值接近 1，判据过不去。
        let data = burst_in_groups(&[0]);
        assert!(!run_coincident(&data, 2, 1).is_empty(), "单组符合下本该找得到");
        assert!(
            run_coincident(&data, 2, 2).is_empty(),
            "要求两组符合时，只点亮一组的暴发必须被挡掉"
        );
    }

    #[test]
    fn requiring_two_groups_keeps_a_burst_seen_by_both() {
        // 两组同时超出 —— 各向同性照射的真信号该有的样子。
        assert!(
            !run_coincident(&burst_in_groups(&[0, 1]), 2, 2).is_empty(),
            "两组同时超出的暴发必须留下"
        );
    }

    #[test]
    fn a_burst_confined_to_a_later_group_is_still_found() {
        // 信号整个落在第 1 组。判据只看第 0 组的话这里什么都搜不到 —— 对 GBM
        // 就是把 BGO 分出去之后 NaI 侧的暴发全丢了。
        let found = run_grouped(&burst_on_quiet_background(1), 2);
        assert!(!found.is_empty(), "落在非第 0 组的暴发必须照样被找到");
    }
}
