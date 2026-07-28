use crate::{constants::DAYS_PER_YEAR, types::candidate::Candidate};
use blink_core::{traits::Event, types::MissionElapsedTime};
use statrs::distribution::{DiscreteCDF, Poisson};
use uom::si::f64::*;

pub struct SearchConfig {
    pub min_duration: Time,
    pub max_duration: Time,
    pub neighbor: Time,
    pub hollow: Time,
    pub false_positive_per_year: f64,
    pub min_number: u32,
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
        }
    }
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

pub fn search_new<E: Event>(
    data: &[E],
    group_number: usize,
    start: MissionElapsedTime<E::Instrument>,
    stop: MissionElapsedTime<E::Instrument>,
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

    loop {
        let mut step = 0;
        let mut numbers: Vec<u32> = vec![0; group_number];
        numbers[data[cursor].group() as usize] = 1;
        let mut mean_stop = mean_stop_snapshot;
        let mut mean_numbers = mean_numbers_snapshot.clone();
        let mut hollow_stop = hollow_stop_snapshot;
        let mut hollow_numbers = hollow_numbers_snapshot.clone();

        loop {
            let total_number = numbers.iter().sum(); // [TODO] Use real total number calculation
            let duration = data[cursor + step].time() - data[cursor].time();
            if total_number >= config.min_number && duration >= config.min_duration {
                let mean_start_time = (data[cursor].time() - config.neighbor / 2.0).max(start);
                let mean_stop_time = (data[cursor + step].time() + config.neighbor / 2.0).min(stop);
                let hollow_start_time = (data[cursor].time() - config.hollow / 2.0).max(start);
                let hollow_stop_time = (data[cursor + step].time() + config.hollow / 2.0).min(stop);
                let pure_mean_duration =
                    (mean_stop_time - mean_start_time) - (hollow_stop_time - hollow_start_time);
                let pure_mean_percent =
                    (duration / pure_mean_duration).get::<uom::si::ratio::ratio>();
                let fps = (0..group_number)
                    .map(|group| {
                        let pure_mean_number = mean_numbers[group] - hollow_numbers[group];
                        let equivalent_background_number =
                            pure_mean_number as f64 * pure_mean_percent;
                        match (equivalent_background_number, numbers[group]) {
                            (0.0, 0) => 1.0,
                            (0.0, _) => 1.0,
                            _ => Poisson::new(equivalent_background_number)
                                .unwrap()
                                .sf(numbers[group] as u64),
                        }
                    })
                    .collect::<Vec<f64>>();
                let fp = fps[0];
                let threshold = config.false_positive_per_year
                    / (uom::si::f64::Time::new::<uom::si::time::second>(3600.0)
                        * 24.0
                        * DAYS_PER_YEAR
                        / duration)
                        .get::<uom::si::ratio::ratio>();
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
            0
        }
        fn keep(&self) -> bool {
            true
        }
    }

    fn at(seconds: &[f64]) -> Vec<TestEvent> {
        seconds.iter().map(|&s| TestEvent { seconds: s }).collect()
    }

    fn run(data: &[TestEvent]) -> Vec<crate::types::candidate::Candidate<TestInstrument>> {
        search_new(
            data,
            1,
            MissionElapsedTime::<TestInstrument>::new(0.0),
            MissionElapsedTime::<TestInstrument>::new(3600.0),
            SearchConfig::default(),
        )
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
        // 本底 200 个事例铺在 20 s 上（10/s），再在 10.0 s 处塞 30 个事例进 100 us
        let mut seconds: Vec<f64> = (0..200).map(|i| i as f64 * 0.1).collect();
        seconds.extend((0..30).map(|i| 10.0 + i as f64 * 3e-6));
        seconds.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let found = run(&at(&seconds));
        assert!(!found.is_empty(), "本底之上 30 个事例挤在 100us 里必须被找到");
    }
}
