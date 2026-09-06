use blink_algorithms::detector_share::{MAX_DETECTOR_FRACTION, max_detector_fraction};
use blink_algorithms::snapshot_stepping::{SearchConfig, search_new};
use blink_core::traits::Event as _;
use blink_core::types::{Attitude, MissionElapsedTime, Position, Signal, Trajectory};
use std::sync::atomic::Ordering;
use uom::si::f64::*;

use super::Chunk;
use crate::types::{Event, FermiGbm};

/// 允许挤在同一个时间戳上的计数占窗内计数的最大比例。
///
/// 带电粒子穿过整星时十几个探头在同一个时间戳上各留一个沉积，而 TGF 的光子
/// 铺在几十到几百微秒上，这个比例是两者之间最干净的一刀。实测 2019-01-01 一
/// 天里 fa≤1e-5 的 41 个候选分成不重叠的两群：时长 54–394 µs、事例散在
/// 14–89 个时间戳上的那 9 个比例为 0.03–0.32，而时长 6–32 µs、8–14 个计数
/// 散在 6–10 个探头上（每探头一个）的那 32 个是 0.45–0.83。全体 47904 个
/// 候选里 96.3% 带后一种形态。
///
/// 阈值取在两群中间，还没有用更多天的数据定过，见本 crate 的
/// `OPEN-QUESTIONS.md`。
const MAX_SIMULTANEOUS_FRACTION: f64 = 0.35;

/// 窗内挤在同一个时间戳上的最大计数占窗内总计数的比例。
///
/// `events` 已按时间排序，时间相同的事例必然相邻，扫一遍取最长的一段即可。
/// 窗口两端都是事例的时刻本身（`Candidate` 的 start/stop 取自 `data[cursor]`
/// 与 `data[cursor + step]`），所以用闭区间比较，边界上的那一簇不会漏掉。
fn simultaneous_fraction(
    events: &[Event],
    start: MissionElapsedTime<FermiGbm>,
    stop: MissionElapsedTime<FermiGbm>,
) -> f64 {
    let low = events.partition_point(|event| event.time() < start);
    let high = events.partition_point(|event| event.time() <= stop);
    let window = &events[low..high];
    if window.is_empty() {
        return 0.0;
    }

    let (mut longest, mut run) = (1usize, 1usize);
    for pair in window.windows(2) {
        run = if pair[1].time() == pair[0].time() {
            run + 1
        } else {
            1
        };
        longest = longest.max(run);
    }

    longest as f64 / window.len() as f64
}

pub(super) fn search(chunk: &Chunk) -> Vec<Signal<Event>> {
    // NaI 与 BGO 分组搜索：两者的响应差得太远，合成一路等于拿 BGO 的信号去
    // 配 NaI 的本底。`group_number` 取本小时实际到齐的类型数，这样只有 BGO
    // 的年份自动退回单组、不必付分组的试验次数惩罚。
    // 先按事例总数把容量要够：十几路加起来一小时有三千多万条，从零翻倍
    // 增长要重新分配二十几次、白搬约 1 GB。keep 之后只会更少，不会更多。
    let mut events: Vec<Event> =
        Vec::with_capacity(chunk.tte_files.iter().map(|file| file.len()).sum());
    for file in &chunk.tte_files {
        let group = chunk
            .groups
            .iter()
            .position(|detector| *detector == file.detector)
            .expect("每个 TTE 的类型都在 groups 里") as u8;
        events.extend(
            file.time()
                .iter()
                .zip(file.pha().iter())
                .map(|(time, channel)| Event {
                    time: MissionElapsedTime::new(*time),
                    channel: *channel,
                    detector: file.detector,
                    group,
                })
                .filter(Event::keep),
        );
    }
    // 十几路事例流拼起来必然是乱的，而 search_new 假定输入按时间有序。
    events.sort();

    // 活时间按各探头 GTI 的并集，事例也按它过滤：GTI 在 SAA 进入处截止，
    // 之后到下一文件之前是死区，本底窗不能伸进去（见 `Chunk::gti_union`）。
    let gti_seconds = chunk.gti_union();
    let inside = |t: f64| gti_seconds.iter().any(|g| t >= g[0] && t <= g[1]);
    let before = events.len();
    events.retain(|e| inside(e.time().met()));
    let n_outside = before - events.len();
    let gti: Vec<[MissionElapsedTime<FermiGbm>; 2]> = gti_seconds
        .iter()
        .map(|g| [MissionElapsedTime::new(g[0]), MissionElapsedTime::new(g[1])])
        .collect();

    let results = search_new(
        &events,
        chunk.groups.len(),
        chunk.span[0],
        chunk.span[1],
        &gti,
        SearchConfig {
            min_duration: Time::new::<uom::si::time::microsecond>(0.0),
            max_duration: Time::new::<uom::si::time::millisecond>(1.0),
            neighbor: Time::new::<uom::si::time::second>(1.0),
            hollow: Time::new::<uom::si::time::millisecond>(10.0),
            false_positive_per_year: 20.0,
            min_number: 8,
            // 三组里要求两组同时越线，用来挡穿过整星的宇宙线
            coincidence: 2,
        },
    );

    let attitudes = Trajectory::<MissionElapsedTime<FermiGbm>, Attitude>::from(&chunk.poshist);
    let positions = Trajectory::<MissionElapsedTime<FermiGbm>, Position>::from(&chunk.poshist);

    let mut n_dropped = 0usize;
    let mut n_without_attitude = 0usize;
    let mut n_single_detector = 0usize;
    let mut n_simultaneous = 0usize;
    let signals = results
        .into_iter()
        .filter_map(|candidate| {
            // 逐组符合挡不住穿过整星的带电粒子——它的定义特征恰恰就是多个
            // 探测器同时响应，组间符合反而在挑它。分开两者要看时间结构。
            // 判据用最显著的那一格，它由 start 偏移 delay 得到；两者都是事例
            // 的时刻，差与和都在 f64 下精确，复原不出偏差。
            let best_start = candidate.start + candidate.delay;
            let best_stop = best_start + candidate.bin_size_best;
            if simultaneous_fraction(&events, best_start, best_stop) > MAX_SIMULTANEOUS_FRACTION {
                n_simultaneous += 1;
                return None;
            }

            // 单路毛刺否决：候选窗里一路探测器占了绝大多数计数就不是暴发。
            // 组间符合已经挡住只在 NaI 里的毛刺，试跑日没有候选超过 0.8；这里与其余
            // 仪器口径一致。
            if max_detector_fraction(&events, candidate.start, candidate.stop, |e| e.detector)
                > MAX_DETECTOR_FRACTION
            {
                n_single_detector += 1;
                return None;
            }
            let peak = candidate.start + candidate.bin_size_best / 2.0;
            let Some(position) = positions.interpolate(peak) else {
                n_dropped += 1;
                return None;
            };
            // 姿态只是元数据，缺了不丢候选，留空并记账
            let attitude = attitudes.interpolate(peak).map(|a| a.state);
            if attitude.is_none() {
                n_without_attitude += 1;
            }
            Some(Signal {
                start: candidate.start,
                stop: candidate.stop,
                bin_size_min: candidate.bin_size_min,
                bin_size_max: candidate.bin_size_max,
                bin_size_best: candidate.bin_size_best,
                delay: candidate.delay,
                count: candidate.count,
                mean: candidate.mean,
                sf: candidate.sf(),
                false_positive_per_year: candidate.false_positive_per_year(),
                attitude,
                position: position.state,
                // GBM 没有反符合探测器
                acd: None,
            })
        })
        .collect::<Vec<_>>();

    chunk
        .dropped_no_ephemeris
        .store(n_dropped, Ordering::Relaxed);
    chunk
        .without_attitude
        .store(n_without_attitude, Ordering::Relaxed);
    chunk
        .dropped_single_detector
        .store(n_single_detector, Ordering::Relaxed);
    chunk
        .dropped_simultaneous
        .store(n_simultaneous, Ordering::Relaxed);
    chunk
        .events_outside_gti
        .store(n_outside, Ordering::Relaxed);

    signals
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Detector;

    fn events(times: &[f64]) -> Vec<Event> {
        times
            .iter()
            .map(|time| Event {
                time: MissionElapsedTime::new(*time),
                channel: 40,
                detector: Detector::Nai,
                group: 0,
            })
            .collect()
    }

    fn fraction(times: &[f64]) -> f64 {
        let events = events(times);
        let (start, stop) = (events[0].time(), events[events.len() - 1].time());
        simultaneous_fraction(&events, start, stop)
    }

    #[test]
    fn photons_spread_in_time_give_a_small_fraction() {
        // 十个各不相同的时刻，最长的一段就是一个事例。
        let times: Vec<f64> = (0..10).map(|i| 1000.0 + i as f64 * 1e-5).collect();
        assert!((fraction(&times) - 0.1).abs() < 1e-12);
    }

    #[test]
    fn a_particle_crossing_the_spacecraft_gives_a_large_fraction() {
        // 一个孤立事例加九个同时刻的沉积，正是宇宙线的形态。
        let mut times = vec![1000.0];
        times.extend(std::iter::repeat_n(1000.000_006, 9));
        assert!(fraction(&times) > MAX_SIMULTANEOUS_FRACTION);
        assert!((fraction(&times) - 0.9).abs() < 1e-12);
    }

    #[test]
    fn the_window_bounds_are_inclusive() {
        // 端点上的那一簇必须算进来：候选窗的两端本来就取自事例的时刻，
        // 用开区间会把整簇漏在窗外。
        let times = [1000.0, 1000.000_002, 1000.000_002, 1000.000_004];
        let events = events(&times);
        let low = MissionElapsedTime::new(1000.000_002);
        let high = MissionElapsedTime::new(1000.000_004);
        assert!((simultaneous_fraction(&events, low, high) - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn an_empty_window_is_not_rejected() {
        let events = events(&[1000.0, 1001.0]);
        let start = MissionElapsedTime::new(1000.5);
        let stop = MissionElapsedTime::new(1000.6);
        assert_eq!(simultaneous_fraction(&events, start, stop), 0.0);
    }
}
