//! 单路探测器占比：挡住一路探测器自己闹出来的"暴发"。

use std::collections::BTreeMap;

use blink_core::traits::Event;
use blink_core::types::MissionElapsedTime;

/// 最显著一格里单路探测器贡献的事例占比上限，四个仪器共用。
///
/// 真暴发照亮整台仪器，实测真候选的单路最大占比都在 0.6 以下：GRID 四路并排
/// 同向，v3 全量真候选中位 0.36、最高 0.56；SVOM 三个 GRD 朝向不同，占比 ≤ 0.6
/// 的候选闪电关联 42%，> 0.8 的 12 个关联 0 个；GBM 试跑日过同时性判据的候选
/// ≤ 0.6。单路毛刺则是 0.9–1.0（GRID 03B 2023-08-08 一路占 100%、PI 6；SVOM
/// 2026-01-25 一路占 100%、60 个计数）。取 0.8：8 个计数里 7 个来自一路也否决。
pub const MAX_DETECTOR_FRACTION: f64 = 0.8;

/// 窗口 `[start, stop]` 里贡献最多的那路探测器占该窗事例数的比例。
///
/// `events` 已按时间排好；`detector` 给出每个事例的探测器标识。两端都是事例本身
/// 的时刻（`Candidate` 的 start 偏移 delay 得到最显著一格），闭区间比较。
pub fn max_detector_fraction<E: Event, K: Ord>(
    events: &[E],
    start: MissionElapsedTime<E::Instrument>,
    stop: MissionElapsedTime<E::Instrument>,
    detector: impl Fn(&E) -> K,
) -> f64 {
    let lo = events.partition_point(|e| e.time() < start);
    let hi = events.partition_point(|e| e.time() <= stop);
    let window = &events[lo..hi];
    if window.is_empty() {
        return 0.0;
    }
    let mut counts: BTreeMap<K, usize> = BTreeMap::new();
    for e in window {
        *counts.entry(detector(e)).or_insert(0) += 1;
    }
    counts.values().copied().max().unwrap_or(0) as f64 / window.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::TestEvent;

    // 用 group 字段充当探测器编号
    fn events(seconds_and_detectors: &[(f64, u8)]) -> Vec<TestEvent> {
        seconds_and_detectors
            .iter()
            .map(|&(seconds, group)| TestEvent { seconds, group })
            .collect()
    }

    fn share(events: &[TestEvent]) -> f64 {
        max_detector_fraction(
            events,
            events[0].time(),
            events[events.len() - 1].time(),
            |e| e.group,
        )
    }

    #[test]
    fn a_burst_shared_by_the_detectors_is_below_the_limit() {
        let v: Vec<(f64, u8)> = (0..12)
            .map(|i| (100.0 + i as f64 * 1e-5, (i % 4) as u8))
            .collect();
        let e = events(&v);
        assert!((share(&e) - 0.25).abs() < 1e-12);
        assert!(share(&e) <= MAX_DETECTOR_FRACTION);
    }

    #[test]
    fn seven_of_eight_counts_on_one_detector_is_over_the_limit() {
        let mut v: Vec<(f64, u8)> = (0..7).map(|i| (100.0 + i as f64 * 1e-5, 2)).collect();
        v.push((100.00005, 0));
        v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let e = events(&v);
        assert!((share(&e) - 0.875).abs() < 1e-12);
        assert!(share(&e) > MAX_DETECTOR_FRACTION);
    }

    #[test]
    fn only_the_window_is_counted_and_its_edges_are_inclusive() {
        // 窗外的一路毛刺不影响窗内的判断；窗口两端的事例都算在内
        let v = [
            (99.0, 5),
            (99.1, 5),
            (100.0, 0),
            (100.0001, 1),
            (100.0002, 2),
            (100.0003, 3),
            (101.0, 5),
        ];
        let e = events(&v);
        let f = max_detector_fraction(&e, e[2].time(), e[5].time(), |x| x.group);
        assert!((f - 0.25).abs() < 1e-12);
        assert_eq!(
            max_detector_fraction(
                &e,
                MissionElapsedTime::new(100.5),
                MissionElapsedTime::new(100.6),
                |x| x.group
            ),
            0.0
        );
    }
}
