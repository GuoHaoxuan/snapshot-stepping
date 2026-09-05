use blink_algorithms::snapshot_stepping::{SearchConfig, search_new};
use blink_core::traits::Event as _;
use blink_core::types::{MissionElapsedTime, Signal};
use std::sync::atomic::Ordering;
use uom::si::f64::*;

use super::Chunk;
use crate::io::posatt::{attitude_trajectory, position_trajectory};
use crate::types::Event;
use crate::types::instrument::{Grid, Satellite};

/// 判为"读出空洞"的门槛：本底窗里最长的空段在本地速率下应有的计数 r·L。
///
/// 天格的读出在高计数率下会成帧丢数：辐射带里（|lat| > 40°，5–75 kc/s）事例
/// 以 5–16 ms 的密集帧到达，帧间 3–12 ms 一个事例都没有，四个探测器同步。
/// 搜索窗落在帧内、本底窗横跨帧和空洞，均值被帧间空洞拉低，帧内的普通计数
/// 就成了 fa=1e-199 的"暴发"——GRID-03B 一天 9213 个候选全是这种。实测
/// 2–15 kc/s 的秒空洞占比中位 0.0–0.7%，安静时段一天搜不出任何候选，所以
/// 只在本底窗里出现统计上不可能的空段时否决：r·L > 14 即 P(0) < 1e-6。
/// 在 300 c/s 要 46 ms 的空段才触发，在 75 kc/s 只要 0.2 ms。
///
/// 更正的做法是把空洞当 GTI 缺口、让本底按真实活时间归一，但 `search_new`
/// 的活时间夹取假设候选所在的一段 GTI 覆盖整个本底窗，毫秒级的洞会把本底窗
/// 夹到一帧里去，统计反而更差。见 `OPEN-QUESTIONS.md`。
const DEAD_GAP_EXPECTED_COUNTS: f64 = 14.0;

/// 本底窗 `[from, to]`（已夹到候选所在的 GTI 段内）里是否有读出空洞。
///
/// `events` 已按时间排好。空段包括窗口两端到最近事例的距离——窗口已夹在
/// GTI 内，端点上没有事例不是过境边界造成的。
fn has_dead_gap<S: Satellite>(events: &[Event<S>], from: f64, to: f64) -> bool {
    let lo = events.partition_point(|e| e.time().met() < from);
    let hi = events.partition_point(|e| e.time().met() <= to);
    let n = hi - lo;
    if n < 2 || to <= from {
        return false;
    }
    let window = &events[lo..hi];
    let mut longest =
        (window[0].time().met() - from).max(to - window[window.len() - 1].time().met());
    for pair in window.windows(2) {
        longest = longest.max(pair[1].time().met() - pair[0].time().met());
    }
    let rate = n as f64 / (to - from);
    rate * longest > DEAD_GAP_EXPECTED_COUNTS
}

pub(super) fn search<S: Satellite>(chunk: &Chunk<S>) -> Vec<Signal<Event<S>>> {
    let gti: Vec<[MissionElapsedTime<Grid<S>>; 2]> = chunk
        .gti
        .iter()
        .map(|g| [MissionElapsedTime::new(g[0]), MissionElapsedTime::new(g[1])])
        .collect();
    let inside = |t: f64| gti.iter().any(|g| t >= g[0].met() && t <= g[1].met());

    // 四路探测器、若干次过境拼在一起再排序。事例准入见 `Event::keep`；实测事例
    // 都落在各自文件的 GTI 内，这里再按 GTI 过滤一次是为了和曝光核算同一口径，
    // 丢掉的数目记进诊断。
    let mut n_outside = 0usize;
    let mut events: Vec<Event<S>> = chunk
        .passes
        .iter()
        .flat_map(|p| p.events::<S>())
        .filter(|e| e.keep())
        .filter(|e| {
            let ok = inside(e.time().met());
            if !ok {
                n_outside += 1;
            }
            ok
        })
        .collect();
    events.sort();

    let results = search_new(
        &events,
        1,
        chunk.span[0],
        chunk.span[1],
        // 活时间就是这几段过境。天格一小时里往往只有十几分钟有数据，两头是
        // 硬边界，本底窗必须按活时间归一，否则边界候选的本底会被压低。
        &gti,
        SearchConfig {
            min_duration: Time::new::<uom::si::time::microsecond>(0.0),
            max_duration: Time::new::<uom::si::time::millisecond>(1.0),
            neighbor: Time::new::<uom::si::time::second>(1.0),
            hollow: Time::new::<uom::si::time::millisecond>(10.0),
            false_positive_per_year: 20.0,
            min_number: 8,
            // 单组：4 个 GAGG 同型，合成一路
            coincidence: 1,
        },
    );

    let attitudes = attitude_trajectory::<S>(&chunk.posatt);
    let positions = position_trajectory::<S>(&chunk.posatt);

    let mut n_dropped = 0usize;
    let mut n_dead_gap = 0usize;
    let half_neighbor = 0.5_f64;
    let signals = results
        .into_iter()
        .filter_map(|candidate| {
            // 读出空洞否决。本底窗夹到候选所在的那段 GTI 里，过境边界不算空洞。
            let (cs, ce) = (candidate.start.met(), candidate.stop.met());
            let seg = chunk.gti.iter().find(|g| cs >= g[0] && cs <= g[1]);
            let (from, to) = match seg {
                Some(g) => (
                    (cs - half_neighbor).max(g[0]),
                    (ce + half_neighbor).min(g[1]),
                ),
                None => (cs - half_neighbor, ce + half_neighbor),
            };
            if has_dead_gap(&events, from, to) {
                n_dead_gap += 1;
                return None;
            }
            let peak = candidate.start + candidate.bin_size_best / 2.0;
            // 2024-02 之后位姿文件里没有位置解（POS_TYPE=0，全 NaN），这样的
            // 候选定不了位。丢可以，静默丢不行——记账见 `diagnostics`。
            let (Some(attitude), Some(position)) =
                (attitudes.interpolate(peak), positions.interpolate(peak))
            else {
                n_dropped += 1;
                return None;
            };
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
                attitude: attitude.state,
                position: position.state,
                // 天格没有反符合探测器
                acd: None,
            })
        })
        .collect::<Vec<_>>();

    chunk
        .dropped_no_ephemeris
        .store(n_dropped, Ordering::Relaxed);
    chunk.events_outside_gti.store(n_outside, Ordering::Relaxed);
    chunk.dropped_dead_gap.store(n_dead_gap, Ordering::Relaxed);
    signals
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Sat03B;

    fn at(times: &[f64]) -> Vec<Event<Sat03B>> {
        times
            .iter()
            .map(|t| Event {
                time: MissionElapsedTime::new(*t),
                channel: 20,
                detector: 0,
                evt_type: 1,
                energy_kev: 100.0,
                overflow: false,
            })
            .collect()
    }

    #[test]
    fn a_uniform_stream_has_no_dead_gap() {
        // 2000 c/s 均匀铺满 1 s：最长空段 0.5 ms，r·L = 1
        let events = at(&(0..2000)
            .map(|i| 100.0 + i as f64 * 5e-4)
            .collect::<Vec<_>>());
        assert!(!has_dead_gap(&events, 100.0, 101.0));
    }

    #[test]
    fn a_frame_gap_at_high_rate_is_a_dead_gap() {
        // 帧结构：每 10 ms 前 5 ms 有 75 kc/s 的事例，后 5 ms 全空——r·L ≈ 37500×0.005 = 187
        let mut times = Vec::new();
        for frame in 0..100 {
            let t0 = 100.0 + frame as f64 * 0.01;
            times.extend((0..375).map(|i| t0 + i as f64 * 5e-3 / 375.0));
        }
        assert!(has_dead_gap(&at(&times), 100.0, 101.0));
    }

    #[test]
    fn a_long_but_poisson_plausible_gap_at_low_rate_is_not_a_dead_gap() {
        // 300 c/s，最长空段 20 ms：r·L = 6，P(0) = 2.5e-3，不能算读出空洞
        let mut times: Vec<f64> = (0..300).map(|i| 100.0 + i as f64 / 300.0).collect();
        times.retain(|t| !(100.50..100.52).contains(t));
        assert!(!has_dead_gap(&at(&times), 100.0, 101.0));
    }

    #[test]
    fn the_window_edges_count_as_gaps() {
        // 窗口前 200 ms 一个事例都没有，之后 5 kc/s：r·L = 4000×0.2 = 800
        let events = at(&(0..4000)
            .map(|i| 100.2 + i as f64 * 2e-4)
            .collect::<Vec<_>>());
        assert!(has_dead_gap(&events, 100.0, 101.0));
    }
}
