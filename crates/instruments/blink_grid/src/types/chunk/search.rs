use blink_algorithms::detector_share::{MAX_DETECTOR_FRACTION, max_detector_fraction};
use blink_algorithms::snapshot_stepping::{SearchConfig, search_new};
use blink_core::traits::Event as _;
use blink_core::types::{MissionElapsedTime, Signal};
use std::sync::atomic::Ordering;
use uom::si::f64::*;

use super::Chunk;
use crate::io::posatt::{attitude_trajectory, interpolate_sampled, position_trajectory};
use crate::types::Event;
use crate::types::instrument::{Grid, Satellite};

/// 判为"读出空洞"的门槛：本底窗里最长的空段在本地速率下应有的计数 r·L。
///
/// 天格的读出在高计数率下会成帧丢数：辐射带里（|lat| > 40°，5–75 kc/s）事例
/// 以 5–16 ms 的密集帧到达，帧间 3–12 ms 一个事例都没有，四个探测器同步。
/// 搜索窗落在帧内、本底窗横跨帧和空洞，均值被帧间空洞拉低，帧内的普通计数
/// 就成了 fa=1e-199 的"暴发"——GRID-03B 一天 9213 个候选全是这种。实测
/// 2–15 kc/s 的秒空洞占比中位 0.0–0.7%，安静时段一天搜不出任何候选，所以
/// 只在本底窗里出现统计上不可能的空段时否决：r·L > 9.2 即 P(0) < 1e-4。
/// 在 1 kc/s 要 9 ms 的空段才触发（一个本底窗里约百次机会，误否决期望 1%），
/// 在 13 kc/s 只要 0.7 ms。全量首跑用 14（P<1e-6）时，辐射带 10–17 kc/s 的
/// 帧边界（空段接近 1 ms）漏了一批，幸存的显著候选一半仍在 |lat| ≥ 40°。
///
/// 速率 r 取本底窗速率与整次过境速率中的大者：成片丢数时窗内速率本身被压到
/// 几个 c/s（GRID-04 有 count=8、mean=0.00、率 2 c/s 的"候选"），按窗内速率
/// 算 r·L 永远到不了门槛，按过境速率算半秒空段就是几百。
///
/// 更正的做法是把空洞当 GTI 缺口、让本底按真实活时间归一，但 `search_new`
/// 的活时间夹取假设候选所在的一段 GTI 覆盖整个本底窗，毫秒级的洞会把本底窗
/// 夹到一帧里去，统计反而更差。见 `OPEN-QUESTIONS.md`。
const DEAD_GAP_EXPECTED_COUNTS: f64 = 9.2;

/// 读出可信的计数率上限（四路合计，c/s）。
///
/// 成帧丢数从这个量级开始：GRID-02 与 04 在 5 kc/s 附近就成片出现空洞，03B
/// 到 15 kc/s 才成片，但 13 kc/s 时已经有另一种形态——四路同时一个 1 ms 的
/// 尖峰（86 个对本底 4 个/ms）后面跟几毫秒的稀疏，合并四路后最长空段只有
/// 0.6 ms，空洞判据抓不到。全量首跑显著候选的一半落在 |lat| ≥ 40° 的高速率
/// 区，本底率中位 5.4 kc/s；加这道门后 |lat| < 40° 的显著候选一个不少
/// （03B 277/277），失去的只有辐射带。取本底窗速率与过境速率中的大者比较。
const RATE_CEILING: f64 = 5000.0;

/// 最显著一格里落在同一个时间戳上的事例占比上限。
///
/// 全量首跑速率门以下的 348 个显著候选（fa ≤ 1e-5）里，261 个有一半以上的计数
/// 挤在同一个时间戳上：四路各一个、时间戳完全相同（跨探测器 dt = 0），最显著
/// 一格只有这四个加零星本底（计数中位 4），谱比本底硬得多（PI 中位比 3.8 对
/// 另一群的 1.7），经度均匀分布。这是穿过整星的带电粒子——四块 GAGG 同时响应、
/// 各留一个大沉积——不是 TGF。另一群 87 个的计数在几百微秒里铺开（时长中位
/// 564 µs），最小间隔在 4.77 µs 的时间戳分辨率量级，经度集中在 120°E–180°
/// （海洋大陆雷暴区占 43%）。两群在占比 0.45–0.50 之间一个候选都没有。门槛与
/// GBM 同取 0.35：一个时间戳上超过三分之一的计数就当粒子否决。四路合成一路
/// 搜索没有组间符合可用，这道门就是天格挡带电粒子的唯一手段。见
/// `OPEN-QUESTIONS.md` 第 10 条。
const MAX_SIMULTANEOUS_FRACTION: f64 = 0.35;

/// 本底窗 `[from, to]`（已夹到候选所在的 GTI 段内）里是否有读出空洞。
///
/// `events` 已按时间排好。空段包括窗口两端到最近事例的距离——窗口已夹在
/// GTI 内，端点上没有事例不是过境边界造成的。
fn has_dead_gap<S: Satellite>(events: &[Event<S>], from: f64, to: f64, pass_rate: f64) -> bool {
    if to <= from {
        return false;
    }
    let lo = events.partition_point(|e| e.time().met() < from);
    let hi = events.partition_point(|e| e.time().met() <= to);
    let window = &events[lo..hi];
    let longest = if window.is_empty() {
        to - from
    } else {
        let mut longest =
            (window[0].time().met() - from).max(to - window[window.len() - 1].time().met());
        for pair in window.windows(2) {
            longest = longest.max(pair[1].time().met() - pair[0].time().met());
        }
        longest
    };
    let rate = (window.len() as f64 / (to - from)).max(pass_rate);
    rate * longest > DEAD_GAP_EXPECTED_COUNTS
}

/// 最显著一格 `[start, stop]` 里同一时间戳上最多有几个事例，占该格事例数的比例。
///
/// `events` 已按时间排好，时间戳相同的事例必然相邻，扫一遍取最长的一段。两端
/// 都是事例本身的时刻（`Candidate` 的 start 与 delay、bin_size_best 都来自事例），
/// 闭区间比较，边界上那一簇不会漏。
fn simultaneous_fraction<S: Satellite>(
    events: &[Event<S>],
    start: MissionElapsedTime<Grid<S>>,
    stop: MissionElapsedTime<Grid<S>>,
) -> f64 {
    let lo = events.partition_point(|e| e.time() < start);
    let hi = events.partition_point(|e| e.time() <= stop);
    let window = &events[lo..hi];
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
    let fitted = crate::io::orbit_fit::trajectory::<S>(&chunk.orbit_fit);

    let mut n_dropped = 0usize;
    let mut n_fitted = 0usize;
    let mut n_dead_gap = 0usize;
    let mut n_high_rate = 0usize;
    let mut n_simultaneous = 0usize;
    let mut n_no_attitude = 0usize;
    let mut n_single_detector = 0usize;
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
            // 过境速率：准入后的事例数 / 过境时长
            let pass_rate = chunk
                .passes
                .iter()
                .find(|p| cs >= p.start && cs <= p.stop)
                .map(|p| {
                    let a = events.partition_point(|e| e.time().met() < p.start);
                    let b = events.partition_point(|e| e.time().met() <= p.stop);
                    (b - a) as f64 / (p.stop - p.start).max(1e-9)
                })
                .unwrap_or(0.0);
            if has_dead_gap(&events, from, to, pass_rate) {
                n_dead_gap += 1;
                return None;
            }
            let window_rate = {
                let a = events.partition_point(|e| e.time().met() < from);
                let b = events.partition_point(|e| e.time().met() <= to);
                (b - a) as f64 / (to - from).max(1e-9)
            };
            if window_rate.max(pass_rate) > RATE_CEILING {
                n_high_rate += 1;
                return None;
            }
            // 带电粒子否决看最显著的那一格：由 start 偏移 delay 得到，两者都是
            // 事例的时刻，差与和在 f64 下精确。
            let best_start = candidate.start + candidate.delay;
            let best_stop = best_start + candidate.bin_size_best;
            if simultaneous_fraction(&events, best_start, best_stop) > MAX_SIMULTANEOUS_FRACTION {
                n_simultaneous += 1;
                return None;
            }
            // 单路毛刺否决。四块 GAGG 并排同向，真暴发四路均分：v3 全量真候选的单路
            // 最大占比中位 0.36、最高 0.56；超过 0.9 的 3 个全是一路探测器自己在闹
            // （03B 2023-08-08T19:53:50 一路占 100%、PI 6；2023-08-10 三分钟内两个
            // 候选同一路占 90–95%）。单组搜索没有组间符合，这是挡单路毛刺的唯一手段。
            if max_detector_fraction(&events, candidate.start, candidate.stop, |e| e.detector)
                > MAX_DETECTOR_FRACTION
            {
                n_single_detector += 1;
                return None;
            }
            let peak = candidate.start + candidate.bin_size_best / 2.0;
            // 2024-02 之后位姿文件里没有位置解（POS_TYPE=0，全 NaN），这样的
            // 候选定不了位。丢可以，静默丢不行——记账见 `diagnostics`。
            // 位姿有位置解就用它；没有（2024-02-09 起）退回拟合轨道表；都没有才丢
            let position = match interpolate_sampled(&positions, peak) {
                Some(p) => p,
                None => match interpolate_sampled(&fitted, peak) {
                    Some(p) => {
                        n_fitted += 1;
                        p
                    }
                    None => {
                        n_dropped += 1;
                        return None;
                    }
                },
            };
            // 姿态解也会整段缺失（v3 全量 82 个候选、3 个显著落在里面）。候选的
            // 实质是时间加位置，姿态只是元数据：缺了照留，留空并记账。
            let attitude = interpolate_sampled(&attitudes, peak);
            if attitude.is_none() {
                n_no_attitude += 1;
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
                position,
                // 天格没有反符合探测器
                acd: None,
            })
        })
        .collect::<Vec<_>>();

    chunk
        .dropped_no_ephemeris
        .store(n_dropped, Ordering::Relaxed);
    chunk
        .positions_from_orbit_fit
        .store(n_fitted, Ordering::Relaxed);
    chunk.events_outside_gti.store(n_outside, Ordering::Relaxed);
    chunk.dropped_dead_gap.store(n_dead_gap, Ordering::Relaxed);
    chunk
        .dropped_high_rate
        .store(n_high_rate, Ordering::Relaxed);
    chunk
        .dropped_simultaneous
        .store(n_simultaneous, Ordering::Relaxed);
    chunk
        .without_attitude
        .store(n_no_attitude, Ordering::Relaxed);
    chunk
        .dropped_single_detector
        .store(n_single_detector, Ordering::Relaxed);
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

    fn on_detectors(times_and_detectors: &[(f64, u8)]) -> Vec<Event<Sat03B>> {
        times_and_detectors
            .iter()
            .map(|(t, d)| Event {
                time: MissionElapsedTime::new(*t),
                channel: 20,
                detector: *d,
                evt_type: 1,
                energy_kev: 100.0,
                overflow: false,
            })
            .collect()
    }

    fn detector_share(events: &[Event<Sat03B>]) -> f64 {
        max_detector_fraction(
            events,
            events[0].time(),
            events[events.len() - 1].time(),
            |e| e.detector,
        )
    }

    #[test]
    fn a_burst_shared_by_the_four_detectors_is_kept() {
        let events = on_detectors(
            &(0..12)
                .map(|i| (100.0 + i as f64 * 1e-5, (i % 4) as u8))
                .collect::<Vec<_>>(),
        );
        assert!((detector_share(&events) - 0.25).abs() < 1e-12);
        assert!(detector_share(&events) <= MAX_DETECTOR_FRACTION);
    }

    #[test]
    fn a_glitch_in_one_detector_is_vetoed_even_with_a_background_count() {
        // 8 个计数里 7 个来自 2 号探测器：7/8 = 0.875
        let mut v: Vec<(f64, u8)> = (0..7).map(|i| (100.0 + i as f64 * 1e-5, 2)).collect();
        v.push((100.00005, 0));
        v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        assert!(detector_share(&on_detectors(&v)) > MAX_DETECTOR_FRACTION);
    }

    fn fraction(times: &[f64]) -> f64 {
        let events = at(times);
        simultaneous_fraction(&events, events[0].time(), events[events.len() - 1].time())
    }

    #[test]
    fn a_particle_crossing_puts_half_the_counts_on_one_timestamp() {
        // 四路各一个、时间戳完全相同，再加 4 个零星本底：4/8 = 0.5
        let times = [
            100.0, 100.0, 100.0, 100.0, 100.0002, 100.0004, 100.0006, 100.0008,
        ];
        assert!(fraction(&times) > MAX_SIMULTANEOUS_FRACTION);
    }

    #[test]
    fn a_burst_spread_over_the_timestamp_grid_is_kept() {
        // 20 个事例铺在 4.77 µs 的时间戳格上，每格最多两个：2/20 = 0.1
        let tick = 4.768e-6;
        let times: Vec<f64> = (0..20).map(|i| 100.0 + (i / 2) as f64 * tick).collect();
        assert!(fraction(&times) < MAX_SIMULTANEOUS_FRACTION);
    }

    #[test]
    fn the_cluster_on_the_window_edge_is_counted() {
        // 簇正好落在窗口末端：闭区间比较才数得到它，4/7
        let times = [
            100.0, 100.0001, 100.0002, 100.0003, 100.0003, 100.0003, 100.0003,
        ];
        assert!((fraction(&times) - 4.0 / 7.0).abs() < 1e-12);
    }

    #[test]
    fn a_uniform_stream_has_no_dead_gap() {
        // 2000 c/s 均匀铺满 1 s：最长空段 0.5 ms，r·L = 1
        let events = at(&(0..2000)
            .map(|i| 100.0 + i as f64 * 5e-4)
            .collect::<Vec<_>>());
        assert!(!has_dead_gap(&events, 100.0, 101.0, 2000.0));
    }

    #[test]
    fn a_frame_gap_at_high_rate_is_a_dead_gap() {
        // 帧结构：每 10 ms 前 5 ms 有 75 kc/s 的事例，后 5 ms 全空——r·L ≈ 37500×0.005 = 187
        let mut times = Vec::new();
        for frame in 0..100 {
            let t0 = 100.0 + frame as f64 * 0.01;
            times.extend((0..375).map(|i| t0 + i as f64 * 5e-3 / 375.0));
        }
        assert!(has_dead_gap(&at(&times), 100.0, 101.0, 37500.0));
    }

    #[test]
    fn a_long_but_poisson_plausible_gap_at_low_rate_is_not_a_dead_gap() {
        // 300 c/s，最长空段 20 ms：r·L = 6，P(0) = 2.5e-3，不能算读出空洞
        let mut times: Vec<f64> = (0..300).map(|i| 100.0 + i as f64 / 300.0).collect();
        times.retain(|t| !(100.50..100.52).contains(t));
        assert!(!has_dead_gap(&at(&times), 100.0, 101.0, 300.0));
    }

    #[test]
    fn the_window_edges_count_as_gaps() {
        // 窗口前 200 ms 一个事例都没有，之后 5 kc/s：r·L = 4000×0.2 = 800
        let events = at(&(0..4000)
            .map(|i| 100.2 + i as f64 * 2e-4)
            .collect::<Vec<_>>());
        assert!(has_dead_gap(&events, 100.0, 101.0, 2000.0));
    }

    #[test]
    fn a_nearly_empty_window_is_judged_by_the_pass_rate() {
        // 窗里只有 8 个事例挤在 1 ms 内（窗内速率 8 c/s），但整次过境是 1 kc/s：
        // 半秒的空段按过境速率是 500 个期望计数——这是成片丢数里的孤帧，不是暴发
        let events = at(&(0..8).map(|i| 100.5 + i as f64 * 1e-4).collect::<Vec<_>>());
        assert!(has_dead_gap(&events, 100.0, 101.0, 1000.0));
        // 同样的 8 个事例，如果整次过境本来就只有 8 c/s，那半秒空段是正常的
        assert!(!has_dead_gap(&events, 100.0, 101.0, 8.0));
    }

    #[test]
    fn a_millisecond_frame_edge_at_belt_rates_is_caught() {
        // 13 kc/s 的帧边界留 0.9 ms 空段：r·L = 11.7，P(0) = 8e-6——用 1e-6 的门槛会漏
        let mut times: Vec<f64> = (0..13000).map(|i| 100.0 + i as f64 / 13000.0).collect();
        times.retain(|t| !(100.5000..100.5009).contains(t));
        assert!(has_dead_gap(&at(&times), 100.0, 101.0, 13000.0));
    }
}
