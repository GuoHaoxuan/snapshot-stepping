use super::Chunk;
use crate::algorithms::acd::acd_counts;
use crate::types::{Event, HxmtHe};
use blink_algorithms::snapshot_stepping::{SearchConfig, search_new};
use blink_core::traits::Event as _;
use blink_core::types::{Attitude, MissionElapsedTime, Position, Signal, Trajectory};
use std::sync::atomic::Ordering;
use uom::si::f64::*;

pub fn search(chunk: &Chunk) -> Vec<Signal<Event>> {
    let events = chunk
        .event_file
        .into_iter()
        .filter(|event| event.keep())
        .collect::<Vec<_>>();

    let results = search_new(
        &events,
        1,
        chunk.span[0],
        chunk.span[1],
        SearchConfig {
            min_duration: Time::new::<uom::si::time::microsecond>(0.0),
            max_duration: Time::new::<uom::si::time::millisecond>(1.0),
            neighbor: Time::new::<uom::si::time::second>(1.0),
            hollow: Time::new::<uom::si::time::millisecond>(10.0),
            false_positive_per_year: 20.0,
            min_number: 8,
            // 单组：18 个 CsI 合成一路，符合判据退化
            coincidence: 1,
        },
    );

    // 饱和排除：用 1B FIFO reset 检测得到的（已扩 ±1s 并求并的）饱和区间，
    // 剔除落在其中的候选。旧的 continuous() 成簇判据已移除——成簇本身不再
    // 作为饱和信号，留作可能的真实发现。
    let saturation_intervals = chunk.saturation_intervals();

    let results = results
        .into_iter()
        .filter(|candidate| {
            let idx = saturation_intervals.partition_point(|iv| iv.1 < candidate.start);
            // 不在任何饱和区间内才保留
            !(idx < saturation_intervals.len() && saturation_intervals[idx].0 <= candidate.start)
        })
        .collect::<Vec<_>>();

    // 两条轨道各建一次。以前写在下面的闭包里，每个候选都要把整小时的采样点
    // （0.25 s 一采，约 14400 个）重新物化一遍。
    let attitudes = Trajectory::<MissionElapsedTime<HxmtHe>, Attitude>::from(&chunk.att_file);
    let positions = Trajectory::<MissionElapsedTime<HxmtHe>, Position>::from(&chunk.orbit_file);

    let mut n_dropped = 0usize;
    let signals = results
        .into_iter()
        .filter_map(|candidate| {
            let peak = candidate.start + candidate.bin_size_best / 2.0;
            // 星历表两端会外推一小截（见 Trajectory::interpolate）。仍取不到
            // 说明表缺了一整段以上，此时候选只能丢 —— 但要计数，不能静默。
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
                // ACD 符合计数必须在此刻统计：候选表落盘后事例流就不在手边了
                acd: Some(acd_counts(
                    &events,
                    candidate.start.met(),
                    candidate.stop.met(),
                )),
            })
        })
        .collect::<Vec<_>>();

    chunk
        .dropped_no_ephemeris
        .store(n_dropped, Ordering::Relaxed);

    signals
}
