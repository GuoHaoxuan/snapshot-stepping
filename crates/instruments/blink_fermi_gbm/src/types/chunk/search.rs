use blink_algorithms::snapshot_stepping::{SearchConfig, search_new};
use blink_core::traits::Event as _;
use blink_core::types::{Attitude, MissionElapsedTime, Position, Signal, Trajectory};
use std::sync::atomic::Ordering;
use uom::si::f64::*;

use super::Chunk;
use crate::types::{Event, FermiGbm};

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

    let results = search_new(
        &events,
        chunk.groups.len(),
        chunk.span[0],
        chunk.span[1],
        SearchConfig {
            min_duration: Time::new::<uom::si::time::microsecond>(0.0),
            max_duration: Time::new::<uom::si::time::millisecond>(1.0),
            neighbor: Time::new::<uom::si::time::second>(1.0),
            hollow: Time::new::<uom::si::time::millisecond>(10.0),
            false_positive_per_year: 20.0,
            min_number: 8,
        },
    );

    let attitudes = Trajectory::<MissionElapsedTime<FermiGbm>, Attitude>::from(&chunk.poshist);
    let positions = Trajectory::<MissionElapsedTime<FermiGbm>, Position>::from(&chunk.poshist);

    let mut n_dropped = 0usize;
    let signals = results
        .into_iter()
        .filter_map(|candidate| {
            let peak = candidate.start + candidate.bin_size_best / 2.0;
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
                // GBM 没有反符合探测器
                acd: None,
            })
        })
        .collect::<Vec<_>>();

    chunk
        .dropped_no_ephemeris
        .store(n_dropped, Ordering::Relaxed);

    signals
}
