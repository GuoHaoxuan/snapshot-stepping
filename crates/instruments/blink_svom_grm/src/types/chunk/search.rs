use crate::types::Chunk;
use crate::types::Event;
use crate::types::SvomGrm;
use blink_algorithms::snapshot_stepping::SearchConfig;
use blink_algorithms::snapshot_stepping::search_new;
use blink_core::types::Attitude;
use blink_core::types::MissionElapsedTime;
use blink_core::types::Position;
use blink_core::types::Signal;
use blink_core::types::Trajectory;
use uom::si::f64::*;

pub(super) fn search(chunk: &Chunk) -> Vec<Signal<Event>> {
    let events = chunk.evt_file.into_iter().collect::<Vec<_>>();
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
        },
    );

    results
        .into_iter()
        .filter_map(|candidate| {
            let peak = candidate.start + candidate.bin_size_best / 2.0;
            let attitude =
                Trajectory::<MissionElapsedTime<SvomGrm>, Attitude>::from(&chunk.att_file)
                    .interpolate(peak)?;
            let position =
                Trajectory::<MissionElapsedTime<SvomGrm>, Position>::from(&chunk.orb_file)
                    .interpolate(peak)?;
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
                // GRM 无逐事例反符合信息
                acd: None,
            })
        })
        .collect::<Vec<_>>()
}
