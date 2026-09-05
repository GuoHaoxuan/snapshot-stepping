use blink_core::error::Error;
use blink_core::types::{Coverage, ExclusionReason, MissionElapsedTime};
use chrono::prelude::*;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::io::{PassFile, PosAttFile};
use crate::types::Event;
use crate::types::instrument::{Grid, Satellite};

mod from_epoch;
mod search;

/// 一小时。天格是逐过境记录的，一天只有 1–21 次过境、0.5–10 小时数据，
/// 所以多数小时根本没有文件（记为 `missing_data`，曝光为零），有文件的小时
/// 里也只有过境那几段是活时间——曝光按 GTI 算，本底也按 GTI 归一。
pub struct Chunk<S: Satellite> {
    pub span: [MissionElapsedTime<Grid<S>>; 2],
    /// 与本小时有交集的过境事例文件
    pub passes: Vec<PassFile>,
    /// 同一天（及跨午夜时前一天）的全部位姿文件
    pub posatt: Vec<PosAttFile>,
    /// 各过境的 GTI 截到本小时内，按时间排好
    pub gti: Vec<[f64; 2]>,
    pub(super) dropped_no_ephemeris: AtomicUsize,
    pub(super) events_outside_gti: AtomicUsize,
    /// 本底窗里有读出空洞而被否决的候选数，见 `search`
    pub(super) dropped_dead_gap: AtomicUsize,
    _satellite: PhantomData<S>,
}

impl<S: Satellite> blink_core::traits::Chunk for Chunk<S> {
    type Event = Event<S>;

    fn from_epoch(epoch: &DateTime<Utc>) -> Result<Self, Error> {
        from_epoch::from_epoch::<S>(epoch)
    }

    fn search(&self) -> Vec<blink_core::types::Signal<Self::Event>> {
        search::search(self)
    }

    fn exclusion(&self) -> Option<ExclusionReason> {
        if self.passes.iter().all(PassFile::is_empty) {
            return Some(ExclusionReason::NoEvents);
        }
        // 与 GRM/GBM 同一道门：探测器内时间回跳会让搜索的窗长判据失效
        if self.passes.iter().any(|p| p.time_reversals() > 0) {
            return Some(ExclusionReason::UnorderedEvents);
        }
        None
    }

    fn coverage(&self) -> Coverage {
        let span_seconds = self.span[1].met() - self.span[0].met();
        let live: f64 = self.gti.iter().map(|g| g[1] - g[0]).sum();
        Coverage {
            span_seconds,
            masked_seconds: (span_seconds - live).max(0.0),
        }
    }

    fn diagnostics(&self) -> Vec<(&'static str, f64)> {
        let mut d = vec![
            ("n_passes", self.passes.len() as f64),
            (
                "n_events",
                self.passes.iter().map(PassFile::len).sum::<usize>() as f64,
            ),
            (
                "time_reversals",
                self.passes
                    .iter()
                    .map(PassFile::time_reversals)
                    .sum::<usize>() as f64,
            ),
            (
                "posatt_rows",
                self.posatt.iter().map(PosAttFile::len).sum::<usize>() as f64,
            ),
            (
                "posatt_positioned_rows",
                self.posatt
                    .iter()
                    .map(PosAttFile::positioned_rows)
                    .sum::<usize>() as f64,
            ),
        ];
        let dropped = self.dropped_no_ephemeris.load(Ordering::Relaxed);
        if dropped > 0 {
            d.push(("dropped_no_ephemeris", dropped as f64));
        }
        let outside = self.events_outside_gti.load(Ordering::Relaxed);
        if outside > 0 {
            d.push(("events_outside_gti", outside as f64));
        }
        let dead = self.dropped_dead_gap.load(Ordering::Relaxed);
        if dead > 0 {
            d.push(("dropped_dead_gap", dead as f64));
        }
        d
    }

    fn last_modified(epoch: &DateTime<Utc>) -> Result<DateTime<Utc>, Error> {
        from_epoch::last_modified::<S>(epoch)
    }
}
