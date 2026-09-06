use blink_core::error::Error;
use blink_core::types::{Coverage, ExclusionReason, MissionElapsedTime};
use chrono::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::io::file::{find_poshist, find_tte};
use crate::io::{PosHistFile, TteFile};
use crate::types::{Detector, Event, FermiGbm};

mod from_epoch;
mod search;

pub struct Chunk {
    pub span: [MissionElapsedTime<FermiGbm>; 2],
    /// 这一小时找得到的所有探头，NaI 在前 BGO 在后。数量随年份变：2017-10
    /// 之前只有 BGO 的逐小时数据，2020 年主目录又只剩 NaI（BGO 从单独的
    /// BGO/ 目录取）。
    pub tte_files: Vec<TteFile>,
    pub poshist: PosHistFile,
    /// 本小时实际到齐的探测器类型，顺序即分组下标。
    pub groups: Vec<Detector>,
    pub(super) dropped_no_ephemeris: AtomicUsize,
    /// 峰值时刻取不到姿态、姿态留空的候选数（位置有，候选照留），见 `search`
    pub(super) without_attitude: AtomicUsize,
    /// 候选窗里单路探测器占比过高（单路毛刺）而被否决的候选数，见 `search`
    pub(super) dropped_single_detector: AtomicUsize,
    /// 因同一时间戳上挤了太多计数而被判为带电粒子、丢掉的候选数。
    pub(super) dropped_simultaneous: AtomicUsize,
    /// 本小时内、落在各探头 GTI 并集之外而搜索前丢掉的事例数，见 `search`
    /// （整点前 120 s 的文件引导段不算在内）
    pub(super) events_outside_gti: AtomicUsize,
}

/// 把一批区间合并成不相交的并集（相隔不到 1 s 的也接上），并截到 `[start, stop]`。
pub(super) fn union_intervals(mut rows: Vec<[f64; 2]>, start: f64, stop: f64) -> Vec<[f64; 2]> {
    rows.retain(|r| r[1] > r[0]);
    for r in rows.iter_mut() {
        r[0] = r[0].max(start);
        r[1] = r[1].min(stop);
    }
    rows.retain(|r| r[1] > r[0]);
    rows.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap());
    let mut out: Vec<[f64; 2]> = Vec::new();
    for r in rows {
        match out.last_mut() {
            Some(last) if r[0] <= last[1] + 1.0 => last[1] = last[1].max(r[1]),
            _ => out.push(r),
        }
    }
    out
}

impl Chunk {
    /// 这一小时的活时间：各探头 TTE 的 GTI 并集。
    ///
    /// 逐小时 TTE 在 SAA 期间降压停数，文件的 GTI 在进入 SAA 处截止；搜索若把
    /// 整小时当活时间，紧挨 SAA 的本底窗会伸进死区、把均值压低（SVOM 上同样
    /// 的机制造出过两类假信号）。合并流里只要有一路在记数就算活着，故取并集。
    pub fn gti_union(&self) -> Vec<[f64; 2]> {
        let rows = self
            .tte_files
            .iter()
            .flat_map(|file| file.gti_rows().map(|(a, b)| [a, b]))
            .collect();
        union_intervals(rows, self.span[0].met(), self.span[1].met())
    }
}

impl blink_core::traits::Chunk for Chunk {
    type Event = Event;

    fn from_epoch(epoch: &DateTime<Utc>) -> Result<Self, Error> {
        from_epoch::from_epoch(epoch)
    }

    fn search(&self) -> Vec<blink_core::types::Signal<Self::Event>> {
        search::search(self)
    }

    fn exclusion(&self) -> Option<ExclusionReason> {
        if self.tte_files.iter().all(TteFile::is_empty) {
            return Some(ExclusionReason::NoEvents);
        }
        // 与 GRM 同一道门：事例流乱序会让搜索的窗长判据失效、计数虚涨。
        if self.tte_files.iter().any(|file| file.time_reversals() > 0) {
            return Some(ExclusionReason::UnorderedEvents);
        }
        None
    }

    /// 曝光取各探头 GTI 的交集长度的下界：以事例最多的那一路为准。
    ///
    /// 14 个探头各有各的 GTI，逐一求交太细而收益不明；SAA 期间 GBM 会降压，
    /// TTE 本身就断，所以 GTI 已经把 SAA 排掉了（poshist 的 FLAGS 位可另做
    /// 交叉验证，见 `PosHistFile::saa_seconds_within`）。
    fn coverage(&self) -> Coverage {
        let (start, stop) = (self.span[0].met(), self.span[1].met());
        let gti_seconds = self
            .tte_files
            .iter()
            .map(|file| file.gti_seconds_within(start, stop))
            .fold(0.0_f64, f64::max);
        let span_seconds = stop - start;

        Coverage {
            span_seconds,
            masked_seconds: (span_seconds - gti_seconds).max(0.0),
        }
    }

    fn diagnostics(&self) -> Vec<(&'static str, f64)> {
        let (start, stop) = (self.span[0].met(), self.span[1].met());
        let mut diagnostics = vec![
            ("n_detectors", self.tte_files.len() as f64),
            ("n_groups", self.groups.len() as f64),
            (
                "n_events",
                self.tte_files.iter().map(TteFile::len).sum::<usize>() as f64,
            ),
            (
                "time_reversals",
                self.tte_files
                    .iter()
                    .map(TteFile::time_reversals)
                    .sum::<usize>() as f64,
            ),
            ("saa_seconds", self.poshist.saa_seconds_within(start, stop)),
        ];
        let dropped = self.dropped_no_ephemeris.load(Ordering::Relaxed);
        if dropped > 0 {
            diagnostics.push(("dropped_no_ephemeris", dropped as f64));
        }
        let without_attitude = self.without_attitude.load(Ordering::Relaxed);
        if without_attitude > 0 {
            diagnostics.push(("without_attitude", without_attitude as f64));
        }
        let single = self.dropped_single_detector.load(Ordering::Relaxed);
        if single > 0 {
            diagnostics.push(("dropped_single_detector", single as f64));
        }
        let simultaneous = self.dropped_simultaneous.load(Ordering::Relaxed);
        if simultaneous > 0 {
            diagnostics.push(("dropped_simultaneous", simultaneous as f64));
        }
        let outside = self.events_outside_gti.load(Ordering::Relaxed);
        if outside > 0 {
            diagnostics.push(("events_outside_gti", outside as f64));
        }
        diagnostics
    }

    fn last_modified(epoch: &DateTime<Utc>) -> Result<DateTime<Utc>, Error> {
        let mut paths: Vec<std::path::PathBuf> = Detector::NAI_NAMES
            .iter()
            .chain(Detector::BGO_NAMES.iter())
            .filter_map(|detector| find_tte(epoch, detector))
            .collect();
        if let Some(poshist) = find_poshist(epoch) {
            paths.push(poshist);
        }

        paths
            .iter()
            .map(|path| {
                let modified = std::fs::metadata(path)?.modified()?;
                Ok::<DateTime<Utc>, Error>(DateTime::<Utc>::from(modified))
            })
            .collect::<Result<Vec<_>, Error>>()?
            .into_iter()
            .max()
            .ok_or_else(|| Error::FileNotFound("no GBM files for this hour".to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::union_intervals;

    #[test]
    fn overlapping_and_touching_rows_merge_and_clip_to_the_span() {
        let rows = vec![[10.0, 100.0], [50.0, 120.0], [120.5, 200.0], [300.0, 400.0], [500.0, 450.0]];
        assert_eq!(
            union_intervals(rows, 20.0, 350.0),
            vec![[20.0, 200.0], [300.0, 350.0]]
        );
    }

    #[test]
    fn rows_apart_by_more_than_a_second_stay_separate() {
        let rows = vec![[0.0, 10.0], [11.5, 20.0]];
        assert_eq!(union_intervals(rows, 0.0, 100.0), vec![[0.0, 10.0], [11.5, 20.0]]);
    }
}
