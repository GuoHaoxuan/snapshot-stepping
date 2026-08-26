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
