use blink_core::types::{Coverage, MissionElapsedTime};

use crate::io::file::{find_att_by_time, find_evt_by_time, find_orb_by_time};
use crate::io::{AttFile, EvtFile, OrbFile};
use crate::types::event::Event;
use crate::types::instrument::SvomGrm;
use blink_core::error::Error;
use chrono::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

mod from_epoch;
mod search;

pub struct Chunk {
    pub span: [MissionElapsedTime<SvomGrm>; 2],
    pub att_file: AttFile,
    pub evt_file: EvtFile,
    pub orb_file: OrbFile,
    /// 本小时因取不到姿态/轨道而丢掉的候选数，由 `search` 写入。
    /// GRM 的 att/orb 是逐小时文件，尾端比事例流早收约 8 s（实测），
    /// 落在那一截里的候选只能丢——但要计数，不能静默。
    pub(super) dropped_no_ephemeris: AtomicUsize,
}

impl blink_core::traits::Chunk for Chunk {
    type Event = Event;

    fn from_epoch(epoch: &chrono::DateTime<chrono::Utc>) -> Result<Self, blink_core::error::Error>
    where
        Self: Sized,
    {
        from_epoch::from_epoch(epoch)
    }

    fn search(&self) -> Vec<blink_core::types::Signal<Self::Event>> {
        search::search(self)
    }

    /// 曝光核算走 L1B 的 GTI。GRM 的小时文件并不覆盖整点到整点：实测
    /// 2025-09-19 全天 GTI 合计只有 86400 s 的 85.3%，单个小时低到 1778 s。
    /// 缺口全部落在南大西洋异常区，即 GTI 已经把 SAA 排除掉了，所以这里
    /// 不需要再叠一层 SAA 掩模。
    fn coverage(&self) -> Coverage {
        let (start, stop) = (self.span[0].met(), self.span[1].met());
        // 相邻小时文件系统性重叠 100 s（TSTART = 整点 − 100 s），故 GTI 要
        // 先截到本小时内再累加，否则重叠段会被两个小时各记一次。
        let gti_seconds = self.evt_file.gti_seconds_within(start, stop);
        let span_seconds = stop - start;

        Coverage {
            span_seconds,
            masked_seconds: (span_seconds - gti_seconds).max(0.0),
        }
    }

    /// 在 `search` 之后取：`dropped_no_ephemeris` 是搜索过程中才产生的。
    fn diagnostics(&self) -> Vec<(&'static str, f64)> {
        let mut diagnostics = Vec::new();
        let dropped = self.dropped_no_ephemeris.load(Ordering::Relaxed);
        if dropped > 0 {
            diagnostics.push(("dropped_no_ephemeris", dropped as f64));
        }
        diagnostics
    }

    fn last_modified(epoch: &DateTime<Utc>) -> Result<DateTime<Utc>, Error> {
        let filenames = [
            find_att_by_time(epoch),
            find_evt_by_time(epoch),
            find_orb_by_time(epoch),
        ];

        let last_modifieds: Vec<DateTime<Utc>> = filenames
            .iter()
            .flatten()
            .map(|filename| {
                let last_modified = std::fs::metadata(filename)?.modified()?;
                let datetime: DateTime<Utc> = last_modified.into();
                Ok::<DateTime<Utc>, Error>(datetime)
            })
            .collect::<Result<Vec<DateTime<Utc>>, Error>>()?;

        let max_last_modified = last_modifieds
            .into_iter()
            .max()
            .ok_or_else(|| Error::FileNotFound("No files found".to_string()))?;

        Ok(max_last_modified)
    }
}
