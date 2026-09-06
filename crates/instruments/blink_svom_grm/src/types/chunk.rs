use blink_core::types::{Coverage, ExclusionReason, MissionElapsedTime};

use crate::io::file::{find_att_by_time, find_evt_by_time, find_orb_by_time};
use crate::io::{AttFile, EvtFile, OrbFile};
use crate::types::event::Event;
use crate::types::instrument::SvomGrm;
use blink_core::error::Error;
use chrono::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

mod from_epoch;
mod search;

/// 判为「整段重复」的时间回跳幅度下限（秒）。
///
/// 实测回跳幅度是双峰的，中间空着：一半在 1–5 个时间量化步（2⁻²⁰ s ≈ 0.95 µs）
/// 以内，另一半在 0.49–1.0 s。阈值取在空当里，并且正好等于搜索的
/// `max_duration` —— 回跳幅度超过搜索窗长才谈得上破坏窗长判据，小于它的
/// 抖动排序一下就没了。
pub const REVERSAL_THRESHOLD: f64 = 1e-3;

pub struct Chunk {
    pub span: [MissionElapsedTime<SvomGrm>; 2],
    pub att_file: AttFile,
    pub evt_file: EvtFile,
    pub orb_file: OrbFile,
    /// 本小时因取不到姿态/轨道而丢掉的候选数，由 `search` 写入。
    /// GRM 的 att/orb 是逐小时文件，尾端比事例流早收约 8 s（实测），
    /// 落在那一截里的候选只能丢——但要计数，不能静默。
    pub(super) dropped_no_ephemeris: AtomicUsize,
    /// 峰值时刻取不到姿态、姿态留空的候选数（位置有，候选照留），见 `search`
    pub(super) without_attitude: AtomicUsize,
    /// 候选窗里单路探测器占比过高（单路毛刺）而被否决的候选数，见 `search`
    pub(super) dropped_single_detector: AtomicUsize,
    /// 本小时因落在 GTI 之外而没有进入搜索的事例数，由 `search` 写入。
    /// 正常一小时是 SAA 缺口两端的两截半秒，合计几千到几万个。
    pub(super) events_outside_gti: AtomicUsize,
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

    /// 搜索前体检。
    fn exclusion(&self) -> Option<ExclusionReason> {
        // 大幅时间回跳＝同一段物理时间被记录了两次。实测 2025-09-19 02 时的
        // 一处：两份事例内容互不相同（逐行去重抓不到），重叠区间的计数率从
        // 1846 c/s 翻到 4036 c/s。排序修不了——率翻倍依旧；不排序更糟，搜索
        // 假设输入有序，乱序会让 max_duration 判据失效，实测把一个真实只有
        // 10 个事例的 1 ms 窗报成 count=315、fa 下溢为 0，畅通无阻地进判选。
        //
        // 小幅回跳不在此列，交给 `search` 排序：见 `REVERSAL_THRESHOLD`。
        if self.evt_file.time_reversals().max_magnitude > REVERSAL_THRESHOLD {
            return Some(ExclusionReason::UnorderedEvents);
        }
        if self.evt_file.into_iter().next().is_none() {
            return Some(ExclusionReason::NoEvents);
        }
        None
    }

    /// 曝光核算走 L1B 的 GTI。GRM 的小时文件并不覆盖整点到整点：实测
    /// 2025-09-19 全天 GTI 合计只有 86400 s 的 85.3%，单个小时低到 1778 s。
    /// 缺口全部落在南大西洋异常区，即 GTI 已经把 SAA 排除掉了，所以这里
    /// 不需要再叠一层 SAA 掩模。事例流在 `search` 里也按同一份 GTI 过滤，
    /// 分子分母口径一致。
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
        // 判没判都记，便于事后审计阈值
        let reversals = self.evt_file.time_reversals();
        let mut diagnostics = vec![
            ("time_reversals", reversals.count as f64),
            ("max_time_reversal", reversals.max_magnitude),
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
        let outside = self.events_outside_gti.load(Ordering::Relaxed);
        if outside > 0 {
            diagnostics.push(("events_outside_gti", outside as f64));
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
