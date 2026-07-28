use crate::algorithms::{config_guard, duplicate_guard};
use crate::io::level_1b::{SciFile, get_sci_filenames};
use crate::io::level_1k::{AttFile, EventFile, OrbitFile};
use crate::types::{Event, HxmtHe};
use blink_core::error::Error;
use blink_core::types::{Coverage, ExclusionReason, MissionElapsedTime};
use chrono::prelude::*;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicUsize, Ordering};

mod check_saturation;
mod from_epoch;
mod search;

pub type Interval = (MissionElapsedTime<HxmtHe>, MissionElapsedTime<HxmtHe>);

/// HE 的三个机箱。缺任何一个，该机箱的 FIFO reset 就无从检测，
/// 饱和掩模是残的。
pub const N_BOXES: usize = 3;

pub struct Chunk {
    pub event_file: EventFile,
    pub sci_files: Vec<(String, SciFile)>, // (box_name, file)
    pub stime_offsets: Vec<(String, f64)>, // (box_name, offset)
    pub orbit_file: OrbitFile,
    pub att_file: AttFile,
    pub span: [MissionElapsedTime<HxmtHe>; 2],
    /// 饱和区间缓存。求它要把三个机箱的 1B 全部重建一遍（本 pipeline 里最贵的
    /// 一步），而 search 和 coverage 都要用，所以只算一次。
    pub(super) saturation_cache: OnceLock<Vec<Interval>>,
    /// 本小时因为取不到姿态/轨道而丢掉的候选数，由 `search` 写入。
    /// 星历表两端已经会外推一小截，落到这里的是表真缺了一整段以上的情形。
    pub(super) dropped_no_ephemeris: AtomicUsize,
    /// 重复行占比缓存。体检和诊断都要它，而它是一遍全表扫描，只算一次。
    pub(super) duplicate_cache: OnceLock<Option<f64>>,
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

    /// 搜索前体检。任一条不过就整仪器跳过这一小时 —— 一处错，整机的结果就是错的。
    fn exclusion(&self) -> Option<ExclusionReason> {
        let channels = self.event_file.channels();
        if channels.is_empty() {
            return Some(ExclusionReason::NoEvents);
        }
        // 1B 缺机箱 → 该机箱的饱和窗未知。此时照常搜索等于拿一副残缺的饱和掩模
        // 去挡候选，漏挡的会当成真信号，是静默污染。
        if self.sci_files.len() < N_BOXES {
            return Some(ExclusionReason::MissingData);
        }
        // 星上阈值被调过的时段，能道解码是错的（见 config_guard）。
        if config_guard::is_nonstandard(channels) {
            return Some(ExclusionReason::NonstandardConfig);
        }
        // 事例表被整行重复写入的时段，显著性凭空虚高（见 duplicate_guard）。
        if self
            .duplicate_fraction()
            .is_some_and(|fraction| fraction > duplicate_guard::DUPLICATE_FRACTION_THRESHOLD)
        {
            return Some(ExclusionReason::DuplicatedEvents);
        }
        None
    }

    /// 注意：算 `masked_seconds` 要把三个机箱的 1B 全时间重建一遍。只对真正
    /// 要搜的小时调它 —— 被排除的小时曝光本来就是 0，没必要付这个代价。
    fn coverage(&self) -> Coverage {
        let (start, stop) = (self.span[0].met(), self.span[1].met());
        let masked_seconds = self
            .saturation_intervals()
            .iter()
            .map(|(interval_start, interval_stop)| {
                // 饱和窗扩过 ±1s，可能探出小时边界，先截到本小时内再累加
                (interval_stop.met().min(stop) - interval_start.met().max(start)).max(0.0)
            })
            .sum();

        Coverage {
            span_seconds: stop - start,
            masked_seconds,
        }
    }

    /// 在 `search` 之后取：`dropped_no_ephemeris` 是搜索过程中才产生的。
    fn diagnostics(&self) -> Vec<(&'static str, f64)> {
        // n_events 一并记：非标准配置判据的"统计量够不够"就是按它判的
        let mut diagnostics = vec![("n_events", self.event_file.channels().len() as f64)];
        if let Some(fraction) = self.config_gap_fraction() {
            diagnostics.push(("config_gap_fraction", fraction));
        }
        if let Some(fraction) = self.duplicate_fraction() {
            diagnostics.push(("duplicate_fraction", fraction));
        }
        let dropped = self.dropped_no_ephemeris.load(Ordering::Relaxed);
        if dropped > 0 {
            diagnostics.push(("dropped_no_ephemeris", dropped as f64));
        }
        diagnostics
    }

    fn last_modified(epoch: &DateTime<Utc>) -> Result<DateTime<Utc>, Error> {
        let sci_last_modifieds: Vec<DateTime<Utc>> = get_sci_filenames(*epoch)
            .iter()
            .map(|(_, filename)| {
                let last_modified = std::fs::metadata(filename)?.modified()?;
                let datetime: DateTime<Utc> = last_modified.into();
                Ok::<DateTime<Utc>, Error>(datetime)
            })
            .collect::<Result<Vec<DateTime<Utc>>, Error>>()?;

        let other_last_modifieds: Vec<DateTime<Utc>> = vec![
            EventFile::last_modified(epoch)?,
            OrbitFile::last_modified(epoch)?,
            AttFile::last_modified(epoch)?,
        ];

        let last_modifieds: Vec<DateTime<Utc>> = sci_last_modifieds
            .into_iter()
            .chain(other_last_modifieds)
            .collect();

        let max_last_modified = last_modifieds
            .into_iter()
            .max()
            .ok_or_else(|| Error::FileNotFound("No files found".to_string()))?;

        Ok(max_last_modified)
    }
}

impl Chunk {
    /// 本小时非标准配置判据的实测值（gap 道段占比），事例表为空时为 None。
    /// 无论判决与否都记录下来，便于事后审计阈值。
    pub fn config_gap_fraction(&self) -> Option<f64> {
        config_guard::gap_fraction(self.event_file.channels())
    }

    /// 本小时重复行判据的实测值。一遍全表扫描，体检和诊断共用同一份。
    pub fn duplicate_fraction(&self) -> Option<f64> {
        *self.duplicate_cache.get_or_init(|| {
            duplicate_guard::duplicate_fraction(
                self.event_file.times(),
                self.event_file.det_ids(),
                self.event_file.channels(),
                self.event_file.pulse_widths(),
                self.event_file.event_types(),
            )
        })
    }
}
