use chrono::{DateTime, TimeDelta, Utc};
use serde::{Deserialize, Serialize};
use uom::si::f64::*;

use crate::{
    traits::{Event, Instrument},
    types::{Attitude, MissionElapsedTime, Position},
};

/// 反符合探测器（ACD）符合计数，候选生成时从事例流现场统计——事后无法从
/// 候选表复原，必须随产物保存。只存原始计数不存比例：保泊松误差，阈值可
/// 复议而不必重跑全量。仅 HXMT HE 填写；无 ACD 的仪器留 `None`。
#[derive(Clone, Serialize, Deserialize)]
pub struct AcdCounts {
    /// 候选窗 [start, stop] 内 kept 事例总数。自带分母：不依赖 `count`
    /// 的分箱语义，占比 n_acd/n 自洽。
    pub n: u32,
    /// 候选窗内任意一块 ACD 着火的事例数
    pub n_acd: u32,
    /// 候选窗内 ≥2 块 ACD 同时着火的事例数（区分真带电粒子与偶然符合）
    pub n_acd_multi: u32,
    /// 邻域基线窗内 kept 事例总数
    pub n_bg: u32,
    /// 邻域基线窗内任意一块 ACD 着火的事例数
    pub n_acd_bg: u32,
}

#[derive(Serialize, Deserialize)]
pub struct Signal<E: Event> {
    pub start: MissionElapsedTime<E::Instrument>,
    pub stop: MissionElapsedTime<E::Instrument>,
    pub bin_size_min: Time,
    pub bin_size_max: Time,
    pub bin_size_best: Time,
    pub delay: Time,
    pub count: u32,
    pub mean: f64,
    pub sf: f64,
    pub false_positive_per_year: f64,
    pub attitude: Attitude,
    pub position: Position,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub acd: Option<AcdCounts>,
}

impl<E: Event> Signal<E> {
    pub fn to_unified(&self) -> UnifiedSignal {
        UnifiedSignal {
            start: self.start.to_utc(),
            stop: self.stop.to_utc(),
            bin_size_min: self.bin_size_min,
            bin_size_max: self.bin_size_max,
            bin_size_best: self.bin_size_best,
            delay: self.delay,
            count: self.count,
            mean: self.mean,
            sf: self.sf,
            false_positive_per_year: self.false_positive_per_year,
            attitude: self.attitude.clone(),
            position: self.position.clone(),
            instrument: <E::Instrument as Instrument>::name().to_string(),
            acd: self.acd.clone(),
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct UnifiedSignal {
    pub start: DateTime<Utc>,
    pub stop: DateTime<Utc>,
    pub bin_size_min: Time,
    pub bin_size_max: Time,
    pub bin_size_best: Time,
    pub delay: Time,
    pub count: u32,
    pub mean: f64,
    pub sf: f64,
    pub false_positive_per_year: f64,
    pub attitude: Attitude,
    pub position: Position,
    pub instrument: String,
    /// `default` 兼容旧 signals.json（无此字段 → None），None 不序列化。
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub acd: Option<AcdCounts>,
}

impl UnifiedSignal {
    pub fn peak_time(&self) -> DateTime<Utc> {
        self.start
            + TimeDelta::nanoseconds(self.delay.get::<uom::si::time::nanosecond>().round() as i64)
            + TimeDelta::nanoseconds(
                (self.bin_size_best / 2.0)
                    .get::<uom::si::time::nanosecond>()
                    .round() as i64,
            )
    }
}
