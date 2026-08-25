use blink_core::types::MissionElapsedTime;
use serde::Serialize;

use crate::types::instrument::SvomGrm;

#[derive(Serialize, Debug, Clone)]
pub struct Event {
    pub time: MissionElapsedTime<SvomGrm>,
    pub channel: i16,
    pub detector_id: u8,
    pub gain_type: u8,
    pub dead_time: f32,
    pub evt_type: u8,
    pub anti_coin: u8,
    pub flag: u8,
}

impl blink_core::traits::Event for Event {
    type Instrument = SvomGrm;
    type ChannelType = i16;

    fn time(&self) -> MissionElapsedTime<Self::Instrument> {
        self.time
    }

    fn channel(&self) -> Self::ChannelType {
        self.channel
    }

    fn group(&self) -> u8 {
        0
    }

    /// 事例准入。判据按 GRM 自己的数据定标，不照抄 HXMT（见
    /// 本 crate 的 `OPEN-QUESTIONS.md`）：
    ///
    /// * `EVT_TYPE == 1` 是量程溢出事例（实测 PI 中位 256、`GAIN_TYPE` 全为
    ///   低增益），能量不可信。
    /// * `PI >= 256` 同样落在 EBOUNDS 顶端的溢出道。
    /// * `PI >= CHANNEL_THRESHOLD` 是低能截止。TGF 是硬谱，而实测未过滤时
    ///   一小时会冒出 15–30 个候选，其窗内事例 PI 中位 8–24（低于全局中位
    ///   28）、溢出占比 16–22%（全局 0.83%）——软到不可能是 TGF。加上本阈值
    ///   后两个独立小时分别降到 2 个（即 GRB 250919A 自身）和 0 个，而 GRB
    ///   仍被抓到。
    ///
    /// 阈值本身尚未用真实 TGF 能谱定标，只是有实测支撑的起点。
    fn keep(&self) -> bool {
        // ch15 ≈ 22.5 keV（EBOUNDS 实测）。取 15 而非 HXMT 的 38：GRM 的
        // ch38 已是 78 keV，切得过狠；且实测 PI≥50 时候选数不降反升——本底
        // 压得过低，少数事例就能触发。甜区在 15–30。
        const CHANNEL_THRESHOLD: i16 = 15;
        // EBOUNDS 顶端两道（257/258）是 10–20 MeV 溢出道。
        const OVERFLOW_CHANNEL: i16 = 256;

        self.evt_type == 0 && self.channel < OVERFLOW_CHANNEL && self.channel >= CHANNEL_THRESHOLD
    }
}

impl Ord for Event {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.time.cmp(&other.time)
    }
}
impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl PartialEq for Event {
    fn eq(&self, other: &Self) -> bool {
        self.time == other.time
    }
}
impl Eq for Event {}
