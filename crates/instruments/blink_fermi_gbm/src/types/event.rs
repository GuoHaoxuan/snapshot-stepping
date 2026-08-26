use blink_core::types::MissionElapsedTime;
use serde::Serialize;

use crate::types::{Detector, FermiGbm};

#[derive(Serialize, Debug, Clone)]
pub struct Event {
    pub time: MissionElapsedTime<FermiGbm>,
    pub channel: i16,
    pub detector: Detector,
    /// 本事例所属的搜索分组下标。由 `Chunk` 在载入时按这一小时**实际存在**的
    /// 探测器类型连续编号——只有 BGO 的时段编成一组（下标恒为 0），NaI 与 BGO
    /// 都在的时段编成两组。这样 Bonferroni 的试验次数惩罚才跟真实组数一致，
    /// 不会给单类型的时段白白乘上 2。
    pub group: u8,
}

impl blink_core::traits::Event for Event {
    type Instrument = FermiGbm;
    type ChannelType = i16;

    fn time(&self) -> MissionElapsedTime<Self::Instrument> {
        self.time
    }

    fn channel(&self) -> Self::ChannelType {
        self.channel
    }

    fn group(&self) -> u8 {
        self.group
    }

    /// 事例准入。
    ///
    /// 两端的道是溢出道：0 收所有低于阈值的、127 收所有高于量程的，能量都不
    /// 可信。中间怎么切还没有用数据定过——GRM 那边实测低能噪声会造出成堆的
    /// 假候选，GBM 是否同样、阈值该落在哪一道，要等第一批候选出来再定标。
    /// 见本 crate 的 `OPEN-QUESTIONS.md`。
    fn keep(&self) -> bool {
        const OVERFLOW_LOW: i16 = 0;
        const OVERFLOW_HIGH: i16 = 127;

        self.channel > OVERFLOW_LOW && self.channel < OVERFLOW_HIGH
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
