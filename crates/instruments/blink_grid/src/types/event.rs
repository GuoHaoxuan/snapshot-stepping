use blink_core::types::MissionElapsedTime;
use serde::Serialize;

use crate::types::instrument::{Grid, Satellite};

/// 低能截止（keV）。
///
/// 实测 v03 产品的谱从 ch5–7（21–29 keV）才起，ch12（56 keV）达峰，往上平滑
/// 下降，低端没有 GRM 那种噪声堆——硬件阈值本身就切在 30 keV 附近。这里取
/// 30 keV 是把阈值放在硬件截止之上、不再额外砍统计的位置。用能量而不用道号，
/// 是因为各星、各版本的道能量不同（03B/04/07 是 97 道，02 是 98 道，同一道号
/// 差几个 keV）。尚未用真实 TGF 定标，见 `OPEN-QUESTIONS.md`。
pub const ENERGY_THRESHOLD_KEV: f32 = 30.0;

#[derive(Serialize, Debug, Clone)]
#[serde(bound(serialize = ""))]
pub struct Event<S: Satellite> {
    pub time: MissionElapsedTime<Grid<S>>,
    /// `PI`，1 起算
    pub channel: i16,
    /// 0–3，对应 `EVENTS0..3`
    pub detector: u8,
    /// `EVT_TYPE`。实测只见 1 和 2：2 的 `PI` 恒等于最高道（97），是溢出/饱和
    /// 一类的标记而非光子；3 从未出现。
    pub evt_type: u8,
    /// 本道下限能量（keV），由所在文件的 `EBOUNDS` 查得
    pub energy_kev: f32,
    /// 是否落在最高道（溢出道）
    pub overflow: bool,
}

impl<S: Satellite> blink_core::traits::Event for Event<S> {
    type Instrument = Grid<S>;
    type ChannelType = i16;

    fn time(&self) -> MissionElapsedTime<Self::Instrument> {
        self.time
    }

    fn channel(&self) -> Self::ChannelType {
        self.channel
    }

    /// 单组：4 个探测器同型，合成一路搜。
    fn group(&self) -> u8 {
        0
    }

    fn keep(&self) -> bool {
        self.evt_type == 1 && !self.overflow && self.energy_kev >= ENERGY_THRESHOLD_KEV
    }
}

impl<S: Satellite> Ord for Event<S> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.time.cmp(&other.time)
    }
}
impl<S: Satellite> PartialOrd for Event<S> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<S: Satellite> PartialEq for Event<S> {
    fn eq(&self, other: &Self) -> bool {
        self.time == other.time
    }
}
impl<S: Satellite> Eq for Event<S> {}
