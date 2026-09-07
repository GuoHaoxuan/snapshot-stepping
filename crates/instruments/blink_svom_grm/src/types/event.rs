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
        // 能阈 ch25（EBOUNDS 实测 42 keV）。GRM 硬件下限是 ch10（15 keV），本底谱在
        // ch11–19 有个鼓包（每道计数是 ch23 以后的 2–3 倍），是平稳的软本底、不成簇，
        // 不会自己造候选，但会稀释信噪。用 72 个闪电证实的 TGF 扫能阈：证实 TGF 的
        // 窗内谱在 ch10–24 都低于本底（窗内/本底 0.19、0.33、0.61），从 ch25 起反超；
        // ch15→ch25 少 15% 的 TGF 计数、少 43% 的本底，信噪中位 12.4→13.9，达到
        // p≤1e-8 的比例 85%→93%；ch25–30 是平台，ch40 以上信号计数掉得快。早先的
        // ch15 是经验值（HXMT 用 ch38，GRM 的 ch38 已是 78 keV 切得过狠；实测 PI≥50
        // 候选数不降反升）。见 scripts/plot_svom_threshold.py。
        const CHANNEL_THRESHOLD: i16 = 25;
        // EBOUNDS 顶端两道（257/258）是 10–20 MeV 溢出道。
        const OVERFLOW_CHANNEL: i16 = 256;

        // ANTI_COIN=1 是星上标定源事例，不是反符合标志：绝对速率恒定在 15–16 c/s
        // 与总计数率无关，能谱是一条 49–57 keV 的线（占这类事例的 54%），SAA 内
        // 反而比 SAA 外低。它们不是天体光子，不该进搜索。ch25 之上占 1.56%，
        // 排除后本底降 1.56%、信号最多损 0.24%。见 evidence/anticoin/。
        self.evt_type == 0
            && self.anti_coin == 0
            && self.channel < OVERFLOW_CHANNEL
            && self.channel >= CHANNEL_THRESHOLD
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
