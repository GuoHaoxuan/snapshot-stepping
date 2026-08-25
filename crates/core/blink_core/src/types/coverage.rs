use crate::error::Error;
use serde::{Deserialize, Serialize};

/// 一个搜索单元（对 HXMT HE 是一小时）被排除、不进入搜索的原因。
///
/// 排除必须是显式记录的：上层把每个单元要么记成 searched、要么记成
/// excluded(reason)，不允许出现"既没搜也没记"的静默缺口 —— 否则事件率
/// 和误报率的分母是错的。
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum ExclusionReason {
    /// 输入文件缺失：归档里没有这个单元，或只有一部分（如三机箱缺其一）。
    MissingData,
    /// 文件在但读不出来：FITS 结构损坏、I/O 失败。
    CorruptData,
    /// 仪器处于非标准配置，既有解码规则不适用，搜出来的候选是假的。
    NonstandardConfig,
    /// 事例表里有整行重复写入：重复把暴发和本底一起放大，显著性凭空虚高。
    DuplicatedEvents,
    /// 事例表内部时间不单调。同一段物理时间被记录了两次（两份事例内容不同，
    /// 逐行去重抓不到），合并后计数率翻倍；而搜索假设输入按时间有序，乱序会
    /// 让窗长判据失效、计数虚涨。排序救不了——率翻倍依然在。
    UnorderedEvents,
    /// 文件可读但没有事例。
    NoEvents,
    /// 平台 UTC 广播冻结/漂移，utc−stime 找不到稳定众数：时间基准无法建立，
    /// 1B 时间重建与饱和掩模不可信。
    UtcFreeze,
}

impl ExclusionReason {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::MissingData => "missing_data",
            Self::CorruptData => "corrupt_data",
            Self::NonstandardConfig => "nonstandard_config",
            Self::DuplicatedEvents => "duplicated_events",
            Self::UnorderedEvents => "unordered_events",
            Self::NoEvents => "no_events",
            Self::UtcFreeze => "utc_freeze",
        }
    }
}

impl std::fmt::Display for ExclusionReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl From<&Error> for ExclusionReason {
    fn from(error: &Error) -> Self {
        match error {
            Error::FileNotFound(_) => Self::MissingData,
            // 目录整个不存在时是 read_dir 抛出的 io::ErrorKind::NotFound
            //（例如 1K 缺当天的产品目录），那是缺数据、不是坏数据。
            Error::Io(e) if e.kind() == std::io::ErrorKind::NotFound => Self::MissingData,
            Error::TimeReferenceInvalid(_) => Self::UtcFreeze,
            _ => Self::CorruptData,
        }
    }
}

/// 一个搜索单元的曝光时间核算，供事件率 / 误报率的分母使用。
///
/// 只记可测的量：单元名义跨度、单元内被算法屏蔽掉的时间。不在这里猜
/// "数据空档"——那需要另立判据，猜出来的曝光比没有更糟。
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Coverage {
    /// 单元名义跨度（秒）。HXMT HE = 3600。
    pub span_seconds: f64,
    /// 单元内被算法屏蔽、不参与搜索的时间（秒）。
    /// HXMT HE = 饱和窗（FIFO reset 空洞扩 ±1s 后求并、再截到单元跨度内）。
    pub masked_seconds: f64,
}

impl Coverage {
    /// 实际进入搜索的时间（秒）= 跨度 − 屏蔽。
    pub fn searched_seconds(&self) -> f64 {
        (self.span_seconds - self.masked_seconds).max(0.0)
    }
}
