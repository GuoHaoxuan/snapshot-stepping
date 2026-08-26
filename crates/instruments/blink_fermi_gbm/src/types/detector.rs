use serde::Serialize;

/// GBM 的两类探测器。
///
/// 二者的响应差得很远，不能合成一路搜：NaI 是 1.27 cm 薄片、能段 4.5–2000 keV；
/// BGO 是 12.7 cm 厚柱、能段 113–50000 keV。TGF 谱硬，光子主要落在 BGO 里，
/// 而实测本底是 12 个 NaI 合计约 15600 c/s、2 个 BGO 合计约 5000 c/s——合并
/// 等于拿 BGO 的信号去配四倍本底。所以两类各成一组，见 `Chunk::search`。
#[derive(Serialize, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Detector {
    /// n0–n9, na, nb 共 12 个
    Nai,
    /// b0, b1 共 2 个
    Bgo,
}

impl Detector {
    /// 归档文件名里的探头代号。
    pub const NAI_NAMES: [&'static str; 12] = [
        "n0", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9", "na", "nb",
    ];
    pub const BGO_NAMES: [&'static str; 2] = ["b0", "b1"];

    pub fn from_name(name: &str) -> Option<Self> {
        match name.as_bytes().first() {
            Some(b'n') => Some(Self::Nai),
            Some(b'b') => Some(Self::Bgo),
            _ => None,
        }
    }

    pub fn names(&self) -> &'static [&'static str] {
        match self {
            Self::Nai => &Self::NAI_NAMES,
            Self::Bgo => &Self::BGO_NAMES,
        }
    }
}
