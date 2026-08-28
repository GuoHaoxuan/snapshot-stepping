use serde::Serialize;

/// GBM 的两类探测器。
///
/// 二者的响应差得很远，不能合成一路搜：NaI 是 1.27 cm 薄片、能段 4.5–2000 keV；
/// BGO 是 12.7 cm 厚柱、能段 113–50000 keV。TGF 谱硬，光子主要落在 BGO 里，
/// 而实测本底是 12 个 NaI 合计约 15600 c/s、2 个 BGO 合计约 5000 c/s——合并
/// 等于拿 BGO 的信号去配四倍本底。所以两类各成一组，见 `Chunk::search`。
/// 分组即按这个枚举来，与 Fermi/GBM 团队的做法一致：12 个 NaI 合成一组，
/// 两个 BGO 各自一组，共三组。实测支持这样分——一天 12.6 万个候选里 72.6%
/// 是宇宙线穿过整星，它在 12 个 NaI 里能凑够 min_number，但两个 BGO 各只
/// 收到一个事例，要求两组符合就能挡住；而真信号在 12 个 NaI 里是累加的，
/// 若按探头拆成 14 组反而会把它摊薄到每组都够不着阈值。
#[derive(Serialize, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Detector {
    /// n0–n9, na, nb 共 12 个，合成一组
    Nai,
    /// b0，自成一组
    Bgo0,
    /// b1，自成一组
    Bgo1,
}

impl Detector {
    /// 归档文件名里的探头代号。
    pub const NAI_NAMES: [&'static str; 12] = [
        "n0", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9", "na", "nb",
    ];
    pub const BGO_NAMES: [&'static str; 2] = ["b0", "b1"];

    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "b0" => Some(Self::Bgo0),
            "b1" => Some(Self::Bgo1),
            n if n.starts_with('n') => Some(Self::Nai),
            _ => None,
        }
    }

    pub fn names(&self) -> &'static [&'static str] {
        match self {
            Self::Nai => &Self::NAI_NAMES,
            Self::Bgo0 => &["b0"],
            Self::Bgo1 => &["b1"],
        }
    }

    /// 分组顺序。缺哪一类就少哪一组，`Chunk` 按实际到齐的类型连续编号。
    pub const ALL: [Self; 3] = [Self::Nai, Self::Bgo0, Self::Bgo1];
}
