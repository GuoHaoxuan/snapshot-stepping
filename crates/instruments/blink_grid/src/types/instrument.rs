use blink_core::traits::Instrument;
use chrono::prelude::*;
use std::{marker::PhantomData, str::FromStr, sync::OnceLock};

/// 天格计划的一颗星。
///
/// 每颗星 4 个 GAGG+SiPM 探测器、同一套文件格式，归档按星分目录
/// （`GRID-03B/fits7/...`），文件名带星的代号（`G03_evt_...`）。
pub trait Satellite:
    Clone + Copy + PartialEq + Eq + PartialOrd + Ord + std::fmt::Debug + Send + Sync + 'static
{
    /// 归档里的目录名，也是输出目录名
    const DIR: &'static str;
    /// 文件名前缀
    const CODE: &'static str;
    /// 事例产品（`fits7`）的第一天。不是发射日：盲搜从有数据的那天起算。
    fn first_day() -> NaiveDate;
}

macro_rules! satellite {
    ($ty:ident, $dir:literal, $code:literal, ($y:literal, $m:literal, $d:literal), $doc:literal) => {
        #[doc = $doc]
        #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
        pub struct $ty;
        impl Satellite for $ty {
            const DIR: &'static str = $dir;
            const CODE: &'static str = $code;
            fn first_day() -> NaiveDate {
                NaiveDate::from_ymd_opt($y, $m, $d).unwrap()
            }
        }
    };
}

satellite!(
    Sat02,
    "GRID-02",
    "G02",
    (2020, 11, 14),
    "GRID-02：事例产品 2020-11-14 .. 2021-03-08，75 天"
);
satellite!(
    Sat03B,
    "GRID-03B",
    "G03",
    (2022, 3, 11),
    "GRID-03B：事例产品 2022-03-11 .. 2024-08-01，789 天"
);
satellite!(
    Sat04,
    "GRID-04",
    "G04",
    (2022, 3, 11),
    "GRID-04：事例产品 2022-03-11 .. 2024-08-22，527 天（2023-03..12 缺）"
);
satellite!(
    Sat07,
    "GRID-07",
    "G07",
    (2024, 1, 1),
    "GRID-07：事例产品 2024-01-01 .. 2024-07-22，195 天"
);

/// 某颗天格星作为一台仪器。
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Grid<S: Satellite>(PhantomData<S>);

impl<S: Satellite> Instrument for Grid<S> {
    type Chunk = crate::types::Chunk<S>;

    /// 事例表 `TIME` 的零点。文件头 `DATE_REF = '2018-01-01T00:00:00.000'`；
    /// 对账：`TIME = 132217726.086` ↔ `DATE_OBS = '2022-03-11T07:08:46.083'`，
    /// 按 UTC 直接数 1530 天 + 25726 s = 132217726 s，完全一致。2017 年以后
    /// 没有新闰秒，`MissionElapsedTime` 的闰秒表在这里恒为零修正。
    fn ref_time() -> &'static DateTime<Utc> {
        static REF_TIME: OnceLock<DateTime<Utc>> = OnceLock::new();
        REF_TIME
            .get_or_init(|| DateTime::<Utc>::from_str("2018-01-01T00:00:00.000000000 UTC").unwrap())
    }

    fn launch_day() -> NaiveDate {
        S::first_day()
    }

    fn name() -> &'static str {
        S::DIR
    }
}

pub type Grid02 = Grid<Sat02>;
pub type Grid03B = Grid<Sat03B>;
pub type Grid04 = Grid<Sat04>;
pub type Grid07 = Grid<Sat07>;
