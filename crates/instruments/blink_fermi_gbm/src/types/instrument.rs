use blink_core::traits::Instrument;
use chrono::prelude::*;
use std::{str::FromStr, sync::OnceLock};

/// Fermi Gamma-ray Burst Monitor (GBM)
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct FermiGbm;

impl Instrument for FermiGbm {
    type Chunk = crate::types::Chunk;

    /// MET 的零点。
    ///
    /// 文件头给的是 `MJDREFI=51910, MJDREFF=0.00074287037, TIMESYS='TT'`，
    /// 而 MJDREFF 那一截正好是 64.184 s = 2001 年的 TT−UTC，所以 MJDREF 精确
    /// 落在 2001-01-01T00:00:00 UTC。TIME 走的是连续 TT 秒，2001 年之后的 5 个
    /// 闰秒（2005/2008/2012/2015/2016）由 `MissionElapsedTime` 的闰秒表补回。
    ///
    /// 实测对账：`TSTART=567993485.00096` 对应 `DATE-OBS='2018-12-31T23:58:00'`，
    /// 按 UTC 直接数是 567993480 s，差的 5 s 正是那 5 个闰秒。
    fn ref_time() -> &'static DateTime<Utc> {
        static REF_TIME: OnceLock<DateTime<Utc>> = OnceLock::new();
        REF_TIME
            .get_or_init(|| DateTime::<Utc>::from_str("2001-01-01T00:00:00.000000000 UTC").unwrap())
    }

    /// 连续 TTE 的起点，不是卫星发射日。GBM 2008 年就在天上，但逐小时的
    /// continuous TTE 要到 2012-11 才常态化；在那之前归档里只有触发前后的
    /// 短文件，盲搜无从谈起。
    fn launch_day() -> NaiveDate {
        NaiveDate::from_ymd_opt(2012, 11, 1).unwrap()
    }

    fn name() -> &'static str {
        "Fermi/GBM"
    }
}
