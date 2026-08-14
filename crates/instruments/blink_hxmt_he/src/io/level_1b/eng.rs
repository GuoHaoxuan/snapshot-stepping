use blink_core::error::Error;

/// 众数至少要有多少行支持才可信。健康小时全部行都是同一个值（数千行）；
/// 平台 UTC 冻结时 UTC 不动、sTime 照走，offset 逐行 −1 递减、互不相同，
/// 永远形不成竞争众数。定 10 行（连续 10 秒健康数据）是为了防极短文件里
/// 两三行偶然同值就当真。不满足即整小时拒绝（UtcFreeze 入账本）。
const MIN_MODE_SUPPORT: usize = 10;

/// 从工程数据文件中读取 stime→UTC 的固定偏移量。
///
/// 工程文件每秒一包，包中 `UTC_Last_Bdc` 和 `sTime_Last_Bdc` 列
/// 给出精确的 UTC↔stime 映射。offset = UTC - stime，在整个小时内恒定。
///
/// 取众数而非首行：平台 UTC 广播存在冻结时段（实测 2024-05-11 G5 风暴期
/// 冻结 1963 s，三箱同值），冻结行的 offset 是错的且逐行漂移；首行恰逢
/// 冻结就会把整小时时间线平移上千秒（Gannon T20 实测 −1042 s）。
pub fn read_stime_offset(filename: &str) -> Result<f64, Error> {
    let mut fptr = fitsio::FitsFile::open(filename)?;
    let hdu = fptr.hdu("HE_Eng")?;

    let utc: Vec<i64> = hdu.read_col(&mut fptr, "UTC_Last_Bdc")?;
    let stime: Vec<i64> = hdu.read_col(&mut fptr, "sTime_Last_Bdc")?;

    if utc.is_empty() || stime.is_empty() {
        return Err(Error::InvalidData("Empty eng data".into()));
    }

    offset_mode(&utc, &stime)
}

/// utc−stime 的众数（并列时取较小值，保证确定性）。
fn offset_mode(utc: &[i64], stime: &[i64]) -> Result<f64, Error> {
    let mut counts: std::collections::BTreeMap<i64, usize> = std::collections::BTreeMap::new();
    for (u, s) in utc.iter().zip(stime.iter()) {
        *counts.entry(u - s).or_insert(0) += 1;
    }
    // BTreeMap 升序遍历 + 严格大于：并列众数取最小 offset，结果确定。
    let (mode, n_mode) = counts
        .iter()
        .fold((0i64, 0usize), |best, (&off, &n)| {
            if n > best.1 { (off, n) } else { best }
        });

    if n_mode < MIN_MODE_SUPPORT {
        return Err(Error::TimeReferenceInvalid(format!(
            "utc-stime offset mode support {n_mode} < {MIN_MODE_SUPPORT} over {} rows",
            utc.len().min(stime.len())
        )));
    }
    Ok(mode as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 健康小时：全部行同值。
    #[test]
    fn healthy_hour_returns_constant() {
        let stime: Vec<i64> = (0..3600).collect();
        let utc: Vec<i64> = stime.iter().map(|s| s + 387527433).collect();
        assert_eq!(offset_mode(&utc, &stime).unwrap(), 387527433.0);
    }

    /// Gannon T20 形态：文件开头 921 行 UTC 冻结（同一 UTC 值、sTime 照走，
    /// offset 逐行 −1），其余行健康。首行 offset = 众数 − 1042，取首行即整
    /// 小时平移；取众数必须给出健康值。
    #[test]
    fn frozen_head_is_ignored() {
        let mode = 387527433i64;
        let frozen_utc = 390081362i64;
        let mut utc = Vec::new();
        let mut stime = Vec::new();
        for i in 0..921 {
            utc.push(frozen_utc);
            stime.push(frozen_utc - mode + 1042 + i); // 首行 offset = mode - 1042
        }
        for i in 0..1873 {
            let s = frozen_utc - mode + 1042 + 921 + i;
            stime.push(s);
            utc.push(s + mode);
        }
        assert_ne!(utc[0] - stime[0], mode); // 首行确实是坏的
        assert_eq!(offset_mode(&utc, &stime).unwrap(), mode as f64);
    }

    /// 整小时全冻结：offset 逐行漂移、无稳定众数 → 必须拒绝而不是给个错值。
    #[test]
    fn fully_frozen_hour_is_rejected() {
        let frozen_utc = 390081362i64;
        let stime: Vec<i64> = (0..3600).map(|i| 2554971 + i).collect();
        let utc: Vec<i64> = vec![frozen_utc; 3600];
        assert!(matches!(
            offset_mode(&utc, &stime),
            Err(Error::TimeReferenceInvalid(_))
        ));
    }
}
