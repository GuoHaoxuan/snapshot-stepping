//! 归档里找文件。
//!
//! GBM 的镜像有三种并存的布局，而且 NaI 与 BGO 的连续数据不在一处：
//!
//! * `{ROOT}/{YYYY}/{MM}/{DD}/current/` —— 2015 年起的主目录。NaI 的逐小时
//!   TTE 只在这里，且要到 2017-10 才开始有；2020 年这里只剩 12 个 NaI。
//! * `{ROOT}/{YYYY}/{YYMMDD}/tte/` —— 2014 年及以前的主目录，里面是按触发
//!   命名的短文件，没有逐小时 TTE。poshist 在这里。
//! * `{ROOT}/BGO/{YYYY}/{YYMMDD}/` —— BGO 的逐小时 TTE 单独存放，2012-11
//!   起一直到 2020，是覆盖最长的一份。
//!
//! 所以查找一个探头时按上面的顺序逐个目录试，取版本号最大的那个文件。

use chrono::prelude::*;
use std::path::PathBuf;
use std::{env, fs, sync::LazyLock};

static ROOT: LazyLock<String> = LazyLock::new(|| {
    env::var("FERMI_GBM_DIR").unwrap_or_else(|_| "/hxmtfs/data/Fermi_GBM".to_string())
});

/// 一天可能落在的目录，按优先级排列。
fn day_dirs(time: &DateTime<Utc>) -> Vec<PathBuf> {
    let root = ROOT.as_str();
    let ymd = time.format("%y%m%d").to_string();
    vec![
        PathBuf::from(format!("{root}/{}/current", time.format("%Y/%m/%d"))),
        PathBuf::from(format!("{root}/{}/{ymd}/tte", time.format("%Y"))),
        PathBuf::from(format!("{root}/BGO/{}/{ymd}", time.format("%Y"))),
    ]
}

/// 目录里以 `prefix` 开头、版本号最大的文件。
fn latest_version(dir: &PathBuf, prefix: &str) -> Option<PathBuf> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .ok()?
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.file_name().to_string_lossy().to_string())
        .filter(|name| name.starts_with(prefix))
        .collect();
    // 文件名尾部形如 `_vNN.fit.gz`，版本号定长，字典序即版本序
    names.sort();
    names.last().map(|name| dir.join(name))
}

fn find(time: &DateTime<Utc>, prefix: &str) -> Option<PathBuf> {
    day_dirs(time)
        .iter()
        .find_map(|dir| latest_version(dir, prefix))
}

/// 某探头某小时的 TTE。`detector` 是归档里的代号，如 `n0` / `b1`。
pub fn find_tte(time: &DateTime<Utc>, detector: &str) -> Option<PathBuf> {
    find(
        time,
        &format!(
            "glg_tte_{detector}_{}_{:02}z_v",
            time.format("%y%m%d"),
            time.hour()
        ),
    )
}

/// 当天的 poshist（整天一个文件，1 Hz）。
pub fn find_poshist(time: &DateTime<Utc>) -> Option<PathBuf> {
    find(time, &format!("glg_poshist_all_{}_v", time.format("%y%m%d")))
}
