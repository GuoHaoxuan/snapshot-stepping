use chrono::prelude::*;
use std::fs;
use std::path::PathBuf;

const DIR_PREFIX: &str = "/gecamfs/SVOM/Archived-DATA/GRM-DATA/L1B/daily/";

/// 归档根。默认是 IHEP 上的挂载点，`SVOM_ARCHIVE` 可覆盖——为的是能拿一小时
/// 的样本文件在本地跑通整条链路（与 `WWLLN_DB_PATH` 同一套做法）。
fn dir_prefix() -> String {
    std::env::var("SVOM_ARCHIVE").unwrap_or_else(|_| DIR_PREFIX.to_string())
}

fn subdir(time: &DateTime<Utc>, sub: &str) -> String {
    let prefix = dir_prefix();
    let separator = if prefix.ends_with('/') { "" } else { "/" };
    format!("{}{}{}/{}/", prefix, separator, time.format("%Y/%m/%d"), sub)
}

fn evt_dir(time: &DateTime<Utc>) -> String {
    subdir(time, "grm_evt")
}

fn att_dir(time: &DateTime<Utc>) -> String {
    subdir(time, "att")
}

fn orb_dir(time: &DateTime<Utc>) -> String {
    subdir(time, "orb")
}

const EVT_TEMPLATE: &str = "svom_grm_evt_%y%m%d_%H_v";
const ATT_TEMPLATE: &str = "svom_att_%y%m%d_%H_v";
const ORB_TEMPLATE: &str = "svom_orb_%y%m%d_%H_v";

fn find_by_time(
    dir: &str,
    template: &str,
    time: &DateTime<Utc>,
) -> Result<PathBuf, std::io::Error> {
    let mut files: Vec<_> = fs::read_dir(dir)?
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.file_name().to_string_lossy().to_string())
        .collect();

    let file_str = time.format(template).to_string();
    files.retain(|f| f.starts_with(&file_str));
    files.sort();

    files
        .last()
        .map(|filename| PathBuf::from(dir).join(filename))
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "No file found for the given time.",
            )
        })
}

pub fn find_evt_by_time(time: &DateTime<Utc>) -> Result<PathBuf, std::io::Error> {
    find_by_time(&evt_dir(time), EVT_TEMPLATE, time)
}

pub fn find_att_by_time(time: &DateTime<Utc>) -> Result<PathBuf, std::io::Error> {
    find_by_time(&att_dir(time), ATT_TEMPLATE, time)
}

pub fn find_orb_by_time(time: &DateTime<Utc>) -> Result<PathBuf, std::io::Error> {
    find_by_time(&orb_dir(time), ORB_TEMPLATE, time)
}
