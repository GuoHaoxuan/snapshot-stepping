//! 归档里找文件。
//!
//! 布局：`{ROOT}/{DIR}/fits7/{YYYY}/{MM}/{DD}/evt_vXX_YY/` 是一天的事例产品，
//! 里面按过境分文件（`G03_evt_YYMMDDhhmm_YYMMDDhhmm_vXX_YY.fits`，一次过境
//! 10–25 分钟，一天 1–21 次）；`fits8/.../posatt_vXX_YY_ZZ/` 是同一批过境的
//! 1 Hz 位姿，文件名的时间段与事例文件一一对应。
//!
//! 一天常有多个版本目录（`evt_v02_02`、`evt_v03_00`、`evt_v03_01` …），版本
//! 高的过境更全（实测 2023-08-12 的 v03_00 只有 1 次过境，v03_01 有 21 次），
//! 所以取字典序最大的那个。归档自带的 `latest_evt.csv` 是过时的（GRID-03B
//! 列了 763 天而实际有 789 天，GRID-07 根本没有），不用它。
//!
//! `fits5` 是逐过境的遥测监视（SiPM 电压/温度），`fits3` 是观测计划，
//! `source_orbit` 是 STK 导出的 10 s 星下点表但只覆盖到 2023 年——都不用。
use chrono::prelude::*;
use std::{env, fs, path::PathBuf};

const DIR_PREFIX: &str = "/gecamfs/Exchange/GSDC/missions/GRID";

fn root() -> String {
    env::var("GRID_ARCHIVE").unwrap_or_else(|_| DIR_PREFIX.to_string())
}

fn day_dir(sat_dir: &str, product: &str, day: NaiveDate) -> PathBuf {
    PathBuf::from(root())
        .join(sat_dir)
        .join(product)
        .join(day.format("%Y/%m/%d").to_string())
}

/// 一天的某个产品里版本最高的子目录，没有这一天则 `None`。
pub fn latest_version_dir(sat_dir: &str, product: &str, day: NaiveDate) -> Option<PathBuf> {
    let dir = day_dir(sat_dir, product, day);
    let mut versions: Vec<String> = fs::read_dir(&dir)
        .ok()?
        .filter_map(|entry| entry.ok())
        .filter(|entry| entry.path().is_dir())
        .map(|entry| entry.file_name().to_string_lossy().to_string())
        .collect();
    // 形如 `evt_v03_01`，位宽固定，字典序即版本序
    versions.sort();
    versions.last().map(|v| dir.join(v))
}

/// 目录里的 `.fits`，按文件名（即按过境起始时间）排序。
pub fn fits_files(dir: &PathBuf) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = fs::read_dir(dir)
        .map(|entries| {
            entries
                .filter_map(|entry| entry.ok())
                .map(|entry| entry.path())
                .filter(|path| path.extension().is_some_and(|ext| ext == "fits"))
                .collect()
        })
        .unwrap_or_default();
    files.sort();
    files
}

/// 事例产品与位姿产品的文件名共有的那段 `YYMMDDhhmm_YYMMDDhhmm`，用来配对。
pub fn pass_token(path: &std::path::Path) -> Option<String> {
    let name = path.file_name()?.to_string_lossy();
    // G03_evt_2203110707_2203110721_v03_02.fits / G03_posatt_2203110707_2203110721_v03_02_00.fits
    let mut parts = name.split('_');
    let _code = parts.next()?;
    let _kind = parts.next()?;
    let start = parts.next()?;
    let stop = parts.next()?;
    if start.len() == 10 && stop.len() == 10 {
        Some(format!("{start}_{stop}"))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evt_and_posatt_names_share_the_pass_token() {
        let evt = PathBuf::from("/x/G03_evt_2203110707_2203110721_v03_02.fits");
        let pos = PathBuf::from("/x/G03_posatt_2203110707_2203110721_v03_02_00.fits");
        assert_eq!(pass_token(&evt).as_deref(), Some("2203110707_2203110721"));
        assert_eq!(pass_token(&evt), pass_token(&pos));
        assert_eq!(pass_token(&PathBuf::from("/x/latest_evt.csv")), None);
    }
}
