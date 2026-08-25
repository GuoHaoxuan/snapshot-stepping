//! 1B 归档的路径与文件名解析。
//!
//! 归档按 APID 分子目录：`$HXMT_1B_DIR/YYYY/YYYYMMDD/<APID>/`。本模块只解
//! 析下面用到的两类；同级还有几个目录，找遥测时会碰到：
//!
//! * `0531` = HE_ANBL（模拟板 HK，`TMYxxx` 原始遥测量）
//! * `0549` = HE_HV（高压）
//! * `0162` = GPS_Status
//!
//! `sci` 那三个（`0642/0922/1686` = HE_Evt_Src）装的是 882 B 的原始 CCSDS
//! 包，要自己解；`eng` 那三个（`0766/1009/1781`）含 `sTime_Last_Bdc`。

use chrono::prelude::*;
use std::path::Path;
use std::{env, sync::LazyLock};

pub static HXMT_1B_DIR: LazyLock<String> = LazyLock::new(|| {
    env::var("HXMT_1B_DIR").unwrap_or_else(|_| "/hxmtfs/data/Archive_tmp/1B".to_string())
});

pub fn find_filename(type_: &str, time: DateTime<Utc>, serial_num: &str) -> Option<String> {
    let code = match (type_, serial_num) {
        ("eng", "A") => "0766",
        ("eng", "B") => "1009",
        ("eng", "C") => "1781",
        ("sci", "A") => "0642",
        ("sci", "B") => "0922",
        ("sci", "C") => "1686",
        _ => panic!("Invalid type or serial number"),
    };

    let mut path = None;
    let mut version = -1;

    let folder_path = Path::new(HXMT_1B_DIR.as_str())
        .join(format!("{}", time.year()))
        .join(format!(
            "{}{:02}{:02}",
            time.year(),
            time.month(),
            time.day()
        ))
        .join(code);

    if !folder_path.exists() {
        return None;
    }

    let prefix = format!(
        "HXMT_1B_{}_{}{:02}{:02}T{:02}",
        code,
        time.year(),
        time.month(),
        time.day(),
        time.hour()
    );

    // 读取目录内容
    if let Ok(entries) = std::fs::read_dir(folder_path) {
        for entry in entries.flatten() {
            let entry_path = entry.path();
            if entry_path.is_file() {
                let filename = entry_path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or("");
                if filename.len() >= 40
                    && filename.starts_with(&prefix)
                    && filename
                        .chars()
                        .nth(39)
                        .map(|c| c.is_ascii_digit())
                        .unwrap_or(false)
                {
                    let ver = filename
                        .chars()
                        .nth(39)
                        .and_then(|c| c.to_digit(10))
                        .map(|d| d as i32)
                        .unwrap_or(-1);

                    if ver > version {
                        version = ver;
                        path = Some(entry_path.to_string_lossy().into_owned());
                    }
                }
            }
        }
    }

    path
}

/// 获取指定小时的科学数据文件路径（仅返回存在的 box）。
pub fn get_sci_filenames(time: DateTime<Utc>) -> Vec<(String, String)> {
    get_filenames("sci", time)
}

/// 获取指定小时的工程数据文件路径（仅返回存在的 box）。
pub fn get_eng_filenames(time: DateTime<Utc>) -> Vec<(String, String)> {
    get_filenames("eng", time)
}

fn get_filenames(type_: &str, time: DateTime<Utc>) -> Vec<(String, String)> {
    let serial_nums = ["A", "B", "C"];
    let mut result = Vec::new();
    for &sn in &serial_nums {
        if let Some(path) = find_filename(type_, time, sn) {
            result.push((sn.to_string(), path));
        }
    }
    result
}
