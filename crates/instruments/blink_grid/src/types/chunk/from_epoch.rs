use blink_core::{error::Error, types::MissionElapsedTime};
use chrono::{TimeDelta, prelude::*};
use std::marker::PhantomData;
use std::path::PathBuf;

use super::Chunk;
use crate::io::file::{fits_files, latest_version_dir};
use crate::io::{PassFile, PosAttFile};
use crate::types::instrument::{Grid, Satellite};

/// 本小时要看的两个日目录：当天，加上（0 时）前一天——过境文件归到起始那天，
/// 跨午夜的那一次会落在前一天的目录里。
fn day_dirs<S: Satellite>(epoch: &DateTime<Utc>, product: &str) -> Vec<PathBuf> {
    let mut days = vec![epoch.date_naive()];
    if epoch.hour() == 0 {
        days.push(epoch.date_naive() - TimeDelta::days(1));
    }
    days.into_iter()
        .filter_map(|day| latest_version_dir(S::DIR, product, day))
        .collect()
}

pub(super) fn from_epoch<S: Satellite>(epoch: &DateTime<Utc>) -> Result<Chunk<S>, Error> {
    let span = [
        MissionElapsedTime::<Grid<S>>::from(*epoch),
        MissionElapsedTime::<Grid<S>>::from(*epoch + TimeDelta::hours(1)),
    ];
    let (start, stop) = (span[0].met(), span[1].met());

    let evt_dirs = day_dirs::<S>(epoch, "fits7");
    if evt_dirs.is_empty() {
        return Err(Error::FileNotFound(format!(
            "{}: no event product for {}",
            S::DIR,
            epoch.date_naive()
        )));
    }

    // 先只读 GTI 挑出与本小时有交集的过境，再把它们的事例读进来
    let mut passes = Vec::new();
    let mut gti = Vec::new();
    for dir in &evt_dirs {
        for path in fits_files(dir) {
            let p = path.to_str().unwrap();
            let (s, e) = PassFile::gti_of_file(p)?;
            if e <= start || s >= stop {
                continue;
            }
            gti.push([s.max(start), e.min(stop)]);
            passes.push(PassFile::from_fits_file(p)?);
        }
    }
    if passes.is_empty() {
        return Err(Error::FileNotFound(format!(
            "{}: no pass overlaps {}",
            S::DIR,
            epoch.format("%Y-%m-%dT%H")
        )));
    }
    gti.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap());

    // 位姿：整天的文件都读（一天 1–21 个小文件），省得逐过境配对。读不出来的
    // 跳过并计数：位姿是辅助数据，归档里有 0 字节的 posatt 文件（GRID-03B
    // 2023-11-20 一个），不能让它把整天的事例都拖成 corrupt_data。
    let mut posatt = Vec::new();
    let mut posatt_unreadable = 0usize;
    for dir in day_dirs::<S>(epoch, "fits8") {
        for path in fits_files(&dir) {
            match PosAttFile::from_fits_file(path.to_str().unwrap()) {
                Ok(file) => posatt.push(file),
                Err(_) => posatt_unreadable += 1,
            }
        }
    }

    Ok(Chunk {
        span,
        passes,
        posatt,
        gti,
        dropped_no_ephemeris: Default::default(),
        events_outside_gti: Default::default(),
        dropped_dead_gap: Default::default(),
        dropped_high_rate: Default::default(),
        dropped_simultaneous: Default::default(),
        without_attitude: Default::default(),
        dropped_single_detector: Default::default(),
        posatt_unreadable,
        _satellite: PhantomData,
    })
}

pub(super) fn last_modified<S: Satellite>(epoch: &DateTime<Utc>) -> Result<DateTime<Utc>, Error> {
    let mut latest: Option<DateTime<Utc>> = None;
    for product in ["fits7", "fits8"] {
        for dir in day_dirs::<S>(epoch, product) {
            for path in fits_files(&dir) {
                let modified: DateTime<Utc> = std::fs::metadata(&path)?.modified()?.into();
                latest = Some(latest.map_or(modified, |l| l.max(modified)));
            }
        }
    }
    latest.ok_or_else(|| {
        Error::FileNotFound(format!("{}: no files for {}", S::DIR, epoch.date_naive()))
    })
}
