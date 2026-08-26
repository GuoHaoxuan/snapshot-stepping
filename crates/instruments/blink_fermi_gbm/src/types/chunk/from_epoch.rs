use blink_core::{error::Error, types::MissionElapsedTime};
use chrono::{TimeDelta, prelude::*};

use super::Chunk;
use crate::io::file::{find_poshist, find_tte};
use crate::io::{PosHistFile, TteFile};
use crate::types::{Detector, FermiGbm};

pub(super) fn from_epoch(epoch: &DateTime<Utc>) -> Result<Chunk, Error> {
    // 没有星历就没法给候选定位，这一小时也就没有意义。
    let poshist_path = find_poshist(epoch)
        .ok_or_else(|| Error::FileNotFound(format!("no GBM poshist for {epoch}")))?;
    let poshist = PosHistFile::from_fits_file(poshist_path.to_str().unwrap())?;

    // 探头的到齐情况随年份变，所以是「找得到几个用几个」而不是要求满编：
    // 2017-10 之前 NaI 没有逐小时数据，2020 年 BGO 只在单独的 BGO/ 目录里。
    let mut tte_files = Vec::new();
    let mut groups = Vec::new();
    for detector in [Detector::Nai, Detector::Bgo] {
        let mut found_any = false;
        for name in detector.names() {
            let Some(path) = find_tte(epoch, name) else {
                continue;
            };
            let file = TteFile::from_fits_file(path.to_str().unwrap(), detector)?;
            if file.is_empty() {
                continue;
            }
            tte_files.push(file);
            found_any = true;
        }
        if found_any {
            groups.push(detector);
        }
    }

    if tte_files.is_empty() {
        return Err(Error::FileNotFound(format!("no GBM TTE for {epoch}")));
    }

    Ok(Chunk {
        span: [
            MissionElapsedTime::<FermiGbm>::from(*epoch),
            MissionElapsedTime::<FermiGbm>::from(*epoch + TimeDelta::hours(1)),
        ],
        tte_files,
        poshist,
        groups,
        dropped_no_ephemeris: Default::default(),
    })
}
