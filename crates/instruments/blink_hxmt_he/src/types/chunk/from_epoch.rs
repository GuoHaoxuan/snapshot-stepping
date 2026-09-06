use crate::{
    io::{
        level_1b::{SciFile, get_eng_filenames, get_sci_filenames, read_stime_offset},
        level_1k::{AttFile, EventFile, OrbitFile},
    },
    types::HxmtHe,
};

use super::Chunk;
use blink_core::{error::Error, types::MissionElapsedTime};
use chrono::{TimeDelta, prelude::*};
use std::sync::OnceLock;
use std::sync::atomic::AtomicUsize;

pub(super) fn from_epoch(epoch: &DateTime<Utc>) -> Result<Chunk, Error> {
    let event_file = EventFile::from_epoch(epoch)?;
    let orbit_file = OrbitFile::from_epoch(epoch)?;
    let att_file = AttFile::from_epoch(epoch)?;

    let sci_pairs = get_sci_filenames(*epoch);
    let eng_pairs = get_eng_filenames(*epoch);

    let mut sci_files = Vec::new();
    let mut stime_offsets = Vec::new();

    for (box_name, sci_path) in &sci_pairs {
        let sci = SciFile::new(sci_path)?;
        // 找对应的 eng 文件。offset 拿不到就必须失败：曾经的 `.ok().unwrap_or(0.0)`
        // 会让该箱锚点全体平移约 3.9 亿秒，饱和掩模整体落空，搜索照常进行——
        // 静默污染。缺文件按 MissingData 入账，读出无稳定众数按 UtcFreeze 入账。
        let offset = match eng_pairs.iter().find(|(bn, _)| bn == box_name) {
            Some((_, eng_path)) => read_stime_offset(eng_path)?,
            None => {
                return Err(Error::FileNotFound(format!(
                    "HE_Eng for box {box_name} at {epoch}"
                )));
            }
        };
        sci_files.push((box_name.clone(), sci));
        stime_offsets.push((box_name.clone(), offset));
    }

    Ok(Chunk {
        event_file,
        sci_files,
        stime_offsets,
        orbit_file,
        att_file,
        span: [
            MissionElapsedTime::<HxmtHe>::from(*epoch),
            MissionElapsedTime::<HxmtHe>::from(*epoch + TimeDelta::hours(1)),
        ],
        saturation_cache: OnceLock::new(),
        dropped_no_ephemeris: AtomicUsize::new(0),
        without_attitude: Default::default(),
        dropped_single_detector: Default::default(),
        duplicate_cache: OnceLock::new(),
    })
}
