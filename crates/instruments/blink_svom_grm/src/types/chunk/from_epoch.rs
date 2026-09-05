use crate::{
    io::{
        AttFile, EvtFile, OrbFile,
        file::{find_att_by_time, find_evt_by_time, find_orb_by_time},
    },
    types::instrument::SvomGrm,
};

use super::Chunk;
use blink_core::{error::Error, types::MissionElapsedTime};
use chrono::{TimeDelta, prelude::*};

pub(super) fn from_epoch(epoch: &DateTime<Utc>) -> Result<Chunk, Error> {
    let att_filename = find_att_by_time(epoch)?;
    let att_file = AttFile::from_fits_file(att_filename.to_str().unwrap())?;
    let evt_filename = find_evt_by_time(epoch)?;
    let evt_file = EvtFile::from_fits_file(evt_filename.to_str().unwrap())?;
    let orb_filename = find_orb_by_time(epoch)?;
    let orb_file = OrbFile::from_fits_file(orb_filename.to_str().unwrap())?;
    Ok(Chunk {
        span: [
            MissionElapsedTime::<SvomGrm>::from(*epoch),
            MissionElapsedTime::<SvomGrm>::from(*epoch + TimeDelta::hours(1)),
        ],
        att_file,
        evt_file,
        orb_file,
        dropped_no_ephemeris: Default::default(),
        dropped_single_detector: Default::default(),
        events_outside_gti: Default::default(),
    })
}
