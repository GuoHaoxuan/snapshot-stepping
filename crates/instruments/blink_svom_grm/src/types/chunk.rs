use blink_core::types::{Coverage, MissionElapsedTime};

use crate::io::file::{find_att_by_time, find_evt_by_time, find_orb_by_time};
use crate::io::{AttFile, EvtFile, OrbFile};
use crate::types::event::Event;
use crate::types::instrument::SvomGrm;
use blink_core::error::Error;
use chrono::prelude::*;

mod from_epoch;
mod search;

pub struct Chunk {
    pub span: [MissionElapsedTime<SvomGrm>; 2],
    pub att_file: AttFile,
    pub evt_file: EvtFile,
    pub orb_file: OrbFile,
}

impl blink_core::traits::Chunk for Chunk {
    type Event = Event;

    fn from_epoch(epoch: &chrono::DateTime<chrono::Utc>) -> Result<Self, blink_core::error::Error>
    where
        Self: Sized,
    {
        from_epoch::from_epoch(epoch)
    }

    fn search(&self) -> Vec<blink_core::types::Signal<Self::Event>> {
        search::search(self)
    }

    fn coverage(&self) -> Coverage {
        Coverage {
            span_seconds: self.span[1].met() - self.span[0].met(),
            // GRM 的搜索目前不做任何时间屏蔽
            masked_seconds: 0.0,
        }
    }

    fn last_modified(epoch: &DateTime<Utc>) -> Result<DateTime<Utc>, Error> {
        let filenames = [
            find_att_by_time(epoch),
            find_evt_by_time(epoch),
            find_orb_by_time(epoch),
        ];

        let last_modifieds: Vec<DateTime<Utc>> = filenames
            .iter()
            .flatten()
            .map(|filename| {
                let last_modified = std::fs::metadata(filename)?.modified()?;
                let datetime: DateTime<Utc> = last_modified.into();
                Ok::<DateTime<Utc>, Error>(datetime)
            })
            .collect::<Result<Vec<DateTime<Utc>>, Error>>()?;

        let max_last_modified = last_modifieds
            .into_iter()
            .max()
            .ok_or_else(|| Error::FileNotFound("No files found".to_string()))?;

        Ok(max_last_modified)
    }
}
