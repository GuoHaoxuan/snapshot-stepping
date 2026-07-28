use super::super::path::get_path;
use crate::types::{Detector, Event, Scintillator};
use blink_core::{error::Error, types::MissionElapsedTime};
use chrono::prelude::*;

pub struct EventFile {
    // HDU 1: Events
    time: Vec<f64>,
    det_id: Vec<u8>,
    channel: Vec<u8>,
    pulse_width: Vec<u8>,
    acd: Vec<[bool; 18]>,
    event_type: Vec<u8>,
    // flag: Vec<u8>,
}

impl EventFile {
    fn get_path(epoch: &DateTime<Utc>) -> Result<String, Error> {
        get_path(epoch, "Evt")
    }

    pub fn last_modified(epoch: &DateTime<Utc>) -> Result<DateTime<Utc>, Error> {
        let path = Self::get_path(epoch)?;
        let metadata = std::fs::metadata(path)?;
        let modified_time = metadata.modified()?;
        let datetime: DateTime<Utc> = modified_time.into();
        Ok(datetime)
    }

    pub fn from_epoch(epoch: &DateTime<Utc>) -> Result<Self, Error> {
        let path = Self::get_path(epoch)?;
        Self::new(&path)
    }

    pub fn times(&self) -> &[f64] {
        &self.time
    }

    pub fn det_ids(&self) -> &[u8] {
        &self.det_id
    }

    pub fn channels(&self) -> &[u8] {
        &self.channel
    }

    pub fn pulse_widths(&self) -> &[u8] {
        &self.pulse_width
    }

    pub fn event_types(&self) -> &[u8] {
        &self.event_type
    }

    fn new(filename: &str) -> Result<Self, Error> {
        let mut fptr = fitsio::FitsFile::open(filename)?;

        // HDU 1: Events
        let events = fptr.hdu("Events")?;
        let time = events.read_col::<f64>(&mut fptr, "Time")?;
        let det_id = events.read_col::<u8>(&mut fptr, "Det_ID")?;
        let channel = events.read_col::<u8>(&mut fptr, "Channel")?;
        let pulse_width = events.read_col::<u8>(&mut fptr, "Pulse_Width")?;

        let acd_raw = events.read_col::<u32>(&mut fptr, "ACD")?;
        let mut acd = Vec::with_capacity(acd_raw.len());
        for &value in &acd_raw {
            let mut array = [false; 18];
            array.iter_mut().enumerate().for_each(|(i, bit)| {
                *bit = ((value >> i) & 1) == 1;
            });
            acd.push(array);
        }

        let event_type = events.read_col::<u8>(&mut fptr, "Event_Type")?;
        // let flag = events.read_col::<u8>(&mut fptr, "FLAG")?;

        Ok(Self {
            time,
            det_id,
            channel,
            pulse_width,
            acd,
            event_type,
            // flag,
        })
    }
}

impl<'a> IntoIterator for &'a EventFile {
    type Item = Event;
    type IntoIter = Iter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        Iter {
            event_file: self,
            index: 0,
        }
    }
}

pub struct Iter<'a> {
    event_file: &'a EventFile,
    index: usize,
}

impl Iterator for Iter<'_> {
    type Item = Event;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.event_file.time.len() {
            let event = Event::new(
                MissionElapsedTime::new(self.event_file.time[self.index]),
                self.event_file.channel[self.index],
                Detector {
                    id: self.event_file.det_id[self.index],
                    scintillator: if self.event_file.pulse_width[self.index] < 75 {
                        Scintillator::Nai
                    } else {
                        Scintillator::Csi
                    },
                },
                self.event_file.event_type[self.index] == 1,
                self.event_file.acd[self.index],
            );
            self.index += 1;
            Some(event)
        } else {
            None
        }
    }
}
