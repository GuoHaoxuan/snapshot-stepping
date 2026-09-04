/*
Filename: svom_grm_evt_250101_00_v00.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PrimaryHDU    1 PrimaryHDU      39   ()
  1  EBOUNDS       1 BinTableHDU     66   259R x 3C   [I, E, E]
  2  GTI           1 BinTableHDU     52   1R x 2C   [D, D]
  3  EVENTS01      1 BinTableHDU     68   1094846R x 7C   [D, I, B, E, B, B, B]
  4  EVENTS02      1 BinTableHDU     68   1028001R x 7C   [D, I, B, E, B, B, B]
  5  EVENTS03      1 BinTableHDU     68   1270337R x 7C   [D, I, B, E, B, B, B]
*/

use std::{cmp::Reverse, collections::BinaryHeap};

mod ebounds_hdu;
mod events_hdu;
mod gti_hdu;

// use ebounds_hdu::EboundsHdu;
use events_hdu::EventsHdu;
pub use events_hdu::TimeReversals;
use gti_hdu::GtiHdu;

use crate::{io::evt::events_hdu::EventsHduIterator, types::Event};

pub struct EvtFile {
    // ebounds: EboundsHdu,
    gti: GtiHdu,
    events01: EventsHdu,
    events02: EventsHdu,
    events03: EventsHdu,
    /// 三路合计的时间回跳统计，载入时算一次（体检和诊断共用）。
    time_reversals: TimeReversals,
}

impl EvtFile {
    pub fn from_fits_file(path: &str) -> Result<Self, fitsio::errors::Error> {
        let mut fptr = fitsio::FitsFile::open(path)?;

        // let ebounds = EboundsHdu::from_fptr(&mut fptr)?;
        let gti = GtiHdu::from_fptr(&mut fptr)?;
        let events01 = EventsHdu::from_fptr(&mut fptr, 1)?;
        let events02 = EventsHdu::from_fptr(&mut fptr, 2)?;
        let events03 = EventsHdu::from_fptr(&mut fptr, 3)?;

        let time_reversals = events01
            .time_reversals()
            .merge(events02.time_reversals())
            .merge(events03.time_reversals());

        Ok(Self {
            // ebounds,
            gti,
            events01,
            events02,
            events03,
            time_reversals,
        })
    }
}

impl EvtFile {
    /// 本文件 GTI 与 `[from, to]` 的交集长度（秒）——曝光核算的分子。
    pub fn gti_seconds_within(&self, from: f64, to: f64) -> f64 {
        self.gti.seconds_within(from, to)
    }

    /// `time` 是否在本文件的 GTI 内，见 `GtiHdu::contains`。
    pub fn gti_contains(&self, time: f64) -> bool {
        self.gti.contains(time)
    }

    /// `[from, to]` 是否整个落在同一个 GTI 段内，见 `GtiHdu::covers`。
    pub fn gti_covers(&self, from: f64, to: f64) -> bool {
        self.gti.covers(from, to)
    }

    /// 三路事例表合计的时间回跳统计，见 `TimeReversals`。
    pub fn time_reversals(&self) -> TimeReversals {
        self.time_reversals
    }
}

impl<'a> IntoIterator for &'a EvtFile {
    type Item = Event;
    type IntoIter = EvtFileIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        let mut file_iters = [
            self.events01.into_iter(),
            self.events02.into_iter(),
            self.events03.into_iter(),
        ];
        let mut buffer = BinaryHeap::new();
        for (index, file_iter) in file_iters.iter_mut().enumerate() {
            if let Some(event) = file_iter.next() {
                buffer.push(Reverse((event, index)));
            }
        }
        EvtFileIter { file_iters, buffer }
    }
}

pub struct EvtFileIter<'a> {
    file_iters: [EventsHduIterator<'a>; 3],
    buffer: BinaryHeap<Reverse<(Event, usize)>>,
}

impl Iterator for EvtFileIter<'_> {
    type Item = Event;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(Reverse((event, index))) = self.buffer.pop() {
            if let Some(next_event) = self.file_iters[index].next() {
                self.buffer.push(Reverse((next_event, index)));
            }
            Some(event)
        } else {
            None
        }
    }
}
