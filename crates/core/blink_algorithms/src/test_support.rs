//! 测试用的最小仪器：给核心算法的单元测试当事例源，各模块共用。

use blink_core::{
    error::Error,
    traits::{Chunk, Event, Instrument},
    types::{Coverage, MissionElapsedTime, Signal},
};
use chrono::{DateTime, NaiveDate, TimeZone, Utc};
use serde::Serialize;
use std::sync::OnceLock;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub(crate) struct TestInstrument;

pub(crate) struct TestChunk;

impl Chunk for TestChunk {
    type Event = TestEvent;

    fn from_epoch(_: &DateTime<Utc>) -> Result<Self, Error> {
        Err(Error::Unknown)
    }
    fn search(&self) -> Vec<Signal<Self::Event>> {
        Vec::new()
    }
    fn last_modified(_: &DateTime<Utc>) -> Result<DateTime<Utc>, Error> {
        Err(Error::Unknown)
    }
    fn coverage(&self) -> Coverage {
        Coverage {
            span_seconds: 0.0,
            masked_seconds: 0.0,
        }
    }
}

impl Instrument for TestInstrument {
    type Chunk = TestChunk;

    fn ref_time() -> &'static DateTime<Utc> {
        static REF_TIME: OnceLock<DateTime<Utc>> = OnceLock::new();
        REF_TIME.get_or_init(|| Utc.with_ymd_and_hms(2012, 1, 1, 0, 0, 0).unwrap())
    }
    fn launch_day() -> NaiveDate {
        NaiveDate::from_ymd_opt(2017, 6, 15).unwrap()
    }
    fn name() -> &'static str {
        "test"
    }
}

#[derive(Serialize, Debug, Clone)]
pub(crate) struct TestEvent {
    pub(crate) seconds: f64,
    pub(crate) group: u8,
}

impl Event for TestEvent {
    type Instrument = TestInstrument;
    type ChannelType = u16;

    fn time(&self) -> MissionElapsedTime<TestInstrument> {
        MissionElapsedTime::new(self.seconds)
    }
    fn channel(&self) -> u16 {
        100
    }
    fn group(&self) -> u8 {
        self.group
    }
    fn keep(&self) -> bool {
        true
    }
}
