use blink_core::types::MissionElapsedTime;

use crate::types::Event;

pub(super) struct EventsHdu {
    id: u8,
    time: Vec<f64>,
    pi: Vec<i16>,
    gain_type: Vec<u8>,
    dead_time: Vec<f32>,
    evt_type: Vec<u8>,
    anti_coin: Vec<u8>,
    flag: Vec<u8>,
}

impl EventsHdu {
    pub fn from_fptr(fptr: &mut fitsio::FitsFile, id: u8) -> Result<Self, fitsio::errors::Error> {
        let events = fptr.hdu(format!("EVENTS0{}", id).as_str())?;

        let time = events.read_col::<f64>(fptr, "TIME")?;
        let pi = events.read_col::<i16>(fptr, "PI")?;
        let gain_type = events.read_col::<u8>(fptr, "GAIN_TYPE")?;
        let dead_time = events.read_col::<f32>(fptr, "DEAD_TIME")?;
        let evt_type = events.read_col::<u8>(fptr, "EVT_TYPE")?;
        let anti_coin = events.read_col::<u8>(fptr, "ANTI_COIN")?;
        let flag = events.read_col::<u8>(fptr, "FLAG")?;

        Ok(Self {
            id,
            time,
            pi,
            gain_type,
            dead_time,
            evt_type,
            anti_coin,
            flag,
        })
    }
}

impl EventsHdu {
    /// 本路事例表里的时间回跳。三路各自有序是 k 路归并的前提，也是搜索
    /// 窗长判据的前提。
    pub(super) fn time_reversals(&self) -> TimeReversals {
        self.time
            .windows(2)
            .filter(|pair| pair[1] < pair[0])
            .fold(TimeReversals::default(), |acc, pair| TimeReversals {
                count: acc.count + 1,
                max_magnitude: acc.max_magnitude.max(pair[0] - pair[1]),
            })
    }
}

/// 事例表内时间回跳的统计量。
///
/// 实测幅度呈双峰，中间是空的：一半落在 1–5 个时间量化步（2⁻²⁰ s ≈ 0.95 µs）
/// 之内，另一半在 0.49–1.0 s，两者之间一处也没有。前者是同一微秒尺度内的
/// 排序抖动，排序即可；后者是整段数据被记录了两次，必须整小时排除。
#[derive(Default, Clone, Copy)]
pub struct TimeReversals {
    /// 回跳处数（含无害的微回跳）。
    pub count: usize,
    /// 最大回跳幅度（秒）。判据看它，不看处数。
    pub max_magnitude: f64,
}

impl TimeReversals {
    pub fn merge(self, other: Self) -> Self {
        Self {
            count: self.count + other.count,
            max_magnitude: self.max_magnitude.max(other.max_magnitude),
        }
    }
}

impl<'a> IntoIterator for &'a EventsHdu {
    type Item = Event;
    type IntoIter = EventsHduIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        EventsHduIterator {
            hdu: self,
            index: 0,
        }
    }
}

pub struct EventsHduIterator<'a> {
    hdu: &'a EventsHdu,
    index: usize,
}

impl Iterator for EventsHduIterator<'_> {
    type Item = Event;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.hdu.time.len() {
            let event = Event {
                time: MissionElapsedTime::new(self.hdu.time[self.index]),
                channel: self.hdu.pi[self.index],
                detector_id: self.hdu.id,
                gain_type: self.hdu.gain_type[self.index],
                dead_time: self.hdu.dead_time[self.index],
                evt_type: self.hdu.evt_type[self.index],
                anti_coin: self.hdu.anti_coin[self.index],
                flag: self.hdu.flag[self.index],
            };
            self.index += 1;
            Some(event)
        } else {
            None
        }
    }
}
