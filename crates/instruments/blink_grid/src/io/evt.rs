/*
G03_evt_2203110707_2203110721_v03_02.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU       4   ()
  1  EBOUNDS       1 BinTableHDU     ..   97R x 3C    [I, D, D]        Channel, E_MIN, E_MAX (keV)
  2  GTI           1 BinTableHDU     ..   1R x 2C     [D, D]           START, STOP (s)
  3  EVENTS0       1 BinTableHDU     ..   NR x 4C     [D, I, B, B]     TIME, PI, DEAD_TIME, EVT_TYPE
  4  EVENTS1  ... 5 EVENTS2 ... 6 EVENTS3
*/

use blink_core::types::MissionElapsedTime;

use crate::types::{Event, Satellite};

/// 一次过境的事例文件。
pub struct PassFile {
    /// GTI，一行：这次过境的起止（s）
    pub start: f64,
    pub stop: f64,
    /// 各道的下限能量（keV），下标 = `PI - 1`
    ebounds_emin: Vec<f32>,
    detectors: [DetectorHdu; 4],
}

struct DetectorHdu {
    time: Vec<f64>,
    pi: Vec<i16>,
    evt_type: Vec<u8>,
}

impl PassFile {
    pub fn from_fits_file(path: &str) -> Result<Self, fitsio::errors::Error> {
        let mut fptr = fitsio::FitsFile::open(path)?;

        let ebounds = fptr.hdu("EBOUNDS")?;
        let emin = ebounds.read_col::<f64>(&mut fptr, "E_MIN")?;
        let ebounds_emin: Vec<f32> = emin.into_iter().map(|e| e as f32).collect();

        let (start, stop) = Self::gti_of(&mut fptr)?;

        let mut detectors = Vec::with_capacity(4);
        for id in 0..4 {
            let events = fptr.hdu(format!("EVENTS{id}").as_str())?;
            detectors.push(DetectorHdu {
                time: events.read_col::<f64>(&mut fptr, "TIME")?,
                pi: events.read_col::<i16>(&mut fptr, "PI")?,
                evt_type: events.read_col::<u8>(&mut fptr, "EVT_TYPE")?,
            });
        }
        let detectors: [DetectorHdu; 4] = detectors.try_into().ok().expect("four detector HDUs");

        Ok(Self {
            start,
            stop,
            ebounds_emin,
            detectors,
        })
    }

    /// 只读 GTI，用来判断这次过境是否落在要搜的那个小时里，不必把事例读进来。
    pub fn gti_of_file(path: &str) -> Result<(f64, f64), fitsio::errors::Error> {
        let mut fptr = fitsio::FitsFile::open(path)?;
        Self::gti_of(&mut fptr)
    }

    fn gti_of(fptr: &mut fitsio::FitsFile) -> Result<(f64, f64), fitsio::errors::Error> {
        let gti = fptr.hdu("GTI")?;
        let start = gti.read_col::<f64>(fptr, "START")?;
        let stop = gti.read_col::<f64>(fptr, "STOP")?;
        // 实测每个文件恰好一行；多行就取包络，缺行当作空
        let s = start.iter().copied().fold(f64::INFINITY, f64::min);
        let e = stop.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if s.is_finite() && e.is_finite() {
            Ok((s, e))
        } else {
            Ok((0.0, 0.0))
        }
    }

    pub fn len(&self) -> usize {
        self.detectors.iter().map(|d| d.time.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// 各探测器各自的时间回跳数之和（搜索假定输入有序，回跳要在体检里挡）。
    pub fn time_reversals(&self) -> usize {
        self.detectors
            .iter()
            .map(|d| d.time.windows(2).filter(|w| w[1] < w[0]).count())
            .sum()
    }

    /// 四个探测器的事例，未排序。
    pub fn events<S: Satellite>(&self) -> impl Iterator<Item = Event<S>> + '_ {
        let n_channels = self.ebounds_emin.len() as i16;
        self.detectors
            .iter()
            .enumerate()
            .flat_map(move |(id, det)| {
                det.time
                    .iter()
                    .zip(det.pi.iter())
                    .zip(det.evt_type.iter())
                    .map(move |((time, pi), evt_type)| {
                        let index = (*pi as isize) - 1;
                        let energy_kev = if index >= 0 && (index as usize) < self.ebounds_emin.len()
                        {
                            self.ebounds_emin[index as usize]
                        } else {
                            0.0
                        };
                        Event {
                            time: MissionElapsedTime::new(*time),
                            channel: *pi,
                            detector: id as u8,
                            evt_type: *evt_type,
                            energy_kev,
                            // 最高道是溢出道（EVT_TYPE=2 的事例全落在这里）；道号越界也按溢出处理
                            overflow: *pi >= n_channels || *pi < 1,
                        }
                    })
            })
    }
}
