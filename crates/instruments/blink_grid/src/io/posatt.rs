/*
G03_posatt_2203110707_2203110721_v03_02_00.fits
  1  ORBIT_ATTITUDE  NR x 19C   TIME(s), Q1..Q4, wx wy wz (deg/s), ATT_TYPE, X/Y/Z_J2000 (m),
                                X/Y/Z_WGS84 (m), Latitude, Longitude (deg), Altitude (m), POS_TYPE
  1 Hz，与事例文件按过境一一对应。
*/

use blink_core::types::{Attitude, MissionElapsedTime, Position, TemporalState, Trajectory};
use uom::si::f64::*;
use uom::si::length::meter;

use crate::types::{Grid, Satellite};

/// 一次过境的位姿文件。
///
/// 位置并不总有：实测 2024-02 之后（GRID-07 几乎全程、03B/04 的 2024 年）
/// `POS_TYPE = 0`，经纬度、WGS84/J2000 坐标全是 NaN，只剩四元数。这样的
/// 过境里搜出的候选定不了位，只能丢弃并计数（`dropped_no_ephemeris`）。
pub struct PosAttFile {
    time: Vec<f64>,
    q1: Vec<f32>,
    q2: Vec<f32>,
    q3: Vec<f32>,
    latitude: Vec<f32>,
    longitude: Vec<f32>,
    altitude_m: Vec<f32>,
}

impl PosAttFile {
    pub fn from_fits_file(path: &str) -> Result<Self, fitsio::errors::Error> {
        let mut fptr = fitsio::FitsFile::open(path)?;
        let hdu = fptr.hdu("ORBIT_ATTITUDE")?;
        Ok(Self {
            time: hdu.read_col::<f64>(&mut fptr, "TIME")?,
            q1: hdu.read_col::<f32>(&mut fptr, "Q1")?,
            q2: hdu.read_col::<f32>(&mut fptr, "Q2")?,
            q3: hdu.read_col::<f32>(&mut fptr, "Q3")?,
            latitude: hdu.read_col::<f32>(&mut fptr, "Latitude")?,
            longitude: hdu.read_col::<f32>(&mut fptr, "Longitude")?,
            altitude_m: hdu.read_col::<f32>(&mut fptr, "Altitude")?,
        })
    }

    pub fn len(&self) -> usize {
        self.time.len()
    }

    pub fn is_empty(&self) -> bool {
        self.time.is_empty()
    }

    /// 有位置解的行数
    pub fn positioned_rows(&self) -> usize {
        self.latitude.iter().filter(|v| v.is_finite()).count()
    }
}

/// 把一批过境的姿态拼成一条轨迹（按时间排序）。
pub fn attitude_trajectory<S: Satellite>(
    files: &[PosAttFile],
) -> Trajectory<MissionElapsedTime<Grid<S>>, Attitude> {
    let mut points: Vec<TemporalState<MissionElapsedTime<Grid<S>>, Attitude>> = files
        .iter()
        .flat_map(|f| {
            f.time
                .iter()
                .zip(f.q1.iter())
                .zip(f.q2.iter())
                .zip(f.q3.iter())
                .map(|(((t, q1), q2), q3)| TemporalState {
                    timestamp: MissionElapsedTime::new(*t),
                    state: Attitude {
                        q1: *q1 as f64,
                        q2: *q2 as f64,
                        q3: *q3 as f64,
                    },
                })
        })
        .collect();
    points.sort_by(|a, b| a.timestamp.cmp(&b.timestamp));
    Trajectory { points }
}

/// 把一批过境里**有位置解**的行拼成一条轨迹；没有解的行跳过（NaN 不能插值）。
pub fn position_trajectory<S: Satellite>(
    files: &[PosAttFile],
) -> Trajectory<MissionElapsedTime<Grid<S>>, Position> {
    let mut points: Vec<TemporalState<MissionElapsedTime<Grid<S>>, Position>> = files
        .iter()
        .flat_map(|f| {
            f.time
                .iter()
                .zip(f.latitude.iter())
                .zip(f.longitude.iter())
                .zip(f.altitude_m.iter())
                .filter(|(((_, lat), lon), alt)| {
                    lat.is_finite() && lon.is_finite() && alt.is_finite()
                })
                .map(|(((t, lat), lon), alt)| TemporalState {
                    timestamp: MissionElapsedTime::new(*t),
                    state: Position {
                        longitude: *lon as f64,
                        latitude: *lat as f64,
                        altitude: Length::new::<meter>(*alt as f64),
                    },
                })
        })
        .collect();
    points.sort_by(|a, b| a.timestamp.cmp(&b.timestamp));
    Trajectory { points }
}
