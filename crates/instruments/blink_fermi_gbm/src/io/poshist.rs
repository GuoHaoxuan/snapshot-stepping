/*
glg_poshist_all_190101_v00.fit.gz
  1  GLAST POS HIST  BinTableHDU  86510R x 19C
  cols: SCLK_UTC, QSJ_1..4, WSJ_1..3, POS_X/Y/Z, VEL_X/Y/Z, SC_LAT, SC_LON, SADA_PY, SADA_NY, FLAGS
*/

use blink_core::types::{Attitude, MissionElapsedTime, Position, TemporalState, Trajectory};

use crate::types::FermiGbm;

/// 地球平均半径（m）。poshist 只给 ECI 位置向量，没有高度列，所以高度按
/// 球近似由 |POS| 减去它得到——量级用（实测 |POS| 中位 6908 km → 537 km），
/// 判据不依赖它，逐事例定位用的是 SC_LAT/SC_LON。
const EARTH_RADIUS_M: f64 = 6_371_000.0;

/// FLAGS 里标记「正处在南大西洋异常区」的位。
///
/// 实测：bit0 恒为 1（不是 SAA），bit1 置位占一天的 13.2%，其覆盖范围是
/// 经度 −94.5°..24.0°、纬度 −25.7°..1.9°，正是 SAA。
const FLAG_BIT_SAA: u16 = 1;

pub struct PosHistFile {
    time: Vec<f64>,
    q1: Vec<f64>,
    q2: Vec<f64>,
    q3: Vec<f64>,
    latitude: Vec<f32>,
    longitude: Vec<f32>,
    altitude: Vec<f64>,
    flags: Vec<i16>,
}

impl PosHistFile {
    pub fn from_fits_file(path: &str) -> Result<Self, fitsio::errors::Error> {
        let mut fptr = fitsio::FitsFile::open(path)?;
        let hdu = fptr.hdu(1)?;

        let time = hdu.read_col::<f64>(&mut fptr, "SCLK_UTC")?;
        let q1 = hdu.read_col::<f64>(&mut fptr, "QSJ_1")?;
        let q2 = hdu.read_col::<f64>(&mut fptr, "QSJ_2")?;
        let q3 = hdu.read_col::<f64>(&mut fptr, "QSJ_3")?;
        let latitude = hdu.read_col::<f32>(&mut fptr, "SC_LAT")?;
        let longitude = hdu.read_col::<f32>(&mut fptr, "SC_LON")?;
        let flags = hdu.read_col::<i16>(&mut fptr, "FLAGS")?;

        let x = hdu.read_col::<f32>(&mut fptr, "POS_X")?;
        let y = hdu.read_col::<f32>(&mut fptr, "POS_Y")?;
        let z = hdu.read_col::<f32>(&mut fptr, "POS_Z")?;
        let altitude = x
            .iter()
            .zip(y.iter())
            .zip(z.iter())
            .map(|((x, y), z)| {
                let (x, y, z) = (*x as f64, *y as f64, *z as f64);
                (x * x + y * y + z * z).sqrt() - EARTH_RADIUS_M
            })
            .collect();

        Ok(Self {
            time,
            q1,
            q2,
            q3,
            latitude,
            longitude,
            altitude,
            flags,
        })
    }

    /// `[from, to]` 内被标记为 SAA 的秒数。poshist 是 1 Hz，一个样本记一秒。
    pub fn saa_seconds_within(&self, from: f64, to: f64) -> f64 {
        self.time
            .iter()
            .zip(self.flags.iter())
            .filter(|(time, _)| **time >= from && **time < to)
            .filter(|(_, flags)| (*flags >> FLAG_BIT_SAA) & 1 == 1)
            .count() as f64
    }
}

/// SC_LON 在文件里是 0..360，全流程统一用 −180..180。
fn wrap_longitude(longitude: f64) -> f64 {
    if longitude > 180.0 {
        longitude - 360.0
    } else {
        longitude
    }
}

impl From<&PosHistFile> for Trajectory<MissionElapsedTime<FermiGbm>, Position> {
    fn from(file: &PosHistFile) -> Self {
        let points = file
            .time
            .iter()
            .zip(file.longitude.iter())
            .zip(file.latitude.iter())
            .zip(file.altitude.iter())
            .map(|(((time, longitude), latitude), altitude)| TemporalState {
                timestamp: MissionElapsedTime::new(*time),
                state: Position {
                    longitude: wrap_longitude(*longitude as f64),
                    latitude: *latitude as f64,
                    altitude: uom::si::f64::Length::new::<uom::si::length::meter>(*altitude),
                },
            })
            .collect();

        Trajectory { points }
    }
}

impl From<&PosHistFile> for Trajectory<MissionElapsedTime<FermiGbm>, Attitude> {
    fn from(file: &PosHistFile) -> Self {
        let points = file
            .time
            .iter()
            .zip(file.q1.iter())
            .zip(file.q2.iter())
            .zip(file.q3.iter())
            .map(|(((time, q1), q2), q3)| TemporalState {
                timestamp: MissionElapsedTime::new(*time),
                state: Attitude {
                    q1: *q1,
                    q2: *q2,
                    q3: *q3,
                },
            })
            .collect();

        Trajectory { points }
    }
}
