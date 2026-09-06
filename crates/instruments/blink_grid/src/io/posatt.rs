/*
G03_posatt_2203110707_2203110721_v03_02_00.fits
  1  ORBIT_ATTITUDE  NR x 19C   TIME(s), Q1..Q4, wx wy wz (deg/s), ATT_TYPE, X/Y/Z_J2000 (m),
                                X/Y/Z_WGS84 (m), Latitude, Longitude (deg), Altitude (m), POS_TYPE
  1 Hz，与事例文件按过境一一对应。
*/

use blink_core::traits::Interpolatable;
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

    /// 把位置解抹掉（姿态保留）：已知位置错误的日子用，之后由拟合轨道表接手。
    pub fn drop_positions(&mut self) {
        for v in self.latitude.iter_mut().chain(self.longitude.iter_mut()).chain(self.altitude_m.iter_mut()) {
            *v = f32::NAN;
        }
    }
}

/// 允许插值跨越的最大采样间隔（s）。
///
/// 位姿文件的采样不是一律 1 Hz：GRID-03B/04/07 是 1 s，GRID-02 全程 10 s，03B
/// 有些过境也是 10 s（2022-10-01 占 7%）。用 3 s 做门槛时 GRID-02 的候选一个都
/// 定不了位（v5 全量 51 个全进了 `dropped_no_ephemeris`）。
///
/// 姿态解和位置解又都会整段缺失：GRID-02 2020-11-21 的姿态 92% 是 NaN，成段
/// 154–177 s；GRID-03B 2022-10-01 有 6 段各 74 s。去掉 NaN 行以后，缺失段两头
/// 的有效行会被 `Trajectory::interpolate` 当成相邻点连成直线——把几分钟的翻滚
/// 或轨道弧当直线，比丢掉更糟。
///
/// 30 s 把两者分开：10 s 的采样在容忍之内，最短的缺失段（74 s）和过境之间的空
/// 档都在外面。30 s 内把轨道弧当直线的位置误差约 1 km（弓高 Rθ²/8，θ = 30 s ×
/// 1.1 mrad/s），对 800 km 的闪电关联半径无关紧要；姿态在 30 s 内线性插值只
/// 是元数据精度的问题，候选本身不依赖它。
const MAX_SAMPLE_GAP_SECONDS: f64 = 30.0;

/// 在采样密度正常的段内插值；落在缺失段里（两侧有效采样相隔超过
/// [`MAX_SAMPLE_GAP_SECONDS`]）返回 `None`，由上层记账。
pub fn interpolate_sampled<S: Satellite, State: Interpolatable + Clone>(
    trajectory: &Trajectory<MissionElapsedTime<Grid<S>>, State>,
    time: MissionElapsedTime<Grid<S>>,
) -> Option<State> {
    let points = &trajectory.points;
    if points.len() < 2 {
        return None;
    }
    // 与 `Trajectory::interpolate` 取同一段：包住 time 的那段，两端之外夹到最近的一段
    let i = points
        .partition_point(|p| p.timestamp < time)
        .saturating_sub(1)
        .min(points.len() - 2);
    let gap = points[i + 1].timestamp.met() - points[i].timestamp.met();
    if gap > MAX_SAMPLE_GAP_SECONDS {
        return None;
    }
    trajectory.interpolate(time).map(|s| s.state)
}

/// 把一批过境里**有姿态解**的行拼成一条轨迹（按时间排序）；NaN 行跳过。
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
                .filter(|(((_, q1), q2), q3)| q1.is_finite() && q2.is_finite() && q3.is_finite())
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Sat02;

    fn file(rows: &[(f64, f32, f32)]) -> PosAttFile {
        // (t, q, lat)：q 同时填三个四元数分量，lat 也当 lon 用，高度固定
        PosAttFile {
            time: rows.iter().map(|r| r.0).collect(),
            q1: rows.iter().map(|r| r.1).collect(),
            q2: rows.iter().map(|r| r.1).collect(),
            q3: rows.iter().map(|r| r.1).collect(),
            latitude: rows.iter().map(|r| r.2).collect(),
            longitude: rows.iter().map(|r| r.2).collect(),
            altitude_m: rows.iter().map(|_| 5.0e5).collect(),
        }
    }

    #[test]
    fn nan_attitude_rows_are_left_out_and_not_bridged() {
        // 0–9 s 有姿态，10–99 s 全 NaN，100–109 s 又有
        let mut rows: Vec<(f64, f32, f32)> = (0..10).map(|i| (i as f64, 0.1, 1.0)).collect();
        rows.extend((10..100).map(|i| (i as f64, f32::NAN, 1.0)));
        rows.extend((100..110).map(|i| (i as f64, 0.2, 1.0)));
        let traj = attitude_trajectory::<Sat02>(&[file(&rows)]);
        assert_eq!(traj.points.len(), 20);
        assert!(interpolate_sampled(&traj, MissionElapsedTime::new(5.5)).is_some());
        // 缺失段里两侧有效采样相隔 91 s，不能连线
        assert!(interpolate_sampled(&traj, MissionElapsedTime::new(50.0)).is_none());
        assert!(interpolate_sampled(&traj, MissionElapsedTime::new(104.5)).is_some());
    }

    #[test]
    fn a_dropped_second_is_still_interpolated() {
        // 1 Hz 里掉了第 5 秒：相邻有效采样相隔 2 s，在容忍之内
        let rows: Vec<(f64, f32, f32)> = (0..10)
            .filter(|i| *i != 5)
            .map(|i| (i as f64, 0.1, i as f32))
            .collect();
        let traj = position_trajectory::<Sat02>(&[file(&rows)]);
        let p = interpolate_sampled(&traj, MissionElapsedTime::new(5.0)).unwrap();
        assert!((p.latitude - 5.0).abs() < 1e-9);
    }

    #[test]
    fn a_ten_second_cadence_is_interpolated() {
        // GRID-02 的位姿全程 10 s 一采
        let rows: Vec<(f64, f32, f32)> = (0..10).map(|i| (10.0 * i as f64, 0.1, i as f32)).collect();
        let traj = position_trajectory::<Sat02>(&[file(&rows)]);
        let p = interpolate_sampled(&traj, MissionElapsedTime::new(45.0)).unwrap();
        assert!((p.latitude - 4.5).abs() < 1e-9);
    }

    #[test]
    fn a_pass_boundary_is_not_bridged() {
        // 两次过境相隔一小时：过境之间的时刻不能拿两头连线
        let a = file(&(0..5).map(|i| (i as f64, 0.1, 1.0)).collect::<Vec<_>>());
        let b = file(
            &(0..5)
                .map(|i| (3600.0 + i as f64, 0.1, 1.0))
                .collect::<Vec<_>>(),
        );
        let traj = position_trajectory::<Sat02>(&[a, b]);
        assert!(interpolate_sampled(&traj, MissionElapsedTime::new(1800.0)).is_none());
        assert!(interpolate_sampled(&traj, MissionElapsedTime::new(3602.5)).is_some());
    }
}
