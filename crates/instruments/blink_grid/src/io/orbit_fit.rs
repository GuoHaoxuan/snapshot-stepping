//! 拟合出来的轨道表：位姿文件没有位置解时的后备。
//!
//! 2024-02-09 起三颗星的位姿文件只剩姿态（POS_TYPE=0），归档里也没有轨道根数。
//! `scripts/cluster/grid_orbit_fit.py` 用太阳同步圆轨道模型（轨道面从有位置解的
//! 日子定标）加本底计数率随磁纬的变化拟合每天的沿轨相位，写出每次过境内 10 s
//! 一行的 `time,lon,lat,alt_m`（time 为 MET 秒）。盲拟合有真值的日子，位置误差
//! 中位 73–174 km、最大 380 km，对 800 km 的闪电关联半径够用。
//!
//! 目录由 `GRID_ORBIT_FIT_DIR` 指定，按 `<星目录名>/<YYYYMMDD>.csv` 存放；没设
//! 或没有文件就没有后备，候选照旧计入 `dropped_no_ephemeris`。

use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use blink_core::types::{MissionElapsedTime, Position, TemporalState, Trajectory};
use chrono::NaiveDate;
use uom::si::f64::*;
use uom::si::length::meter;

use crate::types::{Grid, Satellite};

pub struct OrbitFitFile {
    time: Vec<f64>,
    lon: Vec<f64>,
    lat: Vec<f64>,
    alt_m: Vec<f64>,
}

impl OrbitFitFile {
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self, String> {
        let mut file = Self {
            time: Vec::new(),
            lon: Vec::new(),
            lat: Vec::new(),
            alt_m: Vec::new(),
        };
        for (n, line) in reader.lines().enumerate() {
            let line = line.map_err(|e| e.to_string())?;
            if n == 0 || line.trim().is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split(',').collect();
            if cols.len() < 4 {
                return Err(format!("line {}: expected time,lon,lat,alt_m", n + 1));
            }
            let parse = |s: &str| s.trim().parse::<f64>().map_err(|e| format!("line {}: {e}", n + 1));
            file.time.push(parse(cols[0])?);
            file.lon.push(parse(cols[1])?);
            file.lat.push(parse(cols[2])?);
            file.alt_m.push(parse(cols[3])?);
        }
        Ok(file)
    }

    pub fn from_path(path: &PathBuf) -> Result<Self, String> {
        let f = std::fs::File::open(path).map_err(|e| e.to_string())?;
        Self::from_reader(BufReader::new(f))
    }

    pub fn len(&self) -> usize {
        self.time.len()
    }

    pub fn is_empty(&self) -> bool {
        self.time.is_empty()
    }
}

/// 拟合轨道表的根目录；没设环境变量就是没有后备。
pub fn orbit_fit_dir() -> Option<PathBuf> {
    std::env::var("GRID_ORBIT_FIT_DIR").ok().map(PathBuf::from)
}

/// 某颗星某一天的拟合轨道表；目录没设或文件不存在返回 None。
pub fn load_day(sat_dir: &str, day: NaiveDate) -> Option<Result<OrbitFitFile, String>> {
    let path = orbit_fit_dir()?
        .join(sat_dir)
        .join(format!("{}.csv", day.format("%Y%m%d")));
    if !path.exists() {
        return None;
    }
    Some(OrbitFitFile::from_path(&path))
}

/// 把几天的拟合轨道表拼成一条位置轨迹（按时间排序）。
pub fn trajectory<S: Satellite>(
    files: &[OrbitFitFile],
) -> Trajectory<MissionElapsedTime<Grid<S>>, Position> {
    let mut points: Vec<TemporalState<MissionElapsedTime<Grid<S>>, Position>> = files
        .iter()
        .flat_map(|f| {
            f.time
                .iter()
                .zip(f.lon.iter())
                .zip(f.lat.iter())
                .zip(f.alt_m.iter())
                .map(|(((t, lon), lat), alt)| TemporalState {
                    timestamp: MissionElapsedTime::new(*t),
                    state: Position {
                        longitude: *lon,
                        latitude: *lat,
                        altitude: Length::new::<meter>(*alt),
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
    use crate::io::posatt::interpolate_sampled;
    use crate::types::Sat03B;

    #[test]
    fn a_fitted_table_interpolates_inside_a_pass_and_not_across_passes() {
        let csv = "time,lon,lat,alt_m\n100.0,10.0,-5.0,508000\n110.0,10.6,-4.4,508000\n120.0,11.2,-3.8,508000\n7000.0,50.0,20.0,508000\n7010.0,50.6,20.6,508000\n";
        let file = OrbitFitFile::from_reader(csv.as_bytes()).unwrap();
        assert_eq!(file.len(), 5);
        let traj = trajectory::<Sat03B>(&[file]);
        let p = interpolate_sampled(&traj, MissionElapsedTime::new(105.0)).unwrap();
        assert!((p.longitude - 10.3).abs() < 1e-9 && (p.latitude + 4.7).abs() < 1e-9);
        // 两次过境之间隔了 6880 s，不能连线
        assert!(interpolate_sampled(&traj, MissionElapsedTime::new(3000.0)).is_none());
    }

    #[test]
    fn a_malformed_line_is_an_error() {
        assert!(OrbitFitFile::from_reader("time,lon,lat,alt_m\n1,2\n".as_bytes()).is_err());
        assert!(OrbitFitFile::from_reader("time,lon,lat,alt_m\n1,x,3,4\n".as_bytes()).is_err());
    }
}
