/*
Filename: svom_orb_250101_00_v00.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PrimaryHDU    1 PrimaryHDU      39   ()
  1  ORB           1 BinTableHDU     92   2067R x 16C   [1D, 1E, 1E, 1E, 1E, 1E, 1E, 1E, 1E, 1E, 1E, 1E, 1E, 1E, 1E, 1E]
*/

use blink_core::types::{MissionElapsedTime, Position, TemporalState, Trajectory};

use crate::types::SvomGrm;

pub struct OrbFile {
    orb: OrbHdu,
}

impl OrbFile {
    pub fn from_fits_file(path: &str) -> Result<Self, fitsio::errors::Error> {
        let mut fptr = fitsio::FitsFile::open(path)?;

        let orb = OrbHdu::from_fptr(&mut fptr)?;

        Ok(Self { orb })
    }
}

struct OrbHdu {
    time: Vec<f64>,
    // x_j2000: Vec<f32>,
    // y_j2000: Vec<f32>,
    // z_j2000: Vec<f32>,
    // vx_j2000: Vec<f32>,
    // vy_j2000: Vec<f32>,
    // vz_j2000: Vec<f32>,
    // x_wgs84: Vec<f32>,
    // y_wgs84: Vec<f32>,
    // z_wgs84: Vec<f32>,
    // vx_wgs84: Vec<f32>,
    // vy_wgs84: Vec<f32>,
    // vz_wgs84: Vec<f32>,
    lon: Vec<f32>,
    lat: Vec<f32>,
    alt: Vec<f32>,
}

impl OrbHdu {
    fn from_fptr(fptr: &mut fitsio::FitsFile) -> Result<Self, fitsio::errors::Error> {
        let orb = fptr.hdu("ORB")?;

        let time = orb.read_col::<f64>(fptr, "TIME")?;
        // let x_j2000 = orb.read_col::<f32>(fptr, "X_J2000")?;
        // let y_j2000 = orb.read_col::<f32>(fptr, "Y_J2000")?;
        // let z_j2000 = orb.read_col::<f32>(fptr, "Z_J2000")?;
        // let vx_j2000 = orb.read_col::<f32>(fptr, "VX_J2000")?;
        // let vy_j2000 = orb.read_col::<f32>(fptr, "VY_J2000")?;
        // let vz_j2000 = orb.read_col::<f32>(fptr, "VZ_J2000")?;
        // let x_wgs84 = orb.read_col::<f32>(fptr, "X_WGS84")?;
        // let y_wgs84 = orb.read_col::<f32>(fptr, "Y_WGS84")?;
        // let z_wgs84 = orb.read_col::<f32>(fptr, "Z_WGS84")?;
        // let vx_wgs84 = orb.read_col::<f32>(fptr, "VX_WGS84")?;
        // let vy_wgs84 = orb.read_col::<f32>(fptr, "VY_WGS84")?;
        // let vz_wgs84 = orb.read_col::<f32>(fptr, "VZ_WGS84")?;
        let lon = orb.read_col::<f32>(fptr, "LON")?;
        let lat = orb.read_col::<f32>(fptr, "LAT")?;
        let alt = orb.read_col::<f32>(fptr, "ALT")?;

        Ok(Self {
            time,
            // x_j2000,
            // y_j2000,
            // z_j2000,
            // vx_j2000,
            // vy_j2000,
            // vz_j2000,
            // x_wgs84,
            // y_wgs84,
            // z_wgs84,
            // vx_wgs84,
            // vy_wgs84,
            // vz_wgs84,
            lon,
            lat,
            alt,
        })
    }
}

/// 地球平均半径（km），把早期文件里的地心距折算成高度用。
const EARTH_RADIUS_KM: f64 = 6371.0;

/// 判为「这个值是地心距而不是高度」的下限（km）。SVOM 的轨道高度实测在
/// 611–635 km，地心距则在 7000 km 上下，两者之间隔着一个数量级，取 3000
/// 落在空当里，怎么取都不会含糊。
const GEOCENTRIC_THRESHOLD_KM: f64 = 3000.0;

/// ORB 表的 `ALT` 列，头里写的是 `Altitude of sub-satellite point`、单位 km，
/// 全任务 99.8% 的样本也确实是 611–635 km 的轨道高度。但 2024 年早期有一小批
/// 文件在这一列里存的是**地心距**（约 7000 km）——同一列两种语义。数量级差
/// 十倍，判别没有歧义，这里统一折回高度。
///
/// 换算用的是球近似，纬度带来的椭球差可达 ±10 km；`altitude` 不参与任何判据，
/// 只是随候选记录的量，这个精度够用。
fn normalise_altitude(alt_km: f64) -> f64 {
    if alt_km > GEOCENTRIC_THRESHOLD_KM {
        alt_km - EARTH_RADIUS_KM
    } else {
        alt_km
    }
}

impl From<&OrbFile> for Trajectory<MissionElapsedTime<SvomGrm>, Position> {
    fn from(orb_file: &OrbFile) -> Self {
        let points = orb_file
            .orb
            .time
            .iter()
            .zip(orb_file.orb.lon.iter())
            .zip(orb_file.orb.lat.iter())
            .zip(orb_file.orb.alt.iter())
            .map(|(((t, lon), lat), alt)| TemporalState {
                timestamp: MissionElapsedTime::new(*t),
                state: Position {
                    longitude: *lon as f64,
                    latitude: *lat as f64,
                    altitude: uom::si::f64::Length::new::<uom::si::length::kilometer>(
                        normalise_altitude(*alt as f64),
                    ),
                },
            })
            .collect();

        Trajectory { points }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordinary_altitude_passes_through() {
        // 实测轨道高度落在 611–635 km
        assert_eq!(normalise_altitude(622.25), 622.25);
        assert_eq!(normalise_altitude(611.7), 611.7);
    }

    #[test]
    fn early_geocentric_values_are_folded_back() {
        // 2024 年早期的一批文件在同一列里存地心距
        assert!((normalise_altitude(7000.0) - 629.0).abs() < 1e-9);
    }
}
