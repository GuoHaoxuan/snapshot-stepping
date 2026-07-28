use crate::traits::Interpolatable;
use serde::{Deserialize, Serialize};
use uom::si::f64::*;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Position {
    pub longitude: f64,
    pub latitude: f64,
    pub altitude: Length,
}

/// 把角度差折到 (-180, 180] —— 两点之间的短弧。
fn shortest_arc_degrees(delta: f64) -> f64 {
    delta - 360.0 * (delta / 360.0).round()
}

/// 把经度折回 `[origin, origin + 360)`。
fn wrap_into(longitude: f64, origin: f64) -> f64 {
    origin + (longitude - origin).rem_euclid(360.0)
}

impl Interpolatable for Position {
    /// 经度按短弧插值。
    ///
    /// 经度是周期量：两个采样点跨过取值区间的接缝时，裸的线性插值会把结果甩到
    /// 地球另一面 —— 359.98° 和 0.02° 的中点会算成 180°，而不是 0°。HXMT 1K 的
    /// `Orbit.Lon` 用 [0,360)，实测一天有 14 次跨 0/360 的采样间隔（相邻 0.25 s
    /// 采样点之间 |ΔLon| = 359.96°），落在这些间隔里的候选位置会整整错到对跖点。
    /// 位置错了闪电关联就是错的，所以这里按短弧插，最后折回与输入相同的区间。
    ///
    /// 纬度不做这个处理：它不是周期量，在 ±90° 之间连续变化。
    fn interpolate(&self, other: &Self, ratio: f64) -> Self {
        // 各任务的经度约定不同（HXMT 用 [0,360)，也有用 [-180,180) 的），
        // 按输入落在哪一边决定折回哪个区间，不改动约定本身。
        let origin = if self.longitude < 0.0 || other.longitude < 0.0 {
            -180.0
        } else {
            0.0
        };
        let delta = shortest_arc_degrees(other.longitude - self.longitude);

        Position {
            longitude: wrap_into(self.longitude + delta * ratio, origin),
            latitude: self.latitude + (other.latitude - self.latitude) * ratio,
            altitude: self.altitude + (other.altitude - self.altitude) * ratio,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::length::meter;

    fn at(lon: f64) -> Position {
        Position {
            longitude: lon,
            latitude: 10.0,
            altitude: Length::new::<meter>(550_000.0),
        }
    }

    fn lon_of(a: f64, b: f64, ratio: f64) -> f64 {
        at(a).interpolate(&at(b), ratio).longitude
    }

    #[test]
    fn interpolates_without_crossing_the_seam() {
        assert!((lon_of(100.0, 100.08, 0.5) - 100.04).abs() < 1e-9);
        assert!((lon_of(100.08, 100.0, 0.5) - 100.04).abs() < 1e-9);
    }

    #[test]
    fn crossing_the_seam_takes_the_short_arc() {
        // 这一条就是 bug 本体：裸线性插值会给出 180，即对跖点
        assert!((lon_of(359.98, 0.02, 0.5) - 0.0).abs() < 1e-9);
        assert!((lon_of(0.02, 359.98, 0.5) - 0.0).abs() < 1e-9);
        // 端点必须原样返回
        assert!((lon_of(359.98, 0.02, 0.0) - 359.98).abs() < 1e-9);
        assert!((lon_of(359.98, 0.02, 1.0) - 0.02).abs() < 1e-9);
    }

    #[test]
    fn result_stays_in_the_input_convention() {
        // [0,360) 进 → [0,360) 出
        for ratio in [-1.0, 0.0, 0.5, 1.0, 2.0] {
            let lon = lon_of(359.98, 0.02, ratio);
            assert!((0.0..360.0).contains(&lon), "ratio={ratio} lon={lon}");
        }
        // [-180,180) 进 → [-180,180) 出
        for ratio in [-1.0, 0.0, 0.5, 1.0, 2.0] {
            let lon = lon_of(-179.98, 179.94, ratio);
            assert!((-180.0..180.0).contains(&lon), "ratio={ratio} lon={lon}");
        }
    }

    #[test]
    fn extrapolation_also_takes_the_short_arc() {
        // 外推同样不能跨接缝乱走：ratio=2 是从 359.94 走两段 (2×0.04)，
        // 到 360.02，折回后是 0.02
        assert!((lon_of(359.94, 359.98, 2.0) - 0.02).abs() < 1e-9);
        // 反方向外推越过接缝：ratio=-1 从 0.02 退一段 (0.04) 到 -0.02 → 359.98
        assert!((lon_of(0.02, 0.06, -1.0) - 359.98).abs() < 1e-9);
    }

    #[test]
    fn latitude_and_altitude_are_plain_linear() {
        let a = Position {
            longitude: 0.0,
            latitude: -40.0,
            altitude: Length::new::<meter>(500_000.0),
        };
        let b = Position {
            longitude: 0.0,
            latitude: 40.0,
            altitude: Length::new::<meter>(600_000.0),
        };
        let mid = a.interpolate(&b, 0.5);
        assert!((mid.latitude - 0.0).abs() < 1e-9);
        assert!((mid.altitude.get::<meter>() - 550_000.0).abs() < 1e-6);
    }
}
