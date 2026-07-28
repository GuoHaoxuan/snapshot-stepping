use serde::Serialize;

use crate::{
    traits::{Interpolatable, Temporal},
    types::TemporalState,
};

/// 允许外推的幅度，以所用那一段的长度为单位。
///
/// 两端都外推：时刻早于第一个采样点、或晚于最后一个采样点时，取最靠近的那
/// 一段把直线延出去。但只延有限的一小截 —— 表两端的采样相位差（HXMT 1K 的
/// 轨道/姿态是 0.25 s 一采，末点常落在整点前 0~0.24 s）就靠这一截兜住。
///
/// 超过这个幅度说明表是真的缺了一块，此时返回 `None`：外推一大截等于凭空
/// 造一个位置，比丢掉更糟（LEO 每秒走 7.6 km）。上层须把 `None` 记账，
/// 不要静默丢弃。
const EXTRAPOLATION_LIMIT: f64 = 1.0;

#[derive(Serialize, Debug)]
pub struct Trajectory<Time: Temporal, State: Interpolatable + Clone> {
    pub points: Vec<TemporalState<Time, State>>,
}

impl<Time: Temporal, State: Interpolatable + Clone> Trajectory<Time, State> {
    /// 取 `time` 时刻的状态。段内插值，两端按 [`EXTRAPOLATION_LIMIT`] 外推；
    /// 采样点不足两个、外推过头、或时间戳重复（lerp factor 为 NaN）时返回 `None`。
    pub fn interpolate(&self, time: Time) -> Option<TemporalState<Time, State>> {
        // 定一条直线要两个点
        if self.points.len() < 2 {
            return None;
        }

        // 找包住 time 的那一段；两端之外就夹到最近的一段，让 lerp factor
        // 落到 [0,1] 之外去表达外推。
        let last_segment = self.points.len() - 2;
        let mut i = 0;
        while i < last_segment && self.points[i + 1].timestamp < time {
            i += 1;
        }

        let t0 = self.points[i].timestamp;
        let t1 = self.points[i + 1].timestamp;
        let lerp_factor = time.lerp_factor(t0, t1);

        // NaN（t0 == t1）不满足任何比较，会一并落到这里被拒
        if !(-EXTRAPOLATION_LIMIT..=1.0 + EXTRAPOLATION_LIMIT).contains(&lerp_factor) {
            return None;
        }

        Some(TemporalState {
            timestamp: time,
            state: self.points[i]
                .state
                .interpolate(&self.points[i + 1].state, lerp_factor),
        })
    }

    pub fn window(&self, time: Time, half_width: Time::Duration) -> Self {
        let start_time = time - half_width;
        let end_time = time + half_width;

        Trajectory {
            points: self
                .points
                .iter()
                .filter(|point| point.timestamp >= start_time && point.timestamp <= end_time)
                .cloned()
                .collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Debug, PartialEq)]
    struct Scalar(f64);

    impl Interpolatable for Scalar {
        fn interpolate(&self, other: &Self, factor: f64) -> Self {
            Scalar(self.0 + (other.0 - self.0) * factor)
        }
    }

    impl Temporal for f64 {
        type Duration = f64;

        fn lerp_factor(self, start: Self, end: Self) -> f64 {
            (self - start) / (end - start)
        }
    }

    /// 采样点在 t = 0, 1, 2，值等于 t。
    fn ramp() -> Trajectory<f64, Scalar> {
        Trajectory {
            points: (0..3)
                .map(|i| TemporalState {
                    timestamp: i as f64,
                    state: Scalar(i as f64),
                })
                .collect(),
        }
    }

    fn at(trajectory: &Trajectory<f64, Scalar>, time: f64) -> Option<f64> {
        trajectory.interpolate(time).map(|point| point.state.0)
    }

    #[test]
    fn interpolates_inside_the_table() {
        assert_eq!(at(&ramp(), 0.0), Some(0.0));
        assert_eq!(at(&ramp(), 0.5), Some(0.5));
        assert_eq!(at(&ramp(), 2.0), Some(2.0));
    }

    #[test]
    fn extrapolates_past_the_last_point() {
        // 这一条就是原来被静默丢掉的情形：候选晚于最后一个星历采样点。
        assert_eq!(at(&ramp(), 2.25), Some(2.25));
        assert_eq!(at(&ramp(), 3.0), Some(3.0)); // 正好一整段，仍接受
    }

    #[test]
    fn extrapolates_before_the_first_point() {
        assert_eq!(at(&ramp(), -0.25), Some(-0.25));
        assert_eq!(at(&ramp(), -1.0), Some(-1.0));
    }

    #[test]
    fn refuses_to_extrapolate_far() {
        // 表真的缺了一大块时不能凭空造位置 —— 宁可返回 None 让上层记账。
        assert_eq!(at(&ramp(), 3.01), None);
        assert_eq!(at(&ramp(), -1.01), None);
        assert_eq!(at(&ramp(), 1e6), None);
    }

    #[test]
    fn degenerate_tables_are_rejected_not_panicking() {
        // 空表曾经会在 `points.len() - 1` 上下溢
        let empty: Trajectory<f64, Scalar> = Trajectory { points: Vec::new() };
        assert!(empty.interpolate(0.0).is_none());

        let single = Trajectory {
            points: vec![TemporalState {
                timestamp: 0.0,
                state: Scalar(7.0),
            }],
        };
        assert!(single.interpolate(0.0).is_none());

        // 时间戳重复 → lerp factor 是 NaN，必须拒掉而不是产出 NaN 状态
        let duplicated = Trajectory {
            points: vec![
                TemporalState {
                    timestamp: 5.0,
                    state: Scalar(1.0),
                },
                TemporalState {
                    timestamp: 5.0,
                    state: Scalar(2.0),
                },
            ],
        };
        assert!(duplicated.interpolate(5.0).is_none());
    }
}
