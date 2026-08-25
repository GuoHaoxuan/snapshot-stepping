//! GTI（好时间区间）。
//!
//! GRM 的小时文件并不覆盖整点到整点：实测一天 24 个文件的 GTI 合计只有
//! 86400 s 的 85%，缺口全部落在南大西洋异常区（实测 9 个缺口的中点经纬度
//! 无一例外在 lon ∈ [−88°, 0°]、lat ∈ [−29°, −8°]）。也就是说 L1B 的 GTI
//! 已经把 SAA 排除掉了，曝光核算直接用它即可，不必另立 SAA 判据。

pub(super) struct GtiHdu {
    pub start: Vec<f64>,
    pub stop: Vec<f64>,
}

impl GtiHdu {
    pub fn from_fptr(fptr: &mut fitsio::FitsFile) -> Result<Self, fitsio::errors::Error> {
        let gti = fptr.hdu("GTI")?;

        let start = gti.read_col::<f64>(fptr, "START")?;
        let stop = gti.read_col::<f64>(fptr, "STOP")?;

        Ok(Self { start, stop })
    }

    /// GTI 与 `[from, to]` 的交集长度（秒）。
    pub fn seconds_within(&self, from: f64, to: f64) -> f64 {
        self.start
            .iter()
            .zip(self.stop.iter())
            .map(|(start, stop)| (stop.min(to) - start.max(from)).max(0.0))
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gti(segments: &[(f64, f64)]) -> GtiHdu {
        GtiHdu {
            start: segments.iter().map(|(s, _)| *s).collect(),
            stop: segments.iter().map(|(_, e)| *e).collect(),
        }
    }

    #[test]
    fn segment_inside_window_counts_in_full() {
        assert_eq!(gti(&[(100.0, 400.0)]).seconds_within(0.0, 3600.0), 300.0);
    }

    #[test]
    fn segment_is_clipped_to_the_window() {
        // 相邻小时文件重叠 100 s：GTI 从整点前 100 s 起，越界那截不能计入本小时，
        // 否则重叠段会被相邻两个小时各记一次，曝光分母偏大。
        let hour = gti(&[(-100.0, 3600.0)]);
        assert_eq!(hour.seconds_within(0.0, 3600.0), 3600.0);
        // 尾端同理
        assert_eq!(gti(&[(3000.0, 4000.0)]).seconds_within(0.0, 3600.0), 600.0);
    }

    #[test]
    fn multiple_segments_add_up() {
        // SAA 穿越会把一小时切成两段
        let hour = gti(&[(0.0, 1000.0), (2500.0, 3600.0)]);
        assert_eq!(hour.seconds_within(0.0, 3600.0), 2100.0);
    }

    #[test]
    fn segment_outside_window_contributes_nothing() {
        // 负的交集必须夹到 0，不能把"离窗多远"当成曝光减掉
        assert_eq!(gti(&[(5000.0, 6000.0)]).seconds_within(0.0, 3600.0), 0.0);
        assert_eq!(
            gti(&[(5000.0, 6000.0), (0.0, 500.0)]).seconds_within(0.0, 3600.0),
            500.0
        );
    }

    #[test]
    fn empty_gti_is_zero_exposure() {
        assert_eq!(gti(&[]).seconds_within(0.0, 3600.0), 0.0);
    }
}
