use super::{Chunk, Interval};
use crate::algorithms::saturation::scan_saturation_intervals;
use crate::types::HxmtHe;
use blink_core::types::MissionElapsedTime;

impl Chunk {
    /// 三个机箱的饱和时间段合并（任一机箱饱和即算饱和）。
    ///
    /// 计算要把三个机箱的 1B 全部时间重建一遍，是本 pipeline 里最贵的一步，
    /// 因此结果缓存在 chunk 里 —— search 和 coverage 共用同一份。
    pub fn saturation_intervals(&self) -> &[Interval] {
        self.saturation_cache.get_or_init(|| {
            let mut all_intervals: Vec<Interval> = Vec::new();
            for ((_, sci_file), (_, offset)) in self.sci_files.iter().zip(self.stime_offsets.iter())
            {
                all_intervals.extend(scan_saturation_intervals(sci_file, *offset));
            }

            // 按起始时间排序
            all_intervals.sort_by(|a, b| a.0.cmp(&b.0));

            // 合并有重叠的区间（并集）
            let mut merged: Vec<Interval> = Vec::new();
            for interval in all_intervals {
                if let Some(last) = merged.last_mut()
                    && interval.0 <= last.1
                {
                    if interval.1 > last.1 {
                        last.1 = interval.1;
                    }
                    continue;
                }
                merged.push(interval);
            }

            merged
        })
    }

    /// 判断给定时间点是否处于饱和状态（二分查找）。
    pub fn check_saturation(&self, time: MissionElapsedTime<HxmtHe>) -> bool {
        is_in_intervals(self.saturation_intervals(), time)
    }
}

/// 二分查找判断时间点是否落在某个饱和区间内。
fn is_in_intervals(intervals: &[Interval], time: MissionElapsedTime<HxmtHe>) -> bool {
    let idx = intervals.partition_point(|interval| interval.1 < time);
    idx < intervals.len() && intervals[idx].0 <= time
}
