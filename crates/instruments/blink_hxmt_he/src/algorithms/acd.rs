//! 候选窗口的 ACD 符合计数。
//!
//! 电子必须穿过反符合塑闪才能到达 NaI/CsI，ACD 几乎必然同时着火；伽马对
//! 薄塑闪只有几个百分点的康普顿效率。候选窗内 ACD 着火事例占比因此是
//! REP（电子沉降）与 TGF（伽马）的直接物理判别量，与候选密度（train）
//! 判据正交。搜索时统计并随候选保存；`blink acd-audit` 用同一实现离线
//! 复算认证样本。

use crate::types::Event;
use blink_core::traits::Event as _;
use blink_core::types::AcdCounts;

/// 基线窗与搜索邻域一致：候选两侧各 1 s，紧贴候选的 10 ms hollow 挖掉，
/// 避免瞬变自身拖尾污染基线（对应 SearchConfig 的 neighbor / hollow）。
const NEIGHBOR_S: f64 = 1.0;
const HOLLOW_S: f64 = 0.01;

/// 对时间有序的 kept 事例数组统计候选窗 [start, stop]（闭区间）与两侧
/// 基线窗 [start−1s, start−10ms) ∪ (stop+10ms, stop+1s] 的 ACD 计数。
/// 窗口越过数据边缘时按实际存在的事例截断——搜索按小时分块，基线在小时
/// 边缘天然截断，审计侧同一行为，两边数字可比。
pub fn acd_counts(events: &[Event], start_met: f64, stop_met: f64) -> AcdCounts {
    // (总数, 任意着火数, ≥2 块着火数)；区间端点由调用处的 partition_point 谓词决定
    let count = |i0: usize, i1: usize| -> (u32, u32, u32) {
        let mut n = 0u32;
        let mut n_acd = 0u32;
        let mut n_multi = 0u32;
        for event in &events[i0..i1] {
            n += 1;
            let fired = event.acds.iter().filter(|&&b| b).count();
            if fired >= 1 {
                n_acd += 1;
            }
            if fired >= 2 {
                n_multi += 1;
            }
        }
        (n, n_acd, n_multi)
    };
    let below = |t: f64| events.partition_point(|e| e.time().met() < t);
    let at_or_below = |t: f64| events.partition_point(|e| e.time().met() <= t);

    let (n, n_acd, n_acd_multi) = count(below(start_met), at_or_below(stop_met));
    let (n_left, acd_left, _) = count(below(start_met - NEIGHBOR_S), below(start_met - HOLLOW_S));
    let (n_right, acd_right, _) =
        count(at_or_below(stop_met + HOLLOW_S), at_or_below(stop_met + NEIGHBOR_S));

    AcdCounts {
        n,
        n_acd,
        n_acd_multi,
        n_bg: n_left + n_right,
        n_acd_bg: acd_left + acd_right,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Detector, Scintillator};
    use blink_core::types::MissionElapsedTime;

    fn event(met: f64, fired_paddles: usize) -> Event {
        let mut acds = [false; 18];
        acds.iter_mut().take(fired_paddles).for_each(|b| *b = true);
        Event::new(
            MissionElapsedTime::new(met),
            100,
            Detector {
                id: 0,
                scintillator: Scintillator::Csi,
            },
            false,
            acds,
        )
    }

    #[test]
    fn window_and_baseline_partition() {
        // 候选窗 [10.0, 10.001]：3 事例，其中 1 个单块着火、1 个三块着火。
        // 左基线 [9.0, 9.99)：2 事例 1 着火；右基线 (10.011, 11.001]：1 事例未着火。
        // hollow 内（9.995、10.005）的事例任何窗都不计。
        let events = vec![
            event(9.5, 0),
            event(9.8, 1),
            event(9.995, 5),  // 左 hollow，丢弃
            event(10.0, 0),   // 窗内（闭区间左端）
            event(10.0005, 1),
            event(10.001, 3), // 窗内（闭区间右端）
            event(10.005, 5), // 右 hollow，丢弃
            event(10.5, 0),
        ];
        let counts = acd_counts(&events, 10.0, 10.001);
        assert_eq!(counts.n, 3);
        assert_eq!(counts.n_acd, 2);
        assert_eq!(counts.n_acd_multi, 1);
        assert_eq!(counts.n_bg, 3);
        assert_eq!(counts.n_acd_bg, 1);
    }

    #[test]
    fn empty_and_edge_truncation() {
        let counts = acd_counts(&[], 10.0, 10.001);
        assert_eq!((counts.n, counts.n_bg), (0, 0));
        // 数据只覆盖候选之后：左基线截断为空，不出错
        let events = vec![event(10.0, 1), event(10.7, 2)];
        let counts = acd_counts(&events, 10.0, 10.001);
        assert_eq!(counts.n, 1);
        assert_eq!(counts.n_acd, 1);
        assert_eq!(counts.n_bg, 1);
        assert_eq!(counts.n_acd_bg, 1);
    }
}
