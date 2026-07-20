//! Band-free 能量恢复：给 gap-fill 的 filler 事例确定性地补占位 channel。
//!
//! 与时间摆放同一哲学：不抽样、不引入 RNG。channel 取自**参考箱 in-gap**
//! 分布的等间隔分位（复现谱形，含尾），再用**位反转（van der Corput）**排列
//! 撒到时间有序的 filler 槽上，使 channel 与窗内时间无关（消除排序造成的
//! 假时间-能量漂移）。同一 1B 输入 → 逐字节相同输出。
//!
//! 设计与验证见
//! `docs/superpowers/specs/2026-07-03-eband-gapfill-prototype-design.md`。

/// SEC（秒脉冲）槽位的哨兵 channel：非真实道址，能量恢复时跳过。
pub const CHANNEL_SEC: u16 = u16::MAX;

/// 1B 原始 8-bit 道址 → wrapped 道址（raw < 20 表示 256+raw）。
/// 与 types::Event::channel() 的 pulse-height wrap 语义一致。
pub fn wrap_channel(raw: u8) -> u16 {
    if raw < 20 {
        raw as u16 + 256
    } else {
        raw as u16
    }
}

/// wrapped 道址 → 1B 原始 8-bit 道址（CSV 输出用，保持原始约定）。
pub fn unwrap_channel(ch: u16) -> u8 {
    if ch >= 256 {
        (ch - 256) as u8
    } else {
        ch as u8
    }
}

/// 从已排序样本里取第 ell 个（共 n 个）等间隔分位值：分位 (ell+0.5)/n。
/// n 个分位铺满经验 CDF，复现分布形状（含尾），零抽样噪声。
/// 样本为空时返回 0（调用方保证非空；空只在无任何参考+无标定窗时发生）。
pub fn quantile_value(sorted: &[u16], ell: usize, n: usize) -> u16 {
    if sorted.is_empty() {
        return 0;
    }
    let n = n.max(1);
    let q = (ell as f64 + 0.5) / n as f64;
    let idx = ((q * sorted.len() as f64) as usize).min(sorted.len() - 1);
    sorted[idx]
}

/// 基-2 radical inverse（van der Corput）：把 i 的二进制位倒过来当小数。
fn radical_inverse2(mut i: usize) -> f64 {
    let mut f = 0.0f64;
    let mut b = 0.5f64;
    while i > 0 {
        f += (i & 1) as f64 * b;
        i >>= 1;
        b *= 0.5;
    }
    f
}

/// 低差异（位反转 / van der Corput）排列：返回长度 n 的向量，
/// ranks[k] = 时间第 k 个槽位该拿的 channel 秩（0=最软）。
/// 读时间顺序时秩序列低差异散开，任意时间子段都拿到软硬均匀混合，
/// 从而 channel 与窗内时间去相关。n 为 2 的幂时即经典比特反转。
pub fn lowdisc_ranks(n: usize) -> Vec<usize> {
    if n == 0 {
        return Vec::new();
    }
    let phi: Vec<f64> = (0..n).map(radical_inverse2).collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| phi[a].total_cmp(&phi[b]));
    let mut ranks = vec![0usize; n];
    for (k, &idx) in order.iter().enumerate() {
        ranks[idx] = k;
    }
    ranks
}

use super::detect::{BoxReconstructionData, ReconstructedGap, UnreliableInterval};

/// gap-fill filler 的能量恢复结果：channels 与对应 gap 的
/// `ReconstructedGap.filled_events` 一一对应、同序。channel 为占位分位值
/// （取自参考箱 in-gap 分布），只支持光变，不支持逐事例谱拟合。
#[derive(Debug, Clone)]
pub struct GapFillChannels {
    pub gap_idx: usize,
    pub channels: Vec<u16>,
    /// 与 channels 一一对应的恢复脉宽（与 channel 取自同一参考事例，保留关联）
    pub pulse_widths: Vec<u8>,
    /// 谱来自子窗参考箱 in-gap（正常路径）的 filler 数
    pub n_from_window: usize,
    /// 谱退化到整 gap 参考池的 filler 数（子窗空、gap 内仍有参考）
    pub n_from_whole_gap: usize,
    /// 谱退化到 target 邻标定窗的 filler 数（三箱共饱和；能量/NaI 会偏，次优）
    pub n_from_calib: usize,
    /// 无任何谱来源、channel/pw 记 0 的 filler 数（极端：连邻窗都空）
    pub n_unfilled: usize,
}

/// fallback 标定窗半宽：无任何参考 in-gap 数据时退化用 target 邻窗谱（次优，M3）。
const CALIB_MARGIN: f64 = 0.5;
/// 谱子窗目标宽度：真实 ~30ms reset gap → 1 窗；长 gap 细分以承载窗间谱演化。
const WIN_TARGET: f64 = 0.05;

fn in_unreliable(t: f64, intervals: &[UnreliableInterval]) -> bool {
    intervals.iter().any(|iv| t >= iv.start && t <= iv.stop)
}

/// 按 met 排序的 (met, channel, pulse_width) 三元组，剔除 NaN 时间、SEC 槽、
/// 该 box 的 unreliable 区间。channel 与 pulse_width 来自同一事例（保留关联）。
fn sorted_triples(b: &BoxReconstructionData) -> Vec<(f64, u16, u8)> {
    let mut v: Vec<(f64, u16, u8)> = b
        .events
        .iter()
        .zip(b.channels.iter())
        .zip(b.pulse_widths.iter())
        .filter(|&((&t, &c), _)| !t.is_nan() && c != CHANNEL_SEC && !in_unreliable(t, &b.unreliable))
        .map(|((&t, &c), &w)| (t, c, w))
        .collect();
    v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    v
}

/// [lo, hi) 内的 (channel, pulse_width) 对，按 channel 稳定排序（供分位取样）。
fn pairs_in(triples: &[(f64, u16, u8)], lo: f64, hi: f64) -> Vec<(u16, u8)> {
    let a = triples.partition_point(|&(t, _, _)| t < lo);
    let b = triples.partition_point(|&(t, _, _)| t < hi);
    let mut out: Vec<(u16, u8)> = triples[a..b].iter().map(|&(_, c, w)| (c, w)).collect();
    out.sort_by(|x, y| x.0.cmp(&y.0).then(x.1.cmp(&y.1)));
    out
}

/// 长度为 len 的已排序数组在等间隔分位 (ell+0.5)/n 处的下标。
fn quantile_index(len: usize, ell: usize, n: usize) -> usize {
    if len == 0 {
        return 0;
    }
    let q = (ell as f64 + 0.5) / n.max(1) as f64;
    ((q * len as f64) as usize).min(len - 1)
}

/// 给 `reconstruct_gaps` 的 filler 事例补占位 channel（band-free 后处理）。
/// 不改变 filler 的数量与时刻——只补 channel。
///
/// 每 gap：谱形状取自参考箱 in-gap 分布（尺度无关，无需 k_tot）；按时间子窗
/// （`WIN_TARGET`）取该窗参考 channel 的等间隔分位，用位反转排列撒到时间有序的
/// filler 槽上（channel 与窗内时间去相关）。窗内无参考 → 退化到整 gap 参考池
/// → 再退化到 target 邻标定窗（次优）。全程确定性，无 RNG。
pub fn assign_gap_fill_channels(
    target: &BoxReconstructionData,
    references: &[&BoxReconstructionData],
    gap_results: &[ReconstructedGap],
) -> Vec<GapFillChannels> {
    let ref_triples: Vec<Vec<(f64, u16, u8)>> =
        references.iter().map(|r| sorted_triples(r)).collect();
    let tgt_triples = sorted_triples(target);
    let by_ch = |x: &(u16, u8), y: &(u16, u8)| x.0.cmp(&y.0).then(x.1.cmp(&y.1));

    let mut out = Vec::with_capacity(gap_results.len());
    for gr in gap_results {
        let gap = &target.gaps[gr.gap_idx];
        let (g_lo, g_hi) = (gap.start_met, gap.stop_met);
        let filled = &gr.filled_events;
        let n = filled.len();
        if n == 0 {
            out.push(GapFillChannels {
                gap_idx: gr.gap_idx,
                channels: Vec::new(),
                pulse_widths: Vec::new(),
                n_from_window: 0,
                n_from_whole_gap: 0,
                n_from_calib: 0,
                n_unfilled: 0,
            });
            continue;
        }

        // 一级 fallback：整 gap 参考池
        let mut whole_gap: Vec<(u16, u8)> =
            ref_triples.iter().flat_map(|p| pairs_in(p, g_lo, g_hi)).collect();
        whole_gap.sort_by(by_ch);
        // 二级 fallback：target 邻标定窗（无任何参考 in-gap 时）
        let mut calib: Vec<(u16, u8)> = pairs_in(&tgt_triples, g_lo - CALIB_MARGIN, g_lo);
        calib.extend(pairs_in(&tgt_triples, g_hi, g_hi + CALIB_MARGIN));
        calib.sort_by(by_ch);

        let d = g_hi - g_lo;
        let n_win = ((d / WIN_TARGET).round() as usize).max(1);
        let mut channels = vec![0u16; n];
        let mut pulse_widths = vec![0u8; n];
        let (mut n_from_window, mut n_from_whole_gap, mut n_from_calib, mut n_unfilled) =
            (0usize, 0usize, 0usize, 0usize);

        for wi in 0..n_win {
            let w_lo = g_lo + d * wi as f64 / n_win as f64;
            let w_hi = if wi + 1 == n_win {
                g_hi
            } else {
                g_lo + d * (wi + 1) as f64 / n_win as f64
            };
            let s = filled.partition_point(|&t| t < w_lo);
            let e = if wi + 1 == n_win { n } else { filled.partition_point(|&t| t < w_hi) };
            let n_w = e - s;
            if n_w == 0 {
                continue;
            }
            let mut src: Vec<(u16, u8)> =
                ref_triples.iter().flat_map(|p| pairs_in(p, w_lo, w_hi)).collect();
            src.sort_by(by_ch);
            // fallback 分级：子窗参考 → 整 gap 参考池 → target 邻窗(三箱共饱和)
            let (spectrum, kind): (&[(u16, u8)], usize) = if !src.is_empty() {
                (&src, 0)
            } else if !whole_gap.is_empty() {
                (&whole_gap, 1)
            } else {
                (&calib, 2)
            };
            if spectrum.is_empty() {
                n_unfilled += n_w;
                continue;
            }
            match kind {
                0 => n_from_window += n_w,
                1 => n_from_whole_gap += n_w,
                _ => n_from_calib += n_w,
            }
            // 分位取样：filler 与 channel 都取自参考事例，pulse_width 随之（保关联）
            let ranks = lowdisc_ranks(n_w);
            for k in 0..n_w {
                let idx = quantile_index(spectrum.len(), ranks[k], n_w);
                channels[s + k] = spectrum[idx].0;
                pulse_widths[s + k] = spectrum[idx].1;
            }
        }
        out.push(GapFillChannels {
            gap_idx: gr.gap_idx,
            channels,
            pulse_widths,
            n_from_window,
            n_from_whole_gap,
            n_from_calib,
            n_unfilled,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrap_unwrap_roundtrip() {
        assert_eq!(wrap_channel(19), 275);
        assert_eq!(wrap_channel(20), 20);
        assert_eq!(wrap_channel(0), 256);
        assert_eq!(wrap_channel(255), 255);
        assert_eq!(unwrap_channel(275), 19);
        assert_eq!(unwrap_channel(256), 0);
        assert_eq!(unwrap_channel(44), 44);
        for raw in 0u8..=255 {
            assert_eq!(unwrap_channel(wrap_channel(raw)), raw);
        }
    }

    #[test]
    fn quantile_walks_the_cdf() {
        let sorted = [30u16, 31, 40, 41];
        // n == len：等间隔分位逐一取到每个观测值
        let drawn: Vec<u16> = (0..4).map(|l| quantile_value(&sorted, l, 4)).collect();
        assert_eq!(drawn, vec![30, 31, 40, 41]);
        // n=2：分位 0.25 / 0.75 → 取到 31 与 41
        assert_eq!(quantile_value(&sorted, 0, 2), 31);
        assert_eq!(quantile_value(&sorted, 1, 2), 41);
        // n=1：分位 0.5 → 中位
        assert_eq!(quantile_value(&sorted, 0, 1), 40);
        // 空样本兜底
        assert_eq!(quantile_value(&[], 0, 1), 0);
    }

    #[test]
    fn lowdisc_ranks_bit_reversal() {
        // 手推：8 槽位比特反转 = [0,4,2,6,1,5,3,7]
        assert_eq!(lowdisc_ranks(8), vec![0, 4, 2, 6, 1, 5, 3, 7]);
        assert_eq!(lowdisc_ranks(1), vec![0]);
        assert_eq!(lowdisc_ranks(2), vec![0, 1]);
        assert_eq!(lowdisc_ranks(0), Vec::<usize>::new());
    }

    #[test]
    fn lowdisc_ranks_is_a_permutation() {
        for n in [1usize, 2, 3, 5, 7, 8, 13, 64, 100] {
            let r = lowdisc_ranks(n);
            assert_eq!(r.len(), n);
            let mut seen = r.clone();
            seen.sort_unstable();
            assert_eq!(seen, (0..n).collect::<Vec<_>>(), "n={n} not a permutation");
        }
    }

    #[test]
    fn lowdisc_ranks_deterministic() {
        assert_eq!(lowdisc_ranks(50), lowdisc_ranks(50));
    }

    #[test]
    fn lowdisc_first_half_spans_range() {
        // 低差异性：任意前缀均匀散开——前半段应同时含低秩与高秩。
        let r = lowdisc_ranks(8);
        let first_half = &r[..4];
        assert!(first_half.iter().any(|&x| x < 4), "前半无低秩");
        assert!(first_half.iter().any(|&x| x >= 4), "前半无高秩");
    }

    // ---- assign_gap_fill_channels (band-free 后处理) ----
    use crate::algorithms::saturation::detect::{SaturationInterval, SaturationType};

    fn si(lo: f64, hi: f64) -> SaturationInterval {
        SaturationInterval {
            start_met: lo,
            stop_met: hi,
            gap_seconds: hi - lo,
            prev_pkt_idx: 0,
            next_pkt_idx: 0,
            saturation_type: SaturationType::FifoReset,
        }
    }

    fn make_box(
        events: Vec<f64>,
        channels: Vec<u16>,
        gaps: Vec<SaturationInterval>,
    ) -> BoxReconstructionData {
        let pulse_widths = vec![60u8; channels.len()];
        BoxReconstructionData {
            events,
            channels,
            pulse_widths,
            gaps,
            packets: Vec::new(),
            packet_events: Vec::new(),
            unreliable: Vec::new(),
        }
    }

    fn gap0(filled: Vec<f64>) -> ReconstructedGap {
        let n = filled.len();
        ReconstructedGap {
            gap_idx: 0,
            filled_events: filled,
            n_lost: n,
            has_cross_ref: true,
            filler_weight: 0.0,
        }
    }

    /// [lo, hi) 内均匀铺 n 个等间隔时间戳（升序）。
    fn spread(lo: f64, hi: f64, n: usize) -> Vec<f64> {
        (0..n).map(|i| lo + (i as f64 + 0.5) * (hi - lo) / n as f64).collect()
    }

    #[test]
    fn constant_reference_spectrum_gives_constant_channel() {
        let target = make_box(vec![], vec![], vec![si(1.9, 2.1)]);
        let reference = make_box(spread(1.9, 2.1, 40), vec![100u16; 40], vec![]);
        let out = assign_gap_fill_channels(&target, &[&reference], &[gap0(spread(1.9, 2.1, 5))]);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].channels.len(), 5);
        assert!(out[0].channels.iter().all(|&c| c == 100), "{:?}", out[0].channels);
        // 正常路径:全部来自子窗参考,无 fallback
        assert_eq!(out[0].n_from_window, 5);
        assert_eq!(out[0].n_from_calib, 0);
        assert_eq!(out[0].n_unfilled, 0);
    }

    #[test]
    fn shape_from_reference_ingap_not_target_calib() {
        // target 邻标定窗全软(50)；参考箱 in-gap 全硬(200) → 应取 200
        let mut te = spread(1.4, 1.9, 40);
        te.extend(spread(2.1, 2.6, 40));
        let target = make_box(te, vec![50u16; 80], vec![si(1.9, 2.1)]);
        let reference = make_box(spread(1.9, 2.1, 40), vec![200u16; 40], vec![]);
        let out = assign_gap_fill_channels(&target, &[&reference], &[gap0(spread(1.9, 2.1, 8))]);
        assert!(
            out[0].channels.iter().all(|&c| c == 200),
            "应随参考 in-gap(200): {:?}",
            out[0].channels
        );
    }

    #[test]
    fn no_reference_falls_back_to_target_calib() {
        let mut te = spread(1.4, 1.9, 40);
        te.extend(spread(2.1, 2.6, 40));
        let target = make_box(te, vec![77u16; 80], vec![si(1.9, 2.1)]);
        let out = assign_gap_fill_channels(&target, &[], &[gap0(spread(1.9, 2.1, 5))]);
        assert!(
            out[0].channels.iter().all(|&c| c == 77),
            "应退化到标定窗(77): {:?}",
            out[0].channels
        );
        // 三箱共饱和 fallback 计数可见
        assert_eq!(out[0].n_from_calib, 5, "5 个 filler 应记为 calib fallback");
        assert_eq!(out[0].n_from_window, 0);
    }

    #[test]
    fn single_window_uses_bit_reversal_order() {
        // D=0.04 → 1 窗；参考 8 个不同 channel
        let reference = make_box(
            spread(1.98, 2.02, 8),
            vec![20u16, 40, 60, 80, 100, 120, 140, 160],
            vec![],
        );
        let target = make_box(vec![], vec![], vec![si(1.98, 2.02)]);
        let out = assign_gap_fill_channels(&target, &[&reference], &[gap0(spread(1.98, 2.02, 8))]);
        // quantile_value(src, ell, 8) == src[ell]；槽 k 取秩 lowdisc_ranks(8)=[0,4,2,6,1,5,3,7]
        assert_eq!(out[0].channels, vec![20, 100, 60, 140, 40, 120, 80, 160]);
        assert!(
            out[0].channels.windows(2).any(|w| w[0] > w[1]),
            "不应按时间单调递增（位反转去相关）"
        );
    }

    #[test]
    fn assignment_is_deterministic() {
        let rc: Vec<u16> = (0..40).map(|i| 20 + i as u16 * 6).collect();
        let reference = make_box(spread(1.9, 2.1, 40), rc, vec![]);
        let target = make_box(vec![], vec![], vec![si(1.9, 2.1)]);
        let a = assign_gap_fill_channels(&target, &[&reference], &[gap0(spread(1.9, 2.1, 17))]);
        let b = assign_gap_fill_channels(&target, &[&reference], &[gap0(spread(1.9, 2.1, 17))]);
        assert_eq!(a[0].channels, b[0].channels);
    }

    #[test]
    fn empty_gap_yields_empty_channels() {
        let target = make_box(vec![], vec![], vec![si(1.9, 2.1)]);
        let out = assign_gap_fill_channels(&target, &[], &[gap0(vec![])]);
        assert!(out[0].channels.is_empty());
    }
}
