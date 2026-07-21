//! 掩掉-重建注入验证（spec §11 实验 2）。
//!
//! 在一段**未饱和**的真实数据里，把目标盒的一个时间区间人为当作 FIFO reset gap，
//! 用另两盒 cross-ref 重建它。cross-ref 只用参考盒事件建 shape、用目标盒 gap **外**
//! 的标定窗算 k，**不碰**目标盒 gap 内事件——那些正好留作真值。逐 bin 比"重建 vs
//! 真值"即可实测 U 地板、系统 bias、误差棒覆盖。

use crate::cli::InjectArgs;
use crate::commands::reconstruct::{
    gap_cov_row, gapbins_row, GAPBINS_HEADER, GAPCOV_HEADER,
};
use blink_hxmt_he::algorithms::saturation::{
    assign_gap_fill_channels, detect_fifo_reset_intervals, detect_unreliable_intervals,
    extract_packet_infos, reconstruct_gaps, reconstruct_met_channels,
    reconstruct_met_pulse_widths, reconstruct_met_times, reconstruct_with_wrap_tracking,
    unwrap_channel, BoxReconstructionData, SaturationInterval, SaturationType,
    UnreliableInterval, CHANNEL_SEC,
};
use blink_hxmt_he::io::level_1b::SciFile;

/// True iff `t` falls inside any injected gap interval `[s, e)`.
fn in_any_gap(t: f64, intervals: &[(f64, f64)]) -> bool {
    intervals.iter().any(|&(s, e)| t >= s && t < e)
}

/// Build the reconstructed event-stream rows (spec ①) for injection validation.
///
/// Reference boxes emit an EVT row for every in-window event (their source
/// counts C drive the S·diag(C)·Sᵀ leg). The target box is *masked* inside each
/// injected gap: in-gap events are withheld as truth, so only its out-of-gap
/// events become EVT rows, and each reconstructed filler becomes one FILL_GAP
/// row. Row layout matches `reconstruct`'s events.csv exactly.
fn injected_event_rows(
    boxes: &[(String, BoxReconstructionData)],
    target_idx: usize,
    intervals: &[(f64, f64)],
    fillers: &[(f64, u16, u8)],
    met_min: f64,
    met_max: f64,
) -> Vec<String> {
    let mut rows = Vec::new();
    for (bi, (name, data)) in boxes.iter().enumerate() {
        let is_target = bi == target_idx;
        for (i, &t) in data.events.iter().enumerate() {
            if t < met_min || t > met_max {
                continue;
            }
            // Mask the target box inside injected gaps — those events are truth.
            if is_target && in_any_gap(t, intervals) {
                continue;
            }
            let ch = data.channels[i];
            let raw = if ch == CHANNEL_SEC { 0 } else { unwrap_channel(ch) };
            rows.push(format!(
                "{},EVT,{:.6},{},{},-1,-1",
                name, t, raw, data.pulse_widths[i]
            ));
        }
        if is_target {
            for &(t, ch, pw) in fillers {
                if t < met_min || t > met_max {
                    continue;
                }
                rows.push(format!(
                    "{},FILL_GAP,{:.6},{},{},-1,-1",
                    name,
                    t,
                    unwrap_channel(ch),
                    pw
                ));
            }
        }
    }
    rows
}

/// 在目标盒注入人为 gap。返回：
/// - 注入后的 target（`gaps`/`unreliable` 设为注入区间；`events` 等原样保留，供
///   标定窗（gap 外）与真值提取用）；
/// - 每个注入区间内目标盒的**真值**事件 (met, channel, pulse_width)。
pub fn inject_gaps(
    target: &BoxReconstructionData,
    intervals: &[(f64, f64)],
) -> (BoxReconstructionData, Vec<Vec<(f64, u16, u8)>>) {
    let gaps: Vec<SaturationInterval> = intervals
        .iter()
        .map(|&(s, e)| SaturationInterval {
            start_met: s,
            stop_met: e,
            gap_seconds: e - s,
            prev_pkt_idx: 0,
            next_pkt_idx: 0,
            saturation_type: SaturationType::FifoReset,
        })
        .collect();

    let mut unreliable = target.unreliable.clone();
    for &(s, e) in intervals {
        unreliable.push(UnreliableInterval { start: s, stop: e });
    }

    // 真值：每注入区间 [s,e) 内目标盒的真实事件（cross-ref 完全不用它们）
    let truth: Vec<Vec<(f64, u16, u8)>> = intervals
        .iter()
        .map(|&(s, e)| {
            let a = target.events.partition_point(|&t| t < s);
            let b = target.events.partition_point(|&t| t < e);
            (a..b)
                .map(|i| (target.events[i], target.channels[i], target.pulse_widths[i]))
                .collect()
        })
        .collect();

    let injected = BoxReconstructionData {
        events: target.events.clone(),
        channels: target.channels.clone(),
        pulse_widths: target.pulse_widths.clone(),
        gaps,
        packets: target.packets.clone(),
        packet_events: target.packet_events.clone(),
        unreliable,
    };
    (injected, truth)
}

/// 注入验证命令：在未饱和数据的目标盒注入假 gap，cross-ref 重建，输出真值与 filler。
pub fn cmd_inject(args: &InjectArgs, boxes: &[(String, SciFile, f64)]) {
    eprintln!("Preparing injection data...");
    let mut box_data: Vec<(String, BoxReconstructionData)> = Vec::new();
    for (box_name, sci, offset) in boxes {
        let events = reconstruct_met_times(sci, *offset);
        let channels = reconstruct_met_channels(sci, *offset);
        let pulse_widths = reconstruct_met_pulse_widths(sci, *offset);
        let gaps = detect_fifo_reset_intervals(sci, *offset);
        let packets = extract_packet_infos(sci, *offset);
        let packet_events: Vec<Vec<f64>> = reconstruct_with_wrap_tracking(sci, *offset)
            .into_iter()
            .map(|mut times| {
                times.retain(|t| !t.is_nan());
                times.sort_by(|a, b| a.partial_cmp(b).unwrap());
                times
            })
            .collect();
        let unreliable = detect_unreliable_intervals(&gaps, &packets, &packet_events);
        box_data.push((
            box_name.clone(),
            BoxReconstructionData {
                events, channels, pulse_widths, gaps, packets, packet_events, unreliable,
            },
        ));
    }

    let ti = box_data
        .iter()
        .position(|(n, _)| n.eq_ignore_ascii_case(&args.target))
        .unwrap_or_else(|| panic!("target box {} not found", args.target));

    let met = args.window.trigger_met();
    let intervals: Vec<(f64, f64)> = args
        .at
        .iter()
        .map(|&off| {
            let c = met + off;
            (c - args.width / 2.0, c + args.width / 2.0)
        })
        .collect();

    let (injected, truth) = inject_gaps(&box_data[ti].1, &intervals);
    let refs: Vec<&BoxReconstructionData> = box_data
        .iter()
        .enumerate()
        .filter(|&(j, _)| j != ti)
        .map(|(_, (_, d))| d)
        .collect();

    let (gap_results, _rw) = reconstruct_gaps(&injected, &refs);
    let banded = assign_gap_fill_channels(&injected, &refs, &gap_results);

    // 真值 + 重建 filler → 供 Python 逐 gap/逐 bin 比较（truth 已知，可直接检误差）
    println!("gap_id,kind,met,channel");
    let mut n_truth = 0usize;
    for (gi, tv) in truth.iter().enumerate() {
        for &(t, ch, _pw) in tv {
            let raw = if ch == CHANNEL_SEC { 0 } else { unwrap_channel(ch) };
            println!("{gi},TRUTH,{t:.6},{raw}");
            n_truth += 1;
        }
    }
    let mut n_fill = 0usize;
    // Collect banded fillers (met, channel, pulse_width) for the spec-① stream.
    let mut all_fillers: Vec<(f64, u16, u8)> = Vec::new();
    for (r, b) in gap_results.iter().zip(banded.iter()) {
        for ((&t, &ch), &pw) in r
            .filled_events
            .iter()
            .zip(b.channels.iter())
            .zip(b.pulse_widths.iter())
        {
            println!("{},FILL,{:.6},{}", r.gap_idx, t, unwrap_channel(ch));
            all_fillers.push((t, ch, pw));
            n_fill += 1;
        }
    }

    // ── spec ①/②/③ tables for the new analytic covariance (recovery_cov.py) ──
    // The injected `intervals` become the target box's gaps; ref_idx in the
    // GapCovBlock is the index into `refs`, which is box_data with the target
    // removed — map it back to the box name by skipping the target position ti.
    if let Some(path) = &args.gapcov_out {
        let rows: Vec<String> = gap_results
            .iter()
            .map(|gr| {
                let gap = &injected.gaps[gr.gap_idx];
                gap_cov_row(
                    gr.gap_idx,
                    &box_data[ti].0,
                    gr.has_cross_ref,
                    gap.start_met,
                    gap.stop_met,
                    &gr.cov,
                    |local| {
                        let global = if local < ti { local } else { local + 1 };
                        box_data[global].0.clone()
                    },
                )
            })
            .collect();
        write_table(path, GAPCOV_HEADER, &rows, "gapcov");
    }

    if let Some(path) = &args.gapbins_out {
        let mut rows: Vec<String> = Vec::new();
        for gr in &gap_results {
            for b in &gr.bins {
                rows.push(gapbins_row(gr.gap_idx, b));
            }
        }
        write_table(path, GAPBINS_HEADER, &rows, "gapbins");
    }

    if let Some(path) = &args.events_out {
        let met_min = args.window.met_min();
        let met_max = args.window.met_max();
        let rows = injected_event_rows(&box_data, ti, &intervals, &all_fillers, met_min, met_max);
        write_table(
            path,
            "box,type,met,channel,pulse_width,pkt_idx,evt_idx",
            &rows,
            "events",
        );
    }

    if let Some(path) = &args.truth_out {
        let mut rows: Vec<String> = Vec::new();
        for (gi, tv) in truth.iter().enumerate() {
            let gr = gap_results.iter().find(|r| r.gap_idx == gi);
            let (nfill, xref) = gr.map_or((0, false), |r| (r.filled_events.len(), r.has_cross_ref));
            let (ts, te) = intervals[gi];
            rows.push(format!(
                "{gi},{ts:.6},{te:.6},{},{},{}",
                tv.len(),
                nfill,
                if xref { "crossref" } else { "degenerate" }
            ));
        }
        write_table(
            path,
            "gap_id,t_start,t_stop,n_truth,n_fill,type",
            &rows,
            "truth",
        );
    }

    let n_xref = gap_results.iter().filter(|r| r.has_cross_ref).count();
    eprintln!(
        "injected {} gaps on box {} (width {:.3}s): {} truth, {} fillers, {} cross-ref",
        intervals.len(), box_data[ti].0, args.width, n_truth, n_fill, n_xref,
    );
    for (gi, tv) in truth.iter().enumerate() {
        let gr = gap_results.iter().find(|r| r.gap_idx == gi);
        let (nfill, xref) = gr.map_or((0, false), |r| (r.filled_events.len(), r.has_cross_ref));
        eprintln!(
            "  gap[{gi}] @ T0{:+.3}s: truth={} fill={} {}",
            args.at.get(gi).copied().unwrap_or(f64::NAN),
            tv.len(),
            nfill,
            if xref { "cross-ref" } else { "DEGENERATE" },
        );
    }
}

/// Write `header` + `rows` to `path`, logging success/failure to stderr.
fn write_table(path: &std::path::Path, header: &str, rows: &[String], label: &str) {
    use std::io::Write;
    match std::fs::File::create(path) {
        Ok(mut f) => {
            let body = std::iter::once(header.to_string())
                .chain(rows.iter().cloned())
                .collect::<Vec<_>>()
                .join("\n");
            if let Err(e) = writeln!(f, "{body}") {
                eprintln!("  WARN: writing {label} {}: {e}", path.display());
            } else {
                eprintln!("  {label} table → {} ({} rows)", path.display(), rows.len());
            }
        }
        Err(e) => eprintln!("  WARN: cannot create {label} file {}: {e}", path.display()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn synth() -> BoxReconstructionData {
        BoxReconstructionData {
            events: vec![0.5, 1.02, 1.05, 1.08, 1.5],
            channels: vec![10, 20, 30, 40, 50],
            pulse_widths: vec![60; 5],
            gaps: Vec::new(),
            packets: Vec::new(),
            packet_events: Vec::new(),
            unreliable: Vec::new(),
        }
    }

    #[test]
    fn inject_sets_gap_and_extracts_truth() {
        let (inj, truth) = inject_gaps(&synth(), &[(1.0, 1.1)]);
        // 注入区间成为 gap
        assert_eq!(inj.gaps.len(), 1);
        assert!((inj.gaps[0].start_met - 1.0).abs() < 1e-12);
        assert!((inj.gaps[0].stop_met - 1.1).abs() < 1e-12);
        // 进入 unreliable（标定窗排除注入区间内的目标事件）
        assert!(inj
            .unreliable
            .iter()
            .any(|u| (u.start - 1.0).abs() < 1e-12 && (u.stop - 1.1).abs() < 1e-12));
        // events 原样保留（真值需要）
        assert_eq!(inj.events.len(), 5);
        // 真值 = 区间 [1.0,1.1) 内的 3 个事件（1.02/1.05/1.08）
        assert_eq!(truth.len(), 1);
        assert_eq!(truth[0].len(), 3);
        assert!((truth[0][0].0 - 1.02).abs() < 1e-12);
        assert_eq!(truth[0][0].1, 20);
        assert_eq!(truth[0][2].1, 40);
    }

    #[test]
    fn injected_stream_masks_target_in_gap_and_keeps_refs() {
        // Target A: 0.5 (out), 1.05 (in-gap → masked), 1.5 (out).
        let target = BoxReconstructionData {
            events: vec![0.5, 1.05, 1.5],
            channels: vec![10, 20, 30],
            pulse_widths: vec![60; 3],
            gaps: Vec::new(),
            packets: Vec::new(),
            packet_events: Vec::new(),
            unreliable: Vec::new(),
        };
        // Ref B: all events kept as EVT.
        let refb = BoxReconstructionData {
            events: vec![1.02, 1.08],
            channels: vec![40, 50],
            pulse_widths: vec![60; 2],
            gaps: Vec::new(),
            packets: Vec::new(),
            packet_events: Vec::new(),
            unreliable: Vec::new(),
        };
        let boxes = vec![("A".to_string(), target), ("B".to_string(), refb)];
        let intervals = [(1.0, 1.1)];
        let fillers = vec![(1.03, 20u16, 60u8), (1.07, 30u16, 60u8)];
        let rows = injected_event_rows(&boxes, 0, &intervals, &fillers, 0.0, 10.0);

        // A: 0.5 and 1.5 EVT (NOT 1.05, masked); two FILL_GAP fillers.
        let a_evt: Vec<&String> = rows
            .iter()
            .filter(|r| r.starts_with("A,EVT"))
            .collect();
        assert_eq!(a_evt.len(), 2, "target A keeps only out-of-gap EVT: {rows:?}");
        assert!(a_evt.iter().any(|r| r.contains(",0.500000,")));
        assert!(a_evt.iter().any(|r| r.contains(",1.500000,")));
        assert!(!rows.iter().any(|r| r.starts_with("A,EVT") && r.contains(",1.050000,")));
        let a_fill = rows.iter().filter(|r| r.starts_with("A,FILL_GAP")).count();
        assert_eq!(a_fill, 2, "target A gets one FILL_GAP per filler");
        // B: both events EVT, no fillers.
        let b_evt = rows.iter().filter(|r| r.starts_with("B,EVT")).count();
        assert_eq!(b_evt, 2, "ref B keeps every event as EVT");
        assert!(!rows.iter().any(|r| r.starts_with("B,FILL_GAP")));
    }

    #[test]
    fn in_any_gap_is_half_open() {
        let g = [(1.0, 1.1)];
        assert!(in_any_gap(1.0, &g), "start inclusive");
        assert!(in_any_gap(1.05, &g));
        assert!(!in_any_gap(1.1, &g), "stop exclusive");
        assert!(!in_any_gap(0.9, &g));
    }

    #[test]
    fn multiple_intervals_each_get_truth() {
        let (inj, truth) = inject_gaps(&synth(), &[(0.4, 0.6), (1.4, 1.6)]);
        assert_eq!(inj.gaps.len(), 2);
        assert_eq!(truth.len(), 2);
        assert_eq!(truth[0].len(), 1); // 0.5
        assert_eq!(truth[1].len(), 1); // 1.5
        assert!((truth[0][0].0 - 0.5).abs() < 1e-12);
    }
}
