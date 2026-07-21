use blink_hxmt_he::algorithms::saturation::{
    assign_gap_fill_channels, detect_fifo_reset_intervals, detect_unreliable_intervals,
    extract_packet_infos, reconstruct_gaps, reconstruct_met_channels,
    reconstruct_met_pulse_widths, reconstruct_met_times, reconstruct_with_wrap_tracking,
    unwrap_channel, BoxReconstructionData, CHANNEL_SEC, GapBinInfo, GapBinKind, GapCovBlock,
    GapRefCalib,
};
use blink_hxmt_he::io::level_1b::SciFile;

use crate::cli::ReconstructArgs;

pub fn cmd_reconstruct(
    args: &ReconstructArgs,
    boxes: &[(String, SciFile, f64)],
    filter_box: &Option<String>,
) {
    let met_min = args.window.met_min();
    let met_max = args.window.met_max();

    eprintln!("Preparing reconstruction data...");
    let mut box_data: Vec<(String, BoxReconstructionData)> = Vec::new();
    for (box_name, sci, offset) in boxes {
        let events = reconstruct_met_times(sci, *offset);
        let channels = reconstruct_met_channels(sci, *offset);
        let pulse_widths = reconstruct_met_pulse_widths(sci, *offset);
        assert_eq!(
            events.len(),
            channels.len(),
            "events/channels misaligned for box {box_name}"
        );
        assert_eq!(
            events.len(),
            pulse_widths.len(),
            "events/pulse_widths misaligned for box {box_name}"
        );
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
        eprintln!(
            "  Box {}: {} events, {} gaps, {} unreliable, {} packets",
            box_name, events.len(), gaps.len(), unreliable.len(), packets.len()
        );
        box_data.push((
            box_name.clone(),
            BoxReconstructionData {
                events, channels, pulse_widths, gaps, packets, packet_events, unreliable,
            },
        ));
    }

    let original_events: Vec<(String, Vec<f64>, Vec<u16>, Vec<u8>)> = box_data
        .iter()
        .map(|(name, data)| {
            (name.clone(), data.events.clone(), data.channels.clone(), data.pulse_widths.clone())
        })
        .collect();

    eprintln!("Reconstructing (FIFO reset gaps)...");
    let mut all_filled: Vec<(String, Vec<(f64, u16, u8)>)> = Vec::new();
    // gap 协方差块表行（spec §13），写到 --gapcov-out 单独文件。观测事件回归权重 1
    // （普通泊松）；恢复引入的协方差全部收进每 gap 的块，不再挂逐粒子权重。
    let mut gapcov_rows: Vec<String> = Vec::new();
    // gap 格结构表行（spec ③），写到 --gapbins-out。每 (gap,1ms 格) 一行，只有
    // cross-ref gap 有 S filler↔参考结构（退化 gap 的 bins 为空、不产行）。
    let mut gapbins_rows: Vec<String> = Vec::new();

    for i in 0..box_data.len() {
        let refs_with_idx: Vec<(usize, &BoxReconstructionData)> = box_data
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(j, (_, d))| (j, d))
            .collect();
        let refs: Vec<&BoxReconstructionData> =
            refs_with_idx.iter().map(|(_, d)| *d).collect();

        let (gap_results, _ref_weights) = reconstruct_gaps(&box_data[i].1, &refs);
        let n_gap_filled: usize = gap_results.iter().map(|r| r.n_lost).sum();
        let n_gap_ref = gap_results.iter().filter(|r| r.has_cross_ref).count();
        let banded = assign_gap_fill_channels(&box_data[i].1, &refs, &gap_results);
        let n_calib: usize = banded.iter().map(|b| b.n_from_calib).sum();
        let n_unfill: usize = banded.iter().map(|b| b.n_unfilled).sum();
        if n_calib > 0 || n_unfill > 0 {
            eprintln!(
                "  WARN Box {}: energy/pulse-width fallback — {} fillers from target \
                 calib window (three-box co-saturation; energy & NaI/CsI degraded), \
                 {} left unfilled",
                box_data[i].0, n_calib, n_unfill,
            );
        }

        // 收集本盒各 gap 的协方差块行（与事件流一样受 filter_box 约束）。
        // gr.cov.refs 的 ref_idx 是 references 切片内下标 → refs_with_idx 映射回盒名。
        let box_shown = filter_box
            .as_ref()
            .map_or(true, |fb| box_data[i].0.eq_ignore_ascii_case(fb));
        if box_shown {
            for gr in &gap_results {
                let gap = &box_data[i].1.gaps[gr.gap_idx];
                gapcov_rows.push(gap_cov_row(
                    gr.gap_idx,
                    &box_data[i].0,
                    gr.has_cross_ref,
                    gap.start_met,
                    gap.stop_met,
                    &gr.cov,
                    |local| box_data[refs_with_idx[local].0].0.clone(),
                ));
                // spec ③:每 1ms 格一行（cross-ref gap 才有 bins）
                for b in &gr.bins {
                    gapbins_rows.push(gapbins_row(gr.gap_idx, b));
                }
            }
        }

        let mut gap_events: Vec<(f64, u16, u8)> = gap_results
            .iter()
            .zip(banded.iter())
            .flat_map(|(r, b)| {
                r.filled_events
                    .iter()
                    .copied()
                    .zip(b.channels.iter().copied())
                    .zip(b.pulse_widths.iter().copied())
                    .map(|((t, c), w)| (t, c, w))
            })
            .collect();
        gap_events.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        eprintln!(
            "  Box {}: gaps={} ({} evt, {} ref)",
            box_data[i].0, box_data[i].1.gaps.len(), n_gap_filled, n_gap_ref,
        );

        all_filled.push((box_data[i].0.clone(), gap_events));
    }

    // 事件流（均值）：观测事件 EVT + gap-fill FILL_GAP，不带权重列。
    // 均值 = 数所有行；协方差 = 由块表 + 事件流重算源计数装配（spec §7）。
    println!("box,type,met,channel,pulse_width,pkt_idx,evt_idx");
    for (box_name, _data) in box_data.iter() {
        let (obs_events, obs_channels, obs_pw) = original_events
            .iter()
            .find(|(n, _, _, _)| n == box_name)
            .map(|(_, e, c, w)| (e.as_slice(), c.as_slice(), w.as_slice()))
            .unwrap_or((&[], &[], &[]));
        let gap_events = all_filled
            .iter()
            .find(|(n, _)| n == box_name)
            .map(|(_, f)| f.as_slice())
            .unwrap_or(&[]);

        if let Some(fb) = filter_box {
            if !box_name.eq_ignore_ascii_case(fb) {
                continue;
            }
        }

        let mut n_obs = 0u64;
        let mut n_gap = 0u64;

        for (idx, &t) in obs_events.iter().enumerate() {
            if t >= met_min && t <= met_max {
                let ch = obs_channels[idx];
                let raw = if ch == CHANNEL_SEC { 0 } else { unwrap_channel(ch) };
                println!("{},EVT,{:.6},{},{},-1,-1", box_name, t, raw, obs_pw[idx]);
                n_obs += 1;
            }
        }
        for &(t, ch, pw) in gap_events {
            if t >= met_min && t <= met_max {
                println!("{},FILL_GAP,{:.6},{},{},-1,-1", box_name, t, unwrap_channel(ch), pw);
                n_gap += 1;
            }
        }

        eprintln!(
            "  Box {}: {} observed, {} gap-filled, bin={:.3}s",
            box_name, n_obs, n_gap, args.bin,
        );
    }

    // gap 协方差块表（spec §13）→ 单独文件
    if let Some(path) = &args.gapcov_out {
        use std::io::Write;
        match std::fs::File::create(path) {
            Ok(mut f) => {
                let body = std::iter::once(GAPCOV_HEADER.to_string())
                    .chain(gapcov_rows.iter().cloned())
                    .collect::<Vec<_>>()
                    .join("\n");
                if let Err(e) = writeln!(f, "{body}") {
                    eprintln!("  WARN: writing gapcov {}: {e}", path.display());
                } else {
                    eprintln!(
                        "  gap covariance block table → {} ({} gaps)",
                        path.display(),
                        gapcov_rows.len()
                    );
                }
            }
            Err(e) => eprintln!("  WARN: cannot create gapcov file {}: {e}", path.display()),
        }
    }

    // gap 格结构表（spec ③）→ 单独文件
    if let Some(path) = &args.gapbins_out {
        use std::io::Write;
        match std::fs::File::create(path) {
            Ok(mut f) => {
                let body = std::iter::once(GAPBINS_HEADER.to_string())
                    .chain(gapbins_rows.iter().cloned())
                    .collect::<Vec<_>>()
                    .join("\n");
                if let Err(e) = writeln!(f, "{body}") {
                    eprintln!("  WARN: writing gapbins {}: {e}", path.display());
                } else {
                    eprintln!(
                        "  gap bin structure table → {} ({} bins)",
                        path.display(),
                        gapbins_rows.len()
                    );
                }
            }
            Err(e) => eprintln!("  WARN: cannot create gapbins file {}: {e}", path.display()),
        }
    }
}

/// gap 块表(spec §13)表头,单独文件 `--gapcov-out`。事件流给均值,块表给协方差。
pub(crate) const GAPCOV_HEADER: &str = "gap_id,target_box,type,t_start,t_stop,ref_boxes,k,\
     c_ref_cal,c_a_cal,rho,r_pre,r_post,n_pre,n_post,maskable,sys_bias_flag,sys_bias_scale";

/// gap 格结构表(spec ③)表头,单独文件 `--gapbins-out`。每 (gap,1ms 格) 一行,
/// 下游据此在 1ms 网格精确拼 S 的 filler↔参考系数。
pub(crate) const GAPBINS_HEADER: &str = "gap_id,bin_index,t_lo,n_m,kind,left_bin,right_bin,tau";

/// 把一个 gap 的协方差块格式化成块表 CSV 行(spec §13)。变长字段(参考盒)用分号
/// 分隔;`ref_name` 把 references 切片内下标映射为盒名。cross-ref 段填 refs/k/标定
/// 计数、退化字段留空;degenerate 段反之。
pub(crate) fn gap_cov_row(
    gap_idx: usize,
    target_box: &str,
    has_cross_ref: bool,
    t_start: f64,
    t_stop: f64,
    cov: &GapCovBlock,
    ref_name: impl Fn(usize) -> String,
) -> String {
    let typ = if has_cross_ref { "crossref" } else { "degenerate" };
    let join = |f: &dyn Fn(&GapRefCalib) -> String| {
        cov.refs.iter().map(f).collect::<Vec<_>>().join(";")
    };
    let ref_boxes = join(&|r| ref_name(r.ref_idx));
    let ks = join(&|r| format!("{:.4}", r.k));
    let c_refs = join(&|r| format!("{:.0}", r.c_ref_cal));
    let c_a = cov.refs.first().map(|r| format!("{:.0}", r.c_a_cal)).unwrap_or_default();
    // rho 是 cross-ref 专有的重整因子(spec §7);退化段留空(与 c_a_cal 一致)。
    let rho = if has_cross_ref { format!("{:.4}", cov.rho) } else { String::new() };
    let opt_r = |o: Option<f64>| o.map(|v| format!("{v:.4}")).unwrap_or_default();
    let opt_n = |o: Option<f64>| o.map(|v| format!("{v:.0}")).unwrap_or_default();
    format!(
        "{},{},{},{:.6},{:.6},{},{},{},{},{},{},{},{},{},{},{},{:.4}",
        gap_idx,
        target_box,
        typ,
        t_start,
        t_stop,
        ref_boxes,
        ks,
        c_refs,
        c_a,
        rho,
        opt_r(cov.r_pre),
        opt_r(cov.r_post),
        opt_n(cov.n_pre),
        opt_n(cov.n_post),
        cov.maskable,
        !has_cross_ref,
        cov.sys_bias_scale,
    )
}

/// 把一个 gap 内 1ms 格格式化成 gapbins 表 CSV 行(spec ③)。measured 格记 n_m、
/// 端点/τ 留空;empty 格记端点 bin_index 与 τ、n_m 留空(下游用端点格的 n_m)。
pub(crate) fn gapbins_row(gap_idx: usize, b: &GapBinInfo) -> String {
    let (kind, n_m, left, right, tau) = match b.kind {
        GapBinKind::Measured => (
            "measured",
            b.n_m.to_string(),
            String::new(),
            String::new(),
            String::new(),
        ),
        GapBinKind::Empty => (
            "empty",
            String::new(),
            b.left_bin.map(|v| v.to_string()).unwrap_or_default(),
            b.right_bin.map(|v| v.to_string()).unwrap_or_default(),
            b.tau.map(|v| format!("{v:.4}")).unwrap_or_default(),
        ),
    };
    format!(
        "{},{},{:.6},{},{},{},{},{}",
        gap_idx, b.bin_index, b.t_lo, n_m, kind, left, right, tau,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// cross-ref 行:列出参考盒(ref_idx→盒名)、k、标定窗计数;退化字段留空。
    #[test]
    fn cross_ref_row_lists_refs_and_calib() {
        let cov = GapCovBlock {
            refs: vec![
                GapRefCalib { ref_idx: 0, k: 2.0, c_a_cal: 200.0, c_ref_cal: 100.0 },
                GapRefCalib { ref_idx: 1, k: 1.5, c_a_cal: 200.0, c_ref_cal: 80.0 },
            ],
            rho: 1.07,
            ..Default::default()
        };
        // references 切片内下标 0→"B"、1→"C"
        let names = ["A", "B", "C"];
        let row = gap_cov_row(3, "A", true, 1.0, 1.05, &cov, |i| names[i + 1].to_string());
        let f: Vec<&str> = row.split(',').collect();
        assert_eq!(f.len(), 17, "行应有 17 列(含 rho): {row}");
        assert_eq!(f[0], "3");
        assert_eq!(f[1], "A");
        assert_eq!(f[2], "crossref");
        assert_eq!(f[5], "B;C", "ref_boxes 应按 ref_idx 映射盒名");
        assert_eq!(f[6], "2.0000;1.5000");
        assert_eq!(f[7], "100;80");
        assert_eq!(f[8], "200");
        assert_eq!(f[9], "1.0700", "rho 列在 c_a_cal 之后");
        assert_eq!(f[10], "", "cross-ref 段 r_pre 应空");
        assert_eq!(f[14], "false");
        assert_eq!(f[15], "false", "cross-ref sys_bias_flag=false");
    }

    /// degenerate 行:端点率/事例数/系统偏,参考字段留空,sys_bias_flag=true。
    #[test]
    fn degenerate_row_has_rates_and_sys_flag() {
        let cov = GapCovBlock {
            r_pre: Some(400.0),
            r_post: Some(600.0),
            n_pre: Some(40.0),
            n_post: Some(60.0),
            maskable: false,
            sys_bias_scale: 0.2,
            ..Default::default()
        };
        let row = gap_cov_row(0, "C", false, 2.0, 2.1, &cov, |_| String::new());
        let f: Vec<&str> = row.split(',').collect();
        assert_eq!(f.len(), 17, "行应有 17 列(含 rho): {row}");
        assert_eq!(f[2], "degenerate");
        assert_eq!(f[5], "", "退化段 ref_boxes 应空");
        assert_eq!(f[8], "", "退化段 c_a_cal 应空");
        assert_eq!(f[9], "", "退化段 rho 应空");
        assert_eq!(f[10], "400.0000");
        assert_eq!(f[11], "600.0000");
        assert_eq!(f[12], "40");
        assert_eq!(f[13], "60");
        assert_eq!(f[14], "false");
        assert_eq!(f[15], "true", "degenerate sys_bias_flag=true");
        assert_eq!(f[16], "0.2000");
    }

    /// spec ③:gapbins 行——measured 格记 n_m、kind=measured、端点/τ 留空;
    /// empty 格记 kind=empty、端点 bin_index 与 τ、n_m 留空。
    #[test]
    fn gapbins_measured_and_empty_rows() {
        let measured = GapBinInfo {
            bin_index: 0, t_lo: 1.0, n_m: 2, kind: GapBinKind::Measured,
            left_bin: None, right_bin: None, tau: None,
        };
        let rm = gapbins_row(7, &measured);
        let fm: Vec<&str> = rm.split(',').collect();
        assert_eq!(fm.len(), 8, "gapbins 行应 8 列: {rm}");
        assert_eq!(fm[0], "7", "gap_id");
        assert_eq!(fm[1], "0", "bin_index");
        assert_eq!(fm[2], "1.000000", "t_lo");
        assert_eq!(fm[3], "2", "measured 记 n_m");
        assert_eq!(fm[4], "measured");
        assert_eq!(fm[5], "", "measured left_bin 留空");
        assert_eq!(fm[6], "", "measured right_bin 留空");
        assert_eq!(fm[7], "", "measured tau 留空");

        let empty = GapBinInfo {
            bin_index: 5, t_lo: 1.005, n_m: 0, kind: GapBinKind::Empty,
            left_bin: Some(3), right_bin: Some(9), tau: Some(0.25),
        };
        let re = gapbins_row(7, &empty);
        let fe: Vec<&str> = re.split(',').collect();
        assert_eq!(fe.len(), 8);
        assert_eq!(fe[1], "5");
        assert_eq!(fe[3], "", "empty n_m 留空(下游用端点格 n_m)");
        assert_eq!(fe[4], "empty");
        assert_eq!(fe[5], "3", "left_bin");
        assert_eq!(fe[6], "9", "right_bin");
        assert_eq!(fe[7], "0.2500", "tau");
    }
}
