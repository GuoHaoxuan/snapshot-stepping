use blink_hxmt_he::algorithms::saturation::{
    assign_gap_fill_channels, detect_fifo_reset_intervals, detect_unreliable_intervals,
    extract_packet_infos, reconstruct_gaps, reconstruct_met_channels,
    reconstruct_met_pulse_widths, reconstruct_met_times, reconstruct_with_wrap_tracking,
    unwrap_channel, BoxReconstructionData, CHANNEL_SEC,
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
    let mut all_filled: Vec<(String, Vec<(f64, u16, u8, f64)>)> = Vec::new();
    // 逐盒逐观测事件的粒子权重（particle weight），初值 1.0（自己那一份）；
    // 被别的盒 gap 当参考时累加它的贡献 ρ·k/n_m·(1+Λ)。
    let mut weights: Vec<Vec<f64>> =
        box_data.iter().map(|(_, d)| vec![1.0f64; d.events.len()]).collect();

    for i in 0..box_data.len() {
        let refs_with_idx: Vec<(usize, &BoxReconstructionData)> = box_data
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(j, (_, d))| (j, d))
            .collect();
        let refs: Vec<&BoxReconstructionData> =
            refs_with_idx.iter().map(|(_, d)| *d).collect();

        let (gap_results, ref_weights) = reconstruct_gaps(&box_data[i].1, &refs);
        // 参考贡献累加回各参考盒观测事件（退化 gap 的方差改挂 filler、不碰观测事件）
        for (r, (gj, _)) in refs_with_idx.iter().enumerate() {
            for (ev, &c) in ref_weights[r].iter().enumerate() {
                weights[*gj][ev] += c;
            }
        }
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
        let mut gap_events: Vec<(f64, u16, u8, f64)> = gap_results
            .iter()
            .zip(banded.iter())
            .flat_map(|(r, b)| {
                let fw = r.filler_weight; // 退化 gap 的每-filler 仅方差权重（cross-ref=0）
                r.filled_events
                    .iter()
                    .copied()
                    .zip(b.channels.iter().copied())
                    .zip(b.pulse_widths.iter().copied())
                    .map(move |((t, c), w)| (t, c, w, fw))
            })
            .collect();
        gap_events.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        eprintln!(
            "  Box {}: gaps={} ({} evt, {} ref)",
            box_data[i].0, box_data[i].1.gaps.len(), n_gap_filled, n_gap_ref,
        );

        all_filled.push((box_data[i].0.clone(), gap_events));
    }

    println!("box,type,met,channel,pulse_width,pkt_idx,evt_idx,particle_weight");
    for (box_idx, (box_name, _data)) in box_data.iter().enumerate() {
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
                println!(
                    "{},EVT,{:.6},{},{},-1,-1,{:.6}",
                    box_name, t, raw, obs_pw[idx], weights[box_idx][idx]
                );
                n_obs += 1;
            }
        }
        for &(t, ch, pw, fw) in gap_events {
            if t >= met_min && t <= met_max {
                // filler 是假粒子：cross-ref gap 的 fw=0（方差在参考事件）；
                // 退化 gap 的 fw>0（该 gap 的 σ² 均摊到自己 filler 上，Σfw²=σ²_gap）
                println!(
                    "{},FILL_GAP,{:.6},{},{},-1,-1,{:.6}",
                    box_name, t, unwrap_channel(ch), pw, fw
                );
                n_gap += 1;
            }
        }

        eprintln!(
            "  Box {}: {} observed, {} gap-filled, bin={:.3}s",
            box_name, n_obs, n_gap, args.bin,
        );
    }
}
