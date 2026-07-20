use blink_hxmt_he::algorithms::saturation::{
    assign_gap_fill_channels, detect_fifo_reset_intervals, detect_unreliable_intervals,
    extract_packet_infos, reconstruct_gaps, reconstruct_met_channels,
    reconstruct_met_pulse_widths, reconstruct_met_times, reconstruct_with_wrap_tracking,
    solve_events, unwrap_channel, BoxReconstructionData,
};
use blink_core::types::MissionElapsedTime;
use blink_hxmt_he::io::level_1b::{get_eng_filenames, get_sci_filenames};
use blink_hxmt_he::io::level_1k::EventFile;
use blink_hxmt_he::types::HxmtHe;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write as IoWrite};

use crate::cli::ReportArgs;
use crate::util::{
    epoch_hour_of_met, json_escape, load_boxes, parse_met_or_utc, warn_if_window_crosses_hour,
};

pub fn cmd_report(args: &ReportArgs) -> std::io::Result<()> {
    let trigger_met = parse_met_or_utc(&args.trigger);
    let epoch = epoch_hour_of_met(trigger_met);
    warn_if_window_crosses_hour(trigger_met, args.before, args.after, epoch);

    let trigger_utc = MissionElapsedTime::<HxmtHe>::new(trigger_met).to_utc();
    let met_min = trigger_met - args.before;
    let met_max = trigger_met + args.after;

    eprintln!("Loading 1B files for {}...", epoch.format("%Y-%m-%dT%H"));
    let boxes = load_boxes(epoch);
    eprintln!(
        "  Found {} boxes: {:?}",
        boxes.len(),
        boxes.iter().map(|(n, _, _)| n.as_str()).collect::<Vec<_>>()
    );

    create_dir_all(&args.out)?;

    // ── 1B path discovery (for manifest) ──
    let sci_paths = get_sci_filenames(epoch);
    let eng_paths = get_eng_filenames(epoch);

    // ── Run reconstruction for all boxes (needed for events_obs + events_rec) ──
    eprintln!("Preparing reconstruction data...");
    let mut box_data: Vec<(String, BoxReconstructionData)> = Vec::new();
    for (box_name, sci, offset) in &boxes {
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

    // ── Cross-box gap-fill ──
    eprintln!("Reconstructing (FIFO reset gaps)...");
    let mut filled_per_box: Vec<(String, Vec<(f64, u16, u8)>, usize, usize)> = Vec::new();
    for i in 0..box_data.len() {
        let refs: Vec<&BoxReconstructionData> = box_data
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, (_, d))| d)
            .collect();
        let (gap_results, _ref_weights) = reconstruct_gaps(&box_data[i].1, &refs);
        let n_lost_total: usize = gap_results.iter().map(|r| r.n_lost).sum();
        let n_ref = gap_results.iter().filter(|r| r.has_cross_ref).count();
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
        filled_per_box.push((box_data[i].0.clone(), gap_events, n_lost_total, n_ref));
    }

    // ── 1K loading ──
    eprintln!("Loading 1K EventFile...");
    let evt_1k = EventFile::from_epoch(&epoch).expect("Failed to load 1K EventFile");
    let k1_times = evt_1k.times();
    let k1_dets = evt_1k.det_ids();
    let k1_channels = evt_1k.channels();
    let box_ranges: [(&str, u8, u8); 3] = [("A", 0, 5), ("B", 6, 11), ("C", 12, 17)];

    // ── 1B FITS path discovery for manifest ──
    let mut he_eng_files: Vec<(String, String)> = Vec::new();
    for (box_name, _, _) in &boxes {
        if let Some((_, path)) = eng_paths.iter().find(|(b, _)| b == box_name) {
            he_eng_files.push((box_name.clone(), path.clone()));
        }
    }

    // ── Write per-box CSVs + collect summary ──
    let mut summary: Vec<(String, u64, u64, u64, u64, u64)> = Vec::new();

    for (i, (box_name, sci, offset)) in boxes.iter().enumerate() {
        let box_dir = args.out.join(format!("box_{}", box_name.to_lowercase()));
        create_dir_all(&box_dir)?;

        // events_obs.csv (observed 1B events in window, full detail incl. det_id)
        let obs_path = box_dir.join("events_obs.csv");
        let mut w = BufWriter::new(File::create(&obs_path)?);
        writeln!(w, "met,channel,det_id,pkt_idx,evt_idx,is_second,pulse_width")?;
        let detailed = solve_events(sci, *offset, Some(met_min), Some(met_max));
        let mut n_obs = 0u64;
        for e in &detailed {
            writeln!(
                w, "{:.6},{},{},{},{},{},{}",
                e.met, e.channel, e.det_id, e.pkt_index, e.evt_index,
                if e.is_second { 1 } else { 0 }, e.raw_bytes[1],
            )?;
            if !e.is_second { n_obs += 1; }
        }
        w.flush()?;

        // events_rec.csv (gap-filled reconstructed events in window)
        let rec_path = box_dir.join("events_rec.csv");
        let mut w = BufWriter::new(File::create(&rec_path)?);
        writeln!(w, "met,channel,pulse_width")?;
        let mut n_rec = 0u64;
        for &(t, ch, pw) in &filled_per_box[i].1 {
            if t >= met_min && t <= met_max {
                writeln!(w, "{:.6},{},{}", t, unwrap_channel(ch), pw)?;
                n_rec += 1;
            }
        }
        w.flush()?;

        // events_1k.csv (1K pipeline events for this box's det range, in window).
        // det_id is normalized to box-local (0..5) for consistency with events_obs.csv.
        let k1_path = box_dir.join("events_1k.csv");
        let mut w = BufWriter::new(File::create(&k1_path)?);
        writeln!(w, "met,channel,det_id")?;
        let (_, d_lo, d_hi) = box_ranges.iter().find(|(n, _, _)| n == box_name).unwrap();
        let mut n_1k = 0u64;
        for j in 0..k1_times.len() {
            let (t, d, ch) = (k1_times[j], k1_dets[j], k1_channels[j]);
            if d >= *d_lo && d <= *d_hi && t >= met_min && t <= met_max {
                writeln!(w, "{:.6},{},{}", t, ch, d - d_lo)?;
                n_1k += 1;
            }
        }
        w.flush()?;

        // resets.csv (FIFO resets in window, with cluster_id)
        let resets_path = box_dir.join("resets.csv");
        let mut w = BufWriter::new(File::create(&resets_path)?);
        writeln!(w, "start_met,stop_met,gap_s,prev_pkt_idx,next_pkt_idx,n_lost,cluster_id")?;
        let packets = &box_data[i].1.packets;
        let mut in_window: Vec<(f64, f64, f64, usize, usize, usize)> = Vec::new();
        for iv in &box_data[i].1.gaps {
            if iv.stop_met < met_min || iv.start_met > met_max { continue; }
            let r_true = packets.iter()
                .find(|p| p.pkt_idx == iv.next_pkt_idx)
                .map(|p| 109.0 / p.span().max(1e-9))
                .unwrap_or(15797.0);
            let n_lost = (r_true * iv.gap_seconds).round() as usize;
            in_window.push((iv.start_met, iv.stop_met, iv.gap_seconds,
                            iv.prev_pkt_idx, iv.next_pkt_idx, n_lost));
        }
        // Cluster: resets within 1.0 s of each other share cluster_id.
        in_window.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mut cluster_id: usize = 0;
        let mut last_stop: f64 = f64::NEG_INFINITY;
        let mut n_lost_total = 0u64;
        for (s, e, g, p, nx, nl) in &in_window {
            if *s - last_stop > 1.0 { cluster_id += 1; }
            writeln!(w, "{:.6},{:.6},{:.6},{},{},{},{}",
                     s, e, g, p, nx, nl, cluster_id)?;
            last_stop = *e;
            n_lost_total += *nl as u64;
        }
        w.flush()?;
        let n_resets = in_window.len() as u64;

        eprintln!(
            "  Box {}: obs={} rec={} resets={} lost={} 1k={}",
            box_name, n_obs, n_rec, n_resets, n_lost_total, n_1k
        );
        summary.push((box_name.clone(), n_obs, n_rec, n_resets, n_lost_total, n_1k));
    }

    // ── manifest.json ──
    let manifest_path = args.out.join("manifest.json");
    let mut w = BufWriter::new(File::create(&manifest_path)?);
    writeln!(w, "{{")?;
    writeln!(w, "  \"trigger_raw\": \"{}\",", json_escape(&args.trigger))?;
    writeln!(w, "  \"trigger_met\": {:.6},", trigger_met)?;
    writeln!(w, "  \"trigger_utc\": \"{}\",", trigger_utc.format("%Y-%m-%dT%H:%M:%S%.3f"))?;
    writeln!(w, "  \"before_s\": {:.3},", args.before)?;
    writeln!(w, "  \"after_s\": {:.3},", args.after)?;
    writeln!(w, "  \"epoch\": \"{}\",", epoch.format("%Y-%m-%dT%H"))?;
    writeln!(w, "  \"met_min\": {:.6},", met_min)?;
    writeln!(w, "  \"met_max\": {:.6},", met_max)?;

    writeln!(w, "  \"boxes\": [{}],",
             boxes.iter().map(|(n, _, _)| format!("\"{}\"", n))
                 .collect::<Vec<_>>().join(", "))?;

    writeln!(w, "  \"level_1b_sci\": {{")?;
    for (i, (b, p)) in sci_paths.iter().enumerate() {
        let comma = if i + 1 < sci_paths.len() { "," } else { "" };
        writeln!(w, "    \"{}\": \"{}\"{}", b, json_escape(p), comma)?;
    }
    writeln!(w, "  }},")?;

    writeln!(w, "  \"level_1b_eng\": {{")?;
    for (i, (b, p)) in he_eng_files.iter().enumerate() {
        let comma = if i + 1 < he_eng_files.len() { "," } else { "" };
        writeln!(w, "    \"{}\": \"{}\"{}", b, json_escape(p), comma)?;
    }
    writeln!(w, "  }},")?;

    let evt_1k_path = blink_hxmt_he::io::path::get_path(&epoch, "Evt").ok();
    let orbit_1k_path = blink_hxmt_he::io::path::get_path(&epoch, "Orbit").ok();
    writeln!(w, "  \"level_1k_he_evt\": {},",
             evt_1k_path.as_deref().map(|p| format!("\"{}\"", json_escape(p)))
                 .unwrap_or_else(|| "null".to_string()))?;
    writeln!(w, "  \"level_1k_orbit\": {},",
             orbit_1k_path.as_deref().map(|p| format!("\"{}\"", json_escape(p)))
                 .unwrap_or_else(|| "null".to_string()))?;

    writeln!(w, "  \"summary\": {{")?;
    for (i, (b, n_obs, n_rec, n_resets, n_lost, n_1k)) in summary.iter().enumerate() {
        let comma = if i + 1 < summary.len() { "," } else { "" };
        writeln!(w,
            "    \"{}\": {{\"n_obs\": {}, \"n_rec\": {}, \"n_resets\": {}, \"n_lost\": {}, \"n_1k\": {}}}{}",
            b, n_obs, n_rec, n_resets, n_lost, n_1k, comma)?;
    }
    writeln!(w, "  }}")?;
    writeln!(w, "}}")?;
    w.flush()?;

    eprintln!("Wrote pack to {}", args.out.display());
    Ok(())
}
