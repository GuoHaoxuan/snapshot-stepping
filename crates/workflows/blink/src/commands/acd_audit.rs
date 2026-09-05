//! 离线复算候选的 ACD 符合计数，用于在认证样本上验证 REP/TGF 判别量。
//!
//! 输入任意带 `start`/`stop` 列（UTC ISO 或 MET 数字）的 CSV——
//! `sig_all_v5.csv` 原样可用；输出保留全部输入列，追加
//! `n,n_acd,n_acd_multi,n_bg,n_acd_bg` 五列原始计数，以及 `det_share_max`——
//! 窗内贡献最多的那个探头占窗内事例数的比例（单路毛刺判据的审计量，
//! 这里按候选的 start/stop 整窗算，搜索里按最显著一格算）。
//! 事例选择与窗口定义复用搜索侧同一实现（`blink_hxmt_he::algorithms::acd`），
//! 两边数字可直接互校。需在能访问 1K 档案的机器上运行。

use blink_algorithms::detector_share::max_detector_fraction;
use blink_core::traits::Event as _;
use blink_core::types::MissionElapsedTime;
use blink_hxmt_he::algorithms::acd::acd_counts;
use blink_hxmt_he::io::level_1k::EventFile;
use blink_hxmt_he::types::{Event, HxmtHe};
use chrono::prelude::*;
use std::collections::BTreeMap;
use std::path::Path;

/// 与 util::parse_met_or_utc 同语义，但不逐行打印（审计表几千行）。
fn parse_time(s: &str) -> Option<f64> {
    let s = s.trim();
    if let Ok(met) = s.parse::<f64>() {
        return Some(met);
    }
    let utc = s
        .parse::<DateTime<Utc>>()
        .or_else(|_| format!("{s}Z").parse::<DateTime<Utc>>())
        .ok()?;
    Some(MissionElapsedTime::<HxmtHe>::from(utc).met())
}

fn epoch_hour(met: f64) -> DateTime<Utc> {
    let utc = MissionElapsedTime::<HxmtHe>::new(met).to_utc();
    utc.date_naive()
        .and_hms_opt(utc.hour(), 0, 0)
        .expect("valid hour")
        .and_utc()
}

pub fn cmd_acd_audit(list: &Path, out: &Path, scint: &str) {
    let select: fn(&Event) -> bool = match scint {
        "csi" => |e| e.keep(),
        // NaI 是电子截止层：REP 电子若真穿过 ACD，应在这里留下 ACD 符合升高
        "nai" => {
            |e| e.detector.scintillator == blink_hxmt_he::types::Scintillator::Nai && !e.is_am241
        }
        other => panic!("--scint must be `csi` or `nai`, got `{other}`"),
    };
    let content = std::fs::read_to_string(list).expect("failed to read input csv");
    let mut lines = content.lines();
    let header = lines.next().expect("empty input csv");
    let columns: Vec<&str> = header.split(',').collect();
    let find = |name: &str| {
        columns
            .iter()
            .position(|c| c.trim() == name)
            .unwrap_or_else(|| panic!("input csv has no `{name}` column"))
    };
    let (i_start, i_stop) = (find("start"), find("stop"));

    struct Row<'a> {
        line: &'a str,
        window: Option<(f64, f64)>,
        counts: Option<blink_core::types::AcdCounts>,
        det_share_max: Option<f64>,
    }
    let mut rows: Vec<Row> = lines
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split(',').collect();
            let window = match (fields.get(i_start), fields.get(i_stop)) {
                (Some(a), Some(b)) => parse_time(a).zip(parse_time(b)),
                _ => None,
            };
            Row {
                line,
                window,
                counts: None,
                det_share_max: None,
            }
        })
        .collect();

    // 按候选所在小时分组，每个 1K 小时文件只载入一次
    let mut by_hour: BTreeMap<DateTime<Utc>, Vec<usize>> = BTreeMap::new();
    for (index, row) in rows.iter().enumerate() {
        if let Some((start, _)) = row.window {
            by_hour.entry(epoch_hour(start)).or_default().push(index);
        }
    }

    let n_hours = by_hour.len();
    let mut n_missing_hours = 0usize;
    for (i, (hour, indices)) in by_hour.into_iter().enumerate() {
        let event_file = match EventFile::from_epoch(&hour) {
            Ok(file) => file,
            Err(error) => {
                n_missing_hours += 1;
                eprintln!(
                    "acd-audit: skip hour {} ({} candidates): {error}",
                    hour.format("%Y-%m-%dT%H"),
                    indices.len()
                );
                continue;
            }
        };
        // 与搜索同一事例选择；1K 表按时间有序，防御性排序兜底
        let mut events: Vec<Event> = event_file.into_iter().filter(|e| select(e)).collect();
        events.sort_by(|a, b| a.time().cmp(&b.time()));
        for index in indices {
            let (start, stop) = rows[index].window.expect("grouped rows have windows");
            rows[index].counts = Some(acd_counts(&events, start, stop));
            rows[index].det_share_max = Some(max_detector_fraction(
                &events,
                MissionElapsedTime::new(start),
                MissionElapsedTime::new(stop),
                |e| e.detector.id,
            ));
        }
        if (i + 1) % 50 == 0 || i + 1 == n_hours {
            eprintln!("acd-audit: {}/{n_hours} hours", i + 1);
        }
    }

    let mut output = String::with_capacity(content.len() * 2);
    output.push_str(header);
    output.push_str(",n,n_acd,n_acd_multi,n_bg,n_acd_bg,det_share_max\n");
    let mut n_unresolved = 0usize;
    for row in &rows {
        output.push_str(row.line);
        match &row.counts {
            Some(c) => output.push_str(&format!(
                ",{},{},{},{},{},{:.3}\n",
                c.n,
                c.n_acd,
                c.n_acd_multi,
                c.n_bg,
                c.n_acd_bg,
                row.det_share_max.unwrap_or(f64::NAN)
            )),
            None => {
                n_unresolved += 1;
                output.push_str(",,,,,,\n");
            }
        }
    }
    std::fs::write(out, output).expect("failed to write output csv");
    eprintln!(
        "acd-audit: {} rows audited, {} unresolved ({} missing hours) -> {}",
        rows.len() - n_unresolved,
        n_unresolved,
        n_missing_hours,
        out.display()
    );
}
