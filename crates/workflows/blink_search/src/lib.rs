use blink_core::traits::{Chunk, Instrument};
use blink_workflow::process;
use chrono::prelude::*;
use indicatif::{MultiProgress, ProgressBar};
use std::collections::BTreeMap;
use std::fs;

mod report;

pub use report::{DayReport, HourRecord, HourStatus};

/// 搜索一天。产出两个文件，缺一不可：
///
/// * `YYYYMMDD_signals.json` —— 候选
/// * `YYYYMMDD_hours.json`   —— 逐小时账本：每个小时要么 searched、要么
///   excluded(reason)，外加曝光秒数。没有它，"这天没候选"分不清是真没有
///   还是根本没搜，率估计的分母就是错的。
pub fn search_day<I: Instrument>(day: NaiveDate, multi_progress: &MultiProgress) {
    let spin_bar = multi_progress.add(ProgressBar::new(24));
    spin_bar.set_style(
        indicatif::ProgressStyle::default_spinner()
            .template("{spinner} {msg}")
            .unwrap(),
    );

    spin_bar.set_message("ensure folder exist");
    let year = day.year();
    let month = day.month();
    let output_dir = format!(
        "data/{}/{:04}/{:02}/",
        I::name().replace("/", "_"),
        year,
        month
    );
    std::fs::create_dir_all(&output_dir).expect("failed to create output directory");
    let output_file = format!(
        "{}{:04}{:02}{:02}_signals.json",
        output_dir,
        year,
        month,
        day.day(),
    );
    let report_file = format!(
        "{}{:04}{:02}{:02}_hours.json",
        output_dir,
        year,
        month,
        day.day(),
    );

    spin_bar.set_message("check last modified");
    let last_modified = (0..24)
        .filter_map(|hour| {
            let naive = day.and_hms_opt(hour, 0, 0).expect("invalid time");
            let epoch = Utc.from_utc_datetime(&naive);
            I::Chunk::last_modified(&epoch).ok()
        })
        .max();
    if let Some(last_modified) = last_modified {
        // 两个产出都得存在且不比源数据旧才算这天做完了。少了逐小时账本
        // （老版本只写候选）就得重跑，否则这天永远没有曝光核算。
        let up_to_date = |path: &str| {
            fs::metadata(path)
                .and_then(|metadata| metadata.modified())
                .map(|modified| DateTime::<Utc>::from(modified) >= last_modified)
                .unwrap_or(false)
        };
        if up_to_date(&output_file) && up_to_date(&report_file) {
            return;
        }
    }

    spin_bar.finish_and_clear();

    let progress_bar = multi_progress.add(ProgressBar::new(24));
    progress_bar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] [{bar:40.yellow/red}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    let mut all_signals = Vec::new();
    let mut hours: Vec<HourRecord> = Vec::with_capacity(24);
    for hour in 0..24 {
        let naive = day.and_hms_opt(hour, 0, 0).expect("invalid time");
        let record = match I::Chunk::from_epoch(&Utc.from_utc_datetime(&naive)) {
            // 诊断量一律在这一小时的活干完之后再取：有些量（如搜索中因取不到
            // 星历而丢弃的候选数）是干活过程中才产生的。
            Ok(chunk) => {
                let metrics = |chunk: &I::Chunk| {
                    chunk
                        .diagnostics()
                        .into_iter()
                        .map(|(name, value)| (name.to_string(), value))
                        .collect::<BTreeMap<_, _>>()
                };
                // 先体检再核算：coverage 要重建 1B，别为一个不搜的小时白付。
                match chunk.exclusion() {
                    // 体检没过：不搜。搜了也是错的，产出假候选比没有更糟。
                    Some(reason) => {
                        HourRecord::excluded(hour, reason, None).with_metrics(metrics(&chunk))
                    }
                    None => {
                        let mut signals = chunk
                            .search()
                            .into_iter()
                            .map(|signal| signal.to_unified())
                            .collect::<Vec<_>>();
                        let n_signals = signals.len();
                        all_signals.append(&mut signals);
                        HourRecord::searched(hour, chunk.coverage(), n_signals)
                            .with_metrics(metrics(&chunk))
                    }
                }
            }
            // 载不进来同样是显式排除，不是"这小时没候选"。
            Err(error) => HourRecord::excluded(hour, (&error).into(), Some(error.to_string())),
        };
        hours.push(record);
        progress_bar.inc(1);
    }
    progress_bar.finish_and_clear(); // 使用 finish_and_clear() 以便完成后清除内层进度条

    let spin_bar_writting = multi_progress.add(ProgressBar::new_spinner());
    spin_bar_writting.set_style(
        indicatif::ProgressStyle::default_spinner()
            .template("{spinner} {msg}")
            .unwrap(),
    );
    spin_bar_writting.set_message("writing output files");

    let suffix = format!(".{}.tmp", nanoid::nanoid!(3));
    let write_atomic = |path: &str, contents: &str| {
        let temp = format!("{}{}", path, &suffix);
        std::fs::write(&temp, contents).expect("failed to write output file");
        std::fs::rename(&temp, path).expect("failed to rename output file");
    };

    let json = serde_json::to_string_pretty(&all_signals).expect("failed to serialize signals");
    write_atomic(&output_file, &json);

    spin_bar_writting.set_message("writing hour report");
    let report = DayReport::new(day, I::name(), hours);
    let json = serde_json::to_string_pretty(&report).expect("failed to serialize hour report");
    write_atomic(&report_file, &json);

    // 老版本的自由文本 errors.txt 已被逐小时账本取代；留着只会跟账本打架。
    let _ = fs::remove_file(format!(
        "{}{:04}{:02}{:02}_errors.txt",
        output_dir,
        year,
        month,
        day.day(),
    ));

    spin_bar_writting.finish_and_clear();
}

pub fn search_all<I: Instrument>(total_workers: usize, idx_worker: usize) {
    process::<I, _, _>(None, None, search_day::<I>, total_workers, idx_worker);
}

/// 在 [start, end] 闭区间日期范围内搜索（按天 round-robin 分片）。
///
/// 第 `idx_worker`/`total_workers` 个 worker 只处理 `day_offset % total_workers == idx_worker`
/// 的天。每天结果原子写入 `data/<I>/YYYY/MM/YYYYMMDD_signals.json`（temp + rename），
/// 并按源文件 last_modified 跳过已处理的天，因此可安全地并行、断点重跑。
pub fn search_range<I: Instrument>(
    start: NaiveDate,
    end: NaiveDate,
    total_workers: usize,
    idx_worker: usize,
) {
    process::<I, _, _>(
        Some(start),
        Some(end),
        search_day::<I>,
        total_workers,
        idx_worker,
    );
}
