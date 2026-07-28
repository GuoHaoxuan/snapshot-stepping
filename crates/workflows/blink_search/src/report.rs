//! 逐小时搜索账本。
//!
//! 搜索的原子单元是**小时**，不是天。所以账本也按小时记：一天 24 条，每条
//! 要么 `searched`（附曝光秒数、候选数），要么 `excluded`（附原因）。这样
//! 一天的产出完全自解释，不存在"没搜也没记"的静默缺口 —— 那种缺口会让
//! 事件率、误报率的分母悄悄偏大。

use blink_core::types::{Coverage, ExclusionReason};
use chrono::NaiveDate;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum HourStatus {
    Searched,
    Excluded,
}

/// 一小时的账目。
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct HourRecord {
    pub hour: u32,
    pub status: HourStatus,
    /// 被排除的原因，`searched` 时为空。
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reason: Option<ExclusionReason>,
    /// 排除原因的细节（如载入失败的错误原文）。
    #[serde(skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
    /// 真正进入搜索的时间（秒）。被排除的小时恒为 0。
    pub searched_seconds: f64,
    /// 曝光时间明细。只有真正搜了的小时才有 —— 被排除的小时曝光就是 0，
    /// 而算这份明细要重建 1B，不值当。
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coverage: Option<Coverage>,
    /// 体检判据的实测值（如 `config_gap_fraction`、`n_events`）。
    /// 判没判都记，便于事后审计阈值。
    #[serde(skip_serializing_if = "BTreeMap::is_empty", default)]
    pub metrics: BTreeMap<String, f64>,
    pub n_signals: usize,
}

impl HourRecord {
    pub fn searched(hour: u32, coverage: Coverage, n_signals: usize) -> Self {
        Self {
            hour,
            status: HourStatus::Searched,
            reason: None,
            detail: None,
            searched_seconds: coverage.searched_seconds(),
            coverage: Some(coverage),
            metrics: BTreeMap::new(),
            n_signals,
        }
    }

    pub fn excluded(hour: u32, reason: ExclusionReason, detail: Option<String>) -> Self {
        Self {
            hour,
            status: HourStatus::Excluded,
            reason: Some(reason),
            detail,
            searched_seconds: 0.0,
            coverage: None,
            metrics: BTreeMap::new(),
            n_signals: 0,
        }
    }

    pub fn with_metrics(mut self, metrics: BTreeMap<String, f64>) -> Self {
        self.metrics = metrics;
        self
    }
}

/// 一天的账本。
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct DayReport {
    pub date: NaiveDate,
    pub instrument: String,
    pub searched_hours: usize,
    pub excluded_hours: usize,
    /// 这天真正进入搜索的时间（秒）—— 率估计要用的分母。
    pub searched_seconds: f64,
    /// 逐原因的排除小时数。
    pub excluded_by_reason: BTreeMap<String, usize>,
    pub n_signals: usize,
    pub hours: Vec<HourRecord>,
}

impl DayReport {
    pub fn new(date: NaiveDate, instrument: &str, hours: Vec<HourRecord>) -> Self {
        let searched_hours = hours
            .iter()
            .filter(|record| record.status == HourStatus::Searched)
            .count();
        let mut excluded_by_reason: BTreeMap<String, usize> = BTreeMap::new();
        for reason in hours.iter().filter_map(|record| record.reason) {
            *excluded_by_reason.entry(reason.to_string()).or_default() += 1;
        }

        Self {
            date,
            instrument: instrument.to_string(),
            searched_hours,
            excluded_hours: hours.len() - searched_hours,
            searched_seconds: hours.iter().map(|record| record.searched_seconds).sum(),
            excluded_by_reason,
            n_signals: hours.iter().map(|record| record.n_signals).sum(),
            hours,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn coverage(masked: f64) -> Coverage {
        Coverage {
            span_seconds: 3600.0,
            masked_seconds: masked,
        }
    }

    #[test]
    fn every_hour_is_accounted_for() {
        let hours = vec![
            HourRecord::searched(0, coverage(12.0), 2),
            HourRecord::excluded(1, ExclusionReason::NonstandardConfig, None),
            HourRecord::excluded(2, ExclusionReason::MissingData, Some("no file".into())),
        ];
        let report = DayReport::new(
            NaiveDate::from_ymd_opt(2019, 12, 6).unwrap(),
            "HXMT/HE",
            hours,
        );

        assert_eq!(report.searched_hours, 1);
        assert_eq!(report.excluded_hours, 2);
        assert_eq!(report.searched_seconds, 3588.0);
        assert_eq!(report.n_signals, 2);
        assert_eq!(report.excluded_by_reason["nonstandard_config"], 1);
        assert_eq!(report.excluded_by_reason["missing_data"], 1);
        // 排除的小时不能贡献曝光
        assert!(
            report
                .hours
                .iter()
                .filter(|record| record.status == HourStatus::Excluded)
                .all(|record| record.searched_seconds == 0.0)
        );
    }
}
