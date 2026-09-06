use chrono::prelude::*;
use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

use crate::util::{epoch_hour_of_met, parse_met_or_utc};

#[derive(Parser)]
#[command(about = "HXMT HE analysis toolkit")]
pub struct Cli {
    #[command(subcommand)]
    pub command: TopCommands,
}

#[derive(Subcommand)]
pub enum TopCommands {
    /// Saturation analysis (detect FIFO resets, reconstruct gaps, generate reports)
    Sat {
        #[command(subcommand)]
        command: SatCommands,
    },
    /// TGF search (scan date range for candidate signals)
    Search {
        /// Start date (YYYY-MM-DD)
        from: String,
        /// End date (YYYY-MM-DD)
        to: String,
        /// Total number of parallel workers (days are sharded round-robin)
        #[arg(long, default_value_t = 1)]
        workers: usize,
        /// This worker's index in [0, workers)
        #[arg(long, default_value_t = 0)]
        worker: usize,
        /// Instrument to search
        #[arg(long, value_enum, default_value_t = Instrument::HxmtHe)]
        instrument: Instrument,
    },
    /// WWLLN lightning association + REP train-density enrichment for detected signals
    Wwlln {
        /// Instrument whose candidates to enrich
        #[arg(long, value_enum, default_value_t = Instrument::HxmtHe)]
        instrument: Instrument,
        /// Half-width of the time window around the peak within which a stroke counts as
        /// associated (ms). 5 ms is the TGF criterion; wider windows serve special tests
        /// such as electron beams that arrive tens of ms after the stroke.
        #[arg(long, default_value_t = 5)]
        window_ms: i64,
    },
    /// Recompute per-candidate ACD coincidence counts from 1K events (offline audit;
    /// needs 1K archive access). Input: CSV with `start`/`stop` columns (UTC or MET),
    /// e.g. sig_all_v5.csv as-is. Output: input columns + n,n_acd,n_acd_multi,n_bg,n_acd_bg
    AcdAudit {
        /// Input candidate list CSV
        list: PathBuf,
        /// Output CSV path
        #[arg(long, short = 'o')]
        out: PathBuf,
        /// Event selection: `csi` = search-side keep filter (CsI, ch>=38);
        /// `nai` = NaI non-Am241 events (electron stopping layer, positive control)
        #[arg(long, default_value = "csi")]
        scint: String,
    },
    /// Catalog stage: pool-level REP train removal, then the paper selection
    /// criteria (fa <= 1e-5, or fa <= 1 with lightning association).
    /// Input: tgfs.json from `blink wwlln`. Output: catalog CSV.
    Catalog {
        /// Enriched candidate list from `blink wwlln`
        #[arg(default_value = "tgfs.json")]
        input: PathBuf,
        /// Output catalog CSV path
        #[arg(long, short = 'o')]
        out: PathBuf,
    },
}

#[derive(Subcommand)]
pub enum SatCommands {
    /// Full diagnostic data pack for one burst (events, resets, summary)
    Report(ReportArgs),
    /// Detect FIFO resets in a burst window
    Detect(BurstArgs),
    /// Gap-filled light curve (1B + cross-box reconstruction)
    Reconstruct(ReconstructArgs),
    /// Mask-and-reconstruct injection validation (spec §11): inject fake gaps on a
    /// target box in unsaturated data, cross-ref reconstruct, dump truth vs fill
    Inject(InjectArgs),
    /// Per-event dump from 1B (raw) or 1K pipeline
    Extract(ExtractArgs),
    /// Compare 1B vs 1K event data
    Compare(CompareArgs),
    /// Scan a 1B hour for FIFO resets (no trigger; for offline sweeps)
    Scan(ScanArgs),
    /// Low-level diagnostic dumps
    Dump {
        #[command(subcommand)]
        sub: DumpCommands,
    },
}

/// Shared positional + flags for burst-centric subcommands.
/// EPOCH is derived from TRIGGER (1B archive is per-hour partitioned).
#[derive(Args)]
pub struct BurstWindow {
    /// Trigger time (MET number or UTC datetime, e.g. 2020-04-15T08:48:05.560)
    pub trigger: String,
    /// Seconds before trigger
    #[arg(long)]
    pub before: f64,
    /// Seconds after trigger
    #[arg(long)]
    pub after: f64,
    /// Filter to a single box (a, b, or c). If omitted, all boxes.
    #[arg(long = "box")]
    pub box_filter: Option<String>,
}

impl BurstWindow {
    pub fn trigger_met(&self) -> f64 {
        parse_met_or_utc(&self.trigger)
    }
    pub fn met_min(&self) -> f64 {
        self.trigger_met() - self.before
    }
    pub fn met_max(&self) -> f64 {
        self.trigger_met() + self.after
    }
    pub fn epoch(&self) -> DateTime<Utc> {
        epoch_hour_of_met(self.trigger_met())
    }
}

#[derive(Args)]
pub struct BurstArgs {
    #[command(flatten)]
    pub window: BurstWindow,
}

#[derive(Args)]
pub struct ReportArgs {
    /// Trigger time (MET number or UTC datetime)
    pub trigger: String,
    /// Seconds before trigger
    #[arg(long)]
    pub before: f64,
    /// Seconds after trigger
    #[arg(long)]
    pub after: f64,
    /// Output directory for the data pack
    #[arg(long, short = 'o')]
    pub out: PathBuf,
}

#[derive(Args)]
pub struct ReconstructArgs {
    #[command(flatten)]
    pub window: BurstWindow,
    /// Bin width in seconds
    #[arg(long, default_value_t = 1.0)]
    pub bin: f64,
    /// Optional: write per-gap covariance block table (spec §13) to this file
    #[arg(long)]
    pub gapcov_out: Option<std::path::PathBuf>,
    /// Optional: write per-gap 1ms bin structure table (spec ③ gapbins) to this file
    #[arg(long)]
    pub gapbins_out: Option<std::path::PathBuf>,
}

#[derive(Args)]
pub struct InjectArgs {
    #[command(flatten)]
    pub window: BurstWindow,
    /// Target box to inject fake gaps on (a, b, or c)
    #[arg(long)]
    pub target: String,
    /// Fake-gap centers as second offsets from trigger (comma-separated)
    #[arg(long, value_delimiter = ',')]
    pub at: Vec<f64>,
    /// Width of each injected gap in seconds
    #[arg(long, default_value_t = 0.03)]
    pub width: f64,
    /// Co-saturation sub-interval width (seconds) centered in each gap: reference
    /// boxes are marked unreliable there (simulating them also saturating), which
    /// produces genuine empty cells (co-saturation, not a Poisson void). 0 = off.
    #[arg(long, default_value_t = 0.0)]
    pub cosat_width: f64,
    /// Optional: write the reconstructed event stream (spec ①) to this file.
    /// Reference boxes contribute EVT rows (source counts C); the target box's
    /// in-gap events are masked out (they are the withheld truth) and replaced
    /// by FILL_GAP filler rows.
    #[arg(long)]
    pub events_out: Option<std::path::PathBuf>,
    /// Optional: write per-gap covariance block table (spec §13) to this file.
    #[arg(long)]
    pub gapcov_out: Option<std::path::PathBuf>,
    /// Optional: write per-gap 1ms bin structure table (spec ③) to this file.
    #[arg(long)]
    pub gapbins_out: Option<std::path::PathBuf>,
    /// Optional: write per-gap truth vs fill count summary to this file.
    #[arg(long)]
    pub truth_out: Option<std::path::PathBuf>,
}

#[derive(Args)]
pub struct ExtractArgs {
    #[command(flatten)]
    pub window: BurstWindow,
    /// Source: 1b (raw with MET reconstruction) or 1k (pipeline)
    #[arg(long, default_value = "1b")]
    pub source: String,
}

#[derive(Args)]
pub struct CompareArgs {
    #[command(flatten)]
    pub window: BurstWindow,
    /// Coarse bin width in seconds
    #[arg(long, default_value_t = 1.0)]
    pub coarse_bin: f64,
    /// Fine bin width in seconds
    #[arg(long, default_value_t = 0.1)]
    pub fine_bin: f64,
    /// Max lag in ms for cross-correlation
    #[arg(long, default_value_t = 50)]
    pub max_lag: usize,
    /// Threshold percentage for flagging fine bins
    #[arg(long, default_value_t = 30.0)]
    pub threshold: f64,
    /// Output CSV format
    #[arg(long)]
    pub csv: bool,
}

#[derive(Args)]
pub struct ScanArgs {
    /// Epoch in YYYY-MM-DDTHH format
    #[arg(long)]
    pub epoch: String,
    /// Filter to a single box (a, b, or c). If omitted, all boxes.
    #[arg(long = "box")]
    pub box_filter: Option<String>,
}

#[derive(Subcommand)]
pub enum DumpCommands {
    /// Dump event MET times
    Times(DumpBurstArgs),
    /// Dump packet time ranges
    Packets(DumpBurstArgs),
    /// Dump event details
    Events(DumpBurstArgs),
    /// Histogram of events
    Hist(DumpHistArgs),
    /// Per-packet diagnostics
    Diag(DumpBurstArgs),
    /// Dump ptime/UTC mapping for a packet range
    Ptime(DumpRangeArgs),
    /// Check byte offsets for CRC for a packet range
    CheckOffset(DumpRangeArgs),
}

#[derive(Args)]
pub struct DumpBurstArgs {
    /// Epoch in YYYY-MM-DDTHH format
    #[arg(long)]
    pub epoch: String,
    /// Trigger time (MET number or UTC datetime)
    pub trigger: String,
    /// Seconds before trigger
    #[arg(long, default_value_t = 10.0)]
    pub before: f64,
    /// Seconds after trigger
    #[arg(long, default_value_t = 100.0)]
    pub after: f64,
    /// Filter to a single box (a, b, or c). If omitted, all boxes.
    #[arg(long = "box")]
    pub box_filter: Option<String>,
}

impl DumpBurstArgs {
    pub fn trigger_met(&self) -> f64 { parse_met_or_utc(&self.trigger) }
    pub fn met_min(&self) -> f64 { self.trigger_met() - self.before }
    pub fn met_max(&self) -> f64 { self.trigger_met() + self.after }
}

#[derive(Args)]
pub struct DumpHistArgs {
    #[command(flatten)]
    pub window: DumpBurstArgs,
    /// Bin width in seconds
    #[arg(long, default_value_t = 0.01)]
    pub bin: f64,
}

#[derive(Args)]
pub struct DumpRangeArgs {
    /// Epoch in YYYY-MM-DDTHH format
    #[arg(long)]
    pub epoch: String,
    /// Minimum packet index
    pub pkt_min: usize,
    /// Maximum packet index
    pub pkt_max: usize,
    /// Filter to a single box (a, b, or c). If omitted, all boxes.
    #[arg(long = "box")]
    pub box_filter: Option<String>,
}

/// 可搜索的仪器。搜索管线本身与仪器无关（见 `blink_search::search_day`），
/// 每颗星只提供自己的 `Chunk`/`Event` 实现。
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Instrument {
    /// Insight-HXMT/HE（1B 重建 + 饱和掩模）
    HxmtHe,
    /// SVOM/GRM（L1B 事例 + GTI 曝光）
    SvomGrm,
    /// Fermi/GBM（continuous TTE，NaI 与 BGO 分组）
    FermiGbm,
    /// 天格 GRID-02（逐过境事例，2020-11 .. 2021-03）
    Grid02,
    /// 天格 GRID-03B（2022-03 .. 2024-08）
    Grid03b,
    /// 天格 GRID-04（2022-03 .. 2024-08）
    Grid04,
    /// 天格 GRID-07（2024-01 .. 2024-07）
    Grid07,
}
