use super::rec_sci_data::reconstruct_with_wrap_tracking;
use crate::io::level_1b::SciFile;

/// 饱和类型
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SaturationType {
    /// 整包丢失：FIFOAFullReset 触发，整个 FIFO A 被清空。
    /// 表现为重建时间序列中出现远大于正常包间隔的空洞。
    FifoReset,
}

/// 单个饱和区间
#[derive(Debug, Clone)]
pub struct SaturationInterval {
    /// 空洞起始 MET（前一个包的最晚事件时间）
    pub start_met: f64,
    /// 空洞结束 MET（后一个包的最早事件时间）
    pub stop_met: f64,
    /// 空洞持续时间（秒）
    pub gap_seconds: f64,
    /// 空洞前一个包的索引（排序后）
    pub prev_pkt_idx: usize,
    /// 空洞后一个包的索引（排序后）
    pub next_pkt_idx: usize,
    /// 饱和类型
    pub saturation_type: SaturationType,
}

/// 每个 CCSDS 包的时间摘要
struct PacketTimeSummary {
    /// 原始包索引（在 SciFile.ccsds 中的位置）
    pkt_idx: usize,
    /// 包内最早事件的重建 MET
    min_met: f64,
    /// 包内最晚事件的重建 MET
    max_met: f64,
    /// 包内有效事件数
    n_events: usize,
}

const GAP_FACTOR: f64 = 100.0;

/// FIFO Reset gap 的最大持续时间（秒）。超过此值的 gap 不认为是 FIFO 复位，
/// 而是数据传输中断、SAA 等其他原因。正常 FIFO reset gap 在 8ms~100ms 量级。
const MAX_FIFO_RESET_GAP: f64 = 1.0;

/// MCU 读取速率下限 (events/s)。
/// MCU 以固定速率从 FIFO A 读取：109 events / ~7ms ≈ 15,600 evt/s。
/// 只有当物理事件率超过此值时，FIFO 才可能溢出触发 FIFOAFullReset。
/// 设为 15000 略低于理论值，留一点余量。
const MCU_READ_RATE_FLOOR: f64 = 15000.0;

/// 从单个包的时间跨度和事件数估算平均事件间隔 (秒/事件)。
/// 如果包内时间跨度过小（<1μs）或事件数不足，返回 None。
fn mean_event_interval(summary: &PacketTimeSummary) -> Option<f64> {
    let span = summary.max_met - summary.min_met;
    if span < 1e-6 || summary.n_events < 2 {
        return None;
    }
    Some(span / summary.n_events as f64)
}

/// 从单个包估算事件率 (events/s)。
/// 如果包内时间跨度过小（<1μs）或事件数不足，返回 None。
fn event_rate(summary: &PacketTimeSummary) -> Option<f64> {
    let span = summary.max_met - summary.min_met;
    if span < 1e-6 || summary.n_events < 2 {
        return None;
    }
    Some(summary.n_events as f64 / span)
}

/// 检测整包丢失（FIFO reset）造成的饱和区间。
///
/// 算法：
/// 1. 对每个 CCSDS 包重建所有事件的 MET 时间，提取 (min_met, max_met, n_events)
/// 2. 按 min_met 排序
/// 3. 对每对相邻包：
///    - baseline = 紧邻两包中平均事件间隔较小的那个（事件率较高的包）
///    - local_max_rate = ±5 包窗口内（包含紧邻 2 包，共最多 12 包）的最大事件率
///    - 若 local_max_rate < MCU_READ_RATE_FLOOR → 跳过（源率不到饱和阈值）
///    - 若 gap > baseline × GAP_FACTOR → 标记为 FifoReset
pub fn detect_fifo_reset_intervals(sci_data: &SciFile, offset: f64) -> Vec<SaturationInterval> {
    let packet_times = reconstruct_with_wrap_tracking(sci_data, offset);

    let mut summaries: Vec<PacketTimeSummary> = Vec::new();
    for (pkt_idx, times) in packet_times.iter().enumerate() {
        let valid: Vec<f64> = times.iter().copied().filter(|t| !t.is_nan()).collect();
        if valid.is_empty() {
            continue;
        }
        let min_met = valid.iter().cloned().reduce(f64::min).unwrap();
        let max_met = valid.iter().cloned().reduce(f64::max).unwrap();
        summaries.push(PacketTimeSummary {
            pkt_idx,
            min_met,
            max_met,
            n_events: valid.len(),
        });
    }

    summaries.sort_by(|a, b| a.min_met.partial_cmp(&b.min_met).unwrap());

    let mut intervals = Vec::new();
    for (wi, window) in summaries.windows(2).enumerate() {
        let gap = window[1].min_met - window[0].max_met;
        if gap <= 0.0 {
            continue;
        }

        let iv_prev = mean_event_interval(&window[0]);
        let iv_next = mean_event_interval(&window[1]);
        let baseline = match (iv_prev, iv_next) {
            (Some(a), Some(b)) => a.min(b),
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => continue,
        };

        // 用 ±5 包窗口（含紧邻 2 包）的最大事件率作为本地源率估计：
        // 单包率有涨落，扩展到 12 包窗口取最大值更稳健。
        let lo = wi.saturating_sub(5);
        let hi = (wi + 6).min(summaries.len() - 1);
        let mut local_max_rate = 0.0_f64;
        let mut found = false;
        for k in lo..=hi {
            if let Some(r) = event_rate(&summaries[k]) {
                local_max_rate = local_max_rate.max(r);
                found = true;
            }
        }
        if !found || local_max_rate < MCU_READ_RATE_FLOOR {
            continue;
        }

        if gap > baseline * GAP_FACTOR && gap <= MAX_FIFO_RESET_GAP {
            intervals.push(SaturationInterval {
                start_met: window[0].max_met,
                stop_met: window[1].min_met,
                gap_seconds: gap,
                prev_pkt_idx: window[0].pkt_idx,
                next_pkt_idx: window[1].pkt_idx,
                saturation_type: SaturationType::FifoReset,
            });
        }
    }

    intervals
}

/// 每个包的时间摘要（公开版本，用于重建）
#[derive(Debug, Clone)]
pub struct PacketInfo {
    pub pkt_idx: usize,
    pub min_met: f64,
    pub max_met: f64,
    pub n_events: usize,
}

impl PacketInfo {
    pub fn span(&self) -> f64 {
        self.max_met - self.min_met
    }
}

/// 从 SciFile 提取包时间摘要列表（按 min_met 排序）
pub fn extract_packet_infos(sci_data: &SciFile, offset: f64) -> Vec<PacketInfo> {
    let packet_times = reconstruct_with_wrap_tracking(sci_data, offset);
    let mut infos: Vec<PacketInfo> = packet_times
        .iter()
        .enumerate()
        .filter_map(|(pkt_idx, times)| {
            let valid: Vec<f64> = times.iter().copied().filter(|t| !t.is_nan()).collect();
            if valid.is_empty() {
                return None;
            }
            let min_met = valid.iter().cloned().reduce(f64::min).unwrap();
            let max_met = valid.iter().cloned().reduce(f64::max).unwrap();
            Some(PacketInfo {
                pkt_idx,
                min_met,
                max_met,
                n_events: valid.len(),
            })
        })
        .collect();
    infos.sort_by(|a, b| a.min_met.partial_cmp(&b.min_met).unwrap());
    infos
}

/// 不可信时间区间（FIFO reset gap 或拥塞宽包的时间覆盖）
#[derive(Debug, Clone)]
pub struct UnreliableInterval {
    pub start: f64,
    pub stop: f64,
}

/// 单个 box 的饱和重建数据
#[derive(Debug)]
pub struct BoxReconstructionData {
    /// 原始事件 MET 时间（已排序）
    pub events: Vec<f64>,
    /// 与 events 一一对应的 wrapped channel（SEC 槽 = CHANNEL_SEC）
    pub channels: Vec<u16>,
    /// 与 events 一一对应的脉宽 pulinfo（NaI/CsI 甄别；SEC 槽 = 0）
    pub pulse_widths: Vec<u8>,
    /// FIFO reset 区间
    pub gaps: Vec<SaturationInterval>,
    /// 包信息
    pub packets: Vec<PacketInfo>,
    /// 每个包内的事件时间（索引 = 原始包号，内部已排序）
    pub packet_events: Vec<Vec<f64>>,
    /// 不可信区间（FIFO reset gap），用于交叉参考时排除
    pub unreliable: Vec<UnreliableInterval>,
}

/// 检测不可信时间区间：仅 FIFO reset gap。
///
/// 拥塞宽包和包内泊松异常检测已移除：
/// - 拥塞宽包：实际触发多为 SAA 开关机导致的包跨时异常，非真正拥塞
/// - 包内异常：与静默丢数检测同一判据（泊松 log₁₀(p) < -10），
///   因 λ 在单包时间跨度内不稳定导致大量误报
pub fn detect_unreliable_intervals(
    gaps: &[SaturationInterval],
    _packets: &[PacketInfo],
    _packet_events: &[Vec<f64>],
) -> Vec<UnreliableInterval> {
    let mut intervals: Vec<UnreliableInterval> = Vec::new();

    for g in gaps {
        intervals.push(UnreliableInterval {
            start: g.start_met,
            stop: g.stop_met,
        });
    }

    // 按 start 排序
    intervals.sort_by(|a, b| a.start.partial_cmp(&b.start).unwrap());
    intervals
}

fn is_in_unreliable(t: f64, intervals: &[UnreliableInterval]) -> bool {
    intervals.iter().any(|iv| t >= iv.start && t <= iv.stop)
}

/// 单个参考盒对某 cross-ref gap 的标定描述子(spec §13/§5b)。
/// σ_k 由计数给:(σ_k/k)² = 1/C_a_cal + 1/C_ref_cal;跨-ref 相关用共用分子
/// C_a_cal(评审①):Σ_k[b,b'] = k_b k_b'(δ/C_ref_cal + 1/C_a_cal)。
#[derive(Debug, Clone)]
pub struct GapRefCalib {
    /// 在 `references` 切片里的下标
    pub ref_idx: usize,
    /// 标定系数 k = 目标/参考 率比(标定窗内)
    pub k: f64,
    /// 标定窗内 target 计数(所有 ref 共用分子)
    pub c_a_cal: f64,
    /// 标定窗内该参考盒计数
    pub c_ref_cal: f64,
}

/// gap 的协方差块描述子(spec §13)。cross-ref 用 `refs`(每参考盒的 k 与标定窗计数);
/// degenerate 用端点率 `r_pre/r_post` 与实际事例数 `n_pre/n_post`(自由度)。
/// `sys_bias_scale` 是退化外推偏差的粗代理 |r_post−r_pre|/(r_pre+r_post)(评审④,
/// 下界启发式);cross-ref 段为 0。`maskable`=(None,None) 地板填充,不可信、可屏蔽(评审③)。
#[derive(Debug, Clone, Default)]
pub struct GapCovBlock {
    /// cross-ref:参与形状重建的各参考盒标定(空 = degenerate)
    pub refs: Vec<GapRefCalib>,
    /// degenerate:gap 前端点率(None = 该端无有效率)
    pub r_pre: Option<f64>,
    /// degenerate:gap 后端点率
    pub r_post: Option<f64>,
    /// degenerate:前端点 packet 实际事例数(自由度 = n_pre−1)
    pub n_pre: Option<f64>,
    /// degenerate:后端点 packet 实际事例数
    pub n_post: Option<f64>,
    /// (None,None) 地板填充:率无从估、不入协方差、下游可屏蔽(评审③)
    pub maskable: bool,
    /// 系统偏差量级粗代理(评审④);cross-ref 段为 0
    pub sys_bias_scale: f64,
    /// cross-ref 重整因子 ρ=N_lost/Σshape(spec §7,S 系数 ρk_b/n_m);退化段为 0
    pub rho: f64,
}

/// gap 内单个 1ms 格的 kind(spec ③ gapbins)。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapBinKind {
    /// 有参考盒覆盖:S 用恒等,filler↔参考系数 ρk_b/n_m
    Measured,
    /// 无参考盒,靠端点插值:S 借 left/right 端点格的参考事件
    Empty,
}

/// gap 内单个 1ms 格的结构描述(spec ③ gapbins)。measured 格记有效参考盒数 n_m;
/// empty 格记插值端点 (left_bin,right_bin) 与权重 τ,下游据此在 1ms 网格精确拼 S,
/// 不必重猜插值算子(避免与 Rust 漂移)。
#[derive(Debug, Clone)]
pub struct GapBinInfo {
    /// gap 内 1ms 格序号(0..n_sbins)
    pub bin_index: usize,
    /// 该格下沿的 MET 时间
    pub t_lo: f64,
    /// measured:有效参考盒数;empty 不适用(下游用端点格的 n_m),置 0
    pub n_m: usize,
    /// measured / empty
    pub kind: GapBinKind,
    /// empty:插值左端点 bin_index(measured 为 None)
    pub left_bin: Option<usize>,
    /// empty:插值右端点 bin_index(measured 为 None)
    pub right_bin: Option<usize>,
    /// empty:插值权重 τ=(i−l)/(r−l)(measured 为 None)
    pub tau: Option<f64>,
}

/// 重建后的补全事件
#[derive(Debug, Clone)]
pub struct ReconstructedGap {
    /// 对应的 gap 索引
    pub gap_idx: usize,
    /// 补全的事件 MET 时间
    pub filled_events: Vec<f64>,
    /// N_lost
    pub n_lost: usize,
    /// 是否使用了交叉参考
    pub has_cross_ref: bool,
    /// 每个 filler 的“仅方差”权重 v。cross-ref gap = 0（方差由参考事件承载）；
    /// 退化 gap = √(σ²_gap / N_lost)，σ²_gap 由 pre/post 包率涨落传播，使方差
    /// 局域在 gap 内（挂 filler）而非畸变相邻真实事件。
    pub filler_weight: f64,
    /// spec §13 协方差块描述子(下游装配协方差用)
    pub cov: GapCovBlock,
    /// spec ③:gap 内每 1ms 格结构(measured/empty、n_m、插值端点/τ)。cross-ref
    /// gap 才有;退化 gap 为空(方差走 r項/filler_weight,不经 S filler↔参考)。
    pub bins: Vec<GapBinInfo>,
}

const SHAPE_BIN_WIDTH: f64 = 0.001; // 1ms

/// 对单个 box 的 FIFO reset gap 进行光变曲线重建。
///
/// 算法：
/// 1. 用参考 box 的事件分布构建校准后的形状函数（1ms bin）
/// 2. N_lost = shape 总和（校准后的参考计数直接给出丢失事件数）
/// 3. 按形状分配事件到各 bin
///
/// 当参考 box 不可用（所有 box 同时饱和）时，退化为 post-reset 包率估算 + 均匀分配。
pub fn reconstruct_gaps(
    target: &BoxReconstructionData,
    references: &[&BoxReconstructionData],
) -> (Vec<ReconstructedGap>, Vec<Vec<f64>>) {
    let mut results = Vec::new();
    // per-reference-box, per-event particle-weight contribution from THIS
    // target's gaps (aligned with `references`); accumulated below.
    let mut ref_weights: Vec<Vec<f64>> = references
        .iter()
        .map(|r| vec![0.0f64; r.events.len()])
        .collect();

    for (gap_idx, gap) in target.gaps.iter().enumerate() {
        let gap_start = gap.start_met;
        let gap_stop = gap.stop_met;
        let gap_dur = gap_stop - gap_start;
        if gap_dur <= 0.0 {
            continue;
        }

        // 步骤一：构建形状 bin
        let n_sbins = ((gap_dur / SHAPE_BIN_WIDTH).ceil() as usize).max(1);
        let actual_sbin = gap_dur / n_sbins as f64;
        let mut shape = vec![0.0f64; n_sbins];
        // 每个 bin 记录贡献它的 (ref_idx, 事件起, 事件止, k)，供逐事件权重反算
        let mut bin_refs: Vec<Vec<(usize, usize, usize, f64)>> =
            vec![Vec::new(); n_sbins];

        // 汇总所有参考 box 的事件，统一构建形状函数
        let mut has_ref = false;
        for si in 0..n_sbins {
            let bin_lo = gap_start + si as f64 * actual_sbin;
            let bin_hi = bin_lo + actual_sbin;
            let bin_mid = (bin_lo + bin_hi) / 2.0;

            let mut total_ref_count = 0.0;

            for (ref_idx, ref_data) in references.iter().enumerate() {
                if is_in_unreliable(bin_mid, &ref_data.unreliable) {
                    continue;
                }

                let lo_idx = ref_data.events.partition_point(|&t| t < bin_lo);
                let hi_idx = ref_data.events.partition_point(|&t| t < bin_hi);
                let count = (hi_idx - lo_idx) as f64;

                if count > 0.0 {
                    let k = calibrate_ratio_sorted(
                        &target.events, &ref_data.events,
                        &target.unreliable, &ref_data.unreliable,
                        gap_start, gap_stop, 0.5,
                    );
                    total_ref_count += count * k;
                    bin_refs[si].push((ref_idx, lo_idx, hi_idx, k));
                }
            }

            let n_valid_refs = bin_refs[si].len();
            if n_valid_refs > 0 {
                shape[si] = total_ref_count / n_valid_refs as f64;
                has_ref = true;
            }
        }

        // ── diagnostic: calibration ratio breakdown ──
        for (ri, ref_data) in references.iter().enumerate() {
            let dw = [(gap_start - 0.5, gap_start), (gap_stop, gap_stop + 0.5)];
            let (mut tc, mut rc) = (0usize, 0usize);
            let (mut te, mut re) = (0.0f64, 0.0f64);
            for &(wl, wh) in &dw {
                let t = effective_duration(wl, wh, &target.unreliable);
                if t > 1e-6 {
                    let a = target.events.partition_point(|&x| x < wl);
                    let b = target.events.partition_point(|&x| x < wh);
                    tc += target.events[a..b].iter()
                        .filter(|&&x| !is_in_unreliable(x, &target.unreliable)).count();
                    te += t;
                }
                let r = effective_duration(wl, wh, &ref_data.unreliable);
                if r > 1e-6 {
                    let a = ref_data.events.partition_point(|&x| x < wl);
                    let b = ref_data.events.partition_point(|&x| x < wh);
                    rc += ref_data.events[a..b].iter()
                        .filter(|&&x| !is_in_unreliable(x, &ref_data.unreliable)).count();
                    re += r;
                }
            }
            if re > 1e-6 && rc > 10 && te > 1e-6 {
                let tr = tc as f64 / te;
                let rr = rc as f64 / re;
                eprintln!("  gap[{gap_idx}] ref[{ri}]: k={:.4}  tgt={tc}/{te:.4}s={tr:.0}/s  ref={rc}/{re:.4}s={rr:.0}/s", tr/rr);
            } else {
                eprintln!("  gap[{gap_idx}] ref[{ri}]: k=1.0(default)  tgt={tc}/{te:.4}s  ref={rc}/{re:.4}s");
            }
        }

        // 步骤二：确定 N_lost、退化方差、是否真用了 cross-ref
        let n_lost;
        // None = cross-ref(filler 权重 0，方差在参考事件)；Some(σ²) = 退化
        // (Some(NaN) = (None,None) 地板填充，按 ~100% 不确定处理)
        let filler_sigma2: Option<f64>;
        let used_cross_ref: bool;
        // spec ③/§7:cross-ref 分支填这两项(退化分支保持默认:ρ=0、空格结构)。
        let mut cross_rho = 0.0f64;
        let mut gap_bins: Vec<GapBinInfo> = Vec::new();
        if has_ref {
            let n_filled = shape.iter().filter(|&&v| v > 0.0).count();
            if n_filled * 100 / n_sbins >= 30 {
                // 参考覆盖充分：用 shape 总和作为 N_lost
                let (lambda, interp) = interpolate_empty_bins(&mut shape);
                let total: f64 = shape.iter().sum();
                n_lost = total.round() as usize;
                // 逐参考事件累加权重贡献 w_contrib = ρ · k/n_m · (1+Λ)。
                // 每个 filler 都完整归回它的源参考事件（守恒：Σ contrib = n_lost）。
                let rho = if total > 0.0 { n_lost as f64 / total } else { 0.0 };
                for (si2, refs_in_bin) in bin_refs.iter().enumerate() {
                    let n_valid = refs_in_bin.len();
                    if n_valid == 0 {
                        continue;
                    }
                    let amp = rho * (1.0 + lambda[si2]) / n_valid as f64;
                    for &(ref_idx, lo, hi, k) in refs_in_bin {
                        let contrib = amp * k;
                        for ev in lo..hi {
                            ref_weights[ref_idx][ev] += contrib;
                        }
                    }
                }
                // spec ③:把逐格结构序列化出来(内部已算好,别让下游重推)。interp[si2]
                // None ⟺ 该格 measured(有参考);Some((l,r,τ)) ⟺ 空格,借端点插值。
                cross_rho = rho;
                gap_bins = (0..n_sbins)
                    .map(|si2| {
                        let t_lo = gap_start + si2 as f64 * actual_sbin;
                        match interp[si2] {
                            None => GapBinInfo {
                                bin_index: si2,
                                t_lo,
                                n_m: bin_refs[si2].len(),
                                kind: GapBinKind::Measured,
                                left_bin: None,
                                right_bin: None,
                                tau: None,
                            },
                            Some((l, r, t)) => GapBinInfo {
                                bin_index: si2,
                                t_lo,
                                n_m: 0,
                                kind: GapBinKind::Empty,
                                left_bin: Some(l),
                                right_bin: Some(r),
                                tau: Some(t),
                            },
                        }
                    })
                    .collect();
                filler_sigma2 = None; // 方差由参考事件承载
                used_cross_ref = true;
                eprintln!("gap[{gap_idx}]: {gap_dur:.4}s  n_lost={n_lost}  cross-ref  cov={n_filled}/{n_sbins}");
            } else {
                // 参考覆盖不足：退化为 pre/post 率线性插值（不是 cross-ref）
                fill_shape_fallback(&mut shape, gap, &target.packets);
                n_lost = (shape.iter().sum::<f64>() * actual_sbin).round() as usize;
                filler_sigma2 = Some(degenerate_gap_variance(&target.packets, gap, gap_dur));
                used_cross_ref = false;
                eprintln!("gap[{gap_idx}]: {gap_dur:.4}s  n_lost={n_lost}  FALLBACK  cov={n_filled}/{n_sbins}");
            }
        } else {
            // 无参考：pre/post 率线性插值
            fill_shape_fallback(&mut shape, gap, &target.packets);
            n_lost = (shape.iter().sum::<f64>() * actual_sbin).round() as usize;
            filler_sigma2 = Some(degenerate_gap_variance(&target.packets, gap, gap_dur));
            used_cross_ref = false;
            eprintln!("gap[{gap_idx}]: {gap_dur:.4}s  n_lost={n_lost}  NO-REF");
        }

        if n_lost == 0 {
            continue;
        }

        // 步骤三：按形状分配事件
        let total: f64 = shape.iter().sum();
        if total <= 0.0 {
            continue;
        }

        let mut filled_events = Vec::with_capacity(n_lost);
        for (si, &s) in shape.iter().enumerate() {
            let n_in_bin = (s / total * n_lost as f64).round() as usize;
            if n_in_bin > 0 {
                let bin_lo = gap_start + si as f64 * actual_sbin;
                let bin_hi = bin_lo + actual_sbin;
                let step = (bin_hi - bin_lo) / n_in_bin as f64;
                for j in 0..n_in_bin {
                    filled_events.push(bin_lo + (j as f64 + 0.5) * step);
                }
            }
        }

        // 退化 filler 的“仅方差”权重 v=√(σ²_gap/实际filler数),使 Σv²=σ²_gap。
        // 用 filled_events.len()(逐 bin round 后可 ≠ n_lost),否则下游按实际
        // filler 数求和时不守恒。cross-ref=0；(None,None) 地板填充按 σ²=filler数²
        // (~100% 不确定)→ v=√(filler数),绝不给 0。
        let n_filled_actual = filled_events.len();
        let filler_weight = match filler_sigma2 {
            None => 0.0,
            Some(_) if n_filled_actual == 0 => 0.0,
            Some(s2) if s2.is_nan() => (n_filled_actual as f64).sqrt(),
            Some(s2) => (s2 / n_filled_actual as f64).sqrt(),
        };

        // spec §13 协方差块描述子
        let cov = if used_cross_ref {
            // 收集参与形状重建的参考盒(bin_refs 里出现过),各算 k 与标定窗计数
            let mut seen = vec![false; references.len()];
            let mut refs_calib: Vec<GapRefCalib> = Vec::new();
            for refs_in_bin in &bin_refs {
                for &(ref_idx, _, _, _) in refs_in_bin {
                    if !seen[ref_idx] {
                        seen[ref_idx] = true;
                        let (k, c_a, c_ref) = calibrate_counts(
                            &target.events, &references[ref_idx].events,
                            &target.unreliable, &references[ref_idx].unreliable,
                            gap_start, gap_stop, 0.5,
                        );
                        refs_calib.push(GapRefCalib {
                            ref_idx, k, c_a_cal: c_a as f64, c_ref_cal: c_ref as f64,
                        });
                    }
                }
            }
            refs_calib.sort_by_key(|r| r.ref_idx);
            GapCovBlock { refs: refs_calib, rho: cross_rho, ..Default::default() }
        } else {
            // degenerate:端点率、实际事例数(自由度)、可屏蔽、系统偏代理
            let pre = packet_rate_and_n(&target.packets, gap.prev_pkt_idx);
            let post = packet_rate_and_n(&target.packets, gap.next_pkt_idx);
            let (r_pre, n_pre) = pre.map_or((None, None), |(r, n)| (Some(r), Some(n)));
            let (r_post, n_post) = post.map_or((None, None), |(r, n)| (Some(r), Some(n)));
            let maskable = pre.is_none() && post.is_none();
            // 评审④:退化外推偏差粗代理 = 端点率相对变化;地板段 ~100%;单侧变化未知→0
            let sys_bias_scale = match (r_pre, r_post) {
                (Some(a), Some(b)) if a + b > 0.0 => (b - a).abs() / (a + b),
                _ if maskable => 1.0,
                _ => 0.0,
            };
            GapCovBlock {
                r_pre, r_post, n_pre, n_post, maskable, sys_bias_scale,
                ..Default::default()
            }
        };

        results.push(ReconstructedGap {
            gap_idx,
            filled_events,
            n_lost,
            has_cross_ref: used_cross_ref,
            filler_weight,
            cov,
            bins: gap_bins,
        });
    }

    (results, ref_weights)
}

/// 退化(fallback)gap 的 gap 级方差 σ²_gap,由 pre/post 包率涨落传播:
/// Var(r)=r²/(n_events−1)。**评审②:率的分子与自由度都用各端点 packet 的实际
/// n_events**(reset 前后常是残包,≠名义 109),两端各用各的。两端有率 →
/// (T/2)²(r_pre²/(n_pre−1)+r_post²/(n_post−1)),单侧外推 → T²r²/(n−1)。方差最终挂在
/// gap 内的 filler 上(局域、不畸变相邻真实观测事件),每 filler v=√(σ²_gap/实际filler
/// 数),在分配后算(见调用处)。
/// (None,None):无相邻包率、只能用 MCU 地板率硬填,方差无从估 → 返回 NaN,
/// 由调用方按 ~100% 不确定(σ²=filler数²)处理,绝不给 0。
fn degenerate_gap_variance(
    packets: &[PacketInfo],
    gap: &SaturationInterval,
    gap_dur: f64,
) -> f64 {
    let dof = |n: f64| (n - 1.0).max(1.0);
    let r_pre = packet_rate_and_n(packets, gap.prev_pkt_idx);
    let r_post = packet_rate_and_n(packets, gap.next_pkt_idx);
    match (r_pre, r_post) {
        (Some((rp, np)), Some((rn, nn))) => {
            (gap_dur / 2.0).powi(2) * (rp * rp / dof(np) + rn * rn / dof(nn))
        }
        (Some((r, n)), None) | (None, Some((r, n))) => gap_dur.powi(2) * r * r / dof(n),
        (None, None) => f64::NAN,
    }
}

/// 从包的 span 估算事件率与实际事例数(自由度用)。评审②:分子用该包实际 n_events。
fn packet_rate_and_n(packets: &[PacketInfo], pkt_idx: usize) -> Option<(f64, f64)> {
    packets.iter().find(|p| p.pkt_idx == pkt_idx).and_then(|info| {
        let span = info.span();
        if span > 1e-9 {
            Some((info.n_events as f64 / span, info.n_events as f64))
        } else {
            None
        }
    })
}

/// 从包的 span 估算事件率。
fn packet_rate(packets: &[PacketInfo], pkt_idx: usize) -> Option<f64> {
    packet_rate_and_n(packets, pkt_idx).map(|(r, _)| r)
}

/// 无参考时的 fallback：用 pre/post-reset 包的率线性插值构建 shape。
fn fill_shape_fallback(
    shape: &mut [f64],
    gap: &SaturationInterval,
    packets: &[PacketInfo],
) {
    let r_pre = packet_rate(packets, gap.prev_pkt_idx);
    let r_post = packet_rate(packets, gap.next_pkt_idx);
    let n = shape.len();

    match (r_pre, r_post) {
        (Some(rp), Some(rn)) => {
            for (i, s) in shape.iter_mut().enumerate() {
                let t = (i as f64 + 0.5) / n as f64;
                *s = rp * (1.0 - t) + rn * t;
            }
        }
        (Some(r), None) | (None, Some(r)) => {
            shape.iter_mut().for_each(|v| *v = r);
        }
        (None, None) => {
            shape.iter_mut().for_each(|v| *v = MCU_READ_RATE_FLOOR);
        }
    }
}

/// 计算 target/ref 的事件率比值，用于将参考 box 的计数换算为 target box 的计数。
///
/// 在 gap 前后各 margin 秒的窗口内统计双方的事件率（events/有效秒），
/// 排除各自的 unreliable 区间。
/// 标定窗(gap 前后 ±margin)内 target/ref 的率比 k,连同两侧**事件计数**
/// (spec §5b 算 σ_k 用:(σ_k/k)²=1/C_a_cal+1/C_ref_cal)。返回 (k, C_a_cal, C_ref_cal);
/// 参考计数不足(≤10)或时长为 0 时 k=1.0(默认)。
fn calibrate_counts(
    target_events: &[f64],
    ref_events: &[f64],
    target_unreliable: &[UnreliableInterval],
    ref_unreliable: &[UnreliableInterval],
    gap_start: f64,
    gap_stop: f64,
    margin: f64,
) -> (f64, usize, usize) {
    let windows = [
        (gap_start - margin, gap_start),
        (gap_stop, gap_stop + margin),
    ];

    let mut target_count = 0usize;
    let mut ref_count = 0usize;
    let mut target_effective = 0.0f64;
    let mut ref_effective = 0.0f64;

    for &(win_lo, win_hi) in &windows {
        // target 侧：排除 target 的 unreliable 区间
        let t_eff = effective_duration(win_lo, win_hi, target_unreliable);
        if t_eff > 1e-6 {
            let a = target_events.partition_point(|&t| t < win_lo);
            let b = target_events.partition_point(|&t| t < win_hi);
            // 只计落在可信时段内的事件
            let cnt = target_events[a..b]
                .iter()
                .filter(|&&t| !is_in_unreliable(t, target_unreliable))
                .count();
            target_count += cnt;
            target_effective += t_eff;
        }

        // ref 侧：排除 ref 的 unreliable 区间
        let r_eff = effective_duration(win_lo, win_hi, ref_unreliable);
        if r_eff > 1e-6 {
            let a = ref_events.partition_point(|&t| t < win_lo);
            let b = ref_events.partition_point(|&t| t < win_hi);
            let cnt = ref_events[a..b]
                .iter()
                .filter(|&&t| !is_in_unreliable(t, ref_unreliable))
                .count();
            ref_count += cnt;
            ref_effective += r_eff;
        }
    }

    // 用事件率比值（而非事件数比值），补偿双方有效时长不同
    let k = if ref_effective > 1e-6 && ref_count > 10 && target_effective > 1e-6 {
        let target_rate = target_count as f64 / target_effective;
        let ref_rate = ref_count as f64 / ref_effective;
        target_rate / ref_rate
    } else {
        1.0
    };
    (k, target_count, ref_count)
}

/// 只要 k 的薄封装(bin 级 shape 构建用,不需要计数)。
fn calibrate_ratio_sorted(
    target_events: &[f64],
    ref_events: &[f64],
    target_unreliable: &[UnreliableInterval],
    ref_unreliable: &[UnreliableInterval],
    gap_start: f64,
    gap_stop: f64,
    margin: f64,
) -> f64 {
    calibrate_counts(
        target_events, ref_events, target_unreliable, ref_unreliable, gap_start,
        gap_stop, margin,
    )
    .0
}

/// 计算窗口 [lo, hi] 内排除 unreliable 区间后的有效时长。
fn effective_duration(lo: f64, hi: f64, unreliable: &[UnreliableInterval]) -> f64 {
    let mut excluded = 0.0;
    for iv in unreliable {
        let overlap_lo = iv.start.max(lo);
        let overlap_hi = iv.stop.min(hi);
        if overlap_hi > overlap_lo {
            excluded += overlap_hi - overlap_lo;
        }
    }
    (hi - lo - excluded).max(0.0)
}




/// 空 bin 插值：从最近的有值 bin 做线性插值，边缘用最近有值 bin 常数外推。
#[cfg(test)]
mod weight_tests {
    use super::*;

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

    fn make_box(events: Vec<f64>, gaps: Vec<SaturationInterval>) -> BoxReconstructionData {
        let n = events.len();
        let unreliable = gaps
            .iter()
            .map(|g| UnreliableInterval { start: g.start_met, stop: g.stop_met })
            .collect();
        BoxReconstructionData {
            events,
            channels: vec![100u16; n],
            pulse_widths: vec![60u8; n],
            gaps,
            packets: Vec::new(),
            packet_events: Vec::new(),
            unreliable,
        }
    }

    fn spread(lo: f64, hi: f64, n: usize) -> Vec<f64> {
        (0..n).map(|i| lo + (i as f64 + 0.5) * (hi - lo) / n as f64).collect()
    }

    /// 守恒:每个参考事件的权重贡献之和 == 总 filler 数。每个 filler 都完整
    /// 归回它的源参考事件,不漏不重,所以两者必须相等。
    #[test]
    fn ref_weight_contributions_sum_to_filler_count() {
        // 目标盒 A:仅在标定窗(gap 前后 0.5s)有事件,gap 内饱和、无事件
        let target = make_box(
            [spread(0.5, 1.0, 50), spread(1.1, 1.6, 50)].concat(),
            vec![si(1.0, 1.1)],
        );
        // 两个参考盒:标定窗 + gap 内密集(每盒 gap 内 200 事件 → 全覆盖)
        let ref_events = || {
            [spread(0.5, 1.0, 50), spread(1.0, 1.1, 200), spread(1.1, 1.6, 50)].concat()
        };
        let ref_b = make_box(ref_events(), vec![]);
        let ref_c = make_box(ref_events(), vec![]);
        let refs: Vec<&BoxReconstructionData> = vec![&ref_b, &ref_c];

        let (gaps, ref_weights) = reconstruct_gaps(&target, &refs);

        let n_lost_total: usize = gaps.iter().map(|g| g.n_lost).sum();
        let weight_sum: f64 = ref_weights.iter().flatten().sum();

        assert!(n_lost_total > 0, "should have filled events");
        assert!(
            (weight_sum - n_lost_total as f64).abs() < 1.0,
            "Σ ref weight contributions = {weight_sum}, expected filler count = {n_lost_total}"
        );
    }

    /// 守恒必须在**插值**(Λ>0)下也成立:参考只覆盖 gap 前半,后半靠外推。
    #[test]
    fn ref_weight_conserved_with_interpolation() {
        let target = make_box(
            [spread(0.5, 1.0, 50), spread(1.1, 1.6, 50)].concat(),
            vec![si(1.0, 1.1)],
        );
        // 参考只在 gap 前半 [1.0,1.05) 密集,后半 [1.05,1.1) 空 → 触发插值/外推
        let ref_events = || {
            [spread(0.5, 1.0, 50), spread(1.0, 1.05, 100), spread(1.1, 1.6, 50)].concat()
        };
        let ref_b = make_box(ref_events(), vec![]);
        let ref_c = make_box(ref_events(), vec![]);
        let refs: Vec<&BoxReconstructionData> = vec![&ref_b, &ref_c];

        let (gaps, ref_weights) = reconstruct_gaps(&target, &refs);
        let n_lost: usize = gaps.iter().map(|g| g.n_lost).sum();
        let wsum: f64 = ref_weights.iter().flatten().sum();

        assert!(n_lost > 0, "should have filled events");
        assert!(
            (wsum - n_lost as f64).abs() < 1.0,
            "插值下 Σcontrib = {wsum}, 应等于 n_lost = {n_lost}"
        );
    }

    /// 比守恒更强:验证权重落在**正确的事件**上、且**数值正确**。
    /// 单参考盒、k=1、满覆盖、无插值 → gap 内每个参考事件贡献恰为 ρ=1,
    /// gap 外(标定窗)事件贡献恰为 0。
    #[test]
    fn weight_lands_on_in_gap_refs_with_correct_value() {
        let target = make_box(
            [spread(0.5, 1.0, 50), spread(1.1, 1.6, 50)].concat(),
            vec![si(1.0, 1.1)],
        );
        // gap [1.0,1.1) 内密集 300 事件(每 bin 多个 → 无空 bin、无插值 Λ=0)
        // + 标定窗使 k=1。于是 gap 内每事件贡献 = ρ·k/n_m·(1+Λ) = 1·1/1·1 = 1
        let ref_ev =
            [spread(0.5, 1.0, 50), spread(1.0, 1.1, 300), spread(1.1, 1.6, 50)].concat();
        let ref_b = make_box(ref_ev.clone(), vec![]);
        let (_gaps, rw) = reconstruct_gaps(&target, &[&ref_b]);
        let w = &rw[0];

        for (i, &t) in ref_ev.iter().enumerate() {
            if (1.0..1.1).contains(&t) {
                // gap 内:贡献 = ρ·k/n_m·(1+Λ) = 1·1/1·1 = 1
                assert!(
                    (w[i] - 1.0).abs() < 1e-6,
                    "gap 内 ref[{i}] t={t}: 贡献={} 应=1.0",
                    w[i]
                );
            } else {
                // gap 外(标定窗)的参考事件不参与填充,贡献必须为 0
                assert!(
                    w[i].abs() < 1e-9,
                    "gap 外 ref[{i}] t={t}: 贡献={} 应=0",
                    w[i]
                );
            }
        }
    }

    /// 验证 1/n_m 分摊:两个满覆盖参考盒 → gap 内每个参考事件贡献恰为 0.5。
    #[test]
    fn two_references_halve_the_weight() {
        let target = make_box(
            [spread(0.5, 1.0, 50), spread(1.1, 1.6, 50)].concat(),
            vec![si(1.0, 1.1)],
        );
        let ref_ev =
            || [spread(0.5, 1.0, 50), spread(1.0, 1.1, 300), spread(1.1, 1.6, 50)].concat();
        let rb = make_box(ref_ev(), vec![]);
        let rc = make_box(ref_ev(), vec![]);
        let (_gaps, rw) = reconstruct_gaps(&target, &[&rb, &rc]);
        let ev = ref_ev();
        for (i, &t) in ev.iter().enumerate() {
            if (1.0..1.1).contains(&t) {
                assert!(
                    (rw[0][i] - 0.5).abs() < 1e-6 && (rw[1][i] - 0.5).abs() < 1e-6,
                    "n_m=2 时 gap 内贡献应=0.5, 得 {} / {}",
                    rw[0][i], rw[1][i]
                );
            }
        }
    }

    /// 退化(fallback):无参考盒 → 用 target 自己 pre/post 包率恢复。退化 gap 的
    /// 每个 filler 携带仅方差权重 v,Σv² == σ²_gap;方差挂 filler(gap 内)、不碰
    /// 相邻真实观测事件(ref_weights 全 0)。
    #[test]
    fn degenerate_gap_filler_weight_conserved() {
        // target 有 pre 包 [0.9,1.0)、post 包 [1.1,1.2),中间 gap [1.0,1.1)
        let events = [spread(0.9, 1.0, 109), spread(1.1, 1.2, 109)].concat();
        let mut target = make_box(events, vec![si(1.0, 1.1)]);
        target.packets = vec![
            PacketInfo { pkt_idx: 0, min_met: 0.9, max_met: 1.0, n_events: 109 },
            PacketInfo { pkt_idx: 1, min_met: 1.1, max_met: 1.2, n_events: 109 },
        ];
        target.gaps[0].prev_pkt_idx = 0;
        target.gaps[0].next_pkt_idx = 1;

        // 无参考盒 → NO-REF fallback
        let refs: Vec<&BoxReconstructionData> = vec![];
        let (gaps, _rw) = reconstruct_gaps(&target, &refs);
        let g = &gaps[0];
        assert!(g.n_lost > 0, "退化 gap 应产生 filler");
        assert!(g.filler_weight > 0.0, "退化 gap 的 filler 应带正方差权重");

        // σ²_gap = (T/2)²(r_pre²+r_post²)/(N-1)，r=109/0.1=1090，T=0.1，N-1=108
        let r = 1090.0_f64;
        let t = 0.1_f64;
        let sigma2 = (t / 2.0).powi(2) * (r * r + r * r) / 108.0;
        // Σv² 必须用**实际生成的 filler 数**(filled_events.len(),逐 bin round 后可
        // ≠ n_lost),否则下游按 filled_events 求和时不守恒。
        let sv2 = g.filled_events.len() as f64 * g.filler_weight * g.filler_weight;
        assert!(
            (sv2 - sigma2).abs() < sigma2 * 0.02,
            "Σv²(按实际 filler 数)={sv2} 应≈σ²_gap={sigma2}"
        );
    }

    /// 守护 k:标定系数必须真进入贡献。target 率=2×ref → k=2 → 每个 in-gap
    /// 参考事件 contrib = k/n_m = 2。删掉 `*k` 或写反 k 会让此测试变红。
    #[test]
    fn weight_uses_calibration_k() {
        let target = make_box(
            [spread(0.5, 1.0, 100), spread(1.1, 1.6, 100)].concat(),
            vec![si(1.0, 1.1)],
        );
        let ref_ev =
            [spread(0.5, 1.0, 50), spread(1.0, 1.1, 300), spread(1.1, 1.6, 50)].concat();
        let ref_b = make_box(ref_ev.clone(), vec![]);
        let (_g, rw) = reconstruct_gaps(&target, &[&ref_b]);
        for (i, &t) in ref_ev.iter().enumerate() {
            if (1.0..1.1).contains(&t) {
                assert!(
                    (rw[0][i] - 2.0).abs() < 1e-6,
                    "k=2 时 in-gap 贡献应=2,得 {}",
                    rw[0][i]
                );
            }
        }
    }

    /// 守护跨盒路由:只有盒 A(idx0) 有 gap,B/C 作参考。编排后 A 自己不该收到
    /// 任何贡献(weights[0] 全 1)——若 gj 映射错位、把 B/C 的贡献落到 A 会红。
    #[test]
    fn contributions_route_to_correct_box() {
        let a = make_box(
            [spread(0.5, 1.0, 50), spread(1.1, 1.6, 50)].concat(),
            vec![si(1.0, 1.1)],
        );
        let ref_ev =
            || [spread(0.5, 1.0, 50), spread(1.0, 1.1, 300), spread(1.1, 1.6, 50)].concat();
        let boxes = [a, make_box(ref_ev(), vec![]), make_box(ref_ev(), vec![])];

        let mut weights: Vec<Vec<f64>> =
            boxes.iter().map(|bx| vec![1.0f64; bx.events.len()]).collect();
        for i in 0..boxes.len() {
            let refs_with_idx: Vec<(usize, &BoxReconstructionData)> = boxes
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(j, d)| (j, d))
                .collect();
            let refs: Vec<&BoxReconstructionData> =
                refs_with_idx.iter().map(|(_, d)| *d).collect();
            let (_g, rw) = reconstruct_gaps(&boxes[i], &refs);
            for (r, (gj, _)) in refs_with_idx.iter().enumerate() {
                for (ev, &cc) in rw[r].iter().enumerate() {
                    weights[*gj][ev] += cc;
                }
            }
        }
        assert!(
            weights[0].iter().all(|&w| (w - 1.0).abs() < 1e-9),
            "盒 A 不应收到任何贡献,却有 w≠1"
        );
        assert!(weights[1].iter().any(|&w| w > 1.0001), "B 应有被参考的事件");
        assert!(weights[2].iter().any(|&w| w > 1.0001), "C 应有被参考的事件");
    }

    /// 守护防双计:cross-ref gap 的 filler_weight 必须为 0(方差在参考事件,
    /// 不能同时挂到 filler)。
    #[test]
    fn cross_ref_gap_has_zero_filler_weight() {
        let target = make_box(
            [spread(0.5, 1.0, 50), spread(1.1, 1.6, 50)].concat(),
            vec![si(1.0, 1.1)],
        );
        let ref_ev =
            [spread(0.5, 1.0, 50), spread(1.0, 1.1, 300), spread(1.1, 1.6, 50)].concat();
        let (gaps, _rw) = reconstruct_gaps(&target, &[&make_box(ref_ev, vec![])]);
        assert!(gaps[0].has_cross_ref, "该 gap 应判为 cross-ref");
        assert!(gaps[0].filler_weight.abs() < 1e-12, "cross-ref filler_weight 必须为 0");
    }

    /// Bug:有参考但覆盖<30% 的 gap 走退化路径,不该标 has_cross_ref=true。
    #[test]
    fn low_coverage_gap_not_marked_cross_ref() {
        let mut target = make_box(
            [spread(0.5, 1.0, 50), spread(1.1, 1.6, 50)].concat(),
            vec![si(1.0, 1.1)],
        );
        target.packets = vec![
            PacketInfo { pkt_idx: 0, min_met: 0.9, max_met: 1.0, n_events: 109 },
            PacketInfo { pkt_idx: 1, min_met: 1.1, max_met: 1.2, n_events: 109 },
        ];
        target.gaps[0].prev_pkt_idx = 0;
        target.gaps[0].next_pkt_idx = 1;
        // 参考只在 gap 前 ~20% [1.0,1.02) 有事件 → 覆盖<30% → 退化分支
        let ref_ev =
            [spread(0.5, 1.0, 50), spread(1.0, 1.02, 60), spread(1.1, 1.6, 50)].concat();
        let (gaps, _rw) = reconstruct_gaps(&target, &[&make_box(ref_ev, vec![])]);
        assert!(
            !gaps[0].has_cross_ref,
            "覆盖<30% 的退化 gap 不应标 has_cross_ref=true"
        );
    }

    /// Bug:两侧包都无有效率 (None,None) 时用 MCU 地板率造 filler,但方差不应为 0
    /// (最不确定的纯猜测填充反而零方差、被当精确测量)。
    #[test]
    fn no_packet_rate_degenerate_gap_gets_variance() {
        let mut target = make_box(
            [spread(0.5, 1.0, 50), spread(1.1, 1.6, 50)].concat(),
            vec![si(1.0, 1.1)],
        );
        // 包 span=0 → packet_rate 返回 None → (None,None)
        target.packets = vec![
            PacketInfo { pkt_idx: 0, min_met: 1.0, max_met: 1.0, n_events: 1 },
            PacketInfo { pkt_idx: 1, min_met: 1.0, max_met: 1.0, n_events: 1 },
        ];
        target.gaps[0].prev_pkt_idx = 0;
        target.gaps[0].next_pkt_idx = 1;
        let refs: Vec<&BoxReconstructionData> = vec![];
        let (gaps, _rw) = reconstruct_gaps(&target, &refs);
        assert!(gaps[0].n_lost > 0, "地板率仍会造 filler");
        assert!(
            gaps[0].filler_weight > 0.0,
            "(None,None) 地板填充的 filler 方差不应为 0"
        );
    }

    /// 评审②:率的分子必须是该 packet 的**实际事例数**,不是名义 109。
    /// 残包(reset 前后常见)只有 50 个事例、span=1s → 率=50/s 而非 109/s。
    #[test]
    fn packet_rate_uses_actual_n_events_not_109() {
        let packets =
            vec![PacketInfo { pkt_idx: 0, min_met: 0.0, max_met: 1.0, n_events: 50 }];
        let r = packet_rate(&packets, 0).expect("有 span 应给率");
        assert!((r - 50.0).abs() < 1e-9, "率应=n_events/span=50，得 {r}");
    }

    /// 评审②:退化方差的自由度必须用各端点 packet 的**实际 n_events−1**,不是
    /// 名义 108;且两端率各用各的 n_events(残包 pre=40、post=60)。
    #[test]
    fn degenerate_variance_uses_actual_dof_not_108() {
        let packets = vec![
            PacketInfo { pkt_idx: 0, min_met: 0.9, max_met: 1.0, n_events: 40 },
            PacketInfo { pkt_idx: 1, min_met: 1.1, max_met: 1.2, n_events: 60 },
        ];
        let mut gap = si(1.0, 1.1);
        gap.prev_pkt_idx = 0;
        gap.next_pkt_idx = 1;
        let gap_dur = 0.1_f64;
        let got = degenerate_gap_variance(&packets, &gap, gap_dur);
        // r_pre=40/0.1=400, r_post=60/0.1=600
        // σ²=(T/2)²(r_pre²/(n_pre−1)+r_post²/(n_post−1))=0.05²(400²/39+600²/59)
        let (rp, rn) = (400.0_f64, 600.0_f64);
        let expected = (gap_dur / 2.0).powi(2) * (rp * rp / 39.0 + rn * rn / 59.0);
        assert!(
            (got - expected).abs() < expected * 1e-9,
            "σ²_gap 应用实际自由度:得 {got}, 期望 {expected}"
        );
    }

    /// spec §13:退化 gap 的协方差块记录端点率、实际事例数、可屏蔽标志与系统偏
    /// 代理。r_pre=40/0.1=400、r_post=60/0.1=600、sys_bias_scale=|600−400|/1000=0.2。
    #[test]
    fn degenerate_gap_block_records_rates_and_sys_bias() {
        let events = [spread(0.9, 1.0, 40), spread(1.1, 1.2, 60)].concat();
        let mut target = make_box(events, vec![si(1.0, 1.1)]);
        target.packets = vec![
            PacketInfo { pkt_idx: 0, min_met: 0.9, max_met: 1.0, n_events: 40 },
            PacketInfo { pkt_idx: 1, min_met: 1.1, max_met: 1.2, n_events: 60 },
        ];
        target.gaps[0].prev_pkt_idx = 0;
        target.gaps[0].next_pkt_idx = 1;
        let refs: Vec<&BoxReconstructionData> = vec![];
        let (gaps, _rw) = reconstruct_gaps(&target, &refs);
        let cov = &gaps[0].cov;
        assert!((cov.r_pre.unwrap() - 400.0).abs() < 1e-6, "r_pre={:?}", cov.r_pre);
        assert!((cov.r_post.unwrap() - 600.0).abs() < 1e-6, "r_post={:?}", cov.r_post);
        assert!((cov.n_pre.unwrap() - 40.0).abs() < 1e-9);
        assert!((cov.n_post.unwrap() - 60.0).abs() < 1e-9);
        assert!(!cov.maskable, "有率的退化段不可屏蔽");
        assert!((cov.sys_bias_scale - 0.2).abs() < 1e-6, "sys_bias={}", cov.sys_bias_scale);
        assert!(cov.refs.is_empty(), "退化 gap 无 cross-ref 标定");
    }

    /// spec §13/§5b:cross-ref gap 的块记录每参考盒的 k 与标定窗计数(算 σ_k 用)。
    /// target 率=2×ref → k=2;标定窗 ±0.5s 内 target 200 事例、ref 100 事例。
    #[test]
    fn cross_ref_gap_block_records_k_and_calib_counts() {
        let target = make_box(
            [spread(0.5, 1.0, 100), spread(1.1, 1.6, 100)].concat(),
            vec![si(1.0, 1.1)],
        );
        let ref_ev =
            [spread(0.5, 1.0, 50), spread(1.0, 1.1, 300), spread(1.1, 1.6, 50)].concat();
        let ref_b = make_box(ref_ev, vec![]);
        let (gaps, _rw) = reconstruct_gaps(&target, &[&ref_b]);
        let cov = &gaps[0].cov;
        assert_eq!(cov.refs.len(), 1, "应有一个参考盒");
        let rc = &cov.refs[0];
        assert_eq!(rc.ref_idx, 0);
        assert!((rc.k - 2.0).abs() < 1e-6, "k 应=2,得 {}", rc.k);
        assert!((rc.c_a_cal - 200.0).abs() < 1e-9, "C_a_cal 应=200,得 {}", rc.c_a_cal);
        assert!((rc.c_ref_cal - 100.0).abs() < 1e-9, "C_ref_cal 应=100,得 {}", rc.c_ref_cal);
        assert!(cov.r_pre.is_none(), "cross-ref 段无退化端点率");
    }

    /// spec ③:interpolate_empty_bins 除 lambda 外还须逐格报告插值端点/τ。
    /// shape=[10,0,0,20] → bin1,bin2 为空,端点 l=0,r=3,τ=(i−l)/(r−l)。
    #[test]
    fn interp_reports_endpoints_and_tau() {
        let mut shape = vec![10.0, 0.0, 0.0, 20.0];
        let (_lambda, interp) = interpolate_empty_bins(&mut shape);
        assert!(interp[0].is_none(), "measured 格 interp 应为 None");
        assert!(interp[3].is_none(), "measured 格 interp 应为 None");
        let (l1, r1, t1) = interp[1].unwrap();
        assert_eq!((l1, r1), (0, 3), "端点应为 (0,3)");
        assert!((t1 - 1.0 / 3.0).abs() < 1e-9, "τ 应=(1−0)/(3−0)=1/3,得 {t1}");
        let (l2, r2, t2) = interp[2].unwrap();
        assert_eq!((l2, r2), (0, 3));
        assert!((t2 - 2.0 / 3.0).abs() < 1e-9, "τ 应=2/3,得 {t2}");
    }

    /// spec ③:单侧外推的空格也记端点(left=right=唯一有值格),τ 取 0/1 使
    /// (1−τ)、τ 两支全落到同一端点、权重和=1。
    #[test]
    fn interp_extrapolation_endpoints() {
        // bin0 左空(None,Some(1)) → 右外推;bin2 右空(Some(1),None) → 左外推
        let mut shape = vec![0.0, 5.0, 0.0];
        let (_lambda, interp) = interpolate_empty_bins(&mut shape);
        assert_eq!(interp[0].map(|(l, r, t)| (l, r, t)), Some((1, 1, 1.0)));
        assert_eq!(interp[2].map(|(l, r, t)| (l, r, t)), Some((1, 1, 0.0)));
    }

    /// spec ③:cross-ref gap 全覆盖 → 每 1ms 格 measured、n_m=参考盒数、端点留空;
    /// bins 长度=n_sbins;t_lo 从 gap_start 起。
    #[test]
    fn cross_ref_gap_populates_bin_structure() {
        let target = make_box(
            [spread(0.5, 1.0, 50), spread(1.1, 1.6, 50)].concat(),
            vec![si(1.0, 1.1)],
        );
        // 两参考盒,gap 内 400 事件/盒 → 每 1ms 格都有事件
        let ref_ev =
            || [spread(0.5, 1.0, 50), spread(1.0, 1.1, 400), spread(1.1, 1.6, 50)].concat();
        let rb = make_box(ref_ev(), vec![]);
        let rc = make_box(ref_ev(), vec![]);
        let (gaps, _rw) = reconstruct_gaps(&target, &[&rb, &rc]);
        let g = &gaps[0];
        // n_sbins = ceil(gap_dur/1ms).max(1);gap_dur 用 stop−start(浮点误差使
        // (1.1−1.0)/0.001≈100.0000009 → 101,勿硬编码 100)
        let gap_dur = 1.1_f64 - 1.0_f64;
        let n_sbins = ((gap_dur / SHAPE_BIN_WIDTH).ceil() as usize).max(1);
        assert_eq!(g.bins.len(), n_sbins, "格数应=n_sbins={n_sbins},得 {}", g.bins.len());
        assert!(
            g.bins.iter().all(|b| b.kind == GapBinKind::Measured),
            "全覆盖应全 measured"
        );
        assert!(g.bins.iter().all(|b| b.n_m == 2), "两参考盒 → n_m=2");
        assert!(
            g.bins.iter().all(|b| b.left_bin.is_none() && b.tau.is_none()),
            "measured 格端点/τ 应留空"
        );
        assert!((g.bins[0].t_lo - 1.0).abs() < 1e-9, "首格 t_lo=gap_start");
        assert_eq!(g.bins[0].bin_index, 0);
        assert_eq!(g.bins[n_sbins - 1].bin_index, n_sbins - 1);
        // bin_index 与 t_lo 单调递增
        assert!(g.bins.windows(2).all(|w| w[1].bin_index == w[0].bin_index + 1));
    }

    /// spec ③:cross-ref gap 后半无参考 → 那些格 kind=empty、n_m 不用、带插值端点/τ;
    /// measured 格反之。
    #[test]
    fn cross_ref_gap_marks_empty_bins() {
        let target = make_box(
            [spread(0.5, 1.0, 50), spread(1.1, 1.6, 50)].concat(),
            vec![si(1.0, 1.1)],
        );
        // 参考只在 gap 前半 [1.0,1.05) 有事件(覆盖~50%≥30% → cross-ref),后半空
        let ref_ev =
            [spread(0.5, 1.0, 50), spread(1.0, 1.05, 300), spread(1.1, 1.6, 50)].concat();
        let (gaps, _rw) = reconstruct_gaps(&target, &[&make_box(ref_ev, vec![])]);
        let g = &gaps[0];
        assert!(g.has_cross_ref, "覆盖≥30% 应判 cross-ref");
        let n_empty = g.bins.iter().filter(|b| b.kind == GapBinKind::Empty).count();
        assert!(n_empty > 0, "后半应出现空格");
        for b in g.bins.iter().filter(|b| b.kind == GapBinKind::Empty) {
            assert!(
                b.left_bin.is_some() && b.right_bin.is_some() && b.tau.is_some(),
                "空格必须带插值端点/τ"
            );
        }
        for b in g.bins.iter().filter(|b| b.kind == GapBinKind::Measured) {
            assert!(b.n_m >= 1, "measured 格 n_m≥1");
            assert!(b.left_bin.is_none() && b.tau.is_none(), "measured 格端点/τ 留空");
        }
    }

    /// spec ③/§7:cross-ref 块须记录重整因子 ρ=N_lost/Σshape。k=1 全覆盖 → ρ≈1。
    #[test]
    fn cross_ref_block_records_rho() {
        let target = make_box(
            [spread(0.5, 1.0, 50), spread(1.1, 1.6, 50)].concat(),
            vec![si(1.0, 1.1)],
        );
        let ref_ev =
            [spread(0.5, 1.0, 50), spread(1.0, 1.1, 400), spread(1.1, 1.6, 50)].concat();
        let (gaps, _rw) = reconstruct_gaps(&target, &[&make_box(ref_ev, vec![])]);
        let rho = gaps[0].cov.rho;
        assert!(rho > 0.0, "cross-ref 块应记录正的 ρ,得 {rho}");
        assert!((rho - 1.0).abs() < 0.1, "k=1 全覆盖 ρ 应≈1,得 {rho}");
    }

    /// degenerate gap 无 S filler↔参考结构 → bins 为空、ρ=0。
    #[test]
    fn degenerate_gap_has_empty_bins_and_zero_rho() {
        let events = [spread(0.9, 1.0, 109), spread(1.1, 1.2, 109)].concat();
        let mut target = make_box(events, vec![si(1.0, 1.1)]);
        target.packets = vec![
            PacketInfo { pkt_idx: 0, min_met: 0.9, max_met: 1.0, n_events: 109 },
            PacketInfo { pkt_idx: 1, min_met: 1.1, max_met: 1.2, n_events: 109 },
        ];
        target.gaps[0].prev_pkt_idx = 0;
        target.gaps[0].next_pkt_idx = 1;
        let refs: Vec<&BoxReconstructionData> = vec![];
        let (gaps, _rw) = reconstruct_gaps(&target, &refs);
        assert!(gaps[0].bins.is_empty(), "退化 gap 不产 S 格结构");
        assert!(gaps[0].cov.rho.abs() < 1e-12, "退化 gap ρ=0");
    }

    fn dense_except(lo: f64, hi: f64, per_sec: usize, gap: (f64, f64)) -> Vec<f64> {
        let n = ((hi - lo) * per_sec as f64) as usize;
        (0..n)
            .map(|i| lo + (i as f64 + 0.5) * (hi - lo) / n as f64)
            .filter(|&t| t < gap.0 || t >= gap.1)
            .collect()
    }

    /// 端到端守恒:三盒各有一个 gap、彼此互为参考。模拟编排层(每观测事件
    /// 初值 1、跨所有 target 累加贡献),验证 Σ(观测权重) == 总事件数(观测+填充)。
    /// 这多覆盖了"一个参考事件被多个盒的 gap 引用、贡献跨盒累加"这一层。
    #[test]
    fn total_evt_weight_equals_total_event_count() {
        let ga = (1.0, 1.1);
        let gb = (1.5, 1.6);
        let gc = (2.0, 2.1);
        let a = make_box(dense_except(0.5, 3.5, 2000, ga), vec![si(ga.0, ga.1)]);
        let b = make_box(dense_except(0.5, 3.5, 2000, gb), vec![si(gb.0, gb.1)]);
        let c = make_box(dense_except(0.5, 3.5, 2000, gc), vec![si(gc.0, gc.1)]);
        let boxes = [a, b, c];

        // 编排:每观测事件初值 1.0,累加各 target gap 传来的贡献
        let mut weights: Vec<Vec<f64>> =
            boxes.iter().map(|bx| vec![1.0f64; bx.events.len()]).collect();
        let mut total_fill = 0usize;
        for i in 0..boxes.len() {
            let refs_with_idx: Vec<(usize, &BoxReconstructionData)> = boxes
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(j, d)| (j, d))
                .collect();
            let refs: Vec<&BoxReconstructionData> =
                refs_with_idx.iter().map(|(_, d)| *d).collect();
            let (gaps, ref_w) = reconstruct_gaps(&boxes[i], &refs);
            total_fill += gaps.iter().map(|g| g.n_lost).sum::<usize>();
            for (r, (gj, _)) in refs_with_idx.iter().enumerate() {
                for (ev, &c) in ref_w[r].iter().enumerate() {
                    weights[*gj][ev] += c;
                }
            }
        }

        let total_obs: usize = boxes.iter().map(|bx| bx.events.len()).sum();
        let w_sum: f64 = weights.iter().flatten().sum();
        let expected = (total_obs + total_fill) as f64;
        assert!(total_fill > 0, "should have filled events");
        assert!(
            (w_sum - expected).abs() < 1.0,
            "Σ EVT weights = {w_sum}, expected total events (obs {total_obs} + fill {total_fill}) = {expected}"
        );
    }
}

/// 线性插值填补空 bin，并返回每个（原本有值的）端点 bin 被空 bin 借用的
/// 系数和 Λ_l：空 bin i 以 l 为左端点则 l 得 (1-t)，为右端点则得 t，外推
/// 则整份 1。用于把插值 bin 的 filler 权重反算回端点的参考事件。
fn interpolate_empty_bins(
    shape: &mut [f64],
) -> (Vec<f64>, Vec<Option<(usize, usize, f64)>>) {
    let n = shape.len();
    let mut lambda = vec![0.0f64; n];
    // spec ③:逐格插值端点/τ。None = 该格 measured(有值);Some((l,r,τ)) = 该格空,
    // 借端点 l/r 插值。单侧外推记 (e,e,0/1) 使权重全落到唯一端点。
    let mut interp: Vec<Option<(usize, usize, f64)>> = vec![None; n];
    if n == 0 {
        return (lambda, interp);
    }

    // 预计算每个位置左边和右边最近的有值 bin 索引
    let mut left_filled: Vec<Option<usize>> = vec![None; n];
    let mut right_filled: Vec<Option<usize>> = vec![None; n];

    let mut last = None;
    for i in 0..n {
        if shape[i] > 0.0 {
            last = Some(i);
        }
        left_filled[i] = last;
    }

    last = None;
    for i in (0..n).rev() {
        if shape[i] > 0.0 {
            last = Some(i);
        }
        right_filled[i] = last;
    }

    for i in 0..n {
        if shape[i] > 0.0 {
            continue;
        }
        match (left_filled[i], right_filled[i]) {
            (Some(l), Some(r)) => {
                // 两侧都有值：线性插值
                let t = (i - l) as f64 / (r - l) as f64;
                shape[i] = shape[l] * (1.0 - t) + shape[r] * t;
                lambda[l] += 1.0 - t;
                lambda[r] += t;
                interp[i] = Some((l, r, t));
            }
            (Some(l), None) => {   // 右侧无值：常数外推
                shape[i] = shape[l];
                lambda[l] += 1.0;
                interp[i] = Some((l, l, 0.0));
            }
            (None, Some(r)) => {   // 左侧无值：常数外推
                shape[i] = shape[r];
                lambda[r] += 1.0;
                interp[i] = Some((r, r, 1.0));
            }
            (None, None) => {}                          // 全空，不应到达此处
        }
    }
    (lambda, interp)
}
