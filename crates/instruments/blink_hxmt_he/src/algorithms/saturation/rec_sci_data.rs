use super::crc_check;
use super::detect::detect_fifo_reset_intervals;
use crate::io::level_1b::SciFile;
use crate::types::HxmtHe;
use blink_core::types::MissionElapsedTime;

// ─────────────────────────────────────────────────────────────────────────────
// 常量
// ─────────────────────────────────────────────────────────────────────────────

const PTIME_MOD: u64 = 1 << 19; // 524288

/// 1B→1K 时间校正 (秒)。在 2026-02-26T10 一小时 Box A 重叠窗口
/// （MET 446724004..446727603，3599 s, 9.52M 匹配事件）逐事例匹配：
/// 中位差 +0.0596μs（恰好 1 ulp of double MET at this magnitude），
/// σ = 0.0375μs，|max| = 0.238μs；1B-only = 0，1K-only = 42（均在
/// SEC 端点 1-5 ticks 内）。4.0 s 修正精确到优于 1μs，是 1K 管线
/// 固定延迟而非经验拟合。
const MET_CORRECTION: f64 = 4.0;

// ─────────────────────────────────────────────────────────────────────────────
// CCSDS 包解析
// ─────────────────────────────────────────────────────────────────────────────

enum Pack {
    Event {
        ptime: u64,
        channel: u8,
        raw_bytes: [u8; 8],
    },
    Second {
        stime: u64,
        ptime: u64,
    },
    Error,
}

/// 解析单个 CCSDS 包中所有事例
fn parse_events(ccsds: &[u8]) -> Vec<Pack> {
    let payload = &ccsds[6..878];
    let mut events = Vec::with_capacity(109);

    for chunk in payload.chunks_exact(8) {
        let mut row = [0u64; 8];
        let mut raw_bytes = [0u8; 8];
        for (i, byte) in chunk.iter().enumerate() {
            row[i] = *byte as u64;
            raw_bytes[i] = *byte;
        }

        let pack = if crc_check(&row) == row[7] & 0x0F {
            let ptime =
                ((row[4] & 1) << 18) + (row[5] << 10) + (row[6] << 2) + ((row[7] & 0xC0) >> 6);
            let channel = raw_bytes[0];
            match row[7] & 0x30 {
                0x00 | 0x20 => Pack::Event {
                    ptime,
                    channel,
                    raw_bytes,
                },
                0x10 => {
                    let stime = (row[0] << 24) + (row[1] << 16) + (row[2] << 8) + row[3];
                    Pack::Second { stime, ptime }
                }
                _ => Pack::Error,
            }
        } else {
            Pack::Error
        };
        events.push(pack);
    }
    events
}

fn get_utc_tail(ccsds: &[u8]) -> f64 {
    ccsds[878] as f64
        + (ccsds[879] as f64) * 256.0
        + (ccsds[880] as f64) * 65536.0
        + (ccsds[881] as f64) * 16777216.0
}

// ─────────────────────────────────────────────────────────────────────────────
// 时间重建（待重新设计，见 DESIGN.md）
// ─────────────────────────────────────────────────────────────────────────────

/// 时间重建主函数：输入 CCSDS 包序列，输出每包每事件的 MET。
/// NaN 表示无法确定的事件。
pub fn reconstruct_with_wrap_tracking(sci_data: &SciFile, offset: f64) -> Vec<Vec<f64>> {
    reconstruct_with_wrap_tracking_labeled(sci_data, offset, "")
}

pub fn reconstruct_with_wrap_tracking_labeled(
    sci_data: &SciFile,
    offset: f64,
    _label: &str,
) -> Vec<Vec<f64>> {
    let debug = std::env::var("DEBUG_WRAP").is_ok();
    let n_packets = sci_data.ccsds.len();

    if n_packets == 0 {
        return Vec::new();
    }

    // =====================================================================
    // Step 1: 解析所有包，CRC 过滤
    // =====================================================================
    // 每个 CCSDS 包 109 个 slot，通过 CRC 的为 EVT 或 SEC，不通过的为 Error。
    // 只保留通过 CRC 的事件，记录其 (pkt_idx, evt_idx, ptime, 类型)。

    let mut n_evt_total = 0u64;
    let mut n_sec_total = 0u64;
    let mut n_err_total = 0u64;

    // parsed[pkt_idx] = Vec of parsed events (只包含通过 CRC 的事件)
    struct ParsedEvent {
        ptime: u64,
        stime: Option<u64>,  // SEC 才有
    }

    let mut parsed: Vec<Vec<ParsedEvent>> = Vec::with_capacity(n_packets);

    for (_pkt_idx, ccsds) in sci_data.ccsds.iter().enumerate() {
        let events = parse_events(ccsds);
        let mut pkt_events: Vec<ParsedEvent> = Vec::new();

        for event in &events {
            match event {
                Pack::Event { ptime, .. } => {
                    pkt_events.push(ParsedEvent {
                        ptime: *ptime,
                        stime: None,
                    });
                    n_evt_total += 1;
                }
                Pack::Second { stime, ptime } => {
                    pkt_events.push(ParsedEvent {
                        ptime: *ptime,
                        stime: Some(*stime),
                    });
                    n_sec_total += 1;
                }
                Pack::Error => {
                    n_err_total += 1;
                }
            }
        }

        parsed.push(pkt_events);
    }

    if debug {
        eprintln!(
            "STEP1 CRC filter: {} EVT + {} SEC = {} pass, {} error ({:.1}%)",
            n_evt_total, n_sec_total, n_evt_total + n_sec_total, n_err_total,
            100.0 * n_err_total as f64 / (n_evt_total + n_sec_total + n_err_total) as f64
        );
    }

    // =====================================================================
    // Step 2: 找出所有有效 SEC，过滤幽灵 SEC
    // =====================================================================
    // SEC 提供绝对时间锚点 (stime, ptime)，但 CRC 碰撞会产生幽灵 SEC。
    //
    // 有效 SEC 满足：(ptime - stime × 500000) mod 524288 ≈ 常数（硬件相位）
    // 幽灵 SEC 的相位随机分布在 0~524287。
    //
    // 算法：
    //   Phase 1: 排序+滑动窗口找最大相位簇
    //   Phase 2: 逐对验证 stime-ptime 一致性，踢掉混入的 ghost

    const TICKS_PER_SEC: i64 = 500000; // 1s / 2μs
    const PHASE_TOLERANCE: i64 = 200;  // ±200 ticks (±0.4ms)，覆盖硬件抖动

    struct SecEvent {
        pkt_idx: usize,
        evt_idx: usize,  // 在 parsed[pkt_idx] 中的下标
        stime: u64,
        ptime: u64,
        met: f64,
    }

    // 收集所有 SEC 候选
    let mut all_secs: Vec<SecEvent> = Vec::new();
    for (pkt_idx, pkt) in parsed.iter().enumerate() {
        for (local_idx, evt) in pkt.iter().enumerate() {
            if let Some(stime) = evt.stime {
                all_secs.push(SecEvent {
                    pkt_idx,
                    evt_idx: local_idx,
                    stime,
                    ptime: evt.ptime,
                    met: stime as f64 + offset,
                });
            }
        }
    }

    // Phase 1: 排序+滑动窗口找最大相位簇
    // 计算每个 SEC 的相位
    let phases: Vec<i64> = all_secs.iter()
        .map(|s| ((s.ptime as i64 - s.stime as i64 * TICKS_PER_SEC) % PTIME_MOD as i64 + PTIME_MOD as i64) % PTIME_MOD as i64)
        .collect();

    // 按 phase 排序的下标
    let mut sorted_idx: Vec<usize> = (0..phases.len()).collect();
    sorted_idx.sort_by_key(|&i| phases[i]);

    // 滑动窗口：找包含最多点的窗口（宽度 = 2 × PHASE_TOLERANCE）
    let window_width = 2 * PHASE_TOLERANCE;
    let mut best_start = 0usize;
    let mut best_count = 0usize;
    let mut left = 0usize;

    for right in 0..sorted_idx.len() {
        // 收缩左边界，保持窗口宽度
        while phases[sorted_idx[right]] - phases[sorted_idx[left]] > window_width {
            left += 1;
        }
        let count = right - left + 1;
        if count > best_count {
            best_count = count;
            best_start = left;
        }
    }

    // 标记簇内的 SEC
    let mut in_cluster = vec![false; all_secs.len()];
    for i in best_start..(best_start + best_count) {
        in_cluster[sorted_idx[i]] = true;
    }

    if debug {
        let cluster_phases: Vec<i64> = (best_start..(best_start + best_count))
            .map(|i| phases[sorted_idx[i]])
            .collect();
        eprintln!(
            "STEP2 phase cluster: {} SECs, phase range {}~{} (span={})",
            best_count,
            cluster_phases.first().unwrap_or(&0),
            cluster_phases.last().unwrap_or(&0),
            cluster_phases.last().unwrap_or(&0) - cluster_phases.first().unwrap_or(&0),
        );
    }

    // Phase 2: 逐对验证 stime-ptime 一致性
    // 按打包顺序（pkt_idx, evt_idx）排序簇内的 SEC
    // 这是事件通过 FIFO 的物理顺序
    let mut cluster_indices: Vec<usize> = (0..all_secs.len())
        .filter(|&i| in_cluster[i])
        .collect();
    cluster_indices.sort_by_key(|&i| (all_secs[i].pkt_idx, all_secs[i].evt_idx));

    // Phase 2: stime 升序检查（LIS）
    // FIFO 保序 → 打包顺序中 stime 必须严格递增。
    // Ghost SEC 的 stime 是垃圾值，会打破升序。
    // 用 LIS 找主序列，不在 LIS 中的是 ghost。
    let mut is_valid = vec![false; all_secs.len()];

    if cluster_indices.len() > 1 {
        let vals: Vec<u64> = cluster_indices.iter()
            .map(|&i| all_secs[i].stime)
            .collect();

        let n = vals.len();
        let mut tails: Vec<u64> = Vec::new();
        let mut tail_pos: Vec<usize> = Vec::new();
        let mut parent: Vec<Option<usize>> = vec![None; n];

        for i in 0..n {
            let pos = tails.partition_point(|&t| t < vals[i]);
            if pos == tails.len() {
                tails.push(vals[i]);
                tail_pos.push(i);
            } else {
                tails[pos] = vals[i];
                tail_pos[pos] = i;
            }
            parent[i] = if pos > 0 { Some(tail_pos[pos - 1]) } else { None };
        }

        let mut in_lis = vec![false; n];
        let mut idx = *tail_pos.last().unwrap();
        loop {
            in_lis[idx] = true;
            match parent[idx] {
                Some(p) => idx = p,
                None => break,
            }
        }

        for (k, &ci) in cluster_indices.iter().enumerate() {
            if in_lis[k] {
                is_valid[ci] = true;
            }
        }
    } else if cluster_indices.len() == 1 {
        is_valid[cluster_indices[0]] = true;
    }

    // Phase 3: 绝对量级卫兵 —— stime 物理合理性
    //
    // 前两道检验都不看绝对量级:相位簇是配对自洽检验,随机 stime 以
    // ~400/524288 的概率漏网;stime-LIS 只查相对顺序,而流末尾/开头的
    // 荒谬值不违反单调性(实测:坏 SEC stime 跳 8-13 年、恰在小时末包,
    // LIS 免费收下 → 假 gap 数亿秒 → 每事件数亿 wrap 候选 → OOM)。
    // UTC 上界在 GPS 失效段又被 clamp 架空,三洞对齐即炸。
    //
    // 小时文件的真实 stime 跨度 ≤~4000 s(含前小时积压),偏离簇中位数
    // 超过一天的 stime 物理上不可能 —— 直接拒。对健康数据零影响(裕度
    // >20×)。星上秒计数器若真发生复位,复位后的 SEC 相位本就落在簇外,
    // 行为与既有逻辑一致,本卫兵不改变该情形。
    const STIME_SANITY_TOL: i64 = 86_400; // ±1 天
    let mut n_stime_outlier = 0u32;
    {
        let mut valid_stimes: Vec<u64> = all_secs.iter()
            .enumerate()
            .filter(|&(i, _)| is_valid[i])
            .map(|(_, s)| s.stime)
            .collect();
        if !valid_stimes.is_empty() {
            valid_stimes.sort_unstable();
            let median = valid_stimes[valid_stimes.len() / 2] as i64;
            for (i, sec) in all_secs.iter().enumerate() {
                if is_valid[i] && (sec.stime as i64 - median).abs() > STIME_SANITY_TOL {
                    is_valid[i] = false;
                    n_stime_outlier += 1;
                }
            }
        }
    }
    if debug && n_stime_outlier > 0 {
        eprintln!("STEP2.5 stime sanity guard: {} outlier(s) rejected", n_stime_outlier);
    }

    // 收集有效 SEC
    let valid_secs: Vec<&SecEvent> = all_secs.iter()
        .enumerate()
        .filter(|&(i, _)| is_valid[i])
        .map(|(_, s)| s)
        .collect();

    let n_valid_sec = valid_secs.len();
    let n_ghost_sec = all_secs.len() - n_valid_sec;

    if debug {
        eprintln!(
            "STEP2 SEC validated: {} valid, {} ghost (total {})",
            n_valid_sec, n_ghost_sec, all_secs.len()
        );

        if let (Some(first), Some(last)) = (valid_secs.first(), valid_secs.last()) {
            eprintln!(
                "  stime range: {} ~ {} ({} seconds)",
                first.stime, last.stime, last.stime - first.stime
            );
            eprintln!(
                "  pkt range: {} ~ {}",
                valid_secs.iter().map(|s| s.pkt_idx).min().unwrap(),
                valid_secs.iter().map(|s| s.pkt_idx).max().unwrap(),
            );

            // stime gap 分布
            let mut gap1 = 0u32;
            let mut gap2 = 0u32;
            let mut gap_other = Vec::new();
            for w in valid_secs.windows(2) {
                let gap = w[1].stime as i64 - w[0].stime as i64;
                match gap {
                    1 => gap1 += 1,
                    2 => gap2 += 1,
                    _ => gap_other.push((w[0].stime, w[1].stime, gap, w[0].pkt_idx, w[1].pkt_idx)),
                }
            }
            eprintln!("  stime gaps: {}×1s, {}×2s, {}×other", gap1, gap2, gap_other.len());
            for (s1, s2, gap, p1, p2) in &gap_other {
                let t_rel = *s1 as f64 + offset - 339945422.0;
                eprintln!("    stime {}→{} (gap={}s) pkt {}→{} T+{:.0}", s1, s2, gap, p1, p2, t_rel);
            }
        }
    }

    // =====================================================================
    // Step 3: 对所有相邻 SEC 对，解算中间事件的 MET
    // =====================================================================
    // Δstime=1: elapsed_fwd 唯一 (k=0)，直接对 elapsed_fwd 求 LIS
    // Δstime>1: 保持贪心+LIS（待改进）

    let mut result: Vec<Vec<f64>> = parsed.iter().map(|pkt| vec![f64::NAN; pkt.len()]).collect();
    let mut n_resolved = 0u64;
    let mut n_ghost_deadzone = 0u64;
    let mut n_ghost_order = 0u64;
    let mut n_sec_pairs = 0u64;

    // 有效 SEC 按打包顺序排列的索引
    let mut valid_indices: Vec<usize> = (0..all_secs.len())
        .filter(|&i| is_valid[i])
        .collect();
    valid_indices.sort_by_key(|&i| (all_secs[i].pkt_idx, all_secs[i].evt_idx));

    // 给 SEC 事件本身赋值 MET
    for &vi in &valid_indices {
        let sec = &all_secs[vi];
        result[sec.pkt_idx][sec.evt_idx] = sec.met + MET_CORRECTION;
    }

    // 对每对相邻有效 SEC
    for w in valid_indices.windows(2) {
        let sec1 = &all_secs[w[0]];
        let sec2 = &all_secs[w[1]];
        let ds = sec2.stime as i64 - sec1.stime as i64;

        if ds <= 0 {
            continue;
        }

        // 环境变量控制最大 gap：MAX_SEC_GAP=1 只处理 1s 对
        let max_gap: i64 = std::env::var("MAX_SEC_GAP")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(i64::MAX);
        if ds > max_gap {
            continue;
        }
        n_sec_pairs += 1;

        let met1 = sec1.met + MET_CORRECTION;
        let met2 = sec2.met + MET_CORRECTION;
        let pt1 = sec1.ptime as i64;

        // 两个 SEC 之间 ptime 的总前进量
        let total_ticks = ds * TICKS_PER_SEC;  // ds 秒 × 500000 ticks/s

        // 收集两个 SEC 之间的所有事件
        let pkt_a = sec1.pkt_idx;
        let evt_a = sec1.evt_idx;
        let pkt_b = sec2.pkt_idx;
        let evt_b = sec2.evt_idx;

        struct Candidate {
            pkt_idx: usize,
            local_idx: usize,
            elapsed_fwd: i64,  // mod PTIME_MOD
            utc_max_elapsed: i64,  // UTC tail 约束：elapsed 上界
        }

        let mut candidates: Vec<Candidate> = Vec::new();
        let mut last_pkt_idx: usize = usize::MAX;
        let mut cached_utc_max: i64 = total_ticks;
        for pkt_idx in pkt_a..=pkt_b {
            // 每个包计算一次 UTC tail 约束
            if pkt_idx != last_pkt_idx {
                let utc_tail = get_utc_tail(&sci_data.ccsds[pkt_idx]);
                // elapsed ≤ (utc_tail + 1 - sec1.met) / 2μs
                // +1 因为 UTC tail 是整秒截断
                let utc_ticks = ((utc_tail + 1.0 - sec1.met) / 2e-6) as i64;
                cached_utc_max = utc_ticks.clamp(0, total_ticks);
                last_pkt_idx = pkt_idx;
            }
            let start = if pkt_idx == pkt_a { evt_a + 1 } else { 0 };
            let end = if pkt_idx == pkt_b { evt_b } else { parsed[pkt_idx].len() };
            for local_idx in start..end {
                let pt = parsed[pkt_idx][local_idx].ptime as i64;
                let elapsed_fwd = (pt - pt1).rem_euclid(PTIME_MOD as i64);
                candidates.push(Candidate { pkt_idx, local_idx, elapsed_fwd, utc_max_elapsed: cached_utc_max });
            }
        }

        let pmod = PTIME_MOD as i64;
        let mut alive = vec![false; candidates.len()];
        let mut actual_elapsed = vec![0i64; candidates.len()];

        if ds == 1 {
            // ─── Δstime=1: 直接 LIS ───
            // elapsed_fwd 唯一（k=0），过滤 dead zone 后对 ef 求 LIS
            // LIS 自动排除幽灵事件（随机 ptime 打破升序），不会级联丢失

            // 收集有效候选（非 dead zone）的下标
            let mut valid_idx: Vec<usize> = Vec::new();
            for (i, c) in candidates.iter().enumerate() {
                if c.elapsed_fwd <= total_ticks {
                    actual_elapsed[i] = c.elapsed_fwd;
                    valid_idx.push(i);
                } else {
                    n_ghost_deadzone += 1;
                }
            }

            // 对有效候选的 elapsed_fwd 求 LIS
            if valid_idx.len() > 1 {
                let vals: Vec<i64> = valid_idx.iter()
                    .map(|&i| actual_elapsed[i])
                    .collect();

                let n = vals.len();
                let mut tails: Vec<i64> = Vec::new();
                let mut tail_pos: Vec<usize> = Vec::new();
                let mut parent: Vec<Option<usize>> = vec![None; n];

                for i in 0..n {
                    let pos = tails.partition_point(|&t| t < vals[i]);
                    if pos == tails.len() {
                        tails.push(vals[i]);
                        tail_pos.push(i);
                    } else {
                        tails[pos] = vals[i];
                        tail_pos[pos] = i;
                    }
                    parent[i] = if pos > 0 { Some(tail_pos[pos - 1]) } else { None };
                }

                // 回溯 LIS
                let mut in_lis = vec![false; n];
                let mut idx = *tail_pos.last().unwrap();
                loop {
                    in_lis[idx] = true;
                    match parent[idx] {
                        Some(p) => idx = p,
                        None => break,
                    }
                }

                for (k, &vi) in valid_idx.iter().enumerate() {
                    if in_lis[k] {
                        alive[vi] = true;
                    } else {
                        n_ghost_order += 1;
                    }
                }
            } else if valid_idx.len() == 1 {
                alive[valid_idx[0]] = true;
            }
        } else {
            // ─── Δstime>1 大 gap：括号快速路径（PPS 失锁段） ───
            //
            // 2020 年起平台 GPS 退化，PPS 进入守时模式后锁存相位出现几十 ms 的
            // 离散台阶（实测 ≤100ms，100ms 为守时量子）。失锁段 SEC 被全局相位
            // 簇（正确地）拒绝，有效锚点间出现数百至上千秒的空洞；穷举分组 LIS
            // 需物化 事件数×ds 个候选，单小时可达数百 GB。
            //
            // 快速路径原理：被拒 SEC 的 stime 是链验证过的精确整秒，FIFO 严格
            // 保序 ⇒ 它们把 gap 切成 ≤~1s 的小窗，窗内每个事件的 wrap 整数被
            // 钳死。事件时间仍完全由两端有效锚点 + ptime 计数器给出——被拒 SEC
            // 的锁存漂移只参与整数判定（容差半圈 524ms ≫ 实测漂移），从不进入
            // 时间值。逐链节残差断言 + 终端闭合校验：错一圈必不闭合 → 回退。
            // 回退 = 下方原穷举路径，行为与旧版完全一致。
            const BRACKET_MIN_DS: i64 = 60; // 只对大 gap 启用；小 gap 走原路径
            const RUN_LINK_TOL: i64 = 200; // 链内相位容差（±0.4ms/s，随 dst 放大）
            const RUN_MIN_LEN: usize = 3; // 链最短长度（幽灵 SEC 连环通过 ~1e-6）
            const BRACKET_GUARD: i64 = 200_000; // 跨链/端点残差上限（400ms < 半圈）
            const WINDOW_SLACK: i64 = 2_500; // 窗端漂移余量（5ms）

            let mut bracketed_done = false;
            if ds >= BRACKET_MIN_DS && std::env::var("BLINK_EXHAUSTIVE").is_err() {
                'bracket: {
                    // 1. 收集 gap 内全部 SEC 槽（流序），按逐秒相位一致性切 run
                    #[derive(Clone, Copy)]
                    struct BNode {
                        pkt_idx: usize,
                        local_idx: usize,
                        pt: i64,
                        stime: i64,
                        elapsed: i64,
                    }
                    let mut runs: Vec<Vec<BNode>> = Vec::new();
                    {
                        let mut cur: Vec<BNode> = Vec::new();
                        for pkt_idx in pkt_a..=pkt_b {
                            let start = if pkt_idx == pkt_a { evt_a + 1 } else { 0 };
                            let end = if pkt_idx == pkt_b { evt_b } else { parsed[pkt_idx].len() };
                            for local_idx in start..end {
                                let Some(st) = parsed[pkt_idx][local_idx].stime else {
                                    continue;
                                };
                                let node = BNode {
                                    pkt_idx,
                                    local_idx,
                                    pt: parsed[pkt_idx][local_idx].ptime as i64,
                                    stime: st as i64,
                                    elapsed: 0,
                                };
                                let linked = if let Some(last) = cur.last() {
                                    let dst = node.stime - last.stime;
                                    dst >= 1 && dst <= 30 && {
                                        let dfwd = (node.pt - last.pt).rem_euclid(pmod);
                                        let expect = (dst * TICKS_PER_SEC).rem_euclid(pmod);
                                        let d = (dfwd - expect).rem_euclid(pmod);
                                        d.min(pmod - d) <= RUN_LINK_TOL * dst
                                    }
                                } else {
                                    true
                                };
                                if !linked {
                                    runs.push(std::mem::take(&mut cur));
                                }
                                cur.push(node);
                            }
                        }
                        if !cur.is_empty() {
                            runs.push(cur);
                        }
                    }
                    // 2. 只保留长 run（防幽灵成链），且 stime 全局严格递增
                    let mut nodes: Vec<BNode> = Vec::new();
                    nodes.push(BNode {
                        pkt_idx: pkt_a,
                        local_idx: evt_a,
                        pt: pt1,
                        stime: sec1.stime as i64,
                        elapsed: 0,
                    });
                    for run in &runs {
                        if run.len() < RUN_MIN_LEN {
                            continue;
                        }
                        for n in run {
                            let last_st = nodes.last().unwrap().stime;
                            if n.stime > last_st && n.stime < sec2.stime as i64 {
                                nodes.push(*n);
                            }
                        }
                    }
                    if nodes.len() < 2 {
                        break 'bracket; // 无可用括号 → 穷举
                    }
                    // 3. 链节 elapsed 累计：wrap 整数由 stime 钳定，残差断言
                    for i in 1..nodes.len() {
                        let prev = nodes[i - 1];
                        let dst = nodes[i].stime - prev.stime;
                        let dfwd = (nodes[i].pt - prev.pt).rem_euclid(pmod);
                        let nominal = dst * TICKS_PER_SEC;
                        let w = ((nominal - dfwd) as f64 / pmod as f64).round() as i64;
                        if w < 0 {
                            break 'bracket;
                        }
                        let resid = nominal - (dfwd + w * pmod);
                        if resid.abs() > BRACKET_GUARD {
                            break 'bracket;
                        }
                        nodes[i].elapsed = prev.elapsed + dfwd + w * pmod;
                    }
                    // 4. 终端闭合：最后节点 → sec2，总 ticks 必须对上
                    {
                        let last = *nodes.last().unwrap();
                        let dst = sec2.stime as i64 - last.stime;
                        let dfwd = (sec2.ptime as i64 - last.pt).rem_euclid(pmod);
                        let nominal = dst * TICKS_PER_SEC;
                        let w = ((nominal - dfwd) as f64 / pmod as f64).round() as i64;
                        if w < 0 {
                            break 'bracket;
                        }
                        let resid = nominal - (dfwd + w * pmod);
                        if resid.abs() > BRACKET_GUARD {
                            break 'bracket;
                        }
                        let closure = last.elapsed + dfwd + w * pmod;
                        if (closure - total_ticks).abs() > 2 * BRACKET_GUARD {
                            break 'bracket;
                        }
                        nodes.push(BNode {
                            pkt_idx: pkt_b,
                            local_idx: evt_b,
                            pt: sec2.ptime as i64,
                            stime: sec2.stime as i64,
                            elapsed: closure,
                        });
                    }
                    // 5. 节点 SEC 本身赋 MET（两端锚点已在上游赋值）。
                    //    注意：时间来自 met1/met2 双向平均——节点自己的
                    //    stime+offset 从不作为时间值使用。
                    for n in &nodes[1..nodes.len() - 1] {
                        let met_fwd = met1 + n.elapsed as f64 * 2e-6;
                        let met_bwd = met2 - (total_ticks - n.elapsed) as f64 * 2e-6;
                        result[n.pkt_idx][n.local_idx] = (met_fwd + met_bwd) / 2.0;
                        n_resolved += 1;
                    }
                    // 6. 逐窗解算：窗内 wrap 候选 ≤2-3 个，套用与穷举同款的
                    //    降序分组 LIS；内存 O(窗内事件数)
                    let mut ci = 0usize; // candidates 与 nodes 同为流序，双指针
                    for wpair in nodes.windows(2) {
                        let (a, b) = (wpair[0], wpair[1]);
                        while ci < candidates.len() {
                            let c = &candidates[ci];
                            if (c.pkt_idx, c.local_idx) <= (a.pkt_idx, a.local_idx) {
                                ci += 1;
                            } else {
                                break;
                            }
                        }
                        let cstart = ci;
                        while ci < candidates.len() {
                            let c = &candidates[ci];
                            if (c.pkt_idx, c.local_idx) < (b.pkt_idx, b.local_idx) {
                                ci += 1;
                            } else {
                                break;
                            }
                        }
                        let cend = ci;
                        if cend == cstart {
                            continue;
                        }
                        let a_efwd = (a.pt - pt1).rem_euclid(pmod);
                        let limit = (b.elapsed - a.elapsed) + WINDOW_SLACK;
                        let mut entries: Vec<(usize, i64)> = Vec::new();
                        let mut tails: Vec<i64> = Vec::new();
                        let mut tail_entry: Vec<usize> = Vec::new();
                        let mut lis_parent: Vec<Option<usize>> = Vec::new();
                        for cand_idx in cstart..cend {
                            let c = &candidates[cand_idx];
                            let efwd_local = (c.elapsed_fwd - a_efwd).rem_euclid(pmod);
                            let mut w = (limit - efwd_local).div_euclid(pmod);
                            let push_entry = |e: i64,
                                                  entries: &mut Vec<(usize, i64)>,
                                                  tails: &mut Vec<i64>,
                                                  tail_entry: &mut Vec<usize>,
                                                  lis_parent: &mut Vec<Option<usize>>| {
                                let eidx = entries.len();
                                entries.push((cand_idx, e));
                                let pos = tails.partition_point(|&t| t < e);
                                lis_parent.push(if pos > 0 { Some(tail_entry[pos - 1]) } else { None });
                                if pos == tails.len() {
                                    tails.push(e);
                                    tail_entry.push(eidx);
                                } else {
                                    tails[pos] = e;
                                    tail_entry[pos] = eidx;
                                }
                            };
                            while w >= 0 {
                                let e = efwd_local + w * pmod; // 降序
                                push_entry(e, &mut entries, &mut tails, &mut tail_entry, &mut lis_parent);
                                w -= 1;
                            }
                            // FIFO 写入竞争：事件在 PPS 前数 μs 触发、却排在 SEC 槽
                            // 之后入流，efwd_local ≈ P−ε 会被 limit 误滤。解释为对
                            // 窗口起点的小负偏移（正向候选 P−ε 早已被滤，无歧义），
                            // met 仍由链 elapsed 给出。
                            if efwd_local > pmod - WINDOW_SLACK {
                                let e = efwd_local - pmod; // ∈ (−SLACK, 0)
                                push_entry(e, &mut entries, &mut tails, &mut tail_entry, &mut lis_parent);
                            }
                        }
                        if !tails.is_empty() {
                            let mut idx = *tail_entry.last().unwrap();
                            loop {
                                let (cand_idx, e_local) = entries[idx];
                                let c = &candidates[cand_idx];
                                if result[c.pkt_idx][c.local_idx].is_nan() {
                                    let eg = a.elapsed + e_local;
                                    let met_fwd = met1 + eg as f64 * 2e-6;
                                    let met_bwd = met2 - (total_ticks - eg) as f64 * 2e-6;
                                    result[c.pkt_idx][c.local_idx] = (met_fwd + met_bwd) / 2.0;
                                    n_resolved += 1;
                                }
                                match lis_parent[idx] {
                                    Some(p) => idx = p,
                                    None => break,
                                }
                            }
                        }
                        for cand_idx in cstart..cend {
                            let c = &candidates[cand_idx];
                            if result[c.pkt_idx][c.local_idx].is_nan() {
                                let efwd_local = (c.elapsed_fwd - a_efwd).rem_euclid(pmod);
                                if efwd_local <= limit || efwd_local > pmod - WINDOW_SLACK {
                                    n_ghost_order += 1;
                                } else {
                                    n_ghost_deadzone += 1;
                                }
                            }
                        }
                    }
                    if debug {
                        eprintln!(
                            "  BRACKET gap ds={}: engaged, {} nodes ({} runs)",
                            ds,
                            nodes.len(),
                            runs.len()
                        );
                    }
                    bracketed_done = true;
                }
                if debug && !bracketed_done {
                    eprintln!("  BRACKET gap ds={}: fallback to exhaustive", ds);
                }
            }

            if !bracketed_done {
            // ─── Δstime>1: 分组 LIS ───
            // 每个事件有 ds 个候选 elapsed = elapsed_fwd + w × PTIME_MOD, w ∈ [0, ds)
            // 全局求解：每个事件最多选一个候选，使选出的 elapsed 严格递增
            // 同一事件的候选按降序处理，避免同组候选互相"抬轿"

            let mut entries: Vec<(usize, i64)> = Vec::new(); // (event_idx, elapsed)
            let mut tails: Vec<i64> = Vec::new();
            let mut tail_entry: Vec<usize> = Vec::new();
            let mut lis_parent: Vec<Option<usize>> = Vec::new();

            for (event_idx, c) in candidates.iter().enumerate() {
                let mut cands: Vec<i64> = (0..ds)
                    .map(|w| c.elapsed_fwd + w * pmod)
                    .filter(|&e| e >= 0 && e <= total_ticks && e <= c.utc_max_elapsed)
                    .collect();
                cands.sort_unstable_by(|a, b| b.cmp(a)); // 降序

                for elapsed in cands {
                    let eidx = entries.len();
                    entries.push((event_idx, elapsed));

                    let pos = tails.partition_point(|&t| t < elapsed);
                    lis_parent.push(if pos > 0 { Some(tail_entry[pos - 1]) } else { None });

                    if pos == tails.len() {
                        tails.push(elapsed);
                        tail_entry.push(eidx);
                    } else {
                        tails[pos] = elapsed;
                        tail_entry[pos] = eidx;
                    }
                }
            }

            // 回溯标记 LIS 成员
            if !tails.is_empty() {
                let mut idx = *tail_entry.last().unwrap();
                loop {
                    let (ev, el) = entries[idx];
                    alive[ev] = true;
                    actual_elapsed[ev] = el;
                    match lis_parent[idx] {
                        Some(p) => idx = p,
                        None => break,
                    }
                }
            }

            // 统计未选中事件
            for (i, c) in candidates.iter().enumerate() {
                if !alive[i] {
                    let has_valid = (0..ds).any(|w| {
                        let e = c.elapsed_fwd + w * pmod;
                        e >= 0 && e <= total_ticks && e <= c.utc_max_elapsed
                    });
                    if has_valid {
                        n_ghost_order += 1;
                    } else {
                        n_ghost_deadzone += 1;
                    }
                }
            }
            } // if !bracketed_done（穷举回退路径结束）
        }

        // 赋值 MET
        for (i, c) in candidates.iter().enumerate() {
            if alive[i] {
                let met_fwd = met1 + actual_elapsed[i] as f64 * 2e-6;
                let remaining = total_ticks - actual_elapsed[i];
                let met_bwd = met2 - remaining as f64 * 2e-6;
                result[c.pkt_idx][c.local_idx] = (met_fwd + met_bwd) / 2.0;
                n_resolved += 1;
            }
        }
    }

    // =====================================================================
    // Step 4: 单边锚点外推（首/末 SEC 之外的事件）
    // =====================================================================
    // Step 3 的 windows(2) 只解算相邻 SEC 之间的事件。文件首个有效 SEC 之前、
    // 末个有效 SEC 之后的事件没有 pair 可以夹住，会被遗漏（典型每文件每盒
    // 约 3000 个）。这些事件相对最近的单边锚点距离若 < 1 ptime wrap 周期
    // （1.048576 s），(锚点ptime ± 事件ptime) mod PMOD 唯一确定 elapsed。
    //
    //   leading edge:  event_met = first.met - ((first.pt - evt.pt) mod PMOD) × 2μs
    //   trailing edge: event_met = last.met  + ((evt.pt - last.pt) mod PMOD)  × 2μs
    //
    // 但 mod 折叠本身分不清"1 圈内"和"跨了 N 圈"——(first_pt−pt) mod PMOD 的
    // 值域恒 < PMOD，拿它跟 PMOD−1 比永远通过（深饱和小时实测 190 万个跨圈
    // 前导事例被折进锚点前 1 s，形成假象超密段）。唯一的圈数判据是本包 UTC
    // 尾戳：外推 MET（换回 utc 标尺）偏离 utc_tail 超过容差即视为跨圈，留
    // NaN。容差沿用 SEC 有效性判据的 2 s（覆盖 FIFO 驻留 + 整秒截断）。
    let pmod = PTIME_MOD as i64;
    const EDGE_UTC_TOL: f64 = 2.0;
    let mut n_lead = 0u64;
    let mut n_trail = 0u64;
    let mut n_edge_rejected = 0u64;

    if let (Some(&first_vi), Some(&last_vi)) = (valid_indices.first(), valid_indices.last()) {
        // ----- Leading edge: events before the first valid SEC -----
        let first_sec = &all_secs[first_vi];
        let first_met = first_sec.met + MET_CORRECTION;
        let first_pt = first_sec.ptime as i64;
        let pkt_first = first_sec.pkt_idx;
        let evt_first = first_sec.evt_idx;

        for pkt_idx in 0..=pkt_first {
            let utc_tail = get_utc_tail(&sci_data.ccsds[pkt_idx]);
            let end = if pkt_idx == pkt_first { evt_first } else { parsed[pkt_idx].len() };
            for local_idx in 0..end {
                let evt = &parsed[pkt_idx][local_idx];
                if evt.stime.is_some() { continue; }                      // skip SECs (incl. ghosts)
                if !result[pkt_idx][local_idx].is_nan() { continue; }     // already resolved

                let pt = evt.ptime as i64;
                let elapsed_back = (first_pt - pt).rem_euclid(pmod);
                let met = first_met - elapsed_back as f64 * 2e-6;
                if (met - MET_CORRECTION - utc_tail).abs() > EDGE_UTC_TOL {
                    n_edge_rejected += 1;
                    continue;
                }

                result[pkt_idx][local_idx] = met;
                n_lead += 1;
            }
        }

        // ----- Trailing edge: events after the last valid SEC -----
        let last_sec = &all_secs[last_vi];
        let last_met = last_sec.met + MET_CORRECTION;
        let last_pt = last_sec.ptime as i64;
        let pkt_last = last_sec.pkt_idx;
        let evt_last = last_sec.evt_idx;

        for pkt_idx in pkt_last..parsed.len() {
            let utc_tail = get_utc_tail(&sci_data.ccsds[pkt_idx]);
            let start = if pkt_idx == pkt_last { evt_last + 1 } else { 0 };
            for local_idx in start..parsed[pkt_idx].len() {
                let evt = &parsed[pkt_idx][local_idx];
                if evt.stime.is_some() { continue; }
                if !result[pkt_idx][local_idx].is_nan() { continue; }

                let pt = evt.ptime as i64;
                let elapsed_fwd = (pt - last_pt).rem_euclid(pmod);
                let met = last_met + elapsed_fwd as f64 * 2e-6;
                if (met - MET_CORRECTION - utc_tail).abs() > EDGE_UTC_TOL {
                    n_edge_rejected += 1;
                    continue;
                }

                result[pkt_idx][local_idx] = met;
                n_trail += 1;
            }
        }
    }

    if debug {
        let n_total_events: usize = parsed.iter().map(|p| p.len()).sum();
        let n_nan = result.iter().flat_map(|p| p.iter()).filter(|m| m.is_nan()).count();
        eprintln!(
            "STEP3 1s-SEC pairs: {} pairs, {} events resolved",
            n_sec_pairs, n_resolved
        );
        eprintln!(
            "  ghosts: {} dead-zone + {} order-violation = {} total",
            n_ghost_deadzone, n_ghost_order, n_ghost_deadzone + n_ghost_order
        );
        eprintln!(
            "STEP4 single-anchor edges: {} leading + {} trailing = {} events resolved, {} rejected by UTC gate",
            n_lead, n_trail, n_lead + n_trail, n_edge_rejected
        );
        eprintln!(
            "  coverage: {}/{} events have MET ({:.1}%), {} NaN",
            n_total_events - n_nan, n_total_events,
            100.0 * (n_total_events - n_nan) as f64 / n_total_events as f64,
            n_nan
        );
    }

    result
}

// ─────────────────────────────────────────────────────────────────────────────
// 公共接口（供 CLI 和其他模块调用）
// ─────────────────────────────────────────────────────────────────────────────

/// 提取所有秒事例（Second event）的重建 MET 时间。
pub fn extract_second_event_times(sci_data: &SciFile, offset: f64) -> Vec<f64> {
    let mut second_times: Vec<f64> = Vec::new();

    for ccsds in sci_data.ccsds.iter() {
        let utc_tail = get_utc_tail(ccsds);
        let events = parse_events(ccsds);
        for event in &events {
            if let Pack::Second { stime, .. } = event {
                let met = *stime as f64 + offset;
                if (met - utc_tail).abs() < 2.0 {
                    second_times.push(met);
                }
            }
        }
    }

    second_times
}

/// 时间基准自检：SEC 槽的 (stime + offset) 与所在包 UTC 尾戳的偏差中位数
/// （秒，带符号）。任何时间线位移（如 eng offset 被平台 UTC 冻结污染）都会
/// 把这个量从 ~0 拉到位移量本身。偏差在 ±1 天外的 SEC 槽视为 CRC 碰撞幽灵
/// （stime 是随机垃圾），不参与统计；真 SEC 数量在整小时尺度上占多数，中位
/// 数对残余幽灵稳健。None = 没有可统计的 SEC 槽。
pub fn sec_utc_discrepancy(sci_data: &SciFile, offset: f64) -> Option<f64> {
    const GHOST_CUTOFF: f64 = 86_400.0;
    let mut diffs: Vec<f64> = Vec::new();
    for ccsds in sci_data.ccsds.iter() {
        let utc_tail = get_utc_tail(ccsds);
        for event in &parse_events(ccsds) {
            if let Pack::Second { stime, .. } = event {
                let diff = *stime as f64 + offset - utc_tail;
                if diff.abs() < GHOST_CUTOFF {
                    diffs.push(diff);
                }
            }
        }
    }
    if diffs.is_empty() {
        return None;
    }
    diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(diffs[diffs.len() / 2])
}

/// 重建所有事例的 MET 时间（扁平化）。
pub fn reconstruct_met_times(sci_data: &SciFile, offset: f64) -> Vec<f64> {
    reconstruct_with_wrap_tracking(sci_data, offset)
        .into_iter()
        .flatten()
        .collect()
}

/// 与 reconstruct_met_times 逐元素对齐的 wrapped channel 序列。
///
/// 按 reconstruct_with_wrap_tracking 的每包时间数截断地走 parse_events：
/// 每个 Event/Second 槽产出一个 channel（Error 跳过），与时间 flatten 的
/// 展平顺序、长度完全一致（含 NaN 时间槽）。SEC 槽记为 CHANNEL_SEC。
/// 与 solve_events 的 `time_idx < times.len()` 走法同构，保证对齐。
pub fn reconstruct_met_channels(sci_data: &SciFile, offset: f64) -> Vec<u16> {
    use super::eband::{wrap_channel, CHANNEL_SEC};
    let packet_times = reconstruct_with_wrap_tracking(sci_data, offset);
    let mut channels = Vec::new();
    for (pkt_idx, ccsds) in sci_data.ccsds.iter().enumerate() {
        let n = packet_times[pkt_idx].len();
        let mut slot = 0usize;
        for event in parse_events(ccsds) {
            match event {
                Pack::Event { channel, .. } => {
                    if slot < n {
                        channels.push(wrap_channel(channel));
                        slot += 1;
                    }
                }
                Pack::Second { .. } => {
                    if slot < n {
                        channels.push(CHANNEL_SEC);
                        slot += 1;
                    }
                }
                Pack::Error => {}
            }
        }
    }
    channels
}

/// 与 reconstruct_met_times 逐元素对齐的脉宽（pulinfo = byte[1]）序列。
///
/// pulinfo 是 NaI/CsI 甄别量（实际脉宽 = pulinfo/48 μs；NaI X 射线窗
/// [54,70]，来自 HXMT handbook）。走法与 reconstruct_met_channels 完全一致，
/// 保证与时间/channel 一一对应。SEC 槽记为 0。
pub fn reconstruct_met_pulse_widths(sci_data: &SciFile, offset: f64) -> Vec<u8> {
    let packet_times = reconstruct_with_wrap_tracking(sci_data, offset);
    let mut widths = Vec::new();
    for (pkt_idx, ccsds) in sci_data.ccsds.iter().enumerate() {
        let n = packet_times[pkt_idx].len();
        let mut slot = 0usize;
        for event in parse_events(ccsds) {
            match event {
                Pack::Event { raw_bytes, .. } => {
                    if slot < n {
                        widths.push(raw_bytes[1]);
                        slot += 1;
                    }
                }
                Pack::Second { .. } => {
                    if slot < n {
                        widths.push(0u8);
                        slot += 1;
                    }
                }
                Pack::Error => {}
            }
        }
    }
    widths
}

/// 单个事例的详细信息（用于 dump-events 调试）。
pub struct EventDetail {
    pub pkt_index: usize,
    pub evt_index: usize,
    pub is_second: bool,
    pub channel: u8,
    pub det_id: u8,
    pub met: f64,
    pub raw_bytes: [u8; 8],
}

/// 输出指定 MET 范围内的所有事例详情。
pub fn dump_event_details(
    sci_data: &SciFile,
    offset: f64,
    met_min: f64,
    met_max: f64,
) -> Vec<EventDetail> {
    solve_events(sci_data, offset, Some(met_min), Some(met_max))
}

/// 时间解算：返回所有事件的详细信息。
/// met_min/met_max 为 None 时不做时间窗口过滤。
pub fn solve_events(
    sci_data: &SciFile,
    offset: f64,
    met_min: Option<f64>,
    met_max: Option<f64>,
) -> Vec<EventDetail> {
    let packet_times = reconstruct_with_wrap_tracking(sci_data, offset);
    let mut result = Vec::new();
    let lo = met_min.unwrap_or(f64::NEG_INFINITY);
    let hi = met_max.unwrap_or(f64::INFINITY);

    for (pkt_idx, ccsds) in sci_data.ccsds.iter().enumerate() {
        let events = parse_events(ccsds);
        let times = &packet_times[pkt_idx];

        let mut time_idx = 0;
        for (evt_idx, event) in events.iter().enumerate() {
            match event {
                Pack::Second { .. } => {
                    if time_idx < times.len() {
                        let computed_met = times[time_idx];
                        time_idx += 1;
                        if computed_met >= lo && computed_met <= hi {
                            result.push(EventDetail {
                                pkt_index: pkt_idx,
                                evt_index: evt_idx,
                                is_second: true,
                                channel: 0,
                                det_id: 0,
                                met: computed_met,
                                raw_bytes: [0; 8],
                            });
                        }
                    }
                }
                Pack::Event {
                    channel,
                    raw_bytes,
                    ..
                } => {
                    if time_idx < times.len() {
                        let computed_met = times[time_idx];
                        time_idx += 1;
                        if computed_met >= lo && computed_met <= hi {
                            result.push(EventDetail {
                                pkt_index: pkt_idx,
                                evt_index: evt_idx,
                                is_second: false,
                                channel: *channel,
                                det_id: (raw_bytes[4] >> 1) & 0x07,
                                met: computed_met,
                                raw_bytes: *raw_bytes,
                            });
                        }
                    }
                }
                Pack::Error => {}
            }
        }
    }

    result
}

/// 单个 CCSDS 包的诊断信息。
pub struct PacketDiag {
    pub pkt_index: usize,
    pub n_event: usize,
    pub n_second: usize,
    pub n_error: usize,
    pub n_second_valid: usize,
    pub n_output: usize,
    pub n_dropped: usize,
    pub has_anchor: bool,
    pub utc_tail: f64,
    pub met_min: Option<f64>,
    pub met_max: Option<f64>,
    pub n_0x5a: usize,
}

/// 对每个 CCSDS 包进行诊断，返回统计信息。
pub fn diagnose_packets(sci_data: &SciFile, offset: f64) -> Vec<PacketDiag> {
    let packet_times = reconstruct_with_wrap_tracking(sci_data, offset);
    let mut result = Vec::new();
    let mut has_anchor = false;

    for (pkt_idx, ccsds) in sci_data.ccsds.iter().enumerate() {
        let utc_tail = get_utc_tail(ccsds);
        let events = parse_events(ccsds);
        let times = &packet_times[pkt_idx];

        // Count occurrences of 0x5A in byte[0] of each 8-byte slot (FIFO start marker leaking
        // into payload — indicates MCU read pointer misalignment from SEE).
        let payload = &ccsds[6..878];
        let n_0x5a = payload
            .chunks_exact(8)
            .filter(|chunk| chunk[0] == 0x5A)
            .count();

        let mut n_event = 0usize;
        let mut n_second = 0usize;
        let mut n_error = 0usize;
        let mut n_second_valid = 0usize;

        for event in &events {
            match event {
                Pack::Second { stime, .. } => {
                    n_second += 1;
                    let met = *stime as f64 + offset;
                    if (met - utc_tail).abs() < 2.0 {
                        n_second_valid += 1;
                        has_anchor = true;
                    }
                }
                Pack::Event { .. } => {
                    n_event += 1;
                }
                Pack::Error => {
                    n_error += 1;
                }
            }
        }

        let n_output = times.len();
        let n_dropped = (n_event + n_second).saturating_sub(n_output);

        let pkt_met_min = times.iter().copied().reduce(f64::min);
        let pkt_met_max = times.iter().copied().reduce(f64::max);

        result.push(PacketDiag {
            pkt_index: pkt_idx,
            n_event,
            n_second,
            n_error,
            n_second_valid,
            n_output,
            n_dropped,
            has_anchor,
            utc_tail,
            met_min: pkt_met_min,
            met_max: pkt_met_max,
            n_0x5a,
        });
    }

    result
}

/// 搜索用饱和排除窗向 FIFO reset 空洞前后各扩展的余量（秒）。
///
/// TGF 搜索的本底用候选周围 1s 窗（±0.5s）估计；只要饱和落在该窗内，本底估计
/// 就被污染，附近所有候选的显著性都不可信。扩 ±1s 保证候选本底窗与饱和空洞
/// 完全分离（0.5s 严格下限之上留半秒余量，兼顾 FIFO 检测对饱和起止边缘的漏检）。
const SATURATION_PADDING_SECONDS: f64 = 1.0;

/// 扫描饱和区间，返回 MissionElapsedTime 类型（供 TGF 搜索排除使用）。
///
/// 在 `scan_saturation_intervals_raw` 给出的原始 FIFO reset 空洞基础上，每段
/// 向前后各扩 `SATURATION_PADDING_SECONDS`。相邻空洞扩窗后的合并（求并）由
/// 上层 `Chunk::get_saturation_intervals` 统一处理。
pub fn scan_saturation_intervals(
    sci_data: &SciFile,
    offset: f64,
) -> Vec<(MissionElapsedTime<HxmtHe>, MissionElapsedTime<HxmtHe>)> {
    scan_saturation_intervals_raw(sci_data, offset)
        .into_iter()
        .map(|(start, stop)| {
            (
                MissionElapsedTime::new(start - SATURATION_PADDING_SECONDS),
                MissionElapsedTime::new(stop + SATURATION_PADDING_SECONDS),
            )
        })
        .collect()
}

/// 扫描饱和区间，直接返回原始 FIFO reset 空洞的 MET 秒数（未扩窗，诊断用）。
pub fn scan_saturation_intervals_raw(sci_data: &SciFile, offset: f64) -> Vec<(f64, f64)> {
    detect_fifo_reset_intervals(sci_data, offset)
        .into_iter()
        .map(|iv| (iv.start_met, iv.stop_met))
        .collect()
}

/// 诊断：打印包信息。
pub fn print_diagnose_packets(sci_data: &SciFile, offset: f64, pkt_min: usize, pkt_max: usize) {
    let diags = diagnose_packets(sci_data, offset);
    for d in &diags {
        if d.pkt_index >= pkt_min && d.pkt_index <= pkt_max {
            println!(
                "pkt={} evt={} sec={} err={} sec_valid={} out={} drop={} anchor={} utc={:.0} met=[{}, {}]",
                d.pkt_index, d.n_event, d.n_second, d.n_error, d.n_second_valid,
                d.n_output, d.n_dropped, d.has_anchor, d.utc_tail,
                d.met_min.map_or("?".to_string(), |v| format!("{:.6}", v)),
                d.met_max.map_or("?".to_string(), |v| format!("{:.6}", v)),
            );
        }
    }
}

/// 打印 ptime/utc 诊断信息。
pub fn dump_ptime_utc(sci_data: &SciFile, offset: f64, pkt_min: usize, pkt_max: usize) {
    println!("pkt,evt_idx,type,ptime,stime,utc_tail,met");
    for (pkt_idx, ccsds) in sci_data.ccsds.iter().enumerate() {
        if pkt_idx < pkt_min || pkt_idx > pkt_max {
            continue;
        }
        let utc_tail = get_utc_tail(ccsds);
        let events = parse_events(ccsds);
        for (evt_idx, event) in events.iter().enumerate() {
            match event {
                Pack::Event { ptime, .. } => {
                    println!("{},{},EVT,{},,,{:.0}", pkt_idx, evt_idx, ptime, utc_tail);
                }
                Pack::Second { stime, ptime } => {
                    let met = *stime as f64 + offset;
                    println!("{},{},SEC,{},{},{:.0},{:.6}", pkt_idx, evt_idx, ptime, stime, utc_tail, met);
                }
                Pack::Error => {
                    println!("{},{},ERR,,,,", pkt_idx, evt_idx);
                }
            }
        }
    }
}

/// 诊断：对每个 CCSDS 包尝试 0~7 字节偏移，输出各偏移下的 CRC 通过数。
pub fn check_byte_offsets(sci_data: &SciFile, pkt_min: usize, pkt_max: usize) {
    println!("pkt,utc_tail,off0,off1,off2,off3,off4,off5,off6,off7");
    for (pkt_idx, ccsds) in sci_data.ccsds.iter().enumerate() {
        if pkt_idx < pkt_min || pkt_idx > pkt_max {
            continue;
        }
        let utc_tail = get_utc_tail(ccsds);
        let mut pass_counts = [0u32; 8];

        for byte_offset in 0..8usize {
            let start = 6 + byte_offset;
            let end = 878;
            if start >= end {
                continue;
            }
            let payload = &ccsds[start..end];
            for chunk in payload.chunks_exact(8) {
                let mut row = [0u64; 8];
                for (i, byte) in chunk.iter().enumerate() {
                    row[i] = *byte as u64;
                }
                if crc_check(&row) == row[7] & 0x0F {
                    pass_counts[byte_offset] += 1;
                }
            }
        }

        println!(
            "{},{:.0},{},{},{},{},{},{},{},{}",
            pkt_idx, utc_tail,
            pass_counts[0], pass_counts[1], pass_counts[2], pass_counts[3],
            pass_counts[4], pass_counts[5], pass_counts[6], pass_counts[7],
        );
    }
}
