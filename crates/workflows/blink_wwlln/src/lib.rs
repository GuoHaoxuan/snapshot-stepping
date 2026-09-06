use blink_core::traits::Instrument;
use blink_core::types::{TemporalState, UnifiedSignal};
use blink_lightning::{algorithms::coincidence_prob, database::coverage, database::get_lightnings};
use blink_load::load_all;
use chrono::TimeDelta;
use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicUsize, Ordering};
use uom::si::f64::*;

/// 列车窗口半宽：±600 s 覆盖一次辐射带沉降区穿越（REP 微暴列车持续几分钟
/// 到几十分钟）；单个雷暴单体在星下可见仅 1–2 分钟，又不至于把同一轨道
/// 前后两个独立雷暴缝在一起。
pub const TRAIN_HALF_WINDOW_US: i64 = 600 * 1_000_000;

/// 列车判定阈（已冻结）：WWLLN 关联认证 TGF（2547 个）的 neighbors_10min
/// 分布 99 分位。阈值只由 TGF 样本单方面定标，REP 样本不参与；认证 REP 人群
/// 从 ~60 起步（中位 656），阈值坐在两人群之间的空沟里，取 30 还是 50 结果
/// 几乎不变。原计划的 ACD 逐事例符合交叉验证已做且判据不成立——HE 里的 REP
/// 信号是沉降电子的轫致辐射光子（NaI 截止层几乎无信号，ACD 符合超出 <1%），
/// 与 TGF 伽马对 ACD 同盲——故本阈值以 TGF 侧定标独立成立。
pub const TRAIN_THRESHOLD: u32 = 34;

#[derive(Serialize, Deserialize)]
pub struct LightningInfo {
    pub associated: bool,
    /// 偶合概率。覆盖外为 `None`——那里既没查也算不出，不是 0 也不是 1。
    pub coincidence_probability: Option<f64>,
    /// 候选时刻是否落在 WWLLN 库的覆盖区间内。
    ///
    /// 库止于 2024-12-31，而 SVOM 的候选到 2026-08——覆盖外的候选查出来必然
    /// 是空，`associated` 于是恒为 false。那是"查不到"，不是"没有闪电"，
    /// 下游不能把两者当成一回事。旧的 tgfs.json 没有这一列，读回时按 true
    /// 处理：HXMT 的搜索范围整个落在覆盖内。
    #[serde(default = "covered_by_default")]
    pub in_coverage: bool,
}

fn covered_by_default() -> bool {
    true
}

#[derive(Serialize, Deserialize)]
pub struct TrainInfo {
    /// 全初始候选池（fa < 20 /yr，含亚阈）中，start 落在本候选 start 的
    /// 半开窗 [t−600s, t+600s) 内的其他候选数
    pub neighbors_10min: u32,
    /// neighbors_10min > 34：REP 微暴列车成员，目录阶段池级摘除
    pub is_train: bool,
}

/// tgfs.json 的单条记录。`blink wwlln` 写出（富集不筛选），
/// `blink catalog` 读回做池级清洁与判选。
#[derive(Serialize, Deserialize)]
pub struct Tgf {
    pub signal: UnifiedSignal,
    pub lightning: LightningInfo,
    pub train: TrainInfo,
}

/// 逐候选数出半开窗 [t−600s, t+600s) 内的其他时刻数（排除自身）。
/// 与原型 np.searchsorted 的左侧语义逐点一致。
fn neighbor_counts_us(times_us: &[i64]) -> Vec<u32> {
    let mut sorted = times_us.to_vec();
    sorted.sort_unstable();
    times_us
        .iter()
        .map(|&t| {
            let lo = sorted.partition_point(|&x| x < t - TRAIN_HALF_WINDOW_US);
            let hi = sorted.partition_point(|&x| x < t + TRAIN_HALF_WINDOW_US);
            (hi - lo - 1) as u32
        })
        .collect()
}

/// 列车密度特征。数的是**全量初始候选池**而非仅显著候选：REP 列车过境会
/// 点亮上千个亚阈候选（2025-09-30 单趟 ±10 min 内 4491 个），宽池里列车
/// 藏不住；TGF 一次过境下方只有一个雷暴系统，连亚阈算上也只有几个到
/// 二十几个（中位 8）。
fn train_neighbor_counts(signals: &[UnifiedSignal]) -> Vec<u32> {
    let times: Vec<i64> = signals
        .iter()
        .map(|s| s.start.timestamp_micros())
        .collect();
    neighbor_counts_us(&times)
}

/// 对单个候选做 WWLLN 闪电关联 + 虚警概率。每次调用的两个 `get_lightnings`
/// 查询走线程本地只读连接（见 blink_lightning::database），可安全并行。
fn associate(
    signal: &UnifiedSignal,
    neighbors_10min: u32,
    coverage: (chrono::DateTime<chrono::Utc>, chrono::DateTime<chrono::Utc>),
    window: TimeDelta,
) -> Tgf {
    let peak_time = signal.peak_time();
    let in_coverage = peak_time >= coverage.0 && peak_time <= coverage.1;
    if !in_coverage {
        // 覆盖外就不查了：查也是空，白付两次百万级检索的代价。
        return Tgf {
            signal: signal.clone(),
            lightning: LightningInfo {
                associated: false,
                coincidence_probability: None,
                in_coverage: false,
            },
            train: TrainInfo {
                neighbors_10min,
                is_train: neighbors_10min > TRAIN_THRESHOLD,
            },
        };
    }
    let position = TemporalState {
        timestamp: peak_time,
        state: signal.position.clone(),
    };
    let lightnings = get_lightnings(
        peak_time - TimeDelta::seconds(1),
        peak_time + TimeDelta::seconds(1),
    )
    .into_iter()
    .filter(|lightning| {
        lightning.is_associated(
            &position,
            window,
            Length::new::<uom::si::length::kilometer>(ASSOCIATION_KILOMETERS),
        )
    })
    .collect::<Vec<_>>();

    Tgf {
        signal: signal.clone(),
        lightning: LightningInfo {
            associated: !lightnings.is_empty(),
            in_coverage: true,
            coincidence_probability: Some(coincidence_prob(
                &position,
                window,
                Length::new::<uom::si::length::kilometer>(ASSOCIATION_KILOMETERS),
                TimeDelta::minutes(2),
            )),
        },
        train: TrainInfo {
            neighbors_10min,
            is_train: neighbors_10min > TRAIN_THRESHOLD,
        },
    }
}

/// 关联窗口。±5 ms 与 800 km 是 HXMT 上定标的：800 km 约合 550 km 轨道高度
/// 下的可见地平，5 ms 覆盖光行时差与两边的计时不确定度。SVOM 轨道高约
/// 625 km、GBM 约 535 km，几何上同量级，先沿用同一组值；哪颗星要单独定标，
/// 得先有它自己的认证样本。
/// 默认的时间窗半宽；`run` 可按需放宽（如电子束的足点检验）。
pub const ASSOCIATION_MILLISECONDS: i64 = 5;
const ASSOCIATION_KILOMETERS: f64 = 800.0;

pub fn run<I: Instrument>(window_ms: i64) {
    let window = TimeDelta::milliseconds(window_ms);
    if window_ms != ASSOCIATION_MILLISECONDS {
        eprintln!("filter: association window ±{window_ms} ms (default {ASSOCIATION_MILLISECONDS})");
    }
    let signals = load_all::<I>();
    let total = signals.len();
    eprintln!("filter: {} candidates to associate ({})", total, I::name());

    // 列车密度在全量池上一次算完（排序 + 二分，O(n log n)），并行阶段按
    // 下标取用即可。
    let coverage = coverage();
    eprintln!(
        "filter: WWLLN coverage {} .. {}",
        coverage.0.format("%Y-%m-%d"),
        coverage.1.format("%Y-%m-%d")
    );
    let n_outside = signals
        .iter()
        .filter(|s| s.peak_time() < coverage.0 || s.peak_time() > coverage.1)
        .count();
    if n_outside > 0 {
        eprintln!(
            "filter: {n_outside}/{total} candidates fall outside it and are left unqueried \
             (in_coverage = false, not `no lightning`)"
        );
    }

    let neighbor_counts = train_neighbor_counts(&signals);
    let n_train = neighbor_counts
        .iter()
        .filter(|&&n| n > TRAIN_THRESHOLD)
        .count();
    eprintln!("filter: {n_train}/{total} flagged as train members (neighbors > {TRAIN_THRESHOLD})");

    // 每候选做 2 次 WWLLN 查询（±1s 关联 + ±62s 虚警概率），成本随该时段闪电
    // 密度差几十倍（活跃季 ±62s 窗返回上万条闪电）。静态分块会严重失衡（空段线程
    // 早退、忙段线程拖尾），故用原子取号做工作窃取：每线程反复领下一个待处理下标，
    // 忙闲自动均衡，56 核吃满到最后。结果带原下标收回后排序，保持原顺序。
    let n_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8);
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let signals_ref = &signals;
    let counts_ref = &neighbor_counts;

    let mut collected: Vec<(usize, Tgf)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..n_threads)
            .map(|_| {
                let next = &next;
                let done = &done;
                scope.spawn(move || {
                    let mut local: Vec<(usize, Tgf)> = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= total {
                            break;
                        }
                        local.push((i, associate(&signals_ref[i], counts_ref[i], coverage, window)));
                        let n = done.fetch_add(1, Ordering::Relaxed) + 1;
                        if n % 100_000 == 0 {
                            eprintln!("filter: {n}/{total}");
                        }
                    }
                    local
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect()
    });

    collected.sort_by_key(|(i, _)| *i);
    let tgfs: Vec<Tgf> = collected.into_iter().map(|(_, tgf)| tgf).collect();

    eprintln!("filter: {total}/{total} associated, writing tgfs.json");
    let json = serde_json::to_string_pretty(&tgfs).expect("failed to serialize to json");
    // 原子写：先写临时文件再 rename，避免下游（pipeline 的 cp / git）读到半截 json。
    let tmp = format!("tgfs.json.{}.tmp", nanoid::nanoid!(6));
    std::fs::write(&tmp, json).expect("failed to write tgfs.json tmp");
    std::fs::rename(&tmp, "tgfs.json").expect("failed to rename tgfs.json");
}

#[cfg(test)]
mod tests {
    use super::*;

    const S: i64 = 1_000_000; // 1 s in µs

    #[test]
    fn isolated_candidates_count_zero() {
        assert_eq!(neighbor_counts_us(&[]), Vec::<u32>::new());
        assert_eq!(neighbor_counts_us(&[42 * S]), vec![0]);
        // 相距远超 ±600s 的两个候选互不计数
        assert_eq!(neighbor_counts_us(&[0, 2000 * S]), vec![0, 0]);
    }

    #[test]
    fn train_members_count_each_other() {
        // 1 s 内的三个候选各自数到另外两个，自身不计
        let times = [0, S / 2, S];
        assert_eq!(neighbor_counts_us(&times), vec![2, 2, 2]);
        // 输入乱序不影响结果（计数走排序副本）
        let times = [S, 0, S / 2];
        assert_eq!(neighbor_counts_us(&times), vec![2, 2, 2]);
    }

    #[test]
    fn window_is_half_open_like_searchsorted() {
        // [t−600s, t+600s)：右端开、左端闭，与原型 np.searchsorted 逐点一致。
        // 对 t=0，+600s 处的邻居落在窗外；对 t=600s，0 处的邻居落在窗内。
        let times = [0, 600 * S];
        assert_eq!(neighbor_counts_us(&times), vec![0, 1]);
    }

    #[test]
    fn duplicated_timestamps_are_neighbors_not_self() {
        // 同一时刻的多条候选互为邻居，但都不数自己
        let times = [7 * S, 7 * S, 7 * S];
        assert_eq!(neighbor_counts_us(&times), vec![2, 2, 2]);
    }
}
