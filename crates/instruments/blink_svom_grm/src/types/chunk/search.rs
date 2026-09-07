use crate::types::Chunk;
use crate::types::Event;
use crate::types::SvomGrm;
use blink_algorithms::detector_share::{MAX_DETECTOR_FRACTION, max_detector_fraction};
use blink_algorithms::snapshot_stepping::SearchConfig;
use blink_algorithms::snapshot_stepping::search_new;
use blink_core::types::Attitude;
use blink_core::types::MissionElapsedTime;
use blink_core::types::Position;
use blink_core::types::Signal;
use blink_core::types::Trajectory;
use blink_core::traits::Event as _;
use std::sync::atomic::Ordering;
use uom::si::f64::*;

pub(super) fn search(chunk: &Chunk) -> Vec<Signal<Event>> {
    // 事例准入见 `Event::keep`。不过滤的话一小时会冒出 15–30 个软谱噪声
    // 候选（实测），它们会淹没真信号。
    //
    // 再按 GTI 过滤。曝光核算一直走的是 GTI，事例流却没有——于是一直在搜
    // L1B 判为坏时段的数据：每次 SAA 穿越的缺口两端各有约 0.5 s 事例留在表里，
    // 两侧的本底窗都有一半落在空区或正常区，均值被压低，普通计数就显得极显著。
    // 它们就是全量结果里的「1 Hz 假信号」（缺口离开侧，相位锁在 0.505 s 只是
    // 整秒边界减 0.495 s 的投影）和「高本底」（缺口进入侧）两类。只筛 FLAG 不够：
    // 进入侧那半秒 FLAG=0。丢掉的事例计数，见 `diagnostics`。
    //
    // 过滤只修分子。分母（本底窗的时长）由 `search_new` 按同一份 GTI 夹到
    // 活时间里，否则紧挨缺口的好数据的本底窗仍会伸进空区——SAA 停机前是
    // 130–210 kc/s 的平台后硬切断，那里的候选本底会被压低一半。
    let mut n_outside_gti = 0usize;
    let mut events = chunk
        .evt_file
        .into_iter()
        .filter(|event| event.keep())
        .filter(|event| {
            let inside = chunk.evt_file.gti_contains(event.time().met());
            if !inside {
                n_outside_gti += 1;
            }
            inside
        })
        .collect::<Vec<_>>();
    // 三路各自内部会有几处微米级的时间抖动（1–5 个量化步），k 路归并把它
    // 原样传了出来，而 search_new 假定输入有序。幅度大到有害的回跳已在
    // `exclusion` 挡掉，剩下的排一下就干净了——数组本来就基本有序，代价很低。
    events.sort();
    let gti: Vec<[MissionElapsedTime<SvomGrm>; 2]> = chunk
        .evt_file
        .gti_segments()
        .into_iter()
        .map(|(start, stop)| [MissionElapsedTime::new(start), MissionElapsedTime::new(stop)])
        .collect();
    let results = search_new(
        &events,
        1,
        chunk.span[0],
        chunk.span[1],
        &gti,
        SearchConfig {
            min_duration: Time::new::<uom::si::time::microsecond>(0.0),
            // 窗长上限 1 ms。试过 5 ms（覆盖期全量重跑，见 OPEN-QUESTIONS 第 17 条）：
            // 候选多出 47%，但那些新增候选的闪电关联数是 1、偶然期望 1.29，即纯本底；
            // 82 个已证实 TGF 一个不多一个不少。已证实样本的最佳窗中位只有 0.14 ms，
            // 宽窗量出的 T90 都在 0.43 ms 以内，1 ms 对探测足够。
            max_duration: Time::new::<uom::si::time::millisecond>(1.0),
            neighbor: Time::new::<uom::si::time::second>(1.0),
            hollow: Time::new::<uom::si::time::millisecond>(10.0),
            false_positive_per_year: 20.0,
            min_number: 8,
            // 单组：三个 GRD 合成一路，符合判据退化
            coincidence: 1,
        },
    );

    // 两条轨道各建一次，而不是每个候选重建一遍（1 Hz 采样，一小时约 3700 点）。
    let attitudes = Trajectory::<MissionElapsedTime<SvomGrm>, Attitude>::from(&chunk.att_file);
    let positions = Trajectory::<MissionElapsedTime<SvomGrm>, Position>::from(&chunk.orb_file);

    let mut n_dropped = 0usize;
    let mut n_without_attitude = 0usize;
    let mut n_single_detector = 0usize;
    let signals = results
        .into_iter()
        .filter_map(|candidate| {
            // 单路毛刺否决：候选窗里一路探测器占了绝大多数计数就不是暴发。
            // 三个 GRD 朝向不同，但 v3 全量里占比 > 0.8 的 12 个覆盖内候选闪电关联 0 个，
            // 占比 ≤ 0.6 的关联 42%。按整窗算而不按最显著一格，见 `detector_share`。
            if max_detector_fraction(&events, candidate.start, candidate.stop, |e| e.detector_id)
                > MAX_DETECTOR_FRACTION
            {
                n_single_detector += 1;
                return None;
            }
            let peak = candidate.start + candidate.bin_size_best / 2.0;
            // GRM 的 att/orb 是逐小时文件，尾端比事例流早收约 8 s，那一截里的
            // 候选取不到星历。丢可以，静默丢不行——记账见 `diagnostics`。
            let Some(position) = positions.interpolate(peak) else {
                n_dropped += 1;
                return None;
            };
            // 姿态只是元数据，缺了不丢候选，留空并记账
            let attitude = attitudes.interpolate(peak).map(|a| a.state);
            if attitude.is_none() {
                n_without_attitude += 1;
            }
            Some(Signal {
                start: candidate.start,
                stop: candidate.stop,
                bin_size_min: candidate.bin_size_min,
                bin_size_max: candidate.bin_size_max,
                bin_size_best: candidate.bin_size_best,
                delay: candidate.delay,
                count: candidate.count,
                mean: candidate.mean,
                sf: candidate.sf(),
                false_positive_per_year: candidate.false_positive_per_year(),
                attitude,
                position: position.state,
                // 事例表里有一列 ANTI_COIN，但它的语义未经确认，尚未接入
                // —— 见本 crate 的 `OPEN-QUESTIONS.md`。
                acd: None,
            })
        })
        .collect::<Vec<_>>();

    chunk
        .dropped_no_ephemeris
        .store(n_dropped, Ordering::Relaxed);
    chunk
        .without_attitude
        .store(n_without_attitude, Ordering::Relaxed);
    chunk
        .dropped_single_detector
        .store(n_single_detector, Ordering::Relaxed);
    chunk
        .events_outside_gti
        .store(n_outside_gti, Ordering::Relaxed);

    signals
}
