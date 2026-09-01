//! 目录生成：tgfs.json → 池级清洁（REP 列车摘除）→ 论文判选 → 目录 CSV。
//!
//! 阶段顺序有方法论含义：列车摘除是数据清洁，位于一切统计之前（对应旧搜索
//! 中 continuous() 簇否决所在的位置）——fp-distribution 与阈值论证都画在
//! 清洁池上；判选公式维持论文原样，不因新增判据而改动。被摘候选仍完整保留
//! 在 tgfs.json 中（带 is_train 标记），摘除显式可逆。
//!
//! 逐步计数打印到 stderr——每一步扔掉多少一目了然，无静默过滤。

use blink_wwlln::Tgf;
use std::io::BufReader;
use std::path::Path;

/// 判选 A 层（直接显著）：fa ≤ 1e-5。依据 = 清洁池 fp-distribution 上
/// 噪声/TGF 双幂律交点（v6 实测 4.2e-5）右侧留保守余量；论文原值不动。
const SELECT_FA_DIRECT: f64 = 1e-5;
/// 判选 B 层（闪电救援）：fa ≤ 1 且 WWLLN 关联。期望误关联（Σ coincidence
/// probability）按统计整体报告，不逐候选收紧——那会系统性偏压雷暴活跃区。
const SELECT_FA_ASSOC: f64 = 1.0;

/// 候选归入的判选层。
enum Tier {
    /// A 层：fa ≤ 1e-5，自身显著性足够
    Direct,
    /// B 层：1e-5 < fa ≤ 1，靠闪电关联救援
    Lightning,
}

fn select(tgf: &Tgf) -> Option<Tier> {
    let fa = tgf.signal.false_positive_per_year;
    if fa <= SELECT_FA_DIRECT {
        Some(Tier::Direct)
    } else if fa <= SELECT_FA_ASSOC && tgf.lightning.associated {
        Some(Tier::Lightning)
    } else {
        None
    }
}

/// 与 tgfs.json 内串一致的时间字符串（serde 同一格式，供下游按 start 精确 join）。
fn time_string(t: &chrono::DateTime<chrono::Utc>) -> String {
    serde_json::to_string(t)
        .expect("serialize datetime")
        .trim_matches('"')
        .to_string()
}

pub fn cmd_catalog(input: &Path, out: &Path) {
    let file = std::fs::File::open(input)
        .unwrap_or_else(|e| panic!("failed to open {}: {e}", input.display()));
    eprintln!("catalog: reading {} ...", input.display());
    let mut tgfs: Vec<Tgf> =
        serde_json::from_reader(BufReader::new(file)).expect("failed to parse tgfs json");
    let n_pool = tgfs.len();

    // 池级清洁：REP 微暴列车整体摘除
    tgfs.retain(|t| !t.train.is_train);
    let n_clean = tgfs.len();
    eprintln!(
        "catalog: pool {n_pool} -> train members removed {} -> clean pool {n_clean}",
        n_pool - n_clean
    );

    // 判选（论文公式），防御性按 start 排序保持年表序
    tgfs.sort_by_key(|t| t.signal.start);
    let mut rows = String::from(
        "date,start,stop,duration,count,false_positive_per_year,\
         longitude,latitude,altitude,tier,associated,coincidence_probability,\
         neighbors_10min,n,n_acd,n_acd_multi,n_bg,n_acd_bg\n",
    );
    let (mut n_direct, mut n_lightning) = (0usize, 0usize);
    // 救援带内任何候选只要偶然关联上就会被误救进目录，故期望误救 =
    // 整个带（不只被选中者）的 Σ coincidence_probability。
    let mut expected_misassoc = 0.0f64;
    for tgf in &tgfs {
        let fa = tgf.signal.false_positive_per_year;
        if fa > SELECT_FA_DIRECT && fa <= SELECT_FA_ASSOC {
            // 覆盖外的候选没有偶合概率，也不可能被关联救进来，不计入期望误救。
            expected_misassoc += tgf.lightning.coincidence_probability.unwrap_or(0.0);
        }
        let tier = match select(tgf) {
            Some(Tier::Direct) => {
                n_direct += 1;
                "direct"
            }
            Some(Tier::Lightning) => {
                n_lightning += 1;
                "lightning"
            }
            None => continue,
        };
        let s = &tgf.signal;
        let duration = (s.stop - s.start).num_nanoseconds().expect("sub-ns range") as f64 / 1e9;
        let acd = match &s.acd {
            Some(a) => format!(
                "{},{},{},{},{}",
                a.n, a.n_acd, a.n_acd_multi, a.n_bg, a.n_acd_bg
            ),
            None => ",,,,".to_string(),
        };
        rows.push_str(&format!(
            "{},{},{},{:.9},{},{:e},{:.4},{:.4},{:.1},{},{},{},{},{}\n",
            s.start.format("%Y%m%d"),
            time_string(&s.start),
            time_string(&s.stop),
            duration,
            s.count,
            s.false_positive_per_year,
            s.position.longitude,
            s.position.latitude,
            s.position.altitude.get::<uom::si::length::meter>(),
            tier,
            u8::from(tgf.lightning.associated),
            tgf.lightning
                .coincidence_probability
                .map(|p| format!("{p:e}"))
                .unwrap_or_default(),
            tgf.train.neighbors_10min,
            acd,
        ));
    }
    std::fs::write(out, rows).expect("failed to write catalog csv");
    eprintln!(
        "catalog: selected {} = direct(fa<={SELECT_FA_DIRECT:e}) {n_direct} \
         + lightning(fa<={SELECT_FA_ASSOC:e} & assoc) {n_lightning} -> {}",
        n_direct + n_lightning,
        out.display()
    );
    eprintln!(
        "catalog: expected false rescues (sum coinc over the {:e} < fa <= {:e} band): \
         {expected_misassoc:.1}",
        SELECT_FA_DIRECT, SELECT_FA_ASSOC
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 一条最小完整记录；只有测试关心的字段可变。
    fn record(start: &str, fa: f64, associated: bool, is_train: bool) -> String {
        format!(
            r#"{{
              "signal": {{
                "start": "{start}",
                "stop": "{}",
                "bin_size_min": 0.0001, "bin_size_max": 0.001, "bin_size_best": 0.0005,
                "delay": 0.0, "count": 10, "mean": 0.5, "sf": 1e-10,
                "false_positive_per_year": {fa},
                "attitude": {{ "q1": 0.0, "q2": 0.0, "q3": 0.0 }},
                "position": {{ "longitude": 100.0, "latitude": 1.0, "altitude": 540000.0 }},
                "instrument": "Insight-HXMT/HE",
                "acd": {{ "n": 10, "n_acd": 1, "n_acd_multi": 0, "n_bg": 5000, "n_acd_bg": 400 }}
              }},
              "lightning": {{ "associated": {associated}, "coincidence_probability": 0.001 }},
              "train": {{ "neighbors_10min": {}, "is_train": {is_train} }}
            }}"#,
            start.replace(".000", ".001"),
            if is_train { 500 } else { 3 },
        )
    }

    #[test]
    fn pool_cleaning_then_paper_selection() {
        let json = format!(
            "[{},{},{},{},{}]",
            // 列车成员：无论显著性多高都在池级被摘
            record("2024-10-09T12:00:00.000000000Z", 1e-12, false, true),
            // A 层：直接显著
            record("2024-10-09T13:00:00.000000000Z", 1e-9, false, false),
            // B 层：弱显著但闪电关联
            record("2024-10-09T14:00:00.000000000Z", 0.5, true, false),
            // 弃：弱显著且无关联
            record("2024-10-09T15:00:00.000000000Z", 0.5, false, false),
            // 弃：fa 超出救援上限，关联也不救
            record("2024-10-09T16:00:00.000000000Z", 15.0, true, false),
        );
        let dir = std::env::temp_dir().join(format!("blink_catalog_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let input = dir.join("tgfs.json");
        let out = dir.join("catalog.csv");
        std::fs::write(&input, json).unwrap();

        cmd_catalog(&input, &out);

        let csv = std::fs::read_to_string(&out).unwrap();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 3, "header + 2 selected: {csv}");
        assert!(lines[1].contains("2024-10-09T13:00:00Z") && lines[1].contains(",direct,"));
        assert!(lines[2].contains("2024-10-09T14:00:00Z") && lines[2].contains(",lightning,"));
        // 列车成员即使 fa=1e-12 也不得出现
        assert!(!csv.contains("12:00:00"));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
