use crate::types::Lightning;
use chrono::prelude::*;
use rusqlite::{Connection, params};
use std::env;

thread_local! {
    // 每个线程持有独立的只读连接：SQLite 允许多读者并发，用线程本地连接（而非
    // 全局 Mutex<Connection>）让 filter 的百万级查询能真正并行，而不是串行等锁。
    static LIGHTNING_CONNECTION: Connection = open_connection();
}

fn open_connection() -> Connection {
    let path = env::var("WWLLN_DB_PATH")
        .unwrap_or_else(|_| String::from("/Volumes/Graphite/WWLLN/WWLLN.db"));
    let conn =
        Connection::open_with_flags(path, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY).unwrap();
    // Set a longer busy timeout (e.g., 30 seconds = 30000 ms)
    conn.busy_timeout(std::time::Duration::from_secs(30))
        .unwrap();
    // 用 mmap 直接映射数据库读取，绕过 SQLite 全局页缓存（pcache1）那把互斥锁——
    // filter 用几十个线程各开一个连接跑百万级查询时，该锁的争用是主要瓶颈
    // （perf 显示 ~40% 时间在 pthread_mutex_lock/pcache1）。映射整库（>422GB），
    // 实际驻留由 OS 页缓存按热度管理。
    conn.pragma_update(None, "mmap_size", 549_755_813_888i64)
        .unwrap();
    conn
}

/// 库里最早与最晚一条闪电的时刻。
///
/// 用来把"这一刻库里没有闪电"和"这一刻库根本没有数据"分开——库目前止于
/// 2024-12-31，而 SVOM 的候选到 2026 年，落在覆盖外的候选查出来必然是空，
/// 那不是没有闪电。
///
/// 走 `ORDER BY ... LIMIT 1` 而不是 `min()`/`max()`：后者在这个 422 GB 的表上
/// 进不了索引，实测跑 200 s 还没回来；前者是纯索引扫描，几十毫秒。
pub fn coverage() -> (DateTime<Utc>, DateTime<Utc>) {
    LIGHTNING_CONNECTION.with(|connection| {
        let bound = |order: &str| -> DateTime<Utc> {
            let text: String = connection
                .query_row(
                    &format!("SELECT time FROM lightning ORDER BY time {order} LIMIT 1"),
                    [],
                    |row| row.get(0),
                )
                .expect("WWLLN database is empty");
            NaiveDateTime::parse_from_str(&text, "%Y-%m-%d %H:%M:%S%.6f")
                .expect("unparseable time in WWLLN database")
                .and_utc()
        };
        (bound("ASC"), bound("DESC"))
    })
}

pub fn get_lightnings(time_start: DateTime<Utc>, time_end: DateTime<Utc>) -> Vec<Lightning> {
    let time_start_str = time_start.format("%Y-%m-%d %H:%M:%S%.6f").to_string();
    let time_end_str = time_end.format("%Y-%m-%d %H:%M:%S%.6f").to_string();
    LIGHTNING_CONNECTION.with(|connection| {
    let mut statement = connection
        .prepare(
            "
                SELECT
                    time,
                    lat,
                    lon,
                    resid,
                    nstn,
                    energy,
                    energy_uncertainty,
                    estn
                FROM
                    lightning
                WHERE
                    time BETWEEN ?1 AND ?2
                ORDER BY time ASC
                ",
        )
        .unwrap();
    statement
        .query_map(params![time_start_str, time_end_str], |row| {
            Ok(Lightning {
                time: NaiveDateTime::parse_from_str(
                    &row.get::<_, String>(0).unwrap(),
                    "%Y-%m-%d %H:%M:%S%.6f",
                )
                .unwrap()
                .and_utc(),
                lat: row.get::<_, f64>(1).unwrap(),
                lon: row.get::<_, f64>(2).unwrap(),
                resid: row.get::<_, f64>(3).unwrap(),
                nstn: row.get::<_, i64>(4).unwrap() as u32,
                energy: row.get::<_, Option<f64>>(5).unwrap(),
                energy_uncertainty: row.get::<_, Option<f64>>(6).unwrap(),
                estn: row.get::<_, Option<i64>>(7).unwrap().map(|x| x as u32),
            })
        })
        .unwrap()
        .map(|x| x.unwrap())
        .collect::<Vec<_>>()
    })
}
