//! 天格计划（GRID）立方星的 TGF 搜索接入。
//!
//! 四颗有事例产品的星（GRID-02 / 03B / 04 / 07）探测器相同、归档布局相同，
//! 只有目录名、文件名前缀和数据起始日不同，所以实现只写一份，用类型参数
//! [`types::Satellite`] 区分。归档里核实过的事实与尚未定标的判据见本 crate 的
//! `OPEN-QUESTIONS.md`。
pub mod io;
pub mod types;
