"""饱和恢复光变的下游协方差装配（spec §7）。

读 blink `reconstruct` 的**事件流 CSV**（box,type,met,channel,pulse_width,pkt_idx,evt_idx）
与可选的 **gap 块表**（spec §13），按任意时间 binning / 能段输出均值与协方差。

当前实现（第一步）：均值 N_i + **对角**协方差
  - D_i：观测事件（EVT）计数 → 观测泊松（spec §4 项 I 的无 gap 退化）
  - U_i：填充事件（FILL_GAP）计数 → **不可约丢失涨落**（spec §4 项 IV），保证填充 bin
         误差棒 ≥ √(重建计数)，堵住"平滑重建假装很确定"的漏。
  对角方差 = D + U = bin 总计数。

**未含**（待读块表装配）：cross-ref 的插值非对角 + 标定 k 满协方差（§5），degenerate 秩-2
块（§6），gap 间共模（§4）。故当前 `var_diag` 只够画**逐 bin 误差棒**；功率谱/积分/χ²
需要完整协方差（off-diagonal），走后续的块装配 + §9 的 MC 闭环。
"""
from __future__ import annotations

import csv
from dataclasses import dataclass

import numpy as np


@dataclass
class LcCov:
    """一条恢复光变的均值与（当前仅对角）协方差。"""

    edges: np.ndarray  # bin 边界
    N: np.ndarray      # 均值：bin 内所有事件（EVT + FILL_GAP）
    D: np.ndarray      # 观测泊松对角（EVT 计数）
    U: np.ndarray      # 填充不可约泊松地板（FILL_GAP 计数，spec §4 项 IV）

    @property
    def var_diag(self) -> np.ndarray:
        """对角方差 = D + U。误差棒 = sqrt(var_diag)。仅对角，见模块 docstring。"""
        return self.D + self.U

    @property
    def err(self) -> np.ndarray:
        return np.sqrt(self.var_diag)


def load_events(path: str) -> dict:
    """读事件流 CSV → dict of numpy arrays（box/type/met/channel）。"""
    boxes, types, mets, chans = [], [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            boxes.append(row["box"])
            types.append(row["type"])
            mets.append(float(row["met"]))
            chans.append(int(row["channel"]))
    return {
        "box": np.array(boxes),
        "type": np.array(types),
        "met": np.array(mets, dtype=float),
        "channel": np.array(chans, dtype=int),
    }


def mean_and_diag(events: dict, bin_edges, box=None, chan_range=None) -> LcCov:
    """按 bin_edges 分箱,返回均值 N 与对角分量 D（观测）/U（填充）。
    box=None → 全 HE 总光变（三盒相加）；chan_range=(lo,hi) 闭区间能段选择。"""
    edges = np.asarray(bin_edges, dtype=float)
    met = events["met"]
    sel = np.ones(len(met), dtype=bool)
    if box is not None:
        sel &= events["box"] == box
    if chan_range is not None:
        lo, hi = chan_range
        sel &= (events["channel"] >= lo) & (events["channel"] <= hi)
    is_evt = events["type"] == "EVT"
    is_fill = events["type"] == "FILL_GAP"
    n_all, _ = np.histogram(met[sel], bins=edges)
    d_obs, _ = np.histogram(met[sel & is_evt], bins=edges)
    u_fill, _ = np.histogram(met[sel & is_fill], bins=edges)
    return LcCov(
        edges=edges,
        N=n_all.astype(float),
        D=d_obs.astype(float),
        U=u_fill.astype(float),
    )
