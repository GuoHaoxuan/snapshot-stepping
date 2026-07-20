"""Tests for analysis/recovery_cov.py (spec §7 下游装配:均值 + 对角 D/U)。"""
from __future__ import annotations

import numpy as np

import recovery_cov as R

HEADER = "box,type,met,channel,pulse_width,pkt_idx,evt_idx"


def _write(tmp_path, rows):
    p = tmp_path / "ev.csv"
    p.write_text(HEADER + "\n" + "\n".join(rows) + "\n")
    return str(p)


def test_mean_counts_all_events(tmp_path):
    """均值 = 数 bin 内所有事件(EVT + FILL_GAP);对角 = 观测 D + 填充 U。"""
    rows = [
        "A,EVT,0.1,50,60,-1,-1",
        "A,EVT,0.2,50,60,-1,-1",
        "A,FILL_GAP,1.1,50,60,-1,-1",
    ]
    ev = R.load_events(_write(tmp_path, rows))
    edges = np.array([0.0, 1.0, 2.0])
    lc = R.mean_and_diag(ev, edges)
    assert list(lc.N) == [2, 1], "均值应数所有事件"
    assert list(lc.D) == [2, 0], "bin0 两观测 → D=2;bin1 无观测"
    assert list(lc.U) == [0, 1], "bin1 一填充 → U=1(泊松地板)"
    assert list(lc.var_diag) == [2, 1], "对角方差 = D + U"


def test_filler_floor_separate_from_observed(tmp_path):
    """一个 bin 里 3 观测 + 2 填充:均值 5,对角方差 5(D=3 与 U=2 分开累计)。"""
    rows = [f"A,EVT,0.{i},50,60,-1,-1" for i in range(1, 4)] + [
        "A,FILL_GAP,0.5,50,60,-1,-1",
        "A,FILL_GAP,0.6,50,60,-1,-1",
    ]
    ev = R.load_events(_write(tmp_path, rows))
    lc = R.mean_and_diag(ev, np.array([0.0, 1.0]))
    assert lc.N[0] == 5
    assert lc.D[0] == 3
    assert lc.U[0] == 2
    assert lc.var_diag[0] == 5
    assert np.isclose(lc.err[0], np.sqrt(5.0))


def test_energy_band_filter(tmp_path):
    """chan_range 只数落在能段内的事件(观测与填充都过滤)。"""
    rows = [
        "A,EVT,0.1,30,60,-1,-1",     # 软,出段
        "A,EVT,0.2,100,60,-1,-1",    # 段内
        "A,FILL_GAP,0.3,120,60,-1,-1",  # 段内
        "A,FILL_GAP,0.4,300,60,-1,-1",  # 硬,出段
    ]
    ev = R.load_events(_write(tmp_path, rows))
    lc = R.mean_and_diag(ev, np.array([0.0, 1.0]), chan_range=(80, 200))
    assert lc.N[0] == 2, "只有 channel∈[80,200] 的 2 个事件"
    assert lc.D[0] == 1
    assert lc.U[0] == 1


def test_box_filter(tmp_path):
    """box 过滤只保留指定盒(全 HE 总光变则不过滤)。"""
    rows = [
        "A,EVT,0.1,50,60,-1,-1",
        "B,EVT,0.2,50,60,-1,-1",
        "C,FILL_GAP,0.3,50,60,-1,-1",
    ]
    ev = R.load_events(_write(tmp_path, rows))
    lc = R.mean_and_diag(ev, np.array([0.0, 1.0]), box="B")
    assert lc.N[0] == 1
    assert lc.D[0] == 1
    # 不过滤 → 全 HE 总光变
    lc_all = R.mean_and_diag(ev, np.array([0.0, 1.0]))
    assert lc_all.N[0] == 3
