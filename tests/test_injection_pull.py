"""Tests for analysis/injection_pull.py (M5: repo-reproducible injection pull)."""
from __future__ import annotations

import numpy as np

import injection_pull as ip


def _events(recs):
    return {
        "box": np.array([r[0] for r in recs]),
        "type": np.array([r[1] for r in recs]),
        "met": np.array([r[2] for r in recs], dtype=float),
        "channel": np.array([r[3] for r in recs], dtype=int),
    }


def test_gap_pull_is_fill_minus_truth_over_sigma():
    """A gap [0,1ms): refs B(10)/C(20), k=1, rho=1, n_m=2, 15 fillers; truth=12.
    fill=15; σ²=U(15)+ S·diag(C)·Sᵀ((0.5)²·10+(0.5)²·20=7.5)=22.5; pull=(15−12)/√22.5.
    标定计数取巨大 → k 项≈0。
    """
    recs = (
        [("B", "EVT", 0.0005, 50)] * 10
        + [("C", "EVT", 0.0005, 50)] * 20
        + [("A", "FILL_GAP", 0.0005, 50)] * 15
    )
    ev = _events(recs)
    blk = {
        "gap_id": 0, "target_box": "A", "type": "crossref",
        "t_start": 0.0, "t_stop": 0.001, "ref_boxes": ["B", "C"], "k": [1.0, 1.0],
        "c_ref_cal": [1e12, 1e12], "c_a_cal": 1e12, "rho": 1.0,
        "r_pre": None, "r_post": None, "n_pre": None, "n_post": None,
        "maskable": False, "sys_bias_flag": False, "sys_bias_scale": 0.0,
    }
    bins = {0: {0: {"t_lo": 0.0, "n_m": 2, "kind": "measured",
                    "left_bin": None, "right_bin": None, "tau": None}}}
    truth = {0: 12}
    rows = ip.gap_pulls(ev, [blk], bins, truth, "A", include_u=True)
    assert len(rows) == 1
    gid, t, fill, sigma, pull = rows[0]
    assert gid == 0 and t == 12
    assert np.isclose(fill, 15.0)
    assert np.isclose(sigma, np.sqrt(22.5), rtol=1e-6)
    assert np.isclose(pull, 3.0 / np.sqrt(22.5), rtol=1e-6)


def test_degenerate_gap_skipped():
    """退化 gap 不产 gapbins/无 cross-ref → gap_pulls 跳过（只做 cross-ref 腿）。"""
    ev = _events([("A", "FILL_GAP", 0.5, 50)])
    blk = {
        "gap_id": 0, "target_box": "A", "type": "degenerate",
        "t_start": 0.0, "t_stop": 1.0, "ref_boxes": [], "k": [],
        "c_ref_cal": [], "c_a_cal": None, "rho": 0.0,
        "r_pre": 400.0, "r_post": 600.0, "n_pre": 40.0, "n_post": 60.0,
        "maskable": False, "sys_bias_flag": True, "sys_bias_scale": 0.2,
    }
    rows = ip.gap_pulls(ev, [blk], {0: {}}, {0: 5}, "A")
    assert rows == []


def test_summarize_mean_std_bias():
    # (gid, truth, fill, sigma, pull)
    rows = [(0, 10, 12, 2.0, 1.0), (1, 20, 18, 2.0, -1.0), (2, 10, 10, 1.0, 0.0)]
    s = ip.summarize(rows)
    assert s["n"] == 3
    assert np.isclose(s["pull_mean"], 0.0)
    # bias frac: (2/10, -2/20, 0) = (0.2,-0.1,0) → mean 0.0333 → 3.33%
    assert np.isclose(s["bias_pct"], 100 * (0.2 - 0.1 + 0.0) / 3, rtol=1e-6)
