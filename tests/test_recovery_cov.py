"""Tests for analysis/recovery_cov.py.

下游协方差装配的权威数学(spec §7,三张表 + 1ms 网格解析协方差
S·diag(C)·Sᵀ):恢复光变 N 是独立源 x=(观测源计数 C、标定 k、退化率 r、
不可约丢失 U)的线性像。方差记在 filler 填补的 target 位置,filler↔参考的完全
相关自动出现,总光变逐 bin 方差复现逐粒子 Σw²C(修 pull=1.32)。
"""
from __future__ import annotations

import numpy as np

import recovery_cov as R

HEADER = "box,type,met,channel,pulse_width,pkt_idx,evt_idx"


def _write(tmp_path, rows):
    p = tmp_path / "ev.csv"
    p.write_text(HEADER + "\n" + "\n".join(rows) + "\n")
    return str(p)


def _events(recs):
    """recs: list of (box, type, met, channel) → events dict。"""
    return {
        "box": np.array([r[0] for r in recs]),
        "type": np.array([r[1] for r in recs]),
        "met": np.array([r[2] for r in recs], dtype=float),
        "channel": np.array([r[3] for r in recs], dtype=int),
    }


# ────────────────────────── 均值 / 事件流 ──────────────────────────

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
    assert list(lc.N) == [2, 1]
    assert list(lc.D) == [2, 0]
    assert list(lc.U) == [0, 1]


def test_energy_band_filter_mean(tmp_path):
    rows = [
        "A,EVT,0.1,30,60,-1,-1",
        "A,EVT,0.2,100,60,-1,-1",
        "A,FILL_GAP,0.3,120,60,-1,-1",
        "A,FILL_GAP,0.4,300,60,-1,-1",
    ]
    ev = R.load_events(_write(tmp_path, rows))
    lc = R.mean_and_diag(ev, np.array([0.0, 1.0]), chan_range=(80, 200))
    assert lc.N[0] == 2
    assert lc.D[0] == 1
    assert lc.U[0] == 1


# ────────────────────────── 三张表 I/O ──────────────────────────

BHEADER = (
    "gap_id,target_box,type,t_start,t_stop,ref_boxes,k,c_ref_cal,c_a_cal,rho,"
    "r_pre,r_post,n_pre,n_post,maskable,sys_bias_flag,sys_bias_scale"
)


def _write_blocks(tmp_path, rows):
    p = tmp_path / "gapcov.csv"
    p.write_text(BHEADER + "\n" + "\n".join(rows) + "\n")
    return str(p)


def test_load_blocks_parses_rho_and_both_types(tmp_path):
    rows = [
        "0,A,crossref,1.000000,1.100000,B;C,2.0000;1.5000,100;80,200,0.9500,,,,,false,false,0.0000",
        "1,A,degenerate,2.000000,2.100000,,,,,0.0000,400.0000,600.0000,40,60,false,true,0.2000",
    ]
    blocks = R.load_blocks(_write_blocks(tmp_path, rows))
    assert len(blocks) == 2
    cr, dg = blocks[0], blocks[1]
    assert cr["type"] == "crossref"
    assert cr["ref_boxes"] == ["B", "C"]
    assert cr["k"] == [2.0, 1.5]
    assert cr["c_ref_cal"] == [100.0, 80.0]
    assert cr["c_a_cal"] == 200.0
    assert np.isclose(cr["rho"], 0.95)
    assert dg["type"] == "degenerate"
    assert dg["ref_boxes"] == []
    assert dg["r_pre"] == 400.0 and dg["n_pre"] == 40.0
    assert np.isclose(dg["rho"], 0.0)


BINHEADER = "gap_id,bin_index,t_lo,n_m,kind,left_bin,right_bin,tau"


def _write_bins(tmp_path, rows):
    p = tmp_path / "gapbins.csv"
    p.write_text(BINHEADER + "\n" + "\n".join(rows) + "\n")
    return str(p)


def test_load_bins_measured_and_empty(tmp_path):
    rows = [
        "0,0,0.000000,2,measured,,,",
        "0,1,0.001000,,empty,0,2,0.5000",
        "0,2,0.002000,1,measured,,,",
    ]
    bins = R.load_bins(_write_bins(tmp_path, rows))
    g0 = bins[0]
    assert g0[0]["kind"] == "measured"
    assert g0[0]["n_m"] == 2
    assert g0[0]["left_bin"] is None
    e = g0[1]
    assert e["kind"] == "empty"
    assert e["n_m"] is None
    assert e["left_bin"] == 0 and e["right_bin"] == 2
    assert np.isclose(e["tau"], 0.5)
    assert np.isclose(g0[0]["t_lo"], 0.0)
    assert np.isclose(g0[2]["t_lo"], 0.002)


# ────────────────── ① 总光变对角复现 Σw²C ──────────────────

def _crossref_block(gid, target, t0, t1, refs, ks, rho,
                    c_ref_cal=1e12, c_a_cal=1e12):
    return {
        "gap_id": gid, "target_box": target, "type": "crossref",
        "t_start": t0, "t_stop": t1, "ref_boxes": list(refs), "k": list(ks),
        "c_ref_cal": [c_ref_cal] * len(refs), "c_a_cal": c_a_cal, "rho": rho,
        "r_pre": None, "r_post": None, "n_pre": None, "n_post": None,
        "maskable": False, "sys_bias_flag": False, "sys_bias_scale": 0.0,
    }


def test_total_diagonal_reproduces_sum_w2_c():
    """单盒 gap 的总光变对角 = Σ_b w_b² C_b(w_b = 1 + ρk_b/n_m)。

    A 在 [0,1ms) 饱和,参考 B(10 EVT)、C(20 EVT);ρ=1,k=1,n_m=2。
    → w_B = w_C = 1.5;Σw²C = 1.5²·10 + 1.5²·20 = 67.5。
    标定计数取巨大 → k 项 →0;U 关。只留源腿 S diag(C) Sᵀ。
    """
    recs = [("B", "EVT", 0.0005, 50)] * 10 + [("C", "EVT", 0.0005, 50)] * 20
    recs += [("A", "FILL_GAP", 0.0005, 50)] * 5
    ev = _events(recs)
    blk = _crossref_block(0, "A", 0.0, 0.001, ["B", "C"], [1.0, 1.0], rho=1.0)
    bins = {0: {0: {"t_lo": 0.0, "n_m": 2, "kind": "measured",
                    "left_bin": None, "right_bin": None, "tau": None}}}
    edges = np.array([0.0, 0.001])
    _N, Cov = R.cov_matrix(ev, [blk], bins, edges, box=None, include_u=False)
    w = 1.0 + 1.0 * 1.0 / 2.0
    expected = w * w * 10 + w * w * 20
    assert np.isclose(Cov[0, 0], expected, rtol=1e-6), (Cov[0, 0], expected)


def test_total_diag_matches_particlewise_sum_w2c_multibin():
    """多 1ms 格 + 多参考的一般自检:P(S diag C Sᵀ)Pᵀ 对角 == Σ_i Σ_b w_{b,i}² C_{b,i}。"""
    recs = []
    # B,C 每个 1ms 格若干 EVT;A 三格全饱和
    for i, (nb, nc) in enumerate([(5, 8), (7, 3), (2, 9)]):
        t = 0.0005 + i * 0.001
        recs += [("B", "EVT", t, 50)] * nb + [("C", "EVT", t, 50)] * nc
        recs += [("A", "FILL_GAP", t, 50)]
    ev = _events(recs)
    blk = _crossref_block(0, "A", 0.0, 0.003, ["B", "C"], [2.0, 1.0], rho=0.8)
    bins = {0: {i: {"t_lo": i * 0.001, "n_m": 2, "kind": "measured",
                    "left_bin": None, "right_bin": None, "tau": None}
                for i in range(3)}}
    edges = np.array([0.0, 0.001, 0.002, 0.003])
    _N, Cov = R.cov_matrix(ev, [blk], bins, edges, box=None, include_u=False)
    # 手算逐格 Σ_b w_b² C_b
    wB = 1.0 + 0.8 * 2.0 / 2.0
    wC = 1.0 + 0.8 * 1.0 / 2.0
    counts = [(5, 8), (7, 3), (2, 9)]
    for i, (nb, nc) in enumerate(counts):
        exp = wB * wB * nb + wC * wC * nc
        assert np.isclose(Cov[i, i], exp, rtol=1e-6), (i, Cov[i, i], exp)


# ────────────────── ② filler↔参考跨盒非对角 ──────────────────

def test_cross_box_offdiagonal_sign_and_size():
    """Cov_fine[(A_filler,i),(B_identity,i)] = ρk_B/n_m · C_B > 0(自动出现,符号正确)。"""
    recs = [("B", "EVT", 0.0005, 50)] * 10 + [("C", "EVT", 0.0005, 50)] * 20
    ev = _events(recs)
    blk = _crossref_block(0, "A", 0.0, 0.001, ["B", "C"], [1.0, 1.0], rho=1.0)
    bins = {0: {0: {"t_lo": 0.0, "n_m": 2, "kind": "measured",
                    "left_bin": None, "right_bin": None, "tau": None}}}
    edges = np.array([0.0, 0.001])
    fc = R.assemble_fine(ev, [blk], bins, edges, box=None, include_u=False)
    a = fc.cell("A", 0)
    b = fc.cell("B", 0)
    c = fc.cell("C", 0)
    cov = fc.cov.toarray()
    assert cov[a, b] > 0
    assert np.isclose(cov[a, b], 1.0 * 1.0 / 2.0 * 10)   # 0.5*10 = 5
    assert np.isclose(cov[a, c], 1.0 * 1.0 / 2.0 * 20)   # 0.5*20 = 10
    assert np.isclose(cov[a, b], cov[b, a])              # 对称


# ────────────────── ③ measured vs empty 格系数 ──────────────────

def test_measured_vs_empty_coefficients():
    """measured 格 S=ρk/n_m·恒等;empty 格从左右端点插值 (1−τ)/τ。"""
    # gap A [0,3ms):bin0 measured(ref B), bin1 empty(l=0,r=2,τ=0.5), bin2 measured
    recs = [("B", "EVT", 0.0005, 50)] * 4 + [("B", "EVT", 0.0025, 50)] * 6
    ev = _events(recs)
    blk = _crossref_block(0, "A", 0.0, 0.003, ["B"], [1.0], rho=1.0)
    bins = {0: {
        0: {"t_lo": 0.000, "n_m": 1, "kind": "measured",
            "left_bin": None, "right_bin": None, "tau": None},
        1: {"t_lo": 0.001, "n_m": None, "kind": "empty",
            "left_bin": 0, "right_bin": 2, "tau": 0.5},
        2: {"t_lo": 0.002, "n_m": 1, "kind": "measured",
            "left_bin": None, "right_bin": None, "tau": None},
    }}
    edges = np.array([0.0, 0.001, 0.002, 0.003])
    fc = R.assemble_fine(ev, [blk], bins, edges, box=None, include_u=False)
    cov = fc.cov.toarray()
    a0, a1, a2 = fc.cell("A", 0), fc.cell("A", 1), fc.cell("A", 2)
    b0, b2 = fc.cell("B", 0), fc.cell("B", 2)
    C_b0, C_b2 = 4.0, 6.0
    # measured bin0: coeff = ρk/n_m = 1 → cov[a0,b0] = 1*C_b0
    assert np.isclose(cov[a0, b0], 1.0 * C_b0)
    assert np.isclose(cov[a2, b2], 1.0 * C_b2)
    # empty bin1: 从 b0 权重 (1-τ)=0.5、从 b2 权重 τ=0.5
    assert np.isclose(cov[a1, b0], 0.5 * C_b0)
    assert np.isclose(cov[a1, b2], 0.5 * C_b2)


# ────────────────── ④ U 开关腿 ──────────────────

def test_u_leg_switchable():
    """include_u=False 默认无泊松地板;True 时 filler 格对角加重建计数(泊松地板)。"""
    recs = [("B", "EVT", 0.0005, 50)] * 10
    recs += [("A", "FILL_GAP", 0.0005, 50)] * 3   # 3 filler → U 地板 3
    ev = _events(recs)
    blk = _crossref_block(0, "A", 0.0, 0.001, ["B"], [1.0], rho=1.0)
    bins = {0: {0: {"t_lo": 0.0, "n_m": 1, "kind": "measured",
                    "left_bin": None, "right_bin": None, "tau": None}}}
    edges = np.array([0.0, 0.001])
    _N, cov_off = R.cov_matrix(ev, [blk], bins, edges, box="A", include_u=False)
    _N, cov_on = R.cov_matrix(ev, [blk], bins, edges, box="A", include_u=True)
    assert np.isclose(cov_on[0, 0] - cov_off[0, 0], 3.0), (cov_off, cov_on)


# ────────────────── ⑤ 能段选择 ──────────────────

def test_energy_band_selects_source_counts():
    """chan_range 只数落在能段内的参考计数 → 影响 C,进而影响协方差。"""
    recs = [("B", "EVT", 0.0005, 100)] * 6 + [("B", "EVT", 0.0005, 500)] * 4
    recs += [("A", "FILL_GAP", 0.0005, 100)]
    ev = _events(recs)
    blk = _crossref_block(0, "A", 0.0, 0.001, ["B"], [1.0], rho=1.0)
    bins = {0: {0: {"t_lo": 0.0, "n_m": 1, "kind": "measured",
                    "left_bin": None, "right_bin": None, "tau": None}}}
    edges = np.array([0.0, 0.001])
    # 能段 [80,200] → 只 6 个 B 计数;filler 系数 ρk/n_m = 1 → A 行方差 = 1²·6
    _N, cov = R.cov_matrix(ev, [blk], bins, edges, box="A", chan_range=(80, 200))
    assert np.isclose(cov[0, 0], 6.0)
    _N, cov_all = R.cov_matrix(ev, [blk], bins, edges, box="A")
    assert np.isclose(cov_all[0, 0], 10.0)


# ────────────────── ⑥ 任意 binning(1ms 与 0.1s 一致) ──────────────────

def test_binning_invariance_1ms_vs_coarse():
    """0.1s 目标 == 1ms 目标结果按块聚合(同一 1ms 完整表示 P·Cov·Pᵀ)。"""
    rng = np.random.default_rng(0)
    recs = []
    n_fine = 20
    for i in range(n_fine):
        t = 0.0005 + i * 0.001
        for _ in range(int(rng.integers(3, 9))):
            recs.append(("B", "EVT", t, 50))
        for _ in range(int(rng.integers(3, 9))):
            recs.append(("C", "EVT", t, 50))
        recs.append(("A", "FILL_GAP", t, 50))
    ev = _events(recs)
    blk = _crossref_block(0, "A", 0.0, 0.020, ["B", "C"], [1.5, 1.0], rho=0.9)
    bins = {0: {i: {"t_lo": i * 0.001, "n_m": 2, "kind": "measured",
                    "left_bin": None, "right_bin": None, "tau": None}
                for i in range(n_fine)}}
    fine_edges = 0.0 + np.arange(n_fine + 1) * 0.001
    coarse_edges = np.array([0.0, 0.010, 0.020])
    _N, cov_fine = R.cov_matrix(ev, [blk], bins, fine_edges, box=None)
    _N, cov_coarse = R.cov_matrix(ev, [blk], bins, coarse_edges, box=None)
    # 把 1ms 结果按 10 格块聚合
    P = np.zeros((2, n_fine))
    P[0, :10] = 1.0
    P[1, 10:] = 1.0
    agg = P @ cov_fine @ P.T
    assert np.allclose(agg, cov_coarse, rtol=1e-9, atol=1e-9)


# ────────────────── 退化段 r 项(秩-2 斜坡) ──────────────────

def test_degenerate_r_term_total_variance():
    """退化 gap 的总 gap 方差(整段求和)= σ²_gap = (T/2)²(Var_pre + Var_post)。

    斜坡 Jacobian J_pre=(1−s)Δt、J_post=s·Δt 精确给 Σ J_pre = Σ J_post = T/2。
    """
    n_fine = 10
    T = n_fine * 0.001
    recs = [("A", "FILL_GAP", 0.0005 + i * 0.001, 50) for i in range(n_fine)]
    ev = _events(recs)
    blk = {
        "gap_id": 0, "target_box": "A", "type": "degenerate",
        "t_start": 0.0, "t_stop": T, "ref_boxes": [], "k": [],
        "c_ref_cal": [], "c_a_cal": None, "rho": 0.0,
        "r_pre": 400.0, "r_post": 600.0, "n_pre": 40.0, "n_post": 60.0,
        "maskable": False, "sys_bias_flag": True, "sys_bias_scale": 0.2,
    }
    edges = np.array([0.0, T])   # 单一目标 bin = 整段
    _N, cov = R.cov_matrix(ev, [blk], {}, edges, box="A", include_u=False)
    var_pre = 400.0 ** 2 / 39.0
    var_post = 600.0 ** 2 / 59.0
    sigma2 = (T / 2.0) ** 2 * (var_pre + var_post)
    assert np.isclose(cov[0, 0], sigma2, rtol=1e-9), (cov[0, 0], sigma2)


def test_maskable_degenerate_skipped():
    n_fine = 4
    T = n_fine * 0.001
    recs = [("A", "FILL_GAP", 0.0005 + i * 0.001, 50) for i in range(n_fine)]
    ev = _events(recs)
    blk = {
        "gap_id": 0, "target_box": "A", "type": "degenerate",
        "t_start": 0.0, "t_stop": T, "ref_boxes": [], "k": [], "c_ref_cal": [],
        "c_a_cal": None, "rho": 0.0, "r_pre": None, "r_post": None,
        "n_pre": None, "n_post": None, "maskable": True,
        "sys_bias_flag": True, "sys_bias_scale": 0.0,
    }
    edges = np.array([0.0, T])
    _N, cov = R.cov_matrix(ev, [blk], {}, edges, box="A", include_u=False)
    assert np.isclose(cov[0, 0], 0.0)   # 无源、无 r、U 关 → 全 0


# ────────────────── k 项(标定满协方差,含共用分子) ──────────────────

def test_k_term_shared_numerator_cross_ref():
    """k 项 = Σ_b g_b g_b/C_ref_cal_b + (Σg)²/C_a_cal(共用分子跨-ref 相关)。

    单 1ms 格,measured,ρ=1,n_m=2,k=[1,1];B=10、C=20 EVT。
    g_B = ρk_B/n_m·C_B = 5,g_C = 10,f = 15。
    k 对角(A 行)= 25/C_refB + 100/C_refC + 15²/C_a。
    """
    recs = [("B", "EVT", 0.0005, 50)] * 10 + [("C", "EVT", 0.0005, 50)] * 20
    ev = _events(recs)
    blk = _crossref_block(0, "A", 0.0, 0.001, ["B", "C"], [1.0, 1.0], rho=1.0,
                          c_ref_cal=None)  # 覆盖下面手设
    blk["c_ref_cal"] = [100.0, 80.0]
    blk["c_a_cal"] = 200.0
    bins = {0: {0: {"t_lo": 0.0, "n_m": 2, "kind": "measured",
                    "left_bin": None, "right_bin": None, "tau": None}}}
    edges = np.array([0.0, 0.001])
    fc = R.assemble_fine(ev, [blk], bins, edges, box="A", include_u=False)
    # 源腿 A 行对角 = Σ (ρk/n_m)² C = 0.5²·10 + 0.5²·20 = 7.5
    # k 腿 A 行对角 = 5²/100 + 10²/80 + 15²/200
    g_B, g_C, f = 5.0, 10.0, 15.0
    k_diag = g_B ** 2 / 100.0 + g_C ** 2 / 80.0 + f ** 2 / 200.0
    src_diag = 0.5 ** 2 * 10 + 0.5 ** 2 * 20
    a = fc.cell("A", 0)
    assert np.isclose(fc.cov.toarray()[a, a], src_diag + k_diag), fc.cov.toarray()[a, a]
