"""饱和恢复光变的下游协方差装配(spec §7,三张表 + 1ms 网格解析协方差)。

恢复光变 N 是独立源 x=(观测源计数 C、标定 k、退化率 r、不可约丢失 U)的**线性像**。
在 1ms 网格上组装稀疏灵敏度矩阵 S,协方差为

    Cov(N) = S · diag(C) · Sᵀ  +  Σ_gap( k项 + r项 )  +  U对角(可开关)

方差记在 filler 填补的 **target 盒 / gap 时间位置**;filler↔参考的完全相关是 S 的跨盒
非对角,**自动出现、符号正确**(修复"总光变偏小"的 pull=1.32 病)。任意目标 binning 都
从同一 1ms 完整表示聚合:P · Cov · Pᵀ。**纯解析、确定、可复现(无 MC/RNG)。**

三张表(Rust 产出 / Python 消费):
  ① 事件流 events.csv:box,type,met,channel,...(type∈EVT|FILL_GAP)→ 均值 + 源计数 C
  ② gap 参数表 gapcov.csv(含 rho):每 gap 一行,变长字段分号分隔 → k项/r项/U/系统台账
  ③ gap 格结构表 gapbins.csv:每 (gap,1ms 格) 一行 → S 的 filler↔参考逐格系数(measured
     的 n_m;empty 的插值端点 left_bin/right_bin + τ)。**不在 Python 重推**(避免与 Rust 漂移)。
"""
from __future__ import annotations

import csv
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

SHAPE_BIN_WIDTH = 0.001  # 1ms,与 Rust SHAPE_BIN_WIDTH 一致


# ────────────────────────── ① 事件流 + 均值 ──────────────────────────

@dataclass
class LcCov:
    """一条恢复光变的均值与对角分量。"""

    edges: np.ndarray
    N: np.ndarray      # 均值:bin 内所有事件(EVT + FILL_GAP)
    D: np.ndarray      # 观测泊松对角(EVT 计数)
    U: np.ndarray      # 填充泊松地板(FILL_GAP 计数,spec §4 项 IV)

    @property
    def var_diag(self) -> np.ndarray:
        return self.D + self.U

    @property
    def err(self) -> np.ndarray:
        return np.sqrt(self.var_diag)


def load_events(path: str) -> dict:
    """读事件流 CSV → dict of numpy arrays(box/type/met/channel)。"""
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


def _select(events, box=None, chan_range=None) -> np.ndarray:
    """公共选择掩码:可选 box 与能段过滤。"""
    sel = np.ones(len(events["met"]), dtype=bool)
    if box is not None:
        sel &= events["box"] == box
    if chan_range is not None:
        lo, hi = chan_range
        sel &= (events["channel"] >= lo) & (events["channel"] <= hi)
    return sel


def mean_and_diag(events: dict, bin_edges, box=None, chan_range=None) -> LcCov:
    """按 bin_edges 分箱,返回均值 N 与对角分量 D(观测)/U(填充)。
    box=None → 全 HE 总光变(三盒相加);chan_range=(lo,hi) 闭区间能段选择。"""
    edges = np.asarray(bin_edges, dtype=float)
    met = events["met"]
    sel = _select(events, box, chan_range)
    is_evt = events["type"] == "EVT"
    is_fill = events["type"] == "FILL_GAP"
    n_all, _ = np.histogram(met[sel], bins=edges)
    d_obs, _ = np.histogram(met[sel & is_evt], bins=edges)
    u_fill, _ = np.histogram(met[sel & is_fill], bins=edges)
    return LcCov(edges=edges, N=n_all.astype(float),
                 D=d_obs.astype(float), U=u_fill.astype(float))


# ────────────────────────── ② gapcov + ③ gapbins I/O ──────────────────────────

def load_blocks(path: str) -> list:
    """读 gap 参数表 CSV(含 rho)→ list of dicts。变长字段(ref_boxes/k/c_ref_cal)
    分号分隔;空字段 → None/[]。"""
    def f_opt(s):
        return float(s) if s != "" else None

    def f_list(s):
        return [float(x) for x in s.split(";")] if s != "" else []

    def s_list(s):
        return s.split(";") if s != "" else []

    def as_bool(s):
        return s.strip().lower() == "true"

    blocks = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            blocks.append({
                "gap_id": int(row["gap_id"]),
                "target_box": row["target_box"],
                "type": row["type"],
                "t_start": float(row["t_start"]),
                "t_stop": float(row["t_stop"]),
                "ref_boxes": s_list(row["ref_boxes"]),
                "k": f_list(row["k"]),
                "c_ref_cal": f_list(row["c_ref_cal"]),
                "c_a_cal": f_opt(row["c_a_cal"]),
                "rho": f_opt(row.get("rho", "")) or 0.0,
                "r_pre": f_opt(row["r_pre"]),
                "r_post": f_opt(row["r_post"]),
                "n_pre": f_opt(row["n_pre"]),
                "n_post": f_opt(row["n_post"]),
                "maskable": as_bool(row["maskable"]),
                "sys_bias_flag": as_bool(row["sys_bias_flag"]),
                "sys_bias_scale": float(row["sys_bias_scale"]),
            })
    return blocks


def load_bins(path: str) -> dict:
    """读 gap 格结构表 CSV → dict: gap_id → {bin_index: row}。
    row 键:t_lo, n_m(measured=int / empty=None), kind, left_bin, right_bin, tau。"""
    def i_opt(s):
        return int(s) if s != "" else None

    def f_opt(s):
        return float(s) if s != "" else None

    out: dict = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            gid = int(row["gap_id"])
            out.setdefault(gid, {})[int(row["bin_index"])] = {
                "t_lo": float(row["t_lo"]),
                "n_m": i_opt(row["n_m"]),
                "kind": row["kind"],
                "left_bin": i_opt(row["left_bin"]),
                "right_bin": i_opt(row["right_bin"]),
                "tau": f_opt(row["tau"]),
            }
    return out


def degenerate_sigma2(block: dict) -> float:
    """退化 gap 的 gap 级方差 σ²_gap(spec §6):两端有率 →
    (T/2)²(r_pre²/(n_pre−1)+r_post²/(n_post−1));单侧 → T²r²/(n−1);
    (None,None) → NaN(maskable,不入协方差)。"""
    span = block["t_stop"] - block["t_start"]
    rp, rn = block.get("r_pre"), block.get("r_post")
    np_, nn = block.get("n_pre"), block.get("n_post")

    def var(r, n):
        return r * r / max((n or 1.0) - 1.0, 1.0)

    if rp is not None and rn is not None:
        return (span / 2.0) ** 2 * (var(rp, np_) + var(rn, nn))
    if rp is not None:
        return span * span * var(rp, np_)
    if rn is not None:
        return span * span * var(rn, nn)
    return float("nan")


# ────────────────────────── 1ms 网格协方差装配 ──────────────────────────

@dataclass
class FineCov:
    """1ms 网格上的完整稀疏协方差 + 索引信息(聚合前)。

    cov:scipy.sparse (n_cells×n_cells);cell index = box_pos·n_fine + fine_bin。
    """
    cov: "sp.spmatrix"
    boxes: list
    box_pos: dict
    n_fine: int
    fine_dt: float
    t0: float
    C: np.ndarray
    fine_edges: np.ndarray  # 并集网格边界(非均匀);cov_matrix 的聚合 P 必须用它

    def cell(self, box: str, fine_bin: int) -> int:
        return self.box_pos[box] * self.n_fine + fine_bin


def _var_rn(r, n):
    """退化率方差 r²/(n−1),n−1 下夹到 1。"""
    return r * r / max((n or 1.0) - 1.0, 1.0)


def assemble_fine(events, blocks, bins, bin_edges, box=None, chan_range=None,
                  include_u=False, fine_dt=SHAPE_BIN_WIDTH) -> FineCov:
    """在 1ms 网格上组装 Cov = S diag(C) Sᵀ + k项 + r项 + U(可开关),返回稀疏表示。

    - S 恒等:每观测格 (b,i) → S[(b,i),(b,i)]=1(gap-target 格除外,由 filler 覆盖)。
    - cross-ref filler(measured 格):S[(a,i),(b,i)] = ρ k_b / n_m。
    - cross-ref filler(empty 格):从端点 l,r 插值
        S[(a,i),(b,l)] += (1−τ)ρk_b/n_m_l;S[(a,i),(b,r)] += τ ρk_b/n_m_r。
    - k项:J^k Σ_k J^k = Σ_b g_b g_bᵀ/C_ref_cal_b + (Σ_b g_b)(·)ᵀ/C_a_cal(共用分子相关),
        g_{b,i} = Σ_{i'} S[(a,i),(b,i')] C_{b,i'}(参考 b 对 filler 格 i 的重建贡献)。
    - r项(退化,秩-2 斜坡):J_pre=(1−s)Δt、J_post=s·Δt,
        Cov += Var(r_pre) J_pre J_preᵀ + Var(r_post) J_post J_postᵀ。
    - U(可开关):include_u=True 时对 filler 格加重建计数(泊松地板);默认关(测量 scatter)。
    """
    edges = np.asarray(bin_edges, dtype=float)
    t0 = float(edges[0])
    t1 = float(edges[-1])
    tol = fine_dt * 1e-6

    # ── 并集细网格 = 分析 bin 边界 ∪ 每 cross-ref gap 的 Rust 格边界 ∪ 退化 gap 的 1ms 细分 ──
    # M2:每个细格既落唯一分析 bin(聚合 P 干净 0/1)、又落唯一 Rust 格(源计数按 Rust 格分箱)。
    # gap 跨分析边界时 Rust 格被切:measured 逐子格用各自源计数、empty 按宽度分摊。全程 Python。
    # 不 round(避免全精度 t_lo 查 round 过的边界产生碰撞);用容差去重。
    sel_blocks = [b for b in blocks if box is None or b["target_box"] == box]
    ev_edges = list(edges)  # 分析 bin 边界
    for b in sel_blocks:
        gb = bins.get(b["gap_id"])
        if gb:  # cross-ref:Rust 格边界(t_lo)+ 末边界 t_stop
            ev_edges.extend(row["t_lo"] for row in gb.values())
            ev_edges.append(b["t_stop"])
        else:   # 退化 gap:1ms 细分(r 项秩-2 斜坡需分辨率)
            nsd = max(int(round((b["t_stop"] - b["t_start"]) / fine_dt)), 1)
            ev_edges.extend(b["t_start"] + i * fine_dt for i in range(nsd + 1))
    ev_edges = np.array(sorted(x for x in ev_edges if t0 - tol <= x <= t1 + tol))
    if len(ev_edges) >= 2:
        fine_edges = ev_edges[np.concatenate(([True], np.diff(ev_edges) > tol))]
    if len(ev_edges) < 2 or len(fine_edges) < 2:
        fine_edges = np.array([t0, t1])
    n_fine = len(fine_edges) - 1
    fine_centers = 0.5 * (fine_edges[:-1] + fine_edges[1:])
    fine_widths = np.diff(fine_edges)

    # 盒集合:事件流出现的盒 ∪ 块表 target/参考盒(target 盒可能全饱和无事件)。
    box_set = set(events["box"].tolist())
    for blk in blocks:
        box_set.add(blk["target_box"])
        box_set.update(blk.get("ref_boxes") or [])
    boxes = sorted(box_set)
    box_pos = {b: i for i, b in enumerate(boxes)}
    nb = len(boxes)
    n_cells = nb * n_fine

    def cell(bp, fi):
        return bp * n_fine + fi

    def fine_of(t):
        # 并集网格上定位(+tol 让恰在边上的 t_lo 落到以它为左边界的格,不碰撞)
        i = int(np.searchsorted(fine_edges, t + tol, side="right")) - 1
        return min(max(i, 0), n_fine - 1)

    def fine_range(a, b):
        """Rust 格 [a,b) 覆盖的细格(中心落其中);gap 跨分析边界时 Rust 格被切成多个子格。"""
        return [fi for fi in range(n_fine) if a - tol <= fine_centers[fi] < b - tol]

    # ── 源计数 C_{b,i}(观测 EVT,能段过滤)──
    met = events["met"]
    is_evt = events["type"] == "EVT"
    is_fill = events["type"] == "FILL_GAP"
    if chan_range is not None:
        lo, hi = chan_range
        band = (events["channel"] >= lo) & (events["channel"] <= hi)
    else:
        band = np.ones(len(met), dtype=bool)

    C = np.zeros(n_cells, dtype=float)
    for b, bp in box_pos.items():
        m = (events["box"] == b) & is_evt & band
        cnt, _ = np.histogram(met[m], bins=fine_edges)
        C[bp * n_fine:(bp + 1) * n_fine] = cnt

    # ── 组装 S filler、k项、r项 ──
    s_rows, s_cols, s_vals = [], [], []
    k_rows, k_cols, k_vals = [], [], []
    r_rows, r_cols, r_vals = [], [], []
    gap_cells: set = set()

    def scatter(block, cells, rows, cols, vals):
        """把稠密小块 block[(len,len)] 按 cells 索引散射进 COO 三元组。"""
        for ii in range(len(cells)):
            ci = cells[ii]
            for jj in range(len(cells)):
                v = block[ii, jj]
                if v != 0.0:
                    rows.append(ci)
                    cols.append(cells[jj])
                    vals.append(v)

    for blk in blocks:
        if box is not None and blk["target_box"] != box:
            continue
        tb = blk["target_box"]
        if tb not in box_pos:
            continue
        tbp = box_pos[tb]
        gid = blk["gap_id"]
        is_cross = blk["type"] == "crossref" and blk["ref_boxes"]

        if is_cross:
            rho = blk.get("rho", 0.0) or 0.0
            ref_names = blk["ref_boxes"]
            ks = blk["k"]
            crefs = blk.get("c_ref_cal") or [None] * len(ref_names)
            ca = blk.get("c_a_cal")
            gb = bins.get(gid)
            if gb is None:  # 无格结构表(罕见):整段一个 measured 格,n_m=参考数
                gb = {0: {"t_lo": blk["t_start"], "n_m": len(ref_names),
                          "kind": "measured", "left_bin": None,
                          "right_bin": None, "tau": None}}
            # Rust 格按 bin_index 排序;各自区间 [t_lo, t_hi)(t_hi=下一格 t_lo 或 gap t_stop)
            items = sorted(gb.items())
            t_lo_of = {bi: row["t_lo"] for bi, row in items}
            n_m_of = {bi: row["n_m"] for bi, row in items}

            def t_hi_at(kk, _items=items, _tstop=blk["t_stop"], _tlo=t_lo_of):
                return _tlo[_items[kk + 1][0]] if kk + 1 < len(_items) else _tstop

            def rep_fine(bi, _items=items, _tlo=t_lo_of):
                # 端点 Rust 格的代表细格(格中心所在细格)
                if bi is None or bi not in _tlo:
                    return None
                kk = next(i for i, (b2, _) in enumerate(_items) if b2 == bi)
                return fine_of(0.5 * (_tlo[bi] + t_hi_at(kk)))

            rows_local = []  # (target_cell, g_vec over ref index)
            M = len(ref_names)
            for kk, (bi, row) in enumerate(items):
                lo, hi = t_lo_of[bi], t_hi_at(kk)
                fis = fine_range(lo, hi)  # 通常 1 个;gap 跨分析边界 → 多个子格
                if not fis:
                    continue
                w_rust = max(hi - lo, tol)
                if row["kind"] == "measured":
                    nm = row["n_m"]
                    if not nm or nm <= 0:
                        continue
                    for fi in fis:  # 各子格用各自源计数(切格核心)
                        tc = cell(tbp, fi)
                        gap_cells.add(tc)
                        gvec = np.zeros(M)
                        for m, (rb, kb) in enumerate(zip(ref_names, ks)):
                            if rb not in box_pos:
                                continue
                            rbp = box_pos[rb]
                            coeff = rho * kb / nm
                            col = cell(rbp, fi)
                            s_rows.append(tc); s_cols.append(col); s_vals.append(coeff)
                            gvec[m] = coeff * C[col]
                        rows_local.append((tc, gvec))
                else:  # empty:从端点 l,r 插值;切格时按子格宽度分摊
                    l, r = row["left_bin"], row["right_bin"]
                    tau = row["tau"] if row["tau"] is not None else 0.0
                    fl, fr = rep_fine(l), rep_fine(r)
                    nml, nmr = n_m_of.get(l), n_m_of.get(r)
                    if fl is None or fr is None or not nml or not nmr:
                        continue
                    for fi in fis:
                        tc = cell(tbp, fi)
                        gap_cells.add(tc)
                        frac = fine_widths[fi] / w_rust  # 不切=1;切格按宽度分摊,和=1
                        gvec = np.zeros(M)
                        for m, (rb, kb) in enumerate(zip(ref_names, ks)):
                            if rb not in box_pos:
                                continue
                            rbp = box_pos[rb]
                            cl = (1.0 - tau) * rho * kb / nml * frac
                            cr = tau * rho * kb / nmr * frac
                            col_l = cell(rbp, fl); col_r = cell(rbp, fr)
                            s_rows.append(tc); s_cols.append(col_l); s_vals.append(cl)
                            s_rows.append(tc); s_cols.append(col_r); s_vals.append(cr)
                            gvec[m] = cl * C[col_l] + cr * C[col_r]
                        rows_local.append((tc, gvec))

            # k项:J diag(1/c_ref) Jᵀ + (1/c_a)(J1)(J1)ᵀ
            if rows_local:
                cells_l = [tc for tc, _ in rows_local]
                Jm = np.array([g for _, g in rows_local])  # (nrows × M)
                inv_cref = np.array([
                    (1.0 / c) if (c is not None and c > 0) else 0.0
                    for c in crefs])
                Kblock = (Jm * inv_cref) @ Jm.T
                if ca is not None and ca > 0:
                    f_i = Jm.sum(axis=1)
                    Kblock = Kblock + np.outer(f_i, f_i) / ca
                scatter(Kblock, cells_l, k_rows, k_cols, k_vals)

        else:
            # 退化 gap:秩-2 斜坡 r项(无源 C 依赖 → S 无贡献)
            if blk.get("maskable"):
                # (None,None) 地板:不入协方差(只留 U 腿,若开)
                fis = [fi for fi in range(n_fine)
                       if blk["t_start"] <= fine_centers[fi] < blk["t_stop"]]
                gap_cells.update(cell(tbp, fi) for fi in fis)
                continue
            rp, rn = blk.get("r_pre"), blk.get("r_post")
            np_, nn = blk.get("n_pre"), blk.get("n_post")
            fis = [fi for fi in range(n_fine)
                   if blk["t_start"] <= fine_centers[fi] < blk["t_stop"]]
            if not fis:
                continue
            cells_d = [cell(tbp, fi) for fi in fis]
            gap_cells.update(cells_d)
            T = blk["t_stop"] - blk["t_start"]
            s = (fine_centers[fis] - blk["t_start"]) / T
            w = fine_widths[fis]  # 并集网格:逐格实际宽度(退化区≈fine_dt)
            if rp is not None and rn is not None:
                j_pre = (1.0 - s) * w
                j_post = s * w
                Rblock = (_var_rn(rp, np_) * np.outer(j_pre, j_pre)
                          + _var_rn(rn, nn) * np.outer(j_post, j_post))
            elif rp is not None:
                Rblock = _var_rn(rp, np_) * np.outer(w, w)
            elif rn is not None:
                Rblock = _var_rn(rn, nn) * np.outer(w, w)
            else:
                continue
            scatter(Rblock, cells_d, r_rows, r_cols, r_vals)

    # ── S 恒等(非 gap-target 格)──
    all_cells = np.arange(n_cells)
    if gap_cells:
        ident = np.setdiff1d(all_cells, np.array(sorted(gap_cells)),
                             assume_unique=True)
    else:
        ident = all_cells
    S = sp.coo_matrix(
        (np.concatenate([s_vals, np.ones(len(ident))]) if s_vals or len(ident)
         else np.zeros(0),
         (np.concatenate([s_rows, ident]).astype(int),
          np.concatenate([s_cols, ident]).astype(int))),
        shape=(n_cells, n_cells)).tocsr()

    Cdiag = sp.diags(C)
    cov = S @ Cdiag @ S.T

    if k_vals:
        cov = cov + sp.coo_matrix((k_vals, (k_rows, k_cols)),
                                  shape=(n_cells, n_cells)).tocsr()
    if r_vals:
        cov = cov + sp.coo_matrix((r_vals, (r_rows, r_cols)),
                                  shape=(n_cells, n_cells)).tocsr()

    if include_u:
        u = np.zeros(n_cells)
        for b, bp in box_pos.items():
            m = (events["box"] == b) & is_fill & band
            cnt, _ = np.histogram(met[m], bins=fine_edges)
            u[bp * n_fine:(bp + 1) * n_fine] = cnt
        cov = cov + sp.diags(u)

    return FineCov(cov=cov.tocsr(), boxes=boxes, box_pos=box_pos,
                   n_fine=n_fine, fine_dt=fine_dt, t0=t0, C=C, fine_edges=fine_edges)


def cov_matrix(events, blocks, bins, bin_edges, box=None, chan_range=None,
               include_u=False, fine_dt=SHAPE_BIN_WIDTH):
    """均值 N 与协方差 Cov(任意目标 binning)。

    在 1ms 网格上组装完整稀疏协方差(assemble_fine),再用聚合算子 P 聚合到目标
    binning:Cov_out = P · Cov_fine · Pᵀ(box=None → 三盒行求和的总光变;box=b →
    仅该盒行)。返回 (N, Cov_out),Cov_out 为稠密 (nbin×nbin)。

    include_u=False(默认):测量方差(画误差棒 / χ²)。
    include_u=True:叠加不可约丢失泊松地板(注入验证 / 审稿人闭环,非测量 scatter)。
    """
    edges = np.asarray(bin_edges, dtype=float)
    nbin = len(edges) - 1
    fc = assemble_fine(events, blocks, bins, edges, box=box,
                       chan_range=chan_range, include_u=include_u,
                       fine_dt=fine_dt)
    n_fine = fc.n_fine
    nb = len(fc.boxes)
    n_cells = nb * n_fine

    fine_edges = fc.fine_edges  # 并集网格(非均匀):必须用真实边界,不能假设 t0+arange·dt
    fine_centers = 0.5 * (fine_edges[:-1] + fine_edges[1:])
    tb_of_fine = np.clip(np.digitize(fine_centers, edges) - 1, 0, nbin - 1)

    if box is None:
        p_cols = np.arange(n_cells)
        p_rows = np.tile(tb_of_fine, nb)
    else:
        bp0 = fc.box_pos[box]
        p_cols = np.arange(bp0 * n_fine, (bp0 + 1) * n_fine)
        p_rows = tb_of_fine
    P = sp.coo_matrix((np.ones(len(p_cols)), (p_rows, p_cols)),
                      shape=(nbin, n_cells)).tocsr()

    cov_out = (P @ fc.cov @ P.T).toarray()
    N = mean_and_diag(events, edges, box=box, chan_range=chan_range).N
    return N, cov_out
