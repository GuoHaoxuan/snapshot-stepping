"""注入验证的 pull 组装（M5：把 pull 计算从 scratchpad 收进仓库、可复现）。

读 `blink sat inject` 产出的四张表（events/gapcov/gapbins/truth，每盒一套），逐
cross-ref gap 算 pull = (fill − truth) / σ，其中 σ 由 `recovery_cov.cov_matrix`
**从三表解析算**（include_u=True，"恢复 vs 真值"腿），而不是朴素 √fill。汇总
pull.mean / pull.std / bias。

复现（250919A，需 HXMT_1B_DIR 指向本地 data/1B、blink release 已 build）：
  scripts/injection_validation.sh <out_dir>          # measured 腿 → pull.std≈1.05
  scripts/injection_validation.sh <out_dir> --cosat  # 共饱和 empty 腿 → pull.std≈1.00
或手动：blink sat inject … --events-out/--gapcov-out/--gapbins-out/--truth-out，再
  python analysis/injection_pull.py <out_dir>
"""
from __future__ import annotations

import csv

import numpy as np

import recovery_cov as rc


def load_truth(path: str) -> dict:
    """读 truth 表（gap_id,t_start,t_stop,n_truth,n_fill,type）→ {gap_id: n_truth}。"""
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            out[int(r["gap_id"])] = int(r["n_truth"])
    return out


def gap_pulls(events, blocks, bins, truth, box, include_u=True):
    """逐 gap 的 (gap_id, truth, fill, sigma, pull)。σ 覆盖整个 gap 单 bin，由
    cov_matrix 从三表解析算：cross-ref 用 S·diag(C)·Sᵀ+k 项，degenerate（全宽共饱和、
    无有效参考）用 r 项秩-2 斜坡。无真值、或零方差（maskable 退化：只有 MCU 地板率、
    无从估方差）的 gap 跳过。"""
    rows = []
    for blk in blocks:
        gid = blk["gap_id"]
        if gid not in truth:
            continue
        edges = np.array([blk["t_start"], blk["t_stop"]])
        N, cov = rc.cov_matrix(events, blocks, bins, edges, box=box, include_u=include_u)
        var = float(cov[0, 0])
        if var <= 0:
            continue
        fill = float(N[0])
        sigma = np.sqrt(var)
        rows.append((gid, truth[gid], fill, sigma, (fill - truth[gid]) / sigma))
    return rows


def summarize(rows) -> dict:
    if not rows:
        return {"n": 0}
    pull = np.array([r[4] for r in rows])
    bias = np.array([(r[2] - r[1]) / r[1] for r in rows if r[1] > 0])
    return {
        "n": len(pull),
        "pull_mean": float(pull.mean()),
        "pull_std": float(pull.std(ddof=1)) if len(pull) > 1 else float("nan"),
        "bias_pct": float(100 * bias.mean()) if len(bias) else float("nan"),
    }


def run_dir(ws, boxes=("a", "b", "c"), include_u=True):
    """读 {ws}/{events,gapcov,gapbins,truth}_{box}.csv，汇总三盒 pull。"""
    rows = []
    for box in boxes:
        ev = rc.load_events(f"{ws}/events_{box}.csv")
        bl = rc.load_blocks(f"{ws}/gapcov_{box}.csv")
        bn = rc.load_bins(f"{ws}/gapbins_{box}.csv")
        tr = load_truth(f"{ws}/truth_{box}.csv")
        rows += gap_pulls(ev, bl, bn, tr, box.upper(), include_u=include_u)
    return rows, summarize(rows)


def main():
    import argparse

    ap = argparse.ArgumentParser(description="injection pull from three-table output")
    ap.add_argument("dir", help="dir with {events,gapcov,gapbins,truth}_{a,b,c}.csv")
    ap.add_argument("--boxes", default="a,b,c")
    ap.add_argument("--no-u", action="store_true", help="drop U leg (measurement variance)")
    args = ap.parse_args()
    _rows, s = run_dir(args.dir, tuple(args.boxes.split(",")), include_u=not args.no_u)
    if s["n"] == 0:
        print("no cross-ref gaps found")
        return
    print(
        f"n={s['n']}  pull.mean={s['pull_mean']:+.3f}  pull.std={s['pull_std']:.3f}  "
        f"bias={s['bias_pct']:+.2f}%"
    )


if __name__ == "__main__":
    main()
