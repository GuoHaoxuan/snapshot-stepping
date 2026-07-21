#!/usr/bin/env bash
# Reproduce the injection-pull validation of the analytic error model (spec §11/§12).
#   measured leg (all-measured cells):        pull.std ~ 1.05, bias ~ +0.15%
#   co-saturation empty-cell leg (--cosat):    pull.std ~ 1.00, bias ~ +0.10%
#   degenerate leg (--degenerate, M4):         pre/post ramp extrapolation, r-term sigma
#
# Mechanism: inject fake 30ms gaps into 250919A's flat pre-burst baseline (target
# box masked, its in-gap events kept as truth); reconstruct; compare fill vs truth
# with sigma from the three-table analytic covariance (include_u).
#   --cosat      marks reference boxes unreliable in each gap's MIDDLE 10ms
#                (partial co-saturation) so genuine cross-ref empty cells appear.
#   --degenerate marks reference boxes unreliable across the WHOLE gap (cosat width
#                > gap width) so no valid reference survives → has_cross_ref=false →
#                the degenerate leg reconstructs from the target's own pre/post
#                neighbor-packet rates (needs inject's real prev/next_pkt_idx, M4),
#                and sigma comes from the r-term rank-2 ramp.
#
# Requires: release build (cargo build -p blink --release), HXMT_1B_DIR pointing
# at a local data/1B tree, and a Python with numpy/scipy+recovery_cov importable
# (set PYTHON=.venv/bin/python if your venv is not the default python3).
#
# Usage: HXMT_1B_DIR="$PWD/data/1B" ./scripts/injection_validation.sh <out_dir> [--cosat|--degenerate]
set -euo pipefail

OUT=${1:?usage: injection_validation.sh <out_dir> [--cosat|--degenerate]}
COSAT=""
case "${2:-}" in
  --cosat)      COSAT="--cosat-width 0.01" ;;   # gap 中段 10ms 共饱和 → cross-ref empty 格
  --degenerate) COSAT="--cosat-width 0.031" ;;  # 全宽共饱和(>width 0.03) → 无参考 → 退化腿
  "")           ;;
  *)            echo "unknown mode: ${2}" >&2; exit 2 ;;
esac
PYTHON=${PYTHON:-python3}
BIN=${BLINK_BIN:-./target/release/blink}
TRIG=432865758                            # GRB 250919A T0 MET (2025-09-19T00:29:15 UTC)
AT=$(seq -4.8 0.1 -0.3 | paste -sd, -)    # 46 flat-baseline gaps per box

mkdir -p "$OUT"
for BOX in a b c; do
  "$BIN" sat inject "$TRIG" --before 5 --after 20 --target "$BOX" --at="$AT" \
    --width 0.03 $COSAT \
    --events-out "$OUT/events_$BOX.csv" --gapcov-out "$OUT/gapcov_$BOX.csv" \
    --gapbins-out "$OUT/gapbins_$BOX.csv" --truth-out "$OUT/truth_$BOX.csv" >/dev/null
done
"$PYTHON" analysis/injection_pull.py "$OUT"
