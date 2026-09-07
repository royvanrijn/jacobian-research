#!/usr/bin/env sage -python
"""Certify that the q=6 RR kernel has an actual resolved local-module cover.

``assemble_h92_q6_global_rr.sage`` constructs the ten-dimensional bounded
coefficient ambient and its sole finite (smooth-collision) 8-by-10 block.
This companion supplies the actual H92 E7 provenance: an all-edge resolved
module cover and the corrected marked +/-P1 frame. It then checks that the
displayed ambient already satisfies the E8 and E7 coefficient bounds, so no
unrecorded rows are being silently omitted.

Consequently the rank-eight collision matrix is the complete q=6
vertical-condition matrix for the declared resolved cover, and its
two-dimensional kernel certifies h0(D)=2.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
RR = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-global-rr.json"
E8 = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-e8-p1-branch-module.json"
SMOOTH = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-smooth-po-module.json"
E7 = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-all-edge-module.json"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-resolved-rr-cover.json"
RR_SHA256 = "88466369aea8f838dde85051b174f6f2ea4c6df3edb29828b70ee8e79b9975cb"
E8_SHA256 = "97dbadd5a9e00b95106266aab04e3d37911b9800f8f219b9028a34b32530fb85"
SMOOTH_SHA256 = "a950defe18be876d96b215d287d02aebfe4c0a375b788b1c3bfbf0ff864b839d"
E7_SHA256 = "0127d9f0591e8cd6a4c9370cd9362beacfad5711ee1d336789bcc34d532f34d9"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Replay the finite smooth-chart block through the generic compiler interface.
# The source RR artifact retains the chart derivation; this script supplies its
# exact coefficient matrix as a finite quotient evaluator so that the q=6
# regression exercises the same stacking/kernel code used by later hops.
exec(compile(CORE.read_text(), str(CORE), "exec"))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rr", type=Path, default=RR)
parser.add_argument("--e8", type=Path, default=E8)
parser.add_argument("--smooth", type=Path, default=SMOOTH)
parser.add_argument("--e7", type=Path, default=E7)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

for path, expected in ((args.rr, RR_SHA256), (args.e8, E8_SHA256),
                       (args.smooth, SMOOTH_SHA256), (args.e7, E7_SHA256)):
    defaults = {RR: RR_SHA256, E8: E8_SHA256, SMOOTH: SMOOTH_SHA256, E7: E7_SHA256}
    if path in defaults:
        assert digest(path) == expected
rr = json.loads(args.rr.read_text())
e8 = json.loads(args.e8.read_text())
smooth = json.loads(args.smooth.read_text())
e7 = json.loads(args.e7.read_text())
assert rr["status"] == "PASS_EXACT_GLOBAL_RR_KERNEL"
assert e8["status"] == "PASS_EXACT_E8_MARKED_CHORD_MODULE"
assert smooth["status"] == "PASS_EXACT_SMOOTH_PO_CHORD_MODULE"
assert e7["status"] == "PASS_EXACT_Q6_ACTUAL_E7_ALL_EDGE_MODULE"

ambient = rr["ambient"]
assert ambient["coefficient_form"] == "a=A/h^2, b=B/h, section=a+b*m"
assert ambient["A_powers"] == list(range(1, 9))
assert ambient["B_powers"] == [3, 4]
assert ambient["dimension"] == 10

# At E8 (u=0), the exact module is u*<1,Q>, Q=u^2*m.  These bounds are
# imposed in the ambient itself. At E7 (u=infinity), degree(A)<=8 and
# degree(B)<=4 make a=A/h^2 and b=B/h base-regular; the actual all-edge cover
# identifies <1,m> with the required marked module.
e8_orders = {
    "a": ambient["A_powers"],
    "b": ambient["B_powers"],
}
assert min(e8_orders["a"]) >= 1
assert min(e8_orders["b"]) >= 3
e7_infinity_orders = {
    "a": [8-exponent for exponent in ambient["A_powers"]],
    "b": [4-exponent for exponent in ambient["B_powers"]],
}
assert min(e7_infinity_orders["a"]) >= 0
assert min(e7_infinity_orders["b"]) >= 0
assert e7["exceptional_orders"]["m_lower_bound"] == [1, 1, 3, 2, 0, 2, 2]
assert e7["marked_horizontal_frame"]["corrected_identity"] == "Z*m/t=unit/W"
assert "no additional linear condition row" in e7["compiler_conclusion"]

collision = rr["collision_condition"]
condition_matrix = matrix(QQ, [[QQ(value) for value in row] for row in collision["matrix"]])
assert condition_matrix.nrows() == 8 and condition_matrix.ncols() == 10
assert condition_matrix.rank() == 8
assert rr["kernel"]["dimension"] == 2
assert rr["claims"]["h0_D"] == 2

ambient_basis = tuple(range(ambient["dimension"]))
collision_block = quotient_condition(
    "smooth_P1_O_collision",
    ambient_basis,
    lambda column: tuple(condition_matrix[row, column] for row in range(condition_matrix.nrows())),
    tuple("collision_{}".format(row) for row in range(condition_matrix.nrows())),
    "exact h^2 collision quotient derived in assemble_h92_q6_global_rr.sage",
)
compiler_replay = compile_resolved_conditions(
    ambient_basis, (collision_block,), complete=True, compute_kernel=False
)
assert compiler_replay["condition_matrix"] == condition_matrix
assert (compiler_replay["ambient_dimension"], compiler_replay["condition_rows"],
        compiler_replay["rank"], compiler_replay["kernel_dimension"]) == (10, 8, 8, 2)
assert compiler_replay["h0_certified"]

payload = {
    "schema": "elkies-k3.h92-q6-actual-resolved-rr-cover.v1",
    "status": "PASS_EXACT_Q6_ACTUAL_RESOLVED_RR_COVER",
    "inputs": {
        "global_rr": {"path": str(args.rr.relative_to(ROOT)), "sha256": digest(args.rr)},
        "actual_e7_all_edge_module": {"path": str(args.e7.relative_to(ROOT)), "sha256": digest(args.e7)},
        "e8_module": {"path": str(args.e8.relative_to(ROOT)), "sha256": digest(args.e8)},
        "smooth_collision_module": {"path": str(args.smooth.relative_to(ROOT)), "sha256": digest(args.smooth)},
    },
    "resolved_cover": {
        "E8": "ambient bounds ord_u(a)>=1 and ord_u(b)>=3 give u*<1,Q>",
        "E7": "actual all-edge H92 module cover: m has no exceptional pole, and m=unit/W is the marked -P1 frame; ambient coefficients are regular at u=infinity",
        "smooth_P1_O": "the displayed h^2 collision quotient is the sole finite condition block",
    },
    "vertical_condition_matrix": {
        "ambient_dimension": 10,
        "rows": 8,
        "rank": 8,
        "codimension": 8,
        "kernel_dimension": 2,
        "h0_D": 2,
    },
    "compiler_replay": {
        "core": {"path": str(CORE.relative_to(ROOT)), "sha256": digest(CORE)},
        "condition_blocks": [collision_block["name"]],
        "complete_resolved_chart_cover": compiler_replay["complete_resolved_chart_cover"],
        "rank": int(compiler_replay["rank"]),
        "kernel_dimension": int(compiler_replay["kernel_dimension"]),
        "kernel_basis_materialized": compiler_replay["kernel_basis"] is not None,
        "h0_certified": compiler_replay["h0_certified"],
    },
    "conclusion": "The first H3 q=6 pencil has an actual resolved-chart local-module cover. Its exact rank-eight condition matrix has a two-dimensional kernel.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6ACTUALRRCOVER|ambient=10|rank=8|codimension=8|kernel=2|h0=2|"
    "status=PASS_EXACT_Q6_ACTUAL_RESOLVED_RR_COVER",
    flush=True,
)
