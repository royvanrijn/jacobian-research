#!/usr/bin/env sage -python
"""
status: ACTIVE_PROOF_AID
claim: exact preflight for H3-04 D12/MW5 --q6/orbit42--> A11/MW6
parent: artifacts/local/elkies-k3/q24-d13-to-d12-component-valuation-qq.json
output: artifacts/local/elkies-k3/q24-d12-to-a11-orbit42-divval-preflight.json

This script intentionally DOES NOT construct the RR pencil yet.

It pins the corrected orbit42 profile and the exact D12 parent equation:
  mw = (-1,0,-1,-1,0)
  height = 7
  correction = 3
  P.O = 3
  fibre_twist = 0
  D42 = O + P + V

It then certifies the exact I8* local Weierstrass profile (2,3,14) on the
newly certified D12 equation and exports the exact abstract D12 vertical
root vector that must be transported to the PHYSICAL I8* components.

Any correction=1 / P.O=2 / denominator-degree-2 branch is rejected loudly.

Final equation status is deliberately reserved for the later component-
valuation RR compiler:
  PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR
"""

import argparse
import json
import subprocess
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing, matrix, sage_eval, vector


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents/jacobian-research",
        home / "jacobian-research",
        home / "src/jacobian-research",
        home / "git/jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/local/elkies-k3").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"
SCRIPTS = ROOT / "elkies-k3/scripts"

PARENT = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
BACKWARD = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
BRIDGE = LOCAL / "q24-orbit42-current-equation-bridge.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-d12-to-a11-orbit42-divval-preflight.json"
)

if not PARENT.exists():
    raise SystemExit(f"Missing exact D12 parent: {PARENT}")

parent = json.loads(PARENT.read_text())
if parent.get("status") != "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR":
    raise SystemExit(
        "D12 parent is not the exact passing q24 component-valuation artifact"
    )

# Generate the exact suffix planning manifest if needed.  This is target
# lattice metadata only; the equation parent remains PARENT above.
if not BACKWARD.exists():
    producer = SCRIPTS / "build_h3_r17_backward_exact_lift_manifest.sage"
    if not producer.exists():
        raise SystemExit(f"Missing backward-manifest producer: {producer}")
    subprocess.run(
        ["sage", "-python", str(producer)],
        cwd=str(ROOT),
        check=True,
    )

manifest = json.loads(BACKWARD.read_text())
if manifest.get("status") != "PASS_H3_R17_BACKWARD_EXACT_LIFT_MANIFEST":
    raise SystemExit("Backward exact lift manifest is not passing")

steps = [
    row for row in manifest["forward_steps"]
    if int(row["orbit"]) == 42
    and row["parent"] == "D12/MW5"
    and row["child"] == "A11/MW6"
]
if len(steps) != 1:
    raise ArithmeticError(f"Expected one orbit42 D12->A11 step, got {len(steps)}")
step = steps[0]
h = step["horizontal"]

# -------------------------------------------------------------------------
# Corrected orbit42 profile: hard gates against every obsolete branch.
# -------------------------------------------------------------------------
expected_mw = [-1, 0, -1, -1, 0]

checks = {
    "q": int(step["q"]) == 6,
    "orbit": int(step["orbit"]) == 42,
    "child_root_data": list(map(int, step["child_root_data"])) == [11, 132, 12],
    "child_mw_rank": int(step["child_mw_rank"]) == 6,
    "mw": list(map(int, h["mw_projection"])) == expected_mw,
    "height": QQ(h["height"]) == QQ(7),
    "correction": QQ(h["local_correction"]) == QQ(3),
    "P_dot_O": int(h["P_dot_O"]) == 3,
    "fibre_twist": int(h["fibre_twist"]) == 0,
}
bad = [name for name, ok in checks.items() if not ok]
if bad:
    raise ArithmeticError(
        "CORRECTED_ORBIT42_PROFILE_FAILED: " + ",".join(bad)
    )

# Explicit anti-regression guards.
if QQ(h["local_correction"]) == QQ(1):
    raise ArithmeticError("OBSOLETE_ORBIT42_BRANCH: correction=1")
if int(h["P_dot_O"]) == 2:
    raise ArithmeticError("OBSOLETE_ORBIT42_BRANCH: P.O=2")
if int(h["fibre_twist"]) != 0:
    raise ArithmeticError(
        "OBSOLETE_ORBIT42_BRANCH: expected D42=O+P+V with no fibre twist"
    )

vertical_root = vector(ZZ, h["vertical_root_coefficients"])
if len(vertical_root) != 12:
    raise ArithmeticError(
        f"orbit42 vertical root vector has length {len(vertical_root)}, expected 12"
    )

print(
    "Q42DIVVALQQ|stage=PROFILE|"
    "mw=-1,0,-1,-1,0|height=7|correction=3|PdotO=3|"
    "fibre_twist=0|formula=D42=O+P+V|"
    f"vertical_L1={sum(abs(int(v)) for v in vertical_root)}|"
    f"vertical_support={sum(bool(v) for v in vertical_root)}|"
    "status=PASS_CORRECTED_ORBIT42_PROFILE",
    flush=True,
)

# Cross-check the current-equation bridge when available.  It is an exact
# lattice/marking certificate and must agree with the same corrected profile.
bridge_record = None
if BRIDGE.exists():
    bridge = json.loads(BRIDGE.read_text())
    if bridge.get("status") != "PASS_Q24_ORBIT42_CURRENT_EQUATION_LATTICE_BRIDGE":
        raise ArithmeticError("current-equation orbit42 bridge exists but is not passing")
    d12 = bridge["D12"]
    if list(map(int, d12["orbit42_mw_projection"])) != expected_mw:
        raise ArithmeticError("bridge orbit42 MW class disagrees with backward manifest")
    if QQ(d12["orbit42_height"]) != QQ(7):
        raise ArithmeticError("bridge orbit42 height disagrees")
    if QQ(d12["orbit42_local_correction"]) != QQ(3):
        raise ArithmeticError("bridge orbit42 correction disagrees")
    if int(d12["orbit42_P_dot_O"]) != 3:
        raise ArithmeticError("bridge orbit42 P.O disagrees")
    selected = d12["selected_section_representative"]
    if int(selected["fibre_twist"]) != 0:
        raise ArithmeticError("bridge selected representative has nonzero fibre twist")
    bridge_record = {
        "selected_section_root_coordinates": selected["root_coordinates"],
        "selected_vertical_root_coefficients": selected[
            "vertical_root_coefficients"
        ],
        "selected_dual_pairing": selected["dual_pairing"],
    }

print(
    "Q42DIVVALQQ|stage=MARKING_REGRESSION|"
    f"bridge={int(bridge_record is not None)}|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# Exact parent equation only from the newly certified D12 artifact.
# -------------------------------------------------------------------------
child = parent["child"]
R = PolynomialRing(QQ, "V")
V = R.gen()

A = R([QQ(v) for v in child["minimal_A_coefficients_low_to_high"]])
B = R([QQ(v) for v in child["minimal_B_coefficients_low_to_high"]])
Delta = R([QQ(v) for v in child["minimal_discriminant_coefficients_low_to_high"]])

if not A or not B or not Delta:
    raise ArithmeticError("exact D12 parent has a zero Weierstrass datum")

i8 = [
    item for item in child["finite_fibres"]
    if str(item["kodaira"]) == "I8*" and int(item["degree"]) == 1
]
if len(i8) != 1:
    raise ArithmeticError(
        f"expected one finite rational I8* fibre, found {len(i8)}"
    )

factor = R(sage_eval(str(i8[0]["factor"]), locals={"V": V}))
if factor.degree() != 1:
    raise ArithmeticError("D12 I8* factor is not linear")
alpha = QQ(-factor[0] / factor[1])

T = PolynomialRing(QQ, "t")
t = T.gen()
Al = T(A(alpha + t))
Bl = T(B(alpha + t))
Dl = T(Delta(alpha + t))

orders = (
    int(Al.valuation()),
    int(Bl.valuation()),
    int(Dl.valuation()),
)
if orders != (2, 3, 14):
    raise ArithmeticError(
        f"exact D12 local I8* orders are {orders}, expected (2,3,14)"
    )

print(
    "Q42DIVVALQQ|stage=PARENT|"
    f"Adeg={A.degree()}|Bdeg={B.degree()}|Ddeg={Delta.degree()}|"
    f"I8star_root={alpha}|orders=2,3,14|"
    "status=PASS_EXACT_D12_PARENT",
    flush=True,
)

# -------------------------------------------------------------------------
# Root-frame integrity.  The next compiler must NOT guess a physical
# component orientation from the Kodaira symbol.
# -------------------------------------------------------------------------
artifact_path = ROOT / step["artifact"]
if not artifact_path.exists():
    raise SystemExit(f"Missing orbit42 neighbor artifact: {artifact_path}")

q6 = json.loads(artifact_path.read_text())
if q6.get("status") != "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS":
    raise ArithmeticError("orbit42 neighbor artifact is not passing")

frame_path = ROOT / q6["frame"]
if not frame_path.exists():
    raise SystemExit(f"Missing D12 root-adapted frame: {frame_path}")

G = matrix(ZZ, [
    [ZZ(v) for v in line.split()]
    for line in frame_path.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
if G.dimensions() != (17, 17) or G.det() != 948:
    raise ArithmeticError("D12 parent frame has wrong dimension/determinant")

D12 = G[:12, :12]
if D12.det() != 4:
    raise ArithmeticError(f"D12 root block determinant is {D12.det()}, expected 4")
if any(D12[i, i] != 2 for i in range(12)):
    raise ArithmeticError("D12 root block is not normalized to simple roots")
if any(
    D12[i, j] not in (0, -1)
    for i in range(12)
    for j in range(12)
    if i != j
):
    raise ArithmeticError("D12 root block has non-Dynkin off-diagonal entries")

root_edges = [
    [i, j]
    for i in range(12)
    for j in range(i + 1, 12)
    if D12[i, j] == -1
]
if len(root_edges) != 11:
    raise ArithmeticError(
        f"D12 simple-root graph has {len(root_edges)} edges, expected 11"
    )

print(
    "Q42DIVVALQQ|stage=ROOT_FRAME|"
    f"rank=12|det=4|edges={len(root_edges)}|"
    f"vertical={','.join(map(str, vertical_root))}|"
    "status=PASS_ABSTRACT_D12_ROOT_FRAME",
    flush=True,
)

# No ambient dimension is guessed here.  The next stage must:
#   1. explicitly resolve the exact I8* germ;
#   2. match the physical affine-D12 graph to THIS exact root frame;
#   3. transport vertical_root to physical effective components;
#   4. derive the actual line-bundle/trivialization bounds;
#   5. compile those exact conditions and require kernel=2.
#
# Only after that may elliptic_neighbor_compiler.sage be used to compile
# the child equation.

payload = {
    "schema": "elkies-k3.h3-q24-d12-orbit42-divval-preflight.v1",
    "status": "PASS_Q42_DIVVAL_PREFLIGHT",
    "parent": {
        "artifact": str(PARENT.relative_to(ROOT)),
        "status": parent["status"],
        "A_coefficients_low_to_high": [str(v) for v in A.list()],
        "B_coefficients_low_to_high": [str(v) for v in B.list()],
        "discriminant_coefficients_low_to_high": [str(v) for v in Delta.list()],
        "I8star_factor": str(factor),
        "I8star_root": str(alpha),
        "I8star_local_orders": list(orders),
    },
    "orbit42": {
        "q_search_label": 6,
        "actual_old_fibre_degree": 2,
        "historical_orbit": 42,
        "mw_projection": expected_mw,
        "height": "7",
        "local_correction": "3",
        "P_dot_O": 3,
        "fibre_twist": 0,
        "divisor_formula": "D42 = O + P + V",
        "vertical_root_coefficients_abstract_D12": [
            int(v) for v in vertical_root
        ],
        "child_root_data": [11, 132, 12],
        "child_ade": "A11",
        "child_mw_rank": 6,
    },
    "abstract_D12_marking": {
        "frame": str(frame_path.relative_to(ROOT)),
        "root_gram": [
            [int(v) for v in row]
            for row in D12.rows()
        ],
        "root_edges": root_edges,
        "bridge_regression": bridge_record,
    },
    "acceptance_contract": {
        "rr": {
            "ambient_dimension": "DERIVE_EXACTLY",
            "smooth_or_global_conditions": "DERIVE_EXACTLY",
            "resolved_component_rank": "DERIVE_EXACTLY",
            "kernel_dimension": 2,
            "h0": 2,
        },
        "quartic_degree_allowed": [3, 4],
        "child_root_data": [11, 132, 12],
        "child_ade": "A11",
        "child_root_rank": 11,
        "child_root_determinant": 12,
        "child_euler_number": 24,
        "child_mw_rank": 6,
        "final_status_reserved": (
            "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
        ),
    },
    "next_required": {
        "stage": "PHYSICAL_I8STAR_COMPONENT_MARKING",
        "instruction": (
            "Resolve the exact I8* germ on this D12 equation, build its affine-D12 "
            "component graph, and match the physical effective components to the "
            "stored abstract D12 simple-root frame. Carry both spinor-arm "
            "orientations if graph data alone does not distinguish them. Then "
            "transport the orbit42 vertical vector and derive the resolved RR "
            "line-bundle conditions. Do not use the historical zero-pole section "
            "generation route."
        ),
    },
    "proof_boundary": (
        "Exact equation parent and exact corrected orbit42 lattice/divisor profile "
        "only. No RR ambient dimension, physical I8* component orientation, "
        "kernel dimension, quartic, or child equation is claimed here."
    ),
}

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q42DIVVALQQ_RESULT|"
    "parent=D12|profile=CORRECTED|I8star=PASS|"
    "next=PHYSICAL_I8STAR_COMPONENT_MARKING|"
    "status=PASS_Q42_DIVVAL_PREFLIGHT",
    flush=True,
)
