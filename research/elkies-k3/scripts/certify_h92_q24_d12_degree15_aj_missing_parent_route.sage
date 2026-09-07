#!/usr/bin/env sage -python
"""Certify a degree-14 AJ route to the missing D12 MW coordinate.

In the compatible R3-zero D12 marking, ``close_P24`` has Abel--Jacobi vector
``(-6,0,-7,0,1)``.  Its resolved pole branch is an actual degree-14 curve;
the exact spinor pair supplies ``+/- (1,0,0,0,0)``, saturating the parity gap
left by the eighteen identity points.  One degree-14 trace plus those exact
zero-pole sections constructs the primitive missing parent section
``(0,0,0,0,1)``.  This certificate also
enumerates every closest-root representative of that section and transports
it through the selected orbit64 A11 transition.

The lattice word is exact.  The degree-14 map and its P.O=104 reconstructed
AJ section are presently certified at the pinned good prime; characteristic-
zero lifting of that trace remains to be executed.
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import IntegralLattice, QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-q24-d12-degree14-aj-missing-parent-route.json",
)
args = parser.parse_args()

SPIN = LOCAL / "q24-orbit42-spinor-zero-profiles.json"
Q6 = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
CURVES = LOCAL / "q24-downstream-lift/explicit-curves-a11-span-p100003.json"
IDENTITY = LOCAL / "q24-orbit42-identity-halving-audit.json"
ALTERNATIVES = LOCAL / "q24-a11-missing-direction-alternatives.json"
AJ14 = LOCAL / "q24-close-p24-aj14-plus-section-mod100003.json"
FRAME = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
INPUTS = (SPIN, Q6, CURVES, IDENTITY, ALTERNATIVES, AJ14, FRAME)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


spin = json.loads(SPIN.read_text())
q6 = json.loads(Q6.read_text())
curves = json.loads(CURVES.read_text())
identity = json.loads(IDENTITY.read_text())
alternatives = json.loads(ALTERNATIVES.read_text())
aj14 = json.loads(AJ14.read_text())
profiles = {row["zero"]: row for row in spin["profiles"]}
parent_frame = load_matrix(FRAME)
parent_basis = matrix(ZZ, profiles["R3"]["parent_to_child_basis"])
assert matrix(ZZ, profiles["R3"]["frame"]) == parent_frame
record = next(row for row in q6["neighbors"] if int(row["orbit_index"]) == 64)
transition = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, record["child_root_adapted_basis"])
) * matrix(ZZ, record["neighbor_basis"])
assert abs(parent_basis.det()) == abs(transition.det()) == 1

close_raw = next(row for row in curves["explicit_curve_records"] if row["name"] == "close_P24")
close_parent = vector(ZZ, close_raw["class"]) * parent_basis.inverse().change_ring(ZZ)
close_mw = vector(ZZ, close_parent[-5:])
assert int(close_parent[1]) == 15
assert close_mw == vector(ZZ, (-6, 0, -7, 0, 1))
assert alternatives["status"] == "PASS_EXACT_A11_MISSING_DIRECTION_ALTERNATIVES_AUDIT"
assert alternatives["exact_zero_pole_shell"]["spinor_vectors"] == [
    [-1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0],
]
assert aj14["status"] == "PASS_Q24_CLOSE_P24_AJ14_SECTION_RECONSTRUCTION_MODP"
assert aj14["prime"] == 100003
assert aj14["section"]["P_dot_O"] == 104

identity_shell = [vector(ZZ, row) for row in identity["exact_model_R3_zero"]["identity_vectors"]]
assert len(identity_shell) == 18 and all(row[-1] == 0 for row in identity_shell)
spinor = vector(ZZ, (1, 0, 0, 0, 0))
shell = identity_shell + [spinor]
target = vector(ZZ, (0, 0, 0, 0, 1))
residual = target - close_mw


def shortest_subset_word(generators, wanted):
    best = None
    for size in range(1, 6):
        for indices in itertools.combinations(range(len(generators)), size):
            basis = matrix(ZZ, [list(generators[index][:-1]) for index in indices])
            if basis.rank() != size:
                continue
            try:
                local = basis.solve_left(vector(ZZ, wanted[:-1]))
            except ValueError:
                continue
            if not all(value in ZZ for value in local):
                continue
            coefficients = vector(ZZ, [0] * len(generators))
            for index, value in zip(indices, local):
                coefficients[index] = ZZ(value)
            if sum((coefficients[i] * generators[i] for i in range(len(generators))), vector(ZZ, [0] * 5)) != wanted:
                continue
            score = (size, max(abs(value) for value in coefficients), sum(abs(value) for value in coefficients), tuple(coefficients))
            if best is None or score < best[0]:
                best = (score, coefficients)
        if best is not None:
            break
    return None if best is None else best[1]


shell_coefficients = shortest_subset_word(shell, residual)
if shell_coefficients is None:
    raise ArithmeticError("degree-14 carrier residual is not in the exact parent shell plus spinor")
assert close_mw + sum(
    (shell_coefficients[index] * shell[index] for index in range(len(shell))),
    vector(ZZ, [0] * 5),
) == target

# Enumerate all closest-root representatives of the target parent section.
root_rank = 12
root = parent_frame[:root_rank, :root_rank]
coupling = parent_frame[:root_rank, root_rank:]
tail = parent_frame[root_rank:, root_rank:]
height = tail - coupling.transpose() * root.inverse() * coupling
h = QQ(target * height * target)
base = vector(ZZ, [0] * root_rank + list(target))
dual = vector(QQ, base * parent_frame[:, :root_rank]) * root.inverse()
root_lattice = IntegralLattice(root)
iterator = root_lattice.enumerate_close_vectors(-dual)
minimum = None
representatives = []
for unused in range(100000):
    shift = vector(ZZ, next(iterator))
    lifted = base + vector(ZZ, list(shift) + [0] * 5)
    norm = QQ(lifted * parent_frame * lifted)
    if minimum is None:
        minimum = norm
    elif norm > minimum:
        break
    correction = norm - h
    pole = (norm - 4) / 2
    if pole not in ZZ or pole < 0:
        continue
    section = vector(ZZ, [ZZ(pole) + 1, 1] + list(lifted))
    child = section * transition.inverse().change_ring(ZZ)
    representatives.append(
        {
            "parent_section": entries(section),
            "P_dot_O": int(pole),
            "local_correction": str(correction),
            "A11_degree": int(child[1]),
            "A11_MW_Abel_Jacobi": entries(child[-6:]),
            "A11_child_coordinates": entries(child),
        }
    )

assert h == 12
assert representatives and all(row["P_dot_O"] == 4 for row in representatives)
representatives.sort(key=lambda row: (row["A11_degree"], row["parent_section"]))

payload = {
    "schema": "elkies-k3.h3-q24-d12-degree14-aj-missing-parent-route.v1",
    "status": "PASS_EXACT_Q24_D12_DEGREE14_AJ_MISSING_PARENT_ROUTE",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "trace_curve": {
            "name": "close_P24",
            "raw_NS_intersection_before_basepoint_resolution": int(close_parent[1]),
            "resolved_pole_branch_degree_mod_100003": 14,
            "D12_MW_Abel_Jacobi": entries(close_mw),
            "parent_coordinates": entries(close_parent),
            "reconstructed_AJ_P_dot_O_mod_100003": int(aj14["section"]["P_dot_O"]),
        },
    "parent_bridge_word": {
        "formula": "E5=AJ(close_P24)+sum_i c_i*S_i+c_spin*Spinor",
        "target_D12_MW": entries(target),
        "identity_shell_coefficients": entries(shell_coefficients[:18]),
        "spinor_coefficient_for_vector_1_0_0_0_0": int(shell_coefficients[18]),
        "verified_sum": entries(target),
    },
    "target_section_profile": {
        "height": str(h),
        "P_dot_O": 4,
        "closest_representative_count": len(representatives),
        "representatives_transported_to_A11": representatives,
    },
    "next_exact_computation": {
        "method": "resolved pole-branch restriction and one fibrewise Abel-Jacobi trace",
        "trace_degree": 14,
        "linear_trace_space": "L(15O)",
        "linear_trace_dimension": 15,
        "large_Groebner_required": False,
    },
    "proof_boundary": (
        "Exact NS transport, shell/spinor word, and complete closest-root enumeration, "
        "plus a pinned-good-prime degree-14 trace reconstruction. The AJ trace still "
        "requires characteristic-zero lifting before the parent section is exact."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q24D14PARENT|curve=close_P24|degree=14|target=0,0,0,0,1|trace=L15O|"
    "representatives={}|A11_degrees={}|status={}".format(
        len(representatives),
        ",".join(str(row["A11_degree"]) for row in representatives),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
