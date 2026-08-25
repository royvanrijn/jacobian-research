#!/usr/bin/env sage -python
"""Audit exact explicit curves as Abel--Jacobi carriers in the selected A11 marking.

The selected equation model uses the R3-zero D12 frame.  This checker keeps
that zero choice fixed, transports every already-recorded explicit (-2)-curve
through the orbit64 transition, and tests whether a single carrier plus the
eighteen exact identity-shell points gives the bridge
``M=(1,0,0,0,0,1)``.  It is exact integral lattice arithmetic only.
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-a11-explicit-aj-carriers.json",
)
args = parser.parse_args()

SPIN = LOCAL / "q24-orbit42-spinor-zero-profiles.json"
Q6 = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
CURVES = LOCAL / "q24-downstream-lift/explicit-curves-a11-span-p100003.json"
BRIDGE = LOCAL / "q24-a11-target-coset-bridge.json"
FRAME = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
INPUTS = (SPIN, Q6, CURVES, BRIDGE, FRAME)
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
bridge = json.loads(BRIDGE.read_text())
assert spin["status"] == "PASS_Q24_ORBIT42_EXACT_SPINOR_ZERO_PROFILES"
assert q6["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert bridge["status"] == "PASS_EXACT_A11_TARGET_COSET_BRIDGE"

profiles = {row["zero"]: row for row in spin["profiles"]}
parent_frame = load_matrix(FRAME)
assert matrix(ZZ, profiles["R3"]["frame"]) == parent_frame
parent_basis = matrix(ZZ, profiles["R3"]["parent_to_child_basis"])
record = next(row for row in q6["neighbors"] if int(row["orbit_index"]) == 64)
transition = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, record["child_root_adapted_basis"])
) * matrix(ZZ, record["neighbor_basis"])
assert abs(parent_basis.det()) == abs(transition.det()) == 1
transport_inverse = (
    parent_basis.inverse().change_ring(ZZ)
    * transition.inverse().change_ring(ZZ)
)

shell = [
    vector(ZZ, row)
    for row in bridge["exact_identity_shell"]["MW_vectors_in_equation_order"]
]
shell_module = matrix(ZZ, [entries(row) for row in shell]).row_module()
assert shell_module.rank() == 5
target = vector(ZZ, bridge["selected_bridge"]["mw"])
assert target == vector(ZZ, (1, 0, 0, 0, 0, 1))

rows = []
single_carriers = []
for raw in curves["explicit_curve_records"]:
    if int(raw["square"]) != -2:
        continue
    parent = vector(ZZ, raw["class"]) * parent_basis.inverse().change_ring(ZZ)
    child = parent * transition.inverse().change_ring(ZZ)
    mw = vector(ZZ, child[-6:])
    row = {
        "name": raw["name"],
        "origin": raw["origin"],
        "A11_degree": int(child[1]),
        "A11_MW_Abel_Jacobi": entries(mw),
        "child_coordinates": entries(child),
        "D12_degree": int(parent[1]),
        "D12_MW_Abel_Jacobi": entries(parent[-5:]),
        "parent_coordinates": entries(parent),
    }
    coefficient = None
    if mw[-1] in (-1, 1):
        coefficient = ZZ(target[-1] // mw[-1])
        residual = target - coefficient * mw
        if residual in shell_module:
            word = shell_module.coordinate_vector(residual)
            row["single_carrier_bridge"] = {
                "carrier_coefficient": int(coefficient),
                "shell_basis_coefficients": entries(word),
                "verified_sum": entries(coefficient * mw + word * shell_module.basis_matrix()),
            }
            single_carriers.append(row)
    rows.append(row)

rows.sort(key=lambda row: (row["A11_degree"], row["name"]))
single_carriers.sort(key=lambda row: (row["A11_degree"], row["name"]))
positive_single_carriers = [row for row in single_carriers if row["A11_degree"] > 0]


def integral_word(generators, wanted):
    """Return coefficients in the supplied generators, or ``None``."""

    generator_matrix = matrix(ZZ, [entries(row) for row in generators])
    hnf, transform = generator_matrix.hermite_form(
        include_zero_rows=False, transformation=True
    )
    if hnf.nrows() != hnf.ncols():
        return None
    basis_coefficients = hnf.solve_left(wanted)
    if not all(value in ZZ for value in basis_coefficients):
        return None
    coefficients = vector(ZZ, basis_coefficients) * transform
    assert coefficients * generator_matrix == wanted
    return coefficients


# A group inverse costs no new trace, so score a subset by the sum of the
# degrees of the distinct curves whose AJ points must be computed.  Search
# subsets of size at most three; this is exhaustive for the stated finite
# candidate list and bound, not a claim about all curves on the surface.
positive_carriers = [
    row for row in rows
    if row["A11_degree"] > 0 and row["A11_MW_Abel_Jacobi"][-1] != 0
]
subset_routes = []
for size in range(1, min(3, len(positive_carriers)) + 1):
    for subset in itertools.combinations(positive_carriers, size):
        generators = shell + [vector(ZZ, row["A11_MW_Abel_Jacobi"]) for row in subset]
        coefficients = integral_word(generators, target)
        if coefficients is None:
            continue
        carrier_coefficients = coefficients[len(shell) :]
        if any(value == 0 for value in carrier_coefficients):
            continue
        subset_routes.append(
            {
                "names": [row["name"] for row in subset],
                "degrees": [row["A11_degree"] for row in subset],
                "distinct_trace_degree_sum": sum(row["A11_degree"] for row in subset),
                "shell_coefficients": entries(coefficients[: len(shell)]),
                "carrier_coefficients": entries(carrier_coefficients),
                "verified_sum": entries(target),
            }
        )
subset_routes.sort(
    key=lambda row: (
        row["distinct_trace_degree_sum"], len(row["names"]), tuple(row["names"])
    )
)

payload = {
    "schema": "elkies-k3.h3-a11-explicit-aj-carriers.v1",
    "status": "PASS_EXACT_A11_EXPLICIT_CARRIER_SUBSET_AUDIT",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "selected_marking": {"D12_zero": "R3", "A11_orbit": 64},
    "target_MW": entries(target),
    "shell_rank": int(shell_module.rank()),
    "explicit_curve_count": len(rows),
    "single_carrier_count": len(single_carriers),
    "best_single_carrier": single_carriers[0] if single_carriers else None,
    "single_carriers": single_carriers,
    "positive_single_carrier_count": len(positive_single_carriers),
    "best_positive_single_carrier": (
        positive_single_carriers[0] if positive_single_carriers else None
    ),
    "positive_nonzero_sixth_carrier_count": len(positive_carriers),
    "subset_search_maximum_size": 3,
    "subset_route_count": len(subset_routes),
    "best_positive_subset_route": subset_routes[0] if subset_routes else None,
    "positive_subset_routes": subset_routes,
    "curves": rows,
    "proof_boundary": (
        "Exact transport and integral subgroup membership in the selected R3-zero "
        "orbit64 A11 marking. A retained carrier still requires an equation-level "
        "restriction and Abel--Jacobi trace; absence would apply only to the stored "
        "finite list of explicit curves. Positive A11 degree alone does not "
        "certify that a stored root orientation has an equation-level effective "
        "parametrization."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "A11AJCARRIERS|curves={}|single={}|best={}|degree={}|subset_routes={}|"
    "best_subset={}|status={}".format(
        len(rows),
        len(positive_single_carriers),
        "none" if not positive_single_carriers else positive_single_carriers[0]["name"],
        "none" if not positive_single_carriers else positive_single_carriers[0]["A11_degree"],
        len(subset_routes),
        "none" if not subset_routes else "+".join(subset_routes[0]["names"]),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
