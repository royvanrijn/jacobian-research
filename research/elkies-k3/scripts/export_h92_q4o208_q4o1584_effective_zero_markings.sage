#!/usr/bin/env sage -python
"""Export physical effective-zero markings for the literal q4/orbit1584 edge.

The edge is compiled from the exact divisor

    O(C5) + first_I6_affine + r3 + 2*r4 + r5.

For each inherited degree-one curve, this script rebuilds U from the actual
curve, recovers the complete 3A1+A3+D4 component chamber, and transports all
equation-explicit curves and pinned target classes losslessly.
"""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SOURCE = GENERATED / "elkies-k3-h3-q4o208-physical-3a3-marking.json"
CERT = GENERATED / "elkies-k3-h3-q4o208-physical-q4o1584-lateral-certificate.json"
OUTPUT = GENERATED / "elkies-k3-h3-q4o208-physical-q4o1584-effective-zero-markings.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))

FIBRE = vector(ZZ, [2, 2, 1, 0, 0, 2, 3, 2, 1, 2, 1, -1, -1, 0, 0, 0, 0, 0, 0])
ZERO_NAMES = (
    "old_A11_component_0",
    "old_A11_component_4",
    "old_A11_component_7",
    "second_I6_affine_component",
)


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


source = json.loads(SOURCE.read_text())
cert = json.loads(CERT.read_text())
assert source["status"] == "PASS_EXACT_Q4O208_PHYSICAL_3A3_MARKING"
assert cert["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
assert cert["candidate_id"]["label"] == "q4o1584_lateral"

source_frame = load_matrix(ROOT / source["frame_output"])
source_gram = block_diagonal_matrix(U2, -source_frame)
explicit = {
    name: vector(ZZ, value)
    for name, value in source["equation_explicit_curves_in_child"].items()
}
assert FIBRE * source_gram * FIBRE == 0
assert all(curve * source_gram * curve == -2 for curve in explicit.values())

# This is the exact, already attached horizontal selected by the scorer.
horizontal = explicit["first_I6_affine_component"]
old_zero = explicit["old_A11_component_5"]
vertical = FIBRE - old_zero - horizontal
expected_vertical = vector(ZZ, [2, 0] + [0, 0, 0, 1, 2, 1, 0, 0, 0] + [0] * 8)
assert vertical == expected_vertical

# Physical components already visible on the source equation.  The D4 list is
# the complete affine I0* diagram; the A3 list is a finite chain whose fourth
# component is F minus the chain.  The two A1 classes have their complementary
# component F-C available exactly.
d4_affine = [explicit[name] for name in (
    "first_I6_affine_component",
    "old_A11_component_1",
    "old_A11_component_5",
    "old_A11_component_6",
    "second_old_I6_I4_missing_component",
)]
a3_chain = [explicit[name] for name in (
    "P1229", "old_A11_component_8", "old_A11_component_9",
)]
a3_cycle = a3_chain + [FIBRE - sum(a3_chain, vector(ZZ, [0] * 19))]
a1_pairs = [
    (explicit["first_old_I6_I4_missing_component"], FIBRE - explicit["first_old_I6_I4_missing_component"]),
    (explicit["old_A11_component_3"], FIBRE - explicit["old_A11_component_3"]),
]

outputs = {}
for zero_name in ZERO_NAMES:
    zero = explicit[zero_name]
    assert zero * source_gram * FIBRE == 1
    u_rows = matrix(ZZ, [list(FIBRE), list(zero + FIBRE)])
    kernel = matrix(ZZ, [list(FIBRE * source_gram), list((zero + FIBRE) * source_gram)]).right_kernel_matrix()
    preliminary = u_rows.stack(kernel)
    assert abs(preliminary.det()) == 1
    prelim_gram = preliminary * source_gram * preliminary.transpose()
    assert prelim_gram[:2, :2] == U2 and not any(prelim_gram[:2, 2:].list())
    complement = -prelim_gram[2:, 2:]

    # Select the nonidentity components from each completely or partially
    # visible reducible fibre using the actual zero intersection.
    d4_identity = [curve for curve in d4_affine if zero * source_gram * curve == 1]
    assert len(d4_identity) == 1
    effective_roots = [curve for curve in d4_affine if curve != d4_identity[0]]

    a3_identity = [curve for curve in a3_cycle if zero * source_gram * curve == 1]
    assert len(a3_identity) == 1
    effective_roots.extend(curve for curve in a3_cycle if curve != a3_identity[0])

    for left, right in a1_pairs:
        pair = (left, right)
        identities = [curve for curve in pair if zero * source_gram * curve == 1]
        assert len(identities) == 1
        effective_roots.append(next(curve for curve in pair if curve != identities[0]))
    assert len(effective_roots) == 9

    # The third A1 was not among the 25 inherited curves.  Recover its exact
    # physical sign as the unique root orthogonal to the nine visible simple
    # components and nonnegative on every distinct inherited effective curve.
    half = [vector(ZZ, column) for column in matrix(ZZ, pari(complement).qfminim(2)[2]).columns()]
    root_classes = []
    for tail in half + [-item for item in half]:
        root_classes.append(vector(ZZ, [0, 0] + list(tail)) * preliminary)
    missing = []
    unique_explicit = {tuple(curve): curve for curve in explicit.values()}.values()
    for curve in root_classes:
        if any(curve * source_gram * known for known in effective_roots):
            continue
        if all(curve * source_gram * known >= 0 for known in unique_explicit):
            missing.append(curve)
    assert len({tuple(curve) for curve in missing}) == 1
    effective_roots.append(missing[0])

    # Root-adapted search code represents effective components by -e_i.
    partial = u_rows.stack(matrix(ZZ, [list(-curve) for curve in effective_roots]))
    cartan = -(partial[2:] * source_gram * partial[2:].transpose())
    assert abs(cartan.det()) == 128
    assert ZZ(pari(cartan).qfminim(2)[0]) == 42

    smith, left, right = partial.smith_form()
    assert smith[:, :12] == matrix.identity(ZZ, 12)
    assert not any(smith[:, 12:].list())
    completion = right.inverse().change_ring(ZZ)
    change = block_diagonal_matrix(left.inverse().change_ring(ZZ), matrix.identity(ZZ, 7))
    basis_in_source = change * completion
    assert basis_in_source[:12] == partial
    # Smith completion preserves the displayed rows but need not put its seven
    # new rows in U-perp.  Clear their two U coordinates by determinant-one row
    # operations while leaving the physical roots untouched.
    for index in range(12, 19):
        row = basis_in_source.row(index)
        row -= ZZ(row * source_gram * (zero + FIBRE)) * FIBRE
        row -= ZZ(row * source_gram * FIBRE) * (zero + FIBRE)
        basis_in_source[index] = row
    assert abs(basis_in_source.det()) == 1
    source_in_basis = basis_in_source.inverse().change_ring(ZZ)
    gram = basis_in_source * source_gram * basis_in_source.transpose()
    frame = -gram[2:, 2:]
    assert gram[:2, :2] == U2 and not any(gram[:2, 2:].list())
    assert frame[:10, :10] == cartan

    slug = zero_name.lower()
    frame_path = GENERATED / f"elkies-k3-h3-q4o208-physical-q4o1584-{slug}-frame.txt"
    marking_path = GENERATED / f"elkies-k3-h3-q4o208-physical-q4o1584-{slug}-marking.json"
    frame_path.write_text(
        f"# q4/orbit1584 physical 3A1+A3+D4 frame; zero={zero_name}\n"
        + "\n".join(" ".join(map(str, row)) for row in frame.rows()) + "\n"
    )

    explicit_child = {
        name: entries(curve * source_in_basis)
        for name, curve in explicit.items()
    }
    targets = {
        name: entries(vector(ZZ, value) * source_in_basis)
        for name, value in source["target_fibres_in_root_adapted_hub"].items()
    }
    equation_to_source = matrix(ZZ, source["equation_A11_to_root_adapted_hub_basis"])
    equation_to_child = basis_in_source * equation_to_source
    payload = {
        "schema": "elkies-k3.h3-q4o208-q4o1584-effective-zero-marking.v1",
        "status": "PASS_EXACT_Q4O1584_PHYSICAL_EFFECTIVE_ZERO_MARKING",
        "hub": f"q4o208_q4o1584_physical_3A1_A3_D4_{zero_name}_zero",
        "zero": zero_name,
        "root_data": [10, 42, 128],
        "ade": "3A1+A3+D4",
        "frame_output": str(frame_path.relative_to(ROOT)),
        "frame_sha256": sha256(frame_path),
        "physical_simple_components_in_source": [entries(curve) for curve in effective_roots],
        "basis_in_source": rows(basis_in_source),
        "source_in_basis": rows(source_in_basis),
        "equation_A11_to_root_adapted_hub_basis": rows(equation_to_child),
        "target_fibres_in_root_adapted_hub": targets,
        "equation_explicit_curves_in_child": explicit_child,
        "prefix_operational_score": (source.get("prefix_operational_score") or 0) - 1580,
        "literal_divisor": {
            "fibre_in_source": entries(FIBRE),
            "zero": "old_A11_component_5",
            "horizontal": "first_I6_affine_component",
            "physical_root_coefficients": [0, 0, 0, 1, 2, 1, 0, 0, 0],
            "P_dot_O": 0,
            "expected_RR_ambient": 4,
        },
        "proof_boundary": (
            "Exact physical effective zero, complete 3A1+A3+D4 component chamber, "
            "literal attached horizontal, all inherited equation curves, targets, and "
            "bidirectional unimodular NS transport. Successor edges require independent "
            "finite-horizontal-wall certification."
        ),
        "inputs": {
            "paths": [str(SOURCE.relative_to(ROOT)), str(CERT.relative_to(ROOT))],
            "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (SOURCE, CERT)},
        },
    }
    marking_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    outputs[zero_name] = {
        "frame": str(frame_path.relative_to(ROOT)),
        "marking": str(marking_path.relative_to(ROOT)),
        "frame_max_abs": max(abs(int(value)) for value in basis_in_source.list()),
    }

summary = {
    "schema": "elkies-k3.h3-q4o208-q4o1584-effective-zero-markings.v1",
    "status": "PASS_EXACT_Q4O1584_ALL_PHYSICAL_EFFECTIVE_ZERO_MARKINGS",
    "edge": {
        "q": 4,
        "orbit": 1584,
        "ade": "3A1+A3+D4",
        "mw_rank": 7,
        "literal_divisor": "O(C5)+first_I6_affine+r3+2r4+r5",
    },
    "outputs": outputs,
    "proof_boundary": "Each output is a full physical marked state, not an ADE-only identification.",
}
OUTPUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(f"Q4O1584ZEROS|count={len(outputs)}|status={summary['status']}|output={OUTPUT}")
