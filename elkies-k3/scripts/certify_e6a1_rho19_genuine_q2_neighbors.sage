#!/usr/bin/env sage-python
"""Certify the first genuine q=2 neighbors of the E6+A1 rho-19 K3.

status: ACTIVE_PROOF
claim: complete norm-eight q=2 Weyl census and eighteen nef MW-rank-three frames
inputs: elkies-k3-e6a1-rho19-k3-dissection-v1.json
outputs: elkies-k3-e6a1-rho19-genuine-q2-neighbors-v1.json

The nominal classes ``e+q*f+w`` at norms ``2q`` meet the old zero negatively
and reduce to degree one.  The first genuine degree-two layer is therefore
``2e+2f-w`` with ``w^2=8``.  This script enumerates that layer without forming
the full signed norm shell: dominant Dynkin labels reduce the search to a
rank-two closest-vector enumeration.  It then applies the physical fibre
chamber, classifies every child frame, and gives a complete nefness proof for
the root-rank-fourteen candidates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import product
from pathlib import Path

from sage.all import (
    block_diagonal_matrix,
    CartanMatrix,
    IntegralLattice,
    QQ,
    ZZ,
    floor,
    gcd,
    identity_matrix,
    lcm,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GEN = ROOT / "artifacts/generated-results"
SOURCE = GEN / "elkies-k3-e6a1-rho19-k3-dissection-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-e6a1-rho19-genuine-q2-neighbors-v1.json"

_engine_path = HERE / "exact_neighbor_engine.sage"
exec(compile(_engine_path.read_text(), str(_engine_path), "exec"), globals())


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integer_rows(value):
    return [[int(entry) for entry in row] for row in matrix(ZZ, value).rows()]


def rational_rows_local(value):
    return [[str(entry) for entry in row] for row in matrix(QQ, value).rows()]


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def dot(left, right, gram):
    return vector(left) * gram * vector(right)


def connected_components(cartan):
    cartan = matrix(ZZ, cartan)
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        pending = [min(unseen)]
        unseen.remove(pending[0])
        component = []
        while pending:
            index = pending.pop()
            component.append(index)
            adjacent = [j for j in unseen if cartan[index, j] != 0]
            for j in adjacent:
                unseen.remove(j)
                pending.append(j)
        result.append(tuple(sorted(component)))
    return tuple(result)


def component_label_vectors(cartan, norm_bound):
    """All nonnegative Dynkin labels of dual norm at most norm_bound."""
    inverse = matrix(QQ, cartan).inverse()
    bounds = [
        int(floor((QQ(norm_bound) / inverse[i, i]).sqrt()))
        for i in range(inverse.nrows())
    ]
    labels = []
    for candidate in product(*[range(bound + 1) for bound in bounds]):
        item = vector(ZZ, candidate)
        norm = item * inverse * item
        if norm <= norm_bound:
            labels.append((item, norm))
    return tuple(labels)


def enumerate_dominant_norm_vectors(frame, simple, target_norm):
    """Enumerate Weyl-dominant vectors by labels plus a rank-two CVP."""
    frame = matrix(ZZ, frame)
    simple = matrix(ZZ, simple)
    cartan = simple * frame * simple.transpose()
    components = connected_components(cartan)
    component_labels = []
    for component in components:
        block = cartan.matrix_from_rows_and_columns(component, component)
        records = component_label_vectors(block, target_norm)
        component_labels.append((component, records))

    combined = []
    for choices in product(*[records for _, records in component_labels]):
        total = sum((norm for _, norm in choices), QQ(0))
        if total > target_norm:
            continue
        label = vector(ZZ, simple.nrows())
        for (component, _), (part, _) in zip(component_labels, choices):
            for index, value in zip(component, part):
                label[index] = value
        combined.append((label, total))

    pairing = simple * frame
    smith, left, right = pairing.smith_form()
    if smith != left * pairing * right:
        raise AssertionError("unexpected Smith convention")
    rank = pairing.rank()
    diagonal = tuple(abs(ZZ(smith[i, i])) for i in range(rank))
    kernel = right[:, rank:]
    kernel_gram = kernel.transpose() * frame * kernel
    lattice = IntegralLattice(kernel_gram)

    answers = []
    compatible_labels = 0
    for label, projection_norm in combined:
        rhs = left * label
        if any(rhs[i] % diagonal[i] for i in range(rank)):
            continue
        compatible_labels += 1
        coordinates = vector(ZZ, frame.nrows())
        for i in range(rank):
            coordinates[i] = rhs[i] // diagonal[i]
        particular = right * coordinates
        center = -kernel_gram.inverse() * kernel.transpose() * frame * particular
        for z in lattice.enumerate_close_vectors(center):
            candidate = particular + kernel * vector(ZZ, z)
            norm = candidate * frame * candidate
            if norm > target_norm:
                break
            if norm == target_norm:
                if pairing * candidate != label:
                    raise AssertionError("dominant label reconstruction failed")
                answers.append(candidate)

    answers = sorted({tuple(entries(item)) for item in answers})
    return {
        "vectors": tuple(vector(ZZ, item) for item in answers),
        "cartan": cartan,
        "components": components,
        "component_label_counts": tuple(len(records) for _, records in component_labels),
        "combined_label_count": len(combined),
        "compatible_label_count": compatible_labels,
        "pairing_smith_invariants": diagonal,
        "kernel_gram": kernel_gram,
    }


def highest_root_coefficients(cartan):
    """Return highest-root coefficients in this deterministic simple basis."""
    cartan = matrix(ZZ, cartan)
    # The current deterministic bases have the standard E6,E6,A3 ordering.
    if cartan.nrows() == 6 and cartan.det() == 3:
        return vector(ZZ, [1, 2, 2, 3, 2, 1])
    if cartan.nrows() == 3 and cartan.det() == 4:
        return vector(ZZ, [1, 1, 1])
    raise ValueError("unexpected old fibre component")


def old_fibre_walls(ns, old_fibre):
    walls = [("O", vector(ZZ, [1] + [0] * 18))]
    starts = (2, 8, 14)
    sizes = (6, 6, 3)
    names = ("E6a", "E6b", "A3")
    for name, start, size in zip(names, starts, sizes):
        cartan = -ns[start:start + size, start:start + size]
        components = []
        for i in range(size):
            curve = vector(ZZ, 19)
            curve[start + i] = 1
            components.append(curve)
            walls.append((f"{name}_{i + 1}", curve))
        highest = sum(
            (coefficient * curve for coefficient, curve in zip(
                highest_root_coefficients(cartan), components
            )),
            vector(ZZ, 19),
        )
        walls.append((f"{name}_0", old_fibre - highest))
    return tuple(walls)


ROOT_TYPES = {
    (14, 88, 256): "4A1+2D5",
    (14, 92, 144): "2A2+2D5",
    (14, 136, 48): "2A1+D6+E6",
    (14, 140, 32): "A7+D7",
}


def correction_vector(a, b):
    """Shioda correction for the section a*P0+b*P1."""
    result = vector(QQ, 19)
    residue = a % 3
    if residue:
        node = 0 if residue == 1 else 5
        correction = CartanMatrix(["E", 6]).inverse().row(node)
        result[2:8] = correction
        result[8:14] = correction
    if (a + b) % 2:
        result[14:17] = CartanMatrix(["A", 3]).inverse().row(1)
    return result


def section_class(a, b, phi0, phi1, ns):
    correction = correction_vector(a, b)
    height = QQ(a * a) / 3 + 3 * b * b
    local = -(correction * ns * correction)
    intersection_zero = (height + local - 4) / 2
    if intersection_zero not in ZZ:
        raise AssertionError("section/zero intersection is not integral")
    section = (
        a * phi0 + b * phi1
        + vector(QQ, [1] + [0] * 18)
        + (2 + intersection_zero) * vector(QQ, [0, 1] + [0] * 17)
        - correction
    )
    if not all(value in ZZ for value in section):
        raise AssertionError("Shioda section class is not integral")
    section = vector(ZZ, section)
    old_fibre_local = vector(ZZ, [0, 1] + [0] * 17)
    if section * ns * section != -2 or section * ns * old_fibre_local != 1:
        raise AssertionError("reconstructed section class changed")
    return section


def saturated_mw_data(frame):
    """Compute torsion and the free height lattice after root saturation."""
    frame = matrix(ZZ, frame)
    roots, root_basis, root_data = roots_and_data(frame)
    rank = root_data[0]
    root_smith, _, _ = root_basis.smith_form()
    invariants = tuple(
        abs(ZZ(root_smith[i, i])) for i in range(rank)
    )
    torsion_order = ZZ(1)
    for value in invariants:
        torsion_order *= value

    saturated = root_basis.row_module().saturation().basis_matrix()
    smith, left, right = saturated.smith_form()
    if smith != left * saturated * right:
        raise AssertionError("unexpected saturated-root Smith convention")
    if tuple(abs(ZZ(smith[i, i])) for i in range(rank)) != (1,) * rank:
        raise AssertionError("root saturation is not primitive")
    completion = right.inverse()
    basis = saturated.stack(completion[rank:])
    if abs(basis.det()) != 1:
        raise AssertionError("saturated root completion is not unimodular")
    adapted = basis * frame * basis.transpose()
    root_block = adapted[:rank, :rank]
    coupling = adapted[:rank, rank:]
    tail = adapted[rank:, rank:]
    height = tail - coupling.transpose() * root_block.inverse() * coupling
    scale = lcm(value.denominator() for value in height.list())
    quotient_change = (scale * height).change_ring(ZZ).LLL_gram().transpose()
    if abs(quotient_change.det()) != 1:
        raise AssertionError("MW quotient LLL change is not unimodular")
    full_change = block_diagonal_matrix(identity_matrix(ZZ, rank), quotient_change)
    adapted = full_change * adapted * full_change.transpose()
    root_block = adapted[:rank, :rank]
    coupling = adapted[:rank, rank:]
    tail = adapted[rank:, rank:]
    height = tail - coupling.transpose() * root_block.inverse() * coupling
    return {
        "torsion_order": int(torsion_order),
        "root_smith_invariants": invariants,
        "height": height,
    }


if not SOURCE.exists():
    raise FileNotFoundError(SOURCE)
source = json.loads(SOURCE.read_text())
if not source.get("status", "").startswith("PASS_EXACT"):
    raise ArithmeticError("source dissection is not exact")

ns = matrix(ZZ, source["generic_k3"]["ns_gram"])
frame = matrix(
    ZZ,
    source["generic_k3"]["positive_frame_for_neighbor_search"]["gram"],
)
old_fibre = vector(ZZ, [0, 1] + [0] * 17)
old_zero = vector(ZZ, [1] + [0] * 18)
constraints = matrix(ZZ, [ns * old_fibre, ns * (old_zero + old_fibre)])
frame_basis = constraints.right_kernel().basis_matrix()
transport = matrix(ZZ, [old_fibre, old_zero + old_fibre, *frame_basis.rows()])
if abs(transport.det()) != 1:
    raise ArithmeticError("old U split is not unimodular")
if transport * ns * transport.transpose() != matrix(ZZ, 19, 19, lambda i, j:
        1 if (i, j) in ((0, 1), (1, 0)) else
        (-frame[i - 2, j - 2] if i >= 2 and j >= 2 else 0)):
    raise ArithmeticError("old U split changed")

simple, _, old_cartan_rows = deterministic_simple_roots(frame)
enumeration = enumerate_dominant_norm_vectors(frame, simple, ZZ(8))
dominant_vectors = enumeration["vectors"]
if enumeration["component_label_counts"] != (20, 17, 17):
    raise ArithmeticError(
        f"component label census changed: {enumeration['component_label_counts']}"
    )
if enumeration["combined_label_count"] != 553:
    raise ArithmeticError("combined dominant label census changed")
if enumeration["pairing_smith_invariants"] != (1,) * 14 + (6,):
    raise ArithmeticError("root-pairing Smith data changed")
if (
    enumeration["kernel_gram"][0, 0] != 12
    or enumeration["kernel_gram"][1, 1] != 30
    or abs(enumeration["kernel_gram"][0, 1]) != 18
):
    raise ArithmeticError(f"rank-two kernel changed: {enumeration['kernel_gram']}")
if len(dominant_vectors) != 119:
    raise ArithmeticError("norm-eight Weyl orbit count changed")

# The nominal minimal presentations are obstructed before a neighbor search.
# For D=e+q*f-w with w in the old frame, D.O=1-q.  Removing (q-1)O leaves
# old-fibre degree one.  In particular norm four at q=2 and norm six at q=3
# are section presentations, not genuine q-neighbors.
nominal_obstruction = []
for q in (2, 3):
    nominal_obstruction.append({
        "q": q,
        "required_norm": 2 * q,
        "intersection_with_old_zero": 1 - q,
        "fixed_old_zero_multiplicity": q - 1,
        "degree_after_removal": 1,
    })

walls = old_fibre_walls(ns, old_fibre)
raw_distribution = Counter()
degree_distribution = Counter()
all_records = []
rank_three_records = []
nonprimitive_records = []

# Shioda projections used by the complete section gate.
corr0 = correction_vector(1, 0)
corr1 = correction_vector(0, 1)
p0 = vector(QQ, [0] * 17 + [1, 0])
p1 = vector(QQ, [0] * 18 + [1])
phi0 = p0 - old_zero - 2 * old_fibre + corr0
phi1 = p1 - old_zero - 2 * old_fibre + corr1
height = matrix(QQ, [[QQ(1) / 3, 0], [0, 3]])
if matrix(QQ, [[phi0 * ns * phi0, phi0 * ns * phi1],
               [phi1 * ns * phi0, phi1 * ns * phi1]]) != -height:
    raise ArithmeticError("old Shioda height basis changed")

for orbit, dominant in enumerate(dominant_vectors):
    split_divisor = vector(ZZ, [2, 2] + [-entry for entry in dominant])
    divisor = split_divisor * transport
    if divisor * ns * divisor != 0:
        raise ArithmeticError("norm-eight vector did not give an isotropic class")
    divisor_gcd = gcd(list(divisor))
    if divisor_gcd != 1:
        nonprimitive_records.append({
            "orbit": orbit,
            "dominant_w": entries(dominant),
            "divisor_in_split_basis": entries(split_divisor),
            "divisor_gcd": int(divisor_gcd),
            "primitive_degree": 1,
            "reason": "D=2(e+f-r) for a norm-two root r",
        })
        continue
    if divisor * ns * old_fibre != 2 or divisor * ns * old_zero != 0:
        raise ArithmeticError("genuine q=2 intersection profile changed")

    reduced, reduction = reduce_fixed_components(divisor, walls, ns)
    degree = ZZ(reduced * ns * old_fibre)
    degree_distribution[int(degree)] += 1
    split = primitive_hyperbolic_split(ns, divisor)
    minimized = minimize_child_frame(split["child_frame"])
    root_data = tuple(int(value) for value in minimized["root_data"])
    raw_distribution[root_data] += 1
    record = {
        "orbit": orbit,
        "dominant_w": entries(dominant),
        "divisor_in_split_basis": entries(split_divisor),
        "divisor_in_ns_basis": entries(divisor),
        "divisor_complexity": {
            "max_abs": int(max(abs(value) for value in divisor)),
            "l1": int(sum(abs(value) for value in divisor)),
        },
        "physical_reduction": {
            "sequence": [[name, int(pairing)] for name, pairing in reduction],
            "reduced_divisor_in_ns_basis": entries(reduced),
            "old_fibre_degree": int(degree),
        },
        "child_frame_digest": matrix_digest(split["child_frame"]),
        "child_root_data": list(root_data),
    }

    if root_data[0] == 14:
        if degree != 2 or reduction:
            raise ArithmeticError("rank-three candidate failed physical chamber")
        if root_data not in ROOT_TYPES:
            raise ArithmeticError("unexpected root-rank-fourteen type")
        saturated_data = saturated_mw_data(split["child_frame"])

        # Project the divisor's frame part to the two-dimensional old MW
        # subspace.  Proposition C2 implies that a negative section must lie
        # in (n-delta/2)^T H (n-delta/2)<2.
        pairings = vector(QQ, [divisor * ns * phi0, divisor * ns * phi1])
        delta = -pairings * height.inverse()
        if not all(value in ZZ for value in delta):
            raise ArithmeticError("section-search centre is not integral")
        delta = vector(ZZ, delta)
        candidates = []
        # H is diagonal, so |a-delta_0/2|<sqrt(6) and
        # |b-delta_1/2|<sqrt(2/3).  The fixed radius six is a strict superset.
        for a in range(int(floor(QQ(delta[0]) / 2)) - 6,
                       int(floor(QQ(delta[0]) / 2)) + 7):
            for b in range(int(floor(QQ(delta[1]) / 2)) - 2,
                           int(floor(QQ(delta[1]) / 2)) + 3):
                offset = vector(QQ, [a, b]) - vector(QQ, delta) / 2
                if offset * height * offset >= 2:
                    continue
                section = section_class(a, b, phi0, phi1, ns)
                pairing = ZZ(divisor * ns * section)
                candidates.append({
                    "mw_coordinates": [a, b],
                    "pairing": int(pairing),
                })
        if not candidates or min(item["pairing"] for item in candidates) < 0:
            raise ArithmeticError("rank-three candidate has a negative section")

        mw_height = saturated_data["height"]
        torsion_order = saturated_data["torsion_order"]
        expected_regulator = QQ(36 * torsion_order**2) / QQ(root_data[2])
        if mw_height.det() != expected_regulator:
            raise ArithmeticError("child regulator mutation changed")
        record.update({
            "root_type": ROOT_TYPES[root_data],
            "mw_rank_if_rho_19": 3,
            "root_lattice_primitive": torsion_order == 1,
            "root_smith_invariants": [
                int(value) for value in saturated_data["root_smith_invariants"]
            ],
            "torsion": "trivial" if torsion_order == 1 else f"order {torsion_order}",
            "torsion_order": torsion_order,
            "mw_height_gram": rational_rows_local(mw_height),
            "mw_regulator": str(mw_height.det()),
            "complete_section_gate": {
                "centre_delta": entries(delta),
                "strict_ellipsoid": "(n-delta/2)^T*diag(1/3,3)*(n-delta/2)<2",
                "tested_sections": candidates,
                "minimum_pairing": min(item["pairing"] for item in candidates),
            },
            "negative_bisection_exclusion": (
                "If an irreducible (-2)-curve C of old degree two had D.C<0, "
                "then D-C=lF; Hodge index gives D.C=-1 and (D-C)^2=0, "
                "while C^2=-4l=-2 forces l=1/2, contradicting primitive F."
            ),
            "nef": True,
        })
        rank_three_records.append(record)
    all_records.append(record)

expected_distribution = {
    (15, 156, 36): 26,
    (16, 132, 108): 16,
    (16, 120, 216): 8,
    (15, 190, 32): 8,
    (15, 96, 224): 8,
    (15, 132, 80): 8,
    (14, 140, 32): 6,
    (15, 238, 8): 6,
    (15, 206, 32): 5,
    (16, 156, 81): 4,
    (14, 88, 256): 4,
    (16, 148, 192): 4,
    (14, 92, 144): 4,
    (14, 136, 48): 4,
    (16, 172, 48): 3,
    (15, 76, 576): 2,
}
if dict(raw_distribution) != expected_distribution:
    raise ArithmeticError(
        f"norm-eight child root distribution changed: {dict(raw_distribution)}"
    )
if dict(degree_distribution) != {0: 12, 1: 14, 2: 90}:
    raise ArithmeticError(
        f"physical degree distribution changed: {dict(degree_distribution)}"
    )
if len(rank_three_records) != 18:
    raise ArithmeticError("nef rank-three candidate count changed")

best = min(
    rank_three_records,
    key=lambda item: (
        item["divisor_complexity"]["max_abs"],
        item["divisor_complexity"]["l1"],
        item["orbit"],
    ),
)
fibre_simple_best = min(
    (item for item in rank_three_records if item["root_type"] == "A7+D7"),
    key=lambda item: (
        item["divisor_complexity"]["max_abs"],
        item["divisor_complexity"]["l1"],
        item["orbit"],
    ),
)

payload = {
    "schema": "elkies-k3.e6a1-rho19-genuine-q2-neighbors.v1",
    "status": "PASS_EXACT_COMPLETE_GENUINE_Q2_CENSUS_AND_18_NEF_MW3_FRAMES",
    "inputs": {relative(SOURCE): digest(SOURCE)},
    "old_fibration": {
        "fibre_profile": "2IV*+I4+4I1",
        "root_type": "2E6+A3",
        "root_rank": 15,
        "mw_rank": 2,
        "ns_determinant": 36,
    },
    "zero_section_obstruction": nominal_obstruction,
    "genuine_q2_shape": {
        "split_class": "D=2e+2f-w",
        "norm_equation": "w^2=8",
        "D_squared": 0,
        "D_dot_old_fibre": 2,
        "D_dot_old_zero": 0,
    },
    "dominant_enumeration": {
        "component_label_counts": list(enumeration["component_label_counts"]),
        "combined_label_count": enumeration["combined_label_count"],
        "compatible_label_count": enumeration["compatible_label_count"],
        "pairing_smith_invariants": [
            int(value) for value in enumeration["pairing_smith_invariants"]
        ],
        "rank_two_kernel_gram": integer_rows(enumeration["kernel_gram"]),
        "weyl_orbit_count": len(dominant_vectors),
        "nonprimitive_orbit_count": len(nonprimitive_records),
        "primitive_orbit_count": len(dominant_vectors) - len(nonprimitive_records),
    },
    "physical_chamber_degree_distribution": {
        str(key): value for key, value in sorted(degree_distribution.items())
    },
    "child_root_data_distribution": [
        {"root_data": list(key), "count": value}
        for key, value in sorted(raw_distribution.items())
    ],
    "all_orbit_records": all_records,
    "nonprimitive_orbit_records": nonprimitive_records,
    "nef_rank_three_neighbors": rank_three_records,
    "nef_rank_three_count": len(rank_three_records),
    "root_type_counts_among_nef_rank_three": dict(sorted(Counter(
        item["root_type"] for item in rank_three_records
    ).items())),
    "preferred_equation_compiler_target": {
        "selection_rule": "minimum (max_abs,l1,orbit) in the old NS basis",
        "orbit": best["orbit"],
        "root_type": best["root_type"],
        "divisor_in_ns_basis": best["divisor_in_ns_basis"],
        "divisor_complexity": best["divisor_complexity"],
        "mw_height_gram": best["mw_height_gram"],
        "reason": (
            "The class is exactly P0+P1+A3_2 in the old marking, making it "
            "the lowest-complexity resolved Riemann-Roch target."
        ),
    },
    "secondary_fibre_simple_compiler_target": {
        "selection_rule": "A7+D7, then minimum (max_abs,l1,orbit) in the old NS basis",
        "orbit": fibre_simple_best["orbit"],
        "root_type": fibre_simple_best["root_type"],
        "divisor_in_ns_basis": fibre_simple_best["divisor_in_ns_basis"],
        "divisor_complexity": fibre_simple_best["divisor_complexity"],
        "mw_height_gram": fibre_simple_best["mw_height_gram"],
        "reason": "Only two reducible fibres, at the cost of a less sparse old-marking divisor.",
    },
    "proof_boundary": {
        "proved": (
            "The norm-eight q=2 Weyl census is complete; all eighteen root-rank-14 "
            "classes pass the physical fibre chamber, the complete old-section "
            "gate, and the degree-two horizontal-root exclusion, hence are nef "
            "Jacobian fibrations with MW rank three."
        ),
        "open": (
            "No Riemann-Roch pencil, Weierstrass equation, transported child "
            "sections, or arithmetic specialization rank is claimed here."
        ),
    },
}

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
output_path = arguments.output.resolve()
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "E6A1Q2|orbits=119|nef_MW3=18|types=4|status=PASS_EXACT",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
