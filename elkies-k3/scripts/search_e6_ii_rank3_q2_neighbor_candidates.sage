#!/usr/bin/env sage-python
"""Enumerate the first q=2 lattice-neighbor shell of the E6+II MW3 K3.

status: ACTIVE_PROOF
claim: complete norm-eight Weyl census and six nef A6+D7/MW4 frames
inputs: elkies-k3-e6-ii-rank3-quadratic-base-change-v1.json
outputs: elkies-k3-e6-ii-rank3-q2-neighbor-candidates-v1.json

The source fibration has frame roots 2E6+A2 and MW rank three.  In the
split marking NS=U+M(-1), a zero-neutral degree-two class has the form
D=2e+2f-w with w.M.w=8.  Dominant Dynkin labels and an exact closest-vector
enumeration in the rank-three root-orthogonal kernel give one representative
of every Weyl orbit.  The script then splits the child U and records its
complete norm-two root system.

This is a complete lattice/physical-wall census and a global nefness proof for
the six minimum-root classes.  It does not claim equations for the child
fibrations or that their four geometric sections descend to QQ(r).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import product
from pathlib import Path

from sage.all import CartanMatrix, IntegralLattice, QQ, ZZ, floor, gcd, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GEN = ROOT / "artifacts/generated-results"
SOURCE = GEN / "elkies-k3-e6-ii-rank3-quadratic-base-change-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-e6-ii-rank3-q2-neighbor-candidates-v1.json"

_engine_path = HERE / "exact_neighbor_engine.sage"
exec(compile(_engine_path.read_text(), str(_engine_path), "exec"), globals())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def integer_rows(value):
    return [[int(entry) for entry in row] for row in matrix(ZZ, value).rows()]


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
    frame = matrix(ZZ, frame)
    simple = matrix(ZZ, simple)
    cartan = simple * frame * simple.transpose()
    components = connected_components(cartan)
    component_labels = []
    for component in components:
        block = cartan.matrix_from_rows_and_columns(component, component)
        component_labels.append((component, component_label_vectors(block, target_norm)))

    combined = []
    for choices in product(*[records for _, records in component_labels]):
        total = sum((norm for _, norm in choices), QQ(0))
        if total > target_norm:
            continue
        label = vector(ZZ, simple.nrows())
        for (component, _), (part, _) in zip(component_labels, choices):
            for index, value in zip(component, part):
                label[index] = value
        combined.append(label)

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
    compatible = 0
    for label in combined:
        rhs = left * label
        if any(rhs[i] % diagonal[i] for i in range(rank)):
            continue
        compatible += 1
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
                answers.append(tuple(entries(candidate)))

    answers = tuple(vector(ZZ, item) for item in sorted(set(answers)))
    return {
        "vectors": answers,
        "cartan": cartan,
        "component_label_counts": tuple(len(records) for _, records in component_labels),
        "combined_label_count": len(combined),
        "compatible_label_count": compatible,
        "pairing_smith_invariants": diagonal,
        "kernel_gram": kernel_gram,
    }


def ade_component_name(cartan):
    cartan = matrix(ZZ, cartan)
    rank = cartan.nrows()
    determinant = abs(ZZ(cartan.det()))
    roots = ZZ(pari(cartan).qfminim(2)[0])
    if roots == rank * (rank + 1) and determinant == rank + 1:
        return f"A{rank}"
    if rank >= 4 and roots == 2 * rank * (rank - 1) and determinant == 4:
        return f"D{rank}"
    exceptional = {
        (6, 72, 3): "E6",
        (7, 126, 2): "E7",
        (8, 240, 1): "E8",
    }
    key = (rank, int(roots), int(determinant))
    if key not in exceptional:
        raise ArithmeticError(f"unrecognized ADE component {key}")
    return exceptional[key]


def root_type(frame):
    simple, _, rows = deterministic_simple_roots(frame)
    if not simple.nrows():
        return "rootless"
    cartan = matrix(ZZ, rows)
    names = []
    for component in connected_components(cartan):
        block = cartan.matrix_from_rows_and_columns(component, component)
        names.append(ade_component_name(block))
    counts = Counter(names)
    return "+".join(
        (name if count == 1 else f"{count}{name}")
        for name, count in sorted(counts.items())
    )


def correction_vector(a, b, s):
    """Shioda correction for a*P+b*Q+s*S in the old fibration."""
    result = vector(QQ, 19)
    residues = ((a + b + s) % 3, (a + b - s) % 3, s % 3)
    blocks = ((2, 6, "E"), (8, 6, "E"), (14, 2, "A"))
    for residue, (start, rank, kind) in zip(residues, blocks):
        if not residue:
            continue
        node = 0 if residue == 1 else rank - 1
        result[start:start + rank] = CartanMatrix([kind, rank]).inverse().row(node)
    return result


def section_class(a, b, s, phis):
    correction = correction_vector(a, b, s)
    coefficients = vector(ZZ, [a, b, s])
    height_value = coefficients * old_height * coefficients
    local = -(correction * ns * correction)
    intersection_zero = (height_value + local - 4) / 2
    if intersection_zero not in ZZ:
        raise ArithmeticError("section/zero intersection is not integral")
    section = (
        sum((coefficient * phi for coefficient, phi in zip(coefficients, phis)),
            vector(QQ, 19))
        + old_zero + (2 + intersection_zero) * old_fibre - correction
    )
    if not all(value in ZZ for value in section):
        raise ArithmeticError("reconstructed section class is not integral")
    section = vector(ZZ, section)
    if section * ns * section != -2 or section * ns * old_fibre != 1:
        raise ArithmeticError("reconstructed section class changed")
    return section


if not SOURCE.exists():
    raise FileNotFoundError(SOURCE)
source = json.loads(SOURCE.read_text())
if not source.get("status", "").startswith("PASS_EXACT"):
    raise ArithmeticError("source certificate is not exact")

ns = matrix(ZZ, source["k3"]["neron_severi_gram"])
if ns.det() != 24:
    raise ArithmeticError("source NS changed")
old_fibre = vector(ZZ, [0, 1] + [0] * 17)
old_zero = vector(ZZ, [1] + [0] * 18)
constraints = matrix(ZZ, [ns * old_fibre, ns * (old_zero + old_fibre)])
frame_basis = constraints.right_kernel().basis_matrix()
transport = matrix(ZZ, [old_fibre, old_zero + old_fibre, *frame_basis.rows()])
if abs(transport.det()) != 1:
    raise ArithmeticError("old U split is not unimodular")
frame = -(transport * ns * transport.transpose())[2:, 2:]
if frame.det() != 24:
    raise ArithmeticError("positive frame determinant changed")

simple, _, old_cartan_rows = deterministic_simple_roots(frame)
if tuple(roots_and_data(frame)[2]) != (14, 150, 27):
    raise ArithmeticError("old 2E6+A2 root data changed")
if root_type(frame) != "A2+2E6":
    raise ArithmeticError("old root type changed")

# The exact source certificate proves that P,Q,S are a saturated torsion-free
# MW basis.  Recover their Shioda vectors for the complete old-section gate.
section_basis = []
for index, coefficients in zip((16, 17, 18), ((1, 0, 0), (0, 1, 0), (0, 0, 1))):
    section = vector(ZZ, 19)
    section[index] = 1
    section_basis.append(section - old_zero - 2 * old_fibre + correction_vector(*coefficients))
section_basis = tuple(section_basis)
old_height = matrix(QQ, source["k3"]["height_gram"])
shioda_gram = matrix(QQ, 3, 3, lambda i, j: section_basis[i] * ns * section_basis[j])
if shioda_gram != -old_height:
    raise ArithmeticError("old Shioda basis changed")

enumeration = enumerate_dominant_norm_vectors(frame, simple, ZZ(8))
dominant_vectors = enumeration["vectors"]

split_fibre = vector(ZZ, [1, 0] + [0] * 17)
split_zero = vector(ZZ, [-1, 1] + [0] * 17)
split_walls, _ = component_walls(frame, split_fibre, include_zero=split_zero)
walls = tuple((name, curve * transport) for name, curve in split_walls)

degree_distribution = Counter()
root_distribution = Counter()
root_type_distribution = Counter()
records = []
nonprimitive = []
for orbit, dominant in enumerate(dominant_vectors):
    split_divisor = vector(ZZ, [2, 2] + [-entry for entry in dominant])
    divisor = split_divisor * transport
    if divisor * ns * divisor != 0:
        raise ArithmeticError("norm-eight vector did not give an isotropic class")
    divisor_gcd = gcd(list(divisor))
    if divisor_gcd != 1:
        nonprimitive.append({
            "orbit": orbit,
            "dominant_w": entries(dominant),
            "divisor_gcd": int(divisor_gcd),
        })
        continue
    if divisor * ns * old_fibre != 2 or divisor * ns * old_zero != 0:
        raise ArithmeticError("q=2 intersection profile changed")

    reduced, reduction = reduce_fixed_components(divisor, walls, ns)
    degree = int(reduced * ns * old_fibre)
    degree_distribution[degree] += 1
    split = primitive_hyperbolic_split(ns, divisor)
    root_data = tuple(int(value) for value in roots_and_data(split["child_frame"])[2])
    child_type = root_type(split["child_frame"])
    root_distribution[root_data] += 1
    root_type_distribution[child_type] += 1
    records.append({
        "orbit": orbit,
        "dominant_w": entries(dominant),
        "divisor_in_ns_basis": entries(divisor),
        "divisor_complexity": {
            "max_abs": int(max(abs(value) for value in divisor)),
            "l1": int(sum(abs(value) for value in divisor)),
        },
        "physical_reduction": {
            "sequence": [[name, int(pairing)] for name, pairing in reduction],
            "old_fibre_degree": degree,
        },
        "child_root_data": list(root_data),
        "child_root_type": child_type,
        "geometric_mw_rank_if_nef": 17 - root_data[0],
        "child_frame_digest": matrix_digest(split["child_frame"]),
    })

minimum_root_rank = min(record["child_root_data"][0] for record in records)
best = [record for record in records if record["child_root_data"][0] == minimum_root_rank]
rootless = [record for record in records if record["child_root_data"][0] == 0]
if rootless:
    raise ArithmeticError("rootless q=2 candidate contradicts determinant obstruction")

# Complete nefness for the minimum-root candidates.  Once old components are
# nonnegative, a negative horizontal root has old degree one or two.  Degree
# one roots are sections and lie in the strict height ellipsoid below.  The
# standard primitive-fibre argument excludes degree-two roots.
nef_best = []
for record in best:
    if record["physical_reduction"]["sequence"] or record["physical_reduction"]["old_fibre_degree"] != 2:
        raise ArithmeticError("minimum-root candidate failed the physical chamber")
    divisor = vector(ZZ, record["divisor_in_ns_basis"])
    pairings = vector(QQ, [divisor * ns * phi for phi in section_basis])
    delta = -pairings * old_height.inverse()
    if not all(value in ZZ for value in delta):
        raise ArithmeticError("section-search centre is not integral")
    delta = vector(ZZ, delta)
    tested = []
    centres = [int(floor(QQ(value) / 2)) for value in delta]
    for a in range(centres[0] - 4, centres[0] + 5):
        for b in range(centres[1] - 4, centres[1] + 5):
            for s in range(centres[2] - 4, centres[2] + 5):
                coefficients = vector(ZZ, [a, b, s])
                offset = vector(QQ, coefficients) - vector(QQ, delta) / 2
                if offset * old_height * offset >= 2:
                    continue
                section = section_class(a, b, s, section_basis)
                tested.append({
                    "mw_coordinates": [a, b, s],
                    "pairing": int(divisor * ns * section),
                })
    if not tested or min(item["pairing"] for item in tested) < 0:
        raise ArithmeticError("minimum-root candidate has a negative old section")
    split = primitive_hyperbolic_split(ns, divisor)
    minimized = minimize_child_frame(split["child_frame"])
    if not minimized["root_lattice_primitive"] or minimized["mw_height"] is None:
        raise ArithmeticError("minimum-root candidate has unresolved root glue")
    height = minimized["mw_height"]
    if height.nrows() != 4 or height.det() != QQ(6) / 7:
        raise ArithmeticError("minimum-root MW4 regulator changed")
    record["complete_nefness"] = {
        "physical_components": "nonnegative",
        "section_centre_delta": entries(delta),
        "strict_ellipsoid": "(n-delta/2)^T*H*(n-delta/2)<2",
        "tested_sections": tested,
        "minimum_section_pairing": min(item["pairing"] for item in tested),
        "degree_two_root_exclusion": (
            "If C is a negative old-degree-two root, D-C=lF; then "
            "C^2=-4l=-2 forces l=1/2, contradicting primitive F."
        ),
        "nef": True,
    }
    record["child_mw_height_gram"] = [
        [str(value) for value in row] for row in height.rows()
    ]
    record["child_mw_regulator"] = str(height.det())
    record["child_torsion"] = "trivial"
    nef_best.append(record)

preferred = min(
    nef_best,
    key=lambda item: (
        item["divisor_complexity"]["max_abs"],
        item["divisor_complexity"]["l1"],
        item["orbit"],
    ),
)

payload = {
    "schema": "elkies-k3.e6-ii-rank3-q2-neighbor-candidates.v1",
    "status": "PASS_EXACT_COMPLETE_Q2_CENSUS_AND_6_NEF_MW4_FRAMES",
    "inputs": {str(SOURCE.relative_to(ROOT)): digest(SOURCE)},
    "source_fibration": {
        "root_type": "2E6+A2",
        "root_rank": 14,
        "mw_rank": 3,
        "ns_determinant": 24,
    },
    "genuine_q2_shape": {
        "split_class": "D=2e+2f-w",
        "norm_equation": "w^2=8",
        "D_dot_old_fibre": 2,
        "D_dot_old_zero": 0,
    },
    "dominant_enumeration": {
        "component_label_counts": [int(value) for value in enumeration["component_label_counts"]],
        "combined_label_count": enumeration["combined_label_count"],
        "compatible_label_count": enumeration["compatible_label_count"],
        "pairing_smith_invariants": [
            int(value) for value in enumeration["pairing_smith_invariants"]
        ],
        "root_orthogonal_kernel_gram": integer_rows(enumeration["kernel_gram"]),
        "weyl_orbit_count": len(dominant_vectors),
        "primitive_orbit_count": len(records),
        "nonprimitive_orbit_count": len(nonprimitive),
    },
    "physical_wall_degree_distribution": {
        str(key): value for key, value in sorted(degree_distribution.items())
    },
    "child_root_data_distribution": [
        {"root_data": list(key), "count": value}
        for key, value in sorted(root_distribution.items())
    ],
    "child_root_type_distribution": dict(sorted(root_type_distribution.items())),
    "minimum_child_root_rank": minimum_root_rank,
    "maximum_geometric_mw_rank_if_nef": 17 - minimum_root_rank,
    "minimum_root_rank_records": best,
    "nef_mw4_neighbor_count": len(nef_best),
    "nef_mw4_neighbors": nef_best,
    "preferred_equation_compiler_target": {
        "selection_rule": "minimum (max_abs,l1,orbit) in the old NS basis",
        "orbit": preferred["orbit"],
        "divisor_in_ns_basis": preferred["divisor_in_ns_basis"],
        "divisor_complexity": preferred["divisor_complexity"],
        "root_type": preferred["child_root_type"],
        "mw_height_gram": preferred["child_mw_height_gram"],
    },
    "rootless_candidate_count": len(rootless),
    "nonprimitive_records": nonprimitive,
    "all_primitive_orbit_records": records,
    "proof_boundary": {
        "proved": (
            "The norm-eight q=2 Weyl-orbit census, finite physical-wall reductions, "
            "primitive U-splits, and child norm-two root systems are exact.  The six "
            "minimum-root A6+D7 classes pass the complete old-section ellipsoid and "
            "degree-two-root gates, hence are nef Jacobian fibrations of MW rank four."
        ),
        "not_proved": (
            "Global nefness is not asserted for the other 259 primitive classes.  No "
            "equation-level child fibration or arithmetic descent of its sections is asserted."
        ),
        "rootless_obstruction": (
            "The independent Blichfeldt-Hermite argument in the source certificate "
            "rules out every rootless rank-17 frame in determinant 24, not only q=2."
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
    f"E6IIR3Q2|orbits={len(dominant_vectors)}|primitive={len(records)}|"
    f"min_root_rank={minimum_root_rank}|nef_MW4={len(nef_best)}|"
    "rootless=0|status=PASS_EXACT_Q2_AND_NEF_MW4",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
