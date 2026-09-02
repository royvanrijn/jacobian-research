#!/usr/bin/env sage-python
"""Certify the complete genuine degree 2--4 rootless search for E6 rank four.

status: ACTIVE_COMPUTATION
claim: the determinant-78 E6 rank-four NS has no rootless primitive isotropic
       class in the complete zero-neutral old-degree 2, 3, or 4 shells
inputs: elkies-k3-e6-rank4-linear-chord-incidence-v1.json and the stored frame
outputs: elkies-k3-e6-rank4-rootless-low-degree-search-v1.json

This is an exact finite lattice census, not a global nonexistence theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import product
from pathlib import Path

from sage.all import (
    CartanMatrix, IntegralLattice, QQ, ZZ, floor, gcd, matrix, pari, vector,
    zero_matrix,
)


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GEN = ROOT / "artifacts/generated-results"
SOURCE = GEN / "elkies-k3-e6-rank4-linear-chord-incidence-v1.json"
FRAME_PATH = ROOT / "elkies-k3/data/lattice/e6_rank4_det78_frame.txt"
DEFAULT_OUTPUT = GEN / "elkies-k3-e6-rank4-rootless-low-degree-search-v1.json"

_engine_path = HERE / "exact_neighbor_engine.sage"
exec(compile(_engine_path.read_text(), str(_engine_path), "exec"), globals())


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def integer_rows(value):
    return [[int(entry) for entry in row] for row in matrix(ZZ, value).rows()]


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        pending = [min(unseen)]
        unseen.remove(pending[0])
        component = []
        while pending:
            index = pending.pop()
            component.append(index)
            adjacent = [other for other in unseen if cartan[index, other] != 0]
            for other in adjacent:
                unseen.remove(other)
                pending.append(other)
        result.append(tuple(sorted(component)))
    return tuple(result)


def component_label_vectors(cartan, norm_bound):
    inverse = cartan.inverse()
    bounds = [
        int(floor((QQ(norm_bound) / inverse[index, index]).sqrt()))
        for index in range(inverse.nrows())
    ]
    result = []
    for values in product(*[range(bound + 1) for bound in bounds]):
        label = vector(ZZ, values)
        norm = label * inverse * label
        if norm <= norm_bound:
            result.append((label, norm))
    return tuple(result)


def enumerate_dominant_norm_vectors(frame, simple, target_norm):
    """Enumerate the complete norm shell modulo the old fibre Weyl group."""
    cartan = simple * frame * simple.transpose()
    components = connected_components(cartan)
    component_labels = []
    for component in components:
        block = cartan.matrix_from_rows_and_columns(component, component)
        component_labels.append(
            (component, component_label_vectors(block, target_norm))
        )

    combined = []
    for choices in product(*[records for unused, records in component_labels]):
        total = sum((norm for unused, norm in choices), QQ(0))
        if total > target_norm:
            continue
        label = vector(ZZ, simple.nrows())
        for (component, unused), (part, unused_norm) in zip(
            component_labels, choices
        ):
            for index, value in zip(component, part):
                label[index] = value
        combined.append(label)

    pairing = simple * frame
    smith, left, right = pairing.smith_form()
    if smith != left * pairing * right:
        raise AssertionError("unexpected Smith convention")
    rank = pairing.rank()
    diagonal = tuple(abs(ZZ(smith[index, index])) for index in range(rank))
    kernel = right[:, rank:]
    kernel_gram = kernel.transpose() * frame * kernel
    lattice = IntegralLattice(kernel_gram)

    answers = []
    compatible = 0
    for label in combined:
        rhs = left * label
        if any(rhs[index] % diagonal[index] for index in range(rank)):
            continue
        compatible += 1
        coordinates = vector(ZZ, frame.nrows())
        for index in range(rank):
            coordinates[index] = rhs[index] // diagonal[index]
        particular = right * coordinates
        centre = -kernel_gram.inverse() * kernel.transpose() * frame * particular
        for close in lattice.enumerate_close_vectors(centre):
            candidate = particular + kernel * vector(ZZ, close)
            norm = candidate * frame * candidate
            if norm > target_norm:
                break
            if norm == target_norm:
                if pairing * candidate != label:
                    raise AssertionError("dominant-vector reconstruction failed")
                answers.append(tuple(entries(candidate)))

    answers = tuple(vector(ZZ, item) for item in sorted(set(answers)))
    return {
        "vectors": answers,
        "component_label_counts": tuple(
            len(records) for unused, records in component_labels
        ),
        "combined_label_count": len(combined),
        "compatible_label_count": compatible,
        "pairing_smith_invariants": diagonal,
        "kernel_gram": kernel_gram,
    }


def integral_ns():
    """Rebuild the source's O,F,2E6,A1,P,Q,R1,R2 marking."""
    ns = zero_matrix(ZZ, 19)
    ns[0, 0] = -2
    ns[0, 1] = ns[1, 0] = 1
    ns[2:8, 2:8] = -CartanMatrix(["E", 6])
    ns[8:14, 8:14] = -CartanMatrix(["E", 6])
    ns[14, 14] = -2
    for section_index in range(15, 19):
        ns[section_index, section_index] = -2
        ns[1, section_index] = ns[section_index, 1] = 1
    for section_index in (15, 16):
        for component_index in (2, 8):
            ns[section_index, component_index] = 1
            ns[component_index, section_index] = 1
    for section_index in (17, 18):
        for component_index in (2, 14):
            ns[section_index, component_index] = 1
            ns[component_index, section_index] = 1
    for left, right in ((15, 18), (16, 17)):
        ns[left, right] = ns[right, left] = 1
    return ns


def physical_walls(ns, old_fibre, old_zero):
    walls = [("O", old_zero)]
    highest_e6 = (1, 2, 2, 3, 2, 1)
    for name, start, size, highest in (
        ("E6a", 2, 6, highest_e6),
        ("E6b", 8, 6, highest_e6),
        ("A1", 14, 1, (1,)),
    ):
        components = []
        for index in range(size):
            component = vector(ZZ, 19)
            component[start + index] = 1
            components.append(component)
            walls.append((f"{name}_{index + 1}", component))
        theta = sum(
            (coefficient * component for coefficient, component in zip(
                highest, components
            )),
            vector(ZZ, 19),
        )
        walls.append((f"{name}_0", old_fibre - theta))
    for name, curve in walls:
        if curve * ns * curve != -2:
            raise ArithmeticError(f"physical wall {name} is not a root")
    return tuple(walls)


if not SOURCE.exists():
    raise FileNotFoundError(SOURCE)
source = json.loads(SOURCE.read_text())
if not source.get("status", "").startswith("PASS_EXACT"):
    raise ArithmeticError("rank-four source is not exact")

ns = integral_ns()
if integer_rows(ns) != source["neron_severi"]["integral_gram"]:
    raise ArithmeticError("source integral NS marking changed")
if ns.det() != 78 or ns.elementary_divisors() != [1] * 18 + [78]:
    raise ArithmeticError("determinant-78 integral NS gate failed")

old_fibre = vector(ZZ, [0, 1] + [0] * 17)
old_zero = vector(ZZ, [1] + [0] * 18)
constraints = matrix(ZZ, [ns * old_fibre, ns * (old_zero + old_fibre)])
frame_basis = constraints.right_kernel().basis_matrix()
transport = matrix(ZZ, [old_fibre, old_zero + old_fibre, *frame_basis.rows()])
if abs(transport.det()) != 1:
    raise ArithmeticError("old marked U does not split unimodularly")
split = transport * ns * transport.transpose()
frame = -split[2:, 2:]
if frame != load_gram(FRAME_PATH) or frame.det() != 78:
    raise ArithmeticError("stored positive frame changed")

simple, unused_positive, old_root_data = deterministic_simple_roots(frame)
if roots_and_data(frame)[2] != (13, 146, 18):
    raise ArithmeticError("old 2E6+A1 root system changed")
walls = physical_walls(ns, old_fibre, old_zero)

expected = {
    2: {
        "orbits": 280,
        "nonprimitive": 3,
        "degree_distribution": {0: 8, 1: 28, 2: 241},
    },
    3: {
        "orbits": 6242,
        "nonprimitive": 3,
        "degree_distribution": {0: 27, 1: 204, 2: 1626, 3: 4382},
    },
    4: {
        "orbits": 73601,
        "nonprimitive": 280,
        "degree_distribution": {0: 92, 1: 748, 2: 7992, 3: 19464, 4: 45025},
    },
}

layers = []
for degree in (2, 3, 4):
    norm = 2 * degree * degree
    enumeration = enumerate_dominant_norm_vectors(frame, simple, ZZ(norm))
    dominant = enumeration["vectors"]
    nonprimitive = 0
    degree_distribution = Counter()
    signed_root_count_distribution = Counter()
    rootless = []
    for orbit, horizontal in enumerate(dominant):
        divisor_split = vector(ZZ, [degree, degree] + [-entry for entry in horizontal])
        divisor = divisor_split * transport
        if divisor * ns * divisor != 0:
            raise ArithmeticError("norm shell failed the isotropic equation")
        if gcd(list(divisor)) != 1:
            nonprimitive += 1
            continue
        reduced, reduction = reduce_fixed_components(divisor, walls, ns)
        degree_distribution[int(reduced * ns * old_fibre)] += 1

        child = primitive_hyperbolic_split(ns, divisor)["child_frame"]
        # maxnum=1 limits only returned representatives.  PARI still returns
        # the exact total number of vectors of norm at most two in slot zero.
        signed_root_count = int(pari(child).qfminim(2, 1)[0])
        signed_root_count_distribution[signed_root_count] += 1
        if signed_root_count == 0:
            rootless.append({
                "orbit": orbit,
                "horizontal": entries(horizontal),
                "divisor_in_ns_basis": entries(divisor),
                "physical_reduction": [
                    [name, int(pairing)] for name, pairing in reduction
                ],
            })

    check = expected[degree]
    if len(dominant) != check["orbits"]:
        raise ArithmeticError(f"degree-{degree} Weyl orbit count changed")
    if nonprimitive != check["nonprimitive"]:
        raise ArithmeticError(f"degree-{degree} primitive census changed")
    if dict(degree_distribution) != check["degree_distribution"]:
        raise ArithmeticError(f"degree-{degree} physical degree census changed")
    if rootless:
        raise ArithmeticError(
            f"degree-{degree} rootless hit requires a promoted certificate"
        )
    layers.append({
        "old_degree": degree,
        "zero_neutral_shape": f"D={degree}e+{degree}f-w",
        "required_frame_norm": norm,
        "component_label_counts": list(enumeration["component_label_counts"]),
        "combined_label_count": enumeration["combined_label_count"],
        "compatible_label_count": enumeration["compatible_label_count"],
        "pairing_smith_invariants": list(
            map(int, enumeration["pairing_smith_invariants"])
        ),
        "rank_four_kernel_gram": integer_rows(enumeration["kernel_gram"]),
        "weyl_orbit_count": len(dominant),
        "nonprimitive_orbit_count": nonprimitive,
        "primitive_orbit_count": len(dominant) - nonprimitive,
        "physical_reduced_old_degree_distribution": {
            str(key): value for key, value in sorted(degree_distribution.items())
        },
        "signed_root_count_distribution": {
            str(key): value
            for key, value in sorted(signed_root_count_distribution.items())
        },
        "minimum_signed_root_count": min(signed_root_count_distribution),
        "rootless_count": 0,
    })
    print(
        f"E6R4ROOTLESS|degree={degree}|orbits={len(dominant)}|"
        f"primitive={len(dominant)-nonprimitive}|rootless=0",
        flush=True,
    )

payload = {
    "schema": "elkies-k3.e6-rank4-rootless-low-degree-search.v1",
    "status": "PASS_EXACT_COMPLETE_ZERO_NEUTRAL_DEGREE_2_TO_4_NO_ROOTLESS",
    "inputs": {
        relative(SOURCE): digest(SOURCE),
        relative(FRAME_PATH): digest(FRAME_PATH),
    },
    "source": {
        "ns_determinant": 78,
        "ns_smith_invariants": [1] * 18 + [78],
        "old_root_type": "2E6+A1",
        "old_root_data": [13, 146, 18],
        "old_mw_rank": 4,
        "positive_frame_gram": integer_rows(frame),
    },
    "search": {
        "method": (
            "complete Weyl-dominant Dynkin-label enumeration followed by an "
            "exact rank-four closest-vector reconstruction, primitivity, physical "
            "fixed-component reduction, primitive U splitting, and exact PARI "
            "norm-two enumeration in every child frame"
        ),
        "layers": layers,
        "total_dominant_orbits": sum(layer["weyl_orbit_count"] for layer in layers),
        "total_primitive_orbits": sum(layer["primitive_orbit_count"] for layer in layers),
        "rootless_hits": 0,
    },
    "proof_boundary": {
        "proved": (
            "Every zero-neutral primitive isotropic class D=q(e+f)-w with "
            "q in {2,3,4} has a norm-two root in D-perp/<D>; hence none gives "
            "a rootless MW17 fibration.  The shell enumeration is complete modulo "
            "the full Weyl group of the old 2E6+A1 fibre roots."
        ),
        "open": (
            "No claim is made for degree at least five, a multi-step route, or "
            "global nonexistence of a rootless U embedding in this NS lattice."
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
    "E6R4ROOTLESS|degrees=2,3,4|rootless=0|status=PASS_EXACT_BOUNDED",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
