#!/usr/bin/env sage-python
"""Exhaust a Golay-720 three-support pole-zero MW2 marking modulo p.

The fibre artifact supplies a semistable rank-15 three-A-component stratum.
The source and Smith-quotient pole artifacts supply a complete pole-zero MW
basis.  This script reconstructs the deterministic simple-root components,
derives both local component depths and their cross-correction exactly, then
enumerates every polynomial section and every candidate basis pair on every
stored fibre model.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DEFAULT_FIBRES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-3a5-source-ansatz-mod5-v1.json"
)
DEFAULT_SOURCES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json"
)
DEFAULT_POLES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-source-poles-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-3a5-pole0-pairs-mod5-v1.json"
)

_engine_path = HERE / "exact_neighbor_engine.sage"
exec(compile(_engine_path.read_text(), str(_engine_path), "exec"), globals())


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def component_ranks(root_type):
    ranks = []
    for term in root_type.split("+"):
        match = re.fullmatch(r"(?:(\d+))?A(\d+)", term)
        if match is None:
            raise ValueError("source is not a sum of A components")
        ranks.extend([int(match.group(2))] * int(match.group(1) or 1))
    return ranks


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    answer = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        stack = [first]
        component = []
        while stack:
            left = stack.pop()
            component.append(left)
            adjacent = [
                right
                for right in sorted(unseen)
                if cartan[left, right] != 0
            ]
            for right in adjacent:
                unseen.remove(right)
                stack.append(right)
        answer.append(sorted(component))
    return answer


def ordered_components(cartan, ranks):
    available = connected_components(cartan)
    answer = []
    for rank in ranks:
        match = next(
            (component for component in available if len(component) == rank), None
        )
        if match is None:
            raise ArithmeticError("Cartan components do not match the root type")
        answer.append(match)
        available.remove(match)
    if available:
        raise ArithmeticError("unused Cartan component")
    return answer


def depth_from_correction(root_rank, correction):
    fibre_order = root_rank + 1
    matches = [
        depth
        for depth in range(fibre_order // 2 + 1)
        if QQ(depth * (fibre_order - depth)) / fibre_order == correction
    ]
    if len(matches) != 1:
        raise ArithmeticError(
            f"cannot identify component depth from correction {correction}"
        )
    return matches[0]


def formal_center(A, B, point, precision):
    field = A.base_ring()
    base = A.parent()
    t = base.gen()
    shifted_A = base(A(t + point))
    shifted_B = base(B(t + point))
    node = -field(3) * shifted_B[0] / (field(2) * shifted_A[0])
    series_ring = PowerSeriesRing(field, "s", default_prec=precision + 3)
    center = series_ring(node)
    series_A = series_ring(shifted_A)
    for unused in range(8):
        center = (center + (-series_A / 3) / center) / 2
    if (center**2 + series_A / 3).valuation() < precision + 1:
        raise ArithmeticError("finite formal center did not converge")
    return series_ring, center


def reversed_local(poly, weight, series_ring):
    return series_ring(
        sum(
            poly[index] * series_ring.gen() ** (weight - index)
            for index in range(poly.degree() + 1)
        )
    )


def polynomial_roots(right):
    if not right.is_square():
        return []
    positive = right.sqrt()
    return [positive] if not positive else [positive, -positive]


def serialize(poly):
    return [int(value) for value in poly]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
parser.add_argument("--poles", type=Path, default=DEFAULT_POLES)
parser.add_argument("--quadratic-twist", type=int, default=1)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

fibres_path = arguments.fibres.resolve()
sources_path = arguments.sources.resolve()
poles_path = arguments.poles.resolve()
output_path = arguments.output.resolve()
fibres = json.loads(fibres_path.read_text())
sources = json.loads(sources_path.read_text())
poles = json.loads(poles_path.read_text())
if fibres.get("schema") != (
    "elkies-k3.lattice-foundry-three-support-semistable-source-ansatz-modp.v1"
):
    raise ValueError("unexpected three-support fibre schema")
if not fibres["scan"]["exhausted"]:
    raise ValueError("marked scan requires an exhaustive fibre artifact")
source_id = fibres["source"]["source_id"]
source_entry = next(row for row in sources["sources"] if row["source_id"] == source_id)
source = source_entry["source"]
if source["gram_sha256"] != fibres["source"]["source_gram_sha256"]:
    raise ValueError("fibre/source Gram mismatch")
if source["mw_rank_for_rho_19"] != 2 or source["support_count"] != 3:
    raise ValueError("scanner requires a three-support MW2 source")
ranks = component_ranks(source["root_type"])
if len(ranks) != 3 or sum(ranks) != 15:
    raise ValueError("scanner requires three A components of total rank 15")
expected_supports = [
    f"0:I{ranks[0]+1}", f"1:I{ranks[1]+1}", f"infinity:I{ranks[2]+1}"
]
if fibres["ansatz"]["normalized_reducible_supports"] != expected_supports:
    raise ValueError("fibre support order disagrees with the source root type")

pole_entry = next(
    row
    for row in poles["audits"]
    if row["source_id"] == source_id
    and row["source_gram_sha256"] == source["gram_sha256"]
)
audit = pole_entry["audit"]
if audit["basis_sorted_pole_profile"] != [0, 0] or audit["torsion_order"] != 1:
    raise ValueError("source does not have a primitive pole-zero MW basis")

frame = matrix(QQ, source["gram"])
simple, unused_positive, cartan_rows = deterministic_simple_roots(frame)
cartan = matrix(QQ, cartan_rows)
components = ordered_components(cartan, ranks)
basis_pairings = [vector(QQ, row["simple_root_pairings"]) for row in audit["basis"]]
depth_profiles = [[], []]
component_cross = QQ(0)
for rank, component in zip(ranks, components):
    block = cartan.matrix_from_rows_and_columns(component, component)
    inverse = block.inverse()
    local = [vector(QQ, [pairing[index] for index in component]) for pairing in basis_pairings]
    self_corrections = [value * inverse * value for value in local]
    for basis_index in range(2):
        depth_profiles[basis_index].append(
            depth_from_correction(rank, self_corrections[basis_index])
        )
    component_cross += local[0] * inverse * local[1]

height = matrix(QQ, audit["height_gram"])
basis_coordinates = matrix(
    QQ, [section["free_mw_coordinates"] for section in audit["basis"]]
)
basis_height = basis_coordinates * height * basis_coordinates.transpose()
for index in range(2):
    self_correction = sum(
        QQ(depth * (rank + 1 - depth)) / (rank + 1)
        for rank, depth in zip(ranks, depth_profiles[index])
    )
    if QQ(4) - self_correction != basis_height[index, index]:
        raise ArithmeticError("component depths do not recover the section height")
required_intersection = QQ(2) - component_cross - basis_height[0, 1]
if required_intersection.denominator() != 1 or required_intersection < 0:
    raise ArithmeticError("physical basis requires an invalid intersection degree")
required_intersection = int(required_intersection)

prime = int(fibres["prime"])
field = GF(prime)
twist = field(arguments.quadratic_twist)
if not twist:
    raise ValueError("quadratic twist must be nonzero")
ring = PolynomialRing(field, "t")
t = ring.gen()
precision = max(max(profile) for profile in depth_profiles) + 2

models = []
x_candidates = 0
sections_by_basis = [0, 0]
marked_pairs = 0
pairs_meeting_singular_fibres = 0
for example_index, example in enumerate(fibres["examples"]):
    A = twist**2 * ring(example["A_coefficients_low_to_high"])
    B = twist**3 * ring(example["B_coefficients_low_to_high"])
    discriminant_core = 4 * A**3 + 27 * B**2
    zero_ring, zero_center = formal_center(A, B, field.zero(), precision)
    one_ring, one_center = formal_center(A, B, field.one(), precision)
    infinity_ring = PowerSeriesRing(field, "u", default_prec=precision + 3)
    infinity_A = reversed_local(A, 8, infinity_ring)
    infinity_B = reversed_local(B, 12, infinity_ring)
    infinity_center = infinity_ring(
        -field(3) * infinity_B[0] / (field(2) * infinity_A[0])
    )
    for unused in range(8):
        infinity_center = (
            infinity_center + (-infinity_A / 3) / infinity_center
        ) / 2
    if (infinity_center**2 + infinity_A / 3).valuation() < precision + 1:
        raise ArithmeticError("infinity formal center did not converge")

    candidates = [[], []]
    for coefficients in itertools.product(field, repeat=5):
        X = ring(list(coefficients))
        x_candidates += 1
        for Y in polynomial_roots(X**3 + A * X + B):
            depths = [
                int(
                    min(
                        (zero_ring(X(t)) - zero_center).valuation(),
                        zero_ring(Y(t)).valuation(),
                    )
                ),
                int(
                    min(
                        (one_ring(X(t + 1)) - one_center).valuation(),
                        one_ring(Y(t + 1)).valuation(),
                    )
                ),
                int(
                    min(
                        (
                            reversed_local(X, 4, infinity_ring)
                            - infinity_center
                        ).valuation(),
                        reversed_local(Y, 6, infinity_ring).valuation(),
                    )
                ),
            ]
            for basis_index, profile in enumerate(depth_profiles):
                if depths == profile:
                    candidates[basis_index].append({"X": X, "Y": Y})
                    sections_by_basis[basis_index] += 1
    if not candidates[0] or not candidates[1]:
        continue
    pairs = []
    for left_index, left in enumerate(candidates[0]):
        for right_index, right in enumerate(candidates[1]):
            if left["X"] == right["X"] and left["Y"] == right["Y"]:
                continue
            common = (left["X"] - right["X"]).gcd(left["Y"] - right["Y"])
            if common.gcd(discriminant_core).degree() != 0:
                pairs_meeting_singular_fibres += 1
                continue
            intersection = int(common.degree())
            if intersection != required_intersection:
                continue
            marked_pairs += 1
            pairs.append(
                {
                    "left_section_index": left_index,
                    "right_section_index": right_index,
                    "intersection_on_smooth_fibres": intersection,
                    "component_cross_correction": str(component_cross),
                    "shioda_height_pairing": str(basis_height[0, 1]),
                    "height_determinant": str(basis_height.det()),
                    "mw_index_from_determinant": 1,
                }
            )
    models.append(
        {
            "example_index": example_index,
            "basis_section_candidates": [
                [
                    {
                        "X_coefficients_low_to_high": serialize(row["X"]),
                        "Y_coefficients_low_to_high": serialize(row["Y"]),
                        "component_depths": depth_profiles[basis_index],
                    }
                    for row in candidates[basis_index]
                ]
                for basis_index in range(2)
            ],
            "marked_mw2_pairs": pairs,
        }
    )

output = {
    "schema": "elkies-k3.golay-det720-three-support-pole0-pairs-modp.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH_MARKED_MW2_PAIRS"
        if marked_pairs
        else "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_EMPTY_MARKED_MW2_PAIR_LOCUS"
    ),
    "prime": prime,
    "quadratic_twist": int(twist),
    "quadratic_twist_square_class": "square" if twist.is_square() else "nonsquare",
    "inputs": {
        relative(fibres_path): digest(fibres_path),
        relative(sources_path): digest(sources_path),
        relative(poles_path): digest(poles_path),
    },
    "source": {
        "source_id": source_id,
        "source_gram_sha256": source["gram_sha256"],
        "root_type": source["root_type"],
        "mw_rank": 2,
        "basis_pole_profile": [0, 0],
        "basis_component_depth_profiles": depth_profiles,
        "basis_height_gram": [
            [str(value) for value in row] for row in basis_height.rows()
        ],
        "component_cross_correction": str(component_cross),
        "required_smooth_pair_intersection": required_intersection,
    },
    "scope": {
        "fibre_census_exhaustive": True,
        "stored_fibre_models": len(fibres["examples"]),
        "polynomial_X_candidates_per_model": prime**5,
        "all_polynomial_Y_square_roots_retained": True,
        "all_component-depth-matched_pairs_tested": True,
    },
    "accounting": {
        "models_with_both_section_types": len(models),
        "X_candidates_scanned": x_candidates,
        "basis_section_candidates": sections_by_basis,
        "marked_mw2_pairs": marked_pairs,
        "pairs_meeting_singular_fibres": pairs_meeting_singular_fibres,
    },
    "models": models,
    "proof_boundary": {
        "proved": (
            "On every model in the exhaustive fibre chart, every polynomial "
            "X and every polynomial Y square root is checked. Component profiles "
            "and cross-corrections are recomputed from the exact source Cartan "
            "data, and a retained determinant-matched pair has MW index one."
        ),
        "not_proved": (
            "A finite-field pair does not construct a characteristic-zero family, "
            "Q-rational marking, parameterization, target multisection spectrum, "
            "or neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_golay_det720_three_support_pole0_pairs_modp.sage "
        f"--fibres {relative(fibres_path)} --output {relative(output_path)}"
        + (f" --quadratic-twist {int(twist)}" if twist != 1 else "")
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("Golay-720 three-support pair artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "GOLAY720THREEPAIR|"
    f"source={source_id}|p={prime}|twist={int(twist)}|"
    f"profiles={depth_profiles}|sections={sections_by_basis}|pairs={marked_pairs}|"
    f"status={'PASS' if marked_pairs else 'EMPTY'}",
    flush=True,
)
