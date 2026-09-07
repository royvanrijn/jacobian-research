#!/usr/bin/env sage-python
"""Exhaust the minimum NS0031 ``A1+2A7/MW2`` marking on fibre models.

The exact source ``NS0031-S001`` has displayed height Gram
``[[2,1],[1,41/8]]``.  Its minimum basis is the tail pair ``(1,0)`` and
``(1,-1)`` with pole profile ``(0,1)``.  In the pinned root-component order
their component depths at ``I2,I8,I8`` are respectively ``(1,0,2)`` and
``(0,1,0)``.  Their component cross-correction is zero, so the displayed
height pairing one requires geometric intersection degree two.

The three-support fibre stratum depends only on the normalized fibre orders,
so the script may reuse any exhaustive ``I2+I8+I8`` fibre artifact.  It scans
the complete pole-zero numerator chart on every model, then scans the pole-one
chart only on surviving models, in the requested quadratic-twist class.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIBRES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-source-ansatz-mod5-v1.json"
)
DEFAULT_SOURCE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json"
)
DEFAULT_POLES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-rank2-section-basis-poles-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-a1-2a7-marking-mod5-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    for unused in range(7):
        center = (center + (-series_A / 3) / center) / 2
    if (center**2 + series_A / 3).valuation() < precision + 2:
        raise ArithmeticError("formal center did not converge")
    return series_ring, center, node


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


def connected_components(gram):
    unseen = set(range(gram.nrows()))
    components = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        stack = [first]
        component = []
        while stack:
            left = stack.pop()
            component.append(left)
            adjacent = [right for right in sorted(unseen) if gram[left, right]]
            for right in adjacent:
                unseen.remove(right)
                stack.append(right)
        components.append(sorted(component))
    return sorted(components, key=lambda component: component[0])


def serialize_poly(poly):
    return [int(value) for value in poly]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument("--poles", type=Path, default=DEFAULT_POLES)
parser.add_argument("--quadratic-twist", type=int, default=1)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

fibres_path = arguments.fibres.resolve()
source_path = arguments.source.resolve()
poles_path = arguments.poles.resolve()
output_path = arguments.output.resolve()
fibres = json.loads(fibres_path.read_text())
source_payload = json.loads(source_path.read_text())
poles_payload = json.loads(poles_path.read_text())
if (
    fibres.get("schema")
    != "elkies-k3.lattice-foundry-three-support-semistable-source-ansatz-modp.v1"
):
    raise ValueError("unexpected three-support fibre schema")
if fibres["ansatz"]["normalized_reducible_supports"] != [
    "0:I2",
    "1:I8",
    "infinity:I8",
]:
    raise ValueError("NS0031 scan requires the I2+I8+I8 fibre stratum")
if not fibres["scan"]["exhausted"]:
    raise ValueError("NS0031 scan requires an exhaustive fibre census")

source_entry = next(
    row
    for row in source_payload["sources"]
    if row["ns_id"] == "NS0031"
    and row["source_id"] == "NS0031-S001"
    and row["source"]["root_type"] == "A1+2A7"
)
source = source_entry["source"]
if source["mw_height_gram"] != [["2", "1"], ["1", "41/8"]]:
    raise ValueError("unexpected NS0031 height Gram")
pole_entry = next(
    row
    for row in poles_payload["sources"]
    if row["ns_id"] == "NS0031"
    and row["source_id"] == "NS0031-S001"
    and row["source_gram_sha256"] == source["gram_sha256"]
)
minimum_basis = sorted(
    pole_entry["minimum_basis"],
    key=lambda row: (int(row["section_pole_order"]), tuple(row["tail"])),
)
if [row["tail"] for row in minimum_basis] != [[1, 0], [1, -1]]:
    raise ValueError("unexpected NS0031 minimum tail basis")
if [int(row["section_pole_order"]) for row in minimum_basis] != [0, 1]:
    raise ValueError("unexpected NS0031 basis pole profile")

root_adapted = matrix(QQ, source["root_adapted_gram"])
root_gram = root_adapted[:15, :15]
basis_tails = [vector(QQ, row["tail"]) for row in minimum_basis]
component_data = []
for component in connected_components(root_gram):
    block = root_gram.matrix_from_rows_and_columns(component, component)
    inverse = block.inverse()
    crosses = []
    for tail in basis_tails:
        crosses.append(
            vector(
                QQ,
                [
                    tail[0] * root_adapted[index, 15]
                    + tail[1] * root_adapted[index, 16]
                    for index in component
                ],
            )
        )
    component_data.append(
        (
            len(component),
            crosses[0] * inverse * crosses[0],
            crosses[1] * inverse * crosses[1],
            crosses[0] * inverse * crosses[1],
        )
    )
expected_component_data = [
    (1, QQ(1) / 2, QQ(0), QQ(0)),
    (7, QQ(0), QQ(7) / 8, QQ(0)),
    (7, QQ(3) / 2, QQ(0), QQ(0)),
]
if component_data != expected_component_data:
    raise ArithmeticError("NS0031 component corrections changed")

prime = int(fibres["prime"])
if prime in (2, 3):
    raise ValueError("bad characteristic for the short Weierstrass chart")
field = GF(prime)
twist = field(arguments.quadratic_twist)
if not twist:
    raise ValueError("quadratic twist must be nonzero")
ring = PolynomialRing(field, "t")
t = ring.gen()

records = []
pole_zero_x_scanned = 0
pole_zero_sections = 0
pole_one_denominators = 0
pole_one_x_scanned = 0
pole_one_sections = 0
marked_pairs = 0
pair_histogram = Counter()
pairs_meeting_singular_fibres = 0

for example_index, example in enumerate(fibres["examples"]):
    A = twist**2 * ring(example["A_coefficients_low_to_high"])
    B = twist**3 * ring(example["B_coefficients_low_to_high"])
    discriminant_core = 4 * A**3 + 27 * B**2
    series_zero, center_zero, node_zero = formal_center(A, B, field.zero(), 3)
    series_one, center_one, node_one = formal_center(A, B, field.one(), 2)
    infinity_ring = PowerSeriesRing(field, "u", default_prec=7)
    local_A_infinity = reversed_local(A, 8, infinity_ring)
    local_B_infinity = reversed_local(B, 12, infinity_ring)
    node_infinity = -field(3) * local_B_infinity[0] / (
        field(2) * local_A_infinity[0]
    )
    center_infinity = infinity_ring(node_infinity)
    for unused in range(7):
        center_infinity = (
            center_infinity + (-local_A_infinity / 3) / center_infinity
        ) / 2
    if (center_infinity**2 + local_A_infinity / 3).valuation() < 5:
        raise ArithmeticError("infinity formal center did not converge")

    # Pole-zero basis vector: depths (1,0,2).  Three X-jet equations on the
    # five polynomial coefficients leave a complete affine plane.
    def pole_zero_jets(poly):
        local_zero = series_zero(poly(t))
        local_infinity = reversed_local(poly, 4, infinity_ring)
        return [local_zero[0], local_infinity[0], local_infinity[1]]

    zero_linear = matrix(
        field, [pole_zero_jets(t**degree) for degree in range(5)]
    ).transpose()
    zero_target = vector(
        field, [center_zero[0], center_infinity[0], center_infinity[1]]
    )
    if zero_linear.rank() != 3:
        raise ArithmeticError("NS0031 pole-zero jet system lost rank")
    zero_particular = zero_linear.solve_right(zero_target)
    zero_kernel = zero_linear.right_kernel().basis()
    if len(zero_kernel) != 2:
        raise ArithmeticError("NS0031 pole-zero chart is not two-dimensional")
    zero_sections = []
    for parameters in itertools.product(field, repeat=2):
        coefficients = zero_particular + sum(
            (parameters[index] * zero_kernel[index] for index in range(2)),
            vector(field, 5),
        )
        X = ring(list(coefficients))
        pole_zero_x_scanned += 1
        for Y in polynomial_roots(X**3 + A * X + B):
            depth_zero = int(
                min(
                    (series_zero(X(t)) - center_zero).valuation(),
                    series_zero(Y(t)).valuation(),
                )
            )
            depth_infinity = int(
                min(
                    (reversed_local(X, 4, infinity_ring) - center_infinity).valuation(),
                    reversed_local(Y, 6, infinity_ring).valuation(),
                )
            )
            smooth_one = not (X(1) == node_one and Y(1) == 0)
            if [depth_zero, depth_infinity] != [1, 2] or not smooth_one:
                continue
            zero_sections.append({"X": X, "Y": Y})
    if not zero_sections:
        continue
    pole_zero_sections += len(zero_sections)

    model_one_sections = []
    for c0 in field:
        C = t + c0
        if C(0) == 0 or C(1) == 0:
            continue
        pole_one_denominators += 1

        # Pole-one basis vector: depths (0,1,0).  Only the value at the finite
        # I8 node is prescribed, leaving a six-dimensional numerator chart.
        values = []
        for degree in range(7):
            local = series_one((t**degree)(t + 1)) / series_one(C(t + 1)) ** 2
            values.append(local[0])
        one_linear = matrix(field, [values])
        one_particular = one_linear.solve_right(vector(field, [center_one[0]]))
        one_kernel = one_linear.right_kernel().basis()
        if len(one_kernel) != 6:
            raise ArithmeticError("NS0031 pole-one chart is not six-dimensional")
        for parameters in itertools.product(field, repeat=6):
            coefficients = one_particular + sum(
                (parameters[index] * one_kernel[index] for index in range(6)),
                vector(field, 7),
            )
            X_numerator = ring(list(coefficients))
            pole_one_x_scanned += 1
            if X_numerator.gcd(C).degree() != 0:
                continue
            right = X_numerator**3 + A * X_numerator * C**4 + B * C**6
            for Y_numerator in polynomial_roots(right):
                if Y_numerator.gcd(C).degree() != 0:
                    continue
                x_zero = series_zero(X_numerator(t)) / series_zero(C(t)) ** 2
                y_zero = series_zero(Y_numerator(t)) / series_zero(C(t)) ** 3
                x_one = series_one(X_numerator(t + 1)) / series_one(C(t + 1)) ** 2
                y_one = series_one(Y_numerator(t + 1)) / series_one(C(t + 1)) ** 3
                x_infinity = reversed_local(X_numerator, 6, infinity_ring) / (
                    reversed_local(C, 1, infinity_ring) ** 2
                )
                y_infinity = reversed_local(Y_numerator, 9, infinity_ring) / (
                    reversed_local(C, 1, infinity_ring) ** 3
                )
                depths = [
                    int(min((x_zero - center_zero).valuation(), y_zero.valuation())),
                    int(min((x_one - center_one).valuation(), y_one.valuation())),
                    int(
                        min(
                            (x_infinity - center_infinity).valuation(),
                            y_infinity.valuation(),
                        )
                    ),
                ]
                if depths != [0, 1, 0]:
                    continue
                pole_one_sections += 1
                pairs = []
                for zero_index, P in enumerate(zero_sections):
                    common = (X_numerator - P["X"] * C**2).gcd(
                        Y_numerator - P["Y"] * C**3
                    )
                    if common.gcd(discriminant_core).degree() != 0:
                        pairs_meeting_singular_fibres += 1
                        continue
                    intersection = int(common.degree())
                    pair_histogram[str(intersection)] += 1
                    if intersection != 2:
                        continue
                    marked_pairs += 1
                    pairs.append(
                        {
                            "pole_zero_index": zero_index,
                            "intersection_on_smooth_fibres": 2,
                            "component_cross_correction": "0",
                            "shioda_height_pairing": "1",
                        }
                    )
                model_one_sections.append(
                    {
                        "C_coefficients_low_to_high": serialize_poly(C),
                        "X_numerator_coefficients_low_to_high": serialize_poly(
                            X_numerator
                        ),
                        "Y_numerator_coefficients_low_to_high": serialize_poly(
                            Y_numerator
                        ),
                        "component_depths_at_I2_I8_I8": depths,
                        "marked_pairs": pairs,
                    }
                )

    records.append(
        {
            "example_index": example_index,
            "pole_zero_sections": [
                {
                    "X_coefficients_low_to_high": serialize_poly(row["X"]),
                    "Y_coefficients_low_to_high": serialize_poly(row["Y"]),
                    "component_depths_at_I2_I8_I8": [1, 0, 2],
                }
                for row in zero_sections
            ],
            "pole_one_sections": model_one_sections,
        }
    )

output = {
    "schema": "elkies-k3.lattice-foundry-ns0031-a1-2a7-marking-modp.v1",
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
        relative(source_path): digest(source_path),
        relative(poles_path): digest(poles_path),
    },
    "source": {
        "ns_id": "NS0031",
        "source_id": "NS0031-S001",
        "source_gram_sha256": source["gram_sha256"],
        "root_type": "A1+2A7",
        "mw_height_gram": source["mw_height_gram"],
        "minimum_basis_tails": [row["tail"] for row in minimum_basis],
        "minimum_basis_pole_profile": [0, 1],
        "component_ranks_and_corrections": [
            {
                "root_rank": rank,
                "pole_zero_self_correction": str(zero_self),
                "pole_one_self_correction": str(one_self),
                "cross_correction": str(cross),
            }
            for rank, zero_self, one_self, cross in component_data
        ],
        "required_smooth_pair_intersection": 2,
    },
    "scope": {
        "fibre_census_exhaustive": True,
        "stored_fibre_models": len(fibres["examples"]),
        "pole_zero_X_candidates_per_model": prime**2,
        "normalized_pole_one_denominators_per_surviving_model": prime - 2,
        "pole_one_X_candidates_per_denominator": prime**6,
        "all_polynomial_Y_square_roots_retained": True,
        "pair_intersections_restricted_to_smooth_fibres": True,
    },
    "accounting": {
        "models_with_pole_zero_section": len(records),
        "pole_zero_X_candidates_scanned": pole_zero_x_scanned,
        "pole_zero_sections": pole_zero_sections,
        "pole_one_denominators_scanned": pole_one_denominators,
        "pole_one_X_candidates_scanned": pole_one_x_scanned,
        "pole_one_sections": pole_one_sections,
        "smooth_pair_intersection_histogram": dict(
            sorted(pair_histogram.items(), key=lambda item: int(item[0]))
        ),
        "pairs_meeting_singular_fibres": pairs_meeting_singular_fibres,
        "marked_mw2_pairs": marked_pairs,
    },
    "models": records,
    "proof_boundary": {
        "proved": (
            "For the displayed prime, twist, normalized fibre supports, and every "
            "model in the exhaustive fibre census, both minimum-basis section "
            "charts and all polynomial Y roots are exhausted. Source corrections "
            "and pole orders are recomputed from the pinned exact lattice artifacts."
        ),
        "not_proved": (
            "A finite-field marked pair is not a characteristic-zero family, a "
            "rational parameterization over Q, or a physical neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_lattice_foundry_ns0031_a1_2a7_marking_modp.sage"
        + (
            f" --fibres {relative(fibres_path)}"
            if fibres_path != DEFAULT_FIBRES.resolve()
            else ""
        )
        + (f" --quadratic-twist {int(twist)}" if twist != 1 else "")
        + (
            f" --output {relative(output_path)}"
            if output_path != DEFAULT_OUTPUT.resolve()
            else ""
        )
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0031 A1+2A7 marking artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0031MARKING|"
    f"p={prime}|twist={int(twist)}|modelsP={len(records)}|P={pole_zero_sections}|"
    f"R={pole_one_sections}|pairs={marked_pairs}|"
    f"status={'PASS' if marked_pairs else 'EMPTY'}",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
