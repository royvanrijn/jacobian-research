#!/usr/bin/env sage-python
"""Exhaust the NS0005 pole-one/pole-zero MW basis on stored fibre models.

The exact ``NS0005-S001`` source has root type ``A1+2A7`` and MW height
Gram ``[[5/2,-1/4],[-1/4,25/8]]``.  In the selected root-component ordering,
the pole-one generator has component depths ``(1,2,2)`` at the normalized
``I2,I8,I8`` supports.  The pole-zero generator has depths ``(0,0,1)``.

A pole-one section is written globally as ``x=N/C^2, y=M/C^3`` with a
projective linear denominator C.  Requiring C to be nonzero at all three
reducible supports leaves three normalized denominators over GF(5).  The five
component-jet conditions leave an affine two-dimensional numerator-X chart,
which this script exhausts.  A retained pair must also have smooth-fibre
intersection number three, giving the required Shioda pairing ``-1/4`` after
the exact component cross-correction ``1/4``.
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
DEFAULT_Q_SECTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-infinity-pole0-sections-mod5-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-pole1-pole0-pairs-mod5-v1.json"
)
DEFAULT_SOURCE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json"
)
DEFAULT_POLES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-rank2-section-basis-poles-v1.json"
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


def serialize_poly(poly):
    return [int(value) for value in poly]


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


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
parser.add_argument("--q-sections", type=Path, default=DEFAULT_Q_SECTIONS)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument("--poles", type=Path, default=DEFAULT_POLES)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

fibres_path = arguments.fibres.resolve()
q_sections_path = arguments.q_sections.resolve()
source_path = arguments.source.resolve()
poles_path = arguments.poles.resolve()
output_path = arguments.output.resolve()
fibres = json.loads(fibres_path.read_text())
q_payload = json.loads(q_sections_path.read_text())
source_payload = json.loads(source_path.read_text())
poles_payload = json.loads(poles_path.read_text())
if (
    fibres.get("schema")
    != "elkies-k3.lattice-foundry-three-support-semistable-source-ansatz-modp.v1"
):
    raise ValueError("unexpected NS0005 fibre schema")
if (
    q_payload.get("schema")
    != "elkies-k3.lattice-foundry-three-support-infinity-pole0-sections-modp.v1"
):
    raise ValueError("unexpected NS0005 pole-zero section schema")
if q_payload["input"]["sha256"] != digest(fibres_path):
    raise ValueError("pole-zero section artifact is not attached to the fibre input")
if not fibres["scan"]["exhausted"] or not q_payload["scope"]["fibre_ansatz_scan_exhausted"]:
    raise ValueError("pair scan requires an exhaustive fibre census")
if fibres["source"].get("ns_id") != "NS0005" or fibres["source"].get("root_type") != "A1+2A7":
    raise ValueError("pair scan is specialized to NS0005-S001/A1+2A7")
if fibres["source"].get("mw_height_gram") != [["5/2", "-1/4"], ["-1/4", "25/8"]]:
    raise ValueError("unexpected NS0005 MW height Gram")

source_entry = next(
    row
    for row in source_payload["sources"]
    if row["ns_id"] == "NS0005"
    and row["source_id"] == "NS0005-S001"
    and row["source"]["gram_sha256"] == fibres["source"]["source_gram_sha256"]
)
source = source_entry["source"]
pole_entry = next(
    row
    for row in poles_payload["sources"]
    if row["ns_id"] == "NS0005"
    and row["source_id"] == "NS0005-S001"
    and row["source_gram_sha256"] == source["gram_sha256"]
)
if pole_entry["minimum_basis_sorted_pole_profile"] != [0, 1]:
    raise ValueError("unexpected NS0005 complete-basis pole profile")
if {tuple(row["tail"]) for row in pole_entry["minimum_basis"]} != {(1, 0), (0, 1)}:
    raise ValueError("minimum pole basis does not use the displayed tail basis")

root_adapted = matrix(QQ, source["root_adapted_gram"])
root_gram = root_adapted[:15, :15]
component_corrections = []
for component in connected_components(root_gram):
    block = root_gram.matrix_from_rows_and_columns(component, component)
    p_cross = vector(QQ, [root_adapted[index, 15] for index in component])
    q_cross = vector(QQ, [root_adapted[index, 16] for index in component])
    inverse = block.inverse()
    component_corrections.append(
        (
            len(component),
            p_cross * inverse * p_cross,
            q_cross * inverse * q_cross,
            p_cross * inverse * q_cross,
        )
    )
expected_component_corrections = [
    (1, QQ(1) / 2, QQ(0), QQ(0)),
    (7, QQ(3) / 2, QQ(0), QQ(0)),
    (7, QQ(3) / 2, QQ(7) / 8, QQ(1) / 4),
]
if component_corrections != expected_component_corrections:
    raise ArithmeticError("root-adapted component corrections changed")

prime = int(fibres["prime"])
if prime in (2, 3) or int(q_payload["prime"]) != prime:
    raise ValueError("fibre and section artifacts must use the same good prime")
field = GF(prime)
twist = field(q_payload["quadratic_twist"])
if not twist:
    raise ValueError("quadratic twist must be nonzero")
ring = PolynomialRing(field, "t")
t = ring.gen()

q_models = {
    int(row["example_index"]): row
    for row in q_payload["models"]
    if row["Q_marked_section_count"]
}
if not q_models:
    raise ValueError("no pole-zero section survives on the supplied models")

records = []
total_denominators = 0
total_x_candidates = 0
total_square_sections = 0
total_pairs = 0
pair_intersection_histogram = Counter()
pairs_meeting_singular_fibres = 0
for example_index, q_record in sorted(q_models.items()):
    example = fibres["examples"][example_index]
    A = twist**2 * ring(example["A_coefficients_low_to_high"])
    B = twist**3 * ring(example["B_coefficients_low_to_high"])
    discriminant_core = 4 * A**3 + 27 * B**2
    series_zero, center_zero, node_zero = formal_center(A, B, field.zero(), 2)
    series_one, center_one, node_one = formal_center(A, B, field.one(), 3)
    infinity_ring = PowerSeriesRing(field, "u", default_prec=6)
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
    if (center_infinity**2 + local_A_infinity / 3).valuation() < 4:
        raise ArithmeticError("infinity formal center did not converge")

    q_sections = [
        {
            "X": ring(row["X_coefficients_low_to_high"]),
            "Y": ring(row["Y_coefficients_low_to_high"]),
            "serialized": row,
        }
        for row in q_record["Q_sections"]
    ]
    model_sections = []

    # C=t+c0 is the unique affine normalization with nonzero infinity value.
    # Exclude roots at zero and one because P is prescribed nonidentity there.
    for c0 in field:
        C = t + c0
        if C(0) == 0 or C(1) == 0:
            continue
        total_denominators += 1

        def jet_values(poly):
            local_zero = series_zero(poly(t)) / series_zero(C(t)) ** 2
            local_one = series_one(poly(t + 1)) / series_one(C(t + 1)) ** 2
            local_infinity = reversed_local(poly, 6, infinity_ring) / (
                reversed_local(C, 1, infinity_ring) ** 2
            )
            return [
                local_zero[0],
                local_one[0],
                local_one[1],
                local_infinity[0],
                local_infinity[1],
            ]

        numerator_basis = [t**degree for degree in range(7)]
        linear = matrix(field, [jet_values(poly) for poly in numerator_basis]).transpose()
        target = vector(
            field,
            [
                center_zero[0],
                center_one[0],
                center_one[1],
                center_infinity[0],
                center_infinity[1],
            ],
        )
        if linear.rank() != 5:
            raise ArithmeticError("component-jet numerator system lost rank")
        particular = linear.solve_right(target)
        kernel = linear.right_kernel().basis()
        if len(kernel) != 2:
            raise ArithmeticError("pole-one numerator chart is not two-dimensional")

        for parameters in itertools.product(field, repeat=2):
            coefficients = particular + sum(
                (parameters[index] * kernel[index] for index in range(2)),
                vector(field, 7),
            )
            X_numerator = ring(list(coefficients))
            total_x_candidates += 1
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
                if depths != [1, 2, 2]:
                    continue
                total_square_sections += 1
                pairs = []
                for q_index, Q in enumerate(q_sections):
                    common = (X_numerator - Q["X"] * C**2).gcd(
                        Y_numerator - Q["Y"] * C**3
                    )
                    if common.gcd(discriminant_core).degree() != 0:
                        pairs_meeting_singular_fibres += 1
                        continue
                    intersection = int(common.degree())
                    pair_intersection_histogram[str(intersection)] += 1
                    if intersection != 3:
                        continue
                    pairs.append(
                        {
                            "Q_index": q_index,
                            "intersection_on_smooth_fibres": intersection,
                            "component_cross_correction": "1/4",
                            "shioda_height_pairing": "-1/4",
                        }
                    )
                total_pairs += len(pairs)
                model_sections.append(
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
            "pole_zero_section_count": len(q_sections),
            "pole_one_sections_with_required_components": len(model_sections),
            "pole_one_sections": model_sections,
        }
    )

output = {
    "schema": "elkies-k3.lattice-foundry-ns0005-pole1-pole0-pairs-modp.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_STORED_MODELS_WITH_MARKED_MW2_PAIRS"
        if total_pairs
        else "PASS_EXACT_EXHAUSTIVE_STORED_MODELS_EMPTY_MARKED_MW2_PAIR_CHART"
    ),
    "prime": prime,
    "quadratic_twist": int(twist),
    "quadratic_twist_square_class": q_payload["quadratic_twist_square_class"],
    "inputs": {
        relative(fibres_path): digest(fibres_path),
        relative(q_sections_path): digest(q_sections_path),
        relative(source_path): digest(source_path),
        relative(poles_path): digest(poles_path),
    },
    "source": {
        "ns_id": "NS0005",
        "source_id": fibres["source"]["source_id"],
        "root_type": "A1+2A7",
        "mw_height_gram": fibres["source"]["mw_height_gram"],
        "determinant": int(fibres["source"]["determinant"]),
        "displayed_tail_basis_pole_profile": [1, 0],
        "sorted_basis_pole_profile": [0, 1],
        "source_gram_sha256": source["gram_sha256"],
        "component_ranks_and_corrections": [
            {
                "root_rank": rank,
                "pole_one_self_correction": str(p_self),
                "pole_zero_self_correction": str(q_self),
                "cross_correction": str(cross),
            }
            for rank, p_self, q_self, cross in component_corrections
        ],
    },
    "scope": {
        "fibre_census_exhaustive": True,
        "stored_models_with_pole_zero_generator": len(q_models),
        "normalized_projective_linear_denominators_per_model": prime - 2,
        "pole_one_X_numerator_affine_dimension": 2,
        "pole_one_X_numerator_candidates_per_denominator": prime**2,
        "all_polynomial_Y_square_roots_retained": True,
        "pair_intersections_restricted_to_smooth_fibres": True,
    },
    "accounting": {
        "denominators_scanned": total_denominators,
        "pole_one_X_numerators_scanned": total_x_candidates,
        "pole_one_sections_with_required_components": total_square_sections,
        "smooth_pair_intersection_histogram": dict(
            sorted(pair_intersection_histogram.items(), key=lambda item: int(item[0]))
        ),
        "pairs_meeting_singular_fibres": pairs_meeting_singular_fibres,
        "marked_mw2_pairs": total_pairs,
    },
    "models": records,
    "proof_boundary": {
        "proved": (
            f"On every stored exhaustive GF({prime}) fibre model carrying the pole-zero "
            "generator, all normalized pole-one denominators, the complete affine "
            "component-jet numerator charts, and all polynomial Y roots are "
            "exhausted. The component corrections and basis pole profile are "
            "recomputed from hash-pinned exact lattice artifacts. Retained pairs "
            "have the exact component depths and smooth-fibre intersection required "
            "by the displayed MW height Gram."
        ),
        "not_proved": (
            "A finite-field marked pair is not a characteristic-zero family, a "
            "rational parameterization over Q, or a physical neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_lattice_foundry_ns0005_pole1_pole0_pairs_modp.sage"
        + (
            f" --fibres {relative(fibres_path)} --q-sections {relative(q_sections_path)} "
            f"--output {relative(output_path)}"
            if (
                fibres_path != DEFAULT_FIBRES.resolve()
                or q_sections_path != DEFAULT_Q_SECTIONS.resolve()
                or output_path != DEFAULT_OUTPUT.resolve()
            )
            else ""
        )
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0005 pole-one/pole-zero pair artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0005PAIR|"
    f"models={len(q_models)}|denominators={total_denominators}|"
    f"X={total_x_candidates}|P={total_square_sections}|pairs={total_pairs}|"
    f"status={'PASS' if total_pairs else 'EMPTY'}",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
