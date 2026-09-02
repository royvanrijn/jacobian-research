#!/usr/bin/env sage-python
"""Exhaust the pole-one MW marking on the determinant-500 fibre models.

The promoted ``A3+A4+A9/MW1`` source has one generator of height 5/2 and
pole order one.  Its component corrections at the normalized
``I4,I5,I10`` supports are respectively ``1,0,5/2``, hence its component
depths are ``2,0,5``.  These impose two X-jets at zero and five at infinity.
For every monic linear pole denominator this determines the seven
coefficients of the X numerator uniquely; the script then tests every Y
square root and the exact component depths.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIBRES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-fibre-ansatz-mod5-v1.json"
)
DEFAULT_SOURCES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-prescribed-root-sources-large-a-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-pole1-marking-mod5-v1.json"
)
SOURCE_ID = "K3-04b86146cc6b284b-S0160"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
            for right in tuple(sorted(unseen)):
                if gram[left, right]:
                    unseen.remove(right)
                    stack.append(right)
        components.append(sorted(component))
    return components


def formal_center(A, B, point, precision):
    field = A.base_ring()
    ring = A.parent()
    t = ring.gen()
    shifted_A = ring(A(t + point))
    shifted_B = ring(B(t + point))
    node = -field(3) * shifted_B[0] / (field(2) * shifted_A[0])
    series_ring = PowerSeriesRing(field, "s", default_prec=precision + 3)
    center = series_ring(node)
    series_A = series_ring(shifted_A)
    for unused in range(8):
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


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
parser.add_argument("--quadratic-twist", type=int, default=1)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

fibres_path = arguments.fibres.resolve()
sources_path = arguments.sources.resolve()
output_path = arguments.output.resolve()
fibres = json.loads(fibres_path.read_text())
sources = json.loads(sources_path.read_text())
if fibres["ansatz"]["normalized_reducible_supports"] != [
    "0:I4",
    "1:I5",
    "infinity:I10",
]:
    raise ValueError("marking scan requires normalized I4+I5+I10 fibres")
if not fibres["scan"]["exhausted"]:
    raise ValueError("marking scan requires an exhaustive fibre census")
source_row = next(row for row in sources["sources"] if row["source_id"] == SOURCE_ID)
source = source_row["source"]
if source["root_type"] != "A3+A4+A9" or source["mw_height_gram"] != [["5/2"]]:
    raise ValueError("unexpected determinant-500 source lattice")
if source["pole_audit"]["minimum_nonzero_section_pole_order"] != 1:
    raise ValueError("source no longer has a pole-one generator")

root_adapted = matrix(QQ, source["root_adapted_gram"])
root = root_adapted[:16, :16]
labels = vector(QQ, source["pole_audit"]["basis"][0]["simple_root_pairings"])
corrections_by_rank = {}
for component in connected_components(root):
    block = root.matrix_from_rows_and_columns(component, component)
    block_labels = vector(QQ, [labels[index] for index in component])
    corrections_by_rank[len(component)] = block_labels * block.inverse() * block_labels
if corrections_by_rank != {3: QQ(1), 4: QQ(0), 9: QQ(5) / 2}:
    raise ArithmeticError("component correction profile changed")

prime = int(fibres["prime"])
field = GF(prime)
twist = field(arguments.quadratic_twist)
if not twist:
    raise ValueError("quadratic twist must be nonzero")
ring = PolynomialRing(field, "t")
t = ring.gen()

denominators = 0
x_candidates = 0
sections = 0
records = []
depth_histogram = {}
for example_index, example in enumerate(fibres["examples"]):
    A = twist**2 * ring(example["A_coefficients_low_to_high"])
    B = twist**3 * ring(example["B_coefficients_low_to_high"])
    zero_ring, zero_center, zero_node = formal_center(A, B, field.zero(), 3)
    one_ring, one_center, one_node = formal_center(A, B, field.one(), 2)
    infinity_ring = PowerSeriesRing(field, "u", default_prec=9)
    infinity_A = reversed_local(A, 8, infinity_ring)
    infinity_B = reversed_local(B, 12, infinity_ring)
    infinity_node = -field(3) * infinity_B[0] / (field(2) * infinity_A[0])
    infinity_center = infinity_ring(infinity_node)
    for unused in range(8):
        infinity_center = (
            infinity_center + (-infinity_A / 3) / infinity_center
        ) / 2
    if (infinity_center**2 + infinity_A / 3).valuation() < 7:
        raise ArithmeticError("infinity formal center did not converge")

    model_sections = []
    for c0 in field:
        C = t + c0
        if C(0) == 0 or C(1) == 0:
            continue
        denominators += 1

        def jets(poly):
            at_zero = zero_ring(poly(t)) / zero_ring(C(t)) ** 2
            at_infinity = reversed_local(poly, 6, infinity_ring) / (
                reversed_local(C, 1, infinity_ring) ** 2
            )
            return [at_zero[index] for index in range(2)] + [
                at_infinity[index] for index in range(5)
            ]

        linear = matrix(field, [jets(t**degree) for degree in range(7)]).transpose()
        target = vector(
            field,
            [zero_center[index] for index in range(2)]
            + [infinity_center[index] for index in range(5)],
        )
        if linear.rank() != 7:
            raise ArithmeticError("pole-one X-jet system lost rank")
        coefficients = linear.solve_right(target)
        X_numerator = ring(list(coefficients))
        x_candidates += 1
        if X_numerator.gcd(C).degree() != 0:
            continue
        right = X_numerator**3 + A * X_numerator * C**4 + B * C**6
        for Y_numerator in polynomial_roots(right):
            if Y_numerator.gcd(C).degree() != 0:
                continue
            x_zero = zero_ring(X_numerator(t)) / zero_ring(C(t)) ** 2
            y_zero = zero_ring(Y_numerator(t)) / zero_ring(C(t)) ** 3
            x_one = one_ring(X_numerator(t + 1)) / one_ring(C(t + 1)) ** 2
            y_one = one_ring(Y_numerator(t + 1)) / one_ring(C(t + 1)) ** 3
            x_infinity = reversed_local(X_numerator, 6, infinity_ring) / (
                reversed_local(C, 1, infinity_ring) ** 2
            )
            y_infinity = reversed_local(Y_numerator, 9, infinity_ring) / (
                reversed_local(C, 1, infinity_ring) ** 3
            )
            depths = [
                int(min((x_zero - zero_center).valuation(), y_zero.valuation())),
                int(min((x_one - one_center).valuation(), y_one.valuation())),
                int(
                    min(
                        (x_infinity - infinity_center).valuation(),
                        y_infinity.valuation(),
                    )
                ),
            ]
            key = ",".join(map(str, depths))
            depth_histogram[key] = depth_histogram.get(key, 0) + 1
            if depths != [2, 0, 5]:
                continue
            sections += 1
            model_sections.append(
                {
                    "C_coefficients_low_to_high": serialize_poly(C),
                    "X_numerator_coefficients_low_to_high": serialize_poly(X_numerator),
                    "Y_numerator_coefficients_low_to_high": serialize_poly(Y_numerator),
                    "component_depths_at_I4_I5_I10": depths,
                }
            )
    if model_sections:
        records.append(
            {"example_index": example_index, "pole_one_sections": model_sections}
        )

payload = {
    "schema": "elkies-k3.k3-04b-a3-a4-a9-pole1-marking-modp.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH_MARKED_MW1_SECTION"
        if sections
        else "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_EMPTY_MARKED_MW1_LOCUS"
    ),
    "prime": prime,
    "quadratic_twist": int(twist),
    "quadratic_twist_square_class": "square" if twist.is_square() else "nonsquare",
    "source": {
        "surface_id": "K3-04b86146cc6b284b",
        "source_id": SOURCE_ID,
        "source_gram_sha256": source["gram_sha256"],
        "root_type": "A3+A4+A9",
        "mw_height_gram": [["5/2"]],
        "minimum_basis_pole_profile": [1],
        "component_corrections_at_I4_I5_I10": ["1", "0", "5/2"],
        "component_depths_at_I4_I5_I10": [2, 0, 5],
    },
    "scope": {
        "fibre_census_exhaustive": True,
        "fibre_models": len(fibres["examples"]),
        "normalized_pole_denominators_per_model": prime - 2,
        "X_numerators_per_denominator": 1,
        "all_polynomial_Y_square_roots_retained": True,
    },
    "accounting": {
        "pole_denominators_scanned": denominators,
        "X_numerators_scanned": x_candidates,
        "exact_depth_histogram_for_square_sections": dict(sorted(depth_histogram.items())),
        "marked_mw1_sections": sections,
        "models_with_marked_section": len(records),
    },
    "models": records,
    "inputs": {
        relative(fibres_path): digest(fibres_path),
        relative(sources_path): digest(sources_path),
    },
    "proof_boundary": {
        "proved": (
            "For every model in the exhaustive normalized fibre census and every "
            "monic linear pole denominator away from the reducible supports, the "
            "seven component X-jets determine and test the unique numerator; all "
            "polynomial Y square roots and exact component depths are retained."
        ),
        "not_proved": (
            "A finite-field marked section is not a characteristic-zero lift, a "
            "rational parameterization over Q, or a physical neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_k3_04b_a3_a4_a9_pole1_marking_modp.sage"
        + (
            f" --quadratic-twist {int(twist)}"
            if twist != 1
            else ""
        )
        + (
            f" --output {relative(output_path)}"
            if output_path != DEFAULT_OUTPUT.resolve()
            else ""
        )
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("determinant-500 pole-one marking artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)
print(
    "K304BEMARKING|"
    f"p={prime}|twist={int(twist)}|models={len(fibres['examples'])}|"
    f"sections={sections}|status={'PASS' if sections else 'EMPTY'}",
    flush=True,
)
