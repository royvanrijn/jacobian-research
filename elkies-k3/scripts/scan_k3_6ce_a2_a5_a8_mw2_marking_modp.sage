#!/usr/bin/env sage-python
"""Exhaust the selected pole-[0,1] MW2 marking on I3+I6+I9 models."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
DEFAULT_FIBRES = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a2-a5-a8-mw2-fibre-ansatz-mod5-v1.json"
DEFAULT_SOURCES = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-semistable-mw0-2-sources-large-a-partner1-v1.json"
CLASSIFIER = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a2-a5-a8-mw2-source-isometries-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a2-a5-a8-mw2-marking-mod5-square-v1.json"
SOURCE_ID = "K3-6ce16abb9de3c7c5-S0005"


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def formal_center(A, B, point, precision):
    field = A.base_ring()
    ring = A.parent()
    t = ring.gen()
    shifted_A = ring(A(t + point))
    shifted_B = ring(B(t + point))
    series_ring = PowerSeriesRing(field, "s", default_prec=precision + 3)
    center = series_ring(-field(3) * shifted_B[0] / (field(2) * shifted_A[0]))
    series_A = series_ring(shifted_A)
    for unused in range(8):
        center = (center + (-series_A / 3) / center) / 2
    if (center**2 + series_A / 3).valuation() < precision + 1:
        raise ArithmeticError("finite formal center did not converge")
    return series_ring, center


def reversed_local(poly, weight, series_ring):
    u = series_ring.gen()
    return series_ring(
        sum(poly[index] * u ** (weight - index) for index in range(poly.degree() + 1))
    )


def polynomial_roots(right):
    if not right.is_square():
        return []
    positive = right.sqrt()
    return [positive] if not positive else [positive, -positive]


def serialize(poly):
    return [int(value) for value in poly]


def main():
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
    classifier = json.loads(CLASSIFIER.read_text())
    source = next(
        row["source"] for row in sources["sources"] if row["source_id"] == SOURCE_ID
    )
    if fibres["ansatz"]["normalized_reducible_supports"] != [
        "0:I3", "1:I6", "infinity:I9"
    ]:
        raise ValueError("marking scan requires normalized I3+I6+I9 fibres")
    if not fibres["scan"]["exhausted"]:
        raise ValueError("marking scan requires an exhaustive fibre census")
    basis = source["pole_audit"]["basis"]
    if not (
        source["root_type"] == "A2+A5+A8"
        and source["mw_rank_for_rho_19"] == 2
        and [section["pole_order"] for section in basis] == [0, 1]
    ):
        raise ValueError("selected source marking changed")

    # These exact profiles are independently recomputed in the classifier.
    depth_profiles = [[0, 2, 1], [1, 2, 4]]
    required_intersection = 1
    selected_profile = next(
        member["marking_profile"]
        for class_row in classifier["classes"]
        for member in class_row["members"]
        if member["source_id"] == SOURCE_ID
    )
    if (
        [
            [row[1] for row in selected_profile["support_profile_I_order_left_depth_right_depth"]],
            [row[2] for row in selected_profile["support_profile_I_order_left_depth_right_depth"]],
        ]
        != depth_profiles
        or selected_profile["required_smooth_pair_intersection"] != required_intersection
    ):
        raise ArithmeticError("classifier marking profile changed")
    prime = int(fibres["prime"])
    field = GF(prime)
    twist = field(arguments.quadratic_twist)
    if not twist:
        raise ValueError("quadratic twist must be nonzero")
    ring = PolynomialRing(field, "t")
    t = ring.gen()

    accounting = {
        "models_with_pole_zero_section": 0,
        "pole_zero_X_candidates_scanned": 0,
        "pole_zero_sections": 0,
        "pole_one_denominators_scanned": 0,
        "pole_one_X_candidates_scanned": 0,
        "pole_one_sections": 0,
        "component_matched_pair_candidates": 0,
        "pairs_meeting_singular_fibres": 0,
        "smooth_pair_intersection_histogram": Counter(),
        "marked_mw2_pairs": 0,
    }
    records = []
    for example_index, example in enumerate(fibres["examples"]):
        A = twist**2 * ring(example["A_coefficients_low_to_high"])
        B = twist**3 * ring(example["B_coefficients_low_to_high"])
        discriminant_core = 4 * A**3 + 27 * B**2
        zero_ring, zero_center = formal_center(A, B, field.zero(), 3)
        one_ring, one_center = formal_center(A, B, field.one(), 4)
        infinity_ring = PowerSeriesRing(field, "u", default_prec=8)
        infinity_A = reversed_local(A, 8, infinity_ring)
        infinity_B = reversed_local(B, 12, infinity_ring)
        inf_center = infinity_ring(
            -field(3) * infinity_B[0] / (field(2) * infinity_A[0])
        )
        for unused in range(8):
            inf_center = (inf_center + (-infinity_A / 3) / inf_center) / 2
        if (inf_center**2 + infinity_A / 3).valuation() < 6:
            raise ArithmeticError("infinity formal center did not converge")

        def zero_section_jets(poly):
            at_one = one_ring(poly(t + 1))
            at_infinity = reversed_local(poly, 4, infinity_ring)
            return [at_one[0], at_one[1], at_infinity[0]]

        zero_linear = matrix(
            field, [zero_section_jets(t**degree) for degree in range(5)]
        ).transpose()
        zero_target = vector(field, [one_center[0], one_center[1], inf_center[0]])
        if zero_linear.rank() != 3:
            raise ArithmeticError("pole-zero jet system lost rank")
        zero_particular = zero_linear.solve_right(zero_target)
        zero_kernel = zero_linear.right_kernel().basis()
        if len(zero_kernel) != 2:
            raise ArithmeticError("pole-zero chart dimension changed")
        zero_sections = []
        for parameters in itertools.product(field, repeat=2):
            solution = zero_particular + sum(
                (coefficient * row for coefficient, row in zip(parameters, zero_kernel)),
                vector(field, 5),
            )
            X = ring(list(solution))
            accounting["pole_zero_X_candidates_scanned"] += 1
            for Y in polynomial_roots(X**3 + A * X + B):
                depths = [
                    int(min((zero_ring(X) - zero_center).valuation(), zero_ring(Y).valuation())),
                    int(min((one_ring(X(t + 1)) - one_center).valuation(), one_ring(Y(t + 1)).valuation())),
                    int(min((reversed_local(X, 4, infinity_ring) - inf_center).valuation(), reversed_local(Y, 6, infinity_ring).valuation())),
                ]
                if depths == depth_profiles[0]:
                    zero_sections.append({"X": X, "Y": Y})
        if not zero_sections:
            continue
        accounting["models_with_pole_zero_section"] += 1
        accounting["pole_zero_sections"] += len(zero_sections)

        pole_one_rows = []
        for c0 in field:
            C = t + c0
            if not C(0) or not C(1):
                continue
            accounting["pole_one_denominators_scanned"] += 1

            def pole_one_jets(poly):
                at_zero = zero_ring(poly) / zero_ring(C) ** 2
                at_one = one_ring(poly(t + 1)) / one_ring(C(t + 1)) ** 2
                local_C = reversed_local(C, 1, infinity_ring)
                at_infinity = reversed_local(poly, 6, infinity_ring) / local_C**2
                return [at_zero[0], at_one[0], at_one[1]] + [
                    at_infinity[index] for index in range(4)
                ]

            one_linear = matrix(
                field, [pole_one_jets(t**degree) for degree in range(7)]
            ).transpose()
            one_target = vector(
                field,
                [zero_center[0], one_center[0], one_center[1]]
                + [inf_center[index] for index in range(4)],
            )
            if one_linear.rank() != 7:
                raise ArithmeticError("pole-one jet system lost rank")
            Xn = ring(list(one_linear.solve_right(one_target)))
            accounting["pole_one_X_candidates_scanned"] += 1
            if Xn.gcd(C).degree() != 0:
                continue
            right = Xn**3 + A * Xn * C**4 + B * C**6
            for Yn in polynomial_roots(right):
                if Yn.gcd(C).degree() != 0:
                    continue
                local_C = reversed_local(C, 1, infinity_ring)
                depths = [
                    int(min((zero_ring(Xn) / zero_ring(C)**2 - zero_center).valuation(), (zero_ring(Yn) / zero_ring(C)**3).valuation())),
                    int(min((one_ring(Xn(t + 1)) / one_ring(C(t + 1))**2 - one_center).valuation(), (one_ring(Yn(t + 1)) / one_ring(C(t + 1))**3).valuation())),
                    int(min((reversed_local(Xn, 6, infinity_ring) / local_C**2 - inf_center).valuation(), (reversed_local(Yn, 9, infinity_ring) / local_C**3).valuation())),
                ]
                if depths != depth_profiles[1]:
                    continue
                accounting["pole_one_sections"] += 1
                pairs = []
                for zero_index, P in enumerate(zero_sections):
                    accounting["component_matched_pair_candidates"] += 1
                    common = (Xn - P["X"] * C**2).gcd(Yn - P["Y"] * C**3)
                    if common.gcd(discriminant_core).degree() != 0:
                        accounting["pairs_meeting_singular_fibres"] += 1
                        continue
                    intersection = int(common.degree())
                    accounting["smooth_pair_intersection_histogram"][intersection] += 1
                    if intersection != required_intersection:
                        continue
                    accounting["marked_mw2_pairs"] += 1
                    pairs.append(
                        {
                            "pole_zero_index": zero_index,
                            "intersection_on_smooth_fibres": intersection,
                            "mw_index_from_height_determinant": 1,
                        }
                    )
                pole_one_rows.append(
                    {
                        "C_coefficients_low_to_high": serialize(C),
                        "X_numerator_coefficients_low_to_high": serialize(Xn),
                        "Y_numerator_coefficients_low_to_high": serialize(Yn),
                        "component_depths_at_I3_I6_I9": depths,
                        "marked_pairs": pairs,
                    }
                )
        records.append(
            {
                "example_index": example_index,
                "pole_zero_sections": [
                    {
                        "X_coefficients_low_to_high": serialize(row["X"]),
                        "Y_coefficients_low_to_high": serialize(row["Y"]),
                        "component_depths_at_I3_I6_I9": depth_profiles[0],
                    }
                    for row in zero_sections
                ],
                "pole_one_sections": pole_one_rows,
            }
        )

    histogram = accounting["smooth_pair_intersection_histogram"]
    accounting["smooth_pair_intersection_histogram"] = {
        str(key): value for key, value in sorted(histogram.items())
    }
    status = (
        "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH_MARKED_MW2_PAIRS"
        if accounting["marked_mw2_pairs"]
        else "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_EMPTY_MARKED_MW2_PAIR_LOCUS"
    )
    payload = {
        "schema": "elkies-k3.k3-6ce-a2-a5-a8-mw2-marking-modp.v1",
        "status": status,
        "prime": prime,
        "quadratic_twist": int(twist),
        "quadratic_twist_square_class": "square" if twist.is_square() else "nonsquare",
        "source": {
            "surface_id": "K3-6ce16abb9de3c7c5",
            "source_id": SOURCE_ID,
            "source_gram_sha256": source["gram_sha256"],
            "root_type": source["root_type"],
            "mw_height_gram": source["mw_height_gram"],
            "minimum_basis_pole_profile": [0, 1],
            "component_depths_at_I3_I6_I9": depth_profiles,
            "required_smooth_pair_intersection": required_intersection,
        },
        "scope": {
            "fibre_census_exhaustive": True,
            "fibre_models": len(fibres["examples"]),
            "pole_zero_X_candidates_per_model": prime**2,
            "normalized_pole_one_denominators_per_surviving_model": prime - 2,
            "pole_one_X_candidates_per_denominator": 1,
            "all_polynomial_Y_square_roots": True,
            "all_pair_intersections_tested": True,
        },
        "accounting": accounting,
        "models": records,
        "inputs": {
            relative(fibres_path): digest(fibres_path),
            relative(sources_path): digest(sources_path),
            relative(CLASSIFIER): digest(CLASSIFIER),
        },
        "proof_boundary": {
            "proved": (
                "Both section charts and every smooth pair intersection are exhausted "
                "over the displayed finite field, twist, and normalized fibre census."
            ),
            "not_proved": (
                "A finite-field pair is not a characteristic-zero lift, rational marking, "
                "primitive determinant-384 specialization, or neighbour route."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/scan_k3_6ce_a2_a5_a8_mw2_marking_modp.sage"
            + (f" --quadratic-twist {int(twist)}" if twist != 1 else "")
            + (f" --output {relative(output_path)}" if output_path != DEFAULT_OUTPUT.resolve() else "")
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "K36CEA2A5A8MARKING|"
        f"p={prime}|twist={int(twist)}|models={len(fibres['examples'])}|"
        f"P={accounting['pole_zero_sections']}|R={accounting['pole_one_sections']}|"
        f"pairs={accounting['marked_mw2_pairs']}|"
        f"status={'PASS' if accounting['marked_mw2_pairs'] else 'EMPTY'}"
    )


if __name__ == "__main__":
    main()
