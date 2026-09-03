#!/usr/bin/env python3
"""Derive the leakage-aware R17 visibility-complexity baseline.

This is intentionally a small certificate join, not a multisection search.  It
combines the exact rational-bisection quotient ledger with the exact
target-fitted genus-one pencil certificate for the 38 displayed generators at
the rank-25--28 controls.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FILTRATION_PATH = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-multisection-visibility-filtration-v1.json"
)
QUOTIENT_PATH = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/elkies_2026_deep_cover_exceptional_quotients_v1.json"
)
RELATIONS_PATH = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/elkies_2026_exceptional_specialization_relations_v1.json"
)
OUTPUT_PATH = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-visibility-complexity-v1.json"
)


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def primitive_integer_coefficients(values: list[str | int]) -> list[int]:
    coefficients = [Fraction(value) for value in values]
    denominator = math.lcm(*(value.denominator for value in coefficients))
    integers = [value.numerator * (denominator // value.denominator) for value in coefficients]
    content = math.gcd(*(abs(value) for value in integers))
    assert content > 0
    integers = [value // content for value in integers]
    if next(value for value in reversed(integers) if value) < 0:
        integers = [-value for value in integers]
    return integers


def branch_equation_cost(coefficients: list[int]) -> dict:
    nonzero = [abs(value) for value in coefficients if value]
    assert nonzero
    degree = max(index for index, value in enumerate(coefficients) if value)
    return {
        "metric": (
            "primitive branch polynomial over QQ; tuple ordered as "
            "(support, maximum coefficient bits, total coefficient bits)"
        ),
        "polynomial_degree": degree,
        "coefficient_support": len(nonzero),
        "maximum_absolute_coefficient_bits": max(value.bit_length() for value in nonzero),
        "total_absolute_coefficient_bits": sum(value.bit_length() for value in nonzero),
        "primitive_integer_coefficients_low_to_high": coefficients,
    }


def cost_key(cost: dict) -> tuple[int, int, int]:
    return (
        cost["coefficient_support"],
        cost["maximum_absolute_coefficient_bits"],
        cost["total_absolute_coefficient_bits"],
    )


def build() -> dict:
    filtration = load(FILTRATION_PATH)
    quotient = load(QUOTIENT_PATH)
    relations = load(RELATIONS_PATH)

    filtration_by_parameter = {
        fibre["parameter"]: fibre
        for fibre in filtration["literal_all_genus_filtration"]["fibres"]
    }
    quotient_by_parameter = {fibre["parameter"]: fibre for fibre in quotient["fibres"]}
    relations_by_parameter = {fibre["parameter"]: fibre for fibre in relations["fibres"]}

    rows = []
    summaries = []
    for parameter, post_hoc_fibre in filtration_by_parameter.items():
        quotient_fibre = quotient_by_parameter[parameter]
        relation_fibre = relations_by_parameter[parameter]
        ordered_basis = quotient_fibre["exceptional_quotient"]["ordered_basis"]
        targets = post_hoc_fibre["targets"]
        assert len(ordered_basis) == len(targets)

        relation_rows = relation_fibre[
            "split_rational_bisection_specialization_equations"
        ]
        relation_by_label = {
            item["cover_label"]: item
            for item in relation_rows
        }
        point_by_label = {
            item["cover_label"]: item
            for item in quotient_fibre["points"]
            if item["cover_degree"] == 2
        }
        assert set(relation_by_label) == set(point_by_label)
        direct_relations_by_index = {
            index: [
                item
                for item in relation_rows
                if all(
                    int(value) == (1 if row == index else 0)
                    or int(value) == (-1 if row == index else 0)
                    for row, value in enumerate(item["exceptional_quotient_vector"])
                )
            ]
            for index in range(len(ordered_basis))
        }
        relation_cost_by_label = {
            label: branch_equation_cost(
                primitive_integer_coefficients(
                    item["primitive_integer_q_coefficients_low_to_high"]
                )
            )
            for item in relation_fibre["split_rational_bisection_specialization_equations"]
            for label in [item["cover_label"]]
        }
        rigid_count = 0
        for index, (label, target) in enumerate(zip(ordered_basis, targets)):
            assert label == target["target_label"]
            post_hoc_coefficients = primitive_integer_coefficients(
                target["branch_polynomial_q_coefficients_low_to_high"]
            )
            post_hoc_cost = branch_equation_cost(post_hoc_coefficients)
            pencil_parameter = Fraction(target["pencil_parameter_lambda"])
            post_hoc = {
                "status": "PROVED_POST_HOC_INCIDENCE",
                "mechanism_class": "post-hoc",
                "genus": 1,
                "multisection_degree": 2,
                "linear_system_dimension": 1,
                "incidence_codimension": 1,
                "coefficient_field": {"field": "QQ", "degree_over_QQ": 1},
                "common_trace": "-P2-P5",
                "pencil_parameter_lambda": str(pencil_parameter),
                "pencil_parameter_height_bits": max(
                    abs(pencil_parameter.numerator).bit_length(),
                    pencil_parameter.denominator.bit_length(),
                ),
                "branch_equation_cost": post_hoc_cost,
            }

            rigid = None
            direct_relations = direct_relations_by_index[index]
            if direct_relations:
                rigid_count += 1
                relation = min(
                    direct_relations,
                    key=lambda item: cost_key(
                        relation_cost_by_label[item["cover_label"]]
                    ),
                )
                cover_label = relation["cover_label"]
                point = point_by_label[cover_label]
                rigid_cost = relation_cost_by_label[cover_label]
                rigid = {
                    "status": "PROVED_MINIMUM",
                    "mechanism_class": "rigid",
                    "genus": 0,
                    "multisection_degree": 2,
                    "linear_system_dimension": 0,
                    "incidence_codimension": 0,
                    "coefficient_field": {"field": "QQ", "degree_over_QQ": 1},
                    "cover_label": cover_label,
                    "exceptional_quotient_vector": relation[
                        "exceptional_quotient_vector"
                    ],
                    "minimum_lattice_norm": point["minimum_lattice_norm"],
                    "branch_equation_cost": rigid_cost,
                    "existing_construction_complexity": point[
                        "equation_complexity"
                    ],
                }

            known_costs = [("post-hoc", post_hoc_cost)]
            if rigid is not None:
                known_costs.append(("rigid", rigid["branch_equation_cost"]))
            minimum_cost_class, minimum_cost = min(
                known_costs, key=lambda item: cost_key(item[1])
            )

            row = {
                "parameter": parameter,
                "known_rank_lower_bound": 17 + len(ordered_basis),
                "direction": label,
                "source_public_point_index_one_based": target[
                    "source_public_point_index_one_based"
                ],
                "minimum_known_rational_curve_degree": (
                    {
                        "status": "PROVED_MINIMUM",
                        "degree": 2,
                        "reason_degree_one_is_impossible": (
                            "Degree-one curves are sections and specialize into M_t, "
                            "so cannot represent a nonzero class in L_t/M_t."
                        ),
                    }
                    if rigid is not None
                    else {
                        "status": "UNKNOWN",
                        "degree": None,
                        "search_boundary": (
                            "The rational degree-two atlas is complete; the degree-three "
                            "and degree-four equation atlases are incomplete."
                        ),
                    }
                ),
                "minimum_known_predeclared_trace_pencil_degree": {
                    "status": "UNKNOWN_NO_LEAKAGE_FREE_LEDGER",
                    "degree": None,
                    "reason": (
                        "The available genus-one certificate chooses the pencil member "
                        "using this already known target. It cannot be promoted "
                        "retrospectively to a predeclared predictor."
                    ),
                },
                "minimum_known_incidence_codimension": 0 if rigid is not None else 1,
                "minimum_known_branch_equation_cost": {
                    "witness_class": minimum_cost_class,
                    **minimum_cost,
                },
                "rigid_visibility": rigid,
                "post_hoc_visibility": post_hoc,
                "predictive_visibility_status": (
                    "RIGID_VISIBLE" if rigid is not None else "UNKNOWN"
                ),
            }
            rows.append(row)

        dimension = len(ordered_basis)
        summaries.append(
            {
                "parameter": parameter,
                "known_rank_lower_bound": 17 + dimension,
                "displayed_generator_count": dimension,
                "individually_rigid_visible_generator_count": rigid_count,
                "rational_bisection_visible_span_dimension": quotient_fibre[
                    "minimum_norm_and_equation_complexity_exposure"
                ]["final_rank"],
                "certified_predeclared_pencil_visible_generator_count": 0,
                "post_hoc_genus_one_degree_two_generator_count": dimension,
            }
        )

    summaries.sort(key=lambda item: item["known_rank_lower_bound"])
    rows.sort(
        key=lambda item: (
            item["known_rank_lower_bound"],
            int(item["direction"].removeprefix("Q")),
        )
    )
    assert [item["displayed_generator_count"] for item in summaries] == [8, 9, 10, 11]
    assert [item["rational_bisection_visible_span_dimension"] for item in summaries] == [
        5,
        3,
        2,
        1,
    ]
    assert [item["individually_rigid_visible_generator_count"] for item in summaries] == [
        4,
        0,
        0,
        0,
    ]
    assert len(rows) == 38
    assert sum(row["rigid_visibility"] is not None for row in rows) == 4

    return {
        "schema": "elkies-k3-r17-visibility-complexity-v1",
        "status": "PASS_EXACT_DERIVED_BASELINE",
        "definition": {
            "object": (
                "For a leakage-free atlas A fixed before the target point Q and "
                "control fibre t, VC_A(Q) is the Pareto frontier of witness tuples."
            ),
            "witness_tuple": [
                "mechanism class (categorical: rigid, predeclared-pencil, post-hoc)",
                "geometric genus",
                "multisection degree",
                "linear-system dimension",
                "incidence codimension",
                "branch-equation cost",
                "coefficient-field complexity",
            ],
            "comparison_rule": (
                "Do not compress the tuple to a scalar without declaring an order or "
                "weights before inspecting the controls."
            ),
            "predictive_classes": ["rigid", "predeclared-pencil"],
            "diagnostic_only_class": "post-hoc",
        },
        "equation_cost_convention": {
            "normalization": (
                "Clear rational denominators in the branch polynomial, divide by "
                "integer content, and make the leading coefficient positive."
            ),
            "ordered_tuple": [
                "coefficient_support",
                "maximum_absolute_coefficient_bits",
                "total_absolute_coefficient_bits",
            ],
            "minimum_known_cost_rule": "lexicographic in the displayed ordered tuple",
        },
        "inputs": {
            str(path.relative_to(ROOT)): {"sha256": sha256(path)}
            for path in (FILTRATION_PATH, QUOTIENT_PATH, RELATIONS_PATH)
        },
        "summary_by_rank": summaries,
        "directions": rows,
        "claim_boundary": {
            "exact": (
                "Four named generators have minimum rational-curve degree two; all 38 "
                "have post-hoc genus-one degree-two, codimension-one witnesses."
            ),
            "unknown": (
                "Rational-curve degree remains unknown for 34 named generators, and "
                "predeclared-pencil degree remains uncertified for all 38."
            ),
            "span_warning": (
                "The rational-bisection span dimensions 5,3,2,1 do not mean that "
                "5,3,2,1 ordered basis generators each have one rational-curve witness."
            ),
            "no_exhaustive_campaign": (
                "No norm-20 rational-trisection cosets or M/4M quadrisection cosets "
                "are enumerated by this derivation."
            ),
        },
        "conjecture": (
            "The high-rank controls have increasing rigid visibility complexity, "
            "although unrestricted post-hoc all-genus visibility is uniformly degree two."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    rendered = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.check:
        assert OUTPUT_PATH.read_text() == rendered, f"stale output: {OUTPUT_PATH}"
    else:
        OUTPUT_PATH.write_text(rendered)
    print(
        "R17VISIBILITYCOMPLEXITY|directions=38|rigid_individual=4|"
        "predeclared_certified=0|post_hoc_degree2=38|status=PASS_EXACT_DERIVED_BASELINE"
    )


if __name__ == "__main__":
    main()
