#!/usr/bin/env python3
"""Finite-field Schur scout for cubic target completion of K12 graphs.

The exact bilinear theorem supplies a full-column pivot for the five graph
families z4,z9,z10,z11,z12.  At each parameter point this script first builds
that certified degree-at-most-two column space, then adds the 220 cubic target
monomials only in its quotient.  It tests whether every high-source-degree
retained component can be repaired through target degree three.

This is a bounded discovery experiment.  A survivor needs rational
reconstruction and a complete target-automorphism audit; a search with no
survivor is not a characteristic-zero obstruction.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

import sympy as sp

from audit_k12_coordinate_pair_frontier import build_k12
from search_k12_cubic_graph_bilinear_completions import (
    DESIRED_DEGREE,
    GraphFamily,
    SparseColumnSpace,
    SparsePolynomial,
    as_sparse_rational,
    build_graph_mod_prime,
    high_degree_part,
    linear_graph_families,
    parameter_points,
    polynomial_degree,
    restrict_polynomial,
    signed_mod,
    sparse_multiply,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "k12_schur_cubic_completion_modular_search.json"
)

# Zero-based family pivots and the retained components used by the exact
# constant-minor bilinear certificates.
SELECTED_COMPONENTS = {
    3: 2,
    8: 1,
    9: 2,
    10: 0,
    11: 0,
}

EXPECTED_BILINEAR_RANK = {
    3: 56,
    8: 55,
    9: 55,
    10: 55,
    11: 55,
}


def target_product(
    available: list[SparsePolynomial],
    indices: tuple[int, ...],
    prime: int,
) -> SparsePolynomial:
    product = available[indices[0]]
    for index in indices[1:]:
        product = sparse_multiply(product, available[index], prime)
    return product


def schur_screen_component(
    restricted: list[SparsePolynomial],
    component: int,
    prime: int,
) -> dict[str, object]:
    """Test a component after quotienting by linear/bilinear columns."""

    available = [
        polynomial
        for index, polynomial in enumerate(restricted)
        if index != component
    ]
    assert len(available) == 10
    space = SparseColumnSpace(prime)
    nonzero_by_degree: dict[int, int] = {}

    for degree in (1, 2):
        nonzero = 0
        for indices in itertools.combinations_with_replacement(
            range(10), degree
        ):
            column = high_degree_part(
                target_product(available, indices, prime)
            )
            if column:
                nonzero += 1
                space.add(column)
        nonzero_by_degree[degree] = nonzero

    bilinear_rank = space.rank
    target = high_degree_part(restricted[component])
    bilinear_residual = space.reduce(target)

    nonzero_cubics = 0
    for indices in itertools.combinations_with_replacement(range(10), 3):
        column = high_degree_part(
            target_product(available, indices, prime)
        )
        if column:
            nonzero_cubics += 1
            space.add(column)

    total_rank = space.rank
    cubic_residual = space.reduce(target)
    return {
        "retained_component": component,
        "target_term_count": len(target),
        "formal_linear_count": 10,
        "formal_quadratic_count": 55,
        "formal_cubic_count": 220,
        "nonzero_linear_count": nonzero_by_degree[1],
        "nonzero_quadratic_count": nonzero_by_degree[2],
        "nonzero_cubic_count": nonzero_cubics,
        "bilinear_rank": bilinear_rank,
        "bilinear_residual_term_count": len(bilinear_residual),
        "bilinear_residual_leading_monomial": (
            list(min(bilinear_residual)) if bilinear_residual else None
        ),
        "cubic_quotient_rank": total_rank - bilinear_rank,
        "total_completion_rank": total_rank,
        "cubic_solvable": not cubic_residual,
        "cubic_residual_term_count": len(cubic_residual),
        "cubic_residual_leading_monomial": (
            list(min(cubic_residual)) if cubic_residual else None
        ),
    }


def audit_point(
    family: GraphFamily,
    integer_point: tuple[int, ...],
    prime: int,
    k12_sparse: list[dict[tuple[int, ...], sp.Rational]],
    nonlinear_sparse: list[dict[tuple[int, ...], sp.Rational]],
    collision_image: tuple[sp.Rational, ...],
) -> dict[str, object]:
    modular_point = tuple(signed_mod(value, prime) for value in integer_point)
    graph, coefficients = build_graph_mod_prime(
        family,
        modular_point,
        nonlinear_sparse,
        collision_image,
        prime,
    )
    restricted = [
        restrict_polynomial(polynomial, graph, family.pivot, prime)
        for index, polynomial in enumerate(k12_sparse)
        if index != family.pivot
    ]
    degrees = [polynomial_degree(polynomial) for polynomial in restricted]
    bad = [index for index, degree in enumerate(degrees) if degree > 3]
    selected = SELECTED_COMPONENTS[family.pivot]
    assert selected in bad

    ordered_bad = [selected] + [index for index in bad if index != selected]
    screens = []
    for component in ordered_bad:
        screen = schur_screen_component(restricted, component, prime)
        screens.append(screen)
        if not screen["cubic_solvable"]:
            break
    all_solvable = len(screens) == len(bad) and all(
        bool(screen["cubic_solvable"]) for screen in screens
    )
    selected_screen = screens[0]
    assert selected_screen["bilinear_rank"] == EXPECTED_BILINEAR_RANK[
        family.pivot
    ]
    assert selected_screen["bilinear_residual_term_count"] > 0
    return {
        "parameter_point": list(integer_point),
        "linear_target_coefficients_mod_prime": list(coefficients),
        "graph_degree": polynomial_degree(graph),
        "graph_term_count": len(graph),
        "restricted_maximum_degree": max(degrees),
        "bad_retained_components": bad,
        "screened_components": screens,
        "all_bad_components_cubically_solvable": all_solvable,
        "score": [
            0 if all_solvable else 1,
            selected_screen["cubic_residual_term_count"],
            len(bad),
            len(graph),
        ],
    }


def parse_integer_list(raw: str) -> tuple[int, ...]:
    values = tuple(dict.fromkeys(int(value.strip()) for value in raw.split(",")))
    if not values or 0 in values:
        raise argparse.ArgumentTypeError(
            "values must be a nonempty list excluding zero"
        )
    return values


def parse_primes(raw: str) -> tuple[int, ...]:
    primes = tuple(dict.fromkeys(int(value.strip()) for value in raw.split(",")))
    if not primes or any(not sp.isprime(prime) for prime in primes):
        raise argparse.ArgumentTypeError(
            "--primes requires a nonempty list of primes"
        )
    return primes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primes", type=parse_primes, default=(101, 103))
    parser.add_argument("--support-max", type=int, default=1)
    parser.add_argument("--values", type=parse_integer_list, default=(-2, -1, 1, 2))
    parser.add_argument("--random-samples", type=int, default=250)
    parser.add_argument("--random-seed", type=int, default=20260804)
    parser.add_argument("--keep-closest", type=int, default=5)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.support_max < 0 or args.random_samples < 0 or args.keep_closest < 1:
        parser.error(
            "support/random bounds must be nonnegative and --keep-closest positive"
        )

    variables, k12 = build_k12()
    nonlinear = [
        sp.expand(component - variable)
        for component, variable in zip(k12, variables)
    ]
    k12_sparse = [as_sparse_rational(component, variables) for component in k12]
    nonlinear_sparse = [
        as_sparse_rational(component, variables) for component in nonlinear
    ]
    collision_image = (
        sp.Integer(0),
        sp.Integer(0),
        -sp.Rational(1, 4),
        *([sp.Integer(0)] * 9),
    )
    families = {
        family.pivot: family
        for family in linear_graph_families(k12, variables)
        if family.pivot in SELECTED_COMPONENTS
    }

    records = []
    cross_prime_survivors = []
    total_evaluations = 0
    for pivot in SELECTED_COMPONENTS:
        family = families[pivot]
        points = parameter_points(
            len(family.parameters),
            args.support_max,
            args.values,
            args.random_samples,
            args.random_seed + pivot,
        )
        prime_records = []
        survivor_sets = []
        for prime in args.primes:
            audits = []
            survivor_points = set()
            signature_counts: Counter[tuple[int, int, int]] = Counter()
            cubic_graph_count = 0
            for point in points:
                audit = audit_point(
                    family,
                    point,
                    prime,
                    k12_sparse,
                    nonlinear_sparse,
                    collision_image,
                )
                total_evaluations += 1
                if audit["graph_degree"] == 3:
                    cubic_graph_count += 1
                selected_screen = audit["screened_components"][0]
                signature = (
                    int(selected_screen["bilinear_rank"]),
                    int(selected_screen["cubic_quotient_rank"]),
                    0 if selected_screen["cubic_solvable"] else 1,
                )
                signature_counts[signature] += 1
                if audit["all_bad_components_cubically_solvable"]:
                    survivor_points.add(point)
                    audits.append(audit)
                elif audit["graph_degree"] == 3:
                    audits.append(audit)
            audits.sort(key=lambda audit: (audit["score"], audit["parameter_point"]))
            survivor_sets.append(survivor_points)
            prime_records.append(
                {
                    "prime": prime,
                    "point_count": len(points),
                    "genuinely_cubic_graph_count": cubic_graph_count,
                    "selected_component_rank_signatures": [
                        {
                            "bilinear_rank": signature[0],
                            "cubic_quotient_rank": signature[1],
                            "augmented_rank_increment": signature[2],
                            "point_count": count,
                        }
                        for signature, count in sorted(signature_counts.items())
                    ],
                    "survivor_count": len(survivor_points),
                    "survivors": [
                        audit
                        for audit in audits
                        if audit["all_bad_components_cubically_solvable"]
                    ],
                    "closest_cubic_graph_points": [
                        audit
                        for audit in audits
                        if audit["graph_degree"] == 3
                        and not audit["all_bad_components_cubically_solvable"]
                    ][: args.keep_closest],
                }
            )
        common_survivors = set.intersection(*survivor_sets)
        cross_prime_survivors.extend(
            {
                "source_pivot": pivot + 1,
                "parameter_point": list(point),
            }
            for point in sorted(common_survivors)
        )
        records.append(
            {
                "source_pivot": pivot + 1,
                "parameter_count": len(family.parameters),
                "normalized_linear_target_coefficients": [
                    str(coefficient) for coefficient in family.coefficients
                ],
                "selected_retained_component": SELECTED_COMPONENTS[pivot],
                "point_count_per_prime": len(points),
                "prime_records": prime_records,
                "cross_prime_survivor_count": len(common_survivors),
            }
        )
        print(
            f"pivot z{pivot + 1}: {len(points)} points x "
            f"{len(args.primes)} primes; cross-prime survivors="
            f"{len(common_survivors)}"
        )

    artifact = {
        "format": "k12-schur-cubic-completion-modular-search-v1",
        "status": "bounded finite-field discovery experiment; not a proof",
        "primes": list(args.primes),
        "source_pivots": [pivot + 1 for pivot in SELECTED_COMPONENTS],
        "support_maximum": args.support_max,
        "nonzero_integer_values": list(args.values),
        "random_samples_per_parameter_count": args.random_samples,
        "random_seed": args.random_seed,
        "desired_source_degree": DESIRED_DEGREE,
        "target_completion_degree": 3,
        "target_completion_basis": (
            "all 10 linear, 55 quadratic, and 220 cubic monomials in the "
            "other retained raw outputs"
        ),
        "schur_strategy": (
            "build the exact-theorem bilinear column space first, then add "
            "cubic columns in its finite-field quotient"
        ),
        "total_modular_point_evaluations": total_evaluations,
        "families": records,
        "cross_prime_survivors": cross_prime_survivors,
        "scope": (
            "The search covers only the stated parameter points in the five "
            "constant-bilinear-rank normalized linear graph families. A "
            "survivor requires reconstruction over Q and a complete target-"
            "automorphism, determinant, degree, and collision audit. Absence "
            "of survivors is not a characteristic-zero obstruction."
        ),
    }
    output = args.output
    if not output.is_absolute():
        output = ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2) + "\n")
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    try:
        displayed_output = output.relative_to(ROOT)
    except ValueError:
        displayed_output = output
    print(f"PASS wrote {displayed_output}")
    print(f"SHA256 {digest}")
    if cross_prime_survivors:
        print("DISCOVERY: cross-prime survivors require exact reconstruction")
    else:
        print("NO CROSS-PRIME SURVIVOR IN THE DECLARED BOUNDED SEARCH")


if __name__ == "__main__":
    main()
