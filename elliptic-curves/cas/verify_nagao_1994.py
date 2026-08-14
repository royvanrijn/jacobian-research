#!/usr/bin/env python3
"""Replay Nagao's 1994 rank-13 family and published rank-21 specialization.

The output records exact family identities, exact point membership, current
PARI conductor/local data, and numerical canonical-height determinants.  The
rank lower bounds remain citations to Nagao's paper: this verifier does not
turn a floating determinant into a repository-local exact independence proof.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
from typing import Any, Sequence

from fractions import Fraction

from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK13_CONSTRUCTION,
    RANK13_CONDUCTOR_FACTORIZATION,
    RANK13_PUBLISHED_CONDUCTOR,
    RANK13_PUBLISHED_MODEL,
    RANK13_PUBLISHED_POINTS,
    RANK13_ROOTS,
    RANK21_CONSTRUCTION,
    RANK21_CONDUCTOR_FACTORIZATION,
    RANK21_CONSTRUCTOR_PARAMETER,
    RANK21_PUBLISHED_CONDUCTOR,
    RANK21_PUBLISHED_MODEL,
    RANK21_PUBLISHED_PARAMETER,
    RANK21_PUBLISHED_POINTS,
    RANK21_ROOTS,
    even_discriminant_polynomial,
    factorization_product,
    point_on_extended_weierstrass,
    polynomial_content,
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    rank13_base_changed_discriminant_numerator,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
    rank13_extra_point,
    rank13_known_quartic_points,
    rank13_leading_square,
    rank13_published_quartic_coefficients,
    rank21_short_jacobian_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version


Q = Fraction
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/verify_nagao_1994.py"
)


def _gp_polynomial(coefficients: Sequence[Fraction], variable: str) -> str:
    terms = []
    for degree, coefficient in enumerate(coefficients):
        coefficient = Q(coefficient)
        if coefficient:
            terms.append(
                f"({coefficient.numerator}/{coefficient.denominator})*"
                f"{variable}^{degree}"
            )
    return "+".join(terms) if terms else "0"


def pari_polynomial_is_irreducible(
    coefficients: Sequence[Fraction], *, timeout: float
) -> bool:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    polynomial = _gp_polynomial(coefficients, "x")
    program = f"x='x;P={polynomial};print(polisirreducible(P));quit\n"
    result = subprocess.run(
        [executable, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI polynomial check failed: {result.stderr.strip()}")
    return result.stdout.strip() == "1"


def pari_all_prime(values: Sequence[int], *, timeout: float) -> bool:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    vector = ",".join(str(value) for value in values)
    program = f"print(vecmin(apply(isprime,[{vector}])));quit\n"
    result = subprocess.run(
        [executable, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI primality check failed: {result.stderr.strip()}")
    return result.stdout.strip() == "1"


def polynomial_geometry(timeout: float) -> dict[str, Any]:
    rank13_discriminant = RANK13_CONSTRUCTION.primitive_discriminant_polynomial
    rank21_discriminant = RANK21_CONSTRUCTION.primitive_discriminant_polynomial
    rank13_u_polynomial = even_discriminant_polynomial(RANK13_CONSTRUCTION)
    rank21_u_polynomial = even_discriminant_polynomial(RANK21_CONSTRUCTION)
    base_changed = rank13_base_changed_discriminant_numerator()
    return {
        "rank13": {
            "primitive_quartic_fixed_square_content": str(
                RANK13_CONSTRUCTION.quartic_content
            ),
            "primitive_discriminant_degree_in_T": len(rank13_discriminant) - 1,
            "primitive_discriminant_even": all(
                not rank13_discriminant[index]
                for index in range(1, len(rank13_discriminant), 2)
            ),
            "primitive_discriminant_coefficient_content": str(
                polynomial_content(rank13_discriminant)
            ),
            "degree_in_U_equals_T_squared": len(rank13_u_polynomial) - 1,
            "pari_irreducible_in_T_over_Q": pari_polynomial_is_irreducible(
                rank13_discriminant, timeout=timeout
            ),
            "pari_irreducible_in_U_over_Q": pari_polynomial_is_irreducible(
                rank13_u_polynomial, timeout=timeout
            ),
            "base_changed_cleared_discriminant_degree_in_u": len(base_changed) - 1,
            "base_changed_pari_irreducible_over_Q": pari_polynomial_is_irreducible(
                base_changed, timeout=timeout
            ),
        },
        "rank21": {
            "primitive_quartic_fixed_square_content": str(
                RANK21_CONSTRUCTION.quartic_content
            ),
            "primitive_discriminant_degree_in_T": len(rank21_discriminant) - 1,
            "primitive_discriminant_even": all(
                not rank21_discriminant[index]
                for index in range(1, len(rank21_discriminant), 2)
            ),
            "primitive_discriminant_coefficient_content": str(
                polynomial_content(rank21_discriminant)
            ),
            "degree_in_U_equals_T_squared": len(rank21_u_polynomial) - 1,
            "pari_irreducible_in_T_over_Q": pari_polynomial_is_irreducible(
                rank21_discriminant, timeout=timeout
            ),
            "pari_irreducible_in_U_over_Q": pari_polynomial_is_irreducible(
                rank21_u_polynomial, timeout=timeout
            ),
        },
        "engineering_interpretation": (
            "Each original family has one irreducible nonconstant degree-20 "
            "discriminant factor (degree 10 in T^2), so p-adic power forcing "
            "targets one binary form.  The rank-13 quadratic base change raises "
            "the cleared discriminant degree to 40."
        ),
    }


def _factorization_json(
    factorization: Sequence[tuple[int, int]], *, timeout: float
) -> dict[str, Any]:
    factors = [prime for prime, _ in factorization]
    return {
        "factors": [
            {"prime": str(prime), "exponent": exponent}
            for prime, exponent in factorization
        ],
        "all_factors_pari_isprime": pari_all_prime(factors, timeout=timeout),
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_nagao_1994.json",
    )
    parser.add_argument("--timeout", type=float, default=120.0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if RANK13_CONSTRUCTION.quartic_condition != 0:
        raise AssertionError("Nagao's rank-13 tuple is not a Mestre quartic")
    if RANK21_CONSTRUCTION.quartic_condition != 0:
        raise AssertionError("Nagao's rank-21 tuple is not a Mestre quartic")

    # Pin the paper's printed quartic and thirteenth section at several exact
    # values, then verify all thirteen sections and their covariant images.
    for parameter in (Q(1), Q(2), Q(7, 3), Q(-5, 2)):
        if primitive_quartic_coefficients(
            RANK13_CONSTRUCTION, parameter
        ) != rank13_published_quartic_coefficients(parameter):
            raise AssertionError("Nagao's printed rank-13 quartic changed")
        rank13_extra_point(parameter)
    rank13_quartic_points = rank13_known_quartic_points(rank13_base_parameter(Q(1)))
    rank13_jacobian_points = tuple(
        quartic_point_to_short_jacobian(
            RANK13_CONSTRUCTION, rank13_base_parameter(Q(1)), point
        )
        for point in rank13_quartic_points
    )
    if len(rank13_jacobian_points) != 13:
        raise AssertionError("the rank-13 family did not expose thirteen sections")

    rank13_leading = primitive_quartic_coefficients(
        RANK13_CONSTRUCTION, rank13_base_parameter(Q(1))
    )[-1]
    if rank13_leading_square(Q(1)) ** 2 != rank13_leading:
        raise AssertionError("the quadratic base change did not split infinity")

    rank13_family_curve = minimal_curve_data(
        rank13_base_changed_short_jacobian_coefficients(Q(1)),
        timeout=args.timeout,
    )
    if tuple(rank13_family_curve["minimal_model"]) != RANK13_PUBLISHED_MODEL:
        raise AssertionError("the rank-13 u=1 family replay missed the printed model")
    if rank13_family_curve["conductor"] != RANK13_PUBLISHED_CONDUCTOR:
        raise AssertionError("the rank-13 u=1 conductor changed")
    if not all(
        point_on_extended_weierstrass(RANK13_PUBLISHED_MODEL, point)
        for point in RANK13_PUBLISHED_POINTS
    ):
        raise AssertionError("a printed rank-13 point failed exactly")
    if factorization_product(RANK13_CONDUCTOR_FACTORIZATION) != RANK13_PUBLISHED_CONDUCTOR:
        raise AssertionError("the rank-13 conductor factorization changed")
    rank13_curve = minimal_curve_data(
        tuple(Q(value) for value in RANK13_PUBLISHED_MODEL),
        timeout=args.timeout,
        known_points=RANK13_PUBLISHED_POINTS,
        local_primes=tuple(prime for prime, _ in RANK13_CONDUCTOR_FACTORIZATION),
        stack_bytes=512_000_000,
    )

    rank21_quartic_points = primitive_visible_points(
        RANK21_CONSTRUCTION, RANK21_CONSTRUCTOR_PARAMETER
    )
    rank21_jacobian_points = tuple(
        quartic_point_to_short_jacobian(
            RANK21_CONSTRUCTION, RANK21_CONSTRUCTOR_PARAMETER, point
        )
        for point in rank21_quartic_points
    )
    if len(rank21_jacobian_points) != 12:
        raise AssertionError("the rank-21 tuple did not expose twelve Mestre sections")
    rank21_family_curve = minimal_curve_data(
        rank21_short_jacobian_coefficients(), timeout=args.timeout
    )
    if tuple(rank21_family_curve["minimal_model"]) != RANK21_PUBLISHED_MODEL:
        raise AssertionError("the rank-21 specialization missed the printed model")
    if rank21_family_curve["conductor"] != RANK21_PUBLISHED_CONDUCTOR:
        raise AssertionError("the rank-21 conductor changed")
    if not all(
        point_on_extended_weierstrass(RANK21_PUBLISHED_MODEL, point)
        for point in RANK21_PUBLISHED_POINTS
    ):
        raise AssertionError("a printed rank-21 point failed exactly")
    if factorization_product(RANK21_CONDUCTOR_FACTORIZATION) != RANK21_PUBLISHED_CONDUCTOR:
        raise AssertionError("the rank-21 conductor factorization changed")
    rank21_curve = minimal_curve_data(
        tuple(Q(value) for value in RANK21_PUBLISHED_MODEL),
        timeout=args.timeout,
        known_points=RANK21_PUBLISHED_POINTS,
        local_primes=tuple(prime for prime, _ in RANK21_CONDUCTOR_FACTORIZATION),
        stack_bytes=512_000_000,
    )

    geometry = polynomial_geometry(args.timeout)
    artifact = {
        "schema_version": 1,
        "status": (
            "verified computation: exact Mestre identities, family-to-minimal-model "
            "replay, conductor/local data, and exact printed-point membership; "
            "Nagao's independence results are cited and the height determinants "
            "are numerical replays, not repository-local exact certificates"
        ),
        "primary_source": PRIMARY_SOURCE,
        "family_geometry": geometry,
        "rank13_base_changed_u1": {
            "root_tuple_as_printed": list(RANK13_ROOTS),
            "base_change_parameter_u": "1",
            "resulting_T": str(rank13_base_parameter(Q(1))),
            "leading_coefficient_square_root": str(rank13_leading_square(Q(1))),
            "exact_quartic_sections_checked": len(rank13_quartic_points),
            "exact_short_jacobian_images_checked": len(rank13_jacobian_points),
            "family_minimal_model_matches_paper": True,
            "curve": rank13_curve,
            "conductor_factorization": _factorization_json(
                RANK13_CONDUCTOR_FACTORIZATION, timeout=args.timeout
            ),
            "published_rank_lower_bound": 13,
            "published_points_replayed": len(RANK13_PUBLISHED_POINTS),
            "repository_local_exact_independence_certificate": None,
            "strict_rank21_logN_target_met": False,
        },
        "rank21_record": {
            "root_tuple_as_printed": list(RANK21_ROOTS),
            "published_parameter": str(RANK21_PUBLISHED_PARAMETER),
            "constructor_parameter_reproducing_printed_model": str(
                RANK21_CONSTRUCTOR_PARAMETER
            ),
            "parameter_normalization_note": (
                "The printed model is recovered at constructor T=2*(14721/376); "
                "both values are pinned explicitly."
            ),
            "exact_visible_quartic_points_checked": len(rank21_quartic_points),
            "exact_short_jacobian_images_checked": len(rank21_jacobian_points),
            "family_minimal_model_matches_paper": True,
            "curve": rank21_curve,
            "conductor_factorization": _factorization_json(
                RANK21_CONDUCTOR_FACTORIZATION, timeout=args.timeout
            ),
            "published_rank_lower_bound": 21,
            "published_points_replayed": len(RANK21_PUBLISHED_POINTS),
            "repository_local_exact_independence_certificate": None,
            "strict_rank21_logN_target_met": False,
        },
        "strict_target": {"rank_at_least": 21, "log_conductor_less_than": "182.72"},
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
        },
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reproducing_command": REPRODUCING_COMMAND,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(
        "rank-13 u=1: "
        f"log(N)={rank13_curve['log_conductor']}, "
        f"points={rank13_curve['supplied_points']['on_curve_count']}"
    )
    print(
        "rank-21 record: "
        f"log(N)={rank21_curve['log_conductor']}, "
        f"points={rank21_curve['supplied_points']['on_curve_count']}"
    )


if __name__ == "__main__":
    main()

