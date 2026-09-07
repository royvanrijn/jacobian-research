#!/usr/bin/env sage
"""Certify finite local facts about BNF-free cubic 2-covers by reduction.

For a cover output by ``build_bnf_free_two_covers.py``, a smooth projective
``F_p`` point on its two defining quadrics lifts to ``Q_p``.  At singular
reduction points the script additionally applies a two-variable valuation
Hensel criterion and, when necessary, exhaustively lifts the finite set of
primitive projective residue classes to a requested precision.  Conversely,
an empty reduction at any precision proves there is no ``Q_p`` point.  An
untested large prime or the real place is explicitly left inconclusive.

It is therefore a useful local-obstruction/witness stage before a residual
class is admitted to a full local-Selmer calculation.  It is not a replacement
for all local conditions or a Cassels--Tate computation.
"""

from __future__ import annotations

import argparse
from functools import reduce
from itertools import product
from math import gcd
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, lcm

from build_bnf_free_two_covers import (
    multiply_mod_cubic,
    rational,
    verify_rational_cover_witness,
)


PROTOCOL = "BNFFREECOVERLOCAL"
INPUT_SCHEMA = "elliptic-curves.bnf-free-2cover-equations.v1"


def parse_primes(text: str) -> list[int]:
    try:
        values = sorted({int(part.strip()) for part in text.split(",") if part.strip()})
    except ValueError as exc:
        raise ValueError("--primes must be comma-separated rational primes") from exc
    if not values:
        raise ValueError("--primes must not be empty")
    if any(value < 2 or not ZZ(value).is_prime(proof=False) for value in values):
        raise ValueError("--primes must contain rational primes")
    return values


def cover_quadrics(alpha, coefficients, ring):
    u, v, w, z = ring.gens()
    beta = [u, v, w]
    constant, theta_coefficient, theta_squared_coefficient = multiply_mod_cubic(
        alpha,
        multiply_mod_cubic(beta, beta, coefficients),
        coefficients,
    )
    del constant
    return theta_coefficient + z**2, theta_squared_coefficient


def primitive_integral(polynomial):
    """Return a primitive integral polynomial defining the same quadric."""
    coefficients = polynomial.coefficients()
    denominator = lcm([coefficient.denominator() for coefficient in coefficients])
    variables = polynomial.parent().gens()
    integral = polynomial.parent()(sum(
        ZZ(coefficient * denominator)
        * reduce(
            lambda left, right: left * right,
            (variable**power for variable, power in zip(variables, monomial)),
            polynomial.parent().base_ring()(1),
        )
        for monomial, coefficient in polynomial.dict().items()
    ))
    content = reduce(gcd, (abs(int(value)) for value in integral.coefficients() if value), 0)
    if content == 0:
        raise ValueError("cover equation is identically zero")
    return integral / content


def projective_points_mod_prime(first, second):
    """Enumerate P^3(F_p), retaining singular points if no smooth point occurs."""
    field = first.base_ring()
    variables = first.parent().gens()
    jacobian = [
        [first.derivative(variable) for variable in variables],
        [second.derivative(variable) for variable in variables],
    ]
    point_count = 0
    singular_points = []
    for leading in range(4):
        prefix = [field(0)] * leading + [field(1)]
        for tail in product(field, repeat=3 - leading):
            point = prefix + list(tail)
            if first(*point) != 0 or second(*point) != 0:
                continue
            point_count += 1
            minors = [
                jacobian[0][left](*point) * jacobian[1][right](*point)
                - jacobian[0][right](*point) * jacobian[1][left](*point)
                for left in range(4)
                for right in range(left + 1, 4)
            ]
            if any(minor != 0 for minor in minors):
                return point_count, [int(value) for value in point], ()
            singular_points.append((leading, tuple(int(value) for value in point)))
    return point_count, None, tuple(singular_points)


def valuation(value: int, prime: int):
    if value == 0:
        return None
    value = abs(value)
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


def integer_evaluation(polynomial, point) -> int:
    return int(polynomial(*point))


def hensel_witness(first, second, point, prime):
    """Return a valuation-Hensel certificate for a projective affine chart."""
    variables = first.parent().gens()
    jacobian = [
        [first.derivative(variable) for variable in variables],
        [second.derivative(variable) for variable in variables],
    ]
    values = [integer_evaluation(first, point), integer_evaluation(second, point)]
    value_valuations = [valuation(value, prime) for value in values]
    best = None
    for left in range(4):
        for right in range(left + 1, 4):
            determinant = (
                integer_evaluation(jacobian[0][left], point)
                * integer_evaluation(jacobian[1][right], point)
                - integer_evaluation(jacobian[0][right], point)
                * integer_evaluation(jacobian[1][left], point)
            )
            determinant_valuation = valuation(determinant, prime)
            if determinant_valuation is None:
                continue
            if best is None or determinant_valuation < best[0]:
                best = (determinant_valuation, left, right)
    if best is None:
        return None
    determinant_valuation, left, right = best
    if all(
        item is None or item >= 2 * determinant_valuation + 1
        for item in value_valuations
    ):
        return {
            "free_coordinate_indices": [left, right],
            "jacobian_minor_valuation": determinant_valuation,
            "equation_valuations": [
                "infinity" if item is None else item for item in value_valuations
            ],
        }
    return None


def lift_singular_points(
    first,
    second,
    prime: int,
    singular_points,
    max_precision: int,
    max_states: int,
):
    """Certify a Hensel lift or an empty finite p-adic reduction tree."""
    states = list(singular_points)
    for precision in range(1, max_precision + 1):
        for leading, point in states:
            certificate = hensel_witness(first, second, point, prime)
            if certificate is not None:
                return {
                    "classification": "PROVED_QP_POINT_BY_SINGULAR_HENSEL_LIFT",
                    "lift_precision": precision,
                    "projective_lift_witness": list(point),
                    "hensel_certificate": certificate,
                }
        if precision == max_precision:
            return {
                "classification": "INCONCLUSIVE_SINGULAR_LIFT_PRECISION",
                "lift_precision": precision,
                "surviving_projective_lifts": len(states),
            }
        modulus = prime ** (precision + 1)
        increment = prime**precision
        next_states = []
        for leading, point in states:
            free_indices = [index for index in range(4) if index != leading]
            for digits in product(range(prime), repeat=3):
                lifted = list(point)
                for index, digit in zip(free_indices, digits):
                    lifted[index] += increment * digit
                if (
                    integer_evaluation(first, lifted) % modulus == 0
                    and integer_evaluation(second, lifted) % modulus == 0
                ):
                    next_states.append((leading, tuple(lifted)))
                    if len(next_states) > max_states:
                        return {
                            "classification": "INCONCLUSIVE_SINGULAR_LIFT_STATE_CAP",
                            "lift_precision": precision + 1,
                            "surviving_projective_lifts_lower_bound": len(next_states),
                        }
        if not next_states:
            return {
                "classification": "PROVED_NO_QP_POINT_BY_EMPTY_MOD_P_TO_N_REDUCTION",
                "empty_lift_precision": precision + 1,
            }
        states = next_states
    raise AssertionError("unreachable lift precision state")


def audit_cover(
    alpha,
    coefficients,
    primes,
    max_enumeration_prime,
    max_lift_precision,
    max_lift_states,
    rational_cover_witness=None,
):
    ring = PolynomialRing(QQ, names=("u", "v", "w", "z"))
    first, second = cover_quadrics(alpha, coefficients, ring)
    first = primitive_integral(first)
    second = primitive_integral(second)
    if rational_cover_witness is not None:
        affine_x = verify_rational_cover_witness(
            alpha, coefficients, rational_cover_witness, ring
        )
        return {
            "primitive_integral_quadrics": [str(first), str(second)],
            "global_rational_cover_witness": [
                str(value) for value in rational_cover_witness
            ],
            "global_rational_cover_witness_verified": True,
            "rational_witness_affine_x": (
                None if affine_x is None else str(affine_x)
            ),
            "finite_places": [
                {
                    "rational_prime": prime,
                    "classification": "PROVED_QP_POINT_BY_GLOBAL_Q_WITNESS",
                }
                for prime in primes
            ],
        }
    records = []
    for prime in primes:
        if prime > max_enumeration_prime:
            records.append(
                {
                    "rational_prime": prime,
                    "classification": "SKIPPED_PRIME_EXCEEDS_ENUMERATION_LIMIT",
                }
            )
            continue
        finite_field = GF(prime)
        reduction_ring = PolynomialRing(finite_field, names=("u", "v", "w", "z"))
        first_mod = reduction_ring(first)
        second_mod = reduction_ring(second)
        point_count, smooth_witness, singular_points = projective_points_mod_prime(
            first_mod, second_mod
        )
        if smooth_witness is not None:
            classification = "PROVED_QP_POINT_BY_SMOOTH_FP_LIFT"
        elif point_count == 0:
            classification = "PROVED_NO_QP_POINT_BY_EMPTY_FP_REDUCTION"
        else:
            lift_result = lift_singular_points(
                first,
                second,
                prime,
                singular_points,
                max_lift_precision,
                max_lift_states,
            )
            classification = lift_result.pop("classification")
        record = {
            "rational_prime": prime,
            "classification": classification,
        }
        if smooth_witness is None:
            record["projective_fp_point_count"] = point_count
        else:
            record["projective_fp_point_count_lower_bound"] = point_count
        if smooth_witness is not None:
            record["smooth_projective_fp_witness"] = smooth_witness
        elif point_count:
            record.update(lift_result)
        records.append(record)
    return {
        "primitive_integral_quadrics": [str(first), str(second)],
        "finite_places": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--covers", type=Path, required=True)
    parser.add_argument("--primes", required=True)
    parser.add_argument("--max-enumeration-prime", type=int, default=251)
    parser.add_argument("--max-lift-precision", type=int, default=12)
    parser.add_argument("--max-lift-states", type=int, default=10000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.max_enumeration_prime < 2:
        raise ValueError("--max-enumeration-prime must be at least 2")
    if args.max_lift_precision < 1 or args.max_lift_states < 1:
        raise ValueError("lifting limits must be positive")

    record = json.loads(args.covers.read_text())
    if not isinstance(record, dict) or record.get("schema") != INPUT_SCHEMA:
        raise ValueError("expected BNF-free two-cover equations")
    coefficients = [rational(value) for value in record["field_polynomial_ascending"]]
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("the cover field must be a monic cubic")
    primes = parse_primes(args.primes)

    covers = []
    for cover in record.get("covers", []):
        alpha = [rational(value) for value in cover["alpha_coefficients"]]
        if len(alpha) != 3:
            raise ValueError("a cubic cover alpha needs three power-basis coordinates")
        result = audit_cover(
            alpha,
            coefficients,
            primes,
            args.max_enumeration_prime,
            args.max_lift_precision,
            args.max_lift_states,
            (
                [rational(value) for value in cover["rational_cover_witness"]]
                if cover.get("rational_cover_witness") is not None
                else None
            ),
        )
        result.update(
            {
                "label": str(cover["label"]),
                "alpha_coefficients": [str(value) for value in alpha],
            }
        )
        covers.append(result)

    rational_witness_count = sum(
        cover.get("global_rational_cover_witness_verified") is True
        for cover in covers
    )
    output = {
        "protocol": "BNFFREECOVERLOCAL-v1",
        "status": (
            "GLOBAL_RATIONAL_POINT_POSITIVE_CONTROL_AUDIT"
            if covers and rational_witness_count == len(covers)
            else "SELECTED_FINITE_LOCAL_REDUCTION_AUDIT_ONLY"
        ),
        "tested_rational_primes": primes,
        "max_enumeration_prime": args.max_enumeration_prime,
        "max_lift_precision": args.max_lift_precision,
        "max_lift_states": args.max_lift_states,
        "global_rational_witness_cover_count": rational_witness_count,
        "covers": covers,
        "claim_boundary": [
            "A verified global rational cover witness proves local solubility at every place for that cover.",
            "Finite-reduction witnesses and obstructions retain their one-place meanings.",
            "This audit gives no ambient class-group completeness or Selmer upper bound.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    classifications = {
        item["classification"]
        for cover in covers
        for item in cover["finite_places"]
    }
    print(
        f"{PROTOCOL}|stage=complete|covers={len(covers)}|primes={len(primes)}"
        f"|classifications={','.join(sorted(classifications)) or 'none'}"
        f"|status={output['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
