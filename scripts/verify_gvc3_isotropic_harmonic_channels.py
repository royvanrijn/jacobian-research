#!/usr/bin/env python3
"""Exact obstruction for two/three isotropic harmonic channels below degree 12.

For an even balanced degree d, the coherent harmonic channel attached to an
isotropic covector v and an even positive degree ell is

    rho^((d-ell)/2) * <v,z>^ell.

The checker compiles Reynolds moments from cross-pairing counts.  Modular
Groebner calculations at three primes discover the cutoff on every
three-channel interaction chart; exact Groebner calculations over Q then
certify the same unit ideals.  Collision charts (two proportional isotropic
directions) are treated separately.  No finite-prefix survivor is promoted
to an all-order candidate.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "gvc3_isotropic_harmonic_channels.json"
)

A, B, C, Z = sp.symbols("a b c z")
COEFFICIENTS = (A, B, C)
HARMONIC_DEGREES = (2, 4, 6, 8, 10)
TRIPLES = tuple(itertools.combinations(HARMONIC_DEGREES, 3))
COLLISION_PAIRS = ((0, 1), (0, 2), (1, 2))
PRIMES = (101, 103, 107)


def odd_double_factorial(value: int) -> int:
    answer = 1
    for factor in range(1, value + 1, 2):
        answer *= factor
    return answer


def primitive_chart_polynomial(expression: sp.Expr) -> sp.Expr:
    """Put a=1 and return the primitive integral numerator in b,c."""
    polynomial = sp.Poly(expression.subs(A, 1), B, C, domain=sp.QQ)
    return polynomial.clear_denoms()[1].primitive()[1].as_expr()


def distinct_direction_moment(
    degrees: tuple[int, int, int],
    order: int,
    extra_degrees: tuple[int, int, int] = (0, 0, 0),
) -> sp.Expr:
    """Reynolds moment for three pairwise distinct isotropic directions.

    Pairwise contractions have been absorbed into the three channel
    coefficients.  The optional ``extra_degrees`` inserts a coherent
    multiplier before performing the same Wick pairing count.
    """
    answer = sp.Integer(0)
    for n0 in range(order + 1):
        for n1 in range(order - n0 + 1):
            counts = (n0, n1, order - n0 - n1)
            stubs = tuple(
                degree * count + extra
                for degree, count, extra in zip(
                    degrees, counts, extra_degrees, strict=True
                )
            )
            total = sum(stubs)
            if total % 2 or 2 * max(stubs) > total:
                continue
            edges = (
                (stubs[0] + stubs[1] - stubs[2]) // 2,
                (stubs[0] + stubs[2] - stubs[1]) // 2,
                (stubs[1] + stubs[2] - stubs[0]) // 2,
            )
            multinomial = math.factorial(order)
            for count in counts:
                multinomial //= math.factorial(count)
            numerator = multinomial
            for stub_count in stubs:
                numerator *= math.factorial(stub_count)
            denominator = odd_double_factorial(total + 1)
            for edge_count in edges:
                denominator *= math.factorial(edge_count)
            monomial = sp.prod(
                coefficient**count
                for coefficient, count in zip(COEFFICIENTS, counts, strict=True)
            )
            answer += sp.Rational(numerator, denominator) * monomial
    return sp.factor(answer)


def collision_moment(
    degrees: tuple[int, int, int],
    collision_pair: tuple[int, int],
    order: int,
) -> sp.Expr:
    """Reynolds moment when exactly two channel directions coincide."""
    remaining = next(index for index in range(3) if index not in collision_pair)
    answer = sp.Integer(0)
    for n0 in range(order + 1):
        for n1 in range(order - n0 + 1):
            counts = (n0, n1, order - n0 - n1)
            left_stubs = sum(degrees[index] * counts[index] for index in collision_pair)
            right_stubs = degrees[remaining] * counts[remaining]
            if left_stubs != right_stubs:
                continue
            multinomial = math.factorial(order)
            for count in counts:
                multinomial //= math.factorial(count)
            coefficient = sp.Rational(
                multinomial * math.factorial(left_stubs),
                odd_double_factorial(2 * left_stubs + 1),
            )
            monomial = sp.prod(
                variable**count
                for variable, count in zip(COEFFICIENTS, counts, strict=True)
            )
            answer += coefficient * monomial
    return sp.factor(answer)


def is_unit_groebner(
    equations: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    modulus: int | None,
) -> tuple[bool, str]:
    options = {} if modulus is None else {"modulus": modulus}
    basis = sp.groebner(equations, *variables, order="grevlex", **options)
    unit = len(basis.polys) == 1 and basis.polys[0].as_expr() == 1
    serialized = "\n".join(str(polynomial.as_expr()) for polynomial in basis.polys)
    return unit, hashlib.sha256(serialized.encode()).hexdigest()


def first_unit_cutoff_distinct(
    degrees: tuple[int, int, int], modulus: int | None
) -> tuple[int, str, list[int]]:
    equations = [Z * B * C - 1]
    nonzero_orders: list[int] = []
    for order in range(2, 9):
        moment = distinct_direction_moment(degrees, order)
        if moment:
            equations.append(primitive_chart_polynomial(moment))
            nonzero_orders.append(order)
        unit, basis_hash = is_unit_groebner(equations, (Z, B, C), modulus)
        if unit:
            return order, basis_hash, nonzero_orders
    raise AssertionError(f"distinct-direction chart survived through order 8: {degrees}")


def first_unit_cutoff_collision(
    degrees: tuple[int, int, int],
    collision_pair: tuple[int, int],
    modulus: int | None,
) -> tuple[int, str, list[int]]:
    equations = [Z * A * B - 1]
    nonzero_orders: list[int] = []
    for order in range(2, 9):
        moment = collision_moment(degrees, collision_pair, order)
        if moment:
            polynomial = sp.Poly(moment.subs(C, 1), A, B, domain=sp.QQ)
            equations.append(
                polynomial.clear_denoms()[1].primitive()[1].as_expr()
            )
            nonzero_orders.append(order)
        unit, basis_hash = is_unit_groebner(equations, (Z, A, B), modulus)
        if unit:
            return order, basis_hash, nonzero_orders
    raise AssertionError(
        f"two-direction collision chart survived through order 8: "
        f"{degrees}, {collision_pair}"
    )


def two_channel_return(left: int, right: int) -> dict[str, int]:
    common = math.gcd(left, right)
    left_count = right // common
    right_count = left // common
    order = left_count + right_count
    stubs = left * left_count
    coefficient = sp.Rational(
        math.comb(order, left_count) * math.factorial(stubs),
        odd_double_factorial(2 * stubs + 1),
    )
    assert coefficient != 0
    return {
        "left_count": left_count,
        "right_count": right_count,
        "first_return_order": order,
    }


def near_survivor() -> dict[str, object]:
    """The H2+H4+H6 prefix survivor, killed by its fifth pure moment."""
    degrees = (2, 4, 6)
    c_value = -sp.Rational(143, 60)
    b_relation = 56 * B**2 + 272 * B + 85
    substitutions = {A: 1, C: c_value}
    pure_remainders: dict[str, str] = {}
    for order in (3, 4, 5):
        numerator = sp.together(
            distinct_direction_moment(degrees, order).subs(substitutions)
        ).as_numer_denom()[0]
        pure_remainders[str(order)] = str(sp.factor(sp.rem(numerator, b_relation, B)))
    assert pure_remainders["3"] == "0"
    assert pure_remainders["4"] == "0"
    assert pure_remainders["5"] != "0"
    assert sp.gcd(
        sp.Poly(pure_remainders["5"], B), sp.Poly(b_relation, B)
    ).degree() == 0

    multiplier = distinct_direction_moment(degrees, 2, (2, 0, 0))
    multiplier_value = sp.factor(multiplier.subs(substitutions))
    multiplier_remainder = sp.factor(
        sp.rem(sp.together(multiplier_value).as_numer_denom()[0], b_relation, B)
    )
    assert multiplier_remainder == 8 * B
    return {
        "degrees": list(degrees),
        "normalization": "a=1",
        "c": str(c_value),
        "b_relation": str(b_relation),
        "pure_numerator_remainders_mod_b_relation": pure_remainders,
        "H2_multiplier_order_2_numerator_remainder": str(multiplier_remainder),
        "status": (
            "pure moments 1..4 vanish and the H2 multiplier channel survives "
            "at order 2, but pure moment 5 is nonzero"
        ),
    }


def main() -> None:
    two_channel = {
        f"{left},{right}": two_channel_return(left, right)
        for left, right in itertools.combinations(HARMONIC_DEGREES, 2)
    }

    distinct_results = []
    collision_results = []
    for degrees in TRIPLES:
        modular_cutoffs = {}
        discovered_cutoff = None
        discovered_orders = None
        for prime in PRIMES:
            cutoff, _, modular_orders = first_unit_cutoff_distinct(degrees, prime)
            if discovered_cutoff is None:
                discovered_cutoff = cutoff
                discovered_orders = modular_orders
            assert cutoff == discovered_cutoff
            assert modular_orders == discovered_orders
            modular_cutoffs[str(prime)] = cutoff
        exact_cutoff, exact_hash, nonzero_orders = first_unit_cutoff_distinct(
            degrees, None
        )
        assert exact_cutoff == discovered_cutoff
        assert nonzero_orders == discovered_orders
        distinct_results.append(
            {
                "degrees": list(degrees),
                "nonzero_moment_orders_used": nonzero_orders,
                "exact_Q_unit_cutoff": exact_cutoff,
                "exact_basis_sha256": exact_hash,
                "modular_unit_cutoffs": modular_cutoffs,
            }
        )

        for collision_pair in COLLISION_PAIRS:
            modular_cutoffs = {}
            discovered_cutoff = None
            discovered_orders = None
            for prime in PRIMES:
                cutoff, _, modular_orders = first_unit_cutoff_collision(
                    degrees, collision_pair, prime
                )
                if discovered_cutoff is None:
                    discovered_cutoff = cutoff
                    discovered_orders = modular_orders
                assert cutoff == discovered_cutoff
                assert modular_orders == discovered_orders
                modular_cutoffs[str(prime)] = cutoff
            exact_cutoff, exact_hash, nonzero_orders = first_unit_cutoff_collision(
                degrees, collision_pair, None
            )
            assert exact_cutoff == discovered_cutoff
            assert nonzero_orders == discovered_orders
            collision_results.append(
                {
                    "degrees": list(degrees),
                    "coincident_channel_indices": list(collision_pair),
                    "nonzero_moment_orders_used": nonzero_orders,
                    "exact_Q_unit_cutoff": exact_cutoff,
                    "exact_basis_sha256": exact_hash,
                    "modular_unit_cutoffs": modular_cutoffs,
                }
            )

    assert max(result["exact_Q_unit_cutoff"] for result in distinct_results) == 6
    assert max(result["exact_Q_unit_cutoff"] for result in collision_results) == 8
    assert max(result["first_return_order"] for result in two_channel.values()) == 9

    artifact = {
        "format": "gvc3-isotropic-harmonic-channels-v1",
        "status": "exact characteristic-zero obstruction in the declared coherent-channel architecture",
        "balanced_polynomial_degrees": [4, 6, 8, 10],
        "laplacian_powers": [2, 3, 4, 5],
        "harmonic_decomposition": "P=sum_i c_i*rho^((d-ell_i)/2)*<v_i,z>^ell_i",
        "profile_condition": "each v_i is isotropic and ell_i is positive, even, and distinct",
        "modular_discovery_primes": list(PRIMES),
        "two_channel_first_returns": two_channel,
        "three_distinct_directions": distinct_results,
        "two_distinct_directions": collision_results,
        "one_direction_terminal_reason": (
            "all active coherent states have one positive phase direction; "
            "for a fixed multiplier Q the phase-weight/degree bound kills "
            "Delta^(d*m/2)(Q*P^m) for all sufficiently large m"
        ),
        "first_nine_moment_conclusion": (
            "within the declared two/three-channel architecture, vanishing "
            "of moments 1..9 forces the one-direction terminal stratum"
        ),
        "near_survivor": near_survivor(),
        "candidate_promotion": "none; every multi-direction prefix survivor is eliminated",
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print("PASS two-channel return obstruction for 10 harmonic-degree pairs")
    print("PASS modular discovery and exact-Q promotion for 10 distinct charts")
    print("PASS modular discovery and exact-Q promotion for 30 collision charts")
    print("PASS first nine moments force the one-direction terminal stratum")
    print("PASS H2+H4+H6 near survivor is killed by pure moment five")


if __name__ == "__main__":
    main()
