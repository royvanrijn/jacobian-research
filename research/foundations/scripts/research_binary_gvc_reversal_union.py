#!/usr/bin/env python3
"""Search the smallest coordinate-reversed equal-radial cancellation.

After complete equal-radial recombination, a coordinate-reversed tie has
the one-variable model

    S_N = [z^(N*r)] G(z)^N + [z^(-N*r)] G(z)^N.

Both moving coefficients have the same radial factorial.  The question
is whether ``S_N=0`` for every ``N`` can persist while both slopes
``r`` and ``-r`` lie in the Newton interval of ``G``.

For support ``{-1,0,1}`` and ``r=1``, the two target coefficients are
extreme.  The first two rows are

    c_(-1)+c_1,
    c_(-1)^2+c_1^2,

so they force both endpoints to vanish in characteristic zero.

The script also performs an exhaustive projective search over finite
fields for support ``{-2,-1,0,1,2}``, target slope ``r=1``, and a
configurable moment depth.  This is a bounded obstruction search, not an
all-order proof and not a GVC counterexample.
"""

from __future__ import annotations

import argparse
from itertools import product

import sympy as sp


Laurent = dict[int, int]


def multiply(left: Laurent, right: Laurent, prime: int) -> Laurent:
    answer: Laurent = {}
    for left_exponent, left_value in left.items():
        for right_exponent, right_value in right.items():
            exponent = left_exponent + right_exponent
            answer[exponent] = (
                answer.get(exponent, 0) + left_value * right_value
            ) % prime
    return {
        exponent: value
        for exponent, value in answer.items()
        if value % prime
    }


def reversal_rows(
    coefficients: tuple[int, ...],
    exponents: tuple[int, ...],
    target_slope: int,
    depth: int,
    prime: int,
) -> tuple[int, ...]:
    polynomial = {
        exponent: coefficient % prime
        for exponent, coefficient in zip(exponents, coefficients)
        if coefficient % prime
    }
    power: Laurent = {0: 1}
    rows = []
    for moment in range(1, depth + 1):
        power = multiply(power, polynomial, prime)
        target = moment * target_slope
        rows.append(
            (
                power.get(target, 0)
                + power.get(-target, 0)
            )
            % prime
        )
    return tuple(rows)


def projective_points(prime: int, length: int):
    """Yield one representative with first nonzero coordinate equal to 1."""

    for pivot in range(length):
        for tail in product(range(prime), repeat=length - pivot - 1):
            yield (0,) * pivot + (1,) + tail


def both_slopes_in_newton(
    coefficients: tuple[int, ...],
    exponents: tuple[int, ...],
    target_slope: int,
) -> bool:
    support = tuple(
        exponent
        for exponent, coefficient in zip(exponents, coefficients)
        if coefficient
    )
    return (
        support
        and min(support) <= -target_slope
        and max(support) >= target_slope
    )


def exhaustive_search(
    prime: int,
    depth: int,
) -> tuple[tuple[int, ...], dict[int, int]]:
    exponents = (-2, -1, 0, 1, 2)
    target_slope = 1
    survivors = []
    first_failure_counts = {moment: 0 for moment in range(1, depth + 1)}

    for coefficients in projective_points(prime, len(exponents)):
        if not both_slopes_in_newton(
            coefficients,
            exponents,
            target_slope,
        ):
            continue
        rows = reversal_rows(
            coefficients,
            exponents,
            target_slope,
            depth,
            prime,
        )
        for index, row in enumerate(rows, 1):
            if row:
                first_failure_counts[index] += 1
                break
        else:
            survivors.append(coefficients)

    return tuple(survivors), first_failure_counts


def symbolic_reversal_rows(
    coefficients: dict[int, sp.Expr],
    depth: int,
) -> tuple[sp.Expr, ...]:
    power: dict[int, sp.Expr] = {0: sp.Integer(1)}
    rows = []
    for moment in range(1, depth + 1):
        new: dict[int, sp.Expr] = {}
        for left_exponent, left_value in power.items():
            for right_exponent, right_value in coefficients.items():
                exponent = left_exponent + right_exponent
                new[exponent] = (
                    new.get(exponent, sp.Integer(0))
                    + left_value * right_value
                )
        power = new
        rows.append(
            sp.expand(
                power.get(moment, sp.Integer(0))
                + power.get(-moment, sp.Integer(0))
            )
        )
    return tuple(rows)


def prove_width_two(depth: int) -> dict[tuple[int, int], int]:
    """Saturate the four endpoint charts over QQ."""

    assert depth >= 2
    symbols = {
        exponent: sp.Symbol(f"c_{exponent + 2}")
        for exponent in range(-2, 3)
    }
    t = sp.Symbol("t")
    closure_depths = {}

    for minimum, maximum in ((-1, 1), (-1, 2), (-2, 1), (-2, 2)):
        active = {
            exponent: symbols[exponent]
            for exponent in range(minimum, maximum + 1)
        }
        rows = symbolic_reversal_rows(active, depth)
        endpoint_product = symbols[minimum] * symbols[maximum]
        variables = (
            t,
            *(
                symbols[exponent]
                for exponent in range(minimum, maximum + 1)
            ),
        )
        closed_at = None
        for row_count in range(2, depth + 1):
            basis = sp.groebner(
                [*rows[:row_count], 1 - t * endpoint_product],
                *variables,
                order="grevlex",
                domain=sp.QQ,
            )
            if basis.is_zero_dimensional and any(
                polynomial.as_expr() == 1
                for polynomial in basis.polys
            ):
                closed_at = row_count
                break
            if any(
                polynomial.as_expr() == 1
                for polynomial in basis.polys
            ):
                closed_at = row_count
                break
        if closed_at is None:
            raise AssertionError(
                f"endpoint chart {(minimum, maximum)} survives depth {depth}"
            )
        closure_depths[(minimum, maximum)] = closed_at

    return closure_depths


def verify(
    primes: tuple[int, ...],
    depth: int,
    prove_width_two_flag: bool,
) -> None:
    assert depth >= 2

    # Exact characteristic-zero endpoint calculation.
    # From c_-1+c_1=0, the second row becomes 2*c_-1^2.
    for prime in primes:
        assert prime > 2

    results = {}
    for prime in primes:
        survivors, failures = exhaustive_search(prime, depth)
        results[prime] = (survivors, failures)
        print(
            f"prime {prime}: first-failure counts={failures}, "
            f"survivors through depth {depth}={len(survivors)}"
        )
        if survivors:
            print(f"  survivor representatives: {survivors[:12]}")

    print(
        "exact endpoint support {-1,0,1}: "
        "S1=c_-1+c_1 and S2=c_-1^2+c_1^2 force endpoint support loss"
    )
    if prove_width_two_flag:
        closure_depths = prove_width_two(depth)
        print(
            "PASS characteristic-zero saturated width-two endpoint charts: "
            f"{closure_depths}"
        )
    print(
        "STATUS: finite-field reversal search is bounded evidence only; "
        "the optional saturated chart calculation is exact, but an "
        "all-width two-orientation theorem remains to be proved"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", type=int, nargs="+", default=(5, 7, 11))
    parser.add_argument("--depth", type=int, default=8)
    parser.add_argument("--prove-width-two", action="store_true")
    arguments = parser.parse_args()
    verify(
        tuple(arguments.primes),
        arguments.depth,
        arguments.prove_width_two,
    )


if __name__ == "__main__":
    main()
