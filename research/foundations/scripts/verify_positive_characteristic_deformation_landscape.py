#!/usr/bin/env python3
"""Dependency-free regressions for the positive-characteristic landscape.

The proofs and the precise scope are in
extended-geometry/POSITIVE_CHARACTERISTIC_DEFORMATION_LANDSCAPE.md.
The bounded checks here cover normalized-seed dimensions, p-typical Hasse
recovery, the AGL(1,p) Frobenius-monodromy calibration, and the first
tame/wild multiplicity ledger.
"""

from __future__ import annotations

from math import comb, factorial


PRIMES = (2, 3, 5, 7, 11, 13, 17, 19)
DEGREE_BOUND = 40


def matrix_rank_mod_p(matrix: list[list[int]], prime: int) -> int:
    work = [[entry % prime for entry in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next(
            (
                row
                for row in range(rank, len(work))
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [(entry * inverse) % prime for entry in work[rank]]
        for row in range(len(work)):
            if row == rank:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(work[row], work[rank])
                ]
        rank += 1
    return rank


def p_powers_at_most(bound: int, prime: int) -> tuple[int, ...]:
    result = []
    order = 1
    while order <= bound:
        result.append(order)
        order *= prime
    return tuple(result)


def base_p_digit(value: int, place: int, prime: int) -> int:
    return value // place % prime


def verify_normalized_seed_dimensions() -> None:
    """The equations H(1)=0 and H'(1)=-1 have rank two."""
    for prime in PRIMES:
        for degree in range(3, DEGREE_BOUND + 1):
            exponents = range(2, degree + 1)
            coefficient_matrix = [
                [1 for _ in exponents],
                [exponent for exponent in exponents],
            ]
            assert matrix_rank_mod_p(coefficient_matrix, prime) == 2
            assert len(tuple(exponents)) - 2 == degree - 3


def verify_p_typical_recovery() -> None:
    """Each monomial is recovered from one nonzero base-p digit."""
    for prime in PRIMES:
        for degree in range(1, DEGREE_BOUND + 1):
            channels = p_powers_at_most(degree, prime)
            for exponent in range(1, degree + 1):
                detecting = [
                    order
                    for order in channels
                    if comb(exponent, order) % prime
                ]
                assert detecting
                chosen = next(
                    order
                    for order in channels
                    if base_p_digit(exponent, order, prime)
                )
                digit = base_p_digit(exponent, chosen, prime)
                assert comb(exponent, chosen) % prime == digit
                assert digit * pow(digit, -1, prime) % prime == 1

            for forced in channels:
                # W^(p^i) is invisible to every other positive Hasse order.
                for order in range(1, degree + 1):
                    observed = comb(forced, order) % prime
                    assert (observed != 0) == (order == forced)
                # Its ordinary derivative coefficient is p^i=0 mod p.
                assert forced % prime == 0 or forced == 1


def prime_factors(value: int) -> tuple[int, ...]:
    factors = []
    candidate = 2
    remaining = value
    while candidate * candidate <= remaining:
        if remaining % candidate == 0:
            factors.append(candidate)
            while remaining % candidate == 0:
                remaining //= candidate
        candidate += 1
    if remaining > 1:
        factors.append(remaining)
    return tuple(factors)


def primitive_root(prime: int) -> int:
    if prime == 2:
        return 1
    factors = prime_factors(prime - 1)
    return next(
        candidate
        for candidate in range(2, prime)
        if all(
            pow(candidate, (prime - 1) // factor, prime) != 1
            for factor in factors
        )
    )


Permutation = tuple[int, ...]


def compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[value]] for value in range(len(left)))


def generate_group(generators: tuple[Permutation, ...]) -> set[Permutation]:
    degree = len(generators[0])
    identity = tuple(range(degree))
    group = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            product = compose(generator, current)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return group


def verify_affine_monodromy_groups() -> None:
    for prime in PRIMES:
        multiplier = primitive_root(prime)
        translation = tuple((value + 1) % prime for value in range(prime))
        scaling = tuple(multiplier * value % prime for value in range(prime))
        group = generate_group((translation, scaling))

        assert len(group) == prime * (prime - 1)
        assert len(group) == factorial(prime) if prime in (2, 3) else (
            len(group) < factorial(prime)
        )

        # AGL(1,p) is sharply two-transitive.
        images = {
            (permutation[0], permutation[1])
            for permutation in group
        }
        assert len(images) == prime * (prime - 1)
        assert all(left != right for left, right in images)

        # Freshman's-dream covariance behind W=alpha*Y.
        for coefficient in range(prime):
            assert pow(coefficient, prime, prime) == coefficient
            for value in range(prime):
                assert (
                    pow(value + coefficient, prime, prime)
                    - (value + coefficient)
                ) % prime == (pow(value, prime, prime) - value) % prime


def verify_missing_prime_phase_rows() -> None:
    for prime in PRIMES:
        rows = []
        for multiplicity in range(2, 2 * prime + 3):
            rows.append(
                {
                    "e": multiplicity,
                    "wild": multiplicity % prime == 0,
                }
            )
        assert any(row["wild"] for row in rows)
        assert all(
            row["wild"] == (row["e"] % prime == 0)
            for row in rows
        )
        assert rows[prime - 2] == {"e": prime, "wild": True}
        assert rows[prime - 1] == {"e": prime + 1, "wild": False}


def main() -> None:
    verify_normalized_seed_dimensions()
    verify_p_typical_recovery()
    verify_affine_monodromy_groups()
    verify_missing_prime_phase_rows()
    print("Positive-characteristic deformation landscape audit: PASS")
    print(
        "checked normalized dimensions and p-typical channels through "
        f"degree {DEGREE_BOUND}"
    )
    print(
        "checked AGL(1,p) monodromy and tame/wild rows for primes "
        + ", ".join(map(str, PRIMES))
    )


if __name__ == "__main__":
    main()
