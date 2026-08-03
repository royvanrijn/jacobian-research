#!/usr/bin/env python3
"""Exact diagonal blocks in the first binary GVC ghost shell.

The centered three-level block is the Cartier coefficient of
``(a+b*z+c*z**2)^p`` at exponent ``p`` after deleting the pure
``b**p`` term.  After removing support factors it is a truncated
Catalan polynomial in ``x=a*c/b**2`` and has the universal terminal
factor ``x-1``.

The two-by-two beta block is the Mirimanoff truncated logarithm in the
cross-ratio ``x``.  It has the support factor ``x`` and the universal
Hall-annihilator factor ``x+1``.  Its characteristic-zero lifts have
the persistent cyclotomic factor ``x**2+x+1``, so one-row cross-prime
avoidance is false.

The script factors both blocks over exact finite fields, verifies the
central-trinomial congruence behind the centered factor, verifies that
the augmented beta and centered endpoint blocks are terminal, and
searches a bounded rational cross-ratio window for points surviving
several good primes.  The exact identities are proofs; the cross-prime
window is only a regression.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
from math import comb, factorial, gcd

import sympy as sp


X = sp.symbols("X")


def primes_through(limit: int) -> tuple[int, ...]:
    return tuple(int(p) for p in sp.primerange(5, limit + 1))


def centered_block(prime: int) -> sp.Poly:
    coefficients = [
        (
            factorial(prime - 1)
            // (
                factorial(selection) ** 2
                * factorial(prime - 2 * selection)
            )
        )
        % prime
        for selection in range(1, (prime - 1) // 2 + 1)
    ]
    return sp.Poly(
        sum(
            coefficient * X**index
            for index, coefficient in enumerate(coefficients)
        ),
        X,
        modulus=prime,
    )


def beta_block(prime: int) -> sp.Poly:
    return sp.Poly(
        sum(
            ((-1) ** selection * pow(selection, -1, prime) % prime)
            * X**selection
            for selection in range(1, prime)
        ),
        X,
        modulus=prime,
    )


def exact_beta_diagonal(prime: int) -> sp.Poly:
    return sp.Poly(
        sp.expand(((1 + X) ** prime - 1 - X**prime) / prime),
        X,
        domain=sp.QQ,
    )


def polynomial_value_mod(
    polynomial: sp.Poly,
    value: Fraction,
    prime: int,
) -> int:
    numerator = value.numerator % prime
    denominator = value.denominator % prime
    assert denominator
    residue = numerator * pow(denominator, -1, prime) % prime
    return int(polynomial.eval(residue)) % prime


def three_level_singleton_ghost_value(
    left_count: int,
    right_count: int,
    value: Fraction,
    prime: int,
) -> int:
    """Evaluate the primitive three-level singleton ghost modulo ``p``.

    The primitive affine move is

        (left_count, -(left_count+right_count), right_count).

    At order ``p`` its complete one-dimensional fibre is parametrized by
    ``t``.  After removing the pure singleton and dividing by ``p``, the
    cross-ratio polynomial has the coefficients used below.
    """

    denominator = value.denominator % prime
    assert denominator
    residue = value.numerator % prime
    residue = residue * pow(denominator, -1, prime) % prime
    step = left_count + right_count
    answer = 0
    for multiplicity in range(1, prime // step + 1):
        coefficient = (
            factorial(prime - 1)
            // (
                factorial(left_count * multiplicity)
                * factorial(right_count * multiplicity)
                * factorial(prime - step * multiplicity)
            )
        )
        answer += coefficient * pow(residue, multiplicity, prime)
    return answer % prime


def three_level_singleton_ghost_polynomial(
    left_count: int,
    right_count: int,
    prime: int,
) -> sp.Poly:
    """Return the complete primitive three-level ghost over F_p."""

    step = left_count + right_count
    terms = []
    for multiplicity in range(1, prime // step + 1):
        denominator = (
            factorial(left_count * multiplicity)
            * factorial(right_count * multiplicity)
            * factorial(prime - step * multiplicity)
        )
        coefficient = factorial(prime - 1) // denominator
        logarithmic_coefficient = (
            (-1) ** (step * multiplicity - 1)
            * comb(step * multiplicity, left_count * multiplicity)
            * pow(step * multiplicity, -1, prime)
        ) % prime
        assert coefficient % prime == logarithmic_coefficient
        terms.append(coefficient * X**multiplicity)

    return sp.Poly(sum(terms), X, modulus=prime)


def verify_three_level_cyclotomic_ghosts(
    prime_limit: int,
    step_limit: int,
    cyclotomic_limit: int,
) -> None:
    """Search for bounded root-of-unity factors surviving every prime.

    If an algebraic root of unity of order n is a ghost zero at every
    unramified prime ideal, its cyclotomic polynomial divides the reduced
    ghost at every rational prime not dividing n.  This gives an exact
    bounded algebraic test without choosing embeddings into finite fields.
    """

    survivor_table: dict[tuple[int, int], tuple[int, ...]] = {}
    divisions = 0
    types_tested = 0

    for left_count in range(1, step_limit + 1):
        for right_count in range(left_count, step_limit + 1):
            if gcd(left_count, right_count) != 1:
                continue
            step = left_count + right_count
            primes = tuple(
                prime
                for prime in primes_through(prime_limit)
                if prime > step
            )
            assert primes
            ghosts = {
                prime: three_level_singleton_ghost_polynomial(
                    left_count,
                    right_count,
                    prime,
                )
                for prime in primes
            }
            survivors = []
            for order in range(1, cyclotomic_limit + 1):
                cyclotomic = sp.cyclotomic_poly(order, X)
                usable_primes = tuple(
                    prime for prime in primes if order % prime
                )
                if not usable_primes:
                    continue
                survives = True
                for prime in usable_primes:
                    divisor = sp.Poly(cyclotomic, X, modulus=prime)
                    divisions += 1
                    if not ghosts[prime].rem(divisor).is_zero:
                        survives = False
                        break
                if survives:
                    survivors.append(order)

            survivor_table[(left_count, right_count)] = tuple(survivors)
            expected = (1,) if (left_count, right_count) == (1, 1) else ()
            assert survivor_table[(left_count, right_count)] == expected
            types_tested += 1

    print(
        "PASS primitive three-level cyclotomic ghosts: "
        f"{types_tested} coprime types, root orders through "
        f"{cyclotomic_limit}, {divisions} exact finite-field divisions"
    )
    print(
        "three-level root-of-unity survivors: "
        f"{{(1, 1): {survivor_table[(1, 1)]}}}"
    )
    print(
        "STATUS: bounded cyclotomic-order evidence; this does not "
        "classify arbitrary algebraic ghost roots"
    )


def bounded_irreducible_polynomials(
    degree_limit: int,
    height_limit: int,
) -> tuple[sp.Poly, ...]:
    """Enumerate primitive irreducibles in a finite projective box."""

    polynomials = []
    for degree in range(2, degree_limit + 1):
        for leading in range(1, height_limit + 1):
            for tail in product(
                range(-height_limit, height_limit + 1),
                repeat=degree,
            ):
                if tail[-1] == 0:
                    continue
                coefficients = (leading, *tail)
                content = 0
                for coefficient in coefficients:
                    content = gcd(content, abs(coefficient))
                if content != 1:
                    continue
                polynomial = sp.Poly(
                    sum(
                        coefficient * X ** (degree - index)
                        for index, coefficient in enumerate(coefficients)
                    ),
                    X,
                    domain=sp.ZZ,
                )
                if polynomial.is_irreducible:
                    polynomials.append(polynomial)
    return tuple(polynomials)


def verify_three_level_bounded_algebraic_ghosts(
    prime_limit: int,
    step_limit: int,
    degree_limit: int,
    height_limit: int,
) -> None:
    """Exclude bounded minimal-polynomial factors across all good primes."""

    candidates = bounded_irreducible_polynomials(
        degree_limit,
        height_limit,
    )
    divisions = 0
    types_tested = 0

    for left_count in range(1, step_limit + 1):
        for right_count in range(left_count, step_limit + 1):
            if gcd(left_count, right_count) != 1:
                continue
            step = left_count + right_count
            primes = tuple(
                prime
                for prime in primes_through(prime_limit)
                if prime > step
            )
            assert primes
            ghosts = {
                prime: three_level_singleton_ghost_polynomial(
                    left_count,
                    right_count,
                    prime,
                )
                for prime in primes
            }
            survivors = []
            for candidate in candidates:
                leading = int(candidate.LC())
                discriminant = int(sp.discriminant(candidate.as_expr(), X))
                usable_primes = tuple(
                    prime
                    for prime in primes
                    if leading % prime and discriminant % prime
                )
                assert usable_primes
                survives = True
                for prime in usable_primes:
                    divisor = sp.Poly(
                        candidate.as_expr(),
                        X,
                        modulus=prime,
                    )
                    divisions += 1
                    if not ghosts[prime].rem(divisor).is_zero:
                        survives = False
                        break
                if survives:
                    survivors.append(candidate.as_expr())
            assert not survivors
            types_tested += 1

    print(
        "PASS primitive three-level bounded algebraic ghosts: "
        f"{types_tested} coprime types, {len(candidates)} primitive "
        f"irreducibles of degrees 2..{degree_limit} and coefficient "
        f"height {height_limit}, {divisions} finite-field divisions"
    )
    print(
        "STATUS: bounded minimal-polynomial evidence; no degree/height "
        "box can replace an unrestricted algebraic-root theorem"
    )


def verify_three_level_singleton_ghosts(
    prime_limit: int,
    step_limit: int,
    rational_limit: int,
) -> None:
    """Search all primitive three-level affine ghosts in a bounded box.

    The centered ``(1,-2,1)`` move has the universal nonzero root ``X=1``.
    The search asks whether another coprime pair of endpoint multiplicities
    has a rational root surviving every configured good prime.
    """

    rational_values = {
        Fraction(numerator, denominator)
        for denominator in range(1, rational_limit + 1)
        for numerator in range(-2 * rational_limit, 2 * rational_limit + 1)
    }
    survivor_table: dict[tuple[int, int], tuple[Fraction, ...]] = {}
    types_tested = 0
    evaluations = 0

    for left_count in range(1, step_limit + 1):
        for right_count in range(left_count, step_limit + 1):
            if gcd(left_count, right_count) != 1:
                continue
            step = left_count + right_count
            primes = tuple(
                prime
                for prime in primes_through(prime_limit)
                if prime > step
            )
            assert primes
            survivors = []
            for value in sorted(rational_values):
                usable_primes = tuple(
                    prime
                    for prime in primes
                    if value.denominator % prime
                )
                if not usable_primes:
                    continue
                survives = True
                for prime in usable_primes:
                    evaluations += 1
                    if three_level_singleton_ghost_value(
                        left_count,
                        right_count,
                        value,
                        prime,
                    ):
                        survives = False
                        break
                if survives:
                    survivors.append(value)

            survivor_table[(left_count, right_count)] = tuple(survivors)
            expected = (
                (Fraction(0), Fraction(1))
                if (left_count, right_count) == (1, 1)
                else (Fraction(0),)
            )
            assert survivor_table[(left_count, right_count)] == expected
            types_tested += 1

    print(
        "PASS primitive three-level singleton ghosts: "
        f"{types_tested} coprime types through endpoint count "
        f"{step_limit}, primes through {prime_limit}, "
        f"{evaluations} rational/prime evaluations"
    )
    print(
        "three-level nonzero all-prime survivors: "
        f"{{(1, 1): {survivor_table[(1, 1)][1:]}}}"
    )
    print(
        "STATUS: bounded rational-window evidence; only the adjacent "
        "centered triple has a non-support survivor, and this is not an "
        "all-span proof"
    )


def verify_blocks(prime_limit: int, rational_limit: int) -> None:
    primes = primes_through(prime_limit)
    assert primes
    centered_residual_degrees: dict[int, int] = {}
    beta_residual_degrees: dict[int, int] = {}
    centered_extra_roots: dict[int, tuple[int, ...]] = {}
    beta_extra_roots: dict[int, tuple[int, ...]] = {}
    exact_beta_gcd = None
    augmented_beta_zeros: dict[int, tuple[int, ...]] = {}

    for prime in primes:
        centered = centered_block(prime)
        beta = beta_block(prime)

        # [z^p](1+z+z^2)^p is 1 modulo p^2.  This proves that
        # the centered diagonal vanishes at x=1.
        central_trinomial = sum(
            comb(prime, 2 * selection) * comb(2 * selection, selection)
            for selection in range((prime - 1) // 2 + 1)
        )
        assert central_trinomial % prime**2 == 1
        assert int(centered.eval(1)) % prime == 0

        centered_terminal = sp.Poly(X - 1, X, modulus=prime)
        centered_residual, centered_remainder = sp.div(
            centered,
            centered_terminal,
            domain=sp.GF(prime),
        )
        assert centered_remainder.is_zero
        centered_residual_degrees[prime] = centered_residual.degree()

        # Pairing n with p-n proves the harmonic identity K_p(-1)=0.
        assert int(beta.eval(0)) % prime == 0
        assert int(beta.eval(-1)) % prime == 0
        beta_terminal = sp.Poly(X * (X + 1), X, modulus=prime)
        beta_residual, beta_remainder = sp.div(
            beta,
            beta_terminal,
            domain=sp.GF(prime),
        )
        assert beta_remainder.is_zero
        beta_residual_degrees[prime] = beta_residual.degree()
        exact_beta = exact_beta_diagonal(prime)
        exact_beta_gcd = (
            exact_beta
            if exact_beta_gcd is None
            else sp.gcd(exact_beta_gcd, exact_beta)
        )

        centered_extra_roots[prime] = tuple(
            residue
            for residue in range(prime)
            if residue != 1
            and int(centered_residual.eval(residue)) % prime == 0
        )
        beta_extra_roots[prime] = tuple(
            residue
            for residue in range(prime)
            if residue not in (0, prime - 1)
            and int(beta_residual.eval(residue)) % prime == 0
        )
        augmented_beta_zeros[prime] = tuple(
            residue
            for residue in range(prime)
            if (1 + residue) % prime == 0
            and int(beta.eval(residue)) % prime == 0
        )
        assert augmented_beta_zeros[prime] == (prime - 1,)

    rational_values = {
        Fraction(numerator, denominator)
        for denominator in range(1, rational_limit + 1)
        for numerator in range(-rational_limit, rational_limit + 1)
    }
    centered_survivors = []
    beta_survivors = []
    for value in sorted(rational_values):
        usable_primes = [
            prime for prime in primes if value.denominator % prime
        ]
        if not usable_primes:
            continue
        if all(
            polynomial_value_mod(centered_block(prime), value, prime) == 0
            for prime in usable_primes
        ):
            centered_survivors.append(value)
        if all(
            polynomial_value_mod(beta_block(prime), value, prime) == 0
            for prime in usable_primes
        ):
            beta_survivors.append(value)

    assert centered_survivors == [Fraction(1)]
    assert beta_survivors == [Fraction(-1), Fraction(0)]
    assert exact_beta_gcd is not None
    exact_beta_gcd = sp.Poly(
        exact_beta_gcd.monic(),
        X,
        domain=sp.QQ,
    )
    persistent_beta = sp.Poly(
        X * (X + 1) * (X**2 + X + 1),
        X,
        domain=sp.QQ,
    )
    assert exact_beta_gcd == persistent_beta.monic()

    # The centered Bessel endpoint rows are U and U^2+2V in
    # U=C^p, V=D^p.  Their triangular Jacobian at the first row is zero
    # has determinant two, so every odd-prime block forces U=V=0.
    U, V = sp.symbols("U V")
    centered_endpoint_rows = sp.Matrix([U, U**2 + 2 * V])
    centered_endpoint_jacobian = centered_endpoint_rows.jacobian([U, V])
    assert sp.expand(centered_endpoint_jacobian.det()) == 2

    print(
        "PASS ghost-shell universal factors through prime "
        f"{prime_limit}: centered (X-1), beta X(X+1)"
    )
    print(f"centered residual degrees: {centered_residual_degrees}")
    print(f"beta residual degrees: {beta_residual_degrees}")
    print(f"centered prime-dependent extra roots: {centered_extra_roots}")
    print(f"beta prime-dependent extra roots: {beta_extra_roots}")
    print(
        "exact beta diagonal gcd across tested primes: "
        f"{sp.factor(exact_beta_gcd.as_expr())}"
    )
    print(f"augmented beta block zeros: {augmented_beta_zeros}")
    print("centered Bessel endpoint block determinant: 2")
    print(
        "bounded rational survivors: "
        f"centered={centered_survivors}, beta={beta_survivors}"
    )
    print(
        "STATUS: single-row avoidance is false, but the augmented beta "
        "and centered atom blocks are exactly terminal; the displayed "
        "rational window is only a regression"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=43)
    parser.add_argument("--rational-limit", type=int, default=20)
    parser.add_argument("--three-level-limit", type=int, default=6)
    parser.add_argument("--cyclotomic-limit", type=int, default=80)
    parser.add_argument("--algebraic-degree", type=int, default=3)
    parser.add_argument("--algebraic-height", type=int, default=4)
    arguments = parser.parse_args()
    verify_blocks(arguments.prime_limit, arguments.rational_limit)
    verify_three_level_singleton_ghosts(
        arguments.prime_limit,
        arguments.three_level_limit,
        arguments.rational_limit,
    )
    verify_three_level_cyclotomic_ghosts(
        arguments.prime_limit,
        arguments.three_level_limit,
        arguments.cyclotomic_limit,
    )
    verify_three_level_bounded_algebraic_ghosts(
        arguments.prime_limit,
        arguments.three_level_limit,
        arguments.algebraic_degree,
        arguments.algebraic_height,
    )


if __name__ == "__main__":
    main()
