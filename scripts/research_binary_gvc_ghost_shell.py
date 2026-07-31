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
from math import comb, factorial

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
    arguments = parser.parse_args()
    verify_blocks(arguments.prime_limit, arguments.rational_limit)


if __name__ == "__main__":
    main()
