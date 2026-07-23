#!/usr/bin/env python3
"""Adversarial finite-prime tests in the degree-five weighted family.

For any finite list of certified good primes, choose two parameters in the
degree-five stable-moduli family which are congruent at all of those primes.
The resulting integral/rational models have coefficientwise identical
reductions and hence identical finite-field fiber histograms, although the
characteristic-zero Hessian invariant proves that they are not stably
polynomially left-right equivalent.

This does *not* manufacture a permutation reduction.  The weighted maps have
too many rational source points above the target plane C=0 at every good
prime, so injectivity (and therefore surjectivity) is impossible.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import math
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.chebotarev import rational_good_reduction_certificate  # noqa: E402
from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402


def degree_five_primitive(parameter: int) -> sp.Expr:
    """Return H_lambda from the explicit degree-five stable-moduli family."""
    lam = sp.Integer(parameter)
    return sp.factor(
        w**2
        * (w - 1)
        * (3 * w**2 - (5 * lam + 1) * w + 3 * lam)
        / 60
    )


def degree_five_model(parameter: int) -> WeightedSeedModel:
    primitive = degree_five_primitive(parameter)
    seed = sp.diff(primitive, w)
    c = -seed.subs(w, 1)
    return WeightedSeedModel(seed, c=c)


def exceptional_polynomial(parameter: int) -> int:
    """The explicit excluded-parameter polynomial from the moduli theorem."""
    lam = int(parameter)
    return (
        lam
        * (lam - 1)
        * (25 * lam - 1)
        * (10 * lam**2 - 8 * lam + 1)
        * (100 * lam**2 - 29 * lam + 10)
    )


def hessian_invariant(parameter: int) -> sp.Rational:
    """The stable-equivalence obstruction J(lambda)."""
    lam = sp.Integer(parameter)
    return sp.factor(
        sp.Rational(64, 75)
        * (50 * lam**2 - 40 * lam + 17) ** 3
        / (10 * lam**2 - 8 * lam + 1) ** 2
    )


def modular_terms(expression: sp.Expr, prime: int):
    """Sparse terms with rational coefficients reduced modulo ``prime``."""
    terms = []
    for monomial, coefficient in sp.Poly(expression, x, y, z).terms():
        numerator, denominator = map(int, sp.fraction(coefficient))
        if denominator % prime == 0:
            raise ValueError(f"coefficient denominator is not a unit modulo {prime}")
        terms.append(
            (monomial, numerator * pow(denominator, -1, prime) % prime)
        )
    return tuple(terms)


def reduced_mapping(model: WeightedSeedModel, prime: int):
    return tuple(modular_terms(component, prime) for component in model.mapping())


def evaluate_terms(terms, point, prime: int) -> int:
    return sum(
        coefficient
        * math.prod(
            pow(coordinate, exponent, prime)
            for coordinate, exponent in zip(point, monomial)
        )
        for monomial, coefficient in terms
    ) % prime


def full_fiber_histogram(reduced_map, prime: int) -> dict[int, int]:
    """Count targets by their number of F_p-rational source points."""
    fibers = Counter(
        tuple(evaluate_terms(terms, source, prime) for terms in reduced_map)
        for source in product(range(prime), repeat=3)
    )
    histogram = Counter(fibers.values())
    histogram[0] = prime**3 - len(fibers)
    return dict(sorted(histogram.items()))


def parse_primes(text: str) -> tuple[int, ...]:
    primes = tuple(int(value) for value in text.split(",") if value)
    if not primes or len(set(primes)) != len(primes):
        raise ValueError("provide a nonempty list of distinct primes")
    if any(not sp.isprime(prime) for prime in primes):
        raise ValueError("every modulus must be prime")
    return primes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="11,13,17,19,23")
    parser.add_argument("--base-parameter", type=int, default=2)
    args = parser.parse_args()

    primes = parse_primes(args.primes)
    parameter_1 = args.base_parameter
    if exceptional_polynomial(parameter_1) == 0:
        raise ValueError("the base parameter lies on the explicit exceptional locus")

    modulus = math.prod(primes)
    shift_multiplier = 1
    invariant_1 = hessian_invariant(parameter_1)
    while True:
        parameter_2 = parameter_1 + shift_multiplier * modulus
        if (
            exceptional_polynomial(parameter_2) != 0
            and hessian_invariant(parameter_2) != invariant_1
        ):
            break
        shift_multiplier += 1

    model_1 = degree_five_model(parameter_1)
    model_2 = degree_five_model(parameter_2)
    invariant_2 = hessian_invariant(parameter_2)
    assert invariant_1 != invariant_2

    certificate_1 = rational_good_reduction_certificate(
        model_1.primitive, w, c=model_1.c, b0=model_1.b, a0=model_1.a
    )
    certificate_2 = rational_good_reduction_certificate(
        model_2.primitive, w, c=model_2.c, b0=model_2.b, a0=model_2.a
    )

    print(f"parameters: lambda={parameter_1}, mu={parameter_2}")
    print(f"parameter difference: {parameter_2 - parameter_1}")
    print(f"CRT shift multiplier: {shift_multiplier}")
    print(f"J(lambda)={invariant_1}")
    print(f"J(mu)={invariant_2}")
    print("PASS: unequal Hessian invariants obstruct stable left-right equivalence")

    for prime in primes:
        assert parameter_1 % prime == parameter_2 % prime
        assert certificate_1["bad_integer"] % prime != 0
        assert certificate_2["bad_integer"] % prime != 0

        reduction_1 = reduced_mapping(model_1, prime)
        reduction_2 = reduced_mapping(model_2, prime)
        assert reduction_1 == reduction_2
        histogram = full_fiber_histogram(reduction_1, prime)
        assert histogram.get(0, 0) > 0
        assert any(size > 1 and count for size, count in histogram.items())
        print(
            f"PASS F_{prime}: certified good, coefficientwise identical, "
            f"fiber histogram {histogram}"
        )

    print("PASS: no selected reduction is injective or surjective")
    print("PASS: any prescribed finite good-prime test set can be fooled this way")


if __name__ == "__main__":
    main()
