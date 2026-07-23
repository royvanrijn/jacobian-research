#!/usr/bin/env python3
"""Reduce the bounded hbar^5 period conditions on the 16-term support.

At an exact generic-chart fiber this script constructs:

* all hbar^5 correction-image conditions on a functional supported on the
  16 monomials used by the all-pole certificates;
* all first variations along the complete bounded hbar^3 affine family; and
* all quadratic lower-lift variations.

It then reports an exact independent subset.  At ``(kappa,tau)=(1,1)``,
only 15 conditions are independent, so the desired functional is the
one-dimensional kernel of a 15-by-16 system.  This is a reduction for the
future function-field computation, not itself a parameter-uniform proof.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement

from sympy.polys.domains import QQ
from sympy.polys.matrices.sdm import sdm_irref

from explore_degree_five_quantum_residue import (
    add,
    degree_five_family,
    laurent_monomials,
    pi_power,
    poisson,
    third_order_family,
)
from verify_degree_five_laurent_quantum_obstruction import SUPPORT


def parse_rational(value: str):
    parsed = Fraction(value)
    return QQ(parsed.numerator, parsed.denominator)


def support_vector(poly):
    return {
        index: poly.get(monomial, QQ.zero)
        for index, monomial in enumerate(SUPPORT)
        if poly.get(monomial, QQ.zero)
    }


def analyze(kappa, tau) -> tuple[list[tuple], list[int]]:
    if kappa in (QQ(-2), QQ(-1)):
        raise ValueError("this analyzer uses the generic chart kappa != -2,-1")
    a = -(QQ.one + kappa) / (QQ(2) + kappa)
    S, T = degree_five_family(QQ, a, tau)
    family = third_order_family(S, T, QQ)
    base_s, base_t = family.base

    constraints: list[dict[int, object]] = []
    labels: list[tuple] = []

    for label, monomials, left in (
        ("d5-S", laurent_monomials(21, 1, 0, 3), True),
        ("d5-T", laurent_monomials(17, 0, 0, 3), False),
    ):
        for monomial in monomials:
            image = (
                poisson({monomial: QQ.one}, T)
                if left
                else poisson(S, {monomial: QQ.one})
            )
            vector = support_vector(image)
            if vector:
                constraints.append(vector)
                labels.append((label, monomial))

    for index, (direction_s, direction_t) in enumerate(family.kernel):
        variation = add(
            poisson(direction_s, base_t),
            poisson(base_s, direction_t),
        )
        variation = add(
            variation,
            pi_power(direction_s, T, 3),
            QQ.one / QQ(24),
        )
        variation = add(
            variation,
            pi_power(S, direction_t, 3),
            QQ.one / QQ(24),
        )
        vector = support_vector(variation)
        if vector:
            constraints.append(vector)
            labels.append(("linear", index))

    for left, right in combinations_with_replacement(
        range(len(family.kernel)),
        2,
    ):
        left_s, left_t = family.kernel[left]
        right_s, right_t = family.kernel[right]
        variation = (
            poisson(left_s, left_t)
            if left == right
            else add(
                poisson(left_s, right_t),
                poisson(right_s, left_t),
            )
        )
        vector = support_vector(variation)
        if vector:
            constraints.append(vector)
            labels.append(("quadratic", left, right))

    # The rows below are functional coordinates and the columns are
    # conditions.  Pivot columns therefore select independent conditions.
    transpose = {
        coordinate: {
            condition: vector.get(coordinate, QQ.zero)
            for condition, vector in enumerate(constraints)
            if vector.get(coordinate, QQ.zero)
        }
        for coordinate in range(len(SUPPORT))
    }
    _, pivots, _ = sdm_irref(transpose)
    return labels, pivots


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kappa", default="1")
    parser.add_argument("--tau", default="1")
    args = parser.parse_args()
    kappa = parse_rational(args.kappa)
    tau = parse_rational(args.tau)

    labels, pivots = analyze(kappa, tau)
    print(
        "bounded hbar^5 period constraints at "
        f"(kappa,tau)=({kappa},{tau})"
    )
    print(f"  support dimension: {len(SUPPORT)}")
    print(f"  nonzero conditions: {len(labels)}")
    print(f"  condition types: {dict(Counter(label[0] for label in labels))}")
    print(f"  exact rank: {len(pivots)}")
    print(f"  functional kernel dimension: {len(SUPPORT) - len(pivots)}")
    print("  independent conditions:")
    for pivot in pivots:
        print(f"    {labels[pivot]}")


if __name__ == "__main__":
    main()
