#!/usr/bin/env python3
"""Explore sparse three-real Gaussian supports in circular coordinates.

This is a discovery tool, not a characteristic-zero or all-order proof.
It performs two rigorous finite-cutoff operations:

1. reject a full support if one of the tested moments consists of a single
   coefficient monomial (which cannot vanish on the coefficient torus);
2. test the remaining truncated saturated moment ideal over a finite field.

Every modular survivor needs characteristic-zero analysis.  A modular
exclusion also needs an exact replay before it is promoted to a certificate.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from math import factorial
import subprocess


Monomial = tuple[int, int, int]  # exponents of W, Z, T
Composition = tuple[int, ...]


@dataclass(frozen=True)
class Candidate:
    support: tuple[Monomial, ...]
    moment_terms: tuple[tuple[tuple[int, Composition], ...], ...]


def monomials(degree: int) -> tuple[Monomial, ...]:
    return tuple(
        (w_degree, z_degree, total_degree - w_degree - z_degree)
        for total_degree in range(degree + 1)
        for w_degree in range(total_degree + 1)
        for z_degree in range(total_degree - w_degree + 1)
    )


@lru_cache(maxsize=None)
def compositions(total: int, length: int) -> tuple[Composition, ...]:
    if length == 1:
        return ((total,),)
    return tuple(
        (first,) + tail
        for first in range(total + 1)
        for tail in compositions(total - first, length - 1)
    )


def odd_double_factorial(even_exponent: int) -> int:
    value = 1
    for factor in range(1, even_exponent, 2):
        value *= factor
    return value


def wick_moment_terms(
    support: tuple[Monomial, ...], moment: int
) -> tuple[tuple[int, Composition], ...]:
    """Return coefficient-monomial terms of E[P^moment]."""
    out: list[tuple[int, Composition]] = []
    for multiplicities in compositions(moment, len(support)):
        w_degree = sum(n * monomial[0] for n, monomial in zip(multiplicities, support))
        z_degree = sum(n * monomial[1] for n, monomial in zip(multiplicities, support))
        t_degree = sum(n * monomial[2] for n, monomial in zip(multiplicities, support))
        if w_degree != z_degree or t_degree % 2:
            continue
        multinomial = factorial(moment)
        for multiplicity in multiplicities:
            multinomial //= factorial(multiplicity)
        coefficient = (
            multinomial
            * factorial(w_degree)
            * odd_double_factorial(t_degree)
        )
        out.append((coefficient, multiplicities))
    return tuple(out)


def reflected_support(support: tuple[Monomial, ...]) -> tuple[Monomial, ...]:
    return tuple(sorted((z_degree, w_degree, t_degree) for w_degree, z_degree, t_degree in support))


def canonical_support(support: tuple[Monomial, ...]) -> tuple[Monomial, ...]:
    ordered = tuple(sorted(support))
    return min(ordered, reflected_support(ordered))


def enumerate_candidates(
    degree: int, term_count: int, moment_bound: int
) -> tuple[int, int, tuple[Candidate, ...]]:
    raw_count = 0
    singleton_rejections = 0
    candidates: list[Candidate] = []
    seen: set[tuple[Monomial, ...]] = set()
    for raw_support in combinations(monomials(degree), term_count):
        raw_count += 1
        support = canonical_support(raw_support)
        if support in seen:
            continue
        seen.add(support)
        moment_terms = tuple(
            wick_moment_terms(support, moment)
            for moment in range(1, moment_bound + 1)
        )
        if any(len(terms) == 1 for terms in moment_terms):
            singleton_rejections += 1
            continue
        # If every tested moment vanishes termwise, the support has no active
        # cancellation problem at this cutoff and is reported separately.
        if not any(moment_terms):
            continue
        candidates.append(Candidate(support, moment_terms))
    return raw_count, singleton_rejections, tuple(candidates)


def coefficient_monomial(multiplicities: Composition) -> str:
    # Projective scaling permits c0=1.  The remaining coefficients are x1...
    factors: list[str] = []
    for index, exponent in enumerate(multiplicities[1:], start=1):
        if exponent == 1:
            factors.append(f"x{index}")
        elif exponent > 1:
            factors.append(f"x{index}^{exponent}")
    return "*".join(factors) if factors else "1"


def singular_polynomial(terms: tuple[tuple[int, Composition], ...]) -> str:
    return "+".join(
        f"{coefficient}*{coefficient_monomial(multiplicities)}"
        for coefficient, multiplicities in terms
    )


def modular_screen(
    candidates: tuple[Candidate, ...], prime: int, timeout: int
) -> tuple[Candidate, ...]:
    commands = [f"ring r={prime},(x1,x2,x3,u),dp;"]
    for index, candidate in enumerate(candidates):
        equations = [
            singular_polynomial(terms)
            for terms in candidate.moment_terms
            if terms
        ]
        coefficient_product = "*".join(
            f"x{index}" for index in range(1, len(candidate.support))
        )
        saturation = f"u*{coefficient_product}-1"
        commands.append(f"ideal I={','.join(equations + [saturation])};")
        commands.append("ideal G=std(I);")
        commands.append(
            f'if (reduce(1,G)==0) {{ print("@{index} unit"); }} '
            f'else {{ print("@{index} keep"); }};'
        )
        commands.append("kill I,G;")
    commands.append("quit;")
    completed = subprocess.run(
        ["Singular", "-q"],
        input="\n".join(commands),
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    decisions: dict[int, str] = {}
    for line in completed.stdout.splitlines():
        if not line.startswith("@"):
            continue
        raw_index, decision = line[1:].split()
        decisions[int(raw_index)] = decision
    if len(decisions) != len(candidates):
        raise RuntimeError(
            f"Singular returned {len(decisions)} decisions for "
            f"{len(candidates)} candidates.\n{completed.stderr}"
        )
    return tuple(
        candidate
        for index, candidate in enumerate(candidates)
        if decisions[index] == "keep"
    )


def format_support(support: tuple[Monomial, ...]) -> str:
    return "{" + ", ".join(f"W^{w}Z^{z}T^{t}" for w, z, t in support) + "}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, default=4)
    parser.add_argument("--terms", type=int, default=4)
    parser.add_argument("--moments", type=int, default=8)
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--show", type=int, default=25)
    args = parser.parse_args()

    raw_count, singleton_rejections, candidates = enumerate_candidates(
        args.degree, args.terms, args.moments
    )
    print(f"raw supports: {raw_count}")
    print(f"singleton rejections after W/Z quotient: {singleton_rejections}")
    print(f"active candidates sent to Singular: {len(candidates)}")
    survivors = modular_screen(candidates, args.prime, args.timeout)
    print(
        f"modular survivors at p={args.prime}, "
        f"moments<= {args.moments}: {len(survivors)}"
    )
    for candidate in survivors[: args.show]:
        print(format_support(candidate.support))


if __name__ == "__main__":
    main()
