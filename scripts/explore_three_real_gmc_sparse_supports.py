#!/usr/bin/env python3
"""Explore sparse three-real Gaussian supports in circular coordinates.

This is a discovery tool, not a characteristic-zero or all-order proof.
It performs two rigorous finite-cutoff operations:

1. reject a full support if one of the tested moments consists of a single
   coefficient monomial (which cannot vanish on the coefficient torus);
2. test the remaining truncated moment ideal directly on the coefficient
   torus over a finite field.

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
ChargeParitySignature = tuple[tuple[int, int], ...]


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
    signature = contraction_signature(support)
    for multiplicities in contraction_vectors(signature, moment):
        w_degree = sum(n * monomial[0] for n, monomial in zip(multiplicities, support))
        z_degree = sum(n * monomial[1] for n, monomial in zip(multiplicities, support))
        t_degree = sum(n * monomial[2] for n, monomial in zip(multiplicities, support))
        assert w_degree == z_degree
        assert t_degree % 2 == 0
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


def coefficient_charge(monomial: Monomial) -> int:
    """Weight for W -> tW, Z -> t^-1 Z."""
    return monomial[0] - monomial[1]


def contraction_signature(
    support: tuple[Monomial, ...],
) -> ChargeParitySignature:
    """Charge/parity data, canonicalized under simultaneous charge reversal."""
    signature = tuple(
        (coefficient_charge(monomial), monomial[2] % 2)
        for monomial in support
    )
    opposite = tuple((-charge, parity) for charge, parity in signature)
    return min(signature, opposite)


@lru_cache(maxsize=None)
def contraction_vectors(
    signature: ChargeParitySignature, moment: int
) -> tuple[Composition, ...]:
    """Multiplicity vectors of total size `moment`, generated from the Hilbert basis."""
    zero = (0,) * len(signature)
    basis = contraction_hilbert_basis_from_signature(signature)
    by_degree: list[set[Composition]] = [set() for _ in range(moment + 1)]
    by_degree[0].add(zero)
    for degree in range(moment + 1):
        for vector in tuple(by_degree[degree]):
            for generator in basis:
                next_degree = degree + sum(generator)
                if next_degree <= moment:
                    by_degree[next_degree].add(_vector_add(vector, generator))
    return tuple(sorted(by_degree[moment]))


def _vector_add(left: Composition, right: Composition) -> Composition:
    return tuple(a + b for a, b in zip(left, right))


def _vector_subtract(left: Composition, right: Composition) -> Composition:
    return tuple(a - b for a, b in zip(left, right))


def _vector_le(left: Composition, right: Composition) -> bool:
    return all(a <= b for a, b in zip(left, right))


def _generated_by(target: Composition, generators: tuple[Composition, ...]) -> bool:
    """Whether `target` is in the additive monoid generated by `generators`."""
    zero = (0,) * len(target)
    reachable = {zero}
    frontier = [zero]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = _vector_add(current, generator)
            if not _vector_le(candidate, target) or candidate in reachable:
                continue
            if candidate == target:
                return True
            reachable.add(candidate)
            frontier.append(candidate)
    return False


@lru_cache(maxsize=None)
def charge_hilbert_basis(charges: tuple[int, ...]) -> tuple[Composition, ...]:
    """Hilbert basis of {nu >= 0 : charges dot nu = 0}.

    Pottier's one-row bound is 1 + sum(abs(q)); enumerating through that
    total degree therefore finds every coordinatewise-minimal solution.
    """
    length = len(charges)
    zero = (0,) * length
    bound = 1 + sum(abs(charge) for charge in charges)
    solutions: set[Composition] = {zero}
    basis: list[Composition] = []
    for total in range(1, bound + 1):
        for vector in compositions(total, length):
            if sum(entry * charge for entry, charge in zip(vector, charges)) != 0:
                continue
            decomposable = any(
                part != zero
                and _vector_le(part, vector)
                and _vector_subtract(vector, part) in solutions
                for part in solutions
            )
            if not decomposable:
                basis.append(vector)
            solutions.add(vector)
    return tuple(basis)


@lru_cache(maxsize=None)
def contraction_hilbert_basis_from_signature(
    signature: ChargeParitySignature,
) -> tuple[Composition, ...]:
    """Hilbert basis after imposing even total T-parity.

    Decompose first into charge atoms.  An even-parity sum consists of
    even atoms and pairs of odd atoms.  Minimalizing those candidates gives
    the Hilbert basis of the parity kernel.
    """
    charges = tuple(charge for charge, _ in signature)
    parities = tuple(parity for _, parity in signature)
    charge_basis = charge_hilbert_basis(charges)

    def parity(vector: Composition) -> int:
        return sum(
            entry * entry_parity
            for entry, entry_parity in zip(vector, parities)
        ) % 2

    even_atoms = [vector for vector in charge_basis if parity(vector) == 0]
    odd_atoms = [vector for vector in charge_basis if parity(vector) == 1]
    candidates = set(even_atoms)
    candidates.update(
        _vector_add(left, right)
        for left_index, left in enumerate(odd_atoms)
        for right in odd_atoms[left_index:]
    )

    basis: list[Composition] = []
    for vector in sorted(candidates, key=lambda item: (sum(item), item)):
        if not _generated_by(vector, tuple(basis)):
            basis.append(vector)
    return tuple(basis)


def contraction_hilbert_basis(
    support: tuple[Monomial, ...],
) -> tuple[Composition, ...]:
    return contraction_hilbert_basis_from_signature(contraction_signature(support))


def normalization_anchors(support: tuple[Monomial, ...]) -> tuple[int, ...]:
    """Coefficients fixed to one on an algebraic-closure orbit slice.

    Overall scaling fixes the first coefficient.  If another monomial has
    distinct Gaussian-torus charge, the torus action fixes that coefficient
    as well.  When all charges agree (necessarily charge zero for an active
    non-one-sided support), only projective normalization is available.
    """
    anchors = [0]
    first_charge = coefficient_charge(support[0])
    for index, monomial in enumerate(support[1:], start=1):
        if coefficient_charge(monomial) != first_charge:
            anchors.append(index)
            break
    return tuple(anchors)


def free_coefficient_indices(support: tuple[Monomial, ...]) -> tuple[int, ...]:
    anchors = set(normalization_anchors(support))
    return tuple(index for index in range(len(support)) if index not in anchors)


def invariant_ratio_exponents(
    support: tuple[Monomial, ...], free_index: int
) -> tuple[int, ...]:
    """Laurent exponents of a quotient-torus invariant.

    For anchors a,b with charge difference d=q_b-q_a, the invariant is
    c_s^d c_a^(q_s-q_b) c_b^(q_a-q_s).  On the slice c_a=c_b=1 it becomes
    c_s^d, so the slice is a finite cover of the invariant-ratio quotient.
    """
    anchors = normalization_anchors(support)
    if len(anchors) != 2 or free_index in anchors:
        raise ValueError("invariant ratios require two anchors and a free index")
    anchor_a, anchor_b = anchors
    charge_a = coefficient_charge(support[anchor_a])
    charge_b = coefficient_charge(support[anchor_b])
    charge_s = coefficient_charge(support[free_index])
    difference = charge_b - charge_a
    exponents = [0] * len(support)
    exponents[free_index] = difference
    exponents[anchor_a] = charge_s - charge_b
    exponents[anchor_b] = charge_a - charge_s
    assert sum(exponents) == 0
    assert sum(
        exponent * coefficient_charge(monomial)
        for exponent, monomial in zip(exponents, support)
    ) == 0
    return tuple(exponents)


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
        signature = contraction_signature(support)
        primitive_contractions = contraction_hilbert_basis_from_signature(signature)
        if len(primitive_contractions) == 1:
            # A rank-one contraction semigroup has a unique monomial at the
            # degree of its primitive generator, so it cannot vanish on the
            # coefficient torus.
            singleton_rejections += 1
            continue
        moment_contractions = tuple(
            contraction_vectors(signature, moment)
            for moment in range(1, moment_bound + 1)
        )
        if any(len(vectors) == 1 for vectors in moment_contractions):
            singleton_rejections += 1
            continue
        # If every tested moment vanishes termwise, the support has no active
        # cancellation problem at this cutoff and is reported separately.
        if not any(moment_contractions):
            continue
        moment_terms = tuple(
            wick_moment_terms(support, moment)
            for moment in range(1, moment_bound + 1)
        )
        candidates.append(Candidate(support, moment_terms))
    return raw_count, singleton_rejections, tuple(candidates)


def coefficient_monomial(
    multiplicities: Composition, free_indices: tuple[int, ...]
) -> str:
    # Slice coordinates are x1, x2, ...; anchor coefficients are set to one.
    factors: list[str] = []
    for variable_index, coefficient_index in enumerate(free_indices, start=1):
        exponent = multiplicities[coefficient_index]
        if exponent == 1:
            factors.append(f"x{variable_index}")
        elif exponent > 1:
            factors.append(f"x{variable_index}^{exponent}")
    return "*".join(factors) if factors else "1"


def singular_polynomial(
    terms: tuple[tuple[int, Composition], ...],
    free_indices: tuple[int, ...],
) -> str:
    return "+".join(
        f"{coefficient}*{coefficient_monomial(multiplicities, free_indices)}"
        for coefficient, multiplicities in terms
    )


def laurent_screen(
    candidates: tuple[Candidate, ...], prime: int, timeout: int
) -> tuple[Candidate, ...]:
    """Test the moment ideal in the Laurent coefficient ring.

    On the normalization slice the free coefficients are units.  Singular's
    ideal-quotient saturation computes the contraction

        (I K[x1^+-1,...,xr^+-1]) intersect K[x1,...,xr]

    directly, without adjoining an inverse or Rabinowitsch variable.  Thus
    the standard-basis calculation has only the genuine quotient coordinates.
    """
    commands: list[str] = ['LIB "elim.lib";']
    for free_count in range(4):
        variables = [f"x{index}" for index in range(1, free_count + 1)]
        if not variables:
            variables = ["dummy"]
        commands.append(
            f"ring r{free_count}={prime},({','.join(variables)}),dp;"
        )
    for index, candidate in enumerate(candidates):
        free_indices = free_coefficient_indices(candidate.support)
        commands.append(f"setring r{len(free_indices)};")
        equations = [
            singular_polynomial(terms, free_indices)
            for terms in candidate.moment_terms
            if terms
        ]
        coefficient_product = "*".join(
            f"x{index}" for index in range(1, len(free_indices) + 1)
        ) or "1"
        commands.append(f"ideal I={','.join(equations)};")
        commands.append(
            f"ideal G=sat(I,ideal({coefficient_product}))[1];"
        )
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
    survivors = laurent_screen(candidates, args.prime, args.timeout)
    if args.prime:
        label = f"modular survivors at p={args.prime}"
    else:
        label = "exact rational survivors"
    print(f"{label}, moments<= {args.moments}: {len(survivors)}")
    for candidate in survivors[: args.show]:
        print(format_support(candidate.support))


if __name__ == "__main__":
    main()
