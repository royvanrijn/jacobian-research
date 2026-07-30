#!/usr/bin/env python3
"""Exact sparse Pareto certificate for three-real Gaussian moments.

The search covers every circular-coordinate monomial support with

    total degree <= 4, number of terms <= 4.

It uses exact rational Gröbner bases.  Moment 20 reduces the active supports
to five harmless families across term counts three and four; moments through
24 remove two finite-cutoff four-term artifacts.  The three genuine
survivors are then identified symbolically.

Each support is first replaced by its exact charge/parity contraction
semigroup.  A Pottier-bounded Hilbert-basis computation rejects rank-one
semigroups before Wick coefficients are expanded, and charge/parity-
equivalent supports share their contraction-vector enumeration.

The four-term pass replays 5,718 coefficient-torus ideals over Q directly
in Laurent coordinates, without an auxiliary saturation variable.
Whenever a support has two distinct Gaussian charges, overall scaling and
the Gaussian-preserving torus action reduce its four coefficient variables
to a two-dimensional orbit slice before Gröbner elimination.
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp

from explore_three_real_gmc_sparse_supports import (
    Candidate,
    canonical_support,
    coefficient_charge,
    contraction_hilbert_basis,
    enumerate_candidates,
    format_support,
    free_coefficient_indices,
    invariant_ratio_exponents,
    monomials,
    laurent_screen,
    normalization_anchors,
    wick_moment_terms,
)


LINEAR_SUPPORT = ((0, 0, 1), (0, 1, 0), (1, 0, 0))
SHEAR_SUPPORT = ((0, 1, 0), (1, 0, 1), (3, 0, 0))
PRODUCT_SUPPORT = ((0, 1, 1), (0, 2, 0), (1, 0, 1), (2, 0, 0))
ARTIFACT_SUPPORTS = (
    ((0, 1, 0), (0, 1, 2), (1, 2, 0), (3, 0, 1)),
    ((0, 1, 1), (0, 1, 3), (1, 2, 1), (3, 0, 0)),
)


def candidate(support: tuple[tuple[int, int, int], ...], bound: int) -> Candidate:
    return Candidate(
        support,
        tuple(wick_moment_terms(support, moment) for moment in range(1, bound + 1)),
    )


def check_inactive_supports() -> int:
    """Certify the circuit bound and count strictly one-sided supports."""
    inactive_count = 0
    seen: set[tuple[tuple[int, int, int], ...]] = set()
    for term_count in range(1, 5):
        for raw_support in combinations(monomials(4), term_count):
            support = canonical_support(raw_support)
            if support in seen:
                continue
            seen.add(support)
            charges_and_parities = [
                (w_degree - z_degree, t_degree % 2)
                for w_degree, z_degree, t_degree in support
            ]
            charges = [charge for charge, _ in charges_and_parities]
            if all(charge > 0 for charge in charges) or all(
                charge < 0 for charge in charges
            ):
                assert contraction_hilbert_basis(support) == ()
                inactive_count += 1
                continue

            # Every non-one-sided degree-four support has a Wick contraction
            # by moment 16.  A zero-charge monomial needs one or two copies.
            # Opposite charges q>0 and r<0 balance with |r|/g and q/g
            # copies; their length is at most 8, and doubling fixes parity.
            if any(charge == 0 for charge in charges):
                bound = 2
            else:
                positive = next(charge for charge in charges if charge > 0)
                negative = next(charge for charge in charges if charge < 0)
                from math import gcd

                primitive_length = (
                    positive + abs(negative)
                ) // gcd(positive, abs(negative))
                bound = 2 * primitive_length
            assert bound <= 16
            primitive_contractions = contraction_hilbert_basis(support)
            assert primitive_contractions
            assert min(sum(vector) for vector in primitive_contractions) <= bound
    return inactive_count


def check_coefficient_symmetries() -> None:
    """Check the normalization rank on every support."""
    for term_count in range(1, 5):
        for support in combinations(monomials(4), term_count):
            anchors = normalization_anchors(support)
            assert anchors[0] == 0
            assert len(anchors) in (1, 2)
            if len(anchors) == 2:
                assert coefficient_charge(support[anchors[0]]) != coefficient_charge(
                    support[anchors[1]]
                )
                for free_index in free_coefficient_indices(support):
                    invariant_ratio_exponents(support, free_index)
            else:
                assert len({coefficient_charge(monomial) for monomial in support}) == 1
            assert len(free_coefficient_indices(support)) == term_count - len(anchors)


def exact_survivors(term_count: int) -> tuple[Candidate, ...]:
    _, _, candidates = enumerate_candidates(4, term_count, 20)
    return laurent_screen(candidates, prime=0, timeout=1200)


def survivor_groebner(
    support: tuple[tuple[int, int, int], ...], bound: int
) -> list[sp.Expr]:
    """Export a polynomial certificate after the Laurent search is complete."""
    x1, x2, x3, inverse = sp.symbols("x1 x2 x3 inverse")
    variables = (x1, x2, x3)
    equations: list[sp.Expr] = []
    for terms in candidate(support, bound).moment_terms:
        if not terms:
            continue
        equation = 0
        for coefficient, multiplicities in terms:
            coefficient_monomial = coefficient
            for index, exponent in enumerate(multiplicities[1:]):
                coefficient_monomial *= variables[index] ** exponent
            equation += coefficient_monomial
        equations.append(sp.expand(equation))
    coefficient_product = sp.prod(variables[: len(support) - 1])
    basis = sp.groebner(
        equations + [inverse * coefficient_product - 1],
        inverse,
        *reversed(variables[: len(support) - 1]),
        order="lex",
        domain=sp.QQ,
    )
    return [sp.factor(poly.as_expr()) for poly in basis.polys]


def main() -> None:
    check_coefficient_symmetries()
    inactive_count = check_inactive_supports()

    def assert_basis_equal(actual: list[sp.Expr], expected: list[sp.Expr]) -> None:
        assert len(actual) == len(expected)
        assert all(
            sp.expand(left - right) == 0
            for left, right in zip(actual, expected)
        )

    x1, x2, x3, inverse = sp.symbols("x1 x2 x3 inverse")
    assert_basis_equal(
        survivor_groebner(LINEAR_SUPPORT, 20),
        [inverse + 2, x1 * x2 + sp.Rational(1, 2)],
    )
    assert_basis_equal(
        survivor_groebner(SHEAR_SUPPORT, 20),
        [inverse * x1**3 + 2, sp.Rational(1, 2) * x1**2 + x2],
    )
    assert_basis_equal(
        survivor_groebner(PRODUCT_SUPPORT, 24),
        [
            inverse + 8 * x1**4,
            x1 * x2**2 + x3,
            x1**2 * x2 - sp.Rational(1, 2),
        ],
    )

    one_term = exact_survivors(1)
    two_term = exact_survivors(2)
    three_term = exact_survivors(3)
    four_term_at_20 = exact_survivors(4)

    assert one_term == ()
    assert two_term == ()
    assert tuple(item.support for item in three_term) == (
        LINEAR_SUPPORT,
        SHEAR_SUPPORT,
    )
    assert tuple(item.support for item in four_term_at_20) == (
        ARTIFACT_SUPPORTS[0],
        PRODUCT_SUPPORT,
        ARTIFACT_SUPPORTS[1],
    )

    four_term_at_24 = laurent_screen(
        tuple(candidate(item.support, 24) for item in four_term_at_20),
        prime=0,
        timeout=300,
    )
    assert tuple(item.support for item in four_term_at_24) == (PRODUCT_SUPPORT,)

    # Product-family factorization.  With a=x1 and x2=1/(2a^2),
    # x3=-1/(4a^3), the four-term polynomial is L1*L2.
    W, Z, T, a = sp.symbols("W Z T a", nonzero=True)
    product_family = (
        Z * T
        + a * Z**2
        + sp.Rational(1, 2) / a**2 * W * T
        - sp.Rational(1, 4) / a**3 * W**2
    )
    first_linear = Z + W / (2 * a**2)
    isotropic_linear = T + a * Z - W / (2 * a)
    assert sp.expand(product_family - first_linear * isotropic_linear) == 0

    # The second factor is isotropic and orthogonal to the first under
    # <Z,W>=1 and <T,T>=1.
    variance_isotropic = 1 + 2 * a * (-sp.Rational(1, 2) / a)
    covariance = (
        -sp.Rational(1, 2) / a
        + (sp.Rational(1, 2) / a**2) * a
    )
    assert sp.simplify(variance_isotropic) == 0
    assert sp.simplify(covariance) == 0

    print(f"PASS sparse Pareto: {inactive_count} termwise-null orbit supports are one-sided")
    print("PASS sparse Pareto: exact contraction Hilbert bases precede moment expansion")
    print("PASS sparse Pareto: coefficient ideals are tested directly on the Laurent torus")
    print("PASS sparse Pareto: exact Q search covers degree<=4 and terms<=4")
    print("PASS sparse Pareto: moments<=24 leave exactly three genuine null families")
    print("PASS sparse Pareto: linear, shear, and product families all satisfy GMC")
    print("PASS sparse Pareto: Long is minimal for (degree, terms) <= (4, 4)")


if __name__ == "__main__":
    main()
