#!/usr/bin/env python3
"""Verify central-vertex Hurwitz multiplicity and root-mark selection."""

from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "central_hurwitz_selection_degree6.json"
)


Permutation = tuple[int, ...]
DEGREE = 6
IDENTITY = tuple(range(DEGREE))


def compose(left: Permutation, right: Permutation) -> Permutation:
    """Return left after right."""

    return tuple(left[right[index]] for index in range(DEGREE))


def inverse(permutation: Permutation) -> Permutation:
    result = [0] * DEGREE
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def conjugate(element: Permutation, permutation: Permutation) -> Permutation:
    return compose(compose(element, permutation), inverse(element))


def matching_permutations(
    remaining: tuple[int, ...],
) -> list[Permutation]:
    """All fixed-point-free involutions on the specified six letters."""

    pairings: list[tuple[tuple[int, int], ...]] = []

    def pair(
        items: tuple[int, ...],
        current: tuple[tuple[int, int], ...],
    ) -> None:
        if not items:
            pairings.append(current)
            return
        first = items[0]
        for second in items[1:]:
            pair(
                tuple(
                    item
                    for item in items
                    if item not in (first, second)
                ),
                current + ((first, second),),
            )

    pair(remaining, ())
    result = []
    for pairing in pairings:
        permutation = list(range(DEGREE))
        for first, second in pairing:
            permutation[first], permutation[second] = second, first
        result.append(tuple(permutation))
    return result


infinity_monodromy = (1, 2, 3, 4, 5, 0)
transpositions = []
for first, second in combinations(range(DEGREE), 2):
    permutation = list(range(DEGREE))
    permutation[first], permutation[second] = second, first
    transpositions.append(tuple(permutation))
zero_monodromies = matching_permutations(tuple(range(DEGREE)))
assert len(transpositions) == 15
assert len(zero_monodromies) == 15

# Ordered branch cycles have profiles
#   (6), (2,2,2), (2,1^4), (2,1^4)
# and product one.  Fixing the six-cycle leaves conjugation by its cyclic
# centralizer.
factorizations = {
    (sigma_zero, sigma_first, sigma_second)
    for sigma_zero in zero_monodromies
    for sigma_first in transpositions
    for sigma_second in transpositions
    if compose(
        compose(
            compose(
                infinity_monodromy,
                sigma_zero,
            ),
            sigma_first,
        ),
        sigma_second,
    )
    == IDENTITY
}
assert len(factorizations) == 12

centralizer = [
    tuple((index + shift) % DEGREE for index in range(DEGREE))
    for shift in range(DEGREE)
]


def canonical_factorization(
    factorization: tuple[Permutation, Permutation, Permutation],
) -> tuple[Permutation, Permutation, Permutation]:
    return min(
        tuple(conjugate(element, permutation) for permutation in factorization)
        for element in centralizer
    )


conjugacy_classes = {
    canonical_factorization(factorization)
    for factorization in factorizations
}
assert len(conjugacy_classes) == 2

# Forgetting the order of the two simple branch points uses the B_2 Hurwitz
# move (a,b) -> (b,b^{-1}ab).  It does not merge the two cover classes.
hurwitz_images = set()
for sigma_zero, sigma_first, sigma_second in conjugacy_classes:
    hurwitz_images.add(
        canonical_factorization(
            (
                sigma_zero,
                sigma_second,
                conjugate(inverse(sigma_second), sigma_first),
            )
        )
    )
assert hurwitz_images == conjugacy_classes

# Every polynomial cover with profiles (6) at infinity and (2,2,2) over zero
# is the square of a cubic.  Normalize its labelled roots to 0,1,a.
W, target_value, root_cross_ratio = sp.symbols("W target_value a")
cubic = W * (W - 1) * (W - root_cross_ratio)
critical_value_polynomial = sp.factor(
    sp.discriminant(cubic - target_value, W)
)
critical_value_coefficients = sp.Poly(
    critical_value_polynomial,
    target_value,
).all_coeffs()
assert len(critical_value_coefficients) == 3
leading, linear, constant = critical_value_coefficients
critical_sum = sp.factor(-linear / leading)
critical_product = sp.factor(constant / leading)

# If u,v are the two critical values of the cubic, the simple branch values
# of its square are u^2,v^2.  The following invariant forgets their order and
# common target scaling.
branch_invariant = sp.factor(
    (
        critical_sum**2 - 2 * critical_product
    ) ** 2
    / critical_product**2
)
branch_invariant_at_three = sp.factor(
    branch_invariant.subs(root_cross_ratio, 3)
)
fiber_numerator = sp.factor(
    sp.together(
        branch_invariant - branch_invariant_at_three
    ).as_numer_denom()[0]
)

label_orbit_factor = sp.factor(
    (root_cross_ratio - 3)
    * (root_cross_ratio + 2)
    * (2 * root_cross_ratio - 3)
    * (2 * root_cross_ratio + 1)
    * (3 * root_cross_ratio - 2)
    * (3 * root_cross_ratio - 1)
)
other_cover_factor = (
    9 * root_cross_ratio**6
    - 27 * root_cross_ratio**5
    + 79 * root_cross_ratio**4
    - 113 * root_cross_ratio**3
    + 79 * root_cross_ratio**2
    - 27 * root_cross_ratio
    + 9
)
assert sp.factor(
    fiber_numerator
    / (label_orbit_factor * other_cover_factor)
).is_Rational
assert sp.Poly(
    other_cover_factor,
    root_cross_ratio,
).is_irreducible

# The first six roots are exactly the S_3 orbit obtained by relabelling the
# three zero points.  Fixing the labelled source-root cross-ratio a=3 chooses
# one simple point of this fiber and excludes the second Hurwitz class.
source_relabelling_orbit = {
    sp.Integer(3),
    sp.Rational(1, 3),
    -sp.Integer(2),
    -sp.Rational(1, 2),
    sp.Rational(3, 2),
    sp.Rational(2, 3),
}
assert all(
    label_orbit_factor.subs(root_cross_ratio, value) == 0
    for value in source_relabelling_orbit
)
assert sp.diff(
    fiber_numerator,
    root_cross_ratio,
).subs(root_cross_ratio, 3) != 0
assert other_cover_factor.subs(root_cross_ratio, 3) != 0

artifact = {
    "experiment": "degree-six central Hurwitz component selection",
    "ramification_profiles": [
        [6],
        [2, 2, 2],
        [2, 1, 1, 1, 1],
        [2, 1, 1, 1, 1],
    ],
    "monodromy": {
        "factorizations_with_fixed_six_cycle": len(factorizations),
        "classes_with_ordered_simple_branch_points": len(conjugacy_classes),
        "classes_after_unordering_simple_branch_points": len(hurwitz_images),
    },
    "square_cubic_model": {
        "map": "(W(W-1)(W-a))^2",
        "unordered_scaled_branch_invariant": str(branch_invariant),
        "value_at_a=3": str(branch_invariant_at_three),
        "fiber_factor_at_a=3": str(fiber_numerator),
        "source_relabelling_orbit": [
            str(value)
            for value in sorted(source_relabelling_orbit, key=float)
        ],
        "other_cover_factor": str(other_cover_factor),
        "other_cover_factor_irreducible_over_Q": True,
        "a=3_is_reduced": True,
    },
    "conclusion": (
        "the ambient central Hurwitz problem has two cover classes, but the "
        "labelled source-root cross-ratio selects the polynomial class "
        "etale-locally at a=3"
    ),
    "scope": (
        "central vertex component selection in the labelled degree-six "
        "(2,2,2) chart; construction of the global comparison morphism, "
        "stack inertia, and label descent remain open"
    ),
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print("PASS central Hurwitz enumeration: two algebraic cover classes")
print("PASS unordering simple branch points: classes remain distinct")
print("PASS square-cubic branch invariant: exact degree-twelve fiber")
print("PASS source-root labels: six relabellings plus one other cover orbit")
print("PASS a=3: reduced and separated from the second cover class")
print("DEGREE_SIX_CENTRAL_HURWITZ_SELECTION_PASS")
