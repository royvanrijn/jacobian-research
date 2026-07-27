#!/usr/bin/env python3
"""Exact regression for the two cubic cusp finite-flat local models."""

import itertools
import sympy as sp


t, u, v, T, Z = sp.symbols("t u v T Z")

# The source is k[t,u], finite free of rank three over k[u,v] because t is
# integral with monic equation T^3 + uT - v.
minimal_equation = T**3 + u * T - v
assert sp.Poly(minimal_equation, T).monic()
assert sp.Poly(minimal_equation, T).degree() == 3

source_v = t**3 + u * t
jacobian = sp.diff(source_v, t)
assert sp.expand(jacobian) == u + 3 * t**2

# The discriminant cusp and its full pullback.
discriminant = 4 * u**3 + 27 * v**2
pulled_discriminant = sp.factor(discriminant.subs(v, source_v))
assert pulled_discriminant == (u + 3 * t**2) ** 2 * (4 * u + 3 * t**2)

# The affine companion meets the ramification boundary with contact two.
boundary_equation = u + 3 * t**2
companion_equation = 4 * u + 3 * t**2
companion_on_boundary = sp.expand(companion_equation.subs(u, -3 * t**2))
assert companion_on_boundary == -9 * t**2
assert sp.Poly(companion_on_boundary, t).degree() == 2

# On the ramification curve u=-3t^2, the image is a cusp and the derivative
# of its normalization parametrization vanishes at the origin.
cusp_u = -3 * t**2
cusp_v = sp.expand(source_v.subs(u, cusp_u))
assert cusp_v == -2 * t**3
assert sp.expand(discriminant.subs({u: cusp_u, v: cusp_v})) == 0
assert sp.diff(cusp_u, t).subs(t, 0) == 0
assert sp.diff(cusp_v, t).subs(t, 0) == 0

# The special fiber over the cusp has coordinate algebra k[t]/(t^3), hence
# length three: it is flat, curvilinear, and not a packet-length defect.
special_fiber = sp.Poly(minimal_equation.subs({u: 0, v: 0}), T)
assert special_fiber.as_expr() == T**3
assert special_fiber.degree() == 3


# The ordinary-cusp complement has braid presentation
# <sigma_1,sigma_2 | sigma_1 sigma_2 sigma_1 =
#                         sigma_2 sigma_1 sigma_2>.
# Generic simple ramification sends both meridians to transpositions in S3.
# Up to relabeling, equal transpositions give the 2+1 Kummer cover and
# distinct transpositions give the transitive S3 root cover.
def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def generated_group(
    generators: tuple[tuple[int, ...], ...],
) -> set[tuple[int, ...]]:
    identity = tuple(range(3))
    group = {identity}
    frontier = [identity]
    while frontier:
        element = frontier.pop()
        for generator in generators:
            product = compose(generator, element)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return group


def inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * len(permutation)
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def conjugate(
    permutation: tuple[int, ...], relabeling: tuple[int, ...]
) -> tuple[int, ...]:
    return compose(
        relabeling,
        compose(permutation, inverse(relabeling)),
    )


transpositions = []
for left, right in itertools.combinations(range(3), 2):
    permutation = list(range(3))
    permutation[left], permutation[right] = (
        permutation[right],
        permutation[left],
    )
    transpositions.append(tuple(permutation))

braid_pairs = []
for sigma_1 in transpositions:
    for sigma_2 in transpositions:
        if compose(sigma_1, compose(sigma_2, sigma_1)) == compose(
            sigma_2, compose(sigma_1, sigma_2)
        ):
            braid_pairs.append((sigma_1, sigma_2))

assert len(braid_pairs) == 9
assert sum(left == right for left, right in braid_pairs) == 3
assert {
    len(generated_group((left, right)))
    for left, right in braid_pairs
    if left == right
} == {2}
assert {
    len(generated_group((left, right)))
    for left, right in braid_pairs
    if left != right
} == {6}
relabelings = tuple(itertools.permutations(range(3)))
remaining_pairs = set(braid_pairs)
pair_orbits = []
while remaining_pairs:
    seed = next(iter(remaining_pairs))
    orbit = {
        (conjugate(seed[0], relabeling), conjugate(seed[1], relabeling))
        for relabeling in relabelings
    }
    pair_orbits.append(orbit)
    remaining_pairs -= orbit
assert sorted(len(orbit) for orbit in pair_orbits) == [3, 6]
assert {
    tuple(sorted(left == right for left, right in orbit))
    for orbit in pair_orbits
} == {(True, True, True), (False,) * 6}

# The 2+1 representation is the direct sum of a quadratic Kummer algebra
# and the trivial sheet.  The quadratic equation is monic, hence free; its
# A2 form XY=u^3 is normal.
kummer_equation = Z**2 - discriminant
assert sp.Poly(kummer_equation, Z).monic()
assert sp.Poly(kummer_equation, Z).degree() == 2
alpha = sp.symbols("alpha")
kummer_a2_relation = sp.expand(
    (Z + alpha * v) * (Z - alpha * v) - 4 * u**3
)
assert sp.expand(
    kummer_a2_relation.subs(alpha**2, 27) - kummer_equation
) == 0

print("cubic cusp local models and braid representations: exact checks passed")
