#!/usr/bin/env sage
"""Lift the q=80 rank-19 branch from its singular CM24 seed over GF(11).

The cubic-denominator presentation has a nine-dimensional Zariski tangent
space at the CM point.  Six directions arise from section cancellation and
singular component charts.  The quadratic obstruction calculation in
``verify_q80_cm24_mod7_seed.sage`` selects a tangent that moves the surface.

After that first singular step, lifting one digit while requiring that the
next digit also exist is a linear calculation.  This script performs that
look-ahead Hensel iteration exactly over the integers and verifies all 62
equations at every precision.
"""

import argparse

from sage.all import GF, Matrix, ZZ, vector


parser = argparse.ArgumentParser()
parser.add_argument("--digits", type=int, default=12)
arguments = parser.parse_args()
if arguments.digits < 3:
    parser.error("--digits must be at least 3")

# This loads the exact 54-variable system, its seed, Jacobian, kernels, and
# quadratic obstruction ideal.  It also independently prints the seed audit.
load("elkies-k3/scripts/verify_q80_cm24_mod7_seed.sage")

prime = ZZ(11)
field = GF(prime)
seed_vector = vector(ZZ, seed_values)
generators = parameters.gens()


def evaluate(candidate):
    substitution = dict(zip(generators, candidate))
    return vector(ZZ, [equation.subs(substitution) for equation in equations])


def divided_residual(candidate, modulus):
    values = evaluate(candidate)
    assert all(value % modulus == 0 for value in values)
    return vector(field, [field(-(value // modulus)) for value in values])


def compatible_correction(current, current_modulus):
    """Choose a correction that lifts now and remains liftable next digit."""

    rhs = divided_residual(current, current_modulus)
    obstruction = left_kernel * rhs
    assert not any(obstruction), tuple(obstruction)
    particular = jacobian.solve_right(rhs)
    assert jacobian * particular == rhs

    def next_obstruction(correction):
        trial = current + current_modulus * vector(
            ZZ, [ZZ(value) for value in correction]
        )
        next_rhs = divided_residual(trial, current_modulus * prime)
        return left_kernel * next_rhs

    base_obstruction = next_obstruction(particular)
    columns = []
    for tangent in kernel.rows():
        columns.append(next_obstruction(particular + tangent) - base_obstruction)
    compatibility = Matrix(field, columns).transpose()
    target = -base_obstruction
    assert target in compatibility.column_space(), (
        current_modulus, tuple(base_obstruction)
    )
    kernel_coordinates = compatibility.solve_right(target)
    correction = vector(field, particular)
    for coefficient, tangent in zip(kernel_coordinates, kernel.rows()):
        correction += coefficient * tangent
    assert not any(next_obstruction(correction))
    return correction, compatibility.rank(), compatibility.right_kernel().dimension()


# The radical of the quadratic obstruction is an affine five-space.  Its
# projection to the three surface tangent coordinates has only two degrees of
# freedom; the other three are cancellation gauges.  The following section
# c6=c7=c9=0 gives one representative for each of the 11^2 surface tangents.
# Search these representatives to the requested depth instead of trusting the
# first tangent that happens to cross the quadratic gate.
initial_candidates_tested = 0
candidate = None
initial_kernel_coordinates = None
for surface_s in field:
    for surface_c2 in field:
        c8 = surface_c2 + 3*surface_s
        c5 = surface_s + 3*c8
        point = vector(
            field,
            (
                (surface_s-3)/5,
                surface_c2,
                (-surface_s-1)/4,
                5,
                c5,
                0,
                0,
                c8,
                0,
            ),
        )
        substitution = dict(zip(tangent_field_variables, point))
        assert all(
            polynomial.subs(substitution) == 0
            for polynomial in obstruction_polynomials
        )
        first_correction = vector(field, first_delta)
        for coefficient, tangent in zip(point, kernel.rows()):
            first_correction += coefficient * tangent
        trial = seed_vector + prime * vector(
            ZZ, [ZZ(value) for value in first_correction]
        )
        trial_modulus = prime**2
        initial_candidates_tested += 1
        survived = True
        while trial_modulus < prime**arguments.digits:
            try:
                correction, _, _ = compatible_correction(trial, trial_modulus)
            except (AssertionError, ValueError):
                survived = False
                break
            trial += trial_modulus * vector(
                ZZ, [ZZ(value) for value in correction]
            )
            trial_modulus *= prime
        if survived:
            candidate = trial
            modulus = trial_modulus
            initial_kernel_coordinates = point
            break
    if candidate is not None:
        break

assert candidate is not None, (
    f"none of {initial_candidates_tested} obstruction representatives "
    f"survived to {arguments.digits} digits"
)
print(
    f"Q80RANK19LIFT|stage=branch_search|tested={initial_candidates_tested}|"
    f"digits={arguments.digits}|kernel_coordinates="
    + ",".join(map(str, map(int, initial_kernel_coordinates))),
    flush=True,
)


while modulus < prime**arguments.digits:
    correction, compatibility_rank, compatibility_nullity = compatible_correction(
        candidate, modulus
    )
    candidate += modulus * vector(ZZ, [ZZ(value) for value in correction])
    modulus *= prime
    values = evaluate(candidate)
    assert all(value % modulus == 0 for value in values)
    print(
        f"Q80RANK19LIFT|digits={modulus.valuation(prime)}|modulus={modulus}|"
        f"compatibility_rank={compatibility_rank}|"
        f"compatibility_nullity={compatibility_nullity}",
        flush=True,
    )

surface_names = names[:11]
surface_values = tuple(ZZ(value % modulus) for value in candidate[:11])
final_values = evaluate(candidate)
minimum_valuation = min(
    value.valuation(prime) if value else arguments.digits + 1
    for value in final_values
)
print(
    "Q80RANK19LIFT|surface="
    + ",".join(
        f"{name}:{value}" for name, value in zip(surface_names, surface_values)
    ),
    flush=True,
)
print(
    f"Q80RANK19LIFT|digits={arguments.digits}|"
    f"minimum_residual_valuation={minimum_valuation}|status=PASS",
    flush=True,
)
