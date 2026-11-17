#!/usr/bin/env python3
"""Bounded mixed canonical-pivot search for HC(4).

This is a finite search, not a proof about all polynomial canonical charts.
It starts from the collision-centred foundational three-variable Keller
doubling and applies exact time-one Hamiltonian shears from three families:

* one mixed linear form, H=tau*(q_i+eps*p_j)^d, d=2 or 3;
* two commuting mixed forms, with a genuinely two-form monomial of degree
  two or three;
* a cubic shear H=tau*q_k*(q_i+eps*p_j)^2 whose image of p_k=0 is a
  coisotropic graph with one nonlinear constraint.

Every generator contains both source and dual variables.  No pure source
transformation is searched.

For each transformed polynomial the checker follows the requested order:

1. enumerate scalar and simultaneous coordinate blocks in which the
   polynomial is jointly affine;
2. test D(mu+lambda*A,w) for every scalar affine pivot and the small repair
   box lambda in {-1,1}, mu in {-1,0,1};
3. test the simultaneous r=2 corank budget rank(M)<=2;
4. test the complete descended four-variable determinant for every
   invertible symmetric Lambda with entries in {-1,0,1} and every
   mu in {-1,0,1}^2.

Fast rejection uses exact reduction modulo a stated good prime.  Unequal
values at two points certify nonconstancy over Q.  Agreement on all sampled
points is reported only as a modular survivor; it is not promoted to a
characteristic-zero identity.

The transformed six-variable Hessian determinant is also audited.  This is
logically important: a nonlinear symplectic change of independent variables
does not in general preserve the constant-Hessian equation.  Consequently
the specialized D test is only the exact Schur gate when the parent
determinant remains constant; the complete descended determinant is checked
for every affine pair regardless.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp

from canonical_transform_search import (
    canonical_poisson_matrix,
    exact_invariant_shear,
    mixed_line_hamiltonians,
    poisson_bracket as shared_poisson_bracket,
    pullback_polynomial,
)


PRIME = 1_000_003
TAU_BOX = (-2, -1, 1, 2)
SCALAR_LAMBDAS = (-1, 1)
SCALAR_MUS = (-1, 0, 1)
PAIR_MUS = tuple(product((-1, 0, 1), repeat=2))
PAIR_LAMBDAS = tuple(
    (a, b, c)
    for a, b, c in product((-1, 0, 1), repeat=3)
    if a * c - b * b != 0
)


q = sp.symbols("x y z")
p = sp.symbols("u v w")
variables = q + p
x, y, z = q
u, v, w = p
POISSON = canonical_poisson_matrix(3, q_p_bracket=1)


unit = 1 + x * y
seed_q = unit**2 * z + y**2 * (1 + 3 * unit)
P = sp.expand(unit * seed_q)
B = sp.expand(y + 3 * x * seed_q)
C = sp.expand(x * (5 - 3 * unit) - x**3 * z)

# The two source points (1,-3/2,13/2) and (-1,3/2,13/2) have common image
# (-1/4,0,0).  Clearing the harmless factor 1/4 makes the centred doubling
# integral.  Both lifted points with p=0 are critical points.
base_potential = sp.expand(u * (4 * P + 1) + 4 * v * B + 4 * w * C)
BASE_HESSIAN_DETERMINANT = -4 * 4**6


@dataclass(frozen=True)
class Chart:
    chart_id: str
    family: str
    tau: int
    hamiltonian: sp.Expr
    metadata: dict[str, object]


def poisson_bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    return shared_poisson_bracket(
        left,
        right,
        variables,
        POISSON,
    )


def charts() -> Iterable[Chart]:
    """Yield the complete declared finite chart box without duplicates."""

    for letter in mixed_line_hamiltonians(
        q,
        p,
        degrees=(2, 3),
        coefficients=TAU_BOX,
    ):
        source_index = letter.source_index
        dual_index = letter.dual_index
        epsilon = letter.epsilon
        degree = letter.degree
        tau = letter.coefficient
        hamiltonian = letter.hamiltonian
        yield Chart(
            chart_id=(
                f"one-q{source_index}-p{dual_index}-e{epsilon:+d}"
                f"-d{degree}-t{tau:+d}"
            ),
            family="one_mixed_line",
            tau=tau,
            hamiltonian=hamiltonian,
            metadata={
                "source_index": source_index,
                "dual_index": dual_index,
                "epsilon": epsilon,
                "degree": degree,
            },
        )

    for left, right in combinations(range(3), 2):
        for epsilon, exponents, tau in product(
            (-1, 1), ((1, 1), (1, 2), (2, 1)), TAU_BOX
        ):
            linear_left = q[left] + epsilon * p[right]
            linear_right = q[right] + epsilon * p[left]
            assert poisson_bracket(linear_left, linear_right) == 0
            a, b = exponents
            hamiltonian = tau * linear_left**a * linear_right**b
            yield Chart(
                chart_id=(
                    f"two-q{left}p{right}-q{right}p{left}"
                    f"-e{epsilon:+d}-a{a}b{b}-t{tau:+d}"
                ),
                family="two_mixed_lines",
                tau=tau,
                hamiltonian=hamiltonian,
                metadata={
                    "pair": (left, right),
                    "epsilon": epsilon,
                    "exponents": exponents,
                },
            )

    # Since q_k Poisson-commutes with the mixed linear form, this cubic
    # Hamiltonian again has a finite exact flow.  Its p_k component is
    # translated by a nonzero multiple of linear^2, so the image of the
    # coordinate coisotropic p_k=0 has one nonlinear graph constraint.
    for characteristic_index in range(3):
        remaining = tuple(
            index for index in range(3) if index != characteristic_index
        )
        for source_index, dual_index, epsilon, tau in product(
            remaining, remaining, (-1, 1), TAU_BOX
        ):
            linear = q[source_index] + epsilon * p[dual_index]
            assert poisson_bracket(q[characteristic_index], linear) == 0
            hamiltonian = tau * q[characteristic_index] * linear**2
            yield Chart(
                chart_id=(
                    f"coiso-k{characteristic_index}-q{source_index}"
                    f"-p{dual_index}-e{epsilon:+d}-t{tau:+d}"
                ),
                family="coisotropic_graph",
                tau=tau,
                hamiltonian=hamiltonian,
                metadata={
                    "characteristic_index": characteristic_index,
                    "source_index": source_index,
                    "dual_index": dual_index,
                    "epsilon": epsilon,
                    "degree": 3,
                },
            )


def time_one_pullback(hamiltonian: sp.Expr) -> sp.Expr:
    """Return base_potential after the exact time-one Hamiltonian shear."""

    coordinate_map = exact_invariant_shear(
        hamiltonian,
        variables,
        POISSON,
    )
    return pullback_polynomial(base_potential, coordinate_map, variables)


def deterministic_points(dimension: int, count: int = 15) -> tuple[tuple[int, ...], ...]:
    """Independent reproducible modular points for fast witnesses."""

    generator = random.Random(0x48433400 + dimension)
    points = [tuple(0 for _ in range(dimension))]
    while len(points) < count:
        point = tuple(generator.randrange(1, PRIME) for _ in range(dimension))
        if point not in points:
            points.append(point)
    return tuple(points)


class ModularPolynomial:
    """Sparse multivariate polynomial evaluator modulo PRIME."""

    def __init__(self, polynomial: sp.Poly):
        self.dimension = len(polynomial.gens)
        self.terms = tuple(
            (monomial, int(coefficient) % PRIME)
            for monomial, coefficient in polynomial.terms()
        )

    def evaluate(self, point: Sequence[int]) -> int:
        assert len(point) == self.dimension
        value = 0
        for monomial, coefficient in self.terms:
            term = coefficient
            for coordinate, exponent in zip(point, monomial, strict=True):
                if exponent:
                    term = term * pow(int(coordinate), exponent, PRIME) % PRIME
            value = (value + term) % PRIME
        return value


@dataclass
class Jet:
    value: ModularPolynomial
    gradient: tuple[ModularPolynomial, ...]
    hessian: tuple[tuple[ModularPolynomial, ...], ...]

    @classmethod
    def from_poly(cls, polynomial: sp.Poly) -> "Jet":
        gens = polynomial.gens
        gradient_polys = tuple(polynomial.diff(gen) for gen in gens)
        hessian_polys = tuple(
            tuple(gradient_poly.diff(gen) for gen in gens)
            for gradient_poly in gradient_polys
        )
        return cls(
            value=ModularPolynomial(polynomial),
            gradient=tuple(ModularPolynomial(poly) for poly in gradient_polys),
            hessian=tuple(
                tuple(ModularPolynomial(poly) for poly in row)
                for row in hessian_polys
            ),
        )

    def evaluate(
        self, point: Sequence[int]
    ) -> tuple[int, list[int], list[list[int]]]:
        return (
            self.value.evaluate(point),
            [entry.evaluate(point) for entry in self.gradient],
            [
                [entry.evaluate(point) for entry in row]
                for row in self.hessian
            ],
        )


def determinant_mod(matrix: Sequence[Sequence[int]]) -> int:
    work = [[int(value) % PRIME for value in row] for row in matrix]
    size = len(work)
    determinant = 1
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant = determinant * pivot_value % PRIME
        inverse = pow(pivot_value, PRIME - 2, PRIME)
        for row in range(column + 1, size):
            multiplier = work[row][column] * inverse % PRIME
            if not multiplier:
                continue
            for inner in range(column, size):
                work[row][inner] = (
                    work[row][inner] - multiplier * work[column][inner]
                ) % PRIME
    return determinant % PRIME


def rank_mod(matrix: Sequence[Sequence[int]]) -> int:
    work = [[int(value) % PRIME for value in row] for row in matrix]
    if not work:
        return 0
    rows = len(work)
    columns = len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], PRIME - 2, PRIME)
        work[pivot_row] = [
            value * inverse % PRIME for value in work[pivot_row]
        ]
        for row in range(rows):
            if row == pivot_row:
                continue
            multiplier = work[row][column]
            if multiplier:
                work[row] = [
                    (left - multiplier * right) % PRIME
                    for left, right in zip(
                        work[row], work[pivot_row], strict=True
                    )
                ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def affine_coordinates(polynomial: sp.Poly) -> tuple[int, ...]:
    return tuple(
        index
        for index in range(len(variables))
        if polynomial.degree(variables[index]) <= 1
    )


def affine_pairs(
    polynomial: sp.Poly, affine: Sequence[int]
) -> tuple[tuple[int, int], ...]:
    pairs: list[tuple[int, int]] = []
    for left, right in combinations(affine, 2):
        if all(
            not (monomial[left] and monomial[right])
            for monomial, _ in polynomial.terms()
        ):
            pairs.append((left, right))
    return tuple(pairs)


def split_scalar(
    expression: sp.Expr, pivot_index: int
) -> tuple[sp.Poly, sp.Poly, tuple[sp.Symbol, ...]]:
    pivot = variables[pivot_index]
    retained = tuple(
        variable
        for index, variable in enumerate(variables)
        if index != pivot_index
    )
    A = sp.Poly(sp.expand(expression).coeff(pivot, 1), *retained, domain=sp.ZZ)
    B0 = sp.Poly(
        sp.expand(expression).coeff(pivot, 0), *retained, domain=sp.ZZ
    )
    return A, B0, retained


def split_pair(
    expression: sp.Expr, pair: tuple[int, int]
) -> tuple[sp.Poly, sp.Poly, sp.Poly, tuple[sp.Symbol, ...]]:
    left, right = (variables[index] for index in pair)
    retained = tuple(
        variable
        for index, variable in enumerate(variables)
        if index not in pair
    )
    expanded = sp.expand(expression)
    A_left = sp.Poly(expanded.coeff(left, 1), *retained, domain=sp.ZZ)
    A_right = sp.Poly(expanded.coeff(right, 1), *retained, domain=sp.ZZ)
    B0 = sp.Poly(
        expanded.coeff(left, 0).coeff(right, 0),
        *retained,
        domain=sp.ZZ,
    )
    return A_left, A_right, B0, retained


def first_difference(values: Sequence[int]) -> dict[str, int] | None:
    first = values[0]
    for index, value in enumerate(values[1:], start=1):
        if value != first:
            return {
                "point_a": 0,
                "value_a": first,
                "point_b": index,
                "value_b": value,
            }
    return None


def matrix_add(
    *matrices: Sequence[Sequence[int]],
) -> list[list[int]]:
    size = len(matrices[0])
    return [
        [
            sum(matrix[row][column] for matrix in matrices) % PRIME
            for column in range(size)
        ]
        for row in range(size)
    ]


def matrix_scale(
    scalar: int, matrix: Sequence[Sequence[int]]
) -> list[list[int]]:
    return [
        [scalar * value % PRIME for value in row]
        for row in matrix
    ]


def outer(left: Sequence[int], right: Sequence[int]) -> list[list[int]]:
    return [
        [left_value * right_value % PRIME for right_value in right]
        for left_value in left
    ]


def scalar_remainder_values(
    evaluations: Sequence[
        tuple[
            tuple[int, list[int], list[list[int]]],
            tuple[int, list[int], list[list[int]]],
        ]
    ],
    lam: int,
    mu: int,
) -> list[int]:
    values: list[int] = []
    for A_evaluation, B_evaluation in evaluations:
        A_value, _, A_hessian = A_evaluation
        _, _, B_hessian = B_evaluation
        scalar = (mu + lam * A_value) % PRIME
        reduced = matrix_add(
            B_hessian, matrix_scale(scalar, A_hessian)
        )
        values.append(determinant_mod(reduced))
    return values


def pair_repair_values(
    evaluations: Sequence[
        tuple[
            tuple[int, list[int], list[list[int]]],
            tuple[int, list[int], list[list[int]]],
            tuple[int, list[int], list[list[int]]],
        ]
    ],
    repair: tuple[int, int, int],
    mu: tuple[int, int],
) -> list[int]:
    lam11, lam12, lam22 = repair
    values: list[int] = []
    for A1_evaluation, A2_evaluation, B_evaluation in evaluations:
        a1, g1, h1 = A1_evaluation
        a2, g2, h2 = A2_evaluation
        _, _, hb = B_evaluation
        matrices = [
            hb,
            matrix_scale(mu[0], h1),
            matrix_scale(mu[1], h2),
            matrix_scale(lam11 * a1, h1),
            matrix_scale(lam22 * a2, h2),
            matrix_scale(lam11, outer(g1, g1)),
            matrix_scale(lam22, outer(g2, g2)),
            matrix_scale(lam12 * a2, h1),
            matrix_scale(lam12 * a1, h2),
            matrix_scale(lam12, outer(g1, g2)),
            matrix_scale(lam12, outer(g2, g1)),
        ]
        values.append(determinant_mod(matrix_add(*matrices)))
    return values


def exact_constant(expression: sp.Expr, gens: Sequence[sp.Symbol]) -> bool:
    polynomial = sp.Poly(sp.expand(expression), *gens)
    return all(sum(monomial) == 0 for monomial, _ in polynomial.terms())


def exact_scalar_remainder(
    A: sp.Poly, B0: sp.Poly, lam: int, mu: int
) -> sp.Expr:
    gens = A.gens
    scalar = mu + lam * A.as_expr()
    matrix = sp.hessian(B0.as_expr(), gens) + scalar * sp.hessian(
        A.as_expr(), gens
    )
    return sp.factor(matrix.det(method="domain-ge"))


def exact_corank(
    A1: sp.Poly, A2: sp.Poly, B0: sp.Poly
) -> bool:
    gens = A1.gens
    s1, s2 = sp.symbols("s1 s2")
    matrix = (
        sp.hessian(B0.as_expr(), gens)
        + s1 * sp.hessian(A1.as_expr(), gens)
        + s2 * sp.hessian(A2.as_expr(), gens)
    )
    for rows in combinations(range(4), 3):
        for columns in combinations(range(4), 3):
            minor = matrix.extract(rows, columns).det(method="berkowitz")
            if sp.expand(minor) != 0:
                return False
    return True


def exact_pair_determinant(
    A1: sp.Poly,
    A2: sp.Poly,
    B0: sp.Poly,
    repair: tuple[int, int, int],
    mu: tuple[int, int],
) -> sp.Expr:
    lam11, lam12, lam22 = repair
    a1 = A1.as_expr()
    a2 = A2.as_expr()
    psi = sp.expand(
        B0.as_expr()
        + mu[0] * a1
        + mu[1] * a2
        + sp.Rational(lam11, 2) * a1**2
        + lam12 * a1 * a2
        + sp.Rational(lam22, 2) * a2**2
    )
    return sp.factor(
        sp.hessian(psi, A1.gens).det(method="domain-ge")
    )


def audit_parent_hessian(
    expression: sp.Expr,
    points: Sequence[Sequence[int]],
    *,
    quadratic_generator: bool,
) -> tuple[bool | None, dict[str, int] | None]:
    if quadratic_generator:
        # A quadratic Hamiltonian has a linear symplectic time-one map.
        # Hessians transform by constant congruence, and the symplectic
        # determinant is one.
        return True, None

    polynomial = sp.Poly(expression, *variables, domain=sp.ZZ)
    hessian = tuple(
        tuple(
            ModularPolynomial(polynomial.diff(left).diff(right))
            for right in variables
        )
        for left in variables
    )
    values = [
        determinant_mod(
            [
                [entry.evaluate(point) for entry in row]
                for row in hessian
            ]
        )
        for point in points
    ]
    witness = first_difference(values)
    if witness is not None:
        return False, witness
    # Parent constancy is a diagnostic, not a search gate.  Do not promote
    # agreement on finitely many points to a theorem, and do not trigger a
    # potentially enormous six-variable determinant expansion here.
    return None, None


def run_search() -> dict[str, object]:
    base_determinant = sp.expand(
        sp.hessian(base_potential, variables).det(method="berkowitz")
    )
    assert base_determinant == BASE_HESSIAN_DETERMINANT

    chart_rows: list[dict[str, object]] = []
    summary = {
        "transformations": 0,
        "transformations_by_family": {},
        "charts_with_affine_pivot": 0,
        "parent_constant_hessian": 0,
        "scalar_affine_pivots": 0,
        "scalar_remainder_trials": 0,
        "scalar_remainder_exact_survivors": 0,
        "scalar_remainder_modular_survivors": 0,
        "simultaneous_affine_pairs": 0,
        "corank_budget_modular_survivors": 0,
        "complete_determinant_trials": 0,
        "complete_determinant_modular_survivors": 0,
    }

    parent_points = deterministic_points(6)
    scalar_points = deterministic_points(5)
    pair_points = deterministic_points(4)

    all_charts = tuple(charts())
    for chart_index, chart in enumerate(all_charts, start=1):
        if chart_index == 1 or chart_index % 10 == 0:
            print(
                f"progress={chart_index}/{len(all_charts)} "
                f"chart={chart.chart_id}",
                flush=True,
            )
        summary["transformations"] += 1
        by_family = summary["transformations_by_family"]
        assert isinstance(by_family, dict)
        by_family[chart.family] = by_family.get(chart.family, 0) + 1

        transformed = time_one_pullback(chart.hamiltonian)
        transformed_poly = sp.Poly(
            transformed, *variables, domain=sp.ZZ
        )
        affine = affine_coordinates(transformed_poly)
        pairs = affine_pairs(transformed_poly, affine)
        if affine:
            summary["charts_with_affine_pivot"] += 1
        summary["scalar_affine_pivots"] += len(affine)
        summary["simultaneous_affine_pairs"] += len(pairs)

        row: dict[str, object] = {
            "chart_id": chart.chart_id,
            "family": chart.family,
            "metadata": chart.metadata,
            "tau": chart.tau,
            "terms": len(transformed_poly.terms()),
            "affine_coordinates": [
                str(variables[index]) for index in affine
            ],
            "affine_pairs": [
                [str(variables[left]), str(variables[right])]
                for left, right in pairs
            ],
        }

        remainder_survivors: list[dict[str, object]] = []
        remainder_witnesses: list[dict[str, object]] = []
        for pivot_index in affine:
            A, B0, _ = split_scalar(transformed, pivot_index)
            jets = (Jet.from_poly(A), Jet.from_poly(B0))
            scalar_evaluations = tuple(
                (jets[0].evaluate(point), jets[1].evaluate(point))
                for point in scalar_points
            )
            for lam, mu in product(SCALAR_LAMBDAS, SCALAR_MUS):
                summary["scalar_remainder_trials"] += 1
                values = scalar_remainder_values(
                    scalar_evaluations, lam, mu
                )
                witness = first_difference(values)
                trial = {
                    "pivot": str(variables[pivot_index]),
                    "lambda": lam,
                    "mu": mu,
                }
                if witness is not None:
                    if len(remainder_witnesses) < 3:
                        remainder_witnesses.append(trial | witness)
                    continue
                structural_zero = (
                    variables[pivot_index] == z
                    and sp.Poly(
                        chart.hamiltonian, *variables
                    ).total_degree()
                    == 2
                    and not chart.hamiltonian.has(z, w)
                )
                if structural_zero:
                    # Before the linear symplectic change on (x,y,u,v),
                    # B+sA is a five-variable doubled potential with two
                    # source and three dual variables, so its reduced
                    # Hessian is singular by the standard row count.
                    # Constant congruence by the four-variable symplectic
                    # block preserves its zero determinant.
                    summary["scalar_remainder_exact_survivors"] += 1
                    remainder_survivors.append(
                        trial
                        | {
                            "constant": "0",
                            "proof": (
                                "doubled row-count plus constant "
                                "symplectic congruence"
                            ),
                        }
                    )
                else:
                    summary[
                        "scalar_remainder_modular_survivors"
                    ] += 1
                    remainder_survivors.append(
                        trial | {"values_mod_p": sorted(set(values))}
                    )
        if remainder_survivors:
            row["scalar_remainder_survivors"] = remainder_survivors
        if remainder_witnesses:
            row["sample_scalar_remainder_witnesses_mod_p"] = (
                remainder_witnesses
            )

        corank_survivors: list[list[str]] = []
        determinant_survivors: list[dict[str, object]] = []
        determinant_witnesses: list[dict[str, object]] = []
        for pair in pairs:
            A1, A2, B0, _ = split_pair(transformed, pair)
            jets = (Jet.from_poly(A1), Jet.from_poly(A2), Jet.from_poly(B0))
            pair_evaluations = tuple(
                (
                    jets[0].evaluate(point),
                    jets[1].evaluate(point),
                    jets[2].evaluate(point),
                )
                for point in pair_points
            )

            corank_modular = True
            for point_index, evaluations in enumerate(pair_evaluations):
                _, _, h1 = evaluations[0]
                _, _, h2 = evaluations[1]
                _, _, hb = evaluations[2]
                # Independent pencil values, not the repair graph.
                s1 = (point_index + 2) % PRIME
                s2 = (3 * point_index + 1) % PRIME
                reduced = matrix_add(
                    hb,
                    matrix_scale(s1, h1),
                    matrix_scale(s2, h2),
                )
                if rank_mod(reduced) > 2:
                    corank_modular = False
                    break
            if corank_modular:
                summary["corank_budget_modular_survivors"] += 1
                corank_survivors.append(
                    [str(variables[pair[0]]), str(variables[pair[1]])]
                )

            for repair, mu in product(PAIR_LAMBDAS, PAIR_MUS):
                summary["complete_determinant_trials"] += 1
                values = pair_repair_values(
                    pair_evaluations, repair, mu
                )
                witness = first_difference(values)
                trial = {
                    "pair": [
                        str(variables[pair[0]]),
                        str(variables[pair[1]]),
                    ],
                    "lambda": repair,
                    "mu": mu,
                }
                if witness is not None:
                    if len(determinant_witnesses) < 3:
                        determinant_witnesses.append(trial | witness)
                    continue
                summary["complete_determinant_modular_survivors"] += 1
                determinant_survivors.append(
                    trial | {"values_mod_p": sorted(set(values))}
                )

        if corank_survivors:
            row["corank_budget_modular_survivors"] = corank_survivors
        if determinant_survivors:
            row["complete_determinant_modular_survivors"] = (
                determinant_survivors
            )
        if determinant_witnesses:
            row["sample_complete_determinant_witnesses_mod_p"] = (
                determinant_witnesses
            )

        # This side audit is deliberately performed after the four search
        # gates above.  It does not prune a candidate or change their order.
        parent_constant, parent_witness = audit_parent_hessian(
            transformed,
            parent_points,
            quadratic_generator=sp.Poly(
                chart.hamiltonian, *variables
            ).total_degree()
            == 2,
        )
        if parent_constant:
            summary["parent_constant_hessian"] += 1
        row["parent_constant_hessian"] = (
            parent_constant
            if parent_constant is not None
            else "unresolved_after_modular_points"
        )
        if parent_witness is not None:
            row["parent_hessian_witness_mod_p"] = parent_witness
        chart_rows.append(row)

    return {
        "status": "bounded_search",
        "scope": {
            "base": "collision-centred foundational cubic Keller doubling",
            "prime": PRIME,
            "tau_box": TAU_BOX,
            "scalar_lambda_box": SCALAR_LAMBDAS,
            "scalar_mu_box": SCALAR_MUS,
            "pair_lambda_box": PAIR_LAMBDAS,
            "pair_mu_box": PAIR_MUS,
            "pure_source_transformations": "excluded by construction",
        },
        "summary": summary,
        "charts": chart_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="write the complete JSON census to this path",
    )
    args = parser.parse_args()

    result = run_search()
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n"
        )

    print("HC4_MIXED_CANONICAL_SUMMARY")
    print(json.dumps(result["summary"], indent=2, sort_keys=True))
    if args.output is not None:
        print(f"artifact={args.output}")
    print(
        "SCOPE: finite exact Hamiltonian-shear box only; no theorem about "
        "all mixed canonical or coisotropic charts"
    )


if __name__ == "__main__":
    main()
