#!/usr/bin/env python3
"""Search the full degree-six foundational equivariant box.

The source weights are (1,-1,-2), the target weights are (-2,-1,1), and
the linear part is normalized to (z,y,x).  On x != 0 write

    F = (x^-2 A(v,s), x^-1 B(v,s), x C(v,s)),
    v = x*y, s = x^2*z.

The Keller equation is the determinant of the three-by-three invariant
matrix below.  The order-two orbit of (1,1,1) is
(-1,-1,1), so B(1,1)=C(1,1)=0 is an exact two-point collision.

This script makes two deliberately separate claims.

* It proves over QQ, by an msolve Groebner basis [1], that the complete
  degree-six z-linear subbox contains no such collision.
* It performs a reproducible modular search in the remaining full box.
  The modular pass is a search, not an exhaustion of the coefficient
  scheme.
"""

from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path
import random
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.msolve import run as run_msolve  # noqa: E402


v, s = sp.symbols("v s")


def invariant_support(
    output_weight: int, degree: int
) -> tuple[tuple[int, int], ...]:
    """Return (v,s)-exponents of all allowed source monomials."""

    answer = set()
    for i, j, k in product(range(degree + 1), repeat=3):
        if not 1 <= i + j + k <= degree:
            continue
        if i - j - 2 * k == output_weight:
            answer.add((j, k))
    return tuple(sorted(answer))


def invariant_determinant(A: sp.Expr, B: sp.Expr, C: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.det(
            sp.Matrix(
                [
                    (-2 * A, sp.diff(A, v), sp.diff(A, s)),
                    (-B, sp.diff(B, v), sp.diff(B, s)),
                    (C, sp.diff(C, v), sp.diff(C, s)),
                ]
            )
        )
    )


def coefficient_equations(expression: sp.Expr) -> tuple[sp.Expr, ...]:
    return tuple(
        coefficient
        for _exponent, coefficient in sp.Poly(expression, v, s).terms()
    )


def exact_z_linear_exclusion() -> None:
    """Prove that the complete z-linear degree-six collision ideal is (1)."""

    A = (
        s
        + sp.Symbol("a20") * v**2
        + sp.Symbol("a11") * v * s
        + sp.Symbol("a30") * v**3
        + sp.Symbol("a21") * v**2 * s
        + sp.Symbol("a40") * v**4
    )
    B = (
        v
        + sp.Symbol("b01") * s
        + sp.Symbol("b20") * v**2
        + sp.Symbol("b11") * v * s
        + sp.Symbol("b30") * v**3
        + sp.Symbol("b21") * v**2 * s
    )
    C = (
        1
        + sp.Symbol("c10") * v
        + sp.Symbol("c01") * s
        + sp.Symbol("c20") * v**2
        + sp.Symbol("c11") * v * s
    )
    variables = tuple(
        sorted(
            (A.free_symbols | B.free_symbols | C.free_symbols) - {v, s},
            key=str,
        )
    )
    equations = (
        *coefficient_equations(invariant_determinant(A, B, C) + 1),
        sp.expand(B.subs({v: 1, s: 1})),
        sp.expand(C.subs({v: 1, s: 1})),
    )
    result = run_msolve(
        equations,
        variables,
        prime=0,
        threads=4,
        groebner=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
    assert result.contains_one, result.output


def exact_full_open_exclusion() -> None:
    """Exclude the three principal opens missing from the z-linear box."""

    a = sp.symbols("a0:7")
    b = sp.symbols("b0:6")
    c = sp.symbols("c0:4")
    inverse = sp.Symbol("inverse")
    A = (
        s
        + a[0] * v**2
        + a[1] * v * s
        + a[2] * v**3
        + a[3] * s**2
        + a[4] * v**2 * s
        + a[5] * v**4
        + a[6] * v * s**2
    )
    B = (
        v
        + b[0] * s
        + b[1] * v**2
        + b[2] * v * s
        + b[3] * v**3
        + b[4] * s**2
        + b[5] * v**2 * s
    )
    C = 1 + c[0] * v + c[1] * s + c[2] * v**2 + c[3] * v * s
    variables = (*a, *b, *c, inverse)
    equations = (
        *coefficient_equations(invariant_determinant(A, B, C) + 1),
        sp.expand(B.subs({v: 1, s: 1})),
        sp.expand(C.subs({v: 1, s: 1})),
    )

    # These three opens cover the complement of the z-linear closed box.
    for extra_coefficient in (a[3], a[6], b[4]):
        result = run_msolve(
            (*equations, inverse * extra_coefficient - 1),
            variables,
            prime=0,
            threads=4,
            groebner=True,
            timeout=360,
        )
        assert result.returncode == 0, result.stderr
        assert result.contains_one, result.output


def expanded_profile(
    mapping: tuple[sp.Expr, ...],
    source_variables: tuple[sp.Symbol, ...],
) -> tuple[int, int]:
    expanded = tuple(
        sp.Poly(sp.expand(component), *source_variables)
        for component in mapping
    )
    degree = max(component.total_degree() for component in expanded)
    nonconstant_terms = sum(
        sum(1 for monomial, _coefficient in component.terms() if any(monomial))
        for component in expanded
    )
    return degree, nonconstant_terms


def exact_single_shear_pass() -> None:
    """Check all single elementary shears of degree at most three.

    Degree lowering is checked over the algebraic closure by a common gcd.
    Support lowering is checked at every rational coefficient root returned
    by exact factorization.
    """

    x, y, z, lam = sp.symbols("x y z lambda")
    source_variables = (x, y, z)
    u = 1 + x * y
    foundational = (
        sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y)),
        sp.expand(y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)),
        sp.expand(2 * x - 3 * x**2 * y - x**3 * z),
    )
    assert expanded_profile(foundational, source_variables) == (7, 16)

    def common_gcd(expressions: list[sp.Expr]) -> sp.Poly | None:
        answer = None
        for expression in expressions:
            item = sp.Poly(expression, lam, domain=sp.QQ)
            if item.is_zero:
                continue
            answer = item if answer is None else sp.gcd(answer, item)
            if answer.degree() == 0:
                break
        return answer

    def audit(mapping: tuple[sp.Expr, ...]) -> None:
        polynomials = tuple(
            sp.Poly(sp.expand(component), *source_variables)
            for component in mapping
        )
        above_six = [
            coefficient
            for component in polynomials
            for monomial, coefficient in component.terms()
            if sum(monomial) > 6
        ]
        divisor = common_gcd(above_six)
        assert divisor is not None
        reduced = divisor
        while reduced.eval(0) == 0:
            reduced = sp.exquo(reduced, sp.Poly(lam, lam))
        assert reduced.degree() == 0

        rational_candidates = set()
        for component in polynomials:
            for monomial, coefficient in component.terms():
                if not any(monomial):
                    continue
                for root in sp.roots(coefficient, lam):
                    if root and root.is_Rational:
                        rational_candidates.add(root)
        for value in rational_candidates:
            specialized = tuple(
                sp.expand(component.subs(lam, value)) for component in mapping
            )
            assert expanded_profile(specialized, source_variables)[1] >= 16

    for index, variable in enumerate(source_variables):
        others = tuple(
            item for position, item in enumerate(source_variables)
            if position != index
        )
        for degree in range(1, 4):
            for exponent in range(degree + 1):
                monomial = others[0] ** exponent * others[1] ** (
                    degree - exponent
                )
                audit(
                    tuple(
                        sp.expand(
                            component.subs(
                                variable,
                                variable + lam * monomial,
                                simultaneous=True,
                            )
                        )
                        for component in foundational
                    )
                )

    targets = sp.symbols("T0:3")
    for index in range(3):
        others = tuple(position for position in range(3) if position != index)
        for degree in range(1, 4):
            for exponent in range(degree + 1):
                monomial = targets[others[0]] ** exponent * targets[
                    others[1]
                ] ** (degree - exponent)
                pulled_back = sp.expand(
                    monomial.subs(dict(zip(targets, foundational)))
                )
                transformed = list(foundational)
                transformed[index] = sp.expand(
                    transformed[index] + lam * pulled_back
                )
                audit(tuple(transformed))


Exponent = tuple[int, int]
Sparse = dict[Exponent, int]


def add(*polynomials: Sparse) -> Sparse:
    answer: Sparse = {}
    for item in polynomials:
        for monomial, coefficient in item.items():
            answer[monomial] = answer.get(monomial, 0) + coefficient
    return {monomial: coefficient for monomial, coefficient in answer.items()
            if coefficient}


def scale(polynomial_: Sparse, scalar: int) -> Sparse:
    return {
        monomial: scalar * coefficient
        for monomial, coefficient in polynomial_.items()
        if scalar * coefficient
    }


def multiply(left: Sparse, right: Sparse) -> Sparse:
    answer: Sparse = {}
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            monomial = (i + k, j + ell)
            answer[monomial] = answer.get(monomial, 0) + a * b
    return {monomial: coefficient for monomial, coefficient in answer.items()
            if coefficient}


def derivative(polynomial_: Sparse, variable: int) -> Sparse:
    answer: Sparse = {}
    for (i, j), coefficient in polynomial_.items():
        exponent = i if variable == 0 else j
        if exponent:
            monomial = (i - 1, j) if variable == 0 else (i, j - 1)
            answer[monomial] = exponent * coefficient
    return answer


def bracket(left: Sparse, right: Sparse) -> Sparse:
    return add(
        multiply(derivative(left, 0), derivative(right, 1)),
        scale(multiply(derivative(left, 1), derivative(right, 0)), -1),
    )


def determinant_linear_in_A(A: Sparse, B: Sparse, C: Sparse) -> Sparse:
    return add(
        scale(multiply(A, bracket(B, C)), -2),
        multiply(B, bracket(A, C)),
        multiply(C, bracket(A, B)),
    )


def solve_mod_prime(
    columns: list[Sparse], rhs: Sparse, prime: int
) -> tuple[list[int], int] | None:
    monomials = sorted(set(rhs).union(*(column.keys() for column in columns)))
    matrix = [
        [column.get(monomial, 0) % prime for column in columns]
        + [rhs.get(monomial, 0) % prime]
        for monomial in monomials
    ]
    row = 0
    pivots: list[int] = []
    for column in range(len(columns)):
        pivot = next(
            (
                index
                for index in range(row, len(matrix))
                if matrix[index][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        inverse = pow(matrix[row][column], prime - 2, prime)
        matrix[row] = [(entry * inverse) % prime for entry in matrix[row]]
        for index in range(len(matrix)):
            if index == row or not matrix[index][column]:
                continue
            factor = matrix[index][column]
            matrix[index] = [
                (entry - factor * pivot_entry) % prime
                for entry, pivot_entry in zip(matrix[index], matrix[row])
            ]
        pivots.append(column)
        row += 1
    for index in range(row, len(matrix)):
        if (
            all(matrix[index][column] == 0 for column in range(len(columns)))
            and matrix[index][-1]
        ):
            return None
    solution = [0] * len(columns)
    for index, column in enumerate(pivots):
        solution[column] = matrix[index][-1]
    return solution, len(pivots)


def random_collision_polynomial(
    constant_monomial: Exponent,
    extra_support: tuple[Exponent, ...],
    extra_count: int,
    rng: random.Random,
    height: int,
) -> Sparse | None:
    selected = rng.sample(extra_support, extra_count)
    coefficients: list[int] = []
    for _index in range(extra_count - 1):
        value = 0
        while value == 0:
            value = rng.randint(-height, height)
        coefficients.append(value)
    last = -1 - sum(coefficients)
    if last == 0 or abs(last) > 3 * height:
        return None
    coefficients.append(last)
    return {
        constant_monomial: 1,
        **dict(zip(selected, coefficients)),
    }


def modular_search(trials: int, seed: int) -> int:
    """Search the full box after exploiting linearity in A."""

    A_support = invariant_support(-2, 6)
    B_support = invariant_support(-1, 6)
    C_support = invariant_support(1, 6)
    assert tuple(map(len, (A_support, B_support, C_support))) == (8, 7, 5)
    A_extra = tuple(item for item in A_support if item != (0, 1))
    B_extra = tuple(item for item in B_support if item != (1, 0))
    C_extra = tuple(item for item in C_support if item != (0, 0))

    rng = random.Random(seed)
    prime = 65521
    survivors = 0
    for _trial in range(trials):
        B = random_collision_polynomial(
            (1, 0), B_extra, rng.randint(1, len(B_extra)), rng, 12
        )
        C = random_collision_polynomial(
            (0, 0), C_extra, rng.randint(1, len(C_extra)), rng, 12
        )
        if B is None or C is None:
            continue
        base = determinant_linear_in_A({(0, 1): 1}, B, C)
        rhs = add({(0, 0): -1}, scale(base, -1))
        columns = [
            determinant_linear_in_A({monomial: 1}, B, C)
            for monomial in A_extra
        ]
        if solve_mod_prime(columns, rhs, prime) is not None:
            survivors += 1
    return survivors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20_260_727)
    parser.add_argument(
        "--full-exact",
        action="store_true",
        help="also run the three slower characteristic-zero principal opens",
    )
    args = parser.parse_args()

    full_support = tuple(
        invariant_support(weight, 6) for weight in (-2, -1, 1)
    )
    assert tuple(map(len, full_support)) == (8, 7, 5)
    print("PASS: full degree-six equivariant support is 8+7+5 monomials")

    exact_z_linear_exclusion()
    print("PASS: exact QQ z-linear collision ideal has Groebner basis [1]")

    if args.full_exact:
        exact_full_open_exclusion()
        print("PASS: all three remaining QQ principal opens have basis [1]")

    exact_single_shear_pass()
    print("PASS: all single degree<=3 source/target shears miss both targets")

    survivors = modular_search(args.trials, args.seed)
    print(
        "SEARCH:",
        args.trials,
        "seeded full-box samples over F_65521;",
        survivors,
        "survivors",
    )
    assert survivors == 0


if __name__ == "__main__":
    main()
