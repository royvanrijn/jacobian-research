"""Exact regressions for the one-additional-resolvent-factor theorem."""

from __future__ import annotations

from itertools import product

import sympy as sp


x, y, z, Avar = sp.symbols("x y z A")


def direct_formula(a: int, b: int, u: int, v: int, f2: sp.Expr) -> None:
    """Compare direct differentiation with the closed two-factor formula."""

    f1 = y**2 + y + 1
    gamma = y * Avar + Avar**2
    A = 1 + x * f1
    B = A**b * z + gamma.subs(Avar, A)
    P = A**a * B
    Q = y + x * A ** (a - 1) * B
    s = x / A

    jac_spq = sp.det(sp.Matrix([s, P, Q]).jacobian([x, y, z]))
    d1 = sp.cancel(1 - s * f1)
    d2 = sp.cancel(1 - s * f2)
    direct = sp.cancel(jac_spq * d1**u * d2**v)

    n = a + b - 2
    rho = sp.cancel(f2 / f1)
    formula = sp.cancel(-A ** (n - u - v) * (rho + (1 - rho) * A) ** v)
    assert sp.cancel(direct - formula) == 0


def is_nonzero_constant(expression: sp.Expr) -> bool:
    reduced = sp.cancel(expression)
    return reduced != 0 and not ({Avar, y} & reduced.free_symbols)


def bounded_classification() -> tuple[int, int]:
    """Exhaust a finite box without assuming the theorem's two outcomes."""

    f1 = y**2 + y + 1
    candidates = [
        sum(coefficient * y**j for j, coefficient in enumerate(coefficients))
        for coefficients in product((-1, 0, 1), repeat=3)
    ]
    tested = 0
    solutions = 0

    for a in range(1, 4):
        for b in range(1, 4):
            n = a + b - 2
            for u in range(0, 5):
                for v in range(1, 3):
                    for f2 in candidates:
                        rho = sp.cancel(f2 / f1)
                        expression = sp.cancel(
                            Avar ** (n - u - v)
                            * (rho + (1 - rho) * Avar) ** v
                        )
                        actual = is_nonzero_constant(expression)
                        expected = (f2 == f1 and n == u + v) or (
                            f2 == 0 and n == u
                        )
                        assert actual == expected
                        tested += 1
                        solutions += int(actual)
    return tested, solutions


def multifactor_formula(
    a: int, b: int, factors: list[tuple[sp.Expr, int]]
) -> None:
    """Check the arbitrary finite-product formula (11)."""

    f1 = y**2 + y + 1
    n = a + b - 2
    total_power = sum(power for _, power in factors)
    direct = Avar**n
    formula = Avar ** (n - total_power)
    for polynomial, power in factors:
        direct *= (1 - (Avar - 1) * polynomial / (Avar * f1)) ** power
        rho = sp.cancel(polynomial / f1)
        formula *= (rho + (1 - rho) * Avar) ** power
    assert sp.cancel(direct - formula) == 0


for multi_instance in [
    (1, 4, [(y**2 + y + 1, 1), (y**2 + y + 1, 2), (0, 3)]),
    (2, 3, [(0, 2), (y + 1, 1), (y**2 + y + 1, 1)]),
    (1, 2, [(0, 1), (0, 2), (y**2 + y + 1, 1)]),
]:
    multifactor_formula(*multi_instance)


for instance in [
    (1, 3, 1, 1, y**2 + y + 1),
    (2, 2, 0, 2, y**2 + y + 1),
    (1, 2, 1, 2, 0),
    (2, 3, 2, 1, y + 1),
]:
    direct_formula(*instance)

tested, solutions = bounded_classification()
assert solutions > 0
print(
    f"PASS: two-factor Jacobian formula and {solutions} rigid solutions "
    f"among {tested} bounded cases"
)
