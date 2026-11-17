#!/usr/bin/env python3
"""Verify the algebra in the smooth-quartic reciprocal frontend.

The exhaustive classification of the binary cubic root types is a written
argument in ``HC4_SMOOTH_QUARTIC_RECIPROCAL_FRONTEND.md``.  This checker
replays the matrix identities, the ten displayed boundary representatives,
their common-factor degrees, and the four basepoint-free quadratic-kernel
matrix families.
"""

from __future__ import annotations

import sympy as sp


def symmetric_matrix(prefix: str) -> sp.Matrix:
    entries = sp.symbols(f"{prefix}00 {prefix}01 {prefix}02 {prefix}11 {prefix}12 {prefix}22")
    a00, a01, a02, a11, a12, a22 = entries
    return sp.Matrix(
        [[a00, a01, a02], [a01, a11, a12], [a02, a12, a22]]
    )


def assert_zero_matrix(matrix: sp.Matrix) -> None:
    assert all(sp.expand(entry) == 0 for entry in matrix)


def verify_reciprocal_identities() -> None:
    q, ell, lam, a, mu = sp.symbols("q ell lam a mu")
    determinant_a, pairing = sp.symbols("detA pairing")

    # det(q*A+lambda*e*e^T), used on adj(C)=q*A+lambda*e*e^T.
    matrix = symmetric_matrix("m")
    e0, e1, e2 = sp.symbols("e0 e1 e2")
    vector = sp.Matrix([e0, e1, e2])
    rank_one_update = q * matrix + lam * vector * vector.T
    expected = q**3 * matrix.det() + lam * q**2 * (
        vector.dot(matrix.adjugate() * vector)
    )
    assert sp.expand(rank_one_update.det() - expected) == 0

    # If E=e^T adj(A)e, the two determinant/pairing equations eliminate E
    # to q*ell*(det(A)-q*ell*mu).  This is the scalar cancellation used in
    # the proof after mu=ell-lambda*a.
    first = q * determinant_a + lam * pairing - q**2 * ell**2
    second = mu * pairing - determinant_a * q * a
    elimination = sp.expand(mu * first - lam * second)
    assert sp.expand(
        elimination.subs(mu, ell - lam * a)
        - q * ell * (determinant_a - q * ell * (ell - lam * a))
    ) == 0

    # Formal matrix syzygies proving adj(A)=mu*C+lambda*d*d^T and
    # adj(A)e=q*ell*d once the displayed paired equations hold.
    matrix_a = symmetric_matrix("a")
    matrix_c = symmetric_matrix("c")
    d = sp.Matrix(sp.symbols("d0:3"))
    e = sp.Matrix(sp.symbols("f0:3"))
    identity = sp.eye(3)
    left = matrix_a * (mu * matrix_c + lam * d * d.T) - q * ell * mu * identity
    right = (
        mu * (matrix_a * matrix_c + lam * e * d.T - q * ell * identity)
        + lam * (matrix_a * d - mu * e) * d.T
    )
    assert_zero_matrix(left - right)

    adjugate_on_e = (mu * matrix_c + lam * d * d.T) * e - q * ell * d
    paired_reduction = (
        mu * (matrix_c * e - q * d)
        + lam * d * (d.dot(e) - q * a)
        + q * (mu + lam * a - ell) * d
    )
    assert_zero_matrix(adjugate_on_e - paired_reduction)

    # The first normal coefficient of det(A0+z*A1).
    z = sp.symbols("z")
    matrix_0 = symmetric_matrix("u")
    matrix_1 = symmetric_matrix("v")
    derivative = sp.diff((matrix_0 + z * matrix_1).det(), z).subs(z, 0)
    assert sp.expand(derivative - sp.trace(matrix_0.adjugate() * matrix_1)) == 0
    boundary_direction = sp.Matrix(sp.symbols("g0:3"))
    rank_one_adjugate = lam * boundary_direction * boundary_direction.T
    assert sp.expand(
        sp.trace(rank_one_adjugate * matrix_1)
        - lam * boundary_direction.dot(matrix_1 * boundary_direction)
    ) == 0


def common_factor_degree(entries: tuple[sp.Expr, ...], x: sp.Symbol, y: sp.Symbol) -> int:
    common = sp.Poly(entries[0], x, y)
    for entry in entries[1:]:
        common = sp.gcd(common, sp.Poly(entry, x, y))
    return common.total_degree()


def verify_boundary_representatives() -> None:
    x, y = sp.symbols("x y")
    one_third = sp.Rational(1, 3)
    charts = (
        ("squarefree-conic", one_third * (x**3 + y**3), x * y, 0),
        ("squarefree-line", one_third * (x**3 + y**3), 0, 0),
        ("double-conic", x**2 * y, y**2, 0),
        ("double-defect", x**2 * y, 0, 1),
        ("triple-line", one_third * x**3, y**2, 0),
        ("triple-simple-defect", one_third * x**3, x * y, 1),
        ("triple-double-defect", one_third * x**3, 0, 2),
        ("zero-squarefree-normal", 0, x * y, 2),
        ("zero-double-normal", 0, x**2, 2),
        ("zero-boundary", 0, 0, -1),
    )
    for name, binary_cubic, normal_quadratic, expected_degree in charts:
        direction = (
            sp.diff(binary_cubic, x),
            sp.diff(binary_cubic, y),
            sp.sympify(normal_quadratic),
        )
        if name == "zero-boundary":
            assert direction == (0, 0, 0)
            continue
        assert common_factor_degree(direction, x, y) == expected_degree


def verify_basepoint_free_boundary_matrices() -> None:
    x, y = sp.symbols("x y")
    p, q, r, t = sp.symbols("p q r t")
    h0, h1, h2 = sp.symbols("h0 h1 h2")
    binary_quadratic = h0 * x**2 + h1 * x * y + h2 * y**2

    cases: list[tuple[str, sp.Matrix, sp.Matrix, sp.Expr]] = []

    direction = sp.Matrix([x**2, y**2, x * y])
    matrix = sp.Matrix(
        [
            [p * y**2, -q * x * y, y * (q * y - p * x)],
            [-q * x * y, r * x**2, x * (q * x - r * y)],
            [
                y * (q * y - p * x),
                x * (q * x - r * y),
                p * x**2 - 2 * q * x * y + r * y**2,
            ],
        ]
    )
    cases.append(("squarefree-conic", direction, matrix, p * r - q**2))

    direction = sp.Matrix([x**2, y**2, 0])
    matrix = sp.Matrix(
        [
            [0, 0, -t * y**2],
            [0, 0, t * x**2],
            [-t * y**2, t * x**2, binary_quadratic],
        ]
    )
    cases.append(("squarefree-line", direction, matrix, -t**2))

    direction = sp.Matrix([2 * x * y, x**2, y**2])
    matrix = sp.Matrix(
        [
            [
                p * x**2 / 4 + q * x * y / 2 + r * y**2 / 4,
                -y * (p * x + q * y) / 2,
                -x * (q * x + r * y) / 2,
            ],
            [-y * (p * x + q * y) / 2, p * y**2, q * x * y],
            [-x * (q * x + r * y) / 2, q * x * y, r * x**2],
        ]
    )
    cases.append(("double-conic", direction, matrix, (p * r - q**2) / 4))

    direction = sp.Matrix([x**2, 0, y**2])
    matrix = sp.Matrix(
        [
            [0, -t * y**2, 0],
            [-t * y**2, binary_quadratic, t * x**2],
            [0, t * x**2, 0],
        ]
    )
    cases.append(("triple-line", direction, matrix, -t**2))

    for name, direction, matrix, scalar in cases:
        assert_zero_matrix(matrix * direction)
        assert_zero_matrix(matrix.adjugate() - scalar * direction * direction.T)
        assert scalar != 0
        print(f"PASS boundary matrix: {name}")


def main() -> None:
    verify_reciprocal_identities()
    print("PASS reciprocal determinant and adjugate identities")
    verify_boundary_representatives()
    print("PASS ten residual-line gradient representatives")
    verify_basepoint_free_boundary_matrices()
    print("THEOREM: smooth-quartic reciprocal frontend identities verified")


if __name__ == "__main__":
    main()
