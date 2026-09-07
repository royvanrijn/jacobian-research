#!/usr/bin/env python3
"""Verify the Fermat-symmetry orbit reduction after HC4NHM20.

The squarefree-line normal form is preserved by the automorphisms of the
binary Fermat cubic.  Rotation by a nontrivial cube root sends tau to
lambda*tau; the normalized reflection sends tau to 1/tau.  The fifteen
degenerate polar fibers therefore have only three orbit types.  The first
orbit is tau^3=-1, so the tau=-1 certificates of HC4NHM17 transport to the
two roots of tau^2-tau+1.
"""

from __future__ import annotations

import sympy as sp


def pivot(tau: sp.Expr, p: sp.Expr, q: sp.Expr, r: sp.Expr) -> sp.Expr:
    return (
        (3 * p**2 - q * r) * tau**5
        + (9 * p * r - q**2) * tau**4
        + (18 * r**2 - 6 * p * q) * tau**3
        + (18 * p**2 - 6 * q * r) * tau**2
        + (9 * p * r - q**2) * tau
        + (3 * r**2 - p * q)
    )


def reduce_cube_root(expression: sp.Expr, lam: sp.Symbol) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    relation = sp.Poly(lam**2 + lam + 1, lam)
    numerator = sp.rem(sp.Poly(sp.expand(numerator), lam), relation).as_expr()
    denominator = sp.rem(sp.Poly(sp.expand(denominator), lam), relation).as_expr()
    return sp.cancel(numerator / denominator)


def assert_zero_matrix(matrix: sp.MatrixBase) -> None:
    assert all(sp.cancel(entry) == 0 for entry in matrix)


def verify_packet_covariance() -> None:
    x, y, z = sp.symbols("x y z")
    tau, sigma, p, q, r = sp.symbols("tau sigma p q r", nonzero=True)
    u, v, w = sp.symbols("u v w")
    lam = sp.symbols("lambda")

    cubic = (x**3 + y**3) / 3 + z**2 * (u * x + v * y) + w * z**3
    matrix_0 = sp.Matrix(
        [
            [0, 0, -y**2],
            [0, 0, x**2],
            [-y**2, x**2, p * x**2 + q * x * y + r * y**2],
        ]
    )
    residual_line = y + tau * x + sigma * z

    # Rotation: (x,y,z) -> (lambda*x,y,z), lambda^3=1.
    rotation = sp.diag(lam, 1, 1)
    rotation_inverse = sp.diag(lam**2, 1, 1)
    rotated_matrix = (
        lam
        * rotation_inverse
        * matrix_0.subs({x: lam * x}, simultaneous=True)
        * rotation_inverse.T
    )
    rotated_expected = matrix_0.subs(
        {p: p, q: lam**2 * q, r: lam * r}, simultaneous=True
    )
    assert all(
        reduce_cube_root(entry, lam) == 0
        for entry in rotated_matrix - rotated_expected
    )
    rotated_cubic = cubic.subs({x: lam * x}, simultaneous=True)
    rotated_cubic_expected = cubic.subs({u: lam * u}, simultaneous=True)
    assert reduce_cube_root(rotated_cubic - rotated_cubic_expected, lam) == 0
    assert reduce_cube_root(
        residual_line.subs({x: lam * x}, simultaneous=True)
        - residual_line.subs({tau: lam * tau}, simultaneous=True),
        lam,
    ) == 0

    # Reflection: swap x,y and scale z by tau.  Dividing the transformed
    # residual line by the z-scale is exactly the packet normalization.
    reflection = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, tau]])
    reflection_inverse = sp.Matrix(
        [[0, 1, 0], [1, 0, 0], [0, 0, 1 / tau]]
    )
    reflected_matrix = sp.cancel(reflection.det()) * (
        reflection_inverse
        * matrix_0.subs({x: y, y: x, z: tau * z}, simultaneous=True)
        * reflection_inverse.T
    )
    reflected_expected = matrix_0.subs(
        {p: -r / tau, q: -q / tau, r: -p / tau}, simultaneous=True
    )
    assert_zero_matrix((reflected_matrix - reflected_expected).applyfunc(sp.cancel))
    reflected_cubic = cubic.subs(
        {x: y, y: x, z: tau * z}, simultaneous=True
    )
    reflected_cubic_expected = cubic.subs(
        {u: tau**2 * v, v: tau**2 * u, w: tau**3 * w},
        simultaneous=True,
    )
    assert sp.expand(reflected_cubic - reflected_cubic_expected) == 0
    reflected_line = residual_line.subs(
        {x: y, y: x, z: tau * z}, simultaneous=True
    )
    normalized_reflected_line = residual_line.subs(
        {tau: 1 / tau}, simultaneous=True
    )
    assert sp.expand(reflected_line - tau * normalized_reflected_line) == 0

    # The adjugate-plus-rank-one tensor transforms by congruence.  This is
    # the algebraic identity behind preservation of all Hessian curls.
    a00, a01, a02, a11, a12, a22 = sp.symbols(
        "a00 a01 a02 a11 a12 a22"
    )
    generic_matrix = sp.Matrix(
        [[a00, a01, a02], [a01, a11, a12], [a02, a12, a22]]
    )
    for transform, inverse in (
        (rotation, rotation_inverse),
        (reflection, reflection_inverse),
    ):
        transformed = transform.det() * inverse * generic_matrix * inverse.T
        adjugate_difference = transformed.adjugate() - (
            transform.T * generic_matrix.adjugate() * transform
        )
        if transform == rotation:
            assert all(
                reduce_cube_root(entry, lam) == 0 for entry in adjugate_difference
            )
        else:
            assert_zero_matrix(adjugate_difference.applyfunc(sp.cancel))

    print("PASS: rotations and normalized reflection preserve the full packet")


def verify_pivot_and_orbits() -> None:
    tau, p, q, r, lam = sp.symbols("tau p q r lambda", nonzero=True)

    rotated_pivot = pivot(lam * tau, p, lam**2 * q, lam * r)
    assert reduce_cube_root(rotated_pivot - lam**2 * pivot(tau, p, q, r), lam) == 0
    reflected_pivot = pivot(1 / tau, -r / tau, -q / tau, -p / tau)
    assert sp.cancel(reflected_pivot - pivot(tau, p, q, r) / tau**7) == 0
    print("PASS: Delta is equivariant under tau -> lambda*tau and tau -> 1/tau")

    degeneration = (
        (tau + 1)
        * (tau**2 - tau + 1)
        * (tau**4 - 4 * tau**3 + 10 * tau**2 - 4 * tau + 1)
        * (
            tau**8
            + 4 * tau**7
            + 6 * tau**6
            + 32 * tau**5
            + 83 * tau**4
            + 32 * tau**3
            + 6 * tau**2
            + 4 * tau
            + 1
        )
    )
    orbit_form = (tau**3 + 1) * (
        tau**12 + 44 * tau**9 + 586 * tau**6 + 44 * tau**3 + 1
    )
    assert sp.expand(degeneration - orbit_form) == 0
    assert reduce_cube_root(
        degeneration.subs(tau, lam * tau) - degeneration, lam
    ) == 0
    assert sp.cancel(degeneration.subs(tau, 1 / tau) - degeneration / tau**15) == 0

    s = sp.symbols("s", nonzero=True)
    j = sp.symbols("j")
    reciprocal_quartic = s**4 + 44 * s**3 + 586 * s**2 + 44 * s + 1
    assert sp.cancel(
        reciprocal_quartic / s**2 - ((s + 1 / s) ** 2 + 44 * (s + 1 / s) + 584)
    ) == 0
    quotient_polynomial = j**2 + 44 * j + 584
    assert sp.discriminant(quotient_polynomial, j) == -400
    assert sp.solve(quotient_polynomial, j) == [-22 - 10 * sp.I, -22 + 10 * sp.I]
    print("PASS: the 15 slopes form Fermat-symmetry orbits of sizes 3, 6, and 6")


def verify_transported_quadratic_orbit() -> None:
    p, q, r, lam = sp.symbols("p q r lambda")
    tau = -lam
    p_new, q_new, r_new = p, lam**2 * q, lam * r

    resultant_line = -3 * tau * p_new + (tau - 1) * q_new + 3 * r_new
    residual_polar_line = tau * p_new + r_new
    assert reduce_cube_root(resultant_line - lam * (3 * p + q + 3 * r), lam) == 0
    assert reduce_cube_root(residual_polar_line - lam * (r - p), lam) == 0

    residual_secondary = tau**2 * q_new**2 - 3 * tau * p_new * q_new + 8 * p_new**2
    resultant_secondary = 7 * p_new - 33 * tau**2 * r_new
    assert reduce_cube_root(
        residual_secondary - (q**2 + 3 * p * q + 8 * p**2), lam
    ) == 0
    assert reduce_cube_root(resultant_secondary - (7 * p - 33 * r), lam) == 0
    print("PASS: both tau=-1 line certificates and secondary strata transport")


def main() -> None:
    verify_packet_covariance()
    verify_pivot_and_orbits()
    verify_transported_quadratic_orbit()
    print("THEOREM: only three Fermat-symmetry normal forms occur among the line fibers")


if __name__ == "__main__":
    main()
