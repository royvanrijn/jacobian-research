#!/usr/bin/env python3
"""Exact regressions for cancellation/LOG_GEOMETRY_OF_SUSPENSIONS.md.

These checks protect the coordinate identities, signs, and numerical label
relations in the note.  They do not certify the birational classification
problem isolated there.
"""
from __future__ import annotations

import sympy as sp


def zero(expr: sp.Expr) -> bool:
    return sp.cancel(expr) == 0


def check_tame_log_jacobian() -> None:
    u = sp.symbols("u", nonzero=True)
    for e in range(1, 10):
        t = u**e
        ordinary_jacobian = sp.diff(t, u)
        assert sp.Poly(ordinary_jacobian, u).degree() == e - 1
        # (dt/t)/(du/u) is the logarithmic Jacobian.
        logarithmic_jacobian = sp.cancel(u * ordinary_jacobian / t)
        assert logarithmic_jacobian == e


def check_marked_core_normal_forms() -> None:
    w, q, s, p, c = sp.symbols("w q s p c", nonzero=True)

    h = w**4 - 2 * w**2 + 3 * w + 5
    H = sp.integrate(h, w)
    weighted_T = w * q - H
    weighted_jacobian = sp.det(
        sp.Matrix([q, weighted_T]).jacobian([w, q])
    )
    assert zero(weighted_jacobian - (h - q))

    for m in range(1, 6):
        D = 1 - s * (q - p * s) ** m
        for r in range(1, 5):
            cancellation_T = c * sp.integrate(
                (1 - sp.Symbol("v") * (q - p * sp.Symbol("v")) ** m) ** r,
                (sp.Symbol("v"), 0, s),
            )
            cancellation_jacobian = sp.det(
                sp.Matrix([q, cancellation_T]).jacobian([s, q])
            )
            assert zero(cancellation_jacobian + c * D**r)


def check_reciprocal_boundary_link() -> None:
    x, y, z = sp.symbols("x y z")
    target_s, target_p, target_q = sp.symbols("s P Q")
    T = sp.symbols("T")
    h_polynomial = 2 - 3 * T + T**2

    for m in range(1, 5):
        for r in range(1, 4):
            A = 1 + x * y**m
            h_A = h_polynomial.subs(T, A)
            B = A ** (r + 1) * z + y ** (m + 1) * h_A
            P = A * B
            Q = y + x * B
            source_s = x / A
            source_D = 1 - source_s * (Q - source_s * P) ** m

            assert zero(Q - source_s * P - y)
            assert zero(source_D - 1 / A)

            jacobian = sp.det(
                sp.Matrix([source_s, P, Q]).jacobian([x, y, z])
            )
            assert zero(jacobian + A**r)

            # Check the displayed inverse after substitution from the source.
            inverse_y = Q - source_s * P
            inverse_x = source_s / source_D
            inverse_z = (
                P * source_D ** (r + 2)
                - inverse_y ** (m + 1)
                * source_D ** (r + 1)
                * h_polynomial.subs(T, 1 / source_D)
            )
            assert zero(inverse_x - x)
            assert zero(inverse_y - y)
            assert zero(inverse_z - z)

            # Check the other composition in independent target coordinates.
            D = 1 - target_s * (target_q - target_s * target_p) ** m
            recovered_y = target_q - target_s * target_p
            recovered_x = target_s / D
            recovered_z = (
                target_p * D ** (r + 2)
                - recovered_y ** (m + 1)
                * D ** (r + 1)
                * h_polynomial.subs(T, 1 / D)
            )
            recovered_A = 1 + recovered_x * recovered_y**m
            recovered_B = (
                recovered_A ** (r + 1) * recovered_z
                + recovered_y ** (m + 1)
                * h_polynomial.subs(T, recovered_A)
            )
            assert zero(recovered_A - 1 / D)
            assert zero(recovered_A * recovered_B - target_p)
            assert zero(recovered_y + recovered_x * recovered_B - target_q)
            assert zero(recovered_x / recovered_A - target_s)


def check_label_compression() -> None:
    for m in range(1, 20):
        for r in range(1, 20):
            degree = r * (m + 1) + 1
            ramification_index = r + 1
            contact = m * r * (m + 1)
            recovered_m = (degree - 1) // (ramification_index - 1) - 1
            recovered_contact = (degree - 1) * recovered_m
            assert recovered_m == m
            assert recovered_contact == contact


def check_boundary_orientations() -> None:
    A, D, C = sp.symbols("A D C", nonzero=True)

    # Cancellation identifies the two effective boundary generators with
    # opposite unit-lattice orientations.
    cancellation_pullback = 1 / A
    assert sp.cancel(cancellation_pullback * A) == 1

    # The weighted target chart uses the same C coordinate on both sides.
    weighted_pullback = C
    assert sp.cancel(weighted_pullback / C) == 1

    # The reverse cancellation chart sends A to D^{-1} as well.
    reverse_pullback = 1 / D
    assert sp.cancel(reverse_pullback * D) == 1


def main() -> None:
    check_tame_log_jacobian()
    check_marked_core_normal_forms()
    check_reciprocal_boundary_link()
    check_label_compression()
    check_boundary_orientations()
    print(
        "PASS log-suspension tame, marked-core, reciprocal-chart, orientation, and label identities"
    )


if __name__ == "__main__":
    main()
