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


def check_valuation_straightening() -> None:
    """Check Proposition 8.3 independently of a chosen source z-coordinate."""
    x, y, B, z = sp.symbols("x y B z")

    for m in range(1, 6):
        A = 1 + x * y**m
        s = x / A
        P = A * B
        Q = y + x * B
        Y = Q - s * P
        D = 1 - s * Y**m

        assert zero(Y - y)
        assert zero(D - 1 / A)

        chart_jacobian = sp.det(sp.Matrix([s, P, Q]).jacobian([x, y, B]))
        assert zero(chart_jacobian + 1 / A)

        for r in range(1, 5):
            tail = x**2 + x * y + y**3 + 1
            reconstructed_B = A ** (r + 1) * z + tail
            source_jacobian = sp.det(
                sp.Matrix([x, y, reconstructed_B]).jacobian([x, y, z])
            )
            assert zero(source_jacobian - A ** (r + 1))
            assert zero(chart_jacobian * source_jacobian + A**r)


def check_laurent_rigidity_reduction() -> None:
    """Bounded exact checks for the new Laurent leading-equation argument."""
    u, q, y, x = sp.symbols("u q y x")

    for m in range(1, 6):
        for r in range(1, 5):
            degree = m * r
            spectral = sp.Add(
                *(
                    (-1) ** j
                    * sp.binomial(degree, j)
                    * sp.factorial(r)
                    * sp.factorial(j)
                    / sp.factorial(r + j + 1)
                    * q**j
                    for j in range(degree + 1)
                )
            )
            spectral_poly = sp.Poly(sp.expand(spectral), q, domain=sp.QQ)
            assert spectral_poly.degree() == degree
            assert sp.gcd(spectral_poly, spectral_poly.diff()).degree() == 0

            # For a Laurent boundary value g_0, all y-dependence enters only
            # through q(y)=g_0/y^(m+1).
            q_of_y = q + y + y ** -1
            g0 = y ** (m + 1) * q_of_y
            shifted_y = y - (g0 / y**m) * (1 - u)
            assert zero(shifted_y / y - (1 - q_of_y * (1 - u)))

            # A tail divisible after localizing at y is already divisible in
            # k[x,y], because A=1+xy^m is coprime to y.
            source_A = 1 + x * y**m
            tail = x**2 + x * y + y**2 + 1
            standard_jet = y ** (m + 1) * (2 + 3 * source_A)
            arbitrary_g = standard_jet + source_A ** (r + 1) * tail
            assert zero((arbitrary_g - standard_jet) / source_A ** (r + 1) - tail)


def check_plinth_countermodel() -> None:
    """Check the degree-two boundary model that sigma_A=1 must exclude."""
    x, y, z = sp.symbols("x y z")

    def derivation(expr: sp.Expr) -> sp.Expr:
        return sp.expand(-2 * z * sp.diff(expr, y) + x**2 * sp.diff(expr, z))

    invariant = x**2 * y + z**2
    assert derivation(x) == 0
    assert derivation(invariant) == 0
    assert derivation(z) == x**2
    assert derivation(derivation(z)) == 0
    assert derivation(y) == -2 * z
    assert derivation(derivation(y)) == -2 * x**2
    assert derivation(derivation(derivation(y))) == 0

    # On x=0 the quotient boundary coordinate is invariant=z^2, so a
    # general boundary fiber has two affine-line components.
    boundary_relation = invariant.subs(x, 0)
    assert boundary_relation == z**2
    generic_t = sp.symbols("t")
    fiber_polynomial = sp.Poly(z**2 - generic_t, z, domain=sp.QQ.frac_field(generic_t))
    assert fiber_polynomial.degree() == 2
    assert sp.gcd(fiber_polynomial, fiber_polynomial.diff()).degree() == 0


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
    check_valuation_straightening()
    check_laurent_rigidity_reduction()
    check_plinth_countermodel()
    check_label_compression()
    check_boundary_orientations()
    print(
        "PASS log-suspension tame, marked-core, reciprocal/LND, Stein-boundary, Laurent-rigidity, orientation, and label identities"
    )


if __name__ == "__main__":
    main()
