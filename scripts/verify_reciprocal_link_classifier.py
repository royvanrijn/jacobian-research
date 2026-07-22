#!/usr/bin/env python3
"""Exact regressions for the reciprocal-link classifier and split obstruction."""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from jcsearch.reciprocal import (
    classify_boundary_reconstruction,
    classify_reciprocal_link,
    masuda_plinth_example,
    spectral_polynomial,
    spectral_obstruction,
    standard_cancellation_example,
)


def check_cancellation_certificate() -> None:
    candidate = standard_cancellation_example()
    certificate = classify_reciprocal_link(candidate)
    A = candidate.boundary
    assert dict(certificate.valuations) == {
        "s": -1,
        "Y": 0,
        "P": 1,
        "B": 0,
        "D": -1,
    }
    assert certificate.marked_valuation_pattern
    assert certificate.reciprocal_identity
    assert certificate.straightening_identity
    assert certificate.polynomial_straightening
    assert sp.cancel(certificate.chart_ratio + 1 / A) == 0
    assert sp.cancel(certificate.chart_jacobian + A) == 0
    assert sp.cancel(certificate.straightening_jacobian - A**2) == 0
    assert certificate.lnd_coefficients == (0, 0, 1)
    assert certificate.lnd_content == 1
    assert sp.cancel(certificate.lnd_on_B - A**2) == 0
    assert certificate.localized_coordinate_certificate
    assert certificate.boundary.boundary_parametrization_valid
    assert certificate.boundary.residue_degree == 1
    assert certificate.boundary.stein_degree == 1
    assert not certificate.spectral.excludes_candidate
    assert certificate.spectral.common_degree == 1
    assert certificate.verdict == "passes_marked_cancellation_boundary_prefilter"


def check_masuda_stein_degree() -> None:
    candidate = masuda_plinth_example()
    certificate = classify_boundary_reconstruction(candidate)
    tau = candidate.boundary_parameter
    Z = candidate.reconstruction_symbol
    assert certificate.boundary_parametrization_valid
    assert sp.factor(certificate.residue_polynomial - (Z**2 - tau)) == 0
    assert certificate.residue_degree == 2
    assert certificate.stein_degree == 2


def check_split_spectral_obstruction() -> None:
    tau = sp.symbols("tau", nonzero=True)
    Z = sp.symbols("Z")

    # A primitive degree-f cover Z^f=tau is incompatible with the Keller
    # leading equation for every bounded (m,r) below.  The manuscript proves
    # the all-degree statement because the spectral roots are constants.
    for cover_degree in range(2, 7):
        residue_polynomial = Z**cover_degree - tau
        for m in range(1, 5):
            for r in range(1, 5):
                certificate = spectral_obstruction(
                    residue_polynomial,
                    Z,
                    tau,
                    tau,
                    m,
                    r,
                )
                assert certificate.excludes_candidate
                assert certificate.common_degree == 0

    # An arithmetic gcd can package conjugate constant spectral roots into a
    # higher-degree factor.  It is still excluded geometrically because the
    # ground field of the theorem is algebraically closed.
    q = sp.symbols("q_spectral")
    conjugate_constraint = spectral_polynomial(2, 1, q).as_expr().subs(
        q, Z / tau**3
    )
    conjugate_numerator, _ = sp.fraction(sp.cancel(conjugate_constraint))
    conjugate = spectral_obstruction(
        conjugate_numerator,
        Z,
        tau,
        tau,
        2,
        1,
    )
    assert conjugate.common_degree == 2
    assert conjugate.excludes_candidate


def integrate_unit_interval(expression: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), variable)
    return sp.factor(
        sum(
            coefficient / (degree[0] + 1)
            for degree, coefficient in polynomial.terms()
        )
    )


def check_unsliced_boundary_spectral_equation() -> None:
    """Verify the pole coefficient used before any global slice exists."""
    A, y, b, u = sp.symbols("A y b u", nonzero=True)
    for m in range(1, 5):
        x = (A - 1) / y**m
        bracket = A - x * u * (y + x * b * (1 - u)) ** m
        boundary_bracket = sp.cancel(bracket.subs(A, 0))
        q = b / y ** (m + 1)
        expected_bracket = u * (1 - q * (1 - u)) ** m
        assert sp.cancel(boundary_bracket - expected_bracket) == 0
        for r in range(1, 5):
            pole_coefficient = integrate_unit_interval(boundary_bracket**r, u)
            spectral_value = integrate_unit_interval(
                u**r * (1 - q * (1 - u)) ** (m * r), u
            )
            assert sp.cancel(pole_coefficient - spectral_value) == 0


def main() -> None:
    check_cancellation_certificate()
    print("PASS reciprocal classifier: valuations, chart, LND, and degree-one residue")
    check_masuda_stein_degree()
    print("PASS boundary elimination: Masuda plinth model has exact Stein degree two")
    check_unsliced_boundary_spectral_equation()
    print("PASS unsliced Keller pole coefficient is J_(mr,r)(bar(B)/y^(m+1))")
    check_split_spectral_obstruction()
    print("PASS split-plinth spectral obstruction for f<=6 and m,r<=4")


if __name__ == "__main__":
    main()
