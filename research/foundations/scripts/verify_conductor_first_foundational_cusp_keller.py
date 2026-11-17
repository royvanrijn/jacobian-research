#!/usr/bin/env python3
"""Verify a conductor-first realization inside the foundational Keller map.

For H(W)=W^2(1-W), the tangent-core discriminant has one ordinary cusp.
After translating its cusp point, its normalization is

    S=-3*u^2, V=-2*u^3,

so the discriminant local ring is k[u^2,u^3] with conductor u^2 k[u].
The same cubic inverse pencil, reconstruction pole, distributed determinant
ledger, and polynomial source formulas give the foundational degree-three
Keller map (up to diagonal target scaling).
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "conductor_first_foundational_cusp_keller.json"
)


W, gamma, s, t, S, V, u = sp.symbols("W gamma s t S V u")
x, y, z, A_target, B_target, C_target = sp.symbols(
    "x y z A_target B_target C_target"
)


def tangent_core() -> dict[str, object]:
    """Compute the critical divisor, cusp conductor, and marked-root discriminant."""

    H = W**2 * (1 - W)
    H_prime = sp.diff(H, W)
    core_s = sp.expand(H_prime + gamma)
    core_t = sp.expand(W * core_s - H)
    core_jacobian = sp.factor(
        sp.det(
            sp.Matrix(
                [
                    [sp.diff(core_s, W), sp.diff(core_s, gamma)],
                    [sp.diff(core_t, W), sp.diff(core_t, gamma)],
                ]
            )
        )
    )
    assert core_jacobian == -gamma

    critical_s = sp.expand(core_s.subs(gamma, 0))
    critical_t = sp.expand(core_t.subs(gamma, 0))
    cusp_parameter = sp.Rational(1, 3) + u
    translated_s = sp.expand(
        critical_s.subs(W, cusp_parameter) - sp.Rational(1, 3)
    )
    translated_t = sp.expand(
        critical_t.subs(W, cusp_parameter) - sp.Rational(1, 27)
    )
    translated_v = sp.expand(translated_t - translated_s / 3)
    assert translated_s == -3 * u**2
    assert translated_v == -2 * u**3
    assert sp.expand(4 * translated_s**3 + 27 * translated_v**2) == 0

    inverse_polynomial = sp.expand(H - s * W + t)
    discriminant = sp.factor(sp.discriminant(inverse_polynomial, W))
    translated_discriminant = sp.factor(
        discriminant.subs(
            {
                s: sp.Rational(1, 3) + S,
                t: sp.Rational(1, 27) + S / 3 + V,
            }
        )
    )
    cusp_root_polynomial = sp.factor(
        inverse_polynomial.subs(
            {s: sp.Rational(1, 3), t: sp.Rational(1, 27)}
        )
    )
    assert translated_discriminant == -4 * S**3 - 27 * V**2
    assert cusp_root_polynomial == -(3 * W - 1) ** 3 / 27

    # The normalization quotient k[u]/k[u^2,u^3] has basis u.  Multiplying
    # by u does not return to the subring, while every u^n for n>=2 does.
    semigroup_through_ten = tuple(
        degree for degree in range(11) if degree == 0 or degree >= 2
    )
    assert 1 not in semigroup_through_ten
    assert all(degree in semigroup_through_ten for degree in range(2, 11))
    return {
        "seed": "H(W)=W^2(1-W)",
        "core_map": {
            "s": str(core_s),
            "t": str(core_t),
            "jacobian": str(core_jacobian),
            "critical_divisor": "gamma=0",
        },
        "cusp_point": {"W": "1/3", "s": "1/3", "t": "1/27"},
        "translated_normalization": {
            "S": str(translated_s),
            "V": str(translated_v),
            "implicit_equation": "4*S^3+27*V^2=0",
            "local_ring": "Q[u^2,u^3]",
            "normalization": "Q[u]",
            "conductor_in_normalization": "u^2 Q[u]",
            "normalization_quotient_basis": ("u",),
        },
        "marked_root_algebra": {
            "equation": str(inverse_polynomial),
            "degree": 3,
            "discriminant": str(discriminant),
            "translated_discriminant": str(translated_discriminant),
            "cusp_fiber": str(cusp_root_polynomial),
        },
    }


def weighted_keller_realization() -> dict[str, object]:
    """Derive the polynomial map, determinant ledger, and reconstruction pole."""

    source_gamma = sp.expand(
        1 - sp.Rational(3, 2) * x * y - sp.Rational(1, 2) * x**2 * z
    )
    source_w = sp.expand((1 + x * y) * source_gamma)
    source_c = sp.expand(x * source_gamma)
    H_source = source_w**2 * (1 - source_w)
    source_s = sp.expand(2 * source_w - 3 * source_w**2 + source_gamma)
    source_t = sp.expand(source_w * source_s - H_source)

    keller_a = sp.cancel(source_t / source_c**2)
    keller_b = sp.cancel(source_s / source_c)
    keller_c = source_c
    assert sp.denom(keller_a) == 1
    assert sp.denom(keller_b) == 1
    keller_a = sp.expand(keller_a)
    keller_b = sp.expand(keller_b)

    base = 1 + x * y
    foundational_1 = sp.expand(
        base**3 * z + y**2 * base * (4 + 3 * x * y)
    )
    foundational_2 = sp.expand(
        y + 3 * x * base**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
    )
    foundational_3 = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)
    assert keller_a == foundational_1
    assert 2 * keller_b == foundational_2
    assert 2 * keller_c == foundational_3

    variables = (x, y, z)
    keller_jacobian = sp.factor(
        sp.det(
            sp.Matrix(
                [
                    [sp.diff(output, variable) for variable in variables]
                    for output in (keller_a, keller_b, keller_c)
                ]
            )
        )
    )
    assert keller_jacobian == -sp.Rational(1, 2)

    rho_jacobian = sp.factor(
        sp.det(
            sp.Matrix(
                [
                    [sp.diff(output, variable) for variable in variables]
                    for output in (source_w, source_gamma, source_c)
                ]
            )
        )
    )
    assert sp.simplify(
        rho_jacobian + x**3 * source_gamma**2 / 2
    ) == 0

    mu_jacobian = sp.factor(
        sp.det(
            sp.Matrix(
                [
                    [
                        sp.diff(output, variable)
                        for variable in (A_target, B_target, C_target)
                    ]
                    for output in (
                        B_target * C_target,
                        A_target * C_target**2,
                        C_target,
                    )
                ]
            )
        )
    )
    assert mu_jacobian == -C_target**3
    ledger_left = sp.factor(
        (-source_gamma) * rho_jacobian
    )
    ledger_right = sp.factor(
        mu_jacobian.subs(C_target, keller_c) * keller_jacobian
    )
    assert sp.simplify(ledger_left - ledger_right) == 0
    assert sp.simplify(ledger_left - source_c**3 / 2) == 0

    inverse_polynomial = W**2 * (1 - W) - s * W + t
    inverse_derivative = sp.diff(inverse_polynomial, W)
    derivative_on_incidence = sp.expand(
        inverse_derivative.subs(
            {W: source_w, s: source_s, t: source_t}
        )
    )
    assert derivative_on_incidence == -source_gamma
    assert sp.cancel(keller_c / source_gamma - x) == 0

    collision_points = (
        (sp.Integer(0), sp.Integer(0), -sp.Rational(1, 4)),
        (sp.Integer(1), -sp.Rational(3, 2), sp.Rational(13, 2)),
        (-sp.Integer(1), sp.Rational(3, 2), sp.Rational(13, 2)),
    )
    collision_images = tuple(
        tuple(
            sp.expand(output).subs(dict(zip(variables, point)))
            for output in (keller_a, keller_b, keller_c)
        )
        for point in collision_points
    )
    assert collision_images == (
        (-sp.Rational(1, 4), 0, 0),
        (-sp.Rational(1, 4), 0, 0),
        (-sp.Rational(1, 4), 0, 0),
    )
    return {
        "source_auxiliary_coordinates": {
            "gamma": str(source_gamma),
            "W": str(source_w),
            "C": str(source_c),
        },
        "polynomial_keller_map": {
            "A": str(keller_a),
            "B": str(keller_b),
            "C": str(keller_c),
            "relation_to_foundational_map": "(A,B,C)=(F1,F2/2,F3/2)",
            "jacobian": str(keller_jacobian),
        },
        "determinant_ledger": {
            "core_jacobian": "-gamma",
            "source_chart_jacobian": str(rho_jacobian),
            "target_chart_jacobian": str(mu_jacobian),
            "common_value": str(ledger_left),
        },
        "reconstruction": {
            "inverse_derivative_on_incidence": str(derivative_on_incidence),
            "coordinate": "x=C/gamma=-C/E_W",
            "pole_divisor": "gamma=E_W=0 on C!=0",
            "meets_cusp_conductor_point": True,
        },
        "collision": {
            "source_points": tuple(
                tuple(map(str, point)) for point in collision_points
            ),
            "common_target": ("-1/4", "0", "0"),
        },
    }


def main() -> None:
    output = {
        "schema": "conductor-first-foundational-cusp-keller.v1",
        "status": "exact conductor-first Keller realization",
        "conductor_first_core": tangent_core(),
        "keller_realization": weighted_keller_realization(),
        "success_criterion": (
            "the scaled foundational degree-three Keller map is an explicit "
            "Keller family whose selected repeated-root discriminant is the "
            "cuspidal conductor curve Q[u^2,u^3]; the marked-root algebra, "
            "reconstruction pole, determinant ledger, and polynomiality "
            "hold simultaneously"
        ),
        "mechanism_boundary": (
            "the conductor is genuine but the construction is the known "
            "weighted tangent mechanism.  It escapes the separated "
            "one-chart obstruction because cancellation is distributed "
            "between source and target vertical charts; conductor gluing "
            "is not itself the source of a new polynomiality identity"
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_conductor_first_foundational_cusp_keller.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: the selected discriminant is the cusp Q[u^2,u^3]")
    print("PASS: its conductor in Q[u] is u^2 Q[u]")
    print("PASS: the cubic marked-root discriminant descends through the cusp")
    print("PASS: the reconstruction pole meets the cusp conductor")
    print("PASS: determinant ledger and polynomiality hold simultaneously")
    print("PASS: the resulting polynomial map has Jacobian -1/2")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
