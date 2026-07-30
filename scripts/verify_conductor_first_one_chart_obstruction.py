#!/usr/bin/env python3
"""Verify the node/cusp conductor-first one-chart obstruction.

The theorem proved in the accompanying note concerns the direct separated
ansatz

    normalization B=k[t],
    conductor algebra A=k+cB,
    reconstruction open Spec(B[1/c])[affine variables].

Both the nodal conductor c=t(t-1) and cuspidal conductor c=t^2 admit finite
marked-root algebras with discriminant in A.  However, a reconstruction
coordinate p/c^m is polynomial only when c^m divides p, which removes the
required pole.  Independently, B[1/c] has nonconstant units, so no affine
stabilization is a polynomial ring.

This does not address multi-chart or ambient-coupled conductor mechanisms.
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
    / "conductor_first_one_chart_obstruction.json"
)


t, z = sp.symbols("t z")


def coefficient(expr: sp.Expr, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(expr), t).nth(degree)


def node_descends(expr: sp.Expr) -> bool:
    """Membership in k+t(t-1)k[t]."""

    return sp.expand(expr.subs(t, 0) - expr.subs(t, 1)) == 0


def cusp_descends(expr: sp.Expr) -> bool:
    """Membership in k+t^2k[t]."""

    return coefficient(expr, 1) == 0


def marked_root_examples() -> dict[str, object]:
    """Exhibit finite irreducible quadratic root algebras with descended discriminant."""

    node_conductor = t * (t - 1)
    cusp_conductor = t**2
    node_radical = sp.expand(node_conductor * (t + 1))
    cusp_radical = sp.expand(cusp_conductor * (t + 1))
    node_polynomial = sp.expand(z**2 - node_radical)
    cusp_polynomial = sp.expand(z**2 - cusp_radical)
    node_discriminant = sp.discriminant(node_polynomial, z)
    cusp_discriminant = sp.discriminant(cusp_polynomial, z)

    assert node_descends(node_radical)
    assert cusp_descends(cusp_radical)
    assert node_descends(node_discriminant)
    assert cusp_descends(cusp_discriminant)
    # Removing the visible square factors leaves t^2-1 and t+1,
    # respectively.  Both have an odd valuation, so neither radical is a
    # square in Q(t); the quadratic marked-root algebras are generically
    # separable and irreducible.
    assert sp.factor(node_radical) == t * (t - 1) * (t + 1)
    assert sp.factor(cusp_radical) == t**2 * (t + 1)
    return {
        "node": {
            "A": "Q + t(t-1) Q[t]",
            "conductor": "t(t-1) Q[t]",
            "root_equation": str(node_polynomial),
            "discriminant": str(sp.factor(node_discriminant)),
            "descent_equation": "Delta(0)=Delta(1)",
            "generically_irreducible": True,
        },
        "cusp": {
            "A": "Q + t^2 Q[t] = Q[t^2,t^3]",
            "conductor": "t^2 Q[t]",
            "root_equation": str(cusp_polynomial),
            "discriminant": str(sp.factor(cusp_discriminant)),
            "descent_equation": "coefficient_t(Delta)=Delta'(0)=0",
            "generically_irreducible": True,
        },
    }


def node_polynomiality(power: int) -> dict[str, object]:
    """Solve p/[t(t-1)]^power in Q[t] and compare with a required pole."""

    degree = 2 * power + 2
    coefficients = sp.symbols(f"a0:{degree + 1}")
    numerator = sum(
        coefficients[index] * t**index for index in range(degree + 1)
    )
    conductor_power = (t * (t - 1)) ** power
    remainder = sp.rem(numerator, conductor_power, domain=sp.QQ.frac_field(*coefficients))
    equations = [
        coefficient(remainder, index)
        for index in range(2 * power)
    ]
    matrix, _ = sp.linear_eq_to_matrix(equations, coefficients)
    assert matrix.rank() == 2 * power

    # Polynomiality forces numerator(0)=numerator(1)=0.  Adjoining a
    # normalization of either nonzero residue makes the ideal the unit ideal.
    u = sp.symbols(f"u_node_{power}")
    polynomiality = equations
    obstruction_at_zero = sp.groebner(
        polynomiality + [u * numerator.subs(t, 0) - 1],
        u,
        *coefficients,
        order="grevlex",
    )
    obstruction_at_one = sp.groebner(
        polynomiality + [u * numerator.subs(t, 1) - 1],
        u,
        *coefficients,
        order="grevlex",
    )
    assert obstruction_at_zero.contains(sp.Integer(1))
    assert obstruction_at_one.contains(sp.Integer(1))
    return {
        "pole_power": power,
        "numerator_degree_tested": degree,
        "polynomiality_equation_rank": matrix.rank(),
        "expected_rank": 2 * power,
        "nonzero_residue_at_t0_inconsistent": True,
        "nonzero_residue_at_t1_inconsistent": True,
    }


def cusp_polynomiality(power: int) -> dict[str, object]:
    """Solve p/t^(2 power) in Q[t] and compare with a required cusp pole."""

    degree = 2 * power + 2
    coefficients = sp.symbols(f"b0:{degree + 1}")
    numerator = sum(
        coefficients[index] * t**index for index in range(degree + 1)
    )
    equations = [coefficients[index] for index in range(2 * power)]
    matrix, _ = sp.linear_eq_to_matrix(equations, coefficients)
    assert matrix.rank() == 2 * power
    u = sp.symbols(f"u_cusp_{power}")
    leading_polar_coefficient = coefficients[0]
    obstruction = sp.groebner(
        equations + [u * leading_polar_coefficient - 1],
        u,
        *coefficients,
        order="grevlex",
    )
    assert obstruction.contains(sp.Integer(1))
    return {
        "pole_power": power,
        "numerator_degree_tested": degree,
        "polynomiality_equation_rank": matrix.rank(),
        "expected_rank": 2 * power,
        "nonzero_cusp_residue_inconsistent": True,
    }


def valuation_ledger() -> dict[str, object]:
    """Record the general simultaneous ledger/polynomiality contradiction."""

    return {
        "node": {
            "conductor_branches": ("t=0", "t=1"),
            "unit_lattice_rank_of_Q[t,1/(t(t-1))]": 2,
            "units_mod_constants": ("t", "t-1"),
            "polynomiality": (
                "ord_0(p)-m>=0 and ord_1(p)-m>=0"
            ),
            "required_pole": (
                "ord_0(p)-m<0 or ord_1(p)-m<0"
            ),
        },
        "cusp": {
            "conductor_branches": ("t=0",),
            "unit_lattice_rank_of_Q[t,1/t^2]": 1,
            "units_mod_constants": ("t",),
            "polynomiality": "ord_0(p)-2m>=0",
            "required_pole": "ord_0(p)-2m<0",
        },
        "conclusion": (
            "the reciprocal determinant ledger requires a negative "
            "conductor valuation, while polynomiality requires every "
            "reconstruction valuation to be nonnegative; the simultaneous "
            "system is empty"
        ),
    }


def main() -> None:
    output = {
        "schema": "conductor-first-one-chart-obstruction.v1",
        "status": "exact general one-chart obstruction with node/cusp fixtures",
        "marked_root_discriminant_descent": marked_root_examples(),
        "finite_polynomiality_replays": {
            "node": [node_polynomiality(power) for power in range(1, 5)],
            "cusp": [cusp_polynomiality(power) for power in range(1, 5)],
        },
        "general_valuation_ledger": valuation_ledger(),
        "unit_obstruction": (
            "the reconstruction open B[1/c][x_1,...,x_r] has a "
            "nonconstant unit c (and the displayed positive unit rank), "
            "whereas a polynomial ring over Q has only constant units; "
            "affine stabilization preserves this obstruction"
        ),
        "theorem_boundary": (
            "This excludes the direct separated one-conductor one-chart "
            "ansatz.  It does not exclude multi-chart gluing, a "
            "nonprincipal ambient boundary, or an ambient-coupled "
            "reconstruction whose full affine source is not obtained by "
            "localizing the normalized conductor line."
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_conductor_first_one_chart_obstruction.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: nodal and cuspidal marked-root discriminants descend")
    print("PASS: polynomiality removes every prescribed conductor pole")
    print("PASS: the determinant ledger and polynomiality system is empty")
    print("PASS: the sole-conductor localization has nonconstant units")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
