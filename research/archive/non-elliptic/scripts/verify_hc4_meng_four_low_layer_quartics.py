#!/usr/bin/env python3
"""Select the four immutable low-layer quartics from the 234 HC4MQ1 parts."""

from __future__ import annotations

import contextlib
import io
from pathlib import Path
import runpy

import sympy as sp


PARENT = Path(__file__).with_name(
    "verify_hc4_meng_sparse_quartic_obstruction.py"
)
with contextlib.redirect_stdout(io.StringIO()):
    parent = runpy.run_path(str(PARENT))

quartic_exponents = parent["quartic_exponents"]
principal_unique_survivors = parent["principal_unique_survivors"]
principal_family_survivors = parent["principal_family_survivors"]
collision_monomials_exact = parent["collision_monomials_exact"]
scaled_collision_target = sp.Matrix(parent["scaled_collision_target"])
base_hessian = sp.Matrix(parent["base_hessian"])
tau = parent["tau"]

quartics = []
for support, _ in principal_unique_survivors:
    exponent_matrix = sp.Matrix.hstack(
        *(sp.Matrix(quartic_exponents[index]) for index in support)
    )
    solution = next(
        iter(sp.linsolve((exponent_matrix, scaled_collision_target)))
    )
    quartics.append(
        (
            support,
            tuple(
                sp.factor(
                    solution[position]
                    / collision_monomials_exact[index]
                )
                for position, index in enumerate(support)
            ),
        )
    )
for support, family, principal_gcd in principal_family_survivors:
    roots = sp.solve(principal_gcd.as_expr(), tau)
    assert len(roots) == 1
    quartics.append(
        (
            support,
            tuple(
                sp.factor(coefficient.subs(tau, roots[0]))
                for coefficient in family
            ),
        )
    )
assert len(quartics) == 234

spatial_variables = sp.symbols("x y r s")
base_adjugate = base_hessian.adjugate()
low_layer_quartics = []
for quartic_index, (support, coefficients) in enumerate(quartics):
    quartic = sp.expand(
        sum(
            coefficient
            * sp.prod(
                variable**exponent
                for variable, exponent in zip(
                    spatial_variables,
                    quartic_exponents[monomial],
                    strict=True,
                )
            )
            for monomial, coefficient in zip(
                support, coefficients, strict=True
            )
        )
    )
    degree_two_signature = sp.expand(
        sp.trace(
            base_adjugate * sp.hessian(quartic, spatial_variables)
        )
    )
    if degree_two_signature == 0:
        low_layer_quartics.append(
            (quartic_index, support, coefficients)
        )

expected = [
    (
        3,
        (0, 1, 5, 34),
        (
            sp.Rational(-5632, 1594323),
            sp.Rational(512, 177147),
            sp.Rational(-2048, 177147),
            sp.Rational(-81, 8),
        ),
    ),
    (
        45,
        (0, 32, 33, 34),
        (
            sp.Rational(-512, 531441),
            sp.Rational(3),
            sp.Rational(-12),
            sp.Rational(-297, 8),
        ),
    ),
    (
        87,
        (3, 4, 14, 18),
        (
            sp.Rational(-1, 54),
            sp.Rational(29, 576),
            sp.Rational(8, 9),
            sp.Rational(-3, 16),
        ),
    ),
    (
        111,
        (4, 12, 14, 24),
        (
            sp.Rational(1, 288),
            sp.Rational(32, 27),
            sp.Rational(116, 9),
            sp.Rational(12),
        ),
    ),
]
assert low_layer_quartics == expected

print(
    "PASS: exactly four of the 234 HC4MQ1 quartics have zero immutable "
    "determinant-degree-two signature"
)
print("DETAIL: quartic indices [3, 45, 87, 111]")
