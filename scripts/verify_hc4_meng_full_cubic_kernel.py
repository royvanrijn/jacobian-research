#!/usr/bin/env python3
"""Exclude arbitrary cubic corrections in the sparse Meng quartic chart.

The sparse-quartic certificate leaves 234 rational quartic principal parts.
For each one, this checker constructs the complete 20-column odd determinant
signature over QQ and parameterizes its full kernel, of dimension 12, 13, or
16.  It then computes the determinant polynomial symbolically, groups its
spatial coefficients by degree, and sends the descending degree-six through
degree-two coefficient ideals to Singular over QQ.  Every graded chain
reaches the unit ideal.

Four exceptional quartics have the form u^3*L and a 16-dimensional odd
kernel equal to the H0-harmonic cubic space.  They use the triangular chart

    h3 = a + u*b + u^2*c + u^3*d,

obtained by integrating the harmonic equation along the variable paired with
u by the wave operator d_y d_r + (1/2)d_x d_s.

Requires Singular on PATH.  This is a characteristic-zero certificate.
"""

from __future__ import annotations

import contextlib
import io
from pathlib import Path
import runpy
import shutil
import subprocess

import sympy as sp


assert shutil.which("Singular") is not None, "Singular is required on PATH"

PARENT = Path(__file__).with_name(
    "verify_hc4_meng_sparse_quartic_obstruction.py"
)
with contextlib.redirect_stdout(io.StringIO()):
    parent = runpy.run_path(str(PARENT))

quartic_exponents = parent["quartic_exponents"]
cubic_exponents = parent["cubic_exponents"]
base_hessian = parent["base_hessian"]
collision_monomials_exact = parent["collision_monomials_exact"]
principal_unique_survivors = parent["principal_unique_survivors"]
principal_family_survivors = parent["principal_family_survivors"]
scaled_collision_target = parent["scaled_collision_target"]
monomial_hessian_exact = parent["monomial_hessian_exact"]
tau = parent["tau"]

x, y, r, s = spatial_variables = sp.symbols("x y r s")
base_potential = 2 * y * r + 4 * x * s
degree_points = tuple(
    tuple(sp.Rational(coordinate) for coordinate in point)
    for point in parent["two_cubic_points"][4:9]
)


def reconstruct_quartics():
    quartics = []
    target = sp.Matrix(scaled_collision_target)
    for support, _ in principal_unique_survivors:
        exponent_matrix = sp.Matrix.hstack(
            *(sp.Matrix(quartic_exponents[index]) for index in support)
        )
        solution = next(iter(sp.linsolve((exponent_matrix, target))))
        coefficients = tuple(
            sp.factor(
                solution[position] / collision_monomials_exact[index]
            )
            for position, index in enumerate(support)
        )
        quartics.append((support, coefficients))
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
    return tuple(quartics)


def polynomial_from_support(support, coefficients, exponents):
    return sp.expand(
        sum(
            coefficient
            * sp.prod(
                spatial_variables[index] ** exponents[monomial][index]
                for index in range(4)
            )
            for coefficient, monomial in zip(
                coefficients, support, strict=True
            )
        )
    )


def exact_odd_kernel(support, coefficients):
    quartic_hessians = []
    for point in degree_points:
        hessian = sp.zeros(4)
        for coefficient, monomial in zip(
            coefficients, support, strict=True
        ):
            hessian += (
                coefficient
                * monomial_hessian_exact(
                    quartic_exponents[monomial], point
                )
            )
        quartic_hessians.append(hessian)

    base_adjugate = base_hessian.adjugate()
    signatures = []
    for exponents in cubic_exponents:
        signature = [
            sp.factor(
                sp.trace(
                    quartic_hessian.adjugate()
                    * monomial_hessian_exact(exponents, point)
                )
            )
            for quartic_hessian, point in zip(
                quartic_hessians, degree_points, strict=True
            )
        ]
        signature.extend(
            sp.factor(
                sp.trace(
                    base_adjugate
                    * monomial_hessian_exact(exponents, point)
                )
            )
            for point in degree_points
        )
        signatures.append(signature)
    matrix = sp.Matrix.hstack(
        *(sp.Matrix(signature) for signature in signatures)
    )
    return matrix.rank(), tuple(matrix.nullspace())


def cubic_from_kernel(kernel):
    parameters = sp.symbols(f"z0:{len(kernel)}")
    cubic = sp.Rational(0)
    monomials = tuple(
        sp.prod(
            spatial_variables[index] ** exponents[index]
            for index in range(4)
        )
        for exponents in cubic_exponents
    )
    for parameter, vector in zip(parameters, kernel, strict=True):
        cubic += parameter * sum(
            vector[index] * monomials[index] for index in range(20)
        )
    return parameters, sp.expand(cubic)


def harmonic_cubic_chart(cube_variable):
    pairs = {
        s: (x, y, r, sp.Rational(1, 2), sp.Rational(1)),
        x: (s, y, r, sp.Rational(1, 2), sp.Rational(1)),
        r: (y, x, s, sp.Rational(1), sp.Rational(1, 2)),
        y: (r, x, s, sp.Rational(1), sp.Rational(1, 2)),
    }
    paired, other_left, other_right, paired_weight, other_weight = (
        pairs[cube_variable]
    )
    a_parameters = sp.symbols("a0:10")
    b_parameters = sp.symbols("b0:3")
    c_parameters = sp.symbols("c0:2")
    d_parameter = sp.symbols("d")
    parameters = (
        a_parameters + b_parameters + c_parameters + (d_parameter,)
    )

    a_monomials = [
        paired**left
        * other_left**middle
        * other_right ** (3 - left - middle)
        for left in range(4)
        for middle in range(4 - left)
    ]
    a = sum(
        parameter * monomial
        for parameter, monomial in zip(
            a_parameters, a_monomials, strict=True
        )
    )
    b_zero = (
        b_parameters[0] * other_left**2
        + b_parameters[1] * other_left * other_right
        + b_parameters[2] * other_right**2
    )
    b_derivative = (
        -other_weight
        / paired_weight
        * sp.diff(a, other_left, other_right)
    )
    b = sp.integrate(b_derivative, paired) + b_zero
    c_zero = (
        c_parameters[0] * other_left
        + c_parameters[1] * other_right
    )
    c_derivative = (
        -other_weight
        / (2 * paired_weight)
        * sp.diff(b, other_left, other_right)
    )
    c = sp.integrate(c_derivative, paired) + c_zero
    cubic = sp.expand(
        a
        + cube_variable * b
        + cube_variable**2 * c
        + cube_variable**3 * d_parameter
    )
    assert sp.expand(
        sp.diff(cubic, y, r)
        + sp.Rational(1, 2) * sp.diff(cubic, x, s)
    ) == 0
    return parameters, cubic


def exceptional_cube_variable(quartic):
    for variable in spatial_variables:
        quotient = sp.cancel(quartic / variable**3)
        if sp.Poly(quotient, *spatial_variables).total_degree() == 1:
            return variable
    return None


def spatial_layers(quartic, cubic, parameters):
    determinant = sp.expand(
        sp.hessian(
            base_potential + cubic + quartic,
            spatial_variables,
        ).det(method="berkowitz")
        - 64
    )
    polynomial = sp.Poly(
        determinant,
        *spatial_variables,
        *parameters,
        domain=sp.QQ,
    )
    grouped = {degree: {} for degree in range(2, 7)}
    outside = {}
    for exponents, coefficient in polynomial.terms():
        spatial_exponents = exponents[:4]
        parameter_exponents = exponents[4:]
        spatial_degree = sum(spatial_exponents)
        term = coefficient
        for parameter, exponent in zip(
            parameters, parameter_exponents, strict=True
        ):
            term *= parameter**exponent
        destination = (
            grouped[spatial_degree]
            if spatial_degree in grouped
            else outside
        )
        destination[spatial_exponents] = (
            destination.get(spatial_exponents, 0) + term
        )
    assert all(sp.expand(value) == 0 for value in outside.values())

    layers = {}
    for degree, coefficients in grouped.items():
        layers[degree] = []
        for coefficient in coefficients.values():
            if sp.expand(coefficient) == 0:
                continue
            _, cleared = sp.Poly(
                coefficient, *parameters, domain=sp.QQ
            ).clear_denoms(convert=True)
            layers[degree].append(cleared.as_expr())
    return layers


def singular_ideal(name, polynomials):
    body = ",\n".join(
        str(polynomial).replace("**", "^")
        for polynomial in polynomials
    )
    return f"ideal {name}={body if body else '0'};\n"


def graded_unit_layer(parameters, layers):
    variables = ",".join(str(parameter) for parameter in parameters)
    script = (
        f"ring rr=0,({variables}),dp; option(redSB);\n"
        + singular_ideal("I6", layers[6])
        + 'ideal G=slimgb(I6); print("L6"); print(size(G));'
        + 'if(size(G)==1&&G[1]==1){print("UNIT6");quit;};\n'
    )
    for degree in (5, 4, 3, 2):
        script += (
            singular_ideal(f"I{degree}", layers[degree])
            + f"ideal J{degree}=G,I{degree};"
            + f"G=slimgb(J{degree});"
            + f'print("L{degree}");print(size(G));'
            + f'if(size(G)==1&&G[1]==1)'
            + f'{{print("UNIT{degree}");quit;}};\n'
        )
    script += 'print("NONUNIT");quit;\n'
    result = subprocess.run(
        ["Singular", "-q"],
        input=script,
        text=True,
        capture_output=True,
        timeout=900,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "?" not in result.stdout, result.stdout
    for degree in (6, 5, 4, 3, 2):
        if f"UNIT{degree}" in result.stdout:
            return degree
    raise AssertionError(result.stdout)


quartics = reconstruct_quartics()
rank_counts = {4: 0, 7: 0, 8: 0}
unit_layer_counts = {2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
exceptional_indices = []

for quartic_index, (support, coefficients) in enumerate(quartics):
    quartic = polynomial_from_support(
        support, coefficients, quartic_exponents
    )
    rank, kernel = exact_odd_kernel(support, coefficients)
    assert rank in rank_counts
    rank_counts[rank] += 1

    cube_variable = exceptional_cube_variable(quartic)
    if cube_variable is None:
        parameters, cubic = cubic_from_kernel(kernel)
    else:
        assert rank == 4 and len(kernel) == 16
        exceptional_indices.append(quartic_index)
        parameters, cubic = harmonic_cubic_chart(cube_variable)

    layers = spatial_layers(quartic, cubic, parameters)
    unit_layer = graded_unit_layer(parameters, layers)
    unit_layer_counts[unit_layer] += 1
    print(
        f"PASS quartic {quartic_index:03d}: rank={rank}, "
        f"kernel={len(parameters)}, unit-at-degree={unit_layer}",
        flush=True,
    )

assert rank_counts == {4: 4, 7: 1, 8: 229}
assert exceptional_indices == [0, 81, 158, 231]
assert sum(unit_layer_counts.values()) == 234

print(
    "PASS: full odd-kernel ranks over QQ are rank 8 for 229 quartics, "
    "rank 7 for one quartic, and rank 4 for four quartics"
)
print(
    "PASS: every full cubic-kernel determinant ideal is the unit ideal "
    "over QQ"
)
print(f"DETAIL: first unit determinant layers {unit_layer_counts}")
print(
    "SCOPE: arbitrary cubic corrections are excluded for every collision "
    "quartic supported on at most four monomials"
)
