#!/usr/bin/env python3
"""Nonlinear Kuranishi screen for modular degree-seven/nine survivors.

The strong order-five consistency test is necessary but not sufficient: its
columns remember all coefficients of the quadratic lower-lift dependence,
not one common point of the lower-lift torsor.  This script reads modular
full-plane scan artifacts, projects the genuine quadratic order-five defect
modulo every current correction, and asks Singular whether the resulting
lower-lift scheme is empty. It also computes the radical; recognized degree
eight and nine fibres receive explicit nilpotent-structure certificates.

The output is still modular discovery.  A nonempty fibre must be reconstructed
over characteristic zero before any order-seven PBW work is authorized.
"""

from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import re
import shutil
import subprocess

import sympy as sp
from sympy.polys.domains import GF
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref

from verify_degree_seven_relative_quantization_obstruction import (
    family_presentation as degree_seven_presentation,
    pairing,
)
from verify_degree_eight_relative_quantization_obstruction import (
    family_presentation as degree_eight_presentation,
)
from verify_degree_nine_relative_quantization_obstruction import (
    family_presentation as degree_nine_presentation,
)


FAMILY_PRESENTATIONS = {
    7: degree_seven_presentation,
    8: degree_eight_presentation,
    9: degree_nine_presentation,
}


def projected_order_five_equations(presentation, field):
    """Return the complete quadratic Kuranishi equations modulo D5."""

    columns = presentation["correction_five"]
    constant = presentation["constant"]
    variations = presentation["lower_variations"]
    monomials = sorted(
        set(constant).union(
            *(set(column) for column in columns),
            *(set(column) for column in variations),
        )
    )
    output_index = {
        monomial: index for index, monomial in enumerate(monomials)
    }
    rows = {
        column_index: {
            output_index[monomial]: coefficient
            for monomial, coefficient in column.items()
        }
        for column_index, column in enumerate(columns)
        if column
    }
    reduced, pivots, nonzero = sdm_irref(rows)
    dual, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        len(monomials),
        pivots,
        nonzero,
    )

    variable_count = len(presentation["kernel_pairs"])
    linear = variations[0 : 2 * variable_count : 2]
    diagonal = variations[1 : 2 * variable_count : 2]
    cross = variations[2 * variable_count :]
    assert len(linear) == len(diagonal) == variable_count
    assert len(cross) == variable_count * (variable_count - 1) // 2

    equations = []
    for vector in dual:
        functional = {
            monomials[index]: coefficient
            for index, coefficient in vector.items()
            if coefficient
        }
        zero = (0,) * variable_count
        terms = {zero: pairing(functional, constant, field)}
        for index in range(variable_count):
            exponent = [0] * variable_count
            exponent[index] = 1
            terms[tuple(exponent)] = pairing(
                functional, linear[index], field
            )
            exponent[index] = 2
            terms[tuple(exponent)] = pairing(
                functional, diagonal[index], field
            )
        cross_index = 0
        for left in range(variable_count):
            for right in range(left + 1, variable_count):
                exponent = [0] * variable_count
                exponent[left] = exponent[right] = 1
                terms[tuple(exponent)] = pairing(
                    functional,
                    cross[cross_index],
                    field,
                )
                cross_index += 1
        equations.append(
            {
                exponent: coefficient
                for exponent, coefficient in terms.items()
                if coefficient
            }
        )
    return equations, len(pivots), len(dual)


def render_polynomial(equation, prime: int) -> str:
    terms = []
    for exponent, coefficient in sorted(equation.items()):
        factors = [str(int(coefficient) % prime)]
        for index, degree in enumerate(exponent):
            if degree == 1:
                factors.append(f"z{index}")
            elif degree > 1:
                factors.append(f"z{index}^{degree}")
        terms.append("*".join(factors))
    return "+".join(terms) if terms else "0"


def singular_fibre(equations, prime: int, variable_count: int):
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for the nonlinear screen")
    variables = ",".join(f"z{index}" for index in range(variable_count))
    ideal = ",".join(render_polynomial(equation, prime) for equation in equations)
    program = f"""LIB "primdec.lib";
ring r={prime},({variables}),dp;
option(redSB);
ideal I={ideal};
ideal G=std(I);
ideal R=radical(G);
ideal RG=std(R);
print("GB_SIZE="+string(size(G)));
print("DIMENSION="+string(dim(G)));
print("REDUCED_ONE="+string(reduce(1,G)));
for (int i=1; i<=size(G); i++) {{ print("G_"+string(i)+"="+string(G[i])); }}
print("RADICAL_SIZE="+string(size(RG)));
print("RADICAL_DIMENSION="+string(dim(RG)));
for (int i=1; i<=size(RG); i++) {{ print("R_"+string(i)+"="+string(RG[i])); }}
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if result.stderr.strip() or "   ?" in result.stdout:
        raise AssertionError(result.stdout + result.stderr)
    values = {}
    basis = []
    radical_basis = []
    for line in result.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key.startswith("G_"):
                basis.append(value)
            elif key.startswith("R_"):
                radical_basis.append(value)
            else:
                values[key] = value
    return (
        values,
        basis,
        radical_basis,
        hashlib.sha256(program.encode()).hexdigest(),
    )


def leading_monomial(text: str, variable_count: int):
    """Parse the first (monic) monomial printed by a reduced Singular basis."""

    match = re.match(r"(z\d+(?:\^\d+)?(?:\*z\d+(?:\^\d+)?)*)", text)
    if match is None:
        raise AssertionError(f"cannot parse leading monomial: {text}")
    exponent = [0] * variable_count
    for factor in match.group(1).split("*"):
        variable, _, degree = factor.partition("^")
        exponent[int(variable[1:])] = int(degree or 1)
    return tuple(exponent)


def degree_nine_structure(basis, radical_basis, prime: int, variable_count: int):
    """Recognize the finite rank-six thickening over affine six-space."""

    if variable_count != 12 or len(basis) != 9 or len(radical_basis) != 6:
        return None
    variables = sp.symbols(f"z0:{variable_count}")
    locals_ = {str(variable): variable for variable in variables}
    radical_polynomials = [
        sp.Poly(
            sp.sympify(text.replace("^", "**"), locals=locals_),
            *variables,
            modulus=prime,
        )
        for text in radical_basis
    ]
    radical_leads = [leading_monomial(text, variable_count) for text in radical_basis]
    dependent = {
        exponent.index(1)
        for exponent in radical_leads
        if sum(exponent) == 1
    }
    free = [index for index in range(variable_count) if index not in dependent]
    affine_radical = (
        len(dependent) == 6
        and all(polynomial.total_degree() <= 1 for polynomial in radical_polynomials)
        and all(
            all(
                exponent == 0 or index in free
                for monomial, _ in polynomial.terms()[1:]
                for index, exponent in enumerate(monomial)
            )
            for polynomial in radical_polynomials
        )
    )

    leads = [leading_monomial(text, variable_count) for text in basis]
    depends_only_on_radical_variables = all(
        all(exponent == 0 or index in dependent for index, exponent in enumerate(lead))
        for lead in leads
    )
    standard_count = 0
    if depends_only_on_radical_variables:
        dependent_order = sorted(dependent)
        for exponents in product(range(4), repeat=6):
            ambient = [0] * variable_count
            for index, exponent in zip(dependent_order, exponents, strict=True):
                ambient[index] = exponent
            if not any(
                all(left >= right for left, right in zip(ambient, lead, strict=True))
                for lead in leads
            ):
                standard_count += 1

    polynomials = [
        sp.Poly(
            sp.sympify(text.replace("^", "**"), locals=locals_),
            *variables,
            modulus=prime,
        )
        for text in basis
    ]
    univariate_cubics = []
    for polynomial in polynomials:
        active = {
            index
            for monomial, _ in polynomial.terms()
            for index, exponent in enumerate(monomial)
            if exponent
        }
        if polynomial.total_degree() == 3 and len(active) == 1:
            univariate_cubics.append((polynomial, active.pop()))
    triple_root = False
    triple_variable = None
    if len(univariate_cubics) == 1 and prime != 3:
        polynomial, triple_variable = univariate_cubics[0]
        variable = variables[triple_variable]
        expression = polynomial.monic().as_expr()
        root = (
            -int(expression.coeff(variable, 2)) * pow(3, -1, prime)
        ) % prime
        cube = sp.Poly((variable - root) ** 3, *variables, modulus=prime)
        triple_root = polynomial.monic() == cube
    certified = affine_radical and standard_count == 6 and triple_root
    return {
        "certified": certified,
        "radical": "affine six-space" if affine_radical else "unrecognized",
        "free_reduced_coordinates": [f"z{index}" for index in free],
        "finite_flat_rank_over_radical": standard_count,
        "triple_root_variable": (
            f"z{triple_variable}" if triple_variable is not None else None
        ),
        "univariate_cubic_is_perfect_cube": triple_root,
    }


def square_structure(basis, prime: int, variable_count: int):
    """Recognize the reduced affine-space pattern in degree-eight fibres."""

    if variable_count != 10 or len(basis) != 6:
        return None
    variables = sp.symbols(f"z0:{variable_count}")
    polynomials = [
        sp.Poly(
            sp.sympify(
                text.replace("^", "**"),
                locals={str(variable): variable for variable in variables},
            ),
            *variables,
            modulus=prime,
        )
        for text in basis
    ]
    linear = [polynomial for polynomial in polynomials if polynomial.total_degree() == 1]
    quadratic = [
        polynomial for polynomial in polynomials if polynomial.total_degree() == 2
    ]
    if len(linear) != 3 or len(quadratic) != 3:
        return None

    def active_indices(polynomial):
        return {
            index
            for monomial, _ in polynomial.terms()
            for index, exponent in enumerate(monomial)
            if exponent
        }

    q9 = next(
        (polynomial for polynomial in quadratic if active_indices(polynomial) == {9}),
        None,
    )
    q78 = next(
        (
            polynomial
            for polynomial in quadratic
            if active_indices(polynomial).issubset({7, 8})
            and polynomial.degree(variables[7]) == 2
        ),
        None,
    )
    if q9 is None or q78 is None:
        return None
    bridge = next(
        polynomial for polynomial in quadratic if polynomial not in {q9, q78}
    )

    z7, z8, z9 = variables[7:10]
    q9_expression = q9.as_expr()
    q9_leading = int(q9_expression.coeff(z9, 2)) % prime
    q9_linear = int(q9_expression.coeff(z9, 1)) % prime
    q9_constant = int(q9_expression.coeff(z9, 0)) % prime
    q9_discriminant = (q9_linear**2 - 4 * q9_leading * q9_constant) % prime
    root9 = (-q9_linear * pow(2 * q9_leading, -1, prime)) % prime

    q78_expression = q78.as_expr()
    q78_leading = int(q78_expression.coeff(z7, 2)) % prime
    q78_linear = sp.expand(q78_expression).coeff(z7, 1)
    q78_constant = sp.expand(q78_expression).coeff(z7, 0)
    q78_discriminant = sp.Poly(
        sp.expand(q78_linear**2 - 4 * q78_leading * q78_constant),
        z8,
        modulus=prime,
    )
    root7 = sp.Poly(
        -q78_linear * pow(2 * q78_leading, -1, prime),
        z8,
        modulus=prime,
    ).as_expr()
    bridge_on_reduction = sp.Poly(
        sp.expand(bridge.as_expr().subs({z9: root9, z7: root7})),
        z8,
        modulus=prime,
    )
    dependent_linear_variables = sorted(
        {
            str(variables[index])
            for polynomial in linear
            for index in (1, 3, 5)
            if polynomial.degree(variables[index]) == 1
        }
    )
    certified = (
        q9_discriminant == 0
        and q78_discriminant.is_zero
        and bridge_on_reduction.is_zero
        and dependent_linear_variables == ["z1", "z3", "z5"]
    )
    return {
        "certified": certified,
        "linear_dependent_variables": dependent_linear_variables,
        "square_root_variables": ["z7", "z9"],
        "free_reduced_coordinates": ["z0", "z2", "z4", "z6", "z8"],
        "z9_quadratic_discriminant_zero": q9_discriminant == 0,
        "z7_quadratic_discriminant_zero": q78_discriminant.is_zero,
        "bridge_vanishes_on_reduced_roots": bridge_on_reduction.is_zero,
    }


def screen_point(
    degree: int, prime: int, sigma_value: int, tau_value: int
):
    field = GF(prime)
    family_presentation = FAMILY_PRESENTATIONS[degree]
    presentation = family_presentation(
        field,
        field(sigma_value),
        field(tau_value),
    )
    equations, correction_rank, dual_rank = projected_order_five_equations(
        presentation, field
    )
    values, basis, radical_basis, digest = singular_fibre(
        equations,
        prime,
        len(presentation["kernel_pairs"]),
    )
    nonempty = values["REDUCED_ONE"] != "0"
    return {
        "prime": prime,
        "degree": degree,
        "sigma": sigma_value,
        "tau": tau_value,
        "correction_rank": correction_rank,
        "dual_rank": dual_rank,
        "projected_equations": len(equations),
        "lower_lift_variables": len(presentation["kernel_pairs"]),
        "nonempty": nonempty,
        "dimension": int(values["DIMENSION"]),
        "groebner_basis_size": int(values["GB_SIZE"]),
        "groebner_basis": basis,
        "radical_dimension": int(values["RADICAL_DIMENSION"]),
        "radical_groebner_basis": radical_basis,
        "reduced_square_structure": square_structure(
            basis, prime, len(presentation["kernel_pairs"])
        ),
        "degree_nine_thickening_structure": degree_nine_structure(
            basis,
            radical_basis,
            prime,
            len(presentation["kernel_pairs"]),
        ),
        "singular_program_sha256": digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("scans", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    records = []
    for scan_path in args.scans:
        scan = json.loads(scan_path.read_text())
        prime = scan["prime"]
        degree = scan.get("degree", 7)
        for point in scan["order_five_consistent_points"]:
            records.append(
                screen_point(
                    degree, prime, point["sigma"], point["tau"]
                )
            )
    certificate = {
        "scope": (
            "modular nonlinear Kuranishi screen; no characteristic-zero "
            "component and no order-seven claim"
        ),
        "points_screened": len(records),
        "nonempty_points": [record for record in records if record["nonempty"]],
        "empty_points": [record for record in records if not record["nonempty"]],
    }
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)

    print(
        "PASS: screened {} modular consistency points; nonlinear={}, empty={}".format(
            len(records),
            len(certificate["nonempty_points"]),
            len(certificate["empty_points"]),
        )
    )
    for record in certificate["nonempty_points"]:
        print(
            "NONLINEAR SURVIVOR: degree={} p={} sigma={} tau={} dimension={} gb={}".format(
                record["degree"],
                record["prime"],
                record["sigma"],
                record["tau"],
                record["dimension"],
                record["groebner_basis_size"],
            )
        )
    print("SCOPE: modular discovery only; reconstruct over Q before order seven")


if __name__ == "__main__":
    main()
