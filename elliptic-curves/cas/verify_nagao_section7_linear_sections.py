#!/usr/bin/env python3
"""Classify Nagao section-7 linear sections and prove generic rank at least 12.

The classification is deliberately limited to polynomial sections

``x=m*T+n,  y in Q[T],  degree(y)<=3``.

It is not a computation of the full Mordell--Weil group.  Exact symbolic
addition proves that five of the six extra linear sections are dependent on
eleven visible Mestre sections and the ``+7/27`` section.  Exact reduction
modulo several good primes at ``T=1`` proves those twelve basis sections
independent, and therefore proves generic rank at least 12 over ``Q(T)``.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
from typing import Any

import sympy as sp

from mod_l_reduction_independence import (
    combined_mod_l_rank,
    find_mod_l_reduction_certificate,
    find_no_rational_l_torsion_prime,
)
from nagao_1994 import (
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    short_jacobian_coefficients,
)
from nagao_1994_section7 import (
    PRIMARY_SOURCE,
    SECTION7_CONSTRUCTION,
    SECTION7_JACOBIAN_RELATIONS,
    SECTION7_LINEAR_COMPANION_SECTIONS,
    SECTION7_QUADRATIC_COMPANION_SECTIONS,
    SECTION7_QUADRATIC_JACOBIAN_RELATIONS,
    SECTION7_RELATION_BASIS_VISIBLE_INDICES,
    SECTION7_ROOTS,
    section7_primitive_quartic_coefficients,
)
from verify_nagao_linear_sections import _covariant_map, _elliptic_add


Q = Fraction
SPECIALIZATION_PARAMETER = Q(1)
INDEPENDENCE_MODULUS = 3
EXPECTED_CERTIFICATE_PRIMES = (
    5,
    17,
    37,
    41,
    47,
    67,
    83,
    89,
    97,
    101,
    107,
    127,
)
EXPECTED_TORSION_CERTIFICATE_PRIME = 7
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/verify_nagao_section7_linear_sections.py"
)


def _sympy_rational(value: Fraction) -> sp.Rational:
    value = Q(value)
    return sp.Rational(value.numerator, value.denominator)


def _quartic(parameter: sp.Symbol, x: sp.Symbol) -> sp.Expr:
    coefficients = (
        9 * parameter**6
        - 910748 * parameter**4
        + 23718659440 * parameter**2
        + 557726319412900,
        18
        * (
            354 * parameter**4
            - 17901331 * parameter**2
            - 884640359570
        ),
        -3
        * (
            6 * parameter**4
            - 668642 * parameter**2
            - 52052853547
        ),
        -54 * (118 * parameter**2 + 11538729),
        9 * (parameter**2 + 96714),
    )
    return sp.expand(
        sum(coefficient * x**degree for degree, coefficient in enumerate(coefficients))
    )


def _monic(polynomial: sp.Poly) -> sp.Poly:
    return sp.Poly(polynomial.as_expr() / polynomial.LC(), *polynomial.gens)


def verify_k3_fiber_geometry() -> dict[str, Any]:
    """Verify the singular fibers and the resulting geometric rank bound."""

    parameter, x = sp.symbols("T X")
    polynomial = sp.Poly(_quartic(parameter, x), x)
    e, d, c, b, a = (
        polynomial.coeff_monomial(x**degree) for degree in range(5)
    )
    invariant_i = sp.expand(12 * a * e - 3 * b * d + c**2)
    invariant_j = sp.expand(
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    coefficient_a = sp.expand(-27 * invariant_i)
    coefficient_b = sp.expand(-27 * invariant_j)
    discriminant = sp.Poly(
        sp.expand(-16 * (4 * coefficient_a**3 + 27 * coefficient_b**2)),
        parameter,
    )
    if (
        sp.degree(coefficient_a, parameter) != 8
        or sp.degree(coefficient_b, parameter) != 12
        or discriminant.degree() != 20
    ):
        raise AssertionError("the section-7 K3 degrees changed")
    if sp.gcd(discriminant, discriminant.diff()).degree() != 0:
        raise AssertionError("the finite discriminant polynomial is not squarefree")
    if sp.gcd(discriminant, sp.Poly(invariant_i, parameter)).degree() != 0:
        raise AssertionError("c4 vanishes at a finite discriminant root")

    local = sp.symbols("s")
    infinity_a = sp.Poly(
        sp.expand(local**8 * coefficient_a.subs(parameter, 1 / local)), local
    )
    infinity_b = sp.Poly(
        sp.expand(local**12 * coefficient_b.subs(parameter, 1 / local)), local
    )
    infinity_delta = sp.Poly(
        sp.expand(local**24 * discriminant.as_expr().subs(parameter, 1 / local)),
        local,
    )
    infinity_delta_valuation = min(
        monomial[0]
        for monomial, coefficient in infinity_delta.terms()
        if coefficient
    )
    if infinity_delta_valuation != 4 or infinity_a.eval(0) == 0:
        raise AssertionError("the infinity fiber is not multiplicative I4")
    special_cubic = sp.Poly(
        x**3 + infinity_a.eval(0) * x + infinity_b.eval(0), x
    )
    expected_cubic = sp.Poly((x - 108) ** 2 * (x + 216), x)
    if special_cubic != expected_cubic:
        raise AssertionError("the infinity nodal cubic changed")

    # After X=108+u, the tangent cone is Y^2=324*u^2, so the node is split
    # over Q.  The fiber configuration is 20 I1 plus one split I4.  A minimal
    # Weierstrass model with coefficient degrees 8 and 12 has chi=2, hence is
    # a K3 surface.  Shioda--Tate and rho<=h^{1,1}=20 give
    # rank <= 20-2-(4-1)=15 over Qbar(T).
    return {
        "short_coefficient_degrees": [8, 12],
        "discriminant_degree_on_affine_line": 20,
        "finite_discriminant_squarefree": True,
        "finite_singular_fibers": {"type": "I1", "count": 20},
        "infinity_fiber": {
            "type": "I4",
            "split": True,
            "minimal_discriminant_valuation": infinity_delta_valuation,
            "special_cubic": "(X-108)^2*(X+216)",
            "tangent_square": "324=18^2",
        },
        "total_euler_number": 24,
        "surface_type": "elliptic K3",
        "shioda_tate_component_rank": 3,
        "complex_picard_upper_bound": 20,
        "geometric_generic_rank_upper_bound": 15,
    }


def classify_abscissae() -> dict[str, Any]:
    """Eliminate the ordinate in the linear-abscissa ansatz exactly."""

    parameter, x, slope, intercept = sp.symbols("T X m n")
    y0, y1, y2 = sp.symbols("y0 y1 y2")
    quartic = _quartic(parameter, x)

    # For slope != +/-1, coefficients T^6,T^5,T^4,T^3 recursively give
    # y3,y2,y1,y0.  The other sign negates the entire ordinate.
    y3_formula = 3 * (1 - slope**2)
    y2_formula = -6 * slope * (intercept - 177)
    y1_formula = -(
        435213 * slope**4
        + 9 * slope**2 * intercept**2
        - 3186 * slope**2 * intercept
        + 439041 * slope**2
        - 9 * intercept**2
        + 3186 * intercept
        - 455374
    ) / (3 * (slope - 1) * (slope + 1))
    y0_formula = -slope * (
        870426 * slope**4 * intercept
        - 157480281 * slope**4
        - 1740852 * slope**2 * intercept
        + 305854218 * slope**2
        + 32666 * intercept
        - 90417
    ) / (3 * (slope - 1) ** 2 * (slope + 1) ** 2)
    ordinate = (
        y3_formula * parameter**3
        + y2_formula * parameter**2
        + y1_formula * parameter
        + y0_formula
    )
    residual = sp.Poly(
        sp.together(
            quartic.subs(x, slope * parameter + intercept) - ordinate**2
        ),
        parameter,
    )
    equations = []
    for degree in (2, 1, 0):
        numerator = sp.factor(
            sp.together(residual.coeff_monomial(parameter**degree))
        ).as_numer_denom()[0]
        equations.append(numerator)
    resultants = (
        sp.Poly(sp.resultant(equations[0], equations[1], intercept), slope),
        sp.Poly(sp.resultant(equations[1], equations[2], intercept), slope),
        sp.Poly(sp.resultant(equations[0], equations[2], intercept), slope),
    )
    common = _monic(sp.gcd(sp.gcd(resultants[0], resultants[1]), resultants[2]))
    expected = _monic(
        sp.Poly(
            (27 * slope - 43)
            * (27 * slope - 17)
            * (27 * slope - 7)
            * (27 * slope + 7)
            * (27 * slope + 17)
            * (27 * slope + 43),
            slope,
        )
    )
    if common != expected:
        raise AssertionError("unexpected resultant gcd away from slopes +/-1")

    recovered = []
    for section in SECTION7_LINEAR_COMPANION_SECTIONS:
        candidate_slope = _sympy_rational(section.slope)
        specialized = [
            sp.Poly(equation.subs(slope, candidate_slope), intercept)
            for equation in equations
        ]
        common_intercept = _monic(
            sp.gcd(sp.gcd(specialized[0], specialized[1]), specialized[2])
        )
        expected_intercept = sp.Poly(
            intercept - _sympy_rational(section.intercept), intercept
        )
        if common_intercept != expected_intercept:
            raise AssertionError(f"failed to recover {section.label}")
        recovered.append(
            {
                "label": section.label,
                "slope": str(candidate_slope),
                "intercept": str(_sympy_rational(section.intercept)),
            }
        )

    # At slopes +/-1 the leading cubic term vanishes.  Direct Groebner
    # elimination recovers exactly the six visible Mestre intercepts.
    root_polynomial = sp.Poly(
        sp.prod(intercept - root for root in SECTION7_ROOTS), intercept
    )
    exceptional: dict[str, str] = {}
    for candidate_slope in (sp.Integer(-1), sp.Integer(1)):
        ordinate2 = y2 * parameter**2 + y1 * parameter + y0
        difference = sp.Poly(
            sp.expand(
                quartic.subs(x, candidate_slope * parameter + intercept)
                - ordinate2**2
            ),
            parameter,
        )
        equations2 = [
            difference.coeff_monomial(parameter**degree) for degree in range(5)
        ]
        basis = sp.groebner(
            equations2, y2, y1, y0, intercept, order="lex"
        )
        eliminants = [
            sp.Poly(item.as_expr(), intercept)
            for item in basis.polys
            if not any(
                item.as_expr().has(variable) for variable in (y2, y1, y0)
            )
        ]
        if len(eliminants) != 1 or _monic(eliminants[0]) != _monic(root_polynomial):
            raise AssertionError("slope +/-1 elimination missed a Mestre root")
        exceptional[str(candidate_slope)] = str(
            sp.factor(eliminants[0].as_expr())
        )

    return {
        "nonsingular_slope_resultant_gcd": str(sp.factor(common.as_expr())),
        "recovered_companion_sections": recovered,
        "slope_plus_or_minus_one_intercept_eliminants": exceptional,
        "classification_scope": (
            "all Q-rational polynomial sections with x=m*T+n and degree(y)<=3, "
            "up to negating y; not arbitrary rational sections"
        ),
    }


def classify_quadratic_abscissae() -> dict[str, Any]:
    """Classify ``deg(x)=2, deg(y)<=5`` beyond the linear branch."""

    parameter, x = sp.symbols("T X")
    quadratic, linear, constant, linear_square = sp.symbols("m n k z")
    ordinates = sp.symbols("y0:6")
    quartic = _quartic(parameter, x)
    abscissa = quadratic * parameter**2 + linear * parameter + constant
    ordinate = sum(
        ordinates[degree] * parameter**degree for degree in range(6)
    )
    residual = sp.Poly(
        sp.expand(quartic.subs(x, abscissa) - ordinate**2), parameter
    )

    # On the genuine quadratic branch m!=0, choose y5=3m^2.  The other
    # choice merely negates y.  Coefficients T^9,...,T^5 then determine all
    # remaining ordinate coefficients without further branching.
    substitutions: dict[sp.Symbol, sp.Expr] = {ordinates[5]: 3 * quadratic**2}
    if sp.factor(residual.coeff_monomial(parameter**10).subs(substitutions)) != 0:
        raise AssertionError("the quadratic leading coefficient failed")
    for degree, variable in zip(range(9, 4, -1), reversed(ordinates[:5])):
        equation = sp.factor(
            residual.coeff_monomial(parameter**degree).subs(substitutions)
        )
        solutions = sp.solve(equation, variable)
        if len(solutions) != 1:
            raise AssertionError("quadratic ordinate recursion branched")
        substitutions[variable] = sp.factor(solutions[0])

    equations = []
    for degree in range(4, -1, -1):
        equations.append(
            sp.factor(
                sp.together(
                    residual.coeff_monomial(parameter**degree).subs(
                        substitutions
                    )
                )
            ).as_numer_denom()[0]
        )

    # If n=0 the odd residual coefficients vanish.  The three even
    # coefficients have a two-polynomial Groebner basis giving exactly three
    # rational m-values and one common k-value.
    even_branch = [
        equations[index].subs(linear, 0) for index in (0, 2, 4)
    ]
    even_basis = sp.groebner(even_branch, constant, quadratic, order="lex")
    expected_constant = sp.Poly(5373 * constant - 1389190, constant, quadratic)
    expected_quadratic = sp.Poly(
        (5373 * quadratic - 56)
        * (5373 * quadratic + 22)
        * (5373 * quadratic + 34),
        constant,
        quadratic,
    )
    if len(even_basis.polys) != 2:
        raise AssertionError("unexpected even quadratic Groebner-basis length")
    if {
        _monic(sp.Poly(item.as_expr(), constant, quadratic))
        for item in even_basis.polys
    } != {_monic(expected_constant), _monic(expected_quadratic)}:
        raise AssertionError("the even quadratic section classification changed")

    # For n!=0, divide the odd residual coefficients by n and put z=n^2.
    # Every remaining polynomial is even in n.  A unit Groebner basis proves
    # that this branch has no solutions, even over an algebraic closure.
    nonzero_linear_branch = []
    for degree, equation in zip(range(4, -1, -1), equations):
        if degree in (3, 1):
            equation = sp.cancel(equation / linear)
        polynomial = sp.Poly(equation, linear)
        converted = sp.Integer(0)
        for (power,), coefficient in polynomial.as_dict().items():
            if power % 2:
                raise AssertionError("an odd power remained after extracting n")
            converted += coefficient * linear_square ** (power // 2)
        nonzero_linear_branch.append(sp.factor(converted))
    nonzero_basis = sp.groebner(
        nonzero_linear_branch,
        constant,
        linear_square,
        quadratic,
        order="lex",
        method="f5b",
    )
    if len(nonzero_basis.polys) != 1 or nonzero_basis.polys[0].as_expr() != 1:
        raise AssertionError("the n!=0 quadratic branch unexpectedly survived")

    return {
        "classification_scope": (
            "all Q-rational polynomial sections with degree(x)<=2 and "
            "degree(y)<=5, up to negating y; m=0 reduces to the separately "
            "proved linear classification"
        ),
        "genuine_quadratic_n_nonzero_groebner_basis": ["1"],
        "genuine_quadratic_n_zero_groebner_basis": [
            str(sp.factor(item.as_expr())) for item in even_basis.polys
        ],
        "recovered_quadratic_sections": [
            {
                "label": section.label,
                "quadratic_coefficient": str(section.quadratic_coefficient),
                "linear_coefficient": str(section.linear_coefficient),
                "constant_coefficient": str(section.constant_coefficient),
            }
            for section in SECTION7_QUADRATIC_COMPANION_SECTIONS
        ],
    }


def classify_cubic_abscissae() -> dict[str, Any]:
    """Prove that the genuine cubic-abscissa branch is empty.

    SymPy derives the exact residual ideal over ``ZZ``.  Singular computes its
    reduced Groebner basis over ``QQ``.  The basis forces the leading cubic
    coefficient to vanish, reducing every solution to the already classified
    degree-at-most-two branch.
    """

    parameter, x = sp.symbols("T X")
    cubic, quadratic, linear, constant = sp.symbols("a b c d")
    ordinates = sp.symbols("y0:8")
    quartic = _quartic(parameter, x)
    abscissa = (
        cubic * parameter**3
        + quadratic * parameter**2
        + linear * parameter
        + constant
    )
    ordinate = sum(
        ordinates[degree] * parameter**degree for degree in range(8)
    )
    residual = sp.Poly(
        sp.expand(quartic.subs(x, abscissa) - ordinate**2), parameter
    )
    substitutions: dict[sp.Symbol, sp.Expr] = {ordinates[7]: 3 * cubic**2}
    if sp.factor(residual.coeff_monomial(parameter**14).subs(substitutions)) != 0:
        raise AssertionError("the cubic leading coefficient failed")
    for degree, variable in zip(range(13, 6, -1), reversed(ordinates[:7])):
        equation = sp.factor(
            residual.coeff_monomial(parameter**degree).subs(substitutions)
        )
        solutions = sp.solve(equation, variable)
        if len(solutions) != 1:
            raise AssertionError("cubic ordinate recursion branched")
        substitutions[variable] = sp.factor(solutions[0])
    equations = [
        sp.factor(
            sp.together(
                residual.coeff_monomial(parameter**degree).subs(substitutions)
            )
        ).as_numer_denom()[0]
        for degree in range(6, -1, -1)
    ]
    expressions = [
        str(sp.expand(equation)).replace("**", "^") for equation in equations
    ]
    program = (
        "ring r=0,(d,c,b,a),dp;\n"
        f"ideal I={','.join(expressions)};\n"
        "option(redSB);\n"
        "ideal G=std(I);\n"
        'print("GBSIZE");\n'
        "size(G);\n"
        'print("GB");\n'
        "G;\n"
        "quit;\n"
    )
    try:
        completed = subprocess.run(
            ["Singular", "-q"],
            input=program,
            text=True,
            capture_output=True,
            timeout=30,
            check=False,
        )
    except FileNotFoundError as error:
        raise RuntimeError("Singular is required for the cubic classification") from error
    except subprocess.TimeoutExpired as error:
        raise RuntimeError("the capped cubic Singular proof timed out") from error
    if completed.returncode != 0 or completed.stderr.strip():
        raise RuntimeError(
            "Singular failed in the cubic classification: "
            + completed.stderr.strip()
        )
    lines = tuple(line.strip() for line in completed.stdout.splitlines() if line.strip())
    expected = (
        "GBSIZE",
        "4",
        "GB",
        "G[1]=a",
        "G[2]=c",
        "G[3]=5373d-1389190",
        "G[4]=155113830117b3-12830724b-41888",
    )
    if lines != expected:
        raise AssertionError(f"the exact cubic Groebner basis changed: {lines}")
    singular_version = subprocess.run(
        ["Singular", "--version"],
        text=True,
        capture_output=True,
        timeout=10,
        check=True,
    ).stdout.splitlines()[0]
    return {
        "classification_scope": (
            "the genuine degree-three branch x=a*T^3+b*T^2+c*T+d with "
            "a nonzero and degree(y)<=7; a=0 is the separately proved "
            "degree-at-most-two classification"
        ),
        "residual_equation_count": len(equations),
        "residual_ideal_sha256": hashlib.sha256(
            "\n".join(expressions).encode()
        ).hexdigest(),
        "singular_reduced_groebner_basis": list(expected[3:]),
        "leading_cubic_coefficient_forced_zero": True,
        "new_cubic_sections": 0,
        "singular_version": singular_version,
    }


def classify_quartic_abscissae() -> dict[str, Any]:
    """Prove that the genuine quartic-abscissa branch is empty."""

    parameter, x = sp.symbols("T X")
    coefficient4, coefficient3, coefficient2, coefficient1, coefficient0 = (
        sp.symbols("a b c d e")
    )
    ordinates = sp.symbols("y0:10")
    quartic = _quartic(parameter, x)
    abscissa = (
        coefficient4 * parameter**4
        + coefficient3 * parameter**3
        + coefficient2 * parameter**2
        + coefficient1 * parameter
        + coefficient0
    )
    ordinate = sum(
        ordinates[degree] * parameter**degree for degree in range(10)
    )
    residual = sp.Poly(
        sp.expand(quartic.subs(x, abscissa) - ordinate**2), parameter
    )
    substitutions: dict[sp.Symbol, sp.Expr] = {
        ordinates[9]: 3 * coefficient4**2
    }
    if sp.factor(residual.coeff_monomial(parameter**18).subs(substitutions)) != 0:
        raise AssertionError("the quartic-abscissa leading coefficient failed")
    for degree, variable in zip(range(17, 8, -1), reversed(ordinates[:9])):
        equation = sp.factor(
            residual.coeff_monomial(parameter**degree).subs(substitutions)
        )
        solutions = sp.solve(equation, variable)
        if len(solutions) != 1:
            raise AssertionError("quartic-abscissa ordinate recursion branched")
        substitutions[variable] = sp.factor(solutions[0])
    equations = [
        sp.factor(
            sp.together(
                residual.coeff_monomial(parameter**degree).subs(substitutions)
            )
        ).as_numer_denom()[0]
        for degree in range(8, -1, -1)
    ]
    expressions = [
        str(sp.expand(equation)).replace("**", "^") for equation in equations
    ]
    program = (
        "ring r=0,(e,d,c,b,a),dp;\n"
        f"ideal I={','.join(expressions)};\n"
        "option(redSB);\n"
        "ideal G=slimgb(I);\n"
        'print("GBSIZE");\n'
        "size(G);\n"
        'print("GB");\n'
        "G;\n"
        "quit;\n"
    )
    try:
        completed = subprocess.run(
            ["Singular", "-q"],
            input=program,
            text=True,
            capture_output=True,
            timeout=60,
            check=False,
        )
    except FileNotFoundError as error:
        raise RuntimeError("Singular is required for the quartic classification") from error
    except subprocess.TimeoutExpired as error:
        raise RuntimeError("the capped quartic Singular proof timed out") from error
    if completed.returncode != 0 or completed.stderr.strip():
        raise RuntimeError(
            "Singular failed in the quartic classification: "
            + completed.stderr.strip()
        )
    lines = tuple(line.strip() for line in completed.stdout.splitlines() if line.strip())
    expected = (
        "GBSIZE",
        "5",
        "GB",
        "G[1]=a",
        "G[2]=b",
        "G[3]=d",
        "G[4]=5373e-1389190",
        "G[5]=155113830117c3-12830724c-41888",
    )
    if lines != expected:
        raise AssertionError(f"the exact quartic Groebner basis changed: {lines}")
    return {
        "classification_scope": (
            "the genuine degree-four branch x=a*T^4+b*T^3+c*T^2+d*T+e "
            "with a nonzero and degree(y)<=9; a=0 reduces to the separately "
            "proved degree-at-most-three classification"
        ),
        "residual_equation_count": len(equations),
        "residual_ideal_sha256": hashlib.sha256(
            "\n".join(expressions).encode()
        ).hexdigest(),
        "singular_algorithm": "slimgb over QQ with degree ordering",
        "singular_reduced_groebner_basis": list(expected[3:]),
        "leading_quartic_coefficient_forced_zero": True,
        "new_quartic_sections": 0,
    }


def classify_quintic_abscissae() -> dict[str, Any]:
    """Prove that the genuine quintic-abscissa branch is empty."""

    parameter, x = sp.symbols("T X")
    coefficients = sp.symbols("a0:6")
    ordinates = sp.symbols("y0:12")
    quartic = _quartic(parameter, x)
    abscissa = sum(
        coefficients[degree] * parameter**degree for degree in range(6)
    )
    ordinate = sum(
        ordinates[degree] * parameter**degree for degree in range(12)
    )
    residual = sp.Poly(
        sp.expand(quartic.subs(x, abscissa) - ordinate**2), parameter
    )
    substitutions: dict[sp.Symbol, sp.Expr] = {
        ordinates[11]: 3 * coefficients[5] ** 2
    }
    if sp.factor(residual.coeff_monomial(parameter**22).subs(substitutions)) != 0:
        raise AssertionError("the quintic-abscissa leading coefficient failed")
    for degree, variable in zip(range(21, 10, -1), reversed(ordinates[:11])):
        equation = sp.factor(
            residual.coeff_monomial(parameter**degree).subs(substitutions)
        )
        solutions = sp.solve(equation, variable)
        if len(solutions) != 1:
            raise AssertionError("quintic-abscissa ordinate recursion branched")
        substitutions[variable] = sp.factor(solutions[0])
    equations = [
        sp.factor(
            sp.together(
                residual.coeff_monomial(parameter**degree).subs(substitutions)
            )
        ).as_numer_denom()[0]
        for degree in range(10, -1, -1)
    ]
    expressions = [
        str(sp.expand(equation)).replace("**", "^") for equation in equations
    ]
    program = (
        "ring r=0,(a0,a1,a2,a3,a4,a5),dp;\n"
        f"ideal I={','.join(expressions)};\n"
        "option(redSB);\n"
        "ideal G=slimgb(I);\n"
        'print("GBSIZE");\n'
        "size(G);\n"
        'print("GB");\n'
        "G;\n"
        "quit;\n"
    )
    try:
        completed = subprocess.run(
            ["Singular", "-q"],
            input=program,
            text=True,
            capture_output=True,
            timeout=60,
            check=False,
        )
    except FileNotFoundError as error:
        raise RuntimeError("Singular is required for the quintic classification") from error
    except subprocess.TimeoutExpired as error:
        raise RuntimeError("the capped quintic Singular proof timed out") from error
    if completed.returncode != 0 or completed.stderr.strip():
        raise RuntimeError(
            "Singular failed in the quintic classification: "
            + completed.stderr.strip()
        )
    lines = tuple(line.strip() for line in completed.stdout.splitlines() if line.strip())
    expected = (
        "GBSIZE",
        "6",
        "GB",
        "G[1]=a5",
        "G[2]=a4",
        "G[3]=a3",
        "G[4]=a1",
        "G[5]=5373*a0-1389190",
        "G[6]=155113830117*a2^3-12830724*a2-41888",
    )
    if lines != expected:
        raise AssertionError(f"the exact quintic Groebner basis changed: {lines}")
    return {
        "classification_scope": (
            "the genuine degree-five branch x=sum(a_i*T^i,i=0..5) with "
            "a5 nonzero and degree(y)<=11; a5=0 reduces to the separately "
            "proved degree-at-most-four classification"
        ),
        "residual_equation_count": len(equations),
        "residual_ideal_sha256": hashlib.sha256(
            "\n".join(expressions).encode()
        ).hexdigest(),
        "singular_algorithm": "slimgb over QQ with degree ordering",
        "singular_reduced_groebner_basis": list(expected[3:]),
        "leading_quintic_coefficient_forced_zero": True,
        "new_quintic_sections": 0,
    }


def _visible_symbolic_points(
    parameter: sp.Symbol,
) -> tuple[tuple[sp.Expr, sp.Expr], ...]:
    points = []
    for index in range(12):
        values = [
            (
                sample,
                _sympy_rational(
                    primitive_visible_points(
                        SECTION7_CONSTRUCTION, Q(sample)
                    )[index][1]
                ),
            )
            for sample in range(1, 9)
        ]
        ordinate = sp.expand(sp.interpolate(values, parameter))
        _, check = primitive_visible_points(SECTION7_CONSTRUCTION, Q(9))[index]
        if ordinate.subs(parameter, 9) != _sympy_rational(check):
            raise AssertionError("visible ordinate interpolation failed")
        root = _sympy_rational(SECTION7_CONSTRUCTION.roots[index // 2])
        sign = -1 if index % 2 == 0 else 1
        points.append((root + sign * parameter, ordinate))
    return tuple(points)


def verify_symbolic_relations() -> list[dict[str, Any]]:
    """Prove all displayed generic section identities and dependencies."""

    parameter, x_symbol = sp.symbols("T X")
    quartic = _quartic(parameter, x_symbol)
    visible_quartic = _visible_symbolic_points(parameter)
    visible_jacobian = []
    coefficient_a: sp.Expr | None = None
    coefficient_b: sp.Expr | None = None
    for point in visible_quartic:
        mapped_x, mapped_y, coefficient_a, coefficient_b = _covariant_map(
            quartic, parameter, x_symbol, point
        )
        visible_jacobian.append((mapped_x, mapped_y))
    assert coefficient_a is not None and coefficient_b is not None

    companion_jacobian: dict[str, tuple[sp.Expr, sp.Expr]] = {}
    for section in SECTION7_LINEAR_COMPANION_SECTIONS:
        point_x = (
            _sympy_rational(section.slope) * parameter
            + _sympy_rational(section.intercept)
        )
        point_y = sum(
            _sympy_rational(coefficient) * parameter**degree
            for degree, coefficient in enumerate(section.ordinate_coefficients)
        )
        if sp.expand(quartic.subs(x_symbol, point_x) - point_y**2) != 0:
            raise AssertionError(f"quartic identity failed for {section.label}")
        mapped_x, mapped_y, _, _ = _covariant_map(
            quartic, parameter, x_symbol, (point_x, point_y)
        )
        if (
            sp.cancel(
                mapped_y**2
                - mapped_x**3
                - coefficient_a * mapped_x
                - coefficient_b
            )
            != 0
        ):
            raise AssertionError(f"Jacobian identity failed for {section.label}")
        companion_jacobian[section.label] = (mapped_x, mapped_y)

    quadratic_jacobian: dict[str, tuple[sp.Expr, sp.Expr]] = {}
    for section in SECTION7_QUADRATIC_COMPANION_SECTIONS:
        point_x = (
            _sympy_rational(section.quadratic_coefficient) * parameter**2
            + _sympy_rational(section.linear_coefficient) * parameter
            + _sympy_rational(section.constant_coefficient)
        )
        point_y = sum(
            _sympy_rational(coefficient) * parameter**degree
            for degree, coefficient in enumerate(section.ordinate_coefficients)
        )
        if sp.expand(quartic.subs(x_symbol, point_x) - point_y**2) != 0:
            raise AssertionError(f"quartic identity failed for {section.label}")
        mapped_x, mapped_y, _, _ = _covariant_map(
            quartic, parameter, x_symbol, (point_x, point_y)
        )
        if (
            sp.cancel(
                mapped_y**2
                - mapped_x**3
                - coefficient_a * mapped_x
                - coefficient_b
            )
            != 0
        ):
            raise AssertionError(f"Jacobian identity failed for {section.label}")
        quadratic_jacobian[section.label] = (mapped_x, mapped_y)

    basis = tuple(
        visible_jacobian[index]
        for index in SECTION7_RELATION_BASIS_VISIBLE_INDICES
    ) + (companion_jacobian["plus-7/27"],)
    targets = {
        "visible-11": visible_jacobian[11],
        **{
            label: point
            for label, point in companion_jacobian.items()
            if label != "plus-7/27"
        },
        **quadratic_jacobian,
    }
    records = []
    relations = {
        **SECTION7_JACOBIAN_RELATIONS,
        **SECTION7_QUADRATIC_JACOBIAN_RELATIONS,
    }
    for label, coefficients in relations.items():
        total: tuple[sp.Expr, sp.Expr] | None = None
        for point, coefficient in zip(basis, coefficients):
            if coefficient == 0:
                continue
            if coefficient not in (-1, 1):
                raise AssertionError("the pinned symbolic relations must be unit-valued")
            summand = point if coefficient == 1 else (point[0], -point[1])
            total = _elliptic_add(total, summand, coefficient_a)
        if total is None:
            raise AssertionError(f"relation for {label} collapsed to infinity")
        target = targets[label]
        if (
            sp.cancel(total[0] - target[0]) != 0
            or sp.cancel(total[1] - target[1]) != 0
        ):
            raise AssertionError(f"symbolic relation failed for {label}")
        records.append(
            {
                "section": label,
                "basis_coefficients": list(coefficients),
                "exact_identity_in_Q_of_T": True,
            }
        )
    return records


def exact_specialization_independence() -> dict[str, Any]:
    """Certify twelve independent sections at T=1 using mod-3 reductions."""

    parameter = SPECIALIZATION_PARAMETER
    visible = tuple(
        quartic_point_to_short_jacobian(SECTION7_CONSTRUCTION, parameter, point)
        for point in primitive_visible_points(SECTION7_CONSTRUCTION, parameter)
    )
    points = visible[:11] + (
        SECTION7_LINEAR_COMPANION_SECTIONS[0].jacobian_point(parameter),
    )
    coefficients = short_jacobian_coefficients(SECTION7_CONSTRUCTION, parameter)
    signatures = find_mod_l_reduction_certificate(
        coefficients, points, modulus=INDEPENDENCE_MODULUS, prime_bound=200
    )
    primes = tuple(signature.prime for signature in signatures)
    rank = combined_mod_l_rank(signatures, len(points), INDEPENDENCE_MODULUS)
    torsion_prime = find_no_rational_l_torsion_prime(
        coefficients, modulus=INDEPENDENCE_MODULUS, prime_bound=200
    )
    if primes != EXPECTED_CERTIFICATE_PRIMES or rank != len(points):
        raise AssertionError("the exact mod-3 independence certificate changed")
    if torsion_prime != EXPECTED_TORSION_CERTIFICATE_PRIME:
        raise AssertionError("the rational 3-torsion certificate prime changed")
    return {
        "parameter_T": str(parameter),
        "basis": "visible indices 0,...,10 followed by plus-7/27",
        "basis_size": len(points),
        "modulus": INDEPENDENCE_MODULUS,
        "combined_column_rank": rank,
        "no_rational_3_torsion_certificate_prime": torsion_prime,
        "reduction_signatures": [
            {
                "prime": signature.prime,
                "group_order": signature.group_order,
                "triple_subgroup_order": signature.multiple_subgroup_order,
                "quotient_dimension": signature.quotient_dimension,
                "rows": [list(row) for row in signature.rows],
            }
            for signature in signatures
        ],
        "inference": (
            "full rank modulo products of E(F_p)/3E(F_p), together with "
            "E(Q)[3]=0, proves the specialized points independent by infinite "
            "descent; any generic relation would specialize at T=1, so the "
            "twelve sections are independent over Q(T)"
        ),
        "generic_rank_lower_bound": 12,
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts/generated-results/elliptic_nagao_section7_linear_sections.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    classification = classify_abscissae()
    fiber_geometry = verify_k3_fiber_geometry()
    quadratic_classification = classify_quadratic_abscissae()
    cubic_classification = classify_cubic_abscissae()
    quartic_classification = classify_quartic_abscissae()
    quintic_classification = classify_quintic_abscissae()
    relations = verify_symbolic_relations()
    independence = exact_specialization_independence()
    script_path = Path(__file__).resolve()
    data_path = script_path.with_name("nagao_1994_section7.py")
    mod_l_path = script_path.with_name("mod_l_reduction_independence.py")
    helper_path = script_path.with_name("verify_nagao_linear_sections.py")
    artifact = {
        "schema_version": 1,
        "status": "exact symbolic theorem and finite-reduction certificate",
        "primary_source": PRIMARY_SOURCE,
        "classification": classification,
        "k3_fiber_geometry": fiber_geometry,
        "quadratic_classification": quadratic_classification,
        "cubic_classification": cubic_classification,
        "quartic_abscissa_classification": quartic_classification,
        "quintic_abscissa_classification": quintic_classification,
        "companion_sections": [
            {
                "label": section.label,
                "slope": str(section.slope),
                "intercept": str(section.intercept),
                "ordinate_coefficients_ascending": [
                    str(value) for value in section.ordinate_coefficients
                ],
            }
            for section in SECTION7_LINEAR_COMPANION_SECTIONS
        ],
        "quadratic_companion_sections": [
            {
                "label": section.label,
                "quadratic_coefficient": str(section.quadratic_coefficient),
                "linear_coefficient": str(section.linear_coefficient),
                "constant_coefficient": str(section.constant_coefficient),
                "ordinate_coefficients_ascending": [
                    str(value) for value in section.ordinate_coefficients
                ],
            }
            for section in SECTION7_QUADRATIC_COMPANION_SECTIONS
        ],
        "exact_jacobian_relations": relations,
        "exact_specialization_independence": independence,
        "proved_consequences": {
            "generic_rank_at_least": 12,
            "geometric_generic_rank_at_most": 15,
            "linear_abscissa_sections_classified": 18,
            "polynomial_abscissa_degree_at_most_two_sections_classified": 21,
            "polynomial_abscissa_degree_at_most_three_sections_classified": 21,
            "polynomial_abscissa_degree_at_most_four_sections_classified": 21,
            "polynomial_abscissa_degree_at_most_five_sections_classified": 21,
            "generic_companions_dependent_on_rank12_basis": 8,
        },
        "scope_limits": [
            "not a classification of abscissae of degree six or higher or arbitrary rational sections",
            "not a proof that the generic Mordell--Weil rank equals 12",
            "not a rank-21 or rank-30 target hit",
        ],
        "target_hit": False,
        "software": {"python": platform.python_version(), "sympy": sp.__version__},
        "reproducing_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        "data_sha256": hashlib.sha256(data_path.read_bytes()).hexdigest(),
        "mod_l_engine_sha256": hashlib.sha256(mod_l_path.read_bytes()).hexdigest(),
        "symbolic_helper_sha256": hashlib.sha256(helper_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print("classified 21 polynomial sections with abscissa degree at most five")
    print("proved nine exact generic Mordell--Weil relations")
    print("proved generic rank at least 12 via exact mod-3 reductions at T=1")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
