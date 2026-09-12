#!/usr/bin/env python3
"""Generate every displayed explicit-quintic interface from one JSON spec."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.keller_fiber import quadratic_gauge_map  # noqa: E402

DEFAULT_SPEC = (
    ROOT / "papers/common-arithmetic-fibers/data/explicit-quintic-spec.json"
)


def rational(value: Any) -> sp.Rational:
    q = Fraction(str(value))
    return sp.Rational(q.numerator, q.denominator)


def rational_text(value: Any) -> str:
    q = Fraction(value)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def coefficients(poly: sp.Poly) -> list[str]:
    return [
        rational_text(poly.coeff_monomial(poly.gens[0] ** degree))
        for degree in range(poly.degree() + 1)
    ]


def polynomial_from_coefficients(
    raw_coefficients: Sequence[Any], variable: sp.Symbol
) -> sp.Expr:
    return sp.expand(
        sum(rational(value) * variable**degree for degree, value in enumerate(raw_coefficients))
    )


def sparse_terms(expression: sp.Expr, variables: Sequence[sp.Symbol]) -> list[dict[str, Any]]:
    poly = sp.Poly(sp.expand(expression), *variables, domain=sp.QQ)
    return [
        {
            "coefficient": rational_text(coefficient),
            "exponents": list(exponents),
        }
        for exponents, coefficient in sorted(poly.terms(), key=lambda term: term[0])
    ]


def sparse_serialization(mapping: Sequence[sp.Expr], variables: Sequence[sp.Symbol]) -> str:
    payload = [sparse_terms(component, variables) for component in mapping]
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def lean_rational(value: Any) -> str:
    q = Fraction(value)
    if q.denominator == 1:
        return f"({q.numerator} : ℚ)"
    return f"(({q.numerator} : ℚ) / ({q.denominator} : ℚ))"


def lean_display_polynomial(raw_coefficients: Sequence[Any]) -> str:
    """Render the conventional descending polynomial syntax used in the paper."""
    pieces: list[str] = []
    for degree in reversed(range(len(raw_coefficients))):
        q = Fraction(str(raw_coefficients[degree]))
        if q == 0:
            continue
        magnitude = abs(q)
        if degree == 0:
            body = rational_text(magnitude)
        else:
            coefficient = "" if magnitude == 1 else f"{rational_text(magnitude)} * "
            power = "Polynomial.X" if degree == 1 else f"Polynomial.X ^ {degree}"
            body = coefficient + power
        if not pieces:
            pieces.append(("-" if q < 0 else "") + body)
        else:
            pieces.append((" - " if q < 0 else " + ") + body)
    return "".join(pieces) if pieces else "0"


def tex_rational(value: Any) -> str:
    q = Fraction(value)
    if q.denominator == 1:
        return str(q.numerator)
    return rf"\frac{{{q.numerator}}}{{{q.denominator}}}"


def signed_polynomial_tex(raw_coefficients: Sequence[Any], variable: str) -> str:
    pieces: list[str] = []
    for degree in reversed(range(len(raw_coefficients))):
        q = Fraction(str(raw_coefficients[degree]))
        if q == 0:
            continue
        magnitude = abs(q)
        if degree == 0:
            body = tex_rational(magnitude)
        else:
            coefficient = "" if magnitude == 1 else tex_rational(magnitude)
            power = variable if degree == 1 else rf"{variable}^{{{degree}}}"
            body = coefficient + power
        if not pieces:
            pieces.append(("-" if q < 0 else "") + body)
        else:
            pieces.append((" - " if q < 0 else " + ") + body)
    return "".join(pieces) if pieces else "0"


def vector_tex(values: Sequence[Any]) -> str:
    return "(" + ",".join(tex_rational(value) for value in values) + ")"


def build(specification: Mapping[str, Any]) -> dict[str, Any]:
    inverse_name = specification["variables"]["inverse"]
    source_names = specification["variables"]["source"]
    S = sp.Symbol(inverse_name)
    x, y, z = sp.symbols(" ".join(source_names))

    factors = [
        polynomial_from_coefficients(raw, S)
        for raw in specification["factor_coefficients_ascending"]
    ]
    polynomial = sp.Poly(sp.prod(factors), S, domain=sp.QQ)
    translation = rational(specification["translation"])
    translated = sp.expand(polynomial.as_expr().subs(S, S + translation))
    seed = sp.expand(translated - polynomial.as_expr().subs(S, translation))
    seed_poly = sp.Poly(seed, S, domain=sp.QQ)
    seed_coefficients = coefficients(seed_poly)
    g1 = seed_poly.coeff_monomial(S)
    g2 = seed_poly.coeff_monomial(S**2)
    g3 = seed_poly.coeff_monomial(S**3)
    if polynomial.degree() < 3 or g1 == 0 or g3 == 0:
        raise ValueError("the canonical example is not an admissible quadratic-gauge seed")

    normalized_map = quadratic_gauge_map(seed, S, (x, y, z))
    integral_scaling = tuple(rational(v) for v in specification["integral_output_scaling"])
    jacobian_one_scaling = tuple(
        rational(v) for v in specification["jacobian_one_output_scaling"]
    )
    integral_to_jacobian_one_scaling = tuple(
        sp.cancel(final / displayed)
        for final, displayed in zip(jacobian_one_scaling, integral_scaling)
    )
    integral_map = tuple(
        sp.expand(scale * component)
        for scale, component in zip(integral_scaling, normalized_map)
    )
    jacobian_one_map = tuple(
        sp.expand(scale * component)
        for scale, component in zip(jacobian_one_scaling, normalized_map)
    )

    normalized_target = (
        sp.Integer(1),
        sp.Integer(0),
        sp.cancel(
            -2 * polynomial.as_expr().subs(S, translation) / g1
        ),
    )
    integral_target = tuple(
        sp.cancel(scale * value)
        for scale, value in zip(integral_scaling, normalized_target)
    )
    jacobian_one_target = tuple(
        sp.cancel(scale * value)
        for scale, value in zip(jacobian_one_scaling, normalized_target)
    )
    inverse_polynomial = sp.expand(
        seed
        - g1
        * (normalized_target[1] * S**2 + normalized_target[2])
        / 2
    )
    if inverse_polynomial != translated:
        raise AssertionError("derived inverse polynomial does not recover the translated input")

    t_symbol, q_symbol = sp.symbols("t q")
    q_in_t = t_symbol**2 * z + (g1 / g3) * y**2 * (1 + 3 * t_symbol)
    factorized_integral = (
        integral_scaling[0] * t_symbol * q_symbol,
        integral_scaling[1] * y
        + integral_scaling[1] * 3 * (g3 / g1) * x * q_symbol
        + integral_scaling[1] * 2 * (g2 / g1) * t_symbol * q_symbol
        + sum(
            integral_scaling[1]
            * index
            * (seed_poly.coeff_monomial(S**index) / g1)
            * t_symbol**2
            * x ** (index - 2)
            * q_symbol**index
            for index in range(4, polynomial.degree() + 1)
        ),
        integral_scaling[2] * x * (5 - 3 * t_symbol)
        - integral_scaling[2] * (g3 / g1) * x**3 * z
        - sum(
            integral_scaling[2]
            * (index - 2)
            * (seed_poly.coeff_monomial(S**index) / g1)
            * x**index
            * q_symbol**index
            for index in range(4, polynomial.degree() + 1)
        ),
    )
    reconstructed_integral = tuple(
        sp.expand(
            component.subs(q_symbol, q_in_t).subs(t_symbol, 1 + x * y)
        )
        for component in factorized_integral
    )
    if reconstructed_integral != integral_map:
        raise AssertionError("factorized paper map and shared SymPy compiler disagree")

    display_coefficients = {
        "q_y_squared": sp.cancel(g1 / g3),
        "f2_y": integral_scaling[1],
        "f2_xq": sp.cancel(integral_scaling[1] * 3 * g3 / g1),
        "f2_tq": sp.cancel(integral_scaling[1] * 2 * g2 / g1),
        "f2_high": {
            index: sp.cancel(
                integral_scaling[1]
                * index
                * seed_poly.coeff_monomial(S**index)
                / g1
            )
            for index in range(4, polynomial.degree() + 1)
        },
        "f3_x_factor": integral_scaling[2],
        "f3_x3z": sp.cancel(-integral_scaling[2] * g3 / g1),
        "f3_high": {
            index: sp.cancel(
                -integral_scaling[2]
                * (index - 2)
                * seed_poly.coeff_monomial(S**index)
                / g1
            )
            for index in range(4, polynomial.degree() + 1)
        },
    }

    normalized_jacobian = sp.factor(
        sp.Matrix(normalized_map).jacobian((x, y, z)).det()
    )
    integral_jacobian = sp.factor(
        sp.Matrix(integral_map).jacobian((x, y, z)).det()
    )
    jacobian_one_determinant = sp.factor(
        sp.Matrix(jacobian_one_map).jacobian((x, y, z)).det()
    )
    if (normalized_jacobian, integral_jacobian, jacobian_one_determinant) != (
        -2,
        -722,
        1,
    ):
        raise AssertionError("unexpected Jacobian determinant")

    serialization = sparse_serialization(integral_map, (x, y, z))
    certificate = {
        "base_field": specification["base_field"],
        "canonical_specification": str(DEFAULT_SPEC.relative_to(ROOT)),
        "factor_coefficients_ascending": [
            [rational_text(value) for value in raw]
            for raw in specification["factor_coefficients_ascending"]
        ],
        "integral_map": {
            "expanded_sparse_terms": [
                sparse_terms(component, (x, y, z)) for component in integral_map
            ],
            "jacobian_determinant": rational_text(integral_jacobian),
            "output_scaling": [rational_text(value) for value in integral_scaling],
            "sha256": hashlib.sha256(serialization.encode("utf-8")).hexdigest(),
            "term_counts": [
                len(sparse_terms(component, (x, y, z))) for component in integral_map
            ],
        },
        "inverse_polynomial": {
            "coefficients_ascending": coefficients(
                sp.Poly(inverse_polynomial, S, domain=sp.QQ)
            ),
            "formula": "G(S) - g1 * (B*S^2 + C) / 2",
        },
        "jacobian_one_map": {
            "integral_to_jacobian_one_scaling": [
                rational_text(value)
                for value in integral_to_jacobian_one_scaling
            ],
            "jacobian_determinant": rational_text(jacobian_one_determinant),
            "output_scaling": [
                rational_text(value) for value in jacobian_one_scaling
            ],
            "target": [rational_text(value) for value in jacobian_one_target],
        },
        "lean_correspondence": {
            "generated_module": specification["generated_outputs"]["lean"],
            "proof_module": (
                "formal/finite-etale-keller/FiniteEtaleKeller/"
                "PaperExampleCorrespondence.lean"
            ),
            "theorems": [
                "p5_eq_generatedPaperP5",
                "g5_eq_generatedPaperG5",
                "integralMap_eq_generatedPaperMap",
                "generatedPaperInversePolynomial",
                "generatedPaperTarget_scaling",
                "integralFiberPoint_eval_eq_generatedTarget",
                "generatedPaperJacobianOneTarget_scaling",
                "generatedPaper_map_scalings",
                "generatedPaper_integral_to_jacobianOne",
            ],
        },
        "name": specification["name"],
        "normalized_map": {
            "jacobian_determinant": rational_text(normalized_jacobian),
            "target": [rational_text(value) for value in normalized_target],
        },
        "polynomial": {
            "coefficients_ascending": coefficients(polynomial),
            "degree": polynomial.degree(),
        },
        "provenance": {
            "generator": "scripts/compile_common_arithmetic_fibers_example.py",
            "regeneration_command": (
                ".venv/bin/python "
                "scripts/compile_common_arithmetic_fibers_example.py"
            ),
            "verification_command": (
                ".venv/bin/python "
                "scripts/verify_common_arithmetic_fibers_correspondence.py"
            ),
        },
        "schema": "common-arithmetic-fibers/explicit-example-certificate/v1",
        "seed": {
            "coefficients_ascending": seed_coefficients,
            "linear_coefficient": rational_text(g1),
            "third_coefficient": rational_text(g3),
            "translation": rational_text(translation),
        },
        "tex_macros": specification["generated_outputs"]["tex"],
    }
    return {
        "certificate": certificate,
        "display_coefficients": display_coefficients,
        "factorized_integral": factorized_integral,
        "factors": factors,
        "integral_map": integral_map,
        "integral_scaling": integral_scaling,
        "integral_target": integral_target,
        "integral_to_jacobian_one_scaling": integral_to_jacobian_one_scaling,
        "inverse_polynomial": inverse_polynomial,
        "jacobian_one_scaling": jacobian_one_scaling,
        "jacobian_one_target": jacobian_one_target,
        "normalized_target": normalized_target,
        "polynomial": polynomial,
        "q_in_t": q_in_t,
        "seed": seed_poly,
        "specification": specification,
        "symbols": (S, x, y, z, t_symbol, q_symbol),
    }


def render_lean(data: Mapping[str, Any]) -> str:
    specification = data["specification"]
    namespace = specification["lean_namespace"]
    factor_coefficients = data["certificate"]["factor_coefficients_ascending"]
    p_coefficients = data["certificate"]["polynomial"]["coefficients_ascending"]
    g_coefficients = data["certificate"]["seed"]["coefficients_ascending"]
    integral_scaling = data["integral_scaling"]
    normalized_target = data["normalized_target"]
    integral_target = data["integral_target"]
    jacobian_one_scaling = data["jacobian_one_scaling"]
    integral_to_jacobian_one_scaling = data["integral_to_jacobian_one_scaling"]
    jacobian_one_target = data["jacobian_one_target"]
    display = data["display_coefficients"]
    q_y_squared = Fraction(display["q_y_squared"])
    f2_y = Fraction(display["f2_y"])
    f2_xq = Fraction(display["f2_xq"])
    f2_tq = Fraction(display["f2_tq"])
    f2_high = {
        index: Fraction(value) for index, value in display["f2_high"].items()
    }
    f3_x_factor = Fraction(display["f3_x_factor"])
    f3_x3z = Fraction(display["f3_x3z"])
    f3_high = {
        index: Fraction(value) for index, value in display["f3_high"].items()
    }
    expected_display = (
        q_y_squared,
        f2_y,
        f2_xq,
        f2_tq,
        f2_high,
        f3_x_factor,
        f3_x3z,
        f3_high,
    )
    if expected_display != (
        Fraction(-19),
        Fraction(19),
        Fraction(-3),
        Fraction(38),
        {4: Fraction(-4), 5: Fraction(-5)},
        Fraction(19),
        Fraction(1),
        {4: Fraction(2), 5: Fraction(3)},
    ):
        raise ValueError("the explicit-quintic Lean display shape is no longer applicable")
    return f"""\
/- This file is generated by scripts/compile_common_arithmetic_fibers_example.py. -/
import Mathlib

/-!
# Canonical data for the paper's explicit quintic

Do not edit this file by hand.  Its JSON source is
`papers/common-arithmetic-fibers/data/explicit-quintic-spec.json`.
-/

noncomputable section

open MvPolynomial

namespace {namespace}

abbrev M := MvPolynomial (Fin 3) ℚ

def paperP5 : Polynomial ℚ :=
  ({lean_display_polynomial(factor_coefficients[0])}) *
    ({lean_display_polynomial(factor_coefficients[1])})

def paperG5 : Polynomial ℚ :=
  {lean_display_polynomial(g_coefficients)}

def paperT : M :=
  1 + MvPolynomial.X 0 * MvPolynomial.X 1

def paperQ : M :=
  paperT ^ 2 * MvPolynomial.X 2
    - MvPolynomial.C {rational_text(abs(q_y_squared))} * MvPolynomial.X 1 ^ 2 *
      (1 + MvPolynomial.C 3 * paperT)

def paperIntegralMap : Fin 3 → M :=
  ![paperT * paperQ,
    MvPolynomial.C {rational_text(f2_y)} * MvPolynomial.X 1
      - MvPolynomial.C {rational_text(abs(f2_xq))} * MvPolynomial.X 0 * paperQ
      + MvPolynomial.C {rational_text(f2_tq)} * paperT * paperQ
      - MvPolynomial.C {rational_text(abs(f2_high[4]))} * paperT ^ 2 *
        MvPolynomial.X 0 ^ 2 * paperQ ^ 4
      - MvPolynomial.C {rational_text(abs(f2_high[5]))} * paperT ^ 2 *
        MvPolynomial.X 0 ^ 3 * paperQ ^ 5,
    MvPolynomial.C {rational_text(f3_x_factor)} * MvPolynomial.X 0 *
        (MvPolynomial.C 5 - MvPolynomial.C 3 * paperT)
      + MvPolynomial.X 0 ^ 3 * MvPolynomial.X 2
      + MvPolynomial.C {rational_text(f3_high[4])} *
        (MvPolynomial.X 0 * paperQ) ^ 4
      + MvPolynomial.C {rational_text(f3_high[5])} *
        (MvPolynomial.X 0 * paperQ) ^ 5]

def paperIntegralOutputScaling : Fin 3 → ℚ :=
  ![{", ".join(lean_rational(value) for value in integral_scaling)}]

def paperJacobianOneOutputScaling : Fin 3 → ℚ :=
  ![{", ".join(lean_rational(value) for value in jacobian_one_scaling)}]

def paperIntegralToJacobianOneScaling : Fin 3 → ℚ :=
  ![{", ".join(lean_rational(value) for value in integral_to_jacobian_one_scaling)}]

def paperNormalizedTarget : Fin 3 → ℚ :=
  ![{", ".join(lean_rational(value) for value in normalized_target)}]

def paperIntegralTarget : Fin 3 → ℚ :=
  ![{", ".join(lean_rational(value) for value in integral_target)}]

def paperJacobianOneTarget : Fin 3 → ℚ :=
  ![{", ".join(lean_rational(value) for value in jacobian_one_target)}]

end {namespace}
"""


def render_tex(data: Mapping[str, Any]) -> str:
    p_coefficients = data["certificate"]["polynomial"]["coefficients_ascending"]
    g_coefficients = data["certificate"]["seed"]["coefficients_ascending"]
    factorized = data["factorized_integral"]
    _, x, y, z, t, q = data["symbols"]
    integral_scaling = data["integral_scaling"]
    normalized_target = data["normalized_target"]
    integral_target = data["integral_target"]
    jacobian_one_scaling = data["jacobian_one_scaling"]
    integral_to_jacobian_one_scaling = data["integral_to_jacobian_one_scaling"]
    jacobian_one_target = data["jacobian_one_target"]
    p_factorized = "".join(
        f"({signed_polynomial_tex(raw, 'T')})"
        for raw in data["certificate"]["factor_coefficients_ascending"]
    )
    q_tex = sp.latex(data["q_in_t"])
    map_tex = [sp.latex(component) for component in data["factorized_integral"]]
    g1 = Fraction(data["certificate"]["seed"]["linear_coefficient"])
    normalized_c = Fraction(normalized_target[2])
    inverse_shift = g1 * normalized_c / 2
    inverse_identity = (
        rf"G(S)-\frac{{{tex_rational(g1)}}}{{2}}"
        rf"({tex_rational(normalized_c)})"
        rf"=G(S)-{tex_rational(inverse_shift)}=P_5(S)"
    )
    degree = data["certificate"]["polynomial"]["degree"]
    normalized_jacobian = data["certificate"]["normalized_map"][
        "jacobian_determinant"
    ]
    integral_jacobian = data["certificate"]["integral_map"][
        "jacobian_determinant"
    ]
    jacobian_one_determinant = data["certificate"]["jacobian_one_map"][
        "jacobian_determinant"
    ]
    return f"""\
% Generated by scripts/compile_common_arithmetic_fibers_example.py.
% Canonical source: papers/common-arithmetic-fibers/data/explicit-quintic-spec.json
\\newcommand{{\\CAFExplicitPolynomialFactorized}}{{{p_factorized}}}
\\newcommand{{\\CAFExplicitPolynomialExpanded}}{{{signed_polynomial_tex(p_coefficients, "T")}}}
\\newcommand{{\\CAFExplicitSeed}}{{{signed_polynomial_tex(g_coefficients, "S")}}}
\\newcommand{{\\CAFExplicitT}}{{1+xy}}
\\newcommand{{\\CAFExplicitQ}}{{{q_tex}}}
\\newcommand{{\\CAFExplicitFOne}}{{{map_tex[0]}}}
\\newcommand{{\\CAFExplicitFTwo}}{{{map_tex[1]}}}
\\newcommand{{\\CAFExplicitFThree}}{{{map_tex[2]}}}
\\newcommand{{\\CAFExplicitNormalizedTarget}}{{{vector_tex(normalized_target)}}}
\\newcommand{{\\CAFExplicitIntegralScaling}}{{{vector_tex(integral_scaling)}}}
\\newcommand{{\\CAFExplicitIntegralTarget}}{{{vector_tex(integral_target)}}}
\\newcommand{{\\CAFExplicitJacobianOneScaling}}{{{vector_tex(jacobian_one_scaling)}}}
\\newcommand{{\\CAFExplicitIntegralToJacobianOneScaling}}{{{vector_tex(integral_to_jacobian_one_scaling)}}}
\\newcommand{{\\CAFExplicitJacobianOneTarget}}{{{vector_tex(jacobian_one_target)}}}
\\newcommand{{\\CAFExplicitNormalizedJacobian}}{{{normalized_jacobian}}}
\\newcommand{{\\CAFExplicitIntegralJacobian}}{{{integral_jacobian}}}
\\newcommand{{\\CAFExplicitJacobianOneJacobian}}{{{jacobian_one_determinant}}}
\\newcommand{{\\CAFExplicitDegree}}{{{degree}}}
\\newcommand{{\\CAFExplicitInverseIdentity}}{{{inverse_identity}}}
"""


def python_expression(expression: sp.Expr) -> str:
    return sp.sstr(expression)


def render_sympy(data: Mapping[str, Any]) -> str:
    certificate = data["certificate"]
    factor_coefficients = certificate["factor_coefficients_ascending"]
    integral_scaling = certificate["integral_map"]["output_scaling"]
    map_expressions = [sp.sstr(component) for component in data["factorized_integral"]]
    q_expression = sp.sstr(data["q_in_t"])
    expected_polynomial = sp.sstr(data["polynomial"].as_expr())
    normalized_target = certificate["normalized_map"]["target"]
    integral_target = [rational_text(value) for value in data["integral_target"]]
    integral_jacobian = certificate["integral_map"]["jacobian_determinant"]
    return f'''#!/usr/bin/env python3
"""Generated SymPy input for the common-arithmetic-fibers explicit example."""
import sympy as sp

S, x, y, z = sp.symbols("S x y z")
factor_coefficients_ascending = {factor_coefficients!r}
factors = [
    sum(sp.Rational(c) * S**i for i, c in enumerate(coefficients))
    for coefficients in factor_coefficients_ascending
]
P5 = sp.expand(sp.prod(factors))
G5 = sp.expand(P5 - P5.subs(S, 0))
t = 1 + x*y
q = {q_expression}
F = (
    {map_expressions[0]},
    {map_expressions[1]},
    {map_expressions[2]},
)
normalized_target = tuple(map(sp.Rational, {normalized_target!r}))
integral_scaling = tuple(map(sp.Rational, {integral_scaling!r}))
integral_target = tuple(a*b for a, b in zip(integral_scaling, normalized_target))
inverse_polynomial = sp.expand(
    G5 - sp.Poly(G5, S).coeff_monomial(S)
    * (normalized_target[1]*S**2 + normalized_target[2]) / 2
)
assert P5 == {expected_polynomial}
assert inverse_polynomial == P5
assert integral_target == tuple(map(sp.Rational, {integral_target!r}))
assert sp.factor(sp.Matrix(F).jacobian((x, y, z)).det()) == sp.Rational("{integral_jacobian}")
print("PASS: generated SymPy paper example")
'''


def render_sage(data: Mapping[str, Any]) -> str:
    certificate = data["certificate"]
    factor_coefficients = certificate["factor_coefficients_ascending"]
    map_expressions = [
        sp.sstr(component).replace("**", "^")
        for component in data["factorized_integral"]
    ]
    q_expression = sp.sstr(data["q_in_t"]).replace("**", "^")
    expected_polynomial = sp.sstr(data["polynomial"].as_expr()).replace("**", "^")
    normalized_target = data["certificate"]["normalized_map"]["target"]
    integral_target = [rational_text(value) for value in data["integral_target"]]
    integral_jacobian = data["certificate"]["integral_map"]["jacobian_determinant"]
    return f"""\
# Generated Sage input for the common-arithmetic-fibers explicit example.
Q.<S> = PolynomialRing(QQ)
R.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
factor_coefficients_ascending = {factor_coefficients!r}
factors = [
    sum(QQ(c) * S^i for i, c in enumerate(coefficients))
    for coefficients in factor_coefficients_ascending
]
P5 = prod(factors)
G5 = P5 - P5(0)
t = 1 + x*y
q = {q_expression}
F = vector(R, [
    {map_expressions[0]},
    {map_expressions[1]},
    {map_expressions[2]},
])
normalized_target = vector(QQ, {normalized_target!r})
integral_scaling = diagonal_matrix(QQ, [1, 19, 19])
integral_target = integral_scaling * normalized_target
inverse_polynomial = G5 - G5[1] * (
    normalized_target[1]*S^2 + normalized_target[2]
) / 2
assert P5 == {expected_polynomial}
assert inverse_polynomial == P5
assert integral_target == vector(QQ, {integral_target!r})
assert det(jacobian(F, (x, y, z))) == QQ("{integral_jacobian}")
print("PASS: generated Sage paper example")
"""


def output_texts(data: Mapping[str, Any]) -> dict[Path, str]:
    outputs = data["specification"]["generated_outputs"]
    lean = render_lean(data)
    sage = render_sage(data)
    sympy = render_sympy(data)
    tex = render_tex(data)
    certificate = dict(data["certificate"])
    certificate["generated_views"] = {
        "lean": {
            "path": outputs["lean"],
            "sha256": hashlib.sha256(lean.encode("utf-8")).hexdigest(),
        },
        "sage": {
            "path": outputs["sage"],
            "sha256": hashlib.sha256(sage.encode("utf-8")).hexdigest(),
        },
        "sympy": {
            "path": outputs["sympy"],
            "sha256": hashlib.sha256(sympy.encode("utf-8")).hexdigest(),
        },
        "tex": {
            "path": outputs["tex"],
            "sha256": hashlib.sha256(tex.encode("utf-8")).hexdigest(),
        },
    }
    return {
        ROOT / outputs["certificate"]: (
            json.dumps(certificate, indent=2, sort_keys=True) + "\n"
        ),
        ROOT / outputs["lean"]: lean,
        ROOT / outputs["sage"]: sage,
        ROOT / outputs["sympy"]: sympy,
        ROOT / outputs["tex"]: tex,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if any generated output is absent or stale",
    )
    arguments = parser.parse_args()
    specification = json.loads(arguments.spec.read_text(encoding="utf-8"))
    data = build(specification)
    stale: list[str] = []
    for path, expected in output_texts(data).items():
        if arguments.check:
            if not path.exists() or path.read_text(encoding="utf-8") != expected:
                stale.append(str(path.relative_to(ROOT)))
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(expected, encoding="utf-8")
            print(f"WROTE: {path.relative_to(ROOT)}")
    if stale:
        raise SystemExit("stale generated explicit-example outputs:\n" + "\n".join(stale))
    if arguments.check:
        print("PASS: generated Lean/TeX/SymPy/Sage/JSON views are current")


if __name__ == "__main__":
    main()
