#!/usr/bin/env python3
"""Continue the exact Q-border component in its degree-five coefficient field.

The exact Q-border projection is

    R_20(s1,t1,t2,u) = 0,       s3 = -B/A,

where ``R_20`` has degree five in ``t2``.  Rather than adjoining ``R_20`` as
one more polynomial equation, this driver makes a root ``a`` of ``R_20`` a
Singular coefficient-field element over ``QQ(s1,t1,u)``.  The remaining
calculation is then only in the two fibre variables ``s6,s5``.

This is a research driver.  Its output records exact generic-point arithmetic
on the irreducible residual component; exceptional coefficient denominators
must be extracted and treated separately before the component is closed.
"""

from __future__ import annotations

import argparse
from itertools import combinations, permutations
from fractions import Fraction
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from research_two_pair_sic_bidegree33_t0_stratum_leading import ROOT
from verify_two_pair_sic_bidegree33_boundary_generic_quotient import (
    substitute,
)
from verify_two_pair_sic_bidegree33_corrected_boundary import (
    t0_open_localized_export,
)


RESULTANT_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_t0_stratum_Q_border_resultant_exact.json"
)
LEADING_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_t0_stratum_Q_leading_exact.json"
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=int, default=20)
    parser.add_argument(
        "--through",
        type=int,
        choices=range(5, 9),
        default=5,
        help="last moment to declare and, when at least six, adjoin",
    )
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument(
        "--specialize",
        action="append",
        default=[],
        metavar="VARIABLE=VALUE",
        help="fix any of s1,ell,u to a rational value before field arithmetic",
    )
    parser.add_argument(
        "--engine",
        choices=("auto", "minpoly", "extension", "python", "kummer"),
        default="auto",
        help=(
            "number-field minpoly for a closed fiber, or a polynomial "
            "degree-five extension presentation over remaining parameters"
        ),
    )
    parser.add_argument(
        "--python-stage",
        choices=(
            "pivot",
            "parse",
            "basis",
            "monic",
            "mu6parse",
            "mu6nf",
            "matrix",
            "mu6",
            "norm",
            "mu7nf",
            "matrices7",
            "matrices8",
            "pencil",
            "pencil_factors",
        ),
        default="norm",
    )
    parser.add_argument(
        "--pencil-factor-degree",
        type=int,
        choices=(3, 100),
        default=None,
    )
    parser.add_argument("--resume-checkpoint", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def residual_data() -> tuple[str, str, str]:
    payload = json.loads(RESULTANT_ARTIFACT.read_text(encoding="utf-8"))
    assert payload["prime"] == 0 and payload["stratum"] == "Q"
    residuals = [
        factor
        for factor in payload["factors"]
        if factor["term_count"] == 200
        and factor["total_degree"] == 20
        and factor["multiplicity"] == 2
    ]
    assert len(residuals) == 1
    pivot = payload["linear_subresultant"]
    assert pivot["gcd_residual_factor_A"] == "1"
    assert pivot["gcd_residual_factor_B"] == "1"
    return residuals[0]["factor"], pivot["A"], pivot["B"]


def adapted_polynomial(expression: str) -> tuple[str, int, int, sp.Rational]:
    """Use weight-zero ratios and remove the common ``u`` power/content."""

    s1, t1, t2, u, ell, T = sp.symbols("s1 t1 t2 u ell T")
    parsed = sp.sympify(
        expression.replace("^", "**"),
        locals={"s1": s1, "t1": t1, "t2": t2, "u": u},
    )
    polynomial = sp.Poly(
        sp.expand(
            parsed.subs(
                {
                    t1: u * (s1 - ell),
                    t2: T * u**2,
                },
                simultaneous=True,
            )
        ),
        s1,
        ell,
        T,
        u,
        domain=sp.QQ,
    )
    u_power = min(monomial[3] for monomial, _ in polynomial.terms())
    if u_power:
        polynomial = sp.Poly(
            polynomial.as_expr() / u**u_power,
            s1,
            ell,
            T,
            u,
            domain=sp.QQ,
        )
    content, primitive = polynomial.primitive()
    return (
        sp.sstr(primitive.as_expr()).replace("**", "^"),
        len(primitive.terms()),
        u_power,
        content,
    )


def parameter_to_root(expression: str) -> str:
    return re.sub(r"\bt2\b", "a", expression)


def parse_specializations(items: list[str]) -> dict[str, Fraction]:
    answer: dict[str, Fraction] = {}
    for item in items:
        match = re.fullmatch(r"(s1|ell|u)=(-?\d+(?:/\d+)?)", item)
        if match is None:
            raise ValueError(f"invalid specialization: {item}")
        variable, value = match.groups()
        answer[variable] = Fraction(value)
    return answer


def specialized_polynomial(
    expression: str,
    assignments: dict[str, Fraction],
) -> tuple[str, int, sp.Rational]:
    symbols = {name: sp.symbols(name) for name in ("s1", "ell", "T", "u")}
    parsed = sp.sympify(expression.replace("^", "**"), locals=symbols)
    result = sp.expand(
        parsed.subs(
            {
                symbols[name]: sp.Rational(value.numerator, value.denominator)
                for name, value in assignments.items()
            }
        )
    )
    remaining = [
        symbols[name]
        for name in ("s1", "ell", "T", "u")
        if name not in assignments
    ]
    content, primitive = sp.Poly(
        result,
        *remaining,
        domain=sp.QQ,
    ).primitive()
    return (
        sp.sstr(primitive.as_expr()).replace("**", "^"),
        len(primitive.terms()),
        content,
    )


def singular_fraction(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"({value.numerator}/{value.denominator})"
    )


def singular_bezout_inverse(
    singular: str,
    coefficient_parameters: list[str],
    modulus: str,
    value: str,
    timeout: int,
) -> str:
    """Invert one univariate polynomial modulo another over a function field."""

    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring inverseRing=(0,{",".join(coefficient_parameters)}),(a),dp;
option(redSB);
poly modulus={modulus};
poly value={value};
ideal inputIdeal=modulus,value;
matrix transformation;
ideal G=liftstd(inputIdeal,transformation);
poly inverse=transformation[2,1]/G[1];
ideal modulusBasis=std(ideal(modulus));
print(
  "META "+string(size(G))+" "
  +string(reduce(value*inverse-1,modulusBasis)==0)
);
print("INVERSE "+string(reduce(inverse,modulusBasis)));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=min(timeout, 10),
    )
    assert "\n   ? " not in completed.stdout, (
        completed.stdout[-4000:],
        completed.stderr[-2000:],
    )
    meta = re.search(r"(?m)^META (\d+) ([01])$", completed.stdout)
    inverse = re.search(r"(?m)^INVERSE (.*)$", completed.stdout)
    assert (
        meta is not None
        and meta.groups() == ("1", "1")
        and inverse is not None
    ), completed.stdout[-4000:]
    return inverse.group(1)


def singular_simplify_polynomial(
    singular: str,
    expression: str,
    timeout: int,
) -> str:
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring simplifyRing=0,(s6,s5,s3,a,u),dp;
poly value={expression};
print("VALUE "+string(value));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=min(timeout, 10),
    )
    assert "\n   ? " not in completed.stdout, (
        completed.stdout[-4000:],
        completed.stderr[-2000:],
    )
    value = re.search(r"(?m)^VALUE (.*)$", completed.stdout)
    assert value is not None
    return value.group(1)


class ExtensionArithmetic:
    """Arithmetic in K[a]/(modulus) with a fixed monic modulus."""

    def __init__(self, base, root: sp.Symbol, modulus: sp.Expr):
        self.base = base
        self.root = root
        self.modulus = sp.Poly(modulus, root, domain=base).monic()
        self.degree = self.modulus.degree()
        self.zero = sp.Poly(0, root, domain=base)
        self.one = sp.Poly(1, root, domain=base)
        self.generator = sp.Poly(root, root, domain=base)

    def constant(self, value) -> sp.Poly:
        return sp.Poly.from_dict(
            {(0,): self.base.convert(value)},
            self.root,
            domain=self.base,
        )

    def make(self, value) -> sp.Poly:
        if isinstance(value, sp.Poly):
            polynomial = sp.Poly(value.as_expr(), self.root, domain=self.base)
        else:
            polynomial = sp.Poly(value, self.root, domain=self.base)
        return polynomial.rem(self.modulus)

    def add(self, left: sp.Poly, right: sp.Poly) -> sp.Poly:
        return (left + right).rem(self.modulus)

    def neg(self, value: sp.Poly) -> sp.Poly:
        return -value

    def mul(self, left: sp.Poly, right: sp.Poly) -> sp.Poly:
        return (left * right).rem(self.modulus)

    def inverse(self, value: sp.Poly) -> sp.Poly:
        assert not value.is_zero
        if value == self.one:
            return self.one
        if value == -self.one:
            return -self.one
        if value.degree() == 0:
            coefficient = self.base.convert(value.nth(0))
            return self.constant(self.base.one / coefficient)
        matrix = [
            [self.base.zero for _ in range(self.degree)]
            for _ in range(self.degree)
        ]
        for column in range(self.degree):
            product = self.mul(
                value,
                self.power(self.generator, column),
            )
            for row in range(self.degree):
                matrix[row][column] = self.base.convert(product.nth(row))
        rhs = [
            [self.base.one if row == 0 else self.base.zero]
            for row in range(self.degree)
        ]
        matrix_domain = DomainMatrix(
            matrix,
            (self.degree, self.degree),
            self.base,
        )
        rhs_domain = DomainMatrix(
            rhs,
            (self.degree, 1),
            self.base,
        )
        solution_numerator, solution_denominator = (
            matrix_domain.solve_den(rhs_domain, method="rref")
        )
        solution = solution_numerator.to_list()
        return sp.Poly.from_dict(
            {
                (index,): solution[index][0] / solution_denominator
                for index in range(self.degree)
                if solution[index][0]
            },
            self.root,
            domain=self.base,
        )

    def power(self, value: sp.Poly, exponent: int) -> sp.Poly:
        answer = self.one
        factor = value
        power = exponent
        while power:
            if power & 1:
                answer = self.mul(answer, factor)
            factor = self.mul(factor, factor)
            power //= 2
        return answer


BiPolynomial = dict[tuple[int, int], sp.Poly]


def monomial_key(monomial: tuple[int, int]) -> tuple[int, int]:
    """Degree-reverse-lex key for s6>s5."""

    return sum(monomial), -monomial[1]


def clean_bipolynomial(polynomial: BiPolynomial) -> BiPolynomial:
    return {
        monomial: coefficient
        for monomial, coefficient in polynomial.items()
        if not coefficient.is_zero
    }


def leading_term(
    polynomial: BiPolynomial,
) -> tuple[tuple[int, int], sp.Poly]:
    monomial = max(polynomial, key=monomial_key)
    return monomial, polynomial[monomial]


def shifted_scaled(
    arithmetic: ExtensionArithmetic,
    polynomial: BiPolynomial,
    shift: tuple[int, int],
    scalar: sp.Poly,
) -> BiPolynomial:
    return clean_bipolynomial(
        {
            (monomial[0] + shift[0], monomial[1] + shift[1]): (
                arithmetic.mul(coefficient, scalar)
            )
            for monomial, coefficient in polynomial.items()
        }
    )


def add_bipolynomials(
    arithmetic: ExtensionArithmetic,
    left: BiPolynomial,
    right: BiPolynomial,
    right_sign: int = 1,
) -> BiPolynomial:
    answer = dict(left)
    for monomial, coefficient in right.items():
        if right_sign < 0:
            coefficient = arithmetic.neg(coefficient)
        answer[monomial] = arithmetic.add(
            answer.get(monomial, arithmetic.zero),
            coefficient,
        )
    return clean_bipolynomial(answer)


def monic_bipolynomial(
    arithmetic: ExtensionArithmetic,
    polynomial: BiPolynomial,
) -> BiPolynomial:
    _, coefficient = leading_term(polynomial)
    if coefficient == arithmetic.one:
        return polynomial
    return shifted_scaled(
        arithmetic,
        polynomial,
        (0, 0),
        arithmetic.inverse(coefficient),
    )


def normal_form_bipolynomial(
    arithmetic: ExtensionArithmetic,
    polynomial: BiPolynomial,
    basis: list[BiPolynomial],
) -> BiPolynomial:
    remainder: BiPolynomial = {}
    work = dict(polynomial)
    while work:
        monomial, coefficient = leading_term(work)
        reducer = None
        for candidate in basis:
            leading_monomial, _ = leading_term(candidate)
            if all(
                monomial[index] >= leading_monomial[index]
                for index in range(2)
            ):
                reducer = candidate
                break
        if reducer is None:
            remainder[monomial] = arithmetic.add(
                remainder.get(monomial, arithmetic.zero),
                coefficient,
            )
            del work[monomial]
            continue
        leading_monomial, leading_coefficient = leading_term(reducer)
        scalar = arithmetic.mul(
            coefficient,
            (
                arithmetic.one
                if leading_coefficient == arithmetic.one
                else arithmetic.inverse(leading_coefficient)
            ),
        )
        multiple = shifted_scaled(
            arithmetic,
            reducer,
            (
                monomial[0] - leading_monomial[0],
                monomial[1] - leading_monomial[1],
            ),
            scalar,
        )
        work = add_bipolynomials(arithmetic, work, multiple, -1)
    return clean_bipolynomial(remainder)


def buchberger_basis(
    arithmetic: ExtensionArithmetic,
    generators: list[BiPolynomial],
) -> list[BiPolynomial]:
    basis = [
        monic_bipolynomial(arithmetic, generator)
        for generator in generators
        if generator
    ]
    pairs = [
        (left, right)
        for left in range(len(basis))
        for right in range(left)
    ]
    while pairs:
        left_index, right_index = pairs.pop(0)
        left = basis[left_index]
        right = basis[right_index]
        left_monomial, _ = leading_term(left)
        right_monomial, _ = leading_term(right)
        common = (
            max(left_monomial[0], right_monomial[0]),
            max(left_monomial[1], right_monomial[1]),
        )
        s_polynomial = add_bipolynomials(
            arithmetic,
            shifted_scaled(
                arithmetic,
                left,
                (
                    common[0] - left_monomial[0],
                    common[1] - left_monomial[1],
                ),
                arithmetic.one,
            ),
            shifted_scaled(
                arithmetic,
                right,
                (
                    common[0] - right_monomial[0],
                    common[1] - right_monomial[1],
                ),
                arithmetic.one,
            ),
            -1,
        )
        remainder = normal_form_bipolynomial(
            arithmetic,
            s_polynomial,
            basis,
        )
        if remainder:
            remainder = monic_bipolynomial(arithmetic, remainder)
            new_index = len(basis)
            if new_index >= 12:
                raise RuntimeError(
                    "unexpected border-basis growth: "
                    + str(
                        [
                            leading_term(polynomial)[0]
                            for polynomial in basis
                        ]
                    )
                )
            pairs.extend((new_index, old) for old in range(new_index))
            basis.append(remainder)
    return basis


def pseudo_normal_form_bipolynomial(
    arithmetic: ExtensionArithmetic,
    polynomial: BiPolynomial,
    basis: list[BiPolynomial],
) -> BiPolynomial:
    """Fraction-free normal form, sufficient for Buchberger zero tests."""

    remainder: BiPolynomial = {}
    work = dict(polynomial)
    while work:
        monomial, coefficient = leading_term(work)
        reducer = None
        for candidate in basis:
            leading_monomial, _ = leading_term(candidate)
            if all(
                monomial[index] >= leading_monomial[index]
                for index in range(2)
            ):
                reducer = candidate
                break
        if reducer is None:
            remainder[monomial] = arithmetic.add(
                remainder.get(monomial, arithmetic.zero),
                coefficient,
            )
            del work[monomial]
            continue
        leading_monomial, leading_coefficient = leading_term(reducer)
        work = shifted_scaled(
            arithmetic,
            work,
            (0, 0),
            leading_coefficient,
        )
        remainder = shifted_scaled(
            arithmetic,
            remainder,
            (0, 0),
            leading_coefficient,
        )
        multiple = shifted_scaled(
            arithmetic,
            reducer,
            (
                monomial[0] - leading_monomial[0],
                monomial[1] - leading_monomial[1],
            ),
            coefficient,
        )
        work = add_bipolynomials(arithmetic, work, multiple, -1)
    return clean_bipolynomial(remainder)


def fraction_free_buchberger_basis(
    arithmetic: ExtensionArithmetic,
    generators: list[BiPolynomial],
) -> list[BiPolynomial]:
    basis = [generator for generator in generators if generator]
    pairs = [
        (left, right)
        for left in range(len(basis))
        for right in range(left)
    ]
    while pairs:
        left_index, right_index = pairs.pop(0)
        left = basis[left_index]
        right = basis[right_index]
        left_monomial, left_coefficient = leading_term(left)
        right_monomial, right_coefficient = leading_term(right)
        common = (
            max(left_monomial[0], right_monomial[0]),
            max(left_monomial[1], right_monomial[1]),
        )
        s_polynomial = add_bipolynomials(
            arithmetic,
            shifted_scaled(
                arithmetic,
                left,
                (
                    common[0] - left_monomial[0],
                    common[1] - left_monomial[1],
                ),
                right_coefficient,
            ),
            shifted_scaled(
                arithmetic,
                right,
                (
                    common[0] - right_monomial[0],
                    common[1] - right_monomial[1],
                ),
                left_coefficient,
            ),
            -1,
        )
        remainder = pseudo_normal_form_bipolynomial(
            arithmetic,
            s_polynomial,
            basis,
        )
        if remainder:
            new_index = len(basis)
            if new_index >= 12:
                raise RuntimeError(
                    "unexpected fraction-free border-basis growth: "
                    + str(
                        [
                            leading_term(polynomial)[0]
                            for polynomial in basis
                        ]
                    )
                )
            pairs.extend((new_index, old) for old in range(new_index))
            basis.append(remainder)
    return basis


def parse_extension_bipolynomial(
    expression: str,
    arithmetic: ExtensionArithmetic,
    base_symbols: list[sp.Symbol],
    pivot: sp.Poly,
    extension_symbol: sp.Symbol,
) -> BiPolynomial:
    s6, s5, s3 = sp.symbols("s6 s5 s3")
    locals_map = {
        "s6": s6,
        "s5": s5,
        "s3": s3,
        str(extension_symbol): extension_symbol,
        **{str(symbol): symbol for symbol in base_symbols},
    }
    parsed = sp.sympify(expression.replace("^", "**"), locals=locals_map)
    polynomial = sp.Poly(
        parsed,
        s6,
        s5,
        s3,
        extension_symbol,
        domain=arithmetic.base,
    )
    pivot_powers = {
        exponent: arithmetic.power(pivot, exponent)
        for exponent in range(polynomial.degree(s3) + 1)
    }
    root_powers = {
        exponent: arithmetic.power(arithmetic.generator, exponent)
        for exponent in range(polynomial.degree(extension_symbol) + 1)
    }
    answer: BiPolynomial = {}
    for (s6_power, s5_power, s3_power, a_power), coefficient in (
        polynomial.terms()
    ):
        value = arithmetic.constant(coefficient)
        value = arithmetic.mul(value, pivot_powers[s3_power])
        value = arithmetic.mul(value, root_powers[a_power])
        monomial = (s6_power, s5_power)
        answer[monomial] = arithmetic.add(
            answer.get(monomial, arithmetic.zero),
            value,
        )
    return clean_bipolynomial(answer)


def standard_monomials(
    leading_monomials: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    x_bound = min(
        monomial[0] for monomial in leading_monomials if monomial[1] == 0
    )
    y_bound = min(
        monomial[1] for monomial in leading_monomials if monomial[0] == 0
    )
    answer = [
        (x_power, y_power)
        for x_power in range(x_bound)
        for y_power in range(y_bound)
        if not any(
            x_power >= leading[0] and y_power >= leading[1]
            for leading in leading_monomials
        )
    ]
    return sorted(answer, key=monomial_key)


def determinant_over_extension(
    arithmetic: ExtensionArithmetic,
    matrix: list[list[sp.Poly]],
) -> sp.Poly:
    size = len(matrix)
    answer = arithmetic.zero
    for permutation in permutations(range(size)):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(size)
            for right in range(left + 1, size)
        )
        term = arithmetic.one
        for row, column in enumerate(permutation):
            term = arithmetic.mul(term, matrix[row][column])
        answer = arithmetic.add(
            answer,
            term if inversions % 2 == 0 else arithmetic.neg(term),
        )
    return answer


def multiplication_matrix_from_normal_form(
    arithmetic: ExtensionArithmetic,
    basis: list[BiPolynomial],
    monomials: list[tuple[int, int]],
    normal_form: BiPolynomial,
) -> list[list[sp.Poly]]:
    matrix = [
        [arithmetic.zero for _ in monomials] for _ in monomials
    ]
    monomial_index = {
        monomial: index for index, monomial in enumerate(monomials)
    }
    monomial_normal_forms: dict[tuple[int, int], BiPolynomial] = {
        monomial: {monomial: arithmetic.one}
        for monomial in monomials
    }
    for column, monomial in enumerate(monomials):
        reduced: BiPolynomial = {}
        for multiplier, coefficient in normal_form.items():
            product_monomial = (
                monomial[0] + multiplier[0],
                monomial[1] + multiplier[1],
            )
            if product_monomial not in monomial_normal_forms:
                monomial_normal_forms[product_monomial] = (
                    normal_form_bipolynomial(
                        arithmetic,
                        {product_monomial: arithmetic.one},
                        basis,
                    )
                )
            reduced = add_bipolynomials(
                arithmetic,
                reduced,
                shifted_scaled(
                    arithmetic,
                    monomial_normal_forms[product_monomial],
                    (0, 0),
                    coefficient,
                ),
            )
        assert set(reduced) <= set(monomials)
        for reduced_monomial, coefficient in reduced.items():
            matrix[monomial_index[reduced_monomial]][column] = coefficient
    return matrix


def determinant_pencil_linear_coefficient(
    arithmetic: ExtensionArithmetic,
    left: list[list[sp.Poly]],
    right: list[list[sp.Poly]],
) -> sp.Poly:
    size = len(left)
    answer = arithmetic.zero
    for permutation in permutations(range(size)):
        inversions = sum(
            permutation[first] > permutation[second]
            for first in range(size)
            for second in range(first + 1, size)
        )
        signed_terms = arithmetic.zero
        for selected_row in range(size):
            term = arithmetic.one
            for row, column in enumerate(permutation):
                term = arithmetic.mul(
                    term,
                    (
                        right[row][column]
                        if row == selected_row
                        else left[row][column]
                    ),
                )
            signed_terms = arithmetic.add(signed_terms, term)
        answer = arithmetic.add(
            answer,
            (
                signed_terms
                if inversions % 2 == 0
                else arithmetic.neg(signed_terms)
            ),
        )
    return answer


def determinant_over_base(matrix, base):
    size = len(matrix)
    answer = base.zero
    for permutation in permutations(range(size)):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(size)
            for right in range(left + 1, size)
        )
        term = base.one
        for row, column in enumerate(permutation):
            term *= matrix[row][column]
        answer += term if inversions % 2 == 0 else -term
    return answer


def evaluate_fraction_field_element_at_algebraic_root(
    value,
    fraction_field,
    number_field,
):
    value = fraction_field.convert(value)
    generator = number_field.unit

    def evaluate_polynomial(polynomial):
        answer = number_field.zero
        for monomial, coefficient in polynomial.terms():
            answer += (
                number_field.convert(coefficient)
                * generator ** monomial[0]
            )
        return answer

    return evaluate_polynomial(value.numer) / evaluate_polynomial(value.denom)


def specialize_extension_element_to_number_field(
    value: sp.Poly,
    source: ExtensionArithmetic,
    target: ExtensionArithmetic,
) -> sp.Poly:
    return sp.Poly.from_dict(
        {
            (root_power,): evaluate_fraction_field_element_at_algebraic_root(
                source.base.convert(value.nth(root_power)),
                source.base,
                target.base,
            )
            for root_power in range(source.degree)
            if value.nth(root_power) != 0
        },
        target.root,
        domain=target.base,
    ).rem(target.modulus)


def pencil_factor_tests(
    arithmetic: ExtensionArithmetic,
    left: list[list[sp.Poly]],
    right: list[list[sp.Poly]],
    factors: list[dict[str, object]],
    parameter: sp.Symbol,
    third: list[list[sp.Poly]] | None = None,
) -> list[dict[str, object]]:
    tests = []
    for factor_record in factors:
        factor_expression = sp.sympify(
            str(factor_record["factor"]),
            locals={str(parameter): parameter},
        )
        factor_polynomial = sp.Poly(
            factor_expression,
            parameter,
            domain=sp.QQ,
        )
        number_field = sp.QQ.alg_field_from_poly(
            factor_polynomial,
            alias=f"theta_{factor_polynomial.degree()}",
        )
        root = arithmetic.root
        specialized_modulus = sp.Poly.from_dict(
            {
                (root_power,): (
                    evaluate_fraction_field_element_at_algebraic_root(
                        arithmetic.base.convert(
                            arithmetic.modulus.nth(root_power)
                        ),
                        arithmetic.base,
                        number_field,
                    )
                )
                for root_power in range(arithmetic.degree + 1)
                if arithmetic.modulus.nth(root_power) != 0
            },
            root,
            domain=number_field,
        )
        specialized_arithmetic = ExtensionArithmetic(
            number_field,
            root,
            specialized_modulus.as_expr(),
        )
        specialized_left = [
            [
                specialize_extension_element_to_number_field(
                    entry,
                    arithmetic,
                    specialized_arithmetic,
                )
                for entry in row
            ]
            for row in left
        ]
        specialized_right = [
            [
                specialize_extension_element_to_number_field(
                    entry,
                    arithmetic,
                    specialized_arithmetic,
                )
                for entry in row
            ]
            for row in right
        ]
        specialized_third = (
            None
            if third is None
            else [
                [
                    specialize_extension_element_to_number_field(
                        entry,
                        arithmetic,
                        specialized_arithmetic,
                    )
                    for entry in row
                ]
                for row in third
            ]
        )
        block_columns = [
            [specialized_left[row][column] for row in range(len(left))]
            for column in range(len(left))
        ] + [
            [specialized_right[row][column] for row in range(len(right))]
            for column in range(len(right))
        ]
        nonzero_minor_columns = None
        for selected_columns in combinations(
            range(len(block_columns)),
            len(left),
        ):
            minor = [
                [
                    block_columns[column][row]
                    for column in selected_columns
                ]
                for row in range(len(left))
            ]
            if not determinant_over_extension(
                specialized_arithmetic,
                minor,
            ).is_zero:
                nonzero_minor_columns = list(selected_columns)
                break
        nonzero_mu8_minor_columns = None
        if (
            nonzero_minor_columns is None
            and specialized_third is not None
        ):
            extended_columns = block_columns + [
                [
                    specialized_third[row][column]
                    for row in range(len(third))
                ]
                for column in range(len(third))
            ]
            for selected_columns in combinations(
                range(len(extended_columns)),
                len(left),
            ):
                if max(selected_columns) < len(block_columns):
                    continue
                minor = [
                    [
                        extended_columns[column][row]
                        for column in selected_columns
                    ]
                    for row in range(len(left))
                ]
                if not determinant_over_extension(
                    specialized_arithmetic,
                    minor,
                ).is_zero:
                    nonzero_mu8_minor_columns = list(selected_columns)
                    break
        coefficient = (
            None
            if (
                nonzero_minor_columns is not None
                or nonzero_mu8_minor_columns is not None
            )
            else determinant_pencil_linear_coefficient(
                specialized_arithmetic,
                specialized_left,
                specialized_right,
            )
        )
        tests.append(
            {
                "factor": str(factor_record["factor"]),
                "degree": factor_polynomial.degree(),
                "mu6_norm_multiplicity": factor_record["multiplicity"],
                "joint_rank_full": nonzero_minor_columns is not None,
                "nonzero_block_minor_columns_zero_based": (
                    nonzero_minor_columns
                ),
                "joint_rank_with_mu8_full": (
                    nonzero_minor_columns is not None
                    or nonzero_mu8_minor_columns is not None
                ),
                "nonzero_mu8_block_minor_columns_zero_based": (
                    nonzero_mu8_minor_columns
                ),
                "c1_zero_in_factor_field": (
                    None if coefficient is None else coefficient.is_zero
                ),
                "excluded_by_c1": (
                    None if coefficient is None else not coefficient.is_zero
                ),
                "excluded_through_mu7": (
                    nonzero_minor_columns is not None
                    or (coefficient is not None and not coefficient.is_zero)
                ),
                "excluded_through_mu8": (
                    nonzero_minor_columns is not None
                    or nonzero_mu8_minor_columns is not None
                    or (
                        coefficient is not None
                        and not coefficient.is_zero
                    )
                ),
                "c1_extension_root_degree": (
                    None
                    if coefficient is None
                    else (-1 if coefficient.is_zero else coefficient.degree())
                ),
                "c1_extension_root_terms": (
                    None
                    if coefficient is None
                    else (0 if coefficient.is_zero else len(coefficient.terms()))
                ),
            }
        )
    return tests


def serialize_extension_matrix(
    matrix: list[list[sp.Poly]],
) -> list[list[str]]:
    return [
        [sp.sstr(entry.as_expr()) for entry in row]
        for row in matrix
    ]


def serialize_bipolynomial(polynomial: BiPolynomial) -> dict[str, str]:
    return {
        f"{monomial[0]},{monomial[1]}": sp.sstr(coefficient.as_expr())
        for monomial, coefficient in sorted(polynomial.items())
    }


def resume_pencil_checkpoint(
    checkpoint_path: Path,
    factor_degree: int,
) -> dict[str, object]:
    checkpoint_payload = json.loads(
        checkpoint_path.read_text(encoding="utf-8")
    )
    checkpoint = checkpoint_payload["matrix_checkpoint"]
    parameter = sp.symbols(checkpoint["base_parameter"])
    root = sp.symbols(checkpoint["extension_symbol"])
    base = sp.QQ.frac_field(parameter)
    locals_map = {
        str(parameter): parameter,
        str(root): root,
    }
    arithmetic = ExtensionArithmetic(
        base,
        root,
        sp.sympify(checkpoint["modulus"], locals=locals_map),
    )

    def parse_matrix(serialized):
        return [
            [
                arithmetic.make(sp.sympify(entry, locals=locals_map))
                for entry in row
            ]
            for row in serialized
        ]

    left = parse_matrix(checkpoint["M6"])
    right = parse_matrix(checkpoint["M7"])
    third = (
        parse_matrix(checkpoint["M8"])
        if "M8" in checkpoint
        else None
    )
    norm_artifact = (
        ROOT
        / "artifacts"
        / "generated-results"
        / (
            "two_pair_sic_bidegree33_t0_stratum_Q_"
            "residual_slice_s1_0_ell_0_exact.json"
        )
    )
    norm_payload = json.loads(norm_artifact.read_text(encoding="utf-8"))
    factors = []
    for factor in norm_payload["mu6_norm"]["numerator_factors"]:
        degree = sp.Poly(
            sp.sympify(
                factor["factor"],
                locals={str(parameter): parameter},
            ),
            parameter,
            domain=sp.QQ,
        ).degree()
        if degree == factor_degree:
            factors.append(factor)
    assert len(factors) == 1
    tests = pencil_factor_tests(
        arithmetic,
        left,
        right,
        factors,
        parameter,
        third,
    )
    return {
        "format": "two-pair-sic-bidegree33-t0-Q-residual-pencil-resume-v1",
        "status": (
            "exact characteristic-zero pencil-coefficient test resumed "
            "from exact multiplication matrices"
        ),
        "checkpoint": str(checkpoint_path),
        "pencil_factor_degree": factor_degree,
        "mu6_mu7_factor_tests": tests,
        "all_tested_components_excluded_through_mu7": all(
            test["excluded_through_mu7"] for test in tests
        ),
        "all_tested_components_excluded_through_mu8": all(
            test["excluded_through_mu8"] for test in tests
        ),
    }


def python_extension_certificate(
    residual: str,
    pivot_a: str,
    pivot_b: str,
    pivot_scale: str,
    moment_polynomials: dict[int, str],
    border_polynomial: str,
    remaining_parameters: list[str],
    stage: str,
    projection: str = "quintic",
    pencil_factor_degree: int | None = None,
) -> dict[str, object]:
    """Compute the residual slice by explicit degree-five arithmetic."""

    if projection == "quintic":
        if len(remaining_parameters) != 1:
            raise ValueError(
                "the bounded Python engine requires one base parameter"
            )
        base_symbols = [sp.symbols(name) for name in remaining_parameters]
        root = sp.symbols("a")
        locals_map = {
            "T": root,
            "a": root,
            **{str(symbol): symbol for symbol in base_symbols},
        }
        base = sp.QQ.frac_field(*base_symbols)
        modulus_expression = sp.sympify(
            residual.replace("^", "**"),
            locals=locals_map,
        ).subs(sp.symbols("T"), root)
        projected_moments = moment_polynomials
    else:
        assert projection == "kummer"
        assert remaining_parameters == ["u"]
        parameter = sp.symbols("T")
        root = sp.symbols("u")
        base_symbols = [parameter]
        base = sp.QQ.frac_field(parameter)
        locals_map = {"T": parameter, "u": root}
        residual_expression = sp.sympify(
            residual.replace("^", "**"),
            locals=locals_map,
        )
        residual_polynomial = sp.Poly(
            residual_expression,
            root,
            domain=base,
        )
        common_root_power = min(
            monomial[0] for monomial, _ in residual_polynomial.terms()
        )
        assert common_root_power > 0
        modulus_expression = sp.cancel(
            residual_expression / root**common_root_power
        )
        projected_moments = {
            order: re.sub(r"\ba\b", "T", polynomial)
            for order, polynomial in moment_polynomials.items()
        }
    arithmetic = ExtensionArithmetic(base, root, modulus_expression)
    pivot_a_expression = sp.sympify(
        pivot_a.replace("^", "**"),
        locals=locals_map,
    )
    pivot_b_expression = sp.sympify(
        pivot_b.replace("^", "**"),
        locals=locals_map,
    )
    if projection == "quintic":
        pivot_a_expression = pivot_a_expression.subs(
            sp.symbols("T"),
            root,
        )
        pivot_b_expression = pivot_b_expression.subs(
            sp.symbols("T"),
            root,
        )
    pivot_a_element = arithmetic.make(pivot_a_expression)
    pivot_b_element = arithmetic.make(pivot_b_expression)
    pivot_scale_expression = sp.sympify(
        pivot_scale.replace("^", "**"),
        locals=locals_map,
    )
    pivot = arithmetic.mul(
        arithmetic.constant(-base.from_sympy(pivot_scale_expression)),
        arithmetic.mul(
            pivot_b_element,
            arithmetic.inverse(pivot_a_element),
        ),
    )
    if stage == "pivot":
        return {
            "extension_degree": arithmetic.degree,
            "pivot_A_degree": pivot_a_element.degree(),
            "pivot_B_degree": pivot_b_element.degree(),
            "pivot_A": sp.sstr(pivot_a_element.as_expr()),
            "extension_modulus": sp.sstr(arithmetic.modulus.as_expr()),
            "pivot_degree_in_extension_root": pivot.degree(),
            "pivot_term_count": len(pivot.terms()),
        }
    parsed = {
        order: parse_extension_bipolynomial(
            polynomial,
            arithmetic,
            base_symbols,
            pivot,
            root,
        )
        for order, polynomial in projected_moments.items()
        if order in (3, 4, 5)
    }
    assert not parsed[3]
    if stage == "parse":
        input_leads = {}
        for order in (4, 5):
            monomial, coefficient = leading_term(parsed[order])
            input_leads[f"mu{order}"] = {
                "monomial": f"{monomial[0]},{monomial[1]}",
                "extension_root_degree": coefficient.degree(),
                "extension_root_terms": len(coefficient.terms()),
                "coefficient_lengths": [
                    len(str(coefficient.nth(index)))
                    for index in range(arithmetic.degree)
                ],
            }
        return {
            "extension_degree": arithmetic.degree,
            "mu3_zero_in_extension": True,
            "fiber_term_counts": {
                f"mu{order}": len(parsed[order])
                for order in (4, 5)
            },
            "input_leading_terms": input_leads,
        }
    basis = fraction_free_buchberger_basis(
        arithmetic,
        [parsed[4], parsed[5]],
    )
    minimal_basis: list[BiPolynomial] = []
    for polynomial in sorted(
        basis,
        key=lambda value: monomial_key(leading_term(value)[0]),
    ):
        monomial, _ = leading_term(polynomial)
        if any(
            monomial[0] >= retained_monomial[0]
            and monomial[1] >= retained_monomial[1]
            for retained_monomial, _ in map(
                leading_term,
                minimal_basis,
            )
        ):
            continue
        minimal_basis.append(polynomial)
    basis = minimal_basis
    leading_monomials = [
        leading_term(polynomial)[0] for polynomial in basis
    ]
    monomials = standard_monomials(leading_monomials)
    result: dict[str, object] = {
        "extension_degree": arithmetic.degree,
        "fiber_quotient_length": len(monomials),
        "leading_exponents": [
            f"{monomial[0]},{monomial[1]}"
            for monomial in leading_monomials
        ],
        "standard_monomials": [
            (
                "1"
                if monomial == (0, 0)
                else "*".join(
                    filter(
                        None,
                        (
                            (
                                "s6"
                                if monomial[0] == 1
                                else (
                                    f"s6^{monomial[0]}"
                                    if monomial[0]
                                    else ""
                                )
                            ),
                            (
                                "s5"
                                if monomial[1] == 1
                                else (
                                    f"s5^{monomial[1]}"
                                    if monomial[1]
                                    else ""
                                )
                            ),
                        ),
                    )
                )
            )
            for monomial in monomials
        ],
        "mu3_zero_in_extension": True,
        "border_zero_in_extension": (
            "follows from the certified resultant/pseudo-remainder identity; "
            "not replayed in the bounded vector stage"
        ),
    }
    if stage == "basis":
        return result
    basis = [
        monic_bipolynomial(arithmetic, polynomial)
        for polynomial in basis
    ]
    if stage == "monic":
        result["monic_basis_term_counts"] = [
            len(polynomial) for polynomial in basis
        ]
        return result
    parsed_mu6 = parse_extension_bipolynomial(
        projected_moments[6],
        arithmetic,
        base_symbols,
        pivot,
        root,
    )
    if stage == "mu6parse":
        result["mu6_fiber_term_count"] = len(parsed_mu6)
        return result
    normal_mu6 = normal_form_bipolynomial(
        arithmetic,
        parsed_mu6,
        basis,
    )
    if stage == "mu6nf":
        result["mu6_normal_form_support"] = len(normal_mu6)
        return result
    multiplication_matrix = multiplication_matrix_from_normal_form(
        arithmetic,
        basis,
        monomials,
        normal_mu6,
    )
    if stage == "matrix":
        result["mu6_normal_form_support"] = len(normal_mu6)
        result["mu6_multiplication_matrix_nonzero_entries"] = sum(
            not entry.is_zero
            for row in multiplication_matrix
            for entry in row
        )
        return result
    determinant = determinant_over_extension(
        arithmetic,
        multiplication_matrix,
    )
    assert not determinant.is_zero
    result.update(
        {
            "mu6_normal_form_support": len(normal_mu6),
            "mu6_is_unit_at_generic_slice_point": True,
            "mu6_multiplication_determinant_degree_in_extension_root": (
                determinant.degree()
            ),
        }
    )
    if stage == "mu6":
        return result
    if stage in (
        "mu7nf",
        "matrices7",
        "matrices8",
        "pencil",
        "pencil_factors",
    ):
        if 7 not in projected_moments:
            raise ValueError("the pencil stages require --through 7")
        parsed_mu7 = parse_extension_bipolynomial(
            projected_moments[7],
            arithmetic,
            base_symbols,
            pivot,
            root,
        )
        normal_mu7 = normal_form_bipolynomial(
            arithmetic,
            parsed_mu7,
            basis,
        )
        result["mu7_normal_form_support"] = len(normal_mu7)
        if stage == "mu7nf":
            return result
        multiplication_matrix_mu7 = multiplication_matrix_from_normal_form(
            arithmetic,
            basis,
            monomials,
            normal_mu7,
        )
        if stage in ("matrices7", "matrices8"):
            matrix_checkpoint = {
                "base_parameter": str(base_symbols[0]),
                "extension_symbol": str(root),
                "modulus": sp.sstr(arithmetic.modulus.as_expr()),
                "pivot_s3": sp.sstr(pivot.as_expr()),
                "monic_fiber_basis": [
                    serialize_bipolynomial(polynomial)
                    for polynomial in basis
                ],
                "M6": serialize_extension_matrix(multiplication_matrix),
                "M7": serialize_extension_matrix(multiplication_matrix_mu7),
            }
            if stage == "matrices8":
                if 8 not in projected_moments:
                    raise ValueError("matrices8 requires --through 8")
                parsed_mu8 = parse_extension_bipolynomial(
                    projected_moments[8],
                    arithmetic,
                    base_symbols,
                    pivot,
                    root,
                )
                normal_mu8 = normal_form_bipolynomial(
                    arithmetic,
                    parsed_mu8,
                    basis,
                )
                result["mu8_normal_form_support"] = len(normal_mu8)
                multiplication_matrix_mu8 = (
                    multiplication_matrix_from_normal_form(
                        arithmetic,
                        basis,
                        monomials,
                        normal_mu8,
                    )
                )
                matrix_checkpoint["M8"] = serialize_extension_matrix(
                    multiplication_matrix_mu8
                )
            result["matrix_checkpoint"] = matrix_checkpoint
            return result
        if stage == "pencil_factors":
            norm_artifact = (
                ROOT
                / "artifacts"
                / "generated-results"
                / (
                    "two_pair_sic_bidegree33_t0_stratum_Q_"
                    "residual_slice_s1_0_ell_0_exact.json"
                )
            )
            norm_payload = json.loads(
                norm_artifact.read_text(encoding="utf-8")
            )
            factors = norm_payload["mu6_norm"]["numerator_factors"]
            if pencil_factor_degree is not None:
                factors = [
                    factor
                    for factor in factors
                    if sp.Poly(
                        sp.sympify(
                            factor["factor"],
                            locals={"T": base_symbols[0]},
                        ),
                        base_symbols[0],
                        domain=sp.QQ,
                    ).degree()
                    == pencil_factor_degree
                ]
                assert len(factors) == 1
            result["mu6_mu7_factor_tests"] = pencil_factor_tests(
                arithmetic,
                multiplication_matrix,
                multiplication_matrix_mu7,
                factors,
                base_symbols[0],
            )
            result["all_mu6_norm_components_excluded_through_mu7"] = all(
                test["excluded_through_mu7"]
                for test in result["mu6_mu7_factor_tests"]
            )
            return result
        pencil_coefficient = determinant_pencil_linear_coefficient(
            arithmetic,
            multiplication_matrix,
            multiplication_matrix_mu7,
        )
        assert determinant.degree() == 0
        determinant_base = base.convert(determinant.nth(0))
        determinant_expression = sp.cancel(base.to_sympy(determinant_base))
        determinant_numerator = sp.Poly(
            sp.fraction(determinant_expression)[0],
            *base_symbols,
            domain=sp.QQ,
        ).primitive()[1]
        determinant_factors = sp.factor_list(
            determinant_numerator.as_expr()
        )[1]
        pencil_numerators: list[sp.Poly] = []
        pencil_coefficient_profiles = []
        for root_power in range(arithmetic.degree):
            coefficient = base.convert(pencil_coefficient.nth(root_power))
            coefficient_expression = sp.cancel(base.to_sympy(coefficient))
            numerator_expression = sp.fraction(coefficient_expression)[0]
            numerator = sp.Poly(
                numerator_expression,
                *base_symbols,
                domain=sp.QQ,
            )
            if not numerator.is_zero:
                numerator = numerator.primitive()[1]
                pencil_numerators.append(numerator)
            pencil_coefficient_profiles.append(
                {
                    "extension_root_power": root_power,
                    "zero": numerator.is_zero,
                    "numerator_degree": (
                        -1 if numerator.is_zero else numerator.total_degree()
                    ),
                    "numerator_terms": (
                        0 if numerator.is_zero else len(numerator.terms())
                    ),
                }
            )
        common_gcd = determinant_numerator
        for numerator in pencil_numerators:
            common_gcd = sp.gcd(common_gcd, numerator)
        factor_tests = []
        for factor, multiplicity in determinant_factors:
            factor_polynomial = sp.Poly(
                factor,
                *base_symbols,
                domain=sp.QQ,
            )
            divides_all = all(
                numerator.rem(factor_polynomial).is_zero
                for numerator in pencil_numerators
            )
            factor_tests.append(
                {
                    "factor": sp.sstr(factor),
                    "degree": factor_polynomial.total_degree(),
                    "mu6_determinant_multiplicity": multiplicity,
                    "divides_every_c1_coordinate": divides_all,
                    "excluded_by_c1": not divides_all,
                }
            )
        result["mu6_mu7_pencil"] = {
            "c1_extension_root_degree": pencil_coefficient.degree(),
            "c1_coordinate_profiles": pencil_coefficient_profiles,
            "common_gcd": sp.sstr(common_gcd.as_expr()),
            "common_gcd_degree": common_gcd.total_degree(),
            "mu6_numerator_factor_tests": factor_tests,
        }
        return result
    if determinant.degree() == 0:
        norm = base.convert(determinant.nth(0)) ** arithmetic.degree
    else:
        extension_matrix = [
            [base.zero for _ in range(arithmetic.degree)]
            for _ in range(arithmetic.degree)
        ]
        for column in range(arithmetic.degree):
            product = arithmetic.mul(
                determinant,
                arithmetic.power(arithmetic.generator, column),
            )
            for row in range(arithmetic.degree):
                extension_matrix[row][column] = base.convert(
                    product.nth(row)
                )
        norm = determinant_over_base(extension_matrix, base)
    norm_expression = sp.cancel(base.to_sympy(norm))
    norm_numerator, norm_denominator = sp.fraction(norm_expression)
    numerator_poly = sp.Poly(
        norm_numerator,
        *base_symbols,
        domain=sp.QQ,
    )
    denominator_poly = sp.Poly(
        norm_denominator,
        *base_symbols,
        domain=sp.QQ,
    )
    numerator_factors = sp.factor_list(numerator_poly.as_expr())[1]
    denominator_factors = sp.factor_list(denominator_poly.as_expr())[1]
    result["mu6_norm"] = {
            "numerator": sp.sstr(numerator_poly.as_expr()),
            "numerator_degree": numerator_poly.total_degree(),
            "numerator_terms": len(numerator_poly.terms()),
            "numerator_factors": [
                {
                    "factor": sp.sstr(factor),
                    "multiplicity": multiplicity,
                }
                for factor, multiplicity in numerator_factors
            ],
            "denominator": sp.sstr(denominator_poly.as_expr()),
            "denominator_degree": denominator_poly.total_degree(),
            "denominator_terms": len(denominator_poly.terms()),
            "denominator_factors": [
                {
                    "factor": sp.sstr(factor),
                    "multiplicity": multiplicity,
                }
                for factor, multiplicity in denominator_factors
            ],
        }
    return result


def main() -> None:
    arguments = parse_arguments()
    assert 1 <= arguments.timeout <= 60
    assignments = parse_specializations(arguments.specialize)
    if arguments.resume_checkpoint is not None:
        if arguments.pencil_factor_degree is None:
            raise ValueError(
                "--resume-checkpoint requires --pencil-factor-degree"
            )
        checkpoint_path = arguments.resume_checkpoint
        if not checkpoint_path.is_absolute():
            checkpoint_path = ROOT / checkpoint_path
        payload = resume_pencil_checkpoint(
            checkpoint_path,
            arguments.pencil_factor_degree,
        )
        payload["reproduction_command"] = " ".join(sys.argv)
        if arguments.output is not None:
            output = arguments.output
            if not output.is_absolute():
                output = ROOT / output
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        print(json.dumps(payload, indent=2, sort_keys=True))
        return
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")

    residual_raw, pivot_a_raw, pivot_b_raw = residual_data()
    (
        residual,
        residual_terms,
        residual_u_power,
        residual_content,
    ) = adapted_polynomial(residual_raw)
    (
        pivot_a,
        pivot_a_terms,
        pivot_a_u_power,
        pivot_a_content,
    ) = adapted_polynomial(pivot_a_raw)
    (
        pivot_b,
        pivot_b_terms,
        pivot_b_u_power,
        pivot_b_content,
    ) = adapted_polynomial(pivot_b_raw)
    residual, residual_terms, residual_special_content = specialized_polynomial(
        residual,
        assignments,
    )
    pivot_a, pivot_a_terms, pivot_a_special_content = specialized_polynomial(
        pivot_a,
        assignments,
    )
    pivot_b, pivot_b_terms, pivot_b_special_content = specialized_polynomial(
        pivot_b,
        assignments,
    )
    residual_content *= residual_special_content
    pivot_a_content *= pivot_a_special_content
    pivot_b_content *= pivot_b_special_content
    if arguments.prepare_only:
        print(
            json.dumps(
                {
                    "R20_terms": residual_terms,
                    "A_terms": pivot_a_terms,
                    "B_terms": pivot_b_terms,
                    "R20_u_power": residual_u_power,
                    "A_u_power": pivot_a_u_power,
                    "B_u_power": pivot_b_u_power,
                    "A_content": str(pivot_a_content),
                    "B_content": str(pivot_b_content),
                    "specializations": {
                        name: str(value)
                        for name, value in assignments.items()
                    },
                    "R20_length": len(residual),
                    "A_length": len(pivot_a),
                    "B_length": len(pivot_b),
                    "R20": residual if len(residual) <= 1000 else None,
                    "A": pivot_a if len(pivot_a) <= 1000 else None,
                    "B": pivot_b if len(pivot_b) <= 1000 else None,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return
    export = t0_open_localized_export(
        singular,
        tuple(range(2, arguments.through + 1)),
        0,
        arguments.timeout,
    )
    moment_polynomials = dict(
        zip(
            range(3, arguments.through + 1),
            export["polynomials"][:-1],
            strict=True,
        )
    )
    leading_payload = json.loads(LEADING_ARTIFACT.read_text(encoding="utf-8"))
    border_polynomial = leading_payload["leading_coefficient_lcm"]
    q_replacement = (("s2", "(s1^2*u-(13/3)*u)"),)
    adapted_replacement = (
        ("t1", "(u*(s1-ell))"),
        ("t2", "(a*u^2)"),
    )
    fixed_replacements = tuple(
        (name, singular_fraction(value))
        for name, value in assignments.items()
    )
    prepivot_specialized = {
        order: substitute(
            substitute(polynomial, q_replacement),
            adapted_replacement,
        )
        for order, polynomial in moment_polynomials.items()
    }
    if fixed_replacements:
        prepivot_specialized = {
            order: substitute(polynomial, fixed_replacements)
            for order, polynomial in prepivot_specialized.items()
        }
    # SymPy's parser spends most of its time on thousands of terms which
    # already vanish on this adapted slice (notably powers of s1 and ell in
    # p8).  Let Singular collect those zeros before handing p8 to the exact
    # Kummer arithmetic below.
    if (
        arguments.engine == "kummer"
        and assignments.get("s1") == Fraction(0)
        and assignments.get("ell") == Fraction(0)
        and 8 in prepivot_specialized
    ):
        prepivot_specialized[8] = singular_simplify_polynomial(
            singular,
            prepivot_specialized[8],
            arguments.timeout,
        )
    specialized = {
        order: substitute(polynomial, (("s3", "pivotS"),))
        for order, polynomial in prepivot_specialized.items()
    }
    prepivot_border = substitute(
        border_polynomial,
        adapted_replacement,
    )
    if fixed_replacements:
        prepivot_border = substitute(
            prepivot_border,
            fixed_replacements,
        )
    specialized_border = substitute(
        prepivot_border,
        (("s3", "pivotS"),),
    )
    minpoly = re.sub(r"\bT\b", "a", residual)
    pivot_a_root = re.sub(r"\bT\b", "a", pivot_a)
    pivot_b_root = re.sub(r"\bT\b", "a", pivot_b)
    declarations = "\n".join(
        f"poly p{order}={specialized[order]};"
        for order in range(4, arguments.through + 1)
    )
    reductions = "\n".join(
        (
            f"poly r{order}=reduce(p{order},G); "
            f'print("NF {order} "+string(size(r{order})));'
        )
        for order in range(6, arguments.through + 1)
    )
    unit_program = ""
    if arguments.through >= 6:
        generators = ",".join(
            f"p{order}" for order in range(4, arguments.through + 1)
        )
        unit_program = f"""
ideal U=std(ideal({generators}));
print("UNIT "+string(reduce(1,U)==0)+" "+string(size(U)));
"""
    remaining_parameters = [
        variable
        for variable in ("s1", "ell", "u")
        if variable not in assignments
    ]
    engine = arguments.engine
    if engine == "auto":
        engine = (
            "minpoly"
            if not remaining_parameters
            else ("python" if len(remaining_parameters) == 1 else "extension")
        )
    if engine == "minpoly" and remaining_parameters:
        raise ValueError(
            "Singular minpoly is valid here only after s1,ell,u are fixed"
        )
    if engine == "extension" and not remaining_parameters:
        raise ValueError("use the minpoly engine for a closed number field")
    if engine == "python" and len(remaining_parameters) != 1:
        raise ValueError(
            "the bounded Python engine currently requires one base parameter"
        )
    if engine == "kummer" and not (
        assignments.get("s1") == 0
        and assignments.get("ell") == 0
        and remaining_parameters == ["u"]
    ):
        raise ValueError(
            "the Kummer projection requires s1=0, ell=0 with u free"
        )
    pivot_u_difference = pivot_b_u_power - pivot_a_u_power
    pivot_scalar = Fraction(
        int(pivot_b_content.p),
        int(pivot_b_content.q),
    ) / Fraction(
        int(pivot_a_content.p),
        int(pivot_a_content.q),
    )
    if "u" in assignments:
        pivot_scalar *= assignments["u"] ** pivot_u_difference
        pivot_u_scale = singular_fraction(pivot_scalar)
    elif pivot_u_difference == 0:
        pivot_u_scale = singular_fraction(pivot_scalar)
    elif pivot_u_difference > 0:
        pivot_u_scale = (
            f"({singular_fraction(pivot_scalar)})*u^{pivot_u_difference}"
        )
    else:
        pivot_u_scale = (
            f"({singular_fraction(pivot_scalar)})"
            f"/u^{abs(pivot_u_difference)}"
        )

    if engine in ("python", "kummer"):
        if arguments.through < 6:
            raise ValueError("the Python engine requires --through at least 6")
        python_certificate = python_extension_certificate(
            residual,
            pivot_a,
            pivot_b,
            pivot_u_scale,
            prepivot_specialized,
            prepivot_border,
            remaining_parameters,
            arguments.python_stage,
            "kummer" if engine == "kummer" else "quintic",
            arguments.pencil_factor_degree,
        )
        payload = {
            "format": "two-pair-sic-bidegree33-t0-Q-residual-v2",
            "status": (
                "exact characteristic-zero degree-five quotient arithmetic "
                "on a one-parameter slice of the irreducible Q-border "
                "residual component; the full residual component and the "
                "exceptional norm roots remain open"
            ),
            "engine": (
                "python degree-four Kummer vector arithmetic"
                if engine == "kummer"
                else "python degree-five vector arithmetic"
            ),
            "completed_stage": arguments.python_stage,
            "pencil_factor_degree": arguments.pencil_factor_degree,
            "coefficient_field": (
                "QQ(T)[u]/(u^4*P5(T)+Q3(T))"
                if engine == "kummer"
                else f"QQ({remaining_parameters[0]})[a]/(R20)"
            ),
            "residual_factor_total_degree": 20,
            "residual_factor_t2_degree": 5,
            "dense_pivot": "s3=-B/A",
            "pivot_content_ratio": str(pivot_scalar),
            "specializations": {
                name: str(value) for name, value in assignments.items()
            },
            "remaining_parameter": (
                "T=t2/u^2"
                if engine == "kummer"
                else remaining_parameters[0]
            ),
            "adapted_coordinates": {
                "ell": "s1*u-t1",
                "ratio_ell": "(s1*u-t1)/u",
                "a": "t2/u^2",
            },
            "adapted_term_counts": {
                "R20": residual_terms,
                "A": pivot_a_terms,
                "B": pivot_b_terms,
            },
            **python_certificate,
            "reproduction_command": " ".join(sys.argv),
        }
        if arguments.output is not None:
            output = arguments.output
            if not output.is_absolute():
                output = ROOT / output
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        print(json.dumps(payload, indent=2, sort_keys=True))
        return

    if engine == "minpoly":
        singular_program = f"""
ring residual=(0,a),(s6,s5),dp;
minpoly={minpoly};
option(redSB);
number pivotS=-({pivot_u_scale})*({pivot_b_root})/({pivot_a_root});
{declarations}
poly p3={specialized[3]};
poly border={specialized_border};
print("BASE "+string(size(p3))+" "+string(size(border)));
ideal G=std(ideal(p4,p5));
{reductions}
print(
  "META "+string(size(G))+" "+string(dim(G))+" "+string(vdim(G))
);
print("EXTENSION 5 1");
int basisIndex;
for(basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEAD "+string(leadexp(G[basisIndex])));
}}
{unit_program}
"""
    else:
        pivot_inverse = singular_bezout_inverse(
            singular,
            remaining_parameters,
            minpoly,
            pivot_a_root,
            arguments.timeout,
        )
        extension_declarations = declarations
        extension_reductions = reductions
        extension_unit_program = unit_program.replace(
            "ideal(",
            "ideal(minpolyR,",
            1,
        )
        singular_program = f"""
ring residual=(0,{",".join(remaining_parameters)}),(
  a,s6,s5
),(dp(1),dp(2));
option(redSB);
poly minpolyR={minpoly};
poly pivotA={pivot_a_root};
poly pivotB={pivot_b_root};
ideal Rbasis=std(ideal(minpolyR));
poly pivotInverse={pivot_inverse};
print("BEZOUT 1 "+string(reduce(pivotA*pivotInverse-1,Rbasis)==0));
poly pivotS=-({pivot_u_scale})*pivotB*pivotInverse;
pivotS=reduce(pivotS,Rbasis);
{extension_declarations}
poly p3={specialized[3]};
poly border={specialized_border};
p3=reduce(p3,Rbasis);
border=reduce(border,Rbasis);
print("BASE "+string(size(p3))+" "+string(size(border)));
ideal G=std(ideal(minpolyR,p4,p5));
{extension_reductions}
print(
  "META "+string(size(G))+" "+string(dim(G))+" "+string(vdim(G))
);
print("EXTENSION "+string(deg(minpolyR,a))+" 5");
int basisIndex;
for(basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEAD "+string(leadexp(G[basisIndex])));
}}
{extension_unit_program}
"""

    completed = subprocess.run(
        [singular, "-q"],
        input=singular_program,
        text=True,
        capture_output=True,
        check=True,
        timeout=arguments.timeout,
    )
    assert "\n   ? " not in completed.stdout, (
        completed.stdout[-8000:],
        completed.stderr[-4000:],
    )
    meta = re.search(
        r"(?m)^META (\d+) (-?\d+) (-?\d+)$",
        completed.stdout,
    )
    base = re.search(r"(?m)^BASE (\d+) (\d+)$", completed.stdout)
    extension = re.search(
        r"(?m)^EXTENSION (\d+) (\d+)$",
        completed.stdout,
    )
    bezout = re.search(r"(?m)^BEZOUT (\d+) ([01])$", completed.stdout)
    unit = re.search(
        r"(?m)^UNIT ([01]) (\d+)$",
        completed.stdout,
    )
    assert (
        meta is not None
        and base is not None
        and extension is not None
    ), completed.stdout[-8000:]
    if engine == "extension":
        assert bezout is not None and bezout.groups() == ("1", "1"), (
            completed.stdout[-8000:]
        )
    assert base.groups() == ("0", "0"), (
        "the residual minpoly and pivot must annihilate mu3 and the border",
        completed.stdout[-8000:],
    )
    if arguments.through >= 6:
        assert unit is not None, completed.stdout[-8000:]
    normal_form_counts = {
        f"mu{order}": int(term_count)
        for order, term_count in re.findall(
            r"(?m)^NF (\d+) (\d+)$",
            completed.stdout,
        )
    }
    payload = {
        "format": "two-pair-sic-bidegree33-t0-Q-residual-v1",
        "status": (
            "exact generic-point arithmetic on the characteristic-zero "
            "irreducible Q-border residual component; exceptional "
            "coefficient denominators are not yet classified"
        ),
        "engine": engine,
        "coefficient_field": (
            (
                "QQ[a]/(R20(a))"
                if engine == "minpoly"
                else (
                    f"QQ({','.join(remaining_parameters)})[a]/(R20), "
                    "represented by the polynomial extension basis"
                )
            )
        ),
        "adapted_coordinates": {"ell": "s1*u-t1"},
        "weight_zero_coordinates": {
            "ell": "(s1*u-t1)/u",
            "a": "t2/u^2",
        },
        "adapted_term_counts": {
            "R20": residual_terms,
            "A": pivot_a_terms,
            "B": pivot_b_terms,
        },
        "removed_u_powers": {
            "R20": residual_u_power,
            "A": pivot_a_u_power,
            "B": pivot_b_u_power,
        },
        "pivot_content_ratio": str(pivot_scalar),
        "specializations": {
            name: str(value) for name, value in assignments.items()
        },
        "residual_factor_total_degree": 20,
        "residual_factor_t2_degree": 5,
        "dense_pivot": "s3=-B/A",
        "mu3_zero_in_residual_field": base.group(1) == "0",
        "leading_border_zero_in_residual_field": base.group(2) == "0",
        "fiber_variables": ["s6", "s5"],
        "fiber_groebner_basis_size": int(meta.group(1)),
        "fiber_dimension": int(meta.group(2)),
        "extension_degree": int(extension.group(1)),
        "total_quotient_length_over_base_field": int(meta.group(3)),
        "fiber_quotient_length": (
            int(meta.group(3))
            if engine == "minpoly"
            else int(meta.group(3)) // int(extension.group(1))
        ),
        "leading_exponents": re.findall(
            r"(?m)^LEAD ([0-9,]+)$",
            completed.stdout,
        ),
        "last_moment": arguments.through,
        "normal_form_term_counts": normal_form_counts,
        "unit_ideal_through_last_moment": (
            unit.group(1) == "1" if unit is not None else None
        ),
        "adjoined_basis_size": (
            int(unit.group(2)) if unit is not None else None
        ),
        "reproduction_command": " ".join(sys.argv),
    }
    if arguments.output is not None:
        output = arguments.output
        if not output.is_absolute():
            output = ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
