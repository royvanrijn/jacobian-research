#!/usr/bin/env python3
"""Exact unspecialized arithmetic on the residual Q component.

This research driver replaces the rational-function normalization in
``research_two_pair_sic_bidegree33_t0_Q_residual.py`` by sparse
``python-flint`` arithmetic.  It works in

    QQ(s1,lam,v)[T]/(R20),     lam=(s1*u-t1)/u, v=u^2,

and keeps every coefficient denominator as a product of explicitly recorded
base polynomials.  In particular, no multivariate gcd is taken during the
quotient or border-basis calculation.

The stages are deliberately checkpointed.  ``pivot`` certifies the dense
subresultant pivot, ``moments`` parses the unspecialized corrected moments,
and ``basis`` attempts to construct the length-four fibre algebra.
``pencil`` is reserved for the planned 25 base equations furnished by
det(M6+z*M7); it is not implemented until the basis stage is made bounded.

This is a research calculation, not by itself a theorem.  The pivot-open
result must be saturated by every recorded denominator, and the complementary
pivot-exception locus must be treated separately.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from collections import defaultdict
from fractions import Fraction
import hashlib
from itertools import permutations
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

from flint import (
    fmpq,
    fmpq_mpoly,
    fmpq_mpoly_ctx,
    nmod_poly,
    nmod_mpoly_ctx,
)
import sympy as sp
from sympy.polys.agca.extensions import FiniteExtension, ExtensionElement

from research_two_pair_sic_bidegree33_t0_Q_residual import (
    LEADING_ARTIFACT,
    ROOT,
    residual_data,
)
from verify_two_pair_sic_bidegree33_boundary_generic_quotient import (
    substitute,
)
from verify_two_pair_sic_bidegree33_corrected_boundary import (
    t0_open_localized_export,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_t0_stratum_Q_unspecialized_flint.json"
)
VARIABLES = ("T", "s1", "lam", "v")
CTX = fmpq_mpoly_ctx.get(VARIABLES)
T, S1, LAM, V = CTX.gens()
PRIME = 0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage",
        choices=(
            "prepare",
            "pivot",
            "moments",
            "system",
            "solve",
            "evaluated",
            "remainder5",
            "remainder6",
            "remainders",
            "prepivot5",
            "prepivot6",
            "prepivot7",
            "prepivot",
            "raw_equations",
            "prepivot_equations",
            "prepivot_quadratic",
            "prepivot_cross6",
            "prepivot_cross7",
            "raw8",
            "raw8_solve",
            "raw8_popen",
            "raw8_popen_solve",
            "raw8_exceptional",
            "raw8_exceptional_solve",
            "raw8_direct",
            "raw8_direct_solve",
            "raw8_direct_generic",
            "raw8_direct_generic_lift",
            "raw8_direct_extension",
            "raw8_direct_extension_lift",
            "raw8_direct_extension_gcd",
            "raw8_direct_custom_gcd",
            "raw8_direct_m2_gcd",
            "raw8_direct_m2_exceptional",
            "raw8_direct_m2_open6",
            "raw8_direct_sympy_exceptional",
            "raw8_direct_sympy_open6",
            "raw8_direct_groebner_julia",
            "raw8_direct_popen_groebner_julia",
            "raw8_direct_popen_groebner_julia_change",
            "raw8_popen_groebner_julia",
            "raw8_popen_groebner_julia_change",
            "raw8_exceptional_groebner_julia",
            "raw8_exceptional_groebner_julia_change",
            "basis",
            "pencil",
        ),
        default="pivot",
    )
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument(
        "--moment-through",
        type=int,
        default=7,
        help="expand corrected moments mu4 through this order",
    )
    parser.add_argument(
        "--prime",
        type=int,
        default=0,
        help="0 for characteristic zero, otherwise an odd FLINT word prime",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--s1-value",
        type=int,
        default=None,
        help="specialize s1 before the direct-system solve",
    )
    parser.add_argument(
        "--lam-value",
        type=int,
        default=None,
        help="specialize lam in a split function-field scout",
    )
    parser.add_argument(
        "--v-value",
        type=int,
        default=None,
        help="specialize v in a split function-field scout",
    )
    parser.add_argument(
        "--solver-seconds",
        type=int,
        default=0,
        help="msolve deadline in seconds; 0 means no deadline",
    )
    parser.add_argument(
        "--solver-result",
        type=Path,
        default=None,
        help="persist the full native msolve result outside the temp directory",
    )
    parser.add_argument(
        "--julia-project",
        type=Path,
        default=Path("/tmp/sic33-groebner-env"),
        help="Julia environment containing Groebner.jl and AbstractAlgebra.jl",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    return parser.parse_args()


def configure_coefficient_ring(prime: int) -> None:
    global CTX, T, S1, LAM, V, PRIME
    PRIME = prime
    if prime:
        CTX = nmod_mpoly_ctx.get(VARIABLES, modulus=prime)
    else:
        CTX = fmpq_mpoly_ctx.get(VARIABLES)
    T, S1, LAM, V = CTX.gens()


def coefficient_from_fraction(
    numerator: int,
    denominator: int = 1,
):
    if PRIME:
        return (numerator % PRIME) * pow(denominator, -1, PRIME) % PRIME
    return fmpq(numerator, denominator)


def coefficient_from_fmpq(value: fmpq):
    return coefficient_from_fraction(int(value.p), int(value.q))


def polynomial_profile(value: fmpq_mpoly) -> dict[str, object]:
    return {
        "terms": len(value.to_dict()),
        "degrees": {
            variable: int(degree)
            for variable, degree in zip(
                VARIABLES,
                value.degrees(),
                strict=True,
            )
        },
        "total_degree": int(value.total_degree()),
        "sha256": hashlib.sha256(str(value).encode()).hexdigest(),
    }


def sympy_even_u_to_flint(expression: str) -> fmpq_mpoly:
    """Parse a polynomial and replace every ``u^(2j)`` by ``v^j``."""

    s1, lam, root, u = sp.symbols("s1 lam T u")
    polynomial = sp.Poly(
        sp.sympify(
            expression.replace("^", "**"),
            locals={"s1": s1, "lam": lam, "T": root, "u": u},
        ),
        root,
        s1,
        lam,
        u,
        domain=sp.QQ,
    )
    converted: dict[tuple[int, ...], fmpq] = {}
    for (root_power, s1_power, ell_power, u_power), coefficient in (
        polynomial.terms()
    ):
        if u_power % 2:
            raise ValueError(f"odd u exponent {u_power}")
        converted[(root_power, s1_power, ell_power, u_power // 2)] = (
            coefficient_from_fraction(
                int(coefficient.p),
                int(coefficient.q),
            )
        )
    return CTX.from_dict(converted)


def root_degree(value: fmpq_mpoly) -> int:
    return max((monomial[0] for monomial in value.to_dict()), default=-1)


def root_coefficient(value: fmpq_mpoly, degree: int) -> fmpq_mpoly:
    return CTX.from_dict(
        {
            (0, monomial[1], monomial[2], monomial[3]): coefficient
            for monomial, coefficient in value.to_dict().items()
            if monomial[0] == degree
        }
    )


class DenominatorRegistry:
    """Named base factors used by the fraction-free quotient arithmetic."""

    def __init__(self, leading_coefficient: fmpq_mpoly):
        self.names = ["lc_R20"]
        self.factors = [leading_coefficient]

    def register(self, name: str, factor: fmpq_mpoly) -> int:
        if root_degree(factor) > 0:
            raise ValueError("a coefficient denominator may not contain T")
        for index, existing in enumerate(self.factors):
            if factor == existing:
                return index
            if factor == -existing:
                return index
        self.names.append(name)
        self.factors.append(factor)
        return len(self.factors) - 1

    def pad(self, denominator: tuple[int, ...]) -> tuple[int, ...]:
        return denominator + (0,) * (len(self.factors) - len(denominator))

    def product(self, denominator: tuple[int, ...]) -> fmpq_mpoly:
        result = CTX.constant(1)
        for factor, exponent in zip(
            self.factors,
            self.pad(denominator),
            strict=True,
        ):
            if exponent:
                result *= factor**exponent
        return result

    def describe(self) -> list[dict[str, object]]:
        return [
            {"name": name, **polynomial_profile(factor)}
            for name, factor in zip(self.names, self.factors, strict=True)
        ]


@dataclass(frozen=True)
class QuotientElement:
    numerator: fmpq_mpoly
    denominator: tuple[int, ...]


class QuotientArithmetic:
    """Arithmetic in Frac(QQ[s1,lam,v])[T]/(R20)."""

    def __init__(self, modulus: fmpq_mpoly):
        self.modulus = modulus
        self.degree = root_degree(modulus)
        if self.degree != 5:
            raise ValueError(f"expected a quintic modulus, got {self.degree}")
        self.leading = root_coefficient(modulus, self.degree)
        self.registry = DenominatorRegistry(self.leading)
        self.zero = QuotientElement(CTX.constant(0), ())
        self.one = QuotientElement(CTX.constant(1), ())
        self.generator = QuotientElement(T, ())

    def _cancel_known(self, value: QuotientElement) -> QuotientElement:
        numerator = value.numerator
        denominator = list(self.registry.pad(value.denominator))
        for index, exponent in enumerate(denominator):
            # The first factor is the monomial leading coefficient of R20,
            # so exact division is cheap and controls pseudo-reduction
            # growth.  Trial division by a later multithousand-term norm is
            # more expensive than the arithmetic it is meant to simplify.
            if index:
                continue
            factor = self.registry.factors[index]
            while exponent and not numerator.is_zero():
                quotient, remainder = divmod(numerator, factor)
                if not remainder.is_zero():
                    break
                numerator = quotient
                exponent -= 1
            denominator[index] = exponent
        while denominator and denominator[-1] == 0:
            denominator.pop()
        return QuotientElement(numerator, tuple(denominator))

    def _reduce(
        self,
        numerator: fmpq_mpoly,
        denominator: tuple[int, ...],
    ) -> QuotientElement:
        denominator = list(self.registry.pad(denominator))
        while root_degree(numerator) >= self.degree:
            degree = root_degree(numerator)
            coefficient = root_coefficient(numerator, degree)
            numerator = (
                self.leading * numerator
                - coefficient * T ** (degree - self.degree) * self.modulus
            )
            denominator[0] += 1
        return self._cancel_known(
            QuotientElement(numerator, tuple(denominator))
        )

    def make(
        self,
        numerator: fmpq_mpoly | int,
        denominator: tuple[int, ...] = (),
    ) -> QuotientElement:
        if isinstance(numerator, int):
            numerator = CTX.constant(numerator)
        return self._reduce(numerator, denominator)

    def add(
        self,
        left: QuotientElement,
        right: QuotientElement,
    ) -> QuotientElement:
        left_den = self.registry.pad(left.denominator)
        right_den = self.registry.pad(right.denominator)
        common = tuple(max(a, b) for a, b in zip(left_den, right_den))
        left_scale = tuple(c - a for c, a in zip(common, left_den))
        right_scale = tuple(c - b for c, b in zip(common, right_den))
        numerator = (
            left.numerator * self.registry.product(left_scale)
            + right.numerator * self.registry.product(right_scale)
        )
        return self._cancel_known(QuotientElement(numerator, common))

    def neg(self, value: QuotientElement) -> QuotientElement:
        return QuotientElement(-value.numerator, value.denominator)

    def equal(
        self,
        left: QuotientElement,
        right: QuotientElement,
    ) -> bool:
        return (
            left.numerator * self.registry.product(right.denominator)
            == right.numerator * self.registry.product(left.denominator)
        )

    def sub(
        self,
        left: QuotientElement,
        right: QuotientElement,
    ) -> QuotientElement:
        return self.add(left, self.neg(right))

    def mul(
        self,
        left: QuotientElement,
        right: QuotientElement,
    ) -> QuotientElement:
        left_den = self.registry.pad(left.denominator)
        right_den = self.registry.pad(right.denominator)
        denominator = tuple(a + b for a, b in zip(left_den, right_den))
        return self._reduce(left.numerator * right.numerator, denominator)

    def power(self, value: QuotientElement, exponent: int) -> QuotientElement:
        result = self.one
        factor = value
        while exponent:
            if exponent & 1:
                result = self.mul(result, factor)
            exponent //= 2
            if exponent:
                factor = self.mul(factor, factor)
        return result

    def scalar_mul(
        self,
        value: QuotientElement,
        scalar: fmpq_mpoly,
    ) -> QuotientElement:
        return self._cancel_known(
            QuotientElement(value.numerator * scalar, value.denominator)
        )

    def multiplication_matrix(
        self,
        value: QuotientElement,
    ) -> list[list[QuotientElement]]:
        matrix = [[self.zero for _ in range(self.degree)] for _ in range(self.degree)]
        power = self.one
        for column in range(self.degree):
            product = self.mul(value, power)
            for monomial, coefficient in product.numerator.to_dict().items():
                row = monomial[0]
                base_monomial = (0, monomial[1], monomial[2], monomial[3])
                entry = QuotientElement(
                    CTX.term(coefficient, base_monomial),
                    product.denominator,
                )
                matrix[row][column] = self.add(matrix[row][column], entry)
            power = self.mul(power, self.generator)
        return matrix

    def determinant(
        self,
        matrix: list[list[QuotientElement]],
    ) -> QuotientElement:
        answer = self.zero
        for permutation in permutations(range(len(matrix))):
            inversions = sum(
                permutation[left] > permutation[right]
                for left in range(len(matrix))
                for right in range(left + 1, len(matrix))
            )
            term = self.one
            for row, column in enumerate(permutation):
                term = self.mul(term, matrix[row][column])
            answer = self.add(
                answer,
                term if inversions % 2 == 0 else self.neg(term),
            )
        return answer

    def inverse(
        self,
        value: QuotientElement,
        factor_name: str,
    ) -> QuotientElement:
        """Invert an element by the adjugate of its multiplication matrix."""

        matrix = self.multiplication_matrix(value)
        determinant = self.determinant(matrix)
        if root_degree(determinant.numerator) > 0:
            raise AssertionError("norm retained the extension generator")
        if determinant.numerator.is_zero():
            raise ZeroDivisionError("zero norm")
        denominator_factor = self.registry.register(
            factor_name,
            determinant.numerator,
        )
        coordinates: list[QuotientElement] = []
        for column in range(self.degree):
            minor = [
                [matrix[row][other] for other in range(self.degree) if other != column]
                for row in range(1, self.degree)
            ]
            cofactor = self.determinant(minor)
            if column % 2:
                cofactor = self.neg(cofactor)
            numerator = cofactor.numerator
            # Divide by det = det.num/det.den.
            numerator *= self.registry.product(determinant.denominator)
            denominator = list(self.registry.pad(cofactor.denominator))
            denominator[denominator_factor] += 1
            coordinates.append(
                self._cancel_known(
                    QuotientElement(numerator, tuple(denominator))
                )
            )
        inverse = self.zero
        for degree, coordinate in enumerate(coordinates):
            inverse = self.add(
                inverse,
                self.mul(coordinate, self.power(self.generator, degree)),
            )
        # The matrix inverted multiplication by num/value.den.
        inverse = self.scalar_mul(
            inverse,
            self.registry.product(value.denominator),
        )
        check = self.mul(value, inverse)
        if not self.equal(check, self.one):
            raise AssertionError("adjugate inverse verification failed")
        return inverse


def determinant_polynomial(
    arithmetic: QuotientArithmetic,
    matrix: list[list[QuotientElement]],
) -> list[QuotientElement]:
    """Return coefficients of a determinant whose entries are linear in z."""

    size = len(matrix)
    raise NotImplementedError(size)


def adapted_ratio_polynomial(
    expression: str,
) -> tuple[str, int, int, sp.Rational]:
    """Use ``lam=(s1*u-t1)/u`` and remove common u-power/content."""

    s1, t1, t2, u, lam, root = sp.symbols("s1 t1 t2 u lam T")
    parsed = sp.sympify(
        expression.replace("^", "**"),
        locals={"s1": s1, "t1": t1, "t2": t2, "u": u},
    )
    polynomial = sp.Poly(
        sp.expand(
            parsed.subs(
                {
                    t1: u * (s1 - lam),
                    t2: root * u**2,
                },
                simultaneous=True,
            )
        ),
        s1,
        lam,
        root,
        u,
        domain=sp.QQ,
    )
    u_power = min(monomial[3] for monomial, _ in polynomial.terms())
    if u_power:
        polynomial = sp.Poly(
            polynomial.as_expr() / u**u_power,
            s1,
            lam,
            root,
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


def flint_input_data() -> tuple[
    QuotientArithmetic,
    fmpq_mpoly,
    fmpq_mpoly,
    dict[str, object],
]:
    residual_raw, pivot_a_raw, pivot_b_raw = residual_data()
    residual, residual_terms, residual_u_power, residual_content = (
        adapted_ratio_polynomial(residual_raw)
    )
    pivot_a, pivot_a_terms, pivot_a_u_power, pivot_a_content = (
        adapted_ratio_polynomial(pivot_a_raw)
    )
    pivot_b, pivot_b_terms, pivot_b_u_power, pivot_b_content = (
        adapted_ratio_polynomial(pivot_b_raw)
    )
    if pivot_a_u_power != pivot_b_u_power:
        raise AssertionError("the dense pivot unexpectedly carries an odd u scale")
    modulus = sympy_even_u_to_flint(residual)
    a_value = sympy_even_u_to_flint(pivot_a)
    b_value = sympy_even_u_to_flint(pivot_b)
    arithmetic = QuotientArithmetic(modulus)
    metadata = {
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
        "contents": {
            "R20": str(residual_content),
            "A": str(pivot_a_content),
            "B": str(pivot_b_content),
        },
        "R20": polynomial_profile(modulus),
        "A": polynomial_profile(a_value),
        "B": polynomial_profile(b_value),
    }
    return arithmetic, a_value, b_value, metadata


def compute_pivot(
    arithmetic: QuotientArithmetic,
    a_value: fmpq_mpoly,
    b_value: fmpq_mpoly,
    a_content: str,
    b_content: str,
) -> tuple[QuotientElement, dict[str, object]]:
    content_ratio = Fraction(b_content) / Fraction(a_content)
    a_element = arithmetic.make(a_value)
    b_element = arithmetic.make(b_value)
    started = time.monotonic()
    a_inverse = arithmetic.inverse(a_element, "norm_A")
    pivot = arithmetic.scalar_mul(
        arithmetic.neg(arithmetic.mul(b_element, a_inverse)),
        CTX.constant(
            coefficient_from_fraction(
                content_ratio.numerator,
                content_ratio.denominator,
            )
        ),
    )
    elapsed = time.monotonic() - started
    check = arithmetic.add(
        arithmetic.scalar_mul(
            arithmetic.mul(a_element, pivot),
            CTX.constant(coefficient_from_fraction(int(a_content))),
        ),
        arithmetic.scalar_mul(
            b_element,
            CTX.constant(coefficient_from_fraction(int(b_content))),
        ),
    )
    if not check.numerator.is_zero():
        raise AssertionError("dense pivot equation did not reduce to zero")
    return pivot, {
        "seconds": round(elapsed, 6),
        "content_ratio_B_over_A": str(content_ratio),
        "numerator": polynomial_profile(pivot.numerator),
        "denominator_exponents": list(
            arithmetic.registry.pad(pivot.denominator)
        ),
        "denominator_factors": arithmetic.registry.describe(),
        "verified_A_times_s3_plus_B": True,
    }


RawFibrePolynomial = dict[tuple[int, int, int], fmpq_mpoly]
FibrePolynomial = dict[tuple[int, int], QuotientElement]
RawQuotientPolynomial = dict[tuple[int, int, int], QuotientElement]
YZPolynomial = dict[tuple[int, int], QuotientElement]


def expanded_adapted_moments(
    singular: str,
    timeout: int,
    maximum_order: int = 7,
) -> dict[int, str]:
    """Expand later moments after imposing Q and the weight-zero coordinates."""

    if maximum_order < 4:
        raise ValueError("maximum moment order must be at least four")
    orders = tuple(range(2, maximum_order + 1))
    export = t0_open_localized_export(singular, orders, 0, timeout)
    raw = dict(
        zip(orders[1:], export["polynomials"][:-1], strict=True)
    )
    replacements = (
        ("s2", "(s1^2*u-(13/3)*u)"),
        ("t1", "(u*(s1-lam))"),
        ("t2", "(T*u^2)"),
    )
    declarations = "\n".join(
        f"poly p{order}={substitute(raw[order], replacements)};"
        for order in range(4, maximum_order + 1)
    )
    prints = "\n".join(
        f'print("MOMENT{order} "+string(p{order}));'
        for order in range(4, maximum_order + 1)
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring adapted=0,(s6,s5,s3,T,s1,lam,u),dp;
{declarations}
{prints}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    if "\n   ? " in completed.stdout:
        raise AssertionError(completed.stdout[-8000:])
    moments: dict[int, str] = {}
    for order in range(4, maximum_order + 1):
        marker = re.search(
            rf"(?m)^MOMENT{order} (.*)$",
            completed.stdout,
        )
        if marker is None:
            raise AssertionError(completed.stdout[-8000:])
        moments[order] = marker.group(1)
    return moments


def _expanded_terms(expression: str) -> list[tuple[int, str]]:
    expression = expression.replace(" ", "")
    if expression == "0":
        return []
    starts = [
        index
        for index, character in enumerate(expression)
        if index and character in "+-"
    ]
    pieces = []
    previous = 0
    for boundary in starts + [len(expression)]:
        piece = expression[previous:boundary]
        sign = -1 if piece.startswith("-") else 1
        if piece[:1] in "+-":
            piece = piece[1:]
        pieces.append((sign, piece))
        previous = boundary
    return pieces


def parse_raw_fibre_polynomial(
    expression: str,
) -> tuple[RawFibrePolynomial, dict[str, object]]:
    """Parse one moment after ``s6=u*x6`` and replace ``u^2`` by ``v``."""

    parsed_terms: list[
        tuple[int, dict[str, int], fmpq]
    ] = []
    variable_names = {"s6", "s5", "s3", "T", "s1", "lam", "u"}
    for sign, term in _expanded_terms(expression):
        exponents: dict[str, int] = defaultdict(int)
        coefficient = fmpq(sign)
        for factor in term.split("*"):
            match = re.fullmatch(r"([A-Za-z][A-Za-z0-9]*)(?:\^(\d+))?", factor)
            if match is not None:
                variable, exponent = match.groups()
                if variable not in variable_names:
                    raise ValueError(f"unknown variable {variable}")
                exponents[variable] += int(exponent or 1)
            elif re.fullmatch(r"\d+(?:/\d+)?", factor):
                coefficient *= fmpq(factor)
            else:
                raise ValueError(f"cannot parse factor {factor!r}")
        parsed_terms.append((sign, exponents, coefficient))
    # On the u-open chart the substitution s6=u*x6 makes every coefficient
    # have one u parity.  This retains the smaller base ring
    # QQ[s1,lam,u^2].
    effective_u_powers = [
        exponents["u"] + exponents["s6"]
        for _, exponents, _ in parsed_terms
    ]
    minimum_u = min(
        effective_u_powers,
        default=0,
    )
    by_fibre: dict[
        tuple[int, int, int],
        dict[tuple[int, int, int, int], fmpq],
    ] = defaultdict(dict)
    odd_differences = set()
    for _, exponents, coefficient in parsed_terms:
        u_difference = (
            exponents["u"] + exponents["s6"] - minimum_u
        )
        odd_differences.add(u_difference % 2)
        monomial = (
            exponents["T"],
            exponents["s1"],
            exponents["lam"],
            u_difference // 2,
        )
        fibre = (
            exponents["s6"],
            exponents["s5"],
            exponents["s3"],
        )
        by_fibre[fibre][monomial] = (
            by_fibre[fibre].get(monomial, fmpq(0)) + coefficient
        )
    if odd_differences - {0}:
        raise ValueError(
            f"u exponents are not of one parity: {sorted(odd_differences)}"
        )
    answer = {
        fibre: CTX.from_dict(
            {
                monomial: coefficient_from_fmpq(coefficient)
                for monomial, coefficient in terms.items()
                if coefficient
            }
        )
        for fibre, terms in by_fibre.items()
    }
    answer = {
        fibre: coefficient
        for fibre, coefficient in answer.items()
        if not coefficient.is_zero()
    }
    return answer, {
        "expanded_terms": len(parsed_terms),
        "fibre_monomials": len(answer),
        "removed_u_power": minimum_u,
        "fibre_rescaling": "s6=u*x6 (x6 is subsequently printed as s6)",
        "maximum_s3_degree": max(
            (monomial[2] for monomial in answer),
            default=-1,
        ),
        "maximum_fibre_degree": max(
            (monomial[0] + monomial[1] for monomial in answer),
            default=-1,
        ),
    }


def parse_moments(
    singular: str,
    timeout: int,
    maximum_order: int = 7,
) -> tuple[dict[int, RawFibrePolynomial], dict[str, object]]:
    expanded = expanded_adapted_moments(
        singular,
        timeout,
        maximum_order,
    )
    moments: dict[int, RawFibrePolynomial] = {}
    profiles: dict[str, object] = {}
    for order, expression in expanded.items():
        moments[order], profiles[f"mu{order}"] = (
            parse_raw_fibre_polynomial(expression)
        )
        profiles[f"mu{order}"]["expanded_length"] = len(expression)
    return moments, profiles


def serialize_raw_fibre_polynomial(
    polynomial: RawFibrePolynomial,
) -> str:
    terms = []
    for (s6_power, s5_power, s3_power), coefficient in sorted(
        polynomial.items()
    ):
        factors = [f"({coefficient})"]
        for variable, exponent in (
            ("s6", s6_power),
            ("s5", s5_power),
            ("s3", s3_power),
        ):
            if exponent == 1:
                factors.append(variable)
            elif exponent:
                factors.append(f"{variable}^{exponent}")
        terms.append("*".join(factors))
    return "+".join(terms) or "0"


def direct_residual_system(
    arithmetic: QuotientArithmetic,
    a_value: fmpq_mpoly,
    b_value: fmpq_mpoly,
    moments: dict[int, RawFibrePolynomial],
    s1_value: int | None = None,
) -> dict[str, object]:
    """Return the polynomial system before any function-field division."""

    def specialize(value):
        return (
            value
            if s1_value is None
            else value.subs({"s1": s1_value})
        )

    specialized_moments = {
        order: {
            monomial: specialize(coefficient)
            for monomial, coefficient in polynomial.items()
        }
        for order, polynomial in moments.items()
    }
    polynomials = [
        str(specialize(arithmetic.modulus)),
        (
            f"(6)*({specialize(a_value)})*s3"
            f"+({specialize(b_value)})"
        ),
        *[
            serialize_raw_fibre_polynomial(specialized_moments[order])
            for order in range(4, 8)
        ],
        "v*vinv-1",
    ]
    variables = [
        "s6",
        "s5",
        "s3",
        "T",
        "v",
        "vinv",
        "lam",
    ]
    if s1_value is None:
        variables.insert(4, "s1")
    return {
        "ordinary_variables": variables,
        "polynomials": polynomials,
        "equation_profiles": [
            {
                "index": index,
                "length": len(polynomial),
                "sha256": hashlib.sha256(polynomial.encode()).hexdigest(),
            }
            for index, polynomial in enumerate(polynomials)
        ],
        "scope": (
            "R20=0, 6*A*s3+B=0, mu4=...=mu7=0, and v!=0; "
            "no coefficient-field division or pivot saturation"
        ),
        "s1_value": s1_value,
    }


def run_direct_msolve_unbounded(
    system: dict[str, object],
    prime: int,
    threads: int,
    solver_seconds: int = 0,
    result_output: Path | None = None,
) -> dict[str, object]:
    """Run msolve without a wall-clock deadline; the outer job owns control."""

    msolve = shutil.which("msolve")
    if msolve is None:
        raise RuntimeError("msolve is required")
    with tempfile.TemporaryDirectory(
        prefix="sic33-Q-unspecialized-"
    ) as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        input_path.write_text(
            ",".join(system["ordinary_variables"])
            + "\n"
            + str(prime)
            + "\n"
            + ",\n".join(system["polynomials"])
            + "\n",
            encoding="utf-8",
        )
        print(
            "MSOLVE_START "
            f"prime={prime} threads={threads} "
            f"input_bytes={input_path.stat().st_size}",
            file=sys.stderr,
            flush=True,
        )
        started = time.monotonic()
        try:
            completed = subprocess.run(
                [
                    msolve,
                    "-f",
                    str(input_path),
                    "-o",
                    str(output_path),
                    "-t",
                    str(threads),
                    "-v",
                    "2",
                    "-l",
                    "2",
                ],
                text=True,
                check=False,
                timeout=(solver_seconds or None),
            )
        except subprocess.TimeoutExpired:
            return {
                "status": "timeout",
                "seconds": round(time.monotonic() - started, 3),
                "returncode": None,
                "solver_output": (
                    "inherited stdout/stderr for live job logging"
                ),
                "result_tail": "",
            }
        elapsed = time.monotonic() - started
        result = (
            output_path.read_text(encoding="utf-8")
            if output_path.exists()
            else ""
        )
        if result_output is not None and output_path.exists():
            result_output.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(output_path, result_output)
    if completed.returncode != 0:
        status = f"solver-error-{completed.returncode}"
    elif result.rstrip().endswith("[1]:") or result.rstrip().endswith(
        "[-1]:"
    ):
        status = "unit"
    else:
        status = "nonunit"
    return {
        "status": status,
        "seconds": round(elapsed, 3),
        "returncode": completed.returncode,
        "solver_output": (
            "inherited stdout/stderr for live detached-job logging"
        ),
        "result_tail": result[-4000:],
    }


def run_raw_generic_singular(
    system: dict[str, object],
    prime: int,
    solver_seconds: int,
    lift: bool,
) -> dict[str, object]:
    """Solve the four-variable fibre over ``GF(p)(s1,lam,v)``."""

    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    if system["ordinary_variables"][:7] != [
        "s6",
        "s5",
        "s3",
        "T",
        "s1",
        "lam",
        "v",
    ]:
        raise ValueError("generic raw solve requires the direct system")
    fibre_polynomials = system["polynomials"][:-2]
    ideal = ",\n".join(fibre_polynomials)
    if lift:
        computation = """
matrix transformation;
ideal G=liftstd(I,transformation,"slimgb");
proc polynomialLcm(poly left, poly right)
{
  return(left*right/gcd(left,right));
}
poly transformationDenominator=1;
poly cursor;
number coefficient;
int transformationRow;
int transformationColumn;
for(
  transformationRow=1;
  transformationRow<=nrows(transformation);
  transformationRow++
)
{
  for(
    transformationColumn=1;
    transformationColumn<=ncols(transformation);
    transformationColumn++
  )
  {
    cursor=transformation[transformationRow,transformationColumn];
    while(cursor!=0)
    {
      coefficient=leadcoef(cursor);
      transformationDenominator=polynomialLcm(
        transformationDenominator,
        denominator(coefficient)
      );
      cursor=cursor-lead(cursor);
    }
  }
}
print(
  "QRAW_GENERIC_DENOMINATOR "
  +string(transformationDenominator)
);
"""
    else:
        computation = "ideal G=slimgb(I);"
    program = f"""
ring generic=({prime},s1,lam,v),(s6,s5,s3,T),dp;
option(redSB);
ideal I=
{ideal};
timer=1;
{computation}
int elapsed=timer;
int isUnit=(reduce(1,G)==0);
print(
  "QRAW_GENERIC_META "+string(size(G))+" "+string(dim(G))+" "
  +string(isUnit)+" "+string(elapsed)
);
"""
    print(
        "SINGULAR_GENERIC_START "
        f"prime={prime} lift={int(lift)} input_bytes={len(program)}",
        file=sys.stderr,
        flush=True,
    )
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=False,
            timeout=(solver_seconds or None),
        )
    except subprocess.TimeoutExpired as error:
        return {
            "status": "timeout",
            "seconds": round(time.monotonic() - started, 3),
            "returncode": None,
            "stdout_tail": (
                error.stdout.decode()
                if isinstance(error.stdout, bytes)
                else (error.stdout or "")
            )[-4000:],
            "stderr_tail": (
                error.stderr.decode()
                if isinstance(error.stderr, bytes)
                else (error.stderr or "")
            )[-4000:],
            "program_bytes": len(program),
            "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
        }
    elapsed = time.monotonic() - started
    marker = re.search(
        r"(?m)^QRAW_GENERIC_META (-?\d+) (-?\d+) ([01]) (\d+)$",
        completed.stdout,
    )
    if completed.returncode != 0:
        status = f"solver-error-{completed.returncode}"
    elif marker is None:
        status = "missing-result-marker"
    elif marker.group(3) == "1":
        status = "unit"
    else:
        status = "nonunit"
    result: dict[str, object] = {
        "status": status,
        "seconds": round(elapsed, 3),
        "returncode": completed.returncode,
        "stdout_tail": completed.stdout[-4000:],
        "stderr_tail": completed.stderr[-30000:],
        "program_bytes": len(program),
        "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
    }
    if marker is not None:
        result.update(
            {
                "basis_size": int(marker.group(1)),
                "dimension": int(marker.group(2)),
                "unit_ideal": marker.group(3) == "1",
                "singular_timer_ticks": int(marker.group(4)),
            }
        )
    denominator = re.search(
        r"(?m)^QRAW_GENERIC_DENOMINATOR (.*)$",
        completed.stdout,
    )
    if denominator is not None:
        value = denominator.group(1)
        result.update(
            {
                "transformation_denominator": value,
                "transformation_denominator_length": len(value),
                "transformation_denominator_sha256": hashlib.sha256(
                    value.encode()
                ).hexdigest(),
            }
        )
    return result


def run_raw_extension_singular(
    system: dict[str, object],
    prime: int,
    solver_seconds: int,
    lift: bool,
) -> dict[str, object]:
    """Solve in ``GF(p)(s1,lam,v)[T]/(R20)`` after the dense pivot."""

    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    polynomials = system["polynomials"]
    if len(polynomials) != 9:
        raise ValueError("extension solve requires the direct raw system")
    reduced_declarations = "\n".join(
        f"poly q{order}=reduce(({polynomials[index]}),pivotBasis);"
        for order, index in zip(range(4, 9), range(2, 7), strict=True)
    )
    if lift:
        computation = """
matrix transformation;
ideal G=liftstd(I,transformation,"slimgb");
proc polynomialLcm(poly left, poly right)
{
  return(left*right/gcd(left,right));
}
poly transformationDenominator=1;
poly cursor;
number coefficient;
int transformationRow;
int transformationColumn;
for(
  transformationRow=1;
  transformationRow<=nrows(transformation);
  transformationRow++
)
{
  for(
    transformationColumn=1;
    transformationColumn<=ncols(transformation);
    transformationColumn++
  )
  {
    cursor=transformation[transformationRow,transformationColumn];
    while(cursor!=0)
    {
      coefficient=leadcoef(cursor);
      transformationDenominator=polynomialLcm(
        transformationDenominator,
        denominator(coefficient)
      );
      cursor=cursor-lead(cursor);
    }
  }
}
print(
  "QRAW_EXTENSION_DENOMINATOR "
  +string(transformationDenominator)
);
"""
    else:
        computation = "ideal G=slimgb(I);"
    program = f"""
ring extension=(43,s1,lam,v,T),(s6,s5,s3),dp;
minpoly=({polynomials[0]});
option(redSB);
timer=1;
poly pivot={polynomials[1]};
ideal pivotBasis=std(pivot);
{reduced_declarations}
ideal I=q4,q5,q6,q7,q8;
{computation}
int elapsed=timer;
int isUnit=(reduce(1,G)==0);
print(
  "QRAW_EXTENSION_META "+string(size(G))+" "+string(dim(G))+" "
  +string(isUnit)+" "+string(elapsed)
);
"""
    program = program.replace(
        "ring extension=(43,",
        f"ring extension=({prime},",
        1,
    )
    print(
        "SINGULAR_EXTENSION_START "
        f"prime={prime} lift={int(lift)} input_bytes={len(program)}",
        file=sys.stderr,
        flush=True,
    )
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=False,
            timeout=(solver_seconds or None),
        )
    except subprocess.TimeoutExpired as error:
        return {
            "status": "timeout",
            "seconds": round(time.monotonic() - started, 3),
            "returncode": None,
            "stdout_tail": (
                error.stdout.decode()
                if isinstance(error.stdout, bytes)
                else (error.stdout or "")
            )[-4000:],
            "stderr_tail": (
                error.stderr.decode()
                if isinstance(error.stderr, bytes)
                else (error.stderr or "")
            )[-4000:],
            "program_bytes": len(program),
            "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
        }
    elapsed = time.monotonic() - started
    marker = re.search(
        r"(?m)^QRAW_EXTENSION_META (-?\d+) (-?\d+) ([01]) (\d+)$",
        completed.stdout,
    )
    if completed.returncode != 0:
        status = f"solver-error-{completed.returncode}"
    elif marker is None:
        status = "missing-result-marker"
    elif marker.group(3) == "1":
        status = "unit"
    else:
        status = "nonunit"
    result: dict[str, object] = {
        "status": status,
        "seconds": round(elapsed, 3),
        "returncode": completed.returncode,
        "stdout_tail": completed.stdout[-4000:],
        "stderr_tail": completed.stderr[-30000:],
        "program_bytes": len(program),
        "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
    }
    if marker is not None:
        result.update(
            {
                "basis_size": int(marker.group(1)),
                "dimension": int(marker.group(2)),
                "unit_ideal": marker.group(3) == "1",
                "singular_timer_ticks": int(marker.group(4)),
            }
        )
    denominator = re.search(
        r"(?m)^QRAW_EXTENSION_DENOMINATOR (.*)$",
        completed.stdout,
    )
    if denominator is not None:
        value = denominator.group(1)
        result.update(
            {
                "transformation_denominator": value,
                "transformation_denominator_length": len(value),
                "transformation_denominator_sha256": hashlib.sha256(
                    value.encode()
                ).hexdigest(),
            }
        )
    return result


def run_raw_extension_gcd_singular(
    system: dict[str, object],
    prime: int,
    solver_seconds: int,
) -> dict[str, object]:
    """Test the direct system by two univariate gcds over the residual field."""

    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    polynomials = system["polynomials"]
    if len(polynomials) != 9:
        raise ValueError("extension gcd requires the direct raw system")
    reduced_declarations = "\n".join(
        f"poly q{order}=reduce(({polynomials[index]}),pivotBasis);"
        for order, index in zip(range(4, 9), range(2, 7), strict=True)
    )
    program = f"""
ring extension=({prime},s1,lam,v,T),(s6,s5,s3),dp;
minpoly=({polynomials[0]});
option(redSB);
timer=1;
poly pivot={polynomials[1]};
ideal pivotBasis=std(pivot);
{reduced_declarations}
poly p5=subst(diff(q5,s6),s6,0);
poly z5=subst(q5,s6,0);
poly p6=subst(diff(q6,s6),s6,0);
poly z6=subst(q6,s6,0);
poly p7=subst(diff(q7,s6),s6,0);
poly z7=subst(q7,s6,0);
poly p8=subst(diff(q8,s6),s6,0);
poly z8=subst(q8,s6,0);
poly a4=subst(diff(diff(q4,s6),s6),s6,0)/2;
poly b4=subst(diff(q4,s6),s6,0);
poly c4=subst(q4,s6,0);
poly e4=a4*z5^2-b4*p5*z5+c4*p5^2;
poly e6=p5*z6-p6*z5;
poly e7=p5*z7-p7*z5;
poly e8=p5*z8-p8*z5;
poly exceptionalGcd=gcd(p5,z5);
poly openGcd=gcd(gcd(gcd(e4,e6),e7),e8);
int exceptionalEmpty=(exceptionalGcd!=0 && deg(exceptionalGcd)==0);
int openEmpty=(openGcd!=0 && deg(openGcd)==0);
int isUnit=(exceptionalEmpty && openEmpty);
int elapsed=timer;
print(
  "QRAW_EXTENSION_GCD_META "
  +string(deg(p5))+" "+string(deg(z5))+" "
  +string(deg(e4))+" "+string(deg(e6))+" "
  +string(deg(e7))+" "+string(deg(e8))+" "
  +string(deg(exceptionalGcd))+" "+string(deg(openGcd))+" "
  +string(exceptionalEmpty)+" "+string(openEmpty)+" "
  +string(isUnit)+" "+string(elapsed)
);
"""
    print(
        "SINGULAR_EXTENSION_GCD_START "
        f"prime={prime} input_bytes={len(program)}",
        file=sys.stderr,
        flush=True,
    )
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=False,
            timeout=(solver_seconds or None),
        )
    except subprocess.TimeoutExpired as error:
        return {
            "status": "timeout",
            "seconds": round(time.monotonic() - started, 3),
            "returncode": None,
            "stdout_tail": (
                error.stdout.decode()
                if isinstance(error.stdout, bytes)
                else (error.stdout or "")
            )[-4000:],
            "stderr_tail": (
                error.stderr.decode()
                if isinstance(error.stderr, bytes)
                else (error.stderr or "")
            )[-4000:],
            "program_bytes": len(program),
            "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
        }
    elapsed = time.monotonic() - started
    marker = re.search(
        (
            r"(?m)^QRAW_EXTENSION_GCD_META "
            r"(-?\d+) (-?\d+) (-?\d+) (-?\d+) (-?\d+) (-?\d+) "
            r"(-?\d+) (-?\d+) ([01]) ([01]) ([01]) (\d+)$"
        ),
        completed.stdout,
    )
    if completed.returncode != 0:
        status = f"solver-error-{completed.returncode}"
    elif marker is None:
        status = "missing-result-marker"
    elif marker.group(11) == "1":
        status = "unit"
    else:
        status = "nonunit"
    result: dict[str, object] = {
        "status": status,
        "seconds": round(elapsed, 3),
        "returncode": completed.returncode,
        "stdout_tail": completed.stdout[-4000:],
        "stderr_tail": completed.stderr[-4000:],
        "program_bytes": len(program),
        "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
    }
    if marker is not None:
        result.update(
            {
                "degrees": {
                    "p5": int(marker.group(1)),
                    "z5": int(marker.group(2)),
                    "e4": int(marker.group(3)),
                    "e6": int(marker.group(4)),
                    "e7": int(marker.group(5)),
                    "e8": int(marker.group(6)),
                    "exceptional_gcd": int(marker.group(7)),
                    "open_gcd": int(marker.group(8)),
                },
                "exceptional_empty": marker.group(9) == "1",
                "open_empty": marker.group(10) == "1",
                "unit_ideal": marker.group(11) == "1",
                "singular_timer_ticks": int(marker.group(12)),
            }
        )
    return result


def run_raw_extension_gcd_macaulay2(
    system: dict[str, object],
    arithmetic: QuotientArithmetic,
    pivot: QuotientElement,
    prime: int,
    solver_seconds: int,
) -> dict[str, object]:
    """Run the two univariate gcd tests in a Macaulay2 field tower."""

    macaulay2 = shutil.which("M2")
    if macaulay2 is None:
        raise RuntimeError("Macaulay2 is required")
    polynomials = system["polynomials"]
    if len(polynomials) != 9:
        raise ValueError("Macaulay2 gcd requires the direct raw system")
    declarations = "\n".join(
        f"q{order}S3=({polynomials[index]});"
        for order, index in zip(range(4, 9), range(2, 7), strict=True)
    )
    mapped = "\n".join(
        f"q{order}=phi q{order}S3;"
        for order in range(4, 9)
    )
    pivot_denominator = arithmetic.registry.product(pivot.denominator)
    program = f"""
kk=frac(GF({prime})[s1,lam,v]);
RT=kk[T];
r20=({polynomials[0]});
r20Monic=r20*sub(1/(leadCoefficient r20),RT);
residualAlgebra=RT/ideal(r20Monic);
K=toField residualAlgebra;
pivotNumeratorK=sub(({pivot.numerator}),K);
pivotDenominatorKK=sub(({pivot_denominator}),kk);
s3ValueK=pivotNumeratorK*sub(1/pivotDenominatorKK,K);
S3=K[s6,s5,s3];
{declarations}
S2=K[s6,s5];
s3Value=promote(s3ValueK,S2);
phi=map(S2,S3,{{s6,s5,s3Value}});
{mapped}
U=K[s5];
atS6Zero=map(U,S2,{{0,s5}});
p5=atS6Zero diff(s6,q5);
z5=atS6Zero q5;
p6=atS6Zero diff(s6,q6);
z6=atS6Zero q6;
p7=atS6Zero diff(s6,q7);
z7=atS6Zero q7;
p8=atS6Zero diff(s6,q8);
z8=atS6Zero q8;
a4=(atS6Zero diff(s6,diff(s6,q4)))/2;
b4=atS6Zero diff(s6,q4);
c4=atS6Zero q4;
e4=a4*z5^2-b4*p5*z5+c4*p5^2;
e6=p5*z6-p6*z5;
e7=p5*z7-p7*z5;
e8=p5*z8-p8*z5;
exceptionalGcd=gcd(p5,z5);
openGcd=gcd(gcd(gcd(e4,e6),e7),e8);
exceptionalEmpty=(first degree exceptionalGcd)==0;
openEmpty=(first degree openGcd)==0;
isUnit=exceptionalEmpty and openEmpty;
print(
  "QRAW_M2_GCD_META "
  |toString(first degree p5)|" "
  |toString(first degree z5)|" "
  |toString(first degree e4)|" "
  |toString(first degree e6)|" "
  |toString(first degree e7)|" "
  |toString(first degree e8)|" "
  |toString(first degree exceptionalGcd)|" "
  |toString(first degree openGcd)|" "
  |toString exceptionalEmpty|" "
  |toString openEmpty|" "
  |toString isUnit
  );
"""
    print(
        "M2_EXTENSION_GCD_START "
        f"prime={prime} input_bytes={len(program)}",
        file=sys.stderr,
        flush=True,
    )
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [macaulay2, "--silent", "--no-readline"],
            input=program,
            text=True,
            capture_output=True,
            check=False,
            timeout=(solver_seconds or None),
        )
    except subprocess.TimeoutExpired as error:
        return {
            "status": "timeout",
            "seconds": round(time.monotonic() - started, 3),
            "returncode": None,
            "stdout_tail": (
                error.stdout.decode()
                if isinstance(error.stdout, bytes)
                else (error.stdout or "")
            )[-4000:],
            "stderr_tail": (
                error.stderr.decode()
                if isinstance(error.stderr, bytes)
                else (error.stderr or "")
            )[-4000:],
            "program_bytes": len(program),
            "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
        }
    elapsed = time.monotonic() - started
    marker = re.search(
        (
            r"QRAW_M2_GCD_META "
            r"(-?\d+) (-?\d+) (-?\d+) (-?\d+) (-?\d+) (-?\d+) "
            r"(-?\d+) (-?\d+) (true|false) (true|false) (true|false)"
        ),
        completed.stdout,
    )
    modulus_specialization = irreducible_modulus_specialization(
        arithmetic,
        prime,
    )
    if completed.returncode != 0 or "error:" in completed.stderr:
        status = f"solver-error-{completed.returncode}"
    elif marker is None:
        status = "missing-result-marker"
    elif (
        marker.group(11) == "true"
        and modulus_specialization["irreducible"]
    ):
        status = "unit"
    else:
        status = "nonunit-or-unverified"
    result: dict[str, object] = {
        "status": status,
        "seconds": round(elapsed, 3),
        "returncode": completed.returncode,
        "stdout_tail": completed.stdout[-4000:],
        "stderr_tail": completed.stderr[-30000:],
        "program_bytes": len(program),
        "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
        "modulus_specialization": modulus_specialization,
    }
    if marker is not None:
        result.update(
            {
                "degrees": {
                    "p5": int(marker.group(1)),
                    "z5": int(marker.group(2)),
                    "e4": int(marker.group(3)),
                    "e6": int(marker.group(4)),
                    "e7": int(marker.group(5)),
                    "e8": int(marker.group(6)),
                    "exceptional_gcd": int(marker.group(7)),
                    "open_gcd": int(marker.group(8)),
                },
                "exceptional_empty": marker.group(9) == "true",
                "open_empty": marker.group(10) == "true",
                "unit_ideal_on_generic_base": (
                    marker.group(11) == "true"
                    and modulus_specialization["irreducible"]
                ),
            }
        )
    return result


def run_raw_extension_m2_split(
    system: dict[str, object],
    arithmetic: QuotientArithmetic,
    pivot: QuotientElement,
    prime: int,
    solver_seconds: int,
    mode: str,
) -> dict[str, object]:
    """Run one minimal Macaulay2 gcd branch over the residual field."""

    macaulay2 = shutil.which("M2")
    if macaulay2 is None:
        raise RuntimeError("Macaulay2 is required")
    polynomials = system["polynomials"]
    if len(polynomials) != 9:
        raise ValueError("Macaulay2 split gcd requires the direct system")
    if mode == "exceptional":
        orders = (5,)
        calculation = """
p5=atS6Zero diff(s6,q5);
z5=atS6Zero q5;
branchGcd=gcd(p5,z5);
print(
  "QRAW_M2_SPLIT_META exceptional "
  |toString(first degree p5)|" "
  |toString(first degree z5)|" "
  |toString(first degree branchGcd)
  );
"""
    elif mode == "open6":
        orders = (4, 5, 6)
        calculation = """
p5=atS6Zero diff(s6,q5);
z5=atS6Zero q5;
p6=atS6Zero diff(s6,q6);
z6=atS6Zero q6;
a4=(atS6Zero diff(s6,diff(s6,q4)))/2;
b4=atS6Zero diff(s6,q4);
c4=atS6Zero q4;
e4=a4*z5^2-b4*p5*z5+c4*p5^2;
e6=p5*z6-p6*z5;
branchGcd=gcd(e4,e6);
print(
  "QRAW_M2_SPLIT_META open6 "
  |toString(first degree e4)|" "
  |toString(first degree e6)|" "
  |toString(first degree branchGcd)
  );
"""
    else:
        raise ValueError(f"unknown Macaulay2 split mode {mode}")
    declarations = "\n".join(
        f"q{order}S3=({polynomials[order - 2]});"
        for order in orders
    )
    mapped = "\n".join(
        f"q{order}=phi q{order}S3;"
        for order in orders
    )
    pivot_denominator = arithmetic.registry.product(pivot.denominator)
    program = f"""
kk=frac(GF({prime})[s1,lam,v]);
RT=kk[T];
r20=({polynomials[0]});
r20Monic=r20*sub(1/(leadCoefficient r20),RT);
residualAlgebra=RT/ideal(r20Monic);
K=toField residualAlgebra;
pivotNumeratorK=sub(({pivot.numerator}),K);
pivotDenominatorKK=sub(({pivot_denominator}),kk);
s3ValueK=pivotNumeratorK*sub(1/pivotDenominatorKK,K);
S3=K[s6,s5,s3];
{declarations}
S2=K[s6,s5];
s3Value=promote(s3ValueK,S2);
phi=map(S2,S3,{{s6,s5,s3Value}});
{mapped}
U=K[s5];
atS6Zero=map(U,S2,{{0,s5}});
{calculation}
"""
    print(
        "M2_SPLIT_GCD_START "
        f"mode={mode} prime={prime} input_bytes={len(program)}",
        file=sys.stderr,
        flush=True,
    )
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [macaulay2, "--silent", "--no-readline"],
            input=program,
            text=True,
            capture_output=True,
            check=False,
            timeout=(solver_seconds or None),
        )
    except subprocess.TimeoutExpired as error:
        return {
            "status": "timeout",
            "mode": mode,
            "seconds": round(time.monotonic() - started, 3),
            "returncode": None,
            "stdout_tail": (
                error.stdout.decode()
                if isinstance(error.stdout, bytes)
                else (error.stdout or "")
            )[-4000:],
            "stderr_tail": (
                error.stderr.decode()
                if isinstance(error.stderr, bytes)
                else (error.stderr or "")
            )[-4000:],
            "program_bytes": len(program),
            "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
        }
    elapsed = time.monotonic() - started
    marker = re.search(
        rf"QRAW_M2_SPLIT_META {mode} (-?\d+) (-?\d+) (-?\d+)",
        completed.stdout,
    )
    modulus_specialization = irreducible_modulus_specialization(
        arithmetic,
        prime,
    )
    if completed.returncode != 0 or "error:" in completed.stderr:
        status = f"solver-error-{completed.returncode}"
    elif marker is None:
        status = "missing-result-marker"
    elif (
        int(marker.group(3)) == 0
        and modulus_specialization["irreducible"]
    ):
        status = "coprime"
    else:
        status = "common-factor-or-unverified"
    result: dict[str, object] = {
        "status": status,
        "mode": mode,
        "seconds": round(elapsed, 3),
        "returncode": completed.returncode,
        "stdout_tail": completed.stdout[-4000:],
        "stderr_tail": completed.stderr[-12000:],
        "program_bytes": len(program),
        "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
        "modulus_specialization": modulus_specialization,
    }
    if marker is not None:
        result.update(
            {
                "input_degrees": [
                    int(marker.group(1)),
                    int(marker.group(2)),
                ],
                "gcd_degree": int(marker.group(3)),
                "coprime_on_generic_base": (
                    int(marker.group(3)) == 0
                    and modulus_specialization["irreducible"]
                ),
            }
        )
    return result


def clear_fibre_denominators(
    arithmetic: QuotientArithmetic,
    polynomial: FibrePolynomial,
) -> FibrePolynomial:
    """Multiply a fibre polynomial by one common nonzero base denominator."""

    return clear_quotient_polynomial_denominators(arithmetic, polynomial)


def clear_quotient_polynomial_denominators(
    arithmetic: QuotientArithmetic,
    polynomial: dict,
) -> dict:
    """Clear the recorded base denominators from a sparse polynomial."""

    if not polynomial:
        return {}
    denominators = [
        arithmetic.registry.pad(coefficient.denominator)
        for coefficient in polynomial.values()
    ]
    common = tuple(
        max(denominator[index] for denominator in denominators)
        for index in range(len(arithmetic.registry.factors))
    )
    answer: FibrePolynomial = {}
    for monomial, coefficient in polynomial.items():
        denominator = arithmetic.registry.pad(coefficient.denominator)
        missing = tuple(
            common[index] - denominator[index]
            for index in range(len(common))
        )
        numerator = coefficient.numerator * arithmetic.registry.product(
            missing
        )
        value = arithmetic.make(numerator)
        if not value.numerator.is_zero():
            answer[monomial] = value
    return answer


def quotient_raw_polynomial(
    arithmetic: QuotientArithmetic,
    polynomial: RawFibrePolynomial,
) -> RawQuotientPolynomial:
    """Reduce raw base coefficients modulo the residual quintic."""

    return {
        monomial: arithmetic.make(coefficient)
        for monomial, coefficient in polynomial.items()
    }


def raw_x_coefficient(
    polynomial: RawQuotientPolynomial,
    degree: int,
) -> YZPolynomial:
    return {
        (s5_degree, s3_degree): coefficient
        for (s6_degree, s5_degree, s3_degree), coefficient
        in polynomial.items()
        if s6_degree == degree
    }


def add_yz_polynomials(
    arithmetic: QuotientArithmetic,
    left: YZPolynomial,
    right: YZPolynomial,
    right_sign: int = 1,
) -> YZPolynomial:
    answer = dict(left)
    for monomial, coefficient in right.items():
        if right_sign < 0:
            coefficient = arithmetic.neg(coefficient)
        answer[monomial] = arithmetic.add(
            answer.get(monomial, arithmetic.zero),
            coefficient,
        )
    return {
        monomial: coefficient
        for monomial, coefficient in answer.items()
        if not coefficient.numerator.is_zero()
    }


def multiply_yz_polynomials(
    arithmetic: QuotientArithmetic,
    left: YZPolynomial,
    right: YZPolynomial,
) -> YZPolynomial:
    answer: YZPolynomial = {}
    for (left_y, left_z), left_coefficient in left.items():
        for (right_y, right_z), right_coefficient in right.items():
            monomial = (left_y + right_y, left_z + right_z)
            product = arithmetic.mul(left_coefficient, right_coefficient)
            answer[monomial] = arithmetic.add(
                answer.get(monomial, arithmetic.zero),
                product,
            )
    return clear_quotient_polynomial_denominators(arithmetic, answer)


def scale_yz_polynomial(
    arithmetic: QuotientArithmetic,
    polynomial: YZPolynomial,
    scalar: QuotientElement,
) -> YZPolynomial:
    return {
        monomial: arithmetic.mul(coefficient, scalar)
        for monomial, coefficient in polynomial.items()
    }


def raw_quadratic_pseudo_remainder(
    arithmetic: QuotientArithmetic,
    divisor: RawFibrePolynomial,
    dividend: RawFibrePolynomial,
) -> RawQuotientPolynomial:
    """Pseudo-divide in ``s6`` before substituting the dense ``s3`` pivot."""

    divisor_q = quotient_raw_polynomial(arithmetic, divisor)
    dividend_q = quotient_raw_polynomial(arithmetic, dividend)
    divisor_by_x = {
        degree: raw_x_coefficient(divisor_q, degree)
        for degree in range(3)
    }
    if set(divisor_by_x[2]) != {(0, 0)}:
        raise AssertionError(
            "the raw quadratic leading coefficient depends on s5 or s3"
        )
    leading = divisor_by_x[2][(0, 0)]
    remainder_by_x = {
        degree: raw_x_coefficient(dividend_q, degree)
        for degree in range(
            max((monomial[0] for monomial in dividend_q), default=-1) + 1
        )
    }
    remainder_by_x = {
        degree: polynomial
        for degree, polynomial in remainder_by_x.items()
        if polynomial
    }
    while remainder_by_x and max(remainder_by_x) >= 2:
        top_degree = max(remainder_by_x)
        top = remainder_by_x[top_degree]
        shift = top_degree - 2
        scaled = {
            degree: scale_yz_polynomial(
                arithmetic,
                polynomial,
                leading,
            )
            for degree, polynomial in remainder_by_x.items()
        }
        for divisor_degree, divisor_coefficient in divisor_by_x.items():
            target_degree = divisor_degree + shift
            product = multiply_yz_polynomials(
                arithmetic,
                top,
                divisor_coefficient,
            )
            scaled[target_degree] = add_yz_polynomials(
                arithmetic,
                scaled.get(target_degree, {}),
                product,
                -1,
            )
        remainder_by_x = {
            degree: clear_quotient_polynomial_denominators(
                arithmetic,
                polynomial,
            )
            for degree, polynomial in scaled.items()
            if polynomial
        }
        remainder_by_x = {
            degree: polynomial
            for degree, polynomial in remainder_by_x.items()
            if polynomial
        }
    answer = {
        (s6_degree, s5_degree, s3_degree): coefficient
        for s6_degree, polynomial in remainder_by_x.items()
        for (s5_degree, s3_degree), coefficient in polynomial.items()
    }
    return clear_quotient_polynomial_denominators(arithmetic, answer)


def evaluate_raw_quotient_fraction_free(
    arithmetic: QuotientArithmetic,
    raw: RawQuotientPolynomial,
    a_value: fmpq_mpoly,
    b_value: fmpq_mpoly,
) -> FibrePolynomial:
    """Substitute the pivot only after a raw elimination has been performed."""

    degree = max((monomial[2] for monomial in raw), default=0)
    numerator = arithmetic.make(-b_value)
    denominator = arithmetic.scalar_mul(
        arithmetic.make(a_value),
        CTX.constant(6),
    )
    numerator_powers = [arithmetic.one]
    denominator_powers = [arithmetic.one]
    for _ in range(degree):
        numerator_powers.append(
            arithmetic.mul(numerator_powers[-1], numerator)
        )
        denominator_powers.append(
            arithmetic.mul(denominator_powers[-1], denominator)
        )
    answer: FibrePolynomial = {}
    for (s6_power, s5_power, s3_power), coefficient in raw.items():
        coefficient = arithmetic.mul(
            coefficient,
            numerator_powers[s3_power],
        )
        coefficient = arithmetic.mul(
            coefficient,
            denominator_powers[degree - s3_power],
        )
        monomial = (s6_power, s5_power)
        answer[monomial] = arithmetic.add(
            answer.get(monomial, arithmetic.zero),
            coefficient,
        )
    answer = {
        monomial: coefficient
        for monomial, coefficient in answer.items()
        if not coefficient.numerator.is_zero()
    }
    return clear_fibre_denominators(arithmetic, answer)


def raw_quotient_profile(
    polynomial: RawQuotientPolynomial,
) -> dict[str, object]:
    return {
        "terms": len(polynomial),
        "degrees": {
            "s6": max((monomial[0] for monomial in polynomial), default=-1),
            "s5": max((monomial[1] for monomial in polynomial), default=-1),
            "s3": max((monomial[2] for monomial in polynomial), default=-1),
        },
        "maximum_coefficient_terms": max(
            (
                len(coefficient.numerator.to_dict())
                for coefficient in polynomial.values()
            ),
            default=0,
        ),
    }


def yz_polynomial_profile(
    polynomial: YZPolynomial,
) -> dict[str, object]:
    return {
        "terms": len(polynomial),
        "degrees": {
            "s5": max((monomial[0] for monomial in polynomial), default=-1),
            "s3": max((monomial[1] for monomial in polynomial), default=-1),
        },
        "maximum_coefficient_terms": max(
            (
                len(coefficient.numerator.to_dict())
                for coefficient in polynomial.values()
            ),
            default=0,
        ),
    }


def raw_remainder_univariate_equations(
    arithmetic: QuotientArithmetic,
    mu4: RawFibrePolynomial,
    remainders: dict[int, RawQuotientPolynomial],
) -> dict[str, YZPolynomial]:
    """Eliminate ``s6`` before expanding the dense pivot substitution."""

    if set(remainders) != {5, 6, 7}:
        raise ValueError("the raw equations require remainders 5, 6, 7")
    mu4_q = quotient_raw_polynomial(arithmetic, mu4)
    p = {
        order: raw_x_coefficient(polynomial, 1)
        for order, polynomial in remainders.items()
    }
    q = {
        order: raw_x_coefficient(polynomial, 0)
        for order, polynomial in remainders.items()
    }
    a = raw_x_coefficient(mu4_q, 2)
    f1 = raw_x_coefficient(mu4_q, 1)
    f0 = raw_x_coefficient(mu4_q, 0)
    if set(a) != {(0, 0)}:
        raise AssertionError(
            "the raw quadratic leading coefficient depends on s5 or s3"
        )
    quadratic = add_yz_polynomials(
        arithmetic,
        scale_yz_polynomial(
            arithmetic,
            multiply_yz_polynomials(arithmetic, q[5], q[5]),
            a[(0, 0)],
        ),
        multiply_yz_polynomials(
            arithmetic,
            f1,
            multiply_yz_polynomials(arithmetic, p[5], q[5]),
        ),
        -1,
    )
    quadratic = add_yz_polynomials(
        arithmetic,
        quadratic,
        multiply_yz_polynomials(
            arithmetic,
            f0,
            multiply_yz_polynomials(arithmetic, p[5], p[5]),
        ),
    )
    equations = {
        "quadratic": clear_quotient_polynomial_denominators(
            arithmetic,
            quadratic,
        )
    }
    for order in (6, 7):
        cross = add_yz_polynomials(
            arithmetic,
            multiply_yz_polynomials(arithmetic, p[5], q[order]),
            multiply_yz_polynomials(arithmetic, p[order], q[5]),
            -1,
        )
        equations[f"cross_{order}"] = (
            clear_quotient_polynomial_denominators(
                arithmetic,
                cross,
            )
        )
    return equations


def evaluate_yz_fraction_free(
    arithmetic: QuotientArithmetic,
    polynomial: YZPolynomial,
    a_value: fmpq_mpoly,
    b_value: fmpq_mpoly,
) -> YPolynomial:
    evaluated = evaluate_raw_quotient_fraction_free(
        arithmetic,
        {
            (0, s5_degree, s3_degree): coefficient
            for (s5_degree, s3_degree), coefficient in polynomial.items()
        },
        a_value,
        b_value,
    )
    return {
        s5_degree: coefficient
        for (s6_degree, s5_degree), coefficient in evaluated.items()
        if s6_degree == 0
    }


def raw_cross_equation(
    arithmetic: QuotientArithmetic,
    reference: RawQuotientPolynomial,
    other: RawQuotientPolynomial,
) -> YZPolynomial:
    """Eliminate ``s6`` between two raw linear remainders."""

    reference_p = raw_x_coefficient(reference, 1)
    reference_q = raw_x_coefficient(reference, 0)
    other_p = raw_x_coefficient(other, 1)
    other_q = raw_x_coefficient(other, 0)
    cross = add_yz_polynomials(
        arithmetic,
        multiply_yz_polynomials(arithmetic, reference_p, other_q),
        multiply_yz_polynomials(arithmetic, other_p, reference_q),
        -1,
    )
    return clear_quotient_polynomial_denominators(arithmetic, cross)


def serialize_yz_quotient_polynomial(
    arithmetic: QuotientArithmetic,
    polynomial: YZPolynomial,
) -> str:
    """Serialize a denominator-free raw equation for Singular or msolve."""

    pieces = []
    for (s5_degree, s3_degree), coefficient in sorted(polynomial.items()):
        denominator = arithmetic.registry.pad(coefficient.denominator)
        if any(denominator):
            raise AssertionError(
                "raw system serialization retained a base denominator"
            )
        factors = [f"({coefficient.numerator})"]
        if s5_degree:
            factors.append(
                "s5" if s5_degree == 1 else f"s5^{s5_degree}"
            )
        if s3_degree:
            factors.append(
                "s3" if s3_degree == 1 else f"s3^{s3_degree}"
            )
        pieces.append("*".join(factors))
    return "+".join(pieces) or "0"


def serialize_raw_quotient_polynomial(
    arithmetic: QuotientArithmetic,
    polynomial: RawQuotientPolynomial,
) -> str:
    """Serialize a denominator-free equation in ``s6,s5,s3``."""

    polynomial = clear_quotient_polynomial_denominators(
        arithmetic,
        polynomial,
    )
    pieces = []
    for monomial, coefficient in sorted(polynomial.items()):
        denominator = arithmetic.registry.pad(coefficient.denominator)
        if any(denominator):
            raise AssertionError(
                "raw system serialization retained a base denominator"
            )
        factors = [f"({coefficient.numerator})"]
        for variable, exponent in zip(
            ("s6", "s5", "s3"),
            monomial,
            strict=True,
        ):
            if exponent == 1:
                factors.append(variable)
            elif exponent:
                factors.append(f"{variable}^{exponent}")
        pieces.append("*".join(factors))
    return "+".join(pieces) or "0"


def raw_mu8_system(
    arithmetic: QuotientArithmetic,
    a_value: fmpq_mpoly,
    b_value: fmpq_mpoly,
    moments: dict[int, RawFibrePolynomial],
    p5_open: bool = False,
) -> tuple[dict[str, object], dict[int, RawQuotientPolynomial], dict[str, YZPolynomial]]:
    """Return the compact residual system through ``mu8`` before the pivot."""

    if not set(range(4, 9)) <= set(moments):
        raise ValueError("raw mu8 system requires moments mu4 through mu8")
    remainders = {
        order: raw_quadratic_pseudo_remainder(
            arithmetic,
            moments[4],
            moments[order],
        )
        for order in range(5, 9)
    }
    equations = raw_remainder_univariate_equations(
        arithmetic,
        moments[4],
        {order: remainders[order] for order in (5, 6, 7)},
    )
    equations["cross_8"] = raw_cross_equation(
        arithmetic,
        remainders[5],
        remainders[8],
    )
    equations["p5"] = clear_quotient_polynomial_denominators(
        arithmetic,
        raw_x_coefficient(remainders[5], 1),
    )
    polynomials = [
        str(arithmetic.modulus),
        f"(6)*({a_value})*s3+({b_value})",
        *[
            serialize_yz_quotient_polynomial(
                arithmetic,
                equations[name],
            )
            for name in ("quadratic", "cross_6", "cross_7", "cross_8")
        ],
        "v*vinv-1",
        "(6084*lam^2+4805)*jinv-1",
    ]
    ordinary_variables = [
        "s5",
        "s3",
        "T",
        "s1",
        "lam",
        "v",
        "vinv",
        "jinv",
    ]
    if p5_open:
        polynomials.append(
            "("
            + serialize_yz_quotient_polynomial(
                arithmetic,
                equations["p5"],
            )
            + ")*pinv-1"
        )
        ordinary_variables.append("pinv")
    system = {
        "ordinary_variables": ordinary_variables,
        "polynomials": polynomials,
        "equation_profiles": [
            {
                "index": index,
                "length": len(polynomial),
                "sha256": hashlib.sha256(polynomial.encode()).hexdigest(),
            }
            for index, polynomial in enumerate(polynomials)
        ],
        "scope": (
            "R20=0, dense pivot, the mu5 quadratic resultant, "
            "mu6/mu7/mu8 cross equations, v!=0, and J_Q!=0; "
            + (
                "P5!=0 (hence the common linear remainder has a "
                "well-defined root)"
                if p5_open
                else (
                    "the exceptional case P5=Q5=0 is retained by the "
                    "quadratic equation"
                )
            )
        ),
    }
    return system, remainders, equations


def raw_mu8_exceptional_system(
    arithmetic: QuotientArithmetic,
    a_value: fmpq_mpoly,
    b_value: fmpq_mpoly,
    moments: dict[int, RawFibrePolynomial],
) -> tuple[
    dict[str, object],
    dict[int, RawQuotientPolynomial],
    dict[str, YZPolynomial],
]:
    """Return the complementary ``P5=Q5=0`` residual system."""

    if not set(range(4, 9)) <= set(moments):
        raise ValueError("raw mu8 system requires moments mu4 through mu8")
    reduced = {
        4: clear_quotient_polynomial_denominators(
            arithmetic,
            quotient_raw_polynomial(arithmetic, moments[4]),
        )
    }
    reduced.update(
        {
            order: raw_quadratic_pseudo_remainder(
                arithmetic,
                moments[4],
                moments[order],
            )
            for order in range(5, 9)
        }
    )
    equations = {
        "p5": clear_quotient_polynomial_denominators(
            arithmetic,
            raw_x_coefficient(reduced[5], 1),
        ),
        "q5": clear_quotient_polynomial_denominators(
            arithmetic,
            raw_x_coefficient(reduced[5], 0),
        ),
    }
    polynomials = [
        str(arithmetic.modulus),
        f"(6)*({a_value})*s3+({b_value})",
        serialize_yz_quotient_polynomial(arithmetic, equations["p5"]),
        serialize_yz_quotient_polynomial(arithmetic, equations["q5"]),
        *[
            serialize_raw_quotient_polynomial(
                arithmetic,
                reduced[order],
            )
            for order in (4, 6, 7, 8)
        ],
        "v*vinv-1",
        "(6084*lam^2+4805)*jinv-1",
    ]
    system = {
        "ordinary_variables": [
            "s6",
            "s5",
            "s3",
            "T",
            "s1",
            "lam",
            "v",
            "vinv",
            "jinv",
        ],
        "polynomials": polynomials,
        "equation_profiles": [
            {
                "index": index,
                "length": len(polynomial),
                "sha256": hashlib.sha256(polynomial.encode()).hexdigest(),
            }
            for index, polynomial in enumerate(polynomials)
        ],
        "scope": (
            "R20=0, dense pivot, P5=Q5=0, mu4=0, and the "
            "mu6/mu7/mu8 pseudo-remainders at the common s6 root, "
            "with v!=0 and J_Q!=0"
        ),
    }
    return system, reduced, equations


def raw_mu8_direct_system(
    arithmetic: QuotientArithmetic,
    a_value: fmpq_mpoly,
    b_value: fmpq_mpoly,
    moments: dict[int, RawFibrePolynomial],
    p5_open: bool = False,
) -> tuple[
    dict[str, object],
    dict[int, RawQuotientPolynomial],
    dict[str, YZPolynomial],
]:
    """Keep ``s6`` and impose all compact pseudo-remainders directly."""

    if not set(range(4, 9)) <= set(moments):
        raise ValueError("raw mu8 system requires moments mu4 through mu8")
    reduced = {
        4: clear_quotient_polynomial_denominators(
            arithmetic,
            quotient_raw_polynomial(arithmetic, moments[4]),
        )
    }
    reduced.update(
        {
            order: raw_quadratic_pseudo_remainder(
                arithmetic,
                moments[4],
                moments[order],
            )
            for order in range(5, 9)
        }
    )
    polynomials = [
        str(arithmetic.modulus),
        f"(6)*({a_value})*s3+({b_value})",
        *[
            serialize_raw_quotient_polynomial(
                arithmetic,
                reduced[order],
            )
            for order in range(4, 9)
        ],
        "v*vinv-1",
        "(6084*lam^2+4805)*jinv-1",
    ]
    ordinary_variables = [
        "s6",
        "s5",
        "s3",
        "T",
        "s1",
        "lam",
        "v",
        "vinv",
        "jinv",
    ]
    if p5_open:
        p5 = clear_quotient_polynomial_denominators(
            arithmetic,
            raw_x_coefficient(reduced[5], 1),
        )
        polynomials.append(
            "("
            + serialize_yz_quotient_polynomial(arithmetic, p5)
            + ")*pinv-1"
        )
        ordinary_variables.append("pinv")
    system = {
        "ordinary_variables": ordinary_variables,
        "polynomials": polynomials,
        "equation_profiles": [
            {
                "index": index,
                "length": len(polynomial),
                "sha256": hashlib.sha256(polynomial.encode()).hexdigest(),
            }
            for index, polynomial in enumerate(polynomials)
        ],
        "scope": (
            "R20=0, dense pivot, mu4=0, and the compact "
            "mu5/mu6/mu7/mu8 pseudo-remainders at s6, with v!=0 "
            "and J_Q!=0. Every genuine common moment zero lies in this "
            "system; a vanishing mu4 leading coefficient can only add "
            "extraneous points."
            + (
                " The branch is localized at P5!=0 without expanding "
                "the cross-resultants."
                if p5_open
                else ""
            )
        ),
    }
    return system, reduced, {}


def evaluate_moment_fraction_free(
    arithmetic: QuotientArithmetic,
    raw: RawFibrePolynomial,
    a_value: fmpq_mpoly,
    b_value: fmpq_mpoly,
) -> FibrePolynomial:
    """Substitute ``s3=-B/(6A)`` after clearing one global power of ``6A``."""

    degree = max((monomial[2] for monomial in raw), default=0)
    numerator = arithmetic.make(-b_value)
    denominator = arithmetic.scalar_mul(
        arithmetic.make(a_value),
        CTX.constant(6),
    )
    numerator_powers = [arithmetic.one]
    denominator_powers = [arithmetic.one]
    for _ in range(degree):
        numerator_powers.append(
            arithmetic.mul(numerator_powers[-1], numerator)
        )
        denominator_powers.append(
            arithmetic.mul(denominator_powers[-1], denominator)
        )
    answer: FibrePolynomial = {}
    for (s6_power, s5_power, s3_power), base_coefficient in raw.items():
        coefficient = arithmetic.make(base_coefficient)
        coefficient = arithmetic.mul(
            coefficient,
            numerator_powers[s3_power],
        )
        coefficient = arithmetic.mul(
            coefficient,
            denominator_powers[degree - s3_power],
        )
        monomial = (s6_power, s5_power)
        answer[monomial] = arithmetic.add(
            answer.get(monomial, arithmetic.zero),
            coefficient,
        )
    answer = {
        monomial: coefficient
        for monomial, coefficient in answer.items()
        if not coefficient.numerator.is_zero()
    }
    return clear_fibre_denominators(arithmetic, answer)


def evaluate_quotient_polynomial_fraction_free(
    arithmetic: QuotientArithmetic,
    raw: RawQuotientPolynomial,
    a_value: fmpq_mpoly,
    b_value: fmpq_mpoly,
) -> FibrePolynomial:
    """Substitute the dense ``s3`` pivot in a pre-reduced raw polynomial."""

    degree = max((monomial[2] for monomial in raw), default=0)
    numerator = arithmetic.make(-b_value)
    denominator = arithmetic.scalar_mul(
        arithmetic.make(a_value),
        CTX.constant(6),
    )
    numerator_powers = [arithmetic.one]
    denominator_powers = [arithmetic.one]
    for _ in range(degree):
        numerator_powers.append(
            arithmetic.mul(numerator_powers[-1], numerator)
        )
        denominator_powers.append(
            arithmetic.mul(denominator_powers[-1], denominator)
        )
    answer: FibrePolynomial = {}
    for (s6_power, s5_power, s3_power), raw_coefficient in raw.items():
        coefficient = arithmetic.mul(
            raw_coefficient,
            numerator_powers[s3_power],
        )
        coefficient = arithmetic.mul(
            coefficient,
            denominator_powers[degree - s3_power],
        )
        monomial = (s6_power, s5_power)
        answer[monomial] = arithmetic.add(
            answer.get(monomial, arithmetic.zero),
            coefficient,
        )
    answer = {
        monomial: coefficient
        for monomial, coefficient in answer.items()
        if not coefficient.numerator.is_zero()
    }
    return clear_fibre_denominators(arithmetic, answer)


def quotient_profile(
    arithmetic: QuotientArithmetic,
    value: QuotientElement,
) -> dict[str, object]:
    return {
        "numerator": polynomial_profile(value.numerator),
        "denominator_exponents": list(
            arithmetic.registry.pad(value.denominator)
        ),
    }


def fibre_profile(
    arithmetic: QuotientArithmetic,
    polynomial: FibrePolynomial,
) -> dict[str, object]:
    leading = max(polynomial, key=monomial_key)
    return {
        "fibre_terms": len(polynomial),
        "leading_monomial": list(leading),
        "coefficient_numerator_terms": {
            f"{monomial[0]},{monomial[1]}": len(
                coefficient.numerator.to_dict()
            )
            for monomial, coefficient in sorted(polynomial.items())
        },
        "maximum_coefficient_terms": max(
            len(coefficient.numerator.to_dict())
            for coefficient in polynomial.values()
        ),
        "leading_coefficient": quotient_profile(
            arithmetic,
            polynomial[leading],
        ),
    }


def primitive_fibre_polynomial(
    polynomial: FibrePolynomial,
) -> tuple[FibrePolynomial, fmpq_mpoly]:
    """Remove the exact polynomial gcd of the displayed coefficients."""

    coefficients = [
        coefficient.numerator for coefficient in polynomial.values()
    ]
    content = coefficients[0]
    for coefficient in coefficients[1:]:
        content = content.gcd(coefficient)
        if content.is_constant():
            break
    if content.is_constant():
        return polynomial, CTX.constant(1)
    primitive = {
        monomial: QuotientElement(coefficient.numerator // content, ())
        for monomial, coefficient in polynomial.items()
    }
    return primitive, content


def monomial_key(monomial: tuple[int, int]) -> tuple[int, int]:
    """Degree-reverse-lex order with s6 greater than s5."""

    return sum(monomial), -monomial[1]


def leading_term_fibre(
    polynomial: FibrePolynomial,
) -> tuple[tuple[int, int], QuotientElement]:
    monomial = max(polynomial, key=monomial_key)
    return monomial, polynomial[monomial]


def scale_shift_fibre(
    arithmetic: QuotientArithmetic,
    polynomial: FibrePolynomial,
    scalar: QuotientElement,
    shift: tuple[int, int] = (0, 0),
) -> FibrePolynomial:
    return {
        (monomial[0] + shift[0], monomial[1] + shift[1]): (
            arithmetic.mul(coefficient, scalar)
        )
        for monomial, coefficient in polynomial.items()
        if not coefficient.numerator.is_zero()
    }


def add_fibre(
    arithmetic: QuotientArithmetic,
    left: FibrePolynomial,
    right: FibrePolynomial,
    right_sign: int = 1,
) -> FibrePolynomial:
    answer = dict(left)
    for monomial, coefficient in right.items():
        if right_sign < 0:
            coefficient = arithmetic.neg(coefficient)
        answer[monomial] = arithmetic.add(
            answer.get(monomial, arithmetic.zero),
            coefficient,
        )
    return {
        monomial: coefficient
        for monomial, coefficient in answer.items()
        if not coefficient.numerator.is_zero()
    }


def normalize_fibre_generic(
    arithmetic: QuotientArithmetic,
    polynomial: FibrePolynomial,
    removed_contents: list[dict[str, object]],
    reason: str,
) -> FibrePolynomial:
    if not polynomial:
        return {}
    polynomial = clear_fibre_denominators(arithmetic, polynomial)
    if PRIME:
        return polynomial
    polynomial, content = primitive_fibre_polynomial(polynomial)
    if not content.is_constant():
        removed_contents.append(
            {"reason": reason, **polynomial_profile(content)}
        )
    return polynomial


def pseudo_normal_form_fibre(
    arithmetic: QuotientArithmetic,
    polynomial: FibrePolynomial,
    basis: list[FibrePolynomial],
    removed_contents: list[dict[str, object]],
) -> FibrePolynomial:
    """Fraction-free normal form over the generic degree-five field."""

    remainder: FibrePolynomial = {}
    work = polynomial
    reductions = 0
    while work:
        monomial, coefficient = leading_term_fibre(work)
        reducer = None
        for candidate in basis:
            leading, _ = leading_term_fibre(candidate)
            if (
                monomial[0] >= leading[0]
                and monomial[1] >= leading[1]
            ):
                reducer = candidate
                break
        if reducer is None:
            remainder[monomial] = coefficient
            del work[monomial]
            continue
        leading, leading_coefficient = leading_term_fibre(reducer)
        work = scale_shift_fibre(
            arithmetic,
            work,
            leading_coefficient,
        )
        remainder = scale_shift_fibre(
            arithmetic,
            remainder,
            leading_coefficient,
        )
        multiple = scale_shift_fibre(
            arithmetic,
            reducer,
            coefficient,
            (
                monomial[0] - leading[0],
                monomial[1] - leading[1],
            ),
        )
        work = add_fibre(arithmetic, work, multiple, -1)
        reductions += 1
        work = normalize_fibre_generic(
            arithmetic,
            work,
            removed_contents,
            f"pseudo-reduction-{reductions}",
        )
        remainder = normalize_fibre_generic(
            arithmetic,
            remainder,
            removed_contents,
            f"pseudo-remainder-{reductions}",
        )
    return normalize_fibre_generic(
        arithmetic,
        remainder,
        removed_contents,
        "pseudo-normal-form",
    )


def fraction_free_fibre_basis(
    arithmetic: QuotientArithmetic,
    generators: list[FibrePolynomial],
) -> tuple[list[FibrePolynomial], list[dict[str, object]]]:
    basis = list(generators)
    removed_contents: list[dict[str, object]] = []
    pairs = [
        (left, right)
        for left in range(len(basis))
        for right in range(left)
    ]
    while pairs:
        left_index, right_index = pairs.pop(0)
        left = basis[left_index]
        right = basis[right_index]
        left_monomial, left_coefficient = leading_term_fibre(left)
        right_monomial, right_coefficient = leading_term_fibre(right)
        common = (
            max(left_monomial[0], right_monomial[0]),
            max(left_monomial[1], right_monomial[1]),
        )
        s_polynomial = add_fibre(
            arithmetic,
            scale_shift_fibre(
                arithmetic,
                left,
                right_coefficient,
                (
                    common[0] - left_monomial[0],
                    common[1] - left_monomial[1],
                ),
            ),
            scale_shift_fibre(
                arithmetic,
                right,
                left_coefficient,
                (
                    common[0] - right_monomial[0],
                    common[1] - right_monomial[1],
                ),
            ),
            -1,
        )
        s_polynomial = normalize_fibre_generic(
            arithmetic,
            s_polynomial,
            removed_contents,
            "S-polynomial",
        )
        remainder = pseudo_normal_form_fibre(
            arithmetic,
            s_polynomial,
            basis,
            removed_contents,
        )
        if remainder:
            if len(basis) >= 8:
                raise RuntimeError(
                    "unexpected generic border-basis growth: "
                    + str(
                        [
                            leading_term_fibre(polynomial)[0]
                            for polynomial in basis
                        ]
                    )
                )
            new_index = len(basis)
            pairs.extend((new_index, old) for old in range(new_index))
            basis.append(remainder)
    minimal: list[FibrePolynomial] = []
    for polynomial in sorted(
        basis,
        key=lambda value: monomial_key(
            leading_term_fibre(value)[0]
        ),
    ):
        leading, _ = leading_term_fibre(polynomial)
        if any(
            leading[0] >= retained[0] and leading[1] >= retained[1]
            for retained, _ in map(leading_term_fibre, minimal)
        ):
            continue
        minimal.append(polynomial)
    return minimal, removed_contents


YPolynomial = dict[int, QuotientElement]


def x_coefficient(
    polynomial: FibrePolynomial,
    degree: int,
) -> YPolynomial:
    return {
        y_degree: coefficient
        for (x_degree, y_degree), coefficient in polynomial.items()
        if x_degree == degree
    }


def add_y_polynomials(
    arithmetic: QuotientArithmetic,
    left: YPolynomial,
    right: YPolynomial,
    right_sign: int = 1,
) -> YPolynomial:
    answer = dict(left)
    for degree, coefficient in right.items():
        if right_sign < 0:
            coefficient = arithmetic.neg(coefficient)
        answer[degree] = arithmetic.add(
            answer.get(degree, arithmetic.zero),
            coefficient,
        )
    return {
        degree: coefficient
        for degree, coefficient in answer.items()
        if not coefficient.numerator.is_zero()
    }


def multiply_y_polynomials(
    arithmetic: QuotientArithmetic,
    left: YPolynomial,
    right: YPolynomial,
) -> YPolynomial:
    answer: YPolynomial = {}
    for left_degree, left_coefficient in left.items():
        for right_degree, right_coefficient in right.items():
            degree = left_degree + right_degree
            product = arithmetic.mul(
                left_coefficient,
                right_coefficient,
            )
            answer[degree] = arithmetic.add(
                answer.get(degree, arithmetic.zero),
                product,
            )
    fibre = {
        (0, degree): coefficient
        for degree, coefficient in answer.items()
    }
    fibre = clear_fibre_denominators(arithmetic, fibre)
    return {
        degree: coefficient
        for (_, degree), coefficient in fibre.items()
    }


def scale_y_polynomial(
    arithmetic: QuotientArithmetic,
    polynomial: YPolynomial,
    scalar: QuotientElement,
) -> YPolynomial:
    return {
        degree: arithmetic.mul(coefficient, scalar)
        for degree, coefficient in polynomial.items()
    }


def clear_y_denominators(
    arithmetic: QuotientArithmetic,
    polynomial: YPolynomial,
) -> YPolynomial:
    cleared = clear_fibre_denominators(
        arithmetic,
        {
            (0, degree): coefficient
            for degree, coefficient in polynomial.items()
        },
    )
    return {
        degree: coefficient
        for (_, degree), coefficient in cleared.items()
    }


def primitive_y_polynomial(
    arithmetic: QuotientArithmetic,
    polynomial: YPolynomial,
) -> tuple[YPolynomial, fmpq_mpoly]:
    """Clear known denominators and remove displayed polynomial content."""

    cleared = clear_y_denominators(arithmetic, polynomial)
    if not cleared:
        return {}, CTX.constant(1)
    primitive, content = primitive_fibre_polynomial(
        {
            (0, degree): coefficient
            for degree, coefficient in cleared.items()
        }
    )
    return (
        {
            degree: coefficient
            for (_, degree), coefficient in primitive.items()
        },
        content,
    )


def pseudo_remainder_y(
    arithmetic: QuotientArithmetic,
    dividend: YPolynomial,
    divisor: YPolynomial,
) -> tuple[YPolynomial, list[dict[str, object]]]:
    """Return a fraction-free univariate pseudo-remainder in ``s5``."""

    if not divisor:
        raise ZeroDivisionError("zero univariate divisor")
    remainder, initial_content = primitive_y_polynomial(
        arithmetic,
        dividend,
    )
    divisor, divisor_content = primitive_y_polynomial(
        arithmetic,
        divisor,
    )
    contents = [
        polynomial_profile(content)
        for content in (initial_content, divisor_content)
        if not content.is_constant()
    ]
    divisor_degree = max(divisor)
    divisor_leading = divisor[divisor_degree]
    while remainder and max(remainder) >= divisor_degree:
        remainder_degree = max(remainder)
        remainder_leading = remainder[remainder_degree]
        shift = remainder_degree - divisor_degree
        scaled_remainder = scale_y_polynomial(
            arithmetic,
            remainder,
            divisor_leading,
        )
        scaled_divisor = {
            degree + shift: arithmetic.mul(
                coefficient,
                remainder_leading,
            )
            for degree, coefficient in divisor.items()
        }
        remainder = add_y_polynomials(
            arithmetic,
            scaled_remainder,
            scaled_divisor,
            -1,
        )
        remainder, content = primitive_y_polynomial(
            arithmetic,
            remainder,
        )
        if not content.is_constant():
            contents.append(polynomial_profile(content))
    return remainder, contents


def pseudo_gcd_y(
    arithmetic: QuotientArithmetic,
    left: YPolynomial,
    right: YPolynomial,
    label: str,
) -> tuple[YPolynomial, list[dict[str, object]]]:
    """Compute a generic-field gcd degree using a primitive PRS."""

    left, left_content = primitive_y_polynomial(arithmetic, left)
    right, right_content = primitive_y_polynomial(arithmetic, right)
    contents = [
        polynomial_profile(content)
        for content in (left_content, right_content)
        if not content.is_constant()
    ]
    if not left:
        return right, contents
    if not right:
        return left, contents
    if max(left) < max(right):
        left, right = right, left
    step = 0
    while right and max(right) > 0:
        started = time.monotonic()
        remainder, removed = pseudo_remainder_y(
            arithmetic,
            left,
            right,
        )
        contents.extend(removed)
        print(
            "CUSTOM_GCD_STEP "
            f"label={label} step={step} "
            f"degrees={max(left)},{max(right)},"
            f"{max(remainder, default=-1)} "
            f"seconds={time.monotonic() - started:.3f}",
            file=sys.stderr,
            flush=True,
        )
        left, right = right, remainder
        step += 1
    return (right if right else left), contents


def irreducible_modulus_specialization(
    arithmetic: QuotientArithmetic,
    prime: int,
) -> dict[str, object]:
    """Find a degree-preserving irreducible specialization of ``R20``."""

    if not prime:
        raise ValueError("the specialization scout requires a prime")
    for s1_value in range(1, min(prime, 13)):
        for lam_value in range(1, min(prime, 13)):
            if (6084 * lam_value * lam_value + 4805) % prime == 0:
                continue
            for v_value in range(1, min(prime, 13)):
                specialized = arithmetic.modulus.subs(
                    {
                        "s1": s1_value,
                        "lam": lam_value,
                        "v": v_value,
                    }
                )
                coefficients = [0] * (arithmetic.degree + 1)
                valid = True
                for monomial, coefficient in specialized.to_dict().items():
                    if any(monomial[index] for index in (1, 2, 3)):
                        valid = False
                        break
                    coefficients[monomial[0]] = int(coefficient)
                if not valid or not coefficients[-1] % prime:
                    continue
                polynomial = nmod_poly(coefficients, prime)
                unit, factors = polynomial.factor()
                if (
                    len(factors) == 1
                    and factors[0][0].degree() == arithmetic.degree
                    and factors[0][1] == 1
                ):
                    return {
                        "point": {
                            "s1": s1_value,
                            "lam": lam_value,
                            "v": v_value,
                        },
                        "polynomial": str(polynomial),
                        "factorization": str((unit, factors)),
                        "degree": polynomial.degree(),
                        "irreducible": True,
                    }
    return {"irreducible": False}


def run_raw_sympy_split_gcd(
    arithmetic: QuotientArithmetic,
    a_value: fmpq_mpoly,
    b_value: fmpq_mpoly,
    reduced: dict[int, RawQuotientPolynomial],
    prime: int,
    mode: str,
    base_specialization: dict[str, int] | None = None,
) -> dict[str, object]:
    """Use SymPy's generic finite extension for one minimal gcd branch."""

    if mode == "exceptional":
        orders = (5,)
    elif mode == "open6":
        orders = (4, 5, 6)
    else:
        raise ValueError(f"unknown SymPy split mode {mode}")
    started = time.monotonic()
    pivoted: dict[int, FibrePolynomial] = {}
    pivot_profiles: dict[str, object] = {}
    for order in orders:
        order_started = time.monotonic()
        polynomial = evaluate_quotient_polynomial_fraction_free(
            arithmetic,
            reduced[order],
            a_value,
            b_value,
        )
        polynomial, content = primitive_fibre_polynomial(polynomial)
        pivoted[order] = polynomial
        pivot_profiles[f"mu{order}"] = {
            **fibre_profile(arithmetic, polynomial),
            "seconds": round(time.monotonic() - order_started, 6),
            "removed_content": polynomial_profile(content),
        }
        print(
            "SYMPY_GCD_PIVOT "
            f"mode={mode} order={order} "
            f"seconds={time.monotonic() - order_started:.3f}",
            file=sys.stderr,
            flush=True,
        )

    root_symbol, s1_symbol, lam_symbol, v_symbol, y_symbol = sp.symbols(
        "T s1 lam v s5"
    )
    base_specialization = dict(base_specialization or {})
    base_coordinates = (
        ("s1", s1_symbol, 1),
        ("lam", lam_symbol, 2),
        ("v", v_symbol, 3),
    )
    active_coordinates = tuple(
        coordinate
        for coordinate in base_coordinates
        if coordinate[0] not in base_specialization
    )
    active_symbols = tuple(
        coordinate[1] for coordinate in active_coordinates
    )
    if active_symbols:
        base_ring = sp.GF(prime).poly_ring(*active_symbols)
        base_field = sp.GF(prime).frac_field(*active_symbols)
    else:
        base_ring = None
        base_field = sp.GF(prime)

    def as_root_polynomial(value) -> sp.Poly:
        """Transfer a FLINT sparse polynomial without parsing its string.

        The pivoted coefficients are large enough that Python's parser can
        exceed its recursion limit even though the sparse representation is
        quite manageable.  Grouping the FLINT monomial dictionary by the
        root degree also preserves the coefficient field exactly.
        """

        if base_specialization:
            value = value.subs(
                {
                    name: specialized_value % prime
                    for name, specialized_value
                    in base_specialization.items()
                }
            )
        grouped: dict[int, dict[tuple[int, ...], int]] = {}
        for exponents, coefficient in value.to_dict().items():
            root_degree_value = int(exponents[0])
            grouped.setdefault(root_degree_value, {})[
                tuple(
                    int(exponents[index])
                    for _, _, index in active_coordinates
                )
            ] = int(coefficient)
        if active_symbols:
            coefficients = {
                (root_degree_value,): base_field.convert(
                    base_ring.new(base_ring.ring.from_dict(terms)),
                    base_ring,
                )
                for root_degree_value, terms in grouped.items()
            }
        else:
            coefficients = {
                (root_degree_value,): sum(terms.values()) % prime
                for root_degree_value, terms in grouped.items()
            }
        return sp.Poly.from_dict(
            coefficients,
            (root_symbol,),
            domain=base_field,
        )

    modulus = as_root_polynomial(arithmetic.modulus).monic()
    extension = FiniteExtension(modulus)

    def extension_element(value):
        polynomial = as_root_polynomial(value)
        return extension.new(polynomial.rep.rem(extension.mod))

    converted: dict[int, dict[tuple[int, int], object]] = {}
    for order, polynomial in pivoted.items():
        converted[order] = {
            monomial: extension_element(coefficient.numerator)
            for monomial, coefficient in polynomial.items()
        }

    def y_polynomial(order: int, x_degree: int) -> sp.Poly:
        return sp.Poly.from_dict(
            {
                (y_degree,): coefficient
                for (degree, y_degree), coefficient
                in converted[order].items()
                if degree == x_degree
            },
            (y_symbol,),
            domain=extension,
        )

    original_inverse = ExtensionElement.inverse

    def generic_extension_inverse(value):
        bezout, _, gcd_value = value.rep.gcdex(value.ext.mod)
        if gcd_value.degree() > 0:
            raise ZeroDivisionError("nonunit in the residual extension")
        return value.ext.new(
            bezout.exquo_ground(gcd_value.LC())
        )

    ExtensionElement.inverse = generic_extension_inverse
    try:
        p5 = y_polynomial(5, 1)
        z5 = y_polynomial(5, 0)
        if mode == "exceptional":
            left = p5
            right = z5
        else:
            p6 = y_polynomial(6, 1)
            z6 = y_polynomial(6, 0)
            a4 = y_polynomial(4, 2)
            b4 = y_polynomial(4, 1)
            c4 = y_polynomial(4, 0)
            left = a4 * z5**2 - b4 * p5 * z5 + c4 * p5**2
            right = p5 * z6 - p6 * z5
        resultant_started = time.monotonic()
        branch_resultant = None
        if base_specialization:
            left_coefficients = left.rep.to_dict()
            right_coefficients = right.rep.to_dict()
            left_degree = int(left.degree())
            right_degree = int(right.degree())
            size = left_degree + right_degree
            matrix = [
                [extension.zero for _ in range(size)]
                for _ in range(size)
            ]
            left_descending = [
                left_coefficients.get(
                    (left_degree - index,),
                    extension.zero,
                )
                for index in range(left_degree + 1)
            ]
            right_descending = [
                right_coefficients.get(
                    (right_degree - index,),
                    extension.zero,
                )
                for index in range(right_degree + 1)
            ]
            for shift in range(right_degree):
                for index, coefficient in enumerate(left_descending):
                    matrix[shift][shift + index] = coefficient
            for shift in range(left_degree):
                for index, coefficient in enumerate(right_descending):
                    matrix[right_degree + shift][shift + index] = (
                        coefficient
                    )
            # Division-free Laplace dynamic programming evaluates a
            # size-n determinant in O(n*2^n) extension operations.  It
            # avoids both SymPy's broken exact quotient for this generic
            # extension and the factorial cost of a Leibniz expansion.
            minors = {0: extension.one}
            for row in range(size):
                next_minors = {}
                required_size = row + 1
                for mask in range(1 << size):
                    if mask.bit_count() != required_size:
                        continue
                    value = extension.zero
                    position = 0
                    for column in range(size):
                        bit = 1 << column
                        if not mask & bit:
                            continue
                        coefficient = matrix[row][column]
                        if coefficient:
                            term = (
                                minors[mask ^ bit] * coefficient
                            )
                            if (row + position) % 2:
                                term = -term
                            value += term
                        position += 1
                    next_minors[mask] = value
                minors = next_minors
            branch_resultant = minors[(1 << size) - 1]
        resultant_seconds = time.monotonic() - resultant_started
        norm_profile = None
        if branch_resultant is not None:
            branch_norm = modulus.rep.resultant(branch_resultant.rep)
            norm_profile = {
                "numerator": str(branch_norm.numer),
                "denominator": str(branch_norm.denom),
            }
            if len(active_symbols) == 1:
                unit, factors = branch_norm.numer.factor_list()
                norm_profile["factorization"] = {
                    "unit": str(unit),
                    "factors": [
                        {
                            "factor": str(factor),
                            "degree": int(factor.degree()),
                            "multiplicity": int(multiplicity),
                        }
                        for factor, multiplicity in factors
                    ],
                }
        if branch_resultant is not None and branch_resultant:
            branch_gcd = sp.Poly(
                extension.one,
                y_symbol,
                domain=extension,
            )
            gcd_seconds = 0.0
        else:
            gcd_started = time.monotonic()
            branch_gcd = sp.polys.polytools.gcd(left, right)
            gcd_seconds = time.monotonic() - gcd_started
    finally:
        ExtensionElement.inverse = original_inverse
    modulus_specialization = irreducible_modulus_specialization(
        arithmetic,
        prime,
    )
    gcd_degree = int(branch_gcd.degree())
    coprime = (
        gcd_degree == 0
        and modulus_specialization["irreducible"]
    )
    return {
        "status": "coprime" if coprime else "common-factor-or-unverified",
        "mode": mode,
        "input_degrees": [int(left.degree()), int(right.degree())],
        "gcd_degree": gcd_degree,
        "gcd_string_length": len(str(branch_gcd.as_expr())),
        "gcd_seconds": round(gcd_seconds, 6),
        "resultant_is_zero": (
            None if branch_resultant is None else not bool(branch_resultant)
        ),
        "resultant_seconds": round(resultant_seconds, 6),
        "base_specialization": base_specialization,
        "norm_profile": norm_profile,
        "coprime_on_generic_base": coprime,
        "modulus_specialization": modulus_specialization,
        "pivoted": pivot_profiles,
        "seconds": round(time.monotonic() - started, 6),
        "scope": (
            "generic field GF(p)(s1,lam,v)[T]/(R20), with every "
            "fraction-free pivot content and dense-pivot norm inverted"
        ),
    }


def run_raw_groebner_julia(
    system: dict[str, object],
    prime: int,
    solver_seconds: int,
    julia_project: Path,
    threads: int,
    with_change_matrix: bool = False,
) -> dict[str, object]:
    """Run Groebner.jl/F4 on a compact ordinary modular system.

    The ordinary F4 result is only a fast modular scout.  In change-matrix
    mode the exact identity ``matrix * inputs == basis`` is replayed in
    AbstractAlgebra; if a nonzero constant occurs in that verified basis,
    its row is a genuine modular ideal-membership certificate.
    """

    julia = shutil.which("julia")
    if julia is None:
        raise RuntimeError("Julia is required for the Groebner.jl stage")
    variables = list(system["ordinary_variables"])
    variable_tuple = ", ".join(variables)
    declarations = ",\n".join(
        f"    ({polynomial})"
        for polynomial in system["polynomials"]
    )
    solve_expression = (
        """basis, change = groebner_with_change_matrix(
    polynomials;
    ordering=DegRevLex(),
    reduced=false,
    certify=false,
    linalg=:deterministic,
    monoms=:packed
)
change_verified = change * polynomials == basis"""
        if with_change_matrix
        else """basis = groebner(
    polynomials;
    ordering=DegRevLex(),
    reduced=false,
    certify=false,
    linalg=:deterministic,
    monoms=:packed
)
change_verified = false"""
    )
    source = f"""\
using Groebner
using AbstractAlgebra

ring, generators = polynomial_ring(
    GF({prime}),
    {json.dumps(variables)}
)
({variable_tuple}) = generators
polynomials = [
{declarations}
]
started = time()
{solve_expression}
elapsed = time() - started
unit_indices = findall(
    polynomial -> !iszero(polynomial) && total_degree(polynomial) == 0,
    basis
)
unit = !isempty(unit_indices)
println("GROEBNER_VERSION=" * string(Base.pkgversion(Groebner)))
println("ABSTRACTALGEBRA_VERSION=" * string(Base.pkgversion(AbstractAlgebra)))
println("BASIS_LENGTH=" * string(length(basis)))
println("BASIS_MAX_TERMS=" * string(maximum(length, basis)))
println("BASIS_MAX_TOTAL_DEGREE=" * string(maximum(total_degree, basis)))
println("UNIT_FOUND=" * string(unit))
println("CHANGE_MATRIX_VERIFIED=" * string(change_verified))
println("UNIT_INDEX=" * string(isempty(unit_indices) ? 0 : first(unit_indices)))
println("SOLVER_SECONDS=" * string(elapsed))
"""
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="sic33-groebner-julia-") as temporary:
        source_path = Path(temporary) / "solve.jl"
        source_path.write_text(source, encoding="utf-8")
        environment = dict(os.environ)
        environment.setdefault(
            "JULIA_DEPOT_PATH",
            "/tmp/sic33-julia-depot",
        )
        environment["JULIA_NUM_THREADS"] = str(max(1, threads))
        try:
            completed = subprocess.run(
                [
                    julia,
                    f"--project={julia_project}",
                    "--startup-file=no",
                    str(source_path),
                ],
                text=True,
                capture_output=True,
                check=False,
                timeout=solver_seconds or None,
                env=environment,
            )
        except subprocess.TimeoutExpired as error:
            return {
                "status": "timeout",
                "seconds": round(time.monotonic() - started, 6),
                "stdout_tail": (error.stdout or "")[-4000:],
                "stderr_tail": (error.stderr or "")[-4000:],
                "prime": prime,
                "julia_project": str(julia_project),
            }
    markers = {}
    for line in completed.stdout.splitlines():
        if "=" in line and line.split("=", 1)[0] in {
            "GROEBNER_VERSION",
            "ABSTRACTALGEBRA_VERSION",
            "BASIS_LENGTH",
            "BASIS_MAX_TERMS",
            "BASIS_MAX_TOTAL_DEGREE",
            "UNIT_FOUND",
            "CHANGE_MATRIX_VERIFIED",
            "UNIT_INDEX",
            "SOLVER_SECONDS",
        }:
            key, value = line.split("=", 1)
            markers[key.lower()] = value
    unit_found = markers.get("unit_found") == "true"
    change_verified = markers.get("change_matrix_verified") == "true"
    verified_modular_unit = unit_found and change_verified
    return {
        "status": (
            "verified-modular-unit"
            if completed.returncode == 0 and verified_modular_unit
            else (
                "probable-modular-unit"
                if completed.returncode == 0 and unit_found
                else "nonunit-or-failed"
            )
        ),
        "returncode": completed.returncode,
        "probable_unit_mod_p": unit_found,
        "change_matrix_verified": change_verified,
        "verified_modular_unit": verified_modular_unit,
        "with_change_matrix": with_change_matrix,
        "markers": markers,
        "seconds": round(time.monotonic() - started, 6),
        "stdout_tail": completed.stdout[-4000:],
        "stderr_tail": completed.stderr[-4000:],
        "prime": prime,
        "julia_project": str(julia_project),
        "scope": (
            "the compact ordinary Q-component system through mu8 over "
            "GF(p), including the v and J_Q localizers; a modular unit "
            "basis is not by itself a characteristic-zero certificate. "
            "Without a verified change matrix the F4 result is a "
            "probabilistic scout."
        ),
    }


def run_raw_custom_gcd(
    arithmetic: QuotientArithmetic,
    a_value: fmpq_mpoly,
    b_value: fmpq_mpoly,
    reduced: dict[int, RawQuotientPolynomial],
    prime: int,
) -> dict[str, object]:
    """Use compact pre-pivot remainders and a custom univariate PRS."""

    started = time.monotonic()
    pivoted: dict[int, FibrePolynomial] = {}
    pivot_profiles: dict[str, object] = {}
    for order in range(4, 9):
        order_started = time.monotonic()
        polynomial = evaluate_quotient_polynomial_fraction_free(
            arithmetic,
            reduced[order],
            a_value,
            b_value,
        )
        polynomial, content = primitive_fibre_polynomial(polynomial)
        pivoted[order] = polynomial
        pivot_profiles[f"mu{order}"] = {
            **fibre_profile(arithmetic, polynomial),
            "seconds": round(time.monotonic() - order_started, 6),
            "removed_content": polynomial_profile(content),
        }
        print(
            "CUSTOM_GCD_PIVOT "
            f"order={order} "
            f"seconds={time.monotonic() - order_started:.3f}",
            file=sys.stderr,
            flush=True,
        )
    equations = remainder_univariate_equations(
        arithmetic,
        pivoted[4],
        {order: pivoted[order] for order in (5, 6, 7)},
    )
    p5 = x_coefficient(pivoted[5], 1)
    z5 = x_coefficient(pivoted[5], 0)
    p8 = x_coefficient(pivoted[8], 1)
    z8 = x_coefficient(pivoted[8], 0)
    equations["cross_8"] = add_y_polynomials(
        arithmetic,
        multiply_y_polynomials(arithmetic, p5, z8),
        multiply_y_polynomials(arithmetic, p8, z5),
        -1,
    )
    exceptional_gcd, exceptional_contents = pseudo_gcd_y(
        arithmetic,
        p5,
        z5,
        "exceptional",
    )
    open_gcd = equations["quadratic"]
    open_contents: list[dict[str, object]] = []
    for name in ("cross_6", "cross_7", "cross_8"):
        open_gcd, removed = pseudo_gcd_y(
            arithmetic,
            open_gcd,
            equations[name],
            name,
        )
        open_contents.extend(removed)
        if open_gcd and max(open_gcd) == 0:
            break
    exceptional_degree = max(exceptional_gcd, default=-1)
    open_degree = max(open_gcd, default=-1)
    modulus_specialization = irreducible_modulus_specialization(
        arithmetic,
        prime,
    )
    unit_ideal = (
        modulus_specialization["irreducible"]
        and exceptional_degree == 0
        and open_degree == 0
    )
    return {
        "status": "unit" if unit_ideal else "nonunit-or-unverified",
        "unit_ideal_on_generic_base": unit_ideal,
        "modulus_specialization": modulus_specialization,
        "pivoted": pivot_profiles,
        "equations": {
            name: y_polynomial_profile(arithmetic, polynomial)
            for name, polynomial in equations.items()
        },
        "exceptional_gcd": {
            **y_polynomial_profile(arithmetic, exceptional_gcd),
            "removed_contents": exceptional_contents,
        },
        "open_gcd": {
            **y_polynomial_profile(arithmetic, open_gcd),
            "removed_contents": open_contents,
        },
        "denominator_factors": arithmetic.registry.describe(),
        "seconds": round(time.monotonic() - started, 6),
        "scope": (
            "generic base field GF(p)(s1,lam,v)[T]/(R20), with the "
            "dense pivot coefficient and every removed nonzero content "
            "inverted; exceptional base divisors remain to be treated"
        ),
    }


def y_polynomial_profile(
    arithmetic: QuotientArithmetic,
    polynomial: YPolynomial,
) -> dict[str, object]:
    return {
        "degree_s5": max(polynomial, default=-1),
        "terms_s5": len(polynomial),
        "maximum_coefficient_terms": max(
            (
                len(coefficient.numerator.to_dict())
                for coefficient in polynomial.values()
            ),
            default=0,
        ),
        "coefficient_terms": {
            str(degree): len(coefficient.numerator.to_dict())
            for degree, coefficient in sorted(polynomial.items())
        },
    }


def quadratic_pseudo_remainder(
    arithmetic: QuotientArithmetic,
    divisor: FibrePolynomial,
    dividend: FibrePolynomial,
) -> FibrePolynomial:
    """Pseudo-divide ``dividend`` by quadratic ``divisor`` in ``s6``."""

    divisor_by_x = {
        degree: x_coefficient(divisor, degree)
        for degree in range(3)
    }
    if set(divisor_by_x[2]) != {0}:
        raise AssertionError("the quadratic leading coefficient depends on s5")
    leading = divisor_by_x[2][0]
    remainder_by_x = {
        degree: x_coefficient(dividend, degree)
        for degree in range(
            max((monomial[0] for monomial in dividend), default=-1) + 1
        )
    }
    remainder_by_x = {
        degree: polynomial
        for degree, polynomial in remainder_by_x.items()
        if polynomial
    }
    while remainder_by_x and max(remainder_by_x) >= 2:
        top_degree = max(remainder_by_x)
        top = remainder_by_x[top_degree]
        shift = top_degree - 2
        scaled = {
            degree: scale_y_polynomial(
                arithmetic,
                polynomial,
                leading,
            )
            for degree, polynomial in remainder_by_x.items()
        }
        for divisor_degree, divisor_coefficient in divisor_by_x.items():
            target_degree = divisor_degree + shift
            product = multiply_y_polynomials(
                arithmetic,
                top,
                divisor_coefficient,
            )
            scaled[target_degree] = add_y_polynomials(
                arithmetic,
                scaled.get(target_degree, {}),
                product,
                -1,
            )
        remainder_by_x = {
            degree: clear_y_denominators(arithmetic, polynomial)
            for degree, polynomial in scaled.items()
            if polynomial
        }
        remainder_by_x = {
            degree: polynomial
            for degree, polynomial in remainder_by_x.items()
            if polynomial
        }
    answer = {
        (x_degree, y_degree): coefficient
        for x_degree, polynomial in remainder_by_x.items()
        for y_degree, coefficient in polynomial.items()
    }
    return clear_fibre_denominators(arithmetic, answer)


def quadratic_remainder_equations(
    arithmetic: QuotientArithmetic,
    mu4: FibrePolynomial,
    mu5: FibrePolynomial,
    mu6: FibrePolynomial,
    mu7: FibrePolynomial,
) -> tuple[dict[int, FibrePolynomial], dict[str, YPolynomial]]:
    """Eliminate ``s6`` after linear pseudo-reduction modulo ``mu4``."""

    remainders = {
        order: quadratic_pseudo_remainder(
            arithmetic,
            mu4,
            polynomial,
        )
        for order, polynomial in (
            (5, mu5),
            (6, mu6),
            (7, mu7),
        )
    }
    equations = remainder_univariate_equations(
        arithmetic,
        mu4,
        remainders,
    )
    return remainders, equations


def remainder_univariate_equations(
    arithmetic: QuotientArithmetic,
    mu4: FibrePolynomial,
    remainders: dict[int, FibrePolynomial],
) -> dict[str, YPolynomial]:
    """Eliminate ``s6`` from three linear remainders and quadratic ``mu4``."""

    if set(remainders) != {5, 6, 7}:
        raise ValueError("the univariate equations require remainders 5, 6, 7")
    p = {
        order: x_coefficient(polynomial, 1)
        for order, polynomial in remainders.items()
    }
    q = {
        order: x_coefficient(polynomial, 0)
        for order, polynomial in remainders.items()
    }
    a = x_coefficient(mu4, 2)
    f1 = x_coefficient(mu4, 1)
    f0 = x_coefficient(mu4, 0)
    if set(a) != {0}:
        raise AssertionError("the quadratic leading coefficient depends on s5")
    quadratic = add_y_polynomials(
        arithmetic,
        scale_y_polynomial(
            arithmetic,
            multiply_y_polynomials(arithmetic, q[5], q[5]),
            a[0],
        ),
        multiply_y_polynomials(
            arithmetic,
            f1,
            multiply_y_polynomials(arithmetic, p[5], q[5]),
        ),
        -1,
    )
    quadratic = add_y_polynomials(
        arithmetic,
        quadratic,
        multiply_y_polynomials(
            arithmetic,
            f0,
            multiply_y_polynomials(arithmetic, p[5], p[5]),
        ),
    )
    equations = {"quadratic": clear_y_denominators(arithmetic, quadratic)}
    for order in (6, 7):
        cross = add_y_polynomials(
            arithmetic,
            multiply_y_polynomials(arithmetic, p[5], q[order]),
            multiply_y_polynomials(arithmetic, p[order], q[5]),
            -1,
        )
        equations[f"cross_{order}"] = clear_y_denominators(
            arithmetic,
            cross,
        )
    return equations


def generic_triangular_basis(
    arithmetic: QuotientArithmetic,
    mu4: FibrePolynomial,
    mu5: FibrePolynomial,
) -> tuple[list[FibrePolynomial], list[dict[str, object]]]:
    """Use the quadratic/cubic pseudo-remainder sequence in the s6 variable."""

    removed_contents: list[dict[str, object]] = []
    f2 = x_coefficient(mu4, 2)
    g3 = x_coefficient(mu5, 3)
    if set(f2) != {0} or set(g3) != {0}:
        raise AssertionError("unexpected top s6 coefficients")
    a = f2[0]
    d = g3[0]
    f1 = x_coefficient(mu4, 1)
    f0 = x_coefficient(mu4, 0)
    g2 = x_coefficient(mu5, 2)
    g1 = x_coefficient(mu5, 1)
    g0 = x_coefficient(mu5, 0)

    # r2 = a*g - d*s6*f has degree at most two in s6.
    r2_2 = add_y_polynomials(
        arithmetic,
        scale_y_polynomial(arithmetic, g2, a),
        scale_y_polynomial(arithmetic, f1, d),
        -1,
    )
    r2_1 = add_y_polynomials(
        arithmetic,
        scale_y_polynomial(arithmetic, g1, a),
        scale_y_polynomial(arithmetic, f0, d),
        -1,
    )
    r2_0 = scale_y_polynomial(arithmetic, g0, a)
    # r1 = a*r2 - e*f is linear in s6.
    p = add_y_polynomials(
        arithmetic,
        scale_y_polynomial(arithmetic, r2_1, a),
        multiply_y_polynomials(arithmetic, f1, r2_2),
        -1,
    )
    q = add_y_polynomials(
        arithmetic,
        scale_y_polynomial(arithmetic, r2_0, a),
        multiply_y_polynomials(arithmetic, f0, r2_2),
        -1,
    )
    linear_relation = {
        **{(1, degree): coefficient for degree, coefficient in p.items()},
        **{(0, degree): coefficient for degree, coefficient in q.items()},
    }
    linear_relation = normalize_fibre_generic(
        arithmetic,
        linear_relation,
        removed_contents,
        "linear pseudo-remainder",
    )
    p = x_coefficient(linear_relation, 1)
    q = x_coefficient(linear_relation, 0)

    # Res_s6(f, P*s6+Q) = a*Q^2-F1*P*Q+F0*P^2.
    resultant = add_y_polynomials(
        arithmetic,
        scale_y_polynomial(
            arithmetic,
            multiply_y_polynomials(arithmetic, q, q),
            a,
        ),
        multiply_y_polynomials(
            arithmetic,
            f1,
            multiply_y_polynomials(arithmetic, p, q),
        ),
        -1,
    )
    resultant = add_y_polynomials(
        arithmetic,
        resultant,
        multiply_y_polynomials(
            arithmetic,
            f0,
            multiply_y_polynomials(arithmetic, p, p),
        ),
    )
    resultant_relation = normalize_fibre_generic(
        arithmetic,
        {
            (0, degree): coefficient
            for degree, coefficient in resultant.items()
        },
        removed_contents,
        "linear-resultant relation",
    )
    basis = [
        normalize_fibre_generic(
            arithmetic,
            mu4,
            removed_contents,
            "quadratic generator",
        ),
        linear_relation,
        resultant_relation,
    ]
    return basis, removed_contents


def payload_base(stage: str, elapsed: float) -> dict[str, object]:
    return {
        "format": "two-pair-sic-bidegree33-t0-Q-unspecialized-flint-v1",
        "status": (
            (
                "exact characteristic-zero research calculation on the "
                "unspecialized residual Q component"
            )
            if not PRIME
            else (
                f"exact finite-field research calculation modulo {PRIME} "
                "on the residual Q component; not a characteristic-zero "
                "certificate"
            )
        ),
        "stage": stage,
        "coefficient_ring": (
            "QQ[s1,lam,v], lam=(s1*u-t1)/u, v=u^2"
            if not PRIME
            else (
                f"GF({PRIME})[s1,lam,v], "
                "lam=(s1*u-t1)/u, v=u^2"
            )
        ),
        "extension": (
            "QQ(s1,lam,v)[T]/(R20)"
            if not PRIME
            else f"GF({PRIME})(s1,lam,v)[T]/(R20)"
        ),
        "coordinate_map": {
            "geometric_ell": "s1*u-t1",
            "lam": "(s1*u-t1)/u",
            "v": "u^2",
        },
        "prime": PRIME,
        "seconds": round(elapsed, 6),
        "reproduction_command": " ".join(sys.argv),
    }


def main() -> None:
    arguments = parse_arguments()
    if arguments.prime and arguments.prime in (2, 3, 5, 7, 13):
        raise ValueError("choose a prime avoiding the displayed denominators")
    configure_coefficient_ring(arguments.prime)
    started = time.monotonic()
    arithmetic, a_value, b_value, metadata = flint_input_data()
    payload = payload_base(arguments.stage, time.monotonic() - started)
    payload["input"] = metadata
    if arguments.stage == "prepare":
        pass
    else:
        pivot, pivot_metadata = compute_pivot(
            arithmetic,
            a_value,
            b_value,
            metadata["contents"]["A"],
            metadata["contents"]["B"],
        )
        payload["pivot"] = pivot_metadata
        if arguments.stage not in ("pivot",):
            singular = shutil.which("Singular")
            if singular is None:
                raise RuntimeError("Singular is required for moment export")
            moments, moment_profiles = parse_moments(
                singular,
                arguments.timeout,
                arguments.moment_through,
            )
            payload["moments"] = moment_profiles
            if arguments.stage not in ("moments",):
                if arguments.stage in (
                    "raw8",
                    "raw8_solve",
                    "raw8_popen",
                    "raw8_popen_solve",
                    "raw8_exceptional",
                    "raw8_exceptional_solve",
                    "raw8_direct",
                    "raw8_direct_solve",
                    "raw8_direct_generic",
                    "raw8_direct_generic_lift",
                    "raw8_direct_extension",
                    "raw8_direct_extension_lift",
                    "raw8_direct_extension_gcd",
                    "raw8_direct_custom_gcd",
                    "raw8_direct_m2_gcd",
                    "raw8_direct_m2_exceptional",
                    "raw8_direct_m2_open6",
                    "raw8_direct_sympy_exceptional",
                    "raw8_direct_sympy_open6",
                    "raw8_direct_groebner_julia",
                    "raw8_direct_popen_groebner_julia",
                    "raw8_direct_popen_groebner_julia_change",
                    "raw8_popen_groebner_julia",
                    "raw8_popen_groebner_julia_change",
                    "raw8_exceptional_groebner_julia",
                    "raw8_exceptional_groebner_julia_change",
                ):
                    if not PRIME:
                        raise ValueError(
                            "the bounded raw mu8 implementation is modular; "
                            "pass --prime"
                        )
                    if arguments.moment_through < 8:
                        raise ValueError(
                            "raw mu8 stages require --moment-through 8"
                        )
                    raw_started = time.monotonic()
                    if arguments.stage.startswith("raw8_direct"):
                        raw_system, raw_remainders, raw_equations = (
                            raw_mu8_direct_system(
                                arithmetic,
                                a_value,
                                b_value,
                                moments,
                                p5_open=arguments.stage.startswith(
                                    "raw8_direct_popen"
                                ),
                            )
                        )
                    elif arguments.stage.startswith("raw8_exceptional"):
                        raw_system, raw_remainders, raw_equations = (
                            raw_mu8_exceptional_system(
                                arithmetic,
                                a_value,
                                b_value,
                                moments,
                            )
                        )
                    else:
                        raw_system, raw_remainders, raw_equations = (
                            raw_mu8_system(
                                arithmetic,
                                a_value,
                                b_value,
                                moments,
                                p5_open=arguments.stage.startswith(
                                    "raw8_popen"
                                ),
                            )
                        )
                    payload["raw_mu8_system"] = {
                        **{
                            key: value
                            for key, value in raw_system.items()
                            if key != "polynomials"
                        },
                        "seconds": round(
                            time.monotonic() - raw_started,
                            6,
                        ),
                        "remainders": {
                            f"mu{order}": raw_quotient_profile(polynomial)
                            for order, polynomial in raw_remainders.items()
                        },
                        "equations": {
                            name: yz_polynomial_profile(polynomial)
                            for name, polynomial in raw_equations.items()
                        },
                    }
                    if arguments.stage in (
                        "raw8_solve",
                        "raw8_popen_solve",
                        "raw8_exceptional_solve",
                        "raw8_direct_solve",
                    ):
                        payload["raw_mu8_solve"] = (
                            run_direct_msolve_unbounded(
                                raw_system,
                                arguments.prime,
                                arguments.threads,
                                arguments.solver_seconds,
                                (
                                    None
                                    if arguments.solver_result is None
                                    else (
                                        arguments.solver_result
                                        if arguments.solver_result.is_absolute()
                                        else ROOT / arguments.solver_result
                                    )
                                ),
                            )
                        )
                    elif arguments.stage in (
                        "raw8_direct_generic",
                        "raw8_direct_generic_lift",
                    ):
                        payload["raw_mu8_generic_solve"] = (
                            run_raw_generic_singular(
                                raw_system,
                                arguments.prime,
                                arguments.solver_seconds,
                                lift=arguments.stage.endswith("_lift"),
                            )
                        )
                    elif arguments.stage in (
                        "raw8_direct_extension",
                        "raw8_direct_extension_lift",
                    ):
                        payload["raw_mu8_extension_solve"] = (
                            run_raw_extension_singular(
                                raw_system,
                                arguments.prime,
                                arguments.solver_seconds,
                                lift=arguments.stage.endswith("_lift"),
                            )
                        )
                    elif arguments.stage == "raw8_direct_extension_gcd":
                        payload["raw_mu8_extension_gcd"] = (
                            run_raw_extension_gcd_singular(
                                raw_system,
                                arguments.prime,
                                arguments.solver_seconds,
                            )
                        )
                    elif arguments.stage == "raw8_direct_custom_gcd":
                        payload["raw_mu8_custom_gcd"] = run_raw_custom_gcd(
                            arithmetic,
                            a_value,
                            b_value,
                            raw_remainders,
                            arguments.prime,
                        )
                    elif arguments.stage == "raw8_direct_m2_gcd":
                        payload["raw_mu8_m2_gcd"] = (
                            run_raw_extension_gcd_macaulay2(
                                raw_system,
                                arithmetic,
                                pivot,
                                arguments.prime,
                                arguments.solver_seconds,
                            )
                        )
                    elif arguments.stage in (
                        "raw8_direct_m2_exceptional",
                        "raw8_direct_m2_open6",
                    ):
                        mode = (
                            "exceptional"
                            if arguments.stage.endswith("_exceptional")
                            else "open6"
                        )
                        payload["raw_mu8_m2_split"] = (
                            run_raw_extension_m2_split(
                                raw_system,
                                arithmetic,
                                pivot,
                                arguments.prime,
                                arguments.solver_seconds,
                                mode,
                            )
                        )
                    elif arguments.stage in (
                        "raw8_direct_sympy_exceptional",
                        "raw8_direct_sympy_open6",
                    ):
                        mode = (
                            "exceptional"
                            if arguments.stage.endswith("_exceptional")
                            else "open6"
                        )
                        payload["raw_mu8_sympy_split"] = (
                            run_raw_sympy_split_gcd(
                                arithmetic,
                                a_value,
                                b_value,
                                raw_remainders,
                                arguments.prime,
                                mode,
                                {
                                    name: value
                                    for name, value in (
                                        ("s1", arguments.s1_value),
                                        ("lam", arguments.lam_value),
                                        ("v", arguments.v_value),
                                    )
                                    if value is not None
                                },
                            )
                        )
                    elif arguments.stage in (
                        "raw8_direct_groebner_julia",
                        "raw8_direct_popen_groebner_julia",
                        "raw8_direct_popen_groebner_julia_change",
                        "raw8_popen_groebner_julia",
                        "raw8_popen_groebner_julia_change",
                        "raw8_exceptional_groebner_julia",
                        "raw8_exceptional_groebner_julia_change",
                    ):
                        payload["raw_mu8_groebner_julia"] = (
                            run_raw_groebner_julia(
                                raw_system,
                                arguments.prime,
                                arguments.solver_seconds,
                                arguments.julia_project,
                                arguments.threads,
                                with_change_matrix=arguments.stage.endswith(
                                    "_change"
                                ),
                            )
                        )
                    payload["seconds"] = round(
                        time.monotonic() - started,
                        6,
                    )
                    output = arguments.output
                    if not output.is_absolute():
                        output = ROOT / output
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.write_text(
                        json.dumps(payload, indent=2, sort_keys=True)
                        + "\n",
                        encoding="utf-8",
                    )
                    print(json.dumps(payload, indent=2, sort_keys=True))
                    return
                if arguments.stage in ("system", "solve"):
                    system = direct_residual_system(
                        arithmetic,
                        a_value,
                        b_value,
                        moments,
                        arguments.s1_value,
                    )
                    payload["direct_system"] = {
                        key: value
                        for key, value in system.items()
                        if key != "polynomials"
                    }
                    if arguments.stage == "solve":
                        payload["solve"] = run_direct_msolve_unbounded(
                            system,
                            arguments.prime,
                            arguments.threads,
                            arguments.solver_seconds,
                            (
                                None
                                if arguments.solver_result is None
                                else (
                                    arguments.solver_result
                                    if arguments.solver_result.is_absolute()
                                    else ROOT / arguments.solver_result
                                )
                            ),
                        )
                    payload["seconds"] = round(
                        time.monotonic() - started,
                        6,
                    )
                    output = arguments.output
                    if not output.is_absolute():
                        output = ROOT / output
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.write_text(
                        json.dumps(payload, indent=2, sort_keys=True)
                        + "\n",
                        encoding="utf-8",
                    )
                    print(json.dumps(payload, indent=2, sort_keys=True))
                    return
                if arguments.stage in (
                    "prepivot5",
                    "prepivot6",
                    "prepivot7",
                    "prepivot",
                    "raw_equations",
                    "prepivot_equations",
                    "prepivot_quadratic",
                    "prepivot_cross6",
                    "prepivot_cross7",
                ):
                    if not PRIME:
                        raise ValueError(
                            "the bounded pre-pivot remainder implementation "
                            "is modular; pass --prime"
                        )
                    included_orders = {
                        "prepivot5": (5,),
                        "prepivot6": (5, 6),
                        "prepivot7": (7,),
                        "prepivot": (5, 6, 7),
                        "raw_equations": (5, 6, 7),
                        "prepivot_equations": (5, 6, 7),
                        "prepivot_quadratic": (5, 6, 7),
                        "prepivot_cross6": (5, 6, 7),
                        "prepivot_cross7": (5, 6, 7),
                    }[arguments.stage]
                    raw_remainders = {}
                    raw_seconds = {}
                    evaluated_remainders = {}
                    evaluated_seconds = {}
                    for order in included_orders:
                        print(
                            f"pre-pivot pseudo-dividing mu{order} by mu4",
                            file=sys.stderr,
                            flush=True,
                        )
                        remainder_started = time.monotonic()
                        raw_remainders[order] = (
                            raw_quadratic_pseudo_remainder(
                                arithmetic,
                                moments[4],
                                moments[order],
                            )
                        )
                        raw_seconds[f"mu{order}"] = round(
                            time.monotonic() - remainder_started,
                            6,
                        )
                        if arguments.stage not in (
                            "raw_equations",
                            "prepivot_equations",
                            "prepivot_quadratic",
                            "prepivot_cross6",
                            "prepivot_cross7",
                        ):
                            print(
                                f"substituting the pivot in remainder mu{order}",
                                file=sys.stderr,
                                flush=True,
                            )
                            evaluation_started = time.monotonic()
                            evaluated_remainders[order] = (
                                evaluate_raw_quotient_fraction_free(
                                    arithmetic,
                                    raw_remainders[order],
                                    a_value,
                                    b_value,
                                )
                            )
                            evaluated_seconds[f"mu{order}"] = round(
                                time.monotonic() - evaluation_started,
                                6,
                            )
                    equations = {}
                    raw_equations = {}
                    mu4 = None
                    if arguments.stage == "prepivot":
                        print(
                            "evaluating mu4 for univariate elimination",
                            file=sys.stderr,
                            flush=True,
                        )
                        evaluation_started = time.monotonic()
                        mu4 = evaluate_moment_fraction_free(
                            arithmetic,
                            moments[4],
                            a_value,
                            b_value,
                        )
                        evaluated_seconds["mu4"] = round(
                            time.monotonic() - evaluation_started,
                            6,
                        )
                        equations = remainder_univariate_equations(
                            arithmetic,
                            mu4,
                            evaluated_remainders,
                        )
                    if arguments.stage in (
                        "raw_equations",
                        "prepivot_equations",
                        "prepivot_quadratic",
                        "prepivot_cross6",
                        "prepivot_cross7",
                    ):
                        print(
                            "forming univariate equations before pivot substitution",
                            file=sys.stderr,
                            flush=True,
                        )
                        equations_started = time.monotonic()
                        raw_equations = raw_remainder_univariate_equations(
                            arithmetic,
                            moments[4],
                            raw_remainders,
                        )
                        raw_seconds["equations"] = round(
                            time.monotonic() - equations_started,
                            6,
                        )
                    selected_equations = {
                        "prepivot_equations": tuple(raw_equations),
                        "prepivot_quadratic": ("quadratic",),
                        "prepivot_cross6": ("cross_6",),
                        "prepivot_cross7": ("cross_7",),
                    }.get(arguments.stage, ())
                    for name in selected_equations:
                        polynomial = raw_equations[name]
                        print(
                            f"substituting the pivot in {name}",
                            file=sys.stderr,
                            flush=True,
                        )
                        evaluation_started = time.monotonic()
                        equations[name] = evaluate_yz_fraction_free(
                            arithmetic,
                            polynomial,
                            a_value,
                            b_value,
                        )
                        evaluated_seconds[name] = round(
                            time.monotonic() - evaluation_started,
                            6,
                        )
                    payload["prepivot_remainders"] = {
                        "algorithm": (
                            "quadratic pseudo-division in s6 before the "
                            "dense s3 pivot substitution"
                        ),
                        "raw_seconds": raw_seconds,
                        "pivot_substitution_seconds": evaluated_seconds,
                        "raw": {
                            f"mu{order}": raw_quotient_profile(polynomial)
                            for order, polynomial in raw_remainders.items()
                        },
                        "raw_univariate_equations": {
                            name: yz_polynomial_profile(polynomial)
                            for name, polynomial in raw_equations.items()
                        },
                        "evaluated": {
                            f"mu{order}": fibre_profile(
                                arithmetic,
                                polynomial,
                            )
                            for order, polynomial
                            in evaluated_remainders.items()
                        },
                        "mu4": (
                            None
                            if mu4 is None
                            else fibre_profile(arithmetic, mu4)
                        ),
                        "univariate_equations": {
                            name: y_polynomial_profile(
                                arithmetic,
                                polynomial,
                            )
                            for name, polynomial in equations.items()
                        },
                        "exceptional_case": "P5=Q5=0",
                    }
                    payload["seconds"] = round(
                        time.monotonic() - started,
                        6,
                    )
                    output = arguments.output
                    if not output.is_absolute():
                        output = ROOT / output
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.write_text(
                        json.dumps(payload, indent=2, sort_keys=True)
                        + "\n",
                        encoding="utf-8",
                    )
                    print(json.dumps(payload, indent=2, sort_keys=True))
                    return
                evaluated_orders = {
                    "remainder5": (4, 5),
                    "remainder6": (4, 5, 6),
                    "remainders": (4, 5, 6, 7),
                }.get(arguments.stage, (4, 5))
                evaluated = {}
                evaluated_seconds = {}
                for order in evaluated_orders:
                    print(
                        f"evaluating mu{order}",
                        file=sys.stderr,
                        flush=True,
                    )
                    evaluation_started = time.monotonic()
                    evaluated[order] = evaluate_moment_fraction_free(
                        arithmetic,
                        moments[order],
                        a_value,
                        b_value,
                    )
                    evaluated_seconds[f"mu{order}"] = round(
                        time.monotonic() - evaluation_started,
                        6,
                    )
                    print(
                        f"evaluated mu{order} in "
                        f"{evaluated_seconds[f'mu{order}']} seconds",
                        file=sys.stderr,
                        flush=True,
                    )
                contents = {}
                for order in evaluated_orders:
                    if PRIME:
                        contents[order] = CTX.constant(1)
                    else:
                        evaluated[order], contents[order] = (
                            primitive_fibre_polynomial(evaluated[order])
                        )
                payload["evaluated"] = {
                    f"mu{order}": {
                        **fibre_profile(
                            arithmetic,
                            evaluated[order],
                        ),
                        "removed_polynomial_content": polynomial_profile(
                            contents[order]
                        ),
                    }
                    for order in evaluated_orders
                }
                payload["evaluated_seconds"] = evaluated_seconds
                if arguments.stage in (
                    "remainder5",
                    "remainder6",
                    "remainders",
                ):
                    if not PRIME:
                        raise ValueError(
                            "the first bounded remainder implementation "
                            "is modular; pass --prime"
                        )
                    remainder_started = time.monotonic()
                    if arguments.stage == "remainders":
                        remainders, equations = (
                            quadratic_remainder_equations(
                                arithmetic,
                                evaluated[4],
                                evaluated[5],
                                evaluated[6],
                                evaluated[7],
                            )
                        )
                    else:
                        included_orders = (
                            (5,)
                            if arguments.stage == "remainder5"
                            else (5, 6)
                        )
                        remainders = {}
                        for order in included_orders:
                            print(
                                f"pseudo-dividing mu{order} by mu4",
                                file=sys.stderr,
                                flush=True,
                            )
                            remainders[order] = (
                                quadratic_pseudo_remainder(
                                    arithmetic,
                                    evaluated[4],
                                    evaluated[order],
                                )
                            )
                        equations = {}
                    payload["remainders"] = {
                        "algorithm": (
                            "fraction-free pseudo-division by quadratic "
                            "mu4 in the rescaled s6 variable"
                        ),
                        "seconds": round(
                            time.monotonic() - remainder_started,
                            6,
                        ),
                        "linear_remainders": {
                            f"mu{order}": fibre_profile(
                                arithmetic,
                                polynomial,
                            )
                            for order, polynomial in remainders.items()
                        },
                        "univariate_equations": {
                            name: y_polynomial_profile(
                                arithmetic,
                                polynomial,
                            )
                            for name, polynomial in equations.items()
                        },
                        "exceptional_case": "P5=Q5=0",
                    }
                elif arguments.stage not in ("evaluated",):
                    basis_started = time.monotonic()
                    basis, removed_contents = generic_triangular_basis(
                        arithmetic,
                        evaluated[4],
                        evaluated[5],
                    )
                    payload["basis"] = {
                        "algorithm": (
                            "quadratic-cubic pseudo-remainder sequence "
                            "in the rescaled s6 variable"
                        ),
                        "seconds": round(
                            time.monotonic() - basis_started,
                            6,
                        ),
                        "size": len(basis),
                        "leading_monomials": [
                            list(leading_term_fibre(polynomial)[0])
                            for polynomial in basis
                        ],
                        "polynomials": [
                            fibre_profile(arithmetic, polynomial)
                            for polynomial in basis
                        ],
                        "removed_generic_contents": removed_contents,
                    }
                    if arguments.stage not in ("basis",):
                        raise NotImplementedError(
                            f"stage {arguments.stage} will consume the "
                            "generic fibre basis"
                        )
    payload["seconds"] = round(time.monotonic() - started, 6)
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
