#!/usr/bin/env python3
"""Exact unspecialized arithmetic on the residual Q component.

This research driver replaces the rational-function normalization in
``research_two_pair_sic_bidegree33_t0_Q_residual.py`` by sparse
``python-flint`` arithmetic.  It works in

    QQ(s1,ell,v)[T]/(R20),             v = u^2,

and keeps every coefficient denominator as a product of explicitly recorded
base polynomials.  In particular, no multivariate gcd is taken during the
quotient or border-basis calculation.

The stages are deliberately checkpointed.  ``pivot`` certifies the dense
subresultant pivot, ``moments`` parses the unspecialized corrected moments,
``basis`` constructs the length-four fibre algebra, and ``pencil`` emits the
25 base equations furnished by det(M6+z*M7).

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
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from flint import (
    fmpq,
    fmpq_mpoly,
    fmpq_mpoly_ctx,
    nmod_mpoly_ctx,
)
import sympy as sp

from research_two_pair_sic_bidegree33_t0_Q_residual import (
    LEADING_ARTIFACT,
    ROOT,
    adapted_polynomial,
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
VARIABLES = ("T", "s1", "ell", "v")
CTX = fmpq_mpoly_ctx.get(VARIABLES)
T, S1, ELL, V = CTX.gens()
PRIME = 0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage",
        choices=(
            "prepare",
            "pivot",
            "moments",
            "evaluated",
            "basis",
            "pencil",
        ),
        default="pivot",
    )
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument(
        "--prime",
        type=int,
        default=0,
        help="0 for characteristic zero, otherwise an odd FLINT word prime",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    return parser.parse_args()


def configure_coefficient_ring(prime: int) -> None:
    global CTX, T, S1, ELL, V, PRIME
    PRIME = prime
    if prime:
        CTX = nmod_mpoly_ctx.get(VARIABLES, modulus=prime)
    else:
        CTX = fmpq_mpoly_ctx.get(VARIABLES)
    T, S1, ELL, V = CTX.gens()


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

    s1, ell, root, u = sp.symbols("s1 ell T u")
    polynomial = sp.Poly(
        sp.sympify(
            expression.replace("^", "**"),
            locals={"s1": s1, "ell": ell, "T": root, "u": u},
        ),
        root,
        s1,
        ell,
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
    """Arithmetic in Frac(QQ[s1,ell,v])[T]/(R20)."""

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


def flint_input_data() -> tuple[
    QuotientArithmetic,
    fmpq_mpoly,
    fmpq_mpoly,
    dict[str, object],
]:
    residual_raw, pivot_a_raw, pivot_b_raw = residual_data()
    residual, residual_terms, residual_u_power, residual_content = (
        adapted_polynomial(residual_raw)
    )
    pivot_a, pivot_a_terms, pivot_a_u_power, pivot_a_content = (
        adapted_polynomial(pivot_a_raw)
    )
    pivot_b, pivot_b_terms, pivot_b_u_power, pivot_b_content = (
        adapted_polynomial(pivot_b_raw)
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


def expanded_adapted_moments(
    singular: str,
    timeout: int,
) -> dict[int, str]:
    """Expand mu4,...,mu7 after imposing Q and the weight-zero coordinates."""

    orders = tuple(range(2, 8))
    export = t0_open_localized_export(singular, orders, 0, timeout)
    raw = dict(
        zip(orders[1:], export["polynomials"][:-1], strict=True)
    )
    replacements = (
        ("s2", "(s1^2*u-(13/3)*u)"),
        ("t1", "(u*(s1-ell))"),
        ("t2", "(T*u^2)"),
    )
    declarations = "\n".join(
        f"poly p{order}={substitute(raw[order], replacements)};"
        for order in range(4, 8)
    )
    prints = "\n".join(
        f'print("MOMENT{order} "+string(p{order}));'
        for order in range(4, 8)
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring adapted=0,(s6,s5,s3,T,s1,ell,u),dp;
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
    for order in range(4, 8):
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
    variable_names = {"s6", "s5", "s3", "T", "s1", "ell", "u"}
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
    # have one u parity.  This retains the smaller base ring QQ[s1,ell,u^2].
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
            exponents["ell"],
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
) -> tuple[dict[int, RawFibrePolynomial], dict[str, object]]:
    expanded = expanded_adapted_moments(singular, timeout)
    moments: dict[int, RawFibrePolynomial] = {}
    profiles: dict[str, object] = {}
    for order, expression in expanded.items():
        moments[order], profiles[f"mu{order}"] = (
            parse_raw_fibre_polynomial(expression)
        )
        profiles[f"mu{order}"]["expanded_length"] = len(expression)
    return moments, profiles


def clear_fibre_denominators(
    arithmetic: QuotientArithmetic,
    polynomial: FibrePolynomial,
) -> FibrePolynomial:
    """Multiply a fibre polynomial by one common nonzero base denominator."""

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
            "exact characteristic-zero research calculation on the "
            "unspecialized residual Q component; denominator complements "
            "remain separate until explicitly saturated"
        ),
        "stage": stage,
        "coefficient_ring": (
            "QQ[s1,ell,v], v=u^2"
            if not PRIME
            else f"GF({PRIME})[s1,ell,v], v=u^2"
        ),
        "extension": (
            "QQ(s1,ell,v)[T]/(R20)"
            if not PRIME
            else f"GF({PRIME})(s1,ell,v)[T]/(R20)"
        ),
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
            )
            payload["moments"] = moment_profiles
            if arguments.stage not in ("moments",):
                evaluated = {
                    order: evaluate_moment_fraction_free(
                        arithmetic,
                        moments[order],
                        a_value,
                        b_value,
                    )
                    for order in (4, 5)
                }
                contents = {}
                for order in (4, 5):
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
                    for order in (4, 5)
                }
                if arguments.stage not in ("evaluated",):
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
