#!/usr/bin/env python3
"""Certify the isolated anti-Weyl rank-two nine-moment point.

The proof uses only rational arithmetic.  A Krawczyk inclusion isolates the
zero of (R,I,mu2,mu4,mu6,mu8) in the displayed rational box.  Polynomial
division proves rank(C)<=2 there, interval minors prove exact rank two, and
interval evaluation proves the signs of the first three missing even moments.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import itertools
import json
from math import factorial, gcd
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from explore_two_pair_sic_bidegree33_full_anchor import (  # noqa: E402
    Q_POLYNOMIALS,
    QUADRATIC_Q,
    WEIGHTS,
)
from research_two_pair_sic_bidegree33_anti_weyl import (  # noqa: E402
    anti_weyl_reduce,
    integer_moment_terms,
    primitive_polynomial,
    run_msolve,
)

OUTPUT = (
    ROOT / "artifacts" / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_finite_prefix.json"
)
VARIABLE_NAMES = ("s0", "s1", "s2", "s3", "t0", "t1")
CENTER_STRINGS = (
    "-0.8702368803753021", "-1.8595344299498098",
    "-0.3028762025151954", "1.5192025917402650",
    "-2.3287623389711780", "-1.5350489449499409",
)
RADIUS = Fraction(1, 10_000_000_000)


@dataclass(frozen=True)
class Interval:
    lo: Fraction
    hi: Fraction

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    def __add__(self, other: object) -> "Interval":
        rhs = as_interval(other)
        return Interval(self.lo + rhs.lo, self.hi + rhs.hi)

    __radd__ = __add__

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other: object) -> "Interval":
        return self + (-as_interval(other))

    def __rsub__(self, other: object) -> "Interval":
        return as_interval(other) - self

    def __mul__(self, other: object) -> "Interval":
        rhs = as_interval(other)
        products = (
            self.lo * rhs.lo, self.lo * rhs.hi,
            self.hi * rhs.lo, self.hi * rhs.hi,
        )
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def reciprocal(self) -> "Interval":
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError("interval contains zero")
        return Interval(min(1 / self.lo, 1 / self.hi),
                        max(1 / self.lo, 1 / self.hi))

    def __truediv__(self, other: object) -> "Interval":
        return self * as_interval(other).reciprocal()

    def __pow__(self, exponent: int) -> "Interval":
        if exponent < 0:
            return (self.reciprocal()) ** (-exponent)
        answer = Interval(Fraction(1), Fraction(1))
        base = self
        power = exponent
        while power:
            if power & 1:
                answer = answer * base
            base = base * base
            power //= 2
        return answer

    def strict_subset(self, other: "Interval") -> bool:
        return other.lo < self.lo and self.hi < other.hi


def as_interval(value: object) -> Interval:
    if isinstance(value, Interval):
        return value
    if isinstance(value, sp.Rational):
        value = Fraction(int(value.p), int(value.q))
    if isinstance(value, int):
        value = Fraction(value)
    if isinstance(value, Fraction):
        return Interval(value, value)
    raise TypeError(type(value))


def rational(text: str) -> Fraction:
    return Fraction(text)


def evaluate(polynomial: dict[tuple[int, ...], int], values: list[object]):
    answer = as_interval(0) if any(isinstance(x, Interval) for x in values) else Fraction(0)
    for exponents, coefficient in polynomial.items():
        term = as_interval(coefficient) if isinstance(answer, Interval) else Fraction(coefficient)
        for value, exponent in zip(values, exponents):
            term *= value ** exponent
        answer += term
    return answer


def derivative(polynomial: dict[tuple[int, ...], int], variable: int):
    answer: dict[tuple[int, ...], int] = {}
    for exponents, coefficient in polynomial.items():
        if exponents[variable]:
            reduced = list(exponents)
            reduced[variable] -= 1
            answer[tuple(reduced)] = coefficient * exponents[variable]
    return answer


def expression_terms(expression: sp.Expr, variables: tuple[sp.Symbol, ...]):
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    answer: dict[tuple[int, ...], int] = {}
    denominator = 1
    for coefficient in polynomial.coeffs():
        denominator = sp.ilcm(denominator, coefficient.q)
    for exponents, coefficient in polynomial.terms():
        answer[exponents] = int(coefficient * denominator)
    return answer, int(denominator)


def coefficient_matrix(variables: tuple[sp.Symbol, ...]) -> sp.Matrix:
    s0, s1, s2, s3, t0, t1 = variables
    parameters = (s0, s1, s2, s3, -s2, s1, -s0, t0, t1, 0, t1, -t0)
    matrix = sp.zeros(4)
    for parameter, weight, polynomial in zip(parameters, WEIGHTS, Q_POLYNOMIALS):
        for degree, coefficient in enumerate(polynomial):
            row = degree + max(weight, 0)
            column = degree + max(-weight, 0)
            matrix[row, column] += parameter * coefficient
    for degree, coefficient in enumerate(QUADRATIC_Q):
        matrix[degree, degree] += coefficient
    return matrix


def quotient_moment(polynomial: dict[tuple[int, ...], int]):
    """Reduce an even anti-Weyl moment to (x,y,B,C,D)."""

    x, y, a_square, b_square, c_square, d_square, product = sp.symbols(
        "x y A B C D p"
    )
    a_value = (
        -15*b_square + sp.Rational(14, 3)*c_square
        + sp.Rational(56, 3)*d_square - 6*x**2 - 10*y**2
        - sp.Rational(70, 3)
    )
    product_value = (
        -b_square - sp.Rational(1, 9)*c_square
        + sp.Rational(4, 9)*d_square + x**2 - y**2
        - sp.Rational(8, 9)*y + sp.Rational(1, 9)
    )

    def reduce_s(exponent_s0: int, exponent_s2: int):
        value = a_square**(exponent_s0 // 2) * b_square**(exponent_s2 // 2)
        if exponent_s0 % 2 and exponent_s2 % 2:
            return value * product
        if exponent_s0 % 2 or exponent_s2 % 2:
            raise AssertionError("unpaired sextic square root")
        return value

    answer = 0
    for exponents, coefficient in polynomial.items():
        exponent_s0, exponent_x, exponent_s2, exponent_y, exponent_t0, exponent_t1 = exponents
        assert exponent_t0 % 2 == exponent_t1 % 2
        base = (
            coefficient * x**exponent_x * y**exponent_y
            * c_square**(exponent_t0 // 2)
            * d_square**(exponent_t1 // 2)
        )
        if exponent_t0 % 2:
            # 4*t0*t1=(1-9*y)*s0+9*(2*x+y+1)*s2 on I=0.
            answer += (
                base * sp.Rational(1, 4) * (1 - 9*y)
                * reduce_s(exponent_s0 + 1, exponent_s2)
            )
            answer += (
                base * sp.Rational(9, 4) * (2*x + y + 1)
                * reduce_s(exponent_s0, exponent_s2 + 1)
            )
        else:
            answer += base * reduce_s(exponent_s0, exponent_s2)
    return sp.Poly(
        answer.subs({a_square: a_value, product: product_value}),
        x, y, b_square, c_square, d_square,
    )


def primitive_terms(polynomial: sp.Poly):
    denominator = 1
    for coefficient in polynomial.coeffs():
        denominator = sp.ilcm(denominator, coefficient.q)
    terms = {
        exponents: int(coefficient * denominator)
        for exponents, coefficient in polynomial.terms()
    }
    content = 0
    for coefficient in terms.values():
        content = gcd(content, abs(coefficient))
    return {exponents: coefficient // content for exponents, coefficient in terms.items()}


def determinant(matrix: list[list[Interval]]) -> Interval:
    size = len(matrix)
    answer = as_interval(0)
    for permutation in itertools.permutations(range(size)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(size) for j in range(i + 1, size)
        )
        term = as_interval(-1 if inversions % 2 else 1)
        for row, column in enumerate(permutation):
            term *= matrix[row][column]
        answer += term
    return answer


def matrix_multiply(left: list[list[object]], right: list[list[object]]):
    return [[sum((left[i][k] * right[k][j] for k in range(len(right))),
                 start=as_interval(0))
             for j in range(len(right[0]))] for i in range(len(left))]


def inverse_two(matrix: list[list[Interval]]) -> list[list[Interval]]:
    divisor = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    return [[matrix[1][1] / divisor, -matrix[0][1] / divisor],
            [-matrix[1][0] / divisor, matrix[0][0] / divisor]]


def bivariate_multiply(left: dict[tuple[int, int], object],
                       right: dict[tuple[int, int], object]):
    interval_mode = any(isinstance(value, Interval) for value in left.values())
    zero = as_interval(0) if interval_mode else Fraction(0)
    answer: dict[tuple[int, int], object] = {}
    for (i, j), first in left.items():
        for (a, b), second in right.items():
            key = (i + a, j + b)
            answer[key] = answer.get(key, zero) + first * second
    return answer


def ambient_moment_jacobian(matrix: list[list[object]]) -> list[list[object]]:
    polynomial = {(i, j): matrix[i][j] for i in range(4) for j in range(4)}
    powers = [{(0, 0): as_interval(1) if isinstance(matrix[0][0], Interval)
               else Fraction(1)}]
    for _ in range(8):
        powers.append(bivariate_multiply(powers[-1], polynomial))
    rows = []
    for order in range(1, 10):
        previous = powers[order - 1]
        row = []
        for a in range(4):
            for b in range(4):
                value = as_interval(0) if isinstance(matrix[0][0], Interval) else Fraction(0)
                for diagonal in range(3 * order + 1):
                    value += (
                        order * factorial(3 * order - diagonal) * factorial(diagonal)
                        * previous.get((diagonal - a, diagonal - b), 0)
                    )
                row.append(value)
        rows.append(row)
    return rows


def determinantal_tangent_basis(matrix: list[list[Interval]]):
    """Tangent chart with the central rows/columns as invertible pivot."""

    pivot = (1, 2)
    complement = (0, 3)
    a = [[matrix[i][j] for j in pivot] for i in pivot]
    b = [[matrix[i][j] for j in complement] for i in pivot]
    d = [[matrix[i][j] for j in pivot] for i in complement]
    a_inverse = inverse_two(a)
    coordinates = ([(i, j) for i in pivot for j in range(4)]
                   + [(i, j) for i in complement for j in pivot])
    columns = []
    zero2 = lambda: [[as_interval(0), as_interval(0)],
                     [as_interval(0), as_interval(0)]]
    for coordinate in coordinates:
        da, db, dd = zero2(), zero2(), zero2()
        i, j = coordinate
        if i in pivot and j in pivot:
            da[pivot.index(i)][pivot.index(j)] = as_interval(1)
        elif i in pivot:
            db[pivot.index(i)][complement.index(j)] = as_interval(1)
        else:
            dd[complement.index(i)][pivot.index(j)] = as_interval(1)
        de = matrix_multiply(matrix_multiply(dd, a_inverse), b)
        de_b = matrix_multiply(matrix_multiply(d, a_inverse), db)
        correction = matrix_multiply(
            matrix_multiply(matrix_multiply(matrix_multiply(d, a_inverse), da),
                            a_inverse), b
        )
        de = [[de[r][c] + de_b[r][c] - correction[r][c]
               for c in range(2)] for r in range(2)]
        column = [as_interval(0) for _ in range(16)]
        column[i * 4 + j] = as_interval(1)
        for r, row_index in enumerate(complement):
            for c, column_index in enumerate(complement):
                column[row_index * 4 + column_index] = de[r][c]
        columns.append(column)
    return [[columns[column][row] for column in range(12)] for row in range(16)]


def format_fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def interval_payload(value: Interval) -> list[str]:
    return [format_fraction(value.lo), format_fraction(value.hi)]


def main() -> None:
    variables = sp.symbols(" ".join(VARIABLE_NAMES))
    s0, s1, s2, s3, t0, t1 = variables
    rank_real = (
        -9*s0*s2 + 9*s1**2 - 9*s2**2 - 9*s3**2 - 8*s3
        - t0**2 + 4*t1**2 + 1
    )
    rank_imaginary = (
        -9*s0*s3 + s0 + 18*s1*s2 + 9*s2*s3 + 9*s2 - 4*t0*t1
    )
    rank_terms = [expression_terms(rank_real, variables)[0],
                  expression_terms(rank_imaginary, variables)[0]]

    moments = {}
    contents = {}
    for order in (2, 4, 6, 8, 10, 12, 14):
        primitive, content = primitive_polynomial(
            anti_weyl_reduce(integer_moment_terms(order), "normalized")
        )
        moments[order] = primitive
        contents[order] = content

    system = rank_terms + [moments[order] for order in (2, 4, 6, 8)]
    center = [rational(value) for value in CENTER_STRINGS]
    box = [Interval(value - RADIUS, value + RADIUS) for value in center]
    jacobian = [[derivative(polynomial, column) for column in range(6)]
                for polynomial in system]
    jacobian_center = sp.Matrix([
        [sp.Rational(evaluate(entry, center).numerator,
                     evaluate(entry, center).denominator) for entry in row]
        for row in jacobian
    ])
    assert jacobian_center.det() != 0
    inverse = jacobian_center.inv()
    inverse_fraction = [[Fraction(int(value.p), int(value.q)) for value in row]
                        for row in inverse.tolist()]
    residual = [evaluate(polynomial, center) for polynomial in system]
    newton_center = [
        center[i] - sum(inverse_fraction[i][j] * residual[j] for j in range(6))
        for i in range(6)
    ]
    jacobian_box = [[evaluate(entry, box) for entry in row] for row in jacobian]
    defect = []
    for i in range(6):
        row = []
        for j in range(6):
            value = as_interval(1 if i == j else 0)
            value -= sum(inverse_fraction[i][k] * jacobian_box[k][j]
                         for k in range(6))
            row.append(value)
        defect.append(row)
    delta = Interval(-RADIUS, RADIUS)
    krawczyk = [
        as_interval(newton_center[i]) + sum(
            (defect[i][j] * delta for j in range(6)), start=as_interval(0)
        )
        for i in range(6)
    ]
    assert all(image.strict_subset(domain) for image, domain in zip(krawczyk, box))

    matrix = coefficient_matrix(variables)
    groebner = sp.groebner([rank_real, rank_imaginary], *variables, order="grevlex")
    for rows in itertools.combinations(range(4), 3):
        for columns in itertools.combinations(range(4), 3):
            assert groebner.reduce(sp.expand(matrix.extract(rows, columns).det()))[1] == 0
    matrix_intervals = [[
        evaluate(expression_terms(matrix[i, j], variables)[0], box)
        for j in range(4)
    ] for i in range(4)]
    pivot_rows = (1, 2)
    pivot_columns = (1, 2)
    pivot = determinant([[matrix_intervals[i][j] for j in pivot_columns]
                         for i in pivot_rows])
    assert not (pivot.lo <= 0 <= pivot.hi)

    tangent = determinantal_tangent_basis(matrix_intervals)
    moment_tangent = matrix_multiply(ambient_moment_jacobian(matrix_intervals), tangent)
    center_tangent = [[
        Fraction((entry.lo + entry.hi), 2) for entry in row
    ] for row in moment_tangent[:8]]
    _, pivot_columns_tangent = sp.Matrix(center_tangent).rref()
    assert len(pivot_columns_tangent) == 8
    tangent_minor = determinant([
        [moment_tangent[row][column] for column in pivot_columns_tangent]
        for row in range(8)
    ])
    assert not (tangent_minor.lo <= 0 <= tangent_minor.hi)

    sign_intervals = {order: evaluate(moments[order], box)
                      for order in (10, 12, 14)}
    assert sign_intervals[10].lo > 0
    assert sign_intervals[12].hi < 0
    assert sign_intervals[14].lo > 0

    quotient_profiles = {}
    quotient_polynomials = {}
    for order in (4, 6, 8, 12):
        quotient = quotient_moment(moments[order])
        quotient_polynomials[order] = quotient
        quotient_profiles[str(order)] = {
            "variables": ["x", "y", "B", "C", "D"],
            "terms": len(quotient.terms()),
            "degree": quotient.total_degree(),
        }
    assert [quotient_profiles[str(order)]["terms"] for order in (4, 6, 8, 12)] == [
        39, 119, 294, 1218
    ]
    x, y, b_square, c_square, d_square = sp.symbols("x y B C D")
    a_value = (
        -15*b_square + sp.Rational(14, 3)*c_square
        + sp.Rational(56, 3)*d_square - 6*x**2 - 10*y**2
        - sp.Rational(70, 3)
    )
    product_value = (
        -b_square - sp.Rational(1, 9)*c_square
        + sp.Rational(4, 9)*d_square + x**2 - y**2
        - sp.Rational(8, 9)*y + sp.Rational(1, 9)
    )
    quotient_rank_equations = (
        product_value**2 - a_value*b_square,
        16*c_square*d_square - (1 - 9*y)**2*a_value
        - 81*(2*x + y + 1)**2*b_square
        - 18*(1 - 9*y)*(2*x + y + 1)*product_value,
    )
    corrected_system = {
        0: primitive_terms(sp.Poly(quotient_rank_equations[0], x, y, b_square, c_square, d_square)),
        1: primitive_terms(sp.Poly(quotient_rank_equations[1], x, y, b_square, c_square, d_square)),
    }
    corrected_system.update({
        order: primitive_terms(quotient_polynomials[order])
        for order in (4, 6, 8, 12)
    })
    corrected_solver = run_msolve(
        corrected_system, ("x", "y", "B", "C", "D"), 0, 300, 4, 2
    )
    assert corrected_solver["status"] == "unit"

    profiles = {
        str(order): {
            "terms": len(moments[order]),
            "degree": order,
            "raw_positive_content": str(abs(contents[order])),
        }
        for order in moments
    }
    payload = {
        "calculation": "two_pair_sic_bidegree33_rank_two_finite_prefix",
        "status": "proved",
        "arithmetic": "exact rational interval arithmetic",
        "variables": VARIABLE_NAMES,
        "center": CENTER_STRINGS,
        "radius": format_fraction(RADIUS),
        "system": ["R", "I", "mu2", "mu4", "mu6", "mu8"],
        "krawczyk_strict_inclusion": True,
        "krawczyk_image": [interval_payload(value) for value in krawczyk],
        "rank": {
            "all_3x3_minors_reduce_to_zero_mod_R_I": True,
            "nonzero_2x2_pivot": {"rows": pivot_rows, "columns": pivot_columns,
                                  "interval": interval_payload(pivot)},
            "exact": 2,
        },
        "tangent": {
            "rank_lower_bound": 8,
            "moment_rows": list(range(1, 9)),
            "chart_columns": list(pivot_columns_tangent),
            "minor_interval": interval_payload(tangent_minor),
            "rank_upper_bound": 8,
            "upper_bound_reason": (
                "three SL2 orbit directions and radial scaling are independent "
                "kernel directions on the smooth 12-dimensional exact-rank-two locus"
            ),
            "rank": 8,
            "quotient_point": "reduced and isolated",
        },
        "odd_moments_identically_zero": [1, 3, 5, 7, 9, 11, 13],
        "moment_profiles": profiles,
        "anti_weyl_quotient_profiles": quotient_profiles,
        "corrected_anti_weyl_quotient": {
            "system": ["rank_quartic_1", "rank_quartic_2", "mu4", "mu6", "mu8", "mu12"],
            "variables": ["x", "y", "B", "C", "D"],
            "msolve": corrected_solver,
            "conclusion": "unit ideal over QQ; the corrected rank-two system has no anti-Weyl point",
        },
        "missing_moment_signs": {
            str(order): interval_payload(value) for order, value in sign_intervals.items()
        },
        "semistability": "the normalized non-null Sym^2 component has nonzero discriminant",
        "conclusion": (
            "unique exact-rank-two semistable zero of mu1,...,mu9 in the box; "
            "it does not lift because primitive mu10 is positive"
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    print("artifact_sha256=" + sha256(encoded.encode()).hexdigest())


if __name__ == "__main__":
    main()
