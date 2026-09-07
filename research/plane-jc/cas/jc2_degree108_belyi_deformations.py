#!/usr/bin/env python3
"""Exact dessin-first deformation compiler for the degree-(72,108) JC2 row.

This module deliberately does not eliminate the original bivariate coefficient
system.  It consumes the already certified intrinsic quintic graph of the top
Wronskian block, reconstructs its normalized degree-21 Belyi map, and compiles
the four remaining Wronskian equations as successive linear maps.

All coefficient arithmetic is performed in the exact etale algebra

    K = QQ[q]/(H(q))

using ``fractions.Fraction``.  The implementation is dependency-free except
for the optional SymPy arithmetic-monodromy audit in ``compile_report``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Iterable, Sequence

from verify_firstblock_quotient_graph import GRAPH, H


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts" / "generated-results" / "jc2_degree108_belyi_deformations.json"
GRAPH_FILE = Path(__file__).resolve().parent / "firstblock_mu7_quotient_lex_basis.txt"
SINGULAR = "Singular"


# ---------------------------------------------------------------------------
# The exact quintic coefficient algebra.


@dataclass(frozen=True)
class KElement:
    """An element of QQ[q]/(H), stored in the basis 1,q,...,q^4."""

    coeffs: tuple[Fraction, Fraction, Fraction, Fraction, Fraction]

    def __add__(self, other: object) -> "KElement":
        rhs = as_k(other)
        return KElement(tuple(a + b for a, b in zip(self.coeffs, rhs.coeffs)))

    __radd__ = __add__

    def __neg__(self) -> "KElement":
        return KElement(tuple(-a for a in self.coeffs))

    def __sub__(self, other: object) -> "KElement":
        return self + (-as_k(other))

    def __rsub__(self, other: object) -> "KElement":
        return as_k(other) - self

    def __mul__(self, other: object) -> "KElement":
        rhs = as_k(other)
        raw = [Fraction(0) for _ in range(9)]
        for i, a in enumerate(self.coeffs):
            if a:
                for j, b in enumerate(rhs.coeffs):
                    if b:
                        raw[i + j] += a * b
        for degree in range(8, 4, -1):
            lead = raw[degree]
            if not lead:
                continue
            raw[degree] = Fraction(0)
            shift = degree - 5
            for i in range(5):
                raw[shift + i] -= lead * Fraction(H[i], H[5])
        return KElement(tuple(raw[:5]))

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "KElement":
        if exponent < 0:
            return (self.inverse()) ** (-exponent)
        result = K_ONE
        base = self
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent >>= 1
        return result

    def inverse(self) -> "KElement":
        if not self:
            raise ZeroDivisionError("zero is not invertible in K")
        # Extended Euclid in QQ[q].  Lists are ascending coefficient vectors.
        old_r, r = trim(list(self.coeffs)), trim([Fraction(v) for v in H])
        old_s, s = [Fraction(1)], []
        while r:
            quotient, remainder = qq_divmod(old_r, r)
            old_r, r = r, remainder
            old_s, s = s, qq_sub(old_s, qq_mul(quotient, s))
        if len(old_r) != 1:
            raise ZeroDivisionError("element is a zero divisor in K")
        scalar = old_r[0]
        representative = [value / scalar for value in old_s]
        return poly_reduce_k(representative)

    def __truediv__(self, other: object) -> "KElement":
        return self * as_k(other).inverse()

    def __rtruediv__(self, other: object) -> "KElement":
        return as_k(other) / self

    def __bool__(self) -> bool:
        return any(self.coeffs)

    def text(self, variable: str = "q") -> str:
        terms: list[str] = []
        for degree, value in enumerate(self.coeffs):
            if not value:
                continue
            coefficient = fraction_text(value)
            if degree == 0:
                terms.append(f"({coefficient})")
            else:
                terms.append(f"({coefficient})*{variable}^{degree}")
        return "+".join(terms) if terms else "0"


def as_k(value: object) -> KElement:
    if isinstance(value, KElement):
        return value
    if isinstance(value, Fraction):
        scalar = value
    elif isinstance(value, int):
        scalar = Fraction(value)
    else:
        raise TypeError(f"cannot coerce {type(value)!r} to K")
    return KElement((scalar, Fraction(0), Fraction(0), Fraction(0), Fraction(0)))


K_ZERO = as_k(0)
K_ONE = as_k(1)
K_Q = KElement((Fraction(0), Fraction(1), Fraction(0), Fraction(0), Fraction(0)))


def trim(values: list) -> list:
    while values and not values[-1]:
        values.pop()
    return values


def qq_sub(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    size = max(len(left), len(right))
    out = [Fraction(0) for _ in range(size)]
    for i in range(size):
        out[i] = (left[i] if i < len(left) else 0) - (right[i] if i < len(right) else 0)
    return trim(out)


def qq_mul(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    if not left or not right:
        return []
    out = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim(out)


def qq_divmod(
    numerator: list[Fraction], denominator: list[Fraction]
) -> tuple[list[Fraction], list[Fraction]]:
    if not denominator:
        raise ZeroDivisionError
    remainder = numerator[:]
    quotient = [Fraction(0) for _ in range(max(0, len(numerator) - len(denominator) + 1))]
    while len(remainder) >= len(denominator) and remainder:
        shift = len(remainder) - len(denominator)
        coefficient = remainder[-1] / denominator[-1]
        quotient[shift] += coefficient
        for i, value in enumerate(denominator):
            remainder[shift + i] -= coefficient * value
        trim(remainder)
    return trim(quotient), remainder


def poly_reduce_k(values: Sequence[Fraction]) -> KElement:
    raw = list(values) + [Fraction(0)] * max(0, 5 - len(values))
    for degree in range(len(raw) - 1, 4, -1):
        lead = raw[degree]
        if lead:
            shift = degree - 5
            for i in range(5):
                raw[shift + i] -= lead * Fraction(H[i], H[5])
    return KElement(tuple(raw[:5]))


def graph_element(name: str) -> KElement:
    denominator, numerators = GRAPH[name]
    value = K_ZERO
    for degree, numerator in enumerate(numerators, 1):
        value += Fraction(numerator, denominator) * K_Q**degree
    return value


# ---------------------------------------------------------------------------
# Univariate polynomials over K, represented by ascending coefficient lists.


KPoly = list[KElement]


def kp(values: Iterable[object]) -> KPoly:
    return ktrim([as_k(value) for value in values])


def ktrim(values: KPoly) -> KPoly:
    while values and not values[-1]:
        values.pop()
    return values


def kadd(left: KPoly, right: KPoly) -> KPoly:
    size = max(len(left), len(right))
    out = [K_ZERO for _ in range(size)]
    for i in range(size):
        out[i] = (left[i] if i < len(left) else K_ZERO) + (right[i] if i < len(right) else K_ZERO)
    return ktrim(out)


def kscale(polynomial: KPoly, scalar: object) -> KPoly:
    factor = as_k(scalar)
    return ktrim([factor * value for value in polynomial])


def kmul(left: KPoly, right: KPoly) -> KPoly:
    if not left or not right:
        return []
    out = [K_ZERO for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    out[i + j] = out[i + j] + a * b
    return ktrim(out)


def kderivative(polynomial: KPoly) -> KPoly:
    return ktrim([i * polynomial[i] for i in range(1, len(polynomial))])


def kshift(polynomial: KPoly, degree: int) -> KPoly:
    return ([K_ZERO] * degree + polynomial) if polynomial else []


def kdivmod(numerator: KPoly, denominator: KPoly) -> tuple[KPoly, KPoly]:
    if not denominator:
        raise ZeroDivisionError
    remainder = numerator[:]
    quotient = [K_ZERO for _ in range(max(0, len(numerator) - len(denominator) + 1))]
    inverse_lead = denominator[-1].inverse()
    while len(remainder) >= len(denominator) and remainder:
        shift = len(remainder) - len(denominator)
        coefficient = remainder[-1] * inverse_lead
        quotient[shift] = quotient[shift] + coefficient
        for i, value in enumerate(denominator):
            remainder[shift + i] = remainder[shift + i] - coefficient * value
        ktrim(remainder)
    return ktrim(quotient), remainder


def kgcd(left: KPoly, right: KPoly) -> KPoly:
    a, b = left[:], right[:]
    while b:
        _, remainder = kdivmod(a, b)
        a, b = b, remainder
    return kscale(a, a[-1].inverse()) if a else []


def coefficient(polynomial: KPoly, degree: int) -> KElement:
    return polynomial[degree] if degree < len(polynomial) else K_ZERO


# ---------------------------------------------------------------------------
# Multivariate parameter polynomials over K.


@dataclass(frozen=True)
class ParameterPolynomial:
    terms: dict[tuple[int, ...], KElement]
    variables: int

    @staticmethod
    def constant(value: object, variables: int) -> "ParameterPolynomial":
        scalar = as_k(value)
        return ParameterPolynomial({(0,) * variables: scalar} if scalar else {}, variables)

    @staticmethod
    def variable(index: int, variables: int) -> "ParameterPolynomial":
        exponent = [0] * variables
        exponent[index] = 1
        return ParameterPolynomial({tuple(exponent): K_ONE}, variables)

    def __add__(self, other: object) -> "ParameterPolynomial":
        rhs = as_parameter(other, self.variables)
        out = dict(self.terms)
        for monomial, value in rhs.terms.items():
            new_value = out.get(monomial, K_ZERO) + value
            if new_value:
                out[monomial] = new_value
            elif monomial in out:
                del out[monomial]
        return ParameterPolynomial(out, self.variables)

    __radd__ = __add__

    def __neg__(self) -> "ParameterPolynomial":
        return ParameterPolynomial({m: -c for m, c in self.terms.items()}, self.variables)

    def __sub__(self, other: object) -> "ParameterPolynomial":
        return self + (-as_parameter(other, self.variables))

    def __rsub__(self, other: object) -> "ParameterPolynomial":
        return as_parameter(other, self.variables) - self

    def __mul__(self, other: object) -> "ParameterPolynomial":
        rhs = as_parameter(other, self.variables)
        out: dict[tuple[int, ...], KElement] = {}
        for left_monomial, left_value in self.terms.items():
            for right_monomial, right_value in rhs.terms.items():
                monomial = tuple(a + b for a, b in zip(left_monomial, right_monomial))
                value = out.get(monomial, K_ZERO) + left_value * right_value
                if value:
                    out[monomial] = value
                elif monomial in out:
                    del out[monomial]
        return ParameterPolynomial(out, self.variables)

    __rmul__ = __mul__

    def __bool__(self) -> bool:
        return bool(self.terms)

    def total_degree(self) -> int:
        return max((sum(monomial) for monomial in self.terms), default=-1)

    def text(self, names: Sequence[str], field_variable: str = "q") -> str:
        if len(names) != self.variables:
            raise ValueError("parameter-name count mismatch")
        pieces = []
        for monomial in sorted(self.terms, key=lambda item: (sum(item), item)):
            factors = [f"({self.terms[monomial].text(field_variable)})"]
            for name, exponent in zip(names, monomial):
                if exponent:
                    factors.append(name if exponent == 1 else f"{name}^{exponent}")
            pieces.append("*".join(factors))
        return "+".join(pieces) if pieces else "0"


def as_parameter(value: object, variables: int) -> ParameterPolynomial:
    if isinstance(value, ParameterPolynomial):
        if value.variables != variables:
            raise ValueError("parameter-ring mismatch")
        return value
    return ParameterPolynomial.constant(value, variables)


PPoly = list[ParameterPolynomial]


def ptrim(values: PPoly) -> PPoly:
    while values and not values[-1]:
        values.pop()
    return values


def padd(left: PPoly, right: PPoly) -> PPoly:
    variables = left[0].variables if left else right[0].variables
    size = max(len(left), len(right))
    zero = ParameterPolynomial.constant(0, variables)
    return ptrim([
        (left[i] if i < len(left) else zero) + (right[i] if i < len(right) else zero)
        for i in range(size)
    ])


def pscale(polynomial: PPoly, scalar: object) -> PPoly:
    return ptrim([value * scalar for value in polynomial])


def pmul(left: PPoly, right: PPoly) -> PPoly:
    if not left or not right:
        return []
    variables = left[0].variables
    zero = ParameterPolynomial.constant(0, variables)
    out = [zero for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    out[i + j] = out[i + j] + a * b
    return ptrim(out)


def pderivative(polynomial: PPoly) -> PPoly:
    return ptrim([i * polynomial[i] for i in range(1, len(polynomial))])


def plift(polynomial: KPoly, variables: int) -> PPoly:
    return [ParameterPolynomial.constant(value, variables) for value in polynomial]


# ---------------------------------------------------------------------------
# Exact linear algebra over K.


def rref_with_transform(matrix: list[list[KElement]]) -> tuple[list[list[KElement]], list[int], list[list[KElement]]]:
    rows = len(matrix)
    columns = len(matrix[0]) if rows else 0
    work = [row[:] for row in matrix]
    transform = [[K_ONE if i == j else K_ZERO for j in range(rows)] for i in range(rows)]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(columns):
        selected = next((row for row in range(pivot_row, rows) if work[row][column]), None)
        if selected is None:
            continue
        work[pivot_row], work[selected] = work[selected], work[pivot_row]
        transform[pivot_row], transform[selected] = transform[selected], transform[pivot_row]
        inverse = work[pivot_row][column].inverse()
        work[pivot_row] = [inverse * value for value in work[pivot_row]]
        transform[pivot_row] = [inverse * value for value in transform[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [a - factor * b for a, b in zip(work[row], work[pivot_row])]
            transform[row] = [a - factor * b for a, b in zip(transform[row], transform[pivot_row])]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == rows:
            break
    return work, pivots, transform


def kernel_from_rref(rref: list[list[KElement]], pivots: list[int]) -> list[list[KElement]]:
    columns = len(rref[0]) if rref else 0
    free = [column for column in range(columns) if column not in pivots]
    basis = []
    for free_column in free:
        vector = [K_ZERO for _ in range(columns)]
        vector[free_column] = K_ONE
        for row, pivot in enumerate(pivots):
            vector[pivot] = -rref[row][free_column]
        basis.append(vector)
    return basis


def transform_rhs(
    transform: list[list[KElement]], rhs: list[ParameterPolynomial]
) -> list[ParameterPolynomial]:
    variables = rhs[0].variables
    zero = ParameterPolynomial.constant(0, variables)
    return [sum((value * coefficient for coefficient, value in zip(row, rhs)), zero) for row in transform]


def solve_with_parameters(
    rref: list[list[KElement]],
    pivots: list[int],
    transform: list[list[KElement]],
    rhs: list[ParameterPolynomial],
    free_parameters: Sequence[ParameterPolynomial],
) -> tuple[list[ParameterPolynomial], list[ParameterPolynomial]]:
    transformed = transform_rhs(transform, rhs)
    columns = len(rref[0]) if rref else 0
    free_columns = [column for column in range(columns) if column not in pivots]
    if len(free_columns) != len(free_parameters):
        raise ValueError("free-parameter count mismatch")
    zero = ParameterPolynomial.constant(0, rhs[0].variables)
    solution = [zero for _ in range(columns)]
    for column, parameter in zip(free_columns, free_parameters):
        solution[column] = parameter
    for row, pivot in enumerate(pivots):
        solution[pivot] = transformed[row] - sum(
            (solution[column] * rref[row][column] for column in free_columns), zero
        )
    compatibility = [value for value in transformed[len(pivots):] if value]
    return solution, compatibility


def column_polynomial(
    fixed: dict[str, KPoly],
    variable: str,
    degree: int,
    stage: str,
) -> KPoly:
    monomial = [K_ZERO] * degree + [K_ONE]
    A, D = fixed["A"], fixed["D"]
    if stage == "BE":
        if variable == "B":
            return kadd(kscale(kmul(monomial, kderivative(D)), -1), kscale(kmul(D, kderivative(monomial)), 3))
        return kadd(kscale(kmul(A, kderivative(monomial)), -2), kscale(kmul(monomial, kderivative(A)), 2))
    if stage == "CF":
        if variable == "C":
            return kscale(kmul(D, kderivative(monomial)), 3)
        return kadd(kscale(kmul(A, kderivative(monomial)), -2), kmul(monomial, kderivative(A)))
    if stage == "G":
        return kscale(kmul(A, kderivative(monomial)), -2)
    raise ValueError(stage)


def linear_matrix(
    fixed: dict[str, KPoly], columns: Sequence[tuple[str, int]], stage: str, target_degree: int = 19
) -> list[list[KElement]]:
    images = [column_polynomial(fixed, variable, degree, stage) for variable, degree in columns]
    return [[coefficient(image, row) for image in images] for row in range(target_degree + 1)]


def vector_to_pair(
    values: Sequence[ParameterPolynomial], columns: Sequence[tuple[str, int]], names: tuple[str, str]
) -> tuple[PPoly, PPoly]:
    variables = values[0].variables
    zero = ParameterPolynomial.constant(0, variables)
    degrees = {name: max(degree for variable, degree in columns if variable == name) for name in names}
    out = {name: [zero for _ in range(degrees[name] + 1)] for name in names}
    for value, (variable, degree) in zip(values, columns):
        out[variable][degree] = value
    return ptrim(out[names[0]]), ptrim(out[names[1]])


# ---------------------------------------------------------------------------
# Reconstruction, dessin enumeration, and the deformation cascade.


def reconstruct_top() -> tuple[KPoly, KPoly]:
    x = {index: graph_element(f"x{index}") for index in range(2, 7)}
    A = [K_ZERO, K_ONE] + [x[index] for index in range(2, 7)] + [K_Q, K_Q]
    D = [K_ZERO for _ in range(13)]
    for target_degree in range(2, 13):
        known = K_ZERO
        for i in range(2, min(8, target_degree - 1) + 1):
            j = target_degree + 1 - i
            if 2 <= j <= 12:
                known += (2 * j - 3 * i) * A[i] * D[j]
        rhs = K_ONE if target_degree == 2 else K_ZERO
        D[target_degree] = (rhs - known) / (2 * target_degree - 3)
    top = kadd(kscale(kmul(A, kderivative(D)), 2), kscale(kmul(kderivative(A), D), -3))
    assert top == [K_ZERO, K_ZERO, K_ONE]
    return A, ktrim(D)


def belyi_audit(A: KPoly, D: KPoly) -> dict:
    assert len(A) == 9 and len(D) == 13 and A[0] == K_ZERO and D[:2] == [K_ZERO, K_ZERO]
    a = A[1:]
    d = D[2:]
    scale = d[-1] ** 2 / a[-1] ** 3
    third_fiber = kadd(kshift(kmul(d, d), 1), kscale(kmul(kmul(a, a), a), -scale))
    assert len(third_fiber) == 5
    assert kgcd(a, d) == [K_ONE]
    assert kgcd(a, kderivative(a)) == [K_ONE]
    assert kgcd(d, kderivative(d)) == [K_ONE]
    assert kgcd(third_fiber, kderivative(third_fiber)) == [K_ONE]
    assert kgcd(a, third_fiber) == [K_ONE]
    assert kgcd(d, third_fiber) == [K_ONE]
    return {
        "normalized_map": "beta=X*(D/X^2)^2/(c*(A/X)^3)",
        "c": scale.text(),
        "third_fiber_degree": len(third_fiber) - 1,
        "infinity_multiplicity": 17,
        "squarefree_finite_fibers": True,
        "passport": [[2] * 10 + [1], [3] * 7, [17, 1, 1, 1, 1]],
        "third_fiber_hash": hash_kpolys([third_fiber]),
        "third_fiber_coefficients_ascending": [value.text() for value in third_fiber],
    }


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[i]] for i in range(len(left)))


def inverse_permutation(permutation: tuple[int, ...]) -> tuple[int, ...]:
    inverse = [0] * len(permutation)
    for index, image in enumerate(permutation):
        inverse[image] = index
    return tuple(inverse)


def cycle_type(permutation: tuple[int, ...]) -> tuple[int, ...]:
    seen = [False] * len(permutation)
    cycles = []
    for start in range(len(permutation)):
        if seen[start]:
            continue
        point, length = start, 0
        while not seen[point]:
            seen[point] = True
            length += 1
            point = permutation[point]
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))


def permutation_orbit(generators: Sequence[tuple[int, ...]], start: int = 0) -> set[int]:
    orbit, frontier = {start}, [start]
    while frontier:
        point = frontier.pop()
        for generator in generators:
            image = generator[point]
            if image not in orbit:
                orbit.add(image)
                frontier.append(image)
    return orbit


def dessin_automorphism_order(generators: Sequence[tuple[int, ...]]) -> int:
    """Return the centralizer order of a transitive permutation tuple."""

    degree = len(generators[0])
    order = 0
    for image_of_zero in range(degree):
        image = {0: image_of_zero}
        frontier = [0]
        valid = True
        while frontier and valid:
            point = frontier.pop()
            for generator in generators:
                source_image = generator[point]
                target_image = generator[image[point]]
                if source_image in image:
                    valid = image[source_image] == target_image
                    if not valid:
                        break
                else:
                    image[source_image] = target_image
                    frontier.append(source_image)
        if valid and len(image) == degree and len(set(image.values())) == degree:
            order += 1
    return order


def canonical_rotation(center_set: Sequence[int]) -> tuple[int, ...]:
    return min(tuple(sorted((value + rotation) % 17 for value in center_set)) for rotation in range(17))


def enumerate_dessins() -> list[dict]:
    degree = 21
    sigma_infinity = tuple(list(range(1, 17)) + [0, 17, 18, 19, 20])
    representatives: dict[tuple[int, ...], tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]] = {}
    labelled = 0
    for centers in itertools.combinations(range(17), 4):
        forced: dict[int, int] = {}
        valid = True
        for fixed_point, center in zip(range(17, 21), centers):
            for left, right in ((fixed_point, center), ((center + 1) % 17, (center - 1) % 17)):
                if left == right or left in forced or right in forced:
                    valid = False
                    break
                forced[left], forced[right] = right, left
            if not valid:
                break
        if not valid:
            continue
        remaining = [value for value in range(17) if value not in forced]
        if len(remaining) != 5:
            continue
        found: list[tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]] = []
        for fixed_zero in remaining:
            rest = [value for value in remaining if value != fixed_zero]
            a, b, c, d = rest
            for pairs in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))):
                sigma_zero = list(range(degree))
                for left, right in forced.items():
                    sigma_zero[left] = right
                for left, right in pairs:
                    sigma_zero[left], sigma_zero[right] = right, left
                sigma_zero = tuple(sigma_zero)
                sigma_one = compose(inverse_permutation(sigma_zero), inverse_permutation(sigma_infinity))
                if cycle_type(sigma_zero) == (2,) * 10 + (1,) and cycle_type(sigma_one) == (3,) * 7:
                    found.append((sigma_zero, sigma_one, sigma_infinity))
        if found:
            assert len(found) == 1
            labelled += len(found)
            key = canonical_rotation(centers)
            representatives.setdefault(key, found[0])
    expected = {(0, 3, 7, 11), (0, 3, 7, 12), (0, 3, 8, 11), (0, 3, 8, 13), (0, 3, 9, 13)}
    assert labelled == 85 and set(representatives) == expected
    result = []
    for center_set in sorted(representatives):
        sigma_zero, sigma_one, sigma_infinity = representatives[center_set]
        assert cycle_type(sigma_zero) == (2,) * 10 + (1,)
        assert cycle_type(sigma_one) == (3,) * 7
        assert cycle_type(sigma_infinity) == (17, 1, 1, 1, 1)
        assert compose(compose(sigma_zero, sigma_one), sigma_infinity) == tuple(range(degree))
        assert len(permutation_orbit((sigma_zero, sigma_one))) == degree
        automorphism_order = dessin_automorphism_order((sigma_zero, sigma_one, sigma_infinity))
        assert automorphism_order == 1
        result.append({
            "center_set": list(center_set),
            "sigma_0": list(sigma_zero),
            "sigma_1": list(sigma_one),
            "sigma_infinity": list(sigma_infinity),
            "transitive": True,
            "automorphism_order": automorphism_order,
        })
    return result


def compile_deformations(A: KPoly, D: KPoly) -> tuple[dict, list[ParameterPolynomial], ParameterPolynomial]:
    fixed = {"A": A, "D": D}
    parameter_count = 5
    parameters = [ParameterPolynomial.variable(index, parameter_count) for index in range(parameter_count)]
    zero = ParameterPolynomial.constant(0, parameter_count)

    be_columns = [("B", degree) for degree in range(1, 9)] + [("E", degree) for degree in range(2, 13)]
    be_matrix = linear_matrix(fixed, be_columns, "BE")
    be_rref, be_pivots, be_transform = rref_with_transform(be_matrix)
    be_kernel = kernel_from_rref(be_rref, be_pivots)
    assert len(be_kernel) == 2

    # Closed-form kernel basis from the density-preserving source modes.
    a = A[1:]
    d = D[2:]
    phi_modes = [kadd(kmul(a, a), [as_k(-1)]), [K_ZERO, K_ONE]]
    analytic_be = []
    for phi in phi_modes:
        B_mode = kadd(
            kmul(phi, kadd(kderivative(A), kscale(a, -1))),
            kscale(kmul(kderivative(phi), A), Fraction(-1, 2)),
        )
        E_mode = kadd(
            kmul(phi, kadd(kderivative(D), kscale(D[1:], Fraction(-3, 2)))),
            kscale(kmul(kderivative(phi), D), Fraction(-3, 4)),
        )
        assert not column_combination(be_matrix, be_columns, {"B": B_mode, "E": E_mode})
        analytic_be.append((B_mode, E_mode))
    analytic_be_vectors = [pair_vector(be_columns, {"B": B_mode, "E": E_mode}) for B_mode, E_mode in analytic_be]
    assert len(rref_with_transform([list(row) for row in zip(*analytic_be_vectors)])[1]) == 2
    be_vector = [parameters[0] * basis[0] + parameters[1] * basis[1] for basis in zip(*be_kernel)]
    B, E = vector_to_pair(be_vector, be_columns, ("B", "E"))

    cf_columns = [("C", degree) for degree in range(0, 9)] + [("F", degree) for degree in range(1, 13)]
    cf_matrix = linear_matrix(fixed, cf_columns, "CF")
    cf_rref, cf_pivots, cf_transform = rref_with_transform(cf_matrix)
    cf_kernel = kernel_from_rref(cf_rref, cf_pivots)
    assert len(cf_kernel) == 3

    analytic_cf = []
    for f_mode in (d, [K_ONE], [K_ZERO, K_ONE]):
        C_mode = kadd(
            kadd(kmul(f_mode, kderivative(A)), kscale(kmul(kderivative(f_mode), A), Fraction(-2, 3))),
            kscale(kmul(f_mode, a), Fraction(-4, 3)),
        )
        F_mode = kadd(
            kadd(kmul(f_mode, kderivative(D)), kscale(kmul(kderivative(f_mode), D), -1)),
            kscale(kmul(f_mode, D[1:]), -2),
        )
        assert not column_combination(cf_matrix, cf_columns, {"C": C_mode, "F": F_mode})
        analytic_cf.append((C_mode, F_mode))
    analytic_cf_vectors = [pair_vector(cf_columns, {"C": C_mode, "F": F_mode}) for C_mode, F_mode in analytic_cf]
    assert len(rref_with_transform([list(row) for row in zip(*analytic_cf_vectors)])[1]) == 3
    rhs_cf_poly = padd(pmul(B, pderivative(E)), pscale(pmul(E, pderivative(B)), -2))
    rhs_cf = [rhs_cf_poly[degree] if degree < len(rhs_cf_poly) else zero for degree in range(20)]
    cf_solution, cf_compatibility = solve_with_parameters(
        cf_rref, cf_pivots, cf_transform, rhs_cf, parameters[2:5]
    )
    assert not cf_compatibility
    C, F = vector_to_pair(cf_solution, cf_columns, ("C", "F"))

    g_columns = [("G", degree) for degree in range(1, 13)]
    g_matrix = linear_matrix(fixed, g_columns, "G")
    g_rref, g_pivots, g_transform = rref_with_transform(g_matrix)
    assert len(kernel_from_rref(g_rref, g_pivots)) == 0
    rhs_g_poly = padd(
        padd(pmul(B, pderivative(F)), pscale(pmul(E, pderivative(C)), -2)),
        pscale(pmul(F, pderivative(B)), -1),
    )
    rhs_g = [rhs_g_poly[degree] if degree < len(rhs_g_poly) else zero for degree in range(20)]
    g_solution, g_compatibility = solve_with_parameters(g_rref, g_pivots, g_transform, rhs_g, [])
    # ``g_solution[j-1]`` is the coefficient of X^j, 1 <= j <= 12.
    # The omitted coefficient g_0 is the one-dimensional target-translation
    # kernel and does not occur in either remaining Wronskian equation.
    Gprime = [(index + 1) * value for index, value in enumerate(g_solution)]

    j0 = padd(pscale(pmul(B, Gprime), -1), pmul(F, pderivative(C)))
    final_equations = [value for value in g_compatibility + j0 if value]
    leading_open = B[8]
    assert leading_open
    matrix_data = {
        "BE": matrix_summary(be_matrix, be_pivots, be_transform, len(be_columns)),
        "CF": matrix_summary(cf_matrix, cf_pivots, cf_transform, len(cf_columns)),
        "G_derivative": matrix_summary(g_matrix, g_pivots, g_transform, len(g_columns)),
    }
    matrix_data["BE"].update({
        "formula": "(B,E) -> -2*A*E'-B*D'+3*D*B'+2*E*A'",
        "source_basis": "B:X^1..X^8; E:X^2..X^12",
        "target_basis": "X^0..X^19",
        "complete_kernel_modes": ["phi=(A/X)^2-1", "phi=X"],
        "analytic_kernel_hash": hash_kpolys([polynomial for pair in analytic_be for polynomial in pair]),
    })
    matrix_data["CF"].update({
        "formula": "(C,F) -> 3*D*C'-2*A*F'+A'*F",
        "source_basis": "C:X^0..X^8; F:X^1..X^12",
        "target_basis": "X^0..X^19",
        "complete_homogeneous_kernel_modes": ["f=D/X^2", "f=1", "f=X"],
        "analytic_kernel_hash": hash_kpolys([polynomial for pair in analytic_cf for polynomial in pair]),
    })
    matrix_data["G_derivative"].update({
        "formula": "G_without_constant -> -2*A*G'",
        "source_basis": "G:X^1..X^12",
        "target_basis": "X^0..X^19",
    })
    report = {
        "parameter_names": [f"p{i}" for i in range(parameter_count)],
        "linear_maps": matrix_data,
        "CF_forcing_compatibilities": len(cf_compatibility),
        "G_cokernel_equations": len(g_compatibility),
        "G_constant_kernel_dimension": 1,
        "J0_equations": sum(bool(value) for value in j0),
        "final_equation_count": len(final_equations),
        "final_equation_degrees": sorted({value.total_degree() for value in final_equations}),
        "final_equation_term_counts": sorted(len(value.terms) for value in final_equations),
        "final_equations_hash": hash_parameter_polynomials(final_equations),
        "degree_B_open": leading_open.text([f"p{i}" for i in range(parameter_count)]),
    }
    return report, final_equations, leading_open


def matrix_summary(
    matrix: list[list[KElement]],
    pivots: Sequence[int],
    transform: list[list[KElement]],
    columns: int,
) -> dict:
    return {
        "shape": [len(matrix), columns],
        "rank": len(pivots),
        "kernel_dimension": columns - len(pivots),
        "cokernel_dimension_in_declared_target": len(matrix) - len(pivots),
        "matrix_hash": hash_kmatrix(matrix),
        "cokernel_functionals_hash": hash_kmatrix(transform[len(pivots):]),
    }


def column_combination(
    matrix: list[list[KElement]],
    columns: Sequence[tuple[str, int]],
    polynomials: dict[str, KPoly],
) -> KPoly:
    vector = pair_vector(columns, polynomials)
    return ktrim([sum((entry * value for entry, value in zip(row, vector)), K_ZERO) for row in matrix])


def pair_vector(
    columns: Sequence[tuple[str, int]], polynomials: dict[str, KPoly]
) -> list[KElement]:
    return [coefficient(polynomials[variable], degree) for variable, degree in columns]


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def hash_kmatrix(matrix: Sequence[Sequence[KElement]]) -> str:
    payload = [[element.text() for element in row] for row in matrix]
    return hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()


def hash_kpolys(polynomials: Sequence[KPoly]) -> str:
    payload = [[element.text() for element in polynomial] for polynomial in polynomials]
    return hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()


def hash_parameter_polynomials(polynomials: Sequence[ParameterPolynomial]) -> str:
    names = [f"p{i}" for i in range(polynomials[0].variables)] if polynomials else []
    payload = [polynomial.text(names) for polynomial in polynomials]
    return hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()


def singular_final_check(
    equations: Sequence[ParameterPolynomial], leading_open: ParameterPolynomial, timeout: int = 300
) -> dict:
    names = [f"p{i}" for i in range(leading_open.variables)]
    field_polynomial = "+".join(f"({coefficient})*q^{degree}" for degree, coefficient in enumerate(H) if coefficient)
    generators = [equation.text(names) for equation in equations]
    generators.append(f"z*({leading_open.text(names)})-1")
    source = "\n".join([
        'ring R=(0,q),(p0,p1,p2,p3,p4,z),dp;',
        f"minpoly={field_polynomial};",
        "option(redSB);",
        "ideal I=" + ",\n".join(generators) + ";",
        "ideal G=std(I);",
        'print("BASIS_SIZE="+string(size(G)));',
        'if(size(G)==1 && (G[1]==1 || G[1]==-1)){print("UNIT_IDEAL=1");}else{print("UNIT_IDEAL=0");}',
        "quit;",
    ])
    with tempfile.TemporaryDirectory(prefix="jc2-belyi-") as temporary:
        path = Path(temporary) / "final.sing"
        path.write_text(source, encoding="utf-8")
        completed = subprocess.run(
            [SINGULAR, "-q", str(path)],
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    if completed.returncode:
        raise RuntimeError(f"Singular failed: {completed.stderr or completed.stdout}")
    unit = "UNIT_IDEAL=1" in completed.stdout
    size_line = next(line for line in completed.stdout.splitlines() if line.startswith("BASIS_SIZE="))
    return {
        "saturation": "<final equations, z*B_8-1>",
        "unit_ideal": unit,
        "basis_size": int(size_line.split("=", 1)[1]),
        "singular_input_hash": hashlib.sha256(source.encode()).hexdigest(),
    }


def arithmetic_audit() -> dict:
    import sympy as sp

    q = sp.symbols("q")
    polynomial = sum(int(coefficient) * q**degree for degree, coefficient in enumerate(H))
    exact_polynomial = sp.Poly(polynomial, q, domain=sp.QQ)
    group, _ = sp.polys.numberfields.galois_group(polynomial, q)
    roots = sp.nroots(polynomial, maxsteps=200)
    real_roots = sum(multiplicity for _, multiplicity in exact_polynomial.intervals())
    return {
        "irreducible_over_Q": bool(exact_polynomial.is_irreducible),
        "galois_group_order": int(group.order()),
        "galois_group": "S5" if group.order() == 120 else str(group),
        "field_of_moduli_degree": 5,
        "signature": [real_roots, (exact_polynomial.degree() - real_roots) // 2],
        "embedding_q_values_approx": [str(root) for root in roots],
        "interpretation": "the five normalized maps form one Galois orbit; X=0 and infinity are sections over the quintic field, so the distinguished-point/polynomial-normalization gate does not eliminate an orbit member",
    }


def compile_report(run_singular: bool = True, singular_timeout: int = 300) -> dict:
    A, D = reconstruct_top()
    dessins = enumerate_dessins()
    deformation_report, equations, leading_open = compile_deformations(A, D)
    report = {
        "scope": "no-vertical-edge (72,108) / (8,28) Laurent residue only",
        "status": "exact finite Belyi reconstruction and linear-deformation compilation",
        "claim_boundary": "This independently compiles the finite branch; it is not a stand-alone proof of the published degree reduction or JC(2).",
        "coefficient_field": {
            "presentation": "QQ[q]/(H)",
            "H_coefficients_ascending": H,
            **arithmetic_audit(),
        },
        "top_reconstruction": {
            "normalization": "A=X+x2*X^2+...+x6*X^6+q*X^7+q*X^8; ord(D)=2; 2*A*D'-3*A'*D=X^2",
            "A_D_hash": hash_kpolys([A, D]),
            "A_coefficients_ascending": [value.text() for value in A],
            "D_coefficients_ascending": [value.text() for value in D],
            "quotient_graph_source": "plane-jc/cas/firstblock_mu7_quotient_lex_basis.txt",
            "quotient_graph_sha256": hashlib.sha256(GRAPH_FILE.read_bytes()).hexdigest(),
            **belyi_audit(A, D),
        },
        "dessins": {
            "labelled_center_sets": 85,
            "rotation_orbits": len(dessins),
            "representatives": dessins,
            "galois_orbit_matching": "five exact embeddings and five combinatorial types; the S5 action is transitive setwise",
        },
        "deformations": deformation_report,
    }
    if run_singular:
        report["terminal_open_check"] = singular_final_check(equations, leading_open, singular_timeout)
    else:
        report["terminal_open_check"] = {"status": "not run"}
    return report


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-singular", action="store_true", help="compile exact maps and linear stages without the terminal small Groebner check")
    parser.add_argument("--singular-timeout", type=int, default=300)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="output path; the full run defaults to the pinned artifact, while --no-singular writes only when this is explicit",
    )
    arguments = parser.parse_args()
    report = compile_report(not arguments.no_singular, arguments.singular_timeout)
    output = arguments.output if arguments.output is not None else (None if arguments.no_singular else OUT)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
