"""Bounded discovery of one-parameter elliptic families from target curves.

Unlike :mod:`ecsearch.conductor_engineering`, this module does not assume a
single candidate family.  It screens a declared collection of family
generators by modular roots of the target j-equation, factors only surviving
equations over QQ, specializes every rational parameter, and distinguishes an
exact j-match from a Q-isomorphism.

The engine is construction-agnostic.  This file supplies adapters for
polynomial Weierstrass models and symmetric six-root Mestre constructions;
additional construction spaces can implement the same small protocol.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from math import isqrt
from typing import Any, Callable, Protocol, Sequence

from .conductor_engineering import (
    polynomial_add,
    polynomial_multiply,
    polynomial_scale,
    small_prime_valuations,
    weierstrass_invariant_polynomials,
)


Q = Fraction
RationalModel = tuple[Fraction, Fraction, Fraction, Fraction, Fraction]


def rational_text(value: Fraction | int) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def evaluate_polynomial(
    coefficients: Sequence[int | Fraction], value: Fraction
) -> Fraction:
    result = Q(0)
    for coefficient in reversed(coefficients):
        result = result * value + Q(coefficient)
    return result


def rational_weierstrass_c_data(
    coefficients: Sequence[int | Fraction],
) -> dict[str, Fraction]:
    """Return exact ``c4``, ``c6`` and discriminant, including singular models."""

    if len(coefficients) != 5:
        raise ValueError("five Weierstrass coefficients are required")
    a1, a2, a3, a4, a6 = map(Q, coefficients)
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = (
        a1 * a1 * a6
        + 4 * a2 * a6
        - a1 * a3 * a4
        + a2 * a3 * a3
        - a4 * a4
    )
    c4 = b2 * b2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    discriminant = (
        -b2 * b2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    )
    return {
        "c4": c4,
        "c6": c6,
        "discriminant": discriminant,
    }


def rational_weierstrass_invariants(
    coefficients: Sequence[int | Fraction],
) -> dict[str, Fraction]:
    """Return exact invariants and j for a nonsingular rational model."""

    result = rational_weierstrass_c_data(coefficients)
    if result["discriminant"] == 0:
        raise ValueError("the Weierstrass model is singular")
    result["j"] = result["c4"] ** 3 / result["discriminant"]
    return result


def _integer_nth_root(value: int, exponent: int) -> int | None:
    if value < 0:
        if exponent % 2 == 0:
            return None
        root = _integer_nth_root(-value, exponent)
        return -root if root is not None else None
    if value in (0, 1):
        return value
    low, high = 0, 1
    while high**exponent < value:
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if middle**exponent < value:
            low = middle
        else:
            high = middle
    return high if high**exponent == value else None


def rational_nth_root(value: Fraction, exponent: int) -> Fraction | None:
    if exponent < 1:
        raise ValueError("root exponent must be positive")
    value = Q(value)
    numerator = _integer_nth_root(value.numerator, exponent)
    denominator = _integer_nth_root(value.denominator, exponent)
    if numerator is None or denominator is None:
        return None
    return Q(numerator, denominator)


def q_isomorphism_scale(
    source: dict[str, Fraction], target: dict[str, Fraction]
) -> Fraction | None:
    """Return invariant scale ``u`` when two models are Q-isomorphic.

    The convention is ``target.c4=u^4*source.c4`` and similarly in weights
    six and twelve.  For the ordinary ``j != 0,1728`` case both ratios are
    checked.  The two special j-values use the nonzero invariant.
    """

    if source["j"] != target["j"]:
        return None
    if source["c4"] == 0:
        scale = rational_nth_root(target["c6"] / source["c6"], 6)
    else:
        scale = rational_nth_root(target["c4"] / source["c4"], 4)
    if scale is None:
        return None
    if source["c4"] and target["c4"] != scale**4 * source["c4"]:
        return None
    if source["c6"] and target["c6"] != scale**6 * source["c6"]:
        return None
    if target["discriminant"] != scale**12 * source["discriminant"]:
        return None
    return scale


def fraction_valuation(value: Fraction, prime: int) -> int:
    value = Q(value)

    def integer_valuation(integer: int) -> int:
        integer = abs(integer)
        exponent = 0
        while integer and integer % prime == 0:
            integer //= prime
            exponent += 1
        return exponent

    return integer_valuation(value.numerator) - integer_valuation(value.denominator)


def _evaluate_mod(coefficients: Sequence[int], value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def _multiply_mod(
    left: Sequence[int], right: Sequence[int], prime: int
) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] = (
                result[left_index + right_index]
                + left_value * right_value
            ) % prime
    return result


def _add_mod(
    left: Sequence[int], right: Sequence[int], prime: int
) -> list[int]:
    result = [0] * max(len(left), len(right))
    for index in range(len(result)):
        result[index] = (
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
        ) % prime
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def _interpolate_mod(
    points: Sequence[tuple[int, int]], prime: int
) -> tuple[int, ...]:
    result = [0]
    for index, (x_value, y_value) in enumerate(points):
        basis = [1]
        denominator = 1
        for other_index, (other_x, _) in enumerate(points):
            if index == other_index:
                continue
            basis = _multiply_mod(basis, (-other_x, 1), prime)
            denominator = denominator * (x_value - other_x) % prime
        scale = y_value * pow(denominator, -1, prime) % prime
        result = _add_mod(
            result, [scale * coefficient % prime for coefficient in basis], prime
        )
    return tuple(result)


def _multiply_fraction(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> list[Fraction]:
    result = [Q(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] += left_value * right_value
    return result


def _add_fraction(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> list[Fraction]:
    result = [Q(0)] * max(len(left), len(right))
    for index in range(len(result)):
        result[index] = (
            left[index] if index < len(left) else Q(0)
        ) + (right[index] if index < len(right) else Q(0))
    return result


def interpolate_fraction_polynomial(
    points: Sequence[tuple[Fraction, Fraction]],
) -> tuple[Fraction, ...]:
    result = [Q(0)]
    for index, (x_value, y_value) in enumerate(points):
        basis = [Q(1)]
        denominator = Q(1)
        for other_index, (other_x, _) in enumerate(points):
            if index == other_index:
                continue
            basis = _multiply_fraction(basis, (-other_x, Q(1)))
            denominator *= x_value - other_x
        result = _add_fraction(
            result,
            [y_value * coefficient / denominator for coefficient in basis],
        )
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def rational_roots(coefficients: Sequence[Fraction]) -> tuple[Fraction, ...]:
    """Factor an exact polynomial over QQ and return its rational roots."""

    try:
        import sympy as sp
    except ImportError as error:  # pragma: no cover - repository dependency
        raise RuntimeError("sympy is required for exact family discovery") from error
    variable = sp.symbols("z")
    polynomial = sp.Poly(
        sum(
            sp.Rational(value.numerator, value.denominator) * variable**index
            for index, value in enumerate(coefficients)
        ),
        variable,
        domain=sp.QQ,
    )
    return tuple(
        sorted(
            Q(int(root.p), int(root.q))
            for root in sp.polys.polytools.ground_roots(polynomial)
        )
    )


@dataclass(frozen=True)
class DiscoveryTarget:
    label: str
    coefficients: RationalModel
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def invariants(self) -> dict[str, Fraction]:
        return rational_weierstrass_invariants(self.coefficients)


class DiscoveryFamily(Protocol):
    identifier: str
    kind: str
    parameter_name: str
    metadata: dict[str, Any]

    def modular_roots(
        self, target_j: Fraction, prime: int
    ) -> tuple[int, ...] | None: ...

    def exact_parameters(self, target_j: Fraction) -> tuple[Fraction, ...]: ...

    def specialize(self, parameter: Fraction) -> RationalModel: ...


@dataclass(frozen=True)
class PolynomialWeierstrassFamily:
    """A family whose five a-invariants are integral polynomials."""

    identifier: str
    coefficient_polynomials: dict[str, tuple[int, ...]]
    metadata: dict[str, Any] = field(default_factory=dict)
    kind: str = "polynomial-weierstrass"
    parameter_name: str = "t"

    @property
    def invariant_polynomials(self) -> dict[str, tuple[int, ...]]:
        return weierstrass_invariant_polynomials(self.coefficient_polynomials)

    def _j_equation(self, target_j: Fraction) -> tuple[int, ...]:
        invariants = self.invariant_polynomials
        numerator = polynomial_scale(
            polynomial_multiply(
                invariants["c4"], invariants["c4"], invariants["c4"]
            ),
            target_j.denominator,
        )
        denominator = polynomial_scale(
            invariants["discriminant"], target_j.numerator
        )
        equation = polynomial_add(numerator, polynomial_scale(denominator, -1))
        return tuple(equation)

    def modular_roots(
        self, target_j: Fraction, prime: int
    ) -> tuple[int, ...] | None:
        if target_j.denominator % prime == 0:
            return None
        equation = self._j_equation(target_j)
        # A rational parameter whose denominator is divisible by ``prime``
        # reduces to the point at infinity.  Only use an affine no-root result
        # as an exclusion when the homogenized equation is nonzero there.
        if equation[-1] % prime == 0:
            return None
        return tuple(
            residue
            for residue in range(prime)
            if _evaluate_mod(equation, residue, prime) == 0
        )

    def exact_parameters(self, target_j: Fraction) -> tuple[Fraction, ...]:
        roots = rational_roots(tuple(map(Q, self._j_equation(target_j))))
        invariants = self.invariant_polynomials
        sign_symmetric = all(
            coefficient == 0
            for polynomial in invariants.values()
            for index, coefficient in enumerate(polynomial)
            if index % 2
        )
        if sign_symmetric:
            roots = tuple(sorted({abs(root) for root in roots}))
        return roots

    def specialize(self, parameter: Fraction) -> RationalModel:
        return tuple(
            evaluate_polynomial(self.coefficient_polynomials[name], parameter)
            for name in ("a1", "a2", "a3", "a4", "a6")
        )  # type: ignore[return-value]


def mestre_family_equation_mod(
    roots: Sequence[int], parameter: int, target_j: Fraction, prime: int
) -> int:
    """Evaluate the six-root Mestre target-j equation modulo an odd prime."""

    parameter %= prime
    if parameter == 0 or target_j.denominator % prime == 0:
        raise ValueError("bad modular evaluation point")
    product_coefficients = [1]
    for root in roots:
        product_coefficients = _multiply_mod(
            product_coefficients, (-(root + parameter), 1), prime
        )
    for root in roots:
        product_coefficients = _multiply_mod(
            product_coefficients, (-(root - parameter), 1), prime
        )
    approximant = [0] * 7
    approximant[6] = 1
    inverse_two = pow(2, -1, prime)
    for index in range(5, -1, -1):
        square = _multiply_mod(approximant, approximant, prime)
        degree = 6 + index
        approximant[index] = (
            product_coefficients[degree] - square[degree]
        ) * inverse_two % prime
    square = _multiply_mod(approximant, approximant, prime)
    inverse_t2 = pow(parameter * parameter % prime, -1, prime)
    quartic = [
        (square[index] - product_coefficients[index]) * inverse_t2 % prime
        for index in range(5)
    ]
    e, d, c, b, a = quartic
    invariant_i = (12 * a * e - 3 * b * d + c * c) % prime
    invariant_j = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d * d
        - 27 * b * b * e
        - 2 * c**3
    ) % prime
    a4 = -27 * invariant_i % prime
    a6 = -27 * invariant_j % prime
    c4 = -48 * a4 % prime
    discriminant = -16 * (4 * a4**3 + 27 * a6**2) % prime
    reduced_j = target_j.numerator * pow(target_j.denominator, -1, prime) % prime
    return (c4**3 - reduced_j * discriminant) % prime


@dataclass(frozen=True)
class SixRootMestreFamily:
    """A generated symmetric-shift family from one fixed six-root tuple."""

    roots: tuple[int, ...]
    exact_model: Callable[[Fraction], Sequence[Fraction]] = field(
        compare=False, repr=False
    )
    metadata: dict[str, Any] = field(default_factory=dict, compare=False)
    kind: str = "six-root-mestre"
    parameter_name: str = "T"

    @property
    def identifier(self) -> str:
        return "six-root-mestre:" + ",".join(map(str, self.roots))

    def modular_roots(
        self, target_j: Fraction, prime: int
    ) -> tuple[int, ...] | None:
        if target_j.denominator % prime == 0:
            return None
        samples = []
        occupied: set[int] = set()
        for parameter in range(1, 14):
            z_value = parameter * parameter % prime
            if z_value in occupied:
                return None
            occupied.add(z_value)
            samples.append(
                (
                    z_value,
                    mestre_family_equation_mod(
                        self.roots, parameter, target_j, prime
                    ),
                )
            )
        polynomial = _interpolate_mod(samples, prime)
        if len(polynomial) != 13 or polynomial[-1] == 0:
            return None
        for parameter in (14, 15, prime - 14):
            expected = mestre_family_equation_mod(
                self.roots, parameter, target_j, prime
            )
            if (
                _evaluate_mod(polynomial, parameter * parameter % prime, prime)
                != expected
            ):
                raise AssertionError("the Mestre j-equation exceeded degree 12")
        return tuple(
            value
            for value in range(prime)
            if _evaluate_mod(polynomial, value, prime) == 0
        )

    def _exact_equation(self, target_j: Fraction) -> tuple[Fraction, ...]:
        points = []
        for parameter in range(1, 14):
            model = tuple(map(Q, self.exact_model(Q(parameter))))
            invariants = rational_weierstrass_c_data(model)
            points.append(
                (
                    Q(parameter * parameter),
                    invariants["c4"] ** 3
                    - target_j * invariants["discriminant"],
                )
            )
        equation = interpolate_fraction_polynomial(points)
        if len(equation) != 13:
            raise AssertionError("the exact Mestre j-equation lost degree 12")
        for parameter in (Q(14), Q(15), Q(-14)):
            invariants = rational_weierstrass_c_data(
                self.specialize(parameter)
            )
            expected = invariants["c4"] ** 3 - target_j * invariants[
                "discriminant"
            ]
            observed = evaluate_polynomial(equation, parameter * parameter)
            if observed != expected:
                raise AssertionError("exact Mestre j interpolation failed")
        return equation

    def exact_parameters(self, target_j: Fraction) -> tuple[Fraction, ...]:
        parameters = []
        for z_value in rational_roots(self._exact_equation(target_j)):
            if z_value <= 0:
                continue
            numerator = isqrt(z_value.numerator)
            denominator = isqrt(z_value.denominator)
            if (
                numerator * numerator == z_value.numerator
                and denominator * denominator == z_value.denominator
            ):
                parameters.append(Q(numerator, denominator))
        return tuple(sorted(set(parameters)))

    def specialize(self, parameter: Fraction) -> RationalModel:
        return tuple(map(Q, self.exact_model(Q(parameter))))  # type: ignore[return-value]


def discover_target_families(
    target: DiscoveryTarget,
    families: Sequence[DiscoveryFamily],
    *,
    modular_primes: Sequence[int],
    fingerprint_prime_bound: int = 97,
) -> dict[str, Any]:
    """Run the modular, exact-factorization, and Q-isomorphism gates."""

    target_invariants = target.invariants
    target_j = target_invariants["j"]
    witness_histogram: dict[str, int] = {}
    exact_survivors = 0
    rational_parameter_survivors = []
    exact_matches = []
    for family in families:
        excluded = False
        local_profile = []
        for prime in modular_primes:
            roots = family.modular_roots(target_j, prime)
            if roots is None:
                continue
            local_profile.append([prime, len(roots)])
            if not roots:
                witness_histogram[str(prime)] = (
                    witness_histogram.get(str(prime), 0) + 1
                )
                excluded = True
                break
        if excluded:
            continue
        exact_survivors += 1
        parameters = family.exact_parameters(target_j)
        if parameters:
            rational_parameter_survivors.append(
                {
                    "family_id": family.identifier,
                    "family_kind": family.kind,
                    "modular_root_counts": local_profile,
                    "parameters": [rational_text(value) for value in parameters],
                    "metadata": family.metadata,
                }
            )
        for parameter in parameters:
            try:
                source_model = family.specialize(parameter)
                source_invariants = rational_weierstrass_invariants(source_model)
            except (ArithmeticError, ValueError, ZeroDivisionError):
                continue
            if source_invariants["j"] != target_j:
                raise AssertionError("an exact equation root failed j verification")
            scale = q_isomorphism_scale(source_invariants, target_invariants)
            repeated_primes = [
                prime
                for prime, exponent in small_prime_valuations(
                    target_invariants["discriminant"].numerator,
                    fingerprint_prime_bound,
                )
                if exponent >= 2
            ]
            exact_matches.append(
                {
                    "family_id": family.identifier,
                    "family_kind": family.kind,
                    "parameter_name": family.parameter_name,
                    "parameter": rational_text(parameter),
                    "metadata": family.metadata,
                    "source_model": [rational_text(value) for value in source_model],
                    "source_discriminant": rational_text(
                        source_invariants["discriminant"]
                    ),
                    "q_isomorphic": scale is not None,
                    "q_isomorphism_invariant_scale": (
                        rational_text(scale) if scale is not None else None
                    ),
                    "target_repeated_prime_valuations": {
                        str(prime): fraction_valuation(
                            target_invariants["discriminant"], prime
                        )
                        for prime in repeated_primes
                    },
                    "source_discriminant_valuations_at_target_repeated_primes": {
                        str(prime): fraction_valuation(
                            source_invariants["discriminant"], prime
                        )
                        for prime in repeated_primes
                    },
                }
            )
    return {
        "target": target.label,
        "target_j": rational_text(target_j),
        "families_tested": len(families),
        "modular_primes": list(modular_primes),
        "modular_no_root_witness_histogram": witness_histogram,
        "exact_factorization_survivor_count": exact_survivors,
        "rational_parameter_survivors": rational_parameter_survivors,
        "exact_j_matches": exact_matches,
        "q_isomorphism_matches": [
            match for match in exact_matches if match["q_isomorphic"]
        ],
    }
