#!/usr/bin/env python3
"""Dependency-free verifier for a proof-carrying arithmetic compilation.

Only Python's standard library is used.  In particular, this script imports
neither SymPy nor any module from ``jcsearch``.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import reduce
import hashlib
import json
from math import comb, gcd, isqrt, prod
from pathlib import Path


EXPECTED_SCHEMA = "proof-carrying-arithmetic-compilation/v1"
EXPECTED_MAP_SERIALIZATION = "expanded-sparse-rational-json/v1"
DEFAULT_CERTIFICATE = (
    Path(__file__).resolve().parents[1]
    / "artifacts/generated-results/arithmetic_keller_quintic.json"
)


def q(value: str | int) -> Fraction:
    return Fraction(value)


def qtext(value: Fraction) -> str:
    value = Fraction(value)
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def trim(polynomial):
    polynomial = list(polynomial)
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def padd(left, right):
    return trim(
        [
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
            for index in range(max(len(left), len(right)))
        ]
    )


def pscale(scalar, polynomial):
    return trim([scalar * coefficient for coefficient in polynomial])


def psub(left, right):
    return padd(left, pscale(-1, right))


def pmul(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return trim(result)


def pderivative(polynomial):
    return trim(
        [degree * polynomial[degree] for degree in range(1, len(polynomial))]
        or [0]
    )


def peval(polynomial, value):
    result = 0
    for coefficient in reversed(polynomial):
        result = result * value + coefficient
    return result


def pdivmod(dividend, divisor):
    dividend = trim(dividend)
    divisor = trim(divisor)
    assert divisor != [0]
    quotient = [0] * max(1, len(dividend) - len(divisor) + 1)
    while dividend != [0] and len(dividend) >= len(divisor):
        shift = len(dividend) - len(divisor)
        scalar = dividend[-1] / divisor[-1]
        quotient[shift] += scalar
        dividend = psub(dividend, [0] * shift + pscale(scalar, divisor))
    return trim(quotient), trim(dividend)


def pgcd(left, right):
    while trim(right) != [0]:
        _, remainder = pdivmod(left, right)
        left, right = right, remainder
    left = trim(left)
    return pscale(1 / left[-1], left) if left != [0] else [0]


def translate(polynomial, amount):
    result = [Fraction(0)] * len(polynomial)
    for degree, coefficient in enumerate(polynomial):
        for exponent in range(degree + 1):
            result[exponent] += (
                coefficient * comb(degree, exponent) * amount ** (degree - exponent)
            )
    return trim(result)


def determinant(matrix):
    matrix = [[Fraction(value) for value in row] for row in matrix]
    result = Fraction(1)
    for column in range(len(matrix)):
        pivot = next(
            (row for row in range(column, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            result *= -1
        pivot_value = matrix[column][column]
        result *= pivot_value
        for row in range(column + 1, len(matrix)):
            scalar = matrix[row][column] / pivot_value
            for index in range(column + 1, len(matrix)):
                matrix[row][index] -= scalar * matrix[column][index]
    return result


def resultant(left, right):
    left = trim(left)
    right = trim(right)
    m, n = len(left) - 1, len(right) - 1
    left_descending = list(reversed(left))
    right_descending = list(reversed(right))
    matrix = []
    for shift in range(n):
        matrix.append(
            [0] * shift + left_descending + [0] * (n - 1 - shift)
        )
    for shift in range(m):
        matrix.append(
            [0] * shift + right_descending + [0] * (m - 1 - shift)
        )
    return determinant(matrix)


def discriminant(polynomial):
    degree = len(trim(polynomial)) - 1
    return (
        (-1) ** (degree * (degree - 1) // 2)
        * resultant(polynomial, pderivative(polynomial))
        / polynomial[-1]
    )


def valuation(value, prime):
    value = Fraction(value)
    if not value:
        return None
    numerator = abs(value.numerator)
    denominator = value.denominator
    result = 0
    while numerator % prime == 0:
        result += 1
        numerator //= prime
    while denominator % prime == 0:
        result -= 1
        denominator //= prime
    return result


def is_prime(value):
    return value >= 2 and all(
        value % divisor for divisor in range(2, isqrt(value) + 1)
    )


def mod_trim(polynomial, prime):
    return trim([int(value) % prime for value in polynomial])


def mod_divmod(dividend, divisor, prime):
    dividend = mod_trim(dividend, prime)
    divisor = mod_trim(divisor, prime)
    quotient = [0] * max(1, len(dividend) - len(divisor) + 1)
    inverse = pow(divisor[-1], -1, prime)
    while dividend != [0] and len(dividend) >= len(divisor):
        shift = len(dividend) - len(divisor)
        scalar = dividend[-1] * inverse % prime
        quotient[shift] = (quotient[shift] + scalar) % prime
        for index, coefficient in enumerate(divisor):
            dividend[index + shift] = (
                dividend[index + shift] - scalar * coefficient
            ) % prime
        dividend = mod_trim(dividend, prime)
    return mod_trim(quotient, prime), dividend


def mod_gcd(left, right, prime):
    while mod_trim(right, prime) != [0]:
        _, remainder = mod_divmod(left, right, prime)
        left, right = right, remainder
    left = mod_trim(left, prime)
    inverse = pow(left[-1], -1, prime)
    return mod_trim([inverse * coefficient for coefficient in left], prime)


def mod_mul(left, right, modulus, prime):
    product_polynomial = pmul(left, right)
    _, remainder = mod_divmod(product_polynomial, modulus, prime)
    return remainder


def mod_pow(base, exponent, modulus, prime):
    result = [1]
    while exponent:
        if exponent & 1:
            result = mod_mul(result, base, modulus, prime)
        base = mod_mul(base, base, modulus, prime)
        exponent //= 2
    return result


def prime_divisors(value):
    result = set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            result.add(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        result.add(value)
    return result


def mod_irreducible(polynomial, prime):
    polynomial = mod_trim(polynomial, prime)
    degree = len(polynomial) - 1
    assert degree >= 1 and polynomial[-1] % prime
    inverse = pow(polynomial[-1], -1, prime)
    polynomial = mod_trim(
        [inverse * coefficient for coefficient in polynomial], prime
    )
    x_polynomial = [0, 1]
    for divisor in prime_divisors(degree):
        power = mod_pow(x_polynomial, prime ** (degree // divisor), polynomial, prime)
        if len(mod_gcd(psub(power, x_polynomial), polynomial, prime)) > 1:
            return False
    _, final_remainder = mod_divmod(
        psub(
            mod_pow(x_polynomial, prime**degree, polynomial, prime),
            x_polynomial,
        ),
        polynomial,
        prime,
    )
    return final_remainder == [0]


def reduce_rational(polynomial, prime):
    return mod_trim(
        [
            coefficient.numerator
            * pow(coefficient.denominator, -1, prime)
            for coefficient in polynomial
        ],
        prime,
    )


def sturm_sequence(polynomial):
    sequence = [trim(polynomial), pderivative(polynomial)]
    while sequence[-1] != [0]:
        _, remainder = pdivmod(sequence[-2], sequence[-1])
        if remainder == [0]:
            break
        sequence.append(pscale(-1, remainder))
    return sequence


def variations(signs):
    signs = [sign for sign in signs if sign]
    return sum(left != right for left, right in zip(signs, signs[1:]))


def sign(value):
    return (value > 0) - (value < 0)


def sturm_variations_at(sequence, point):
    return variations([sign(peval(polynomial, point)) for polynomial in sequence])


def sturm_variations_at_infinity(sequence, positive):
    signs = []
    for polynomial in sequence:
        leading_sign = sign(polynomial[-1])
        if not positive and (len(polynomial) - 1) % 2:
            leading_sign *= -1
        signs.append(leading_sign)
    return variations(signs)


# Sparse multivariate rational polynomials, keyed by (x,y,z)-exponents.
def sadd(left, right):
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = result.get(monomial, 0) + coefficient
        if not result[monomial]:
            del result[monomial]
    return result


def sscale(scalar, polynomial):
    return (
        {
            monomial: scalar * coefficient
            for monomial, coefficient in polynomial.items()
            if scalar * coefficient
        }
        if scalar
        else {}
    )


def smul(left, right):
    result = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(
                a + b for a, b in zip(left_monomial, right_monomial)
            )
            result[monomial] = (
                result.get(monomial, 0)
                + left_coefficient * right_coefficient
            )
            if not result[monomial]:
                del result[monomial]
    return result


def spow(polynomial, exponent):
    result = {(0, 0, 0): Fraction(1)}
    while exponent:
        if exponent & 1:
            result = smul(result, polynomial)
        polynomial = smul(polynomial, polynomial)
        exponent //= 2
    return result


def sderivative(polynomial, variable):
    result = {}
    for monomial, coefficient in polynomial.items():
        if monomial[variable]:
            derived = list(monomial)
            factor = derived[variable]
            derived[variable] -= 1
            result[tuple(derived)] = coefficient * factor
    return result


def canonical_map(mapping):
    return [
        [
            [list(monomial), qtext(coefficient)]
            for monomial, coefficient in sorted(coordinate.items())
        ]
        for coordinate in mapping
    ]


def reconstruct_map(seed, stable_parameter=None):
    one = {(0, 0, 0): Fraction(1)}
    x = {(1, 0, 0): Fraction(1)}
    y = {(0, 1, 0): Fraction(1)}
    z = {(0, 0, 1): Fraction(1)}
    t = sadd(one, smul(x, y))
    g1, g2, g3 = seed[1], seed[2], seed[3]
    qpoly = sadd(
        smul(spow(t, 2), z),
        sscale(
            g1 / g3,
            smul(spow(y, 2), sadd(one, sscale(3, t))),
        ),
    )
    first = smul(t, qpoly)
    second = sadd(
        y,
        sadd(
            sscale(3 * g3 / g1, smul(x, qpoly)),
            sscale(2 * g2 / g1, smul(t, qpoly)),
        ),
    )
    third = sadd(
        smul(x, sadd(sscale(5, one), sscale(-3, t))),
        sscale(-g3 / g1, smul(spow(x, 3), z)),
    )
    if stable_parameter is not None and len(seed) == 4:
        exponent = stable_parameter + 4
        second = sadd(
            second,
            sscale(
                3 * g3 / g1,
                sadd(
                    smul(
                        smul(spow(t, exponent - 1), x),
                        spow(qpoly, exponent),
                    ),
                    sscale(-1, smul(smul(spow(t, 2), x), spow(qpoly, 3))),
                ),
            ),
        )
        third = sadd(
            third,
            sscale(
                -g3 / g1,
                sadd(
                    smul(
                        smul(spow(t, exponent - 3), spow(x, 3)),
                        spow(qpoly, exponent),
                    ),
                    sscale(-1, smul(spow(x, 3), spow(qpoly, 3))),
                ),
            ),
        )
    else:
        shift = stable_parameter or 0
        for degree in range(4, len(seed)):
            second = sadd(
                second,
                sscale(
                    degree * seed[degree] / g1,
                    smul(
                        smul(spow(t, shift + 2), spow(x, degree - 2)),
                        spow(qpoly, degree + shift),
                    ),
                ),
            )
            third = sadd(
                third,
                sscale(
                    -(degree - 2) * seed[degree] / g1,
                    smul(
                        smul(spow(t, shift), spow(x, degree)),
                        spow(qpoly, degree + shift),
                    ),
                ),
            )
    return [first, sscale(Fraction(-1, 2), second), third]


def verify_stable_multiplicity_record(record, degree, seed):
    if record is None:
        return None
    parameter = record["family_parameter"]
    assert isinstance(parameter, int) and not isinstance(parameter, bool)
    assert parameter >= 0
    if degree == 3:
        exponent = parameter + 4
        assert record == {
            "family_parameter": parameter,
            "gauge_exponent": exponent,
            "separation_invariant": "geometric_boundary_target_components",
            "separation_value": exponent,
            "fitting_support": [],
            "boundary_prime_count": None,
            "boundary_ramification_index": None,
        }
    else:
        assert all(seed[index] for index in range(4, degree + 1))
        common_divisor = gcd(degree - 3, parameter + 2)
        assert record == {
            "family_parameter": parameter,
            "gauge_exponent": parameter,
            "separation_invariant": "normalized_fitting_newton_area",
            "separation_value": 2 * degree - 3 + (degree - 2) * parameter,
            "fitting_support": [[0, 0], [1, 2]]
            + [
                [index + parameter, index - 1]
                for index in range(4, degree + 1)
            ],
            "boundary_prime_count": common_divisor,
            "boundary_ramification_index": (
                degree - 3
            )
            // common_divisor,
        }
    return parameter


def verify_local_actions(certificate, global_polynomial, local_models):
    for action in certificate["local_action_certificates"]:
        prime = int(action["prime"])
        factors = [
            [q(value) for value in factor["coefficients_ascending"]]
            for factor in action["factors"]
        ]
        if action["kind"] == "local_factorization":
            assert reduce(pmul, factors, [Fraction(1)]) == local_models[prime]
            for left_index, left in enumerate(factors):
                for right in factors[left_index + 1 :]:
                    assert valuation(resultant(left, right), prime) == 0
            for factor, record in zip(factors, action["factors"]):
                if record["witness"] == "eisenstein":
                    assert factor[-1] == 1
                    assert all(
                        coefficient.denominator == 1
                        and coefficient.numerator % prime == 0
                        for coefficient in factor[:-1]
                    )
                    assert factor[0].numerator % (prime**2) != 0
                elif record["witness"] == "unramified_irreducible_reduction":
                    assert mod_irreducible(reduce_rational(factor, prime), prime)
                    assert valuation(discriminant(factor), prime) == 0
                else:
                    raise AssertionError("unknown local factor witness")
                assert (
                    int(record["ramification_index"])
                    * int(record["residue_degree"])
                    == len(factor) - 1
                )
        elif action["kind"] == "frobenius_factorization":
            reduced_product = reduce(
                lambda left, right: mod_trim(pmul(left, right), prime),
                factors,
                [1],
            )
            assert reduced_product == reduce_rational(global_polynomial, prime)
            assert all(mod_irreducible(factor, prime) for factor in factors)
            assert [len(factor) - 1 for factor in factors] == action["cycle_degrees"]
            for left_index, left in enumerate(factors):
                for right in factors[left_index + 1 :]:
                    assert mod_gcd(left, right, prime) == [1]
        else:
            raise AssertionError("unknown local action certificate kind")


def verify(path: Path) -> None:
    certificate = json.loads(path.read_text(encoding="utf-8"))
    assert certificate["schema"] == EXPECTED_SCHEMA
    degree = int(certificate["degree"])
    global_polynomial = [
        q(value)
        for value in certificate["global_polynomial"][
            "monic_coefficients_ascending"
        ]
    ]
    assert len(global_polynomial) == degree + 1
    assert global_polynomial[-1] == 1

    crt = certificate["coefficient_crt"]
    base_denominator = int(crt["base_denominator"])
    crt_modulus = int(crt["crt_modulus"])
    multiplier = int(crt["multiplier"])
    common_denominator = int(crt["common_denominator"])
    residues = [int(value) for value in crt["coefficient_residues_ascending"]]
    assert common_denominator == base_denominator * (1 + multiplier * crt_modulus)
    assert len(residues) == degree
    assert all(
        (common_denominator * coefficient).denominator == 1
        and int(common_denominator * coefficient) % crt_modulus == residue
        for coefficient, residue in zip(global_polynomial[:-1], residues)
    )
    for coefficient, interval in zip(
        global_polynomial[:-1], crt["real_coefficient_intervals"]
    ):
        assert q(interval[0]) < coefficient < q(interval[1])

    local_models = {}
    local_moduli = {
        int(record["prime"]): int(record["modulus"])
        for record in crt["local_moduli"]
    }
    for record in certificate["local_models_and_precision_claims"]:
        prime = int(record["prime"])
        assert is_prime(prime)
        model = [q(value) for value in record["model_coefficients_ascending"]]
        local_models[prime] = model
        claimed_discriminant = q(record["discriminant"])
        assert discriminant(model) == claimed_discriminant
        discriminant_valuation = valuation(claimed_discriminant, prime)
        assert discriminant_valuation == int(record["discriminant_valuation"])
        universal_precision = 2 * discriminant_valuation + 1
        claimed_precision = int(record["claimed_precision"])
        assert universal_precision == int(record["universal_precision"])
        assert claimed_precision >= universal_precision
        denominator_valuation = valuation(Fraction(base_denominator), prime)
        assert denominator_valuation is not None
        assert local_moduli[prime] == prime ** (
            claimed_precision + denominator_valuation
        )
        for coefficient, center, residue in zip(
            global_polynomial[:-1], model[:-1], residues
        ):
            coefficient_valuation = valuation(coefficient - center, prime)
            assert (
                coefficient_valuation is None
                or coefficient_valuation >= claimed_precision
            )
            scaled_center = base_denominator * center
            assert scaled_center.denominator == 1
            assert int(scaled_center) % local_moduli[prime] == (
                residue % local_moduli[prime]
            )
    assert prod(local_moduli.values()) == crt_modulus

    primitive = certificate["global_polynomial"][
        "primitive_integer_coefficients_ascending"
    ]
    primitive_denominator = int(
        certificate["global_polynomial"]["common_denominator"]
    )
    assert primitive == [
        int(coefficient * primitive_denominator)
        for coefficient in global_polynomial
    ]
    assert gcd(*[abs(value) for value in primitive]) == 1
    assert pgcd(global_polynomial, pderivative(global_polynomial)) == [1]

    sturm = sturm_sequence(global_polynomial)
    total_real_roots = (
        sturm_variations_at_infinity(sturm, False)
        - sturm_variations_at_infinity(sturm, True)
    )
    root_data = certificate["real_roots"]
    assert total_real_roots == int(root_data["count"])
    previous_right = None
    isolated_count = 0
    for interval in root_data["isolating_intervals"]:
        left, right = q(interval["left"]), q(interval["right"])
        multiplicity = int(interval["multiplicity"])
        assert left < right
        assert previous_right is None or previous_right < left
        previous_right = right
        assert peval(global_polynomial, left)
        assert peval(global_polynomial, right)
        count = sturm_variations_at(sturm, left) - sturm_variations_at(sturm, right)
        assert count == multiplicity
        isolated_count += multiplicity
    assert isolated_count == total_real_roots

    witness_prime = int(certificate["irreducibility_witness"]["prime"])
    assert is_prime(witness_prime)
    assert all(
        coefficient.denominator % witness_prime
        for coefficient in global_polynomial
    )
    assert mod_irreducible(
        reduce_rational(global_polynomial, witness_prime), witness_prime
    )

    translation_value = q(certificate["selected_translation"])
    translated = translate(global_polynomial, translation_value)
    derivative_at_translation = peval(
        pderivative(global_polynomial), translation_value
    )
    assert derivative_at_translation
    assert len(translated) > 3 and translated[3]
    seed = list(translated)
    seed[0] = 0
    target = [q(value) for value in certificate["target"]]
    assert target == [
        Fraction(1),
        Fraction(0),
        -2 * peval(global_polynomial, translation_value)
        / derivative_at_translation,
    ]
    identity = certificate["inverse_polynomial_identity"]
    assert [q(value) for value in identity["seed_coefficients_ascending"]] == seed
    assert q(identity["linear_coefficient"]) == seed[1]
    inverse = list(seed)
    inverse[0] -= seed[1] * target[2] / 2
    assert inverse == translated
    assert [
        q(value) for value in identity["inverse_coefficients_ascending"]
    ] == inverse
    assert [
        q(value)
        for value in identity["translated_input_coefficients_ascending"]
    ] == translated

    map_record = certificate["keller_map"]
    stable_parameter = verify_stable_multiplicity_record(
        map_record.get("stable_multiplicity"), degree, seed
    )
    mapping = reconstruct_map(seed, stable_parameter)
    jacobian = [[sderivative(coordinate, index) for index in range(3)] for coordinate in mapping]
    jacobian_determinant = {}
    for permutation in (
        (0, 1, 2),
        (0, 2, 1),
        (1, 0, 2),
        (1, 2, 0),
        (2, 0, 1),
        (2, 1, 0),
    ):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(3)
            for right in range(left + 1, 3)
        )
        term = reduce(
            smul,
            [jacobian[row][permutation[row]] for row in range(3)],
        )
        jacobian_determinant = sadd(
            jacobian_determinant,
            sscale(-1 if inversions % 2 else 1, term),
        )
    assert jacobian_determinant == {(0, 0, 0): Fraction(1)}
    degrees = [
        max(sum(monomial) for monomial in coordinate)
        for coordinate in mapping
    ]
    assert degrees == map_record["coordinate_degrees"]
    assert degree == int(map_record["geometric_degree"])
    hash_record = map_record["expanded_map_hash"]
    assert hash_record["serialization"] == EXPECTED_MAP_SERIALIZATION
    expanded = canonical_map(mapping)
    payload = json.dumps(expanded, ensure_ascii=True, separators=(",", ":")).encode(
        "ascii"
    )
    assert hashlib.sha256(payload).hexdigest() == hash_record["digest"]
    assert [len(coordinate) for coordinate in mapping] == hash_record["term_counts"]

    verify_local_actions(certificate, global_polynomial, local_models)
    print(f"PASS: {path}")
    print("PASS: local stability radii, coefficient CRT, and real root isolation")
    print("PASS: irreducibility and all supplied local action certificates")
    print("PASS: inverse identity, Jacobian one, and expanded map SHA-256")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", nargs="?", type=Path, default=DEFAULT_CERTIFICATE)
    verify(parser.parse_args().certificate)


if __name__ == "__main__":
    main()
