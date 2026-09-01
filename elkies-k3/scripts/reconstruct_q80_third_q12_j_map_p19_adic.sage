#!/usr/bin/env sage -python
"""Interpolate and reconstruct the fixed-gauge third-q12 j-map."""

import argparse
import hashlib
import json
import math
from pathlib import Path

from sage.all import CRT_list, GF, Matrix, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--manifest",
    type=Path,
    default=RESULTS / "q80-third-q12-p19-adic-precision1024-equation-sample-manifest.json",
)
parser.add_argument(
    "--transport",
    type=Path,
    default=RESULTS / "q80-third-q12-long-jacobians-exact-quadratic-gauge.json",
)
parser.add_argument(
    "--source",
    type=Path,
    default=RESULTS / "q80-third-q12-exact-pencil-p19-adic-precision1028.json",
)
parser.add_argument(
    "--operands",
    type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=RESULTS / "q80-third-q12-j-map-p19-adic-reconstructed-qq.json",
)
parser.add_argument("--check", action="store_true")
parser.add_argument(
    "--interpolation-only",
    action="store_true",
    help="stop after p-adic interpolation, held-outs, and the numerator-cube check",
)
args = parser.parse_args()
for name in ("manifest", "transport", "source", "operands", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


manifest = json.loads(args.manifest.read_text())
transport = json.loads(args.transport.read_text())
source = json.loads(args.source.read_text())
operands = json.loads(args.operands.read_text())
if manifest.get("status") != "PASS_EXACT_THIRD_Q12_WEIERSTRASS_P19_ADIC_SAMPLE_BATCH":
    raise ValueError("high-precision sample batch is not certified")
if transport.get("status") != "PASS_EXACT_TRANSPORTED_THIRD_Q12_LONG_JACOBIANS_COMMON_QUADRATIC_GAUGE":
    raise ValueError("finite-prime exact-gauge transport is not certified")
if operands.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_CLOSURE_OPERANDS_P19_HENSEL":
    raise ValueError("exact closure operands are not certified")

prime = 19
digits = int(manifest["specialization"]["digits"])
modulus = prime**digits
omega_square_modular = int(source["quadratic_field"]["omega_square_modulus"]) % modulus
ALL_TRANSPORTED_PRIMES = sorted(int(value) for value in transport["transported_models"])
HELD_OUT_TRANSPORTED_PRIMES = []
RECONSTRUCTION_TRANSPORTED_PRIMES = list(ALL_TRANSPORTED_PRIMES)
ZERO = (0, 0)
ONE = (1, 0)


def c(value):
    return int(value[0]) % modulus, int(value[1]) % modulus


def add(left, right):
    return (left[0] + right[0]) % modulus, (left[1] + right[1]) % modulus


def neg(value):
    return -value[0] % modulus, -value[1] % modulus


def sub(left, right):
    return add(left, neg(right))


def mul(left, right):
    return (
        (left[0] * right[0] + omega_square_modular * left[1] * right[1]) % modulus,
        (left[0] * right[1] + left[1] * right[0]) % modulus,
    )


def power(value, exponent):
    result = ONE
    base = value
    while exponent:
        if exponent & 1:
            result = mul(result, base)
        base = mul(base, base)
        exponent >>= 1
    return result


def is_unit(value):
    return (value[0] * value[0] - omega_square_modular * value[1] * value[1]) % prime != 0


def inverse(value):
    norm = (value[0] * value[0] - omega_square_modular * value[1] * value[1]) % modulus
    if norm % prime == 0:
        raise ZeroDivisionError("non-unit quadratic coefficient")
    inverse_norm = pow(norm, -1, modulus)
    return value[0] * inverse_norm % modulus, -value[1] * inverse_norm % modulus


def divide(left, right):
    return mul(left, inverse(right))


samples = []
seen = set()
for record in manifest["samples"]:
    path = ROOT / record["path"]
    if sha256(path) != record["sha256"]:
        raise ArithmeticError(f"sample hash changed: {path}")
    payload = json.loads(path.read_text())
    base = tuple(payload["specialization"]["base_U_coefficients_1_omega"])
    if base in seen:
        raise ArithmeticError(f"duplicate sample base {base}")
    seen.add(base)
    samples.append((path, payload))
samples.sort(key=lambda item: tuple(reversed(item[1]["specialization"]["base_U_coefficients_1_omega"])))
if len(samples) < 20:
    raise ArithmeticError("20 residue-distinct samples are required: 17 training and 3 held out")
training = samples[:17]
held_out = samples[17:20]


def solve_overdetermined(rows, unknown_count):
    matrix = [[c(value) for value in row] for row in rows]
    rank = 0
    for column in range(unknown_count):
        pivot = next(
            (index for index in range(rank, len(matrix)) if is_unit(matrix[index][column])),
            None,
        )
        if pivot is None:
            raise ArithmeticError(f"interpolation lost a unit pivot at column {column}")
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = inverse(matrix[rank][column])
        matrix[rank] = [mul(value, scale) for value in matrix[rank]]
        for row_index in range(len(matrix)):
            if row_index == rank:
                continue
            factor = matrix[row_index][column]
            if factor == ZERO:
                continue
            matrix[row_index] = [
                sub(matrix[row_index][index], mul(factor, matrix[rank][index]))
                for index in range(unknown_count + 1)
            ]
        rank += 1
    for row in matrix[rank:]:
        if any(value != ZERO for value in row):
            raise ArithmeticError("overdetermined interpolation equations are inconsistent")
    return [matrix[index][-1] for index in range(unknown_count)]


def evaluate_modular(coefficients, value):
    result = ZERO
    for coefficient in reversed(coefficients):
        result = add(mul(result, value), tuple(coefficient))
    return result


def modular_poly_trim(value):
    result = [list(coefficient) for coefficient in value]
    while len(result) > 1 and result[-1] == [0, 0]:
        result.pop()
    return result


def modular_poly_add(left, right):
    result = []
    for index in range(max(len(left), len(right))):
        left_value = tuple(left[index]) if index < len(left) else ZERO
        right_value = tuple(right[index]) if index < len(right) else ZERO
        result.append(list(add(left_value, right_value)))
    return modular_poly_trim(result)


def modular_poly_neg(value):
    return [list(neg(tuple(coefficient))) for coefficient in value]


def modular_poly_sub(left, right):
    return modular_poly_add(left, modular_poly_neg(right))


def modular_poly_scale(value, scalar):
    scalar_pair = (int(scalar) % modulus, 0)
    return [list(mul(tuple(coefficient), scalar_pair)) for coefficient in value]


def modular_poly_mul(left, right):
    result = [ZERO for unused in range(len(left) + len(right) - 1)]
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] = add(
                result[left_index + right_index], mul(tuple(left_value), tuple(right_value))
            )
    return modular_poly_trim([list(value) for value in result])


def modular_poly_power(value, exponent):
    result = [[1, 0]]
    base = value
    while exponent:
        if exponent & 1:
            result = modular_poly_mul(result, base)
        base = modular_poly_mul(base, base)
        exponent >>= 1
    return result


def modular_poly_monic(value):
    value = modular_poly_trim(value)
    inverse_leading = inverse(tuple(value[-1]))
    return [list(mul(tuple(coefficient), inverse_leading)) for coefficient in value]


def modular_poly_divmod(dividend, divisor):
    dividend = modular_poly_trim(dividend)
    divisor = modular_poly_trim(divisor)
    if divisor == [[0, 0]]:
        raise ZeroDivisionError
    if len(dividend) < len(divisor):
        return [[0, 0]], dividend
    quotient = [[0, 0] for unused in range(len(dividend) - len(divisor) + 1)]
    remainder = [list(value) for value in dividend]
    inverse_leading = inverse(tuple(divisor[-1]))
    while remainder != [[0, 0]] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = mul(tuple(remainder[-1]), inverse_leading)
        quotient[shift] = list(coefficient)
        subtractor = [[0, 0] for unused in range(shift)] + [
            list(mul(coefficient, tuple(value))) for value in divisor
        ]
        remainder = modular_poly_sub(remainder, subtractor)
    return modular_poly_trim(quotient), modular_poly_trim(remainder)


def modular_poly_exact_divide(dividend, divisor):
    quotient, remainder = modular_poly_divmod(dividend, divisor)
    if remainder != [[0, 0]]:
        raise ArithmeticError("modular polynomial division is not exact")
    return quotient


def modular_poly_gcd(left, right):
    left = modular_poly_trim(left)
    right = modular_poly_trim(right)
    while right != [[0, 0]]:
        unused, remainder = modular_poly_divmod(left, right)
        left, right = right, remainder
    return modular_poly_monic(left)


def modular_poly_derivative(value):
    if len(value) <= 1:
        return [[0, 0]]
    return [
        list(mul(tuple(value[index]), (index % modulus, 0)))
        for index in range(1, len(value))
    ]


def modular_squarefree_decomposition(value):
    repeated = modular_poly_gcd(value, modular_poly_derivative(value))
    current = modular_poly_exact_divide(value, repeated)
    remaining = repeated
    factors = {}
    multiplicity = 1
    while current != [[1, 0]]:
        overlap = modular_poly_gcd(current, remaining)
        factor = modular_poly_exact_divide(current, overlap)
        if factor != [[1, 0]]:
            factors[multiplicity] = modular_poly_monic(factor)
        current = overlap
        remaining = modular_poly_exact_divide(remaining, overlap)
        multiplicity += 1
        if multiplicity > 24:
            raise ArithmeticError("modular squarefree decomposition did not terminate")
    return factors


# Interpolate only the five small-support long coefficients.  This is the
# cheap invariant route: compose c4 and Delta modulo 19^digits without first
# reconstructing any characteristic-zero Weierstrass coefficient or map.
names = ("a1", "a2", "a3", "a4", "a6")
degree_bounds = dict(zip(names, (2, 4, 4, 6, 8)))
long_coefficients = {}
for coefficient_index, name in enumerate(names):
    coefficient_degree = degree_bounds[name]
    unknown_count = 2 * coefficient_degree + 1
    rows = []
    for unused_path, payload in training:
        u_value = c(payload["specialization"]["base_U_coefficients_1_omega"])
        y_value = c(
            payload["weierstrass"]["a1_a2_a3_a4_a6_mod_19_power_1_omega"][
                coefficient_index
            ]
        )
        powers = [power(u_value, exponent) for exponent in range(coefficient_degree + 1)]
        rows.append(
            powers
            + [
                neg(mul(y_value, powers[exponent]))
                for exponent in range(coefficient_degree)
            ]
            + [mul(y_value, powers[coefficient_degree])]
        )
    solution = solve_overdetermined(rows, unknown_count)
    long_coefficients[name] = {
        "numerator": [list(value) for value in solution[: coefficient_degree + 1]],
        "denominator": [list(value) for value in solution[coefficient_degree + 1 :]]
        + [[1, 0]],
    }

H = long_coefficients["a1"]["denominator"]
denominator_exponents = {"a1": 1, "a2": 2, "a3": 2, "a4": 3, "a6": 4}
for name, exponent in denominator_exponents.items():
    if long_coefficients[name]["denominator"] != modular_poly_power(H, exponent):
        raise ArithmeticError(f"{name}: denominator is not H^{exponent}")
A1, A2, A3, A4, A6 = [long_coefficients[name]["numerator"] for name in names]
B2 = modular_poly_add(modular_poly_power(A1, 2), modular_poly_scale(A2, 4))
B4 = modular_poly_add(modular_poly_mul(A1, A3), modular_poly_scale(A4, 2))
B6 = modular_poly_add(modular_poly_power(A3, 2), modular_poly_scale(A6, 4))
B8 = modular_poly_add(
    modular_poly_add(
        modular_poly_mul(modular_poly_power(A1, 2), A6),
        modular_poly_scale(modular_poly_mul(A2, A6), 4),
    ),
    modular_poly_add(
        modular_poly_neg(modular_poly_mul(modular_poly_mul(A1, A3), A4)),
        modular_poly_sub(
            modular_poly_mul(A2, modular_poly_power(A3, 2)),
            modular_poly_power(A4, 2),
        ),
    ),
)
C4 = modular_poly_sub(
    modular_poly_power(B2, 2),
    modular_poly_scale(modular_poly_mul(B4, H), 24),
)
Delta = modular_poly_add(
    modular_poly_add(
        modular_poly_neg(modular_poly_mul(modular_poly_power(B2, 2), B8)),
        modular_poly_scale(modular_poly_mul(modular_poly_power(B4, 3), H), -8),
    ),
    modular_poly_add(
        modular_poly_scale(
            modular_poly_mul(modular_poly_power(B6, 2), modular_poly_power(H, 2)),
            -27,
        ),
        modular_poly_scale(
            modular_poly_mul(
                modular_poly_mul(modular_poly_mul(B2, B4), B6), H
            ),
            9,
        ),
    ),
)
j_numerator = modular_poly_power(C4, 3)
j_denominator = modular_poly_mul(Delta, modular_poly_power(H, 2))
if len(j_numerator) != 25 or len(j_denominator) != 25:
    raise ArithmeticError(
        f"unexpected composed j degrees {len(j_numerator)-1},{len(j_denominator)-1}"
    )
normalization = inverse(tuple(j_denominator[-1]))
interpolated = {
    "degrees_numerator_denominator": [24, 24],
    "numerator": [list(mul(tuple(value), normalization)) for value in j_numerator],
    "denominator": [list(mul(tuple(value), normalization)) for value in j_denominator],
}
degree = 24
modular_denominator_factors = modular_squarefree_decomposition(
    interpolated["denominator"]
)
modular_factor_degrees = {
    multiplicity: len(factor) - 1
    for multiplicity, factor in modular_denominator_factors.items()
}
if modular_factor_degrees != {1: 8, 2: 3, 4: 1, 6: 1}:
    raise ArithmeticError(
        f"unexpected modular denominator multiplicities: {modular_factor_degrees}"
    )


modular_numerator_leading = tuple(interpolated["numerator"][-1])
modular_normalized_numerator = [
    list(divide(tuple(coefficient), modular_numerator_leading))
    for coefficient in interpolated["numerator"]
]
modular_cube_root = [[0, 0] for unused in range(9)]
modular_cube_root[-1] = [1, 0]
inverse_three = pow(3, -1, modulus)
for index in range(7, -1, -1):
    partial_cube = modular_poly_power(modular_cube_root, 3)
    target_degree = 16 + index
    current = tuple(partial_cube[target_degree])
    difference = sub(tuple(modular_normalized_numerator[target_degree]), current)
    modular_cube_root[index] = [
        difference[0] * inverse_three % modulus,
        difference[1] * inverse_three % modulus,
    ]
if [
    list(mul(modular_numerator_leading, tuple(coefficient)))
    for coefficient in modular_poly_power(modular_cube_root, 3)
] != interpolated["numerator"]:
    raise ArithmeticError("the interpolated p-adic j numerator is not a scalar times a cube")

expected_p19_j = transport["transported_models"]["19"]["j"]
for exact_key, finite_key in (
    ("numerator", "numerator_coefficients_low_to_high_1_omega"),
    ("denominator", "denominator_coefficients_low_to_high_1_omega"),
):
    reduction = [[value % prime for value in pair] for pair in interpolated[exact_key]]
    if reduction != expected_p19_j[finite_key]:
        raise ArithmeticError(f"interpolated {exact_key} does not replay the transported p=19 j-map")
for path, payload in held_out:
    u_value = c(payload["specialization"]["base_U_coefficients_1_omega"])
    expected = c(payload["weierstrass"]["j_mod_19_power_1_omega"])
    numerator_value = evaluate_modular(interpolated["numerator"], u_value)
    denominator_value = evaluate_modular(interpolated["denominator"], u_value)
    if not is_unit(denominator_value) or divide(numerator_value, denominator_value) != expected:
        raise ArithmeticError(f"held-out j replay failed at {path}")
if args.interpolation_only:
    print(
        f"Q80THIRDQ12JINTERPOLATE|digits={digits}|samples=20|heldout=3|"
        "degrees=24,24|cube_degree=8|status=PASS_P_ADIC_INTERPOLATION_ONLY"
    )
    raise SystemExit(0)


def rational_mod(value, local_prime):
    numerator = int(value.numerator()) % local_prime
    denominator = int(value.denominator()) % local_prime
    if not denominator:
        raise ZeroDivisionError
    return numerator * pow(denominator, -1, local_prime) % local_prime


def validate_candidate_at_finite_primes(candidate, selected_primes=None):
    selected_primes = set(selected_primes or ALL_TRANSPORTED_PRIMES)
    for local_prime_text, model in transport["transported_models"].items():
        local_prime = int(local_prime_text)
        if local_prime not in selected_primes:
            continue
        expected = model["j"]
        for exact_key, finite_key in (
            ("numerator", "numerator_coefficients_low_to_high_1_omega"),
            ("denominator", "denominator_coefficients_low_to_high_1_omega"),
        ):
            reduction = [
                [rational_mod(QQ(coordinate), local_prime) for coordinate in pair]
                for pair in candidate[exact_key]
            ]
            if reduction != expected[finite_key]:
                return False
    return True


def reconstruct_projectively(record):
    residues = []
    positions = []
    for key in ("numerator", "denominator"):
        for coefficient_index, pair in enumerate(record[key]):
            if key == "denominator" and coefficient_index == len(record[key]) - 1:
                if pair != [1, 0]:
                    raise ArithmeticError("j denominator is not normalized")
                continue
            for omega_index, value in enumerate(pair):
                residues.append(ZZ(value) % modulus)
                positions.append((key, coefficient_index, omega_index))
    dimension = len(residues) + 1
    lattice = Matrix(ZZ, dimension, dimension)
    for index in range(dimension - 1):
        lattice[index, index] = modulus
    for index, value in enumerate(residues):
        lattice[dimension - 1, index] = value
    lattice[dimension - 1, dimension - 1] = 1
    reduced = lattice.LLL(delta=0.99)
    diagnostics = []
    for row in sorted(reduced.rows(), key=lambda value: value.dot_product(value)):
        row = list(row)
        if not row[-1] or row[-1] % prime == 0:
            continue
        common = math.gcd(*(abs(int(value)) for value in row))
        if common > 1:
            row = [value // common for value in row]
        if row[-1] < 0:
            row = [-value for value in row]
        candidate = {
            "degrees_numerator_denominator": [degree, degree],
            "numerator": [[QQ(0), QQ(0)] for unused in record["numerator"]],
            "denominator": [[QQ(0), QQ(0)] for unused in record["denominator"]],
        }
        candidate["denominator"][-1] = [QQ(1), QQ(0)]
        for index, (key, coefficient_index, omega_index) in enumerate(positions):
            candidate[key][coefficient_index][omega_index] = QQ(row[index]) / QQ(row[-1])
        maximum_bits = max(abs(ZZ(value)).nbits() for value in row)
        diagnostics.append(maximum_bits)
        if validate_candidate_at_finite_primes(candidate, RECONSTRUCTION_TRANSPORTED_PRIMES):
            return candidate, {
                "method": "projective LLL with seven-prime acceptance",
                "lattice_dimension": dimension,
                "maximum_primitive_coordinate_bits": maximum_bits,
                "random_lattice_boundary_bits": int(
                    math.ceil(ZZ(modulus).nbits() * (dimension - 1) / dimension)
                ),
                "short_rows_tested": len(diagnostics),
                "validated_primes": RECONSTRUCTION_TRANSPORTED_PRIMES,
            }
    raise ArithmeticError(
        "j: no projective LLL row validates at all aligned primes; "
        f"short-row bits={diagnostics[:8]}"
    )


def reconstruct_pair_block(modular_pairs, validator, label, lattice_modulus=modulus):
    lattice_modulus = ZZ(lattice_modulus)
    residues = [ZZ(value) % lattice_modulus for pair in modular_pairs for value in pair]
    dimension = len(residues) + 1
    lattice = Matrix(ZZ, dimension, dimension)
    for index in range(dimension - 1):
        lattice[index, index] = lattice_modulus
    for index, value in enumerate(residues):
        lattice[dimension - 1, index] = value
    lattice[dimension - 1, dimension - 1] = 1
    reduced = lattice.LLL(delta=0.99)
    diagnostics = []
    for row in sorted(reduced.rows(), key=lambda value: value.dot_product(value)):
        row = list(row)
        if not row[-1] or row[-1] % prime == 0:
            continue
        common = math.gcd(*(abs(int(value)) for value in row))
        if common > 1:
            row = [value // common for value in row]
        if row[-1] < 0:
            row = [-value for value in row]
        pairs = [
            [QQ(row[2 * index]) / QQ(row[-1]), QQ(row[2 * index + 1]) / QQ(row[-1])]
            for index in range(len(modular_pairs))
        ]
        maximum_bits = max(abs(ZZ(value)).nbits() for value in row)
        diagnostics.append(maximum_bits)
        if validator(pairs):
            return pairs, {
                "method": (
                    "structured projective LLL over p-adic plus CRT modulus"
                    if lattice_modulus != modulus
                    else "structured projective LLL over p-adic modulus"
                ),
                "lattice_dimension": dimension,
                "maximum_primitive_coordinate_bits": maximum_bits,
                "random_lattice_boundary_bits": int(
                    math.ceil(lattice_modulus.nbits() * (dimension - 1) / dimension)
                ),
                "short_rows_tested": len(diagnostics),
                "validated_primes": RECONSTRUCTION_TRANSPORTED_PRIMES,
            }
    raise ArithmeticError(
        f"{label}: no structured projective LLL row validates; "
        f"short-row bits={diagnostics[:8]}"
    )


def reconstruct_scalar_asymmetric(residue, validator, label, lattice_modulus=modulus):
    lattice_modulus = ZZ(lattice_modulus)
    old_remainder, remainder = lattice_modulus, ZZ(residue) % lattice_modulus
    old_cofactor, cofactor = ZZ(0), ZZ(1)
    diagnostics = []
    while remainder:
        if cofactor and cofactor % prime:
            candidate = QQ(remainder) / QQ(cofactor)
            numerator_bits = abs(ZZ(candidate.numerator())).nbits()
            denominator_bits = ZZ(candidate.denominator()).nbits()
            diagnostics.append([numerator_bits, denominator_bits])
            if validator(candidate):
                return candidate, {
                    "method": "extended-Euclidean asymmetric rational reconstruction",
                    "numerator_bits": numerator_bits,
                    "denominator_bits": denominator_bits,
                    "bit_sum": numerator_bits + denominator_bits,
                    "euclidean_candidates_tested": len(diagnostics),
                    "validated_primes": RECONSTRUCTION_TRANSPORTED_PRIMES,
                }
        quotient = old_remainder // remainder
        old_remainder, remainder = remainder, old_remainder - quotient * remainder
        old_cofactor, cofactor = cofactor, old_cofactor - quotient * cofactor
    raise ArithmeticError(
        f"{label}: no asymmetric convergent validates; tail={diagnostics[-8:]}"
    )


def rational_record(record):
    return QQ(ZZ(record["numerator"])) / QQ(ZZ(record["denominator"]))


q1 = rational_record(operands["biquadratic_field"]["q1"])
q2 = rational_record(operands["biquadratic_field"]["q2"])
omega_square_exact = QQ(16) * q1 * q2
quadratic_discriminant = QQ(4) * omega_square_exact


def eadd(left, right):
    return [left[0] + right[0], left[1] + right[1]]


def eneg(value):
    return [-value[0], -value[1]]


def esub(left, right):
    return eadd(left, eneg(right))


def emul(left, right):
    return [
        left[0] * right[0] + omega_square_exact * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    ]


def einv(value):
    norm = value[0] * value[0] - omega_square_exact * value[1] * value[1]
    if not norm:
        raise ZeroDivisionError
    return [value[0] / norm, -value[1] / norm]


def ediv(left, right):
    return emul(left, einv(right))


EZERO = [QQ(0), QQ(0)]
EONE = [QQ(1), QQ(0)]


def poly_trim(value):
    result = [list(coefficient) for coefficient in value]
    while len(result) > 1 and result[-1] == EZERO:
        result.pop()
    return result


def poly_add(left, right):
    result = [list(EZERO) for unused in range(max(len(left), len(right)))]
    for index in range(len(result)):
        a = left[index] if index < len(left) else EZERO
        b = right[index] if index < len(right) else EZERO
        result[index] = eadd(a, b)
    return poly_trim(result)


def poly_neg(value):
    return [eneg(coefficient) for coefficient in value]


def poly_sub(left, right):
    return poly_add(left, poly_neg(right))


def poly_mul(left, right):
    result = [list(EZERO) for unused in range(len(left) + len(right) - 1)]
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] = eadd(
                result[left_index + right_index], emul(left_value, right_value)
            )
    return poly_trim(result)


def poly_power(value, exponent):
    result = [list(EONE)]
    base = value
    while exponent:
        if exponent & 1:
            result = poly_mul(result, base)
        base = poly_mul(base, base)
        exponent >>= 1
    return result


def poly_monic(value):
    inverse_lead = einv(poly_trim(value)[-1])
    return [emul(coefficient, inverse_lead) for coefficient in poly_trim(value)]


def poly_divmod(dividend, divisor):
    dividend = poly_trim(dividend)
    divisor = poly_trim(divisor)
    if divisor == [EZERO]:
        raise ZeroDivisionError
    if len(dividend) < len(divisor):
        return [list(EZERO)], dividend
    quotient = [list(EZERO) for unused in range(len(dividend) - len(divisor) + 1)]
    remainder = [list(coefficient) for coefficient in dividend]
    inverse_lead = einv(divisor[-1])
    while remainder != [EZERO] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = emul(remainder[-1], inverse_lead)
        quotient[shift] = coefficient
        subtractor = [list(EZERO) for unused in range(shift)] + [
            emul(coefficient, value) for value in divisor
        ]
        remainder = poly_sub(remainder, subtractor)
    return poly_trim(quotient), poly_trim(remainder)


def poly_exact_divide(dividend, divisor):
    quotient, remainder = poly_divmod(dividend, divisor)
    if remainder != [EZERO]:
        raise ArithmeticError("polynomial division is not exact")
    return quotient


def poly_gcd(left, right):
    left = poly_trim(left)
    right = poly_trim(right)
    while right != [EZERO]:
        unused, remainder = poly_divmod(left, right)
        left, right = right, remainder
    return poly_monic(left)


def poly_derivative(value):
    if len(value) <= 1:
        return [list(EZERO)]
    return [[QQ(index) * a, QQ(index) * b] for index, (a, b) in enumerate(value)][1:]


def finite_cube_root_from_numerator(coefficients, local_prime):
    local_omega_square = rational_mod(omega_square_exact, local_prime)

    def local_add(left, right):
        return [
            (left[0] + right[0]) % local_prime,
            (left[1] + right[1]) % local_prime,
        ]

    def local_sub(left, right):
        return [
            (left[0] - right[0]) % local_prime,
            (left[1] - right[1]) % local_prime,
        ]

    def local_mul(left, right):
        return [
            (left[0] * right[0] + local_omega_square * left[1] * right[1])
            % local_prime,
            (left[0] * right[1] + left[1] * right[0]) % local_prime,
        ]

    def local_inverse(value):
        norm = (value[0] * value[0] - local_omega_square * value[1] * value[1]) % local_prime
        inverse_norm = pow(norm, -1, local_prime)
        return [value[0] * inverse_norm % local_prime, -value[1] * inverse_norm % local_prime]

    def local_poly_mul(left, right):
        result = [[0, 0] for unused in range(len(left) + len(right) - 1)]
        for left_index, left_value in enumerate(left):
            for right_index, right_value in enumerate(right):
                result[left_index + right_index] = local_add(
                    result[left_index + right_index], local_mul(left_value, right_value)
                )
        return result

    def local_poly_cube(value):
        return local_poly_mul(local_poly_mul(value, value), value)

    leading_coefficient = coefficients[-1]
    inverse_leading = local_inverse(leading_coefficient)
    normalized = [local_mul(value, inverse_leading) for value in coefficients]
    result = [[0, 0] for unused in range(9)]
    result[-1] = [1, 0]
    inverse_local_three = pow(3, -1, local_prime)
    for index in range(7, -1, -1):
        partial = local_poly_cube(result)
        target_degree = 16 + index
        difference = local_sub(normalized[target_degree], partial[target_degree])
        result[index] = [
            difference[0] * inverse_local_three % local_prime,
            difference[1] * inverse_local_three % local_prime,
        ]
    if [local_mul(leading_coefficient, value) for value in local_poly_cube(result)] != coefficients:
        raise ArithmeticError(f"finite j numerator is not a cube at p={local_prime}")
    return result


finite_denominator_factor_cache = {}


def finite_denominator_factors(local_prime):
    if local_prime in finite_denominator_factor_cache:
        return finite_denominator_factor_cache[local_prime]
    finite = GF(local_prime)
    z_ring = PolynomialRing(finite, "z")
    z = z_ring.gen()
    local_omega_square = finite(omega_square_exact)
    quadratic = GF(
        local_prime**2,
        "w",
        modulus=z**2 - local_omega_square,
    )
    w = quadratic.gen()
    polynomial_ring = PolynomialRing(quadratic, "U")
    record = transport["transported_models"][str(local_prime)]["j"]
    denominator_polynomial = polynomial_ring(
        [
            quadratic(a) + quadratic(b) * w
            for a, b in record["denominator_coefficients_low_to_high_1_omega"]
        ]
    )
    grouped = {}
    for factor, multiplicity in denominator_polynomial.factor():
        grouped[multiplicity] = grouped.get(multiplicity, polynomial_ring.one()) * factor.monic()
    encoded = {}
    for multiplicity, factor in grouped.items():
        encoded[multiplicity] = []
        for coefficient in factor.monic().list():
            coordinates = list(coefficient.polynomial())
            coordinates += [finite(0)] * (2 - len(coordinates))
            encoded[multiplicity].append([int(coordinates[0]), int(coordinates[1])])
    if {key: len(value) - 1 for key, value in encoded.items()} != {1: 8, 2: 3, 4: 1, 6: 1}:
        raise ArithmeticError(f"unexpected finite denominator factors at p={local_prime}")
    finite_denominator_factor_cache[local_prime] = encoded
    return encoded


CRT_AUXILIARY_PRIMES = [
    value
    for value in RECONSTRUCTION_TRANSPORTED_PRIMES
    if value != prime
]
reconstruction_modulus = ZZ(modulus)
for local_prime in CRT_AUXILIARY_PRIMES:
    reconstruction_modulus *= local_prime


def crt_extend_pair(padic_pair, finite_pair_getter):
    result = []
    moduli = [ZZ(modulus)] + [ZZ(value) for value in CRT_AUXILIARY_PRIMES]
    for omega_index, padic_residue in enumerate(padic_pair):
        residues = [ZZ(padic_residue)] + [
            ZZ(finite_pair_getter(local_prime)[omega_index])
            for local_prime in CRT_AUXILIARY_PRIMES
        ]
        result.append(int(CRT_list(residues, moduli)))
    return result


crt_cube_root_lower = []
for coefficient_index, padic_pair in enumerate(modular_cube_root[:-1]):
    def finite_cube_coefficient(local_prime, index=coefficient_index):
        model = transport["transported_models"][str(local_prime)]
        return finite_cube_root_from_numerator(
            model["j"]["numerator_coefficients_low_to_high_1_omega"], local_prime
        )[index]

    crt_cube_root_lower.append(crt_extend_pair(padic_pair, finite_cube_coefficient))


def validate_cube_root(candidate_lower, selected_primes=RECONSTRUCTION_TRANSPORTED_PRIMES):
    selected_primes = set(selected_primes)
    candidate = candidate_lower + [[QQ(1), QQ(0)]]
    for local_prime_text, model in transport["transported_models"].items():
        local_prime = int(local_prime_text)
        if local_prime not in selected_primes:
            continue
        try:
            reduction = [
                [rational_mod(value, local_prime) for value in pair]
                for pair in candidate
            ]
        except ZeroDivisionError:
            return False
        expected = finite_cube_root_from_numerator(
            model["j"]["numerator_coefficients_low_to_high_1_omega"], local_prime
        )
        if reduction != expected:
            return False
    return True


try:
    exact_cube_root_lower, cube_reconstruction = reconstruct_pair_block(
        crt_cube_root_lower,
        validate_cube_root,
        "monic degree-8 c4 factor",
        lattice_modulus=reconstruction_modulus,
    )
except ArithmeticError as joint_error:
    exact_cube_root_lower = []
    coefficient_diagnostics = []
    for coefficient_index, modular_pair in enumerate(crt_cube_root_lower):
        def validate_cube_root_coefficient(candidate, index=coefficient_index):
            for local_prime_text, model in transport["transported_models"].items():
                local_prime = int(local_prime_text)
                if local_prime not in RECONSTRUCTION_TRANSPORTED_PRIMES:
                    continue
                expected = finite_cube_root_from_numerator(
                    model["j"]["numerator_coefficients_low_to_high_1_omega"],
                    local_prime,
                )[index]
                reduction = [rational_mod(value, local_prime) for value in candidate[0]]
                if reduction != expected:
                    return False
            return True

        try:
            reconstructed, diagnostic = reconstruct_pair_block(
                [modular_pair],
                validate_cube_root_coefficient,
                f"c4 factor coefficient {coefficient_index}",
                lattice_modulus=reconstruction_modulus,
            )
            exact_pair = reconstructed[0]
            diagnostic["fallback_level"] = "coefficient_pair"
        except ArithmeticError as pair_error:
            exact_pair = []
            scalar_diagnostics = []
            for omega_index, residue in enumerate(modular_pair):
                def validate_scalar(
                    candidate,
                    index=coefficient_index,
                    coordinate=omega_index,
                ):
                    for local_prime_text, model in transport["transported_models"].items():
                        local_prime = int(local_prime_text)
                        if local_prime not in RECONSTRUCTION_TRANSPORTED_PRIMES:
                            continue
                        expected = finite_cube_root_from_numerator(
                            model["j"]["numerator_coefficients_low_to_high_1_omega"],
                            local_prime,
                        )[index][coordinate]
                        try:
                            reduction = rational_mod(candidate, local_prime)
                        except ZeroDivisionError:
                            return False
                        if reduction != expected:
                            return False
                    return True

                exact_coordinate, scalar_diagnostic = reconstruct_scalar_asymmetric(
                    residue,
                    validate_scalar,
                    f"c4 factor coefficient {coefficient_index} coordinate {omega_index}",
                    lattice_modulus=reconstruction_modulus,
                )
                exact_pair.append(exact_coordinate)
                scalar_diagnostics.append(scalar_diagnostic)
            diagnostic = {
                "fallback_level": "asymmetric_scalar",
                "pair_failure": str(pair_error),
                "coordinates": scalar_diagnostics,
            }
        exact_cube_root_lower.append(exact_pair)
        coefficient_diagnostics.append(diagnostic)
    if not validate_cube_root(exact_cube_root_lower):
        raise ArithmeticError("separate c4-factor coefficient reconstruction failed full replay")
    cube_reconstruction = {
        "method": "separate coefficient reconstruction after joint boundary failure",
        "joint_failure": str(joint_error),
        "coefficients": coefficient_diagnostics,
    }
cube_root = exact_cube_root_lower + [[QQ(1), QQ(0)]]
if HELD_OUT_TRANSPORTED_PRIMES and not validate_cube_root(
    exact_cube_root_lower, HELD_OUT_TRANSPORTED_PRIMES
):
    raise ArithmeticError("reconstructed degree-8 cube root failed held-out-prime replay")


factor_multiplicities = (1, 2, 4, 6)
def finite_j_leading(local_prime):
    return transport["transported_models"][str(local_prime)]["j"][
        "numerator_coefficients_low_to_high_1_omega"
    ][-1]


modular_invariant_bundle = [
    crt_extend_pair(interpolated["numerator"][-1], finite_j_leading)
]
for multiplicity in factor_multiplicities:
    for coefficient_index, padic_pair in enumerate(
        modular_denominator_factors[multiplicity][:-1]
    ):
        def finite_factor_coefficient(
            local_prime,
            factor_multiplicity=multiplicity,
            index=coefficient_index,
        ):
            return finite_denominator_factors(local_prime)[factor_multiplicity][index]

        modular_invariant_bundle.append(
            crt_extend_pair(padic_pair, finite_factor_coefficient)
        )


def decode_invariant_bundle(candidate):
    candidate_leading = candidate[0]
    position = 1
    candidate_factors = {}
    for multiplicity in factor_multiplicities:
        factor_degree = modular_factor_degrees[multiplicity]
        candidate_factors[multiplicity] = candidate[position : position + factor_degree] + [
            [QQ(1), QQ(0)]
        ]
        position += factor_degree
    if position != len(candidate):
        raise ArithmeticError("invariant factor bundle has an unexpected length")
    candidate_denominator = [list(EONE)]
    for multiplicity in factor_multiplicities:
        candidate_denominator = poly_mul(
            candidate_denominator,
            poly_power(candidate_factors[multiplicity], multiplicity),
        )
    return candidate_leading, candidate_factors, candidate_denominator


def validate_invariant_bundle(candidate):
    candidate_leading, unused_factors, candidate_denominator = decode_invariant_bundle(
        candidate
    )
    candidate_numerator = [
        emul(candidate_leading, coefficient) for coefficient in poly_power(cube_root, 3)
    ]
    return validate_candidate_at_finite_primes(
        {
            "numerator": candidate_numerator,
            "denominator": candidate_denominator,
        },
        RECONSTRUCTION_TRANSPORTED_PRIMES,
    )


invariant_bundle, scale_reconstruction = reconstruct_pair_block(
    modular_invariant_bundle,
    validate_invariant_bundle,
    "j leading scalar and squarefree denominator factors",
    lattice_modulus=reconstruction_modulus,
)
leading, reconstructed_denominator_factors, denominator = decode_invariant_bundle(
    invariant_bundle
)
numerator = [emul(leading, coefficient) for coefficient in poly_power(cube_root, 3)]
exact = {
    "degrees_numerator_denominator": [degree, degree],
    "numerator": numerator,
    "denominator": denominator,
}
if not validate_candidate_at_finite_primes(exact, RECONSTRUCTION_TRANSPORTED_PRIMES):
    raise ArithmeticError("structured exact j candidate failed reconstruction-prime replay")
if HELD_OUT_TRANSPORTED_PRIMES and not validate_candidate_at_finite_primes(
    exact, HELD_OUT_TRANSPORTED_PRIMES
):
    raise ArithmeticError("structured exact j candidate failed held-out-prime replay")
reconstruction = {
    "method": "j=c4^3/Delta structured reconstruction",
    "combined_reconstruction_modulus_bits": int(reconstruction_modulus.nbits()),
    "CRT_auxiliary_primes": CRT_AUXILIARY_PRIMES,
    "monic_degree8_cube_root": cube_reconstruction,
    "numerator_leading_scalar_and_squarefree_denominator_factors": scale_reconstruction,
    "unstructured_dimension_99_fallback_used": False,
}


# Verify the invariant cube shape N = leading(N) * C(U)^3 with C monic of degree 8.
if [emul(leading, coefficient) for coefficient in poly_power(cube_root, 3)] != numerator:
    raise ArithmeticError("the reconstructed j numerator is not a scalar times a polynomial cube")


# The monic denominator factors were reconstructed separately.  Their exact
# powered product is the denominator, and their mod-19 reductions are the
# pairwise-coprime squarefree factors returned by the modular gcd tower.
# Good reduction therefore certifies pairwise coprimality over QQ(omega)
# without repeating a very large rational-function gcd computation.
squarefree_factors = reconstructed_denominator_factors
expected_factor_degrees = {1: 8, 2: 3, 4: 1, 6: 1}
actual_factor_degrees = {
    multiplicity: len(factor) - 1 for multiplicity, factor in squarefree_factors.items()
}
if actual_factor_degrees != expected_factor_degrees:
    raise ArithmeticError(
        f"unexpected exact denominator multiplicities: {actual_factor_degrees}"
    )
replayed_denominator = [list(EONE)]
for multiplicity, factor in squarefree_factors.items():
    replayed_denominator = poly_mul(
        replayed_denominator, poly_power(factor, multiplicity)
    )
if replayed_denominator != denominator:
    raise ArithmeticError("exact squarefree-factor product does not replay the denominator")


def encode_pair(pair):
    return {
        "coefficients_1_omega": [str(value) for value in pair],
        "trace": str(QQ(2) * pair[0]),
        "anti_invariant_coefficient": str(pair[1]),
    }


def encode_polynomial(coefficients):
    return [encode_pair(pair) for pair in coefficients]


output = {
    "schema": "elkies-k3.q80-third-q12-j-map-p19-adic-reconstructed-qq.v1",
    "status": "PASS_CANDIDATE_EXACT_THIRD_Q12_J_MAP_RECONSTRUCTION_QQ",
    "specialization": {
        "u": "-2",
        "coefficient_field_basis": ["1", "omega"],
        "omega_square": str(omega_square_exact),
        "quadratic_polynomial": f"X^2-({omega_square_exact})",
        "quadratic_discriminant": str(quadratic_discriminant),
    },
    "gauge": {
        "base_coordinate": "the exact common U coordinate pinned by the transported models",
        "base_PGL2_change_applied": False,
        "weierstrass_scaling_dependence": False,
        "section_sign_dependence": False,
        "denominator_normalization": "monic of degree 24 in U",
    },
    "j_map": {
        "degrees_numerator_denominator": [degree, degree],
        "numerator_coefficients_low_to_high_U": encode_polynomial(exact["numerator"]),
        "denominator_coefficients_low_to_high_U": encode_polynomial(exact["denominator"]),
    },
    "exact_invariant_structure": {
        "numerator_leading_scalar": encode_pair(leading),
        "monic_c4_factor_degree8_coefficients_low_to_high_U": encode_polynomial(cube_root),
        "identity": "numerator = numerator_leading_scalar * monic_c4_factor^3",
        "denominator_squarefree_factor_degrees_by_multiplicity": {
            str(key): value for key, value in sorted(actual_factor_degrees.items())
        },
        "denominator_squarefree_factors": {
            str(key): encode_polynomial(value)
            for key, value in sorted(squarefree_factors.items())
        },
        "fibre_multiplicity_shape": "I6 + I4 + 3 I2 + 8 simple residual discriminant points",
    },
    "reconstruction": reconstruction,
    "validation": {
        "p19_digits": digits,
        "training_samples": len(training),
        "held_out_samples": len(held_out),
        "independent_aligned_primes": sorted(
            int(value) for value in transport["transported_models"] if value != "19"
        ),
        "transported_primes_used_for_reconstruction_acceptance": RECONSTRUCTION_TRANSPORTED_PRIMES,
        "held_out_transported_primes": HELD_OUT_TRANSPORTED_PRIMES,
        "exact_numerator_cube_identity": True,
        "exact_denominator_squarefree_decomposition": True,
        "literal_exact_pencil_or_weierstrass_replay": False,
    },
    "inputs": {
        "manifest": {"path": str(args.manifest.relative_to(ROOT)), "sha256": sha256(args.manifest)},
        "transport": {"path": str(args.transport.relative_to(ROOT)), "sha256": sha256(args.transport)},
        "source": {"path": str(args.source.relative_to(ROOT)), "sha256": sha256(args.source)},
        "operands": {"path": str(args.operands.relative_to(ROOT)), "sha256": sha256(args.operands)},
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "unique modular interpolation of the five long coefficients at their certified small degree bounds through 19^digits",
            "symbolic modular composition of c4, Delta, and the degree-(24,24) j-map in the pinned U gauge",
            "three held-out certified p-adic sample replays",
            "one structured projective LLL candidate using the p-adic residue and all six aligned auxiliary primes",
            "the displayed characteristic-zero numerator cube and denominator squarefree identities",
        ],
        "not_proved": [
            "replay at a new aligned prime not consumed by the CRT reconstruction",
            "literal characteristic-zero substitution into the exact genus-one pencil or a reconstructed Jacobian",
            "a characteristic-zero fibre theorem or transported component marking",
            "an equation over QQ rather than the displayed quadratic coefficient field",
        ],
    },
    "reproduce": "sage -python elkies-k3/scripts/reconstruct_q80_third_q12_j_map_p19_adic.sage",
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"reconstructed j-map artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    f"Q80THIRDQ12JRECONSTRUCT|digits={digits}|samples={len(training) + len(held_out)}|"
    "primes=19,61,67,83,89,103,131|"
    "status=PASS_CANDIDATE_EXACT_THIRD_Q12_J_MAP_RECONSTRUCTION_QQ"
)
