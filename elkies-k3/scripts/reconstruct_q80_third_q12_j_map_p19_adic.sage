#!/usr/bin/env sage -python
"""Interpolate and reconstruct the fixed-gauge third-q12 j-map."""

import argparse
import hashlib
import json
import math
from pathlib import Path

from sage.all import CRT_list, GF, Matrix, PolynomialRing, QQ, ZZ, inverse_mod


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
parser.add_argument(
    "--holdout-prime",
    action="append",
    default=[],
    type=int,
    help=(
        "transported prime to exclude from CRT/LLL and reserve for literal replay; "
        "repeat for multiple holdouts"
    ),
)
parser.add_argument(
    "--reconstruction-granularity",
    choices=(
        "auto",
        "bundle",
        "component",
        "quadratic-projective",
        "pair",
        "scalar",
    ),
    default="auto",
    help=(
        "shared-denominator bundle reconstruction, separate projective blocks for "
        "the rational and omega components, one projective factor over the quadratic "
        "field with a two-coordinate scale, one algebraic coefficient pair at a "
        "time, or independent rational coordinates"
    ),
)
parser.add_argument(
    "--base-normalization",
    choices=("pinned", "i6-i4", "i6-i4-i2trace"),
    default="pinned",
    help=(
        "reconstruct in the pinned U coordinate or in the intrinsic coordinate "
        "z=L6(U)/L4(U), which sends the I6 and I4 fibres to zero and infinity; "
        "i6-i4-i2trace also removes the remaining scaling with the cubic I2 trace"
    ),
)
parser.add_argument(
    "--intrinsic-basis",
    choices=("monomial", "evaluations", "joint-evaluations"),
    default="monomial",
    help=(
        "reconstruct intrinsic factor coefficients directly, reconstruct each exact "
        "value separately at residue-distinct p-adic sample nodes, or reconstruct all "
        "node values in one joint projective lattice before exact interpolation"
    ),
)
parser.add_argument(
    "--c4-pivot",
    type=int,
    default=8,
    help=(
        "coefficient index 0..8 used to normalize the projective degree-eight "
        "c4 factor before reconstruction; 8 is the usual monic chart"
    ),
)
args = parser.parse_args()
for name in ("manifest", "transport", "source", "operands", "output"):
    setattr(args, name, getattr(args, name).resolve())
if not 0 <= args.c4_pivot <= 8:
    raise ValueError("--c4-pivot must lie between 0 and 8")
if args.c4_pivot != 8 and args.intrinsic_basis != "monomial":
    raise ValueError("non-leading c4 pivots currently require --intrinsic-basis monomial")


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
HELD_OUT_TRANSPORTED_PRIMES = sorted(set(args.holdout_prime))
unknown_holdouts = sorted(set(HELD_OUT_TRANSPORTED_PRIMES) - set(ALL_TRANSPORTED_PRIMES))
if unknown_holdouts:
    raise ValueError(f"holdout primes are absent from the transported models: {unknown_holdouts}")
if prime in HELD_OUT_TRANSPORTED_PRIMES:
    raise ValueError("p=19 supplies the p-adic residue and cannot be a transported-prime holdout")
RECONSTRUCTION_TRANSPORTED_PRIMES = [
    value for value in ALL_TRANSPORTED_PRIMES if value not in HELD_OUT_TRANSPORTED_PRIMES
]
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


raw_modular_denominator_factors = modular_squarefree_decomposition(
    interpolated["denominator"]
)
raw_factor_degrees = {
    multiplicity: len(factor) - 1
    for multiplicity, factor in raw_modular_denominator_factors.items()
}
if raw_factor_degrees != {1: 8, 2: 3, 4: 1, 6: 1}:
    raise ArithmeticError(
        f"unexpected raw modular denominator multiplicities: {raw_factor_degrees}"
    )


def modular_poly_pair_scale(value, scalar):
    return [list(mul(tuple(coefficient), tuple(scalar))) for coefficient in value]


def modular_poly_linear_fraction_transform(value, source_degree, r6, r4):
    """Return (z-1)^source_degree * value((r6-z*r4)/(z-1))."""
    numerator_linear = [list(r6), list(neg(tuple(r4)))]
    denominator_linear = [[-1 % modulus, 0], [1, 0]]
    result = [[0, 0]]
    for index, coefficient in enumerate(value):
        term = modular_poly_mul(
            modular_poly_power(numerator_linear, index),
            modular_poly_power(denominator_linear, source_degree - index),
        )
        result = modular_poly_add(
            result, modular_poly_pair_scale(term, coefficient)
        )
    return modular_poly_trim(result)


def modular_poly_variable_scale(value, scalar):
    return [
        list(mul(tuple(coefficient), power(tuple(scalar), index)))
        for index, coefficient in enumerate(value)
    ]


base_normalization = {
    "kind": "pinned",
    "coordinate": "U",
    "transformation": None,
}
expected_factor_degrees = {1: 8, 2: 3, 4: 1, 6: 1}
if args.base_normalization in ("i6-i4", "i6-i4-i2trace"):
    r6 = raw_modular_denominator_factors[6][0]
    r4 = raw_modular_denominator_factors[4][0]
    transformed_numerator = modular_poly_linear_fraction_transform(
        interpolated["numerator"], degree, r6, r4
    )
    transformed_denominator = modular_poly_linear_fraction_transform(
        interpolated["denominator"], degree, r6, r4
    )
    transformed_normalization = inverse(tuple(transformed_denominator[-1]))
    interpolated = {
        "degrees_numerator_denominator": [
            len(transformed_numerator) - 1,
            len(transformed_denominator) - 1,
        ],
        "numerator": [
            list(mul(tuple(value), transformed_normalization))
            for value in transformed_numerator
        ],
        "denominator": [
            list(mul(tuple(value), transformed_normalization))
            for value in transformed_denominator
        ],
    }
    if interpolated["degrees_numerator_denominator"] != [24, 20]:
        raise ArithmeticError(
            "the intrinsic I6/I4 normalization did not produce j degrees 24/20"
        )
    i2_trace_scale = None
    if args.base_normalization == "i6-i4-i2trace":
        first_normalized_factors = modular_squarefree_decomposition(
            interpolated["denominator"]
        )
        i2_trace_scale = first_normalized_factors[2][2]
        if not is_unit(tuple(i2_trace_scale)):
            raise ArithmeticError("the cubic I2 trace scale is not a p-adic unit")
        rescaled_numerator = modular_poly_variable_scale(
            interpolated["numerator"], i2_trace_scale
        )
        rescaled_denominator = modular_poly_variable_scale(
            interpolated["denominator"], i2_trace_scale
        )
        rescaled_normalization = inverse(tuple(rescaled_denominator[-1]))
        interpolated = {
            "degrees_numerator_denominator": [24, 20],
            "numerator": [
                list(mul(tuple(value), rescaled_normalization))
                for value in rescaled_numerator
            ],
            "denominator": [
                list(mul(tuple(value), rescaled_normalization))
                for value in rescaled_denominator
            ],
        }
    base_normalization = {
        "kind": "intrinsic_i6_i4",
        "coordinate": "w" if i2_trace_scale is not None else "z",
        "transformation": (
            "w=(L6(U)/L4(U))/a2, where a2 is the z^2 coefficient of the monic cubic I2 factor"
            if i2_trace_scale is not None
            else "z=L6(U)/L4(U)=(U+r6)/(U+r4)"
        ),
        "p19_power_r6_r4_coefficients_1_omega": [r6, r4],
        "p19_power_i2_trace_scale_coefficients_1_omega": i2_trace_scale,
        "distinguished_fibres": {
            ("w=0" if i2_trace_scale is not None else "z=0"): "I6",
            ("w=infinity" if i2_trace_scale is not None else "z=infinity"): "I4",
        },
    }
    expected_factor_degrees = {1: 8, 2: 3, 6: 1}

denominator_degree = len(interpolated["denominator"]) - 1
modular_denominator_factors = modular_squarefree_decomposition(
    interpolated["denominator"]
)
modular_factor_degrees = {
    multiplicity: len(factor) - 1
    for multiplicity, factor in modular_denominator_factors.items()
}
if modular_factor_degrees != expected_factor_degrees:
    raise ArithmeticError(
        f"unexpected normalized modular denominator multiplicities: {modular_factor_degrees}"
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
    raise ArithmeticError("the normalized p-adic j numerator is not a scalar times a cube")

c4_pivot_inverse = inverse(tuple(modular_cube_root[args.c4_pivot]))
modular_cube_root_in_pivot_chart = [
    list(mul(tuple(coefficient), c4_pivot_inverse))
    for coefficient in modular_cube_root
]
if modular_cube_root_in_pivot_chart[args.c4_pivot] != [1, 0]:
    raise ArithmeticError("the selected p-adic c4 pivot did not normalize to one")
c4_free_indices = [index for index in range(9) if index != args.c4_pivot]


def rational_mod(value, local_prime):
    numerator = int(value.numerator()) % local_prime
    denominator = int(value.denominator()) % local_prime
    if not denominator:
        raise ZeroDivisionError
    return numerator * pow(denominator, -1, local_prime) % local_prime


finite_normalized_j_cache = {}


def finite_normalized_j(local_prime):
    if local_prime in finite_normalized_j_cache:
        return finite_normalized_j_cache[local_prime]
    record = transport["transported_models"][str(local_prime)]["j"]
    if args.base_normalization == "pinned":
        finite_normalized_j_cache[local_prime] = record
        return record
    finite = GF(local_prime)
    x_ring = PolynomialRing(finite, "x")
    x = x_ring.gen()
    local_omega_square = finite(omega_square_exact)
    quadratic = GF(local_prime**2, "w", modulus=x**2 - local_omega_square)
    w = quadratic.gen()
    z_ring = PolynomialRing(quadratic, "z")
    z = z_ring.gen()

    def decode(coefficients):
        return z_ring([quadratic(a) + quadratic(b) * w for a, b in coefficients])

    numerator_polynomial = decode(
        record["numerator_coefficients_low_to_high_1_omega"]
    )
    denominator_polynomial = decode(
        record["denominator_coefficients_low_to_high_1_omega"]
    )
    grouped = {}
    for factor, multiplicity in denominator_polynomial.factor():
        grouped[multiplicity] = grouped.get(multiplicity, z_ring.one()) * factor.monic()
    r6 = grouped[6].monic()[0]
    r4 = grouped[4].monic()[0]
    numerator_linear = quadratic(r6) - quadratic(r4) * z
    denominator_linear = z - 1

    def transform(polynomial):
        return sum(
            coefficient * numerator_linear**index * denominator_linear ** (degree - index)
            for index, coefficient in enumerate(polynomial.list())
        )

    transformed_numerator = transform(numerator_polynomial)
    transformed_denominator = transform(denominator_polynomial)
    scale = transformed_denominator.leading_coefficient() ** -1
    transformed_numerator *= scale
    transformed_denominator *= scale
    if args.base_normalization == "i6-i4-i2trace":
        transformed_grouped = {}
        for factor, multiplicity in transformed_denominator.factor():
            transformed_grouped[multiplicity] = (
                transformed_grouped.get(multiplicity, z_ring.one()) * factor.monic()
            )
        i2_trace_scale = transformed_grouped[2].monic()[2]
        transformed_numerator = z_ring(
            [
                coefficient * i2_trace_scale**index
                for index, coefficient in enumerate(transformed_numerator.list())
            ]
        )
        transformed_denominator = z_ring(
            [
                coefficient * i2_trace_scale**index
                for index, coefficient in enumerate(transformed_denominator.list())
            ]
        )
        scale = transformed_denominator.leading_coefficient() ** -1
        transformed_numerator *= scale
        transformed_denominator *= scale

    def encode(polynomial):
        encoded = []
        for coefficient in polynomial.list():
            coordinates = list(coefficient.polynomial())
            coordinates += [finite(0)] * (2 - len(coordinates))
            encoded.append([int(coordinates[0]), int(coordinates[1])])
        return encoded

    result = {
        "numerator_coefficients_low_to_high_1_omega": encode(transformed_numerator),
        "denominator_coefficients_low_to_high_1_omega": encode(transformed_denominator),
    }
    if [len(result["numerator_coefficients_low_to_high_1_omega"]) - 1,
        len(result["denominator_coefficients_low_to_high_1_omega"]) - 1] != [24, 20]:
        raise ArithmeticError(
            f"finite I6/I4 normalization has unexpected degrees at p={local_prime}"
        )
    finite_normalized_j_cache[local_prime] = result
    return result


def validate_candidate_at_finite_primes(candidate, selected_primes=None):
    selected_primes = set(selected_primes or ALL_TRANSPORTED_PRIMES)
    for local_prime_text, model in transport["transported_models"].items():
        local_prime = int(local_prime_text)
        if local_prime not in selected_primes:
            continue
        expected = finite_normalized_j(local_prime)
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
            "degrees_numerator_denominator": [degree, denominator_degree],
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


def reconstruct_quadratic_projective_factor(
    modular_free_pairs,
    validator,
    label,
    lattice_modulus=modulus,
):
    """Reconstruct one projective polynomial over QQ(omega).

    The selected pivot coefficient supplies a two-coordinate scale ``u+v*omega``.
    For every other coefficient ``a+b*omega`` the congruences are

        x = a*u + D*b*v,   y = b*u + a*v  (mod M),

    where ``D=omega^2``.  Thus the lattice has sixteen modulus rows and two
    coupled scale rows, rather than imposing one rational scale on all
    quadratic coordinates or reconstructing the two components independently.
    """
    lattice_modulus = ZZ(lattice_modulus)
    if len(modular_free_pairs) != 8:
        raise ValueError("a degree-eight projective factor needs eight free pairs")
    omega_numerator = ZZ(omega_square_exact.numerator())
    omega_denominator = ZZ(omega_square_exact.denominator())
    if math.gcd(int(omega_denominator), int(lattice_modulus)) != 1:
        raise ZeroDivisionError("quadratic discriminant denominator meets lattice modulus")
    omega_square_modulus = (
        omega_numerator * inverse_mod(omega_denominator, lattice_modulus)
    ) % lattice_modulus
    dimension = 18
    lattice = Matrix(ZZ, dimension, dimension)
    for index in range(16):
        lattice[index, index] = lattice_modulus
    scale_one_row = 16
    scale_omega_row = 17
    for pair_index, pair in enumerate(modular_free_pairs):
        a = ZZ(pair[0]) % lattice_modulus
        b = ZZ(pair[1]) % lattice_modulus
        lattice[scale_one_row, 2 * pair_index] = a
        lattice[scale_one_row, 2 * pair_index + 1] = b
        lattice[scale_omega_row, 2 * pair_index] = omega_square_modulus * b
        lattice[scale_omega_row, 2 * pair_index + 1] = a
    lattice[scale_one_row, 16] = 1
    lattice[scale_omega_row, 17] = 1
    reduced = lattice.LLL(delta=0.99)
    diagnostics = []
    for row in sorted(reduced.rows(), key=lambda value: value.dot_product(value)):
        row = list(row)
        pivot_pair = [QQ(row[16]), QQ(row[17])]
        if pivot_pair == EZERO:
            continue
        common = math.gcd(*(abs(int(value)) for value in row))
        if common > 1:
            row = [value // common for value in row]
            pivot_pair = [QQ(row[16]), QQ(row[17])]
        first_nonzero = next((value for value in row if value), ZZ(1))
        if first_nonzero < 0:
            row = [-value for value in row]
            pivot_pair = [QQ(row[16]), QQ(row[17])]
        try:
            candidate_free = [
                ediv([QQ(row[2 * index]), QQ(row[2 * index + 1])], pivot_pair)
                for index in range(8)
            ]
        except ZeroDivisionError:
            continue
        maximum_bits = max(abs(ZZ(value)).nbits() for value in row)
        diagnostics.append(maximum_bits)
        if validator(candidate_free):
            return candidate_free, {
                "method": "projective LLL over QQ(omega) with a coupled two-coordinate scale",
                "lattice_dimension": dimension,
                "congruence_count": 16,
                "maximum_primitive_coordinate_bits": maximum_bits,
                "random_lattice_boundary_bits": int(
                    math.ceil(lattice_modulus.nbits() * 16 / dimension)
                ),
                "short_rows_tested": len(diagnostics),
                "validated_primes": RECONSTRUCTION_TRANSPORTED_PRIMES,
            }
    raise ArithmeticError(
        f"{label}: no quadratic-projective LLL row validates; "
        f"short-row bits={diagnostics[:8]}"
    )


def reconstruct_scalar_block(residues, validator, label, lattice_modulus=modulus):
    lattice_modulus = ZZ(lattice_modulus)
    residues = [ZZ(value) % lattice_modulus for value in residues]
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
        candidate = [QQ(value) / QQ(row[-1]) for value in row[:-1]]
        maximum_bits = max(abs(ZZ(value)).nbits() for value in row)
        diagnostics.append(maximum_bits)
        if validator(candidate):
            return candidate, {
                "method": "projective LLL for one quadratic-basis component",
                "lattice_dimension": dimension,
                "maximum_primitive_coordinate_bits": maximum_bits,
                "random_lattice_boundary_bits": int(
                    math.ceil(lattice_modulus.nbits() * (dimension - 1) / dimension)
                ),
                "short_rows_tested": len(diagnostics),
                "validated_primes": RECONSTRUCTION_TRANSPORTED_PRIMES,
            }
    raise ArithmeticError(
        f"{label}: no component-wise projective LLL row validates; "
        f"short-row bits={diagnostics[:8]}"
    )


def reconstruct_scalar_asymmetric(residue, validator, label, lattice_modulus=modulus):
    lattice_modulus = ZZ(lattice_modulus)
    old_remainder, remainder = lattice_modulus, ZZ(residue) % lattice_modulus
    old_cofactor, cofactor = ZZ(0), ZZ(1)
    diagnostics = []
    candidates = []
    while remainder:
        if cofactor and cofactor % prime:
            candidate = QQ(remainder) / QQ(cofactor)
            numerator_bits = abs(ZZ(candidate.numerator())).nbits()
            denominator_bits = ZZ(candidate.denominator()).nbits()
            diagnostics.append([numerator_bits, denominator_bits])
            if validator(candidate):
                candidates.append(
                    (
                        max(numerator_bits, denominator_bits),
                        numerator_bits + denominator_bits,
                        candidate,
                        numerator_bits,
                        denominator_bits,
                    )
                )
        quotient = old_remainder // remainder
        old_remainder, remainder = remainder, old_remainder - quotient * remainder
        old_cofactor, cofactor = cofactor, old_cofactor - quotient * cofactor
    if not candidates:
        raise ArithmeticError(
            f"{label}: no asymmetric convergent validates; tail={diagnostics[-8:]}"
        )
    unused_maximum, unused_sum, candidate, numerator_bits, denominator_bits = min(candidates)
    return candidate, {
        "method": "minimum-height extended-Euclidean asymmetric rational reconstruction",
        "numerator_bits": numerator_bits,
        "denominator_bits": denominator_bits,
        "bit_sum": numerator_bits + denominator_bits,
        "euclidean_candidates_tested": len(diagnostics),
        "valid_convergents": len(candidates),
        "validated_primes": RECONSTRUCTION_TRANSPORTED_PRIMES,
    }


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


def finite_c4_factor_in_pivot_chart(local_prime):
    coefficients = finite_cube_root_from_numerator(
        finite_normalized_j(local_prime)[
            "numerator_coefficients_low_to_high_1_omega"
        ],
        local_prime,
    )
    local_omega_square = rational_mod(omega_square_exact, local_prime)

    def local_mul(left, right):
        return [
            (left[0] * right[0] + local_omega_square * left[1] * right[1])
            % local_prime,
            (left[0] * right[1] + left[1] * right[0]) % local_prime,
        ]

    pivot = coefficients[args.c4_pivot]
    norm = (pivot[0] * pivot[0] - local_omega_square * pivot[1] * pivot[1]) % local_prime
    if not norm:
        raise ZeroDivisionError(
            f"c4 coefficient {args.c4_pivot} is not a unit at p={local_prime}"
        )
    inverse_norm = pow(norm, -1, local_prime)
    pivot_inverse = [
        pivot[0] * inverse_norm % local_prime,
        -pivot[1] * inverse_norm % local_prime,
    ]
    normalized = [local_mul(value, pivot_inverse) for value in coefficients]
    if normalized[args.c4_pivot] != [1, 0]:
        raise ArithmeticError("finite c4 pivot did not normalize to one")
    return normalized


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
    record = finite_normalized_j(local_prime)
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
    if {key: len(value) - 1 for key, value in encoded.items()} != expected_factor_degrees:
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
for coefficient_index in c4_free_indices:
    padic_pair = modular_cube_root_in_pivot_chart[coefficient_index]

    def finite_cube_coefficient(local_prime, index=coefficient_index):
        return finite_c4_factor_in_pivot_chart(local_prime)[index]

    crt_cube_root_lower.append(crt_extend_pair(padic_pair, finite_cube_coefficient))


def decode_c4_pivot_chart(candidate_free):
    if len(candidate_free) != 8:
        raise ArithmeticError("a degree-eight projective c4 factor needs eight free coefficients")
    candidate_in_chart = []
    position = 0
    for coefficient_index in range(9):
        if coefficient_index == args.c4_pivot:
            candidate_in_chart.append([QQ(1), QQ(0)])
        else:
            candidate_in_chart.append(list(candidate_free[position]))
            position += 1
    leading = candidate_in_chart[-1]
    if leading == EZERO:
        raise ZeroDivisionError("the reconstructed c4 factor lost degree eight")
    return [ediv(coefficient, leading) for coefficient in candidate_in_chart]


def validate_cube_root(candidate_free, selected_primes=RECONSTRUCTION_TRANSPORTED_PRIMES):
    selected_primes = set(selected_primes)
    try:
        candidate = decode_c4_pivot_chart(candidate_free)
    except (ArithmeticError, ZeroDivisionError):
        return False
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
            finite_normalized_j(local_prime)[
                "numerator_coefficients_low_to_high_1_omega"
            ],
            local_prime,
        )
        if reduction != expected:
            return False
    return True


def exact_pair_power(value, exponent):
    result = list(EONE)
    base = list(value)
    while exponent:
        if exponent & 1:
            result = emul(result, base)
        base = emul(base, base)
        exponent >>= 1
    return result


def exact_pair_linear_solve(rows, unknown_count):
    matrix = [[list(value) for value in row] for row in rows]
    rank = 0
    for column in range(unknown_count):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column] != EZERO),
            None,
        )
        if pivot is None:
            raise ArithmeticError(f"exact evaluation interpolation lost pivot {column}")
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse_pivot = einv(matrix[rank][column])
        matrix[rank] = [emul(value, inverse_pivot) for value in matrix[rank]]
        for row_index in range(len(matrix)):
            if row_index == rank:
                continue
            factor = matrix[row_index][column]
            if factor == EZERO:
                continue
            matrix[row_index] = [
                esub(matrix[row_index][index], emul(factor, matrix[rank][index]))
                for index in range(unknown_count + 1)
            ]
        rank += 1
    return [matrix[index][-1] for index in range(unknown_count)]


def finite_pair_evaluate(coefficients, exact_node, local_prime):
    local_omega_square = rational_mod(omega_square_exact, local_prime)

    def local_mul(left, right):
        return [
            (left[0] * right[0] + local_omega_square * left[1] * right[1])
            % local_prime,
            (left[0] * right[1] + left[1] * right[0]) % local_prime,
        ]

    node = [rational_mod(value, local_prime) for value in exact_node]
    result = [0, 0]
    for coefficient in reversed(coefficients):
        product = local_mul(result, node)
        result = [
            (product[0] + coefficient[0]) % local_prime,
            (product[1] + coefficient[1]) % local_prime,
        ]
    return result


exact_cube_root_lower = None
cube_reconstruction = None
if args.intrinsic_basis in ("evaluations", "joint-evaluations"):
    if args.base_normalization != "pinned":
        raise ValueError("evaluation-basis reconstruction currently requires the exact pinned U nodes")
    evaluation_nodes = [
        [QQ(value) for value in payload["specialization"]["base_U_coefficients_1_omega"]]
        for unused_path, payload in training[:8]
    ]
    evaluation_records = []
    for evaluation_index, exact_node in enumerate(evaluation_nodes):
        padic_node = c([int(value) for value in exact_node])
        padic_value = evaluate_modular(modular_cube_root, padic_node)

        def finite_evaluation(local_prime, node=exact_node):
            coefficients = finite_cube_root_from_numerator(
                finite_normalized_j(local_prime)[
                    "numerator_coefficients_low_to_high_1_omega"
                ],
                local_prime,
            )
            return finite_pair_evaluate(coefficients, node, local_prime)

        crt_value = crt_extend_pair(padic_value, finite_evaluation)
        evaluation_records.append((exact_node, crt_value, finite_evaluation))

    reconstructed_values = []
    evaluation_diagnostics = []
    if args.intrinsic_basis == "joint-evaluations":
        if args.reconstruction_granularity == "component":
            reconstructed_components = []
            component_diagnostics = []
            for omega_index in range(2):
                component_residues = [
                    crt_value[omega_index]
                    for unused_node, crt_value, unused_finite in evaluation_records
                ]

                def validate_evaluation_component(candidate, coordinate=omega_index):
                    for candidate_value, (unused_node, unused_crt, finite_evaluation) in zip(
                        candidate, evaluation_records
                    ):
                        for local_prime in RECONSTRUCTION_TRANSPORTED_PRIMES:
                            try:
                                reduction = rational_mod(candidate_value, local_prime)
                            except ZeroDivisionError:
                                return False
                            if reduction != finite_evaluation(local_prime)[coordinate]:
                                return False
                    return True

                component, diagnostic = reconstruct_scalar_block(
                    component_residues,
                    validate_evaluation_component,
                    f"joint c4 factor evaluations quadratic-basis component {omega_index}",
                    lattice_modulus=reconstruction_modulus,
                )
                reconstructed_components.append(component)
                component_diagnostics.append(diagnostic)
            reconstructed_values = [
                [reconstructed_components[0][index], reconstructed_components[1][index]]
                for index in range(len(evaluation_records))
            ]
            evaluation_diagnostics = {
                "method": (
                    "two joint projective evaluation lattices for the rational and "
                    "omega components"
                ),
                "components": component_diagnostics,
            }
        else:
            crt_values = [
                crt_value for unused_node, crt_value, unused_finite in evaluation_records
            ]

            def validate_joint_evaluations(candidate):
                for candidate_value, (unused_node, unused_crt, finite_evaluation) in zip(
                    candidate, evaluation_records
                ):
                    for local_prime in RECONSTRUCTION_TRANSPORTED_PRIMES:
                        try:
                            reduction = [
                                rational_mod(value, local_prime) for value in candidate_value
                            ]
                        except ZeroDivisionError:
                            return False
                        if reduction != finite_evaluation(local_prime):
                            return False
                interpolation_rows = []
                for node, value in zip(evaluation_nodes, candidate):
                    powers = [exact_pair_power(node, exponent) for exponent in range(9)]
                    interpolation_rows.append(powers[:8] + [esub(value, powers[8])])
                try:
                    candidate_lower = exact_pair_linear_solve(interpolation_rows, 8)
                except (ArithmeticError, ZeroDivisionError):
                    return False
                return validate_cube_root(candidate_lower)

            reconstructed_values, joint_diagnostic = reconstruct_pair_block(
                crt_values,
                validate_joint_evaluations,
                "joint residue-distinct c4 factor evaluations",
                lattice_modulus=reconstruction_modulus,
            )
            evaluation_diagnostics = {
                "method": "one joint projective lattice for all quadratic node values",
                "joint_lattice": joint_diagnostic,
            }
    else:
        for evaluation_index, (exact_node, crt_value, finite_evaluation) in enumerate(
            evaluation_records
        ):

            def validate_evaluation(candidate):
                for local_prime in RECONSTRUCTION_TRANSPORTED_PRIMES:
                    expected = finite_evaluation(local_prime)
                    try:
                        reduction = [
                            rational_mod(value, local_prime) for value in candidate[0]
                        ]
                    except ZeroDivisionError:
                        return False
                    if reduction != expected:
                        return False
                return True

            if args.reconstruction_granularity == "scalar":
                reconstructed_value = []
                scalar_diagnostics = []
                for omega_index, residue in enumerate(crt_value):
                    def validate_evaluation_scalar(candidate):
                        return math.gcd(
                            int(candidate.denominator()), int(reconstruction_modulus)
                        ) == 1

                    coordinate, scalar_diagnostic = reconstruct_scalar_asymmetric(
                        residue,
                        validate_evaluation_scalar,
                        f"c4 factor evaluation {evaluation_index} coordinate {omega_index}",
                        lattice_modulus=reconstruction_modulus,
                    )
                    reconstructed_value.append(coordinate)
                    scalar_diagnostics.append(scalar_diagnostic)
                diagnostic = {
                    "method": "independent scalar evaluation reconstruction",
                    "coordinates": scalar_diagnostics,
                }
            else:
                reconstructed, diagnostic = reconstruct_pair_block(
                    [crt_value],
                    validate_evaluation,
                    f"c4 factor evaluation {evaluation_index}",
                    lattice_modulus=reconstruction_modulus,
                )
                reconstructed_value = reconstructed[0]
            reconstructed_values.append(reconstructed_value)
            evaluation_diagnostics.append(diagnostic)
    interpolation_rows = []
    for node, value in zip(evaluation_nodes, reconstructed_values):
        powers = [exact_pair_power(node, exponent) for exponent in range(9)]
        interpolation_rows.append(powers[:8] + [esub(value, powers[8])])
    exact_cube_root_lower = exact_pair_linear_solve(interpolation_rows, 8)
    if not validate_cube_root(exact_cube_root_lower):
        raise ArithmeticError("evaluation-basis c4-factor reconstruction failed full replay")
    cube_reconstruction = {
        "method": (
            "joint residue-distinct intrinsic evaluations followed by exact interpolation"
            if args.intrinsic_basis == "joint-evaluations"
            else "residue-distinct intrinsic evaluations followed by exact interpolation"
        ),
        "evaluation_nodes_coefficients_1_omega": [
            [[str(value) for value in pair] for pair in evaluation_nodes]
        ][0],
        "evaluations": evaluation_diagnostics,
    }

if (
    exact_cube_root_lower is None
    and args.reconstruction_granularity == "quadratic-projective"
):
    exact_cube_root_lower, quadratic_projective_diagnostic = (
        reconstruct_quadratic_projective_factor(
            crt_cube_root_lower,
            validate_cube_root,
            f"degree-eight c4 factor in coefficient-{args.c4_pivot} chart",
            lattice_modulus=reconstruction_modulus,
        )
    )
    cube_reconstruction = {
        "method": "one projective factor over the exact quadratic coefficient field",
        "quadratic_projective_lattice": quadratic_projective_diagnostic,
    }

if exact_cube_root_lower is None and args.reconstruction_granularity == "component":
    reconstructed_components = []
    component_diagnostics = []
    for omega_index in range(2):
        component_residues = [pair[omega_index] for pair in crt_cube_root_lower]

        def validate_cube_component(candidate, coordinate=omega_index):
            for local_prime in RECONSTRUCTION_TRANSPORTED_PRIMES:
                expected = finite_c4_factor_in_pivot_chart(local_prime)
                try:
                    reduction = [rational_mod(value, local_prime) for value in candidate]
                except ZeroDivisionError:
                    return False
                if reduction != [expected[index][coordinate] for index in c4_free_indices]:
                    return False
            return True

        component, diagnostic = reconstruct_scalar_block(
            component_residues,
            validate_cube_component,
            f"c4 factor quadratic-basis component {omega_index}",
            lattice_modulus=reconstruction_modulus,
        )
        reconstructed_components.append(component)
        component_diagnostics.append(diagnostic)
    exact_cube_root_lower = [
        [reconstructed_components[0][index], reconstructed_components[1][index]]
        for index in range(len(crt_cube_root_lower))
    ]
    if not validate_cube_root(exact_cube_root_lower):
        raise ArithmeticError("component-wise c4-factor reconstruction failed full replay")
    cube_reconstruction = {
        "method": "separate projective reconstruction of rational and omega polynomial components",
        "components": component_diagnostics,
    }

joint_error = None
if exact_cube_root_lower is None and args.reconstruction_granularity in ("auto", "bundle"):
    try:
        exact_cube_root_lower, cube_reconstruction = reconstruct_pair_block(
            crt_cube_root_lower,
            validate_cube_root,
            f"degree-8 c4 factor in coefficient-{args.c4_pivot} chart",
            lattice_modulus=reconstruction_modulus,
        )
    except ArithmeticError as error:
        if args.reconstruction_granularity == "bundle":
            raise
        joint_error = error
elif exact_cube_root_lower is None:
    joint_error = ArithmeticError(
        f"bundle reconstruction disabled by --reconstruction-granularity="
        f"{args.reconstruction_granularity}"
    )
if exact_cube_root_lower is None:
    exact_cube_root_lower = []
    coefficient_diagnostics = []
    for coefficient_index, modular_pair in enumerate(crt_cube_root_lower):
        actual_coefficient_index = c4_free_indices[coefficient_index]

        def validate_cube_root_coefficient(candidate, index=actual_coefficient_index):
            for local_prime_text, model in transport["transported_models"].items():
                local_prime = int(local_prime_text)
                if local_prime not in RECONSTRUCTION_TRANSPORTED_PRIMES:
                    continue
                expected = finite_c4_factor_in_pivot_chart(local_prime)[index]
                reduction = [rational_mod(value, local_prime) for value in candidate[0]]
                if reduction != expected:
                    return False
            return True

        pair_error = None
        if args.reconstruction_granularity != "scalar":
            try:
                reconstructed, diagnostic = reconstruct_pair_block(
                    [modular_pair],
                    validate_cube_root_coefficient,
                    f"c4 factor coefficient {actual_coefficient_index}",
                    lattice_modulus=reconstruction_modulus,
                )
                exact_pair = reconstructed[0]
                diagnostic["fallback_level"] = "coefficient_pair"
            except ArithmeticError as error:
                if args.reconstruction_granularity == "pair":
                    raise
                pair_error = error
        else:
            pair_error = ArithmeticError(
                "coefficient-pair reconstruction disabled by scalar mode"
            )
        if pair_error is not None:
            exact_pair = []
            scalar_diagnostics = []
            for omega_index, residue in enumerate(modular_pair):
                def validate_scalar(
                    candidate,
                    index=actual_coefficient_index,
                    coordinate=omega_index,
                ):
                    # Every convergent already satisfies the one combined CRT
                    # congruence.  Only denominator invertibility is needed here;
                    # validate_cube_root performs the full finite-prime replay
                    # after all coordinates have been selected.
                    return math.gcd(
                        int(candidate.denominator()), int(reconstruction_modulus)
                    ) == 1

                exact_coordinate, scalar_diagnostic = reconstruct_scalar_asymmetric(
                    residue,
                    validate_scalar,
                    f"c4 factor coefficient {actual_coefficient_index} coordinate {omega_index}",
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
        "method": "separate intrinsic-coefficient reconstruction",
        "requested_granularity": args.reconstruction_granularity,
        "joint_failure": str(joint_error),
        "coefficients": coefficient_diagnostics,
    }
cube_root = decode_c4_pivot_chart(exact_cube_root_lower)
if HELD_OUT_TRANSPORTED_PRIMES and not validate_cube_root(
    exact_cube_root_lower, HELD_OUT_TRANSPORTED_PRIMES
):
    raise ArithmeticError("reconstructed degree-8 cube root failed held-out-prime replay")


factor_multiplicities = tuple(sorted(expected_factor_degrees))
def finite_j_leading(local_prime):
    return finite_normalized_j(local_prime)[
        "numerator_coefficients_low_to_high_1_omega"
    ][-1]


modular_invariant_bundle = [crt_extend_pair(interpolated["numerator"][-1], finite_j_leading)]
modular_invariant_specs = [("numerator_leading", None, None)]
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
        modular_invariant_specs.append(
            ("denominator_factor", multiplicity, coefficient_index)
        )


def finite_invariant_pair(spec, local_prime):
    kind, multiplicity, coefficient_index = spec
    if kind == "numerator_leading":
        return finite_j_leading(local_prime)
    return finite_denominator_factors(local_prime)[multiplicity][coefficient_index]


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


invariant_bundle = None
invariant_joint_error = None
if args.reconstruction_granularity == "component":
    reconstructed_components = []
    component_diagnostics = []
    for omega_index in range(2):
        component_residues = [pair[omega_index] for pair in modular_invariant_bundle]

        def validate_invariant_component(candidate, coordinate=omega_index):
            for local_prime in RECONSTRUCTION_TRANSPORTED_PRIMES:
                expected = [
                    finite_invariant_pair(spec, local_prime)[coordinate]
                    for spec in modular_invariant_specs
                ]
                try:
                    reduction = [rational_mod(value, local_prime) for value in candidate]
                except ZeroDivisionError:
                    return False
                if reduction != expected:
                    return False
            return True

        component, diagnostic = reconstruct_scalar_block(
            component_residues,
            validate_invariant_component,
            f"intrinsic invariant quadratic-basis component {omega_index}",
            lattice_modulus=reconstruction_modulus,
        )
        reconstructed_components.append(component)
        component_diagnostics.append(diagnostic)
    invariant_bundle = [
        [reconstructed_components[0][index], reconstructed_components[1][index]]
        for index in range(len(modular_invariant_bundle))
    ]
    if not validate_invariant_bundle(invariant_bundle):
        raise ArithmeticError("component-wise intrinsic-invariant reconstruction failed full replay")
    scale_reconstruction = {
        "method": "separate projective reconstruction of rational and omega invariant components",
        "components": component_diagnostics,
    }
elif args.reconstruction_granularity in ("auto", "bundle"):
    try:
        invariant_bundle, scale_reconstruction = reconstruct_pair_block(
            modular_invariant_bundle,
            validate_invariant_bundle,
            "j leading scalar and squarefree denominator factors",
            lattice_modulus=reconstruction_modulus,
        )
    except ArithmeticError as error:
        if args.reconstruction_granularity == "bundle":
            raise
        invariant_joint_error = error
else:
    invariant_joint_error = ArithmeticError(
        f"bundle reconstruction disabled by --reconstruction-granularity="
        f"{args.reconstruction_granularity}"
    )
if invariant_bundle is None:
    invariant_bundle = []
    invariant_diagnostics = []
    for bundle_index, (modular_pair, spec) in enumerate(
        zip(modular_invariant_bundle, modular_invariant_specs)
    ):
        def validate_invariant_pair(candidate, local_spec=spec):
            for local_prime in RECONSTRUCTION_TRANSPORTED_PRIMES:
                expected = finite_invariant_pair(local_spec, local_prime)
                try:
                    reduction = [rational_mod(value, local_prime) for value in candidate[0]]
                except ZeroDivisionError:
                    return False
                if reduction != expected:
                    return False
            return True

        pair_error = None
        if args.reconstruction_granularity != "scalar":
            try:
                reconstructed, diagnostic = reconstruct_pair_block(
                    [modular_pair],
                    validate_invariant_pair,
                    f"intrinsic invariant coefficient {bundle_index}",
                    lattice_modulus=reconstruction_modulus,
                )
                exact_pair = reconstructed[0]
                diagnostic["fallback_level"] = "coefficient_pair"
            except ArithmeticError as error:
                if args.reconstruction_granularity == "pair":
                    raise
                pair_error = error
        else:
            pair_error = ArithmeticError(
                "coefficient-pair reconstruction disabled by scalar mode"
            )
        if pair_error is not None:
            exact_pair = []
            scalar_diagnostics = []
            for omega_index, residue in enumerate(modular_pair):
                def validate_invariant_scalar(
                    candidate,
                    local_spec=spec,
                    coordinate=omega_index,
                ):
                    # As above, the combined residue already encodes every
                    # reconstruction prime.  The full intrinsic bundle is
                    # replayed after these scalar choices are assembled.
                    return math.gcd(
                        int(candidate.denominator()), int(reconstruction_modulus)
                    ) == 1

                exact_coordinate, scalar_diagnostic = reconstruct_scalar_asymmetric(
                    residue,
                    validate_invariant_scalar,
                    f"intrinsic invariant coefficient {bundle_index} coordinate {omega_index}",
                    lattice_modulus=reconstruction_modulus,
                )
                exact_pair.append(exact_coordinate)
                scalar_diagnostics.append(scalar_diagnostic)
            diagnostic = {
                "fallback_level": "asymmetric_scalar",
                "pair_failure": str(pair_error),
                "coordinates": scalar_diagnostics,
            }
        invariant_bundle.append(exact_pair)
        invariant_diagnostics.append(diagnostic)
    if not validate_invariant_bundle(invariant_bundle):
        raise ArithmeticError("separate intrinsic-invariant reconstruction failed full replay")
    scale_reconstruction = {
        "method": "separate intrinsic-coefficient reconstruction",
        "requested_granularity": args.reconstruction_granularity,
        "joint_failure": str(invariant_joint_error),
        "coefficients": invariant_diagnostics,
    }
leading, reconstructed_denominator_factors, denominator = decode_invariant_bundle(
    invariant_bundle
)
numerator = [emul(leading, coefficient) for coefficient in poly_power(cube_root, 3)]
exact = {
    "degrees_numerator_denominator": [degree, denominator_degree],
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
    "requested_granularity": args.reconstruction_granularity,
    "intrinsic_basis": args.intrinsic_basis,
    "c4_projective_pivot_coefficient_index": args.c4_pivot,
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


proved_claims = [
    "unique modular interpolation of the five long coefficients at their certified small degree bounds through 19^digits",
    "symbolic modular composition of c4, Delta, and the degree-(24,24) j-map in the pinned U gauge",
    "three held-out certified p-adic sample replays",
    (
        "one structured projective reconstruction candidate using the p-adic residue and "
        f"{len(CRT_AUXILIARY_PRIMES)} aligned auxiliary primes"
    ),
    "the displayed characteristic-zero numerator cube and denominator squarefree identities",
]
not_proved_claims = [
    "literal characteristic-zero substitution into the exact genus-one pencil or a reconstructed Jacobian",
    "a characteristic-zero fibre theorem or transported component marking",
    "an equation over QQ rather than the displayed quadratic coefficient field",
]
if HELD_OUT_TRANSPORTED_PRIMES:
    proved_claims.append(
        "literal replay at transported primes excluded from CRT/LLL: "
        + ",".join(str(value) for value in HELD_OUT_TRANSPORTED_PRIMES)
    )
else:
    not_proved_claims.insert(
        0, "replay at a new aligned prime not consumed by the CRT reconstruction"
    )
reproduce_parts = [
    "sage -python elkies-k3/scripts/reconstruct_q80_third_q12_j_map_p19_adic.sage",
    f"--manifest {args.manifest.relative_to(ROOT)}",
    f"--transport {args.transport.relative_to(ROOT)}",
    f"--source {args.source.relative_to(ROOT)}",
    f"--operands {args.operands.relative_to(ROOT)}",
    f"--output {args.output.relative_to(ROOT)}",
    f"--reconstruction-granularity {args.reconstruction_granularity}",
    f"--base-normalization {args.base_normalization}",
    f"--intrinsic-basis {args.intrinsic_basis}",
    f"--c4-pivot {args.c4_pivot}",
]
for heldout_prime in HELD_OUT_TRANSPORTED_PRIMES:
    reproduce_parts.append(f"--holdout-prime {heldout_prime}")
output_base_variable = base_normalization["coordinate"]


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
        "base_coordinate": (
            f"the intrinsic {output_base_variable} coordinate defined in base_normalization"
            if args.base_normalization != "pinned"
            else "the exact common U coordinate pinned by the transported models"
        ),
        "base_normalization": base_normalization,
        "base_PGL2_change_applied": args.base_normalization != "pinned",
        "weierstrass_scaling_dependence": False,
        "section_sign_dependence": False,
        "denominator_normalization": (
            f"monic of degree 20 in {output_base_variable} with the I4 fibre at infinity"
            if args.base_normalization != "pinned"
            else "monic of degree 24 in U"
        ),
    },
    "j_map": {
        "degrees_numerator_denominator": [degree, denominator_degree],
        "base_variable": output_base_variable,
        (
            f"numerator_coefficients_low_to_high_{output_base_variable}"
            if args.base_normalization != "pinned"
            else "numerator_coefficients_low_to_high_U"
        ): encode_polynomial(exact["numerator"]),
        (
            f"denominator_coefficients_low_to_high_{output_base_variable}"
            if args.base_normalization != "pinned"
            else "denominator_coefficients_low_to_high_U"
        ): encode_polynomial(exact["denominator"]),
    },
    "exact_invariant_structure": {
        "numerator_leading_scalar": encode_pair(leading),
        (
            f"monic_c4_factor_degree8_coefficients_low_to_high_{output_base_variable}"
            if args.base_normalization != "pinned"
            else "monic_c4_factor_degree8_coefficients_low_to_high_U"
        ): encode_polynomial(cube_root),
        "identity": "numerator = numerator_leading_scalar * monic_c4_factor^3",
        "denominator_squarefree_factor_degrees_by_multiplicity": {
            str(key): value for key, value in sorted(actual_factor_degrees.items())
        },
        "denominator_squarefree_factors": {
            str(key): encode_polynomial(value)
            for key, value in sorted(squarefree_factors.items())
        },
        "fibre_multiplicity_shape": (
            f"I6 at {output_base_variable}=0 + I4 at {output_base_variable}=infinity + 3 I2 + 8 simple residual discriminant points"
            if args.base_normalization != "pinned"
            else "I6 + I4 + 3 I2 + 8 simple residual discriminant points"
        ),
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
        "proved": proved_claims,
        "not_proved": not_proved_claims,
    },
    "reproduce": " ".join(reproduce_parts),
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
    f"reconstruction_primes={','.join(str(value) for value in RECONSTRUCTION_TRANSPORTED_PRIMES)}|"
    f"holdout_primes={','.join(str(value) for value in HELD_OUT_TRANSPORTED_PRIMES) or 'none'}|"
    f"granularity={args.reconstruction_granularity}|"
    "status=PASS_CANDIDATE_EXACT_THIRD_Q12_J_MAP_RECONSTRUCTION_QQ"
)
