#!/usr/bin/env python3
"""Verify the lifted cubic integral basis over (Z/19^n)[omega](V)."""

import argparse
import hashlib
import json
import shlex
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_SOURCE = RESULTS / "q80-third-q12-exact-pencil-p19-adic-precision64.json"
DEFAULT_LIFT = RESULTS / "q80-third-q12-discriminant-factors-p19-adic-precision5.json"
DEFAULT_OUTPUT = RESULTS / "q80-third-q12-integral-basis-mod19-power.json"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument("--lift", type=Path, default=DEFAULT_LIFT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument(
    "--verification-digits",
    type=int,
    help=(
        "reduce the certified source and lift to this many p-adic digits for the "
        "generic divisibility proof; defaults to the full lift precision"
    ),
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
args.source = args.source.resolve()
args.lift = args.lift.resolve()
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


source = json.loads(args.source.read_text())
lift = json.loads(args.lift.read_text())
if source.get("status") != "PASS_EXACT_THIRD_Q12_PENCIL_REDUCTION_MOD_19_POWER":
    raise ValueError("exact p-adic source is not certified")
if lift.get("status") != "PASS_EXACT_THIRD_Q12_DISCRIMINANT_FACTOR_HENSEL_LIFT_P19":
    raise ValueError("p-adic discriminant/repeated-root lift is not certified")
prime = 19
lift_digits = int(lift["specialization"]["digits"])
digits = lift_digits if args.verification_digits is None else int(args.verification_digits)
if digits < 1 or digits > lift_digits:
    raise ValueError("verification digits must lie between one and the lift precision")
modulus = prime**digits
omega_square = int(source["quadratic_field"]["omega_square_modulus"]) % modulus
inverse_two = pow(2, -1, modulus)

# Constants are pairs a+b*omega in the unramified quadratic ring.
ZERO_C = (0, 0)
ONE_C = (1, 0)


def c_add(left, right):
    return ((left[0] + right[0]) % modulus, (left[1] + right[1]) % modulus)


def c_neg(value):
    return ((-value[0]) % modulus, (-value[1]) % modulus)


def c_sub(left, right):
    return c_add(left, c_neg(right))


def c_mul(left, right):
    return (
        (left[0] * right[0] + omega_square * left[1] * right[1]) % modulus,
        (left[0] * right[1] + left[1] * right[0]) % modulus,
    )


def c_scale(value, scalar):
    return (value[0] * scalar % modulus, value[1] * scalar % modulus)


def c_is_zero(value):
    return value == ZERO_C


def coordinates(value):
    return int(value[0]) % modulus, int(value[1]) % modulus


def c_inv(value):
    norm = (value[0] * value[0] - omega_square * value[1] * value[1]) % modulus
    if norm % prime == 0:
        raise ZeroDivisionError("non-unit quadratic coefficient")
    inverse_norm = pow(norm, -1, modulus)
    return value[0] * inverse_norm % modulus, -value[1] * inverse_norm % modulus


# Polynomials in the exact base U, represented low to high.
def u_trim(value):
    value = list(value)
    while value and c_is_zero(value[-1]):
        value.pop()
    return tuple(value)


ZERO_U = ()
ONE_U = (ONE_C,)


def u_add(left, right):
    result = [ZERO_C] * max(len(left), len(right))
    for index in range(len(result)):
        result[index] = c_add(
            left[index] if index < len(left) else ZERO_C,
            right[index] if index < len(right) else ZERO_C,
        )
    return u_trim(result)


def u_neg(value):
    return tuple(c_neg(coefficient) for coefficient in value)


def u_sub(left, right):
    return u_add(left, u_neg(right))


def u_mul(left, right):
    if not left or not right:
        return ZERO_U
    result = [ZERO_C] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result[i + j] = c_add(
                result[i + j], c_mul(left_coefficient, right_coefficient)
            )
    return u_trim(result)


def u_scale(value, scalar):
    return u_trim(c_scale(coefficient, scalar) for coefficient in value)


def u_is_primitive(value):
    return any(coefficient[0] % prime or coefficient[1] % prime for coefficient in value)


def u_pow(value, exponent):
    result = ONE_U
    base = value
    while exponent:
        if exponent & 1:
            result = u_mul(result, base)
        base = u_mul(base, base)
        exponent >>= 1
    return result


def u_div_monic(dividend, divisor):
    dividend = list(u_trim(dividend))
    divisor = u_trim(divisor)
    if not divisor or divisor[-1] != ONE_C:
        raise ArithmeticError("U divisor is not literally monic")
    if len(dividend) < len(divisor):
        return ZERO_U, tuple(dividend)
    quotient = [ZERO_C] * (len(dividend) - len(divisor) + 1)
    while len(dividend) >= len(divisor):
        coefficient = dividend[-1]
        shift = len(dividend) - len(divisor)
        quotient[shift] = coefficient
        for index in range(len(divisor)):
            dividend[shift + index] = c_sub(
                dividend[shift + index], c_mul(coefficient, divisor[index])
            )
        dividend = list(u_trim(dividend))
    return u_trim(quotient), tuple(dividend)


def record_u_polynomial(values):
    return u_trim(coordinates(value) for value in values)


# Put every input coefficient over one fixed square-free-enough common
# denominator H(U).  A rational function is then (numerator, H-exponent).
# This prevents the expression swell caused by multiplying unrelated copies
# of identical denominators during W-polynomial reduction.
input_denominators = []
for container, names in (
    (lift["factorization"], ("L", "Q", "D")),
    (lift["integral_basis_candidate"], ("A", "B")),
):
    for name in names:
        for coefficient in container[name]["coefficients_low_to_high_W"]:
            denominator = record_u_polynomial(
                coefficient["denominator_coefficients_low_to_high_U_1_omega"]
            )
            if denominator != ONE_U and denominator not in input_denominators:
                if denominator[-1] != ONE_C:
                    raise ArithmeticError("compressed input denominator is not monic")
                input_denominators.append(denominator)

leading_u_input = {}
for u_degree, w_degree, x_degree, coefficient in source["pencil"][
    "terms_V_W_old_x_coefficient_1_omega"
]:
    if x_degree == 3 and w_degree == 0:
        leading_u_input[u_degree] = coordinates(coefficient)
leading_u_exact = u_trim(
    leading_u_input.get(index, ZERO_C)
    for index in range(max(leading_u_input) + 1)
)
leading_scalar = leading_u_exact[-1]
leading_scalar_inverse = c_inv(leading_scalar)
leading_monic = u_trim(
    c_mul(coefficient, leading_scalar_inverse) for coefficient in leading_u_exact
)
if leading_monic != ONE_U and leading_monic not in input_denominators:
    input_denominators.append(leading_monic)

GLOBAL_H = ONE_U
for denominator in input_denominators:
    GLOBAL_H = u_mul(GLOBAL_H, denominator)


def rat_atomic(numerator=ZERO_U, denominator=ONE_U):
    numerator = u_trim(numerator)
    denominator = u_trim(denominator)
    if not denominator or not u_is_primitive(denominator):
        raise ZeroDivisionError("non-primitive U denominator")
    if not numerator:
        return ZERO_U, 0
    if denominator == ONE_U:
        return numerator, 0
    quotient, remainder = u_div_monic(GLOBAL_H, denominator)
    if remainder:
        raise ArithmeticError("input denominator does not divide global H")
    return u_mul(numerator, quotient), 1


ZERO_R = (ZERO_U, 0)
ONE_R = (ONE_U, 0)


def r_add(left, right):
    if r_is_zero(left):
        return right
    if r_is_zero(right):
        return left
    exponent = max(left[1], right[1])
    left_numerator = u_mul(left[0], u_pow(GLOBAL_H, exponent - left[1]))
    right_numerator = u_mul(right[0], u_pow(GLOBAL_H, exponent - right[1]))
    return u_add(left_numerator, right_numerator), exponent


def r_neg(value):
    return u_neg(value[0]), value[1]


def r_sub(left, right):
    return r_add(left, r_neg(right))


def r_mul(left, right):
    if r_is_zero(left) or r_is_zero(right):
        return ZERO_R
    if left == ONE_R:
        return right
    if right == ONE_R:
        return left
    return u_mul(left[0], right[0]), left[1] + right[1]


def r_scale(value, scalar):
    return u_scale(value[0], scalar), value[1]


def r_is_zero(value):
    return not value[0]


# Polynomials in old W over the rational-function ring.
def w_trim(value):
    value = list(value)
    while value and r_is_zero(value[-1]):
        value.pop()
    return tuple(value)


ZERO_W = ()
ONE_W = (ONE_R,)


def w_add(left, right):
    result = [ZERO_R] * max(len(left), len(right))
    for index in range(len(result)):
        result[index] = r_add(
            left[index] if index < len(left) else ZERO_R,
            right[index] if index < len(right) else ZERO_R,
        )
    return w_trim(result)


def w_neg(value):
    return tuple(r_neg(coefficient) for coefficient in value)


def w_sub(left, right):
    return w_add(left, w_neg(right))


def w_mul(left, right):
    if not left or not right:
        return ZERO_W
    result = [ZERO_R] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result[i + j] = r_add(
                result[i + j], r_mul(left_coefficient, right_coefficient)
            )
    return w_trim(result)


def w_scale(value, scalar):
    return w_trim(r_scale(coefficient, scalar) for coefficient in value)


def w_shift(value, shift):
    return (ZERO_R,) * shift + tuple(value) if value else ZERO_W


def w_monic_remainder(value, divisor):
    value = list(w_trim(value))
    divisor = w_trim(divisor)
    if not divisor or divisor[-1] != ONE_R:
        raise ArithmeticError("W divisor is not literally monic")
    divisor_degree = len(divisor) - 1
    while len(value) - 1 >= divisor_degree:
        coefficient = value[-1]
        shift = len(value) - len(divisor)
        for index in range(len(divisor)):
            value[shift + index] = r_sub(
                value[shift + index], r_mul(coefficient, divisor[index])
            )
        value = list(w_trim(value))
    return tuple(value)


def w_pow(value, exponent):
    result = ONE_W
    base = value
    while exponent:
        if exponent & 1:
            result = w_mul(result, base)
        base = w_mul(base, base)
        exponent >>= 1
    return result


def w_mul_mod(left, right, divisor):
    return w_monic_remainder(w_mul(left, right), divisor)


def rational_record(record):
    numerator = tuple(
        coordinates(value)
        for value in record["numerator_coefficients_low_to_high_U_1_omega"]
    )
    denominator = tuple(
        coordinates(value)
        for value in record["denominator_coefficients_low_to_high_U_1_omega"]
    )
    return rat_atomic(numerator, denominator)


def factor_record(record):
    return w_trim(rational_record(value) for value in record["coefficients_low_to_high_W"])


L = factor_record(lift["factorization"]["L"])
Q = factor_record(lift["factorization"]["Q"])
A = factor_record(lift["integral_basis_candidate"]["A"])
B = factor_record(lift["integral_basis_candidate"]["B"])
conductor = w_mul(L, Q)
if [len(value) - 1 for value in (L, Q, A, B)] != [1, 4, 4, 4]:
    raise ArithmeticError("lifted conductor/integral-basis degrees changed")
if w_monic_remainder(B, L):
    raise ArithmeticError("B is not divisible by L modulo 19^digits")

# Build the monic cubic coefficients b,c,d from the exact source.
raw = [[{} for unused in range(10)] for unused in range(4)]
for u_degree, w_degree, x_degree, coefficient in source["pencil"][
    "terms_V_W_old_x_coefficient_1_omega"
]:
    raw[x_degree][w_degree][u_degree] = coordinates(coefficient)


def raw_w_polynomial(x_degree):
    values = []
    for by_u in raw[x_degree]:
        if not by_u:
            values.append(ZERO_R)
            continue
        maximum = max(by_u)
        polynomial = tuple(by_u.get(index, ZERO_C) for index in range(maximum + 1))
        values.append(rat_atomic(polynomial))
    return w_trim(values)


raw_coefficients = [raw_w_polynomial(index) for index in range(4)]
if len(raw_coefficients[3]) != 1:
    raise ArithmeticError("cubic leading coefficient is not W-constant")
def divide_by_leading(value):
    inverse_leading = rat_atomic((leading_scalar_inverse,), leading_monic)
    return w_trim(r_mul(coefficient, inverse_leading) for coefficient in value)


d, c, b = [divide_by_leading(raw_coefficients[index]) for index in range(3)]

# Multiplication matrix of n=z^2+A*z+B on 1,z,z^2.
m = [
    [B, w_neg(d), w_mul(d, w_sub(b, A))],
    [A, w_sub(B, c), w_sub(w_mul(c, w_sub(b, A)), d)],
    [ONE_W, w_sub(A, b), w_sub(w_add(B, w_mul(b, b)), w_add(c, w_mul(A, b)))],
]
q1 = conductor
q2 = w_pow(conductor, 2)
q3 = w_pow(conductor, 3)
trace = w_add(w_add(m[0][0], m[1][1]), m[2][2])
trace_remainder = w_monic_remainder(trace, q1)
trace_mod_q2 = w_monic_remainder(trace, q2)
trace_square = w_mul_mod(trace_mod_q2, trace_mod_q2, q2)
trace_m2 = ZERO_W
for row in range(3):
    for column in range(3):
        trace_m2 = w_add(trace_m2, w_mul_mod(m[row][column], m[column][row], q2))
trace_m2 = w_monic_remainder(trace_m2, q2)
second_remainder = w_monic_remainder(
    w_scale(w_sub(trace_square, trace_m2), inverse_two), q2
)


def triple(left, middle, right):
    return w_mul_mod(w_mul_mod(left, middle, q3), right, q3)


determinant_remainder = ZERO_W
for sign, indices in (
    (1, (0, 1, 2)),
    (1, (1, 2, 0)),
    (1, (2, 0, 1)),
    (-1, (2, 1, 0)),
    (-1, (1, 0, 2)),
    (-1, (0, 2, 1)),
):
    term = triple(m[0][indices[0]], m[1][indices[1]], m[2][indices[2]])
    determinant_remainder = w_add(
        determinant_remainder, term if sign == 1 else w_neg(term)
    )
determinant_remainder = w_monic_remainder(determinant_remainder, q3)

remainders = {
    "trace_mod_LQ": trace_remainder,
    "second_symmetric_mod_LQ_squared": second_remainder,
    "determinant_mod_LQ_cubed": determinant_remainder,
}
if any(value for value in remainders.values()):
    raise ArithmeticError(
        "candidate integral-basis characteristic coefficient is not divisible"
    )

output = {
    "schema": "elkies-k3.q80-third-q12-integral-basis-mod19-power.v1",
    "status": "PASS_EXACT_THIRD_Q12_GENERIC_INTEGRAL_BASIS_MOD19_POWER",
    "specialization": {
        "u": "-2",
        "prime": prime,
        "digits": digits,
        "lift_digits": lift_digits,
        "modulus": modulus,
    },
    "integral_basis": {
        "basis": ["1", "z", "e"],
        "e_formula": "(z^2+A*z+B)/(L*Q)",
        "degrees_W_L_Q_A_B": [1, 4, 4, 4, 4],
        "B_divisible_by_L": True,
        "characteristic_coefficient_divisibility": {
            "trace_by_LQ": True,
            "second_symmetric_by_LQ_squared": True,
            "determinant_by_LQ_cubed": True,
        },
        "generic_in_U": True,
    },
    "arithmetic": {
        "constant_ring": "(Z/19^digits)[omega]/(omega^2-D)",
        "base": "fraction field localization in U at primitive denominators",
        "old_base": "monic polynomial arithmetic in W",
        "rational_function_representation": "numerator over a power of one fixed global H(U)",
        "global_H_degree": len(GLOBAL_H) - 1,
        "global_H_distinct_input_factors": len(input_denominators),
        "all_denominators_primitive_mod_19": True,
        "zero_remainders_checked_coefficientwise": True,
    },
    "inputs": {
        "source": {"path": str(args.source.relative_to(ROOT)), "sha256": sha256(args.source)},
        "lift": {"path": str(args.lift.relative_to(ROOT)), "sha256": sha256(args.lift)},
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            f"generic characteristic-coefficient divisibility for the candidate basis through {digits} p-adic digits",
            f"integrality of e modulo 19^{digits} over the localized generic coefficient ring",
            "the exact degree and conductor shapes required by p-adic Riemann--Roch",
        ],
        "not_proved": [
            "characteristic-zero integrality or rational reconstruction of the basis",
            "p-adic Riemann--Roch generators, Jacobian, or maps",
        ],
    },
    "reproduce": shlex.join(
        [
            "python3",
            "elkies-k3/scripts/verify_q80_third_q12_integral_basis_mod19_power.py",
            "--source",
            str(args.source.relative_to(ROOT)),
            "--lift",
            str(args.lift.relative_to(ROOT)),
            "--output",
            str(args.output.relative_to(ROOT)),
        ]
        + (
            ["--verification-digits", str(digits)]
            if args.verification_digits is not None
            else []
        )
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"integral-basis verifier artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    f"Q80THIRDQ12INTEGRALBASIS|p=19|digits={digits}|degrees=1,4,4,4,4|"
    "trace=0|second=0|det=0|status=PASS_EXACT_THIRD_Q12_GENERIC_INTEGRAL_BASIS_MOD19_POWER"
)
