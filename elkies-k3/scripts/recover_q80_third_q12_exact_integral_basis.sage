#!/usr/bin/env sage -python
"""Recover and certify the exact generic cubic integral basis for third q12.

The exact discriminant conductor is C=L*Q.  On C the cubic has a repeated
root rho, recovered without a new factorization from

    rho = (9*a*d-b*c)/(2*(b^2-3*a*c)) mod C.

Here a*x^3+b*x^2+c*x+d is the (not necessarily monic) moving cubic.  The
candidate integral generator is

    e = (x^2+A*x+B)/C,
    A = b/a+rho,  B = -(b/a)*rho-2*rho^2  (mod C).

Every assertion below is in characteristic zero.  The old 19-adic lift is
used only as an independent replay of the final compact coefficients.
"""

import argparse
import hashlib
import json
import multiprocessing
import re
import shlex
import sys
from pathlib import Path

from sage.all import Matrix, NumberField, PolynomialRing, QQ, ZZ, vector


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_OPERANDS = RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json"
DEFAULT_PENCIL = RESULTS / "q80-third-q12-um2-biquadratic-resolved-pencil-qq.json"
DEFAULT_FACTORIZATION = RESULTS / "elkies-k3-q80-third-q12-exact-generic-quartic-factorization-v1.json"
DEFAULT_PADIC_LIFT = ROOT / "artifacts/local/elkies-k3/q80-third-q12-discriminant-factors-p19-adic-precision12288.json"
DEFAULT_OUTPUT = RESULTS / "elkies-k3-q80-third-q12-exact-integral-basis-v1.json"
COEFFICIENT = re.compile(r"^(-?\d+)/(\d+)\*theta\^2 ([+-]) (\d+)/(\d+)$")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path):
    return json.loads(path.read_text())


def record_rational(record):
    return QQ(ZZ(record["numerator"])) / QQ(ZZ(record["denominator"]))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--operands", type=Path, default=DEFAULT_OPERANDS)
parser.add_argument("--pencil", type=Path, default=DEFAULT_PENCIL)
parser.add_argument("--factorization", type=Path, default=DEFAULT_FACTORIZATION)
parser.add_argument("--padic-lift", type=Path, default=DEFAULT_PADIC_LIFT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
parser.add_argument("--workers", type=int, default=4)
parser.add_argument(
    "--direct-xgcd",
    action="store_true",
    help="use the slower direct extended gcd over QQ(delta)(V) instead of exact interpolation",
)
parser.add_argument(
    "--no-padic-replay",
    action="store_true",
    help="skip the independent coefficient replay (exact gates still run)",
)
args = parser.parse_args()
for name in ("operands", "pencil", "factorization", "padic_lift", "output"):
    setattr(args, name, getattr(args, name).resolve())

operands = load_json(args.operands)
pencil = load_json(args.pencil)
factorization = load_json(args.factorization)
padic_lift = load_json(args.padic_lift)
if pencil.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_RESOLVED_PENCIL":
    raise ValueError("the exact moving pencil is not certified")
if factorization.get("status") != "PASS_EXACT_GENERIC_L3_Q2_D_FACTORIZATION":
    raise ValueError("the exact generic conductor factorization is not certified")
if padic_lift.get("status") != "PASS_EXACT_THIRD_Q12_DISCRIMINANT_FACTOR_HENSEL_LIFT_P19":
    raise ValueError("the independent p-adic lift is not certified")
for label, path in (("operands", args.operands), ("pencil", args.pencil)):
    if factorization["inputs"][label]["sha256"] != sha256(path):
        raise ValueError(f"the exact factorization has stale {label} input")

q1 = record_rational(operands["biquadratic_field"]["q1"])
q2 = record_rational(operands["biquadratic_field"]["q2"])
product = q1 * q2
product_numerator_root = ZZ(product.numerator()).isqrt()
if product_numerator_root**2 != product.numerator():
    raise ArithmeticError("q1*q2 numerator is not a square")
delta_square = QQ(product.denominator())
if str(delta_square) != factorization["quadratic_field"]["delta_square"]:
    raise ArithmeticError("descent field changed")
omega_to_delta = QQ(4 * product_numerator_root) / QQ(product.denominator())

delta_polynomial_ring = PolynomialRing(QQ, "delta_polynomial")
delta_polynomial = delta_polynomial_ring.gen()
K = NumberField(delta_polynomial**2 - delta_square, "delta")
delta = K.gen()
V_ring = PolynomialRing(K, "V")
V = V_ring.gen()
F = V_ring.fraction_field()
W_ring = PolynomialRing(F, "W")
W = W_ring.gen()


def exact_field_pair(pair):
    """Decode [rational-record, rational-record] as a+b*delta."""
    return K(record_rational(pair[0])) + K(record_rational(pair[1])) * delta


print("Q80Q12EXACTIB|stage=assemble_cubic", flush=True)
raw = [W_ring.zero() for _ in range(4)]
for v_degree, w_degree, x_degree, encoded in pencil["moving_equation"][
    "terms_T_W_x_coefficient_1_r"
]:
    if len(encoded) != 1:
        raise ArithmeticError("unexpected moving-pencil coefficient list")
    match = COEFFICIENT.fullmatch(encoded[0])
    if match is None:
        raise ArithmeticError("unexpected exact coefficient encoding")
    theta2 = QQ(ZZ(match[1])) / QQ(ZZ(match[2]))
    sign = 1 if match[3] == "+" else -1
    constant = QQ(sign * ZZ(match[4])) / QQ(ZZ(match[5]))
    coefficient = K(constant + theta2 * (q1 + q2))
    coefficient += K(2 * theta2 * product_numerator_root / product.denominator()) * delta
    raw[x_degree] += F(coefficient * V**v_degree) * W**w_degree

d_raw, c_raw, b_raw, a_raw = raw
if a_raw.degree() != 0:
    raise ArithmeticError("cubic leading coefficient is not W-constant")
a = F(a_raw[0])
if not a:
    raise ArithmeticError("cubic leading coefficient vanishes")
b = b_raw
c = c_raw
d = d_raw
b_monic = b / a
c_monic = c / a
d_monic = d / a

print("Q80Q12EXACTIB|stage=assemble_conductor", flush=True)
qhat = W_ring.zero()
qhat_records = factorization["generic_quartic_first_jet"][
    "Q_numerator_coefficients_low_to_high_W_then_V_1_delta"
]
for w_degree, v_records in enumerate(qhat_records):
    qhat += F(
        sum((exact_field_pair(pair) * V**v_degree for v_degree, pair in enumerate(v_records)), K.zero())
    ) * W**w_degree
linear_constant = record_rational(factorization["linear_factor_reconstruction"]["constant"])
L = W + F(K(linear_constant))
C_raw = L * qhat
C = C_raw.monic()
if C.degree() != 5:
    raise ArithmeticError("generic conductor is not degree five in W")

def direct_integral_basis():
    print("Q80Q12EXACTIB|stage=invert_repeated_root_denominator", flush=True)
    root_numerator = (9 * a * d - b * c) % C
    root_denominator = (2 * (b**2 - 3 * a * c)) % C
    gcd_value, inverse_value, unused = root_denominator.xgcd(C)
    if gcd_value.degree() != 0:
        raise ArithmeticError("repeated-root denominator is not invertible modulo LQ")
    inverse_value /= gcd_value[0]
    root = (root_numerator * inverse_value) % C
    return root, (b_monic + root) % C, (-b_monic * root - 2 * root**2) % C


def specialize_w(polynomial, value, target_ring):
    return target_ring([K(F(coefficient)(value)) for coefficient in polynomial.list()])


PAIR_ZERO = (QQ.zero(), QQ.zero())
PAIR_ONE = (QQ.one(), QQ.zero())


def pair(value):
    coefficients = list(K(value))
    coefficients += [QQ.zero()] * (2 - len(coefficients))
    return tuple(coefficients[:2])


def pair_add(left, right):
    return left[0] + right[0], left[1] + right[1]


def pair_neg(value):
    return -value[0], -value[1]


def pair_mul(left, right):
    return (
        left[0] * right[0] + left[1] * right[1] * delta_square,
        left[0] * right[1] + left[1] * right[0],
    )


def pair_inverse(value):
    norm = value[0] ** 2 - value[1] ** 2 * delta_square
    if not norm:
        raise ZeroDivisionError("zero quadratic-field norm")
    return value[0] / norm, -value[1] / norm


def pair_scale(value, scalar):
    return value[0] * scalar, value[1] * scalar


def small_trim(value):
    value = list(value)
    while value and value[-1] == PAIR_ZERO:
        value.pop()
    return value or [PAIR_ZERO]


def small_add(left, right):
    result = [PAIR_ZERO] * max(len(left), len(right))
    for index in range(len(result)):
        result[index] = pair_add(
            left[index] if index < len(left) else PAIR_ZERO,
            right[index] if index < len(right) else PAIR_ZERO,
        )
    return small_trim(result)


def small_neg(value):
    return small_trim([pair_neg(coefficient) for coefficient in value])


def small_scale(value, scalar):
    return small_trim([pair_scale(coefficient, scalar) for coefficient in value])


def small_pair_scale(value, scalar):
    return small_trim([pair_mul(coefficient, scalar) for coefficient in value])


def small_mul(left, right):
    result = [PAIR_ZERO] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result[i + j] = pair_add(
                result[i + j], pair_mul(left_coefficient, right_coefficient)
            )
    return small_trim(result)


def small_divmod(dividend, divisor):
    dividend = list(small_trim(dividend))
    divisor = small_trim(divisor)
    if divisor == [PAIR_ZERO]:
        raise ZeroDivisionError("zero polynomial divisor")
    if len(dividend) < len(divisor):
        return [PAIR_ZERO], dividend
    quotient = [PAIR_ZERO] * (len(dividend) - len(divisor) + 1)
    inverse_lead = pair_inverse(divisor[-1])
    while dividend != [PAIR_ZERO] and len(dividend) >= len(divisor):
        coefficient = pair_mul(dividend[-1], inverse_lead)
        shift = len(dividend) - len(divisor)
        quotient[shift] = pair_add(quotient[shift], coefficient)
        for index, divisor_coefficient in enumerate(divisor):
            position = shift + index
            dividend[position] = pair_add(
                dividend[position], pair_neg(pair_mul(coefficient, divisor_coefficient))
            )
        dividend = small_trim(dividend)
    return small_trim(quotient), small_trim(dividend)


def small_mod(value, modulus):
    return small_divmod(value, modulus)[1]


def small_inverse_mod(value, modulus):
    old_remainder, remainder = small_trim(modulus), small_mod(value, modulus)
    old_coefficient, coefficient = [PAIR_ZERO], [PAIR_ONE]
    if remainder != [PAIR_ZERO]:
        scale = pair_inverse(remainder[-1])
        remainder = small_pair_scale(remainder, scale)
        coefficient = small_pair_scale(coefficient, scale)
    while remainder != [PAIR_ZERO]:
        quotient, next_remainder = small_divmod(old_remainder, remainder)
        next_coefficient = small_add(
            old_coefficient, small_neg(small_mul(quotient, coefficient))
        )
        if next_remainder != [PAIR_ZERO]:
            scale = pair_inverse(next_remainder[-1])
            next_remainder = small_pair_scale(next_remainder, scale)
            next_coefficient = small_pair_scale(next_coefficient, scale)
        old_remainder, remainder = remainder, next_remainder
        old_coefficient, coefficient = coefficient, next_coefficient
    if len(old_remainder) != 1 or old_remainder[0] == PAIR_ZERO:
        raise ZeroDivisionError("polynomial is not invertible modulo conductor")
    answer = small_pair_scale(old_coefficient, pair_inverse(old_remainder[0]))
    return small_mod(answer, modulus)


def small_from_w(polynomial, value):
    return small_trim([pair(F(coefficient)(value)) for coefficient in polynomial.list()])


def small_to_k(value):
    return [K(coefficient[0]) + K(coefficient[1]) * delta for coefficient in value]


def specialized_integral_basis(value):
    """Compute one exact fibre using bare pairs in QQ(delta)[W]."""
    aa = pair(a(value))
    if aa == PAIR_ZERO:
        raise ZeroDivisionError("cubic leading coefficient vanishes at sample")
    bb = small_from_w(b, value)
    cc = small_from_w(c, value)
    dd = small_from_w(d, value)
    conductor = small_from_w(C_raw, value)
    if len(conductor) != 6:
        raise ZeroDivisionError("conductor degree drops at sample")
    conductor = small_pair_scale(conductor, pair_inverse(conductor[-1]))
    numerator = small_mod(
        small_add(
            small_pair_scale(dd, pair_scale(aa, 9)),
            small_neg(small_mul(bb, cc)),
        ),
        conductor,
    )
    denominator = small_mod(
        small_scale(
            small_add(
                small_mul(bb, bb),
                small_neg(small_pair_scale(cc, pair_scale(aa, 3))),
            ),
            2,
        ),
        conductor,
    )
    root = small_mod(small_mul(numerator, small_inverse_mod(denominator, conductor)), conductor)
    bb_monic = small_pair_scale(bb, pair_inverse(aa))
    integral_a = small_mod(small_add(bb_monic, root), conductor)
    integral_b = small_mod(
        small_neg(
            small_add(small_mul(bb_monic, root), small_scale(small_mul(root, root), 2))
        ),
        conductor,
    )
    integral_a += [PAIR_ZERO] * (5 - len(integral_a))
    integral_b += [PAIR_ZERO] * (5 - len(integral_b))
    return small_to_k(integral_a), small_to_k(integral_b)


def specialized_worker(candidate):
    try:
        integral_a, integral_b = specialized_integral_basis(K(candidate))
        return (
            candidate,
            [[pair(value) for value in integral_a], [pair(value) for value in integral_b]],
            None,
        )
    except (ArithmeticError, ZeroDivisionError) as error:
        return candidate, None, str(error)


def interpolate_function(samples, values, numerator_degree, denominator_degree):
    """Interpolate P(V)/Q(V), with Q monic of the requested degree."""
    count = numerator_degree + denominator_degree + 1
    samples = samples[:count]
    values = values[:count]
    rows = []
    right = []
    for sample, value in zip(samples, values):
        sample = K(sample)
        rows.append(
            [sample**index for index in range(numerator_degree + 1)]
            + [-value * sample**index for index in range(denominator_degree)]
        )
        right.append(value * sample**denominator_degree)
    solution = Matrix(K, rows).solve_right(vector(K, right))
    numerator = V_ring(solution[: numerator_degree + 1])
    denominator = V_ring(list(solution[numerator_degree + 1 :]) + [K.one()])
    return F(numerator) / F(denominator)


def interpolated_integral_basis():
    degree_profiles = {
        name: [
            list(record["degrees_numerator_denominator"])
            for record in padic_lift["integral_basis_candidate"][name]["coefficients_low_to_high_W"]
        ]
        for name in ("A", "B")
    }
    maximum_samples = max(sum(profile) + 1 for profiles in degree_profiles.values() for profile in profiles)
    sample_results = {}
    candidates = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6]
    print(
        f"Q80Q12EXACTIB|stage=specialize_for_interpolation|needed={maximum_samples}|"
        f"workers={args.workers}",
        flush=True,
    )
    pending = candidates[:maximum_samples]
    reserve = candidates[maximum_samples:]
    context = multiprocessing.get_context("fork")
    while pending and len(sample_results) < maximum_samples:
        with context.Pool(processes=min(args.workers, len(pending))) as pool:
            for candidate, result, error in pool.imap_unordered(specialized_worker, pending):
                if result is not None:
                    sample_results[candidate] = result
                    print(
                        f"Q80Q12EXACTIB|stage=specialized|V={candidate}|"
                        f"count={len(sample_results)}",
                        flush=True,
                    )
                else:
                    print(f"Q80Q12EXACTIB|stage=singular_sample|V={candidate}", flush=True)
        missing = maximum_samples - len(sample_results)
        pending, reserve = reserve[:missing], reserve[missing:]
    if len(sample_results) != maximum_samples:
        raise ArithmeticError("insufficient nonsingular exact interpolation fibres")
    selected = [candidate for candidate in candidates if candidate in sample_results][:maximum_samples]
    samples = [K(candidate) for candidate in selected]
    a_values = [[] for _ in range(5)]
    b_values = [[] for _ in range(5)]
    for candidate in selected:
        result_a, result_b = sample_results[candidate]
        for index in range(5):
            a_values[index].append(K(result_a[index][0]) + K(result_a[index][1]) * delta)
            b_values[index].append(K(result_b[index][0]) + K(result_b[index][1]) * delta)
    print("Q80Q12EXACTIB|stage=interpolate", flush=True)
    recovered_a = W_ring(
        [
            interpolate_function(samples, a_values[index], *degree_profiles["A"][index])
            for index in range(5)
        ]
    )
    recovered_b = W_ring(
        [
            interpolate_function(samples, b_values[index], *degree_profiles["B"][index])
            for index in range(5)
        ]
    )
    root = (recovered_a - b_monic) % C
    if recovered_a != (b_monic + root) % C:
        raise ArithmeticError("interpolated A does not equal b+rho modulo C")
    if recovered_b != (-b_monic * root - 2 * root**2) % C:
        raise ArithmeticError("interpolated B does not equal -b*rho-2*rho^2 modulo C")
    return root, recovered_a, recovered_b, [int(value) for value in samples]


if args.direct_xgcd:
    rho, A, B = direct_integral_basis()
    interpolation_samples = None
    recovery_method = "direct extended gcd over QQ(delta)(V)[W]"
else:
    rho, A, B, interpolation_samples = interpolated_integral_basis()
    recovery_method = "exact specialization/interpolation with full generic replay"
if rho.degree() > 4:
    raise ArithmeticError("repeated root has degree above four modulo LQ")

print("Q80Q12EXACTIB|stage=verify_repeated_root", flush=True)
derivative_remainder = (3 * a * rho**2 + 2 * b * rho + c) % C
cubic_remainder = (a * rho**3 + b * rho**2 + c * rho + d) % C
if derivative_remainder or cubic_remainder:
    raise ArithmeticError("closed-form rho is not a repeated cubic root modulo LQ")

if A.degree() > 4 or B.degree() > 4:
    raise ArithmeticError("integral-basis numerator has degree above four")
if B % L:
    raise ArithmeticError("B is not divisible by the linear conductor factor")


def matrix_characteristic_remainders():
    """Return the three exact divisibility remainders for n=x^2+A*x+B."""
    one = W_ring.one()
    zero = W_ring.zero()
    m = [
        [B, -d_monic, d_monic * (b_monic - A)],
        [A, B - c_monic, c_monic * (b_monic - A) - d_monic],
        [one, A - b_monic, B + b_monic**2 - c_monic - A * b_monic],
    ]
    trace = m[0][0] + m[1][1] + m[2][2]
    trace_m2 = zero
    for row in range(3):
        for column in range(3):
            trace_m2 += m[row][column] * m[column][row]
    second = (trace**2 - trace_m2) / 2
    determinant = zero
    for sign, indices in (
        (1, (0, 1, 2)),
        (1, (1, 2, 0)),
        (1, (2, 0, 1)),
        (-1, (2, 1, 0)),
        (-1, (1, 0, 2)),
        (-1, (0, 2, 1)),
    ):
        determinant += sign * m[0][indices[0]] * m[1][indices[1]] * m[2][indices[2]]
    return trace % C, second % (C**2), determinant % (C**3)


print("Q80Q12EXACTIB|stage=verify_integrality", flush=True)
trace_remainder, second_remainder, determinant_remainder = matrix_characteristic_remainders()
if trace_remainder or second_remainder or determinant_remainder:
    raise ArithmeticError("integral-basis characteristic divisibility fails")


def normalized_fraction(value):
    value = F(value)
    numerator = V_ring(value.numerator())
    denominator = V_ring(value.denominator())
    leading = denominator.leading_coefficient()
    return numerator / leading, denominator / leading


def field_coordinates(value):
    values = list(K(value))
    values += [QQ.zero()] * (2 - len(values))
    return values[:2]


def rational_output(value):
    value = QQ(value)
    return {"numerator": str(value.numerator()), "denominator": str(value.denominator())}


def field_output(value):
    return [rational_output(coordinate) for coordinate in field_coordinates(value)]


def function_output(value):
    numerator, denominator = normalized_fraction(value)
    return {
        "numerator_coefficients_low_to_high_V_1_delta": [field_output(x) for x in numerator.list()],
        "denominator_coefficients_low_to_high_V_1_delta": [field_output(x) for x in denominator.list()],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
    }


def w_polynomial_output(value):
    value = W_ring(value)
    return {
        "degree_W": int(value.degree()),
        "coefficients_low_to_high_W": [function_output(x) for x in value.list()],
    }


def rational_mod(value, modulus):
    value = QQ(value)
    return int(ZZ(value.numerator()) * ZZ(value.denominator()).inverse_mod(modulus) % modulus)


def replay_one_function(exact, record, modulus):
    numerator, denominator = normalized_fraction(exact)
    exact_parts = (numerator, denominator)
    record_keys = (
        "numerator_coefficients_low_to_high_U_1_omega",
        "denominator_coefficients_low_to_high_U_1_omega",
    )
    expected_degrees = [int(numerator.degree()), int(denominator.degree())]
    if expected_degrees != list(record["degrees_numerator_denominator"]):
        raise ArithmeticError("exact/p-adic rational-function degrees disagree")
    checked = 0
    for polynomial, key in zip(exact_parts, record_keys):
        saved = record[key]
        coefficients = polynomial.list()
        if len(coefficients) != len(saved):
            raise ArithmeticError("exact/p-adic coefficient supports disagree")
        for exact_coefficient, saved_pair in zip(coefficients, saved):
            constant, delta_coefficient = field_coordinates(exact_coefficient)
            omega_pair = (constant, delta_coefficient / omega_to_delta)
            for exact_coordinate, saved_coordinate in zip(omega_pair, saved_pair):
                if rational_mod(exact_coordinate, modulus) != ZZ(saved_coordinate) % modulus:
                    raise ArithmeticError("exact integral basis fails the full p-adic replay")
                checked += 1
    return checked


padic_coordinates_replayed = 0
if not args.no_padic_replay:
    print("Q80Q12EXACTIB|stage=replay_padic", flush=True)
    modulus = ZZ(padic_lift["specialization"]["modulus"])
    for name, exact in (("A", A), ("B", B)):
        records = padic_lift["integral_basis_candidate"][name]["coefficients_low_to_high_W"]
        if len(records) != exact.degree() + 1:
            raise ArithmeticError("exact/p-adic W supports disagree")
        for coefficient, record in zip(exact.list(), records):
            padic_coordinates_replayed += replay_one_function(coefficient, record, modulus)

print("Q80Q12EXACTIB|stage=serialize", flush=True)
payload = {
    "schema": "elkies-k3.q80-third-q12-exact-integral-basis.v1",
    "status": "PASS_EXACT_GENERIC_REPEATED_ROOT_AND_INTEGRAL_BASIS",
    "quadratic_field": {
        "generator": "delta",
        "delta_square": str(delta_square),
        "omega_in_delta_basis": {"coefficient": rational_output(omega_to_delta)},
    },
    "conductor": {
        "formula": "C=(W+r)*Qhat/H, with the scalar H omitted in the polynomial modulus",
        "degree_W": 5,
        "L_constant": rational_output(linear_constant),
    },
    "repeated_root": {
        "formula": "rho=(9*a*d-b*c)/(2*(b^2-3*a*c)) mod C",
        "recovery_method": recovery_method,
        "interpolation_samples": interpolation_samples,
        "rho": w_polynomial_output(rho),
        "cubic_remainder_zero": True,
        "derivative_remainder_zero": True,
        "denominator_invertible_mod_C": True,
    },
    "integral_basis": {
        "basis": ["1", "x", "e"],
        "e_formula": "e=(x^2+A*x+B)/(L*Q)",
        "A": w_polynomial_output(A),
        "B": w_polynomial_output(B),
        "degrees_W_A_B": [int(A.degree()), int(B.degree())],
        "B_divisible_by_L": True,
        "characteristic_coefficient_divisibility": {
            "trace_by_LQ": True,
            "second_symmetric_by_LQ_squared": True,
            "determinant_by_LQ_cubed": True,
        },
    },
    "independent_replay": {
        "prime": 19,
        "digits": int(padic_lift["specialization"]["digits"]),
        "coordinates_replayed": padic_coordinates_replayed,
        "skipped": bool(args.no_padic_replay),
    },
    "inputs": {
        label: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for label, path in (
            ("operands", args.operands),
            ("pencil", args.pencil),
            ("factorization", args.factorization),
            ("padic_lift", args.padic_lift),
        )
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "the exact generic repeated root of the moving cubic modulo the full degree-five conductor LQ",
            "the exact characteristic-zero integral basis 1,x,(x^2+A*x+B)/(LQ)",
            "trace, second-symmetric, and determinant divisibility by LQ, (LQ)^2, and (LQ)^3",
            "agreement with the independent 12,288-digit 19-adic lift when replay is enabled",
        ],
        "not_proved": [
            "Riemann--Roch pencil generators in the normalized order",
            "the characteristic-zero Jacobian, Weierstrass equation, sections, or maps",
        ],
    },
    "reproduce": shlex.join(
        [
            str(Path(sys.executable)),
            str(Path(__file__).resolve().relative_to(ROOT)),
            "--output",
            str(args.output.relative_to(ROOT)),
        ]
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"exact integral-basis artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    "Q80Q12EXACTIB|degree_C=5|degree_A=4|degree_B=4|"
    f"padic_coordinates={padic_coordinates_replayed}|"
    "status=PASS_EXACT_GENERIC_REPEATED_ROOT_AND_INTEGRAL_BASIS",
    flush=True,
)
