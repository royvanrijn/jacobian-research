#!/usr/bin/env sage -python
"""Factor one exact base specialization of the third-q12 cubic discriminant.

status: ACTIVE_COMPILER
claim: feasibility probe for direct exact L^3*Q^2*D conductor recovery
inputs: exact connected pencil and certified biquadratic closure operands
outputs: optional q80-third-q12 exact discriminant specialization artifact

This works in the exact quadratic descent field and specializes only the new
base variable V.  A successful factorization is an exact fibrewise compiler
probe; it is not the generic factorization over QQ(omega)(V).
"""

import argparse
import hashlib
import json
from math import comb
import os
from pathlib import Path
import pickle
import re
import sys
import time

from gmpy2 import gcd, isqrt, lcm, mpq, mpz
from sage.all import GF, NumberField, PolynomialRing, QQ, ZZ, inverse_mod


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_OPERANDS = (
    RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json"
)
DEFAULT_PENCIL = RESULTS / "q80-third-q12-um2-biquadratic-resolved-pencil-qq.json"
DEFAULT_FACTOR_LIFT = (
    ROOT
    / "artifacts/local/elkies-k3/"
    / "q80-third-q12-discriminant-factors-p19-adic-precision12288.json"
)
DEFAULT_SUBRESULTANT_CHECKPOINT = (
    ROOT
    / "artifacts/local/elkies-k3/"
    / "q80-third-q12-exact-quartic-subresultant-checkpoint-v1.pickle"
)
DEFAULT_SPECIALIZED_FACTORIZATION = (
    RESULTS / "elkies-k3-q80-third-q12-exact-discriminant-specialization-v1.json"
)
DEFAULT_GENERIC_LINEAR = (
    RESULTS / "elkies-k3-q80-third-q12-exact-generic-linear-conductor-v1.json"
)
DEFAULT_H_CANDIDATE = (
    RESULTS / "elkies-k3-q80-third-q12-quartic-denominator-candidate-v1.json"
)
DEFAULT_GENERIC_QUARTIC_JET = (
    ROOT
    / "artifacts/local/elkies-k3/"
    / "q80-third-q12-exact-generic-quartic-jet-v1.json"
)
DEFAULT_GENERIC_QUARTIC_FACTORIZATION = (
    RESULTS / "elkies-k3-q80-third-q12-exact-generic-quartic-factorization-v1.json"
)
SUBRESULTANT_CHECKPOINT_SCHEMA = (
    "elkies-k3-q80-third-q12-exact-quartic-subresultant-checkpoint-v1"
)
COEFFICIENT = re.compile(
    r"^(-?\d+)/(\d+)\*theta\^2 ([+-]) (\d+)/(\d+)$"
)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational(record):
    return mpq(mpz(record["numerator"]), mpz(record["denominator"]))


def rational_record(value):
    value = mpq(value)
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def rational_height_bits(value):
    value = mpq(value)
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def field_record(value):
    coefficients = list(value.polynomial())
    coefficients += [QQ.zero()] * (2 - len(coefficients))
    return [rational_record(coefficients[0]), rational_record(coefficients[1])]


def polynomial_record(value):
    if hasattr(value, "list"):
        coefficients = value.list()
    else:
        generator = value.parent().gen()
        coefficients = [
            value.monomial_coefficient(generator**degree)
            for degree in range(value.degree() + 1)
        ]
    return [field_record(coefficient) for coefficient in coefficients]


def pair_record(value):
    return [rational_record(value[0]), rational_record(value[1])]


def pair_polynomial_record(value):
    return [pair_record(coefficient) for coefficient in value]


def pair_from_record(value):
    return rational(value[0]), rational(value[1])


def pair_polynomial_from_record(value):
    return p_trim([pair_from_record(coefficient) for coefficient in value])


ZERO = (mpz(0), mpz(0))
ONE = (mpz(1), mpz(0))


def k_add(left, right):
    return left[0] + right[0], left[1] + right[1]


def k_neg(value):
    return -value[0], -value[1]


def k_mul(left, right, omega_square):
    return (
        left[0] * right[0] + left[1] * right[1] * omega_square,
        left[0] * right[1] + left[1] * right[0],
    )


def k_scale(value, scalar):
    return value[0] * scalar, value[1] * scalar


def k_inverse(value, omega_square):
    norm = value[0] ** 2 - value[1] ** 2 * omega_square
    if not norm:
        raise ZeroDivisionError("quadratic coefficient has zero norm")
    return mpq(value[0]) / mpq(norm), mpq(-value[1]) / mpq(norm)


def p_trim(value):
    value = list(value)
    while value and value[-1] == ZERO:
        value.pop()
    return value or [ZERO]


def p_add(left, right):
    answer = [ZERO] * max(len(left), len(right))
    for index in range(len(answer)):
        answer[index] = k_add(
            left[index] if index < len(left) else ZERO,
            right[index] if index < len(right) else ZERO,
        )
    return p_trim(answer)


def p_neg(value):
    return p_trim([k_neg(coefficient) for coefficient in value])


def p_scale(value, scalar):
    return p_trim([k_scale(coefficient, scalar) for coefficient in value])


def p_k_scale(value, scalar, omega_square):
    return p_trim([k_mul(coefficient, scalar, omega_square) for coefficient in value])


def p_mul(left, right, omega_square):
    answer = [ZERO] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] = k_add(
                answer[left_index + right_index],
                k_mul(left_value, right_value, omega_square),
            )
    return p_trim(answer)


def p_pow(value, exponent, omega_square):
    answer = [ONE]
    power = value
    while exponent:
        if exponent & 1:
            answer = p_mul(answer, power, omega_square)
        exponent //= 2
        if exponent:
            power = p_mul(power, power, omega_square)
    return answer


def k_pow(value, exponent, field_square):
    answer = ONE
    power = value
    while exponent:
        if exponent & 1:
            answer = k_mul(answer, power, field_square)
        exponent //= 2
        if exponent:
            power = k_mul(power, power, field_square)
    return answer


def s_add(left, right):
    return [p_add(left[index], right[index]) for index in range(4)]


def s_neg(value):
    return [p_neg(coefficient) for coefficient in value]


def s_scale(value, scalar):
    return [p_scale(coefficient, scalar) for coefficient in value]


def s_mul(left, right, field_square):
    answer = [[ZERO] for _ in range(4)]
    for left_degree in range(4):
        for right_degree in range(4 - left_degree):
            answer[left_degree + right_degree] = p_add(
                answer[left_degree + right_degree],
                p_mul(left[left_degree], right[right_degree], field_square),
            )
    return answer


def s_pow(value, exponent, field_square):
    answer = [[ONE], [ZERO], [ZERO], [ZERO]]
    power = value
    while exponent:
        if exponent & 1:
            answer = s_mul(answer, power, field_square)
        exponent //= 2
        if exponent:
            power = s_mul(power, power, field_square)
    return answer


def j_add(left, right):
    return [p_add(left[0], right[0]), p_add(left[1], right[1])]


def j_neg(value):
    return [p_neg(value[0]), p_neg(value[1])]


def j_scale(value, scalar):
    return [p_scale(value[0], scalar), p_scale(value[1], scalar)]


def j_mul(left, right, field_square):
    return [
        p_mul(left[0], right[0], field_square),
        p_add(
            p_mul(left[0], right[1], field_square),
            p_mul(left[1], right[0], field_square),
        ),
    ]


def j_pow(value, exponent, field_square):
    answer = [[ONE], [ZERO]]
    power = value
    while exponent:
        if exponent & 1:
            answer = j_mul(answer, power, field_square)
        exponent //= 2
        if exponent:
            power = j_mul(power, power, field_square)
    return answer


def rational_modulus(value, modulus):
    value = mpq(value)
    denominator = ZZ(value.denominator)
    if gcd(mpz(denominator), mpz(modulus)) != 1:
        raise ZeroDivisionError("rational denominator is not invertible modulo modulus")
    return ZZ(value.numerator) * inverse_mod(denominator, modulus) % modulus


def b_trim(value):
    value = [p_trim(coefficient) for coefficient in value]
    while value and value[-1] == [ZERO]:
        value.pop()
    return value or [[ZERO]]


def b_add(left, right):
    answer = [[ZERO] for _ in range(max(len(left), len(right)))]
    for index in range(len(answer)):
        answer[index] = p_add(
            left[index] if index < len(left) else [ZERO],
            right[index] if index < len(right) else [ZERO],
        )
    return b_trim(answer)


def b_neg(value):
    return b_trim([p_neg(coefficient) for coefficient in value])


def b_scale(value, scalar):
    return b_trim([p_scale(coefficient, scalar) for coefficient in value])


def b_mul(left, right, field_square):
    answer = [[ZERO] for _ in range(len(left) + len(right) - 1)]
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] = p_add(
                answer[left_index + right_index],
                p_mul(left_value, right_value, field_square),
            )
    return b_trim(answer)


def b_pow(value, exponent, field_square):
    answer = [[ONE]]
    power = value
    while exponent:
        if exponent & 1:
            answer = b_mul(answer, power, field_square)
        exponent //= 2
        if exponent:
            power = b_mul(power, power, field_square)
    return answer


def b_divide_monic_linear(value, constant, field_square):
    """Divide in W by W+constant, with coefficients in K[V]."""
    value = b_trim(value)
    if len(value) < 2:
        return [[ZERO]], value[0]
    quotient = [[ZERO] for _ in range(len(value) - 1)]
    quotient[-1] = value[-1]
    for degree in range(len(value) - 2, 0, -1):
        quotient[degree - 1] = p_add(
            value[degree],
            p_neg(p_k_scale(quotient[degree], constant, field_square)),
        )
    remainder = p_add(
        value[0],
        p_neg(p_k_scale(quotient[0], constant, field_square)),
    )
    return b_trim(quotient), p_trim(remainder)


def b_pseudo_divmod(dividend, divisor, field_square):
    """Fraction-free W pseudo-division over the domain K[V]."""
    remainder = b_trim(dividend)
    divisor = b_trim(divisor)
    if divisor == [[ZERO]]:
        raise ZeroDivisionError
    divisor_degree = len(divisor) - 1
    divisor_lead = divisor[-1]
    quotient = [[ZERO]]
    steps = 0
    while remainder != [[ZERO]] and len(remainder) - 1 >= divisor_degree:
        shift = len(remainder) - 1 - divisor_degree
        remainder_lead = remainder[-1]
        remainder = b_mul(remainder, [divisor_lead], field_square)
        quotient = b_mul(quotient, [divisor_lead], field_square)
        while len(quotient) <= shift:
            quotient.append([ZERO])
        quotient[shift] = p_add(quotient[shift], remainder_lead)
        subtraction = [[ZERO] for _ in range(shift)] + [
            p_mul(remainder_lead, coefficient, field_square)
            for coefficient in divisor
        ]
        remainder = b_add(remainder, b_neg(subtraction))
        steps += 1
        print(
            "Q80Q12GENERICPSEUDO|"
            f"step={steps}|remainder_W_degree="
            f"{-1 if remainder == [[ZERO]] else len(remainder)-1}",
            flush=True,
        )
    return b_trim(quotient), b_trim(remainder), divisor_lead, steps


def p_divide_monic_linear(value, constant, field_square):
    """Divide by W+constant, returning quotient and exact remainder."""
    value = p_trim(value)
    if len(value) < 2:
        return [ZERO], value[0]
    quotient = [ZERO] * (len(value) - 1)
    quotient[-1] = value[-1]
    for degree in range(len(value) - 2, 0, -1):
        quotient[degree - 1] = k_add(
            value[degree],
            k_neg(k_mul(constant, quotient[degree], field_square)),
        )
    remainder = k_add(
        value[0], k_neg(k_mul(constant, quotient[0], field_square))
    )
    return p_trim(quotient), remainder


def p_derivative(value):
    if len(value) <= 1:
        return [ZERO]
    return p_trim([k_scale(value[index], index) for index in range(1, len(value))])


def p_primitive_rational(value):
    """Remove one rational scalar from a K[W] polynomial.

    This does not attempt algebraic-integer content in K.  It is enough to
    prevent the rational denominators introduced by field division from
    accumulating across a low-degree Euclidean remainder sequence.
    """
    value = p_trim(value)
    denominator = mpz(1)
    for coefficient in value:
        for coordinate in coefficient:
            denominator = lcm(denominator, coordinate.denominator)
    integers = [
        (
            coefficient[0].numerator * (denominator // coefficient[0].denominator),
            coefficient[1].numerator * (denominator // coefficient[1].denominator),
        )
        for coefficient in value
    ]
    content = mpz(0)
    for coefficient in integers:
        content = gcd(content, abs(coefficient[0]))
        content = gcd(content, abs(coefficient[1]))
    if not content:
        return [ZERO]
    result = [(a // content, b // content) for a, b in integers]
    first = next(coordinate for coefficient in reversed(result) for coordinate in coefficient if coordinate)
    if first < 0:
        result = [(-a, -b) for a, b in result]
    return p_trim(result)


def p_divmod_field(dividend, divisor, field_square):
    dividend = p_trim(dividend)
    divisor = p_trim(divisor)
    if divisor == [ZERO]:
        raise ZeroDivisionError
    if len(dividend) < len(divisor):
        return [ZERO], dividend
    quotient = [ZERO] * (len(dividend) - len(divisor) + 1)
    remainder = list(dividend)
    inverse_lead = k_inverse(divisor[-1], field_square)
    while remainder != [ZERO] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = k_mul(remainder[-1], inverse_lead, field_square)
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            position = shift + index
            remainder[position] = k_add(
                remainder[position], k_neg(k_mul(coefficient, value, field_square))
            )
        remainder = p_trim(remainder)
    return p_trim(quotient), remainder


def p_inverse_mod(value, modulus, field_square):
    """Return the inverse of value modulo modulus over the quadratic field."""
    old_remainder = p_trim(modulus)
    remainder = p_divmod_field(p_trim(value), old_remainder, field_square)[1]
    old_coefficient = [ZERO]
    coefficient = [ONE]
    while remainder != [ZERO]:
        quotient, next_remainder = p_divmod_field(
            old_remainder, remainder, field_square
        )
        old_remainder, remainder = remainder, next_remainder
        old_coefficient, coefficient = coefficient, p_add(
            old_coefficient,
            p_neg(p_mul(quotient, coefficient, field_square)),
        )
    if len(old_remainder) != 1 or old_remainder[0] == ZERO:
        raise ArithmeticError("polynomials are not coprime in modular inverse")
    inverse_gcd = k_inverse(old_remainder[0], field_square)
    candidate = p_k_scale(old_coefficient, inverse_gcd, field_square)
    _, candidate = p_divmod_field(candidate, modulus, field_square)
    product = p_mul(value, candidate, field_square)
    _, product_remainder = p_divmod_field(product, modulus, field_square)
    if product_remainder != [ONE]:
        raise ArithmeticError("polynomial modular inverse fails exact replay")
    return candidate


def p_pseudo_remainder(dividend, divisor, field_square):
    """Primitive pseudo-remainder using no quadratic-field inversions."""
    remainder = p_primitive_rational(dividend)
    divisor = p_primitive_rational(divisor)
    if divisor == [ZERO]:
        raise ZeroDivisionError
    divisor_degree = len(divisor) - 1
    divisor_lead = divisor[-1]
    while remainder != [ZERO] and len(remainder) - 1 >= divisor_degree:
        shift = len(remainder) - 1 - divisor_degree
        remainder_lead = remainder[-1]
        scaled = [k_mul(value, divisor_lead, field_square) for value in remainder]
        for index, value in enumerate(divisor):
            position = shift + index
            scaled[position] = k_add(
                scaled[position],
                k_neg(k_mul(remainder_lead, value, field_square)),
            )
        remainder = p_primitive_rational(p_trim(scaled))
    return remainder


def p_pseudo_remainder_raw(dividend, divisor, field_square):
    """Canonical pseudo-remainder used by Brown's subresultant PRS."""
    remainder = p_trim(dividend)
    divisor = p_trim(divisor)
    if divisor == [ZERO]:
        raise ZeroDivisionError
    divisor_degree = len(divisor) - 1
    divisor_lead = divisor[-1]
    remaining_power = len(remainder) - 1 - divisor_degree + 1
    while remainder != [ZERO] and len(remainder) - 1 >= divisor_degree:
        shift = len(remainder) - 1 - divisor_degree
        remaining_power -= 1
        remainder_lead = remainder[-1]
        scaled = [k_mul(value, divisor_lead, field_square) for value in remainder]
        for index, value in enumerate(divisor):
            position = shift + index
            scaled[position] = k_add(
                scaled[position],
                k_neg(k_mul(remainder_lead, value, field_square)),
            )
        remainder = p_trim(scaled)
    if remainder != [ZERO] and remaining_power:
        remainder = p_k_scale(
            remainder,
            k_pow(divisor_lead, remaining_power, field_square),
            field_square,
        )
    return remainder


def p_k_divide(value, divisor, field_square):
    inverse = k_inverse(divisor, field_square)
    return p_k_scale(value, inverse, field_square)


def p_monic(value, field_square):
    value = p_trim(value)
    inverse_lead = k_inverse(value[-1], field_square)
    return p_k_scale(value, inverse_lead, field_square)


def p_maximum_coordinate_bits(value):
    return max(
        rational_height_bits(coordinate)
        for coefficient in value
        for coordinate in coefficient
    )


def p_custom_gcd(left, right, field_square):
    left = p_primitive_rational(left)
    right = p_primitive_rational(right)
    steps = []
    while right != [ZERO]:
        started = time.monotonic()
        remainder = p_pseudo_remainder(left, right, field_square)
        record = {
            "input_degrees": [len(left) - 1, len(right) - 1],
            "remainder_degree": -1 if remainder == [ZERO] else len(remainder) - 1,
            "remainder_maximum_coordinate_bits": (
                0 if remainder == [ZERO] else p_maximum_coordinate_bits(remainder)
            ),
            "seconds": time.monotonic() - started,
        }
        steps.append(record)
        print(
            "Q80Q12CUSTOMPRS|"
            f"degrees={record['input_degrees']}|remainder={record['remainder_degree']}|"
            f"bits={record['remainder_maximum_coordinate_bits']}|"
            f"seconds={record['seconds']:.3f}",
            flush=True,
        )
        left, right = right, remainder
    return p_monic(left, field_square), steps


def write_subresultant_checkpoint(path, identity, state):
    """Atomically persist one exact Brown-PRS state.

    The checkpoint is a local, trusted Python pickle because its million-bit
    GMP integers are prohibitively wasteful in decimal JSON.  A small JSON
    companion records the identity and progress without loading the pickle.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": SUBRESULTANT_CHECKPOINT_SCHEMA,
        "identity": identity,
        "state": state,
    }
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("wb") as handle:
        pickle.dump(payload, handle, protocol=5)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    metadata = {
        "schema": SUBRESULTANT_CHECKPOINT_SCHEMA,
        "identity": identity,
        "status": state["status"],
        "completed_remainder_degrees": [
            record["remainder_degree"] for record in state["steps"]
        ],
        "next_remainder_degree": (
            -1 if state["h"] == [ZERO] else len(state["h"]) - 1
        ),
        "checkpoint_bytes": path.stat().st_size,
    }
    metadata_path = path.with_suffix(path.suffix + ".json")
    metadata_temporary = metadata_path.with_name(
        f".{metadata_path.name}.tmp-{os.getpid()}"
    )
    metadata_temporary.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    os.replace(metadata_temporary, metadata_path)
    print(
        "Q80Q12SUBRESCHECKPOINT|"
        f"status={state['status']}|next={metadata['next_remainder_degree']}|"
        f"bytes={metadata['checkpoint_bytes']}|path={path}",
        flush=True,
    )


def read_subresultant_checkpoint(path, identity):
    # Pickle is intentionally restricted to an explicitly selected local
    # checkpoint produced by this script.  Never use an untrusted file here.
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if payload.get("schema") != SUBRESULTANT_CHECKPOINT_SCHEMA:
        raise ArithmeticError("subresultant checkpoint schema mismatch")
    if payload.get("identity") != identity:
        raise ArithmeticError("subresultant checkpoint belongs to different inputs")
    state = payload.get("state")
    required = {"status", "f", "g", "m", "h", "lc", "c", "steps"}
    if not isinstance(state, dict) or not required.issubset(state):
        raise ArithmeticError("subresultant checkpoint state is incomplete")
    return state


def p_subresultant_gcd(
    left,
    right,
    field_square,
    checkpoint_path=None,
    checkpoint_identity=None,
    resume=True,
):
    """Brown subresultant PRS over the exact quadratic coefficient field."""
    if checkpoint_path is not None and checkpoint_identity is None:
        raise ValueError("checkpoint identity is required with a checkpoint path")
    if checkpoint_path is not None and resume and checkpoint_path.exists():
        state = read_subresultant_checkpoint(checkpoint_path, checkpoint_identity)
        f = state["f"]
        g = state["g"]
        m = state["m"]
        h = state["h"]
        lc = state["lc"]
        c = state["c"]
        steps = state["steps"]
        print(
            "Q80Q12SUBRESRESUME|"
            f"status={state['status']}|completed={len(steps)}|"
            f"next={-1 if h == [ZERO] else len(h)-1}|path={checkpoint_path}",
            flush=True,
        )
        if state["status"] == "COMPLETE":
            return p_monic(g, field_square), steps
    else:
        f = p_primitive_rational(left)
        g = p_primitive_rational(right)
        n = len(f) - 1
        m = len(g) - 1
        if n < m:
            f, g = g, f
            n, m = m, n
        d = n - m
        sign = -1 if (d + 1) % 2 else 1
        h = p_scale(p_pseudo_remainder_raw(f, g, field_square), sign)
        lc = g[-1]
        c = k_neg(k_pow(lc, d, field_square))
        steps = []
        if checkpoint_path is not None:
            write_subresultant_checkpoint(
                checkpoint_path,
                checkpoint_identity,
                {
                    "status": "IN_PROGRESS",
                    "f": f,
                    "g": g,
                    "m": m,
                    "h": h,
                    "lc": lc,
                    "c": c,
                    "steps": steps,
                },
            )
    while h != [ZERO]:
        k = len(h) - 1
        record = {
            "input_degrees": [len(f) - 1, len(g) - 1],
            "remainder_degree": k,
            "remainder_maximum_coordinate_bits": p_maximum_coordinate_bits(h),
        }
        steps.append(record)
        print(
            "Q80Q12SUBRES|"
            f"degrees={record['input_degrees']}|remainder={k}|"
            f"bits={record['remainder_maximum_coordinate_bits']}",
            flush=True,
        )
        old_m = m
        f, g, m = g, h, k
        d = old_m - k
        b = k_neg(k_mul(lc, k_pow(c, d, field_square), field_square))
        started = time.monotonic()
        raw = p_pseudo_remainder_raw(f, g, field_square)
        h = p_k_divide(raw, b, field_square) if raw != [ZERO] else [ZERO]
        if h != [ZERO]:
            print(
                "Q80Q12SUBRES_SCALE|"
                f"degree={len(h)-1}|bits={p_maximum_coordinate_bits(h)}|"
                f"seconds={time.monotonic()-started:.3f}",
                flush=True,
            )
        lc = g[-1]
        if d > 1:
            numerator = k_pow(k_neg(lc), d, field_square)
            denominator = k_pow(c, d - 1, field_square)
            c = k_mul(numerator, k_inverse(denominator, field_square), field_square)
        else:
            c = k_neg(lc)
        if checkpoint_path is not None:
            write_subresultant_checkpoint(
                checkpoint_path,
                checkpoint_identity,
                {
                    "status": "COMPLETE" if h == [ZERO] else "IN_PROGRESS",
                    "f": f,
                    "g": g,
                    "m": m,
                    "h": h,
                    "lc": lc,
                    "c": c,
                    "steps": steps,
                },
            )
    return p_monic(g, field_square), steps


def p_modular_gcd_degree(value, field_square, prime):
    constants = GF(prime)
    z_ring = PolynomialRing(constants, "z_mod")
    z_mod = z_ring.gen()
    finite = GF(prime**2, "delta_mod", modulus=z_mod**2 - constants(field_square))
    delta_mod = finite.gen()
    w_ring = PolynomialRing(finite, "W_mod")

    def reduce_coordinate(coordinate):
        return finite(
            constants(int(coordinate.numerator % prime))
            / constants(int(coordinate.denominator % prime))
        )

    polynomial = w_ring(
        [
            reduce_coordinate(coefficient[0])
            + reduce_coordinate(coefficient[1]) * delta_mod
            for coefficient in value
        ]
    )
    return int(polynomial.gcd(polynomial.derivative()).degree())


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--operands", type=Path, default=DEFAULT_OPERANDS)
parser.add_argument("--pencil", type=Path, default=DEFAULT_PENCIL)
parser.add_argument("--factor-lift", type=Path, default=DEFAULT_FACTOR_LIFT)
parser.add_argument(
    "--specialized-factorization",
    type=Path,
    default=DEFAULT_SPECIALIZED_FACTORIZATION,
)
parser.add_argument("--H-candidate", type=Path, default=DEFAULT_H_CANDIDATE)
parser.add_argument(
    "--generic-quartic-jet-artifact",
    type=Path,
    default=DEFAULT_GENERIC_QUARTIC_JET,
)
parser.add_argument("--base-value", default="0")
parser.add_argument("--output", type=Path)
parser.add_argument("--check", action="store_true")
parser.add_argument("--attempt-quartic-gcd", action="store_true")
parser.add_argument("--attempt-custom-prs", action="store_true")
parser.add_argument("--attempt-subresultant-prs", action="store_true")
parser.add_argument("--attempt-generic-quartic-jet", action="store_true")
parser.add_argument("--attempt-generic-quartic-division", action="store_true")
parser.add_argument(
    "--subresultant-checkpoint",
    type=Path,
    default=DEFAULT_SUBRESULTANT_CHECKPOINT,
    help="trusted local binary checkpoint written atomically after every Brown-PRS step",
)
parser.add_argument(
    "--restart-subresultant-prs",
    action="store_true",
    help="ignore and atomically replace an existing subresultant checkpoint",
)
parser.add_argument("--certify-generic-linear", action="store_true")
args = parser.parse_args()
selected_modes = [
    args.attempt_quartic_gcd,
    args.attempt_custom_prs,
    args.attempt_subresultant_prs,
    args.attempt_generic_quartic_jet,
    args.attempt_generic_quartic_division,
    args.certify_generic_linear,
]
if sum(selected_modes) > 1:
    parser.error("select exactly one factorization/certification mode")
if args.check and args.output is None:
    if args.certify_generic_linear:
        args.output = DEFAULT_GENERIC_LINEAR
    elif args.attempt_custom_prs or args.attempt_subresultant_prs:
        args.output = DEFAULT_SPECIALIZED_FACTORIZATION
    elif args.attempt_generic_quartic_jet:
        args.output = DEFAULT_GENERIC_QUARTIC_JET
    elif args.attempt_generic_quartic_division:
        args.output = DEFAULT_GENERIC_QUARTIC_FACTORIZATION
    else:
        parser.error("--check requires --output for the selected mode")
args.operands = args.operands.resolve()
args.pencil = args.pencil.resolve()
args.factor_lift = args.factor_lift.resolve()
args.specialized_factorization = args.specialized_factorization.resolve()
args.H_candidate = args.H_candidate.resolve()
args.generic_quartic_jet_artifact = args.generic_quartic_jet_artifact.resolve()
args.subresultant_checkpoint = args.subresultant_checkpoint.resolve()
if args.output:
    args.output = args.output.resolve()

operands = json.loads(args.operands.read_text())
pencil = json.loads(args.pencil.read_text())
factor_lift = json.loads(args.factor_lift.read_text())
q1 = rational(operands["biquadratic_field"]["q1"])
q2 = rational(operands["biquadratic_field"]["q2"])
base_value = mpq(args.base_value)
omega_square = 16 * q1 * q2
product = q1 * q2
product_numerator_root = isqrt(product.numerator)
if product_numerator_root**2 != product.numerator:
    raise ArithmeticError("reduced q1*q2 numerator is not a square")
delta_square = product.denominator
omega_to_delta = mpq(
    4 * product_numerator_root,
    product.denominator,
)

cubic_coefficients = [[ZERO] for _ in range(4)]
parsed_terms = []
terms = pencil["moving_equation"]["terms_T_W_x_coefficient_1_r"]
for v_degree, w_degree, x_degree, encoded in terms:
    assert len(encoded) == 1
    match = COEFFICIENT.fullmatch(encoded[0])
    if match is None:
        raise ArithmeticError("unexpected exact coefficient encoding")
    theta2 = mpq(mpz(match[1]), mpz(match[2]))
    sign = 1 if match[3] == "+" else -1
    constant = mpq(sign * mpz(match[4]), mpz(match[5]))
    # Work initially in the denominator-integral basis delta^2=D, where
    # omega=4*sqrt(N)*delta/D for reduced q1*q2=N/D.
    exact_value = (
        constant + theta2 * (q1 + q2),
        2 * theta2 * product_numerator_root / product.denominator,
    )
    parsed_terms.append((v_degree, w_degree, x_degree, exact_value))
    value = k_scale(exact_value, base_value**v_degree)
    while len(cubic_coefficients[x_degree]) <= w_degree:
        cubic_coefficients[x_degree].append(ZERO)
    cubic_coefficients[x_degree][w_degree] = k_add(
        cubic_coefficients[x_degree][w_degree], value
    )

cubic_coefficients = list(map(p_trim, cubic_coefficients))
if len(cubic_coefficients[3]) != 1 or cubic_coefficients[3][0] == ZERO:
    raise ArithmeticError("specialized cubic leading coefficient is not constant")

# Clear one rational denominator and content before any polynomial products.
# The general cubic discriminant is homogeneous of degree four, so this
# projective scaling changes only its constant factor and preserves L,Q,D.
common_denominator = mpz(1)
for polynomial in cubic_coefficients:
    for coefficient in polynomial:
        for coordinate in coefficient:
            common_denominator = lcm(common_denominator, coordinate.denominator)
integer_coordinates = []
for polynomial in cubic_coefficients:
    for coefficient in polynomial:
        for coordinate in coefficient:
            integer_coordinates.append(
                coordinate.numerator * (common_denominator // coordinate.denominator)
            )
content = mpz(0)
for coordinate in integer_coordinates:
    content = gcd(content, abs(coordinate))
if not content:
    raise ArithmeticError("specialized cubic is zero")
cubic_coefficients = [
    [
        tuple(
            coordinate.numerator
            * (common_denominator // coordinate.denominator)
            // content
            for coordinate in coefficient
        )
        for coefficient in polynomial
    ]
    for polynomial in cubic_coefficients
]
a, b, c, d = (
    cubic_coefficients[3],
    cubic_coefficients[2],
    cubic_coefficients[1],
    cubic_coefficients[0],
)
discriminant = p_add(
    p_add(
        p_add(
            p_mul(p_pow(b, 2, delta_square), p_pow(c, 2, delta_square), delta_square),
            p_neg(p_scale(p_mul(a, p_pow(c, 3, delta_square), delta_square), 4)),
        ),
        p_neg(p_scale(p_mul(p_pow(b, 3, delta_square), d, delta_square), 4)),
    ),
    p_add(
        p_neg(
            p_scale(
                p_mul(p_pow(a, 2, delta_square), p_pow(d, 2, delta_square), delta_square),
                27,
            )
        ),
        p_scale(
            p_mul(p_mul(p_mul(a, b, delta_square), c, delta_square), d, delta_square),
            18,
        ),
    ),
)

# The high-precision lift makes the unique linear factor much cheaper to
# recover than a number-field gcd.  At V=0 it is W+r with rational r.  Accept
# the rational reconstruction only after three literal exact divisions of
# the independently constructed characteristic-zero discriminant.
if base_value != 0:
    raise NotImplementedError("the current exact factor-lift seed is implemented at V=0")
modulus = ZZ(factor_lift["specialization"]["modulus"])
linear_record = factor_lift["factorization"]["L"][
    "coefficients_low_to_high_W"
][0]
if linear_record["degrees_numerator_denominator"] != [0, 0]:
    raise ArithmeticError("lifted L constant is not base-independent")
linear_numerator = linear_record[
    "numerator_coefficients_low_to_high_U_1_omega"
][0]
linear_denominator = linear_record[
    "denominator_coefficients_low_to_high_U_1_omega"
][0]
if linear_numerator[1] or linear_denominator[1]:
    raise ArithmeticError("lifted L constant is not rational")
linear_residue = (
    ZZ(linear_numerator[0]) * inverse_mod(ZZ(linear_denominator[0]), modulus)
) % modulus
linear_rational = linear_residue.rational_reconstruction(modulus)
linear_constant = (
    mpq(
        mpz(str(linear_rational.numerator())),
        mpz(str(linear_rational.denominator())),
    ),
    mpq(0),
)
residual_discriminant = discriminant
for multiplicity in range(3):
    residual_discriminant, remainder = p_divide_monic_linear(
        residual_discriminant, linear_constant, delta_square
    )
    if remainder != ZERO:
        raise ArithmeticError(
            f"reconstructed linear factor fails exact division {multiplicity + 1}"
        )
_, fourth_remainder = p_divide_monic_linear(
    residual_discriminant, linear_constant, delta_square
)
if fourth_remainder == ZERO:
    raise ArithmeticError("reconstructed linear factor has multiplicity above three")

generic_linear = None
if args.certify_generic_linear:
    # Expand the raw cubic coefficients in t=W+r through t^3, retaining V as
    # a polynomial.  The discriminant is homogeneous of degree four in cubic
    # coefficients, so its V-degree is at most 8.  Literal vanishing of the
    # first three t coefficients proves (W+r)^3 divides Delta over K[V,W].
    shifted_cubic = [
        [[ZERO] for _ in range(4)] for _ in range(4)
    ]
    generic_denominator = mpz(1)
    for _, _, _, exact_value in parsed_terms:
        for coordinate in exact_value:
            generic_denominator = lcm(generic_denominator, coordinate.denominator)
    generic_integer_terms = []
    generic_content = mpz(0)
    for v_degree, w_degree, x_degree, exact_value in parsed_terms:
        integer_value = tuple(
            coordinate.numerator
            * (generic_denominator // coordinate.denominator)
            for coordinate in exact_value
        )
        generic_integer_terms.append((v_degree, w_degree, x_degree, integer_value))
        for coordinate in integer_value:
            generic_content = gcd(generic_content, abs(coordinate))
    generic_integer_terms = [
        (
            v_degree,
            w_degree,
            x_degree,
            tuple(coordinate // generic_content for coordinate in integer_value),
        )
        for v_degree, w_degree, x_degree, integer_value in generic_integer_terms
    ]
    linear_numerator_exact = mpz(str(linear_rational.numerator()))
    linear_denominator_exact = mpz(str(linear_rational.denominator()))
    maximum_w_degree = max(item[1] for item in generic_integer_terms)
    assert maximum_w_degree == 9
    # Use S=d*W+n for r=n/d and multiply the whole cubic by d^9.
    # A term c*W^w contributes
    # c*d^(9-w)*(S-n)^w, so every shifted coefficient is integral.
    for v_degree, w_degree, x_degree, exact_value in generic_integer_terms:
        for t_degree in range(min(3, w_degree) + 1):
            multiplier = (
                comb(w_degree, t_degree)
                * (-linear_numerator_exact) ** (w_degree - t_degree)
                * linear_denominator_exact ** (maximum_w_degree - w_degree)
            )
            contribution = k_scale(exact_value, multiplier)
            coefficient = shifted_cubic[x_degree][t_degree]
            while len(coefficient) <= v_degree:
                coefficient.append(ZERO)
            coefficient[v_degree] = k_add(coefficient[v_degree], contribution)
            shifted_cubic[x_degree][t_degree] = p_trim(coefficient)
    series_a, series_b, series_c, series_d = (
        shifted_cubic[3],
        shifted_cubic[2],
        shifted_cubic[1],
        shifted_cubic[0],
    )
    discriminant_series = s_add(
        s_add(
            s_add(
                s_mul(
                    s_pow(series_b, 2, delta_square),
                    s_pow(series_c, 2, delta_square),
                    delta_square,
                ),
                s_neg(s_scale(s_mul(series_a, s_pow(series_c, 3, delta_square), delta_square), 4)),
            ),
            s_neg(
                s_scale(
                    s_mul(s_pow(series_b, 3, delta_square), series_d, delta_square),
                    4,
                )
            ),
        ),
        s_add(
            s_neg(
                s_scale(
                    s_mul(
                        s_pow(series_a, 2, delta_square),
                        s_pow(series_d, 2, delta_square),
                        delta_square,
                    ),
                    27,
                )
            ),
            s_scale(
                s_mul(
                    s_mul(
                        s_mul(series_a, series_b, delta_square),
                        series_c,
                        delta_square,
                    ),
                    series_d,
                    delta_square,
                ),
                18,
            ),
        ),
    )
    if any(coefficient != [ZERO] for coefficient in discriminant_series[:3]):
        raise ArithmeticError("the reconstructed linear factor is not generically triple")
    t3 = p_trim(discriminant_series[3])
    if t3 == [ZERO] or len(t3) - 1 > 8:
        raise ArithmeticError("the generic t^3 discriminant coefficient is invalid")
    canonical_t3 = json.dumps(
        [[rational_record(value[0]), rational_record(value[1])] for value in t3],
        sort_keys=True,
    )
    generic_linear = {
        "base_degree_bound": 8,
        "vanishing_t_degrees": [0, 1, 2],
        "first_nonzero_t_degree": 3,
        "t3_base_degree": len(t3) - 1,
        "t3_nonzero_coefficient_count": sum(value != ZERO for value in t3),
        "t3_maximum_coordinate_bits": max(
            rational_height_bits(coordinate)
            for coefficient in t3
            for coordinate in coefficient
        ),
        "t3_sha256": hashlib.sha256(canonical_t3.encode()).hexdigest(),
        "integral_compilation": {
            "coefficient_common_denominator_bits": generic_denominator.bit_length(),
            "coefficient_integer_content_bits": generic_content.bit_length(),
            "shift": "S=denominator(r)*W+numerator(r)",
            "cubic_scale": "denominator(r)^9",
        },
        "conclusion": "(W+r)^3 divides the exact generic discriminant in K[V,W]",
    }

factor_data = None
degree_exponents = [[1, 3]]
custom_prs = None
if args.attempt_custom_prs:
    Q_pair, custom_prs_steps = p_custom_gcd(
        residual_discriminant,
        p_derivative(residual_discriminant),
        delta_square,
    )
    if len(Q_pair) - 1 != 4:
        raise ArithmeticError(
            f"custom residual gcd has degree {len(Q_pair) - 1}, expected four"
        )
    Q_pair_squared = p_mul(Q_pair, Q_pair, delta_square)
    D_pair, pair_remainder = p_divmod_field(
        p_monic(residual_discriminant, delta_square),
        Q_pair_squared,
        delta_square,
    )
    if pair_remainder != [ZERO]:
        raise ArithmeticError("custom exact Q^2 division has nonzero remainder")
    D_pair = p_monic(D_pair, delta_square)
    if len(D_pair) - 1 != 4:
        raise ArithmeticError("custom exact residual factor D is not quartic")
    custom_prs = {
        "steps": custom_prs_steps,
        "Q_coefficients_low_to_high_W_1_delta": pair_polynomial_record(Q_pair),
        "D_coefficients_low_to_high_W_1_delta": pair_polynomial_record(D_pair),
        "Q_maximum_coordinate_bits": p_maximum_coordinate_bits(Q_pair),
        "D_maximum_coordinate_bits": p_maximum_coordinate_bits(D_pair),
        "identity": "monic_specialized_discriminant=L^3*Q^2*D at V=0",
    }
    degree_exponents = [[4, 1], [4, 2], [1, 3]]
subresultant_prs = None
if args.attempt_subresultant_prs:
    test_Q = [(1, 1), (2, 0), (0, 1), (3, 0), ONE]
    test_D = [(5, 0), (1, 0), (4, 0), (0, 1), ONE]
    test_value = p_mul(p_mul(test_Q, test_Q, delta_square), test_D, delta_square)
    test_gcd, unused_test_steps = p_subresultant_gcd(
        test_value, p_derivative(test_value), delta_square
    )
    if len(test_gcd) - 1 != 4:
        raise ArithmeticError(
            f"internal subresultant self-test returned degree {len(test_gcd)-1}"
        )
    modular_gcd_degrees = {
        prime: p_modular_gcd_degree(residual_discriminant, delta_square, prime)
        for prime in (163, 191, 199)
    }
    print(f"Q80Q12SUBRESMOD|degrees={modular_gcd_degrees}", flush=True)
    if set(modular_gcd_degrees.values()) != {4}:
        raise ArithmeticError(
            f"exact pair residual has unexpected modular gcd degrees {modular_gcd_degrees}"
        )
    subresultant_checkpoint_identity = {
        "algorithm": "Brown subresultant PRS over QQ(delta), delta^2=D",
        "base_value": str(base_value),
        "delta_square": str(delta_square),
        "linear_factor_constant": str(linear_constant[0]),
        "linear_stripped_residual_degree": len(residual_discriminant) - 1,
        "linear_stripped_residual_maximum_coordinate_bits": (
            p_maximum_coordinate_bits(residual_discriminant)
        ),
        "operands_sha256": sha256(args.operands),
        "pencil_sha256": sha256(args.pencil),
        "factor_lift_sha256": sha256(args.factor_lift),
    }
    Q_pair, subresultant_steps = p_subresultant_gcd(
        residual_discriminant,
        p_derivative(residual_discriminant),
        delta_square,
        checkpoint_path=args.subresultant_checkpoint,
        checkpoint_identity=subresultant_checkpoint_identity,
        resume=not args.restart_subresultant_prs,
    )
    if len(Q_pair) - 1 != 4:
        raise ArithmeticError(
            f"subresultant residual gcd has degree {len(Q_pair) - 1}, expected four"
        )
    Q_pair_squared = p_mul(Q_pair, Q_pair, delta_square)
    D_pair, pair_remainder = p_divmod_field(
        p_monic(residual_discriminant, delta_square),
        Q_pair_squared,
        delta_square,
    )
    if pair_remainder != [ZERO]:
        raise ArithmeticError("subresultant exact Q^2 division has nonzero remainder")
    D_pair = p_monic(D_pair, delta_square)
    if len(D_pair) - 1 != 4:
        raise ArithmeticError("subresultant exact residual factor D is not quartic")
    subresultant_prs = {
        "steps": subresultant_steps,
        "Q_coefficients_low_to_high_W_1_delta": pair_polynomial_record(Q_pair),
        "D_coefficients_low_to_high_W_1_delta": pair_polynomial_record(D_pair),
        "Q_maximum_coordinate_bits": p_maximum_coordinate_bits(Q_pair),
        "D_maximum_coordinate_bits": p_maximum_coordinate_bits(D_pair),
        "identity": "monic_specialized_discriminant=L^3*Q^2*D at V=0",
    }
    degree_exponents = [[4, 1], [4, 2], [1, 3]]
generic_quartic_jet = None
if args.attempt_generic_quartic_jet:
    if base_value != 0:
        raise NotImplementedError("the generic quartic jet is anchored at V=0")
    specialized = json.loads(args.specialized_factorization.read_text())
    H_artifact = json.loads(args.H_candidate.read_text())
    if specialized.get("status") != "PASS_EXACT_SPECIALIZED_L3_Q2_D_FACTORIZATION":
        raise ArithmeticError("specialized factorization artifact is not certified")
    if specialized["base_value"] != rational_record(mpq(0)):
        raise ArithmeticError("specialized factorization is not anchored at V=0")
    if specialized["quadratic_field"]["delta_square"] != str(delta_square):
        raise ArithmeticError("specialized factorization uses a different descent field")
    q0 = pair_polynomial_from_record(
        specialized["monic_factors_by_exponent"]["2"]
    )
    d0 = pair_polynomial_from_record(
        specialized["monic_factors_by_exponent"]["1"]
    )
    if len(q0) - 1 != 4 or len(d0) - 1 != 4:
        raise ArithmeticError("specialized exact factors are not quartics")

    # Clear one global rational denominator/content from the exact moving
    # cubic, then retain only its V^0 and V^1 coefficients.  This is the exact
    # dual-number fibre over K[V]/(V^2), not finite differencing.
    jet_denominator = mpz(1)
    for _, _, _, exact_value in parsed_terms:
        for coordinate in exact_value:
            jet_denominator = lcm(jet_denominator, coordinate.denominator)
    jet_integer_terms = []
    jet_content = mpz(0)
    for v_degree, w_degree, x_degree, exact_value in parsed_terms:
        integer_value = tuple(
            coordinate.numerator
            * (jet_denominator // coordinate.denominator)
            for coordinate in exact_value
        )
        jet_integer_terms.append((v_degree, w_degree, x_degree, integer_value))
        for coordinate in integer_value:
            jet_content = gcd(jet_content, abs(coordinate))
    jet_integer_terms = [
        (
            v_degree,
            w_degree,
            x_degree,
            tuple(coordinate // jet_content for coordinate in integer_value),
        )
        for v_degree, w_degree, x_degree, integer_value in jet_integer_terms
    ]
    cubic_jets = [[[ZERO], [ZERO]] for _ in range(4)]
    for v_degree, w_degree, x_degree, exact_value in jet_integer_terms:
        if v_degree > 1:
            continue
        coefficient = cubic_jets[x_degree][v_degree]
        while len(coefficient) <= w_degree:
            coefficient.append(ZERO)
        coefficient[w_degree] = k_add(coefficient[w_degree], exact_value)
        cubic_jets[x_degree][v_degree] = p_trim(coefficient)
    jet_a, jet_b, jet_c, jet_d = (
        cubic_jets[3],
        cubic_jets[2],
        cubic_jets[1],
        cubic_jets[0],
    )
    discriminant_jet = j_add(
        j_add(
            j_add(
                j_mul(
                    j_pow(jet_b, 2, delta_square),
                    j_pow(jet_c, 2, delta_square),
                    delta_square,
                ),
                j_neg(
                    j_scale(
                        j_mul(jet_a, j_pow(jet_c, 3, delta_square), delta_square),
                        4,
                    )
                ),
            ),
            j_neg(
                j_scale(
                    j_mul(j_pow(jet_b, 3, delta_square), jet_d, delta_square),
                    4,
                )
            ),
        ),
        j_add(
            j_neg(
                j_scale(
                    j_mul(
                        j_pow(jet_a, 2, delta_square),
                        j_pow(jet_d, 2, delta_square),
                        delta_square,
                    ),
                    27,
                )
            ),
            j_scale(
                j_mul(
                    j_mul(
                        j_mul(jet_a, jet_b, delta_square),
                        jet_c,
                        delta_square,
                    ),
                    jet_d,
                    delta_square,
                ),
                18,
            ),
        ),
    )
    residual_jet = []
    for jet_coefficient in discriminant_jet:
        stripped = jet_coefficient
        for _ in range(3):
            stripped, remainder = p_divide_monic_linear(
                stripped, linear_constant, delta_square
            )
            if remainder != ZERO:
                raise ArithmeticError("generic discriminant jet is not divisible by L^3")
        residual_jet.append(stripped)
    residual0, residual1 = residual_jet
    residual_degree = len(residual0) - 1
    if residual_degree != 12 or len(residual1) - 1 > residual_degree:
        raise ArithmeticError("linear-stripped discriminant jet has wrong W-degree")
    leading0 = residual0[-1]
    leading1 = (
        residual1[residual_degree]
        if len(residual1) > residual_degree
        else ZERO
    )
    inverse_leading0 = k_inverse(leading0, delta_square)
    r0 = p_k_scale(residual0, inverse_leading0, delta_square)
    r1 = p_add(
        p_k_scale(residual1, inverse_leading0, delta_square),
        p_neg(
            p_k_scale(
                r0,
                k_mul(leading1, inverse_leading0, delta_square),
                delta_square,
            )
        ),
    )
    expected_r0 = p_mul(p_mul(q0, q0, delta_square), d0, delta_square)
    if r0 != expected_r0:
        raise ArithmeticError("exact V=0 jet does not replay certified Q0^2*D0")
    quotient_r1_q0, remainder_r1_q0 = p_divmod_field(
        r1, q0, delta_square
    )
    if remainder_r1_q0 != [ZERO]:
        raise ArithmeticError("Q0 does not divide the exact first discriminant jet")
    _, tangent_numerator = p_divmod_field(
        quotient_r1_q0, q0, delta_square
    )
    inverse_d0 = p_inverse_mod(d0, q0, delta_square)
    tangent_product = p_mul(tangent_numerator, inverse_d0, delta_square)
    _, q1 = p_divmod_field(tangent_product, q0, delta_square)
    q1 = p_scale(q1, mpq(1, 2))
    if len(q1) - 1 > 3:
        raise ArithmeticError("monic quartic tangent has degree above three")
    differentiated_numerator = p_add(
        r1,
        p_neg(
            p_scale(
                p_mul(
                    p_mul(q0, q1, delta_square),
                    d0,
                    delta_square,
                ),
                2,
            )
        ),
    )
    q0_squared = p_mul(q0, q0, delta_square)
    d1, differentiated_remainder = p_divmod_field(
        differentiated_numerator, q0_squared, delta_square
    )
    if differentiated_remainder != [ZERO] or len(d1) - 1 > 3:
        raise ArithmeticError("recovered quartic tangent fails differentiated identity")

    H_candidate = H_artifact["candidate"]
    if H_candidate["delta_square"] != str(delta_square):
        raise ArithmeticError("H candidate uses a different descent field")
    h0 = (
        rational(H_candidate["h0_rational"]),
        rational(H_candidate["h0_delta"]),
    )
    q0_coefficients = q0 + [ZERO] * (5 - len(q0))
    q1_coefficients = q1 + [ZERO] * (4 - len(q1))
    q_hat = []
    for w_degree in range(4):
        constant = k_mul(h0, q0_coefficients[w_degree], delta_square)
        slope = k_add(
            q0_coefficients[w_degree],
            k_mul(h0, q1_coefficients[w_degree], delta_square),
        )
        q_hat.append(p_trim([constant, slope]))
    q_hat.append([h0, ONE])

    # The p-adic lift records these same numerator coefficients in the omega
    # basis.  The exact V=0 anchor removes the former projective ambiguity.
    q_lift_records = factor_lift["factorization"]["Q"][
        "coefficients_low_to_high_W"
    ]
    p_adic_mismatches = []
    for w_degree in range(4):
        lifted = q_lift_records[w_degree][
            "numerator_coefficients_low_to_high_U_1_omega"
        ]
        if len(lifted) != 2:
            raise ArithmeticError("p-adic quartic numerator is not linear in V")
        for v_degree in range(2):
            exact_delta = q_hat[w_degree][v_degree]
            exact_omega = (
                exact_delta[0],
                exact_delta[1] / mpq(omega_to_delta),
            )
            for coordinate in range(2):
                exact_residue = rational_modulus(
                    exact_omega[coordinate], modulus
                )
                lifted_residue = ZZ(lifted[v_degree][coordinate]) % modulus
                if exact_residue != lifted_residue:
                    p_adic_mismatches.append(
                        [w_degree, v_degree, coordinate]
                    )
    if p_adic_mismatches:
        raise ArithmeticError(
            "jet-derived quartic numerator fails p-adic replay at "
            f"{p_adic_mismatches}"
        )
    generic_quartic_jet = {
        "method": "exact first-order Hensel lift at V=0",
        "identity": "R1=2*Q0*Q1*D0+Q0^2*D1",
        "R1_maximum_coordinate_bits": p_maximum_coordinate_bits(r1),
        "Q1_coefficients_low_to_high_W_1_delta": pair_polynomial_record(q1),
        "Q1_maximum_coordinate_bits": p_maximum_coordinate_bits(q1),
        "D1_coefficients_low_to_high_W_1_delta": pair_polynomial_record(d1),
        "D1_maximum_coordinate_bits": p_maximum_coordinate_bits(d1),
        "Q_numerator_coefficients_low_to_high_W_then_V_1_delta": [
            pair_polynomial_record(coefficient) for coefficient in q_hat
        ],
        "Q_numerator_maximum_coordinate_bits": max(
            p_maximum_coordinate_bits(coefficient) for coefficient in q_hat
        ),
        "H_coefficients_low_to_high_V_1_delta": pair_polynomial_record(
            [h0, ONE]
        ),
        "p_adic_modulus_bits": int(modulus.nbits()),
        "p_adic_numerator_coordinates_replayed": 16,
        "p_adic_mismatches": p_adic_mismatches,
        "claim_boundary": (
            "The exact first-order deformation and all sixteen p-adic numerator "
            "coordinates replay. Full generic Q^2 divisibility remains a separate gate."
        ),
    }
generic_quartic_division = None
if args.attempt_generic_quartic_division:
    jet_artifact = json.loads(args.generic_quartic_jet_artifact.read_text())
    if jet_artifact.get("status") != "PASS_EXACT_GENERIC_QUARTIC_FIRST_JET_P19_REPLAY":
        raise ArithmeticError("generic quartic jet artifact is not certified")
    jet_inputs = jet_artifact["inputs"]
    for label, current_path in (
        ("operands", args.operands),
        ("pencil", args.pencil),
        ("factor_lift", args.factor_lift),
        ("specialized_factorization", args.specialized_factorization),
        ("H_candidate", args.H_candidate),
    ):
        if jet_inputs[label]["sha256"] != sha256(current_path):
            raise ArithmeticError(f"generic quartic jet has stale {label} input")
    generic_quartic_jet = jet_artifact["generic_quartic_first_jet"]
    q_hat = [
        pair_polynomial_from_record(coefficient)
        for coefficient in generic_quartic_jet[
            "Q_numerator_coefficients_low_to_high_W_then_V_1_delta"
        ]
    ]
    if len(q_hat) != 5 or q_hat[-1][-1] != ONE:
        raise ArithmeticError("saved generic quartic numerator has wrong shape")
    h0 = q_hat[-1][0]
    jet_denominator = mpz(1)
    for _, _, _, exact_value in parsed_terms:
        for coordinate in exact_value:
            jet_denominator = lcm(jet_denominator, coordinate.denominator)
    jet_integer_terms = []
    jet_content = mpz(0)
    for v_degree, w_degree, x_degree, exact_value in parsed_terms:
        integer_value = tuple(
            coordinate.numerator
            * (jet_denominator // coordinate.denominator)
            for coordinate in exact_value
        )
        jet_integer_terms.append((v_degree, w_degree, x_degree, integer_value))
        for coordinate in integer_value:
            jet_content = gcd(jet_content, abs(coordinate))
    jet_integer_terms = [
        (
            v_degree,
            w_degree,
            x_degree,
            tuple(coordinate // jet_content for coordinate in integer_value),
        )
        for v_degree, w_degree, x_degree, integer_value in jet_integer_terms
    ]
    print("Q80Q12GENERICDIV|stage=assemble_full_cubic", flush=True)
    cubic_bivariate = [[[ZERO]] for _ in range(4)]
    for v_degree, w_degree, x_degree, exact_value in jet_integer_terms:
        coefficient = cubic_bivariate[x_degree]
        while len(coefficient) <= w_degree:
            coefficient.append([ZERO])
        v_polynomial = coefficient[w_degree]
        while len(v_polynomial) <= v_degree:
            v_polynomial.append(ZERO)
        v_polynomial[v_degree] = k_add(v_polynomial[v_degree], exact_value)
        coefficient[w_degree] = p_trim(v_polynomial)
        cubic_bivariate[x_degree] = b_trim(coefficient)
    full_a, full_b, full_c, full_d = cubic_bivariate[3], cubic_bivariate[2], cubic_bivariate[1], cubic_bivariate[0]
    print("Q80Q12GENERICDIV|stage=build_full_discriminant", flush=True)
    full_discriminant = b_add(
        b_add(
            b_add(
                b_mul(
                    b_pow(full_b, 2, delta_square),
                    b_pow(full_c, 2, delta_square),
                    delta_square,
                ),
                b_neg(
                    b_scale(
                        b_mul(full_a, b_pow(full_c, 3, delta_square), delta_square),
                        4,
                    )
                ),
            ),
            b_neg(
                b_scale(
                    b_mul(b_pow(full_b, 3, delta_square), full_d, delta_square),
                    4,
                )
            ),
        ),
        b_add(
            b_neg(
                b_scale(
                    b_mul(
                        b_pow(full_a, 2, delta_square),
                        b_pow(full_d, 2, delta_square),
                        delta_square,
                    ),
                    27,
                )
            ),
            b_scale(
                b_mul(
                    b_mul(
                        b_mul(full_a, full_b, delta_square),
                        full_c,
                        delta_square,
                    ),
                    full_d,
                    delta_square,
                ),
                18,
            ),
        ),
    )
    print(
        "Q80Q12GENERICDIV|stage=strip_L3|"
        f"discriminant_W_degree={len(full_discriminant)-1}",
        flush=True,
    )
    full_residual = full_discriminant
    for multiplicity in range(3):
        full_residual, full_remainder = b_divide_monic_linear(
            full_residual, linear_constant, delta_square
        )
        if full_remainder != [ZERO]:
            raise ArithmeticError(
                f"full generic discriminant fails L division {multiplicity+1}"
            )
    H_polynomial = [h0, ONE]
    H_squared = p_mul(H_polynomial, H_polynomial, delta_square)
    scaled_full_residual = [
        p_mul(coefficient, H_squared, delta_square)
        for coefficient in full_residual
    ]
    q_hat_squared = b_mul(q_hat, q_hat, delta_square)
    print(
        "Q80Q12GENERICDIV|stage=pseudo_divide_Q2|"
        f"residual_W_degree={len(full_residual)-1}|"
        f"Q2_W_degree={len(q_hat_squared)-1}",
        flush=True,
    )
    pseudo_quotient, pseudo_remainder, pseudo_lead, pseudo_steps = b_pseudo_divmod(
        scaled_full_residual, q_hat_squared, delta_square
    )
    if pseudo_remainder != [[ZERO]]:
        raise ArithmeticError(
            "jet-derived generic quartic square fails full exact pseudo-division"
        )
    if len(pseudo_quotient) - 1 != 4:
        raise ArithmeticError("generic complementary factor is not quartic in W")
    denominator_H_exponent = 2 * pseudo_steps
    cancelled_H_exponent = 0
    simplified_quotient = pseudo_quotient
    while cancelled_H_exponent < denominator_H_exponent:
        divided_coefficients = []
        for coefficient in simplified_quotient:
            divided, remainder = p_divmod_field(
                coefficient, H_polynomial, delta_square
            )
            if remainder != [ZERO]:
                divided_coefficients = None
                break
            divided_coefficients.append(divided)
        if divided_coefficients is None:
            break
        simplified_quotient = b_trim(divided_coefficients)
        cancelled_H_exponent += 1
    denominator_H_exponent -= cancelled_H_exponent
    generic_quartic_division = {
        "identity": "H(V)^2*Delta/L^3 is divisible by Qhat(V,W)^2 in K(V)[W]",
        "full_discriminant_W_degree": len(full_discriminant) - 1,
        "full_discriminant_V_degree": max(len(value) - 1 for value in full_discriminant),
        "linear_stripped_W_degree": len(full_residual) - 1,
        "Qhat_squared_W_degree": len(q_hat_squared) - 1,
        "pseudo_division_steps": pseudo_steps,
        "pseudo_remainder_zero": True,
        "complementary_factor_W_degree": len(simplified_quotient) - 1,
        "complementary_factor_coefficients_low_to_high_W_then_V_1_delta": [
            pair_polynomial_record(coefficient)
            for coefficient in simplified_quotient
        ],
        "complementary_factor_common_denominator": (
            f"H(V)^{denominator_H_exponent}"
        ),
        "cancelled_H_exponent": cancelled_H_exponent,
        "complementary_factor_numerator_maximum_coordinate_bits": max(
            p_maximum_coordinate_bits(coefficient)
            for coefficient in simplified_quotient
        ),
        "conclusion": "exact generic quartic-square divisibility over the descent field",
        "claim_boundary": (
            "The generic discriminant factor Q^2 and a complementary W-quartic are exact. "
            "Jacobian invariants, minimization, and birational maps remain open."
        ),
    }
    degree_exponents = [[4, 1], [4, 2], [1, 3]]
if args.attempt_quartic_gcd:
    # Convert only the L-stripped degree-12 residual to Singular.  This gcd is
    # deliberately opt-in: both PARI and Singular exceed the bounded probe
    # window on the current million-bit coefficients.
    generator_ring = PolynomialRing(QQ, "z")
    z = generator_ring.gen()
    descent_field = NumberField(z**2 - QQ(delta_square), "delta")
    delta = descent_field.gen()
    w_ring = PolynomialRing(descent_field, "W", implementation="singular")
    W = w_ring.gen()

    def sage_field(value):
        return descent_field(QQ(value[0])) + descent_field(QQ(value[1])) * delta

    residual_sage = sum(
        (
            sage_field(value) * W**degree
            for degree, value in enumerate(residual_discriminant)
        ),
        w_ring.zero(),
    )
    L = W + descent_field(QQ(linear_rational))
    Q = residual_sage.gcd(residual_sage.derivative(W)).monic()
    D, remainder = residual_sage.monic().quo_rem(Q**2)
    if remainder:
        raise ArithmeticError("Q^2 does not divide the L-stripped monic discriminant")
    D = D.monic()
    factor_data = [(D, 1), (Q, 2), (L, 3)]
    degree_exponents = [[factor.degree(), exponent] for factor, exponent in factor_data]
    if degree_exponents != [[4, 1], [4, 2], [1, 3]]:
        raise ArithmeticError(
            f"unexpected exact specialized discriminant factors: {degree_exponents}"
        )
    if Q**2 * D != residual_sage.monic():
        raise ArithmeticError("exact squarefree-chain factorization does not replay")

pair_factorization = (
    subresultant_prs if subresultant_prs is not None else custom_prs
)
payload = {
    "schema": (
        "elkies-k3-q80-third-q12-exact-generic-quartic-factorization-v1"
        if generic_quartic_division is not None
        else (
            "elkies-k3-q80-third-q12-exact-generic-quartic-jet-v1"
            if generic_quartic_jet is not None
            else (
                "elkies-k3-q80-third-q12-exact-discriminant-specialization-v1"
                if factor_data is not None or pair_factorization is not None
                else (
                    "elkies-k3-q80-third-q12-exact-generic-linear-conductor-v1"
                    if generic_linear is not None
                    else "elkies-k3-q80-third-q12-exact-linear-conductor-specialization-v1"
                )
            )
        )
    ),
    "status": (
        "PASS_EXACT_GENERIC_L3_Q2_D_FACTORIZATION"
        if generic_quartic_division is not None
        else (
            "PASS_EXACT_GENERIC_QUARTIC_FIRST_JET_P19_REPLAY"
            if generic_quartic_jet is not None
            else (
                "PASS_EXACT_SPECIALIZED_L3_Q2_D_FACTORIZATION"
                if factor_data is not None or pair_factorization is not None
                else (
                    "PASS_EXACT_GENERIC_LINEAR_CONDUCTOR_MULTIPLICITY_THREE"
                    if generic_linear is not None
                    else "PASS_EXACT_SPECIALIZED_LINEAR_CONDUCTOR_MULTIPLICITY_THREE"
                )
            )
        )
    ),
    "base_value": rational_record(base_value),
    "quadratic_field": {
        "generator": "delta",
        "delta_square": str(delta_square),
        "omega_square": rational_record(omega_square),
        "omega_in_delta_basis": {
            "coefficient": rational_record(
                4 * product_numerator_root / product.denominator
            )
        },
    },
    "projective_scaling": {
        "common_denominator_bits": common_denominator.bit_length(),
        "integer_content_bits": content.bit_length(),
    },
    "discriminant_degree": len(discriminant) - 1,
    "discriminant_maximum_coordinate_bits": max(
        rational_height_bits(coordinate)
        for coefficient in discriminant
        for coordinate in coefficient
    ),
    "linear_stripped_residual_degree": len(residual_discriminant) - 1,
    "linear_stripped_residual_maximum_coordinate_bits": max(
        rational_height_bits(coordinate)
        for coefficient in residual_discriminant
        for coordinate in coefficient
    ),
    "linear_factor_reconstruction": {
        "constant": rational_record(mpq(linear_constant[0])),
        "numerator_bits": int(abs(linear_rational.numerator()).nbits()),
        "denominator_bits": int(linear_rational.denominator().nbits()),
        "modulus_bits": int(modulus.nbits()),
        "exact_multiplicity": 3,
    },
    "generic_linear_factor": generic_linear,
    "factor_degree_exponents_recovered": degree_exponents,
    "custom_primitive_remainder_sequence": custom_prs,
    "subresultant_remainder_sequence": subresultant_prs,
    "generic_quartic_first_jet": generic_quartic_jet,
    "generic_quartic_factorization": generic_quartic_division,
    "monic_factors_by_exponent": (
        {
            str(exponent): polynomial_record(factor)
            for factor, exponent in factor_data
        }
        if factor_data is not None
        else (
            {
                "1": pair_factorization["D_coefficients_low_to_high_W_1_delta"],
                "2": pair_factorization["Q_coefficients_low_to_high_W_1_delta"],
                "3": [pair_record(linear_constant), pair_record(ONE)],
            }
            if pair_factorization is not None
            else {"3": [pair_record(linear_constant), pair_record(ONE)]}
        )
    ),
    "inputs": {
        "operands": {
            "path": str(args.operands.relative_to(ROOT)),
            "sha256": sha256(args.operands),
        },
        "pencil": {
            "path": str(args.pencil.relative_to(ROOT)),
            "sha256": sha256(args.pencil),
        },
        "factor_lift": {
            "path": str(args.factor_lift.relative_to(ROOT)),
            "sha256": sha256(args.factor_lift),
        },
        "specialized_factorization": (
            {
                "path": str(args.specialized_factorization.relative_to(ROOT)),
                "sha256": sha256(args.specialized_factorization),
            }
            if generic_quartic_jet is not None
            else None
        ),
        "H_candidate": (
            {
                "path": str(args.H_candidate.relative_to(ROOT)),
                "sha256": sha256(args.H_candidate),
            }
            if generic_quartic_jet is not None
            else None
        ),
        "generic_quartic_jet_artifact": (
            {
                "path": str(args.generic_quartic_jet_artifact.relative_to(ROOT)),
                "sha256": sha256(args.generic_quartic_jet_artifact),
            }
            if generic_quartic_division is not None
            else None
        ),
        "worker": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
    },
    "claim_boundary": (
        "The generic L^3*Q^2 factorization and complementary W-quartic are exact "
        "over the characteristic-zero descent field. No Jacobian or birational map "
        "is claimed."
        if generic_quartic_division is not None
        else (
            (
                "The exact first-order V-adic quartic deformation is recovered and "
                "replays the p-adic numerator data; full generic Q^2 divisibility is open. "
                if generic_quartic_jet is not None
                else (
                    "The linear factor and multiplicity three are exact over the full "
                    "characteristic-zero base. "
                    if generic_linear is not None
                    else "The linear factor and multiplicity are exact only at the displayed base value. "
                )
            )
            + (
                "The quartic squarefree factors are also exact. "
                if factor_data is not None or pair_factorization is not None
                else "The quartic Q and D factors remain unrecovered. "
            )
            + "No generic factorization, Jacobian, or birational map is claimed."
        )
    ),
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_q80_third_q12_exact_discriminant_specialization.sage "
        f"--base-value {base_value}"
        + (" --attempt-quartic-gcd" if args.attempt_quartic_gcd else "")
        + (" --attempt-custom-prs" if args.attempt_custom_prs else "")
        + (" --attempt-subresultant-prs" if args.attempt_subresultant_prs else "")
        + (" --attempt-generic-quartic-jet" if args.attempt_generic_quartic_jet else "")
        + (
            " --attempt-generic-quartic-division"
            if args.attempt_generic_quartic_division
            else ""
        )
        + (
            f" --subresultant-checkpoint {args.subresultant_checkpoint.relative_to(ROOT)}"
            if args.attempt_subresultant_prs
            else ""
        )
        + (
            " --restart-subresultant-prs"
            if args.attempt_subresultant_prs and args.restart_subresultant_prs
            else ""
        )
        + (" --certify-generic-linear" if args.certify_generic_linear else "")
        + (
            f" --output {args.output.relative_to(ROOT)}"
            if args.output is not None
            else ""
        )
    ),
}
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
digest = hashlib.sha256(encoded.encode()).hexdigest()
if args.output:
    if args.check:
        if args.output.read_text() != encoded:
            raise SystemExit(f"stale exact specialization artifact: {args.output}")
        print(
            f"Q80Q12EXACTDISC|V={base_value}|recovered={degree_exponents}|"
            f"sha256={digest}|status=PASS_CHECK"
        )
    else:
        args.output.write_text(encoded)
        print(
            f"Q80Q12EXACTDISC|V={base_value}|recovered={degree_exponents}|"
            f"artifact={args.output}|sha256={digest}|status=PASS_WRITE"
        )
else:
    print(
        f"Q80Q12EXACTDISC|V={base_value}|recovered={degree_exponents}|"
        f"sha256={digest}|status=PASS"
    )
