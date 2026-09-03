#!/usr/bin/env sage -python
"""Reconstruct and directly certify the generic third-q12 quartic factor.

status: ACTIVE_COMPILER
claim: coefficientwise CRT reconstruction followed by exact Q^2 division
inputs: exact pencil/descent data, the p=19 power lift, H, and modular samples
outputs: optional exact quartic reconstruction certificate

Each modular sample supplies the four degree-at-most-one numerators of

    Q = W^4 + sum(N_i(V)/H(V) * W^i, i=0..3)

in the fixed basis ``(1,delta)``.  The p=19^12288 lift is always included as
one CRT modulus.  Primes 163, 191, and 199 are permanently withheld: they are
rebuilt from the exact characteristic-zero pencil and used only as blind
checks.  Independently, the last modular inputs are withheld from CRT and
used to reject scalar rational reconstructions that are not yet stable.

No characteristic-zero gcd or factorization is used.  Once all sixteen
rational coordinates reconstruct, the worker rebuilds the homogeneous cubic
discriminant in exact pair arithmetic, removes the already known L^3 by
synthetic division, clears H by writing P=H*Q, divides exactly by P twice,
and finally checks the coefficientwise multiplication identity.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys

from sage.all import CRT_list, GF, PolynomialRing, QQ, ZZ, inverse_mod


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
DEFAULT_PENCIL = RESULTS / "q80-third-q12-um2-biquadratic-resolved-pencil-qq.json"
DEFAULT_OPERANDS = (
    RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json"
)
DEFAULT_FACTOR_LIFT = (
    LOCAL / "q80-third-q12-discriminant-factors-p19-adic-precision12288.json"
)
DEFAULT_H = RESULTS / "elkies-k3-q80-third-q12-quartic-denominator-candidate-v1.json"
DEFAULT_LINEAR = RESULTS / "elkies-k3-q80-third-q12-exact-generic-linear-conductor-v1.json"
DEFAULT_EXACT_SPECIALIZATION = (
    RESULTS / "elkies-k3-q80-third-q12-exact-discriminant-specialization-v1.json"
)
DEFAULT_EXACT_JET = (
    LOCAL / "q80-third-q12-exact-generic-quartic-jet-v1.json"
)
DEFAULT_OUTPUT = RESULTS / "elkies-k3-q80-third-q12-quartic-crt-qq-v1.json"
BLIND_PRIMES = (163, 191, 199)
COEFFICIENT = re.compile(r"^(-?\d+)/(\d+)\*theta\^2 ([+-]) (\d+)/(\d+)$")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative_or_absolute(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def rational(record):
    return QQ(ZZ(record["numerator"])) / ZZ(record["denominator"])


def rational_record(value):
    value = QQ(value)
    return {"numerator": str(value.numerator()), "denominator": str(value.denominator())}


def reduce_rational(value, modulus):
    value = QQ(value)
    denominator = ZZ(value.denominator())
    if math.gcd(int(denominator), int(modulus)) != 1:
        raise ZeroDivisionError(f"denominator is not invertible modulo {modulus}")
    return int(ZZ(value.numerator()) * inverse_mod(denominator, ZZ(modulus)) % modulus)


def collect_sha256_values(value):
    answer = set()
    if isinstance(value, dict):
        for key, item in value.items():
            if key == "sha256" and isinstance(item, str):
                answer.add(item)
            answer.update(collect_sha256_values(item))
    elif isinstance(value, list):
        for item in value:
            answer.update(collect_sha256_values(item))
    return answer


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("inputs", nargs="*", type=Path, help="modular quartic numerator samples")
parser.add_argument("--pencil", type=Path, default=DEFAULT_PENCIL)
parser.add_argument("--operands", type=Path, default=DEFAULT_OPERANDS)
parser.add_argument("--factor-lift", type=Path, default=DEFAULT_FACTOR_LIFT)
parser.add_argument("--H-candidate", type=Path, default=DEFAULT_H)
parser.add_argument("--linear-certificate", type=Path, default=DEFAULT_LINEAR)
parser.add_argument(
    "--exact-specialization", type=Path, default=DEFAULT_EXACT_SPECIALIZATION
)
parser.add_argument(
    "--exact-jet",
    type=Path,
    help="exact first-order jet certificate; bypasses CRT reconstruction when supplied",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
parser.add_argument("--write-artifact", action="store_true")
parser.add_argument(
    "--stability-heldout-count",
    type=int,
    default=8,
    help="reserve the last N modular inputs from CRT for scalar stability checks",
)
parser.add_argument(
    "--diagnostic-incomplete",
    action="store_true",
    help="print unresolved CRT coordinates instead of failing; never promotes an artifact",
)
args = parser.parse_args()
if args.check and args.write_artifact:
    parser.error("--check and --write-artifact are mutually exclusive")
if args.stability_heldout_count < 0:
    parser.error("--stability-heldout-count must be nonnegative")
for name in (
    "pencil",
    "operands",
    "factor_lift",
    "H_candidate",
    "linear_certificate",
    "exact_specialization",
    "exact_jet",
    "output",
):
    path = getattr(args, name)
    if path is None:
        continue
    if not path.is_absolute():
        path = ROOT / path
    setattr(args, name, path.resolve())
args.inputs = [(path if path.is_absolute() else ROOT / path).resolve() for path in args.inputs]

pencil = json.loads(args.pencil.read_text())
operands = json.loads(args.operands.read_text())
factor_lift = json.loads(args.factor_lift.read_text())
H_artifact = json.loads(args.H_candidate.read_text())
linear_artifact = json.loads(args.linear_certificate.read_text())
exact_specialization = json.loads(args.exact_specialization.read_text())

q1 = rational(operands["biquadratic_field"]["q1"])
q2 = rational(operands["biquadratic_field"]["q2"])
product = q1 * q2
product_root = ZZ(product.numerator()).isqrt()
if product_root**2 != product.numerator():
    raise ArithmeticError("q1*q2 numerator is not a square")
delta_square = ZZ(product.denominator())
omega_to_delta = QQ(4 * product_root) / delta_square

h0_rational = rational(H_artifact["candidate"]["h0_rational"])
h0_delta = rational(H_artifact["candidate"]["h0_delta"])
expected_hashes = {
    sha256(args.pencil),
    sha256(args.operands),
    sha256(args.H_candidate),
    sha256(args.linear_certificate),
}

if exact_specialization.get("schema") != (
    "elkies-k3-q80-third-q12-exact-discriminant-specialization-v1"
):
    raise ArithmeticError("unexpected exact V=0 specialization schema")
if exact_specialization.get("status") != "PASS_EXACT_SPECIALIZED_L3_Q2_D_FACTORIZATION":
    raise ArithmeticError("exact V=0 specialization is not certified")
if rational(exact_specialization["base_value"]) != 0:
    raise ArithmeticError("exact specialization anchor is not at V=0")
if ZZ(exact_specialization["quadratic_field"]["delta_square"]) != delta_square:
    raise ArithmeticError("exact specialization uses a different delta field")
exact_pair_factorization = exact_specialization.get("subresultant_remainder_sequence")
if not isinstance(exact_pair_factorization, dict):
    exact_pair_factorization = exact_specialization.get("custom_primitive_remainder_sequence")
if not isinstance(exact_pair_factorization, dict):
    raise ArithmeticError("exact specialization has no pair-arithmetic factorization")
exact_Q_at_zero = exact_pair_factorization["Q_coefficients_low_to_high_W_1_delta"]
if len(exact_Q_at_zero) != 5:
    raise ArithmeticError("exact specialized exponent-two factor is not quartic")


def exact_pair(record):
    if not isinstance(record, list) or len(record) != 2:
        raise ArithmeticError("exact quadratic coefficient is not a pair")
    return rational(record[0]), rational(record[1])


def multiply_exact_pairs(left, right):
    return (
        left[0] * right[0] + delta_square * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


if exact_pair(exact_Q_at_zero[-1]) != (QQ.one(), QQ.zero()):
    raise ArithmeticError("exact specialized quartic is not monic")
H_at_zero = (h0_rational, h0_delta)
exact_intercepts = [
    multiply_exact_pairs(exact_pair(record), H_at_zero)
    for record in exact_Q_at_zero[:4]
]


def parse_basis(record, path):
    basis = record.get(
        "coefficient_field_basis",
        record.get("coefficient_basis", record.get("basis")),
    )
    declared_square = record.get("delta_square")
    label = basis
    if isinstance(basis, dict):
        declared_square = basis.get("delta_square", declared_square)
        label = basis.get("name", basis.get("label", basis.get("coordinates")))
    if declared_square is None and isinstance(record.get("field"), dict):
        declared_square = record["field"].get("delta_square")
    normalized_label = re.sub(r"[^a-z0-9]", "", str(label).lower())
    if "delta" not in normalized_label or "1" not in normalized_label:
        raise ArithmeticError(f"{path}: coefficient basis is not explicitly (1,delta)")
    modular = record.get("modular_result")
    if declared_square is None and isinstance(modular, dict):
        declared_square = modular.get("delta_square", delta_square)
    if declared_square is None or ZZ(declared_square) != delta_square:
        raise ArithmeticError(f"{path}: delta-square declaration changed")


def parse_coordinates(record, modulus, path):
    modular = record.get("modular_result", {})
    raw = modular.get("interpolated_N_coefficients_low_to_high_W_then_V_1_delta")
    if not isinstance(raw, list) or len(raw) != 4:
        raise ArithmeticError(f"{path}: expected four W numerator blocks")
    result = []
    for w_degree, v_rows in enumerate(raw):
        if not isinstance(v_rows, list) or len(v_rows) != 2:
            raise ArithmeticError(f"{path}: W^{w_degree} block is not 2x2")
        result.append([])
        for v_degree, pair in enumerate(v_rows):
            if not isinstance(pair, list) or len(pair) != 2:
                raise ArithmeticError(
                    f"{path}: W^{w_degree} V^{v_degree} coordinate is not a pair"
                )
            result[-1].append([int(ZZ(value) % modulus) for value in pair])
    return result


sample_records = []
seen_primes = set()
for path in args.inputs:
    record = json.loads(path.read_text())
    if record.get("schema") != "elkies-k3.q80-third-q12-quartic-modp-sample.v1":
        raise ArithmeticError(f"{path}: unexpected modular sample schema")
    modular_result = record.get("modular_result")
    if not isinstance(modular_result, dict):
        raise ArithmeticError(f"{path}: modular sample has no accepted result")
    prime = ZZ(modular_result["prime"])
    if not prime.is_prime():
        raise ArithmeticError(f"{path}: sample modulus is not prime")
    if int(prime) in BLIND_PRIMES:
        raise ArithmeticError(f"{path}: prime {prime} is reserved for blind replay")
    if int(prime) == 19:
        raise ArithmeticError(f"{path}: p=19 is already supplied by the power lift")
    if int(prime) in seen_primes:
        raise ArithmeticError(f"duplicate modular sample prime {prime}")
    seen_primes.add(int(prime))
    if record.get("status") != "PASS_Q80_THIRD_Q12_QUARTIC_MODP_INTERPOLATION":
        raise ArithmeticError(f"{path}: modular sample does not have PASS status")
    parse_basis(record, path)
    if record.get("coefficient_field_basis") != ["1", "delta"]:
        raise ArithmeticError(f"{path}: modular sample is not in the fixed delta basis")
    if ZZ(modular_result.get("delta_square_mod_prime", -1)) != delta_square % prime:
        raise ArithmeticError(f"{path}: modular delta square does not match exact field")
    supplied_hashes = collect_sha256_values(record.get("inputs", {}))
    if not expected_hashes.issubset(supplied_hashes):
        missing = sorted(expected_hashes - supplied_hashes)
        raise ArithmeticError(f"{path}: missing exact input hashes {missing}")
    if GF(prime)(delta_square).is_square():
        raise ArithmeticError(f"{path}: delta field splits at prime {prime}")
    coordinates = parse_coordinates(record, prime, path)
    flattened = [
        coordinate
        for coefficient in coordinates
        for v_pair in coefficient
        for coordinate in v_pair
    ]
    if flattened != [int(ZZ(value) % prime) for value in modular_result["residue_vector"]]:
        raise ArithmeticError(f"{path}: nested coordinates and residue vector disagree")
    sample_records.append(
        {"path": path, "record": record, "prime": prime, "coordinates": coordinates}
    )

if sample_records and args.stability_heldout_count == 0:
    parser.error("at least one stability-heldout modular sample is required")
stability_count = min(args.stability_heldout_count, len(sample_records))
if stability_count:
    training_records = sample_records[:-stability_count]
    stability_records = sample_records[-stability_count:]
else:
    training_records = sample_records
    stability_records = []

# The p-adic lift is the first, high-precision CRT image.  Its coefficients
# are stored in the omega basis and are converted to the fixed delta basis.
p19_modulus = ZZ(factor_lift["specialization"]["modulus"])
if ZZ(factor_lift["specialization"]["prime"]) != 19:
    raise ArithmeticError("factor lift is not the pinned p=19 lift")
omega_to_delta_modulus = (
    ZZ(omega_to_delta.numerator())
    * inverse_mod(ZZ(omega_to_delta.denominator()), p19_modulus)
) % p19_modulus
q_records = factor_lift["factorization"]["Q"]["coefficients_low_to_high_W"]
if len(q_records) != 5 or q_records[-1]["degrees_numerator_denominator"] != [0, 0]:
    raise ArithmeticError("unexpected p=19 monic quartic schema")
p19_coordinates = []
for w_degree, record in enumerate(q_records[:4]):
    if record["degrees_numerator_denominator"] != [1, 1]:
        raise ArithmeticError(f"p=19 Q coefficient W^{w_degree} lost degree (1,1)")
    numerator = record["numerator_coefficients_low_to_high_U_1_omega"]
    if len(numerator) != 2 or any(len(pair) != 2 for pair in numerator):
        raise ArithmeticError("p=19 numerator lost its 2x2 coordinate shape")
    p19_coordinates.append(
        [
            [
                int(ZZ(pair[0]) % p19_modulus),
                int(ZZ(pair[1]) * omega_to_delta_modulus % p19_modulus),
            ]
            for pair in numerator
        ]
    )

moduli = [p19_modulus] + [item["prime"] for item in training_records]
crt_modulus = math.prod(moduli)
coordinate_results = []
candidate_coordinates = [[[None, None] for _ in range(2)] for _ in range(4)]
stability_rejection_count = 0
for w_degree in range(4):
    for v_degree in range(2):
        for field_coordinate in range(2):
            if v_degree == 0:
                anchored = exact_intercepts[w_degree][field_coordinate]
                expected_images = [
                    (p19_modulus, p19_coordinates[w_degree][0][field_coordinate])
                ] + [
                    (
                        item["prime"],
                        item["coordinates"][w_degree][0][field_coordinate],
                    )
                    for item in sample_records
                ]
                for image_modulus, expected in expected_images:
                    if reduce_rational(anchored, image_modulus) != expected:
                        raise ArithmeticError(
                            "exact V=0 anchor fails modular replay at "
                            f"modulus={image_modulus}, slot={(w_degree, 0, field_coordinate)}"
                        )
                candidate_coordinates[w_degree][0][field_coordinate] = anchored
                coordinate_results.append(
                    {
                        "index_W_V_field": [w_degree, 0, field_coordinate],
                        "resolved": True,
                        "value": rational_record(anchored),
                        "source": "exact_V0_Q_times_H0",
                        "all_modular_reductions_replayed": True,
                    }
                )
                continue
            residues = [p19_coordinates[w_degree][v_degree][field_coordinate]] + [
                item["coordinates"][w_degree][v_degree][field_coordinate]
                for item in training_records
            ]
            residue = ZZ(CRT_list(residues, moduli))
            try:
                reconstructed = residue.rational_reconstruction(ZZ(crt_modulus))
                failed_stability_primes = []
                for item in stability_records:
                    try:
                        actual = reduce_rational(reconstructed, item["prime"])
                    except ZeroDivisionError:
                        actual = None
                    expected = item["coordinates"][w_degree][v_degree][field_coordinate]
                    if actual != expected:
                        failed_stability_primes.append(int(item["prime"]))
                if failed_stability_primes:
                    stability_rejection_count += 1
                    result = {
                        "index_W_V_field": [w_degree, v_degree, field_coordinate],
                        "resolved": False,
                        "reason": "rational reconstruction fails stability-heldout replay",
                        "failed_stability_primes": failed_stability_primes,
                    }
                else:
                    candidate_coordinates[w_degree][v_degree][field_coordinate] = reconstructed
                    result = {
                        "index_W_V_field": [w_degree, v_degree, field_coordinate],
                        "resolved": True,
                        "value": rational_record(reconstructed),
                        "stability_heldout_replay": True,
                    }
            except ArithmeticError:
                result = {
                    "index_W_V_field": [w_degree, v_degree, field_coordinate],
                    "resolved": False,
                    "reason": "no rational reconstruction at current CRT modulus",
                }
            coordinate_results.append(result)

unresolved = [item for item in coordinate_results if not item["resolved"]]

exact_jet = None
if args.exact_jet is not None:
    exact_jet = json.loads(args.exact_jet.read_text())
    if exact_jet.get("status") != "PASS_EXACT_GENERIC_QUARTIC_FIRST_JET_P19_REPLAY":
        raise ArithmeticError("exact jet artifact does not have certified PASS status")
    exact_jet_inputs = exact_jet.get("inputs", {})
    for label, current_path in (
        ("pencil", args.pencil),
        ("operands", args.operands),
        ("factor_lift", args.factor_lift),
        ("specialized_factorization", args.exact_specialization),
        ("H_candidate", args.H_candidate),
    ):
        if exact_jet_inputs.get(label, {}).get("sha256") != sha256(current_path):
            raise ArithmeticError(f"exact jet has stale {label} input")
    jet_data = exact_jet.get("generic_quartic_first_jet", {})
    raw_coordinates = jet_data.get(
        "Q_numerator_coefficients_low_to_high_W_then_V_1_delta"
    )
    if not isinstance(raw_coordinates, list) or len(raw_coordinates) != 5:
        raise ArithmeticError("exact jet quartic numerator has wrong W shape")
    candidate_coordinates = []
    for w_degree, v_rows in enumerate(raw_coordinates[:4]):
        if not isinstance(v_rows, list) or len(v_rows) != 2:
            raise ArithmeticError(f"exact jet W^{w_degree} numerator is not linear")
        candidate_coordinates.append([])
        for pair in v_rows:
            if not isinstance(pair, list) or len(pair) != 2:
                raise ArithmeticError("exact jet coefficient is not a quadratic pair")
            candidate_coordinates[-1].append([rational(value) for value in pair])
    leading = raw_coordinates[4]
    if (
        len(leading) != 2
        or [rational(value) for value in leading[0]] != [h0_rational, h0_delta]
        or [rational(value) for value in leading[1]] != [QQ.one(), QQ.zero()]
    ):
        raise ArithmeticError("exact jet leading coefficient does not equal H(V)")
    coordinate_results = [
        {
            "index_W_V_field": [w_degree, v_degree, field_coordinate],
            "resolved": True,
            "value": rational_record(
                candidate_coordinates[w_degree][v_degree][field_coordinate]
            ),
            "source": "exact_first_order_jet",
        }
        for w_degree in range(4)
        for v_degree in range(2)
        for field_coordinate in range(2)
    ]
    unresolved = []


def verify_sample_reductions(coordinates):
    for item in sample_records:
        prime = item["prime"]
        for w_degree in range(4):
            for v_degree in range(2):
                for field_coordinate in range(2):
                    actual = reduce_rational(
                        coordinates[w_degree][v_degree][field_coordinate], prime
                    )
                    expected = item["coordinates"][w_degree][v_degree][field_coordinate]
                    if actual != expected:
                        raise ArithmeticError(
                            f"candidate fails CRT-input replay at p={prime}, "
                            f"slot={(w_degree, v_degree, field_coordinate)}"
                        )
    for w_degree in range(4):
        for v_degree in range(2):
            for field_coordinate in range(2):
                actual = reduce_rational(
                    coordinates[w_degree][v_degree][field_coordinate], p19_modulus
                )
                expected = p19_coordinates[w_degree][v_degree][field_coordinate]
                if actual != expected:
                    raise ArithmeticError(
                        "candidate fails p=19-power replay at "
                        f"slot={(w_degree, v_degree, field_coordinate)}"
                    )


# Parse the exact pencil only once.  These values are used independently in
# the finite blind replays and in the characteristic-zero identity below.
parsed_terms = []
for v_degree, w_degree, x_degree, encoded in pencil["moving_equation"][
    "terms_T_W_x_coefficient_1_r"
]:
    match = COEFFICIENT.fullmatch(encoded[0])
    if match is None:
        raise ArithmeticError("unexpected exact pencil coefficient encoding")
    theta2 = QQ(ZZ(match[1])) / ZZ(match[2])
    sign = 1 if match[3] == "+" else -1
    constant = sign * QQ(ZZ(match[4])) / ZZ(match[5])
    parsed_terms.append(
        (
            int(v_degree),
            int(w_degree),
            int(x_degree),
            (
                constant + theta2 * (q1 + q2),
                QQ(2) * theta2 * product_root / delta_square,
            ),
        )
    )


def blind_replays(coordinates):
    records = []
    for prime_integer in BLIND_PRIMES:
        prime = ZZ(prime_integer)
        constants = GF(prime)
        if constants(delta_square).is_square():
            raise ArithmeticError(f"delta field splits at blind prime {prime}")
        z_ring = PolynomialRing(constants, "z")
        z = z_ring.gen()
        finite = GF(prime**2, "delta", modulus=z**2 - constants(delta_square))
        delta = finite.gen()
        v_ring = PolynomialRing(finite, "V")
        V = v_ring.gen()
        v_field = v_ring.fraction_field()
        w_ring = PolynomialRing(v_field, "W")
        W = w_ring.gen()
        x_ring = PolynomialRing(w_ring, "old_x")
        cubic_coefficients = [w_ring.zero() for _ in range(4)]
        for v_degree, w_degree, x_degree, pair in parsed_terms:
            coefficient = finite(pair[0]) + finite(pair[1]) * delta
            cubic_coefficients[x_degree] += (
                coefficient * v_field(V) ** v_degree * W**w_degree
            )
        cubic = x_ring(
            [value / cubic_coefficients[3] for value in cubic_coefficients]
        )
        b, c, d = cubic[2], cubic[1], cubic[0]
        discriminant = b**2 * c**2 - 4 * c**3 - 4 * b**3 * d - 27 * d**2 + 18 * b * c * d
        factors = discriminant.factor()
        shape = sorted((int(factor.degree()), int(exponent)) for factor, exponent in factors)
        if shape != [(1, 3), (4, 1), (4, 2)]:
            raise ArithmeticError(f"blind factor shape changed at p={prime}: {shape}")
        finite_Q = next(
            factor.monic() for factor, exponent in factors if int(exponent) == 2
        )
        finite_H = V + finite(h0_rational) + finite(h0_delta) * delta
        expected_coordinates = []
        for w_degree in range(4):
            numerator = finite_Q[w_degree] * finite_H
            if numerator.denominator().degree() != 0:
                raise ArithmeticError(f"blind Q*H is not polynomial at p={prime}")
            numerator = numerator.numerator() / numerator.denominator()[0]
            coefficients = list(numerator) + [finite.zero()] * (2 - len(numerator.list()))
            if len(coefficients) != 2:
                raise ArithmeticError(f"blind Q*H changed V degree at p={prime}")
            block = []
            for coefficient in coefficients:
                pair = list(coefficient.polynomial()) + [constants.zero()] * 2
                block.append([int(pair[0]), int(pair[1])])
            expected_coordinates.append(block)
        for w_degree in range(4):
            for v_degree in range(2):
                for field_coordinate in range(2):
                    actual = reduce_rational(
                        coordinates[w_degree][v_degree][field_coordinate], prime
                    )
                    expected = expected_coordinates[w_degree][v_degree][field_coordinate]
                    if actual != expected:
                        raise ArithmeticError(
                            f"candidate fails blind replay at p={prime}, "
                            f"slot={(w_degree, v_degree, field_coordinate)}"
                        )
        records.append(
            {
                "prime": int(prime),
                "factor_degree_exponents": [list(value) for value in shape],
                "all_16_numerator_coordinates_replayed": True,
            }
        )
    return records


# Exact K[V,W] arithmetic.  A K[V] coefficient is a pair (a,b), denoting
# a+b*delta.  A K[V,W] polynomial is a low-to-high W list of such pairs.
RV = PolynomialRing(QQ, "V")
V_exact = RV.gen()
KZERO = (RV.zero(), RV.zero())
KONE = (RV.one(), RV.zero())


def k_add(left, right):
    return left[0] + right[0], left[1] + right[1]


def k_neg(value):
    return -value[0], -value[1]


def k_mul(left, right):
    return (
        left[0] * right[0] + delta_square * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def k_scale(value, scalar):
    return RV(scalar) * value[0], RV(scalar) * value[1]


def w_trim(value):
    result = list(value)
    while result and result[-1] == KZERO:
        result.pop()
    return result or [KZERO]


def w_add(left, right):
    answer = [KZERO] * max(len(left), len(right))
    for index in range(len(answer)):
        answer[index] = k_add(
            left[index] if index < len(left) else KZERO,
            right[index] if index < len(right) else KZERO,
        )
    return w_trim(answer)


def w_neg(value):
    return w_trim([k_neg(coefficient) for coefficient in value])


def w_scale(value, scalar):
    return w_trim([k_scale(coefficient, scalar) for coefficient in value])


def w_mul(left, right):
    answer = [KZERO] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] = k_add(
                answer[left_index + right_index], k_mul(left_value, right_value)
            )
    return w_trim(answer)


def w_pow(value, exponent):
    answer = [KONE]
    power = value
    while exponent:
        if exponent & 1:
            answer = w_mul(answer, power)
        exponent //= 2
        if exponent:
            power = w_mul(power, power)
    return answer


def divide_monic_linear(value, constant):
    value = w_trim(value)
    if len(value) < 2:
        return [KZERO], value[0]
    quotient = [KZERO] * (len(value) - 1)
    quotient[-1] = value[-1]
    for degree in range(len(value) - 2, 0, -1):
        quotient[degree - 1] = k_add(
            value[degree], k_neg(k_mul(constant, quotient[degree]))
        )
    remainder = k_add(value[0], k_neg(k_mul(constant, quotient[0])))
    return w_trim(quotient), remainder


def exact_divide_by_H(value, H):
    """Divide one K[V] pair by H, using conjugation and no polynomial gcd."""
    conjugate = (H[0], -H[1])
    norm_pair = k_mul(H, conjugate)
    if norm_pair[1] != 0 or norm_pair[0] == 0:
        raise ArithmeticError("H has invalid rational norm")
    numerator = k_mul(value, conjugate)
    quotients = []
    for coordinate in numerator:
        quotient, remainder = coordinate.quo_rem(norm_pair[0])
        if remainder:
            raise ArithmeticError("leading coefficient is not exactly divisible by H")
        quotients.append(quotient)
    result = tuple(quotients)
    if k_mul(result, H) != value:
        raise ArithmeticError("exact H division failed multiplication replay")
    return result


def exact_divide_by_P(dividend, divisor, H):
    dividend = w_trim(dividend)
    divisor = w_trim(divisor)
    if divisor[-1] != H:
        raise ArithmeticError("cleared quartic leading coefficient is not H")
    quotient = [KZERO] * max(1, len(dividend) - len(divisor) + 1)
    remainder = list(dividend)
    while remainder != [KZERO] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = exact_divide_by_H(remainder[-1], H)
        quotient[shift] = coefficient
        for index, divisor_coefficient in enumerate(divisor):
            position = shift + index
            remainder[position] = k_add(
                remainder[position], k_neg(k_mul(coefficient, divisor_coefficient))
            )
        remainder = w_trim(remainder)
    quotient = w_trim(quotient)
    if w_add(w_mul(divisor, quotient), remainder) != dividend:
        raise ArithmeticError("P division failed multiplication replay")
    return quotient, remainder


def polynomial_digest(value):
    digest = hashlib.sha256()
    for w_degree, coefficient in enumerate(w_trim(value)):
        for field_coordinate, polynomial in enumerate(coefficient):
            for v_degree, scalar in enumerate(polynomial.list()):
                digest.update(
                    (
                        f"{w_degree}:{field_coordinate}:{v_degree}:"
                        f"{scalar.numerator()}/{scalar.denominator()}\n"
                    ).encode()
                )
    return digest.hexdigest()


def maximum_bits(value):
    return max(
        max(abs(ZZ(scalar.numerator())).nbits(), ZZ(scalar.denominator()).nbits())
        for coefficient in w_trim(value)
        for polynomial in coefficient
        for scalar in polynomial.list()
    )


def exact_division_certificate(coordinates):
    cubic_coefficients = [[KZERO] for _ in range(4)]
    for v_degree, w_degree, x_degree, pair in parsed_terms:
        while len(cubic_coefficients[x_degree]) <= w_degree:
            cubic_coefficients[x_degree].append(KZERO)
        contribution = (
            RV(pair[0]) * V_exact**v_degree,
            RV(pair[1]) * V_exact**v_degree,
        )
        cubic_coefficients[x_degree][w_degree] = k_add(
            cubic_coefficients[x_degree][w_degree], contribution
        )
    cubic_coefficients = list(map(w_trim, cubic_coefficients))
    a, b, c, d = (
        cubic_coefficients[3],
        cubic_coefficients[2],
        cubic_coefficients[1],
        cubic_coefficients[0],
    )
    discriminant = w_add(
        w_add(
            w_add(w_mul(w_pow(b, 2), w_pow(c, 2)), w_neg(w_scale(w_mul(a, w_pow(c, 3)), 4))),
            w_neg(w_scale(w_mul(w_pow(b, 3), d), 4)),
        ),
        w_add(
            w_neg(w_scale(w_mul(w_pow(a, 2), w_pow(d, 2)), 27)),
            w_scale(w_mul(w_mul(w_mul(a, b), c), d), 18),
        ),
    )

    p19_linear = factor_lift["factorization"]["L"]["coefficients_low_to_high_W"][0]
    numerator = p19_linear["numerator_coefficients_low_to_high_U_1_omega"][0]
    denominator = p19_linear["denominator_coefficients_low_to_high_U_1_omega"][0]
    if numerator[1] or denominator[1]:
        raise ArithmeticError("known L constant is not rational")
    linear_residue = ZZ(numerator[0]) * inverse_mod(ZZ(denominator[0]), p19_modulus) % p19_modulus
    linear_rational = rational(
        linear_artifact["linear_factor_reconstruction"]["constant"]
    )
    if reduce_rational(linear_rational, p19_modulus) != linear_residue:
        raise ArithmeticError("exact L certificate and p=19 lift disagree")
    linear_constant = (RV(linear_rational), RV.zero())

    residual = discriminant
    for multiplicity in range(3):
        residual, remainder = divide_monic_linear(residual, linear_constant)
        if remainder != KZERO:
            raise ArithmeticError(f"known L fails exact division {multiplicity + 1}")
    unused, fourth_remainder = divide_monic_linear(residual, linear_constant)
    if fourth_remainder == KZERO:
        raise ArithmeticError("known L has multiplicity above three")
    if len(residual) - 1 != 12:
        raise ArithmeticError(f"L-stripped discriminant has degree {len(residual)-1}, not 12")

    H = (V_exact + RV(h0_rational), RV(h0_delta))
    if H[0].degree() != 1 or H[0].leading_coefficient() != 1:
        raise ArithmeticError("H is not monic linear")
    numerators = []
    for w_degree in range(4):
        numerators.append(
            (
                RV(
                    [
                        coordinates[w_degree][v_degree][0]
                        for v_degree in range(2)
                    ]
                ),
                RV(
                    [
                        coordinates[w_degree][v_degree][1]
                        for v_degree in range(2)
                    ]
                ),
            )
        )
    conjugate = (H[0], -H[1])
    norm_H = k_mul(H, conjugate)[0]
    numerator_remainders = []
    for numerator in numerators:
        multiplied = k_mul(numerator, conjugate)
        numerator_remainders.append(
            [coordinate.quo_rem(norm_H)[1] for coordinate in multiplied]
        )
    if not any(any(remainder for remainder in pair) for pair in numerator_remainders):
        raise ArithmeticError("all Q numerators are divisible by H; normalization is imprimitive")

    P = numerators + [H]
    # Q=P/H.  Therefore Q^2|R is certified integrally by the cleared
    # identity H^2*R=P^2*D.  Dividing R itself by P would incorrectly demand
    # that its W-leading coefficient already contain H twice.
    H_squared = k_mul(H, H)
    cleared_residual = w_mul(residual, [H_squared])
    first_quotient, first_remainder = exact_divide_by_P(cleared_residual, P, H)
    if first_remainder != [KZERO] or len(first_quotient) - 1 != 8:
        raise ArithmeticError("first cleared-quartic division did not give degree eight")
    second_quotient, second_remainder = exact_divide_by_P(first_quotient, P, H)
    if second_remainder != [KZERO] or len(second_quotient) - 1 != 4:
        raise ArithmeticError("second cleared-quartic division did not give degree four")
    if w_mul(w_mul(P, P), second_quotient) != cleared_residual:
        raise ArithmeticError("direct H^2*R=P^2*D multiplication identity failed")

    return {
        "field": f"QQ(delta), delta^2={delta_square}",
        "known_linear_factor": f"W+({linear_rational})",
        "known_linear_exact_multiplicity": 3,
        "degree_sequence_W": [
            len(discriminant) - 1,
            len(residual) - 1,
            len(first_quotient) - 1,
            len(second_quotient) - 1,
        ],
        "division_remainders_zero": [True, True, True, True, True],
        "cleared_quartic_primitive": True,
        "clearing_identity": "Q=P/H, so H^2*R=P^2*D",
        "direct_multiplication_identity": "H^2*R=P^2*D, with P=H*Q",
        "direct_multiplication_verified": True,
        "maximum_coordinate_bits": {
            "L_stripped_discriminant": maximum_bits(residual),
            "H2_cleared_L_stripped_discriminant": maximum_bits(cleared_residual),
            "cleared_quartic": maximum_bits(P),
            "remaining_degree_four_factor": maximum_bits(second_quotient),
        },
        "sha256": {
            "L_stripped_discriminant": polynomial_digest(residual),
            "H2_cleared_L_stripped_discriminant": polynomial_digest(cleared_residual),
            "cleared_quartic": polynomial_digest(P),
            "remaining_degree_four_factor": polynomial_digest(second_quotient),
        },
        "conclusion": "Q^2 divides the exact L^3-stripped discriminant in QQ(delta)(V)[W]",
    }


coordinates_resolved = not unresolved
blind_records = None
exact_certificate = None
candidate_rejection = None
if coordinates_resolved:
    verify_sample_reductions(candidate_coordinates)
    try:
        blind_records = blind_replays(candidate_coordinates)
        exact_certificate = exact_division_certificate(candidate_coordinates)
    except (ArithmeticError, ZeroDivisionError) as error:
        candidate_rejection = {
            "stage": "blind_replay_or_exact_direct_division",
            "reason": str(error),
        }

complete = coordinates_resolved and candidate_rejection is None

status = (
    (
        "PASS_EXACT_Q80_THIRD_Q12_QUARTIC_JET_AND_DIRECT_DIVISION"
        if exact_jet is not None
        else "PASS_EXACT_Q80_THIRD_Q12_QUARTIC_CRT_AND_DIRECT_DIVISION"
    )
    if complete
    else "INCOMPLETE_Q80_THIRD_Q12_QUARTIC_CRT_RECONSTRUCTION"
)
payload = {
    "schema": (
        "elkies-k3-q80-third-q12-quartic-jet-direct-division-v1"
        if exact_jet is not None
        else "elkies-k3-q80-third-q12-quartic-crt-qq-v1"
    ),
    "status": status,
    "reconstruction": {
        "normalization": "Q=W^4+sum(N_i(V)/H(V)*W^i); H monic linear",
        "coordinate_source": (
            "exact_first_order_jet" if exact_jet is not None else "CRT"
        ),
        "coordinate_order": "low_to_high_W_then_V_then_(1,delta)",
        "coordinate_count": 16,
        "exact_specialization_anchor_coordinate_count": 8,
        "CRT_slope_coordinate_count": 8,
        "resolved_coordinate_count": 16 - len(unresolved),
        "unresolved_coordinate_count": len(unresolved),
        "crt_primes": [19] + [int(item["prime"]) for item in training_records],
        "all_modular_sample_primes": [int(item["prime"]) for item in sample_records],
        "training_primes": [19] + [int(item["prime"]) for item in training_records],
        "stability_heldout_primes": [int(item["prime"]) for item in stability_records],
        "stability_heldout_requested_count": args.stability_heldout_count,
        "stability_rejection_count": stability_rejection_count,
        "p19_power_digits": int(factor_lift["specialization"]["digits"]),
        "crt_modulus_bits": int(ZZ(crt_modulus).nbits()),
        "coordinate_results": coordinate_results,
        "all_CRT_input_reductions_replayed": coordinates_resolved,
        "candidate_rejection": candidate_rejection,
    },
    "blind_exact_pencil_replays": blind_records,
    "exact_direct_division_certificate": exact_certificate,
    "claim_boundary": (
        "This certifies only the reconstructed exponent-two quartic factor by "
        "exact division and multiplication. The remaining factor, Jacobian, minimal "
        "model, maps, and downstream neighbour are not certified here."
        if complete
        else (
            "No characteristic-zero quartic is claimed: the CRT coordinates are unresolved "
            "or the candidate failed blind replay/exact direct division."
        )
    ),
    "inputs": {
        "pencil": {"path": str(args.pencil.relative_to(ROOT)), "sha256": sha256(args.pencil)},
        "operands": {"path": str(args.operands.relative_to(ROOT)), "sha256": sha256(args.operands)},
        "factor_lift": {
            "path": str(args.factor_lift.relative_to(ROOT)),
            "sha256": sha256(args.factor_lift),
        },
        "H_candidate": {
            "path": str(args.H_candidate.relative_to(ROOT)),
            "sha256": sha256(args.H_candidate),
        },
        "linear_certificate": {
            "path": str(args.linear_certificate.relative_to(ROOT)),
            "sha256": sha256(args.linear_certificate),
        },
        "exact_specialization": {
            "path": str(args.exact_specialization.relative_to(ROOT)),
            "sha256": sha256(args.exact_specialization),
        },
        "exact_jet": (
            {
                "path": relative_or_absolute(args.exact_jet),
                "sha256": sha256(args.exact_jet),
            }
            if args.exact_jet is not None
            else None
        ),
        "modular_samples": [
            {"path": relative_or_absolute(item["path"]), "sha256": sha256(item["path"])}
            for item in sample_records
        ],
        "worker": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/reconstruct_and_certify_q80_third_q12_quartic_crt_qq.sage "
        + " ".join(relative_or_absolute(path) for path in args.inputs)
        + (
            f" --exact-jet {relative_or_absolute(args.exact_jet)}"
            if args.exact_jet is not None
            else ""
        )
        + f" --stability-heldout-count {args.stability_heldout_count}"
        + (" --check" if complete else " --diagnostic-incomplete")
    ).strip(),
}
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
digest = hashlib.sha256(encoded.encode()).hexdigest()

if not complete and not args.diagnostic_incomplete:
    detail = (
        f"{len(unresolved)}/16 coordinates unresolved"
        if unresolved
        else f"candidate rejected: {candidate_rejection['reason']}"
    )
    raise SystemExit(
        f"incomplete CRT reconstruction: {detail}; "
        "rerun with more modular samples or --diagnostic-incomplete"
    )
if not complete and (args.check or args.write_artifact):
    raise SystemExit("an incomplete diagnostic cannot be checked or written as a certificate")
if args.write_artifact:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded)
    print(f"Q80Q12QCRT|artifact={args.output}|sha256={digest}|status=PASS_WRITE")
elif args.check:
    if args.output.read_text() != encoded:
        raise SystemExit(f"stale quartic CRT certificate: {args.output}")
    print(f"Q80Q12QCRT|artifact={args.output}|sha256={digest}|status=PASS_CHECK")
else:
    print(encoded, end="")
