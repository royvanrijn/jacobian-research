#!/usr/bin/env sage -python
"""Sample the normalized third-q12 quartic factor at good inert primes.

The exact cubic pencil has discriminant of the expected shape

    Delta = L^3 * Q^2 * D,       deg_W(Q)=deg_W(D)=4,

and the candidate normalization writes

    Q = W^4 + sum(N_i(V)/H(V) * W^i, i=0..3),
    deg_V(N_i) <= 1.

This worker uses univariate specializations in V.  A gcd is used only as a
finite-field extraction oracle.  Every extracted quartic is checked by direct
Q^2 division, its affine numerators are interpolated from two samples, and
the result is replayed at one or more held-out V-values.

The output is modular reconstruction input, not a characteristic-zero proof.
Primes 163, 191, and 199 are reserved for blind replay unless explicitly
enabled with ``--allow-heldout-prime``.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

from sage.all import GF, PolynomialRing, QQ, ZZ


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_PENCIL = RESULTS / "q80-third-q12-um2-biquadratic-resolved-pencil-qq.json"
DEFAULT_OPERANDS = (
    RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json"
)
DEFAULT_H = RESULTS / "elkies-k3-q80-third-q12-quartic-denominator-candidate-v1.json"
DEFAULT_LINEAR = (
    RESULTS / "elkies-k3-q80-third-q12-exact-generic-linear-conductor-v1.json"
)
RESERVED_HELDOUT_PRIMES = (163, 191, 199)
COEFFICIENT = re.compile(r"^(-?\d+)/(\d+)\*theta\^2 ([+-]) (\d+)/(\d+)$")


class BadPrime(Exception):
    """A requested prime does not support the pinned modular normalization."""


class BadSample(Exception):
    """One V-specialization is exceptional at an otherwise usable prime."""


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational(record):
    return QQ(ZZ(record["numerator"])) / ZZ(record["denominator"])


def relative_or_absolute(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def parse_integer_list(values, label):
    answer = []
    for value in values:
        for part in value.split(","):
            part = part.strip()
            if not part:
                continue
            try:
                answer.append(int(part))
            except ValueError as error:
                raise ValueError(f"invalid {label} value {part!r}") from error
    return answer


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--prime",
    action="append",
    required=True,
    help="prime, comma-separated primes, or a repeatable option",
)
parser.add_argument("--pencil", type=Path, default=DEFAULT_PENCIL)
parser.add_argument("--operands", type=Path, default=DEFAULT_OPERANDS)
parser.add_argument("--H-candidate", dest="H_candidate", type=Path, default=DEFAULT_H)
parser.add_argument("--linear-certificate", type=Path, default=DEFAULT_LINEAR)
parser.add_argument("--output-directory", type=Path, required=True)
parser.add_argument(
    "--training-V",
    action="append",
    default=None,
    help="two preferred training values (default: 0,1); repeatable or comma-separated",
)
parser.add_argument(
    "--heldout-V",
    action="append",
    default=None,
    help="preferred held-out values (default: 2); repeatable or comma-separated",
)
parser.add_argument("--allow-heldout-prime", action="store_true")
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
for key in ("pencil", "operands", "H_candidate", "linear_certificate", "output_directory"):
    setattr(args, key, getattr(args, key).resolve())

primes = parse_integer_list(args.prime, "prime")
if not primes:
    raise ValueError("at least one prime is required")
primes = list(dict.fromkeys(primes))
training_requested = parse_integer_list(args.training_V or ["0,1"], "training V")
heldout_requested = parse_integer_list(args.heldout_V or ["2"], "held-out V")
if len(training_requested) != 2 or training_requested[0] == training_requested[1]:
    raise ValueError("exactly two distinct preferred training V-values are required")
if not heldout_requested:
    raise ValueError("at least one preferred held-out V-value is required")

pencil = json.loads(args.pencil.read_text())
operands = json.loads(args.operands.read_text())
H_artifact = json.loads(args.H_candidate.read_text())
linear_artifact = json.loads(args.linear_certificate.read_text())
if pencil.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_RESOLVED_PENCIL":
    raise ValueError("exact resolved pencil is not certified")
if operands.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_CLOSURE_OPERANDS_P19_HENSEL":
    raise ValueError("exact closure operands are not certified")
if H_artifact.get("status") != "PASS_CANDIDATE_THREE_HELDOUTS_EXACT_Q_OPEN":
    raise ValueError("quartic denominator candidate is not certified as a candidate")
if linear_artifact.get("status") != "PASS_EXACT_GENERIC_LINEAR_CONDUCTOR_MULTIPLICITY_THREE":
    raise ValueError("exact generic linear conductor is not certified")

q1 = rational(operands["biquadratic_field"]["q1"])
q2 = rational(operands["biquadratic_field"]["q2"])
product = q1 * q2
product_root = ZZ(product.numerator()).isqrt()
if product_root**2 != product.numerator():
    raise ArithmeticError("q1*q2 numerator is not a square")
delta_square = ZZ(product.denominator())
if ZZ(H_artifact["candidate"]["delta_square"]) != delta_square:
    raise ArithmeticError("H candidate and closure operands use different delta fields")
if ZZ(linear_artifact["quadratic_field"]["delta_square"]) != delta_square:
    raise ArithmeticError("linear certificate and closure operands use different delta fields")

omega_to_delta = QQ(4 * product_root) / product.denominator()
if omega_to_delta**2 * delta_square != 16 * product:
    raise ArithmeticError("omega-to-delta basis conversion is inconsistent")
h0_rational = rational(H_artifact["candidate"]["h0_rational"])
h0_delta = rational(H_artifact["candidate"]["h0_delta"])
linear_constant = rational(linear_artifact["linear_factor_reconstruction"]["constant"])

# Parse the exact pencil only once, before looping over primes.  In the fixed
# delta basis theta^2=q1+q2+(omega_to_delta/2)*delta.
terms = []
for v_degree, w_degree, x_degree, encoded in pencil["moving_equation"][
    "terms_T_W_x_coefficient_1_r"
]:
    if len(encoded) != 1:
        raise ArithmeticError("unexpected exact pencil coefficient record")
    match = COEFFICIENT.fullmatch(encoded[0])
    if match is None:
        raise ArithmeticError("unexpected exact pencil coefficient encoding")
    theta2_coefficient = QQ(ZZ(match[1])) / ZZ(match[2])
    sign = 1 if match[3] == "+" else -1
    constant = sign * QQ(ZZ(match[4])) / ZZ(match[5])
    terms.append(
        (
            int(v_degree),
            int(w_degree),
            int(x_degree),
            constant + theta2_coefficient * (q1 + q2),
            theta2_coefficient * omega_to_delta / 2,
        )
    )
if len(terms) != 63:
    raise ArithmeticError("exact moving-equation support changed")

input_records = {
    "pencil": {"path": relative_or_absolute(args.pencil), "sha256": sha256(args.pencil)},
    "operands": {"path": relative_or_absolute(args.operands), "sha256": sha256(args.operands)},
    "H_candidate": {
        "path": relative_or_absolute(args.H_candidate),
        "sha256": sha256(args.H_candidate),
    },
    "linear_certificate": {
        "path": relative_or_absolute(args.linear_certificate),
        "sha256": sha256(args.linear_certificate),
    },
    "worker": {
        "path": relative_or_absolute(Path(__file__).resolve()),
        "sha256": sha256(Path(__file__).resolve()),
    },
}


def coordinate_pair(value, finite):
    coefficients = list(finite(value).polynomial())
    coefficients += [finite.base_ring().zero()] * (2 - len(coefficients))
    if len(coefficients) != 2:
        raise ArithmeticError("quadratic-field coefficient lost its two-coordinate basis")
    return [int(coefficients[0]), int(coefficients[1])]


def polynomial_pairs(polynomial, degree, finite):
    return [coordinate_pair(polynomial[index], finite) for index in range(degree + 1)]


def process_prime(prime):
    if prime in RESERVED_HELDOUT_PRIMES and not args.allow_heldout_prime:
        raise BadPrime("reserved blind-replay prime; pass --allow-heldout-prime to consume it")
    if prime <= 2 or not ZZ(prime).is_prime():
        raise BadPrime("requested modulus is not an odd prime")

    constants = GF(prime)

    def reduce_rational(value, label):
        value = QQ(value)
        if value.denominator() % prime == 0:
            raise BadPrime(f"prime divides exact denominator of {label}")
        return constants(value.numerator()) / constants(value.denominator())

    delta_square_mod = reduce_rational(QQ(delta_square), "delta_square")
    if not delta_square_mod or delta_square_mod.is_square():
        raise BadPrime("delta field is not inert and unramified")
    modulus_ring = PolynomialRing(constants, "z_delta")
    z_delta = modulus_ring.gen()
    finite = GF(
        prime**2,
        "delta",
        modulus=z_delta**2 - delta_square_mod,
    )
    delta = finite.gen()
    if delta**2 != finite(delta_square_mod):
        raise ArithmeticError("finite delta generator has the wrong square")

    h0 = finite(reduce_rational(h0_rational, "H rational constant")) + finite(
        reduce_rational(h0_delta, "H delta constant")
    ) * delta
    r_value = finite(reduce_rational(linear_constant, "linear conductor constant"))
    finite_terms = []
    for v_degree, w_degree, x_degree, rational_part, delta_part in terms:
        coefficient = finite(reduce_rational(rational_part, "pencil rational coordinate"))
        coefficient += finite(reduce_rational(delta_part, "pencil delta coordinate")) * delta
        finite_terms.append((v_degree, w_degree, x_degree, coefficient))

    w_ring = PolynomialRing(finite, "W")
    W = w_ring.gen()
    L = W + r_value
    rejected_samples = []
    cache = {}

    def extract(value):
        value = int(value) % prime
        if value in cache:
            cached = cache[value]
            if isinstance(cached, Exception):
                raise cached
            return cached
        V_value = constants(value)
        try:
            H_value = finite(V_value) + h0
            if not H_value:
                raise BadSample("H(V) vanishes")
            cubic_coefficients = [w_ring.zero() for unused in range(4)]
            for v_degree, w_degree, x_degree, coefficient in finite_terms:
                cubic_coefficients[x_degree] += (
                    coefficient * finite(V_value) ** v_degree * W**w_degree
                )
            leading = cubic_coefficients[3]
            if not leading or leading.degree() != 0:
                raise BadSample("cubic leading coefficient is zero or W-dependent")
            inverse_leading = leading[0] ** (-1)
            b = cubic_coefficients[2] * inverse_leading
            c = cubic_coefficients[1] * inverse_leading
            d = cubic_coefficients[0] * inverse_leading
            discriminant = b**2 * c**2 - 4 * c**3 - 4 * b**3 * d - 27 * d**2 + 18 * b * c * d
            if not discriminant or discriminant.degree() != 15:
                raise BadSample(
                    f"discriminant degree is {discriminant.degree() if discriminant else -1}, expected 15"
                )
            discriminant = discriminant.monic()
            residual, remainder = discriminant.quo_rem(L**3)
            if remainder:
                raise BadSample("exact L^3 does not divide the modular discriminant")
            if residual.degree() != 12:
                raise BadSample("L-stripped residual is not degree 12")
            unused_quotient, fourth_remainder = residual.quo_rem(L)
            if not fourth_remainder:
                raise BadSample("linear conductor has multiplicity above three")
            Q = residual.gcd(residual.derivative()).monic()
            if Q.degree() != 4:
                raise BadSample(f"square-factor gcd has degree {Q.degree()}, expected 4")
            D, square_remainder = residual.quo_rem(Q**2)
            if square_remainder or D.degree() != 4:
                raise BadSample("Q^2 division does not have a quartic quotient")
            D = D.monic()
            if Q.gcd(D).degree() != 0 or D.gcd(D.derivative()).degree() != 0:
                raise BadSample("quartic factors are not coprime and squarefree")
            if Q**2 * D != residual:
                raise ArithmeticError("direct modular factor multiplication failed")
            numerators = [H_value * Q[index] for index in range(4)]
            record = {
                "V": value,
                "H_1_delta": coordinate_pair(H_value, finite),
                "Q_coefficients_low_to_high_W_1_delta": polynomial_pairs(Q, 4, finite),
                "D_coefficients_low_to_high_W_1_delta": polynomial_pairs(D, 4, finite),
                "N_values_low_to_high_W_1_delta": [
                    coordinate_pair(numerator, finite) for numerator in numerators
                ],
                "checks": {
                    "discriminant_degree": 15,
                    "linear_multiplicity": 3,
                    "linear_stripped_degree": 12,
                    "gcd_degree": 4,
                    "Q_squared_division_quotient_degree": 4,
                    "Q_D_coprime": True,
                    "D_squarefree": True,
                    "direct_factor_multiplication": True,
                },
                "_Q": Q,
                "_N": numerators,
                "_H": H_value,
            }
        except BadSample as error:
            cache[value] = error
            rejected_samples.append({"V": value, "reason": str(error)})
            raise
        cache[value] = record
        return record

    fallback_values = [
        value
        for value in range(min(prime, 512))
        if value not in set(training_requested + heldout_requested)
    ]
    used = set()

    def choose(preferred):
        candidates = [preferred] + fallback_values
        for candidate in candidates:
            reduced = candidate % prime
            if reduced in used:
                continue
            try:
                record = extract(reduced)
            except BadSample:
                continue
            used.add(reduced)
            return record
        raise BadPrime("could not find enough valid V-specializations")

    training = [choose(value) for value in training_requested]
    heldout = [choose(value) for value in heldout_requested]
    v0 = finite(constants(training[0]["V"]))
    v1 = finite(constants(training[1]["V"]))
    if v0 == v1:
        raise ArithmeticError("training specializations collided modulo the prime")
    affine_numerators = []
    for w_degree in range(4):
        y0 = training[0]["_N"][w_degree]
        y1 = training[1]["_N"][w_degree]
        slope = (y1 - y0) / (v1 - v0)
        intercept = y0 - slope * v0
        affine_numerators.append((intercept, slope))

    for record in heldout:
        value = finite(constants(record["V"]))
        predicted = W**4
        for w_degree, (intercept, slope) in enumerate(affine_numerators):
            predicted += ((intercept + slope * value) / record["_H"]) * W**w_degree
        if predicted != record["_Q"]:
            raise BadPrime(f"affine numerator interpolation fails at held-out V={record['V']}")

    numerator_records = [
        [coordinate_pair(value, finite) for value in coefficient]
        for coefficient in affine_numerators
    ]
    residue_vector = [
        coordinate
        for coefficient in numerator_records
        for v_pair in coefficient
        for coordinate in v_pair
    ]
    if len(residue_vector) != 16:
        raise ArithmeticError("quartic numerator residue vector does not have 16 slots")
    labels = [
        f"N{w_degree}.V{v_degree}.{basis}"
        for w_degree in range(4)
        for v_degree in range(2)
        for basis in ("constant", "delta")
    ]

    def public_sample(record):
        return {key: value for key, value in record.items() if not key.startswith("_")}

    return {
        "prime": prime,
        "delta_square_mod_prime": int(delta_square_mod),
        "field_modulus": f"delta^2-{int(delta_square_mod)}",
        "linear_factor_constant": int(r_value),
        "H_constant_1_delta": coordinate_pair(h0, finite),
        "requested_training_V": [value % prime for value in training_requested],
        "requested_heldout_V": [value % prime for value in heldout_requested],
        "training_samples": [public_sample(record) for record in training],
        "heldout_samples": [public_sample(record) for record in heldout],
        "rejected_sample_diagnostics": rejected_samples,
        "interpolated_N_coefficients_low_to_high_W_then_V_1_delta": numerator_records,
        "residue_slot_labels": labels,
        "residue_vector": residue_vector,
        "checks": {
            "quadratic_field_inert_unramified": True,
            "training_sample_count": len(training),
            "heldout_sample_count": len(heldout),
            "all_samples_have_L3_Q2_D_shape": True,
            "all_heldout_affine_interpolation_replays": True,
            "slot_count": len(residue_vector),
        },
    }


def artifact_path(prime):
    return args.output_directory / f"q80-third-q12-quartic-modp-{prime}.json"


def reproduce_command(prime):
    command = (
        "sage -python elkies-k3/scripts/sample_q80_third_q12_quartic_modp.sage "
        f"--prime {prime} "
        f"--pencil {relative_or_absolute(args.pencil)} "
        f"--operands {relative_or_absolute(args.operands)} "
        f"--H-candidate {relative_or_absolute(args.H_candidate)} "
        f"--linear-certificate {relative_or_absolute(args.linear_certificate)} "
        f"--output-directory {relative_or_absolute(args.output_directory)} "
        f"--training-V {','.join(map(str, training_requested))} "
        f"--heldout-V {','.join(map(str, heldout_requested))}"
    )
    if prime in RESERVED_HELDOUT_PRIMES and args.allow_heldout_prime:
        command += " --allow-heldout-prime"
    return command


args.output_directory.mkdir(parents=True, exist_ok=True)
failure_count = 0
for prime in primes:
    try:
        modular = process_prime(prime)
        status = "PASS_Q80_THIRD_Q12_QUARTIC_MODP_INTERPOLATION"
        bad_prime_diagnostics = []
    except BadPrime as error:
        modular = None
        status = "REJECT_BAD_OR_RESERVED_PRIME"
        bad_prime_diagnostics = [str(error)]
        failure_count += 1
    payload = {
        "schema": "elkies-k3.q80-third-q12-quartic-modp-sample.v1",
        "status": status,
        "requested_prime": prime,
        "coefficient_field_basis": ["1", "delta"],
        "delta_square": str(delta_square),
        "normalization": {
            "quartic": "Q=W^4+sum(N_i(V)/H(V)*W^i,i=0..3)",
            "H": "V+h0_rational+h0_delta*delta",
            "N_degree_V_bound": 1,
            "residue_order": "W degree, then V degree, then (constant,delta)",
        },
        "modular_result": modular,
        "bad_prime_diagnostics": bad_prime_diagnostics,
        "reserved_blind_primes": list(RESERVED_HELDOUT_PRIMES),
        "reserved_prime_consumed": bool(
            prime in RESERVED_HELDOUT_PRIMES and args.allow_heldout_prime
        ),
        "claim_boundary": (
            "This is a modular extraction and interpolation record. The gcd is only a "
            "finite-field oracle. Characteristic-zero H, Q, D, and Q^2 divisibility "
            "remain open until coefficientwise CRT reconstruction passes direct exact division."
        ),
        "inputs": input_records,
        "reproduce": reproduce_command(prime),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = artifact_path(prime)
    digest = hashlib.sha256(encoded.encode()).hexdigest()
    if args.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit(f"stale modular quartic sample artifact: {output}")
        action = "PASS_CHECK"
    else:
        temporary = output.with_suffix(output.suffix + ".tmp")
        temporary.write_text(encoded)
        temporary.replace(output)
        action = "PASS_WRITE"
    print(
        f"Q80Q12QMODP|prime={prime}|status={status}|artifact={output}|"
        f"sha256={digest}|action={action}",
        flush=True,
    )

if failure_count:
    print(
        f"Q80Q12QMODP|requested={len(primes)}|rejected={failure_count}|"
        "status=COMPLETED_WITH_REJECTED_PRIMES",
        file=sys.stderr,
    )
