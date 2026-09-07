#!/usr/bin/env sage -python
"""Audit the common denominator candidate of the third-q12 quartic factor.

status: ACTIVE_COMPILER
claim: p-adic reconstruction plus three untouched-prime replays for H(V)
inputs: exact pencil, exact descent field, and the 19^12288 factor lift
outputs: elkies-k3-q80-third-q12-quartic-denominator-candidate-v1.json

This deliberately does not promote H or Q to characteristic zero: the exact
Q^2 divisibility test remains open.  It isolates and preserves the strongest
validated normalization for the next modular/subresultant implementation.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

from sage.all import GF, PolynomialRing, QQ, ZZ, inverse_mod


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_PENCIL = RESULTS / "q80-third-q12-um2-biquadratic-resolved-pencil-qq.json"
DEFAULT_OPERANDS = (
    RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json"
)
DEFAULT_LIFT = (
    ROOT
    / "artifacts/local/elkies-k3/"
    / "q80-third-q12-discriminant-factors-p19-adic-precision12288.json"
)
DEFAULT_OUTPUT = (
    RESULTS / "elkies-k3-q80-third-q12-quartic-denominator-candidate-v1.json"
)
COEFFICIENT = re.compile(
    r"^(-?\d+)/(\d+)\*theta\^2 ([+-]) (\d+)/(\d+)$"
)
HELD_OUT_PRIMES = (163, 191, 199)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rational(record):
    return QQ(ZZ(record["numerator"])) / ZZ(record["denominator"])


def rational_record(value):
    value = QQ(value)
    return {"numerator": str(value.numerator()), "denominator": str(value.denominator())}


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pencil", type=Path, default=DEFAULT_PENCIL)
parser.add_argument("--operands", type=Path, default=DEFAULT_OPERANDS)
parser.add_argument("--factor-lift", type=Path, default=DEFAULT_LIFT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
parser.add_argument("--write-artifact", action="store_true")
args = parser.parse_args()
for name in ("pencil", "operands", "factor_lift", "output"):
    setattr(args, name, getattr(args, name).resolve())

pencil = json.loads(args.pencil.read_text())
operands = json.loads(args.operands.read_text())
lift = json.loads(args.factor_lift.read_text())
q1 = read_rational(operands["biquadratic_field"]["q1"])
q2 = read_rational(operands["biquadratic_field"]["q2"])
product = q1 * q2
product_root = product.numerator().sqrt()
if product_root**2 != product.numerator():
    raise ArithmeticError("q1*q2 numerator is not a square")

modulus = ZZ(lift["specialization"]["modulus"])
q_record = lift["factorization"]["Q"]["coefficients_low_to_high_W"]
common_denominators = [
    record["denominator_coefficients_low_to_high_U_1_omega"]
    for record in q_record[:4]
]
if not all(value == common_denominators[0] for value in common_denominators[1:]):
    raise ArithmeticError("p-adic Q coefficients do not have one common denominator")
H_record = common_denominators[0]
if len(H_record) != 2 or H_record[1] != [1, 0]:
    raise ArithmeticError("p-adic common denominator is not monic linear in V")
h0_omega_residues = [ZZ(value) % modulus for value in H_record[0]]
omega_to_delta = QQ(4 * product_root) / product.denominator()
omega_to_delta_modulus = (
    ZZ(omega_to_delta.numerator())
    * inverse_mod(ZZ(omega_to_delta.denominator()), modulus)
) % modulus
h0_rational = h0_omega_residues[0].rational_reconstruction(modulus)
h0_delta = (
    h0_omega_residues[1] * omega_to_delta_modulus % modulus
).rational_reconstruction(modulus)
h0_omega = h0_delta / omega_to_delta

terms = []
for v_degree, w_degree, x_degree, encoded in pencil["moving_equation"][
    "terms_T_W_x_coefficient_1_r"
]:
    match = COEFFICIENT.fullmatch(encoded[0])
    if match is None:
        raise ArithmeticError("unexpected exact pencil coefficient encoding")
    theta2 = QQ(ZZ(match[1])) / ZZ(match[2])
    sign = 1 if match[3] == "+" else -1
    constant = sign * QQ(ZZ(match[4])) / ZZ(match[5])
    terms.append(
        (
            v_degree,
            w_degree,
            x_degree,
            constant + theta2 * (q1 + q2),
            theta2 / 2,
        )
    )

prime_replays = []
for prime in HELD_OUT_PRIMES:
    constants = GF(prime)
    omega_square = constants(16 * q1 * q2)
    if omega_square.is_square():
        raise ArithmeticError(f"descent field splits at held-out prime {prime}")
    modulus_ring = PolynomialRing(constants, "z")
    z = modulus_ring.gen()
    finite = GF(prime**2, "omega", modulus=z**2 - omega_square)
    omega = finite.gen()
    v_ring = PolynomialRing(finite, "V")
    V = v_ring.gen()
    v_field = v_ring.fraction_field()
    w_ring = PolynomialRing(v_field, "W")
    W = w_ring.gen()
    x_ring = PolynomialRing(w_ring, "old_x")
    coefficients = [w_ring.zero() for _ in range(4)]
    for v_degree, w_degree, x_degree, rational_part, omega_part in terms:
        coefficient = finite(rational_part) + finite(omega_part) * omega
        coefficients[x_degree] += coefficient * v_field(V)**v_degree * W**w_degree
    cubic = x_ring([value / coefficients[3] for value in coefficients])
    b, c, d = cubic[2], cubic[1], cubic[0]
    discriminant = b**2*c**2 - 4*c**3 - 4*b**3*d - 27*d**2 + 18*b*c*d
    factors = discriminant.factor()
    degree_exponents = sorted((int(factor.degree()), int(exponent)) for factor, exponent in factors)
    if degree_exponents != [(1, 3), (4, 1), (4, 2)]:
        raise ArithmeticError(f"held-out factor shape changed at p={prime}")
    Q = next(factor.monic() for factor, exponent in factors if int(exponent) == 2)
    denominators = [coefficient.denominator().monic() for coefficient in Q.list()[:4]]
    if not all(value == denominators[0] for value in denominators[1:]):
        raise ArithmeticError(f"finite Q denominators disagree at p={prime}")
    expected = V + finite(h0_rational) + finite(h0_omega) * omega
    if denominators[0] != expected:
        raise ArithmeticError(f"reconstructed H fails exact-pencil replay at p={prime}")
    prime_replays.append(
        {
            "prime": prime,
            "factor_degree_exponents": [list(value) for value in degree_exponents],
            "H_constant_1_omega": [int(finite(h0_rational)), int(finite(h0_omega))],
            "all_four_Q_denominators_equal": True,
            "candidate_replay": True,
        }
    )

payload = {
    "schema": "elkies-k3-q80-third-q12-quartic-denominator-candidate-v1",
    "status": "PASS_CANDIDATE_THREE_HELDOUTS_EXACT_Q_OPEN",
    "candidate": {
        "formula": "H(V)=V+h0_rational+h0_delta*delta",
        "delta_square": str(product.denominator()),
        "h0_rational": rational_record(h0_rational),
        "h0_delta": rational_record(h0_delta),
        "h0_omega": rational_record(h0_omega),
        "coordinate_height_bits": {
            "rational": [int(h0_rational.numerator().nbits()), int(h0_rational.denominator().nbits())],
            "delta": [int(h0_delta.numerator().nbits()), int(h0_delta.denominator().nbits())],
            "omega": [int(h0_omega.numerator().nbits()), int(h0_omega.denominator().nbits())],
        },
        "p_adic_modulus_bits": int(modulus.nbits()),
        "common_to_all_four_nonleading_Q_coefficients_mod_19_power": True,
    },
    "held_out_exact_pencil_replays": prime_replays,
    "claim_boundary": (
        "H is a p-adically reconstructed candidate validated at three untouched "
        "finite primes. Exact characteristic-zero Q^2 divisibility, Q, D, and the "
        "Jacobian remain open."
    ),
    "inputs": {
        "pencil": {"path": str(args.pencil.relative_to(ROOT)), "sha256": sha256(args.pencil)},
        "operands": {"path": str(args.operands.relative_to(ROOT)), "sha256": sha256(args.operands)},
        "factor_lift": {
            "path": str(args.factor_lift.relative_to(ROOT)),
            "sha256": sha256(args.factor_lift),
        },
        "worker": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/audit_q80_third_q12_quartic_denominator_candidate.sage --check"
    ),
}
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
digest = hashlib.sha256(encoded.encode()).hexdigest()
if args.write_artifact:
    args.output.write_text(encoded)
    print(f"Q80Q12QDEN|artifact={args.output}|sha256={digest}|status=PASS_WRITE")
elif args.check:
    if args.output.read_text() != encoded:
        raise SystemExit(f"stale quartic-denominator candidate artifact: {args.output}")
    print(f"Q80Q12QDEN|artifact={args.output}|sha256={digest}|status=PASS_CHECK")
else:
    print(encoded, end="")
