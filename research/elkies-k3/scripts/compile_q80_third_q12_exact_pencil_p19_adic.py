#!/usr/bin/env python3
"""Compile the exact QQ(omega) q12 pencil modulo a high power of 19."""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_PENCIL = RESULTS / "q80-third-q12-um2-biquadratic-resolved-pencil-qq.json"
DEFAULT_OPERANDS = RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json"
DEFAULT_OUTPUT = RESULTS / "q80-third-q12-exact-pencil-p19-adic-precision64.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pencil", type=Path, default=DEFAULT_PENCIL)
parser.add_argument("--operands", type=Path, default=DEFAULT_OPERANDS)
parser.add_argument("--precision", type=int, default=64)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
args.pencil = args.pencil.resolve()
args.operands = args.operands.resolve()
args.output = args.output.resolve()
if args.precision < 2:
    raise ValueError("precision must be at least two p-adic digits")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


pencil = json.loads(args.pencil.read_text())
operands = json.loads(args.operands.read_text())
if pencil.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_RESOLVED_PENCIL":
    raise ValueError("exact resolved pencil is not certified")
if operands.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_CLOSURE_OPERANDS_P19_HENSEL":
    raise ValueError("exact closure operands are not certified")

prime = 19
modulus = prime**args.precision


def rational_parts(value):
    if "/" in value:
        numerator, denominator = value.split("/")
        return int(numerator), int(denominator)
    return int(value), 1


def rational_record_parts(record):
    return int(record["numerator"]), int(record["denominator"])


def rational_mod(parts):
    numerator, denominator = parts
    denominator_mod = denominator % modulus
    if denominator_mod % prime == 0:
        raise ZeroDivisionError("exact rational has bad p=19 reduction")
    return numerator % modulus * pow(denominator_mod, -1, modulus) % modulus


q1 = rational_mod(rational_record_parts(operands["biquadratic_field"]["q1"]))
q2 = rational_mod(rational_record_parts(operands["biquadratic_field"]["q2"]))
q_sum = (q1 + q2) % modulus
omega_square = 16 * q1 * q2 % modulus
if pow(omega_square % prime, (prime - 1) // 2, prime) != prime - 1:
    raise ArithmeticError("omega^2 is not an unramified quadratic unit at p=19")
inverse_two = pow(2, -1, modulus)

terms = []
for t_degree, w_degree, x_degree, encoded in pencil["moving_equation"][
    "terms_T_W_x_coefficient_1_r"
]:
    if len(encoded) != 1:
        raise ArithmeticError("exact coefficient is not one theta expression")
    match = re.fullmatch(r"(.+)\*theta\^2 ([+-]) (.+)", encoded[0])
    if match is None:
        raise ArithmeticError("exact coefficient does not descend to QQ(theta^2)")
    alpha = rational_mod(rational_parts(match.group(1)))
    beta_parts = rational_parts(match.group(3))
    if match.group(2) == "-":
        beta_parts = (-beta_parts[0], beta_parts[1])
    beta = rational_mod(beta_parts)
    # theta^2=q1+q2+omega/2.
    constant = (beta + alpha * q_sum) % modulus
    omega = alpha * inverse_two % modulus
    terms.append([t_degree, w_degree, x_degree, [constant, omega]])
if len(terms) != 63:
    raise ArithmeticError("exact moving-equation support changed")

output = {
    "schema": "elkies-k3.q80-third-q12-exact-pencil-p19-adic.v1",
    "status": "PASS_EXACT_THIRD_Q12_PENCIL_REDUCTION_MOD_19_POWER",
    "specialization": {
        "u": "-2",
        "prime": prime,
        "precision_digits": args.precision,
        "modulus": modulus,
    },
    "quadratic_field": {
        "basis": ["1", "omega"],
        "omega": "4*a*b",
        "omega_square_modulus": omega_square,
        "omega_square_mod_19": omega_square % prime,
        "unramified": True,
    },
    "pencil": {
        "variables": ["V", "W", "old_x"],
        "degrees_V_W_old_x": [2, 9, 3],
        "terms_V_W_old_x_coefficient_1_omega": terms,
        "coefficient_count": len(terms),
    },
    "inputs": {
        "pencil": {"path": str(args.pencil.relative_to(ROOT)), "sha256": sha256(args.pencil)},
        "operands": {"path": str(args.operands.relative_to(ROOT)), "sha256": sha256(args.operands)},
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "literal reduction of every exact QQ(omega) pencil coefficient modulo 19^precision",
            "all exact rational denominators are p=19 units",
            "the quadratic coefficient field is unramified at p=19",
        ],
        "not_proved": [
            "a p-adic Jacobian or birational maps",
            "Hensel nonsingularity of a coupled map-identity system",
            "a characteristic-zero child reconstruction",
        ],
    },
    "reproduce": (
        "python3 elkies-k3/scripts/compile_q80_third_q12_exact_pencil_p19_adic.py "
        f"--precision {args.precision} --output {args.output.relative_to(ROOT)}"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"p-adic exact-pencil artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    f"Q80THIRDQ12PADICPENCIL|p=19|precision={args.precision}|terms=63|"
    "field=Q19(omega)|status=PASS_EXACT_THIRD_Q12_PENCIL_REDUCTION_MOD_19_POWER"
)
