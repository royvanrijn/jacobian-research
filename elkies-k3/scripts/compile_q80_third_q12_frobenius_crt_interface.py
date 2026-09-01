#!/usr/bin/env python3
"""Compile the branch-mixed local q12 slots into a diagnostic CRT ledger."""

import argparse
import hashlib
import json
import sys
from pathlib import Path


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--prime-input", action="append", default=[], metavar="PRIME=PATH",
    help=(
        "certified Frobenius-invariant input; repeat to replace the default "
        "residue-prime set"
    ),
)
parser.add_argument(
    "--p19", type=Path, default=None,
    help="legacy override for the default p=19 input",
)
parser.add_argument(
    "--p61", type=Path, default=None,
    help="legacy override for the default p=61 input",
)
parser.add_argument(
    "--p67", type=Path, default=None,
    help="legacy override for the default p=67 input",
)
parser.add_argument(
    "--exact-operands", type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json",
)
parser.add_argument(
    "--held-out-audit", type=Path,
    default=RESULTS / "q80-fixed-u-minus2-p71-heldout-good-reduction.json",
)
parser.add_argument(
    "--output", type=Path,
    default=RESULTS / "q80-third-q12-um2-frobenius-crt-interface.json",
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
for name in ("p19", "p61", "p67"):
    if getattr(args, name) is not None:
        setattr(args, name, getattr(args, name).resolve())
for name in ("exact_operands", "held_out_audit", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


default_inputs = {
    prime: RESULTS / f"q80-third-q12-p{prime}-frobenius-invariants.json"
    for prime in (19, 61, 67, 83, 89, 103, 131)
}
for prime in (19, 61, 67):
    override = getattr(args, f"p{prime}")
    if override is not None:
        default_inputs[prime] = override


def parse_prime_input(spec):
    try:
        prime_text, path_text = spec.split("=", 1)
        prime = int(prime_text)
    except (TypeError, ValueError) as error:
        raise ValueError(f"invalid --prime-input {spec!r}; expected PRIME=PATH") from error
    if prime < 3 or any(prime % divisor == 0 for divisor in range(2, int(prime**0.5) + 1)):
        raise ValueError(f"--prime-input modulus is not an odd prime: {prime}")
    return prime, Path(path_text).resolve()


if args.prime_input:
    inputs = {}
    for spec in args.prime_input:
        prime, path = parse_prime_input(spec)
        if prime in inputs:
            raise ValueError(f"duplicate reconstruction prime: {prime}")
        inputs[prime] = path
else:
    inputs = {prime: path.resolve() for prime, path in default_inputs.items()}
primes = tuple(sorted(inputs))
if 19 not in inputs or len(primes) < 3:
    raise ValueError("the residue set must include p=19 and at least two controls")

payloads = {prime: json.loads(path.read_text()) for prime, path in inputs.items()}
for prime, payload in payloads.items():
    allowed_status = (
        "PASS_EXACT_FROBENIUS_INVARIANT_THIRD_Q12_ENCODING_MOD19_QUADRATIC"
        if prime == 19 else
        "PASS_EXACT_FROBENIUS_INVARIANT_THIRD_Q12_ENCODING_COMMON_PRODUCER"
    )
    if payload.get("status") != allowed_status:
        raise ValueError(f"p={prime}: uncertified Frobenius input")
    if int(payload["specialization"]["prime"]) != prime:
        raise ValueError(f"p={prime}: specialization mismatch")

held_out = json.loads(args.held_out_audit.read_text())
exact_operands = json.loads(args.exact_operands.read_text())
if exact_operands.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_CLOSURE_OPERANDS_P19_HENSEL":
    raise ValueError("exact biquadratic operand certificate is not valid")
parameters = held_out.get("parameters", [])
if len(parameters) != 1 or parameters[0].get("u") != "-2":
    raise ValueError("held-out audit is not the u=-2 specialization")
modular = parameters[0].get("modular", [])
if len(modular) != 1 or modular[0].get("prime") != 71 or modular[0].get("status") != "PASS_GOOD_REDUCTION_AUDIT":
    raise ValueError("p=71 is not certified as a held-out good prime")


def rational_mod(record, prime):
    numerator = int(record["numerator"]) % prime
    denominator = int(record["denominator"]) % prime
    if denominator == 0:
        raise ArithmeticError(f"operand denominator vanishes at p={prime}")
    return numerator * pow(denominator, -1, prime) % prime


def quadratic_character(value, prime):
    result = pow(value % prime, (prime - 1) // 2, prime)
    return -1 if result == prime - 1 else result


q1_record = exact_operands["biquadratic_field"]["q1"]
q2_record = exact_operands["biquadratic_field"]["q2"]
branch_diagnostic = {}
for prime in primes:
    q1_mod = rational_mod(q1_record, prime)
    q2_mod = rational_mod(q2_record, prime)
    characters = [
        quadratic_character(q1_mod, prime),
        quadratic_character(q2_mod, prime),
        quadratic_character(q1_mod * q2_mod, prime),
    ]
    if characters[2] != -1 or sorted(characters[:2]) != [-1, 1]:
        raise ArithmeticError(
            f"p={prime}: expected exactly one inert operand and inert product, got {characters}"
        )
    branch_diagnostic[str(prime)] = {
        "characters_q1_q2_q1q2": characters,
        "local_quadratic_operand": "q1" if characters[0] == -1 else "q2",
    }


def flatten_coefficients(value, path=()):
    if isinstance(value, dict):
        coefficient_keys = {
            "trace", "anti_invariant_coefficient", "coefficient_discriminant", "norm"
        }
        if set(value) == coefficient_keys:
            return {path: value}
        result = {}
        for key in sorted(value):
            result.update(flatten_coefficients(value[key], path + (key,)))
        return result
    if isinstance(value, list):
        result = {}
        for index, item in enumerate(value):
            result.update(flatten_coefficients(item, path + (index,)))
        return result
    return {}


flat = {
    prime: flatten_coefficients(payload["encoded_coefficients"])
    for prime, payload in payloads.items()
}
paths = sorted(flat[19], key=lambda path: tuple(map(str, path)))
if len(paths) != 1947 or any(set(flat[prime]) != set(paths) for prime in primes):
    raise ArithmeticError("ordered coefficient slots do not align across reconstruction primes")


def crt(residues):
    value = 0
    modulus = 1
    for prime in primes:
        target = int(residues[prime]) % prime
        correction = (target - value) * pow(modulus, -1, prime) % prime
        value += modulus * correction
        modulus *= prime
    if any(value % prime != int(residues[prime]) % prime for prime in residues):
        raise ArithmeticError("CRT replay failed")
    centered = value if value <= modulus // 2 else value - modulus
    return value, centered, modulus


slot_records = []
common_modulus = 1
for prime in primes:
    common_modulus *= prime
for path in paths:
    record = {
        "path": list(path),
        "anti_invariant_local_only": {
            str(prime): flat[prime][path]["anti_invariant_coefficient"]
            for prime in primes
        },
        "generator_free": {},
    }
    for name in ("trace", "norm", "coefficient_discriminant"):
        residues = {prime: flat[prime][path][name] for prime in primes}
        value, centered, modulus = crt(residues)
        if modulus != common_modulus:
            raise ArithmeticError("inconsistent CRT modulus")
        record["generator_free"][name] = {
            "residues": {str(prime): residues[prime] for prime in primes},
            "crt_nonnegative": value,
            "crt_centered": centered,
        }
    slot_records.append(record)

output = {
    "schema": "elkies-k3.q80-third-q12-multibranch-residue-ledger.v3",
    "status": "PASS_EXACT_MULTIBRANCH_FROBENIUS_RESIDUE_LEDGER",
    "specialization": {"u": "-2"},
    "residue_primes": list(primes),
    "held_out_good_primes": [71],
    "crt_modulus": common_modulus,
    "ordered_coefficient_slots": len(slot_records),
    "generator_free_fields": ["trace", "norm", "coefficient_discriminant"],
    "local_only_fields": ["anti_invariant_coefficient"],
    "branch_diagnostic": {
        "global_field_degree": 4,
        "basis": ["1", "a", "b", "a*b"],
        "relations": ["a^2=q1", "b^2=q2"],
        "prime_records": branch_diagnostic,
        "conclusion": (
            "the local quadratic orbit alternates between q1 and q2; the CRT values "
            "below mix distinct global conjugation quotients and are not rational-"
            "reconstruction candidates"
        ),
    },
    "slots": slot_records,
    "inputs": {
        **{
            str(prime): {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
            for prime, path in inputs.items()
        },
        "held_out_audit": {
            "path": str(args.held_out_audit.relative_to(ROOT)),
            "sha256": sha256(args.held_out_audit),
        },
        "exact_operands": {
            "path": str(args.exact_operands.relative_to(ROOT)),
            "sha256": sha256(args.exact_operands),
        },
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            f"{len(primes)}-prime formal CRT accumulation for all 1,947 ordered local invariant slots",
            "literal residue replay at every listed reconstruction prime",
            "the local quadratic operand alternates between the two independent global square classes q1 and q2",
            "p=71 reserved and certified as a good held-out prime",
        ],
        "not_proved": [
            "that the branch-mixed CRT integers are reductions of common rational coefficients",
            "rational reconstruction from this retired interface",
            "a characteristic-zero child equation or Mordell--Weil rank",
        ],
    },
    "reproduce": "python3 elkies-k3/scripts/compile_q80_third_q12_frobenius_crt_interface.py",
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"CRT interface artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    f"Q80THIRDQ12RESIDUES|primes={','.join(map(str, primes))}|heldout=71|slots={len(slot_records)}|"
    f"modulus={common_modulus}|status=PASS_EXACT_MULTIBRANCH_FROBENIUS_RESIDUE_LEDGER"
)
