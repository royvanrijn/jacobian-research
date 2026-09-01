#!/usr/bin/env python3
"""Compile generator-free q12 coefficient slots into a CRT interface."""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--prime-input", action="append", default=[], metavar="PRIME=PATH",
    help=(
        "certified Frobenius-invariant input; repeat to replace the default "
        "reconstruction-prime set"
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
for name in ("held_out_audit", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


default_inputs = {
    prime: RESULTS / f"q80-third-q12-p{prime}-frobenius-invariants.json"
    for prime in (19, 61, 67, 83)
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
    raise ValueError("the reconstruction set must include p=19 and at least two controls")

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
parameters = held_out.get("parameters", [])
if len(parameters) != 1 or parameters[0].get("u") != "-2":
    raise ValueError("held-out audit is not the u=-2 specialization")
modular = parameters[0].get("modular", [])
if len(modular) != 1 or modular[0].get("prime") != 71 or modular[0].get("status") != "PASS_GOOD_REDUCTION_AUDIT":
    raise ValueError("p=71 is not certified as a held-out good prime")


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
    "schema": "elkies-k3.q80-third-q12-frobenius-crt-interface.v2",
    "status": "PASS_EXACT_THIRD_Q12_FROBENIUS_CRT_INTERFACE",
    "specialization": {"u": "-2"},
    "reconstruction_primes": list(primes),
    "held_out_good_primes": [71],
    "crt_modulus": common_modulus,
    "ordered_coefficient_slots": len(slot_records),
    "generator_free_fields": ["trace", "norm", "coefficient_discriminant"],
    "local_only_fields": ["anti_invariant_coefficient"],
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
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            f"{len(primes)}-prime CRT accumulation for all 1,947 ordered generator-free coefficient slots",
            "literal residue replay at every listed reconstruction prime",
            "p=71 reserved and certified as a good held-out prime",
        ],
        "not_proved": [
            "rational reconstruction from the current small CRT modulus",
            "a common characteristic-zero quadratic field",
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
    f"Q80THIRDQ12CRT|primes={','.join(map(str, primes))}|heldout=71|slots={len(slot_records)}|"
    f"modulus={common_modulus}|status=PASS_EXACT_THIRD_Q12_FROBENIUS_CRT_INTERFACE"
)
