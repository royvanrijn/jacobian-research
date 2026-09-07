#!/usr/bin/env sage -python
"""CRT-lift the marked first-neighbor H92 section coefficient vectors."""

from sage.all import ZZ, crt

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "artifacts/generated-results/h92-p2-modular"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-p2-intermediate-crt.json"


records = [json.loads(path.read_text()) for path in sorted(INPUT.glob("intermediate-*.json"))]
if len(records) < 2:
    raise SystemExit("need at least two modular coefficient records")
primes = [ZZ(record["prime"]) for record in records]
modulus = ZZ.prod(primes)


def lift_vector(vectors):
    lifted = []
    for entries in zip(*vectors, strict=True):
        residue = ZZ(entries[0])
        current_modulus = primes[0]
        for entry, prime in zip(entries[1:], primes[1:], strict=True):
            residue = crt(residue, ZZ(entry), current_modulus, prime)
            current_modulus *= prime
        try:
            lifted.append(str(residue.rational_reconstruction(modulus)))
        except (ArithmeticError, ValueError):
            lifted.append(None)
    return lifted


payload = {
    "schema": "elkies-k3.h92-p2-intermediate-crt.v1",
    "primes": [int(prime) for prime in primes],
    "modulus_bits": int(modulus.nbits()),
    "x": {
        "numerator": lift_vector([record["x"][0] for record in records]),
        "denominator": lift_vector([record["x"][1] for record in records]),
    },
    "y": {
        "numerator": lift_vector([record["y"][0] for record in records]),
        "denominator": lift_vector([record["y"][1] for record in records]),
    },
}
payload["all_coefficients_reconstructed"] = not any(
    value is None
    for coordinate in ("x", "y")
    for part in ("numerator", "denominator")
    for value in payload[coordinate][part]
)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92P2CRT|records={}|modulus_bits={}|complete={}".format(
        len(records), payload["modulus_bits"], int(payload["all_coefficients_reconstructed"])
    ),
    flush=True,
)
