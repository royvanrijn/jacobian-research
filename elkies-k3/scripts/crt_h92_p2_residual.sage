#!/usr/bin/env sage -python
"""CRT-lift the degree-10 residual point that trivializes the H92 torsor."""

from sage.all import ZZ, crt

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUTS = sorted((ROOT / "artifacts/generated-results").glob("h92-p2-residual-mod-*.json"))
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-p2-residual-crt.json"
all_records = [json.loads(path.read_text()) for path in INPUTS]
expected_lengths = {
    coordinate: tuple(len(all_records[0][coordinate][part]) for part in (0, 1))
    for coordinate in ("residual_z", "residual_x")
}
records = [
    record for record in all_records
    if all(
        tuple(len(record[coordinate][part]) for part in (0, 1)) == expected_lengths[coordinate]
        for coordinate in expected_lengths
    )
]
if len(records) < 2:
    raise SystemExit("need two modular residual records")
primes = [ZZ(record["prime"]) for record in records]
modulus = ZZ.prod(primes)


def lift(vectors):
    answer = []
    length = len(vectors[0])
    if any(len(vector) != length for vector in vectors):
        raise ValueError("inconsistent modular coefficient vector lengths")
    for index in range(length):
        residue = ZZ(vectors[0][index])
        current = primes[0]
        for record_index in range(1, len(vectors)):
            residue = crt(residue, ZZ(vectors[record_index][index]), current, primes[record_index])
            current *= primes[record_index]
        try:
            answer.append(str(residue.rational_reconstruction(modulus)))
        except (ArithmeticError, ValueError):
            answer.append(None)
    return answer


payload = {
    "schema": "elkies-k3.h92-p2-torsor-residual-crt.v1",
    "primes": [int(prime) for prime in primes],
    "discarded_degenerate_records": len(all_records) - len(records),
    "modulus_bits": int(modulus.nbits()),
    "residual_z": {
        "numerator": lift([record["residual_z"][0] for record in records]),
        "denominator": lift([record["residual_z"][1] for record in records]),
    },
    "residual_x": {
        "numerator": lift([record["residual_x"][0] for record in records]),
        "denominator": lift([record["residual_x"][1] for record in records]),
    },
}
payload["complete"] = not any(
    coefficient is None
    for coordinate in ("residual_z", "residual_x")
    for part in ("numerator", "denominator")
    for coefficient in payload[coordinate][part]
)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"H92P2RESIDUALCRT|records={len(records)}|bits={modulus.nbits()}|"
    f"complete={int(payload['complete'])}",
    flush=True,
)
