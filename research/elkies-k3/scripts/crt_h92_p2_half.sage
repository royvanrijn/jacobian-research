#!/usr/bin/env sage -python
"""CRT lift of the rational H92 P2 half extracted modulo good primes."""

from pathlib import Path
import json
import os

from sage.all import ZZ, crt


ROOT = Path(__file__).resolve().parents[2]
paths = sorted((ROOT / "artifacts/generated-results").glob(
    os.environ.get("H92P2_HALF_GLOB", "h92-p2-half-mod-*.json")
))
records = [json.loads(path.read_text()) for path in paths]
records = [record for record in records if all(coordinate in record for coordinate in ("x", "y"))]
records = [record for record in records if all(
    len(record[coordinate][part]) == expected
    for coordinate, parts in (("x", (47, 43)), ("y", (70, 64)))
    for part, expected in zip(("numerator", "denominator"), parts)
)]
if len(records) < 2:
    raise SystemExit("need at least two nondegenerate half-section records")
primes = [ZZ(record["prime"]) for record in records]
modulus = ZZ.prod(primes)


def lift(vectors):
    result = []
    for coefficients in zip(*vectors):
        residue, current = ZZ(coefficients[0]), primes[0]
        for value, prime in zip(coefficients[1:], primes[1:]):
            residue = crt(residue, ZZ(value), current, prime)
            current *= prime
        try:
            result.append(str(residue.rational_reconstruction(modulus)))
        except (ArithmeticError, ValueError):
            result.append(None)
    return result


payload = {
    "schema": "elkies-k3.h92-p2-half-crt.v1",
    "primes": [int(prime) for prime in primes],
    "modulus_bits": int(modulus.nbits()),
    "x": {part: lift([record["x"][part] for record in records]) for part in ("numerator", "denominator")},
    "y": {part: lift([record["y"][part] for record in records]) for part in ("numerator", "denominator")},
}
payload["complete"] = not any(
    coefficient is None
    for coordinate in ("x", "y")
    for part in ("numerator", "denominator")
    for coefficient in payload[coordinate][part]
)
output = ROOT / "artifacts/generated-results" / os.environ.get(
    "H92P2_HALF_CRT_OUTPUT", "elkies-k3-h92-p2-half-crt.json"
)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"H92P2HALFCRT|records={len(records)}|bits={modulus.nbits()}|complete={int(payload['complete'])}")
