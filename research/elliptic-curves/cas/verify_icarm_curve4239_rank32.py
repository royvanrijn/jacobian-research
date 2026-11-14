#!/usr/bin/env python3
"""Exact finite-quotient replay of Elkies--Klagsbrun--Li's public rank>=32 curve.

This is an external reproduction, never a new discovery or an exact-rank claim.
No conductor factorization or BSD/GRH upper bound is used.
"""
import argparse
from dataclasses import asdict
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from mod2_reduction_independence import (
    combined_mod2_rank, find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime, mod2_reduction_signature,
    short_curve_has_no_rational_2_torsion_modular_certificate,
)

SOURCE = ROOT/"elliptic-curves/data/external/icarm_curve4239_20261113.json"
DEFAULT = ROOT/"artifacts/generated-results/elliptic-curves/icarm_curve4239_rank32_v1.json"


def prepare():
    data = json.loads(SOURCE.read_text())
    a1, a2, a3, a4, a6 = map(Q, data["ainvs"])
    points = [tuple(map(Q, row)) for row in data["points"]]
    assert len(points) == 32
    assert all(y*y+a1*x*y+a3*y == x**3+a2*x*x+a4*x+a6 for x, y in points)
    b2, b4, b6 = a1*a1+4*a2, a1*a3+2*a4, a3*a3+4*a6
    c4 = b2*b2-24*b4
    c6 = -b2**3+36*b2*b4-216*b6
    A, B = -27*c4, -54*c6
    short = [(36*x+3*b2, 108*(2*y+a1*x+a3)) for x, y in points]
    assert all(y*y == x**3+A*x+B for x, y in short)
    return [Q(0), Q(0), Q(0), A, B], short


def generate(prime_bound):
    started = time.process_time()
    curve, points = prepare()
    signatures = find_mod2_reduction_certificate(curve, points, prime_bound=prime_bound)
    rank = combined_mod2_rank(signatures, len(points))
    assert rank == 32, f"bounded certificate incomplete: binary rank {rank}; no dependence claim"
    torsion_prime = find_two_torsion_certificate_prime(curve)
    return {"schema": "icarm4239-rank32-replay-v1", "status": "PASS",
        "source_url": "https://elliptic-rank.icarm.cloud/curve/4239.json",
        "source_sha256": sha256(SOURCE.read_bytes()).hexdigest(),
        "short_ainvs": [str(x) for x in curve],
        "point_count": 32, "binary_rank": rank,
        "no_rational_2torsion_prime": torsion_prime,
        "signatures": [asdict(s) for s in signatures],
        "prime_bound": prime_bound, "cpu_seconds": round(time.process_time()-started, 6),
        "conclusion": "Unconditional rank >=32 by finite-quotient infinite descent.",
        "boundary": "Public external curve and points. No new curve, upper bound, saturation, conductor or original discovery is claimed."}


def verify(result):
    curve, points = prepare()
    assert result["source_sha256"] == sha256(SOURCE.read_bytes()).hexdigest()
    assert result["short_ainvs"] == [str(x) for x in curve]
    signatures = []
    for row in result["signatures"]:
        actual = mod2_reduction_signature(curve, points, row["prime"])
        assert json.loads(json.dumps(asdict(actual))) == json.loads(json.dumps(row))
        signatures.append(actual)
    assert combined_mod2_rank(signatures, len(points)) == result["binary_rank"] == 32
    assert short_curve_has_no_rational_2_torsion_modular_certificate(
        curve, result["no_rational_2torsion_prime"])
    return True


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=DEFAULT)
    p.add_argument("--check", action="store_true")
    p.add_argument("--prime-bound", type=int, default=1000)
    args = p.parse_args()
    if args.check:
        result = json.loads(args.output.read_text())
        verify(result)
    else:
        result = generate(args.prime_bound)
        verify(result)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
    print(json.dumps({"status": "PASS", "rank_lower_bound": 32,
        "certificate_primes": [s["prime"] for s in result["signatures"]],
        "matrix_rows": sum(len(s["rows"]) for s in result["signatures"]),
        "no_2torsion_prime": result["no_rational_2torsion_prime"]}))
