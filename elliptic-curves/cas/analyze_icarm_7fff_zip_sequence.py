#!/usr/bin/env python3
"""Exact audit of ICARM curves 281, 282, 285 and 286.

The checker verifies the pinned equations and points, recomputes the standard
integral invariants, proves trivial rational torsion by good reduction, and
emits conservative cross-curve fingerprints.  For the rank-at-least-21 curves
285 and 286 it additionally replays global minimalization and every local Tate
reduction with PARI/GP.  It does not identify an algebraic family.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, isqrt, log, prod
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_7fff_zip_public_source_281_282_285_286.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_7fff_zip_independence_analysis_v1.json"
EXPECTED_IDS = (281, 282, 285, 286)
GLOBAL_REDUCTION_IDS = (285, 286)
GP_TIMEOUT_SECONDS = 30.0


def invariants(ainvs: list[str]) -> dict[str, int]:
    a1, a2, a3, a4, a6 = map(int, ainvs)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    c4 = b2 * b2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    delta = -b2 * b2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    assert c4**3 - c6**2 == 1728 * delta
    return {"b2": b2, "c4": c4, "c6": c6, "discriminant": delta}


def on_curve(ainvs: list[str], point: list[str]) -> bool:
    a1, a2, a3, a4, a6 = map(Fraction, ainvs)
    x, y = map(Fraction, point)
    return y * y + a1 * x * y + a3 * y == x**3 + a2 * x * x + a4 * x + a6


def count_points_mod_p(ainvs: list[str], p: int) -> int:
    a1, a2, a3, a4, a6 = (int(a) % p for a in ainvs)
    total = 1
    for x in range(p):
        for y in range(p):
            if (y*y + a1*x*y + a3*y - x**3 - a2*x*x - a4*x - a6) % p == 0:
                total += 1
    return total


def torsion_certificate(curve: dict, delta: int) -> dict:
    orders = []
    running_gcd = 0
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
        if delta % p == 0:
            continue
        order = count_points_mod_p(curve["ainvs"], p)
        running_gcd = gcd(running_gcd, order)
        orders.append({"p": p, "order": order, "running_gcd": running_gcd})
        if running_gcd == 1:
            break
    assert running_gcd == 1
    return {"good_reductions": orders, "gcd": running_gcd, "conclusion": "E(Q)_tors is trivial"}


def discriminant_signature(curve: dict, delta: int) -> dict:
    remainder = abs(delta)
    valuations = []
    for text in curve["bad_primes"]:
        p = int(text)
        exponent = 0
        while remainder % p == 0:
            remainder //= p
            exponent += 1
        assert exponent > 0
        valuations.append([str(p), exponent])
    assert remainder == 1
    return {
        "valuations_on_reported_bad_primes": valuations,
        "repeated_bad_prime_count": sum(e > 1 for _, e in valuations),
        "exponent_one_bad_prime_count": sum(e == 1 for _, e in valuations),
        "support_size": len(valuations),
    }


def denominator_fingerprint(points: list[list[str]]) -> dict:
    roots = []
    for x_text, _ in points:
        denominator = Fraction(x_text).denominator
        root = isqrt(denominator)
        assert root * root == denominator
        roots.append(root)
    roots.sort()
    return {
        "integral_x_count": roots.count(1),
        "x_denominator_square_roots": [str(x) for x in roots],
        "median_root": str(roots[len(roots) // 2]),
        "maximum_root": str(roots[-1]),
    }


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def kodaira_symbol(code: int) -> str:
    """Decode PARI's integer Kodaira convention."""
    if code == 1:
        return "I0"
    if code == 2:
        return "II"
    if code == 3:
        return "III"
    if code == 4:
        return "IV"
    if code >= 5:
        return f"I{code - 4}"
    if code == -1:
        return "I0*"
    if code == -2:
        return "IV*"
    if code == -3:
        return "III*"
    if code == -4:
        return "II*"
    if code <= -5:
        return f"I{-code - 4}*"
    raise ValueError(f"unknown PARI Kodaira code {code}")


def pari_global_reduction(curve: dict, discriminant: int) -> dict:
    """Replay global minimality and Tate data for one pinned integral model."""
    gp = shutil.which("gp")
    if gp is None:
        raise FileNotFoundError("PARI/GP executable 'gp' is required")
    model = ",".join(curve["ainvs"])
    program = f"""
default(realprecision,100);
E=ellinit([{model}]);v=0;M=ellminimalmodel(E,&v);G=ellglobalred(M);
print("PARI_VERSION\\t",version());
print("MIN_MODEL\\t",M.a1,"\\t",M.a2,"\\t",M.a3,"\\t",M.a4,"\\t",M.a6);
print("MIN_CHANGE\\t",v[1],"\\t",v[2],"\\t",v[3],"\\t",v[4]);
print("MIN_DISC\\t",M.disc);
print("CONDUCTOR\\t",G[1]);
print("TAMAGAWA_PRODUCT\\t",G[3]);
print("ROOT_NUMBER\\t",ellrootno(M));
for(i=1,matsize(G[4])[1],p=G[4][i,1];L=elllocalred(M,p);print("LOCAL\\t",p,"\\t",L[1],"\\t",L[2],"\\t",L[3][1],"\\t",L[3][2],"\\t",L[3][3],"\\t",L[3][4],"\\t",L[4]));
"""
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=GP_TIMEOUT_SECONDS,
    )
    combined = completed.stdout + completed.stderr
    if "***" in combined:
        raise RuntimeError(combined)

    scalar: dict[str, list[str]] = {}
    local_rows: list[list[str]] = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if fields[0] == "LOCAL":
            local_rows.append(fields[1:])
        elif len(fields) > 1:
            scalar[fields[0]] = fields[1:]

    minimal_model = scalar["MIN_MODEL"]
    minimal_change = scalar["MIN_CHANGE"]
    minimal_discriminant = int(scalar["MIN_DISC"][0])
    conductor = int(scalar["CONDUCTOR"][0])
    assert minimal_model == curve["ainvs"]
    assert minimal_change == ["1", "0", "0", "0"]
    assert minimal_discriminant == discriminant
    assert conductor == int(curve["conductor"])

    local_reductions = []
    for row in local_rows:
        prime, conductor_exponent, kodaira_code, u, r, s, t, tamagawa = map(int, row)
        local_reductions.append({
            "prime": str(prime),
            "minimal_discriminant_valuation": valuation(abs(minimal_discriminant), prime),
            "conductor_exponent": conductor_exponent,
            "kodaira_code": kodaira_code,
            "kodaira_symbol": kodaira_symbol(kodaira_code),
            "local_minimal_change": [u, r, s, t],
            "tamagawa_number": tamagawa,
        })
    assert [entry["prime"] for entry in local_reductions] == curve["bad_primes"]
    reconstructed_conductor = prod(
        int(entry["prime"]) ** entry["conductor_exponent"]
        for entry in local_reductions
    )
    assert reconstructed_conductor == conductor
    assert prod(entry["tamagawa_number"] for entry in local_reductions) == int(
        scalar["TAMAGAWA_PRODUCT"][0]
    )

    version = scalar["PARI_VERSION"][0].strip("[]").replace(", ", ".")
    return {
        "engine": f"PARI/GP {version}",
        "global_minimal_model": minimal_model,
        "global_minimal_change": [int(value) for value in minimal_change],
        "source_model_is_global_minimal": True,
        "minimal_discriminant": str(minimal_discriminant),
        "conductor": str(conductor),
        "reported_conductor_matches": True,
        "conductor_reconstructed_from_local_exponents": str(reconstructed_conductor),
        "tamagawa_product": int(scalar["TAMAGAWA_PRODUCT"][0]),
        "root_number": int(scalar["ROOT_NUMBER"][0]),
        "local_reductions": local_reductions,
    }


def independence_certificate(curve: dict, inv: dict[str, int]) -> dict:
    """Certify independence after the standard integral short transformation."""
    a1, _a2, a3, _a4, _a6 = map(int, curve["ainvs"])
    b2 = inv["b2"]
    short = (Fraction(0), Fraction(0), Fraction(0), Fraction(-27 * inv["c4"]), Fraction(-54 * inv["c6"]))
    points = []
    for x_text, y_text in curve["points"]:
        x, y = Fraction(x_text), Fraction(y_text)
        points.append((36 * x + 3 * b2, 108 * (2 * y + a1 * x + a3)))
    signatures = find_mod2_reduction_certificate(short, points, prime_bound=1000)
    rank = combined_mod2_rank(signatures, len(points))
    assert rank == len(points)
    return {
        "method": "exact images in products of E(F_p)/2E(F_p)",
        "short_model": [str(value) for value in short],
        "combined_binary_rank": rank,
        "primes": [row.prime for row in signatures],
        "rows": [
            {
                "prime": row.prime,
                "group_order": row.group_order,
                "quotient_dimension": row.quotient_dimension,
                "matrix_rows": [list(bits) for bits in row.rows],
            }
            for row in signatures
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    raw = args.input.read_bytes()
    source = json.loads(raw)
    curves = sorted(source["curves"], key=lambda c: c["id"])
    assert tuple(c["id"] for c in curves) == EXPECTED_IDS
    assert {c["submitter"] for c in curves} == {"7fff-zip"}

    audited = []
    for curve in curves:
        inv = invariants(curve["ainvs"])
        assert inv["discriminant"] == int(curve["discriminant"])
        assert [str(inv["c4"]), str(inv["c6"])] == curve["curve_key"].split(":")
        assert inv["discriminant"] != 0
        assert len(curve["points"]) == curve["rank_lower_bound"]
        assert all(on_curve(curve["ainvs"], point) for point in curve["points"])
        record = {
            "id": curve["id"],
            "rank_lower_bound_source_claim": curve["rank_lower_bound"],
            "created_at_utc": curve["created_at"],
            "point_membership_verified": len(curve["points"]),
            "independence_replayed": True,
            "invariants": {key: str(value) for key, value in inv.items()},
            "log_conductor": log(int(curve["conductor"])),
            "discriminant_signature": discriminant_signature(curve, inv["discriminant"]),
            "denominator_fingerprint": denominator_fingerprint(curve["points"]),
            "torsion_certificate": torsion_certificate(curve, inv["discriminant"]),
            "independence_certificate": independence_certificate(curve, inv),
        }
        if curve["id"] in GLOBAL_REDUCTION_IDS:
            record["pari_global_and_local_reduction"] = pari_global_reduction(
                curve, inv["discriminant"]
            )
        audited.append(record)

    pairwise = []
    for i, left in enumerate(audited):
        for right in audited[i + 1:]:
            lroots = set(left["denominator_fingerprint"]["x_denominator_square_roots"])
            rroots = set(right["denominator_fingerprint"]["x_denominator_square_roots"])
            lc4 = int(left["invariants"]["c4"])
            lc6 = int(left["invariants"]["c6"])
            ld = int(left["invariants"]["discriminant"])
            rc4 = int(right["invariants"]["c4"])
            rc6 = int(right["invariants"]["c6"])
            rd = int(right["invariants"]["discriminant"])
            same_j = lc4**3 * rd == rc4**3 * ld
            pairwise.append({
                "ids": [left["id"], right["id"]],
                "same_j_invariant": same_j,
                "quadratic_twists_or_Qbar_isomorphic": same_j,
                "shared_denominator_roots": sorted(lroots & rroots, key=int),
            })
    assert not any(pair["same_j_invariant"] for pair in pairwise)

    output = {
        "schema": "icarm-7fff-zip-sequence-analysis-v1",
        "input_sha256": sha256(raw).hexdigest(),
        "scope": (
            "exact equation/point/invariant/torsion audit, PARI global and local "
            "reduction for curves 285/286, plus descriptive fingerprints"
        ),
        "curves": audited,
        "pairwise": pairwise,
        "interpretation": {
            "proved": [
                "all 81 submitted rational points lie on their stated curves",
                "all four curves are nonsingular and have trivial rational torsion",
                "the submitted point sets are independent by exact finite-reduction certificates",
                "no pair has the same j-invariant, hence no pair is a quadratic-twist pair",
                "curves 285 and 286 are already global minimal models and their reported conductors are reconstructed from complete local Tate data",
            ],
            "observed": [
                "every discriminant has repeated powers on several small bad primes",
                "the point-denominator sets have no substantial four-curve common core",
                "curve 285 has only six bad primes, the sparsest support of the sequence",
            ],
            "not_proved": [
                "membership in a single parametrized family or use of a specific construction",
            ],
        },
    }
    rendered = json.dumps(output, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit(f"stale or missing artifact: {args.output}")
        print(f"PASS {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        print(f"wrote {args.output}")
    print("verified 4 curves, 81 points, 4 trivial-torsion certificates, 6 unequal j-pairs")


if __name__ == "__main__":
    main()
