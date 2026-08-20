#!/usr/bin/env python3

from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve273 import (
    A, B, POINTS, SHORT_POINTS,
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    short_coefficients,
    on_curve,
)

from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)

from search_extra_points import (
    gp_rational,
    gp_vector,
    run_gp,
)


PROTOCOL = "R30ICARM"

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts/local/elliptic-curves/curve273-rank30"
OUT.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Exact point membership
# ------------------------------------------------------------

print(
    f"{PROTOCOL}|stage=input"
    f"|points={len(POINTS)}",
    flush=True,
)

assert len(POINTS) == 30

for i, P in enumerate(POINTS, 1):
    if not on_curve(P):
        raise AssertionError(f"point {i} is off curve")

print(
    f"{PROTOCOL}|stage=membership"
    f"|checked=30|status=PASS",
    flush=True,
)


# ------------------------------------------------------------
# Exact finite-reduction independence certificate
# ------------------------------------------------------------

short = short_coefficients()

torsion_prime = find_two_torsion_certificate_prime(
    short,
    prime_bound=500,
)

print(
    f"{PROTOCOL}|stage=two_torsion"
    f"|certificate_prime={torsion_prime}"
    f"|status=NO_RATIONAL_2_TORSION",
    flush=True,
)

signatures = find_mod2_reduction_certificate(
    short,
    SHORT_POINTS,
    prime_bound=2000,
)

rank = combined_mod2_rank(
    signatures,
    len(SHORT_POINTS),
)

print(
    f"{PROTOCOL}|stage=mod2"
    f"|rank={rank}"
    f"|target=30"
    f"|primes={','.join(str(s.prime) for s in signatures)}",
    flush=True,
)

for sig in signatures:
    print(
        f"{PROTOCOL}|stage=prime"
        f"|p={sig.prime}"
        f"|group_order={sig.group_order}"
        f"|quotient_dim={sig.quotient_dimension}",
        flush=True,
    )

if rank != 30:
    raise RuntimeError(
        "bounded mod-2 certificate did not reach rank 30; "
        "this does NOT prove dependence -- raise prime bound"
    )

print(
    f"{PROTOCOL}|stage=certificate"
    f"|rank_lower_bound=30"
    f"|status=EXACT_UNCONDITIONAL",
    flush=True,
)


# ------------------------------------------------------------
# Neron-Tate height Gram
# ------------------------------------------------------------

curve = ",".join(
    gp_rational(x)
    for x in GENERAL_WEIERSTRASS_COEFFICIENTS
)

points = ",".join(
    gp_vector(P)
    for P in POINTS
)

program = "\n".join([
    "default(realprecision,120);",
    f"E=ellinit([{curve}]);",
    f"P=[{points}];",
    "H=ellheightmatrix(E,P);",
    'print("HEIGHT_BEGIN");',
    "for(i=1,matsize(H)[1],print(Vec(H[i,])));",
    'print("HEIGHT_END");',
    'print("REGULATOR|",matdet(H));',
    "quit",
]) + "\n"

gp_output, wall = run_gp(
    program,
    timeout=300.0,
    stack_bytes=2_000_000_000,
)

lines = [
    line.strip()
    for line in gp_output.splitlines()
    if line.strip()
]

start = lines.index("HEIGHT_BEGIN") + 1
end = lines.index("HEIGHT_END")

height_rows = lines[start:end]

if len(height_rows) != 30:
    raise AssertionError(
        f"expected 30 height rows, got {len(height_rows)}"
    )

regulator_line = next(
    line for line in lines
    if line.startswith("REGULATOR|")
)

regulator = regulator_line.split("|", 1)[1]

height_path = OUT / "height-gram.txt"
height_path.write_text(
    "\n".join(height_rows) + "\n"
)

points_path = OUT / "points.txt"
with points_path.open("w") as f:
    for x, y in POINTS:
        f.write(f"{x}\t{y}\n")

certificate = {
    "status": "exact_unconditional_rank_at_least_30",
    "curve": {
        "ainvs": [str(x) for x in GENERAL_WEIERSTRASS_COEFFICIENTS],
    },
    "point_count": 30,
    "membership_checks": 30,
    "two_torsion_certificate_prime": torsion_prime,
    "mod2_rank": rank,
    "certificate_primes": [
        sig.prime for sig in signatures
    ],
    "signatures": [
        {
            "prime": sig.prime,
            "group_order": sig.group_order,
            "doubled_subgroup_order": sig.doubled_subgroup_order,
            "quotient_dimension": sig.quotient_dimension,
            "rows": [list(r) for r in sig.rows],
        }
        for sig in signatures
    ],
    "regulator_numeric": regulator,
    "height_matrix_path": str(height_path),
    "pari_height_wall_seconds": wall,
}

cert_path = OUT / "rank30-certificate.json"
cert_path.write_text(
    json.dumps(certificate, indent=2) + "\n"
)

print(
    f"{PROTOCOL}|stage=height"
    f"|rows=30"
    f"|regulator={regulator}"
    f"|seconds={wall:.3f}",
    flush=True,
)

print(
    f"{PROTOCOL}|stage=done"
    f"|status=PASS"
    f"|certificate={cert_path}",
    flush=True,
)
