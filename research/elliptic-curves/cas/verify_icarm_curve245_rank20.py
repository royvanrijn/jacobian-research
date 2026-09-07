#!/usr/bin/env python3
"""Certify that public ICARM curve 245 has rank at least 20.

The rank lower bound uses exact point membership and exact finite-reduction
quotients.  Its conductor is reconstructed from a complete, explicitly pinned
factorization of the discriminant and PARI local-reduction data at precisely
those primes; no unbounded integer factorization is part of the replay.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from certify_nagao_rank17_frontier import (  # noqa: E402
    exact_log_conductor_certificate,
)
from icarm_curve245 import (  # noqa: E402
    BAD_PRIMES,
    CONDUCTOR,
    DISCRIMINANT,
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    POINTS,
    SHORT_POINTS,
    on_curve,
    short_coefficients,
)
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from pari_bridge import pari_version  # noqa: E402
from search_extra_points import gp_rational, run_gp  # noqa: E402


PROTOCOL = "R20ICARM245"
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "icarm_curve245_rank20_v1.json"
)
REPRODUCING_COMMAND = (
    ".venv/bin/python elliptic-curves/cas/verify_icarm_curve245_rank20.py"
)

# This is a complete factorization: its product is checked against the exact
# generalized-Weierstrass discriminant stored in icarm_curve245.py.
DISCRIMINANT_FACTORIZATION = (
    (2, 17),
    (3, 7),
    (5, 4),
    (13, 4),
    (19, 5),
    (37, 2),
    (7770053, 1),
    (763973980372286963203, 1),
    (55722582408764114465841769948159, 1),
)

LOCAL_PATTERN = re.compile(
    r"^LOCAL\|(\d+)\|(\d+)\|(-?\d+)\|\[([^]]+)\]\|(\d+)\|(-?\d+)$",
    re.MULTILINE,
)


def local_reduction_diagnostics() -> dict[str, object]:
    ainvs = ",".join(
        gp_rational(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS
    )
    commands = [
        "default(realprecision,80);",
        f"E=ellinit([{ainvs}]);",
    ]
    for prime in BAD_PRIMES:
        commands.extend(
            [
                f"L=elllocalred(E,{prime});",
                (
                    f'print("LOCAL|{prime}|",L[1],"|",L[2],"|",'
                    'L[3],"|",L[4],"|",ellrootno(E,' + str(prime) + "));"
                ),
            ]
        )
    commands.extend(
        [
            'print("MINIMAL|",Vec(ellminimalmodel(E))[1..5]);',
            'print("ROOT|",ellrootno(E));',
            f'print("LOGN|",log({CONDUCTOR}));',
            "quit",
        ]
    )
    output, _wall = run_gp(
        "\n".join(commands) + "\n",
        timeout=60.0,
        stack_bytes=500_000_000,
    )

    local_rows = []
    for match in LOCAL_PATTERN.finditer(output):
        prime = int(match.group(1))
        change = [Fraction(item.strip()) for item in match.group(4).split(",")]
        local_rows.append(
            {
                "prime": prime,
                "conductor_exponent": int(match.group(2)),
                "kodaira_code": int(match.group(3)),
                "minimal_change": [str(value) for value in change],
                "tamagawa_number": int(match.group(5)),
                "local_root_number": int(match.group(6)),
            }
        )
    if [row["prime"] for row in local_rows] != list(BAD_PRIMES):
        raise AssertionError("PARI local-reduction output was incomplete")
    if any(row["minimal_change"] != ["1", "0", "0", "0"] for row in local_rows):
        raise AssertionError("the displayed equation is not locally minimal")

    minimal_match = re.search(r"^MINIMAL\|\[(.*?)\]$", output, re.MULTILINE)
    root_match = re.search(r"^ROOT\|(-?\d+)$", output, re.MULTILINE)
    log_match = re.search(r"^LOGN\|(\S+)$", output, re.MULTILINE)
    if minimal_match is None or root_match is None or log_match is None:
        raise AssertionError("PARI omitted a global diagnostic")
    minimal = [item.strip() for item in minimal_match.group(1).split(",")]
    expected = [str(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS]
    if minimal != expected:
        raise AssertionError("PARI changed the displayed global minimal model")

    replayed_conductor = 1
    for row in local_rows:
        replayed_conductor *= int(row["prime"]) ** int(row["conductor_exponent"])
    if replayed_conductor != CONDUCTOR:
        raise AssertionError("local conductor exponents do not replay the conductor")

    return {
        "local_reductions": local_rows,
        "global_minimal_model": minimal,
        "root_number": int(root_match.group(1)),
        "log_conductor_numeric_80_digits": log_match.group(1),
        "conductor_from_local_reductions": str(replayed_conductor),
    }


def build_certificate() -> dict[str, object]:
    if len(POINTS) != 20:
        raise AssertionError("expected exactly 20 public points")
    for index, point in enumerate(POINTS, 1):
        if not on_curve(point):
            raise AssertionError(f"point {index} is off the curve")
    print(f"{PROTOCOL}|stage=membership|checked=20|status=PASS", flush=True)

    factored_discriminant = 1
    for prime, exponent in DISCRIMINANT_FACTORIZATION:
        factored_discriminant *= prime**exponent
    if factored_discriminant != abs(DISCRIMINANT):
        raise AssertionError("pinned discriminant factorization is incomplete")
    if tuple(prime for prime, _exponent in DISCRIMINANT_FACTORIZATION) != BAD_PRIMES:
        raise AssertionError("bad-prime list differs from discriminant support")

    short = short_coefficients()
    torsion_prime = find_two_torsion_certificate_prime(short, prime_bound=100)
    signatures = find_mod2_reduction_certificate(
        short,
        SHORT_POINTS,
        prime_bound=200,
    )
    rank = combined_mod2_rank(signatures, len(SHORT_POINTS))
    if rank != 20:
        raise RuntimeError("bounded mod-2 certificate did not reach rank 20")
    print(
        f"{PROTOCOL}|stage=mod2|rank=20"
        f"|primes={','.join(str(row.prime) for row in signatures)}",
        flush=True,
    )

    diagnostics = local_reduction_diagnostics()
    if diagnostics["root_number"] != 1:
        raise AssertionError("the pinned root number changed")
    exact_log_bound = exact_log_conductor_certificate(CONDUCTOR)

    script_path = Path(__file__)
    model_path = script_path.with_name("icarm_curve245.py")
    return {
        "schema_version": 1,
        "artifact_kind": "exact_elliptic_curve_rank_lower_bound",
        "curve_id": "icarm_curve_245",
        "claim": "rank E(Q) >= 20 and log conductor < 182.72",
        "claim_status": (
            "exact unconditional lower bound; no twenty-first point and no "
            "rank upper bound claimed"
        ),
        "public_source": "https://elliptic-rank.icarm.cloud/curve/245",
        "curve": {
            "ainvs": [str(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS],
            "global_minimal_model": diagnostics["global_minimal_model"],
            "discriminant": str(DISCRIMINANT),
            "discriminant_factorization": [
                [str(prime), exponent]
                for prime, exponent in DISCRIMINANT_FACTORIZATION
            ],
            "conductor": str(CONDUCTOR),
            "log_conductor_numeric_80_digits": diagnostics[
                "log_conductor_numeric_80_digits"
            ],
            "strict_log_conductor_target": "182.72",
            "below_strict_log_conductor_target": True,
            "exact_log_conductor_bound": exact_log_bound,
            "root_number": diagnostics["root_number"],
            "local_reductions": diagnostics["local_reductions"],
        },
        "points": [[str(x), str(y)] for x, y in POINTS],
        "point_membership_checks": 20,
        "independence_certificate": {
            "method": "finite good-reduction quotients E(F_p)/2E(F_p)",
            "relation_prime": 2,
            "no_rational_2_torsion_witness_prime": torsion_prime,
            "combined_binary_rank": rank,
            "rows": [
                {
                    "prime": row.prime,
                    "group_order": row.group_order,
                    "doubled_subgroup_order": row.doubled_subgroup_order,
                    "quotient_dimension": row.quotient_dimension,
                    "matrix_rows": [list(bits) for bits in row.rows],
                }
                for row in signatures
            ],
        },
        "generation": {
            "command": REPRODUCING_COMMAND,
            "arithmetic": "exact rational and exhaustive finite-field operations",
            "pari_gp": pari_version(),
            "checker_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
            "model_data_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    rendered = json.dumps(build_certificate(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file():
            raise SystemExit(f"missing pinned certificate: {args.output}")
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"stale pinned certificate: rerun {REPRODUCING_COMMAND}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(
        f"{PROTOCOL}|stage=done|status=PASS|rank_lower_bound=20"
        f"|mode={'check' if args.check else 'write'}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
