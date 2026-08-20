#!/usr/bin/env python3
"""Certify the public ICARM curve 273 has Mordell--Weil rank at least 30.

The rank lower bound uses exact point membership and finite good-reduction
quotients only. The PARI height matrix is retained as a diagnostic under the
ignored local-artifact tree; it is not part of the independence proof.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve273 import (  # noqa: E402
    POINTS,
    SHORT_POINTS,
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    on_curve,
    short_coefficients,
)
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from pari_bridge import pari_version  # noqa: E402
from search_extra_points import gp_rational, gp_vector, run_gp  # noqa: E402


PROTOCOL = "R30ICARM"
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "icarm_curve273_rank30_v1.json"
)
DEFAULT_HEIGHT_DIR = ROOT / "artifacts/local/elliptic-curves/curve273-rank30"
REPRODUCING_COMMAND = (
    ".venv/bin/python elliptic-curves/cas/verify_icarm_curve273_rank30.py"
)


def parse_gp_vector(value: str) -> list[str]:
    value = value.strip()
    if not value.startswith("[") or not value.endswith("]"):
        raise ValueError(f"expected GP vector, got {value!r}")
    body = value[1:-1].strip()
    return [] if not body else [item.strip() for item in body.split(",")]


def run_pari_diagnostics(height_dir: Path) -> dict[str, object]:
    curve = ",".join(gp_rational(x) for x in GENERAL_WEIERSTRASS_COEFFICIENTS)
    points = ",".join(gp_vector(point) for point in POINTS)
    program = "\n".join(
        [
            "default(realprecision,120);",
            f"E=ellinit([{curve}]);",
            "M=ellminimalmodel(E);",
            'print("AINVS|",Vec(M)[1..5]);',
            'print("CONDUCTOR|",ellglobalred(E)[1]);',
            'print("ROOT_NUMBER|",ellrootno(E));',
            'print("TORSION_ORDER|",elltors(E)[1]);',
            'print("DISCRIMINANT|",E.disc);',
            f"P=[{points}];",
            "H=ellheightmatrix(E,P);",
            'print("HEIGHT_BEGIN");',
            "for(i=1,matsize(H)[1],print(Vec(H[i,])));",
            'print("HEIGHT_END");',
            'print("REGULATOR|",matdet(H));',
            "quit",
        ]
    ) + "\n"
    gp_output, _wall = run_gp(
        program,
        timeout=300.0,
        stack_bytes=2_000_000_000,
    )
    lines = [line.strip() for line in gp_output.splitlines() if line.strip()]
    start = lines.index("HEIGHT_BEGIN") + 1
    end = lines.index("HEIGHT_END")
    height_rows = lines[start:end]
    if len(height_rows) != 30:
        raise AssertionError(f"expected 30 height rows, got {len(height_rows)}")

    fields = {}
    for label in (
        "AINVS",
        "CONDUCTOR",
        "ROOT_NUMBER",
        "TORSION_ORDER",
        "DISCRIMINANT",
        "REGULATOR",
    ):
        line = next(item for item in lines if item.startswith(f"{label}|"))
        fields[label] = line.split("|", 1)[1]

    height_dir.mkdir(parents=True, exist_ok=True)
    (height_dir / "height-gram.txt").write_text(
        "\n".join(height_rows) + "\n", encoding="utf-8"
    )
    with (height_dir / "points.txt").open("w", encoding="utf-8") as handle:
        for x_value, y_value in POINTS:
            handle.write(f"{x_value}\t{y_value}\n")

    minimal_ainvs = parse_gp_vector(str(fields["AINVS"]))
    input_ainvs = [str(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS]
    if minimal_ainvs != input_ainvs:
        raise AssertionError("the supplied integral model is not global minimal")

    return {
        "global_minimal_model": minimal_ainvs,
        "minimal_model_same": True,
        "conductor": str(fields["CONDUCTOR"]),
        "root_number": int(str(fields["ROOT_NUMBER"])),
        "torsion_order": int(str(fields["TORSION_ORDER"])),
        "discriminant": str(fields["DISCRIMINANT"]),
        "regulator_numeric_120_digits": str(fields["REGULATOR"]),
    }


def build_certificate(height_dir: Path) -> dict[str, object]:
    print(f"{PROTOCOL}|stage=input|points={len(POINTS)}", flush=True)
    if len(POINTS) != 30:
        raise AssertionError("expected exactly 30 published points")
    for index, point in enumerate(POINTS, 1):
        if not on_curve(point):
            raise AssertionError(f"point {index} is off curve")
    print(f"{PROTOCOL}|stage=membership|checked=30|status=PASS", flush=True)

    short = short_coefficients()
    torsion_prime = find_two_torsion_certificate_prime(short, prime_bound=500)
    signatures = find_mod2_reduction_certificate(
        short,
        SHORT_POINTS,
        prime_bound=2000,
    )
    rank = combined_mod2_rank(signatures, len(SHORT_POINTS))
    if rank != 30:
        raise RuntimeError(
            "bounded mod-2 certificate did not reach rank 30; this does not "
            "prove dependence"
        )
    print(
        f"{PROTOCOL}|stage=mod2|rank={rank}|target=30"
        f"|primes={','.join(str(signature.prime) for signature in signatures)}",
        flush=True,
    )

    diagnostics = run_pari_diagnostics(height_dir)
    if diagnostics["torsion_order"] != 1:
        raise AssertionError("PARI did not report trivial rational torsion")
    if diagnostics["root_number"] != 1:
        raise AssertionError("the pinned root number changed")
    if int(str(diagnostics["discriminant"])) == 0:
        raise AssertionError("the curve is singular")

    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    model_path = Path(__file__).with_name("icarm_curve273.py")
    model_hash = hashlib.sha256(model_path.read_bytes()).hexdigest()
    return {
        "schema_version": 1,
        "artifact_kind": "exact_elliptic_curve_rank_lower_bound",
        "curve_id": "icarm_curve_273",
        "claim": "rank E(Q) >= 30",
        "claim_status": "exact unconditional lower bound; no exact-rank claim",
        "public_source": "https://web.math.pmf.unizg.hr/~duje/tors/rk30.html",
        "curve": {
            "ainvs": [str(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS],
            "global_minimal_model": diagnostics["global_minimal_model"],
            "minimal_model_same": diagnostics["minimal_model_same"],
            "discriminant": diagnostics["discriminant"],
            "conductor": diagnostics["conductor"],
            "root_number": diagnostics["root_number"],
            "torsion_order": diagnostics["torsion_order"],
        },
        "points": [[str(x_value), str(y_value)] for x_value, y_value in POINTS],
        "point_membership_checks": 30,
        "independence_certificate": {
            "method": "finite good-reduction quotients E(F_p)/2E(F_p)",
            "relation_prime": 2,
            "no_rational_2_torsion_witness_prime": torsion_prime,
            "combined_binary_rank": rank,
            "rows": [
                {
                    "prime": signature.prime,
                    "group_order": signature.group_order,
                    "doubled_subgroup_order": signature.doubled_subgroup_order,
                    "quotient_dimension": signature.quotient_dimension,
                    "matrix_rows": [list(row) for row in signature.rows],
                }
                for signature in signatures
            ],
        },
        "height_diagnostic": {
            "regulator_numeric_120_digits": diagnostics[
                "regulator_numeric_120_digits"
            ],
            "used_for_rank_claim": False,
        },
        "generation": {
            "command": REPRODUCING_COMMAND,
            "arithmetic": "exact rational and exhaustive finite-field group operations",
            "pari_gp": pari_version(),
            "checker_sha256": script_hash,
            "model_data_sha256": model_hash,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--height-output-dir", type=Path, default=DEFAULT_HEIGHT_DIR)
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify that --output already equals the deterministic certificate",
    )
    args = parser.parse_args()
    certificate = build_certificate(args.height_output_dir)
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file():
            raise SystemExit(f"missing pinned certificate: {args.output}")
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(
                f"stale pinned certificate: rerun {REPRODUCING_COMMAND}"
            )
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(
        f"{PROTOCOL}|stage=done|status=PASS|rank_lower_bound=30"
        f"|mode={'check' if args.check else 'write'}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
