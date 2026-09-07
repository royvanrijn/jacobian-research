#!/usr/bin/env python3
"""Certify the known R17 Mordell--Weil quotient modulo 2 without a BNF.

For each selected specialization, map every supplied rational point to its
Kummer class ``x(P) - theta`` in the etale cubic algebra.  At odd auxiliary
primes where an integral two-division polynomial is squarefree and every
class is a unit, residue square characters give an exact homomorphism from
global squareclasses to an F_2-vector space.  A full-rank image therefore
certifies independence of the supplied Kummer classes.  This is a
lower-bound/labeling audit only: it neither computes the full 2-Selmer group
nor proves that the displayed local fingerprint is injective on unknown
classes.  It deliberately avoids maximal-order and class-group computation.

Run this file with Sage's Python, not ordinary CPython.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]

from run_elkies_2026_relative_2selmer_open import (  # noqa: E402
    GENERIC_RANK,
    f2_rank,
    load_authoritative_cases,
    selected_manifest_cases,
)
from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, prime_range  # noqa: E402
from sage.version import version as sage_version  # noqa: E402


INPUT_SCHEMA = "elliptic-curves.elkies-2026-relative-2selmer-suite-input.v1"
OUTPUT_SCHEMA = "elliptic-curves.elkies-2026-known-kummer-quotients.v1"
PROTOCOL = "ELKIESR17KNOWNKUMMER"


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def audit_case(case: Any, prime_bound: int) -> dict[str, Any]:
    started = time.monotonic()
    a1, a2, a3, a4, a6 = [QQ(value) for value in case.model]
    curve = EllipticCurve(QQ, [a1, a2, a3, a4, a6])
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    if not all(value.denominator() == 1 for value in (b2, 8 * b4, 16 * b6)):
        raise ArithmeticError("integral input did not produce an integral cubic")
    integer_ring = PolynomialRing(ZZ, "z")
    z = integer_ring.gen()
    # If theta is a root of the completed-square cubic, zeta=4*theta
    # satisfies this monic integral polynomial.  The element below is
    # 4*(x(P)-theta), and 4 is a global square.
    field_polynomial = (
        z**3 + ZZ(b2) * z**2 + ZZ(8 * b4) * z + ZZ(16 * b6)
    )
    polynomial_discriminant = abs(ZZ(field_polynomial.discriminant()))

    labelled_points = [
        (f"P{index + 1}", "specialized_generic_R17", point)
        for index, point in enumerate(case.generic_points)
    ]
    labelled_points.extend(
        (f"Q{index + 1}", "known_exceptional", point)
        for index, point in enumerate(case.exceptional_points)
    )
    x_coordinates = []
    for _label, _role, point in labelled_points:
        x_coordinate, y_coordinate = map(QQ, point)
        curve(x_coordinate, y_coordinate)
        x_coordinates.append(x_coordinate)
    rows: list[list[int]] = [[] for _ in x_coordinates]
    selected_primes = []
    for auxiliary_prime in prime_range(3, prime_bound + 1):
        q = int(auxiliary_prime)
        if polynomial_discriminant % q == 0 or any(
            x.denominator() % q == 0 for x in x_coordinates
        ):
            continue
        residue_field = GF(q)
        residue_ring = PolynomialRing(residue_field, "u")
        residue_polynomial = residue_ring(field_polynomial)
        if not residue_polynomial.is_squarefree():
            continue
        factors = [factor for factor, _exponent in residue_polynomial.factor()]
        extra_rows = []
        all_units = True
        for x_coordinate in x_coordinates:
            residue_x = residue_field(4 * x_coordinate)
            point_row = []
            for factor in factors:
                if factor.degree() == 1:
                    root = -factor[0] / factor[1]
                    residue_alpha = residue_x - root
                else:
                    extension = GF(
                        q ** factor.degree(), name="root", modulus=factor
                    )
                    residue_alpha = extension(residue_x) - extension.gen()
                if residue_alpha == 0:
                    all_units = False
                    break
                point_row.append(0 if residue_alpha.is_square() else 1)
            if not all_units:
                break
            extra_rows.append(point_row)
        if not all_units:
            continue
        candidate_rows = [
            old_row + extra_row for old_row, extra_row in zip(rows, extra_rows)
        ]
        rank_before = f2_rank(rows)
        rank_after = f2_rank(candidate_rows)
        if rank_after == rank_before:
            continue
        rows = candidate_rows
        selected_primes.append(
            {
                "prime": q,
                "prime_ideal_count": len(factors),
                "residue_factors": [str(factor) for factor in factors],
                "columns_added": len(extra_rows[0]),
                "rank_before": rank_before,
                "rank_after": rank_after,
                "generic_rank_after": f2_rank(rows[:GENERIC_RANK]),
            }
        )
        if rank_after == len(x_coordinates):
            break

    generic_rows = rows[:GENERIC_RANK]
    exceptional_rows = rows[GENERIC_RANK:]
    generic_rank = f2_rank(generic_rows)
    total_rank = f2_rank(rows)
    exceptional_quotient_dimension = total_rank - generic_rank
    expected_exceptional_dimension = len(exceptional_rows)
    full_known_rank_certified = total_rank == len(rows)
    quotient_basis_certified = (
        generic_rank == GENERIC_RANK
        and exceptional_quotient_dimension == expected_exceptional_dimension
    )

    point_records = []
    for index, ((label, role, _point), row) in enumerate(zip(labelled_points, rows)):
        quotient_coordinates = None
        if role == "known_exceptional" and quotient_basis_certified:
            exceptional_index = index - GENERIC_RANK
            quotient_coordinates = [
                int(position == exceptional_index)
                for position in range(expected_exceptional_dimension)
            ]
        point_records.append(
            {
                "label": label,
                "role": role,
                "local_squareclass_row": row,
                "exceptional_quotient_coordinates": quotient_coordinates,
            }
        )

    return {
        "case_id": case.case_id,
        "parameter": case.parameter,
        "role": case.role,
        "global_minimal_model": [str(value) for value in case.model],
        "integral_two_division_polynomial_for_zeta_equals_4theta": str(field_polynomial),
        "integral_two_division_polynomial_discriminant": str(polynomial_discriminant),
        "generic_point_count": len(generic_rows),
        "exceptional_point_count": len(exceptional_rows),
        "known_point_count": len(rows),
        "fingerprint_dimension": len(rows[0]) if rows else 0,
        "generic_kummer_rank": generic_rank,
        "known_kummer_rank": total_rank,
        "generic_mod_2_independence_certified": generic_rank == GENERIC_RANK,
        "all_known_mod_2_independence_certified": full_known_rank_certified,
        "known_exceptional_quotient_dimension": exceptional_quotient_dimension,
        "expected_exceptional_quotient_dimension": expected_exceptional_dimension,
        "exceptional_quotient_basis_certified": quotient_basis_certified,
        "known_realized_exceptional_quotient_class_count_including_zero": (
            2**expected_exceptional_dimension if quotient_basis_certified else None
        ),
        "selected_auxiliary_primes": selected_primes,
        "largest_selected_auxiliary_prime": (
            selected_primes[-1]["prime"] if selected_primes else None
        ),
        "points": point_records,
        "runtime_seconds": time.monotonic() - started,
        "status": (
            "PASS_EXACT_KNOWN_KUMMER_QUOTIENT"
            if full_known_rank_certified
            else "INCOMPLETE_AUXILIARY_FINGERPRINT"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--controls-only", action="store_true")
    parser.add_argument("--prime-bound", type=int, default=5000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.prime_bound < 3:
        parser.error("--prime-bound must be at least 3")
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(args.output)

    manifest = json.loads(args.manifest.read_text())
    if manifest.get("schema") != INPUT_SCHEMA:
        raise SystemExit("unexpected input manifest")
    selected = selected_manifest_cases(
        manifest, set(args.case), args.controls_only
    )
    authoritative = load_authoritative_cases()

    suite_started = time.monotonic()
    runs = []
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def write_output(status: str) -> None:
        passed = sum(run["status"].startswith("PASS") for run in runs)
        output = {
            "schema": OUTPUT_SCHEMA,
            "status": status,
            "input_manifest": {
                "path": str(args.manifest),
                "sha256": file_sha256(args.manifest),
            },
            "program": {
                "path": str(Path(__file__).resolve().relative_to(ROOT)),
                "sha256": file_sha256(Path(__file__)),
            },
            "method": (
                "exact odd-prime residue squareclasses of 4*x(P)-zeta, where "
                "zeta=4*theta, using squarefree factors of the integral "
                "two-division polynomial modulo each auxiliary prime"
            ),
            "parameters": {
                "prime_bound": args.prime_bound,
                "selected_case_ids": [row["case_id"] for row in selected],
            },
            "software": {
                "sage_version": str(sage_version),
                "license": "open_source",
            },
            "claim_boundary": [
                "Full row rank is an exact certificate that the supplied global Kummer classes are independent modulo squares.",
                "For controls, exceptional quotient coordinates are relative only to the known subgroup generated by MW17 and the supplied exceptional points.",
                "This audit does not compute or upper-bound the full 2-Selmer group, does not enumerate unknown quotient classes, and is not blind recovery of exceptional points.",
                "Failure to reach full rank before the prime bound would be an incomplete certificate, not a dependence result.",
            ],
            "run_count": len(runs),
            "selected_run_count": len(selected),
            "passed_count": passed,
            "runtime_seconds": time.monotonic() - suite_started,
            "runs": runs,
        }
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        temporary.replace(args.output)

    for manifest_case in selected:
        case_id = manifest_case["case_id"]
        print(f"{PROTOCOL}|case={case_id}|stage=audit|status=start", flush=True)
        record = audit_case(authoritative[case_id], args.prime_bound)
        runs.append(record)
        write_output("IN_PROGRESS_CHECKPOINT")
        print(
            f"{PROTOCOL}|case={case_id}|stage=audit|status={record['status']}"
            f"|generic_rank={record['generic_kummer_rank']}"
            f"|known_rank={record['known_kummer_rank']}"
            f"|largest_prime={record['largest_selected_auxiliary_prime']}"
            f"|seconds={record['runtime_seconds']:.6f}",
            flush=True,
        )

    passed = sum(run["status"].startswith("PASS") for run in runs)
    final_status = (
        "PASS_ALL_EXACT_KNOWN_KUMMER_QUOTIENTS"
        if passed == len(runs)
        else "INCOMPLETE_ONE_OR_MORE_KUMMER_FINGERPRINTS"
    )
    write_output(final_status)
    print(
        f"{PROTOCOL}|stage=complete|status={final_status}"
        f"|passed={passed}|selected={len(runs)}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
