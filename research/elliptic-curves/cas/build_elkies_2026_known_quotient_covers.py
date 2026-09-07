#!/usr/bin/env python3
"""Build explicit 2-covers for every known exceptional R17 quotient class.

For each rank-21 and rank-25--28 control, enumerate every nonzero subset of
the certified exceptional mod-2 basis.  Add the corresponding rational points
on the elliptic curve, map their sum to the monic completed-square cubic, and
use ``alpha = X(P)-theta`` to construct its intersection of two quadrics.
The representative has the exact rational cover point ``[1:0:0:1]``.

These are all classes in the *known exceptional subgroup*.  The script does
not assert that this subgroup is the full relative 2-Selmer quotient and does
not classify any genuinely unknown or unrealized Selmer class.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
DEFAULT_MANIFEST = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_relative_2selmer_suite_inputs_v1.json"
)
DEFAULT_KUMMER_AUDIT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_known_kummer_quotients_controls_v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_known_exceptional_quotient_covers_v1.json"
)
DEFAULT_FULL_OUTPUT = (
    ROOT
    / "artifacts/local/elliptic-curves"
    / "elkies-2026-known-exceptional-quotient-covers-v1"
    / "all-classes.json"
)
OUTPUT_SCHEMA = "elliptic-curves.elkies-2026-known-exceptional-quotient-covers.v1"
PROTOCOL = "ELKIESR17KNOWNCOVERS"

sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]

from build_bnf_free_two_covers import (  # noqa: E402
    cover_for,
    verify_rational_cover_witness,
)
from build_q12o5867_bnf_free_signature import (  # noqa: E402
    evaluate_cubic,
    monic_cubic_coefficients,
    point_on_monic_cubic,
)
from run_elkies_2026_relative_2selmer_open import (  # noqa: E402
    GENERIC_RANK,
    load_authoritative_cases,
)

from sage.all import EllipticCurve, PolynomialRing, QQ  # noqa: E402
from sage.version import version as sage_version  # noqa: E402


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rational(value: Fraction) -> QQ:
    return QQ(value.numerator) / QQ(value.denominator)


def build_case(case, audit_record, ring) -> dict:
    started = time.monotonic()
    exceptional_count = len(case.exceptional_points)
    if not (
        audit_record["generic_kummer_rank"] == GENERIC_RANK
        and audit_record["known_exceptional_quotient_dimension"]
        == exceptional_count
        and audit_record["exceptional_quotient_basis_certified"] is True
    ):
        raise ArithmeticError("the exact known-Kummer quotient audit is incomplete")

    model = tuple(Fraction(value) for value in case.model)
    coefficients_fraction = monic_cubic_coefficients(model)
    coefficients = [rational(value) for value in coefficients_fraction]
    curve = EllipticCurve(QQ, [QQ(value) for value in case.model])
    exceptional_points = [
        curve(QQ(x_coordinate), QQ(y_coordinate))
        for x_coordinate, y_coordinate in case.exceptional_points
    ]
    witness = [QQ(1), QQ(0), QQ(0), QQ(1)]
    classes = []
    for mask in range(1, 1 << exceptional_count):
        point = curve.zero()
        bits = []
        labels = []
        for index, exceptional_point in enumerate(exceptional_points):
            bit = (mask >> index) & 1
            bits.append(bit)
            if bit:
                point += exceptional_point
                labels.append(f"Q{index + 1}")
        if point.is_zero():
            raise ArithmeticError("a nonzero exceptional quotient mask summed to zero")
        source_point = [str(point[0]), str(point[1])]
        x_coordinate, z_coordinate = point_on_monic_cubic(model, source_point)
        norm = evaluate_cubic(coefficients_fraction, x_coordinate)
        if norm != z_coordinate**2:
            raise ArithmeticError("a quotient point failed completed-square transport")
        alpha = [rational(x_coordinate), QQ(-1), QQ(0)]
        cover = cover_for(alpha, coefficients, ring)
        affine_x = verify_rational_cover_witness(
            alpha, coefficients, witness, ring
        )
        if affine_x != rational(x_coordinate):
            raise ArithmeticError("the rational cover witness has the wrong x-map")
        classes.append(
            {
                "known_quotient_class_integer": mask,
                "known_quotient_bits": bits,
                "exceptional_point_labels_in_sum": labels,
                "elliptic_sum_point": source_point,
                "monic_cubic_point": [str(x_coordinate), str(z_coordinate)],
                "alpha_coefficients": [str(value) for value in alpha],
                "norm": str(norm),
                "norm_square_root": str(z_coordinate),
                "quadrics": cover,
                "rational_cover_witness": [str(value) for value in witness],
                "rational_cover_witness_verified": True,
                "known_exceptional_subgroup_realizes_class": True,
            }
        )

    expected = (1 << exceptional_count) - 1
    if len(classes) != expected:
        raise ArithmeticError("known quotient enumeration is incomplete")
    return {
        "case_id": case.case_id,
        "parameter": case.parameter,
        "generic_kummer_dimension": GENERIC_RANK,
        "known_exceptional_quotient_dimension": exceptional_count,
        "nonzero_known_exceptional_quotient_class_count": expected,
        "field_polynomial_ascending": [
            str(value) for value in coefficients_fraction
        ],
        "classes": classes,
        "runtime_seconds": time.monotonic() - started,
        "status": "PASS_ALL_KNOWN_EXCEPTIONAL_QUOTIENT_COVERS",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--kummer-audit", type=Path, default=DEFAULT_KUMMER_AUDIT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--full-output", type=Path, default=DEFAULT_FULL_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(args.output)
    if args.full_output.exists() and not args.overwrite:
        raise FileExistsError(args.full_output)

    audit = json.loads(args.kummer_audit.read_text())
    if audit.get("status") != "PASS_ALL_EXACT_KNOWN_KUMMER_QUOTIENTS":
        raise ValueError("the supplied Kummer audit is incomplete")
    audit_by_case = {record["case_id"]: record for record in audit["runs"]}
    authoritative = load_authoritative_cases()
    controls = [
        authoritative[case_id]
        for case_id in (
            "control-r21-t3_8",
            "control-r25",
            "control-r26",
            "control-r27",
            "control-r28",
        )
    ]
    ring = PolynomialRing(QQ, names=("u", "v", "w", "z"))
    started = time.monotonic()
    runs = []
    for case in controls:
        print(f"{PROTOCOL}|case={case.case_id}|stage=covers|status=start", flush=True)
        record = build_case(case, audit_by_case[case.case_id], ring)
        runs.append(record)
        print(
            f"{PROTOCOL}|case={case.case_id}|stage=covers|status=complete"
            f"|classes={record['nonzero_known_exceptional_quotient_class_count']}"
            f"|seconds={record['runtime_seconds']:.6f}",
            flush=True,
        )

    common = {
        "schema": OUTPUT_SCHEMA,
        "status": "PASS_ALL_CONTROL_KNOWN_EXCEPTIONAL_QUOTIENT_COVERS",
        "program": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": file_sha256(Path(__file__)),
        },
        "input_manifest": {
            "path": str(args.manifest),
            "sha256": file_sha256(args.manifest),
        },
        "known_kummer_audit": {
            "path": str(args.kummer_audit),
            "sha256": file_sha256(args.kummer_audit),
        },
        "software": {"sage_version": str(sage_version), "license": "open_source"},
        "descent_equation": (
            "X-theta = alpha*(u+v*theta+w*theta^2)^2; the two-cover is "
            "the theta and theta^2 coefficient equations"
        ),
        "total_nonzero_known_quotient_classes": sum(
            record["nonzero_known_exceptional_quotient_class_count"]
            for record in runs
        ),
        "runtime_seconds": time.monotonic() - started,
        "claim_boundary": [
            "Every stored cover is an exact intersection of quadrics with a verified rational point.",
            "Enumeration is exhaustive only inside each certified known exceptional quotient subgroup.",
            "The full relative 2-Selmer quotients remain unknown, so this artifact identifies no unknown or unrealized Selmer class.",
            "No bounded-search miss or Tate-Shafarevich conclusion is present.",
        ],
    }
    full_output = {**common, "storage": "full_all_class_records", "runs": runs}
    args.full_output.parent.mkdir(parents=True, exist_ok=True)
    args.full_output.write_text(
        json.dumps(full_output, indent=2, sort_keys=True) + "\n"
    )
    compact_runs = []
    for record in runs:
        class_text = json.dumps(
            record["classes"], sort_keys=True, separators=(",", ":")
        )
        compact_runs.append(
            {
                **{key: value for key, value in record.items() if key != "classes"},
                "all_class_records_sha256": sha256(class_text.encode()).hexdigest(),
                "basis_class_records": [
                    class_record
                    for class_record in record["classes"]
                    if class_record["known_quotient_class_integer"]
                    & (class_record["known_quotient_class_integer"] - 1)
                    == 0
                ],
            }
        )
    output = {
        **common,
        "storage": "compact_basis_records_plus_hashed_full_enumeration",
        "full_all_class_output": {
            "path": str(args.full_output),
            "sha256": file_sha256(args.full_output),
        },
        "runs": compact_runs,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|status={output['status']}"
        f"|classes={output['total_nonzero_known_quotient_classes']}"
        f"|seconds={output['runtime_seconds']:.6f}|output={args.output}"
        f"|full_output={args.full_output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
