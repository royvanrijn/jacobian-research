#!/usr/bin/env sage-python
"""Certify the completed and timed-out 2-descent jobs for four 074d9 twists."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
LOCAL_ROOT = ROOT / "artifacts/local/elkies-k3/r17-074d9-twist-2descent"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-twist-2descent-audit-v1.json"
)
JOBS = (
    ("074d9-orbit-04b07", "04b07", 19),
    ("074d9-orbit-11a44", "11a44", 19),
    ("074d9-orbit-11279", "11279", 19),
    ("074d9-orbit-080fa", "080fa", 31),
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def parse_output(path: Path, label: str, prime: int):
    root = ET.fromstring(path.read_text())
    headers = root.find("headers")
    if headers is None:
        raise ValueError(f"{path}: missing Magma headers")
    lines = [element.text or "" for element in root.findall("./results/line")]
    marker = f"R17074D9_2DESCENT|label={label}|prime={prime}"
    if marker not in lines:
        raise ValueError(f"{path}: output marker mismatch")
    warning = headers.findtext("warning")
    fields = {
        "magma_version": headers.findtext("version"),
        "server_max_time_seconds": int(headers.findtext("max_time")),
        "server_seed": int(headers.findtext("seed")),
        "server_time_seconds": (
            None if headers.findtext("time") is None else float(headers.findtext("time"))
        ),
        "warning": warning,
        "lines": lines,
    }
    factor_line = next(
        (line for line in lines if line.startswith("TWO_DIVISION_FACTOR_DEGREES")),
        None,
    )
    if factor_line != "TWO_DIVISION_FACTOR_DEGREES [ 3 ]":
        raise ArithmeticError(f"{path}: 2-division cubic is not certified irreducible")
    upper_line = next(
        (line for line in lines if line.startswith("MW_RANK_UPPER_FROM_TWO_SELMER")),
        None,
    )
    if upper_line is None:
        if warning is None or "time limit" not in warning:
            raise ArithmeticError(f"{path}: incomplete output without timeout warning")
        fields["outcome"] = "TIMEOUT_AFTER_IRREDUCIBLE_TWO_DIVISION_CUBIC"
        fields["MW_rank_upper_bound_over_displayed_finite_function_field"] = None
    else:
        fields["outcome"] = "PASS_COMPLETED_TWO_SELMER_GROUP"
        fields["MW_rank_upper_bound_over_displayed_finite_function_field"] = int(
            upper_line.rsplit(" ", 1)[1]
        )
        invariants_line = next(
            line for line in lines if line.startswith("TWO_SELMER_INVARIANTS")
        )
        fields["two_selmer_dimension"] = invariants_line.count("2")
    return fields


def build_payload():
    rows = []
    inputs = {}
    for label, tag, discovery_prime in JOBS:
        attempts = []
        for prime, role in ((discovery_prime, "discovery_only"), (131, "good_reduction")):
            directory = LOCAL_ROOT / tag / f"p{prime}"
            certificate_path = directory / "input-certificate.json"
            job_path = directory / "two-descent.m"
            output_path = directory / "two-descent.output.xml"
            certificate = json.loads(certificate_path.read_text())
            if certificate["label"] != label or int(certificate["prime"]) != prime:
                raise ValueError("stale 2-descent input certificate")
            if certificate["magma_job"]["sha256"] != digest(job_path):
                raise ArithmeticError("2-descent job hash mismatch")
            if role == "good_reduction" and not certificate["good_reduction"]:
                raise ArithmeticError("p=131 lost good-reduction status")
            if role == "discovery_only" and certificate["good_reduction"]:
                raise ArithmeticError("discovery prime unexpectedly became globally good")
            parsed = parse_output(output_path, label, prime)
            attempts.append(
                {
                    "prime": prime,
                    "role": role,
                    "good_reduction": bool(certificate["good_reduction"]),
                    **parsed,
                    "input_certificate": {
                        "path": relative(certificate_path),
                        "sha256": digest(certificate_path),
                    },
                    "magma_job": {"path": relative(job_path), "sha256": digest(job_path)},
                    "raw_output": {
                        "path": relative(output_path),
                        "sha256": digest(output_path),
                    },
                }
            )
            for path in (certificate_path, job_path, output_path):
                inputs[relative(path)] = digest(path)
        rows.append(
            {
                "label": label,
                "characteristic_zero_function_field_rank_status": "UNKNOWN",
                "attempts": attempts,
            }
        )
    return {
        "schema": "elkies-k3.r17-074d9-twist-2descent-audit.v1",
        "status": "INCOMPLETE_GOOD_REDUCTION_2DESCENTS_TIMED_OUT",
        "claim": (
            "All four finite-function-field 2-descents complete at smaller "
            "discovery primes, while all four globally good p=131 attempts exceed "
            "the public Magma server's 60-second limit."
        ),
        "twists": rows,
        "proof_boundary": (
            "The nongood-prime Selmer dimensions are discovery diagnostics and do "
            "not bound the characteristic-zero ranks. The p=131 jobs establish no "
            "Selmer bound because they timed out. Every characteristic-zero rank "
            "therefore remains UNKNOWN."
        ),
        "generation": {
            "command": (
                "sage -python elkies-k3/scripts/"
                "certify_r17_074d9_twist_2descent_audit.sage"
            ),
            "magma_calculator_endpoint": (
                "https://magma.maths.usyd.edu.au/xml/calculator.xml"
            ),
        },
        "inputs": inputs,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("stale 074d9 2-descent audit")
        terminal = "PASS"
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        terminal = "WROTE"
    print(
        f"R17074D92DESCENTAUDIT|twists={len(payload['twists'])}"
        f"|status={terminal}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
