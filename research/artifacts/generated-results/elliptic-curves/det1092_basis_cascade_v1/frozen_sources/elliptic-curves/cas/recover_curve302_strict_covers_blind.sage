#!/usr/bin/env sage-python
"""Blind, checkpointed point recovery for the eight strict degree-four controls.

The ``--worker`` entry point is the arithmetic boundary: it reads one sealed
JSON file and is intentionally unable to import the curve, its displayed
points, Kummer words, alpha, or a control witness.  It only verifies and
searches the two supplied quadrics.  ``--capture`` prepares the fixed eight
post-hoc controls, starts one worker per cover, and only then maps any returned
cover point back to E302.

This is a calibration experiment.  A miss is a bounded miss, never a
solubility, Sha, or rank conclusion.
"""

from __future__ import annotations

import argparse
import importlib.machinery
import json
from hashlib import sha256
from math import gcd
from pathlib import Path
import subprocess
import sys
import time

from sage.all import QQ, ZZ, PolynomialRing, matrix


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
SOURCE = CAS / "benchmark_curve302_strict_covers_hidden.sage"
OUTPUT = ART / "curve302_strict_cover_blind_recovery_v2.json"
WORK = ROOT / "artifacts" / "local" / "elliptic-curves" / "curve302-strict-cover-blind-recovery-v2"
POLICY = {
    "fixed_primary_tokens": [f"case-{index:02d}" for index in range(8, 16)],
    "bounded_projective_slice_height": 300,
    "worker_timeout_seconds": 20,
    "selection": "all eight strict-coset controls; no witness-dependent ordering or stopping",
}


def read(path: Path):
    return json.loads(path.read_text())


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def put_new(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def source_module():
    return importlib.machinery.SourceFileLoader("curve302_hidden_source", str(SOURCE)).load_module()


def primitive(values):
    answer = 0
    for value in values:
        answer = gcd(answer, abs(int(value)))
    return answer == 1


def coefficient_bits(quadrics):
    ring = PolynomialRing(QQ, names=("u", "v", "w", "z"))
    return max(
        abs(int(coefficient)).bit_length()
        for text in quadrics
        for coefficient in ring(text).coefficients()
        if coefficient
    )


def check_and_slice(quadrics, bound):
    """Exact raw-coordinate search, with every coordinate used as a pivot.

    The strict-cover construction has Q1=A(u,v,w)+z^2 and Q2=B(u,v,w).
    For each choice of the eliminated coordinate of B, enumerate the other
    two integral coordinates and test the quadratic discriminant, then test
    whether -A is a square.  This is deliberately not a heuristic or a
    witness-guided slice: it exhausts this declared primitive box.
    """
    ring = PolynomialRing(QQ, names=("u", "v", "w", "z"))
    u, v, w, z = ring.gens()
    first, second = [ring(text) for text in quadrics]
    assert first.monomial_coefficient(z * z) == 1
    assert all(first.monomial_coefficient(z * x) == 0 for x in (u, v, w))
    assert all(second.monomial_coefficient(z * x) == 0 for x in (u, v, w, z))
    xyz = (u, v, w)
    seen, hits, conic_hits, tested = set(), [], 0, 0
    for eliminated in xyz:
        others = tuple(value for value in xyz if value != eliminated)
        leading = ZZ(second.monomial_coefficient(eliminated * eliminated))
        assert leading
        cross = [ZZ(second.monomial_coefficient(eliminated * value)) for value in others]
        diagonal = [ZZ(second.monomial_coefficient(value * value)) for value in others]
        mixed = ZZ(second.monomial_coefficient(others[0] * others[1]))
        for left in range(-bound, bound + 1):
            for right in range(-bound, bound + 1):
                if not left and not right:
                    continue
                tested += 1
                linear = cross[0] * left + cross[1] * right
                constant = diagonal[0] * left * left + mixed * left * right + diagonal[1] * right * right
                discriminant = linear * linear - 4 * leading * constant
                if discriminant < 0 or not discriminant.is_square():
                    continue
                root = discriminant.sqrt()
                for numerator in (-linear - root, -linear + root):
                    denominator = 2 * leading
                    if numerator % denominator:
                        continue
                    values = {others[0]: ZZ(left), others[1]: ZZ(right), eliminated: numerator // denominator}
                    point3 = tuple(values[value] for value in xyz)
                    if not primitive(point3):
                        continue
                    conic_hits += 1
                    square = -ZZ(first(*(point3 + (0,))))
                    if square < 0 or not square.is_square():
                        continue
                    ordinate = square.sqrt()
                    for signed in ({ordinate} if not ordinate else {ordinate, -ordinate}):
                        point = point3 + (signed,)
                        assert first(*point) == second(*point) == 0
                        if point not in seen:
                            seen.add(point)
                            hits.append([str(value) for value in point])
    return {
        "algorithm": "all_three_raw_ternary_quadratic_slices",
        "coordinate_height": bound,
        "primitive_pairs_tested": tested,
        "conic_points_in_box": conic_hits,
        "exact_cover_points": hits,
        "completion": "COMPLETE_IN_DECLARED_BOX",
    }


def worker(sealed: Path) -> None:
    """Run without loading the curve-302 source or the benchmark panel."""
    job = read(sealed)
    assert set(job) == {"schema", "token", "integral_quadrics", "policy"}
    assert job["schema"] == "elliptic-curves.sealed-strict-cover-recovery-job.v2"
    assert set(job["policy"]) == {"bounded_projective_slice_height"}
    assert len(job["integral_quadrics"]) == 2
    started = time.monotonic()
    result = {
        "schema": "elliptic-curves.sealed-strict-cover-recovery-result.v2",
        "token": job["token"],
        "sealed_input_sha256": digest(sealed),
        "status": "COMPLETED",
        "input_coefficient_bits": coefficient_bits(job["integral_quadrics"]),
        "methods": [],
    }
    try:
        result["methods"].append(
            check_and_slice(job["integral_quadrics"], int(job["policy"]["bounded_projective_slice_height"]))
        )
    except Exception as error:
        result["status"] = "WORKER_ERROR"
        result["error"] = repr(error)
    result["wall_seconds"] = time.monotonic() - started
    put_new(sealed.with_name(sealed.stem + "-result.json"), result)
    print(job["token"], result["status"], flush=True)


def sealed_case(case):
    return {
        "schema": "elliptic-curves.sealed-strict-cover-recovery-job.v2",
        "token": case["token"],
        "integral_quadrics": case["integral_quadrics"],
        "policy": {"bounded_projective_slice_height": POLICY["bounded_projective_slice_height"]},
    }


def capture():
    source = source_module()
    panel = source.build_panels()
    cases = [case for case in panel["cases"] if case["token"] in POLICY["fixed_primary_tokens"]]
    assert [case["token"] for case in cases] == POLICY["fixed_primary_tokens"]
    WORK.mkdir(parents=True, exist_ok=True)
    rows = []
    for case in cases:
        sealed = WORK / f"{case['token']}.json"
        result_path = sealed.with_name(sealed.stem + "-result.json")
        expected = sealed_case(case)
        if sealed.exists():
            assert read(sealed) == expected
        else:
            put_new(sealed, expected)
        if not result_path.exists():
            try:
                subprocess.run(
                    ["sage", "-python", str(Path(__file__)), "--worker", "--sealed", str(sealed)],
                    check=False,
                    timeout=POLICY["worker_timeout_seconds"],
                )
            except subprocess.TimeoutExpired:
                put_new(result_path, {
                    "schema": "elliptic-curves.sealed-strict-cover-recovery-result.v2",
                    "token": case["token"], "sealed_input_sha256": digest(sealed),
                    "status": "WORKER_TIMEOUT", "methods": [],
                })
        result = read(result_path)
        assert result["token"] == case["token"] and result["sealed_input_sha256"] == digest(sealed)
        points = [point for method in result.get("methods", []) for point in method.get("exact_cover_points", [])]
        verified = []
        ring = panel["ring"]
        first, second = [ring(text) for text in case["integral_quadrics"]]
        for point in points:
            vector = [QQ(value) for value in point]
            assert first(*vector) == second(*vector) == 0
            verified.append(point)
        rows.append({
            "token": case["token"], "sealed_input": str(sealed.relative_to(ROOT)),
            "sealed_input_sha256": digest(sealed), "result": str(result_path.relative_to(ROOT)),
            "result_sha256": digest(result_path), "worker_status": result["status"],
            "input_coefficient_bits": result.get("input_coefficient_bits"),
            "exact_blind_cover_points": verified,
        })
    found = [row for row in rows if row["exact_blind_cover_points"]]
    output = {
        "schema": "elliptic-curves.curve302-strict-cover-blind-recovery.v2",
        "status": "POINT_RECOVERED" if found else "BOUNDED_NO_POINT",
        "protocol": "eight post-hoc known-soluble strict-coset controls; workers receive quadrics only",
        "policy": POLICY,
        "bindings": {str(path.relative_to(ROOT)): digest(path) for path in (SOURCE, Path(__file__))},
        "cases": rows,
        "recovered_control_count": len(found),
        "boundary": "A completed raw-coordinate box miss says nothing about global solubility or rank.",
    }
    put_new(OUTPUT, output)
    print(output["status"], len(found), flush=True)


def check():
    stored = read(OUTPUT)
    assert stored["bindings"] == {str(path.relative_to(ROOT)): digest(path) for path in (SOURCE, Path(__file__))}
    assert stored["policy"] == POLICY and len(stored["cases"]) == 8
    for row in stored["cases"]:
        sealed = ROOT / row["sealed_input"]
        result = ROOT / row["result"]
        assert digest(sealed) == row["sealed_input_sha256"]
        assert digest(result) == row["result_sha256"]
    print("PASS exact replay of blind strict-cover bounded worker", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capture", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--sealed", type=Path)
    args = parser.parse_args()
    assert sum((args.capture, args.check, args.worker)) == 1
    if args.worker:
        assert args.sealed is not None
        worker(args.sealed)
    elif args.capture:
        capture()
    else:
        check()


if __name__ == "__main__":
    main()
