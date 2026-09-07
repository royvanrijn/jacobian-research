#!/usr/bin/env sage-python
"""Representative-normalised, witness-hidden PointsQI benchmark for E302.

Controls are run first.  The missing R+H panel is unavailable unless a
nonzero H control is solved by the isolated solver.  The preparer knows the
retrospective public control words only to construct square-equivalent Kummer
representatives; it never writes a rational cover witness.  The solver process
receives a neutral token and two quadrics only.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path
import subprocess

from sage.all import QQ, EllipticCurve
from sage.version import version as SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
PROTOCOL = CAS / "CURVE302_STRICT_COVER_BENCHMARK_V2_PROTOCOL.json"
BASE_PATH = CAS / "benchmark_curve302_strict_covers_hidden.sage"
BASE = SourceFileLoader("curve302_hidden_cover_base", str(BASE_PATH)).load_module()
CONTROL_OUTPUT = ART / "curve302_strict_cover_control_benchmark_v2.json"
MISSING_OUTPUT = ART / "curve302_strict_cover_missing_panel_v2.json"
WORK = ROOT / "artifacts" / "local" / "elliptic-curves" / "curve302-strict-cover-benchmark-v2"


def read(path: Path):
    return json.loads(path.read_text())


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def public_sum(word, curve, public):
    return sum((bit * point for bit, point in zip(word, public)), curve(0))


def canonical_cases(cohort: str):
    """Build cover equations, but deliberately retain no cover witness."""
    panel = BASE.build_panels()
    selected = panel["cases"][:8] if cohort == "recovered_control" else panel["cases"][8:]
    E = EllipticCurve(QQ, list(map(QQ, BASE.curve.GENERAL_WEIERSTRASS_COEFFICIENTS)))
    public = [E(point) for point in BASE.curve.POINTS]
    ring = panel["ring"]
    u, v, w, z = ring.gens()
    cases = []
    for source in selected:
        word = source["public_word_mod_2"]
        point = public_sum(word, E, public)
        if point.is_zero():
            alpha, root = [QQ(1), QQ(0), QQ(0)], QQ(1)
            normalisation = "ZERO_CLASS"
        else:
            X, Y = 4 * QQ(point[0]), 8 * QQ(point[1]) + 4 * QQ(point[0]) + 4
            assert Y * Y == panel["f"](X)
            alpha, root = [X, QQ(-1), QQ(0)], Y
            normalisation = "NORMALIZED_ONE_FACTOR_KUMMER_REPRESENTATIVE"
        products = BASE.multiply_mod_cubic(
            alpha,
            BASE.multiply_mod_cubic([u, v, w], [u, v, w], panel["coefficients"]),
            panel["coefficients"],
        )
        first = BASE.primitive_integral(products[1] + z**2)
        second = BASE.primitive_integral(products[2])
        cases.append({
            "token": f"{cohort}-{source['member']:02d}",
            "cohort": cohort,
            "member": source["member"],
            "public_word_mod_2": word,
            "representative_normalisation": normalisation,
            "alpha_coefficients": [str(value) for value in alpha],
            "norm_square_root": str(root),
            "integral_quadrics": [str(first), str(second)],
        })
    assert len(cases) == 8
    return panel, cases


def sealed_job(case, protocol):
    return {
        "schema": "elliptic-curves.sealed-two-quadric-job.v1",
        "token": case["token"],
        "integral_quadrics": case["integral_quadrics"],
        "solver": {
            "engine": protocol["solver"]["engine"],
            "model_preparation": protocol["solver"]["model_preparation"],
            "pointsqi_height": protocol["solver"]["pointsqi_height"],
            "remote_request_timeout_seconds": protocol["solver"]["remote_request_timeout_seconds"],
        },
    }


def run_case(case, protocol, panel, stage):
    directory = WORK / stage
    directory.mkdir(parents=True, exist_ok=True)
    sealed = directory / f"{case['token']}.json"
    result = sealed.with_name(sealed.stem + "-result.json")
    job = sealed_job(case, protocol)
    if sealed.exists():
        assert read(sealed) == job
    else:
        BASE.put_new(sealed, job)
    if not result.exists():
        try:
            subprocess.run(
                ["sage", "-python", str(BASE_PATH), "--solver", "--sealed", str(sealed)],
                check=False, timeout=protocol["solver"]["subprocess_timeout_seconds"],
            )
        except subprocess.TimeoutExpired:
            BASE.put_new(result, {
                "schema": "elliptic-curves.sealed-two-quadric-solver-result.v1",
                "token": case["token"], "sealed_input_sha256": digest(sealed),
                "engine": protocol["solver"]["engine"], "status": "SUBPROCESS_TIMEOUT",
                "wall_seconds": protocol["solver"]["subprocess_timeout_seconds"],
            })
    solved = read(result)
    assert solved["token"] == case["token"] and solved["sealed_input_sha256"] == digest(sealed)
    checked = BASE.map_solution(case, solved, panel)
    return {
        **{key: case[key] for key in ("token", "cohort", "member", "public_word_mod_2", "representative_normalisation", "alpha_coefficients", "norm_square_root", "integral_quadrics")},
        "sealed_solver_input": str(sealed.relative_to(ROOT)),
        "sealed_solver_input_sha256": digest(sealed),
        "solver_result": str(result.relative_to(ROOT)),
        "solver_result_sha256": digest(result),
        **checked,
    }


def bindings():
    return {
        str(path.relative_to(ROOT)): digest(path)
        for path in (BASE.FILTER, BASE.SPAN / "result.json", BASE.PARENT, BASE.PUBLIC_SOURCE,
                     BASE.BUILDER, BASE_PATH, PROTOCOL, Path(__file__))
    }


def emit(output, stage, records, protocol):
    solved_nonzero = [row for row in records if row["member"] != 0 and row["exact_rational_cover_points"]]
    status = "PASS" if solved_nonzero else "PARTIAL_NO_ADVANCE"
    payload = {
        "schema": f"elliptic-curves.curve302-strict-cover-{stage}.v2",
        "status": status,
        "baseline": "M24", "software": {"sage": SAGE_VERSION},
        "protocol": str(PROTOCOL.relative_to(ROOT)), "protocol_sha256": digest(PROTOCOL),
        "bindings": bindings(), "cases": records,
        "comparison": {"representatives": 8, "nonzero_exact_cover_solution_count": len(solved_nonzero)},
        "witness_hiding": "The sealed job schema contains only token, quadrics and solver limits. It contains no alpha, public word, point coordinate, coset label or rational cover witness.",
        "boundary": protocol["boundaries"],
    }
    BASE.put_new(output, payload)
    print(stage, status, payload["comparison"], flush=True)


def controls():
    protocol = read(PROTOCOL)
    panel, cases = canonical_cases("recovered_control")
    records = [run_case(case, protocol, panel, "controls") for case in cases]
    emit(CONTROL_OUTPUT, "control-benchmark", records, protocol)


def missing():
    control = read(CONTROL_OUTPUT)
    assert control["status"] == "PASS", "the R+H panel is closed until controls pass"
    protocol = read(PROTOCOL)
    panel, cases = canonical_cases("missing_strict_coset")
    records = [run_case(case, protocol, panel, "missing") for case in cases]
    emit(MISSING_OUTPUT, "missing-panel", records, protocol)


def check(output, cohort, stage):
    stored = read(output)
    protocol = read(PROTOCOL)
    panel, cases = canonical_cases(cohort)
    assert stored["bindings"] == bindings() and stored["protocol_sha256"] == digest(PROTOCOL)
    assert len(stored["cases"]) == len(cases) == 8
    for case, row in zip(cases, stored["cases"]):
        assert all(row[key] == case[key] for key in case)
        sealed, result = ROOT / row["sealed_solver_input"], ROOT / row["solver_result"]
        assert digest(sealed) == row["sealed_solver_input_sha256"] and digest(result) == row["solver_result_sha256"]
        assert BASE.map_solution(case, read(result), panel) == {key: row[key] for key in ("solver_status", "exact_rational_cover_points")}
    assert stored["boundary"] == protocol["boundaries"]
    print("PASS", stage, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("controls", "missing", "check-controls", "check-missing"))
    args = parser.parse_args()
    if args.mode == "controls":
        controls()
    elif args.mode == "missing":
        missing()
    elif args.mode == "check-controls":
        check(CONTROL_OUTPUT, "recovered_control", "controls")
    else:
        check(MISSING_OUTPUT, "missing_strict_coset", "missing")


if __name__ == "__main__":
    main()
