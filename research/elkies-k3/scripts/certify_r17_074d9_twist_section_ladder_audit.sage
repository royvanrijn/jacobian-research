#!/usr/bin/env sage-python
"""Index and check the bounded P.O=0,1,2 section-ladder campaign."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
LOCAL_ROOT = ROOT / "artifacts/local/elkies-k3/r17-074d9-twist-section-ladder"
SECTIONS = ROOT / "artifacts/generated-results/elkies-k3-r17-074d9-record-twist-sections-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-074d9-twist-section-ladder-audit-v1.json"
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


def transformed_known_section(section, export):
    prime = int(export["prime"])
    field = GF(prime)
    rational_ring = PolynomialRing(QQ, "t")
    t = rational_ring.gen()
    X = rational_ring([QQ(value) for value in section["X_coefficients_low_to_high"]])
    Y = rational_ring([QQ(value) for value in section["Y_coefficients_low_to_high"]])
    chart = export["charts"][0]["original_fibre"]
    if chart == "original_infinity":
        X_chart, Y_chart = X, Y
    elif chart.startswith("t="):
        c = QQ(int(chart[2:]))
        X_chart = rational_ring(sum(X[index] * (c*t+1)**index * t**(6-index) for index in range(7)))
        Y_chart = rational_ring(sum(Y[index] * (c*t+1)**index * t**(9-index) for index in range(10)))
    else:
        raise ValueError("unrecognized section chart")
    X_mod = [int(field(value)) for value in X_chart.list()]
    Y_mod = [int(field(value)) for value in Y_chart.list()]
    X_mod += [0] * (7 - len(X_mod))
    Y_mod += [0] * (10 - len(Y_mod))
    leading = [X_mod[6], Y_mod[9]]
    block = next(
        row for row in export["systems"]
        if row["chart_index"] == 0 and row["leading_x_y"] == leading
    )
    return {
        "chart": chart,
        "block_index": int(block["block_index"]),
        "leading_x_y": leading,
        "X_coefficients_low_to_high_mod_p": X_mod,
        "Y_coefficients_low_to_high_mod_p": Y_mod,
    }


def build_payload():
    section_data = json.loads(SECTIONS.read_text())
    section_by_label = {row["label"]: row for row in section_data["records"]}
    rows = []
    inputs = {relative(SECTIONS): digest(SECTIONS)}
    for label, tag, prime in JOBS:
        levels = []
        po0_export = None
        for intersection in (0, 1, 2):
            export_path = LOCAL_ROOT / tag / f"p{prime}" / f"intersection-{intersection}" / "export.json"
            export = json.loads(export_path.read_text())
            if (
                export.get("status") != "PASS_EXACT_COMPLETE_MODP_DISCOVERY_SYSTEM_EXPORT"
                or export["label"] != label
                or int(export["prime"]) != prime
                or int(export["intersection_P_dot_O"]) != intersection
            ):
                raise ValueError("stale section-ladder export")
            for system in export["systems"]:
                system_path = ROOT / system["path"]
                if digest(system_path) != system["sha256"]:
                    raise ArithmeticError("section-ladder system hash mismatch")
            inputs[relative(export_path)] = digest(export_path)
            levels.append(
                {
                    "P_dot_O": intersection,
                    "chart_count": int(export["chart_count"]),
                    "block_count": len(export["systems"]),
                    "variable_count_per_block": int(export["variable_count_per_block"]),
                    "equation_count_per_block": int(export["equation_count_per_block"]),
                    "solver_status": "NOT_RUN_TO_COMPLETION",
                    "export": {"path": relative(export_path), "sha256": digest(export_path)},
                }
            )
            if intersection == 0:
                po0_export = export
        rows.append(
            {
                "label": label,
                "prime": prime,
                "role": "discovery_only",
                "known_section_reduction": transformed_known_section(
                    section_by_label[label], po0_export
                ),
                "levels": levels,
            }
        )

    solved_summary_path = LOCAL_ROOT / "msolve-04b07-po-0.json"
    solved_summary = json.loads(solved_summary_path.read_text())
    if solved_summary.get("status") != "PASS_ALL_LEADING_IDEALS_CLASSIFIED":
        raise ValueError("04b07 P.O=0 solve is incomplete")
    if solved_summary["counts"] != {
        "ERROR": 0,
        "NONUNIT_LEADING_IDEAL": 8,
        "TIMEOUT": 0,
        "UNIT_IDEAL": 4,
    }:
        raise ArithmeticError("04b07 P.O=0 solver counts changed")
    inputs[relative(solved_summary_path)] = digest(solved_summary_path)
    rows[0]["levels"][0]["solver_status"] = "PASS_ALL_LEADING_IDEALS_CLASSIFIED"
    rows[0]["levels"][0]["solver_counts"] = solved_summary["counts"]
    rows[0]["levels"][0]["solver_summary"] = {
        "path": relative(solved_summary_path),
        "sha256": digest(solved_summary_path),
    }

    return {
        "schema": "elkies-k3.r17-074d9-twist-section-ladder-audit.v1",
        "status": "INCOMPLETE_BOUNDED_SECTION_LADDER_EXPORTED_THROUGH_P_DOT_O_2",
        "claim": (
            "Exact complete mod-p discovery systems are exported at P.O=0,1,2 "
            "for all four twists. Only the 04b07 P.O=0 leading ideals have been "
            "completely classified; no new characteristic-zero section is certified."
        ),
        "twists": rows,
        "proof_boundary": (
            "The small primes are not globally good surface reductions. The exports "
            "are exact coefficient-wise discovery schemes, not rank bounds. Unsolved "
            "levels and nonunit ideals cannot be interpreted as absence or existence "
            "of characteristic-zero sections."
        ),
        "inputs": inputs,
        "generation": {
            "command": (
                "sage -python elkies-k3/scripts/"
                "certify_r17_074d9_twist_section_ladder_audit.sage"
            )
        },
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
            raise SystemExit("stale twist section-ladder audit")
        terminal = "PASS"
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        terminal = "WROTE"
    print(
        f"R17074D9SECTIONLADDERAUDIT|twists={len(payload['twists'])}"
        f"|status={terminal}|output={args.output}", flush=True
    )


if __name__ == "__main__":
    main()
