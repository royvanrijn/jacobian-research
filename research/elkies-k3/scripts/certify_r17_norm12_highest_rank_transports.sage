#!/usr/bin/env sage-python
"""Certify the three highest-rank unresolved public-fibre transports.

Curves 11, 391, and 423 are the rank-at-least-28 rows whose displayed
quotients were still unknown in the 69-fibre calibration table.  This replay
specializes the saturated generic MW17 bases on the native 08234 and 07ca9
charts, recovers their coordinates in the pinned public point subgroups, and
verifies every relation by the exact elliptic-curve group law.

No cover inventory is evaluated here.  In particular, this certificate
computes the displayed quotient transports but makes no visibility claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import PolynomialRing, QQ
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SWEEP = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v1.json"
PUBLIC = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
HELPER = ROOT / "elkies-k3/scripts/certify_r17_norm12_native_icarm_quotient_audit.sage"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-highest-rank-transports-v1.json"

CHARTS = (
    {
        "source_chart": "norm12-orbit-08234",
        "representative": "norm12-orbit-08234",
        "curve_ids": (11, 423),
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08234-direct-fibration-v1.json",
    },
    {
        "source_chart": "norm12-orbit-07ca9",
        "representative": "norm12-orbit-07ca9",
        "curve_ids": (391,),
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit07ca9-direct-fibration-v1.json",
    },
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def build():
    helper = runpy.run_path(str(HELPER))
    special_fibre = helper["special_fibre"]
    sweep = json.loads(SWEEP.read_text())
    public = json.loads(PUBLIC.read_text())
    if public["status"] != "PASS_PINNED_PUBLIC_POINT_PROJECTION_FOR_69_RECOGNIZED_FIBRES":
        raise ArithmeticError("the pinned 69-fibre public projection changed")
    hit_records = {
        int(record["curve_id"]): record
        for record in sweep["rational_j_hits_and_twists"]
    }
    public_records = {int(record["id"]): record for record in public["records"]}
    ring = PolynomialRing(QQ, "u")
    fibres = []
    inputs = {
        relative(SWEEP): digest(SWEEP),
        relative(PUBLIC): digest(PUBLIC),
        relative(HELPER): digest(HELPER),
    }
    for config in CHARTS:
        direct_path = config["direct"]
        direct = json.loads(direct_path.read_text())
        if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
            raise ArithmeticError(f"{config['source_chart']} basis is not saturated")
        inputs[relative(direct_path)] = digest(direct_path)
        for curve_id in config["curve_ids"]:
            print(
                f"R17HIGHRANKTRANSPORT|curve={curve_id}|stage=exact_specialization",
                flush=True,
            )
            fibres.append(
                special_fibre(
                    config,
                    hit_records[curve_id],
                    public_records[curve_id],
                    direct,
                    None,
                    ring,
                )
            )

    if [row["curve_id"] for row in fibres] != [11, 423, 391]:
        raise ArithmeticError("highest-rank transport order changed")
    for row in fibres:
        if row["snapshot_rank_lower_bound"] != 28:
            raise ArithmeticError("a selected row is no longer rank at least 28")
        quotient = row["displayed_exceptional_quotient"]
        if quotient["free_rank"] != 11 or quotient["quotient"] != "Z^11":
            raise ArithmeticError("a highest-rank displayed quotient is not Z^11")
        if row["alternate_q80_cover_audit"]["status"] != "NOT_RUN_NO_FROZEN_NATIVE_COVER_INVENTORY":
            raise ArithmeticError("this transport-only replay opened a cover inventory")

    return {
        "schema": "elkies-k3.r17-norm12-highest-rank-transports.v1",
        "status": "PASS_EXACT_THREE_HIGHEST_RANK_DISPLAYED_QUOTIENT_TRANSPORTS",
        "summary": {
            "curve_ids": [11, 423, 391],
            "snapshot_rank_lower_bound_each": 28,
            "native_families": ["norm12-orbit-08234", "norm12-orbit-07ca9"],
            "displayed_quotient_each": "Z^11",
            "new_exact_quotient_rows": 3,
            "remaining_unknown_quotient_rows_in_69_fibre_table": 54,
        },
        "fibres": fibres,
        "claim_boundary": {
            "proved": [
                "independence of each displayed 28-point subgroup by exact finite reductions",
                "exact specialization of a saturated generic MW17 basis on each fibre",
                "primitivity of the generic subgroup in each displayed public subgroup",
                "the displayed-subgroup quotient Z^11 for curves 11, 391, and 423",
            ],
            "not_proved": [
                "that a displayed subgroup is the full Mordell-Weil group",
                "an exact-rank upper bound",
                "any native-cover visibility count for these fibres",
                "a quotient transport for any of the other 54 unresolved public fibres",
            ],
        },
        "inputs": inputs,
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact QQ elliptic-curve group law",
                "exact Smith normal form",
                "canonical heights used only to propose integer relations",
            ],
        },
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves/cas sage -python elkies-k3/scripts/"
            "certify_r17_norm12_highest_rank_transports.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    serialized = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored highest-rank transports differ from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "R17HIGHRANKTRANSPORT|curves=11,423,391|quotients=Z^11,Z^11,Z^11|"
        f"status=PASS|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
