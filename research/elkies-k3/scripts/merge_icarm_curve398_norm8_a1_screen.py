#!/usr/bin/env python3
"""Compact the checkpointed curve-398 norm-eight modular screen chain."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
TABLE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
SCREEN = ROOT / "elkies-k3/scripts/screen_icarm_curve398_norm8_a1_fibrations.sage"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-curve398-11952-norm8-a1-modular-screen-v1.json"
PRIMES = (
    1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061,
    1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123,
    1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213,
    1217, 1223, 1229, 1231,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-directory", type=Path, default=LOCAL)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    ledgers = []
    prior_path = None
    survivors = None
    exclusions = []
    for prime in PRIMES:
        path = args.local_directory / f"curve398-11952-norm8-a1-mod{prime}.json"
        document = json.loads(path.read_text())
        if document.get("status") != "PASS_COMPLETE_DECLARED_CHUNK_MODULAR_SCREEN":
            raise ArithmeticError(f"prime {prime}: incomplete modular ledger")
        search = document["search"]
        if search["prime"] != prime:
            raise ArithmeticError(f"prime {prime}: mislabeled ledger")
        processed = search["processed_priority_ranks"]
        if survivors is None:
            if processed != list(range(1, search["priority_table_class_count"] + 1)):
                raise ArithmeticError("first ledger is not the complete priority table")
        elif processed != survivors:
            raise ArithmeticError(f"prime {prime}: survivor chain is discontinuous")
        if prior_path is not None and search["candidate_rank_filter_sha256"] != digest(prior_path):
            raise ArithmeticError(f"prime {prime}: predecessor hash changed")
        excluded = [
            record["priority_rank"]
            for record in document["records"]
            if record["status"] == "PASS_MODULAR_NO_CURVE398_PARAMETER"
        ]
        exclusions.append({"prime": prime, "excluded_priority_ranks": excluded})
        survivors = search["survivor_priority_ranks"]
        ledgers.append(
            {
                "prime": prime,
                "processed_count": search["processed_count"],
                "excluded_count": len(excluded),
                "survivor_count": len(survivors),
                "checkpoint_sha256": digest(path),
            }
        )
        prior_path = path

    if survivors != [16875, 63669]:
        raise ArithmeticError(f"unexpected final survivor set: {survivors}")
    excluded_union = {rank for row in exclusions for rank in row["excluded_priority_ranks"]}
    if len(excluded_union) != 63915 or excluded_union | set(survivors) != set(range(1, 63918)):
        raise ArithmeticError("compact exclusions do not partition the complete layer")

    payload = {
        "schema": "elkies-k3.curve398-11952-norm8-a1-modular-screen.v1",
        "status": "PASS_COMPLETE_MODULAR_SCREEN_TWO_SURVIVORS",
        "source_chart": "norm12-orbit-11952",
        "target": "ICARM curve 398",
        "search": {
            "priority_table_class_count": 63917,
            "prime_chain": list(PRIMES),
            "checkpoint_summaries": ledgers,
            "excluded_count": len(excluded_union),
            "survivor_count": len(survivors),
            "survivor_priority_ranks": survivors,
            "exclusions_by_first_witness_prime": exclusions,
        },
        "proof_boundary": (
            "Each excluded rank has a displayed prime at which its exact projective curve-398 j-equation has no root. "
            "The two survivors remain candidates until characteristic-zero factorization; completeness applies only to the committed 63,917-class norm-eight layer on source chart 11952."
        ),
        "inputs": {relative(path): digest(path) for path in (MODEL, TABLE, SCREEN)},
        "reproducing_command": (
            "rerun the prime chain with screen_icarm_curve398_norm8_a1_fibrations.sage, then "
            "python3 elkies-k3/scripts/merge_icarm_curve398_norm8_a1_screen.py --check"
        ),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != rendered:
            raise ArithmeticError("stored compact curve-398 modular screen differs from checkpoints")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(
        f"CURVE398A1MERGE|classes=63917|excluded=63915|survivors={','.join(map(str, survivors))}|"
        f"status={payload['status']}|output={relative(args.output)}"
    )


if __name__ == "__main__":
    main()
