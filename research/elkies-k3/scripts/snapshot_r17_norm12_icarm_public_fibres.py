#!/usr/bin/env python3
"""Pin the public point data for the 69 recognized norm-12 ICARM fibres.

The exact six-class sweep deliberately stores only the arithmetic recognition
data.  This companion snapshot retains the smallest public-database projection
needed for quotient and local-arithmetic replays.  Point lists are truncated to
the rank lower bound pinned by the sweep, so later additions to the live ICARM
record cannot silently change this input.  The default check revalidates the
committed projection and its fixed digest offline; ``--live-source`` is the
explicit refresh/drift-audit path.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlopen


ROOT = Path(__file__).resolve().parents[2]
SWEEP = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
ENDPOINT = "https://elliptic-rank.icarm.cloud/database.json"
EXPECTED_PROJECTION_SHA256 = (
    "7b1e89c01812a04fd8bf8d6683e622e65487e54e08616f1e1e6f2e936649f7b2"
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_sha256(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def build(
    *,
    live_source: bool = False,
    database_path: Path | None = None,
    stored_projection: Path = OUTPUT,
) -> dict[str, object]:
    sweep = json.loads(SWEEP.read_text())
    hits = sweep["rational_j_hits_and_twists"]
    pinned = {int(record["curve_id"]): record for record in hits}
    if len(pinned) != 69:
        raise ArithmeticError("the exact sweep no longer contains 69 distinct hits")

    if database_path is not None:
        database = json.loads(database_path.read_text())
        source_records = database["curves"]
        offline = False
    elif live_source:
        with urlopen(ENDPOINT, timeout=60) as response:
            database = json.load(response)
        source_records = database["curves"]
        offline = False
    else:
        if not stored_projection.is_file():
            raise ArithmeticError(
                "offline replay requires the stored projection; use --database "
                "or --live-source only to reconstruct it"
            )
        stored = json.loads(stored_projection.read_text())
        if stored.get("schema") != "elkies-k3.r17-norm12-icarm-public-fibres.v1":
            raise ArithmeticError("stored public-fibre projection has an unknown schema")
        source_records = stored.get("records")
        if not isinstance(source_records, list):
            raise ArithmeticError("stored public-fibre projection is malformed")
        offline = True
    live = {int(record["id"]): record for record in source_records}

    records: list[dict[str, object]] = []
    for curve_id in sorted(pinned):
        hit = pinned[curve_id]
        record = live.get(curve_id)
        if record is None:
            raise ArithmeticError(f"ICARM curve {curve_id} disappeared")
        rank = int(hit["snapshot_rank_lower_bound"])
        if offline:
            if (
                int(record["snapshot_rank_lower_bound"]) != rank
                or len(record["points"]) != rank
            ):
                raise ArithmeticError(
                    f"stored ICARM curve {curve_id} has the wrong pinned point prefix"
                )
        elif int(record["rank_lower_bound"]) < rank or len(record["points"]) < rank:
            raise ArithmeticError(f"ICARM curve {curve_id} lost pinned point data")
        records.append(
            {
                "id": curve_id,
                "representative": hit["representative"],
                "representative_frame_class": hit["representative_frame_class"],
                "representative_parameter": hit["representative_parameter"],
                "snapshot_rank_lower_bound": rank,
                "curve_key": record["curve_key"],
                "ainvs": record["ainvs"],
                "points": record["points"][:rank],
                "bad_primes": record["bad_primes"],
                "conductor": record["conductor"],
                "discriminant": record["discriminant"],
                "created_at": record["created_at"],
            }
        )

    source_projection = [
        {
            key: record[key]
            for key in (
                "id",
                "curve_key",
                "ainvs",
                "points",
                "snapshot_rank_lower_bound",
                "bad_primes",
                "conductor",
                "discriminant",
                "created_at",
            )
        }
        for record in records
    ]
    projection_sha256 = canonical_sha256(source_projection)
    # This digest pins the 2026-09-04 projection used by the audit.  A mismatch
    # is a source-data change, not something the replay may accept implicitly.
    if projection_sha256 != EXPECTED_PROJECTION_SHA256:
        raise ArithmeticError(
            "the recognized public-fibre projection changed: " + projection_sha256
        )
    return {
        "schema": "elkies-k3.r17-norm12-icarm-public-fibres.v1",
        "status": "PASS_PINNED_PUBLIC_POINT_PROJECTION_FOR_69_RECOGNIZED_FIBRES",
        "source": {
            "url": ENDPOINT,
            "live_database_curve_count_observed_at_snapshot": 556,
            "projection_sha256": projection_sha256,
            "projection_rule": (
                "For every curve recognized by the exact six-class sweep, retain "
                "the listed fields and the first snapshot_rank_lower_bound points."
            ),
        },
        "summary": {
            "recognized_fibres": len(records),
            "displayed_independent_points": sum(len(record["points"]) for record in records),
            "all_records_have_the_pinned_number_of_points": True,
        },
        "records": records,
        "inputs": {relative(SWEEP): digest(SWEEP)},
        "claim_boundary": {
            "proved": [
                "the live public endpoint supplied the pinned point prefix for all 69 exact-sweep hits"
            ],
            "not_proved": [
                "that any displayed subgroup is the full Mordell--Weil group",
                "that the public rank lower bounds are exact ranks",
            ],
        },
        "reproducing_command": (
            ".venv/bin/python elkies-k3/scripts/"
            "snapshot_r17_norm12_icarm_public_fibres.py"
        ),
        "software_assumptions": {"python": sys.version.split()[0]},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--database", type=Path)
    parser.add_argument(
        "--live-source",
        action="store_true",
        help="refresh or audit against the current mutable public endpoint",
    )
    args = parser.parse_args()
    if args.database is not None and args.live_source:
        parser.error("--database and --live-source are mutually exclusive")
    output = args.output.resolve()
    payload = json.dumps(
        build(
            live_source=args.live_source,
            database_path=args.database.resolve() if args.database else None,
            stored_projection=output,
        ),
        indent=2,
        sort_keys=True,
    ) + "\n"
    if args.check:
        if not output.exists() or output.read_text() != payload:
            raise ArithmeticError("stored public-fibre projection differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(payload)
    print(
        "R17ICARMSNAPSHOT|fibres=69|points=1545|projection_sha256="
        "7b1e89c01812a04fd8bf8d6683e622e65487e54e08616f1e1e6f2e936649f7b2|"
        f"output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
