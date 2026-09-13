#!/usr/bin/env python3
"""Classify registered software locks without replaying research computations.

``MATH_STATUS.json`` remains the authority.  This audit makes the operational
role of every registered path visible: source program, raw replay input,
compact generated certificate, optional local checkpoint, or supporting
record.  The categories are path roles only; they do not establish that an
input is sufficient for a replay or that a certificate proves a theorem.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).resolve().parents[1]
PROGRAM_SUFFIXES = {
    ".gp", ".ipynb", ".jl", ".m", ".pari", ".py", ".r", ".sage", ".sh",
}
ROLE_ORDER = (
    "source-program",
    "raw-replay-input",
    "compact-certificate",
    "local-checkpoint",
    "supporting-record",
)


def role(lock: str) -> str:
    """Return the registry-path role for one software-lock value."""
    path = PurePosixPath(lock)
    parts = path.parts
    if parts[:2] == ("artifacts", "generated-results"):
        return "compact-certificate"
    if parts[:2] == ("artifacts", "local"):
        return "local-checkpoint"
    if "data" in parts:
        return "raw-replay-input"
    if path.suffix.lower() in PROGRAM_SUFFIXES:
        return "source-program"
    return "supporting-record"


def rows() -> list[dict[str, str]]:
    status = json.loads((ROOT / "MATH_STATUS.json").read_text())
    return [
        {"claim": entry["id"], "lock": lock, "role": role(lock)}
        for entry in status["entries"]
        for lock in entry["software_lock"]
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="emit individual lock rows")
    args = parser.parse_args()
    audit_rows = rows()
    unique = {row["lock"]: row["role"] for row in audit_rows}
    if len(unique) != len({row["lock"] for row in audit_rows}):
        raise SystemExit("software-lock role conflict for a shared path")
    missing = [lock for lock in unique if not (ROOT / lock).is_file()]
    if missing:
        raise SystemExit("missing software locks:\n" + "\n".join(sorted(missing)))
    counts = Counter(unique.values())
    summary = {
        "references": len(audit_rows),
        "unique_paths": len(unique),
        "unique_roles": {name: counts[name] for name in ROLE_ORDER},
    }
    if args.json:
        print(json.dumps({"summary": summary, "rows": audit_rows}, indent=2))
        return
    categories = "; ".join(f"{counts[name]} {name}" for name in ROLE_ORDER)
    print(
        f"PASS software-lock roles: {len(unique)} unique paths across "
        f"{len(audit_rows)} references ({categories}); all registered locks exist"
    )


if __name__ == "__main__":
    main()
