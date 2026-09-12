#!/usr/bin/env python3
"""List or run a canonical H3 equation-lift success-path stage."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
LEDGER = Path(__file__).with_name("ledger.json")


def all_records(ledger: dict) -> list[dict]:
    return ledger["route"] + ledger["shortcut_audits"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", nargs="?")
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    ledger = json.loads(LEDGER.read_text())
    records = all_records(ledger)
    if args.list:
        for record in records:
            runnable = "run" if record.get("command") else "--"
            print(
                f"{record['id']:<42} {record['status']:<30} {runnable}"
            )
        return
    if not args.stage:
        parser.error("provide a stage id or use --list")

    matches = [record for record in records if record["id"] == args.stage]
    if len(matches) != 1:
        parser.error(f"unknown stage: {args.stage}")
    record = matches[0]
    command = record.get("command")
    if not command:
        raise SystemExit(
            f"{args.stage} has no successful canonical command yet "
            f"(status={record['status']})"
        )

    print("SUCCESS_PATH_COMMAND|" + shlex.join(command), flush=True)
    if args.dry_run:
        return
    subprocess.run(command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()

