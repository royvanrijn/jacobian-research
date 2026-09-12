#!/usr/bin/env python3
"""Check preserved navigation bytes and the content-preserving replay relocation."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "archive/repository-cleanup-2026-09-12"


def audit() -> None:
    manifest = json.loads((SNAPSHOT / "MANIFEST.json").read_text())
    records = manifest["files"]
    assert len({r["original_path"] for r in records}) == len(records)
    assert len({r["preserved_path"] for r in records}) == len(records)
    for record in records:
        path = ROOT / record["preserved_path"]
        assert path.is_file(), f"missing preserved document: {path}"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"], (
            f"preserved navigation bytes changed: {path}"
        )
    old_entry = SNAPSHOT / "OP-EC-NEXT.before.json"
    assert hashlib.sha256(old_entry.read_bytes()).hexdigest() == manifest["previous_scope_sha256"]
    for baseline in manifest.get("partial_review_baselines", []):
        path = ROOT / baseline["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == baseline["sha256"], (
            f"partial-review status baseline changed: {path}"
        )

    record = next(r for r in records if r["original_path"] == "research/REPRODUCE.md")
    original = (ROOT / record["preserved_path"]).read_text()
    destination = ROOT / "replay/CATALOGUE.md"

    def rebase(match):
        target = match.group(1)
        if "://" in target or target.startswith(("#", "mailto:")):
            return match.group(0)
        if "." not in target and "/" not in target:
            return match.group(0)  # e.g. mathematical R[... ](-3), not a path
        path, sep, anchor = target.partition("#")
        return "](" + os.path.relpath((ROOT / path).resolve(), destination.parent) + (sep + anchor if sep else "") + ")"

    expected = re.sub(r"\]\(([^)]+)\)", rebase, original)
    notice = (
        "> Full command reference, moved from `REPRODUCE.md` on 2026-09-12. Commands still run from `research/`.\n"
        "> Some sections describe superseded campaigns; consult the [short replay guide](../REPRODUCE.md),\n"
        "> [current status](../STATUS.md) and [algorithmic lessons](../knowledge/ALGORITHMS.md) before executing one.\n\n"
    )
    expected = expected.replace("# Reproducing the results\n\n", "# Full replay catalogue\n\n" + notice, 1)
    # Current consumer markers may be refreshed after a scope review, while the
    # preserved command/proof text remains fixed. audit_status validates markers.
    def without_markers(text):
        return re.sub(r"<!-- status-consumer: [^>]+ -->", "", text)

    assert without_markers(destination.read_text()) == without_markers(expected), (
        "replay relocation changed more than its heading, notice and relative links; "
        "preserve the original catalogue and record new commands in canonical notes"
    )
    def code_blocks(text):
        blocks, current, fence = [], [], None
        for line in text.splitlines(keepends=True):
            match = re.match(r"^(`{3,}|~{3,})", line)
            if match and fence is None:
                fence, current = match[1][0], [line]
            elif match and match[1][0] == fence:
                current.append(line)
                blocks.append("".join(current))
                current, fence = [], None
            elif fence:
                current.append(line)
        return blocks

    assert code_blocks(original) == code_blocks(destination.read_text()), "a replay command block changed during relocation"
    print(f"PASS cleanup preservation: {len(records)} byte-identical snapshots, prior scope, "
          f"{len(code_blocks(original))} unchanged replay command blocks")


if __name__ == "__main__":
    audit()
