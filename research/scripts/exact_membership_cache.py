"""Content-addressed cache for expensive exact Singular membership checks."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA = 1


def program_sha256(program: str) -> str:
    return hashlib.sha256(program.encode("utf-8")).hexdigest()


def cached_basis_size(
    path: Path, program: str, metadata: dict[str, Any]
) -> int | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None
    expected = {
        "schema": SCHEMA,
        "program_sha256": program_sha256(program),
        "membership": True,
        **metadata,
    }
    if any(payload.get(key) != value for key, value in expected.items()):
        return None
    basis_size = payload.get("basis_size")
    return basis_size if isinstance(basis_size, int) and basis_size > 0 else None


def write_exact_membership_cache(
    path: Path,
    program: str,
    metadata: dict[str, Any],
    basis_size: int,
    singular_version: str,
) -> None:
    payload = {
        "schema": SCHEMA,
        "program_sha256": program_sha256(program),
        "membership": True,
        **metadata,
        "basis_size": basis_size,
        "singular_version": singular_version.strip(),
        "verified_at": datetime.now(timezone.utc).isoformat(),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)
