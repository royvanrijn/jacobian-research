"""Enumerate maintained documentation without walking local computation caches."""

from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parent


def document_paths(repository: Path = REPOSITORY) -> list[Path]:
    """Include tracked and new nonignored Markdown, including retained archives.

    Tracked files remain included even under an ignore rule (pinned evidence).
    Ignored, untracked run trees and dependency environments are not source docs.
    Git handles pruning before traversal and NUL delimiters preserve path names.
    """
    result = subprocess.run(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard",
         "--", "*.md"],
        cwd=repository, check=True, capture_output=True,
    )
    return sorted({repository / name.decode() for name in result.stdout.split(b"\0")
                   if name and (repository / name.decode()).is_file()})
