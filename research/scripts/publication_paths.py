"""Paths for current paper editions and separately retained public deposits."""

from pathlib import Path


def result_paper_directories(results: Path) -> list[Path]:
    """Find current papers directly under results or in its published group."""
    group = results / "published"
    papers = [path for path in results.glob("*") if path.is_dir() and path != group]
    papers.extend(path for path in group.glob("*") if path.is_dir())
    return papers


def is_deposited_path(parts: tuple[str, ...]) -> bool:
    """Recognize archive path components after the collection name."""
    if parts[:2] == ("results", "published"):
        parts = parts[2:]
    return "published" in parts
