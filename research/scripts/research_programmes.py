"""Programme scope is independent of the mathematical state of a claim."""

from pathlib import Path
from functools import lru_cache
import json

ACTIVE_AREAS = {'elliptic-curves', 'elkies-k3'}
ARCHIVE_PREFIX = 'archive/non-elliptic/'


def is_active(record: dict) -> bool:
    return record.get('programme_status', 'active') == 'active'


def validate_programme(record: dict) -> None:
    assert record.get('programme_status', 'active') in {'active', 'archived'}, (
        f"{record.get('id')}: invalid programme status"
    )


def programme_path_allowed(path: str, record: dict) -> bool:
    parts = Path(path).parts
    return not parts or parts[0] != 'archive' or (
        not is_active(record) and path.startswith(ARCHIVE_PREFIX)
    )


@lru_cache(maxsize=1)
def archived_roots() -> frozenset[str]:
    manifest = Path(__file__).resolve().parents[1] / (ARCHIVE_PREFIX + 'MANIFEST.json')
    return frozenset(json.loads(manifest.read_text())['root_moves']) if manifest.is_file() else frozenset()


def relocate_reference(path: str) -> str:
    """Map a former research-relative path to its byte-preserved archived source."""
    parts = Path(path)
    roots = archived_roots()
    if path in roots or any(str(parent) in roots for parent in parts.parents):
        return ARCHIVE_PREFIX + path
    return path
