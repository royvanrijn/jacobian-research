"""Resolve retained bulk evidence declared in the repository manifest."""

from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import re


@lru_cache
def artifact_manifest(repository: Path) -> tuple[Path | None, dict[str, str]]:
    path = repository / "external-artifacts.json"
    if not path.is_file():
        return None, {}
    data = json.loads(path.read_text(encoding="utf-8"))
    root = Path(data["external_root"])
    if data["schema_version"] != 1 or not root.is_absolute():
        raise ValueError("invalid external artifact manifest")
    files = {}
    for entry in data["files"]:
        name, digest = entry["path"], entry["sha256"]
        relative = PurePosixPath(name)
        if (
            relative.is_absolute() or PureWindowsPath(name).drive or ".." in relative.parts
            or not relative.parts or relative.as_posix() != name or name in files
            or not re.fullmatch(r"[0-9a-f]{64}", digest)
        ):
            raise ValueError(f"invalid external artifact entry: {name}")
        files[name] = digest
    return root.resolve(), files


def repository_path(path: Path, repository: Path) -> Path:
    local = Path(os.path.abspath(path))
    if not local.is_relative_to(repository):
        raise ValueError(f"source escapes repository: {path}")
    resolved = path.resolve()
    if resolved.is_relative_to(repository):
        return resolved
    root, files = artifact_manifest(repository)
    name = local.relative_to(repository).as_posix()
    if root is not None and name in files and resolved == root / name and local.is_symlink():
        return local
    raise ValueError(f"source escapes repository: {path}")


def sha256_file(path: Path) -> str:
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()
