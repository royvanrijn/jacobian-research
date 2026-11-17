"""Discover retained local proof-source dependencies without editing files or pins."""

from __future__ import annotations

from collections import deque
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import re


_TEX_SUFFIXES = {".tex", ".sty", ".cls", ".ltx"}


def _within(path: Path, repository: Path) -> bool:
    return path.is_relative_to(repository)


def _source_path(value: str, root: Path, repository: Path) -> Path:
    if not isinstance(value, str) or not value or PurePosixPath(value).is_absolute() or PureWindowsPath(value).drive:
        raise ValueError(f"invalid repository-relative source path: {value}")
    path = (root / value).resolve()
    if not _within(path, repository):
        raise ValueError(f"source escapes repository: {value}")
    return path


def _tex_without_comments(text: str) -> str:
    lines = []
    for line in text.splitlines(keepends=True):
        end = len(line)
        for index, char in enumerate(line):
            if char != "%":
                continue
            backslashes, previous = 0, index - 1
            while previous >= 0 and line[previous] == "\\":
                backslashes += 1
                previous -= 1
            if backslashes % 2 == 0:
                end = index
                break
        lines.append(line[:end] + ("\n" if end < len(line) else ""))
    return "".join(lines)


def _lean_without_comments_and_strings(text: str) -> str:
    """Preserve lines while ignoring nested block comments and quoted strings."""
    output, index, depth, quoted = [], 0, 0, False
    while index < len(text):
        pair, char = text[index:index + 2], text[index]
        if depth:
            if pair == "/-":
                depth += 1
                output.extend("  ")
                index += 2
            elif pair == "-/":
                depth -= 1
                output.extend("  ")
                index += 2
            else:
                output.append("\n" if char == "\n" else " ")
                index += 1
        elif quoted:
            if char == "\\" and index + 1 < len(text):
                output.extend("  ")
                index += 2
            else:
                output.append("\n" if char == "\n" else " ")
                if char == '"':
                    quoted = False
                index += 1
        elif pair == "/-":
            depth = 1
            output.extend("  ")
            index += 2
        elif pair == "--":
            end = text.find("\n", index)
            end = len(text) if end < 0 else end
            output.extend(" " * (end - index))
            index = end
        elif char == '"':
            quoted = True
            output.append(" ")
            index += 1
        else:
            output.append(char)
            index += 1
    return "".join(output)


def _lean_project(path: Path, repository: Path) -> Path | None:
    for directory in path.parents:
        if not _within(directory, repository):
            break
        if any((directory / name).is_file() for name in ("lakefile.toml", "lakefile.lean", "lean-toolchain")):
            return directory
    return None


def source_closure(paths: list[str], root: Path, repository: Path) -> list[dict]:
    """Return seeds and reachable retained TeX/Lean sources, with source locators.

    TeX paths must be literal input/include names. Both the entrypoint directory
    and the including file's directory are considered, conservatively retaining
    both when ambiguous. Lean imports resolve only within the nearest local Lake
    project (including literal srcDir settings), never downloaded packages.
    Missing files remain the graph freshness check's responsibility.
    """
    root, repository = Path(root).resolve(), Path(repository).resolve()
    if not _within(root, repository):
        raise ValueError("research root escapes repository")
    found: dict[Path, str] = {}
    pending = deque()
    visited = set()

    def enqueue(path: Path, locator: str, tex_root: Path | None = None) -> None:
        path = path.resolve()
        if not _within(path, repository) or not path.is_file():
            return
        found.setdefault(path, locator)
        context = tex_root if path.suffix.lower() in _TEX_SUFFIXES else None
        key = (path, context)
        if key not in visited:
            visited.add(key)
            pending.append(key)

    for value in paths:
        path = _source_path(value, root, repository)
        enqueue(path, "Explicit retained source", path.parent if path.suffix.lower() in _TEX_SUFFIXES else None)

    while pending:
        path, tex_root = pending.popleft()
        if path.suffix.lower() in _TEX_SUFFIXES:
            text = _tex_without_comments(path.read_text(encoding="utf-8"))
            pattern = r"(?<!\\)\\(?:input|include)\s*(?:\{([^{}]*)\}|([^\s{}%]+))"
            for match in re.finditer(pattern, text):
                value = (match.group(1) or match.group(2)).strip()
                if not value or any(char in value for char in "\\#$~") or PurePosixPath(value).is_absolute() or PureWindowsPath(value).drive:
                    continue
                include = Path(value)
                if not include.suffix:
                    include = include.with_suffix(".tex")
                line = text.count("\n", 0, match.start()) + 1
                for directory in dict.fromkeys((tex_root or path.parent, path.parent)):
                    enqueue(directory / include, f"Literal TeX include from {os.path.relpath(path, root)}:{line}", tex_root)
        elif path.suffix.lower() == ".lean":
            project = _lean_project(path, repository)
            if project is None:
                continue
            source_roots = [project]
            for name in ("lakefile.toml", "lakefile.lean", "lean-toolchain", "lake-manifest.json"):
                config = project / name
                enqueue(config, "Retained local Lean project configuration")
                if config.is_file() and _within(config.resolve(), repository) and name.startswith("lakefile"):
                    for value in re.findall(r'\bsrcDir\s*(?::=|=)\s*"([^"]+)"', config.read_text(encoding="utf-8")):
                        candidate = (project / value).resolve()
                        if _within(candidate, repository):
                            source_roots.append(candidate)
            text = _lean_without_comments_and_strings(path.read_text(encoding="utf-8"))
            for match in re.finditer(r"(?m)^\s*(?:(?:public|private|meta)\s+)*import[ \t]+([^\n]+)", text):
                for module in match.group(1).split():
                    if module == "all" or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_'.]*", module):
                        continue
                    relative = Path(*module.split(".")).with_suffix(".lean")
                    line = text.count("\n", 0, match.start()) + 1
                    for directory in source_roots:
                        enqueue(directory / relative, f"Local Lean import {module} from {os.path.relpath(path, root)}:{line}")
    return [{"path": Path(os.path.relpath(path, root)).as_posix(), "locator": found[path]}
            for path in sorted(found)]


def missing_source_closure(graph) -> list[dict]:
    """Report unpinned local includes/imports of established source leaves.

    Contextual sources on goals/obligations, non-established declarations and
    derived nodes are excluded. No hashes, evidence or catalogue bytes change.
    """
    missing = []
    for owner, node in graph.nodes.items():
        if graph.incoming.get(owner):
            continue
        if node.get("kind") in {"goal", "gap", "mathematical_goal"} or "obligation" in node.get("kind", ""):
            continue
        if node.get("status", node.get("recorded_state")) not in {"proved", "falsified", "refuted"}:
            continue
        paths = [source["path"] for source in node.get("sources", [])]
        pinned = {_source_path(value, graph.root, graph.repository) for value in paths}
        for item in source_closure(paths, graph.root, graph.repository):
            if _source_path(item["path"], graph.root, graph.repository) not in pinned:
                missing.append({"owner": owner, "path": item["path"]})
    return sorted(missing, key=lambda item: (item["owner"], item["path"]))

