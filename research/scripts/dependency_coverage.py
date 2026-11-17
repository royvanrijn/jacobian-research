"""Inventory retained research sources and detect omissions from the one catalogue.

Coverage records map source declarations to graph statements. They do not prove
those statements or certify complete expansion of their mathematical imports.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

ROLES = {"proof", "manuscript", "formal-support", "imported-source",
         "scope", "evidence", "reference", "documentation"}
SUFFIXES = {".md", ".tex", ".lean", ".typ", ".py", ".json", ".pdf", ".txt", ".toml", ".html", ".gz", ".bib", ".bbl", ".zip", ".yml", ".yaml"}
GENERATED = {"README.md", "STATUS.md", "SCOPES.md", "DEPENDENCIES.md", "COVERAGE.md", "REFERENCES.md", "BOUNDS.md", "OPPORTUNITIES.md"}
GUIDES = {"DEPENDENCY_GUIDE.md", "REPRODUCE.md", "MATH_STATUS.json"}
TOPICS = {"moments", "factorial", "jacobian", "hessian", "dixmier",
          "arithmetic", "foundations", "meridional-rank"}


class CoverageError(ValueError):
    """A retained source or declaration has no current inventory mapping."""


def inventory_paths(root: Path, repository: Path) -> list[str]:
    """Ask Git for retained current mathematical and scope documents.

    Generated navigation, tooling documentation and immutable deposited archives
    are excluded explicitly. Current paper revisions remain in the inventory.
    """
    root, repository = Path(root).resolve(), Path(repository).resolve()
    names = subprocess.check_output(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard",
         "--", "research", "results", "unpublished"], cwd=repository).decode().split("\0")
    selected = set()
    for name in names:
        path = repository / name
        parts = Path(name).parts
        if not name or not path.is_file() or (path.suffix.lower() not in SUFFIXES
                                              and path.name not in {"lean-toolchain", "Makefile"}):
            continue
        if "published" in parts or ".lake" in parts or "__pycache__" in parts:
            continue
        if parts[0] == "research":
            relative = path.relative_to(root)
            if len(relative.parts) == 1 and relative.name in GENERATED | GUIDES:
                continue
            if len(relative.parts) == 2 and relative.parts[0] in TOPICS and relative.name == "README.md":
                continue
            if relative.parts[0] in {"scripts", "tests"}:
                continue
        selected.add(Path(os.path.relpath(path, root)).as_posix())
    return sorted(selected)


def _read(root: Path, repository: Path, value: str) -> Path:
    path = (root / value).resolve()
    if Path(value).is_absolute() or not path.is_relative_to(repository.resolve()):
        raise CoverageError("coverage path escapes repository: " + value)
    if not path.is_file():
        raise CoverageError("missing coverage source: " + value)
    return path


def declarations(path: Path) -> list[dict]:
    """Index explicit named statements, retaining their exact source locators."""
    if path.suffix.lower() not in {".md", ".tex", ".lean"}:
        return []
    text = path.read_text(encoding="utf-8")
    found = []
    if path.suffix.lower() == ".tex":
        # An escaped percent is text; ordinary comments are excluded.
        text = re.sub(r"(?m)(?<!\\)%.*$", "", text)
        pattern = r"\\begin\{(theorem|lemma|proposition|corollary|claim|conjecture|question|problem)\*?\}(.*?)(?=\\end\{\1\*?\})"
        for match in re.finditer(pattern, text, re.S):
            line = text.count("\n", 0, match.start()) + 1
            body = match.group(2)
            label = re.search(r"\\label\{([^}]+)\}", body)
            title = re.match(r"\s*\[([^]]+)\]", body)
            locator = label.group(1) if label else match.group(1) + "@" + str(line)
            found.append({"locator": locator, "kind": match.group(1), "line": line,
                          "title": title.group(1) if title else locator})
    elif path.suffix.lower() == ".lean":
        # Source modules are inventoried as formal support. An indexed theorem is
        # not thereby promoted to a whole-paper assurance.
        from dependency_sources import _lean_without_comments_and_strings
        cleaned = _lean_without_comments_and_strings(text)
        pattern = r"(?m)^\s*(?:(?:private|protected|noncomputable|unsafe)\s+)*(theorem|lemma|axiom|def|abbrev|structure|class)\s+([A-Za-z_][\w'.]*)"
        for match in re.finditer(pattern, cleaned):
            line = cleaned.count("\n", 0, match.start()) + 1
            found.append({"locator": match.group(2), "kind": match.group(1), "line": line,
                          "title": match.group(2)})
    elif path.suffix.lower() == ".md":
        text = re.sub(r"(?ms)^\s*```.*?^\s*```[^\n]*", "", text)
        pattern = r"(?mi)^\s*(?:>\s*)*(?:#{1,6}\s+|\*\*)(Theorem|Lemma|Proposition|Corollary|Claim|Conjecture|Open question|Problem|Obligation|Gap)\b([^\n]*)"
        for match in re.finditer(pattern, text):
            title = (match.group(1) + match.group(2)).strip().strip("*# ")
            line = text.count("\n", 0, match.start()) + 1
            found.append({"locator": title, "kind": match.group(1).lower(),
                          "line": line, "title": title})
    # Distinct same-name declarations need distinct locators within one file.
    counts = Counter(item["locator"] for item in found)
    for item in found:
        if counts[item["locator"]] > 1:
            item["locator"] += "@" + str(item["line"])
    return found


def coverage_state(data: dict, graph, paths: list[str] | None = None) -> dict:
    """Return omissions and freshness, without changing files or recorded states."""
    root, repository = graph.root, graph.repository
    wanted = set(inventory_paths(root, repository) if paths is None else paths)
    coverage = data.get("graph", {}).get("coverage", {})
    errors, rows = [], []
    if coverage.get("schema_version") != 1:
        errors.append("missing or unsupported graph.coverage schema")
    records = coverage.get("documents", [])
    if not isinstance(records, list):
        records = []
        errors.append("coverage.documents must be a list")
    indexed = {}
    for record in records:
        value = record.get("path")
        if not isinstance(value, str):
            errors.append("coverage record lacks path")
            continue
        if value in indexed:
            errors.append("duplicate coverage document: " + value)
        indexed[value] = record
    for value in sorted(wanted - indexed.keys()):
        errors.append("unmapped retained document: " + value)
    for value in sorted(indexed.keys() - wanted):
        errors.append("coverage document outside current inventory: " + value)
    known = set(graph.nodes)
    for value in sorted(wanted & indexed.keys()):
        record = indexed[value]
        try:
            path = _read(root, repository, value)
        except CoverageError as error:
            errors.append(str(error))
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != record.get("sha256"):
            errors.append("changed coverage source: " + value)
        role = record.get("role")
        if role not in ROLES:
            errors.append("invalid coverage role: " + value)
        identifiers = []
        for key in ("claims", "dependencies", "gaps"):
            entries = record.get(key, [])
            if not isinstance(entries, list) or not all(isinstance(item, str) for item in entries):
                errors.append("invalid coverage " + key + ": " + value)
                entries = []
            identifiers.extend(entries)
        if not identifiers and role not in {"reference", "documentation"}:
            errors.append("source has no graph mapping: " + value)
        if not str(record.get("explanation", record.get("reason", ""))).strip():
            errors.append("coverage classification lacks explanation: " + value)
        for node_id in identifiers:
            if node_id not in known:
                errors.append("unknown coverage node " + node_id + ": " + value)
        statements = record.get("statements", [])
        if not isinstance(statements, list):
            statements = []
            errors.append("invalid statement mapping: " + value)
        statement_map = {}
        for statement in statements:
            locator = statement.get("locator")
            if locator in statement_map:
                errors.append("duplicate statement locator: " + value + " :: " + str(locator))
            statement_map[locator] = statement
            targets = statement.get("nodes", [])
            if not isinstance(targets, list) or not targets or any(target not in known for target in targets):
                errors.append("unmapped or unknown statement targets: " + value + " :: " + str(locator))
        actual_statements = declarations(path) if role not in {"reference", "documentation"} else []
        for declaration in actual_statements:
            if declaration["locator"] not in statement_map:
                errors.append("unmapped source declaration: " + value + " :: " + declaration["locator"])
        rows.append({**record, "current_sha256": actual,
                     "declaration_count": len(actual_statements),
                     "fresh": actual == record.get("sha256")})
    return {"documents": rows, "expected_documents": len(wanted), "errors": errors,
            "complete": not errors, "declarations": sum(row["declaration_count"] for row in rows)}


def check_coverage(data: dict, graph, paths: list[str] | None = None) -> dict:
    state = coverage_state(data, graph, paths)
    if state["errors"]:
        shown = state["errors"][:25]
        suffix = "\n... " + str(len(state["errors"]) - 25) + " further issues" if len(state["errors"]) > 25 else ""
        raise CoverageError("research coverage is incomplete or stale:\n" + "\n".join(shown) + suffix)
    return state


def render_coverage(data: dict, graph) -> str:
    state = coverage_state(data, graph)
    lines = ["# Research source coverage", "",
             "<!-- Generated from MATH_STATUS.json; run make render-navigation. -->", "",
             "This index maps the retained current research and paper sources to graph statements, "
             "bridges, support and gaps. It detects omitted documents and explicit named declarations. "
             "Source mapping completeness is separate from mathematical dependency completeness: "
             "a mapped written theorem can still have unverified imports or an open proof obligation.", "",
             "Generated navigation, tooling guides and immutable deposited archives are excluded from "
             "this current-source inventory. Archived publication integrity is checked separately.", "",
             f"**Inventory:** {len(state['documents'])}/{state['expected_documents']} documents; "
             f"{state['declarations']} explicit source declarations indexed. "
             + ("Current mappings pass." if state["complete"] else str(len(state["errors"])) + " coverage issues."), "",
             "[Dependency graph](DEPENDENCIES.md) · [Exact scopes](SCOPES.md) · "
             "[Maintenance guide](DEPENDENCY_GUIDE.md)", ""]
    if state["errors"]:
        lines += ["## Coverage issues", ""]
        lines += ["- " + error for error in state["errors"]]
        lines += [""]
    lines += ["| Source | Role | Statements | Mapped claims / gaps | Assessment |",
              "|---|---|---|---|---|"]
    cell = lambda value: str(value).replace("|", "\\|").replace("\n", " ")
    for row in state["documents"]:
        nodes = list(dict.fromkeys([*row.get("claims", []), *row.get("gaps", [])]))
        lines.append("| [" + cell(row["path"]) + "](" + row["path"] + ") | " + row.get("role", "")
                     + " | " + str(row["declaration_count"]) + " | " + cell(", ".join(nodes))
                     + " | " + cell(row.get("explanation", row.get("reason", ""))) + " |")
    return "\n".join(lines) + "\n"
