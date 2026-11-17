"""Resolve explicit theorem references in TeX proofs across local input files.

This module reads sources and metadata only. It neither proves implications nor
changes the catalogue. Ambiguous references are reported for review.
"""
from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import re

CLAIM_ENVS = {"theorem", "lemma", "proposition", "corollary", "claim", "assertion",
              "thm", "lem", "prop", "cor"}
CLAIM_LABEL = re.compile(r"^(?:thm|lem|prop|cor|claim|assert)(?::|[-_.])", re.I)
NONCLAIM_LABEL = re.compile(r"^(?:eq|sec|subsec|fig|tab|alg|app)(?::|[-_.])", re.I)
INPUT = re.compile(r"(?<!\\)\\(?:input|include)\s*(?:\{([^{}]*)\}|([^\s{}%]+))")
REF = re.compile(r"\\(?:eqref|ref|cref|Cref|autoref)\s*\{([^}]+)\}")
EVENT = re.compile(
    r"\\begin\{(?P<begin>[^{}]+)\}(?:\[(?P<title>[^\]]*)\])?"
    r"|\\end\{(?P<end>[^{}]+)\}"
    r"|\\label\s*\{(?P<label>[^}]+)\}"
    r"|\\(?P<ref>eqref|ref|cref|Cref|autoref)\s*\{(?P<refs>[^}]+)\}"
)


def _blank(value: str) -> str:
    return "".join("\n" if char == "\n" else " " for char in value)


def _comments(text: str) -> str:
    """Keep offsets and line numbers while hiding comments and code examples."""
    text = re.sub(r"(?s)\\begin\{(verbatim\*?|lstlisting|minted|comment)\}.*?"
                  r"\\end\{\1\}", lambda match: _blank(match[0]), text)
    text = re.sub(r"\\verb\*?(?P<delimiter>[^\w\s]).*?(?P=delimiter)",
                  lambda match: _blank(match[0]), text)
    result = []
    for line in text.splitlines(keepends=True):
        cut = len(line)
        for index, char in enumerate(line):
            if char != "%":
                continue
            backslashes, previous = 0, index - 1
            while previous >= 0 and line[previous] == "\\":
                backslashes += 1
                previous -= 1
            if backslashes % 2 == 0:
                cut = index
                break
        result.append(line[:cut] + _blank(line[cut:]))
    return "".join(result)


def _relative(path: Path, root: Path) -> str:
    return Path(os.path.relpath(path, root)).as_posix()


@dataclass
class _Piece:
    start: int
    end: int
    path: Path
    offset: int


class _Document:
    def __init__(self, entrypoint: Path, root: Path, repository: Path, issues: list[dict]):
        self.entrypoint, self.root, self.repository = entrypoint, root, repository
        self.issues = issues
        self.parts: list[str] = []
        self.pieces: list[_Piece] = []
        self.sources: dict[Path, str] = {}
        self.length = 0
        self.files: set[Path] = set()
        self._expand(entrypoint, [])
        self.text = "".join(self.parts)
        self.starts = [piece.start for piece in self.pieces]

    def _issue(self, kind: str, path: Path, line: int, **detail) -> None:
        self.issues.append({"kind": kind, "entrypoint": _relative(self.entrypoint, self.root),
                            "source_path": _relative(path, self.root), "line": line, **detail})

    def _append(self, value: str, path: Path, offset: int) -> None:
        if value:
            self.parts.append(value)
            self.pieces.append(_Piece(self.length, self.length + len(value), path, offset))
            self.length += len(value)

    def _expand(self, path: Path, stack: list[Path]) -> None:
        if path in stack:
            self._issue("input_cycle", path, 1, chain=[_relative(p, self.root) for p in [*stack, path]])
            return
        self.files.add(path)
        raw = path.read_text(encoding="utf-8")
        self.sources[path] = raw
        text = _comments(raw)
        position = 0
        for match in INPUT.finditer(text):
            self._append(text[position:match.start()], path, position)
            value = (match[1] or match[2]).strip()
            line = text.count("\n", 0, match.start()) + 1
            if (not value or any(char in value for char in "\\#$~")
                    or PurePosixPath(value).is_absolute() or PureWindowsPath(value).drive):
                self._issue("unsupported_input", path, line, reference=value)
            else:
                included = Path(value)
                if not included.suffix:
                    included = included.with_suffix(".tex")
                candidates = {(directory / included).resolve()
                              for directory in (self.entrypoint.parent, path.parent)}
                candidates = {candidate for candidate in candidates
                              if candidate.is_relative_to(self.repository) and candidate.is_file()}
                if len(candidates) == 1:
                    self._expand(next(iter(candidates)), [*stack, path])
                elif not candidates:
                    self._issue("missing_input", path, line, reference=value)
                else:
                    self._issue("ambiguous_input", path, line, reference=value,
                                candidates=sorted(_relative(p, self.root) for p in candidates))
            position = match.end()
        self._append(text[position:], path, position)

    def location(self, position: int) -> tuple[Path, int]:
        index = max(0, bisect_right(self.starts, position) - 1)
        piece = self.pieces[index]
        offset = piece.offset + position - piece.start
        return piece.path, self.sources[piece.path].count("\n", 0, offset) + 1


def _declarations(document: _Document) -> tuple[dict[str, list[dict | None]], list[dict]]:
    """Attach a proof to an explicit title or the preceding theorem declaration."""
    labels: dict[str, list[dict | None]] = {}
    proofs, stack, pending = [], [], None
    for event in EVENT.finditer(document.text):
        if event["begin"]:
            env = event["begin"].removesuffix("*")
            item = {"env": env, "start": event.start(), "labels": []}
            if env in CLAIM_ENVS:
                item["claim"] = True
                pending = None
            elif env == "proof":
                item.update({"proof": True, "owner": pending, "refs": [],
                             "title_refs": [ref.strip() for match in REF.finditer(event["title"] or "")
                                            for ref in match[1].split(",")]})
                pending = None
            stack.append(item)
        elif event["end"]:
            env = event["end"].removesuffix("*")
            index = next((index for index in range(len(stack) - 1, -1, -1)
                          if stack[index]["env"] == env), None)
            if index is None:
                continue
            item = stack[index]
            del stack[index:]
            item["end"] = event.end()
            if item.get("claim"):
                pending = item
            elif item.get("proof"):
                proofs.append(item)
        elif event["label"]:
            name = event["label"].strip()
            claim = next((item for item in reversed(stack) if item.get("claim")), None)
            top_is_claim = bool(stack and stack[-1].get("claim"))
            attached = (claim if claim and not NONCLAIM_LABEL.match(name)
                        and (top_is_claim or CLAIM_LABEL.match(name)) else None)
            labels.setdefault(name, []).append(attached)
            if attached is not None:
                attached["labels"].append(name)
        elif event["refs"]:
            # A reference in a nested theorem statement belongs to that statement,
            # rather than becoming a dependency of an enclosing proof.
            innermost = next((item for item in reversed(stack)
                              if item.get("proof") or item.get("claim")), None)
            if innermost is not None and innermost.get("proof"):
                for name in event["refs"].split(","):
                    innermost["refs"].append((name.strip(), event.start(), event["ref"]))
    return labels, proofs


def _paths(data: dict, graph) -> list[Path]:
    values = set()
    for node in graph.nodes.values():
        values.update(source["path"] for source in node.get("sources", []))
        if node.get("canonical_source"):
            values.add(node["canonical_source"])
    for rule in graph.rules.values():
        if rule.get("proof"):
            values.add(rule["proof"]["path"])
    coverage = data.get("graph", {}).get("coverage", {})
    documents = coverage.get("documents", []) if isinstance(coverage, dict) else coverage
    for document in documents:
        if isinstance(document, dict) and document.get("path"):
            values.add(document["path"])
    paths = []
    for value in values:
        if PurePosixPath(value).is_absolute() or PureWindowsPath(value).drive:
            continue
        path = (graph.root / value).resolve()
        # Deposited snapshots and primary-source mirrors are evidence leaves,
        # not another current paper whose entire internal proof is inventoried.
        if "published" in path.parts or "sources" in path.parts:
            continue
        if path.is_relative_to(graph.repository) and path.suffix == ".tex" and path.is_file():
            paths.append(path)
    return sorted(set(paths))


def _node_for_label(graph, path: Path, label: str, source_index: dict) -> tuple[str | None, list[str]]:
    candidates = []
    pattern = re.compile(r"(?<![\w:.-])" + re.escape(label) + r"(?![\w:.-])")
    for node_id, locator in source_index.get(path, []):
        if not pattern.search(locator):
            continue
        score = 100 if (node_id.startswith(("MOM-", "FACT-", "MRC-", "ARITH-"))
                        or re.search(r"-(?:THM|LEM|PROP|COR)-", node_id)) else 0
        score += 20 if locator.strip() == label else 10 if locator.startswith(label) else 1
        score -= len(set(re.findall(r"(?:thm|lem|prop|cor):[\w.-]+", locator))) - 1
        candidates.append((score, node_id))
    if not candidates:
        return None, []
    best = max(score for score, _ in candidates)
    selected = sorted({node_id for score, node_id in candidates if score == best})
    return (selected[0] if len(selected) == 1 else None), selected


def tex_proof_link_report(data: dict, graph) -> dict:
    """Return explicit source-backed links plus unresolved/ambiguous input issues.

    Each documentclass entrypoint has a separate label namespace. A source that
    is not reachable from any registered entrypoint is inspected alone and
    reported as an orphan fragment. No catalogue, source, or digest is changed.
    """
    issues: list[dict] = []
    source_index = {}
    for node_id, node in graph.nodes.items():
        if node.get("kind") not in {"source", "gap", "goal"}:
            for source in node.get("sources", []):
                source_path = (graph.root / source["path"]).resolve()
                source_index.setdefault(source_path, []).append((node_id, source.get("locator", "")))
    resolution_cache = {}
    paths = _paths(data, graph)
    entrypoints = [path for path in paths
                   if re.search(r"\\documentclass(?:\[[^\]]*\])?\s*\{",
                                _comments(path.read_text(encoding="utf-8")))]
    documents, reached = [], set()
    for path in entrypoints:
        document = _Document(path, graph.root, graph.repository, issues)
        documents.append(document)
        reached.update(document.files)
    for path in paths:
        if path not in reached:
            document = _Document(path, graph.root, graph.repository, issues)
            document._issue("orphan_tex_fragment", path, 1,
                            detail="No registered documentclass entrypoint includes this fragment.")
            documents.append(document)
    links, seen = [], set()
    for document in documents:
        labels, proofs = _declarations(document)

        def resolve(label: str, position: int, role: str) -> tuple[dict, str] | None:
            definitions = labels.get(label, [])
            path, line = document.location(position)
            if not definitions:
                if not NONCLAIM_LABEL.match(label):
                    document._issue("unresolved_reference", path, line, reference=label, role=role)
                return None
            declarations = [declaration for declaration in definitions if declaration is not None]
            if not declarations:  # Equation, section, figure, etc.
                return None
            if len(definitions) != 1:
                document._issue("ambiguous_label", path, line, reference=label, role=role)
                return None
            declaration = declarations[0]
            declaration_path, _ = document.location(declaration["start"])
            key = declaration_path, label
            if key not in resolution_cache:
                resolution_cache[key] = _node_for_label(graph, declaration_path, label, source_index)
            node_id, candidates = resolution_cache[key]
            if node_id is None:
                document._issue("ambiguous_node" if candidates else "unmapped_theorem",
                                path, line, reference=label, role=role, candidates=candidates,
                                dependency_path=_relative(declaration_path, graph.root))
                return None
            return declaration, node_id

        for proof in proofs:
            title_refs = proof["title_refs"]
            if title_refs:
                owners = [resolve(label, proof["start"], "conclusion") for label in title_refs]
                owners = [owner for owner in owners if owner is not None]
                if len(owners) != 1:
                    path, line = document.location(proof["start"])
                    document._issue("ambiguous_proof_owner", path, line, references=title_refs)
                    continue
                declaration, conclusion = owners[0]
                conclusion_label = next(label for label in declaration["labels"] if label in title_refs)
            else:
                declaration = proof["owner"]
                if declaration is None or not declaration["labels"]:
                    # Unlabelled auxiliary proofs have no stable catalogue target.
                    continue
                conclusion_label = declaration["labels"][0]
                resolved = resolve(conclusion_label, proof["start"], "conclusion")
                if resolved is None:
                    continue
                declaration, conclusion = resolved
            proof_path, proof_line = document.location(proof["start"])
            conclusion_path, _ = document.location(declaration["start"])
            for label, position, command in proof["refs"]:
                resolved = resolve(label, position, "dependency")
                if resolved is None:
                    continue
                dependency_declaration, dependency = resolved
                if dependency == conclusion:
                    continue
                dependency_path, _ = document.location(dependency_declaration["start"])
                path, line = document.location(position)
                identity = (conclusion, dependency, path, line, label)
                if identity in seen:
                    continue
                seen.add(identity)
                links.append({
                    "entrypoint": _relative(document.entrypoint, graph.root),
                    "source_path": _relative(path, graph.root),
                    "locator": f"proof of {conclusion_label}; reference {label} at line {line}",
                    "conclusion": conclusion, "conclusion_path": _relative(conclusion_path, graph.root),
                    "dependency": dependency, "dependency_path": _relative(dependency_path, graph.root),
                    "reference": label,
                    "proof_locator": f"proof of {conclusion_label} beginning at "
                                     f"{_relative(proof_path, graph.root)}:{proof_line}",
                })
    return {"links": sorted(links, key=lambda item: (item["conclusion"], item["dependency"],
                                                    item["source_path"], item["locator"])),
            "issues": issues}


def tex_proof_links(data: dict, graph) -> list[dict]:
    """Return links only; use tex_proof_link_report for every review issue."""
    return tex_proof_link_report(data, graph)["links"]
