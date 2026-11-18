#!/usr/bin/env python3
"""Current research scopes, proof dependencies and publication integrity; no proof replay."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

from zipfile import ZipFile

from dependency_graph import DependencyGraph, GraphValidationError
from dependency_sources import missing_source_closure
from dependency_coverage import CoverageError, check_coverage, coverage_state, render_coverage, declarations
from dependency_bounds import bounds_records, bound_comparisons
from dependency_references import check_references, reference_state, render_references
from dependency_tex_links import tex_proof_link_report

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parent
STATUS = ROOT / "MATH_STATUS.json"
STATES = {"proved", "partial", "falsified"}
ASSURANCES = ("independent_replay", "formal_verification", "external_review")


def load() -> dict:
    return json.loads(STATUS.read_text())


def reference(value: str, base: Path | None = None, *, directory: bool = False) -> Path:
    assert isinstance(value, str) and value, "empty source reference"
    path = Path(value)
    assert not path.is_absolute(), f"absolute source reference: {value}"
    resolved = ((ROOT if base is None else base) / path).resolve()
    assert resolved.is_relative_to(REPOSITORY.resolve()), f"reference escapes repository: {value}"
    assert resolved.is_dir() if directory else resolved.is_file(), f"missing reference: {value}"
    return resolved


def validate(data: dict) -> tuple[int, int]:
    assert data["schema_version"] in {2, 3}
    assert data["authority"] == "MATH_STATUS.json"
    topics = data["topics"]
    assert isinstance(topics, list) and topics, "no current research topics"
    topic_ids, result_ids = set(), set()
    for topic in topics:
        assert topic["id"] not in topic_ids, "duplicate topic"
        topic_ids.add(topic["id"])
        for key in ("id", "title", "summary", "status"):
            assert isinstance(topic[key], str) and topic[key].strip(), f"empty topic {key}"
        assert re.fullmatch(r"[a-z0-9-]+", topic["id"]), "invalid topic identifier"
        assert topic["readme"] == topic["id"] + "/README.md", "topic guide must be its own README"
        reference(topic["readme"])
        if "verification_note" in topic:
            assert isinstance(topic["verification_note"], str) and topic["verification_note"].strip()
        if "replay_record" in topic:
            reference(topic["replay_record"])
        assert isinstance(topic["results"], list) and isinstance(topic["open_questions"], list)
        for result in topic["results"]:
            assert result["id"] not in result_ids, f"duplicate result: {result['id']}"
            result_ids.add(result["id"])
            assert result["state"] in STATES, f"invalid state: {result['id']}"
            for key in ("title", "scope", "proof_type"):
                assert isinstance(result[key], str) and result[key].strip(), f"empty {key}: {result['id']}"
            reference(result["canonical_source"])
            for path in result["evidence"]:
                reference(path)
            for flag in ASSURANCES:
                assert isinstance(result[flag], bool), f"invalid assurance: {result['id']}/{flag}"
            checker = result["checker"]
            if checker is None:
                assert result["artifact_hash"] is None, f"hash without checker: {result['id']}"
            else:
                checker_path = reference(checker)
                expected = "sha256:" + hashlib.sha256(checker_path.read_bytes()).hexdigest()
                assert result["artifact_hash"] == expected, f"changed checker source: {result['id']}"
                if checker_path.suffix == ".py":
                    ast.parse(checker_path.read_text(), filename=str(checker_path))
            for support in result.get("additional_checkers", []):
                support_path = reference(support["path"])
                actual = hashlib.sha256(support_path.read_bytes()).hexdigest()
                assert actual == support["sha256"], f"changed supporting checker: {result['id']}"
                if support_path.suffix == ".py":
                    ast.parse(support_path.read_text(), filename=str(support_path))
        for question in topic["open_questions"]:
            assert isinstance(question["question"], str) and question["question"].strip()
            reference(question["source"])
    graph = dependency_graph(data)
    assert graph is not None, "schema 2+ requires an explicit dependency graph"
    if data["schema_version"] >= 3:
        assert data.get("graph", {}).get("coverage", {}).get("schema_version") == 1, (
            "schema 3 requires the current research source inventory")
    generated = {(ROOT / name).resolve() for name in
                 ["STATUS.md", "SCOPES.md", "DEPENDENCIES.md", "COVERAGE.md", "REFERENCES.md", "BOUNDS.md", "README.md",
                  *(topic["readme"] for topic in topics)]}
    for node in graph.nodes.values():
        if node["supplemental"]:
            assert node["kind"] in {"source", "claim", "gap", "goal"}, (
                f"unknown graph node kind: {node['id']}")
        for source in node.get("sources", []):
            assert (ROOT / source["path"]).resolve() not in generated, (
                f"generated navigation cannot establish graph evidence: {node['id']}")
    for rule in graph.rules.values():
        if rule.get("proof"):
            assert (ROOT / rule["proof"]["path"]).resolve() not in generated, (
                f"generated navigation cannot prove an implication: {rule['id']}")
    for topic in topics:
        assert isinstance(topic.get("scope_notes"), list) and topic["scope_notes"], (
            f"missing scope audit: {topic['id']}")
        assert isinstance(topic.get("gap_ids"), list), f"missing gap mapping: {topic['id']}"
        for node_id in topic["gap_ids"]:
            assert node_id in graph.nodes, f"unknown topic gap: {topic['id']}/{node_id}"
        for question in topic["open_questions"]:
            assert question.get("kind") in {"mathematical", "formalization", "verification", "workflow"}, (
                f"unclassified open question: {topic['id']}")
            assert question.get("node") in topic["gap_ids"], f"unmapped open question: {topic['id']}"
        for result in topic["results"]:
            assert result["id"] in graph.incoming, f"unmapped principal result: {result['id']}"
            assert result.get("dependency_coverage") in {"complete", "partial"}, (
                f"unspecified dependency coverage: {result['id']}")
    bounds_records(graph)
    return len(topics), len(result_ids)


def dependency_graph(data: dict) -> DependencyGraph | None:
    """Fixtures without a graph can still render; schema 2 validation requires it."""
    return DependencyGraph(data, ROOT, REPOSITORY) if data.get("graph") else None


def topic_progress(topic: dict, evaluation: dict) -> list[str]:
    if not evaluation or not topic.get("gap_ids"):
        return []
    resolved = [node_id for node_id in topic["gap_ids"]
                if evaluation[node_id]["state"] in {"proved", "falsified"}]
    unresolved = len(topic["gap_ids"]) - len(resolved)
    lines = ["**Computed obligations:** " + str(unresolved) + " unresolved; "
             + str(len(resolved)) + " resolved within their exact scope.", ""]
    closed_questions = [question["node"] for question in topic["open_questions"]
                        if question.get("node") in resolved]
    if closed_questions:
        lines += ["**Assessment update needed:** " + ", ".join(closed_questions)
                  + " now has a resolved graph state. Read the computed states below "
                    "when the authored assessment still uses open-problem wording.", ""]
    return lines


def tracked_question(question: dict, evaluation: dict) -> str:
    node_id = question.get("node")
    state = evaluation.get(node_id, {}).get("state", "unmapped")
    resolved = state in {"proved", "falsified"}
    label = state + (" (discharged)" if state == "proved" else "")
    return ("**" + question.get("kind", "mathematical") + " / " + label + ":** "
            + ("Authored question: " if resolved else "") + "["
            + question["question"] + "](" + question["source"] + ")"
            + (" — " + node_id if node_id else ""))


def render(data: dict) -> str:
    graph = dependency_graph(data)
    evaluation = graph.evaluate() if graph else {}
    lines = ["# Current mathematical status", "",
             "<!-- Generated by scripts/research.py render from MATH_STATUS.json. -->", "",
             "This is a map of the principal results and current open questions. "
             "A result's state applies to its stated scope. Publication, internal checking, "
             "independent replay, formal verification and external review remain separate.", "",
             "[Published papers](../README.md) · [Unpublished manuscripts](../unpublished/README.md) · "
             "[Research topics](README.md) · [Scope and gaps](SCOPES.md) · "
             "[Implications](DEPENDENCIES.md) · [Source coverage](COVERAGE.md) · [Bounds](BOUNDS.md) · [Cited sources](REFERENCES.md)", "",
             "The authored state describes the exact statement below. Graph status is computed "
             "from declared proof routes and current source pins; it is an internal evidence "
             "classification, not an independent proof audit. A proved scoped result can coexist "
             "with an open broader conjecture. Dependency coverage can remain partial.", ""]
    for topic in data["topics"]:
        lines += ["## " + topic["title"], "", topic["summary"], "",
                  "**Authored assessment:** " + topic["status"], "",
                  "[Topic guide](" + topic["readme"] + ")", ""]
        lines += topic_progress(topic, evaluation)
        if topic["results"]:
            lines += ["| Principal result | Authored state | Graph status | Scope |", "|---|---|---|---|"]
            for result in topic["results"]:
                cell = lambda text: text.replace("|", "\\|").replace("\n", " ")
                lines.append("| [" + cell(result["title"]) + "](" + result["canonical_source"] + ") | "
                             + result["state"] + " | " + evaluation.get(result["id"], {}).get("state", "unmapped")
                             + " | " + cell(result["scope"]) + " |")
            lines.append("")
        if topic["open_questions"]:
            lines += ["Tracked obligations (computed states):", ""]
            for question in topic["open_questions"]:
                lines.append("- " + tracked_question(question, evaluation))
            lines.append("")
    lines += ["The [structured record](MATH_STATUS.json) names each retained source, checker, "
              "evidence file and assurance flag. This navigation check does not rerun mathematics.", ""]
    return "\n".join(lines)


def markdown_cell(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


def topic_link(topic: dict, source: str) -> str:
    return Path(os.path.relpath(ROOT / source, (ROOT / topic["readme"]).parent)).as_posix()


def render_topic(topic: dict, evaluation: dict | None = None) -> str:
    evaluation = evaluation or {}
    lines = ["# " + topic["title"], "",
             "<!-- Generated from research/MATH_STATUS.json; run make render-navigation. -->", "",
             topic["summary"], "", "**Authored assessment:** " + topic["status"], "",
             "| Current result | Evidence / publication | Graph status | Main source |", "|---|---|---|---|"]
    for result in topic["results"]:
        source = result["canonical_source"]
        if source.startswith("../results/"):
            state = "Published preprints"
            source = "../results/README.md"
        elif source.startswith("../unpublished/"):
            state = "Unpublished; independent review pending"
            readme = str(Path(source).parent / "README.md")
            if (ROOT / readme).is_file():
                source = readme
        else:
            state = {"proved": "Internal result", "partial": "Partial result",
                     "falsified": "Counterexample"}[result["state"]]
        lines.append("| " + markdown_cell(result["title"]) + " | " + state
                     + " | " + evaluation.get(result["id"], {}).get("state", "unmapped")
                     + " | [Source](" + topic_link(topic, source) + ") |")
    lines.append("")
    lines += topic_progress(topic, evaluation)
    if topic.get("scope_notes"):
        lines += ["", "Scope limits:", ""]
        lines += ["- " + note for note in topic["scope_notes"]]
    if topic.get("verification_note"):
        note = topic["verification_note"]
        if topic.get("replay_record"):
            note += " [Replay record](" + topic_link(topic, topic["replay_record"]) + ")."
        lines += ["", "**Verification:** " + note]
    if topic["open_questions"]:
        lines += ["", "Tracked obligations (computed states):", ""]
        for question in topic["open_questions"]:
            relative_question = {**question, "source": topic_link(topic, question["source"])}
            lines.append("- " + tracked_question(relative_question, evaluation))
    lines += ["", "[Exact scopes and verification limits](" + topic_link(topic, "STATUS.md") + "). "
              "Publication, formal verification and independent review are separate statuses. "
              "[Dependency routes](" + topic_link(topic, "DEPENDENCIES.md") + ") and "
              "[scope audit](" + topic_link(topic, "SCOPES.md") + ") expose remaining gaps.", ""]
    return "\n".join(lines)


def render_research_index(data: dict) -> str:
    graph = dependency_graph(data)
    evaluation = graph.evaluate() if graph else {}
    lines = ["# Research", "",
             "<!-- Generated from MATH_STATUS.json; run make render-navigation. -->", "",
             "Current AI-assisted mathematics, organized by topic. "
             "[Published papers](../results/README.md) and "
             "[unpublished manuscripts](../unpublished/README.md) have separate catalogues.", "",
             "| Topic | Scope | Authored assessment | Computed obligations |", "|---|---|---|---|"]
    for topic in data["topics"]:
        gaps = [evaluation[node_id] for node_id in topic.get("gap_ids", []) if node_id in evaluation]
        resolved = sum(record["state"] in {"proved", "falsified"} for record in gaps)
        progress = f"{len(gaps) - resolved} unresolved; {resolved} resolved"
        if any(evaluation.get(question.get("node"), {}).get("state") in {"proved", "falsified"}
               for question in topic["open_questions"]):
            progress += "; assessment update needed"
        lines.append("| [" + markdown_cell(topic["title"]) + "](" + topic["readme"] + ") | "
                     + markdown_cell(topic["summary"]) + " | " + markdown_cell(topic["status"])
                     + " | " + progress + " |")
    lines += ["", "The [current result catalogue](STATUS.md) is generated from "
              "[MATH_STATUS.json](MATH_STATUS.json). It records main statements, proof sources "
              "and verification limits. Internal proof status does not imply publication.",
              "", "[Scope and missing gaps](SCOPES.md) · [Dependency DAG](DEPENDENCIES.md) · "
              "[How to maintain implications](DEPENDENCY_GUIDE.md) · [Source coverage](COVERAGE.md) · [Bounds](BOUNDS.md) · [Cited sources](REFERENCES.md)",
              "", "[Verification commands](REPRODUCE.md)", ""]
    return "\n".join(lines)


def navigation_outputs(data: dict) -> dict[Path, str]:
    outputs = {ROOT / "STATUS.md": render(data), ROOT / "README.md": render_research_index(data)}
    graph = dependency_graph(data)
    evaluation = graph.evaluate() if graph else {}
    for topic in data["topics"]:
        assert re.fullmatch(r"[a-z0-9-]+", topic["id"]), "invalid topic identifier"
        assert topic["readme"] == topic["id"] + "/README.md", "topic guide must be its own README"
        outputs[ROOT / topic["readme"]] = render_topic(topic, evaluation)
    if graph:
        from dependency_views import render_dependency_report, render_scope_report
        outputs[ROOT / "DEPENDENCIES.md"] = render_dependency_report(data, graph)
        outputs[ROOT / "SCOPES.md"] = render_scope_report(data, graph)
        if data.get("graph", {}).get("coverage"):
            outputs[ROOT / "COVERAGE.md"] = render_coverage(data, graph)
        outputs[ROOT / "BOUNDS.md"] = render_bounds(data, graph)
        if data.get("graph", {}).get("external_references"):
            outputs[ROOT / "REFERENCES.md"] = render_references(data, graph)
    return outputs


def check_navigation(data: dict) -> None:
    for path, expected in navigation_outputs(data).items():
        assert path.is_file() and path.read_text() == expected, (
            f"{path.relative_to(ROOT)} is stale; run make render-navigation")


def check_extracted_sources(meta: dict, directory: Path, sealed_names: set[str]) -> None:
    """Bind extracted publication sources to the exact deposited archive."""
    if meta["source"] in sealed_names:
        return
    archive_name = meta.get("source_bundle")
    assert archive_name in sealed_names, f"published source lacks a sealed archive: {directory.name}"
    members = meta.get("source_bundle_members")
    assert isinstance(members, dict) and meta["source"] in members, (
        f"published source is not bound to an archive member: {directory.name}")
    with ZipFile(reference(archive_name, directory)) as archive:
        for local, member in members.items():
            assert member in archive.namelist(), f"missing deposited source member: {member}"
            assert reference(local, directory).read_bytes() == archive.read(member), (
                f"extracted publication source changed: {directory.name}/{local}")


def check_deposited_edition(path: Path) -> None:
    """Validate the immutable public edition separately from a local revision."""
    meta = json.loads(path.read_text())
    assert all(isinstance(meta[key], str) and meta[key].strip()
               for key in ("title", "summary", "status", "source", "pdf")), str(path)
    assert meta["author"]["name"]
    assert re.fullmatch(r"\d{4}-\d{4}-\d{4}-[\dX]{4}", meta["author"]["orcid"])
    assert meta["status"] == "published-preprint", f"archive is not a deposited edition: {path}"
    assert re.fullmatch(r"10\.5281/zenodo\.\d+", meta["doi"])
    assert meta["zenodo"].startswith("https://zenodo.org/records/")
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", meta["publication_date"])
    reference("README.md", path.parent)
    reference(meta["source"], path.parent)
    assert meta["verification"]["scope"] and meta["verification"]["limits"]
    for evidence in meta["verification"]["evidence"]:
        reference(evidence, path.parent)
    for field in ("unpublished_revision", "current_edition"):
        if meta.get(field):
            reference(meta[field], path.parent, directory=True)
    seal = json.loads(reference("deposited-files.json", path.parent).read_text())
    assert seal["record"] == meta["zenodo"], f"deposit record mismatch: {path.parent}"
    sealed_names = set()
    for item in seal["files"]:
        contents = reference(item["file"], path.parent).read_bytes()
        assert hashlib.sha256(contents).hexdigest() == item["sha256"], f"deposit SHA mismatch: {item['file']}"
        assert hashlib.md5(contents).hexdigest() == item["md5"], f"deposit MD5 mismatch: {item['file']}"
        sealed_names.add(item["file"])
    assert meta["pdf"] in sealed_names, f"published PDF is not the sealed deposited file: {path.parent}"
    check_extracted_sources(meta, path.parent, sealed_names)


def check_portfolio() -> tuple[int, int]:
    counts = []
    for collection in ("results", "unpublished"):
        directory = REPOSITORY / collection
        assert all((folder / "metadata.json").is_file() for folder in directory.iterdir() if folder.is_dir()), (
            f"uncatalogued paper directory in {collection}")
        folders = sorted(directory.glob("*/metadata.json"))
        assert folders, f"no {collection} metadata"
        counts.append(len(folders))
        for path in folders:
            meta = json.loads(path.read_text())
            assert all(isinstance(meta[key], str) and meta[key].strip()
                       for key in ("title", "summary", "status", "source", "pdf")), str(path)
            assert meta["author"]["name"]
            assert re.fullmatch(r"\d{4}-\d{4}-\d{4}-[\dX]{4}", meta["author"]["orcid"])
            reference("README.md", path.parent)
            reference(meta["source"], path.parent)
            reference(meta["pdf"], path.parent)
            for field in ("source_bundle", "abstract_source"):
                if meta.get(field):
                    reference(meta[field], path.parent)
            assert meta["published"] is False, f"current revision mislabeled as deposited: {path}"
            assert meta["status"] in {"draft", "revision", "ready-for-review", "parked", "partial-draft"}
            assert meta["next_steps"] and all(isinstance(step, str) and step.strip() for step in meta["next_steps"])
            assert meta["verification"] and meta["build"]
            if collection == "results":
                assert meta.get("published_edition") == "published/", f"missing published archive: {path.parent}"
                archive = reference("published/", path.parent, directory=True)
                check_deposited_edition(reference("metadata.json", archive))
            elif meta.get("published_edition") and "://" not in meta["published_edition"]:
                reference(meta["published_edition"], path.parent, directory=meta["published_edition"].endswith("/"))
    return tuple(counts)


def markdown_documents() -> list[Path]:
    names = subprocess.check_output(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard", "--", "*.md"],
        cwd=REPOSITORY).decode().split("\0")
    return sorted({REPOSITORY / name for name in names if name and (REPOSITORY / name).is_file()})


def retained_json_records() -> list[Path]:
    """Ask Git for retained records without walking ignored computation caches."""
    names = subprocess.check_output(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard",
         "--", "research/**/*.json"], cwd=REPOSITORY).decode().split("\0")
    return sorted({REPOSITORY / name for name in names
                   if name and (REPOSITORY / name).is_file()})


def check_current_inputs(documents: list[Path] | None = None) -> int:
    """Validate explicit current manifests; old receipt digests remain snapshots."""
    checked = 0
    for document in retained_json_records() if documents is None else documents:
        name = document.name
        if not (name == "sources.json" or name.endswith(
                ("context.json", "verification.json", "receipt.json"))):
            continue
        data = json.loads(document.read_text())
        declarations = []
        for key in ("retained_sha256", "retained_proof_and_source_sha256"):
            for value, digest in data.get(key, {}).items():
                declarations.append((value, digest, ROOT))
        for item in data.get("local_sources", []):
            declarations.append((item["path"], item["sha256"], ROOT))
        replay = data.get("current_replay", {})
        if isinstance(replay, dict):
            for value, digest in replay.get("files", {}).items():
                declarations.append((value, digest, document.parent))
        for value, digest, base in declarations:
            path = reference(value, base)
            assert re.fullmatch(r"[0-9a-f]{64}", digest), f"invalid input digest: {document}"
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            assert actual == digest, f"changed current input: {document.relative_to(ROOT)} -> {value}"
            checked += 1
    return checked


def source_markdown_mirrors(ledgers: list[Path] | None = None) -> set[Path]:
    """Pinned upstream documents retain upstream-relative links and exact bytes."""
    mirrors = set()
    for ledger in ([p for p in retained_json_records() if p.name == "source_ledger.json"]
                   if ledgers is None else ledgers):
        data = json.loads(ledger.read_text())
        for item in data.get("cached_primary_sources", []):
            path = reference(item["file"], ledger.parent)
            assert item["url"].startswith("https://"), f"source mirror lacks upstream URL: {path}"
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            assert actual == item["sha256"], f"changed pinned upstream Markdown source: {path}"
            if path.suffix == ".md":
                mirrors.add(path)
    return mirrors


def check_links() -> tuple[int, int]:
    missing, checked = [], 0
    mirrors = source_markdown_mirrors()
    documents = [path for path in markdown_documents() if path.resolve() not in mirrors]
    fence = re.escape(chr(96) * 3)
    for path in documents:
        text = re.sub(r"(?ms)^ *" + fence + r".*?^ *" + fence + r"[^\n]*", "", path.read_text())
        text = re.sub(r"(?ms)^ *~~~.*?^ *~~~[^\n]*", "", text)
        text = re.sub(r"(?m)^(?: {4}|\t).*", "", text)
        tick = re.escape(chr(96))
        text = re.sub(r"(?s)(" + tick + r"+).*?\1", "", text)
        text = re.sub(r"(?s)\\\(.*?\\\)|\\\[.*?\\\]", "", text)
        text = re.sub(r"(?s)(?<!\\)" + re.escape(chr(36)) + r"{1,2}.*?(?<!\\)" + re.escape(chr(36)) + r"{1,2}", "", text)
        for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", text):
            target = target.strip()
            if target.startswith(("#", "mailto:", "data:")) or "://" in target:
                continue
            target = target.split(' "', 1)[0].strip("<>").split("#", 1)[0]
            if not target or ("." not in target and "/" not in target):
                continue
            checked += 1
            if not (path.parent / target).exists():
                missing.append(f"{path.relative_to(REPOSITORY)} -> {target}")
    assert not missing, "missing Markdown references:\n" + "\n".join(missing)
    return checked, len(documents)


def check_tex_proof_links(data: dict, graph: DependencyGraph) -> dict:
    report = tex_proof_link_report(data, graph)
    assert not report["issues"], "unresolved current TeX proof references: " + json.dumps(report["issues"][:20])
    stored = data.get("graph", {}).get("proof_links", {}).get("links", [])
    indexed = {(item["source_path"], item["locator"], item["conclusion"], item["dependency"]): item
               for item in stored}
    for link in report["links"]:
        key = (link["source_path"], link["locator"], link["conclusion"], link["dependency"])
        assert key in indexed, "unregistered TeX proof implication: " + json.dumps(link)
        routes = indexed[key].get("routes", [])
        assert routes, "TeX proof implication lacks an attached route: " + json.dumps(link)
        for rule_id in routes:
            assert rule_id in graph.rules and link["dependency"] in graph.rules[rule_id]["premises"], (
                "TeX prerequisite bypassed: " + rule_id + " needs " + link["dependency"])
    return report


def render_bounds(data: dict, graph: DependencyGraph) -> str:
    rows = bounds_records(graph)
    lines = ["# Quantitative claims and bounds", "",
             "<!-- Generated from MATH_STATUS.json; run make render-navigation. -->", "",
             "Bounds retain their exact domains and assumptions. Numeric comparisons use identical "
             "quantities, domains and parameter restrictions; conditional bounds stay separate. "
             "Symbolic criteria are displayed without guessing implication or asymptotic order. "
             "Bounds on an explicitly excluded domain remain visible as intermediate proof steps. "
             "A current proof of exclusion removes those rows from comparisons.", "",
             "[Dependency graph](DEPENDENCIES.md) · [Source coverage](COVERAGE.md) · [Exact scopes](SCOPES.md)", "",
             "| Claim | State | Quantity | Domain | Bound | Assumptions |", "|---|---|---|---|---|---|"]
    for row in rows:
        value = str(row["value"]) if "value" in row else row["expression"]
        domain = row["domain"]
        if row.get("parameters"):
            domain += "; " + json.dumps(row["parameters"], ensure_ascii=False, sort_keys=True)
        if row.get("domain_excluded_by"):
            label = "domain excluded by " if row["domain_excluded"] else "domain exclusion unavailable: "
            domain += " (" + label + row["domain_excluded_by"] + ")"
        lines.append("| " + markdown_cell(row["id"]) + " | " + row["state"] + " | "
                     + markdown_cell(row["quantity"]) + " | " + markdown_cell(domain) + " | "
                     + row["sense"] + ": " + markdown_cell(value) + " | "
                     + markdown_cell(", ".join(row["assumptions"])) + " |")
    comparisons = bound_comparisons(graph)
    opportunities = [item for row in comparisons for item in row["opportunities"]]
    lines += ["", "## Comparable improvements", ""]
    if opportunities:
        for item in opportunities:
            lines.append("- " + item["weaker_claim"] + ": " + str(item["old_value"]) + " -> "
                         + str(item["stronger_value"]) + " via " + ", ".join(item["stronger_claims"])
                         + "; structural consumers: " + (", ".join(item["affected"]) or "none") + ".")
    else:
        lines.append("No weaker registered numeric bound has an available stronger bound in an identical domain.")
    for row in comparisons:
        if row["conflict"]:
            lines += ["", "**Conflicting available bounds:** " + row["quantity"] + " in " + row["domain"] + "."]
    return "\n".join(lines) + "\n"


def graph_selection(graph: DependencyGraph, *, focus: str | None = None,
                    topic: str | None = None, depth: int | None = None,
                    direction: str = "both") -> set[str]:
    selected = set(graph.nodes)
    if topic:
        selected = {node_id for node_id, node in graph.nodes.items()
                    if node.get("topic") == topic or topic in node.get("topics", [])}
    if focus:
        if focus not in graph.nodes:
            raise GraphValidationError("unknown focus node: " + focus)
        seen, pending = {focus}, [(focus, 0)]
        while pending:
            node_id, distance = pending.pop(0)
            if depth is not None and distance >= depth:
                continue
            neighbors = set()
            if direction in {"both", "ancestors"}:
                for rule_id in graph.incoming[node_id]:
                    neighbors.update(graph.rules[rule_id]["premises"])
            if direction in {"both", "descendants"}:
                neighbors.update(graph.outgoing[node_id])
            for neighbor in sorted(neighbors - seen):
                seen.add(neighbor)
                pending.append((neighbor, distance + 1))
        selected &= seen
    return selected


def export_selected_dot(graph: DependencyGraph, selected: set[str]) -> str:
    if selected == set(graph.nodes):
        return graph.to_dot()
    quote = json.dumps
    evaluated = graph.evaluate()
    lines = ["digraph research {", "  rankdir=LR;"]
    for node_id in graph.order:
        if node_id in selected:
            lines.append("  " + quote(node_id) + " [shape=box,label="
                         + quote(node_id + "\n" + evaluated[node_id]["state"]) + "];")
    for rule_id, rule in graph.rules.items():
        visible = [premise for premise in rule["premises"] if premise in selected]
        if rule["conclusion"] not in selected or not visible:
            continue
        missing = len(rule["premises"]) - len(visible)
        label = rule_id + " (AND)" + (f" +{missing} omitted required premises" if missing else "")
        lines.append("  " + quote(rule_id) + " [shape=diamond,label=" + quote(label) + "];")
        lines += ["  " + quote(premise) + " -> " + quote(rule_id) + ";" for premise in visible]
        lines.append("  " + quote(rule_id) + " -> " + quote(rule["conclusion"]) + ";")
    lines.append("}")
    return "\n".join(lines) + "\n"


def print_explanation(record: dict) -> None:
    print(f"{record['id']}: {record['state']} — {record['title']}")
    print(record["scope"])
    print("Authored state: " + str(record["recorded_state"])
          + "; dependency coverage: " + record["dependency_coverage"])
    if record["assumptions"]:
        print("Assumptions: " + ", ".join(record["assumptions"]))
    seen = set()

    def provenance(current: dict, indent: str = "") -> None:
        identity = (current["id"], tuple(current["assumptions"]),
                    tuple((item["type"], item["id"]) for item in current["support"]))
        if identity in seen:
            print(indent + current["id"] + " (this proof and assumptions shown above)")
            return
        seen.add(identity)
        if indent:
            print(indent + current["id"] + ": " + current["state"])
        for source in current.get("sources", []):
            print(indent + "  Source: " + source["path"] + " :: " + source["locator"])
            print(indent + "    sha256:" + source["sha256"])
        for route in current["routes"]:
            print(indent + "  " + route["id"] + ": " + route["state"]
                  + " <= " + " AND ".join(route["premises"]))
            if route.get("proof"):
                proof = route["proof"]
                print(indent + "    " + proof["path"] + " :: " + proof["locator"])
                print(indent + "    sha256:" + proof["sha256"])
        for premise in current.get("premise_provenance", []):
            provenance(premise, indent + "    ")

    provenance(record)
    for blocker in record["blockers"]:
        print("  Blocker: " + blocker["id"] + " (" + blocker["state"] + ") via "
              + " -> ".join(blocker["path"]))
    if record.get("stale_sources"):
        print("Stale evidence: " + ", ".join(item["path"] for item in record["stale_sources"]))
    print("Internal graph status; assurance and publication flags are unchanged.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("render", help="render status, scope audit, DAG and topic guides")
    commands.add_parser("check", help="validate sources, graph, freshness, portfolio and local links")
    commands.add_parser("list", help="list the current research topics")
    show = commands.add_parser("show", help="read one topic, principal result or graph node")
    show.add_argument("id")
    search = commands.add_parser("search", help="search current topic, result and graph metadata")
    search.add_argument("query")
    explain = commands.add_parser("explain", help="show proof routes, assumptions and missing premises")
    explain.add_argument("id")
    explain.add_argument("--json", action="store_true")
    impact = commands.add_parser("impact", help="list descendants affected by a node or source")
    impact.add_argument("id", nargs="?")
    impact.add_argument("--source", help="research-relative evidence path")
    impact.add_argument("--json", action="store_true")
    simulate = commands.add_parser("simulate", help="read-only hypothetical propagation")
    for option in ("prove", "reopen", "falsify", "stale"):
        simulate.add_argument("--" + option, action="append", default=[], metavar="ID")
    simulate.add_argument("--json", action="store_true")
    graph_command = commands.add_parser("graph", help="export the complete implication DAG")
    graph_command.add_argument("--format", choices=("dot", "json"), default="dot")
    graph_command.add_argument("--focus", help="show a claim and its neighborhood")
    graph_command.add_argument("--topic")
    graph_command.add_argument("--depth", type=int, help="maximum edge distance from --focus")
    graph_command.add_argument("--direction", choices=("both", "ancestors", "descendants"), default="both")
    coverage = commands.add_parser("coverage", help="inspect source and declaration mapping completeness")
    coverage.add_argument("--json", action="store_true")
    references = commands.add_parser("references", help="inspect cited papers and external input verification")
    references.add_argument("query", nargs="?")
    references.add_argument("--json", action="store_true")
    source = commands.add_parser("source", help="inspect source declarations, mapping and downstream impact")
    source.add_argument("path", help="research-relative path")
    source.add_argument("--json", action="store_true")
    opportunities = commands.add_parser("opportunities", help="collect changed inputs, missing links and comparable bounds")
    opportunities.add_argument("--json", action="store_true")
    bounds = commands.add_parser("bounds", help="inspect scoped bounds and comparable improvements")
    bounds.add_argument("--opportunities", action="store_true")
    bounds.add_argument("--json", action="store_true")
    gaps = commands.add_parser("gaps", help="list scoped mathematical, formal and workflow gaps")
    gaps.add_argument("--topic")
    gaps.add_argument("--kind", choices=("mathematical", "formalization", "verification", "workflow"))
    gaps.add_argument("--json", action="store_true")
    args = parser.parse_args()
    data = load()
    try:
        graph = dependency_graph(data)
        if args.command in {"render", "check"}:
            topics, results = validate(data)
            if args.command == "render":
                for path, content in navigation_outputs(data).items():
                    path.write_text(content)
                print(f"Rendered status, scope audit, DAG and {topics} topic guides ({results} results)")
            else:
                covered = references = None
                if data["schema_version"] >= 3 or data.get("graph", {}).get("coverage"):
                    covered = check_coverage(data, graph)
                if data["schema_version"] >= 3 or data.get("graph", {}).get("external_references"):
                    references = check_references(data, graph)
                linked = None
                if data["schema_version"] >= 3:
                    linked = check_tex_proof_links(data, graph)
                conflicts = [row for row in bound_comparisons(graph) if row["conflict"] and row["state"] == "proved"]
                assert not conflicts, "conflicting current bounds: " + json.dumps(conflicts)
                missing = missing_source_closure(graph)
                assert not missing, "unpinned local proof includes/imports: " + json.dumps(missing)
                stale = graph.stale_evidence()
                if stale:
                    paths = sorted({item["path"] for item in stale})
                    parser.error("changed graph evidence in " + str(len(paths)) + " files / "
                                 + str(len(stale)) + " declarations: " + ", ".join(paths)
                                 + ". Use impact --source PATH and inspect the proof before accepting new pins.")
                check_navigation(data)
                pinned_inputs = check_current_inputs()
                published, drafts = check_portfolio()
                links, documents = check_links()
                print(f"PASS current research: {topics} topics, {results} principal results; sources, evidence and checker hashes")
                print(f"PASS dependency DAG: {len(graph.nodes)} nodes, {len(graph.rules)} implications; quotient acyclic, pinned and current")
                if covered:
                    print(f"PASS research inventory: {len(covered['documents'])} sources, {covered['declarations']} indexed declarations; no omissions")
                if linked:
                    print(f"PASS explicit TeX proof links: {len(linked['links'])}; cross-input prerequisites attached")
                if references:
                    print(f"PASS cited resources: {len(references['papers'])} entries; "
                          f"{references['bibliography_items']} bibliography items accounted for")
                print(f"PASS portfolio: {published} released papers with sealed published archives, {drafts} unpublished manuscripts")
                print(f"PASS current input manifests: {pinned_inputs} pinned references; cached primary sources checked")
                print(f"PASS local links: {links} references in {documents} maintained Markdown files; no proof replay")
        elif args.command == "list":
            for topic in data["topics"]:
                print(topic["id"] + ": " + topic["title"] + " — " + topic["status"])
        elif args.command == "show":
            records = [*data["topics"], *(result for topic in data["topics"] for result in topic["results"]),
                       *data.get("graph", {}).get("nodes", [])]
            selected = next((record for record in records if record["id"] == args.id), None)
            if selected is None:
                parser.error("unknown current topic or result: " + args.id)
            print(json.dumps(selected, indent=2, ensure_ascii=False))
        elif args.command == "search":
            words = args.query.lower().split()
            for record in [*data["topics"], *data.get("graph", {}).get("nodes", [])]:
                if all(word in json.dumps(record, ensure_ascii=False).lower() for word in words):
                    print(record["id"] + ": " + record.get("title", record.get("label", "")))
        elif args.command == "graph":
            assert graph is not None, "no dependency graph"
            if args.depth is not None and (args.depth < 0 or not args.focus):
                parser.error("--depth requires --focus and a nonnegative value")
            if args.topic and args.topic not in {topic["id"] for topic in data["topics"]}:
                parser.error("unknown topic: " + args.topic)
            selected = graph_selection(graph, focus=args.focus, topic=args.topic, depth=args.depth,
                                       direction=args.direction)
            if args.format == "dot":
                print(export_selected_dot(graph, selected), end="")
            else:
                evaluated = graph.evaluate()
                rules = [rule for rule in graph.rules.values() if rule["conclusion"] in selected
                         and any(premise in selected for premise in rule["premises"])]
                print(json.dumps({"nodes": {node_id: evaluated[node_id] for node_id in graph.order if node_id in selected},
                                  "definitions": {node_id: graph.nodes[node_id] for node_id in selected},
                                  "implications": rules,
                                  "quotient": graph.quotient() if hasattr(graph, "quotient") else None,
                                  "selection_is_complete": selected == set(graph.nodes),
                                  "omitted_premises": sorted({premise for rule in rules for premise in rule["premises"]
                                                             if premise not in selected})},
                                 indent=2, ensure_ascii=False))
        elif args.command == "explain":
            assert graph is not None, "no dependency graph"
            record = graph.explain(args.id)
            if args.json:
                print(json.dumps(record, indent=2, ensure_ascii=False))
            else:
                print_explanation(record)
        elif args.command == "impact":
            assert graph is not None, "no dependency graph"
            if bool(args.id) == bool(args.source):
                parser.error("impact requires either ID or --source PATH")
            affected = graph.source_impact(args.source) if args.source else graph.impact(args.id)
            evaluation = graph.evaluate()
            records = [evaluation[node_id] for node_id in affected]
            if args.json:
                print(json.dumps(records, indent=2, ensure_ascii=False))
            else:
                print("Structural descendants (including blocked alternative routes):")
                for record in records:
                    print(record["id"] + ": " + record["state"] + " — " + record["title"])
        elif args.command == "simulate":
            assert graph is not None, "no dependency graph"
            overrides = {}
            for option, state in (("prove", "proved"), ("reopen", "open"),
                                  ("falsify", "falsified"), ("stale", "stale")):
                for node_id in getattr(args, option):
                    if node_id in overrides:
                        parser.error("multiple overrides for " + node_id)
                    overrides[node_id] = state
            if not overrides:
                parser.error("simulate requires at least one --prove/--reopen/--falsify/--stale ID")
            result = graph.simulate(overrides)
            if args.json:
                print(json.dumps(result, indent=2, ensure_ascii=False))
            else:
                print("HYPOTHETICAL — assumes overrides; catalogue and assurance flags unchanged.")
                for change in result["changes"]:
                    print(json.dumps(change, ensure_ascii=False))
                if not result["changes"]:
                    print("No effective state changes.")
        elif args.command == "coverage":
            assert graph is not None, "no dependency graph"
            state = coverage_state(data, graph)
            if args.json:
                print(json.dumps(state, indent=2, ensure_ascii=False))
            else:
                print(f"{len(state['documents'])}/{state['expected_documents']} source mappings; "
                      f"{state['declarations']} indexed declarations; {len(state['errors'])} issues")
                for error in state["errors"]:
                    print(error)
        elif args.command == "source":
            assert graph is not None, "no dependency graph"
            path = graph._path(args.path)
            assert path.is_file(), "missing source: " + args.path
            mapping = next((row for row in data.get("graph", {}).get("coverage", {}).get("documents", [])
                            if graph._path(row["path"]) == path), None)
            impact = graph.source_impact(args.path)
            evaluation = graph.evaluate()
            info = {"path": args.path, "current_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    "mapping": mapping, "declarations": declarations(path),
                    "affected": [evaluation[node_id] for node_id in impact]}
            if args.json:
                print(json.dumps(info, indent=2, ensure_ascii=False))
            else:
                print(args.path + ": " + ("mapped" if mapping else "UNMAPPED") + "; "
                      + str(len(info["declarations"])) + " named declarations; "
                      + str(len(impact)) + " directly supported claims/descendants")
                for row in info["affected"]:
                    print(row["id"] + ": " + row["state"])
        elif args.command == "opportunities":
            assert graph is not None, "no dependency graph"
            evaluated = graph.evaluate()
            covered = coverage_state(data, graph)
            cited = reference_state(data, graph)
            linked = tex_proof_link_report(data, graph)
            comparisons = [row for row in bound_comparisons(graph, evaluated)
                           if row["opportunities"] or row["conflict"]]
            stale = graph.stale_evidence(evaluated)
            report = {"changed_evidence": stale, "coverage_issues": covered["errors"],
                      "reference_issues": cited["errors"], "proof_link_issues": linked["issues"],
                      "bound_improvements_or_conflicts": comparisons,
                      "pending_external_applications": [
                          {"id": paper["id"], "title": paper["title"], "verification": paper["verification"]}
                          for paper in cited["papers"] if any(use["role"] == "application_unverified"
                                                             for use in paper.get("uses", []))]}
            if args.json:
                print(json.dumps(report, indent=2, ensure_ascii=False))
            else:
                for key, items in report.items():
                    print(key.replace("_", " ") + ": " + str(len(items)))
                for row in comparisons:
                    print(json.dumps(row, ensure_ascii=False))
                print("Pending citations need application review; they are not automatically mathematical blockers.")
        elif args.command == "references":
            assert graph is not None, "no dependency graph"
            state = reference_state(data, graph)
            papers = state["papers"]
            if args.query:
                words = args.query.casefold().split()
                papers = [paper for paper in papers
                          if all(word in json.dumps(paper, ensure_ascii=False).casefold() for word in words)]
            if args.json:
                print(json.dumps({"papers": papers, "issues": state["errors"]}, indent=2, ensure_ascii=False))
            else:
                for paper in papers:
                    print(paper["id"] + " [" + paper["verification"] + "]: " + paper["citation"])
                for error in state["errors"]:
                    print("ISSUE: " + error)
        elif args.command == "bounds":
            assert graph is not None, "no dependency graph"
            records = bound_comparisons(graph) if args.opportunities else bounds_records(graph)
            if args.opportunities:
                records = [row for row in records if row["opportunities"] or row["conflict"]]
            if args.json:
                print(json.dumps(records, indent=2, ensure_ascii=False))
            else:
                for record in records:
                    print(json.dumps(record, ensure_ascii=False))
                if not records:
                    print("No matching bound records or comparable improvement opportunities.")
        elif args.command == "gaps":
            assert graph is not None, "no dependency graph"
            if args.topic and args.topic not in {topic["id"] for topic in data["topics"]}:
                parser.error("unknown topic: " + args.topic)
            evaluation = graph.evaluate()
            rows = []
            for topic in data["topics"]:
                if args.topic and topic["id"] != args.topic:
                    continue
                for node_id in topic.get("gap_ids", []):
                    node = graph.nodes[node_id]
                    kind = node.get("gap_kind", "mathematical")
                    if args.kind and kind != args.kind:
                        continue
                    rows.append({**evaluation[node_id], "topic": topic["id"], "node_kind": node["kind"], "kind": kind})
            if args.json:
                print(json.dumps(rows, indent=2, ensure_ascii=False))
            else:
                for row in rows:
                    print(f"{row['topic']}/{row['id']} [{row['kind']}, {row['state']}]: {row['scope']}")
    except (GraphValidationError, CoverageError, AssertionError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
