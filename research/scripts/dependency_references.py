"""Catalogue cited resources separately from mathematical theorem premises."""

import hashlib
from pathlib import Path
import re

VERIFICATIONS = {"source_read", "citation_only", "application_pending"}
USE_ROLES = {"context", "mathematical_input", "application_unverified"}


def bibliography_items(root, paths):
    items = set()
    for value in paths:
        path = root / value
        if path.suffix not in {".tex", ".bib", ".bbl"}:
            continue
        text = re.sub(r"(?m)(?<!\\)%.*$", "", path.read_text())
        if path.suffix == ".bib":
            for match in re.finditer(r"@([A-Za-z]+)\s*[({]\s*([^,\s]+)\s*,", text):
                if match.group(1).lower() not in {"comment", "preamble", "string"}:
                    items.add((value, match.group(2)))
        else:
            for match in re.finditer(r"\\bibitem(?:\s*\[[^]]*\])?\s*\{([^}]+)\}", text):
                items.add((value, match.group(1)))
    return items



def normalized_url(value):
    value = value.strip().rstrip(".,;:")
    while value.endswith(")") and value.count(")") > value.count("("):
        value = value[:-1]
    return value.rstrip("/").replace("http://", "https://", 1)


def markdown_external_links(graph):
    """Inspect current authored prose, preserving balanced DOI parentheses."""
    records = graph.data.get("graph", {}).get("coverage", {}).get("documents", [])
    found = set()
    for record in records:
        if Path(record["path"]).suffix != ".md" or record["role"] in {"documentation", "imported-source", "reference"}:
            continue
        path = graph._path(record["path"])
        if not path.is_file():
            continue
        text = path.read_text()
        fence = re.escape(chr(96) * 3)
        text = re.sub(r"(?ms)^\s*" + fence + r".*?^\s*" + fence + r"[^\n]*", "", text)
        for value in re.findall(r"https?://[^\s<>\[\]{}]+", text):
            url = normalized_url(value.strip(chr(96) + "*"))
            if "github.com/royvanrijn/" in url or "orcid.org/" in url:
                continue
            found.add((record["path"], url))
    return found


def reference_state(data, graph, paths=None):
    """Check identity, current retained pins and complete explicit bibitem coverage."""
    from dependency_coverage import inventory_paths
    catalogue = data.get("graph", {}).get("external_references", {})
    errors, papers, identifiers, accounted = [], catalogue.get("papers", []), set(), set()
    if catalogue.get("schema_version") != 1:
        errors.append("missing external reference catalogue")
    for paper in papers:
        paper_id = paper.get("id")
        if not isinstance(paper_id, str) or not paper_id or paper_id in identifiers:
            errors.append("empty/duplicate external reference id: " + str(paper_id))
        identifiers.add(paper_id)
        for key in ("title", "citation"):
            if not isinstance(paper.get(key), str) or not paper[key].strip():
                errors.append("missing external " + key + ": " + str(paper_id))
        if paper.get("verification") not in VERIFICATIONS:
            errors.append("invalid external verification: " + str(paper_id))
        for source in paper.get("local_sources", []):
            try:
                graph._source(source, paper_id)
                path = graph._path(source["path"])
                current = hashlib.sha256(path.read_bytes()).hexdigest()
                if current != source["sha256"]:
                    errors.append("changed external source: " + source["path"])
            except (OSError, ValueError) as error:
                errors.append(str(error))
        for use in paper.get("uses", []):
            path = graph._path(use.get("path", ""))
            if not path.is_file() or not str(use.get("locator", "")).strip():
                errors.append("missing external use source/locator: " + str(paper_id))
            if use.get("role") not in USE_ROLES:
                errors.append("invalid external use role: " + str(paper_id))
            for node_id in use.get("nodes", []):
                if node_id not in graph.nodes:
                    errors.append("unknown external use node: " + node_id)
            match = re.match(r"bib(?:item|entry)\s+(\S+)", use.get("locator", ""))
            if match:
                accounted.add((use["path"], match.group(1)))
    for use in catalogue.get("internal_citations", []):
        match = re.match(r"bib(?:item|entry)\s+(\S+)", use.get("locator", ""))
        if match:
            accounted.add((use["path"], match.group(1)))
    wanted = bibliography_items(graph.root, inventory_paths(graph.root, graph.repository) if paths is None else paths)
    for value, key in sorted(wanted - accounted):
        errors.append("unmapped bibliography item: " + value + " :: " + key)
    known_urls = {normalized_url(url) for paper in papers for url in paper.get("urls", [])}
    for value, url in sorted(markdown_external_links(graph)):
        if url not in known_urls:
            errors.append("unmapped external URL: " + value + " :: " + url)
    retained = {graph._path(source["path"]) for paper in papers
                for source in paper.get("local_sources", [])}
    for node in graph.nodes.values():
        if node.get("origin") != "external":
            continue
        for source in node.get("sources", []):
            path = graph._path(source["path"])
            if "/sources/" in source["path"] and path.suffix.lower() in {".pdf", ".tex", ".typ", ".html", ".txt"}:
                if path not in retained:
                    errors.append("external theorem primary/source file lacks cited-resource mapping: "
                                  + node["id"] + " :: " + source["path"])
    return {"papers": papers, "errors": errors, "complete": not errors,
            "bibliography_items": len(wanted), "internal_citations": catalogue.get("internal_citations", [])}


def check_references(data, graph, paths=None):
    state = reference_state(data, graph, paths)
    assert not state["errors"], "external reference catalogue incomplete/stale:\n" + "\n".join(state["errors"][:30])
    return state


def render_references(data, graph):
    state = reference_state(data, graph)
    cell = lambda value: str(value).replace("|", "\\|").replace("\n", " ")
    lines = ["# Cited research and external inputs", "",
             "<!-- Generated from MATH_STATUS.json; run make render-navigation. -->", "",
             "Citation identity, reading a source, verifying a particular application, and proving "
             "a dependent theorem are separate. Context citations supply no mathematical premise. "
             "A source_read entry retains the stated ledger scope; it does not certify every use.", "",
             f"**Catalogue:** {len(state['papers'])} cited resources; {state['bibliography_items']} current "
             "TeX bibliography entries accounted for, including own-paper citations tracked separately.", "",
             "[Graph](DEPENDENCIES.md) · [Source coverage](COVERAGE.md) · [Scope](SCOPES.md)", ""]
    if state["errors"]:
        lines += ["## Reference issues", ""]
        lines += ["- " + error for error in state["errors"]]
        lines += [""]
    lines += ["| Cited resource | Kind | Evidence status | Retained primary source | Uses / mapped inputs |",
              "|---|---|---|---|---|"]
    for paper in sorted(state["papers"], key=lambda row: row["title"].casefold()):
        urls = paper.get("urls", [])
        title = "[" + cell(paper["title"]) + "](" + urls[0] + ")" if urls else cell(paper["title"])
        sources = "; ".join("[" + cell(source["path"]) + "](" + source["path"] + ") :: "
                            + cell(source["locator"]) for source in paper.get("local_sources", []))
        uses = paper.get("uses", [])
        inputs = sorted({node_id for use in uses for node_id in use.get("nodes", [])})
        lines.append("| " + title + " | " + paper.get("resource_type", "paper") + " | "
                     + paper["verification"] + " | " + (sources or "not retained; citation/application boundary recorded")
                     + " | " + str(len(uses)) + " uses; " + cell(", ".join(inputs)) + " |")
    return "\n".join(lines) + "\n"
