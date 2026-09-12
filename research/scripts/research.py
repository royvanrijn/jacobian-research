#!/usr/bin/env python3
"""Find prior research and render navigation. Standard library; never runs a checker.

Search reads the working tree and MATH_STATUS directly, not a cached search index.
Generated catalogues are navigation only; the typed claim ledger is authoritative.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re
import sys

from research_documents import ROOT, document_paths
from research_programmes import ACTIVE_AREAS, is_active, validate_programme

LESSONS = ROOT / "knowledge/lessons.json"
AREAS = {
    "core": "Keller core, cancellation and arithmetic",
    "geometry": "GVC, Gaussian moments, SIC and extended geometry",
    "hessian": "HC4, Hessian and Schur programmes",
    "plane-jc": "Plane Jacobian programme",
    "elliptic-curves": "Elliptic curves and rank jumps",
    "elkies-k3": "K3 constructions and lattice algorithms",
    "formal": "Formal verification",
    "papers": "Papers",
}
GENERATED = {"STATUS.md", "knowledge/ALGORITHMS.md", "knowledge/FAILED_ROUTES.md",
             "knowledge/WORK_LEDGER.md", "knowledge/LEGACY_WORK_REVIEW.md", "knowledge/PARTIAL_REVIEW.md",
             "knowledge/CHECKER_REVIEW.md"}


def ledger() -> dict:
    return json.loads((ROOT / "MATH_STATUS.json").read_text())


def lessons() -> list[dict]:
    return json.loads(LESSONS.read_text())["lessons"]


def area(path: str) -> str:
    if path.startswith('archive/non-elliptic/'):
        return area(path.removeprefix('archive/non-elliptic/'))
    if path == 'replay/CATALOGUE.md':
        return 'elliptic-curves'  # Two retained EC canonical derivations share this catalogue.
    parts = Path(path).parts
    first = parts[0] if parts else ""
    if first == "archive" and len(parts) > 1:
        return area("/".join(parts[1:]))
    if first in AREAS:
        return first
    if first == "extended-geometry":
        return "geometry"
    if first.startswith(("HC4", "HC5", "MENG_", "SCHUR_", "PROJECTIVE_GRADIENT", "JC2_HC4")):
        return "hessian"
    return "core"


def source_documents() -> list[Path]:
    """Exclude generated views/data and frozen source copies, not historical notes."""
    paths = []
    for path in document_paths():
        if ROOT not in path.parents:
            continue
        rel = path.relative_to(ROOT)
        if rel.parts[0] in {"index", "artifacts", "output"} or rel.parts[:2] == ("knowledge", "work"):
            continue
        if rel.as_posix() in GENERATED or rel.parts[:3] == ("elliptic-curves", "data", "research_curves"):
            continue
        if rel.name == "INVENTORY.md":
            continue
        paths.append(path)
    return paths


def title(text: str, fallback: str) -> str:
    return next((line.lstrip("# ").strip() for line in text.splitlines()
                 if line.startswith("# ")), fallback)


def claim_fingerprint(entry: dict) -> str:
    payload = "\0".join(entry[k] for k in ("state", "title", "scope"))
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def stale_claims(item: dict, entries: list[dict]) -> list[str]:
    known = {e["id"]: e for e in entries}
    return [i for i in item["claims"]
            if item["reviewed_claims"].get(i) != claim_fingerprint(known[i])]


def validate_lessons(items: list[dict], entries: list[dict], *, check_reviews: bool = False) -> None:
    known = {e["id"] for e in entries}
    identifiers = [item["id"] for item in items]
    assert len(identifiers) == len(set(identifiers)), "duplicate lesson ID"
    required = {"id", "area", "title", "when", "use", "avoid", "boundary",
                "reopen_when", "claims", "sources", "implementation", "keywords", "reviewed_claims"}
    for item in items:
        assert required <= set(item) <= required | {'programme_status'}, f"{item.get('id')}: invalid lesson schema"
        validate_programme(item)
        assert item["area"] in AREAS, f"{item['id']}: unknown area"
        for field in ("id", "title", "when", "use", "avoid", "boundary", "reopen_when"):
            assert isinstance(item[field], str) and item[field], f"{item['id']}: missing {field}"
        for field in ("claims", "sources", "implementation", "keywords"):
            assert isinstance(item[field], list) and all(isinstance(v, str) and v for v in item[field])
            assert len(item[field]) == len(set(item[field])), f"{item['id']}: repeated {field}"
        assert item["sources"], f"{item['id']}: lesson needs a source"
        assert set(item["claims"]) <= known, f"{item['id']}: unknown claim ID"
        assert isinstance(item["reviewed_claims"], dict) and set(item["reviewed_claims"]) == set(item["claims"]), (
            f"{item['id']}: every supporting claim needs a reviewed scope fingerprint"
        )
        assert all(isinstance(h, str) and re.fullmatch(r"[0-9a-f]{16}", h)
                   for h in item["reviewed_claims"].values()), f"{item['id']}: invalid review fingerprint"
        if check_reviews:
            stale = stale_claims(item, entries)
            assert not stale, f"{item['id']}: review changed claim scopes before updating fingerprints: {', '.join(stale)}"
        for value in item["sources"] + item["implementation"]:
            p = Path(value.split("#", 1)[0])
            assert not p.is_absolute() and ".." not in p.parts, f"{item['id']}: invalid path {p}"
            assert (ROOT / p).is_file(), f"{item['id']}: missing source {p}"


def cell(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


def render_algorithms(items: list[dict], entries: list[dict]) -> str:
    by_id = {e["id"]: e for e in entries}
    lines = ["# Algorithmic lessons", "",
             "<!-- Generated by scripts/research.py render from knowledge/lessons.json. -->", "",
             "Use these methods before designing another calculation. These are scoped engineering "
             "lessons, not additional theorem claims. Follow each linked source and the live "
             "[claim ledger](../MATH_STATUS.json). [Search and preflight](../KNOWLEDGE_BASE.md).", ""]
    for key, label in AREAS.items():
        group = [item for item in items if item["area"] == key and is_active(item)]
        if not group:
            continue
        lines.extend([f"## {label}", ""])
        for item in group:
            stale = stale_claims(item, entries)
            lines.extend([f"### {item['id']}: {item['title']}", "",
                          f"**When:** {item['when']}", "",
                          f"**Use:** {item['use']}", "",
                          f"**Avoid repeating:** {item['avoid']}", "",
                          f"**Boundary:** {item['boundary']}", "",
                          f"**Revisit when:** {item['reopen_when']}", ""])
            if stale:
                lines.extend(["**Review needed:** supporting claim scopes changed: " + ", ".join(stale) +
                              ". Consult their current sources before reusing this lesson.", ""])
            refs = [f"[{Path(p.split('#')[0]).stem}](../{p})" for p in item["sources"]]
            code = [f"[{Path(p).name}](../{p})" for p in item["implementation"]]
            claims = [f"`{i}` ({by_id[i]['state']})" for i in item["claims"]]
            lines.append("Sources: " + "; ".join(refs) + ".")
            if code:
                lines.extend(["", "Implementation: " + "; ".join(code) + "."])
            if claims:
                lines.extend(["", "Recorded claims: " + "; ".join(claims) + "."])
            lines.append("")
    return "\n".join(lines)


def render_outputs(entries: list[dict], items: list[dict], *, work_data: tuple | None = None,
                   resource_rows: list[dict] | None = None, partial_data: dict | None = None) -> dict[Path, str]:
    documents = source_documents()
    by_source = defaultdict(list)
    for entry in entries:
        by_source[entry["canonical_source"]].append(entry)
    groups = {key: [e for e in entries if area(e["canonical_source"]) == key] for key in AREAS}
    active_entries = [e for e in entries if is_active(e)]
    intro = ["# Elliptic-curve research catalogue", "",
             "<!-- Generated by scripts/research.py render; do not edit. -->", "",
             f"The active programme contains **{len(active_entries)} claims**. "
             f"All {len(entries)} registered claims remain indexed, including the "
             "[archived programmes](../archive/non-elliptic/README.md). "
             "State, proof assurance, scope and dependency edges come from "
             "[MATH_STATUS.json](../MATH_STATUS.json); a navigation row does not strengthen a claim.", "",
             "[Current queue](../STATUS.md) · [Unknowns and work](../knowledge/WORK_LEDGER.md) · "
             "[Structured records](resources.md) · [Methods and failed routes](../KNOWLEDGE_BASE.md) · "
             "[Replay guide](../REPRODUCE.md) · [History](../archive/README.md)", "",
             "| Area | Claims | Canonical sources |", "|---|---:|---:|"]
    output = {}
    for key, label in AREAS.items():
        group = sorted(groups[key], key=lambda e: (e["canonical_source"], e["id"]))
        if any(is_active(e) for e in group):
            intro.append(f"| [{label}]({key}.md) | {len(group)} | {len({e['canonical_source'] for e in group})} |")
        lines = [f"# {label}", "", "<!-- Generated by scripts/research.py render; do not edit. -->", "",
                 "[All research](README.md) · [Current queue](../STATUS.md) · "
                 "[Algorithmic lessons](../knowledge/ALGORITHMS.md)", "",
                 "Each title links to its canonical source. `proved` applies only to the scope "
                 "of that claim, including bounded experiments. `parked` can mean resolved or "
                 "deferred: consult the full scope and replacement edges.", "",
                 "From the repository root, inspect any ID without executing its checker:", "",
                 "```sh", "python3 research/scripts/research.py show ID", "```", "",
                 "| ID | Recorded state | Result / canonical source |", "|---|---|---|"]
        if group and not any(is_active(e) for e in group):
            lines[4:4] = ["**Archived programme.** These preserved claims are outside the active elliptic-curve programme. "
                          "Their mathematical states are unchanged. [Archive guide](../archive/non-elliptic/README.md).", ""]
        for e in group:
            replacements = "; replaced by " + ", ".join(e["replaced_by"]) if e["replaced_by"] else ""
            programme = "archived programme; " if not is_active(e) else ""
            lines.append(f"| `{e['id']}` | {programme}{e['state']}{replacements} | [{cell(e['title'])}](../{e['canonical_source']}) |")
        refs = [p for p in documents if area(p.relative_to(ROOT).as_posix()) == key
                and p.relative_to(ROOT).as_posix() not in by_source
                and "archive" not in p.relative_to(ROOT).parts]
        if refs:
            lines.extend(["", "## Other source notes and navigation", "",
                          "These documents are not canonical sources of a registered claim. "
                          "They include workflows, bounded reports and unregistered work; their "
                          "presence does not establish a theorem or an active campaign.", ""])
            for p in refs:
                name = p.relative_to(ROOT).as_posix()
                lines.append(f"- [{cell(title(p.read_text(), p.stem))}](../{name})")
        output[ROOT / f"index/{key}.md"] = "\n".join(lines) + "\n"
    intro.extend(["", "## Find previous work before computing", "", "```sh",
                  'python3 research/scripts/research.py search "pointed quartic"',
                  'python3 research/scripts/research.py search "q8" --history',
                  "python3 research/scripts/research.py show EC-CURVE302-RECOVERED-MW17-PARENT",
                  "python3 research/scripts/research.py routes --area elkies-k3", "```", "",
                  "Search defaults to elliptic curves, supporting K3 work and shared navigation. "
                  "`--history` includes other programmes and byte-preserved cleanup snapshots. Results label their "
                  "authority; an old note is never treated as a current status update.", "",
                  f"[Document inventory](documents.tsv) covers {len(documents)} maintained Markdown "
                  "documents, with source hashes, roles and claim IDs. Generated curve pages, "
                  "artifact snapshots and generated navigation are accessed through their own "
                  "inventories. Missing local computation outputs are not regenerated.", ""])
    output[ROOT / "index/README.md"] = "\n".join(intro)
    records = ["path\trole\tarea\tclaim_ids\tsha256"]
    for p in documents:
        name = p.relative_to(ROOT).as_posix()
        role = "historical" if "archive" in p.relative_to(ROOT).parts else ("canonical" if name in by_source else "reference")
        records.append("\t".join([name, role, area(name), ",".join(e["id"] for e in by_source[name]),
                                   hashlib.sha256(p.read_bytes()).hexdigest()]))
    output[ROOT / "index/documents.tsv"] = "\n".join(records) + "\n"
    output[ROOT / "knowledge/ALGORITHMS.md"] = render_algorithms(items, entries)
    routes = ["# Scoped failed and retired routes", "",
              "<!-- Generated by scripts/research.py render from MATH_STATUS.json. -->", "",
              "These are the claim ledger's existing exclusions and supersessions, not new "
              "prohibitions. Read each reason and witness at its exact scope. A resource stop "
              "or finite miss is not a universal obstruction. [Algorithmic lessons](ALGORITHMS.md) "
              "retain additional implementation decisions.", ""]
    known = {e["id"]: e for e in entries}
    for e in entries:
        if not is_active(e) or not e.get("forbidden_attack_classes"):
            continue
        routes.extend([f"## {e['id']}: {e['title']}", "",
                       f"Recorded programme state: **{e['state']}**. "
                       f"[Canonical scope](../{e['canonical_source']}).", ""])
        if e["replaced_by"]:
            routes.extend(["Replaced by: " + ", ".join(f"`{i}`" for i in e["replaced_by"]) + ".", ""])
        for a in e["forbidden_attack_classes"]:
            refs = ", ".join(f"[{i}](../{known[i]['canonical_source']}) ({known[i]['state']})"
                             for i in a["witnesses"])
            routes.extend([f"- **{a['attack']}** {a['reason']} Witnesses: {refs}.", ""])
    output[ROOT / "knowledge/FAILED_ROUTES.md"] = "\n".join(routes)
    if work_data is not None:
        import research_work
        output.update(research_work.render(*work_data, entries, ROOT))
    if resource_rows is not None:
        import research_resources
        output.update(research_resources.render(resource_rows, ROOT))
    if partial_data is not None:
        import research_partial_reviews
        output.update(research_partial_reviews.render(partial_data, entries, ROOT))
    return output


def query_terms(query: str) -> list[str]:
    return re.findall(r"[a-z0-9]+", query.lower())


def match_score(query: str, heading: str, body: str) -> int:
    terms = query_terms(query)
    if not terms:
        return 0
    h, b = heading.lower(), body.lower()
    if not all(t in h or t in b for t in terms):
        return 0
    return 1 + sum(12 if t in h else 1 for t in terms) + (30 if query.lower() in h else 0)


def excerpt(body: str, query: str, limit: int = 380) -> str:
    flat = re.sub(r"\s+", " ", body).strip()
    terms = query_terms(query)
    position = min((flat.lower().find(t) for t in terms if t in flat.lower()), default=0)
    start = max(0, position - 70)
    return ("…" if start else "") + flat[start:start + limit] + ("…" if start + limit < len(flat) else "")


def search(query: str, entries: list[dict], items: list[dict], *, selected_area: str | None = None,
           history: bool = False, work_data: tuple | None = None,
           resource_rows: list[dict] | None = None) -> list[dict]:
    results = []

    def add(kind, identifier, heading, body, source, group, **extra):
        if not history and extra.get('programme_status') == 'archived':
            return
        if selected_area and group != selected_area:
            return
        score = match_score(query, identifier + " " + heading, body)
        if not score:
            return
        score += {"lesson": 10, "claim": 8}.get(kind, 0)
        if query.lower() == identifier.lower():
            score += 10000
        results.append({"kind": kind, "id": identifier, "title": heading, "source": source,
                        "area": group, "score": score, "excerpt": excerpt(body, query), **extra})

    for e in entries:
        searchable = e["scope"] + " " + str(e.get("checker") or "") + " " + " ".join(e.get("software_lock", []))
        add("claim", e["id"], e["title"], searchable, e["canonical_source"], area(e["canonical_source"]),
            state=e["state"], replaced_by=e["replaced_by"], programme_status=e.get('programme_status', 'active'))
    for item in items:
        searchable = " ".join(f"{k.replace('_', ' ').capitalize()}: {item[k]}" for k in
                              ("when", "use", "avoid", "boundary", "reopen_when"))
        searchable += " " + " ".join(item["sources"] + item["implementation"] + item["keywords"])
        add("lesson", item["id"], item["title"], searchable,
            item["sources"][0], item["area"], review_needed=stale_claims(item, entries),
            programme_status=item.get('programme_status', 'active'))
    if work_data is not None:
        import research_work
        data, legacy = work_data
        for row in research_work.records(data, legacy, entries, history=history):
            body = " ".join(row[k] for k in ("next_step", "done_when", "prerequisites"))
            add("work", row["id"], row["title"], body, row["sources"][0], row["area"],
                disposition=row["disposition"], review_needed=row["review_needed"],
                programme_status=row.get('programme_status', 'active'))
        if history:
            for row in legacy["items"]:
                add("historical-work", row["id"], row["text"][:135],
                    row["text"] + " " + row["rationale"], legacy["source"], "core",
                    disposition=row["disposition"], current_actions=row["actions"])
    for row in resource_rows or []:
        add("resource", row["id"], row["title"], row["summary"], row["source"], row["area"],
            resource_kind=row["kind"], related_claims=row["related_claims"], boundary=row["boundary"],
            temporal_role=row["temporal_role"], programme_status=row.get('programme_status', 'active'))
    canonical = {e["canonical_source"] for e in entries}
    extra_sources = [ROOT / name for name in sorted(canonical)
                     if Path(name).suffix.lower() in {".tex", ".py", ".sage", ".lean", ".sing", ".txt"}
                     and (ROOT / name).is_file()]
    for p in [*source_documents(), *extra_sources]:
        name = p.relative_to(ROOT).as_posix()
        if not history and (name.startswith('archive/non-elliptic/') or (
                area(name) not in ACTIVE_AREAS and '/' in name)):
            continue
        if not history and name == 'replay/CATALOGUE.md':
            continue  # Mixed historical command text; its EC claims still surface above.
        text = p.read_text()
        kind = "historical" if "archive" in p.relative_to(ROOT).parts else "source" if name in canonical else "reference"
        add(kind, name, title(text, p.stem), text, name, area(name))
    if history:
        for manifest in sorted((ROOT / "archive").glob("repository-cleanup-*/MANIFEST.json")):
            for record in json.loads(manifest.read_text())["files"]:
                p = ROOT / record["preserved_path"]
                text = p.read_text()
                add("historical-snapshot", record["preserved_path"], record["original_path"], text,
                    record["preserved_path"], area(record["original_path"].removeprefix("research/")))
    # A verbatim phrase in a retired handoff must not outrank a relevant current
    # claim or method. Exact ID lookup still wins within the authoritative tier.
    priority = {"lesson": 0, "claim": 0, "source": 1, "work": 1, "resource": 2,
                "reference": 2, "historical-work": 3, "historical": 3, "historical-snapshot": 4}
    return sorted(results, key=lambda row: (priority[row["kind"]], -row["score"], row["kind"], row["id"]))


def show(identifier: str, entries: list[dict], items: list[dict], *, work_data: tuple | None = None,
         resource_rows: list[dict] | None = None, partial_data: dict | None = None) -> dict:
    by_id = {e["id"]: e for e in entries}
    if identifier in by_id:
        e = by_id[identifier]
        import research_partial_reviews
        return {"entry": e, "current_scope_fingerprint": claim_fingerprint(e),
                "partial_review": research_partial_reviews.detail(identifier, partial_data, entries, ROOT) if partial_data else None,
                "checker_absence_review": research_partial_reviews.detail(identifier, partial_data, entries, ROOT,
                    collection='checker_absence_reviews') if partial_data else None,
                "related_lessons": [{**i, "review_needed": stale_claims(i, entries)} for i in items
                                    if identifier in i["claims"] or e["canonical_source"] in i["sources"]],
                "used_by": [x["id"] for x in entries if identifier in x["dependencies"]],
                "narrowed_by": [x["id"] for x in entries if identifier in x.get("narrows_problems", [])],
                "related_work": [a["id"] for a in work_data[0]["actions"] if identifier in a["claims"]] if work_data else [],
                "note": "Recorded scope only. Checker paths are references, never executed by this command."}
    for item in items:
        if item["id"] == identifier:
            return {"lesson": item, "review_needed": stale_claims(item, entries),
                    "claims": [by_id[i] for i in item["claims"]],
                    "current_scope_fingerprints": {i: claim_fingerprint(by_id[i]) for i in item["claims"]}}
    if work_data is not None:
        import research_work
        result = research_work.show(identifier, *work_data, entries)
        if result is not None:
            return result
    for row in resource_rows or []:
        if row["id"] == identifier:
            return {k: v for k, v in row.items() if k != "summary"} | {
                "current_claims": [by_id[i] for i in row["related_claims"]]}
    raise ValueError(f"Unknown claim, lesson, work item or resource: {identifier}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    find = sub.add_parser("search", help="search claims, methods, proposed work, structured records and source text")
    find.add_argument("query")
    find.add_argument("--area", choices=AREAS)
    find.add_argument("--history", action="store_true", help="include archived programmes and preserved cleanup snapshots")
    find.add_argument("--limit", type=int, default=12)
    find.add_argument("--json", action="store_true")
    find.add_argument("--kind", choices=("claim", "lesson", "work", "resource", "source", "reference",
                                       "historical", "historical-work", "historical-snapshot"))
    detail = sub.add_parser("show", help="read full scope, evidence, dependencies and related lessons")
    detail.add_argument("id")
    route = sub.add_parser("routes", help="read scoped failed routes from the authority, without reruns")
    route.add_argument("--area", choices=AREAS)
    route.add_argument("--history", action="store_true")
    work = sub.add_parser("work", help="list unscheduled proposals or literal partial/parked records")
    work.add_argument("--area", choices=AREAS)
    work.add_argument("--history", action="store_true")
    work.add_argument("--kind", choices=("actions", "partial", "parked"), default="actions")
    work.add_argument("--json", action="store_true", help="include full scopes or completion gates")
    resource = sub.add_parser("resources", help="read maintained process, support and curve records")
    from research_resources import KINDS
    resource.add_argument("--kind", choices=KINDS)
    resource.add_argument("--unknown-only", action="store_true", help="only explicit UNKNOWN conductor records")
    resource.add_argument("--limit", type=int, default=12)
    resource.add_argument("--json", action="store_true")
    resource.add_argument("--history", action="store_true")
    sub.add_parser("render", help="regenerate navigation only")
    check = sub.add_parser("check", help="check lesson references and generated navigation freshness")
    check.add_argument("--require-partial-review", action="store_true", help="fail unless every current partial result has a fresh source-level review")
    args = parser.parse_args()
    entries, items = ledger()["entries"], lessons()
    validate_lessons(items, entries, check_reviews=args.command == "check")
    import research_work
    import research_resources
    import research_partial_reviews
    work_data = research_work.load_work(ROOT)
    research_work.validate(*work_data, entries, ROOT, check_reviews=args.command == "check",
                           require_coverage=args.command in {"render", "check"})
    resource_rows = research_resources.records(ROOT, entries)
    partial_data = research_partial_reviews.load(ROOT)
    research_partial_reviews.validate(partial_data, entries, ROOT, check_reviews=args.command == "check",
                                      require_complete=args.command == "check" and args.require_partial_review,
                                      require_checker_review=args.command == "check")
    if args.command in {"render", "check"}:
        outputs = render_outputs(entries, items, work_data=work_data, resource_rows=resource_rows, partial_data=partial_data)
        stale = []
        for path, expected in outputs.items():
            if args.command == "check":
                if not path.exists() or path.read_text() != expected:
                    stale.append(str(path.relative_to(ROOT)))
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(expected)
        if stale:
            print("Stale research navigation: " + ", ".join(stale) +
                  "; run python3 research/scripts/research.py render", file=sys.stderr)
            return 1
        print(f"PASS research navigation: {len(entries)} claims, {len(items)} sourced lessons, "
              f"{len(work_data[0]['actions'])} proposed actions, {len(resource_rows)} structured records, "
              f"{len(outputs)} generated views; no calculations")
        current_partial = {e['id'] for e in entries if is_active(e) and e['state'] == 'partial'}
        recorded = current_partial & {r['id'] for r in partial_data['reviews']}
        reviewed = current_partial & {r['id'] for r in partial_data['reviews']
                                      if not research_partial_reviews.stale(r, entries, ROOT)}
        print(f"Source-level partial review: {len(reviewed)}/{len(current_partial)}; "
              f"{len(current_partial - recorded)} still unreviewed; "
              f"{len(recorded - reviewed)} need reconciliation")
        current_absent = {e['id'] for e in entries if is_active(e) and e['checker'] is None}
        reviewed_absent = current_absent & {r['id'] for r in partial_data.get('checker_absence_reviews', [])
                                           if not research_partial_reviews.stale(r, entries, ROOT)}
        print(f"Checker-absence review: {len(reviewed_absent)}/{len(current_absent)}; "
              "classification only, no proof-assurance upgrade")
        unavailable = {path for review in partial_data['reviews']
                       for path in research_partial_reviews.unavailable_local_evidence(review, ROOT)}
        print(f"Local evidence receipts: {len(unavailable)} files unavailable in this checkout; "
              "source review is not a portable replay certificate")
        inherited = work_data[1]['items']
        completed = sum('resolution' in item for item in inherited)
        print(f"Inherited checklist: {completed}/{len(inherited)} snapshot-specific completion records; "
              f"{len(inherited) - completed} unfinished")
        print(f"Programme scope: {sum(is_active(e) for e in entries)} active EC/K3 claims; "
              f"{sum(not is_active(e) for e in entries)} archived claims. "
              "Archiving is not mathematical or checklist completion.")
    elif args.command == "show":
        try:
            print(json.dumps(show(args.id, entries, items, work_data=work_data, resource_rows=resource_rows,
                                  partial_data=partial_data), indent=2, ensure_ascii=False))
        except ValueError as error:
            parser.error(str(error))
    elif args.command == "routes":
        rows = [{"problem": e["id"], "state": e["state"], "source": e["canonical_source"],
                 "replaced_by": e["replaced_by"], "excluded_routes": e["forbidden_attack_classes"]}
                for e in entries if e.get("forbidden_attack_classes")
                and (args.history or is_active(e))
                and (not args.area or area(e["canonical_source"]) == args.area)]
        print(json.dumps(rows, indent=2, ensure_ascii=False))
    elif args.command == "work":
        rows = research_work.records(*work_data, entries, args.kind, args.area, history=args.history)
        if args.json:
            print(json.dumps(rows, indent=2, ensure_ascii=False))
        else:
            print(f"{len(rows)} {args.kind}. Actions are unscheduled; recorded states retain their exact scope.")
            for row in rows:
                print(f"[{row['disposition']}] {row['id']} — {row['title']}")
                if row.get("replaced_by"):
                    print("  Replaced by: " + ", ".join(row["replaced_by"]))
    elif args.command == "resources":
        if args.limit < 1:
            parser.error("--limit must be positive")
        if args.unknown_only and args.kind not in {None, "curve"}:
            parser.error("--unknown-only refers to curve conductor status; use --kind curve")
        rows = research_resources.filtered(resource_rows, args.kind, args.unknown_only, history=args.history)
        if args.json:
            print(json.dumps({"total": len(rows), "results": [
                {k: v for k, v in row.items() if k != "record"} for row in rows[:args.limit]]}, indent=2, ensure_ascii=False))
        else:
            print(f"{len(rows)} structured reference records; showing {min(len(rows), args.limit)}. Mathematical authority remains MATH_STATUS.json.")
            for row in rows[:args.limit]:
                print(f"[{row['temporal_role']} / {row['kind']}] {row['id']} — {row['title']}")
    else:
        if args.limit < 1:
            parser.error("--limit must be positive")
        rows = search(args.query, entries, items, selected_area=args.area, history=args.history,
                      work_data=work_data, resource_rows=resource_rows)
        if args.kind:
            rows = [r for r in rows if r["kind"] == args.kind]
        if args.json:
            print(json.dumps({"total": len(rows), "results": rows[:args.limit]}, indent=2, ensure_ascii=False))
        else:
            print(f"{len(rows)} matches; showing {min(len(rows), args.limit)}. State applies only to the recorded scope.")
            for r in rows[:args.limit]:
                state = f" / {r['state']}" if "state" in r else ""
                kind = r['kind']
                if r.get('programme_status') == 'archived':
                    kind = 'archived programme / ' + kind
                if kind == 'resource':
                    kind = f"{r['temporal_role']} resource / {r['resource_kind']}"
                print(f"\n[{kind}{state}] {r['id']}\n{r['title']}\n  research/{r['source']}\n  {r['excerpt']}")
                if r.get('temporal_role') == 'historical':
                    print('  Historical process record; read current claims before interpreting its frontier.')
                if r.get('related_claims'):
                    print('  Related current claim IDs: ' + ', '.join(r['related_claims']))
                if r.get("replaced_by"):
                    print("  Replaced by: " + ", ".join(r["replaced_by"]))
                if r.get("review_needed"):
                    print("  REVIEW NEEDED: supporting claim scopes changed: " + ", ".join(r["review_needed"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
