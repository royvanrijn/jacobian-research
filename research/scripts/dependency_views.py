"""Generate human-readable views of the catalogue's implication graph.

No status is authored here. All node, implication, scope and assurance values
come from MATH_STATUS.json and the graph engine's current evaluation.
"""

from __future__ import annotations

from collections import Counter
from html import escape
from textwrap import wrap


def _cell(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _link(label: str, target: str) -> str:
    return "[" + _cell(label).replace("[", "\\[").replace("]", "\\]") + "](" + target + ")"


def _topic_nodes(topic: dict, graph) -> set[str]:
    selected = {result["id"] for result in topic["results"]}
    selected.update(node_id for node_id, node in graph.nodes.items()
                    if topic["id"] in node.get("topics", []) or node.get("topic") == topic["id"])
    selected.update(topic.get("gap_ids", []))
    return selected.intersection(graph.nodes)


def _blocker_text(record: dict) -> str:
    values = []
    seen = set()
    for blocker in record["blockers"]:
        identity = (blocker["type"], blocker["id"], blocker["state"])
        if identity not in seen:
            seen.add(identity)
            values.append(blocker["id"] + " (" + blocker["state"] + ")")
    return "; ".join(values) if values else "none in declared routes"


def _assurances(result: dict) -> str:
    names = (("independent_replay", "Independent replay"),
             ("formal_verification", "Formal verification"),
             ("external_review", "External review"))
    return "; ".join(label + ": " + ("recorded" if result.get(key) else "not recorded")
                     for key, label in names)


def _notes(topic: dict) -> list[str]:
    values = topic.get("scope_notes", [])
    if isinstance(values, str):
        values = [values]
    return [str(value) for value in values]


def _obligation_state(state: str) -> str:
    return "proved / discharged" if state == "proved" else state


def _progress(topic: dict, evaluated: dict) -> list[str]:
    records = [evaluated[node_id] for node_id in topic.get("gap_ids", []) if node_id in evaluated]
    if not records:
        return []
    counts = Counter(record["state"] for record in records)
    return ["**Computed obligation progress:** " + "; ".join(
        str(count) + " " + _obligation_state(state) for state, count in sorted(counts.items())) + ".", ""]


def _scope_notes(topic: dict, evaluated: dict) -> list[str]:
    discharged = any(evaluated.get(node_id, {}).get("state") == "proved"
                     for node_id in topic.get("gap_ids", []))
    lines = []
    for note in _notes(topic):
        lines += [("**Authored scope assessment:** " if discharged else "") + note, ""]
    return lines


def _question_lines(topic: dict, graph, evaluated: dict) -> list[str]:
    lines = []
    for question in topic.get("open_questions", []):
        node_id = question.get("node")
        node = graph.nodes.get(node_id, {})
        state = evaluated.get(node_id, {}).get("state", "unmapped")
        kind = node.get("gap_kind", question.get("kind", "mathematical"))
        label = (node_id + " — " if node_id else "") + kind + "; " + _obligation_state(state)
        assessment = "Authored question: " if state == "proved" else ""
        lines += ["- **" + label + ".** " + assessment + _link(question["question"], question["source"])]
    return lines


def _mermaid_id(value: str) -> str:
    # Keep long exact identifiers readable without making the diagram wider
    # than the page. Hyphens remain in the label; only line breaks are added.
    return "<br/>".join(escape(part, quote=True) for part in wrap(value, width=30))


def _mermaid(topic: dict, graph, evaluated: dict) -> list[str]:
    focus = _topic_nodes(topic, graph)
    # Include the ancestors needed to interpret every focused route, even when
    # an imported statement is shared with another topic.
    selected = set(focus)
    pending = list(focus)
    while pending:
        node_id = pending.pop()
        for rule_id in graph.incoming[node_id]:
            for premise in graph.rules[rule_id]["premises"]:
                if premise not in selected:
                    selected.add(premise)
                    pending.append(premise)
    hidden = {node_id for node_id in selected if graph.nodes[node_id]["kind"] == "source"
              and evaluated[node_id]["state"] == "proved" and not graph.incoming[node_id]
              and node_id not in topic.get("gap_ids", [])}
    shown = [node_id for node_id in graph.order if node_id in selected - hidden]
    relevant_rules = [rule_id for rule_id, rule in graph.rules.items()
                      if rule["conclusion"] in selected - hidden]
    if not shown:
        return ["No graph nodes are mapped to this topic.", ""]
    names = {node_id: "n" + str(index) for index, node_id in enumerate(shown)}
    rule_names = {rule_id: "r" + str(index) for index, rule_id in enumerate(relevant_rules)}
    lines = ["```mermaid", "flowchart LR"]
    for node_id in shown:
        record, node = evaluated[node_id], graph.nodes[node_id]
        label = _mermaid_id(node_id) + "<br/>" + record["state"]
        if len(graph.incoming[node_id]) > 1:
            label += "<br/>OR routes"
        if node["kind"] == "goal":
            lines.append(f'  {names[node_id]}(["{label}"])')
        else:
            lines.append(f'  {names[node_id]}["{label}"]')
        lines.append(f"  class {names[node_id]} {record['state']}")
    for rule_id in relevant_rules:
        rule = graph.rules[rule_id]
        route = next(route for route in evaluated[rule["conclusion"]]["routes"] if route["id"] == rule_id)
        omitted = sum(premise in hidden for premise in rule["premises"])
        label = "AND<br/>" + _mermaid_id(rule_id)
        if omitted:
            label += f"<br/>+{omitted} pinned source" + ("s" if omitted != 1 else "")
        label += "<br/>" + route["state"]
        lines.append(f'  {rule_names[rule_id]}{{"{label}"}}')
        state_class = {"blocked": "open", "unverified": "open", "refuted": "falsified"}.get(route["state"], route["state"])
        lines.append(f"  class {rule_names[rule_id]} {state_class}")
        for premise in rule["premises"]:
            if premise in names:
                lines.append(f"  {names[premise]} --> {rule_names[rule_id]}")
        lines.append(f"  {rule_names[rule_id]} --> {names[rule['conclusion']]}")
    lines += ["  classDef proved fill:#e6f4ea,stroke:#2e7d32",
              "  classDef conditional fill:#fff8e1,stroke:#a66b00",
              "  classDef open fill:#fff4e5,stroke:#bf6a00",
              "  classDef stale fill:#fce8e6,stroke:#c62828",
              "  classDef falsified fill:#ffebee,stroke:#8e0000",
              "  classDef unmapped fill:#f1f3f4,stroke:#607d8b", "```", ""]
    if hidden:
        lines += ["Current unconditional proof-source leaves are omitted from this diagram; "
                  "each diamond states how many source premises are hidden. All premises, "
                  "including hidden sources, are required. Use `explain` for their pins and locators.", ""]
    return lines


def render_dependency_report(data: dict, graph) -> str:
    """Render research/DEPENDENCIES.md from the sole status authority."""
    evaluated = graph.evaluate()
    counts = Counter(record["state"] for record in evaluated.values())
    lines = ["# Mathematical dependencies", "",
             "<!-- Generated by scripts/research.py render from MATH_STATUS.json. -->", "",
             "This is a generated view of the " + _link("current catalogue", "MATH_STATUS.json") + ". "
             "An AND diamond requires every premise; alternative incoming diamonds are OR routes. "
             "A checked route becomes available only when its premises and pinned proof sources are current. "
             "A falsified premise blocks its route; it does not falsify the conclusion.", "",
             "Graph state follows the declared dependencies. Dependency coverage separately records "
             "whether that model is complete, partial or unmapped. A source hash detects changed bytes; "
             "it does not check a mathematical argument. Authored state, formal verification, replay, "
             "external review and publication remain separate.", "",
             _link("Exact scopes and recorded assurances", "SCOPES.md") + " · "
             + _link("Result catalogue", "STATUS.md") + " · " + _link("Research topics", "README.md"), "",
             f"The graph contains {len(graph.nodes)} claim/source nodes and {len(graph.rules)} implication rules. "
             + "; ".join(f"{state}: {count}" for state, count in sorted(counts.items())) + ".", "",
             "Inspect the current graph from the repository root:", "", "```sh",
             "python3 research/scripts/research.py gaps",
             "python3 research/scripts/research.py explain ID",
             "python3 research/scripts/research.py impact ID",
             "python3 research/scripts/research.py simulate --prove ID",
             "python3 research/scripts/research.py graph --format dot", "```", "",
             "Replace `ID` with a node identifier below. Simulation is read only and explicitly "
             "hypothetical: assuming a node proved leaves its identifier among the conditional "
             "assumptions of dependent conclusions. It does not edit the catalogue or upgrade assurances.", ""]
    for topic in data["topics"]:
        lines += ["## " + topic["title"], "", topic["summary"], ""]
        lines += _progress(topic, evaluated)
        lines += _scope_notes(topic, evaluated)
        lines += ["| Principal claim | Graph state | Dependency coverage | Current blockers |",
                  "|---|---|---|---|"]
        for result in topic["results"]:
            record = evaluated[result["id"]]
            lines.append("| " + _link(result["id"], result["canonical_source"]) + " | " + record["state"]
                         + " | " + record["dependency_coverage"] + " | " + _cell(_blocker_text(record)) + " |")
        lines += [""]
        interesting = [node_id for node_id in graph.order if node_id in _topic_nodes(topic, graph)
                       and graph.nodes[node_id]["kind"] in {"gap", "goal"}]
        if interesting:
            lines += ["| Goal or obligation | Kind | Graph state | Coverage | Exact scope | Current blockers |",
                      "|---|---|---|---|---|---|"]
            for node_id in interesting:
                record, node = evaluated[node_id], graph.nodes[node_id]
                state = _obligation_state(record["state"]) if node["kind"] == "gap" else record["state"]
                kind = node.get("gap_kind", "mathematical") if node["kind"] == "gap" else node["kind"]
                lines.append("| " + _cell(node_id) + " | " + kind + " | " + state + " | "
                             + record["dependency_coverage"] + " | " + _cell(record["scope"])
                             + " | " + _cell(_blocker_text(record)) + " |")
            lines += [""]
        lines += _mermaid(topic, graph, evaluated)
    stale = graph.stale_evidence(evaluated)
    if stale:
        lines += ["## Changed or missing evidence", "",
                  "These declarations require inspection and revalidation. Updating a digest alone "
                  "does not repair the supporting argument.", "",
                  "| Owner | Source | Freshness |", "|---|---|---|"]
        lines += ["| " + _cell(item["owner"]) + " | " + _link(item["path"], item["path"])
                  + " | " + item["reason"] + " |" for item in stale]
        lines += [""]
    return "\n".join(lines)


def render_scope_report(data: dict, graph) -> str:
    """Render research/SCOPES.md without interpreting a claim beyond its scope."""
    evaluated = graph.evaluate()
    lines = ["# Claim scopes and assurance", "",
             "<!-- Generated by scripts/research.py render from MATH_STATUS.json. -->", "",
             "Each statement below retains its exact authored scope. The graph evaluates only "
             "declared proof routes, with current source pins; its state does not replace the "
             "authored statement or certify that every dependency has been identified. "
             "Dependency coverage is recorded independently.", "",
             "Assurance flags retain their recorded scope. A formal certificate for a component "
             "does not formalize an entire paper, and a passing graph check does not establish "
             "independent replay or external review. Publication remains in the paper metadata.", "",
             _link("Implication diagrams and graph commands", "DEPENDENCIES.md") + " · "
             + _link("Structured authority", "MATH_STATUS.json") + " · " + _link("Current results", "STATUS.md"), ""]
    for topic in data["topics"]:
        lines += ["## " + topic["title"], "", "**Topic scope:** " + topic["summary"], ""]
        lines += _progress(topic, evaluated)
        lines += _scope_notes(topic, evaluated)
        for result in topic["results"]:
            record = evaluated[result["id"]]
            lines += ["### " + result["id"], "",
                      _link(result["title"], result["canonical_source"]), "",
                      "**Exact scope:** " + result["scope"], "",
                      "**Authored state:** " + result["state"] + ". **Graph state:** " + record["state"]
                      + ". **Dependency coverage:** " + record["dependency_coverage"] + ".", "",
                      "**Recorded assurances:** " + _assurances(result) + ".", "",
                      "**Current blockers:** " + _blocker_text(record) + ".", ""]
            if record["assumptions"]:
                lines += ["**Conditional assumptions:** " + "; ".join(record["assumptions"]) + ".", ""]
            if result.get("verification_note"):
                lines += ["**Verification limits:** " + result["verification_note"], ""]
        gaps = topic.get("gap_ids", [])
        supplemental = [node_id for node_id in graph.order if node_id in _topic_nodes(topic, graph)
                        and graph.nodes[node_id]["kind"] not in {"result", "source"} and node_id not in gaps]
        if supplemental:
            lines += ["### Supporting statements and targets", "",
                      "| Identifier | Kind | Graph state / coverage | Exact scope | Pinned source |",
                      "|---|---|---|---|---|"]
            for node_id in supplemental:
                record, node = evaluated[node_id], graph.nodes[node_id]
                references = "; ".join(_link(source["path"], source["path"]) + ": " + _cell(source["locator"])
                                       for source in node.get("sources", []))
                references = references or ("see incoming proof routes" if graph.incoming[node_id] else "no pinned source")
                lines += ["| " + _cell(node_id) + " | " + node["kind"] + " | " + record["state"]
                          + " / " + record["dependency_coverage"] + " | " + _cell(record["scope"])
                          + " | " + references + " |"]
            lines += [""]
        if gaps:
            lines += ["### Topic obligations", ""]
            for node_id in gaps:
                if node_id in evaluated:
                    record = evaluated[node_id]
                    kind = graph.nodes[node_id].get("gap_kind", "mathematical")
                    lines += ["- **" + node_id + "** (" + kind + "; " + _obligation_state(record["state"])
                              + "): " + record["scope"]]
            lines += [""]
        if topic.get("open_questions"):
            lines += ["Tracked obligations from the catalogue:", ""]
            lines += _question_lines(topic, graph, evaluated)
            lines += [""]
    return "\n".join(lines)
