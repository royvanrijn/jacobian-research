"""Conservative mathematical implication DAG, using only the standard library.

The catalogue's authored ``state`` is not a proof edge. A checked implication
has AND premises; different implications into a claim are alternative routes.
File digests detect changed evidence and never establish mathematical validity.
Evaluation and counterfactual simulation do not modify the source catalogue or
its publication, review, replay, or formal-verification flags.
"""

from __future__ import annotations

from collections import defaultdict, deque
from copy import deepcopy
import hashlib
from pathlib import Path, PurePosixPath, PureWindowsPath
import re


NODE_STATES = {"proved", "open", "unverified", "partial", "falsified", "refuted"}
RULE_STATES = {"proved", "unverified", "refuted"}
EFFECTIVE_STATES = {"proved", "conditional", "open", "stale", "falsified", "unmapped"}
AVAILABLE = {"proved", "conditional"}


class GraphValidationError(ValueError):
    """The declared graph is malformed or circular."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GraphValidationError(message)


def _nonempty(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


class DependencyGraph:
    """Read an explicit ``graph`` section in the existing status catalogue.

    Paths in graph evidence records are relative to the research directory.
    ``../unpublished/...`` and ``../results/...`` are permitted inside the
    repository. The constructor validates graph structure, but changed/missing
    evidence is reported as stale by :meth:`evaluate`, not a schema error.
    """

    def __init__(self, data: dict, root: Path, repository: Path | None = None):
        self.data = deepcopy(data)
        self.root = Path(root).resolve()
        self.repository = Path(repository or self.root.parent).resolve()
        self.graph = self.data.get("graph", {})
        self.nodes: dict[str, dict] = {}
        self.rules: dict[str, dict] = {}
        self.incoming: dict[str, list[str]] = defaultdict(list)
        self.outgoing: dict[str, set[str]] = defaultdict(set)
        self.order: list[str] = []
        self.validate()

    def _path(self, value: str) -> Path:
        _require(_nonempty(value), "empty graph source path")
        _require(not PurePosixPath(value).is_absolute()
                 and not PureWindowsPath(value).drive,
                 f"absolute graph source path: {value}")
        path = (self.root / value).resolve()
        _require(path.is_relative_to(self.repository),
                 f"graph source escapes repository: {value}")
        return path

    def _source(self, source: dict, owner: str) -> None:
        _require(isinstance(source, dict), f"invalid source record: {owner}")
        self._path(source.get("path", ""))
        _require(bool(re.fullmatch(r"[0-9a-f]{64}", str(source.get("sha256", "")))),
                 f"invalid source digest: {owner}")
        _require(_nonempty(source.get("locator")), f"missing exact source locator: {owner}")

    @staticmethod
    def _assumptions(item: dict, owner: str) -> None:
        values = item.get("assumptions", [])
        _require(isinstance(values, list) and all(_nonempty(value) for value in values),
                 f"invalid assumptions: {owner}")
        _require(len(values) == len(set(values)), f"duplicate assumptions: {owner}")
        _require(item.get("dependency_coverage", "unmapped") in {"complete", "partial", "unmapped"},
                 f"invalid dependency coverage: {owner}")

    def validate(self) -> None:
        """Reject malformed references, duplicates and every dependency cycle."""
        self.nodes.clear()
        self.rules.clear()
        self.incoming.clear()
        self.outgoing.clear()
        self.order.clear()
        _require(isinstance(self.graph, dict), "graph must be an object")
        if self.graph:
            _require(self.graph.get("schema_version") == 1, "unsupported graph schema_version")
        for topic in self.data.get("topics", []):
            for result in topic.get("results", []):
                node_id = result.get("id")
                _require(_nonempty(node_id), "empty result node id")
                _require(node_id not in self.nodes, f"duplicate node id: {node_id}")
                self.nodes[node_id] = {**result, "kind": "result", "topic": topic["id"],
                                       "recorded_state": result.get("state"), "supplemental": False}
                self._assumptions(result, node_id)
        supplemental = self.graph.get("nodes", [])
        implications = self.graph.get("implications", [])
        _require(isinstance(supplemental, list), "graph.nodes must be a list")
        _require(isinstance(implications, list), "graph.implications must be a list")
        for node in supplemental:
            _require(isinstance(node, dict), "invalid graph node")
            node_id = node.get("id")
            _require(_nonempty(node_id), "empty graph node id")
            _require(node_id not in self.nodes, f"duplicate node id: {node_id}")
            _require(_nonempty(node.get("kind")), f"missing node kind: {node_id}")
            _require(_nonempty(node.get("title", node.get("label"))), f"missing node title: {node_id}")
            _require(_nonempty(node.get("scope")), f"missing exact node scope: {node_id}")
            _require(node.get("status") in NODE_STATES, f"invalid node status: {node_id}")
            self._assumptions(node, node_id)
            sources = node.get("sources", [])
            _require(isinstance(sources, list), f"invalid node sources: {node_id}")
            _require(node.get("status") not in {"proved", "falsified", "refuted"} or bool(sources),
                     f"established node needs pinned sources: {node_id}")
            for source in sources:
                self._source(source, node_id)
            self.nodes[node_id] = {**node, "recorded_state": node["status"], "supplemental": True}
        for rule in implications:
            _require(isinstance(rule, dict), "invalid implication record")
            rule_id = rule.get("id")
            _require(_nonempty(rule_id), "empty implication id")
            _require(rule_id not in self.rules, f"duplicate implication id: {rule_id}")
            _require(rule_id not in self.nodes, f"implication id collides with node: {rule_id}")
            premises = rule.get("premises")
            _require(isinstance(premises, list) and bool(premises),
                     f"implication needs nonempty AND premises: {rule_id}")
            _require(all(_nonempty(node_id) for node_id in premises), f"invalid premises: {rule_id}")
            _require(len(premises) == len(set(premises)), f"duplicate premise: {rule_id}")
            conclusion = rule.get("conclusion")
            _require(_nonempty(conclusion), f"missing conclusion: {rule_id}")
            _require(rule.get("status") in RULE_STATES, f"invalid implication status: {rule_id}")
            self._assumptions(rule, rule_id)
            for node_id in [*premises, conclusion]:
                _require(node_id in self.nodes, f"dangling node {node_id}: {rule_id}")
            _require(conclusion not in premises, f"self dependency: {rule_id}")
            if rule.get("status") == "proved":
                self._source(rule.get("proof"), rule_id)
            elif "proof" in rule:
                self._source(rule["proof"], rule_id)
            self.rules[rule_id] = rule
            self.incoming[conclusion].append(rule_id)
            for premise in premises:
                self.outgoing[premise].add(conclusion)
        # Validate all routes, including unverified/refuted ones. A later status
        # edit must not activate a cycle that was hidden during validation.
        indegree = {node_id: 0 for node_id in self.nodes}
        for conclusions in self.outgoing.values():
            for conclusion in conclusions:
                indegree[conclusion] += 1
        pending = deque(sorted(node_id for node_id, degree in indegree.items() if degree == 0))
        while pending:
            node_id = pending.popleft()
            self.order.append(node_id)
            for conclusion in sorted(self.outgoing[node_id]):
                indegree[conclusion] -= 1
                if indegree[conclusion] == 0:
                    pending.append(conclusion)
        _require(len(self.order) == len(self.nodes),
                 "dependency cycle: " + ", ".join(sorted(node_id for node_id, degree in indegree.items() if degree)))

    def _freshness(self, sources: list[dict], cache: dict[Path, str | None]) -> list[dict]:
        stale = []
        for source in sources:
            path = self._path(source["path"])
            if path not in cache:
                try:
                    cache[path] = hashlib.sha256(path.read_bytes()).hexdigest()
                except OSError:
                    cache[path] = None
            actual = cache[path]
            if actual != source["sha256"]:
                stale.append({"path": source["path"], "expected_sha256": source["sha256"],
                              "actual_sha256": actual,
                              "reason": "missing" if actual is None else "changed"})
        return stale

    def evaluate(self, overrides: dict[str, str] | None = None) -> dict[str, dict]:
        """Calculate current states in topological order without editing records.

        Overrides are explicitly hypothetical. Assuming an open claim true
        leaves its identifier in the propagated assumption set; descendants
        remain conditional. A falsified premise disables a route and never
        proves the negation of its conclusion.
        """
        overrides = overrides or {}
        for node_id, state in overrides.items():
            _require(node_id in self.nodes, f"unknown simulation node: {node_id}")
            _require(state in {"proved", "open", "falsified", "stale"},
                     f"invalid simulation state: {node_id}={state}")
        evaluated = {}
        cache: dict[Path, str | None] = {}
        for node_id in self.order:
            node = self.nodes[node_id]
            sources = node.get("sources", []) if node["supplemental"] else []
            stale_sources = self._freshness(sources, cache)
            result = {"id": node_id, "kind": node["kind"],
                      "title": node.get("title", node.get("label", node_id)),
                      "scope": node.get("scope", ""), "recorded_state": node["recorded_state"],
                      "dependency_coverage": node.get("dependency_coverage", "unmapped"),
                      "state": "open" if node["supplemental"] else "unmapped",
                      "assumptions": [], "support": [], "routes": [],
                      "stale_sources": stale_sources, "blockers": [], "hypothetical": False}
            if node["supplemental"]:
                if stale_sources and node["status"] in {"proved", "falsified", "refuted"}:
                    result["state"] = "stale"
                elif node["status"] in {"falsified", "refuted"}:
                    result["state"] = "falsified"
                elif node["status"] == "proved" and not self.incoming[node_id]:
                    # Only actual leaves supply autonomous evidence. A proved
                    # declaration on an intermediate node must not bypass the
                    # premises of its recorded derivation. An independent proof
                    # can be represented by a separate source leaf and OR rule.
                    result["assumptions"] = sorted(node.get("assumptions", []))
                    result["state"] = "conditional" if result["assumptions"] else "proved"
                    result["support"] = [{"type": "leaf", "id": node_id, "sources": deepcopy(sources)}]
            elif node["recorded_state"] == "falsified":
                # A catalogue counterexample is retained as recorded, but is not
                # an available positive premise and cannot establish anything.
                result["state"] = "falsified"
            candidates = []
            for rule_id in self.incoming[node_id]:
                rule = self.rules[rule_id]
                stale_proof = self._freshness([rule["proof"]], cache) if "proof" in rule else []
                unavailable = [premise for premise in rule["premises"]
                               if evaluated[premise]["state"] not in AVAILABLE]
                assumptions = set(node.get("assumptions", [])) | set(rule.get("assumptions", []))
                for premise in rule["premises"]:
                    assumptions.update(evaluated[premise]["assumptions"])
                if rule["status"] == "refuted":
                    route_state = "refuted"
                elif rule["status"] != "proved":
                    route_state = "unverified"
                elif stale_proof:
                    route_state = "stale"
                elif unavailable:
                    route_state = "stale" if any(evaluated[p]["state"] == "stale" for p in unavailable) else "blocked"
                else:
                    route_state = "conditional" if assumptions else "proved"
                route = {"id": rule_id, "state": route_state, "premises": list(rule["premises"]),
                         "assumptions": sorted(assumptions), "unavailable_premises": unavailable,
                         "stale_sources": stale_proof, "proof": deepcopy(rule.get("proof")),
                         "hypothetical": any(evaluated[p]["hypothetical"] for p in rule["premises"])}
                result["routes"].append(route)
                if route_state in AVAILABLE:
                    candidates.append(route)
            # An intrinsic, current, unconditional leaf is another valid route.
            if result["state"] in AVAILABLE:
                candidates.append({"id": None, "state": result["state"],
                                   "assumptions": result["assumptions"]})
            if candidates and result["state"] != "falsified":
                best = min(candidates, key=lambda route: (len(route["assumptions"]), route["assumptions"], route["id"] or ""))
                result["state"] = best["state"]
                result["assumptions"] = best["assumptions"]
                result["support"] = ([{"type": "implication", "id": best["id"]}]
                                     if best["id"] else result["support"])
            elif self.incoming[node_id] and result["state"] != "falsified":
                result["state"] = ("stale" if result["state"] == "stale" or any(
                    route["state"] == "stale" for route in result["routes"]) else "open")
            if node_id in overrides:
                state = overrides[node_id]
                result["hypothetical"] = True
                result["state"] = "conditional" if state == "proved" else state
                result["assumptions"] = [node_id] if state == "proved" else []
                result["support"] = [{"type": "hypothesis", "id": node_id, "state": state}]
            elif result["support"] and result["support"][0]["type"] == "implication":
                selected = self.rules[result["support"][0]["id"]]
                result["hypothetical"] = any(evaluated[p]["hypothetical"] for p in selected["premises"])
            elif result["state"] not in AVAILABLE:
                result["hypothetical"] = any(route["hypothetical"] for route in result["routes"])
            evaluated[node_id] = result
        for node_id in self.order:
            evaluated[node_id]["blockers"] = self.blockers(node_id, evaluated)
        return evaluated

    def blockers(self, node_id: str, evaluation: dict[str, dict] | None = None) -> list[dict]:
        """Return unresolved leaves/rules for each alternative route, with paths.

        This is a route explanation, not a claim that all listed blockers are
        jointly necessary. One repaired alternative may suffice.
        """
        _require(node_id in self.nodes, f"unknown graph node: {node_id}")
        if evaluation is None:
            return self.evaluate()[node_id]["blockers"]
        found = []

        def visit(current: str, path: list[str]) -> None:
            record = evaluation[current]
            if record["state"] in AVAILABLE:
                for assumption in record["assumptions"]:
                    found.append({"type": "assumption", "id": assumption, "state": "open", "path": path})
                return
            routes = record["routes"]
            if not routes:
                found.append({"type": "node", "id": current, "state": record["state"], "path": path,
                              "stale_sources": record["stale_sources"]})
            for route in routes:
                route_path = [*path, route["id"]]
                if route["state"] in {"unverified", "refuted"} or route["stale_sources"]:
                    found.append({"type": "implication", "id": route["id"], "state": route["state"],
                                  "path": route_path, "stale_sources": route["stale_sources"]})
                for premise in route["unavailable_premises"]:
                    visit(premise, [*route_path, premise])
                for assumption in route["assumptions"]:
                    found.append({"type": "assumption", "id": assumption, "state": "open", "path": route_path})

        visit(node_id, [node_id])
        return found

    def impact(self, node_id: str) -> list[str]:
        """List all structural descendants, including currently blocked routes."""
        _require(node_id in self.nodes, f"unknown graph node: {node_id}")
        seen = set()
        pending = deque(self.outgoing[node_id])
        while pending:
            descendant = pending.popleft()
            if descendant not in seen:
                seen.add(descendant)
                pending.extend(self.outgoing[descendant])
        return [descendant for descendant in self.order if descendant in seen]

    def source_impact(self, source_path: str) -> list[str]:
        """Find directly supported claims and descendants for a source change."""
        path = self._path(source_path)
        impacted = set()
        for node_id, node in self.nodes.items():
            if any(self._path(source["path"]) == path for source in node.get("sources", [])):
                impacted.add(node_id)
                impacted.update(self.impact(node_id))
        for rule in self.rules.values():
            if "proof" in rule and self._path(rule["proof"]["path"]) == path:
                impacted.add(rule["conclusion"])
                impacted.update(self.impact(rule["conclusion"]))
        return [node_id for node_id in self.order if node_id in impacted]

    def explain(self, node_id: str, evaluation: dict[str, dict] | None = None) -> dict:
        """Include recursive proof provenance and all incoming alternatives."""
        _require(node_id in self.nodes, f"unknown graph node: {node_id}")
        evaluated = evaluation if evaluation is not None else self.evaluate()

        def provenance(current: str) -> dict:
            record = deepcopy(evaluated[current])
            record["sources"] = deepcopy(self.nodes[current].get("sources", []))
            record["premise_provenance"] = []
            for support in record["support"]:
                if support["type"] == "implication":
                    rule = self.rules[support["id"]]
                    record["proof"] = deepcopy(rule["proof"])
                    record["premise_provenance"] = [provenance(premise) for premise in rule["premises"]]
            return record

        return {**provenance(node_id), "impact": self.impact(node_id)}

    def simulate(self, overrides: dict[str, str]) -> dict:
        """Compare a hypothetical scenario with current evaluation, read only."""
        baseline = self.evaluate()
        scenario = self.evaluate(overrides)
        changes = []
        for node_id in self.order:
            before, after = baseline[node_id], scenario[node_id]
            if any(before[key] != after[key] for key in ("state", "assumptions", "support", "hypothetical")):
                changes.append({"id": node_id, "before": before["state"], "after": after["state"],
                                "assumptions": after["assumptions"], "hypothetical": after["hypothetical"]})
        return {"hypothetical": True, "overrides": dict(overrides), "changes": changes, "nodes": scenario}

    def stale_evidence(self, evaluation: dict[str, dict] | None = None) -> list[dict]:
        """Report every stale declaration, even when an alternate route works."""
        evaluated = evaluation if evaluation is not None else self.evaluate()
        stale = []
        for node_id, record in evaluated.items():
            stale.extend({"owner": node_id, **source} for source in record["stale_sources"])
            for route in record["routes"]:
                stale.extend({"owner": route["id"], **source} for source in route["stale_sources"])
        return stale

    def to_dot(self, evaluation: dict[str, dict] | None = None) -> str:
        """Render implication diamonds so AND rules are not mistaken for OR."""
        import json

        evaluated = evaluation if evaluation is not None else self.evaluate()
        quote = json.dumps
        lines = ["digraph research {", "  rankdir=LR;"]
        for node_id in self.order:
            record = evaluated[node_id]
            label = f"{node_id}\n{record['state']}"
            lines.append(f"  {quote(node_id)} [label={quote(label)}, shape=box];")
        for rule_id, rule in self.rules.items():
            lines.append(f"  {quote(rule_id)} [label={quote(rule_id + ' (AND)')}, shape=diamond];")
            for premise in rule["premises"]:
                lines.append(f"  {quote(premise)} -> {quote(rule_id)};")
            lines.append(f"  {quote(rule_id)} -> {quote(rule['conclusion'])};")
        lines.append("}")
        return "\n".join(lines) + "\n"
