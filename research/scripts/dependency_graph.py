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
from pathlib import Path, PurePosixPath, PureWindowsPath
import re

from external_artifacts import repository_path, sha256_file


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
        self.equivalences: dict[str, dict] = {}
        self.component_of: dict[str, str] = {}
        self.components: dict[str, list[str]] = {}
        self.component_order: list[str] = []
        self.quotient_edges: dict[str, set[str]] = defaultdict(set)
        self.validate()

    def _path(self, value: str) -> Path:
        _require(_nonempty(value), "empty graph source path")
        _require(not PurePosixPath(value).is_absolute()
                 and not PureWindowsPath(value).drive,
                 f"absolute graph source path: {value}")
        try:
            return repository_path(self.root / value, self.repository)
        except ValueError as error:
            raise GraphValidationError(f"graph source escapes repository: {value}") from error

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
        """Reject malformed references and cycles outside proved-direction groups."""
        self.nodes.clear()
        self.rules.clear()
        self.incoming.clear()
        self.outgoing.clear()
        self.order.clear()
        self.equivalences.clear()
        self.component_of.clear()
        self.components.clear()
        self.component_order.clear()
        self.quotient_edges.clear()
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
        equivalences = self.graph.get("equivalences", [])
        _require(isinstance(equivalences, list), "graph.equivalences must be a list")
        directions = []
        for group in equivalences:
            _require(isinstance(group, dict), "invalid equivalence group")
            group_id = group.get("id")
            _require(_nonempty(group_id), "empty equivalence group id")
            _require(group_id not in self.equivalences and group_id not in self.nodes,
                     f"duplicate equivalence group id: {group_id}")
            members = group.get("members")
            _require(isinstance(members, list) and len(members) >= 2
                     and all(_nonempty(member) for member in members),
                     f"equivalence group needs at least two members: {group_id}")
            _require(len(members) == len(set(members)), f"duplicate equivalence member: {group_id}")
            for member in members:
                _require(member in self.nodes, f"dangling equivalence member {member}: {group_id}")
                _require(member not in self.component_of, f"overlapping equivalence groups: {member}")
                self.component_of[member] = group_id
            group_directions = group.get("directions")
            _require(isinstance(group_directions, list) and bool(group_directions),
                     f"missing equivalence directions: {group_id}")
            for direction in group_directions:
                _require(isinstance(direction, dict), f"invalid equivalence direction: {group_id}")
                _require(direction.get("conclusion") in members,
                         f"equivalence conclusion outside group: {group_id}")
                premises = direction.get("premises", [])
                _require(isinstance(premises, list)
                         and sum(premise in members for premise in premises) == 1,
                         f"equivalence direction needs exactly one member premise: {group_id}")
                directions.append({**direction, "equivalence": group_id})
            self.equivalences[group_id] = group
            self.components[group_id] = sorted(members)
        for node_id in self.nodes:
            if node_id not in self.component_of:
                self.component_of[node_id] = node_id
                self.components[node_id] = [node_id]
        for rule in implications:
            _require(isinstance(rule, dict) and "equivalence" not in rule,
                     "ordinary implications cannot claim an equivalence group")
        for rule in [*implications, *directions]:
            _require(isinstance(rule, dict), "invalid implication record")
            rule_id = rule.get("id")
            _require(_nonempty(rule_id), "empty implication id")
            _require(rule_id not in self.rules, f"duplicate implication id: {rule_id}")
            _require(rule_id not in self.nodes and rule_id not in self.equivalences,
                     f"implication id collides with node/group: {rule_id}")
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
        # Every internal edge must be an explicitly documented directional
        # equivalence, not an ordinary circular proof disguised by the group.
        for rule in self.rules.values():
            conclusion_component = self.component_of[rule["conclusion"]]
            for premise in rule["premises"]:
                premise_component = self.component_of[premise]
                if premise_component == conclusion_component:
                    _require(rule.get("equivalence") == conclusion_component,
                             f"internal dependency must be an equivalence direction: {rule['id']}")
                else:
                    self.quotient_edges[premise_component].add(conclusion_component)
        # Each group really has directional paths in both directions. External
        # context premises remain on the individual rules and are never waived.
        for group_id, group in self.equivalences.items():
            adjacency = defaultdict(set)
            for rule in self.rules.values():
                if rule.get("equivalence") == group_id:
                    member = next(p for p in rule["premises"] if p in group["members"])
                    adjacency[member].add(rule["conclusion"])
            for member in group["members"]:
                reached, pending = {member}, [member]
                while pending:
                    for target in adjacency[pending.pop()]:
                        if target not in reached:
                            reached.add(target)
                            pending.append(target)
                _require(set(group["members"]) <= reached,
                         f"equivalence directions are not strongly connected: {group_id}")
        # The quotient is a DAG, including unverified/refuted routes. A later
        # status change cannot activate an undeclared cross-component cycle.
        indegree = {component: 0 for component in self.components}
        for targets in self.quotient_edges.values():
            for target in targets:
                indegree[target] += 1
        pending = deque(sorted(component for component, degree in indegree.items() if degree == 0))
        while pending:
            component = pending.popleft()
            self.component_order.append(component)
            self.order.extend(self.components[component])
            for target in sorted(self.quotient_edges[component]):
                indegree[target] -= 1
                if indegree[target] == 0:
                    pending.append(target)
        _require(len(self.component_order) == len(self.components),
                 "dependency cycle outside equivalence groups: " + ", ".join(
                     sorted(component for component, degree in indegree.items() if degree)))

    def _freshness(self, sources: list[dict], cache: dict[Path, str | None]) -> list[dict]:
        stale = []
        for source in sources:
            path = self._path(source["path"])
            if path not in cache:
                try:
                    cache[path] = sha256_file(path)
                except OSError:
                    cache[path] = None
            actual = cache[path]
            if actual != source["sha256"]:
                stale.append({"path": source["path"], "expected_sha256": source["sha256"],
                              "actual_sha256": actual,
                              "reason": "missing" if actual is None else "changed"})
        return stale

    def evaluate(self, overrides: dict[str, str] | None = None, *,
                 distinct_blockers: bool = False) -> dict[str, dict]:
        """Evaluate the quotient DAG, grounding each equivalence component.

        Directional equivalences are consequences, never independent evidence.
        A component with no external proof seed remains open. Fresh evaluation
        starts from source leaves every time, so stale or withdrawn evidence
        cannot survive by circulating through an equivalence. Assumptions are
        unioned and retained through hypothetical and conditional directions.
        Navigation can keep only the first witness path for each blocker.
        """
        overrides = overrides or {}
        for node_id, state in overrides.items():
            _require(node_id in self.nodes, f"unknown simulation node: {node_id}")
            _require(state in {"proved", "open", "falsified", "stale"},
                     f"invalid simulation state: {node_id}={state}")
        evaluated, proof_keys, cache = {}, {}, {}
        locked, intrinsic_stale = set(), set()
        for node_id in self.order:
            node = self.nodes[node_id]
            sources = node.get("sources", []) if node["supplemental"] else []
            stale_sources = self._freshness(sources, cache)
            record = {"id": node_id, "kind": node["kind"],
                      "title": node.get("title", node.get("label", node_id)),
                      "scope": node.get("scope", ""), "recorded_state": node["recorded_state"],
                      "dependency_coverage": node.get("dependency_coverage", "unmapped"),
                      "component": self.component_of[node_id],
                      "state": "open" if node["supplemental"] else "unmapped",
                      "assumptions": [], "support": [], "routes": [],
                      "stale_sources": stale_sources, "blockers": [], "hypothetical": False,
                      "proof_depth": None, "proof_certificate": None}
            if node["supplemental"]:
                if stale_sources and node["status"] in {"proved", "falsified", "refuted"}:
                    record["state"] = "stale"
                    intrinsic_stale.add(node_id)
                elif node["status"] in {"falsified", "refuted"}:
                    record["state"] = "falsified"
                    locked.add(node_id)
                elif node["status"] == "proved" and not self.incoming[node_id]:
                    record["assumptions"] = sorted(node.get("assumptions", []))
                    record["state"] = "conditional" if record["assumptions"] else "proved"
                    record["support"] = [{"type": "leaf", "id": node_id, "sources": deepcopy(sources)}]
                    record["proof_depth"] = 0
            elif node["recorded_state"] == "falsified":
                record["state"] = "falsified"
                locked.add(node_id)
            if node_id in overrides:
                state = overrides[node_id]
                record.update(state="conditional" if state == "proved" else state,
                              assumptions=[node_id] if state == "proved" else [],
                              support=[{"type": "hypothesis", "id": node_id, "state": state}],
                              hypothetical=True, proof_depth=0 if state == "proved" else None)
                locked.add(node_id)
            if record["state"] in AVAILABLE:
                proof_keys[node_id] = (len(record["assumptions"]), tuple(record["assumptions"]), 0, "")
                record["proof_certificate"] = {
                    "id": node_id, "state": record["state"],
                    "assumptions": list(record["assumptions"]),
                    "support": [{key: value for key, value in support.items() if key != "sources"}
                                for support in record["support"]],
                    "hypothetical": record["hypothetical"], "proof_depth": 0, "premises": [],
                }
            evaluated[node_id] = record
        rule_stale = {rule_id: self._freshness([rule["proof"]], cache) if "proof" in rule else []
                      for rule_id, rule in self.rules.items()}

        def route_for(rule_id: str) -> dict:
            rule = self.rules[rule_id]
            record = evaluated[rule["conclusion"]]
            unavailable = [p for p in rule["premises"] if evaluated[p]["state"] not in AVAILABLE]
            assumptions = set(self.nodes[rule["conclusion"]].get("assumptions", []))
            assumptions.update(rule.get("assumptions", []))
            for premise in rule["premises"]:
                assumptions.update(evaluated[premise]["assumptions"])
            if rule["status"] == "refuted":
                state = "refuted"
            elif rule["status"] != "proved":
                state = "unverified"
            elif rule_stale[rule_id]:
                state = "stale"
            elif unavailable:
                state = "stale" if any(evaluated[p]["state"] == "stale" for p in unavailable) else "blocked"
            else:
                state = "conditional" if assumptions else "proved"
            return {"id": rule_id, "state": state, "premises": list(rule["premises"]),
                    "assumptions": sorted(assumptions), "unavailable_premises": unavailable,
                    "stale_sources": deepcopy(rule_stale[rule_id]), "proof": deepcopy(rule.get("proof")),
                    "equivalence": rule.get("equivalence"),
                    "hypothetical": any(evaluated[p]["hypothetical"] for p in rule["premises"])}

        for component in self.component_order:
            members = self.components[component]
            changed = True
            while changed:
                changed = False
                for node_id in members:
                    if node_id in locked:
                        continue
                    record = evaluated[node_id]
                    for rule_id in self.incoming[node_id]:
                        route = route_for(rule_id)
                        if route["state"] not in AVAILABLE:
                            continue
                        depth = 1 + max(evaluated[p]["proof_depth"] for p in route["premises"])
                        key = (len(route["assumptions"]), tuple(route["assumptions"]), depth, rule_id)
                        if node_id in proof_keys and key >= proof_keys[node_id]:
                            continue
                        # Along the selected proof, assumptions cannot decrease
                        # towards a conclusion and depth strictly increases when
                        # assumptions match. Monotone key improvement therefore
                        # prevents a selected A<-B<-A provenance cycle.
                        proof_keys[node_id] = key
                        record.update(state=route["state"], assumptions=route["assumptions"],
                                      support=[{"type": "implication", "id": rule_id}],
                                      hypothetical=route["hypothetical"], proof_depth=depth)
                        # Freeze the exact premise proofs used here. A premise
                        # may later prefer a different conditional route; that
                        # must not rewrite this already grounded certificate or
                        # silently change its required assumption union.
                        record["proof_certificate"] = {
                            "id": node_id, "state": route["state"],
                            "assumptions": list(route["assumptions"]),
                            "support": [{"type": "implication", "id": rule_id}],
                            "hypothetical": route["hypothetical"], "proof_depth": depth,
                            "premises": [evaluated[p]["proof_certificate"] for p in route["premises"]],
                        }
                        changed = True
            # Propagate loss of current evidence through ungrounded routes too.
            changed = True
            while changed:
                changed = False
                for node_id in members:
                    record = evaluated[node_id]
                    if node_id in locked or record["state"] in AVAILABLE:
                        continue
                    routes = [route_for(rule_id) for rule_id in self.incoming[node_id]]
                    state = ("stale" if node_id in intrinsic_stale or any(
                        route["state"] == "stale" for route in routes)
                        else "open" if routes or record["kind"] != "result" else "unmapped")
                    hypothetical = any(route["hypothetical"] for route in routes)
                    if (record["state"], record["hypothetical"]) != (state, hypothetical):
                        record.update(state=state, hypothetical=hypothetical)
                        changed = True
        for node_id in self.order:
            record = evaluated[node_id]
            record["routes"] = [route_for(rule_id) for rule_id in self.incoming[node_id]]
            certificate = record["proof_certificate"]
            if certificate and certificate["support"][0]["type"] == "implication":
                selected_id = certificate["support"][0]["id"]
                for route in record["routes"]:
                    if route["id"] != selected_id:
                        continue
                    if route["assumptions"] != certificate["assumptions"]:
                        route["preferred_premise_assumptions"] = list(route["assumptions"])
                    route.update(state=certificate["state"],
                                 assumptions=list(certificate["assumptions"]),
                                 hypothetical=certificate["hypothetical"], selected=True,
                                 selected_proof_depth=certificate["proof_depth"])
        for node_id in self.order:
            evaluated[node_id]["blockers"] = self.blockers(
                node_id, evaluated, distinct=distinct_blockers)
        return evaluated

    def blockers(self, node_id: str, evaluation: dict[str, dict] | None = None, *,
                 distinct: bool = False) -> list[dict]:
        """Return unresolved leaves/rules for each alternative route, with paths.

        This is a route explanation, not a claim that all listed blockers are
        jointly necessary. One repaired alternative may suffice.
        A falsified premise terminates its route. A falsified target has no
        unresolved proof obligations. In distinct mode, keep the same first
        witness that would appear for each blocker identity in the full list.
        """
        _require(node_id in self.nodes, f"unknown graph node: {node_id}")
        if evaluation is None:
            evaluation = self.evaluate(distinct_blockers=distinct)
            return evaluation[node_id]["blockers"]
        if evaluation[node_id]["state"] == "falsified":
            return []
        found = []
        seen_blockers = set()
        seen_visits = set()

        def add(blocker: dict) -> None:
            if distinct:
                identity = (blocker["type"], blocker["id"], blocker["state"])
                if identity in seen_blockers:
                    return
                seen_blockers.add(identity)
            found.append(blocker)

        def visit(current: str, path: list[str], ancestors: frozenset[str] = frozenset()) -> None:
            if current in ancestors:
                add({"type": "equivalence", "id": self.component_of[current],
                     "state": "open", "path": path,
                     "reason": "equivalence has no independent available grounding"})
                return
            if distinct:
                visit_key = (current, ancestors)
                if visit_key in seen_visits:
                    return
                seen_visits.add(visit_key)
            ancestors = ancestors | {current}
            record = evaluation[current]
            if record["state"] in AVAILABLE:
                for assumption in record["assumptions"]:
                    add({"type": "assumption", "id": assumption, "state": "open", "path": path})
                return
            routes = record["routes"]
            if record["state"] == "falsified" or not routes:
                add({"type": "node", "id": current, "state": record["state"], "path": path,
                     "stale_sources": record["stale_sources"]})
                return
            for route in routes:
                route_path = [*path, route["id"]]
                if route["state"] in {"unverified", "refuted"} or route["stale_sources"]:
                    add({"type": "implication", "id": route["id"], "state": route["state"],
                         "path": route_path, "stale_sources": route["stale_sources"]})
                for premise in route["unavailable_premises"]:
                    # An external premise belongs to an earlier quotient component,
                    # so it cannot lead back to ancestors in this component.
                    if distinct and self.component_of[premise] != self.component_of[current]:
                        visit(premise, [*route_path, premise])
                    else:
                        visit(premise, [*route_path, premise], ancestors)
                for assumption in route["assumptions"]:
                    add({"type": "assumption", "id": assumption, "state": "open", "path": route_path})

        visit(node_id, [node_id])
        return found

    def impact(self, node_id: str) -> list[str]:
        """List all structural descendants, including currently blocked routes."""
        _require(node_id in self.nodes, f"unknown graph node: {node_id}")
        seen = {node_id}
        pending = deque(self.outgoing[node_id])
        while pending:
            descendant = pending.popleft()
            if descendant not in seen:
                seen.add(descendant)
                pending.extend(self.outgoing[descendant])
        return [descendant for descendant in self.order if descendant in seen and descendant != node_id]

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

        def provenance(current: str, certificate: dict | None = None,
                       ancestors: frozenset[str] = frozenset()) -> dict:
            _require(current not in ancestors, f"circular selected proof provenance: {current}")
            ancestors = ancestors | {current}
            preferred = evaluated[current]
            # The expanded premise tree is itself the readable certificate; do
            # not duplicate the nested compact certificate at every tree node.
            record = deepcopy({key: value for key, value in preferred.items()
                               if key != "proof_certificate"})
            certificate = certificate or preferred.get("proof_certificate")
            if certificate:
                for key in ("state", "assumptions", "support", "hypothetical", "proof_depth"):
                    record[key] = deepcopy(certificate[key])
                record["provenance_basis"] = "selected grounded proof certificate"
                preferred_support = [{key: value for key, value in support.items() if key != "sources"}
                                     for support in preferred["support"]]
                if record["assumptions"] != preferred["assumptions"] or record["support"] != preferred_support:
                    record["preferred_assumptions"] = list(preferred["assumptions"])
                    record["preferred_support"] = deepcopy(preferred["support"])
            record["sources"] = deepcopy(self.nodes[current].get("sources", []))
            record["premise_provenance"] = []
            for support in record["support"]:
                if support["type"] == "leaf":
                    support["sources"] = deepcopy(record["sources"])
                elif support["type"] == "implication":
                    rule = self.rules[support["id"]]
                    record["proof"] = deepcopy(rule["proof"])
                    _require(certificate is not None, f"missing selected proof certificate: {current}")
                    premises = certificate["premises"]
                    _require([premise["id"] for premise in premises] == rule["premises"],
                             f"selected certificate premises differ from rule: {support['id']}")
                    record["premise_provenance"] = [
                        provenance(premise["id"], premise, ancestors) for premise in premises]
                    required = set(self.nodes[current].get("assumptions", []))
                    required.update(rule.get("assumptions", []))
                    for premise in premises:
                        required.update(premise["assumptions"])
                    _require(set(record["assumptions"]) == required,
                             f"selected proof certificate drops or adds assumptions: {current}")
                    _require(record["proof_depth"] == 1 + max(
                        premise["proof_depth"] for premise in premises),
                        f"invalid selected proof certificate depth: {current}")
                    for route in record["routes"]:
                        route.pop("selected", None)
                        if route["id"] == support["id"]:
                            if route["assumptions"] != record["assumptions"]:
                                route["preferred_premise_assumptions"] = list(route["assumptions"])
                            route.update(state=record["state"],
                                         assumptions=list(record["assumptions"]),
                                         hypothetical=record["hypothetical"], selected=True,
                                         selected_proof_depth=record["proof_depth"])
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

    def quotient(self) -> dict:
        """Return the structural acyclic quotient, retaining member identities."""
        return {"components": [{"id": component, "members": list(self.components[component]),
                                "equivalence": component in self.equivalences}
                               for component in self.component_order],
                "edges": [{"premise": premise, "conclusion": conclusion}
                          for premise in self.component_order
                          for conclusion in sorted(self.quotient_edges[premise])]}

    def to_dot(self, evaluation: dict[str, dict] | None = None) -> str:
        """Render implication diamonds so AND rules are not mistaken for OR."""
        import json

        evaluated = evaluation if evaluation is not None else self.evaluate()
        quote = json.dumps
        lines = ["digraph research {", "  rankdir=LR;"]
        for index, (group_id, group) in enumerate(self.equivalences.items()):
            lines += [f"  subgraph cluster_equivalence_{index} {{",
                      f"    label={quote(group.get('title', group_id) + ' (independent grounding required)')};"]
            lines += [f"    {quote(member)};" for member in group["members"]]
            lines.append("  }")
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
