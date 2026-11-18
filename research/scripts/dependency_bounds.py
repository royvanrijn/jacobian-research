"""Inspect quantitative claims without silently changing their domains."""

from collections import defaultdict
import json
import math

SENSES = {"upper", "lower", "exact", "range", "criterion"}


def _scope_key(bound):
    return json.dumps([bound["quantity"], bound["domain"], bound.get("parameters", {}),
                       bound.get("units")], sort_keys=True)


def _strictly_improves(replacement, bound):
    sense = bound["sense"]
    if (sense not in {"upper", "lower"} or replacement["sense"] not in {sense, "exact"}
            or "value" not in bound or "value" not in replacement
            or _scope_key(replacement) != _scope_key(bound)):
        return False
    return (replacement["value"] < bound["value"] if sense == "upper"
            else replacement["value"] > bound["value"])


def bounds_records(graph, evaluation=None):
    evaluated = graph.evaluate() if evaluation is None else evaluation
    rows = []
    for node_id, node in graph.nodes.items():
        for bound in node.get("bounds", []):
            assert isinstance(bound, dict), "invalid bound: " + node_id
            assert bound.get("sense") in SENSES, "invalid bound sense: " + node_id
            assert isinstance(bound.get("quantity"), str) and bound["quantity"].strip(), "missing bound quantity: " + node_id
            assert isinstance(bound.get("domain"), str) and bound["domain"].strip(), "missing exact bound domain: " + node_id
            numeric = bound.get("value")
            if "value" in bound:
                assert (isinstance(numeric, (int, float)) and not isinstance(numeric, bool)
                        and (isinstance(numeric, int) or math.isfinite(numeric))), (
                            "bound value must be a finite number: " + node_id)
            else:
                assert isinstance(bound.get("expression"), str) and bound["expression"].strip(), (
                    "bound needs finite value or explicit expression: " + node_id)
            assert isinstance(bound.get("parameters", {}), dict), "invalid bound parameters: " + node_id
            exclusion = bound.get("domain_excluded_by")
            if exclusion is not None:
                assert isinstance(exclusion, str) and exclusion in graph.nodes and exclusion != node_id, (
                    "invalid bound domain exclusion: " + node_id)
            domain_excluded = exclusion is not None and evaluated[exclusion]["state"] == "proved"
            rows.append({**bound, "id": node_id, "title": node.get("title", node.get("label", node_id)),
                         "domain_excluded": domain_excluded, "state": evaluated[node_id]["state"],
                         "assumptions": evaluated[node_id]["assumptions"],
                         "dependency_coverage": evaluated[node_id]["dependency_coverage"]})
    by_claim = defaultdict(list)
    for row in rows:
        by_claim[row["id"]].append(row)
    for row in rows:
        replacement = row.get("superseded_by")
        if replacement is not None:
            assert isinstance(replacement, str) and replacement != row["id"], (
                "invalid bound replacement: " + row["id"])
            assert any(_strictly_improves(new, row) for new in by_claim.get(replacement, [])), (
                "bound replacement must strictly improve the same numeric assertion: " + row["id"])
    return rows


def bound_comparisons(graph, evaluation=None):
    """Compare only numeric assertions with identical declared domains/parameters.

    Conditional and unconditional claims are kept separate. Symbolic formulas,
    restricted families and different dimensions are never compared by guessing.
    A recorded replacement is current only if it attains the best available
    value in the same group.
    """
    groups = defaultdict(list)
    for row in bounds_records(graph, evaluation):
        if row["domain_excluded"] or row["state"] not in {"proved", "conditional"}:
            continue
        if not isinstance(row.get("value"), (int, float)):
            continue
        if row["sense"] not in {"upper", "lower", "exact"}:
            continue
        key = _scope_key(row), row["state"], tuple(row["assumptions"])
        groups[key].append(row)
    comparisons = []
    for rows in groups.values():
        upper = [row for row in rows if row["sense"] in {"upper", "exact"}]
        lower = [row for row in rows if row["sense"] in {"lower", "exact"}]
        best_upper = min((row["value"] for row in upper), default=None)
        best_lower = max((row["value"] for row in lower), default=None)
        opportunities, recorded_updates = [], []
        for sense, candidates, best in (("upper", upper, best_upper), ("lower", lower, best_lower)):
            stronger_claims = sorted({row["id"] for row in candidates if row["value"] == best})
            for row in candidates:
                if row["value"] == best:
                    continue
                replacement = row.get("superseded_by")
                update = {"weaker_claim": row["id"], "sense": sense, "old_value": row["value"],
                          "stronger_value": best}
                if replacement in stronger_claims:
                    recorded_updates.append({**update, "replacement": replacement})
                else:
                    opportunities.append({**update, "stronger_claims": stronger_claims,
                                          "affected": graph.impact(row["id"])})
        comparisons.append({"quantity": rows[0]["quantity"], "domain": rows[0]["domain"],
                            "parameters": rows[0].get("parameters", {}), "units": rows[0].get("units"),
                            "state": rows[0]["state"],
                            "assumptions": rows[0]["assumptions"], "best_upper": best_upper,
                            "best_lower": best_lower,
                            "conflict": best_upper is not None and best_lower is not None and best_lower > best_upper,
                            "claims": sorted({row["id"] for row in rows}),
                            "opportunities": opportunities, "recorded_updates": recorded_updates})
    return comparisons
