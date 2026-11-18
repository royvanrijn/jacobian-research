"""Inspect quantitative claims without silently changing their domains."""

from collections import defaultdict
import json
import math

SENSES = {"upper", "lower", "exact", "range", "criterion"}


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
    return rows


def bound_comparisons(graph, evaluation=None):
    """Compare only numeric assertions with identical declared domains/parameters.

    Conditional and unconditional claims are kept separate. Symbolic formulas,
    restricted families and different dimensions are never compared by guessing.
    """
    groups = defaultdict(list)
    for row in bounds_records(graph, evaluation):
        if row["domain_excluded"] or row["state"] not in {"proved", "conditional"}:
            continue
        if not isinstance(row.get("value"), (int, float)):
            continue
        if row["sense"] not in {"upper", "lower", "exact"}:
            continue
        key = json.dumps([row["quantity"], row["domain"], row.get("parameters", {}),
                          row.get("units"), row["state"], row["assumptions"]], sort_keys=True)
        groups[key].append(row)
    comparisons = []
    for rows in groups.values():
        upper = [row for row in rows if row["sense"] in {"upper", "exact"}]
        lower = [row for row in rows if row["sense"] in {"lower", "exact"}]
        best_upper = min((row["value"] for row in upper), default=None)
        best_lower = max((row["value"] for row in lower), default=None)
        opportunities = []
        if best_upper is not None:
            for row in upper:
                if row["value"] > best_upper:
                    opportunities.append({"weaker_claim": row["id"], "sense": "upper",
                                          "old_value": row["value"], "stronger_value": best_upper,
                                          "stronger_claims": sorted({item["id"] for item in upper if item["value"] == best_upper}),
                                          "affected": graph.impact(row["id"])})
        if best_lower is not None:
            for row in lower:
                if row["value"] < best_lower:
                    opportunities.append({"weaker_claim": row["id"], "sense": "lower",
                                          "old_value": row["value"], "stronger_value": best_lower,
                                          "stronger_claims": sorted({item["id"] for item in lower if item["value"] == best_lower}),
                                          "affected": graph.impact(row["id"])})
        comparisons.append({"quantity": rows[0]["quantity"], "domain": rows[0]["domain"],
                            "parameters": rows[0].get("parameters", {}), "units": rows[0].get("units"),
                            "state": rows[0]["state"],
                            "assumptions": rows[0]["assumptions"], "best_upper": best_upper,
                            "best_lower": best_lower,
                            "conflict": best_upper is not None and best_lower is not None and best_lower > best_upper,
                            "claims": sorted({row["id"] for row in rows}),
                            "opportunities": opportunities})
    return comparisons
