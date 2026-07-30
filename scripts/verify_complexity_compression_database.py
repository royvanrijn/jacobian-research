#!/usr/bin/env python3
"""Validate and summarize the Programme 8 complexity database.

The checker is deliberately dependency-free.  It validates the typed census,
recomputes elementary complexity data from stored sparse F=X+H artifacts, and
checks that every declared monotone transformation is actually monotone in
its selected objectives and equal on its machine-checkable preserved fields.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = ROOT / "extended-geometry" / "complexity_compression_database.json"
DEFAULT_REPORT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "programme8_complexity_compression_report.json"
)

VALUE_STATUSES = {
    "exact",
    "theorem",
    "formula",
    "upper_bound",
    "diagnostic",
    "uncomputed",
    "not_defined",
    "not_applicable",
    "not_established",
}
FORMAL_STATUSES = {"complete", "partial", "none"}


def fail(message: str) -> None:
    raise AssertionError(message)


def load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def resolve_path(value: Any, dotted_path: str) -> Any:
    current = value
    for part in dotted_path.split("."):
        if isinstance(current, list):
            try:
                current = current[int(part)]
            except (ValueError, IndexError) as exc:
                fail(f"invalid list path {dotted_path!r} at {part!r}: {exc}")
        elif isinstance(current, dict):
            if part not in current:
                fail(f"missing path {dotted_path!r} at {part!r}")
            current = current[part]
        else:
            fail(f"path {dotted_path!r} descends through scalar {current!r}")
    return current


def validate_metric(metric: dict[str, Any], label: str) -> None:
    if not isinstance(metric, dict):
        fail(f"{label} must be an object")
    if "value" not in metric:
        fail(f"{label} has no value")
    if metric.get("status") not in VALUE_STATUSES:
        fail(f"{label} has invalid status {metric.get('status')!r}")
    if not metric.get("source"):
        fail(f"{label} has no source/explanation")
    missing_statuses = {
        "uncomputed",
        "not_defined",
        "not_applicable",
        "not_established",
    }
    if metric["value"] is None and metric["status"] not in missing_statuses:
        fail(f"{label} has null value with affirmative status {metric['status']!r}")
    if metric["value"] is not None and metric["status"] in missing_statuses:
        fail(f"{label} has non-null value with missing status {metric['status']!r}")


def validate_construction(
    construction: dict[str, Any], required_measures: list[str]
) -> None:
    cid = construction.get("id")
    if not cid:
        fail("construction without id")
    measures = construction.get("measures")
    if not isinstance(measures, dict):
        fail(f"{cid}: measures must be an object")
    missing = sorted(set(required_measures) - set(measures))
    extra = sorted(set(measures) - set(required_measures))
    if missing or extra:
        fail(f"{cid}: measure mismatch; missing={missing}, extra={extra}")

    for name in (
        "source_dimension",
        "coordinate_degree_vector",
        "coefficient_height",
        "geometric_degree",
        "boundary_prime_count",
        "puncture_count",
        "stable_moduli_dimension",
    ):
        validate_metric(measures[name], f"{cid}.{name}")

    support = measures["support_vector"]
    validate_metric(support, f"{cid}.support_vector")
    if "total" not in support:
        fail(f"{cid}.support_vector has no total")
    if isinstance(support["value"], list):
        expected_total = sum(support["value"])
        if support["total"] != expected_total:
            fail(
                f"{cid}.support_vector total {support['total']} "
                f"does not equal {expected_total}"
            )
    elif support["value"] is None and support["total"] is not None:
        fail(f"{cid}.support_vector has null vector but non-null total")

    monodromy = measures["monodromy"]
    validate_metric(monodromy, f"{cid}.monodromy")
    if "action" not in monodromy:
        fail(f"{cid}.monodromy has no action field")

    operator_rows = measures["rank_and_nilpotency_indices"]
    if not isinstance(operator_rows, list) or not operator_rows:
        fail(f"{cid}.rank_and_nilpotency_indices must be a nonempty list")
    for index, row in enumerate(operator_rows):
        label = f"{cid}.rank_and_nilpotency_indices[{index}]"
        if not row.get("operator"):
            fail(f"{label} has no operator")
        validate_metric(row.get("rank"), f"{label}.rank")
        validate_metric(row.get("nilpotency_index"), f"{label}.nilpotency_index")

    formal = measures["formal_verification_status"]
    if formal.get("status") not in FORMAL_STATUSES:
        fail(f"{cid}: invalid formal status {formal.get('status')!r}")
    if not formal.get("scope") or not formal.get("source"):
        fail(f"{cid}: formal status requires scope and source")


def term_monomial(term: dict[str, Any]) -> tuple[tuple[int, int], ...]:
    return tuple((int(variable), int(exponent)) for variable, exponent in term["monomial"])


def derive_sparse_map_measures(artifact: dict[str, Any]) -> dict[str, Any]:
    dimension = int(artifact["dimension"])
    components = artifact["H"]
    if len(components) != dimension:
        fail(
            f"sparse artifact dimension {dimension} does not match "
            f"{len(components)} H components"
        )

    degrees: list[int] = []
    supports: list[int] = []
    heights: list[int] = []
    correction_terms = 0

    for coordinate, terms in enumerate(components):
        polynomial: dict[tuple[tuple[int, int], ...], Fraction] = {
            ((coordinate, 1),): Fraction(1)
        }
        correction_terms += len(terms)
        for term in terms:
            monomial = term_monomial(term)
            polynomial[monomial] = polynomial.get(monomial, Fraction(0)) + Fraction(
                term["coefficient"]
            )
        polynomial = {
            monomial: coefficient
            for monomial, coefficient in polynomial.items()
            if coefficient
        }
        degrees.append(
            max(sum(exponent for _, exponent in monomial) for monomial in polynomial)
        )
        supports.append(len(polynomial))
        heights.extend(
            max(abs(coefficient.numerator), coefficient.denominator)
            for coefficient in polynomial.values()
        )

    return {
        "source_dimension": dimension,
        "coordinate_degree_vector": degrees,
        "support_vector": supports,
        "support_total": sum(supports),
        "coefficient_height": max(heights),
        "correction_component_monomial_terms": correction_terms,
    }


def validate_machine_source(construction: dict[str, Any]) -> None:
    if construction.get("class") != "cubic_homogeneous_keller_map":
        return
    relative = construction.get("machine_source")
    if not relative:
        fail(f"{construction['id']}: cubic map has no machine_source")
    path = ROOT / relative
    artifact = load_json(path)
    derived = derive_sparse_map_measures(artifact)
    measures = construction["measures"]
    expected = {
        "source_dimension": measures["source_dimension"]["value"],
        "coordinate_degree_vector": measures["coordinate_degree_vector"]["value"],
        "support_vector": measures["support_vector"]["value"],
        "support_total": measures["support_vector"]["total"],
        "coefficient_height": measures["coefficient_height"]["value"],
        "correction_component_monomial_terms": construction["auxiliary_measures"][
            "correction_component_monomial_terms"
        ],
    }
    if derived != expected:
        fail(
            f"{construction['id']}: machine-source mismatch\n"
            f"derived={derived}\nexpected={expected}"
        )

    statistics = artifact.get("statistics", {})
    rank = statistics.get("generic_rank_JH_over_QQ_x")
    index = statistics.get("nilpotency_index_JH")
    if rank is not None:
        database_rank = measures["rank_and_nilpotency_indices"][0]["rank"]["value"]
        if rank != database_rank:
            fail(f"{construction['id']}: rank {database_rank} != artifact {rank}")
    if index is not None:
        database_index = measures["rank_and_nilpotency_indices"][0][
            "nilpotency_index"
        ]["value"]
        if index != database_index:
            fail(f"{construction['id']}: index {database_index} != artifact {index}")


def validate_quartic_scaling_rows(
    constructions: dict[str, dict[str, Any]]
) -> None:
    artifact = load_json(
        ROOT / "artifacts/generated-results/quartic_lr_sparsity_search.json"
    )["rational_scalings"]
    balanced = constructions["quadratic-gauge-quartic-sparse-balanced"]
    coefficient = constructions[
        "quadratic-gauge-quartic-coefficient-optimized-grid"
    ]
    expected_balanced = {
        "coefficient_height": balanced["measures"]["coefficient_height"]["value"],
        "collision_height": 19856,
    }
    actual_balanced = {
        "coefficient_height": artifact["best_balanced"]["coefficient_height"],
        "collision_height": artifact["best_balanced"]["collision_height"],
    }
    if actual_balanced != expected_balanced:
        fail(
            "balanced quartic row disagrees with quartic scaling artifact: "
            f"{actual_balanced} != {expected_balanced}"
        )

    expected_coefficient = {
        "alpha": coefficient["parameters"]["alpha"],
        "beta": coefficient["parameters"]["beta"],
        "coefficient_height": coefficient["measures"]["coefficient_height"]["value"],
        "collision_height": coefficient["auxiliary_measures"][
            "collision_coordinate_height"
        ],
    }
    actual_coefficient = {
        key: artifact["best_coefficient"][key] for key in expected_coefficient
    }
    if actual_coefficient != expected_coefficient:
        fail(
            "coefficient-optimized quartic row disagrees with scaling artifact: "
            f"{actual_coefficient} != {expected_coefficient}"
        )
    if (
        artifact["candidate_count"]
        != coefficient["auxiliary_measures"]["grid_candidate_count"]
    ):
        fail("quartic scaling candidate count mismatch")


def validate_status_references(database: dict[str, Any]) -> None:
    status_entries = {entry["id"] for entry in load_json(ROOT / "MATH_STATUS.json")["entries"]}
    for construction in database["constructions"]:
        for source in construction.get("sources", []):
            prefix = "MATH_STATUS.json#"
            if source.startswith(prefix):
                status_id = source[len(prefix) :]
                if status_id not in status_entries:
                    fail(
                        f"{construction['id']}: unknown MATH_STATUS id {status_id!r}"
                    )


def validate_transformations(
    database: dict[str, Any], constructions: dict[str, dict[str, Any]]
) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    seen: set[str] = set()
    for transformation in database["transformations"]:
        tid = transformation.get("id")
        if not tid or tid in seen:
            fail(f"duplicate or missing transformation id {tid!r}")
        seen.add(tid)
        before = constructions.get(transformation.get("from"))
        after = constructions.get(transformation.get("to"))
        if before is None or after is None:
            fail(f"{tid}: unknown from/to construction")

        for path in transformation.get("preserved_measure_paths", []):
            before_value = resolve_path(before["measures"], path)
            after_value = resolve_path(after["measures"], path)
            if before_value != after_value:
                fail(
                    f"{tid}: preserved path {path} changes "
                    f"{before_value!r} -> {after_value!r}"
                )

        objectives = transformation.get("objectives", [])
        if not objectives:
            fail(f"{tid}: no selected objective")
        strict = False
        objective_summaries = []
        for objective in objectives:
            path = objective["measure_path"]
            before_value = resolve_path(before["measures"], path)
            after_value = resolve_path(after["measures"], path)
            if before_value != objective["before"] or after_value != objective["after"]:
                fail(
                    f"{tid}: declared objective {path} does not match construction rows"
                )
            if not isinstance(before_value, (int, float)) or not isinstance(
                after_value, (int, float)
            ):
                fail(f"{tid}: objective {path} is not numeric")
            direction = objective["direction"]
            if direction == "minimize":
                if after_value > before_value:
                    fail(f"{tid}: minimizing {path} worsens")
                strict |= after_value < before_value
            elif direction == "maximize":
                if after_value < before_value:
                    fail(f"{tid}: maximizing {path} worsens")
                strict |= after_value > before_value
            else:
                fail(f"{tid}: invalid direction {direction!r}")
            objective_summaries.append(
                {
                    "measure_path": path,
                    "direction": direction,
                    "before": before_value,
                    "after": after_value,
                }
            )
        if not strict:
            fail(f"{tid}: selected objectives contain no strict improvement")
        summaries.append(
            {
                "id": tid,
                "kind": transformation["kind"],
                "status": transformation["status"],
                "from": transformation["from"],
                "to": transformation["to"],
                "objectives": objective_summaries,
                "preserved_measure_paths": transformation.get(
                    "preserved_measure_paths", []
                ),
            }
        )
    return summaries


def numeric_complexity_vector(construction: dict[str, Any]) -> dict[str, int | float]:
    measures = construction["measures"]
    result: dict[str, int | float] = {}
    candidates = {
        "source_dimension": measures["source_dimension"]["value"],
        "support_total": measures["support_vector"]["total"],
        "coefficient_height": measures["coefficient_height"]["value"],
        "rank": measures["rank_and_nilpotency_indices"][0]["rank"]["value"],
        "nilpotency_index": measures["rank_and_nilpotency_indices"][0][
            "nilpotency_index"
        ]["value"],
    }
    for key, value in candidates.items():
        if isinstance(value, (int, float)):
            result[key] = value
    degree_vector = measures["coordinate_degree_vector"]["value"]
    if isinstance(degree_vector, list) and all(
        isinstance(value, (int, float)) for value in degree_vector
    ):
        result["maximum_coordinate_degree"] = max(degree_vector)
    return result


def comparison_frontiers(
    database: dict[str, Any],
) -> list[dict[str, Any]]:
    classes: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for construction in database["constructions"]:
        if comparison_class := construction.get("comparison_class"):
            classes[comparison_class].append(construction)

    reports = []
    for comparison_class, members in sorted(classes.items()):
        vectors = {member["id"]: numeric_complexity_vector(member) for member in members}
        common_metrics = sorted(
            set.intersection(*(set(vector) for vector in vectors.values()))
        )
        nondominated = []
        domination_edges = []
        for candidate_id, candidate in vectors.items():
            dominated = False
            for other_id, other in vectors.items():
                if candidate_id == other_id:
                    continue
                weak = all(other[key] <= candidate[key] for key in common_metrics)
                strict = any(other[key] < candidate[key] for key in common_metrics)
                if weak and strict:
                    dominated = True
                    domination_edges.append(
                        {"dominator": other_id, "dominated": candidate_id}
                    )
            if not dominated:
                nondominated.append(candidate_id)
        reports.append(
            {
                "comparison_class": comparison_class,
                "common_numeric_complexity_metrics": common_metrics,
                "vectors": vectors,
                "nondominated": sorted(nondominated),
                "domination_edges": sorted(
                    domination_edges,
                    key=lambda edge: (edge["dominator"], edge["dominated"]),
                ),
            }
        )
    return reports


def coverage_report(database: dict[str, Any]) -> dict[str, Any]:
    coverage: dict[str, Counter[str]] = {
        measure: Counter() for measure in database["required_measures"]
    }
    for construction in database["constructions"]:
        measures = construction["measures"]
        for measure in database["required_measures"]:
            if measure == "rank_and_nilpotency_indices":
                rows = measures[measure]
                rank_statuses = {row["rank"]["status"] for row in rows}
                index_statuses = {row["nilpotency_index"]["status"] for row in rows}
                status = (
                    "complete_operator_data"
                    if rank_statuses <= {"exact", "theorem", "formula"}
                    and index_statuses <= {"exact", "theorem", "formula"}
                    else "contains_gap_or_not_applicable"
                )
            elif measure == "formal_verification_status":
                status = measures[measure]["status"]
            else:
                status = measures[measure]["status"]
            coverage[measure][status] += 1
    return {
        measure: dict(sorted(counter.items()))
        for measure, counter in coverage.items()
    }


def make_report(
    database_path: Path,
    database: dict[str, Any],
    transformation_summaries: list[dict[str, Any]],
) -> dict[str, Any]:
    raw = database_path.read_bytes()
    return {
        "format": "programme8-complexity-compression-report-v1",
        "database": str(database_path.relative_to(ROOT)),
        "database_sha256": hashlib.sha256(raw).hexdigest(),
        "construction_count": len(database["constructions"]),
        "transformation_count": len(database["transformations"]),
        "search_hypothesis_count": len(database["search_hypotheses"]),
        "coverage_by_measure_and_status": coverage_report(database),
        "certified_monotone_relations": transformation_summaries,
        "comparison_frontiers": comparison_frontiers(database),
        "guardrails": [
            "Null values remain explicit gaps and are never copied from a nearby construction.",
            "Diagnostic sampled nilpotency indices are excluded from exact objectives.",
            "Construction replacements are not relabeled as equivalences.",
            "Pareto comparison uses only numeric complexity metrics common to every member of the declared comparison class.",
        ],
    }


def write_json_atomic(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rendered = json.dumps(value, indent=2, sort_keys=True) + "\n"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent, text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(rendered)
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument(
        "--write-report",
        nargs="?",
        const=DEFAULT_REPORT,
        type=Path,
        help="Write the deterministic generated report (default path if omitted).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    database_path = args.database.resolve()
    database = load_json(database_path)
    if database.get("schema_version") != 1:
        fail(f"unsupported schema version {database.get('schema_version')!r}")
    required_measures = database.get("required_measures")
    if not isinstance(required_measures, list) or not required_measures:
        fail("required_measures must be a nonempty list")

    constructions: dict[str, dict[str, Any]] = {}
    for construction in database.get("constructions", []):
        validate_construction(construction, required_measures)
        cid = construction["id"]
        if cid in constructions:
            fail(f"duplicate construction id {cid!r}")
        constructions[cid] = construction
        validate_machine_source(construction)

    validate_quartic_scaling_rows(constructions)
    validate_status_references(database)
    transformation_summaries = validate_transformations(database, constructions)
    report = make_report(database_path, database, transformation_summaries)
    if args.write_report:
        write_json_atomic(args.write_report.resolve(), report)

    print(
        "programme8 database: "
        f"{len(constructions)} constructions, "
        f"{len(transformation_summaries)} monotone relations, "
        f"{len(database['search_hypotheses'])} search hypotheses"
    )
    print("machine-source replays: exact")
    print("transformation monotonicity: exact")
    if args.write_report:
        print(f"wrote {args.write_report.resolve().relative_to(ROOT)}")


if __name__ == "__main__":
    main()
