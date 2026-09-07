"""Recover a censoring-aware Fermigier specialization corpus.

The historical search files use two parameter coordinates.  In Fermigier's
published/adapter coordinate ``u`` the tuple-model coordinate is ``T=2*u``.
The tuple model depends on ``T**2``, so this module identifies ``T`` and
``-T`` before joining any search records.

Only the two explicitly supplied exact certificates become positive labels.
Every other row is censored: a bounded point search, a numerical height rank,
or a legacy ``rank>=`` log line is retained as an observation, never promoted
to a mathematical rank statement.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from hashlib import sha256
import csv
import gzip
import io
import json
from pathlib import Path
import re
from typing import Any, Iterable, Iterator


SCHEMA = "elliptic-curves.fermigier-labelled-corpus-row.v1"
FAMILY_ID = "fermigier-mestre-v1"


@dataclass(frozen=True)
class TabularSource:
    name: str
    relative_path: str
    kind: str


@dataclass(frozen=True)
class JsonSource:
    name: str
    relative_path: str
    collections: tuple[str, ...]
    required: bool = True
    lane_filter: str | None = None


@dataclass(frozen=True)
class LegacyLogSource:
    name: str
    relative_path: str
    default_height: int


TABULAR_SOURCES = (
    TabularSource(
        "hot_neighborhood_population",
        "artifacts/local/elliptic-curves/fermigier_hot_neighborhood.tsv",
        "pair_population",
    ),
    TabularSource(
        "global_score_top5000",
        "artifacts/local/elliptic-curves/fermigier_global_score_top5000.tsv",
        "rank_score_no_header",
    ),
    TabularSource(
        "hot_neighborhood_score_top5000",
        "artifacts/local/elliptic-curves/fermigier_hot_neighborhood_top5000.tsv",
        "rank_score_no_header",
    ),
    TabularSource(
        "multibound_score_top5000",
        "artifacts/local/elliptic-curves/fermigier_multibound_top5000.tsv",
        "dict_table",
    ),
    TabularSource(
        "broad_ensemble",
        "artifacts/local/elliptic-curves/fermigier_broad_ensemble.tsv",
        "dict_table",
    ),
    TabularSource(
        "family_residual_score_h100000_b200_top20000",
        "artifacts/local/elliptic-curves/fermigier_family_residual_score_h100000_b200_top20000.tsv",
        "comment_header_dict_table",
    ),
)

HOT_NEIGHBORHOOD_PROVENANCE = (
    "artifacts/local/elliptic-curves/fermigier_hot_neighborhood_provenance.tsv"
)

SCORE_METADATA_SOURCES = (
    "artifacts/local/elliptic-curves/fermigier_global_score_top5000.log",
    "artifacts/local/elliptic-curves/fermigier_multibound.log",
    "artifacts/local/elliptic-curves/fermigier_hot_neighborhood_score.log",
    "artifacts/local/elliptic-curves/fermigier_broad_ensemble.log",
    "artifacts/local/elliptic-curves/fermigier_family_residual_score_h100000_b200_top20000.log",
    "elliptic-curves/ecsearch/fermigier_score_sweep.cpp",
)

CORPUS_IMPLEMENTATION_SOURCES = (
    "elliptic-curves/ecsearch/fermigier_labelled_corpus.py",
    "elliptic-curves/scripts/build_fermigier_labelled_corpus.py",
)


JSON_SOURCES = (
    JsonSource(
        "global_scan",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_global.json",
        (
            "calibrations.*",
            "selection.conductor_tranche.*",
            "conductor_replays.*",
            "point_triage.stages.*.results.*",
        ),
    ),
    JsonSource(
        "batch_rank_triage",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_batch_rank_triage.json",
        ("results.*",),
    ),
    JsonSource(
        "crt_lattice_pilot",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_crt_lattice_pilot.json",
        ("candidates.*",),
    ),
    JsonSource(
        "discovered_local_conditions",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_discovered_local_conditions.json",
        ("candidates.*",),
    ),
    JsonSource(
        "multiple_root_crt",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_multiple_root_crt.json",
        ("candidates.*",),
    ),
    JsonSource(
        "multiple_root_frontier",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_multiple_root_frontier.json",
        ("candidates.*",),
    ),
    JsonSource(
        "multiple_root_height_h50000",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_multiple_root_height_h50000.json",
        ("finalists.*",),
    ),
    JsonSource(
        "record_residue_deep_tranche",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_record_residue_deep_tranche.json",
        ("candidates.*",),
    ),
    JsonSource(
        "record_residue_class",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_record_residue_class.json",
        ("candidates.*",),
    ),
    JsonSource(
        "score_cutoffs",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_score_cutoffs.json",
        ("curves.*",),
    ),
    JsonSource(
        "extra_points",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_extra_points.json",
        ("results.*",),
    ),
    JsonSource(
        "checkpoint_3115_3_h1000000",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_3115_3_h1000000_checkpoint.json",
        ("results.*",),
    ),
    JsonSource(
        "high_power_crt_gauss",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_high_power_crt_gauss.json",
        ("selected.*",),
    ),
    JsonSource(
        "record_rescore_h5000",
        "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_record_rescore_h5000.json",
        ("finalists.*",),
    ),
    JsonSource(
        "rank22_accidental_slices",
        "artifacts/generated-results/elliptic-curves/elliptic_fermigier_rank22_accidental_slices.json",
        ("candidate_conductor_screen.*", "multi_slice_intersections.*"),
    ),
    JsonSource(
        "conductor_first_anchor_pilot",
        "artifacts/generated-results/elliptic-curves/conductor_first_family_anchor_pilot_v1.json",
        ("candidates.*",),
        lane_filter="fermigier",
    ),
)


LEGACY_LOG_SOURCES = (
    LegacyLogSource(
        "global_triage_log",
        "artifacts/local/elliptic-curves/fermigier_global_triage.log",
        100_000,
    ),
    LegacyLogSource(
        "wide_funnel_log",
        "artifacts/local/elliptic-curves/fermigier_wide_funnel_200.log",
        5_000,
    ),
    LegacyLogSource(
        "broad_h5k_log",
        "artifacts/local/elliptic-curves/fermigier_broad_h5k.log",
        5_000,
    ),
    LegacyLogSource(
        "stratified_h20k_log",
        "artifacts/local/elliptic-curves/fermigier_stratified_h20k.log",
        20_000,
    ),
)


POSITIVE_CERTIFICATES = (
    {
        "id": "fermigier-E22",
        "adapter_u": "19754/39",
        "certified_rank_lower_bound": 22,
        "exceptional_quotient_rank_lower_bound": 10,
        "relative_path": "artifacts/generated-results/elliptic-curves/elliptic_fermigier_rank22_points.json",
        "rank_field": "certified_rank_lower_bound_in_this_artifact",
    },
    {
        "id": "fermigier-rank20-near-miss",
        "adapter_u": "28917/20",
        "certified_rank_lower_bound": 20,
        "exceptional_quotient_rank_lower_bound": 8,
        "relative_path": "artifacts/generated-results/elliptic-curves/elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json",
        "rank_field": "global_arithmetic.rank_lower_bound",
    },
)


@dataclass
class AccumulatedRow:
    cohorts: set[str] = field(default_factory=set)
    cheap_features: dict[str, dict[str, Any]] = field(default_factory=dict)
    observations: list[dict[str, Any]] = field(default_factory=list)
    observation_keys: set[str] = field(default_factory=set)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_parameter(value: str | Fraction, *, coordinate: str) -> tuple[int, int]:
    parameter = abs(Fraction(value))
    if coordinate == "u":
        parameter *= 2
    elif coordinate != "T":
        raise ValueError(f"unknown Fermigier coordinate {coordinate!r}")
    return parameter.numerator, parameter.denominator


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def parameter_record(parameter: tuple[int, int]) -> dict[str, Any]:
    normalized_t = Fraction(*parameter)
    adapter_u = normalized_t / 2
    square = normalized_t * normalized_t
    return {
        "adapter_u": fraction_text(adapter_u),
        "normalized_T": fraction_text(normalized_t),
        "normalized_T_squared": fraction_text(square),
        "projective_pair_T": [normalized_t.numerator, normalized_t.denominator],
        "sign_quotient": True,
        "coordinate_relation": "T=2*u; T and -T define the same tuple-model coefficients",
    }


def split_bucket(parameter: tuple[int, int], salt: str) -> str:
    digest = sha256(f"{salt}|{parameter[0]}/{parameter[1]}".encode()).digest()
    value = int.from_bytes(digest[:8], "big") % 100
    if value < 70:
        return "train"
    if value < 85:
        return "validation"
    return "internal_test"


def iter_path(value: Any, expression: str) -> Iterator[tuple[str, Any]]:
    parts = expression.split(".") if expression else []

    def visit(current: Any, remaining: list[str], prefix: list[str]) -> Iterator[tuple[str, Any]]:
        if not remaining:
            yield ".".join(prefix), current
            return
        head, *tail = remaining
        if head == "*":
            if isinstance(current, list):
                for index, child in enumerate(current):
                    yield from visit(child, tail, [*prefix, str(index)])
            elif isinstance(current, dict):
                for key in sorted(current):
                    yield from visit(current[key], tail, [*prefix, str(key)])
            return
        if isinstance(current, dict) and head in current:
            yield from visit(current[head], tail, [*prefix, head])

    yield from visit(value, parts, [])


def _scalar_snapshot(value: Any) -> dict[str, Any]:
    """Keep compact score/outcome scalars while dropping point coordinates."""

    exact_names = {
        "candidate_id",
        "label",
        "role",
        "height",
        "projective_height",
        "quartic_height_bound",
        "stable_numerical_rank",
        "numerical_rank",
        "exact_pool_point_count",
        "signed_quartic_points_found",
        "distinct_quartic_x_values",
        "new_x_values_beyond_visible_sections",
        "root_number",
        "log_conductor",
        "below_strict_log_conductor_target",
    }
    tokens = ("score", "rank_position", "timeout", "point_count")
    answer: dict[str, Any] = {}

    def visit(current: Any, prefix: tuple[str, ...]) -> None:
        if isinstance(current, dict):
            for key, child in current.items():
                if key in {"points", "selected_points", "basis", "minimal_model", "new_points"}:
                    continue
                visit(child, (*prefix, str(key)))
        elif isinstance(current, list):
            # Precision runs duplicate the same numerical rank.  Other long
            # arrays are provenance payloads, not row-level cheap features.
            if prefix and prefix[-1] == "precision_runs":
                for index, child in enumerate(current):
                    visit(child, (*prefix, str(index)))
        elif current is None or isinstance(current, (str, int, float, bool)):
            leaf = prefix[-1] if prefix else "value"
            if leaf in exact_names or leaf == "status" or any(token in leaf for token in tokens):
                answer[".".join(prefix)] = current

    visit(value, ())
    return dict(sorted(answer.items()))


def _find_parameter(record: dict[str, Any]) -> tuple[int, int] | None:
    for key in (
        "literal_shift_T",
        "parameter_t",
        "t",
        "parameter",
        "normalized_parameter_T",
        "family_parameter",
    ):
        value = record.get(key)
        if isinstance(value, (str, int)):
            return canonical_parameter(str(value), coordinate="T")
    value = record.get("adapter_u")
    if isinstance(value, (str, int)):
        return canonical_parameter(str(value), coordinate="u")
    return None


def _add_observation(row: AccumulatedRow, observation: dict[str, Any]) -> None:
    key = json.dumps(observation, sort_keys=True, separators=(",", ":"))
    if key not in row.observation_keys:
        row.observation_keys.add(key)
        row.observations.append(observation)


def _read_tabular_source(
    root: Path,
    source: TabularSource,
    rows: dict[tuple[int, int], AccumulatedRow],
) -> int:
    path = root / source.relative_path
    count = 0
    if source.kind == "pair_population":
        with path.open(newline="") as stream:
            for fields in csv.reader(stream, delimiter="\t"):
                if len(fields) != 2:
                    raise ValueError(f"malformed pair row in {source.relative_path}: {fields}")
                parameter = canonical_parameter(Fraction(int(fields[0]), int(fields[1])), coordinate="u")
                rows.setdefault(parameter, AccumulatedRow()).cohorts.add(source.name)
                count += 1
        return count
    if source.kind == "rank_score_no_header":
        with path.open(newline="") as stream:
            for fields in csv.reader(stream, delimiter="\t"):
                if len(fields) != 4:
                    raise ValueError(f"malformed score row in {source.relative_path}: {fields}")
                rank, numerator, denominator, score = fields
                parameter = canonical_parameter(Fraction(int(numerator), int(denominator)), coordinate="u")
                row = rows.setdefault(parameter, AccumulatedRow())
                row.cohorts.add(source.name)
                row.cheap_features[source.name] = {"rank": int(rank), "score": float(score)}
                count += 1
        return count
    if source.kind in {"dict_table", "comment_header_dict_table"}:
        with path.open(newline="") as stream:
            if source.kind == "comment_header_dict_table":
                header = stream.readline()
                if not header.startswith("# "):
                    raise ValueError(f"missing comment header in {source.relative_path}")
                fieldnames = header[2:].rstrip("\n").split("\t")
                records = csv.DictReader(stream, delimiter="\t", fieldnames=fieldnames)
            else:
                records = csv.DictReader(stream, delimiter="\t")
            for record in records:
                parameter = canonical_parameter(
                    Fraction(int(record["numerator"]), int(record["denominator"])),
                    coordinate="u",
                )
                row = rows.setdefault(parameter, AccumulatedRow())
                row.cohorts.add(source.name)
                row.cheap_features[source.name] = {
                    key: _parse_number(value)
                    for key, value in record.items()
                    if key not in {"numerator", "denominator"}
                }
                count += 1
        return count
    raise ValueError(f"unknown table kind {source.kind!r}")


def _parse_number(value: str) -> int | float | str:
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


LOG_RESULT = re.compile(
    r"^#\d+\s+u=(?P<u>\S+).*?rank>=(?P<rank>\?|\d+)\s+"
    r"qpts=(?P<qpts>\?|\d+)(?:\s+status=(?P<status>[A-Za-z_-]+))?"
    r"(?:\s+(?P<seconds>\d+(?:\.\d+)?)s)?"
)
STAGE = re.compile(r"^===\s+\S+\s+===\s+H=(?P<height>\d+)")


def parse_legacy_log(path: Path, default_height: int) -> Iterator[tuple[tuple[int, int], dict[str, Any]]]:
    height = default_height
    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        stage = STAGE.search(line)
        if stage:
            height = int(stage.group("height"))
            continue
        match = LOG_RESULT.search(line)
        if not match:
            continue
        parameter = canonical_parameter(match.group("u"), coordinate="u")
        rank = match.group("rank")
        qpts = match.group("qpts")
        explicit_status = match.group("status")
        status = explicit_status or ("completed" if qpts != "?" else "unknown")
        yield parameter, {
            "line": line_number,
            "height_bound": height,
            "status": status,
            "reported_rank_floor_uncertified": (
                None if rank == "?" or qpts == "?" else int(rank)
            ),
            "reported_quartic_point_count": None if qpts == "?" else int(qpts),
            "wall_seconds": None if match.group("seconds") is None else float(match.group("seconds")),
            "semantics": "legacy verbatim fields; no exact independence certificate recovered",
        }


def _read_json_source(
    root: Path,
    source: JsonSource,
    rows: dict[tuple[int, int], AccumulatedRow],
) -> int:
    path = root / source.relative_path
    if not path.exists():
        if source.required:
            raise FileNotFoundError(path)
        return 0
    document = json.loads(path.read_text())
    count = 0
    for expression in source.collections:
        for json_path, record in iter_path(document, expression):
            if not isinstance(record, dict):
                continue
            if source.lane_filter is not None and record.get("lane") != source.lane_filter:
                continue
            parameter = _find_parameter(record)
            if parameter is None:
                continue
            row = rows.setdefault(parameter, AccumulatedRow())
            row.cohorts.add(source.name)
            _add_observation(
                row,
                {
                    "source": source.name,
                    "json_path": json_path,
                    "recorded_fields": _scalar_snapshot(record),
                    "semantics": "source fields retained without promotion to a rank certificate",
                },
            )
            count += 1
    return count


def _lookup(document: dict[str, Any], expression: str) -> Any:
    current: Any = document
    for part in expression.split("."):
        current = current[part]
    return current


def build_corpus(
    root: Path,
    *,
    output: Path,
    summary_output: Path,
    split_salt: str = "fermigier-labelled-corpus-v1",
) -> dict[str, Any]:
    rows: dict[tuple[int, int], AccumulatedRow] = {}
    inputs: dict[str, dict[str, Any]] = {}
    source_counts: dict[str, int] = {}

    for source in TABULAR_SOURCES:
        path = root / source.relative_path
        if not path.is_file():
            raise FileNotFoundError(path)
        source_counts[source.name] = _read_tabular_source(root, source, rows)
        inputs[source.relative_path] = {"sha256": file_sha256(path), "kind": source.kind}

    provenance_path = root / HOT_NEIGHBORHOOD_PROVENANCE
    if not provenance_path.is_file():
        raise FileNotFoundError(provenance_path)
    provenance_parameters: set[tuple[int, int]] = set()
    provenance_record_count = 0
    with provenance_path.open(newline="") as stream:
        for record in csv.DictReader(stream, delimiter="\t"):
            parameter = canonical_parameter(
                Fraction(int(record["numerator"]), int(record["denominator"])),
                coordinate="u",
            )
            if "hot_neighborhood_population" not in rows.get(parameter, AccumulatedRow()).cohorts:
                raise AssertionError(
                    f"provenance parameter {parameter} is absent from the hot-neighborhood population"
                )
            provenance_parameters.add(parameter)
            provenance_record_count += 1
    expected_hot_parameters = {
        parameter
        for parameter, row in rows.items()
        if "hot_neighborhood_population" in row.cohorts
    }
    if provenance_parameters != expected_hot_parameters:
        raise AssertionError("hot-neighborhood provenance does not cover the population exactly")
    source_counts["hot_neighborhood_provenance"] = provenance_record_count
    inputs[HOT_NEIGHBORHOOD_PROVENANCE] = {
        "sha256": file_sha256(provenance_path),
        "kind": "companion_parameter_provenance",
        "record_count": provenance_record_count,
        "distinct_parameter_count": len(provenance_parameters),
        "join_key": "canonical adapter u, converted to nonnegative T=2*u",
    }
    for relative_path in SCORE_METADATA_SOURCES:
        path = root / relative_path
        if not path.is_file():
            raise FileNotFoundError(path)
        inputs[relative_path] = {
            "sha256": file_sha256(path),
            "kind": (
                "score_implementation"
                if path.suffix == ".cpp"
                else "companion_score_run_metadata"
            ),
        }
    for relative_path in CORPUS_IMPLEMENTATION_SOURCES:
        path = root / relative_path
        inputs[relative_path] = {
            "sha256": file_sha256(path),
            "kind": "corpus_implementation",
        }

    for source in JSON_SOURCES:
        path = root / source.relative_path
        source_counts[source.name] = _read_json_source(root, source, rows)
        if path.exists():
            inputs[source.relative_path] = {
                "sha256": file_sha256(path),
                "collections": list(source.collections),
            }

    for source in LEGACY_LOG_SOURCES:
        path = root / source.relative_path
        if not path.is_file():
            raise FileNotFoundError(path)
        count = 0
        for parameter, observation in parse_legacy_log(path, source.default_height):
            row = rows.setdefault(parameter, AccumulatedRow())
            row.cohorts.add(source.name)
            _add_observation(row, {"source": source.name, **observation})
            count += 1
        source_counts[source.name] = count
        inputs[source.relative_path] = {"sha256": file_sha256(path), "kind": "legacy_log"}

    positives: dict[tuple[int, int], dict[str, Any]] = {}
    for certificate in POSITIVE_CERTIFICATES:
        path = root / certificate["relative_path"]
        document = json.loads(path.read_text())
        found_rank = int(_lookup(document, certificate["rank_field"]))
        expected_rank = int(certificate["certified_rank_lower_bound"])
        if found_rank != expected_rank:
            raise AssertionError(
                f"{certificate['id']} certificate rank changed: {found_rank} != {expected_rank}"
            )
        parameter = canonical_parameter(certificate["adapter_u"], coordinate="u")
        row = rows.setdefault(parameter, AccumulatedRow())
        row.cohorts.add("certified_positive")
        positives[parameter] = {
            "state": "certified_positive",
            "generic_rank": 12,
            "certified_rank_lower_bound": expected_rank,
            "exceptional_quotient_rank_lower_bound": int(
                certificate["exceptional_quotient_rank_lower_bound"]
            ),
            "certificate": {
                "id": certificate["id"],
                "path": certificate["relative_path"],
                "sha256": file_sha256(path),
            },
        }
        inputs[certificate["relative_path"]] = {
            "sha256": file_sha256(path),
            "kind": "exact_positive_certificate",
        }

    split_counts: dict[str, int] = {}
    cohort_counts: dict[str, int] = {}
    feature_counts: dict[str, int] = {}
    observed_count = 0
    legacy_rank_observed_count = 0
    json_rank_observed_count = 0
    legacy_completed_count = 0
    legacy_error_count = 0
    legacy_zero_point_count = 0
    legacy_positive_point_count = 0
    legacy_rank_histogram: dict[str, int] = {}
    numerical_rank_histogram: dict[str, int] = {}
    label_counts = {"certified_positive": 0, "censored_control": 0}

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as stream:
                for parameter in sorted(
                    rows, key=lambda item: (max(item[0], item[1]), item[1], item[0])
                ):
                    accumulated = rows[parameter]
                    label = positives.get(
                        parameter,
                        {
                            "state": "censored_control",
                            "generic_rank": 12,
                            "certified_rank_lower_bound": None,
                            "exceptional_quotient_rank_lower_bound": None,
                            "reason": "no exact exceptional-rank certificate is attached",
                        },
                    )
                    split = (
                        "positive_holdout"
                        if label["state"] == "certified_positive"
                        else split_bucket(parameter, split_salt)
                    )
                    observations = sorted(
                        accumulated.observations,
                        key=lambda value: (
                            str(value.get("source", "")),
                            int(value.get("height_bound", -1)),
                            str(value.get("json_path", "")),
                            int(value.get("line", -1)),
                        ),
                    )
                    legacy_ranks = [
                        value["reported_rank_floor_uncertified"]
                        for value in observations
                        if value.get("reported_rank_floor_uncertified") is not None
                    ]
                    json_ranks = [
                        rank
                        for value in observations
                        for key, rank in value.get("recorded_fields", {}).items()
                        if key.endswith("stable_numerical_rank") and isinstance(rank, int)
                    ]
                    legacy_searches = [
                        value
                        for value in observations
                        if "reported_quartic_point_count" in value
                    ]
                    point_counts = [
                        value["reported_quartic_point_count"]
                        for value in legacy_searches
                        if value["reported_quartic_point_count"] is not None
                    ]
                    height_bounds = [
                        int(value["height_bound"])
                        for value in observations
                        if value.get("height_bound") is not None
                    ]
                    statuses = sorted(
                        {
                            str(value["status"])
                            for value in observations
                            if value.get("status") is not None
                        }
                    )
                    corpus_row = {
                        "schema": SCHEMA,
                        "family_id": FAMILY_ID,
                        "parameter": parameter_record(parameter),
                        "split": split,
                        "label": label,
                        "cohorts": sorted(accumulated.cohorts),
                        "cheap_features": dict(sorted(accumulated.cheap_features.items())),
                        "search_observations": observations,
                        "outcome_summary": {
                            "observation_count": len(observations),
                            "maximum_legacy_reported_rank_floor_uncertified": max(legacy_ranks, default=None),
                            "maximum_stable_numerical_rank": max(json_ranks, default=None),
                            "maximum_recorded_search_height": max(height_bounds, default=None),
                            "maximum_reported_quartic_point_count": max(point_counts, default=None),
                            "observed_statuses": statuses,
                            "rank_claim_boundary": "only label.certificate can establish exceptional rank",
                        },
                    }
                    stream.write(json.dumps(corpus_row, sort_keys=True, separators=(",", ":")) + "\n")
                    label_counts[label["state"]] += 1
                    split_counts[split] = split_counts.get(split, 0) + 1
                    for cohort in accumulated.cohorts:
                        cohort_counts[cohort] = cohort_counts.get(cohort, 0) + 1
                    for feature in accumulated.cheap_features:
                        feature_counts[feature] = feature_counts.get(feature, 0) + 1
                    if observations:
                        observed_count += 1
                    if legacy_ranks:
                        legacy_rank_observed_count += 1
                        rank_key = str(max(legacy_ranks))
                        legacy_rank_histogram[rank_key] = legacy_rank_histogram.get(rank_key, 0) + 1
                    if json_ranks:
                        json_rank_observed_count += 1
                        rank_key = str(max(json_ranks))
                        numerical_rank_histogram[rank_key] = numerical_rank_histogram.get(rank_key, 0) + 1
                    if any(value.get("status") == "completed" for value in legacy_searches):
                        legacy_completed_count += 1
                    if any(value.get("status") == "error" for value in legacy_searches):
                        legacy_error_count += 1
                    if point_counts and max(point_counts) == 0:
                        legacy_zero_point_count += 1
                    if point_counts and max(point_counts) > 0:
                        legacy_positive_point_count += 1

    summary = {
        "schema": "elliptic-curves.fermigier-labelled-corpus-summary.v1",
        "status": "EXPERIMENTAL_CENSORING_AWARE_CORPUS_TWO_CERTIFIED_POSITIVES",
        "family": {
            "id": FAMILY_ID,
            "generic_rank": 12,
            "canonical_coordinate": "nonnegative normalized T=abs(2*u)",
            "deduplication_key": "reduced normalized T after the exact sign quotient",
            "limitation": "this is not a global Q-isomorphism or twist-class audit",
        },
        "output": {
            "path": str(output.resolve()),
            "sha256": file_sha256(output),
            "compression": "deterministic gzip (mtime=0) over canonical compact JSONL",
            "row_count": len(rows),
        },
        "labels": label_counts,
        "splits": dict(sorted(split_counts.items())),
        "cohort_row_counts_after_deduplication": dict(sorted(cohort_counts.items())),
        "feature_row_counts_after_deduplication": dict(sorted(feature_counts.items())),
        "outcomes": {
            "rows_with_any_search_observation": observed_count,
            "rows_with_legacy_reported_rank_floor_uncertified": legacy_rank_observed_count,
            "rows_with_structured_stable_numerical_rank": json_rank_observed_count,
            "rows_with_completed_legacy_search": legacy_completed_count,
            "rows_with_legacy_error": legacy_error_count,
            "rows_with_zero_reported_quartic_points": legacy_zero_point_count,
            "rows_with_positive_reported_quartic_points": legacy_positive_point_count,
            "maximum_legacy_rank_histogram_uncertified": dict(
                sorted(legacy_rank_histogram.items(), key=lambda item: int(item[0]))
            ),
            "maximum_structured_numerical_rank_histogram": dict(
                sorted(numerical_rank_histogram.items(), key=lambda item: int(item[0]))
            ),
        },
        "source_record_counts_before_deduplication": dict(sorted(source_counts.items())),
        "inputs": dict(sorted(inputs.items())),
        "proof_boundary": [
            "Only the two certificate-backed rows are positive labels.",
            "All other rows are censored controls, never rank-zero or no-jump claims.",
            "Legacy rank>= fields and structured stable numerical ranks are observations only.",
            "A point-search miss is bounded by its recorded search and is not a nonexistence proof.",
            "Some legacy local TSV and log inputs lack their historical generators; the new family-residual table has complete code/run lineage.",
            "The complete 60,815,684-fibre global box remains a population formula/replay, not 60,815,684 materialized corpus labels.",
        ],
    }
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary
