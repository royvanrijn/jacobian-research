#!/usr/bin/env python3
"""Freeze prior-parameter exclusions and raw scopes near Fermigier rank 20.

This is a read-only audit plus exact population enumeration.  It makes no
conductor, point-search, or rank call.  Every legacy Fermigier ``t`` value is
converted from the literal symmetric shift to the canonical adapter
``u=t/2``; imported adapter values are retained directly.  Signs are
quotiented because the canonical family coefficients are even in ``u``.

The two new raw populations are disjoint by denominator:

* dense window: primitive ``1 <= b <= 1200`` and
  ``|20a-28917b| <= 40b``;
* deep strip: primitive ``1201 <= b <= 20000`` and
  ``a=nearest_half_up(28917b/20)+delta``, ``|delta| <= 64``.

All imported and prior parameters are removed before the evaluated manifests
are hashed.  The imported high-power CRT seed is audited and excluded here;
its separate root-ball lane must pass an exact conductor/radical gate.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from math import gcd
import os
from pathlib import Path
import platform
import sys
from typing import Any, Iterable

from search_mestre_root_tuple_scale import sha256_file
from search_mestre_root_tuple_scale_max100 import stable_json_digest


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

ANCHOR = Q(28_917, 20)
DENSE_DENOMINATOR_BOUND = 1_200
DEEP_FIRST_DENOMINATOR = 1_201
DEEP_DENOMINATOR_BOUND = 20_000
DEEP_OFFSET_BOUND = 64
HIGH_POWER_SEED = Q(673_709, 29_965)
EXPECTED_INPUT_INVENTORY_SHA256 = (
    "c7a325864d062cf11324880fe7268446d2c48c68cdb0bdc91173d0a5e45393e2"
)
DEFAULT_OUTPUT = Path(
    "artifacts/local/elliptic-curves/"
    "fermigier_rank20_adapter_neighborhood_audit.json"
)

LEGACY_ARTIFACT_PATHS = (
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_1666_9.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_3115_3_h1000000_checkpoint.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_batch_rank_triage.json",
    "artifacts/generated-results/elliptic-curves/elliptic_fermigier_benchmark.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_crt_lattice_pilot.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_discovered_local_conditions.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_extra_points.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_global.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_multiple_root_crt.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_multiple_root_frontier.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_multiple_root_height_h50000.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_power_pairs.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_published_pair_fiber_products.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_published_pair_fiber_products_h1000000.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_published_pair_fiber_products_h1000000_charts.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_published_pair_fiber_products_h50000.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank20_28917_20_explicit_formula_delta22.json",
    "artifacts/generated-results/elliptic-curves/elliptic_fermigier_rank22_accidental_slices.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank22_auxiliary_orbits.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank22_auxiliary_orbits_l1_7.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank22_missing_preimage_slices.json",
    "artifacts/generated-results/elliptic-curves/elliptic_fermigier_rank22_points.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank22_record_group_directions.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank22_record_group_quads.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank22_record_group_triples.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank22_record_group_triples_remainder.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_record_rescore_h5000.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_record_residue_class.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_record_residue_deep_tranche.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_score_cutoffs.json",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank22_record_group_quads_stream.jsonl",
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank22_record_group_triples_remainder_stream.jsonl",
)

IMPORTED_AND_REFERENCE_PATHS = (
    "artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json",
    "artifacts/generated-results/elliptic-curves/fermigier_rank_certificates_v1.json",
    "artifacts/generated-results/elliptic-curves/fermigier_crt_seed_v1.json",
    "artifacts/generated-results/elliptic-curves/crt_lattice_calibration_v1.json",
    "elliptic-curves/families/fermigier_mestre_rank12.json",
    "elliptic-curves/notes/FERMIGIER_REPRODUCTION.md",
    "elliptic-curves/notes/CRT_LATTICE_PIPELINE.md",
    "elliptic-curves/ecsearch/fermigier.py",
    "elliptic-curves/ecsearch/fermigier_near_miss.py",
    "elliptic-curves/ecsearch/fermigier_seed.py",
    "elliptic-curves/ecsearch/fermigier_rank.py",
    "elliptic-curves/ecsearch/rank_certification.py",
    "elliptic-curves/ecsearch/crt_lattice.py",
    "elliptic-curves/ecsearch/local_data.py",
    "elliptic-curves/scripts/run_fermigier_rank20_near_miss.py",
    "elliptic-curves/scripts/run_fermigier_rank_certificates.py",
    "elliptic-curves/scripts/run_fermigier_crt_seed.py",
    "elliptic-curves/scripts/verify_fermigier_rank20_near_miss.py",
    "elliptic-curves/scripts/verify_fermigier_rank_certificates.py",
    "elliptic-curves/scripts/verify_fermigier_crt_seed.py",
    "elliptic-curves/cas/explicit_formula_fermigier_rank20_28917_20_delta22.py",
)


def fraction_text(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def canonical_adapter(value: Fraction) -> Fraction:
    return abs(Q(value))


def parse_fraction(value: Any) -> Fraction | None:
    if isinstance(value, bool) or isinstance(value, float):
        return None
    if isinstance(value, int):
        return Q(value)
    if not isinstance(value, str):
        return None
    try:
        return Q(value)
    except (ValueError, ZeroDivisionError):
        return None


def parameter_key_kind(key: str) -> str | None:
    lower = key.lower()
    if lower == "t" or "parameter" in lower:
        if "adapter" in lower or lower.endswith("_u") or lower == "published_parameter":
            return "adapter"
        return "literal_shift"
    if "literal_shift" in lower:
        return "literal_shift"
    return None


def extract_legacy_parameters(
    value: Any,
    *,
    source: str,
    path: str = "",
) -> list[tuple[Fraction, str]]:
    answer: list[tuple[Fraction, str]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            child_path = f"{path}.{key}" if path else key
            kind = parameter_key_kind(key)
            parsed = parse_fraction(child)
            if kind is not None and parsed is not None:
                adapter = parsed if kind == "adapter" else parsed / 2
                answer.append(
                    (canonical_adapter(adapter), f"{source}:{child_path}:{kind}")
                )
            answer.extend(
                extract_legacy_parameters(child, source=source, path=child_path)
            )
    elif isinstance(value, list):
        for index, child in enumerate(value):
            answer.extend(
                extract_legacy_parameters(
                    child,
                    source=source,
                    path=f"{path}[{index}]",
                )
            )
    return answer


def load_json_or_lines(path: Path) -> Iterable[Any]:
    if path.suffix == ".jsonl":
        with path.open() as stream:
            for line in stream:
                if line.strip():
                    yield json.loads(line)
    else:
        yield json.loads(path.read_text())


def input_inventory(root: Path) -> tuple[list[dict[str, Any]], str]:
    records = []
    digest = hashlib.sha256()
    for relative in (*LEGACY_ARTIFACT_PATHS, *IMPORTED_AND_REFERENCE_PATHS):
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(f"missing audited Fermigier input: {relative}")
        file_hash = sha256_file(path)
        digest.update(f"{relative}\0{file_hash}\n".encode())
        records.append(
            {"path": relative, "sha256": file_hash, "size_bytes": path.stat().st_size}
        )
    return records, digest.hexdigest()


def prior_parameter_sources(root: Path) -> dict[Fraction, set[str]]:
    sources: dict[Fraction, set[str]] = {}

    def add(value: Fraction, source: str) -> None:
        parameter = canonical_adapter(value)
        sources.setdefault(parameter, set()).add(source)

    for relative in LEGACY_ARTIFACT_PATHS:
        path = root / relative
        for document_index, document in enumerate(load_json_or_lines(path)):
            source = relative
            if path.suffix == ".jsonl":
                source += f":line{document_index + 1}"
            for parameter, location in extract_legacy_parameters(
                document, source=source
            ):
                add(parameter, location)

    near_miss_path = root / IMPORTED_AND_REFERENCE_PATHS[0]
    near_miss = json.loads(near_miss_path.read_text())
    add(
        Q(near_miss["family"]["adapter_parameter"]),
        f"{IMPORTED_AND_REFERENCE_PATHS[0]}:family.adapter_parameter:adapter",
    )
    add(
        Q(near_miss["family"]["literal_shift"]) / 2,
        f"{IMPORTED_AND_REFERENCE_PATHS[0]}:family.literal_shift:literal_shift",
    )

    rank_certificate_path = root / IMPORTED_AND_REFERENCE_PATHS[1]
    rank_certificate = json.loads(rank_certificate_path.read_text())
    add(
        Q(rank_certificate["generic_sections"]["adapter_parameter"]),
        f"{IMPORTED_AND_REFERENCE_PATHS[1]}:generic_sections.adapter_parameter:adapter",
    )
    add(
        Q(rank_certificate["generic_sections"]["literal_shift"]) / 2,
        f"{IMPORTED_AND_REFERENCE_PATHS[1]}:generic_sections.literal_shift:literal_shift",
    )

    seed_path = root / IMPORTED_AND_REFERENCE_PATHS[2]
    seed = json.loads(seed_path.read_text())
    for index, combination in enumerate(seed["search"]["combinations"]):
        add(
            Q(combination["numerator"], combination["denominator"]),
            f"{IMPORTED_AND_REFERENCE_PATHS[2]}:search.combinations[{index}]:adapter",
        )
    add(
        Q(seed["best_seed"]["numerator"], seed["best_seed"]["denominator"]),
        f"{IMPORTED_AND_REFERENCE_PATHS[2]}:best_seed:adapter",
    )

    family_path = root / "elliptic-curves/families/fermigier_mestre_rank12.json"
    family = json.loads(family_path.read_text())
    for label in ("benchmark_specialization", "rank20_near_miss"):
        record = family[label]
        add(
            Q(record["adapter_parameter_u"]),
            f"{family_path.relative_to(root)}:{label}.adapter_parameter_u:adapter",
        )
        add(
            Q(record["literal_shift_s"]) / 2,
            f"{family_path.relative_to(root)}:{label}.literal_shift_s:literal_shift",
        )

    reproduction = (root / "elliptic-curves/notes/FERMIGIER_REPRODUCTION.md").read_text()
    crt_note = (root / "elliptic-curves/notes/CRT_LATTICE_PIPELINE.md").read_text()
    for literal in ("u=28917/20", "u=19754/39", "s=39508/39"):
        if literal not in reproduction:
            raise AssertionError(f"the reproduction note lost {literal}")
    if "u=673709/29965" not in crt_note:
        raise AssertionError("the CRT note lost its imported Fermigier seed")
    add(ANCHOR, "elliptic-curves/notes/FERMIGIER_REPRODUCTION.md:u=28917/20")
    add(Q(19_754, 39), "elliptic-curves/notes/FERMIGIER_REPRODUCTION.md:u=19754/39")
    add(HIGH_POWER_SEED, "elliptic-curves/notes/CRT_LATTICE_PIPELINE.md:u=673709/29965")
    return sources


def dense_scope() -> Iterable[tuple[int, int]]:
    for denominator in range(1, DENSE_DENOMINATOR_BOUND + 1):
        first = (28_877 * denominator + 19) // 20
        last = (28_957 * denominator) // 20
        for numerator in range(first, last + 1):
            if gcd(numerator, denominator) == 1:
                yield numerator, denominator


def deep_scope() -> Iterable[tuple[int, int]]:
    for denominator in range(DEEP_FIRST_DENOMINATOR, DEEP_DENOMINATOR_BOUND + 1):
        center = (28_917 * denominator + 10) // 20
        for offset in range(-DEEP_OFFSET_BOUND, DEEP_OFFSET_BOUND + 1):
            numerator = center + offset
            if gcd(numerator, denominator) == 1:
                yield numerator, denominator


def scope_record(
    label: str,
    pairs: Iterable[tuple[int, int]],
    exclusions: set[Fraction],
) -> tuple[dict[str, Any], list[Fraction]]:
    raw_digest = hashlib.sha256()
    evaluated_digest = hashlib.sha256()
    excluded_digest = hashlib.sha256()
    raw_count = 0
    evaluated_count = 0
    excluded = []
    previous: Fraction | None = None
    for numerator, denominator in pairs:
        parameter = Q(numerator, denominator)
        if parameter.denominator != denominator or parameter <= 0:
            raise AssertionError("a raw adapter pair was not canonical positive")
        # Both iterators are ordered by denominator then numerator, not by value;
        # uniqueness follows from canonical denominators.
        if previous == parameter:
            raise AssertionError("a duplicate adjacent adapter parameter escaped")
        previous = parameter
        text = f"{numerator}/{denominator}\n"
        raw_digest.update(text.encode())
        raw_count += 1
        if parameter in exclusions:
            excluded_digest.update(text.encode())
            excluded.append(parameter)
        else:
            evaluated_digest.update(text.encode())
            evaluated_count += 1
    return (
        {
            "label": label,
            "raw_primitive_parameter_count": raw_count,
            "raw_parameter_manifest_sha256": raw_digest.hexdigest(),
            "prior_parameter_exclusion_count": len(excluded),
            "prior_parameter_exclusion_manifest_sha256": excluded_digest.hexdigest(),
            "evaluated_parameter_count": evaluated_count,
            "evaluated_parameter_manifest_sha256": evaluated_digest.hexdigest(),
            "excluded_parameters": [fraction_text(value) for value in excluded],
        },
        excluded,
    )


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.output.exists():
        raise SystemExit("refusing to overwrite the Fermigier neighborhood audit")
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    inventory, inventory_digest = input_inventory(root)
    if inventory_digest != EXPECTED_INPUT_INVENTORY_SHA256:
        raise AssertionError("the audited Fermigier input inventory changed")

    sources = prior_parameter_sources(root)
    exclusions = set(sources)
    if ANCHOR not in exclusions or HIGH_POWER_SEED not in exclusions:
        raise AssertionError("an imported anchor/CRT seed escaped the exclusion set")
    ordered_exclusions = sorted(exclusions, key=lambda value: (value.denominator, value.numerator))
    exclusion_digest = hashlib.sha256(
        "".join(f"{fraction_text(value)}\n" for value in ordered_exclusions).encode()
    ).hexdigest()

    dense, dense_excluded = scope_record("dense-window", dense_scope(), exclusions)
    deep, deep_excluded = scope_record("deep-strip", deep_scope(), exclusions)
    if set(dense_excluded) & set(deep_excluded):
        raise AssertionError("the disjoint denominator scopes acquired an overlap")
    raw_union_digest = hashlib.sha256(
        (
            f"dense-window:{dense['raw_parameter_manifest_sha256']}\n"
            f"deep-strip:{deep['raw_parameter_manifest_sha256']}\n"
        ).encode()
    ).hexdigest()
    evaluated_union_digest = hashlib.sha256(
        (
            f"dense-window:{dense['evaluated_parameter_manifest_sha256']}\n"
            f"deep-strip:{deep['evaluated_parameter_manifest_sha256']}\n"
        ).encode()
    ).hexdigest()

    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "complete read-only prior-parameter audit and exact raw population freeze",
        "claim_level": "bounded search population definition; no conductor or rank claim",
        "coordinate_normalization": {
            "canonical_coordinate": "adapter u=s/2",
            "legacy_conversion": "every legacy literal-shift t becomes |t|/2",
            "sign_quotient": "u and -u define the same canonical family coefficients",
            "factor_two_discrepancy_remains_unresolved": True,
        },
        "anchor": {
            "adapter_parameter": fraction_text(ANCHOR),
            "literal_shift": fraction_text(2 * ANCHOR),
            "excluded": True,
            "exact_rank_lower_bound": 20,
            "log_conductor": "159.93482522552545339840138409496203860033602762598444576466041725048663568492918",
            "root_number": 1,
            "conditional_delta_11_over_5_upper_under_grh": "21.0335328229846198389",
            "fixed_fiber_search_conditionally_closed": True,
        },
        "input_inventory": {
            "file_count": len(inventory),
            "combined_path_and_sha256": inventory_digest,
            "expected_combined_path_and_sha256": EXPECTED_INPUT_INVENTORY_SHA256,
            "files": inventory,
            "includes_every_preexisting_elliptic_fermigier_json_or_jsonl": True,
            "includes_imported_ecsearch_notes_family_and_crt_inputs": True,
            "crt_lattice_calibration_is_audited_but_outside_the_Fermigier_family": True,
        },
        "prior_parameter_exclusion": {
            "canonical_adapter_parameter_count": len(ordered_exclusions),
            "canonical_adapter_parameter_sha256": exclusion_digest,
            "parameters": [
                {
                    "adapter_u": fraction_text(value),
                    "source_count": len(sources[value]),
                    "sources": sorted(sources[value]),
                }
                for value in ordered_exclusions
            ],
        },
        "new_raw_scopes": {
            "dense_window": {
                **dense,
                "definition": (
                    "gcd(a,b)=1, 1<=b<=1200, |20a-28917b|<=40b"
                ),
            },
            "deep_strip": {
                **deep,
                "definition": (
                    "gcd(a,b)=1, 1201<=b<=20000, "
                    "a=nearest_half_up(28917b/20)+delta, |delta|<=64"
                ),
            },
            "raw_populations_disjoint_by_denominator": True,
            "raw_union_primitive_parameter_count": (
                dense["raw_primitive_parameter_count"]
                + deep["raw_primitive_parameter_count"]
            ),
            "raw_union_manifest_of_manifest_sha256": raw_union_digest,
            "evaluated_union_parameter_count": (
                dense["evaluated_parameter_count"]
                + deep["evaluated_parameter_count"]
            ),
            "evaluated_union_manifest_of_manifest_sha256": evaluated_union_digest,
        },
        "high_power_seed_lane": {
            "imported_seed_adapter_u": fraction_text(HIGH_POWER_SEED),
            "imported_seed_denominator": HIGH_POWER_SEED.denominator,
            "outside_broad_scope_because_denominator_exceeds_20000": True,
            "excluded_from_broad_scope_and_reserved_for_disjoint_CRT_root_ball": True,
            "required_promotion_gate": (
                "exact homogeneous-discriminant radical/cofactor features and exact "
                "minimal conductor must close before any point search"
            ),
        },
        "search_boundary": {
            "expensive_external_calls": 0,
            "conductor_calls": 0,
            "point_or_rank_calls": 0,
            "selection_outcomes_used": [],
            "unsearched": (
                "all adapter parameters outside the two declared broad scopes; all "
                "future discovery non-survivors; the separately owned CRT root-ball"
            ),
        },
        "provenance": {
            "script": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "python": platform.python_version(),
            "external_process_calls": 0,
            "same_stage_retries": 0,
            "owned_processes_remaining": 0,
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "normalization": artifact["coordinate_normalization"],
            "anchor": artifact["anchor"],
            "inventory": artifact["input_inventory"],
            "exclusions": artifact["prior_parameter_exclusion"],
            "scopes": artifact["new_raw_scopes"],
            "seed": artifact["high_power_seed_lane"],
            "boundary": artifact["search_boundary"],
        }
    )
    exclusive_write(args.output, artifact)
    print(
        "audit complete "
        f"prior={len(ordered_exclusions)} "
        f"raw={artifact['new_raw_scopes']['raw_union_primitive_parameter_count']} "
        f"evaluated={artifact['new_raw_scopes']['evaluated_union_parameter_count']} "
        f"digest={evaluated_union_digest} output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
