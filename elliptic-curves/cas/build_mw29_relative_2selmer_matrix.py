#!/usr/bin/env python3
"""Build the exact local-condition matrix after quotienting known MW classes.

This is a backend-independent proof gate.  A number-field worker supplies an
ambient norm-squareclass envelope, certified Kummer coordinates for the known
points, and exact allowed subspaces at the relevant places.  This script first
quotients the known Mordell--Weil subspace and only then forms local equations.

An incomplete relation collection is useful input and produces a reproducible
fingerprint, but it cannot certify a Selmer upper bound.  Residual-zero status
is emitted only when the manifest explicitly certifies that the global
envelope and every required local condition are exhaustive.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Mapping


CAS = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS))

from residual_selmer_quotient import (  # noqa: E402
    F2Error,
    build_relative_local_condition_matrix,
)


INPUT_SCHEMA = "elliptic-curves.mw29-relative-2selmer-ambient.v1"
OUTPUT_SCHEMA = "elliptic-curves.mw29-relative-2selmer-matrix.v1"
PROTOCOL = "MW29REL2MAT-v1"
REQUIRED_UPPER_BOUND_FLAGS = (
    "global_ambient_upper_envelope_certified",
    "norm_condition_incorporated",
    "known_mw_kummer_coordinates_certified",
    "supplied_local_conditions_certified",
    "supplied_subspaces_are_necessary_selmer_conditions",
)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rows(records: object) -> tuple[list[str], list[list[int]]]:
    if not isinstance(records, list):
        raise F2Error("known_mw_rows must be a list")
    labels: list[str] = []
    rows: list[list[int]] = []
    for index, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise F2Error(f"known_mw_rows[{index}] is not an object")
        if "label" not in record or "row" not in record:
            raise F2Error(f"known_mw_rows[{index}] needs label and row")
        row = record["row"]
        if not isinstance(row, list):
            raise F2Error(f"known_mw_rows[{index}].row is not a list")
        labels.append(str(record["label"]))
        values = [int(value) for value in row]
        if any(value not in (0, 1) for value in values):
            raise F2Error(f"known_mw_rows[{index}].row contains a non-binary entry")
        rows.append(values)
    if len(set(labels)) != len(labels):
        raise F2Error("known Mordell--Weil labels are not unique")
    return labels, rows


def build_certificate(
    manifest: Mapping[str, object], *, maximum_cut_size: int,
    maximum_cut_subsets: int = 1_000_000,
) -> Mapping[str, object]:
    if manifest.get("schema") != INPUT_SCHEMA:
        raise F2Error(f"input schema must be {INPUT_SCHEMA}")
    if "ambient_norm_square_dimension" not in manifest:
        raise F2Error("manifest is missing ambient_norm_square_dimension")
    ambient_dimension = int(manifest["ambient_norm_square_dimension"])
    labels, known_rows = _rows(manifest.get("known_mw_rows"))
    places = manifest.get("places")
    if not isinstance(places, list):
        raise F2Error("places must be a list")
    place_labels = [
        str(place.get("place", f"place-{index}"))
        if isinstance(place, Mapping)
        else f"place-{index}"
        for index, place in enumerate(places)
    ]
    if len(set(place_labels)) != len(place_labels):
        raise F2Error("place labels are not unique")

    certification = manifest.get("certification")
    if not isinstance(certification, Mapping):
        raise F2Error("certification must be an object")
    flags: dict[str, bool] = {
        name: certification.get(name) is True
        for name in REQUIRED_UPPER_BOUND_FLAGS
    }
    flags["all_required_local_conditions_complete"] = (
        certification.get("all_required_local_conditions_complete") is True
    )
    flags["global_ambient_exact"] = certification.get("global_ambient_exact") is True
    flags["residual_dimension_parity_certified"] = certification.get(
        "residual_dimension_parity_certified"
    ) is True
    method = str(certification.get("method", "")).strip()
    hypothesis_value = certification.get("hypothesis")
    hypothesis = None if hypothesis_value is None else str(hypothesis_value).strip()
    upper_bound_certified = all(
        flags[name] for name in REQUIRED_UPPER_BOUND_FLAGS
    ) and bool(method)
    complete_local_coverage = flags["all_required_local_conditions_complete"]

    matrix = build_relative_local_condition_matrix(
        ambient_dimension=ambient_dimension,
        known_mw_rows=known_rows,
        places=places,
        maximum_cut_size=maximum_cut_size,
        maximum_cut_subsets=maximum_cut_subsets,
    )
    known_rank = int(matrix["known_mw_kummer_dimension"])
    raw_residual_upper = int(matrix["unexplained_selmer_excess_kernel_dimension"])
    target_rank = int(manifest.get("known_mw_target_rank", 29))
    if target_rank < 0:
        raise F2Error("known_mw_target_rank must be nonnegative")

    parity_value = manifest.get("residual_selmer_dimension_parity")
    if flags["residual_dimension_parity_certified"]:
        if parity_value is None or int(parity_value) not in (0, 1):
            raise F2Error("a certified residual parity must be 0 or 1")
        residual_parity = int(parity_value)
        sharpened_residual_upper = raw_residual_upper - (
            (raw_residual_upper - residual_parity) & 1
        )
        if sharpened_residual_upper < 0:
            raise F2Error("residual parity contradicts the known Selmer lower bound")
    else:
        residual_parity = None
        sharpened_residual_upper = raw_residual_upper
    exact_presentation = (
        upper_bound_certified
        and flags["global_ambient_exact"]
        and complete_local_coverage
    )
    if exact_presentation and residual_parity is not None and (
        raw_residual_upper & 1
    ) != residual_parity:
        raise F2Error("exact residual kernel contradicts the certified parity")
    residual_zero = sharpened_residual_upper == 0

    if known_rank != target_rank:
        status = "INCOMPLETE_KNOWN_MW_QUOTIENT"
    elif upper_bound_certified and residual_zero:
        status = (
            "CERTIFIED_RELATIVE_2SELMER_QUOTIENT_ZERO_UNDER_HYPOTHESIS"
            if hypothesis
            else "CERTIFIED_RELATIVE_2SELMER_QUOTIENT_ZERO"
        )
    elif upper_bound_certified and exact_presentation:
        status = (
            "CERTIFIED_RELATIVE_2SELMER_KERNEL_COMPUTED_UNDER_HYPOTHESIS"
            if hypothesis
            else "CERTIFIED_RELATIVE_2SELMER_KERNEL_COMPUTED"
        )
    elif upper_bound_certified:
        status = (
            "CERTIFIED_RELATIVE_2SELMER_UPPER_BOUND_UNDER_HYPOTHESIS"
            if hypothesis
            else "CERTIFIED_RELATIVE_2SELMER_UPPER_BOUND"
        )
    else:
        status = "INCOMPLETE_RELATIVE_2SELMER_MATRIX"

    missing_flags = [
        name for name in REQUIRED_UPPER_BOUND_FLAGS if not flags[name]
    ]
    return {
        "schema": OUTPUT_SCHEMA,
        "protocol": PROTOCOL,
        "case_id": manifest.get("case_id"),
        "known_mw_target_rank": target_rank,
        "known_mw_labels": labels,
        "certification": {
            **flags,
            "method": method or None,
            "hypothesis": hypothesis or None,
            "relative_upper_bound_certified": upper_bound_certified,
            "complete_local_coverage": complete_local_coverage,
            "missing_or_false_requirements": missing_flags
            + ([] if method else ["method"]),
        },
        "relative_local_matrix": matrix,
        "relative_selmer_bound": {
            "raw_upper_bound_from_supplied_local_matrix": raw_residual_upper,
            "certified_residual_dimension_parity": residual_parity,
            "parity_sharpened_upper_bound": sharpened_residual_upper,
            "global_presentation_exact": flags["global_ambient_exact"],
            "all_required_local_conditions_complete": complete_local_coverage,
            "residual_kernel_exact": exact_presentation,
        },
        "status": status,
        "claim_boundary": (
            "A zero kernel proves Sel_2(E)/im(MW/2MW)=0 from any certified "
            "subset of local conditions when relative_upper_bound_certified "
            "is true and the certified known Mordell--Weil rank equals "
            "known_mw_target_rank. A nonzero kernel is the exact residual "
            "Selmer group only when both global_ambient_exact and "
            "complete_local_coverage are true. A certified parity may sharpen "
            "an upper bound by one. "
            "Auxiliary fingerprints are not eligible local subspaces."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--maximum-cut-size", type=int, default=8)
    parser.add_argument("--maximum-cut-subsets", type=int, default=1_000_000)
    args = parser.parse_args()
    manifest = json.loads(args.input.read_text())
    if not isinstance(manifest, Mapping):
        raise F2Error("top-level JSON value must be an object")
    result = dict(
        build_certificate(
            manifest,
            maximum_cut_size=args.maximum_cut_size,
            maximum_cut_subsets=args.maximum_cut_subsets,
        )
    )
    result["input"] = str(args.input)
    result["input_sha256"] = file_sha256(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    matrix = result["relative_local_matrix"]
    print(
        f"{PROTOCOL}|status={result['status']}"
        f"|known={matrix['known_mw_kummer_dimension']}"
        f"|ambient_quotient={matrix['mw_quotient_ambient_dimension']}"
        f"|kernel={matrix['unexplained_selmer_excess_kernel_dimension']}"
        f"|output={args.output}"
    )


if __name__ == "__main__":
    main()
