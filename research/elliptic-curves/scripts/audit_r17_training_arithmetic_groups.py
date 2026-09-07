#!/usr/bin/env python3
"""Freeze and audit exact R17 isomorphism/twist groups before score reuse."""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
LOCAL = ROOT / "artifacts/local/elliptic-curves"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
POPULATION = LOCAL / "r17-training-population.jsonl"
SELECTED = LOCAL / "r17-training-selected.jsonl"
HOLDOUT = LOCAL / "r17-prospective-holdout.jsonl"
RANKER = ART / "r17_bisection_gain_ranker_quarantined_replay_v1.json"
PROSPECTIVE = ART / "r17_bisection_gain_ranker_prospective_holdout_v1.json"
COMPACT_INPUT = ART / "r17_training_arithmetic_group_inputs_v1.json.gz"
OUTPUT = ART / "r17_training_arithmetic_group_audit_v1.json"

SPLIT_SALT = "r17-training-v1"
EXPECTED_SOURCE_HASHES = {
    "artifacts/local/elliptic-curves/r17-training-population.jsonl": (
        "fec17217ae1349c1231711ac6bfb0d59b00f2214730948e0f1408b23f6b2ede0"
    ),
    "artifacts/local/elliptic-curves/r17-training-selected.jsonl": (
        "11ebd6e58f6bf82d51e688791d9f3585ffe2c3a8fc17137b195c907fa3a6c484"
    ),
    "artifacts/local/elliptic-curves/r17-prospective-holdout.jsonl": (
        "acbd9c534aac5e71090120326eb573275c8fbef0d25c19ca8172845eccc2cce1"
    ),
}
EXPECTED_MODEL_SHA256 = (
    "6f4d8d9d85a5b6880a1f6c1381c7ed8988972e8c6ba551109d745dc0013c9869"
)
EXPECTED_RANKER_SHA256 = (
    "f6f07cf5a19c71adb73ed193fcb98af6640b08749689646debff61201e9e4d36"
)
EXPECTED_PROSPECTIVE_SHA256 = (
    "70491a965cfbdf88c00f610a29809cadde6e78635404b7525b9103ae48b56c2c"
)
QUARANTINED_CONTROLS = ((-2, 377), (-308, 251), (2456, 135), (-9529, 5471))


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def split_bucket(parameter: tuple[int, int]) -> str:
    value = int.from_bytes(
        sha256(f"{SPLIT_SALT}|{parameter[0]}/{parameter[1]}".encode()).digest()[:8],
        "big",
    ) % 100
    if value < 70:
        return "train"
    if value < 85:
        return "validation"
    return "internal_test"


def freeze_compact_input() -> dict[str, Any]:
    for path in (POPULATION, SELECTED, HOLDOUT):
        expected = EXPECTED_SOURCE_HASHES[relative(path)]
        if not path.is_file() or digest(path) != expected:
            raise ArithmeticError(f"source-hash mismatch: {relative(path)}")
    population = read_jsonl(POPULATION)
    selected = {tuple(map(int, row["projective_pair"])) for row in read_jsonl(SELECTED)}
    holdout = {tuple(map(int, row["projective_pair"])) for row in read_jsonl(HOLDOUT)}
    population_parameters = {
        tuple(map(int, row["projective_pair"])) for row in population
    }
    if len(population_parameters) != len(population):
        raise ArithmeticError("the development population repeats a parameter")
    if not selected <= population_parameters or not holdout <= population_parameters:
        raise ArithmeticError("selected or holdout parameters left the frozen population")
    if selected & holdout:
        raise ArithmeticError("the prospective holdout intersects the labelled selection")
    compact = {
        "schema": "elliptic-curves.r17-training-arithmetic-group-inputs.v1",
        "status": "FROZEN_OUTCOME_FREE_PARAMETER_SPLIT_INPUTS",
        "split_salt": SPLIT_SALT,
        "source_hashes": EXPECTED_SOURCE_HASHES,
        "rows": [
            [
                int(row["projective_pair"][0]),
                int(row["projective_pair"][1]),
                row["split"],
                tuple(map(int, row["projective_pair"])) in selected,
                tuple(map(int, row["projective_pair"])) in holdout,
            ]
            for row in population
        ],
        "quarantined_controls": [list(row) for row in QUARANTINED_CONTROLS],
        "proof_boundary": (
            "Only parameters, deterministic split labels, and selection/holdout "
            "membership are frozen; no bisection, rank, or point-search outcome is read."
        ),
    }
    return compact


def write_compact_input(document: dict[str, Any]) -> None:
    COMPACT_INPUT.parent.mkdir(parents=True, exist_ok=True)
    COMPACT_INPUT.write_bytes(gzip.compress(canonical_bytes(document), mtime=0))


def load_compact_input() -> dict[str, Any]:
    with gzip.open(COMPACT_INPUT, "rt") as source:
        return json.load(source)


def homogeneous_value(coefficients: Iterable[int], a: int, b: int, degree: int) -> int:
    return sum(
        coefficient * a**power * b ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


def exact_j_key(
    a_coefficients: tuple[int, ...],
    b_coefficients: tuple[int, ...],
    parameter: tuple[int, int],
) -> tuple[int, int]:
    a, b = parameter
    coefficient_a = homogeneous_value(a_coefficients, a, b, 8)
    coefficient_b = homogeneous_value(b_coefficients, a, b, 12)
    denominator = 4 * coefficient_a**3 + 27 * coefficient_b**2
    if denominator == 0:
        raise ArithmeticError(f"singular R17 fibre at {a}/{b}")
    value = Fraction(6912 * coefficient_a**3, denominator)
    return value.numerator, value.denominator


def sequence_hash(rows: Iterable[Any]) -> str:
    answer = sha256()
    for row in rows:
        answer.update(json.dumps(row, sort_keys=True, separators=(",", ":")).encode())
        answer.update(b"\n")
    return answer.hexdigest()


def build() -> dict[str, Any]:
    compact = load_compact_input()
    if compact.get("schema") != "elliptic-curves.r17-training-arithmetic-group-inputs.v1":
        raise ArithmeticError("the compact grouping input has the wrong schema")
    if compact.get("status") != "FROZEN_OUTCOME_FREE_PARAMETER_SPLIT_INPUTS":
        raise ArithmeticError("the compact grouping input is not frozen")
    if compact.get("source_hashes") != EXPECTED_SOURCE_HASHES:
        raise ArithmeticError("the compact input names different raw sources")
    if compact.get("split_salt") != SPLIT_SALT:
        raise ArithmeticError("the compact input uses another split salt")
    if digest(MODEL) != EXPECTED_MODEL_SHA256:
        raise ArithmeticError("the exact R17 model changed")
    if digest(RANKER) != EXPECTED_RANKER_SHA256:
        raise ArithmeticError("the frozen learned-score artifact changed")
    if digest(PROSPECTIVE) != EXPECTED_PROSPECTIVE_SHA256:
        raise ArithmeticError("the frozen prospective score evaluation changed")

    model = json.loads(MODEL.read_text())
    a_coefficients = tuple(map(int, model["A_coefficients_low_to_high"]))
    b_coefficients = tuple(map(int, model["B_coefficients_low_to_high"]))
    rows = compact["rows"]
    if len(rows) != 100_000:
        raise ArithmeticError("the frozen development population is not 100,000 rows")

    parameter_rows: dict[tuple[int, int], tuple[str, bool, bool]] = {}
    j_groups: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    exact_rows = []
    for raw_a, raw_b, split, selected, holdout in rows:
        parameter = int(raw_a), int(raw_b)
        if gcd(abs(parameter[0]), parameter[1]) != 1 or parameter[1] <= 0:
            raise ArithmeticError("a compact population parameter is not normalized")
        if split != split_bucket(parameter):
            raise ArithmeticError("a parameter is assigned to the wrong deterministic split")
        if type(selected) is not bool or type(holdout) is not bool:
            raise ArithmeticError("selection flags must be boolean")
        if parameter in parameter_rows:
            raise ArithmeticError("the compact population repeats a parameter")
        parameter_rows[parameter] = split, selected, holdout
        j_key = exact_j_key(a_coefficients, b_coefficients, parameter)
        j_groups[j_key].append(parameter)
        exact_rows.append([*parameter, split, selected, holdout, *j_key])

    same_j_groups = [group for group in j_groups.values() if len(group) > 1]
    cross_split_groups = [
        group
        for group in same_j_groups
        if len({parameter_rows[parameter][0] for parameter in group}) > 1
    ]
    selected_parameters = {
        parameter for parameter, (_split, selected, _holdout) in parameter_rows.items()
        if selected
    }
    holdout_parameters = {
        parameter for parameter, (_split, _selected, holdout) in parameter_rows.items()
        if holdout
    }
    if len(selected_parameters) != 4_922 or len(holdout_parameters) != 5_000:
        raise ArithmeticError("the selected or holdout membership count changed")
    if selected_parameters & holdout_parameters:
        raise ArithmeticError("the labelled selection and prospective holdout overlap")

    selected_j = {
        exact_j_key(a_coefficients, b_coefficients, parameter)
        for parameter in selected_parameters
    }
    holdout_j = {
        exact_j_key(a_coefficients, b_coefficients, parameter)
        for parameter in holdout_parameters
    }
    controls = [tuple(map(int, row)) for row in compact["quarantined_controls"]]
    if tuple(controls) != QUARANTINED_CONTROLS:
        raise ArithmeticError("the quarantined control list changed")
    control_j = {
        exact_j_key(a_coefficients, b_coefficients, parameter)
        for parameter in controls
    }
    population_j = set(j_groups)
    if len(control_j) != len(controls):
        raise ArithmeticError("two quarantined controls share an exact twist class")

    pass_gate = (
        not same_j_groups
        and not cross_split_groups
        and not (selected_j & holdout_j)
        and not (control_j & population_j)
        and not (control_j & holdout_j)
    )
    if not pass_gate:
        raise ArithmeticError("exact isomorphism/twist leakage gate failed")

    definition = {
        "grouping_rule": (
            "group by the exact reduced rational j-invariant; this conservatively "
            "places every Q-isomorphism and every rational twist class together, "
            "including the exceptional j=0 and j=1728 cases"
        ),
        "development_population": {
            "row_count": len(rows),
            "exact_j_group_count": len(j_groups),
            "repeated_exact_j_group_count": len(same_j_groups),
            "cross_split_exact_j_group_count": len(cross_split_groups),
            "exact_parameter_split_j_sequence_sha256": sequence_hash(exact_rows),
        },
        "labelled_selection": {
            "row_count": len(selected_parameters),
            "exact_j_group_count": len(selected_j),
        },
        "prospective_holdout": {
            "row_count": len(holdout_parameters),
            "exact_j_group_count": len(holdout_j),
            "parameter_overlap_with_labelled_selection": 0,
            "twist_class_overlap_with_labelled_selection": len(selected_j & holdout_j),
        },
        "quarantined_controls": {
            "row_count": len(controls),
            "exact_j_group_count": len(control_j),
            "twist_class_overlap_with_development_population": len(
                control_j & population_j
            ),
            "twist_class_overlap_with_prospective_holdout": len(control_j & holdout_j),
        },
        "gate": {
            "status": "PASS_EXACT_ISOMORPHISM_TWIST_GROUPING",
            "historical_parameter_hash_split_closed": True,
            "learned_score_reuse_authorized": True,
            "authorized_score_artifact": relative(RANKER),
            "authorized_score_artifact_sha256": EXPECTED_RANKER_SHA256,
            "authorization_scope": (
                "reuse of the immutable v1 learned score; any changed population, "
                "split, holdout, control set, family model, or fitted score requires "
                "a new outcome-free exact grouping audit"
            ),
        },
        "temporal_boundary": (
            "The exact grouping audit is post-fit but reads only pre-outcome family "
            "parameters, split labels, and the exact family equation. It validates "
            "the immutable split before score reuse and does not authorize retuning."
        ),
    }
    body = {
        "schema": "elliptic-curves.r17-training-arithmetic-group-audit.v1",
        "status": "PASS_EXACT_GROUPING_BEFORE_FURTHER_LEARNED_SCORE_REUSE",
        "definition": definition,
        "definition_sha256": canonical_hash(definition),
        "inputs": {
            relative(COMPACT_INPUT): digest(COMPACT_INPUT),
            relative(MODEL): digest(MODEL),
            relative(RANKER): digest(RANKER),
            relative(PROSPECTIVE): digest(PROSPECTIVE),
        },
        "generation": {
            "script": relative(Path(__file__)),
            "command": (
                "python3 elliptic-curves/scripts/"
                "audit_r17_training_arithmetic_groups.py"
            ),
        },
        "claim_boundary": [
            "Exact equal-j grouping is a leakage audit, not a rank or Selmer calculation.",
            "The historical fitted coefficients and opened evaluations remain immutable.",
            "Authorization applies only to the pinned v1 score and exact input universe.",
        ],
    }
    return body


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--freeze-inputs", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.freeze_inputs:
        if args.check:
            raise SystemExit("--freeze-inputs and --check are mutually exclusive")
        write_compact_input(freeze_compact_input())
    document = build()
    encoded = canonical_bytes(document)
    if args.check:
        if not args.output.is_file() or args.output.read_bytes() != encoded:
            raise SystemExit(f"stale or missing arithmetic-group audit: {args.output}")
        print(f"R17GROUPS|status=PASS|sha256={sha256(encoded).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(encoded)
    print(f"R17GROUPS|status=WROTE|sha256={sha256(encoded).hexdigest()}")


if __name__ == "__main__":
    main()
