#!/usr/bin/env python3
"""Audit the exact P.O=0 singleton-twist screen on the top twenty candidates."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
OUTPUT = GENERATED / (
    "elkies-k3-r17-norm12-11952-singleton-twist-po0-top20-v1.json"
)

DISCOVERY = GENERATED / (
    "elkies-k3-r17-norm12-11952-singleton-twist-census-top200-v1.json"
)
HOLDOUTS = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-singleton-twist-census-top20-p499-821-v1.json",
    GENERATED
    / "elkies-k3-r17-norm12-11952-singleton-twist-census-top20-p823-1151-v1.json",
)
KNOWN_REGRESSION = GENERATED / (
    "elkies-k3-r17-norm12-11952-singleton-twist-po0-"
    "orbit1c3d5-known-p23-hensel512-v1.json"
)

# One exhaustive usable-prime shell for every member of the discovery top 20.
# The order is the discovery-census order and is therefore part of the audit.
PRIMARY_SHELLS = (
    ("119aa", 23),
    ("131dc", 23),
    ("1f77b", 23),
    ("1c2c0", 23),
    ("1a20a", 23),
    ("0cd49", 23),
    ("1e4d6", 37),
    ("1a49d", 23),
    ("18979", 23),
    ("0ee78", 37),
    ("08d9d", 23),
    ("16abe", 17),
    ("01247", 29),
    ("0943c", 23),
    ("16cbc", 31),
    ("1a803", 37),
    ("12e7f", 37),
    ("113b1", 41),
    ("1c3d5", 23),
    ("0889c", 23),
)

BRUTE_SCHEMA = "elkies-k3.twist-polynomial-section-bruteforce.v1"
BRUTE_STATUS = "PASS_EXHAUSTIVE_FINITE_FIELD_ENUMERATION_OF_EXPORTED_BLOCKS"
LIFT_SCHEMA = "elkies-k3.r17-norm12-direct-singleton-po0-bruteforce-hensel.v1"
LIFT_STATUS = "PASS_EXACT_BRUTEFORCE_SEED_HENSEL_AUDIT"
CENSUS_SCHEMA = "elkies-k3.r17-norm12-direct-singleton-twist-census.v1"
CENSUS_STATUS = "PASS_BOUNDED_HEURISTIC_DIRECT_SINGLETON_TWIST_CENSUS"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load(path: Path):
    return json.loads(path.read_text())


def brute_path(orbit: str, prime: int) -> Path:
    return GENERATED / (
        "elkies-k3-r17-norm12-11952-singleton-twist-po0-"
        f"orbit{orbit}-p{prime}-v1.json"
    )


def lift_path(orbit: str, prime: int) -> Path:
    return GENERATED / (
        "elkies-k3-r17-norm12-11952-singleton-twist-po0-"
        f"orbit{orbit}-p{prime}-hensel-v1.json"
    )


def verify_input_pins(payload, label: str) -> list[Path]:
    paths = []
    for name, expected in payload.get("inputs", {}).items():
        path = ROOT / name
        if not path.is_file():
            raise FileNotFoundError(f"{label}: missing hash-pinned input {name}")
        if digest(path) != expected:
            raise ArithmeticError(f"{label}: input digest changed for {name}")
        paths.append(path)
    return paths


def census_labels(payload) -> list[str]:
    return [str(record["label"]) for record in payload["top"]]


def obstruction_exponent(message: str) -> int:
    match = re.fullmatch(r"no lift modulo (?:p|\d+)\^(\d+)", message)
    if match is None:
        raise ArithmeticError(f"unexpected exact obstruction message: {message}")
    return int(match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    inputs = [DISCOVERY, *HOLDOUTS]
    discovery = load(DISCOVERY)
    if (
        discovery.get("schema") != CENSUS_SCHEMA
        or discovery.get("status") != CENSUS_STATUS
        or discovery.get("candidate_count") != 39147
    ):
        raise ArithmeticError("discovery census metadata changed")
    expected_labels = [f"alternate-orbit-{orbit}" for orbit, unused in PRIMARY_SHELLS]
    if census_labels(discovery)[:20] != expected_labels:
        raise ArithmeticError("primary shell list no longer equals the discovery top 20")

    holdout_records = []
    for path in HOLDOUTS:
        payload = load(path)
        if (
            payload.get("schema") != CENSUS_SCHEMA
            or payload.get("status") != CENSUS_STATUS
            or payload.get("candidate_count") != 20
            or set(census_labels(payload)) != set(expected_labels)
            or len(payload.get("usable_primes", [])) != 48
        ):
            raise ArithmeticError(f"holdout census changed: {relative(path)}")
        holdout_records.append(
            {
                "artifact": relative(path),
                "usable_prime_count": len(payload["usable_primes"]),
                "prime_minimum": min(payload["usable_primes"]),
                "prime_maximum": max(payload["usable_primes"]),
                "highest_mean_score": max(
                    float(record["mean_block_score"]) for record in payload["top"]
                ),
                "highest_weakest_block_score": max(
                    float(record["weakest_block_score"]) for record in payload["top"]
                ),
            }
        )

    shell_records = []
    total_tested = 0
    total_solutions = 0
    total_isolated_extras = 0
    obstruction_histogram = Counter()
    for discovery_rank, (orbit, prime) in enumerate(PRIMARY_SHELLS, start=1):
        path = brute_path(orbit, prime)
        payload = load(path)
        inputs.append(path)
        inputs.extend(verify_input_pins(payload, f"{orbit} p={prime} brute shell"))
        candidate = payload.get("candidate", {})
        if (
            payload.get("schema") != BRUTE_SCHEMA
            or payload.get("status") != BRUTE_STATUS
            or payload.get("prime") != prime
            or candidate.get("key") != f"alternate-orbit-{orbit}"
            or candidate.get("kind") != "direct_singleton"
            or candidate.get("chi") != 3
            or candidate.get("x_degree_bound") != 6
            or candidate.get("y_degree_bound") != 9
        ):
            raise ArithmeticError(f"{orbit} p={prime}: brute shell metadata changed")

        solutions = payload["solutions"]
        enumeration = payload["enumeration"]
        known = [int(index) for index in payload["known_section_match_indices"]]
        if (
            len(known) != 1
            or enumeration["x_polynomials_tested"] <= 0
            or enumeration["representative_sign_solutions"] != len(solutions)
            or any(index < 0 or index >= len(solutions) for index in known)
        ):
            raise ArithmeticError(f"{orbit} p={prime}: shell completeness changed")
        ranks = [int(solution["full_shell_tangent_rank"]) for solution in solutions]
        recomputed_histogram = {
            str(rank): ranks.count(rank) for rank in sorted(set(ranks))
        }
        if recomputed_histogram != payload["full_shell_tangent_rank_histogram"]:
            raise ArithmeticError(f"{orbit} p={prime}: tangent-rank histogram changed")
        isolated_extras = [
            index for index, rank in enumerate(ranks) if rank == 8 and index not in known
        ]

        obstruction_records = []
        if isolated_extras:
            hensel_path = lift_path(orbit, prime)
            hensel = load(hensel_path)
            inputs.append(hensel_path)
            inputs.extend(verify_input_pins(hensel, f"{orbit} p={prime} Hensel audit"))
            if (
                hensel.get("schema") != LIFT_SCHEMA
                or hensel.get("status") != LIFT_STATUS
                or hensel.get("prime") != prime
                or hensel.get("candidate") != candidate
                or hensel.get("inputs", {}).get(relative(path)) != digest(path)
            ):
                raise ArithmeticError(f"{orbit} p={prime}: Hensel metadata changed")
            lifts = {int(record["solution_index"]): record for record in hensel["lifts"]}
            if not set(isolated_extras).issubset(lifts):
                raise ArithmeticError(f"{orbit} p={prime}: an isolated branch was not lifted")
            for index in isolated_extras:
                lift = lifts[index]
                obstruction = lift.get("exact_local_obstruction")
                if not obstruction:
                    raise ArithmeticError(
                        f"{orbit} p={prime} solution {index}: no exact local obstruction"
                    )
                exponent = obstruction_exponent(obstruction)
                obstruction_histogram[exponent] += 1
                obstruction_records.append(
                    {
                        "solution_index": index,
                        "first_impossible_exponent": exponent,
                    }
                )

        record = {
            "discovery_rank": discovery_rank,
            "label": f"alternate-orbit-{orbit}",
            "prime": prime,
            "x_polynomials_tested": int(enumeration["x_polynomials_tested"]),
            "representative_sign_solution_count": len(solutions),
            "known_solution_index": known[0],
            "isolated_extra_solution_indices": isolated_extras,
            "exact_local_obstructions": obstruction_records,
        }
        shell_records.append(record)
        total_tested += record["x_polynomials_tested"]
        total_solutions += len(solutions)
        total_isolated_extras += len(isolated_extras)

    regression = load(KNOWN_REGRESSION)
    inputs.append(KNOWN_REGRESSION)
    inputs.extend(verify_input_pins(regression, "known-section Hensel regression"))
    regression_lift = regression.get("lifts", [{}])[0]
    reconstruction = regression_lift.get("exact_rational_reconstruction", {})
    if (
        regression.get("schema") != LIFT_SCHEMA
        or regression.get("status") != LIFT_STATUS
        or regression.get("prime") != 23
        or regression.get("candidate", {}).get("key") != "alternate-orbit-1c3d5"
        or regression.get("hensel_depth") != 512
        or regression.get("selected_solution_indices") != [4]
        or regression.get("exact_local_obstruction_count") != 0
        or regression.get("unique_hensel_lift_count") != 1
        or regression.get("exact_rational_section_count") != 1
        or len(regression.get("exact_sections", [])) != 1
        or regression_lift.get("solution_index") != 4
        or regression_lift.get("exact_local_obstruction") is not None
        or len(regression_lift.get("hensel_levels", [])) != 511
        or not all(level["compatible"] for level in regression_lift["hensel_levels"])
        or reconstruction.get("literal_curve_substitution") is not True
    ):
        raise ArithmeticError("known-section depth-512 Hensel regression changed")

    direct_inputs = [
        ROOT / "elkies-k3/scripts/bruteforce_twist_polynomial_sections_modp.cpp",
        ROOT / "elkies-k3/scripts/run_twist_polynomial_sections_bruteforce.py",
        ROOT / "elkies-k3/scripts/lift_r17_norm12_direct_singleton_po0_bruteforce.sage",
    ]
    inputs.extend(direct_inputs)
    unique_inputs = sorted(set(path.resolve() for path in inputs))
    result = {
        "schema": "elkies-k3.r17-norm12-11952-singleton-po0-top20-audit.v1",
        "status": "PASS_EXACT_BOUNDED_TOP20_PO0_NO_SECOND_TWIST_SECTION",
        "search_box": {
            "character_count": len(PRIMARY_SHELLS),
            "chi": 3,
            "intersection_with_zero_section": 0,
            "x_degree_bound": 6,
            "y_degree_bound": 9,
            "total_x_polynomials_tested_across_distinct_finite_field_shells": total_tested,
            "representative_sign_solution_count": total_solutions,
            "isolated_extra_modular_branch_count": total_isolated_extras,
            "exact_local_obstruction_count": sum(obstruction_histogram.values()),
            "first_impossible_exponent_histogram": {
                str(exponent): obstruction_histogram[exponent]
                for exponent in sorted(obstruction_histogram)
            },
            "records": shell_records,
        },
        "known_section_regression": {
            "artifact": relative(KNOWN_REGRESSION),
            "label": "alternate-orbit-1c3d5",
            "prime": 23,
            "hensel_depth": 512,
            "solution_index": 4,
            "literal_curve_substitution": True,
        },
        "heuristic_holdouts": holdout_records,
        "inputs": {relative(path): digest(path) for path in unique_inputs},
        "proof_boundary": (
            "The finite-field enumerations are exhaustive for the chi=3, P.O=0 "
            "polynomial-section boxes of the twenty discovery-ranked characters, "
            "and every isolated extra modular branch in the chosen shells has an "
            "exact obstruction at a stated prime power. Singular modular branches, "
            "the other 39127 direct singleton characters, sections with P.O>0, and "
            "non-singleton characters are not excluded. The two holdout censuses are "
            "target-selection heuristics only. Hence this result does not prove a "
            "second twist section or generic rank at least 19."
        ),
        "reproducing_command": (
            ".venv/bin/python elkies-k3/scripts/"
            "audit_r17_norm12_11952_singleton_po0_top20.py --check"
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored top-20 P.O=0 audit differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "R17SINGLETONPO0TOP20"
        f"|characters={len(PRIMARY_SHELLS)}"
        f"|tested={total_tested}"
        f"|solutions={total_solutions}"
        f"|isolated_extras={total_isolated_extras}"
        f"|obstructed={sum(obstruction_histogram.values())}"
        f"|status={result['status']}"
    )


if __name__ == "__main__":
    main()
