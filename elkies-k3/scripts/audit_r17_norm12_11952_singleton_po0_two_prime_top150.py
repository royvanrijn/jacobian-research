#!/usr/bin/env python3
"""Audit two distinct exhaustive P.O=0 shells on the top 150 singleton twists."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from run_r17_norm12_11952_singleton_po0_top150 import (  # noqa: E402
    verify_manifest,
)


GENERATED = ROOT / "artifacts/generated-results"
PRIMARY = GENERATED / (
    "elkies-k3-r17-norm12-11952-singleton-twist-po0-top150-campaign-v1.json"
)
SECONDARY = GENERATED / (
    "elkies-k3-r17-norm12-11952-singleton-twist-po0-"
    "top150-distinct-second-campaign-v1.json"
)
KNOWN_REGRESSION = GENERATED / (
    "elkies-k3-r17-norm12-11952-singleton-twist-po0-"
    "orbit1c3d5-known-p23-hensel512-v1.json"
)
OUTPUT = GENERATED / (
    "elkies-k3-r17-norm12-11952-singleton-twist-po0-"
    "two-prime-top150-audit-v1.json"
)
CAMPAIGN_STATUS = "PASS_BOUNDED_TOPN_ALL_ISOLATED_EXTRA_BRANCHES_LOCALLY_OBSTRUCTED"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load(path: Path):
    return json.loads(path.read_text())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    # This rechecks every hash pin, census rank, exact shell count, tangent-rank
    # classification, Hensel obstruction, and aggregate in both campaigns.
    verify_manifest(PRIMARY)
    verify_manifest(SECONDARY)
    primary = load(PRIMARY)
    secondary = load(SECONDARY)
    if (
        primary["status"] != CAMPAIGN_STATUS
        or secondary["status"] != CAMPAIGN_STATUS
        or primary["candidate_limit"] != 150
        or secondary["candidate_limit"] != 150
        or secondary.get("excluded_prime_manifest") != relative(PRIMARY)
        or secondary["inputs"].get(relative(PRIMARY)) != digest(PRIMARY)
    ):
        raise ArithmeticError("two-prime campaign metadata changed")

    pair_histogram = Counter()
    both_with_extras = 0
    neither_with_extras = 0
    records = []
    for first, second in zip(primary["records"], secondary["records"], strict=True):
        if (
            first["discovery_rank"] != second["discovery_rank"]
            or first["label"] != second["label"]
            or first["prime"] == second["prime"]
            or second["excluded_prime"] != first["prime"]
        ):
            raise ArithmeticError("campaigns do not give two distinct shells per character")
        first_extras = len(first["isolated_extra_solution_indices"])
        second_extras = len(second["isolated_extra_solution_indices"])
        pair_histogram[(int(first["prime"]), int(second["prime"]))] += 1
        both_with_extras += bool(first_extras and second_extras)
        neither_with_extras += not first_extras and not second_extras
        records.append(
            {
                "discovery_rank": first["discovery_rank"],
                "label": first["label"],
                "first_prime": first["prime"],
                "second_prime": second["prime"],
                "first_isolated_extra_count": first_extras,
                "second_isolated_extra_count": second_extras,
            }
        )

    regression = load(KNOWN_REGRESSION)
    exact = regression.get("exact_sections", [])
    if (
        regression.get("status") != "PASS_EXACT_BRUTEFORCE_SEED_HENSEL_AUDIT"
        or regression.get("hensel_depth") != 512
        or regression.get("exact_rational_section_count") != 1
        or len(exact) != 1
        or exact[0].get("literal_curve_substitution") is not True
    ):
        raise ArithmeticError("depth-512 known-section positive control changed")

    first_summary = primary["summary"]
    second_summary = secondary["summary"]
    obstruction_histogram = Counter(
        {
            int(exponent): int(count)
            for exponent, count in first_summary[
                "first_impossible_exponent_histogram"
            ].items()
        }
    )
    obstruction_histogram.update(
        {
            int(exponent): int(count)
            for exponent, count in second_summary[
                "first_impossible_exponent_histogram"
            ].items()
        }
    )
    result = {
        "schema": "elkies-k3.r17-norm12-11952-singleton-po0-two-prime-top150.v1",
        "status": "PASS_EXACT_BOUNDED_TWO_PRIME_TOP150_NO_ISOLATED_SECOND_SECTION",
        "character_count": 150,
        "complete_shell_count": 300,
        "distinct_prime_pair_count": len(pair_histogram),
        "prime_pair_histogram": {
            f"{first},{second}": pair_histogram[(first, second)]
            for first, second in sorted(pair_histogram)
        },
        "x_polynomials_tested_across_distinct_shells": (
            first_summary["x_polynomials_tested_across_distinct_shells"]
            + second_summary["x_polynomials_tested_across_distinct_shells"]
        ),
        "representative_sign_solution_count": (
            first_summary["representative_sign_solution_count"]
            + second_summary["representative_sign_solution_count"]
        ),
        "isolated_extra_modular_branch_count": (
            first_summary["isolated_extra_branch_count"]
            + second_summary["isolated_extra_branch_count"]
        ),
        "exact_local_obstruction_count": (
            first_summary["exact_local_obstruction_count"]
            + second_summary["exact_local_obstruction_count"]
        ),
        "first_impossible_exponent_histogram": {
            str(exponent): obstruction_histogram[exponent]
            for exponent in sorted(obstruction_histogram)
        },
        "characters_with_isolated_extras_at_both_primes": both_with_extras,
        "characters_with_no_isolated_extra_at_either_prime": neither_with_extras,
        "unresolved_isolated_branch_count": 0,
        "exact_new_section_count": 0,
        "known_section_positive_control": {
            "artifact": relative(KNOWN_REGRESSION),
            "label": "alternate-orbit-1c3d5",
            "prime": 23,
            "hensel_depth": 512,
            "literal_curve_substitution": True,
        },
        "records": records,
        "inputs": {
            relative(path): digest(path)
            for path in (
                Path(__file__).resolve(),
                SCRIPTS / "run_r17_norm12_11952_singleton_po0_top150.py",
                PRIMARY,
                SECONDARY,
                KNOWN_REGRESSION,
            )
        },
        "proof_boundary": (
            "For each of the 150 ranked characters, the two finite-field primes "
            "are distinct and each exported P.O=0 polynomial-section shell is "
            "complete. Every isolated extra point in all 300 shells has an exact "
            "prime-power nonlifting obstruction. A characteristic-zero section "
            "could still reduce into the known-section or singular locus at both "
            "chosen primes. Sections with P.O>0, the remaining 38997 singleton "
            "characters, and non-singleton characters are outside this result."
        ),
        "reproducing_command": (
            ".venv/bin/python elkies-k3/scripts/"
            "audit_r17_norm12_11952_singleton_po0_two_prime_top150.py --check"
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if args.check:
        if not output.is_file() or output.read_text() != serialized:
            raise ArithmeticError("stored two-prime top-150 audit differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "R17SINGLETONPO0TWO_PRIME_TOP150"
        f"|characters={result['character_count']}"
        f"|shells={result['complete_shell_count']}"
        f"|tested={result['x_polynomials_tested_across_distinct_shells']}"
        f"|isolated_extras={result['isolated_extra_modular_branch_count']}"
        f"|obstructed={result['exact_local_obstruction_count']}"
        f"|status={result['status']}"
    )


if __name__ == "__main__":
    main()
