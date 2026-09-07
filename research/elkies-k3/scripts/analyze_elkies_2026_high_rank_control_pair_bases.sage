#!/usr/bin/env sage-python
"""Construct paired bases selected by the four high-rank control fibres.

status: ACTIVE_PROOF
claim: exact control-selected pair catalogue and finite-quotient incidence table
inputs: published R17 model/sections, complete bisection batch, positive controls
outputs: artifacts/generated-results/elkies-2026-high-rank-control-pair-bases.json
supersedes: none

For each published control ``t0`` this script computes

    S(t0) = {i : q_i(t0) is a rational square},

constructs every pair in ``binomial(S(t0),2)``, materializes its two lifted
points on the control fibre, and measures their incidence with the fixed
generic-17 plus public exceptional complement in the same exact mod-2 finite
quotient ensemble.  A zero increment is an incidence statement in that finite
quotient, not a proof of rational dependence or an exact-rank upper bound.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

from sage.all import QQ


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from elliptic_candidate_record import (  # noqa: E402
    finite_quotient_signature,
    matrix_rank_and_pivots_mod_prime,
    primes_up_to,
    source_point_to_target,
)
from elkies_rank25 import POINTS as RANK25_POINTS  # noqa: E402
from elkies_rank26 import POINTS as RANK26_POINTS  # noqa: E402
from elkies_rank27 import POINTS as RANK27_POINTS  # noqa: E402
from elkies_rank28 import POINTS as RANK28_POINTS  # noqa: E402


MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-high-rank-control-pair-bases.json"
PUBLIC_POINTS = {
    25: RANK25_POINTS,
    26: RANK26_POINTS,
    27: RANK27_POINTS,
    28: RANK28_POINTS,
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def to_fraction(value) -> Fraction:
    value = QQ(value)
    return Fraction(int(value.numerator()), int(value.denominator()))


def evaluate(coefficients, value):
    answer = QQ(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + QQ(coefficient)
    return answer


def lifted_affine_point(record, t_value, root):
    lifted = record["lifted_section"]
    return (
        evaluate(lifted["x0_coefficients"], t_value)
        + evaluate(lifted["x1_coefficients"], t_value) * root,
        evaluate(lifted["y0_coefficients"], t_value)
        + evaluate(lifted["y1_coefficients"], t_value) * root,
    )


def rank_from_rows(rows, column_count):
    return int(matrix_rank_and_pivots_mod_prime(rows, column_count, 2)[0])


def incidence_ranks(short_model, points, primes, cuts):
    rows = []
    signatures = []
    for prime in primes:
        try:
            signature = finite_quotient_signature(short_model, points, int(prime), 2)
        except ValueError:
            # A newly materialized point may have a denominator at a prime
            # that was good for the pinned public basis. Such a prime supplies
            # no reduction row for this enlarged list and is omitted exactly.
            continue
        signatures.append(signature)
        rows.extend(signature["rows"])
    if not signatures:
        raise ArithmeticError("no finite-quotient prime survived point denominators")
    ranks = {label: rank_from_rows([row[:cut] for row in rows], cut) for label, cut in cuts.items()}
    return ranks, signatures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--sections", type=Path, default=SECTIONS)
    parser.add_argument("--bisections", type=Path, default=BISECTIONS)
    parser.add_argument("--controls", type=Path, default=CONTROLS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bisection_document = json.loads(args.bisections.read_text())
    control_document = json.loads(args.controls.read_text())
    if len(bisection_document.get("bisections", ())) != 39120:
        raise ArithmeticError("the complete bisection catalogue changed")
    if control_document.get("status") != "PASS_EXACT_ELKIES_2026_HIGH_RANK_POSITIVE_CONTROLS":
        raise ArithmeticError("the four exact positive controls are missing")
    records = bisection_document["bisections"]
    by_mask = {int(record["lattice_orbit_mask"]): record for record in records}
    data = load_q12o5867_data(args.model, args.sections)

    rows = []
    pair_controls = defaultdict(list)
    mask_occurrences = Counter()
    for control in control_document["fibres"]:
        numerator, denominator = map(int, control["projective_parameter"])
        t_value = QQ(numerator) / QQ(denominator)
        split = []
        roots = {}
        for record in records:
            q_value = evaluate(record["residual_chord"]["q_coefficients"], t_value)
            if q_value.is_square():
                mask = int(record["lattice_orbit_mask"])
                split.append(mask)
                roots[mask] = q_value.sqrt()
                mask_occurrences[mask] += 1
        split.sort()

        specialization = evaluate_projective_specialization(data, numerator, denominator)
        minimal_model, minimal_change, _metadata = global_minimal_model_with_change(
            specialization.model
        )
        expected_minimal = tuple(Fraction(value) for value in control["minimal_model"])
        if minimal_model != expected_minimal:
            raise ArithmeticError("a control minimal model changed")
        short_model, short_change = short_certificate_model(minimal_model)
        generic_short = tuple(
            source_point_to_target(
                source_point_to_target(point, minimal_change), short_change
            )
            for point in specialization.points
        )
        published_rank = int(control["published_rank_lower_bound"])
        selected_indices = [
            int(value) - 1
            for value in control["public_positive_control"][
                "selected_public_point_indices_one_based"
            ]
        ]
        public_complement = tuple(PUBLIC_POINTS[published_rank][index] for index in selected_indices)
        public_short = tuple(source_point_to_target(point, short_change) for point in public_complement)
        baseline = generic_short + public_short
        certificate = control["public_positive_control"][
            "combined_generic_plus_complement_independence"
        ]
        primes = tuple(
            dict.fromkeys(
                [int(value) for value in certificate["certificate_primes"]]
                + [value for value in primes_up_to(300) if value != 2]
            )
        )
        if len(baseline) != published_rank:
            raise ArithmeticError("the public complement size changed")

        pair_rows = []
        for left, right in combinations(split, 2):
            key = f"{left}:{right}"
            pair_controls[key].append(control["label"])
            lifted_short = []
            for mask in (left, right):
                affine_x, affine_y = lifted_affine_point(by_mask[mask], t_value, roots[mask])
                projective_point = (
                    to_fraction(affine_x * denominator**4),
                    to_fraction(affine_y * denominator**6),
                )
                minimal_point = source_point_to_target(projective_point, minimal_change)
                lifted_short.append(source_point_to_target(minimal_point, short_change))
            points = baseline + tuple(lifted_short)
            cuts = {
                "generic_17": 17,
                "generic_plus_public": published_rank,
                "generic_plus_public_plus_first_pair_point": published_rank + 1,
                "generic_plus_public_plus_pair": published_rank + 2,
            }
            ranks, signatures = incidence_ranks(short_model, points, primes, cuts)
            if ranks["generic_17"] != 17 or ranks["generic_plus_public"] != published_rank:
                raise ArithmeticError("the pinned control incidence baseline failed")
            generic_pair_points = generic_short + tuple(lifted_short)
            generic_pair_ranks, _ = incidence_ranks(
                short_model,
                generic_pair_points,
                primes,
                {"generic_17": 17, "generic_plus_pair": 19},
            )
            pair_rows.append(
                {
                    "pair_key": key,
                    "orbit_masks": [left, right],
                    "orbit_hex": [f"0x{left:05x}", f"0x{right:05x}"],
                    "base_point_t_u_v": [
                        rational_text(t_value),
                        rational_text(roots[left]),
                        rational_text(roots[right]),
                    ],
                    "lifted_short_points": [
                        [str(point[0]), str(point[1])] for point in lifted_short
                    ],
                    "finite_quotient_primes": [int(signature["prime"]) for signature in signatures],
                    "finite_quotient_signature_sha256": sha256(
                        json.dumps(signatures, sort_keys=True, separators=(",", ":")).encode()
                    ).hexdigest(),
                    "incidence_ranks": {
                        **ranks,
                        **generic_pair_ranks,
                        "pair_increment_beyond_generic_17": generic_pair_ranks["generic_plus_pair"] - 17,
                        "pair_increment_beyond_public_control": (
                            ranks["generic_plus_public_plus_pair"] - published_rank
                        ),
                    },
                    "interpretation": (
                        "These are ranks in the pinned exact mod-2 finite-quotient ensemble. "
                        "A zero public-control increment does not prove rational dependence."
                    ),
                }
            )
        rows.append(
            {
                "label": control["label"],
                "parameter": rational_text(t_value),
                "published_rank_lower_bound": published_rank,
                "S_size": len(split),
                "S_masks": split,
                "S_hex": [f"0x{mask:05x}" for mask in split],
                "pair_count": len(pair_rows),
                "pairs": pair_rows,
            }
        )

    repeated_pairs = {
        key: labels for key, labels in pair_controls.items() if len(labels) >= 2
    }
    repeated_masks = {
        str(mask): count for mask, count in mask_occurrences.items() if count >= 2
    }
    result = {
        "schema": "elkies-k3.elkies-2026-high-rank-control-pair-bases.v1",
        "status": "PASS_EXACT_HIGH_RANK_CONTROL_PAIR_BASE_INCIDENCE",
        "inputs": {
            display_path(Path(__file__).resolve()): digest(Path(__file__).resolve()),
            display_path(args.model): digest(args.model),
            display_path(args.sections): digest(args.sections),
            display_path(args.bisections): digest(args.bisections),
            display_path(args.controls): digest(args.controls),
        },
        "controls": rows,
        "summary": {
            "S_sizes": [row["S_size"] for row in rows],
            "pair_counts": [row["pair_count"] for row in rows],
            "total_control_selected_pairs": sum(row["pair_count"] for row in rows),
            "masks_occurring_at_two_or_more_controls": repeated_masks,
            "pairs_occurring_at_two_or_more_controls": repeated_pairs,
            "repeated_pair_count": len(repeated_pairs),
        },
        "proof_boundary": (
            "The S(t0) square tests, lifted points, and finite-quotient matrices are exact. "
            "The incidence ranks do not give rank upper bounds: a point invisible modulo 2 "
            "at these primes may still be independent over Q. In this four-control batch no "
            "mask, and hence no pair, occurs at two controls."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        "ELKIES2026CONTROLPAIRS|"
        f"S_sizes={','.join(str(row['S_size']) for row in rows)}|"
        f"pairs={sum(row['pair_count'] for row in rows)}|repeated_pairs={len(repeated_pairs)}|"
        f"status={result['status']}|output={display_path(args.output)}"
    )


if __name__ == "__main__":
    main()
