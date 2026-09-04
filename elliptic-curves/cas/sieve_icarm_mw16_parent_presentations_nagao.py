#!/usr/bin/env python3
"""Bounded Nagao/local prefilter on all nine MW16 atlas presentations.

The nine coordinate presentations represent five exact Jacobian fibrations.
Each presentation receives the same projective-height box and frozen local
prime blocks.  Multiple presentations are retained because an affine base
change does not preserve the bounded height box, but they remain nested search
charts rather than independent observations.

This stage evaluates no sections and performs no point search.  Its finalists
may advance only to the bounded half-lattice recovery stage.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
from time import perf_counter


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
AUDIT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_parent_presentation_audit_v1.json"
COMMON = ROOT / "elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_parent_nagao_prefilter_h300_v1.json"
DEFAULT_PRIME_BLOCKS = (
    (19, 41, 43, 61, 71, 73, 79, 83),
    (89, 107, 113, 127, 131, 137, 139, 151),
    (157, 163, 167, 173, 179, 181, 191, 193, 197),
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_sha256(value) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def parent_model(common, parent, input_path: Path):
    pencil = parent["pencil"]
    source_hash = canonical_sha256(
        {
            "parent_id": parent["parent_id"],
            "A": pencil["A_coefficients_low_to_high"],
            "B": pencil["B_coefficients_low_to_high"],
        }
    )
    return common.FamilyModel(
        source=input_path.resolve(),
        source_sha256=source_hash,
        a_coefficients=tuple(
            Fraction(value) for value in pencil["A_coefficients_low_to_high"]
        ),
        b_coefficients=tuple(
            Fraction(value) for value in pencil["B_coefficients_low_to_high"]
        ),
        a_degree=8,
        b_degree=12,
        coordinate=f"{parent['parent_id']}:lambda",
        coefficient_source_keys=(
            "pencil.A_coefficients_low_to_high",
            "pencil.B_coefficients_low_to_high",
        ),
    )


def score_parameter(common, parameter: Fraction, blocks):
    candidate = common.Candidate(parameter.numerator, parameter.denominator, max(abs(parameter.numerator), parameter.denominator))
    inverse_cache = {}
    for block in blocks:
        candidate = common.score_block(candidate, block, inverse_cache)
    return common.candidate_record(candidate)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--audit", type=Path, default=AUDIT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=300)
    parser.add_argument("--height-bucket-width", type=int, default=50)
    parser.add_argument("--keep-per-bucket", default="16,8,4")
    parser.add_argument("--finalists-per-presentation", type=int, default=12)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    keep = tuple(int(value) for value in args.keep_per_bucket.split(","))
    if args.height_bound < 1 or args.height_bucket_width < 1:
        raise SystemExit("height bounds must be positive")
    if len(keep) != len(DEFAULT_PRIME_BLOCKS) or any(value < 1 for value in keep):
        raise SystemExit("--keep-per-bucket requires three positive entries")
    if args.finalists_per_presentation < 1:
        raise SystemExit("--finalists-per-presentation must be positive")

    inputs = json.loads(args.input.read_text())
    audit = json.loads(args.audit.read_text())
    if inputs.get("status") != "PASS_EXACT_COMPLEMENT_BLIND_NINE_PARENT_INPUTS":
        raise ArithmeticError("MW16 parent inputs are not passing")
    if audit.get("status") != "PASS_EXACT_NINE_PRESENTATIONS_FIVE_FIBRATIONS":
        raise ArithmeticError("MW16 presentation audit is not passing")
    if audit["operational_rule"]["coordinate_search_trial_count"] != 9 or audit["operational_rule"]["independent_observation_count"] != 5:
        raise ArithmeticError("presentation/fibration partition changed")

    common = SourceFileLoader("mw16_nagao_common", str(COMMON)).load_module()
    presentations = []
    started = perf_counter()
    for parent in inputs["parents"]:
        parent_started = perf_counter()
        model = parent_model(common, parent, args.input)
        table_blocks, rejected = common.build_residue_tables(
            model, DEFAULT_PRIME_BLOCKS
        )
        survivors, stages = common.run_staged_sieve(
            numerator_bound=args.height_bound,
            denominator_bound=args.height_bound,
            table_blocks=table_blocks,
            keep_per_bucket=keep,
            bucket_width=args.height_bucket_width,
        )
        finalists = sorted(survivors, key=common.candidate_sort_key)[
            : args.finalists_per_presentation
        ]
        control = score_parameter(
            common, Fraction(parent["target_parameter"]), table_blocks
        )
        presentations.append(
            {
                "parent_id": parent["parent_id"],
                "curve_id": int(parent["curve_id"]),
                "priority_rank": int(parent["priority_rank"]),
                "pencil_sha256": model.source_sha256,
                "usable_prime_blocks": [
                    list(block.keys()) for block in table_blocks
                ],
                "rejected_primes": list(rejected),
                "stages": stages,
                "known_target_fibre_local_control": control,
                "finalists": [
                    common.candidate_record(candidate) for candidate in finalists
                ],
                "runtime_seconds": perf_counter() - parent_started,
            }
        )
        print(
            f"MW16NAGAO|parent={parent['parent_id']}|"
            f"population={stages[0]['population_scored']}|"
            f"survivors={len(survivors)}|finalists={len(finalists)}|status=PASS",
            flush=True,
        )

    payload = {
        "schema": "elliptic-curves.icarm-mw16-parent-nagao-prefilter.v1",
        "status": "PASS_BOUNDED_NINE_PRESENTATION_NAGAO_PREFILTER",
        "geometry": {
            "exact_fibration_count": 5,
            "coordinate_presentation_count": 9,
            "statistical_unit": "target-curve/fibration cluster",
            "presentation_role": (
                "nested coordinate-height search charts; never independent outcomes"
            ),
        },
        "nagao_contribution": {
            "formula": "((2-a_p)/(p+1-a_p))*log(p)",
            "integer_scale": common.SCORE_SCALE,
            "bad_reduction_contribution": 0,
        },
        "search": {
            "projective_height_box": args.height_bound,
            "numerator_interval": [-args.height_bound, args.height_bound],
            "denominator_interval": [1, args.height_bound],
            "primitive_pairs_only": True,
            "includes_zero_and_infinity": True,
            "height_bucket_width": args.height_bucket_width,
            "requested_prime_blocks": [
                list(block) for block in DEFAULT_PRIME_BLOCKS
            ],
            "keep_per_bucket": list(keep),
            "finalists_per_presentation": args.finalists_per_presentation,
            "total_finalist_rows_before_exact_curve_deduplication": sum(
                len(row["finalists"]) for row in presentations
            ),
        },
        "presentations": presentations,
        "next_gate": {
            "stage": "bounded_half_lattice_jump_recovery",
            "required_inputs": [
                "exact specialization of the saturated generic MW16 basis",
                "exact target-fibre minimal model",
                "state-bound recomputation of the complete intended chart order",
            ],
            "forbidden_direct_transition": "Nagao finalist -> unrestricted or giant point search",
            "after_positive_recovery": (
                "run a complete same-minimal-curve residual 2-Selmer calculation "
                "before any expensive continuation toward rank 32"
            ),
        },
        "inputs": {
            relative(path): digest(path)
            for path in (args.input, args.audit, COMMON, Path(__file__))
        },
        "runtime_seconds": perf_counter() - started,
        "claim_boundary": [
            "Nagao scores and staged survival are scheduling heuristics, not rank evidence.",
            "No section is evaluated and no point search or Selmer calculation occurs in this stage.",
            "A low score cannot veto an unsearched fibre.",
            "The known target scores are local positive-control measurements, not independent observations across repeated presentations.",
            "Finalists are authorized only for bounded half-lattice recovery, never direct expensive point search.",
        ],
        "reproducing_command": (
            "python3 elliptic-curves/cas/sieve_icarm_mw16_parent_presentations_nagao.py"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        stored = json.loads(args.output.read_text())
        for document in (stored, payload):
            document.pop("runtime_seconds", None)
            for row in document["presentations"]:
                row.pop("runtime_seconds", None)
                for stage in row["stages"]:
                    stage.pop("runtime_seconds", None)
                    stage.pop("parameters_per_second", None)
        if stored != payload:
            raise ArithmeticError("stored Nagao prefilter differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        f"MW16NAGAO|presentations=9|fibrations=5|"
        f"finalists={payload['search']['total_finalist_rows_before_exact_curve_deduplication']}|"
        f"output={relative(args.output)}|status={payload['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
