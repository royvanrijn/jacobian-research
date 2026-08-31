#!/usr/bin/env python3
"""Build an unconditional, descent-first Magma job for Elkies's rank-28 fibre.

The generated job computes the complete 2-Selmer group before constructing a
single residual two-cover.  The repository's exact finite-reduction
certificate proves that the seventeen specialized generic sections have
independent Kummer images, so their quotient dimension is obtained by
subtracting 17 from the complete Selmer dimension.  A residual dimension below
15 rejects rank 32 and exits.  On a pass, the eleven certified public
complement directions are removed as well, so ``TwoDescent`` materializes only
the classes unexplained by the full known rank-28 subgroup.  This program
contains no point search.

Magma is instructed to use ``Bound := -1`` and no GRH class-group bound.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
DEFAULT_CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
EXPECTED_CONTROL_SCHEMA = "elliptic-curves.elkies-2026-high-rank-positive-controls.v2"
EXPECTED_CONTROL_STATUS = "PASS_EXACT_ELKIES_2026_HIGH_RANK_POSITIVE_CONTROLS"
PARAMETER = (-9529, 5471)
GENERIC_RANK = 17
KNOWN_QUOTIENT_GAIN = 11
TARGET_RANK = 32
REQUIRED_RESIDUAL_DIMENSION = TARGET_RANK - GENERIC_RANK
KNOWN_RANK28_LOWER_BOUND = GENERIC_RANK + KNOWN_QUOTIENT_GAIN
REQUIRED_DIRECTIONS_BEYOND_KNOWN_RANK28 = TARGET_RANK - KNOWN_RANK28_LOWER_BOUND

sys.path[:0] = [str(ELLIPTIC_ROOT), str(CAS)]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from elliptic_candidate_record import (  # noqa: E402
    is_on_weierstrass_curve,
    source_point_to_target,
    verify_finite_quotient_certificate,
)
from elkies_rank28 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS as PUBLIC_MODEL,
    POINTS as PUBLIC_POINTS,
)


@dataclass(frozen=True)
class RelativeDescentInput:
    model: tuple[Fraction, ...]
    generic_points: tuple[tuple[Fraction, Fraction], ...]
    public_complement: tuple[tuple[Fraction, Fraction], ...]
    controls_sha256: str
    generic_point_sequence_sha256: str
    combined_point_sequence_sha256: str


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _pinned_hash(inputs: dict[str, str], suffix: str) -> str:
    matches = [value for key, value in inputs.items() if key.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(f"positive-control artifact does not uniquely pin {suffix}")
    return matches[0]


def _rank28_row(controls: dict[str, Any]) -> dict[str, Any]:
    matches = [row for row in controls["fibres"] if row["parameter"] == "-9529/5471"]
    if len(matches) != 1:
        raise ValueError("positive-control artifact has no unique rank-28 fibre")
    return matches[0]


def load_relative_input(
    controls_path: Path = DEFAULT_CONTROLS,
) -> RelativeDescentInput:
    """Reconstruct and replay every exact input needed by the Magma job."""

    controls = json.loads(controls_path.read_text())
    if controls.get("schema") != EXPECTED_CONTROL_SCHEMA:
        raise ValueError("unexpected positive-control schema")
    if controls.get("status") != EXPECTED_CONTROL_STATUS:
        raise ValueError("the exact positive-control suite is not passing")
    inputs = controls.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("positive-control artifact has no pinned inputs")
    current_hashes = {
        "elkies_2026_published_r17_model.json": file_sha256(MODEL),
        "elkies_2026_published_r17_sections.json": file_sha256(SECTIONS),
        "elliptic-curves/cas/elkies_rank28.py": file_sha256(CAS / "elkies_rank28.py"),
    }
    for suffix, current_hash in current_hashes.items():
        if _pinned_hash(inputs, suffix) != current_hash:
            raise ValueError(f"positive-control input changed: {suffix}")

    data = load_q12o5867_data(MODEL, SECTIONS)
    specialization = evaluate_projective_specialization(data, *PARAMETER)
    minimal_model, minimal_change, _ = global_minimal_model_with_change(
        specialization.model
    )
    if minimal_model != tuple(PUBLIC_MODEL):
        raise ArithmeticError("compact rank-28 specialization missed the public model")
    generic_points = tuple(
        source_point_to_target(point, minimal_change) for point in specialization.points
    )
    if len(generic_points) != GENERIC_RANK or any(
        not is_on_weierstrass_curve(minimal_model, point) for point in generic_points
    ):
        raise ArithmeticError("the transported generic subgroup failed exact replay")

    rank28 = _rank28_row(controls)
    if rank28["minimal_model"] != [str(value) for value in minimal_model]:
        raise ArithmeticError("the control artifact is bound to another minimal model")
    if rank28["locally_certified_rank_lower_bound"] != 28:
        raise ArithmeticError("the public rank-28 lower bound is no longer certified")

    short_model, short_change = short_certificate_model(minimal_model)
    generic_short = tuple(
        source_point_to_target(point, short_change) for point in generic_points
    )
    generic_certificate = rank28["generic_sections"]["finite_quotient_independence"]
    verify_finite_quotient_certificate(
        short_model, generic_short, generic_certificate
    )
    if (
        generic_certificate.get("relation_prime") != 2
        or generic_certificate.get("combined_rank_over_relation_field") != GENERIC_RANK
        or generic_certificate.get("certified_independent") is not True
    ):
        raise ArithmeticError("the generic Kummer image is not certified 17-dimensional")

    positive_control = rank28["public_positive_control"]
    indices = tuple(
        int(index) - 1
        for index in positive_control["selected_public_point_indices_one_based"]
    )
    if len(indices) != KNOWN_QUOTIENT_GAIN or len(set(indices)) != len(indices):
        raise ArithmeticError("the public complement is not eleven distinct directions")
    public_complement = tuple(PUBLIC_POINTS[index] for index in indices)
    combined_short = generic_short + tuple(
        source_point_to_target(point, short_change) for point in public_complement
    )
    combined_certificate = positive_control[
        "combined_generic_plus_complement_independence"
    ]
    verify_finite_quotient_certificate(
        short_model, combined_short, combined_certificate
    )
    if (
        positive_control["quotient_gain_beyond_generic_rank_17"]
        != KNOWN_QUOTIENT_GAIN
        or combined_certificate.get("combined_rank_over_relation_field") != 28
        or combined_certificate.get("certified_independent") is not True
    ):
        raise ArithmeticError("the public eleven-direction quotient gain failed replay")

    return RelativeDescentInput(
        model=tuple(minimal_model),
        generic_points=generic_points,
        public_complement=public_complement,
        controls_sha256=file_sha256(controls_path),
        generic_point_sequence_sha256=str(
            generic_certificate["point_sequence_sha256"]
        ),
        combined_point_sequence_sha256=str(
            combined_certificate["point_sequence_sha256"]
        ),
    )


def magma_q(value: Fraction | int | str) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def _magma_points(
    name: str, points: Sequence[tuple[Fraction, Fraction]]
) -> str:
    rows = ",\n    ".join(
        f"E![{magma_q(x)}, {magma_q(y)}, 1]" for x, y in points
    )
    return f"{name} := [\n    {rows}\n];"


def build_magma(source: RelativeDescentInput) -> str:
    model = ", ".join(magma_q(value) for value in source.model)
    generic = _magma_points("generic", source.generic_points)
    public_complement = _magma_points(
        "public_complement", source.public_complement
    )
    return f'''// Generated by build_elkies_2026_rank28_relative_descent_magma.py
// Protocol ELKIESR28REL v1. No SetClassGroupBounds("GRH") is used.
// Complete basis-level Selmer computation precedes all residual-cover work.

Q := Rationals();
E := EllipticCurve([{model}]);
{generic}
{public_complement}

printf "ELKIESR28REL|version=1|stage=input|magma=%o|parameter=-9529/5471|generic=17|known_quotient_floor=11|known_rank_lower_bound=28|target_rank=32|required_beyond_known_rank28=4|controls_sha256={source.controls_sha256}|generic_sha256={source.generic_point_sequence_sha256}|combined_sha256={source.combined_point_sequence_sha256}\\n",
    GetVersion();
assert #generic eq {GENERIC_RANK};
assert #public_complement eq {KNOWN_QUOTIENT_GAIN};
assert &and[P in E : P in generic];
assert &and[P in E : P in public_complement];

// The input builder has independently replayed a full-rank finite-reduction
// certificate for these seventeen Kummer images and for their eleven public
// complements.  TwoSelmerGroup now supplies the missing complete class-group
// and local-solubility layer. Bound=-1 requests unconditional class groups.
T2, T2map := TwoTorsionSubgroup(E);
assert #T2 eq 1;
printf "ELKIESR28REL|stage=two_selmer|status=start|bound=-1\\n";
started := Realtime();
S2, AtoS := TwoSelmerGroup(E : Bound := -1);
seconds := Realtime(started);
selmer_invariants := Invariants(S2);
assert &and[n eq 2 : n in selmer_invariants];
total_selmer_dim := #selmer_invariants;
residual_dim := total_selmer_dim - {GENERIC_RANK};
unexplained_dim := total_selmer_dim - {KNOWN_RANK28_LOWER_BOUND};
assert residual_dim ge {KNOWN_QUOTIENT_GAIN};
assert unexplained_dim ge 0;
printf "ELKIESR28REL|stage=two_selmer|status=complete|seconds=%o|invariants=%o|total_selmer_dim=%o|residual_dim=%o|required_residual_dim={REQUIRED_RESIDUAL_DIMENSION}|unexplained_dim=%o|required_unexplained_dim={REQUIRED_DIRECTIONS_BEYOND_KNOWN_RANK28}\\n",
    seconds, selmer_invariants, total_selmer_dim, residual_dim, unexplained_dim;

if unexplained_dim lt {REQUIRED_DIRECTIONS_BEYOND_KNOWN_RANK28} then
    printf "ELKIESR28REL|classification=REJECT_RANK32_BY_RESIDUAL_2_SELMER|total_selmer_dim=%o|residual_dim=%o|unexplained_dim=%o|expensive_search_authorized=false\\n",
        total_selmer_dim, residual_dim, unexplained_dim;
    quit;
end if;

printf "ELKIESR28REL|classification=PASS_RANK32_RESIDUAL_2_SELMER_GATE|total_selmer_dim=%o|residual_dim=%o|unexplained_dim=%o|expensive_search_authorized=true\\n",
    total_selmer_dim, residual_dim, unexplained_dim;

// This stage is unreachable until the exact dimension gate passes.  It
// quotients by all 28 certified points and materializes only the nonzero
// classes still unexplained by the known rank-28 subgroup.  It never searches
// any cover for rational points.
printf "ELKIESR28REL|stage=unexplained_two_covers|status=start|remove_known_rank=28\\n";
cover_started := Realtime();
covers := TwoDescent(
    E :
    RemoveGens := SequenceToSet(generic cat public_complement),
    RemoveTorsion := true,
    WithMaps := false
);
cover_seconds := Realtime(cover_started);
assert #covers + 1 eq 2^unexplained_dim;
printf "ELKIESR28REL|stage=unexplained_two_covers|status=complete|seconds=%o|nonzero_classes=%o|unexplained_dim=%o|removed_known_rank=28\\n",
    cover_seconds, #covers, unexplained_dim;
'''


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls", type=Path, default=DEFAULT_CONTROLS)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = load_relative_input(args.controls)
    program = build_magma(source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(program)
    print(f"WROTE {args.output}")
    print(f"PROGRAM_SHA256={sha256(program.encode()).hexdigest()}")
    print("PROTOCOL=ELKIESR28REL_v1")
    print("POINT_SEARCH=false")


if __name__ == "__main__":
    main()
