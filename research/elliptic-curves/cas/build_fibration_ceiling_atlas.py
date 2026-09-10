#!/usr/bin/env python3
"""Bind the generic-rank ceiling to the exact frame-census state.

This is a production preflight, not a fibration or point search.  A surface
may enter a rootless-neighbour construction queue only after both its Picard
ceiling and its rootless J2 frame census are certified.  In particular, an
existing rootless MW17 equation proves attainability of the ceiling but does
not by itself enumerate the other ceiling-attaining lattice types.
"""

# <!-- status-consumer: EC-DET1092-NISHIYAMA-AUXILIARY-20260910 dbfaeafbd9c1fbcc -->

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
X948_J2 = GENERATED / "elkies-k3-rootless-j2-niemeier-first.json"
X948_J1 = GENERATED / "elkies-k3-rootless-j1-uniform-bound-v1.json"
X1092_PICARD = GENERATED / "elliptic-curves/curve302_parent_geometric_picard19_v1.json"
X1092_PARENT = GENERATED / "elliptic-curves/curve302_recovered_mw17_parent_v1.json"
X1092_AUXILIARY = GENERATED / "elliptic-curves/det1092_nishiyama_auxiliary_v1.json"
DEFAULT_OUTPUT = GENERATED / "elliptic-curves/fibration-generic-rank-ceiling-v1.json"


def read(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: dict) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def build() -> dict:
    x948_j2, x948_j1 = read(X948_J2), read(X948_J1)
    x1092_picard, x1092_parent = read(X1092_PICARD), read(X1092_PARENT)
    x1092_auxiliary = read(X1092_AUXILIARY)

    assert x948_j2["status"] == "PASS_COMPLETE_ROOTLESS_J2_CLASSIFICATION"
    x948_classes = x948_j2["rootless_classes"]
    assert len(x948_classes) == 2
    assert all(row["determinant"] == 948 and row["minimum"] == 4 for row in x948_classes)
    assert x948_j1["status"] == "PASS_EXACT_ROOTLESS_J1_UNIFORM_BOUND_NOT_CLASSIFICATION"
    assert x948_j1["conclusion"]["complete_rootless_j2_class_count"] == 2

    assert x1092_picard["status"] == "PASS"
    assert x1092_picard["geometric_Picard_rank"] == 19
    assert x1092_picard["geometric_generic_MW_rank"] == 17
    # The parent packet is deliberately an explicit-input packet.  Its exact
    # rank and geometric closure are certified by X1092_PICARD above.
    assert x1092_parent["status"] == "EXPLICIT_MODEL_AND_BASIS_INPUTS"
    assert x1092_parent["height_determinant"] == "1092"
    assert x1092_auxiliary["status"] == "PASS_EXACT_AUXILIARY_AND_UNIMODULAR_GLUE"
    assert x1092_auxiliary["unimodular_glue"]["complement_isometric_to_recovered_frame"]

    return {
        "schema": "elliptic-curves.fibration-generic-rank-ceiling.v1",
        "status": "PASS_EXACT_CEILINGS_X948_J2_COMPLETE_X1092_J2_PENDING",
        "inputs": {relative(path): digest(path) for path in (X948_J2, X948_J1, X1092_PICARD, X1092_PARENT, X1092_AUXILIARY)},
        "shioda_tate": {
            "formula": "rank(MW) = rho(Xbar) - 2 - rank(reducible-fibre root lattice)",
            "consequence": "A Picard-rank-19 Jacobian K3 has generic MW rank at most 17; equality requires a rootless fibration.",
        },
        "surfaces": [
            {
                "surface": "X948",
                "geometric_picard_rank": 19,
                "generic_MW_ceiling": 17,
                "ceiling_attaining_root_rank": 0,
                "j2_frame_census": "COMPLETE",
                "rootless_j2_lattice_types": [
                    {
                        "class_index": row["class_index"],
                        "gram_sha256": row["gram_sha256"],
                        "height_determinant": row["determinant"],
                        "minimum": row["minimum"],
                        "published_R17": row["matches_published_R17"],
                        "alternate_Q80": row["matches_alternate_Q80"],
                    }
                    for row in x948_classes
                ],
                "j1_orbit_count": {
                    "lower_bound": x948_j1["conclusion"]["rootless_j1_class_count_lower_bound"],
                    "upper_bound": x948_j1["conclusion"]["rootless_j1_class_count_upper_bound"],
                    "boundary": "This is a finite J1 bound, not an exact surface-automorphism classification.",
                },
                "rootless_neighbour_equation_search": "ADMISSIBLE_AFTER_EXACT_MARKED_U_AND_EQUATION_GATES",
            },
            {
                "surface": "X1092",
                "geometric_picard_rank": 19,
                "generic_MW_ceiling": 17,
                "ceiling_attaining_root_rank": 0,
                "realized_rootless_j2_lattice_type": {
                    "height_determinant": int(x1092_parent["height_determinant"]),
                    "height_gram": x1092_parent["generic_height_gram"],
                    "reducible_fibre_root_rank": 0,
                    "meaning": "The recovered curve302 parent attains MW17, so the ceiling is attained.",
                },
                "j2_frame_census": "PENDING",
                "nishiyama_auxiliary": {
                    "status": x1092_auxiliary["status"],
                    "height_determinant": x1092_auxiliary["auxiliary"]["determinant"],
                    "discriminant_generator_q": x1092_auxiliary["auxiliary"]["discriminant_generator_q"],
                },
                "rootless_j2_lattice_types": "UNKNOWN_BEYOND_THE_REALIZED_DETERMINANT1092_TYPE",
                "rootless_neighbour_equation_search": "BLOCKED_PENDING_COMPLETE_ROOTLESS_J2_CENSUS",
            },
        ],
        "production_gate": {
            "forbidden": "Launch no rootless-neighbour equation or specialization search from a surface whose j2_frame_census is not COMPLETE.",
            "next_required_proof": "For X1092, enumerate every primitive embedding of the certified rank-seven Nishiyama auxiliary through all Niemeier lattices, and deduplicate rootless rank-17 complements by integral isometry.",
            "MW18_boundary": "Neither active surface can support MW18. A MW18 parent requires a different K3 with geometric Picard rank at least 20 and its own arithmetic descent and fibration-realization gates.",
        },
        "claim_boundary": "This records exact generic-rank ceilings and the present completeness of their frame classifications. It neither classifies X1092 rootless frames nor constructs any new fibration, equation, rational point, specialization rank, or conductor record.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = canonical(build())
    if args.check:
        assert args.output.read_bytes() == payload
        print(f"PASS byte-identical {relative(args.output)}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(payload)
    print(f"WROTE {relative(args.output)}")


if __name__ == "__main__":
    main()
