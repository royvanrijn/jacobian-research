#!/usr/bin/env sage-python
"""Canonicalize the exact N(3A8) fixed-coordinate shell scan.

Quotient every accepted coordinate seed by the complete order-twelve
glue-preserving residual group. Every row-module equality, stabilizer element,
and canonical representative is computed over Z by Hermite normal form.
Compute literal complement stabilizers, exact ternary gates, and
(T,NS)-first surface/auxiliary/frame deduplication.
"""

from __future__ import annotations

import argparse
import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
RESIDUAL = (
    ROOT / "artifacts/generated-results/elkies-k3-3a8-residual-group-v1.json"
)
PROBE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-3a8-fixed-coordinate-shell-probe-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-3a8-fixed-coordinate-shells-v1.json"
)
SHARED_CANON_SOURCE = (
    Path(__file__).resolve().parent / "canonicalize_3d8_fixed_coordinate_shells.sage"
)
COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage"
)
SURFACE_COMMON_SOURCE = (
    Path(__file__).resolve().parent / "canonicalize_8a3_fixed_coordinate_shells.sage"
)
SHARED = runpy.run_path(
    str(SHARED_CANON_SOURCE), run_name="_3a8_fixed_shell_canonicalization_common"
)

canonicalize = SHARED["canonicalize"]
digest = SHARED["digest"]

CONFIG = {
    "ambient_label": "3A8",
    "backend_id": "ROOTED-3A8",
    "residual_schema": "elkies-k3.3a8-residual-group.v1",
    "residual_status": "PASS_EXACT_3A8_GLUE_AND_RESIDUAL_GROUP",
    "probe_schema": "elkies-k3.3a8-fixed-coordinate-shell-probe.v1",
    "probe_status": (
        "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_3A8_COORDINATE_SHELL_SCAN"
    ),
    "expected_group_order": 12,
    "expected_seed_count": 189,
    "orbit_id_prefix": "3A8-CF",
}


def build(catalog, residual, probe):
    backend = canonicalize(catalog, residual, probe, CONFIG)
    assert (
        backend["accounting"]["residual_group_embedding_orbits"],
        backend["accounting"][
            "k3_compatible_residual_group_embedding_orbits"
        ],
        backend["accounting"]["surface_classes_after_T_NS_first_dedup"],
        backend["accounting"][
            "partner_auxiliary_isometry_classes_after_surface_dedup"
        ],
        backend["accounting"][
            "frame_isometry_classes_after_surface_dedup"
        ],
    ) == (189, 189, 25, 30, 64)
    assert backend["accounting"]["residual_group_orbit_size_distribution"] == {
        "6": 189
    }
    assert backend["accounting"]["post_dedup_frame_mw_rank_distribution"] == {
        "12": 55,
        "13": 9,
    }
    return {
        "schema": "elkies-k3.3a8-fixed-coordinate-shells.v1",
        "status": "PASS_EXACT_DECLARED_3A8_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        "proof_scope": {
            "proved": (
                "complete order-twelve residual lift group; every 7-of-r "
                "coordinate summand in all four eligible fixed-class LLL bases; "
                "exact determinant, length, MW12--17, mod-two, residual quotient, "
                "ternary, and T/NS-first gates"
            ),
            "not_proved": (
                "all primitive rank-seven fixed-lattice or ambient sublattices, "
                "full Weyl embedding orbits, ternary class enumeration beyond "
                "genera, or determinant-band completeness"
            ),
        },
        "parameters": probe["parameters"],
        "backends": [backend],
        "accounting": {
            "backends": 1,
            "coordinate_subsets_tested": backend["source_probe_accounting"][
                "coordinate_subsets_tested"
            ],
            "accepted_seeds_before_residual_dedup": backend["accounting"][
                "high_mw_mod2_accepted_seeds_before_residual_dedup"
            ],
            "residual_group_embedding_orbits": backend["accounting"][
                "residual_group_embedding_orbits"
            ],
            "k3_compatible_residual_group_embedding_orbits": backend[
                "accounting"
            ]["k3_compatible_residual_group_embedding_orbits"],
            "surface_classes_before_global_cross_backend_dedup": backend[
                "accounting"
            ]["surface_classes_after_T_NS_first_dedup"],
            "partner_auxiliary_classes_before_global_cross_backend_dedup": (
                backend["accounting"][
                    "partner_auxiliary_isometry_classes_after_surface_dedup"
                ]
            ),
            "frame_classes_before_global_cross_backend_dedup": backend[
                "accounting"
            ]["frame_isometry_classes_after_surface_dedup"],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--residual", type=Path, default=RESIDUAL)
    parser.add_argument("--probe", type=Path, default=PROBE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build(
        json.loads(arguments.catalog.read_text()),
        json.loads(arguments.residual.read_text()),
        json.loads(arguments.probe.read_text()),
    )
    payload["inputs"] = {
        str(arguments.catalog.resolve().relative_to(ROOT)): digest(arguments.catalog),
        str(arguments.residual.resolve().relative_to(ROOT)): digest(arguments.residual),
        str(arguments.probe.resolve().relative_to(ROOT)): digest(arguments.probe),
        str(SHARED_CANON_SOURCE.resolve().relative_to(ROOT)): digest(
            SHARED_CANON_SOURCE
        ),
        str(COMMON_SOURCE.resolve().relative_to(ROOT)): digest(COMMON_SOURCE),
        str(SURFACE_COMMON_SOURCE.resolve().relative_to(ROOT)): digest(
            SURFACE_COMMON_SOURCE
        ),
    }
    payload["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/canonicalize_3a8_fixed_coordinate_shells.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("3A8 canonical shell artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "3A8CANON|seeds={}|orbits={}|surfaces={}|partners={}|frames={}|status=PASS_EXACT".format(
            payload["accounting"]["accepted_seeds_before_residual_dedup"],
            payload["accounting"]["residual_group_embedding_orbits"],
            payload["accounting"][
                "surface_classes_before_global_cross_backend_dedup"
            ],
            payload["accounting"][
                "partner_auxiliary_classes_before_global_cross_backend_dedup"
            ],
            payload["accounting"][
                "frame_classes_before_global_cross_backend_dedup"
            ],
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
