#!/usr/bin/env python3
"""Prepare subgroup-only inputs for held-out half-lattice experiments.

This one-time boundary builder reads the public point fixtures but writes only
the chosen starting subgroups.  The blind search reads the resulting JSON and
cannot access any held-out point coordinate.  Complements are reopened only by
the verification script.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
OUTPUT = ROOT / "elliptic-curves/data/half_lattice_heldout_subgroup_inputs_v1.json"
sys.path.insert(0, str(CAS))

import icarm_curve245  # noqa: E402
import icarm_curve273  # noqa: E402
import icarm_curve302  # noqa: E402


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def point_record(point) -> list[str]:
    return [str(point[0]), str(point[1])]


def configurations(point_count: int) -> list[tuple[str, list[int]]]:
    rows = [
        (f"prefix-d{dimension}", list(range(dimension)))
        for dimension in (12, 15, 17, 18)
    ]
    # A deterministic spread across the displayed basis supplies a same-rank
    # primitive adverse comparison without consulting heights or held-out hits.
    permutation = list(range(0, point_count, 2)) + list(range(1, point_count, 2))
    rows.append(("interleaved-d17", permutation[:17]))
    return rows


def curve_record(label: str, module, *, include_all_configurations: bool) -> dict:
    all_configs = configurations(len(module.SHORT_POINTS))
    if not include_all_configurations:
        all_configs = [row for row in all_configs if row[0] in {"prefix-d12", "prefix-d17", "prefix-d18"}]
    return {
        "label": label,
        "short_model": [str(value) for value in module.short_coefficients()],
        "public_displayed_point_count": len(module.SHORT_POINTS),
        "configurations": [
            {
                "id": identifier,
                "dimension": len(indices),
                "included_public_indices_one_based": [index + 1 for index in indices],
                "starting_subgroup_points": [
                    point_record(module.SHORT_POINTS[index]) for index in indices
                ],
                "selection_policy": (
                    "ordered public prefix"
                    if identifier.startswith("prefix")
                    else "even public indices followed by odd public indices"
                ),
            }
            for identifier, indices in all_configs
        ],
    }


def main() -> None:
    payload = {
        "schema": "elliptic-curves.half-lattice-heldout-subgroup-inputs.v1",
        "status": "FROZEN_SUBGROUP_ONLY_INPUTS_NO_HELDOUT_COORDINATES",
        "boundary": {
            "builder_read_public_fixtures": True,
            "output_contains_heldout_coordinates": False,
            "search_must_not_import_public_curve_modules": True,
            "configuration_selection_used_heights_or_search_hits": False,
        },
        "source_hashes": {
            str((CAS / "icarm_curve245.py").relative_to(ROOT)): digest(CAS / "icarm_curve245.py"),
            str((CAS / "icarm_curve273.py").relative_to(ROOT)): digest(CAS / "icarm_curve273.py"),
            str((CAS / "icarm_curve302.py").relative_to(ROOT)): digest(CAS / "icarm_curve302.py"),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "curves": [
            curve_record("curve273", icarm_curve273, include_all_configurations=True),
            curve_record("curve302", icarm_curve302, include_all_configurations=True),
            curve_record("curve245-adverse-control", icarm_curve245, include_all_configurations=False),
        ],
        "claim_boundary": (
            "Each selected coordinate subgroup is primitive in the displayed free public-point lattice because it is spanned by a subset of its certified basis. "
            "This does not assert primitivity in the full Mordell-Weil group."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"HALFLATTICEINPUTS|status=PASS|configurations=13|output={OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
