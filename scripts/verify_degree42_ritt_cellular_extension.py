#!/usr/bin/env python3
"""Verify the first non-split degree-42 cellular coefficient extension."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from jcsearch.hessian_ritt_cellular import braid_totalization  # noqa: E402
from jcsearch.ritt_complex import MoveType, symmetric_braid_complex  # noqa: E402
from research_degree42_cellular_extension import (  # noqa: E402
    CACHE,
    audit,
)


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_cellular_extension.json"
)


def matrix(data) -> sp.Matrix:
    return sp.Matrix(data)


def check_module_actions(result: dict[str, object]) -> None:
    """Check that every displayed pair defines a Q[tau,zeta]-module."""

    for module_name in ("sector", "spectator", "total"):
        actions = result[f"{module_name}_actions"]
        tau = matrix(actions["tau"])
        zeta = matrix(actions["zeta"])
        assert tau * zeta == zeta * tau
        assert tau**3 == sp.zeros(tau.rows)
        assert zeta**3 == sp.zeros(zeta.rows)


def main() -> None:
    assert CACHE.is_file(), (
        "missing source-ideal cache; run "
        "scripts/research_degree42_cellular_extension.py "
        "--order 2 --rebuild-source"
    )
    order_two = audit(2)
    order_three = audit(3)
    order_four = audit(4)
    check_module_actions(order_two)
    check_module_actions(order_three)
    check_module_actions(order_four)

    assert order_two["ring_lengths"] == {
        "A6": 5,
        "A_boundary": 4,
        "B": 3,
    }
    assert order_two["module_dimensions"] == {
        "sector": 1,
        "spectator": 1,
        "total": 2,
    }
    assert order_two["splits_over_Q_tau_zeta"] is True
    assert order_two["extension_certificate"] == {
        "adapted_coupling": {"tau": [[0]], "zeta": [[0]]},
        "coboundary_rank": 0,
        "augmented_rank": 0,
        "obstruction_functional": None,
        "obstruction_value": None,
    }

    assert order_three["ring_lengths"] == {
        "A6": 14,
        "A_boundary": 9,
        "B": 6,
    }
    assert order_three["module_dimensions"] == {
        "sector": 5,
        "spectator": 3,
        "total": 8,
    }
    assert order_three["splits_over_Q_tau_zeta"] is False
    certificate = order_three["extension_certificate"]
    assert certificate["adapted_coupling"] == {
        "tau": sp.zeros(5, 3).tolist(),
        "zeta": [
            [0, 1, 10],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
        ],
    }
    assert certificate["coboundary_rank"] == 7
    assert certificate["augmented_rank"] == 8
    assert certificate["obstruction_functional"] == {
        "zeta[0,1]": 1,
        "zeta[1,1]": 10,
    }
    assert certificate["obstruction_value"] == "1"

    assert order_four["ring_lengths"] == {
        "A6": 29,
        "A_boundary": 16,
        "B": 10,
    }
    assert order_four["module_dimensions"] == {
        "sector": 13,
        "spectator": 6,
        "total": 19,
    }
    assert order_four["splits_over_Q_tau_zeta"] is False
    certificate_four = order_four["extension_certificate"]
    assert certificate_four["coboundary_rank"] == 52
    assert certificate_four["augmented_rank"] == 53
    assert certificate_four["obstruction_functional"] == {
        "tau[1,0]": 12,
        "tau[2,0]": 100,
        "tau[4,1]": 12,
        "tau[5,1]": 100,
        "tau[7,0]": -100,
        "tau[8,1]": -100,
        "tau[9,1]": -300,
        "tau[11,3]": 80,
        "zeta[1,1]": -3,
        "zeta[2,1]": -30,
    }
    assert certificate_four["obstruction_value"] == "-3"

    braid = symmetric_braid_complex((2, 3, 7), MoveType.CHEBYSHEV)
    order_two_totalization = braid_totalization(
        braid,
        base_dimension=3,
        defect_dimensions=(2,),
        defect_names=("nonsplit-total-order-2",),
        name="degree-42 order-two nested totalization",
    )
    order_three_totalization = braid_totalization(
        braid,
        base_dimension=6,
        defect_dimensions=(8,),
        defect_names=("nonsplit-total-order-3",),
        name="degree-42 order-three nested totalization",
    )
    order_four_totalization = braid_totalization(
        braid,
        base_dimension=10,
        defect_dimensions=(19,),
        defect_names=("nonsplit-total-order-4",),
        name="degree-42 order-four nested totalization",
    )
    assert order_two_totalization.complex.cohomology_dimensions == (3, 2, 0)
    assert order_three_totalization.complex.cohomology_dimensions == (6, 8, 0)
    assert order_four_totalization.complex.cohomology_dimensions == (
        10,
        19,
        0,
    )

    source_hash = hashlib.sha256(CACHE.read_bytes()).hexdigest()
    result = {
        "schema": "degree42-ritt-cellular-extension.v1",
        "status": "exact finite-jet module computation",
        "source_ideal_cache": str(CACHE.relative_to(ROOT)),
        "source_ideal_cache_sha256": source_hash,
        "order_two": order_two,
        "order_three": order_three,
        "order_four": order_four,
        "cellular_totalization": {
            "order_two_cohomology_dimensions": [3, 2, 0],
            "order_three_cohomology_dimensions": [6, 8, 0],
            "order_four_cohomology_dimensions": [10, 19, 0],
            "H1_interpretation": (
                "the order-three and order-four H1 are non-split total "
                "modules K/I_6, "
                "not a direct sum of its sector and spectator subquotients"
            ),
        },
        "conclusion": (
            "the sector-spectator extension splits at first conormal order "
            "and is non-split at orders three and four; its first coupling "
            "is carried by zeta multiplication"
        ),
        "theorem_boundary": {
            "proved": [
                "exact order-two through order-four Q[tau,zeta]-module actions",
                "exact non-splitting certificates at orders three and four",
                "compatibility with the cellular H1/H2 dimensions",
            ],
            "not_proved": [
                "non-splitting of the untruncated completed extension",
                "identification with the full derived cotangent transitivity class",
                "filtered H2-vanishing beyond the displayed jets",
            ],
        },
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_degree42_ritt_cellular_extension.py"
        ),
    }
    ARTIFACT.write_text(
        json.dumps(
            result,
            indent=2,
            default=lambda value: (
                int(value) if value.is_Integer else str(value)
            ),
        )
        + "\n"
    )
    print("PASS: order-two sector/spectator extension splits")
    print("PASS: order-three sector/spectator extension is non-split")
    print("PASS: order-four sector/spectator extension remains non-split")
    print("PASS: explicit obstruction functionals evaluate to one and minus three")
    print("PASS: cellular H1/H2 dimensions through order four agree")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
