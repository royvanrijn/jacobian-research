#!/usr/bin/env python3
"""Replay the pinned CRT--lattice calibration with exact checks and PARI/GP."""

from __future__ import annotations

import ast
import json
import shutil
import subprocess
import sys
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.calibration import build_calibration  # noqa: E402
from ecsearch.crt_lattice import p_adic_valuation  # noqa: E402
from ecsearch.local_data import calibration_family_local_data  # noqa: E402


ARTIFACT = (
    REPOSITORY_ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic-curves"
    / "crt_lattice_calibration_v1.json"
)


def main() -> None:
    gp = shutil.which("gp")
    if gp is None:
        raise SystemExit("PARI/GP executable 'gp' is required")
    rank_two = subprocess.run(
        [gp, "-q", "-f"],
        input="setrand(1);E=ellinit([0,0,0,-25,25]);print(ellrank(E,2));\n",
        text=True,
        capture_output=True,
        check=True,
    )
    assert "***" not in rank_two.stdout + rank_two.stderr
    assert ast.literal_eval(rank_two.stdout.strip()) == [
        2,
        2,
        0,
        [[5, 5], [-5, 5]],
    ]
    expected = json.loads(ARTIFACT.read_text())
    actual = build_calibration(
        maximum_height=expected["search"]["maximum_height"]
    )
    assert actual == expected, "pinned calibration is stale"
    best = actual["best_candidate"]
    assert (best["numerator"], best["denominator"]) == (-110627, 84367)
    assert best["crt_modulus"] == 143227016087
    assert best["binary_discriminant_factor"] == best["crt_modulus"]
    assert best["pari_gp"]["global_minimal_change"] == [1, 0, 0, 0]
    assert best["pari_gp"]["torsion_order"] == 1
    assert best["pari_gp"]["known_point_order"] == 0
    assert best["pari_gp"]["rank_bounds"] == [3, 3]
    conductor_exponents = dict(best["pari_gp"]["conductor_factorization"])
    assert all(conductor_exponents[prime] == 1 for prime in (23, 47, 73))
    assert [
        p_adic_valuation(best["discriminant"], prime)
        for prime in (23, 47, 73)
    ] == [3, 2, 2]
    assert all(
        constraint["split_test_legendre_6"] == 1
        for constraint in actual["search"]["constraints"]
    )
    for constraint, root in zip(
        actual["search"]["constraints"],
        actual["search"]["combinations"][0]["roots"],
        strict=True,
    ):
        local = calibration_family_local_data(root, constraint["prime"])
        assert local.reduction == "split_multiplicative"
        assert local.trace is None
        assert local.local_euler_coefficient == 1
    print(
        "PASS CRT--lattice calibration: exact prime powers in the minimal "
        "discriminant, conductor exponent one at 23,47,73, PARI rank [3,3]"
    )


if __name__ == "__main__":
    main()
