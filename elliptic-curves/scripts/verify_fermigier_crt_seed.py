#!/usr/bin/env python3
"""Replay the pinned Fermigier high-family CRT seed."""

from __future__ import annotations

import json
import sys
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.fermigier_seed import build_fermigier_seed  # noqa: E402


ARTIFACT = (
    REPOSITORY_ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic-curves"
    / "fermigier_crt_seed_v1.json"
)


def main() -> None:
    expected = json.loads(ARTIFACT.read_text())
    actual = build_fermigier_seed(
        maximum_height=expected["search"]["maximum_height"]
    )
    assert actual == expected, "pinned Fermigier CRT seed is stale"
    best = actual["best_seed"]
    assert (best["numerator"], best["denominator"]) == (673709, 29965)
    assert best["crt_modulus"] == 2551312982089
    assert best["height"] == 673709
    assert actual["search"]["root_combinations_tested"] == 8
    for constraint in actual["search"]["constraints"]:
        assert len(constraint["roots_mod_prime_power"]) == 2
        for reduction in constraint["root_reductions"]:
            assert reduction["derivative_mod_prime"] != 0
            assert reduction["reduction"] == "split_multiplicative"
            assert reduction["local_euler_coefficient"] == 1
    assert all(
        value == 2
        for value in actual["search"]["combinations"][0][
            "forced_valuations"
        ].values()
    )
    for prime in ("89", "131", "137"):
        local = best["pari_gp"]["local_data"][prime]
        assert local == {
            "conductor_exponent": 1,
            "kodaira_code": 6,
            "kodaira_symbol": "I_2",
            "tamagawa_number": 2,
            "minimal_discriminant_valuation": 2,
            "local_euler_coefficient": 1,
        }
    assert actual["limitations"]["global_conductor"].startswith("not computed")
    print(
        "PASS Fermigier CRT seed: exhaustive 8-class search, height 673709, "
        "split I_2 and conductor exponent one at 89,131,137; global conductor open"
    )


if __name__ == "__main__":
    main()
