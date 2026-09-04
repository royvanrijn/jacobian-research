#!/usr/bin/env sage-python
"""Cross-class amendment runner for the frozen R17 jump ladder.

Protocol v1 incorrectly asserted that every determinant-948 native MW17
lattice has exactly 43 maximum-depth parity classes.  The first case, curve
478, has that property and completed with response 6.  Before curve 498 ran a
generic-census assertion showed that the property does not transfer across
native classes.  This wrapper changes only the initial set definition to the
first 43 classes ordered by exact generic norm and mask.  It then invokes the
otherwise unchanged v1 blind runner under the v2 protocol.
"""

from __future__ import annotations

from importlib.machinery import SourceFileLoader
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
V1_RUNNER = ROOT / "elliptic-curves/cas/run_r17_refresh_jump_ladder_blind.sage"
V2_PROTOCOL = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_protocol_v2.json"
V2_OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_blind_v2.json"


def complete_generic_census_cross_class(module, gram):
    oracle = module.CosetOracle(gram)
    rows = []
    maximum_error = 0.0
    for mask in range(1 << 17):
        residue = tuple((mask >> index) & 1 for index in range(17))
        norm, representative, error = oracle.solve(residue)
        maximum_error = max(maximum_error, error)
        rows.append((norm, mask, representative))
    rows.sort(key=lambda row: (-row[0], row[1]))
    if len(rows) != 1 << 17 or len(rows[:301]) != 301:
        raise ArithmeticError("the complete cross-class generic census is incomplete")
    return rows, maximum_error


def main():
    runner = SourceFileLoader("r17_jump_ladder_v1_runner", str(V1_RUNNER)).load_module()
    runner.PROTOCOL = V2_PROTOCOL
    runner.OUTPUT = V2_OUTPUT
    runner.complete_generic_census = lambda legacy, gram: complete_generic_census_cross_class(
        legacy, gram
    )
    runner.main()


if __name__ == "__main__":
    main()
