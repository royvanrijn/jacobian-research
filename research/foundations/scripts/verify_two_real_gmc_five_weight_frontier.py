#!/usr/bin/env python3
"""Replay the proved characteristic-zero part of the five-weight census.

The exploratory compiler also retains four slow presentations for future
promotion.  This checker deliberately omits those presentations, recomputes
the 31 rational unit ideals that constitute the proved frontier, and then
uses the exact Z/W chart reflection to audit the published counts of
18 supports and 92 ordinary charts.
"""

from __future__ import annotations

import shutil
from concurrent.futures import ThreadPoolExecutor

from verify_two_real_gmc_five_weight import (
    chart_systems,
    reflected_support,
    run_basis,
)


UNPROMOTED = {
    ((-3, -1, 0, 1, 2), (0, 0, 0, 0, 0)),
    ((-3, -1, 0, 1, 2), (0, 0, 0, 1, 0)),
    ((-2, -1, 0, 1, 2), (0, 0, 0, 0, 0)),
    # This mixed chart is supplied by reflection of the exact (1,0) chart.
    ((-2, -1, 0, 1, 2), (0, 0, 0, 1, 0)),
}


def main() -> None:
    singular = shutil.which("Singular")
    assert singular is not None, "this checker requires Singular"
    systems = chart_systems()
    direct_systems = [
        system
        for system in systems
        if (system["support"], system["chart"]) not in UNPROMOTED
    ]
    assert len(direct_systems) == 31

    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(
            executor.map(
                lambda system: run_basis(
                    singular,
                    system,
                    cutoff=8,
                    characteristic=0,
                    timeout=900,
                ),
                direct_systems,
            )
        )
    assert all(result["is_unit"] == 1 for result in results)

    direct_keys = {
        (system["support"], system["chart"]) for system in direct_systems
    }

    def is_exact(system: dict[str, object]) -> bool:
        key = (system["support"], system["chart"])
        if key in direct_keys:
            return True
        support = system["support"]
        return (
            reflected_support(support) == support
            and (support, tuple(reversed(system["chart"]))) in direct_keys
        )

    ordinary_charts = sum(
        int(system["raw_charts_covered"])
        for system in systems
        if is_exact(system)
    )
    excluded_supports = 0
    for support in sorted({system["support"] for system in systems}):
        support_systems = [
            system for system in systems if system["support"] == support
        ]
        if all(is_exact(system) for system in support_systems):
            excluded_supports += (
                1 if reflected_support(support) == support else 2
            )

    assert ordinary_charts == 92
    assert excluded_supports == 18
    print(
        "PASS five-weight GMC frontier: 31 rational unit presentations "
        "exclude 18 supports / 92 ordinary charts"
    )


if __name__ == "__main__":
    main()
