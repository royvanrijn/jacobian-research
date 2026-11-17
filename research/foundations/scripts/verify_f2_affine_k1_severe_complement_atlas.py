#!/usr/bin/env sage
"""Verify the severe k=1 complement groups and cubic degree-six actions.

Run with SageMath and its optional ``sirocco`` package.  Unlike the nearby
research driver, this file is an assertion-based theorem checker: it
recomputes every implicit curve and Zariski--van Kampen group from the exact
rational witness before checking the pinned complement/action atlas.
"""

from __future__ import annotations

from sage.all import QQ  # noqa: F401 -- initialize Sage before zariski_vankampen
from sage.schemes.curves.zariski_vankampen import fundamental_group

from research_f2_affine_k1_severe_complements import (
    cubic_degree_six_actions,
    generated_group,
    implicit_polynomial,
    witnesses,
)


CYCLIC = {
    "A2+A3+A1",
    "A4+2A1",
    "D4+A1",
    "D4+A2",
    "D5+A1",
    "D6",
    "E7-I3",
    "E7-A4",
}

NONCYCLIC_CUBIC_COUNTS = {
    "A4+A3": (0, 0),
    "A6+A1": (0, 0),
    "A4+A2+A1": (18, 1),
    "D5+A2": (0, 0),
}


def main() -> None:
    rows = witnesses()
    assert {row.label for row in rows} == CYCLIC | set(NONCYCLIC_CUBIC_COUNTS)

    for row in rows:
        polynomial = implicit_polynomial(row.parameters)
        simplified = fundamental_group(polynomial, simplified=True, projective=False)
        assert tuple(simplified.abelian_invariants()) == (0,), row.label

        if row.label in CYCLIC:
            assert simplified.ngens() == 1, row.label
            assert tuple(simplified.relations()) == (), row.label
            continue

        raw = fundamental_group(polynomial, simplified=False, projective=False)
        solutions, representatives = cubic_degree_six_actions(raw)
        expected_solutions, expected_classes = NONCYCLIC_CUBIC_COUNTS[row.label]
        assert len(solutions) == expected_solutions, row.label
        assert len(representatives) == expected_classes, row.label
        if row.label == "A4+A2+A1":
            assert {len(generated_group(images)) for images in solutions} == {360}

    print(
        "PASS: all eight severe cyclic rows have complement group Z; "
        "A4+A3, A6+A1, and D5+A2 have no cubic degree-six action; "
        "A4+A2+A1 has one class with image order 360"
    )


if __name__ == "__main__":
    main()
