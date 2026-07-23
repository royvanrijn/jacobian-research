#!/usr/bin/env python3
"""Exact closure certificates for the cubic two-real-variable GMC frontier.

The previously verified frontier leaves three five-weight chart orbits, all
six-weight charts, and the full seven-weight support.  This script builds the
centered circular-coordinate moment systems directly and asks msolve for a
reduced Groebner basis over QQ after a Rabinowitsch localization.

Every successful calculation must return the literal characteristic-zero
basis ``[1]``.  Finite-field output is not used as a proof.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch import msolve
from verify_two_real_gmc_five_weight import (
    generic_moments,
    primitive_numerator,
    reflected_support,
    support_data,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "cubic_gaussian_null_cone_closure.json"
)


def support_representatives(size: int) -> list[tuple[int, ...]]:
    """Return reflection representatives among supports of the given size."""

    return sorted(
        {
            min(support, reflected_support(support))
            for support in itertools.combinations(range(-3, 4), size)
        }
    )


def chart_choices(
    components: dict[int, tuple[sp.Symbol, ...]],
    support: tuple[int, ...],
) -> list[tuple[int, ...]]:
    return list(
        itertools.product(
            *(range(len(components[weight])) for weight in support)
        )
    )


def build_system(
    support: tuple[int, ...],
    chart: tuple[int, ...],
    cutoff: int,
    normalization: str,
    eliminate_second_moment: bool,
) -> dict[str, object]:
    """Build one normalized, localized characteristic-zero moment system."""

    components, coefficients, terms = support_data(support)
    moments = generic_moments(coefficients, terms, cutoff)
    selected = tuple(
        components[weight][choice]
        for weight, choice in zip(support, chart)
    )

    if normalization == "two_negative":
        normalized = {
            coefficient: 1
            for weight, coefficient in zip(support, selected)
            if weight < 0
        }
        normalized = dict(list(normalized.items())[:2])
    elif normalization == "mixed":
        negative = next(
            coefficient
            for weight, coefficient in zip(support, selected)
            if weight < 0
        )
        positive = next(
            coefficient
            for weight, coefficient in zip(support, selected)
            if weight > 0
        )
        normalized = {negative: 1, positive: 1}
    elif normalization == "all_selected":
        normalized = {coefficient: 1 for coefficient in selected}
    else:
        raise ValueError(f"unknown normalization {normalization!r}")
    assert len(normalized) == 2

    variables = tuple(
        coefficient
        for coefficient in coefficients
        if coefficient not in normalized
    )
    normalized_moments = {
        order: sp.expand(moment.subs(normalized))
        for order, moment in moments.items()
    }
    localization_factor = sp.prod(
        coefficient.subs(normalized)
        for coefficient in selected
        if coefficient not in normalized
    )
    elimination_record = None

    if eliminate_second_moment:
        for variable in reversed(variables):
            polynomial = sp.Poly(normalized_moments[2], variable)
            if (
                polynomial.degree() != 1
                or polynomial.LC().free_symbols
            ):
                continue
            solution = sp.solve(normalized_moments[2], variable)[0]
            remaining = tuple(
                candidate
                for candidate in variables
                if candidate != variable
            )
            normalized_moments = {
                order: (
                    sp.Integer(0)
                    if order == 2
                    else primitive_numerator(moment.subs(variable, solution), remaining)
                )
                for order, moment in normalized_moments.items()
            }
            localization_factor = primitive_numerator(
                localization_factor.subs(variable, solution),
                remaining,
            )
            variables = remaining
            elimination_record = {
                "variable": str(variable),
                "solution": str(solution),
            }
            break
    rho = sp.Symbol("rho")
    equations = [
        normalized_moments[order]
        for order in range(2, cutoff + 1)
        if normalized_moments[order] != 0
    ]
    equations.append(sp.expand(rho * localization_factor - 1))
    solver_variables = variables + (rho,)
    input_text = msolve.input_text(equations, solver_variables, characteristic=0)

    return {
        "support": support,
        "chart": chart,
        "cutoff": cutoff,
        "normalization_kind": normalization,
        "normalization": {
            str(coefficient): value
            for coefficient, value in normalized.items()
        },
        "selected_coefficients": tuple(map(str, selected)),
        "variables": solver_variables,
        "equations": equations,
        "localization_factor": localization_factor,
        "second_moment_elimination": elimination_record,
        "input_sha256": hashlib.sha256(input_text.encode()).hexdigest(),
    }


def closure_systems() -> list[dict[str, object]]:
    """Return the 31 finite characteristic-zero closure systems."""

    systems: list[dict[str, object]] = []

    # The earlier two-weight argument used the whole circuit-moment sequence.
    # These eleven presentations certify the claimed finite cutoff ten.
    two_weight_supports = sorted(
        {
            min(support, reflected_support(support))
            for support in itertools.product(range(-3, 0), range(1, 4))
        }
    )
    for support in two_weight_supports:
        components, _, _ = support_data(support)
        for chart in chart_choices(components, support):
            system = build_system(
                support,
                chart,
                cutoff=10,
                normalization="all_selected",
                eliminate_second_moment=False,
            )
            system["stratum"] = "two_weight_finite"
            systems.append(system)

    # The three chart orbits not promoted over QQ in the earlier five-weight
    # census.  Reflection makes them account for ten ordinary charts.
    five_weight_cases = (
        ((-3, -1, 0, 1, 2), (0, 0, 0, 0, 0)),
        ((-3, -1, 0, 1, 2), (0, 0, 0, 1, 0)),
        ((-2, -1, 0, 1, 2), (0, 0, 0, 0, 0)),
    )
    for support, chart in five_weight_cases:
        system = build_system(
            support,
            chart,
            cutoff=8,
            normalization="two_negative",
            eliminate_second_moment=True,
        )
        system["stratum"] = "five_weight_promotion"
        systems.append(system)

    # Four reflection-representative supports and fourteen centered chart
    # presentations account for all 44 ordinary six-weight charts.
    for support in support_representatives(6):
        components, _, _ = support_data(support)
        for chart in chart_choices(components, support):
            system = build_system(
                support,
                chart,
                cutoff=8,
                normalization="mixed",
                eliminate_second_moment=False,
            )
            system["stratum"] = "six_weight"
            systems.append(system)

    # The full support is self-reflective.  The two mixed +/-1 coefficient
    # charts are exchanged by Z <-> W, leaving three chart orbits.  Centering
    # also identifies the two original weight-zero coefficient charts.
    seven_support = tuple(range(-3, 4))
    seven_charts = (
        (0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 1, 0, 0),
        (0, 0, 1, 0, 1, 0, 0),
    )
    for chart in seven_charts:
        system = build_system(
            seven_support,
            chart,
            cutoff=10,
            normalization="two_negative",
            eliminate_second_moment=False,
        )
        system["stratum"] = "seven_weight"
        systems.append(system)

    assert len(systems) == 31
    assert sum(system["stratum"] == "two_weight_finite" for system in systems) == 11
    assert sum(system["stratum"] == "five_weight_promotion" for system in systems) == 3
    assert sum(system["stratum"] == "six_weight" for system in systems) == 14
    assert sum(system["stratum"] == "seven_weight" for system in systems) == 3
    return systems


def verify_system(system: dict[str, object]) -> dict[str, object]:
    """Run msolve over QQ and require its literal reduced basis [1]."""

    started = time.monotonic()
    result = msolve.run(
        system["equations"],
        system["variables"],
        prime=0,
        threads=4,
        groebner=True,
        timeout=900,
    )
    elapsed = round(time.monotonic() - started, 3)
    exact_header = "#field characteristic: 0" in result.output
    exact_unit = (
        result.returncode == 0
        and result.empty
        and exact_header
        and "#length of basis:      1 element" in result.output
        and result.output.rstrip().endswith("[1]:")
    )
    if not exact_unit:
        raise RuntimeError(
            "characteristic-zero unit certificate failed for "
            f"{system['support']} {system['chart']}: "
            f"returncode={result.returncode}, "
            f"stdout={result.output[-1000:]!r}, stderr={result.stderr[-1000:]!r}"
        )
    return {
        "stratum": system["stratum"],
        "support": list(system["support"]),
        "chart": list(system["chart"]),
        "cutoff": system["cutoff"],
        "normalization_kind": system["normalization_kind"],
        "normalization": system["normalization"],
        "selected_coefficients": list(system["selected_coefficients"]),
        "variables": [str(variable) for variable in system["variables"]],
        "equation_count": len(system["equations"]),
        "localization_factor": str(system["localization_factor"]),
        "second_moment_elimination": system["second_moment_elimination"],
        "input_sha256": system["input_sha256"],
        "field": "QQ",
        "reduced_basis": ["1"],
        "elapsed_seconds": elapsed,
    }


def main() -> None:
    assert msolve.available(), "this certificate requires msolve"
    records = []
    for system in closure_systems():
        record = verify_system(system)
        records.append(record)
        print(
            "PASS cubic null cone:",
            record["stratum"],
            record["support"],
            record["chart"],
            f"M<={record['cutoff']}",
            f"{record['elapsed_seconds']}s",
            flush=True,
        )

    payload = {
        "format": "cubic-gaussian-null-cone-closure-v1",
        "field": "QQ",
        "certificate": "literal reduced Groebner basis [1]",
        "finite_closure_presentations": len(records),
        "two_weight_finite_presentations": 11,
        "five_weight_promotions": 3,
        "five_weight_ordinary_charts_newly_closed": 10,
        "six_weight_presentations": 14,
        "six_weight_ordinary_charts_closed": 44,
        "seven_weight_chart_orbits": 3,
        "seven_weight_ordinary_charts_closed": 8,
        "moment_cutoff": {
            "five_weight": 8,
            "six_weight": 8,
            "seven_weight": 10,
        },
        "theorem": (
            "Together with the existing exact at-most-four-weight and "
            "five-weight certificates, every mixed-sign cubic coefficient "
            "support is excluded by moments through order ten.  Equivalently, "
            "the cubic torus null-cone ideal is contained in the radical of "
            "<M1,...,M10>."
        ),
        "records": records,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS cubic null cone: wrote", OUTPUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
