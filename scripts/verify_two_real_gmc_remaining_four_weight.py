#!/usr/bin/env python3
"""Exclude the remaining three four-weight cubic GMC(2) supports over QQ.

The support census and the symmetric-support quotient certificates leave

    {-2,-1,0,1}, {-2,-1,1,2}, {-1,0,1,2}.

For the zero-weight support, the first moment centers its radial component:
``g+h*Z*W`` becomes ``h*(Z*W-1)``.  The two coefficient choices in that
component therefore describe the same nonvanishing locus.  Four exact chart
representatives cover the support and, after swapping ``Z,W``, its
reflection.

For the balanced support, the second moment eliminates the coefficient of
``Z^2``.  Three representatives cover its four charts up to the same
reflection.  Singular proves that every localized rational moment ideal
through order six has reduced standard basis ``[1]``.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import sympy as sp

from verify_two_real_gmc_frontier import W, Z, expectation
from verify_two_real_gmc_symmetric_chart import (
    primitive_numerator,
    singular_expression,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_real_gmc_remaining_four_weight.json"
)


def moments(polynomial: sp.Expr, maximum: int = 6) -> dict[int, sp.Expr]:
    """Return exact pure moments through ``maximum``."""

    result: dict[int, sp.Expr] = {}
    power = sp.Integer(1)
    for order in range(1, maximum + 1):
        power = sp.expand(power * polynomial)
        result[order] = sp.factor(expectation(power))
    return result


def centered_zero_weight_systems() -> list[dict[str, object]]:
    """Build four representatives for E2 and its reflected support."""

    a, b, c, d, h, rho = sp.symbols("a b c d h rho")
    polynomial = (
        W**2
        + b * W
        + a * Z * W**2
        + h * (Z * W - 1)
        + c * Z
        + d * Z**2 * W
    )
    pure_moments = moments(polynomial)
    assert pure_moments[1] == 0

    systems: list[dict[str, object]] = []
    for negative_kind, normalization in (
        ("constant", {b: 1}),
        ("radial", {a: 1}),
    ):
        for positive_kind, selected in (("constant", c), ("radial", d)):
            variables = tuple(
                variable
                for variable in (a, b, h, c, d)
                if variable not in normalization
            )
            equations = tuple(
                primitive_numerator(
                    pure_moments[order].subs(normalization),
                    variables,
                )
                for order in range(2, 7)
            )
            localization_factor = primitive_numerator(h * selected, variables)
            systems.append(
                {
                    "key": f"centered_{negative_kind}_{positive_kind}",
                    "support_orbit": (
                        (-2, -1, 0, 1),
                        (-1, 0, 1, 2),
                    ),
                    "normalization": (
                        "W^2 coefficient and selected weight -1 coefficient "
                        "normalized to one"
                    ),
                    "selected_negative_coefficient": negative_kind,
                    "selected_positive_coefficient": positive_kind,
                    "variables": variables + (rho,),
                    "ideal_generators": equations
                    + (sp.expand(rho * localization_factor - 1),),
                    "localization_factor": localization_factor,
                    "moments_imposed": list(range(1, 7)),
                    # Two zero-weight coefficient charts coincide after
                    # centering, and reflection supplies two more charts.
                    "original_charts_covered": 4,
                }
            )
    return systems


def balanced_systems() -> list[dict[str, object]]:
    """Build three reflection representatives for E3."""

    a, b, c, d, e, rho = sp.symbols("a b c d e rho")
    polynomial = W**2 + b * W + a * Z * W**2 + c * Z + d * Z**2 * W + e * Z**2
    pure_moments = moments(polynomial)
    assert pure_moments[1] == 0

    specs = (
        ("constant_constant", {b: 1}, c, 1),
        ("constant_radial", {b: 1}, d, 2),
        ("radial_radial", {a: 1}, d, 1),
    )
    systems: list[dict[str, object]] = []
    for key, normalization, selected, charts_covered in specs:
        variables = tuple(
            variable
            for variable in (a, b, c, d)
            if variable not in normalization
        )
        second_moment = sp.expand(pure_moments[2].subs(normalization))
        assert sp.degree(second_moment, e) == 1
        e_solution = sp.solve(second_moment, e)[0]
        assert sp.expand(second_moment.subs(e, e_solution)) == 0
        equations = tuple(
            primitive_numerator(
                pure_moments[order]
                .subs(normalization)
                .subs(e, e_solution),
                variables,
            )
            for order in range(3, 7)
        )
        localization_factor = primitive_numerator(
            selected.subs(normalization) * e_solution,
            variables,
        )
        systems.append(
            {
                "key": f"balanced_{key}",
                "support_orbit": ((-2, -1, 1, 2),),
                "normalization": (
                    "W^2 coefficient and selected weight -1 coefficient "
                    "normalized to one"
                ),
                "variables": variables + (rho,),
                "ideal_generators": equations
                + (sp.expand(rho * localization_factor - 1),),
                "localization_factor": localization_factor,
                "second_moment_elimination": f"e = {e_solution}",
                "moments_imposed": list(range(1, 7)),
                "original_charts_covered": charts_covered,
            }
        )
    return systems


def singular_program(system: dict[str, object]) -> str:
    """Emit one exact rational unit-ideal calculation."""

    variables = system["variables"]
    generators = system["ideal_generators"]
    assert isinstance(variables, tuple)
    assert isinstance(generators, tuple)
    variable_names = ",".join(str(variable) for variable in variables)
    ideal_text = ",\n".join(
        singular_expression(generator) for generator in generators
    )
    return f"""
ring rational_ring=0,({variable_names}),dp;
option(redSB);
ideal I={ideal_text};
ideal G=slimgb(I);
int is_unit=0;
if (size(G)==1)
{{
  if (G[1]==1) {{ is_unit=1; }}
}}
print("BASIS_SIZE "+string(size(G)));
print("IS_UNIT "+string(is_unit));
"""


def run_system(
    singular: str,
    system: dict[str, object],
) -> dict[str, object]:
    """Run and validate one exact chart certificate."""

    completed = subprocess.run(
        (singular, "-q"),
        input=singular_program(system),
        text=True,
        capture_output=True,
        check=True,
        timeout=900,
    )
    if "?" in completed.stdout or "?" in completed.stderr:
        raise RuntimeError(
            f"Singular failed on {system['key']}:\n"
            + completed.stderr
            + completed.stdout[-4000:]
        )
    size_match = re.search(r"(?m)^BASIS_SIZE\s+(\d+)\s*$", completed.stdout)
    unit_match = re.search(r"(?m)^IS_UNIT\s+(\d+)\s*$", completed.stdout)
    if size_match is None or unit_match is None:
        raise RuntimeError(
            f"missing Singular markers for {system['key']}:\n"
            + completed.stdout[-4000:]
        )
    basis_size = int(size_match.group(1))
    is_unit = int(unit_match.group(1))
    assert basis_size == 1
    assert is_unit == 1
    return {
        "representative": system["key"],
        "support_orbit": [list(support) for support in system["support_orbit"]],
        "normalization": system["normalization"],
        "localization_factor": str(system["localization_factor"]),
        "second_moment_elimination": system.get("second_moment_elimination"),
        "moments_imposed": system["moments_imposed"],
        "original_charts_covered": system["original_charts_covered"],
        "exact_reduced_basis_over_QQ": ["1"],
    }


def main() -> None:
    singular = shutil.which("Singular")
    assert singular is not None, "the remaining four-weight certificate requires Singular"
    systems = centered_zero_weight_systems() + balanced_systems()
    assert len(systems) == 7
    with ThreadPoolExecutor(max_workers=len(systems)) as executor:
        records = list(
            executor.map(
                lambda system: run_system(singular, system),
                systems,
            )
        )
    assert sum(int(record["original_charts_covered"]) for record in records) == 20

    payload = {
        "format": "two-real-gmc-remaining-four-weight-v1",
        "supports_excluded": [
            [-2, -1, 0, 1],
            [-2, -1, 1, 2],
            [-1, 0, 1, 2],
        ],
        "support_count": 3,
        "input_chart_count": 20,
        "representative_ideal_count": len(records),
        "maximum_moment_order": 6,
        "representative_certificates": records,
        "conclusion": (
            "no total-degree-three pure-moment-zero polynomial has exactly "
            "four nonzero rotational weights containing both signs"
        ),
        "combined_cubic_frontier": (
            "a cubic GMC(2) counterexample must have at least five nonzero "
            "rotational weights"
        ),
        "scope": "the five-, six-, and seven-weight cubic strata remain open",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for record in records:
        print(
            "PASS remaining four-weight GMC:",
            record["representative"],
            "has exact rational basis [1]",
        )
    print(
        "PASS remaining four-weight GMC: excluded 3 supports / 20 charts; wrote",
        OUTPUT.relative_to(ROOT),
    )


if __name__ == "__main__":
    main()
