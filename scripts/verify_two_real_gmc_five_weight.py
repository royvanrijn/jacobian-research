#!/usr/bin/env python3
"""Resumable exact census of the five-weight cubic GMC(2) stratum.

For a cubic in circular coordinates, the possible rotational weights are
``-3,...,3``.  There are 21 five-weight supports and 102 ordinary
nonvanishing coefficient charts.  Reflection ``Z <-> W`` and centering
the weight-zero component reduce these to 12 support representatives and
35 chart representatives.

Moments are generated directly from zero-weight multinomial compositions.
For each chart, computations at ``PRIME`` are used only to discover the
first promising cutoff among moments 8, 10, and 12.  A modular unit ideal is
promoted only by an independent reduced-basis calculation over QQ.  Results
are written one chart at a time, so interrupted runs resume safely.

This script proves nonexistence on a support only when every representative
chart in its reflection orbit has an exact rational basis ``[1]``.  A chart
that remains nonunit modulo the discovery prime is merely a finite-moment
survivor, not an all-order GMC witness.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations, product
from pathlib import Path

import sympy as sp

from verify_two_real_gmc_symmetric_chart import (
    primitive_numerator,
    singular_expression,
)


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_real_gmc_five_weight_charts"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_real_gmc_five_weight.json"
)
SYSTEMS_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_real_gmc_five_weight_systems.json"
)
PRIME = 1_000_003
CUTOFFS = (8, 10, 12)


def reflected_support(support: tuple[int, ...]) -> tuple[int, ...]:
    """Return the support obtained by interchanging Z and W."""

    return tuple(sorted(-weight for weight in support))


def five_weight_support_representatives() -> list[tuple[int, ...]]:
    """Return the 12 reflection representatives of the 21 supports."""

    supports = list(combinations(range(-3, 4), 5))
    assert len(supports) == 21
    representatives = sorted(
        {min(support, reflected_support(support)) for support in supports}
    )
    assert len(representatives) == 12
    return representatives


def coefficient_name(weight: int, radial_power: int) -> str:
    """Use short deterministic names in Singular rings."""

    if weight < 0:
        return f"n{-weight}_{radial_power}"
    if weight > 0:
        return f"p{weight}_{radial_power}"
    return "h0"


def support_data(
    support: tuple[int, ...],
) -> tuple[
    dict[int, tuple[sp.Symbol, ...]],
    tuple[sp.Symbol, ...],
    tuple[tuple[int, int, int, int], ...],
]:
    """Return components, coefficients, and centered monomial terms.

    A term is ``(z_degree, w_degree, sign, coefficient_index)``.  If weight
    zero occurs, its first moment centers ``c_0+c_1*ZW`` to ``h0*(ZW-1)``;
    the two displayed monomials therefore share one coefficient.
    """

    components: dict[int, tuple[sp.Symbol, ...]] = {}
    coefficients: list[sp.Symbol] = []
    raw_terms: list[tuple[int, int, int, sp.Symbol]] = []
    for weight in support:
        if weight == 0:
            coefficient = sp.Symbol("h0")
            components[weight] = (coefficient,)
            coefficients.append(coefficient)
            raw_terms.extend(
                (
                    (0, 0, -1, coefficient),
                    (1, 1, 1, coefficient),
                )
            )
            continue
        radial_degree = (3 - abs(weight)) // 2
        component: list[sp.Symbol] = []
        for radial_power in range(radial_degree + 1):
            coefficient = sp.Symbol(coefficient_name(weight, radial_power))
            component.append(coefficient)
            coefficients.append(coefficient)
            if weight > 0:
                z_degree = weight + radial_power
                w_degree = radial_power
            else:
                z_degree = radial_power
                w_degree = -weight + radial_power
            raw_terms.append(
                (z_degree, w_degree, 1, coefficient)
            )
        components[weight] = tuple(component)
    coefficient_tuple = tuple(coefficients)
    index = {coefficient: position for position, coefficient in enumerate(coefficients)}
    terms = tuple(
        (z_degree, w_degree, sign, index[coefficient])
        for z_degree, w_degree, sign, coefficient in raw_terms
    )
    return components, coefficient_tuple, terms


def moment_from_compositions(
    order: int,
    coefficients: tuple[sp.Symbol, ...],
    terms: tuple[tuple[int, int, int, int], ...],
) -> sp.Expr:
    """Compute one exact moment without expanding powers in Z and W."""

    accumulated: dict[tuple[int, ...], int] = {}
    coefficient_exponents = [0] * len(coefficients)
    order_factorial = math.factorial(order)

    def visit(
        term_index: int,
        remaining: int,
        weight: int,
        z_degree: int,
        denominator: int,
        sign: int,
    ) -> None:
        if term_index == len(terms) - 1:
            count = remaining
            term_z, term_w, term_sign, coefficient_index = terms[term_index]
            final_weight = weight + count * (term_z - term_w)
            if final_weight != 0:
                return
            coefficient_exponents[coefficient_index] += count
            final_z_degree = z_degree + count * term_z
            integer = (
                sign
                * (term_sign**count)
                * order_factorial
                * math.factorial(final_z_degree)
                // (denominator * math.factorial(count))
            )
            key = tuple(coefficient_exponents)
            accumulated[key] = accumulated.get(key, 0) + integer
            coefficient_exponents[coefficient_index] -= count
            return

        term_z, term_w, term_sign, coefficient_index = terms[term_index]
        for count in range(remaining + 1):
            coefficient_exponents[coefficient_index] += count
            visit(
                term_index + 1,
                remaining - count,
                weight + count * (term_z - term_w),
                z_degree + count * term_z,
                denominator * math.factorial(count),
                sign * (term_sign**count),
            )
            coefficient_exponents[coefficient_index] -= count

    visit(0, order, 0, 0, 1, 1)
    accumulated = {
        monomial: coefficient
        for monomial, coefficient in accumulated.items()
        if coefficient != 0
    }
    if not accumulated:
        return sp.Integer(0)
    return sp.Poly.from_dict(
        accumulated,
        *coefficients,
        domain=sp.ZZ,
    ).as_expr()


def generic_moments(
    coefficients: tuple[sp.Symbol, ...],
    terms: tuple[tuple[int, int, int, int], ...],
    maximum: int = 12,
) -> dict[int, sp.Expr]:
    """Generate the centered generic moment list for one support."""

    return {
        order: moment_from_compositions(order, coefficients, terms)
        for order in range(1, maximum + 1)
    }


def chart_systems() -> list[dict[str, object]]:
    """Build the 35 centered/reflected representative systems."""

    systems: list[dict[str, object]] = []
    for support in five_weight_support_representatives():
        components, coefficients, terms = support_data(support)
        moments = generic_moments(coefficients, terms, max(CUTOFFS))
        assert moments[1] == 0
        choices = tuple(range(len(components[weight])) for weight in support)
        for chart in product(*choices):
            selected = tuple(
                components[weight][choice]
                for weight, choice in zip(support, chart)
            )
            # Every chosen representative has at least two negative weights,
            # so these two distinct torus characters give the scale/torus
            # slice used in the previous three- and four-weight censuses.
            assert support[0] < support[1] < 0
            normalization = {selected[0]: 1, selected[1]: 1}
            variables = tuple(
                coefficient
                for coefficient in coefficients
                if coefficient not in normalization
            )
            normalized_moments = {
                order: sp.expand(moment.subs(normalization))
                for order, moment in moments.items()
                if moment != 0
            }
            localization_factor = sp.prod(
                coefficient.subs(normalization)
                for coefficient in selected[2:]
            )
            # A large majority of charts have a second moment that is
            # affine-linear in an outer positive-weight coefficient with a
            # nonzero constant leading coefficient.  Eliminating that
            # variable is globally valid over QQ (no extra chart split) and
            # turns several otherwise expensive rational unit calculations
            # into small ones.
            second_moment_elimination: tuple[sp.Symbol, sp.Expr] | None = None
            for variable in reversed(variables):
                univariate = sp.Poly(normalized_moments[2], variable)
                if univariate.degree() != 1 or univariate.LC().free_symbols:
                    continue
                solution = sp.solve(normalized_moments[2], variable)[0]
                remaining_variables = tuple(
                    candidate
                    for candidate in variables
                    if candidate != variable
                )
                normalized_moments = {
                    order: (
                        sp.Integer(0)
                        if order == 2
                        else primitive_numerator(
                            moment.subs(variable, solution),
                            remaining_variables,
                        )
                    )
                    for order, moment in normalized_moments.items()
                }
                localization_factor = primitive_numerator(
                    localization_factor.subs(variable, solution),
                    remaining_variables,
                )
                variables = remaining_variables
                second_moment_elimination = (variable, solution)
                break
            rho = sp.Symbol("rho")
            reflected = reflected_support(support)
            orbit_size = 1 if reflected == support else 2
            centered_chart_factor = 2 if 0 in support else 1
            key = (
                "s_"
                + "_".join(str(weight).replace("-", "m") for weight in support)
                + "__c_"
                + "_".join(str(choice) for choice in chart)
            )
            systems.append(
                {
                    "key": key,
                    "support": support,
                    "reflected_support": reflected,
                    "support_orbit_size": orbit_size,
                    "chart": chart,
                    "selected_coefficients": selected,
                    "normalization": normalization,
                    "variables": variables + (rho,),
                    "moments": normalized_moments,
                    "localization_factor": localization_factor,
                    "localization_equation": sp.expand(
                        rho * localization_factor - 1
                    ),
                    "second_moment_elimination": second_moment_elimination,
                    "raw_charts_covered": orbit_size * centered_chart_factor,
                }
            )
    assert len(systems) == 35
    assert sum(int(system["raw_charts_covered"]) for system in systems) == 102
    return systems


def system_record(system: dict[str, object]) -> dict[str, object]:
    """Serialize one normalized moment system before Gröbner computation."""

    normalization = system["normalization"]
    moments = system["moments"]
    assert isinstance(normalization, dict)
    assert isinstance(moments, dict)
    return {
        "key": system["key"],
        "support": list(system["support"]),
        "reflected_support": list(system["reflected_support"]),
        "chart": list(system["chart"]),
        "selected_coefficients": [
            str(coefficient) for coefficient in system["selected_coefficients"]
        ],
        "normalization": {
            str(coefficient): int(value)
            for coefficient, value in normalization.items()
        },
        "variables": [str(variable) for variable in system["variables"]],
        "localization_factor": str(system["localization_factor"]),
        "localization_equation": str(system["localization_equation"]),
        "second_moment_elimination": (
            None
            if system["second_moment_elimination"] is None
            else {
                "variable": str(system["second_moment_elimination"][0]),
                "solution": str(system["second_moment_elimination"][1]),
            }
        ),
        "raw_charts_covered": system["raw_charts_covered"],
        "moment_term_counts": {
            str(order): len(
                sp.Poly(moment, *system["variables"][:-1]).terms()
            )
            for order, moment in moments.items()
        },
        "moments": {
            str(order): str(moment)
            for order, moment in moments.items()
            if order <= max(CUTOFFS)
        },
    }


def singular_program(
    system: dict[str, object],
    cutoff: int,
    characteristic: int,
) -> str:
    """Emit one modular or rational reduced-basis calculation."""

    variables = system["variables"]
    moments = system["moments"]
    assert isinstance(variables, tuple)
    assert isinstance(moments, dict)
    generators = [system["localization_equation"]]
    generators.extend(
        moments[order]
        for order in sorted(moments)
        if order <= cutoff and moments[order] != 0
    )
    generator_text = ",\n".join(
        singular_expression(generator) for generator in generators
    )
    variable_text = ",".join(str(variable) for variable in variables)
    return f"""
ring coefficient_ring={characteristic},({variable_text}),dp;
option(redSB);
ideal I={generator_text};
ideal G=slimgb(I);
int is_unit=0;
if (size(G)==1)
{{
  if (G[1]==1) {{ is_unit=1; }}
}}
print("BASIS_SIZE "+string(size(G)));
print("DIMENSION "+string(dim(G)));
print("IS_UNIT "+string(is_unit));
if (dim(G)==0)
{{
  print("VECTOR_DIMENSION "+string(vdim(G)));
}}
"""


def parse_basis_markers(output: str) -> dict[str, int | None]:
    """Parse the small certificate summary emitted by Singular."""

    markers: dict[str, int | None] = {}
    for name in ("BASIS_SIZE", "DIMENSION", "IS_UNIT"):
        match = re.search(rf"(?m)^{name}\s+(-?\d+)\s*$", output)
        if match is None:
            raise RuntimeError(f"Singular did not report {name}:\n{output[-4000:]}")
        markers[name.lower()] = int(match.group(1))
    vector_match = re.search(
        r"(?m)^VECTOR_DIMENSION\s+(-?\d+)\s*$",
        output,
    )
    markers["vector_dimension"] = (
        int(vector_match.group(1)) if vector_match is not None else None
    )
    return markers


def run_basis(
    singular: str,
    system: dict[str, object],
    cutoff: int,
    characteristic: int,
    timeout: int,
) -> dict[str, object]:
    """Run one finite moment ideal and return its basis metadata."""

    started = time.monotonic()
    completed = subprocess.run(
        (singular, "-q"),
        input=singular_program(system, cutoff, characteristic),
        text=True,
        capture_output=True,
        timeout=timeout,
        check=True,
    )
    if "?" in completed.stdout or "?" in completed.stderr:
        raise RuntimeError(
            f"Singular error for {system['key']}:\n"
            + completed.stderr
            + completed.stdout[-4000:]
        )
    result = parse_basis_markers(completed.stdout)
    result.update(
        {
            "cutoff": cutoff,
            "characteristic": characteristic,
            "elapsed_seconds": round(time.monotonic() - started, 3),
        }
    )
    return result


def result_path(system: dict[str, object]) -> Path:
    return RESULT_DIR / f"{system['key']}.json"


def atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    """Replace one result record only after it is completely written."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=path.name + ".",
        suffix=".tmp",
        delete=False,
    ) as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def base_result(system: dict[str, object]) -> dict[str, object]:
    return {
        "format": "two-real-gmc-five-weight-chart-v1",
        "key": system["key"],
        "support": list(system["support"]),
        "reflected_support": list(system["reflected_support"]),
        "chart": list(system["chart"]),
        "selected_coefficients": [
            str(coefficient) for coefficient in system["selected_coefficients"]
        ],
        "normalization": {
            str(coefficient): int(value)
            for coefficient, value in system["normalization"].items()
        },
        "localization_factor": str(system["localization_factor"]),
        "second_moment_elimination": (
            None
            if system["second_moment_elimination"] is None
            else {
                "variable": str(system["second_moment_elimination"][0]),
                "solution": str(system["second_moment_elimination"][1]),
            }
        ),
        "raw_charts_covered": system["raw_charts_covered"],
        "discovery_prime": PRIME,
        "cutoffs": list(CUTOFFS),
    }


def process_system(
    singular: str,
    system: dict[str, object],
    modular_timeout: int,
    rational_timeout: int,
    force: bool,
) -> dict[str, object]:
    """Discover modularly and promote a chart to an exact QQ certificate."""

    path = result_path(system)
    if path.exists() and not force:
        existing = json.loads(path.read_text())
        if existing.get("status") == "exact_unit":
            return existing

    record = base_result(system)
    modular_runs: list[dict[str, object]] = []
    try:
        for cutoff in CUTOFFS:
            modular = run_basis(
                singular,
                system,
                cutoff,
                PRIME,
                modular_timeout,
            )
            modular_runs.append(modular)
            if modular["is_unit"] == 1:
                rational = run_basis(
                    singular,
                    system,
                    cutoff,
                    0,
                    rational_timeout,
                )
                record["modular_discovery"] = modular_runs
                record["rational_promotion"] = rational
                if rational["is_unit"] == 1:
                    record["status"] = "exact_unit"
                    record["exact_certificate"] = {
                        "field": "QQ",
                        "moments_through": cutoff,
                        "reduced_basis": ["1"],
                    }
                else:
                    record["status"] = "bad_prime_discrepancy"
                atomic_write_json(path, record)
                return record
        record["modular_discovery"] = modular_runs
        record["status"] = "finite_moment_survivor_mod_prime"
    except subprocess.TimeoutExpired as error:
        record["modular_discovery"] = modular_runs
        record["status"] = "timeout"
        record["error"] = (
            f"timeout after {error.timeout} seconds at "
            f"{'QQ' if error.cmd and '0' in str(error.cmd) else 'discovery'}"
        )
    except Exception as error:  # retain a resumable diagnostic record
        record["modular_discovery"] = modular_runs
        record["status"] = "error"
        record["error"] = f"{type(error).__name__}: {error}"
    atomic_write_json(path, record)
    return record


def aggregate(systems: list[dict[str, object]]) -> dict[str, object]:
    """Aggregate completed per-chart files without overstating survivors."""

    records = []
    for system in systems:
        path = result_path(system)
        if path.exists():
            records.append(json.loads(path.read_text()))
    by_key = {record["key"]: record for record in records}
    support_records = []
    exactly_excluded_supports = 0
    exactly_excluded_raw_charts = 0
    for support in five_weight_support_representatives():
        relevant = [
            system for system in systems if system["support"] == support
        ]
        completed = [
            by_key[system["key"]]
            for system in relevant
            if system["key"] in by_key
        ]
        all_exact = (
            len(completed) == len(relevant)
            and all(record["status"] == "exact_unit" for record in completed)
        )
        orbit_size = int(relevant[0]["support_orbit_size"])
        if all_exact:
            exactly_excluded_supports += orbit_size
            exactly_excluded_raw_charts += sum(
                int(record["raw_charts_covered"]) for record in completed
            )
        support_records.append(
            {
                "representative": list(support),
                "reflected_support": list(reflected_support(support)),
                "representative_charts": len(relevant),
                "completed_chart_records": len(completed),
                "status_histogram": {
                    status: sum(
                        record["status"] == status for record in completed
                    )
                    for status in sorted(
                        {record["status"] for record in completed}
                    )
                },
                "exactly_excluded": all_exact,
            }
        )
    status_histogram = {
        status: sum(record["status"] == status for record in records)
        for status in sorted({record["status"] for record in records})
    }
    return {
        "format": "two-real-gmc-five-weight-v1",
        "scope": (
            "all 21 exactly-five-weight cubic supports; six- and seven-weight "
            "cubics are outside this artifact"
        ),
        "ordinary_supports": 21,
        "ordinary_nonvanishing_charts": 102,
        "reflection_support_representatives": 12,
        "centered_reflection_chart_representatives": 35,
        "completed_chart_records": len(records),
        "status_histogram": status_histogram,
        "exactly_excluded_supports": exactly_excluded_supports,
        "exactly_excluded_ordinary_charts": exactly_excluded_raw_charts,
        "support_records": support_records,
        "interpretation": (
            "exact_unit is a characteristic-zero exclusion; "
            "finite_moment_survivor_mod_prime is only a discovery result and "
            "is not an all-order witness"
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--modular-timeout", type=int, default=600)
    parser.add_argument("--rational-timeout", type=int, default=1800)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--only-support",
        action="append",
        default=[],
        help="comma-separated support, for example -2,-1,0,1,2",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    singular = shutil.which("Singular")
    assert singular is not None, "this census requires Singular"
    systems = chart_systems()
    systems_payload = {
        "format": "two-real-gmc-five-weight-systems-v1",
        "ordinary_supports": 21,
        "ordinary_nonvanishing_charts": 102,
        "reflection_support_representatives": 12,
        "centered_reflection_chart_representatives": 35,
        "systems": [system_record(system) for system in systems],
    }
    atomic_write_json(SYSTEMS_OUTPUT, systems_payload)

    if args.only_support:
        requested = {
            tuple(int(part) for part in text.split(","))
            for text in args.only_support
        }
        systems = [
            system for system in systems if system["support"] in requested
        ]
        if not systems:
            raise SystemExit("no representative system matches --only-support")

    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                process_system,
                singular,
                system,
                args.modular_timeout,
                args.rational_timeout,
                args.force,
            ): system
            for system in systems
        }
        for future in as_completed(futures):
            system = futures[future]
            record = future.result()
            print(
                record["status"],
                system["support"],
                system["chart"],
                flush=True,
            )

    all_systems = chart_systems()
    payload = aggregate(all_systems)
    atomic_write_json(OUTPUT, payload)
    print(
        "SUMMARY",
        payload["status_histogram"],
        "exact supports",
        payload["exactly_excluded_supports"],
        "/ 21; wrote",
        OUTPUT.relative_to(ROOT),
    )


if __name__ == "__main__":
    main()
