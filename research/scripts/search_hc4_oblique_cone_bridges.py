#!/usr/bin/env python3
"""Exact HC(4) search on non-coordinate cone/bridge families.

For a nonzero rational slope lambda, put

    u = x2 + lambda*x3.

The degree-d correction is chosen in Q[x0,x1,u], while the fixed collision
carrier is -x1*x0^(d-1).  Its four-variable Hessian determinant therefore
vanishes identically.  Each top parameter expands to a tied collection of
ordinary x2/x3 monomials.  Lower-degree collision-invisible monomials
involving x3 provide complementary fourth-variable bridges.

For every selected support, the script forms the complete coefficient ideal
of det(Hess(psi))-1 and solves it exactly over the requested prime fields
using Singular.  The support sampling is bounded and is not a proof of
HC(4).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from fractions import Fraction
from pathlib import Path
from time import perf_counter

try:
    from search_hc4_finite_field_dense_supports import (
        coefficient_system,
        deterministic_order,
        singular_equations,
        solve_support,
    )
    from search_hc4_finite_field_potentials import (
        Direction,
        direction_record,
        is_prime,
        weak_compositions,
    )
except ModuleNotFoundError:
    from scripts.search_hc4_finite_field_dense_supports import (
        coefficient_system,
        deterministic_order,
        singular_equations,
        solve_support,
    )
    from scripts.search_hc4_finite_field_potentials import (
        Direction,
        direction_record,
        is_prime,
        weak_compositions,
    )


def oblique_top_directions(
    degree_bound: int, slope: Fraction
) -> list[Direction]:
    """Collision-invisible degree-d monomials in adapted x0,x1,u."""

    del slope
    directions: list[Direction] = []
    visible_seeds = {
        (degree_bound, 0, 0),
        (degree_bound - 1, 1, 0),
        (degree_bound - 1, 0, 1),
    }
    for exponent_x0, exponent_x1, exponent_u in weak_compositions(
        degree_bound, 3
    ):
        seed = (exponent_x0, exponent_x1, exponent_u)
        if seed in visible_seeds:
            continue
        directions.append(
            (
                (
                    (
                        exponent_x0,
                        exponent_x1,
                        exponent_u,
                        0,
                    ),
                    Fraction(1),
                ),
            )
        )
    return directions


def complementary_bridge_directions(
    degree_bound: int,
) -> list[Direction]:
    """Lower-degree monomials using x3 and preserving the normalized collision."""

    directions: list[Direction] = []
    for degree in range(3, degree_bound):
        visible_x3 = (degree - 1, 0, 0, 1)
        for exponent in weak_compositions(degree, 4):
            if exponent[3] == 0 or exponent == visible_x3:
                continue
            directions.append(((exponent, Fraction(1)),))
    return directions


def family_directions(
    degree_bound: int, slope: Fraction
) -> tuple[list[Direction], int]:
    top = oblique_top_directions(degree_bound, slope)
    bridges = complementary_bridge_directions(degree_bound)
    return top + bridges, len(top)


def oblique_base_terms(
    degree_bound: int, slope: Fraction
) -> list[tuple[tuple[int, int, int, int], Fraction]]:
    """Base in adapted coordinates (x0,x1,u,v), with x2=u-slope*v."""

    return [
        ((1, 1, 0, 0), Fraction(1)),
        ((0, 0, 1, 1), Fraction(1)),
        ((0, 0, 0, 2), -slope),
        ((degree_bound - 1, 1, 0, 0), Fraction(-1)),
    ]


def choose_oblique_support(
    degree_bound: int,
    slope: Fraction,
    support_size: int,
    trial: int,
    seed: str,
    top_count_available: int,
    total_direction_count: int,
) -> tuple[int, ...]:
    top_count = support_size // 2
    bridge_count = support_size - top_count
    if top_count_available < top_count:
        raise ValueError("not enough oblique top directions")
    bridge_indices = list(
        range(top_count_available, total_direction_count)
    )
    if len(bridge_indices) < bridge_count:
        raise ValueError("not enough complementary bridge directions")
    slope_text = (
        str(slope.numerator)
        if slope.denominator == 1
        else f"{slope.numerator}/{slope.denominator}"
    )
    key = (
        f"{seed}:d={degree_bound}:slope={slope_text}:"
        f"k={support_size}:trial={trial}"
    )
    chosen = deterministic_order(
        list(range(top_count_available)), key + ":top"
    )[:top_count]
    chosen += deterministic_order(
        bridge_indices, key + ":bridge"
    )[:bridge_count]
    assert len(chosen) == support_size
    return tuple(sorted(chosen))


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def parse_fraction(value: str) -> Fraction:
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError(
            f"invalid rational slope: {value}"
        ) from error


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--degrees", type=int, nargs="+", default=[5, 6, 7, 8]
    )
    parser.add_argument(
        "--primes", type=int, nargs="+", default=[11, 13]
    )
    parser.add_argument(
        "--support-sizes",
        type=int,
        nargs="+",
        default=[6, 8, 10, 12],
    )
    parser.add_argument(
        "--slopes",
        type=parse_fraction,
        nargs="+",
        default=[Fraction(-1), Fraction(1), Fraction(2)],
    )
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--seed", default="hc4-oblique-v1")
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    if arguments.trials < 1:
        raise SystemExit("--trials must be positive")
    if arguments.timeout <= 0:
        raise SystemExit("--timeout must be positive")
    if any(degree < 3 for degree in arguments.degrees):
        raise SystemExit("every degree must be at least 3")
    if any(size < 2 or size % 2 for size in arguments.support_sizes):
        raise SystemExit("support sizes must be positive even integers")
    if any(slope == 0 for slope in arguments.slopes):
        raise SystemExit("oblique slopes must be nonzero")
    if any(not is_prime(prime) for prime in arguments.primes):
        raise SystemExit("--primes accepts prime numbers only")
    if any(
        prime <= max(arguments.degrees) for prime in arguments.primes
    ):
        raise SystemExit("every prime must exceed every searched degree")
    for slope in arguments.slopes:
        if any(
            slope.denominator % prime == 0
            for prime in arguments.primes
        ):
            raise SystemExit(
                "slope denominators must be invertible modulo every prime"
            )
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required")

    records: list[dict[str, object]] = []
    total_started = perf_counter()
    for degree_bound in arguments.degrees:
        for slope in arguments.slopes:
            directions, top_direction_count = family_directions(
                degree_bound, slope
            )
            seen_supports: set[tuple[int, ...]] = set()
            for support_size in arguments.support_sizes:
                for trial in range(arguments.trials):
                    support = choose_oblique_support(
                        degree_bound,
                        slope,
                        support_size,
                        trial,
                        arguments.seed,
                        top_direction_count,
                        len(directions),
                    )
                    if support in seen_supports:
                        continue
                    seen_supports.add(support)
                    expansion_started = perf_counter()
                    (
                        coefficient_equations,
                        determinant_term_count,
                        coefficient_equation_count,
                    ) = coefficient_system(
                        degree_bound,
                        directions,
                        support,
                        base_potential_terms=oblique_base_terms(
                            degree_bound, slope
                        ),
                    )
                    expansion_seconds = round(
                        perf_counter() - expansion_started, 6
                    )
                    for prime in arguments.primes:
                        modular_equations = singular_equations(
                            coefficient_equations, prime
                        )
                        solver_record = solve_support(
                            singular,
                            prime,
                            support_size,
                            modular_equations,
                            arguments.timeout,
                        )
                        top_selected = sum(
                            index < top_direction_count
                            for index in support
                        )
                        record: dict[str, object] = {
                            "degree_bound": degree_bound,
                            "slope": fraction_text(slope),
                            "prime": prime,
                            "support_size": support_size,
                            "trial": trial,
                            "support_indices": list(support),
                            "top_direction_count": top_selected,
                            "bridge_direction_count": (
                                support_size - top_selected
                            ),
                            "determinant_term_count": (
                                determinant_term_count
                            ),
                            "coefficient_equation_count": (
                                coefficient_equation_count
                            ),
                            "modular_equation_count": len(
                                modular_equations
                            ),
                            **solver_record,
                        }
                        if solver_record["solver_status"] == "nonunit":
                            record["support_directions"] = [
                                direction_record(directions[index])
                                for index in support
                            ]
                        records.append(record)
                        print(
                            "HC4_OBLIQUE"
                            f" d={degree_bound}"
                            f" slope={fraction_text(slope)}"
                            f" p={prime}"
                            f" k={support_size}"
                            f" trial={trial}"
                            f" terms={determinant_term_count}"
                            f" equations={len(modular_equations)}"
                            f" expansion_s={expansion_seconds}"
                            f" solver={solver_record['solver_status']}"
                            f" solver_s={solver_record['solver_seconds']}"
                        )

    status_counts: dict[str, int] = {}
    for record in records:
        status = str(record["solver_status"])
        status_counts[status] = status_counts.get(status, 0) + 1
    payload = {
        "status": (
            "bounded oblique-cone finite-field experiment; not a proof"
        ),
        "seed": arguments.seed,
        "degrees": arguments.degrees,
        "primes": arguments.primes,
        "support_sizes": arguments.support_sizes,
        "slopes": [fraction_text(slope) for slope in arguments.slopes],
        "trials": arguments.trials,
        "timeout_seconds": arguments.timeout,
        "normalization": {
            "quadratic_part": "x0*x1 + x2*x3",
            "collision": "grad(psi)(0) = grad(psi)(1,0,0,0)",
            "constant_hessian_determinant": 1,
            "top_cone": "u=x2+slope*x3",
            "adapted_coordinates": (
                "(x0,x1,u,v)=(x0,x1,x2+slope*x3,x3)"
            ),
        },
        "records": records,
        "status_counts": status_counts,
    }
    deterministic_payload = json.loads(json.dumps(payload))
    for record in deterministic_payload["records"]:
        record.pop("solver_seconds", None)
    canonical = json.dumps(
        deterministic_payload, sort_keys=True, separators=(",", ":")
    )
    payload["deterministic_content_sha256"] = hashlib.sha256(
        canonical.encode()
    ).hexdigest()
    print(
        "HC4_OBLIQUE_SUMMARY"
        f" records={len(records)}"
        f" statuses={json.dumps(status_counts, sort_keys=True)}"
        f" total_seconds={round(perf_counter() - total_started, 6)}"
    )
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )
        print(f"WROTE {arguments.output}")


if __name__ == "__main__":
    main()
