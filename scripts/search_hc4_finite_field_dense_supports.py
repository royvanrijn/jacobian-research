#!/usr/bin/env python3
"""Exact modular coefficient solving on sampled dense HC(4) supports.

This experiment continues ``search_hc4_finite_field_potentials.py`` beyond
its exhaustive one/two-direction range.  It deterministically chooses
supports of 6--12 vectors in the full collision-kernel basis, forms every
coefficient of

    det(Hess(psi)) - 1,

and asks Singular whether those coefficient equations have a solution over
GF(p).  Field equations ``a_i^p-a_i`` are included, so a unit Groebner basis
is an exact exclusion of the selected support over GF(p), while a nonunit
basis is an exact modular candidate-support certificate.

The support sampling is bounded and is not evidence for unrestricted HC(4).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from fractions import Fraction
from pathlib import Path
from time import perf_counter

import sympy as sp

try:
    from search_hc4_finite_field_potentials import (
        Direction,
        base_terms,
        collision_kernel_directions,
        direction_record,
        is_prime,
    )
except ModuleNotFoundError:
    from scripts.search_hc4_finite_field_potentials import (
        Direction,
        base_terms,
        collision_kernel_directions,
        direction_record,
        is_prime,
    )


def monomial(
    variables: tuple[sp.Symbol, ...], exponents: tuple[int, ...]
) -> sp.Expr:
    return sp.prod(
        variable**exponent
        for variable, exponent in zip(variables, exponents, strict=True)
    )


def deterministic_order(
    population: list[int], key: str
) -> list[int]:
    return sorted(
        population,
        key=lambda index: hashlib.sha256(
            f"{key}:{index}".encode()
        ).digest(),
    )


_AXIS_ELIGIBLE_CACHE: dict[int, list[int]] = {}


def axis_eligible_indices(
    directions: list[Direction], degree_bound: int
) -> list[int]:
    """Directions that can alter both forced base defects on the x0-axis."""

    if degree_bound in _AXIS_ELIGIBLE_CACHE:
        return _AXIS_ELIGIBLE_CACHE[degree_bound]
    x_variables = sp.symbols("x0:4")
    parameter = sp.Symbol("axis_parameter")
    base = sum(
        sp.Rational(coefficient.numerator, coefficient.denominator)
        * monomial(x_variables, exponent)
        for exponent, coefficient in base_terms(degree_bound)
    )
    target_powers = {degree_bound - 2, 2 * degree_bound - 4}
    substitutions = {
        x_variables[1]: 0,
        x_variables[2]: 0,
        x_variables[3]: 0,
    }
    eligible: list[int] = []
    for index, direction in enumerate(directions):
        potential = base + sum(
            parameter
            * sp.Rational(
                coefficient.numerator, coefficient.denominator
            )
            * monomial(x_variables, exponent)
            for exponent, coefficient in direction
        )
        restricted_defect = sp.Poly(
            sp.expand(
                sp.hessian(potential, x_variables)
                .subs(substitutions)
                .det(method="berkowitz")
                - 1
            ),
            x_variables[0],
            parameter,
        )
        affected = {
            exponents[0]
            for exponents, _ in restricted_defect.terms()
            if exponents[0] in target_powers and exponents[1] > 0
        }
        if affected == target_powers:
            eligible.append(index)
    _AXIS_ELIGIBLE_CACHE[degree_bound] = eligible
    return eligible


def choose_support(
    directions: list[Direction],
    degree_bound: int,
    support_size: int,
    strategy: str,
    trial: int,
    seed: str,
) -> tuple[int, ...]:
    all_indices = list(range(len(directions)))
    homogeneous = [
        index
        for index, direction in enumerate(directions)
        if len(direction) == 1
        and sum(direction[0][0]) == degree_bound
    ]
    carriers = [
        index
        for index, direction in enumerate(directions)
        if len(direction) == 2
    ]
    singletons = [
        index
        for index, direction in enumerate(directions)
        if len(direction) == 1
    ]
    key = (
        f"{seed}:d={degree_bound}:k={support_size}:"
        f"strategy={strategy}:trial={trial}"
    )
    if strategy == "uniform":
        chosen = deterministic_order(all_indices, key)[:support_size]
    elif strategy == "homogeneous":
        if len(homogeneous) < support_size:
            raise ValueError(
                f"degree {degree_bound} has only {len(homogeneous)} "
                "homogeneous directions"
            )
        chosen = deterministic_order(homogeneous, key)[:support_size]
    elif strategy == "mixed":
        carrier_count = min(
            len(carriers), max(1, support_size // 3)
        )
        chosen = deterministic_order(carriers, key + ":carrier")[
            :carrier_count
        ]
        chosen += deterministic_order(singletons, key + ":singleton")[
            : support_size - carrier_count
        ]
    elif strategy == "axis":
        axis_eligible = axis_eligible_indices(
            directions, degree_bound
        )
        axis_count = min(len(axis_eligible), support_size)
        chosen = deterministic_order(axis_eligible, key + ":axis")[
            :axis_count
        ]
        remaining = [
            index for index in all_indices if index not in chosen
        ]
        chosen += deterministic_order(remaining, key + ":fill")[
            : support_size - axis_count
        ]
    elif strategy in ("cone2", "cone3"):
        omitted_coordinate = int(strategy[-1])
        top_cone = [
            index
            for index, direction in enumerate(directions)
            if len(direction) == 1
            and sum(direction[0][0]) == degree_bound
            and direction[0][0][omitted_coordinate] == 0
        ]
        lower_bridges = [
            index
            for index, direction in enumerate(directions)
            if len(direction) == 1
            and sum(direction[0][0]) < degree_bound
            and direction[0][0][omitted_coordinate] > 0
        ]
        top_count = max(3, support_size // 2)
        top_count = min(
            top_count,
            len(top_cone),
            support_size - 1,
        )
        bridge_count = support_size - top_count
        if len(lower_bridges) < bridge_count:
            raise ValueError(
                f"degree {degree_bound} has only "
                f"{len(lower_bridges)} lower {strategy} bridges"
            )
        chosen = deterministic_order(top_cone, key + ":top")[
            :top_count
        ]
        chosen += deterministic_order(
            lower_bridges, key + ":bridge"
        )[:bridge_count]
    else:
        raise ValueError(f"unknown support strategy: {strategy}")
    assert len(chosen) == support_size
    assert len(set(chosen)) == support_size
    return tuple(sorted(chosen))


def coefficient_system(
    degree_bound: int,
    directions: list[Direction],
    support: tuple[int, ...],
    base_potential_terms: list[tuple[tuple[int, int, int, int], Fraction]]
    | None = None,
) -> tuple[
    list[list[tuple[tuple[int, ...], Fraction]]],
    int,
    int,
]:
    """Return x-coefficient equations as sparse parameter polynomials."""

    x_variables = sp.symbols("x0:4")
    parameters = sp.symbols(f"a0:{len(support)}")
    potential = sp.Integer(0)
    if base_potential_terms is None:
        base_potential_terms = base_terms(degree_bound)
    for exponent, coefficient in base_potential_terms:
        potential += (
            sp.Rational(coefficient.numerator, coefficient.denominator)
            * monomial(x_variables, exponent)
        )
    for parameter, direction_index in zip(
        parameters, support, strict=True
    ):
        for exponent, coefficient in directions[direction_index]:
            potential += (
                parameter
                * sp.Rational(
                    coefficient.numerator, coefficient.denominator
                )
                * monomial(x_variables, exponent)
            )

    determinant_defect = sp.expand(
        sp.hessian(potential, x_variables).det(method="berkowitz") - 1
    )
    polynomial = sp.Poly(
        determinant_defect,
        *x_variables,
        *parameters,
        domain=sp.QQ,
    )
    grouped: dict[
        tuple[int, ...],
        list[tuple[tuple[int, ...], Fraction]],
    ] = {}
    for exponents, coefficient in polynomial.terms():
        grouped.setdefault(exponents[:4], []).append(
            (
                exponents[4:],
                Fraction(
                    int(coefficient.p),
                    int(coefficient.q),
                ),
            )
        )
    return (
        [grouped[key] for key in sorted(grouped)],
        len(polynomial.terms()),
        len(grouped),
    )


def coefficient_mod(coefficient: Fraction, prime: int) -> int:
    return (
        coefficient.numerator
        * pow(coefficient.denominator, -1, prime)
        % prime
    )


def singular_parameter_term(
    exponents: tuple[int, ...],
    coefficient: int,
) -> str:
    factors: list[str] = []
    if coefficient != 1 or not any(exponents):
        factors.append(str(coefficient))
    for index, exponent in enumerate(exponents):
        if exponent == 1:
            factors.append(f"a{index}")
        elif exponent > 1:
            factors.append(f"a{index}^{exponent}")
    return "*".join(factors) if factors else "1"


def singular_equations(
    equations: list[list[tuple[tuple[int, ...], Fraction]]],
    prime: int,
) -> list[str]:
    result: list[str] = []
    for equation in equations:
        terms = [
            singular_parameter_term(
                exponents, coefficient_mod(coefficient, prime)
            )
            for exponents, coefficient in equation
            if coefficient_mod(coefficient, prime) != 0
        ]
        if terms:
            result.append("+".join(terms))
    return result


def singular_program(
    prime: int,
    parameter_count: int,
    equations: list[str],
) -> str:
    parameters = ",".join(
        f"a{index}" for index in range(parameter_count)
    )
    field_equations = [
        f"a{index}^{prime}-a{index}"
        for index in range(parameter_count)
    ]
    ideal_entries = equations + field_equations
    return "\n".join(
        [
            f"ring R={prime},({parameters}),dp;",
            "option(redSB);",
            f"ideal I={','.join(ideal_entries)};",
            "ideal G=std(I);",
            'print("HC4_BASIS_SIZE="+string(size(G)));',
            'print("HC4_FIRST="+string(G[1]));',
            'print("HC4_VDIM="+string(vdim(G)));',
            "if(G[1]<>1)",
            "{",
            '  print("HC4_GROEBNER_BEGIN");',
            "  G;",
            '  print("HC4_GROEBNER_END");',
            "}",
            "quit;",
        ]
    )


def parse_marker(output: str, marker: str) -> str | None:
    prefix = marker + "="
    for line in output.splitlines():
        if line.startswith(prefix):
            return line[len(prefix) :].strip()
    return None


def solve_support(
    singular: str,
    prime: int,
    support_size: int,
    equations: list[str],
    timeout_seconds: float,
) -> dict[str, object]:
    program = singular_program(
        prime, support_size, equations
    )
    started = perf_counter()
    try:
        process = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {
            "solver_status": "timeout",
            "solver_seconds": round(perf_counter() - started, 6),
        }
    elapsed = round(perf_counter() - started, 6)
    first = parse_marker(process.stdout, "HC4_FIRST")
    basis_size = parse_marker(process.stdout, "HC4_BASIS_SIZE")
    vector_dimension = parse_marker(process.stdout, "HC4_VDIM")
    if process.returncode != 0 or first is None:
        return {
            "solver_status": "error",
            "solver_seconds": elapsed,
            "return_code": process.returncode,
            "stderr_tail": process.stderr[-2000:],
            "stdout_tail": process.stdout[-2000:],
        }
    record: dict[str, object] = {
        "solver_status": "unit" if first == "1" else "nonunit",
        "solver_seconds": elapsed,
        "basis_size": int(basis_size) if basis_size else None,
        "quotient_dimension": (
            int(vector_dimension) if vector_dimension else None
        ),
    }
    if first != "1":
        start_marker = "HC4_GROEBNER_BEGIN\n"
        end_marker = "\nHC4_GROEBNER_END"
        if start_marker in process.stdout and end_marker in process.stdout:
            record["groebner_basis"] = process.stdout.split(
                start_marker, 1
            )[1].split(end_marker, 1)[0]
    return record


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
        "--strategies",
        choices=(
            "uniform",
            "homogeneous",
            "mixed",
            "axis",
            "cone2",
            "cone3",
        ),
        nargs="+",
        default=["uniform", "homogeneous", "mixed"],
    )
    parser.add_argument("--trials", type=int, default=2)
    parser.add_argument("--seed", default="hc4-dense-v1")
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
    if any(size < 3 for size in arguments.support_sizes):
        raise SystemExit("every support size must be at least 3")
    if any(not is_prime(prime) for prime in arguments.primes):
        raise SystemExit("--primes accepts prime numbers only")
    if any(
        prime <= max(arguments.degrees) for prime in arguments.primes
    ):
        raise SystemExit("every prime must exceed every searched degree")
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required")

    records: list[dict[str, object]] = []
    total_started = perf_counter()
    for degree_bound in arguments.degrees:
        directions = collision_kernel_directions(degree_bound)
        seen_supports: set[tuple[int, ...]] = set()
        for support_size in arguments.support_sizes:
            for strategy in arguments.strategies:
                for trial in range(arguments.trials):
                    support = choose_support(
                        directions,
                        degree_bound,
                        support_size,
                        strategy,
                        trial,
                        arguments.seed,
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
                        degree_bound, directions, support
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
                        record: dict[str, object] = {
                            "degree_bound": degree_bound,
                            "prime": prime,
                            "support_size": support_size,
                            "strategy": strategy,
                            "trial": trial,
                            "support_indices": list(support),
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
                            "HC4_DENSE"
                            f" d={degree_bound}"
                            f" p={prime}"
                            f" k={support_size}"
                            f" strategy={strategy}"
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
            "bounded sampled-support finite-field experiment; not a proof"
        ),
        "seed": arguments.seed,
        "degrees": arguments.degrees,
        "primes": arguments.primes,
        "support_sizes": arguments.support_sizes,
        "strategies": arguments.strategies,
        "trials": arguments.trials,
        "timeout_seconds": arguments.timeout,
        "normalization": {
            "quadratic_part": "x0*x1 + x2*x3",
            "collision": "grad(psi)(0) = grad(psi)(1,0,0,0)",
            "constant_hessian_determinant": 1,
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
        "HC4_DENSE_SUMMARY"
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
