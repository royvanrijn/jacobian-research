#!/usr/bin/env python3
"""Degree-eight search for four coherent harmonic channels.

The family is

    P = sum_{i=0}^3 a_i rho^((8-ell_i)/2) <v_i,z>^ell_i,
    (ell_0,ell_1,ell_2,ell_3) = (2,4,6,8),

with isotropic directions ``v_i``.  Pairwise-distinct directions are
normalized to infinity, zero, one, and ``lam``.  The first nonzero moment is
linear in ``a2``.  Its coefficient ``B`` gives two exact charts: ``B != 0``
is saturated in the compact, uneliminated moment ideal, while ``B = 0``
forces an explicit value of ``a3`` and is checked after substitution.

All thirteen nonterminal direction-collision partitions are also checked.
They and the ``B=0`` boundary are certified exactly.  The four-distinct
``B!=0`` chart is currently a three-prime modular unit calculation; it is not
promoted to a characteristic-zero theorem unless ``--exact-open`` succeeds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from research_gvc3_many_coherent_channels import (
    exact_unit,
    modular_cutoff,
    moment,
    primitive_polynomial,
    singular_expression,
)

OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "gvc3_four_coherent_channels.json"
)
DEGREES = (2, 4, 6, 8)
DISTINCT_GROUPS = ((0,), (1,), (2,), (3,))
PRIMES = (101, 103, 107)

A0, A1, A2, A3 = sp.symbols("a0:4")
COEFFICIENTS = (A0, A1, A2, A3)
LAM = sp.Symbol("lam")


def set_partitions(size: int):
    answer: list[tuple[tuple[int, ...], ...]] = []

    def visit(index: int, blocks: list[list[int]]) -> None:
        if index == size:
            answer.append(tuple(tuple(block) for block in blocks))
            return
        for block in blocks:
            block.append(index)
            visit(index + 1, blocks)
            block.pop()
        blocks.append([index])
        visit(index + 1, blocks)
        blocks.pop()

    visit(0, [])
    return tuple(answer)


def compile_moments(
    groups: tuple[tuple[int, ...], ...],
    max_order: int,
    parameters: tuple[sp.Symbol, ...],
    variables: tuple[sp.Symbol, ...],
) -> dict[int, sp.Expr]:
    equations: dict[int, sp.Expr] = {}
    for order in range(2, max_order + 1):
        expression = moment(DEGREES, groups, order, COEFFICIENTS, parameters)
        polynomial = primitive_polynomial(expression.subs(A0, 1), variables)
        if polynomial != 0:
            equations[order] = polynomial
    return equations


def strip_inverted_factors(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    factors: tuple[sp.Expr, ...],
) -> sp.Expr:
    answer = primitive_polynomial(expression, variables)
    for factor in factors:
        divisor = sp.Poly(factor, *variables)
        while sp.rem(sp.Poly(answer, *variables), divisor) == 0:
            answer = sp.cancel(answer / factor)
    return primitive_polynomial(answer, variables)


def wrapped_exact_expression(expression: sp.Expr) -> str:
    return (
        str(sp.expand(expression))
        .replace("**", "^")
        .replace(" + ", "\n+")
        .replace(" - ", "\n-")
    )


def singular_saturation(
    equations: dict[int, sp.Expr],
    saturation: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    characteristic: int,
    timeout: int,
) -> dict[str, object]:
    declarations: list[str] = []
    names: list[str] = []
    for order, equation in equations.items():
        name = f"f{order}"
        names.append(name)
        serialized = (
            wrapped_exact_expression(equation)
            if characteristic == 0
            else singular_expression(equation, variables, characteristic)
        )
        declarations.append(f"poly {name}=\n{serialized};")
    serialized_saturation = (
        wrapped_exact_expression(saturation)
        if characteristic == 0
        else singular_expression(saturation, variables, characteristic)
    )
    source = "\n".join(
        [
            'LIB "elim.lib";',
            "option(redSB);",
            f"ring R={characteristic},({','.join(map(str, variables))}),dp;",
            *declarations,
            f"poly jsat=\n{serialized_saturation};",
            f"ideal I={','.join(names)};",
            "ideal J=jsat;",
            "timer=1;",
            "list L=sat_with_exp(I,J);",
            "ideal G=L[1];",
            "int saturation_exponent=L[2];",
            "int elapsed=timer;",
            "int unit=0;",
            "if ((size(G)==1)&&(G[1]==1)) { unit=1; }",
            'print("RESULT_BEGIN");',
            'print("unit="+string(unit));',
            'print("basis_size="+string(size(G)));',
            'print("dimension="+string(dim(G)));',
            'print("saturation_exponent="+string(saturation_exponent));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("BASIS_BEGIN"); print(G); print("BASIS_END");',
            'print("RESULT_END");',
            "quit;",
        ]
    )
    started = time.monotonic()
    record: dict[str, object] = {
        "characteristic": characteristic,
        "input_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "input_bytes": len(source.encode()),
    }
    try:
        completed = subprocess.run(
            ["Singular", "-q"],
            input=source,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        record.update(
            {
                "status": "timeout",
                "seconds": time.monotonic() - started,
                "stdout_tail": (error.stdout or "")[-1000:],
                "stderr": (error.stderr or "")[-1000:],
            }
        )
        return record
    record["seconds"] = time.monotonic() - started
    record["returncode"] = completed.returncode
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in (
            "unit",
            "basis_size",
            "dimension",
            "saturation_exponent",
            "elapsed_ticks",
        ):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    if "BASIS_BEGIN" in completed.stdout and "BASIS_END" in completed.stdout:
        basis = (
            completed.stdout.split("BASIS_BEGIN", 1)[1]
            .split("BASIS_END", 1)[0]
            .strip()
        )
        record["basis_sha256"] = hashlib.sha256(basis.encode()).hexdigest()
        record["basis"] = basis
    record["output_sha256"] = hashlib.sha256(completed.stdout.encode()).hexdigest()
    record["stdout_tail"] = completed.stdout[-1000:]
    record["stderr"] = completed.stderr[-1000:]
    return record


def distinct_direction_data(exact_open: bool) -> dict[str, object]:
    variables = (A1, A2, A3, LAM)
    equations = compile_moments(DISTINCT_GROUPS, 7, (LAM,), variables)
    first = sp.Poly(equations[3], A2)
    if first.degree() != 1:
        raise AssertionError("the third moment must be linear in a2")
    coefficient = sp.expand(first.coeff_monomial(A2))
    remainder = sp.expand(first.coeff_monomial(1))

    open_saturation = A1 * A2 * A3 * LAM * (LAM - 1) * coefficient
    with ThreadPoolExecutor(max_workers=len(PRIMES)) as executor:
        modular = list(
            executor.map(
                lambda prime: singular_saturation(
                    equations, open_saturation, variables, prime, 600
                ),
                PRIMES,
            )
        )
    for result in modular:
        assert result.get("status") == "completed" and result.get("unit") == 1
    exact = None
    if exact_open:
        exact = singular_saturation(
            equations, open_saturation, variables, 0, 1200
        )
        if exact.get("status") != "completed" or exact.get("unit") != 1:
            raise RuntimeError(
                "the optional characteristic-zero B!=0 saturation failed; "
                "do not promote the modular chart"
            )

    # On B=0 the equation A+B*a2=0 forces A=0.  Since a1 is
    # inverted, the remaining linear factor of A solves for a3.
    reduced_remainder = sp.cancel(remainder / A1)
    remainder_in_a3 = sp.Poly(reduced_remainder, A3)
    if remainder_in_a3.degree() != 1:
        raise AssertionError("the B=0 boundary must solve linearly for a3")
    a3_solution = sp.cancel(
        -remainder_in_a3.coeff_monomial(1)
        / remainder_in_a3.coeff_monomial(A3)
    )
    boundary_variables = (A1, A2, LAM)
    boundary_equations: dict[int, sp.Expr] = {}
    raw_boundary = {3: coefficient}
    raw_boundary.update(
        {order: equation for order, equation in equations.items() if 4 <= order <= 6}
    )
    for order, equation in raw_boundary.items():
        numerator = sp.together(equation.subs(A3, a3_solution)).as_numer_denom()[0]
        boundary_equations[order] = strip_inverted_factors(
            numerator,
            boundary_variables,
            (A1, A2, LAM, LAM - 1),
        )
    boundary_saturation = A1 * A2 * LAM * (LAM - 1)
    boundary_modular = [
        modular_cutoff(
            boundary_equations,
            boundary_saturation,
            boundary_variables,
            prime,
            "Singular",
            120,
        )
        for prime in PRIMES
    ]
    assert all(result.get("unit") == 1 for result in boundary_modular)
    boundary_exact = exact_unit(
        boundary_equations,
        boundary_saturation,
        boundary_variables,
        6,
        300,
        4,
    )
    assert boundary_exact["unit"] == 1

    return {
        "normalization": "directions infinity, zero, one, lam; a0=1",
        "configuration_saturation": "lam*(lam-1)",
        "moment_term_counts": {
            str(order): len(sp.Poly(equation, *variables).terms())
            for order, equation in equations.items()
        },
        "moment_sha256": {
            str(order): hashlib.sha256(str(equation).encode()).hexdigest()
            for order, equation in equations.items()
        },
        "third_moment": "A+B*a2",
        "A_terms": len(sp.Poly(remainder, *variables).terms()),
        "B_terms": len(sp.Poly(coefficient, *variables).terms()),
        "B_nonzero_chart": {
            "cutoff": 7,
            "modular_saturations": modular,
            "exact_Q_saturation": exact,
        },
        "B_zero_chart": {
            "a3_solution": str(a3_solution),
            "cutoff": 6,
            "term_counts": {
                str(order): len(sp.Poly(equation, *boundary_variables).terms())
                for order, equation in boundary_equations.items()
            },
            "modular_results": boundary_modular,
            "exact_Q_result": boundary_exact,
        },
    }


def collision_data() -> list[dict[str, object]]:
    variables = (A1, A2, A3)
    saturation = A1 * A2 * A3
    results: list[dict[str, object]] = []
    for groups in set_partitions(4):
        if len(groups) in (1, 4):
            continue
        equations = compile_moments(groups, 6, (), variables)
        cutoff = 0
        modular_results: list[dict[str, object]] = []
        for candidate in sorted(equations):
            first = modular_cutoff(
                {
                    order: equation
                    for order, equation in equations.items()
                    if order <= candidate
                },
                saturation,
                variables,
                PRIMES[0],
                "Singular",
                120,
            )
            if first.get("unit") == 1:
                cutoff = candidate
                modular_results.append(first)
                break
        if cutoff == 0:
            raise AssertionError(f"collision partition survived through six: {groups}")
        selected = {
            order: equation
            for order, equation in equations.items()
            if order <= cutoff
        }
        for prime in PRIMES[1:]:
            replay = modular_cutoff(
                selected, saturation, variables, prime, "Singular", 120
            )
            assert replay.get("unit") == 1
            modular_results.append(replay)
        exact = exact_unit(equations, saturation, variables, cutoff, 120, 2)
        assert exact["unit"] == 1
        results.append(
            {
                "groups": [list(block) for block in groups],
                "directions": len(groups),
                "cutoff": cutoff,
                "nonzero_moment_orders": sorted(
                    order for order in equations if order <= cutoff
                ),
                "modular_results": modular_results,
                "exact_Q_result": exact,
            }
        )
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--exact-open",
        action="store_true",
        help="attempt the expensive characteristic-zero B!=0 saturation",
    )
    arguments = parser.parse_args()

    distinct = distinct_direction_data(arguments.exact_open)
    collisions = collision_data()
    assert len(collisions) == 13
    assert max(result["cutoff"] for result in collisions) == 6

    exact_open = distinct["B_nonzero_chart"]["exact_Q_saturation"]
    promoted = bool(exact_open and exact_open.get("unit") == 1)
    artifact = {
        "format": "gvc3-four-coherent-channels-v1",
        "status": (
            "exact characteristic-zero obstruction for the declared four-channel degree-eight family"
            if promoted
            else "bounded modular discovery on the four-distinct B!=0 chart; not a characteristic-zero theorem"
        ),
        "balanced_degree": 8,
        "laplacian_power": 4,
        "harmonic_degrees": list(DEGREES),
        "family": "P=sum_i a_i*rho^((8-ell_i)/2)*<v_i,z>^ell_i",
        "profile_condition": "one isotropic coherent state in each of H2,H4,H6,H8",
        "moment_sequence_promotion_gate": {
            "type": "finite proper-hypergeometric sum",
            "wick_coefficient": (
                "W(s)=[u^s] exp(sum_(i<j) g_ij*u_i*u_j)="
                "sum_(deg(e)=s) prod_(i<j) g_ij^e_ij/e_ij!"
            ),
            "moment": (
                "mu_m=sum_(|n|=m) multinomial(m;n)*prod_i(a_i^n_i)*"
                "prod_i((ell_i*n_i)!)/(sum_i ell_i*n_i+1)!!*W(ell*n)"
            ),
            "reason": "every summand has rational shift quotients in the occupation and edge indices",
        },
        "modular_discovery_primes": list(PRIMES),
        "four_distinct_directions": distinct,
        "direction_collision_partitions": collisions,
        "coefficient_boundaries": (
            "a zero coefficient reduces to the exact two/three-channel theorem GVC3IHC"
        ),
        "one_direction_terminal_reason": (
            "phase weight is one-sided, so every fixed multiplier is eventually killed"
        ),
        "conclusion": (
            "moments through seven force the one-direction terminal stratum; "
            "there is no Delta^4 witness in this four-channel degree-eight family"
            if promoted
            else "all collision and B=0 boundaries are excluded exactly; the four-distinct B!=0 chart is a three-prime modular unit candidate awaiting exact promotion"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    if promoted:
        print("PASS pairwise-distinct B!=0 chart: exact Q unit through moment 7")
    else:
        print("PASS pairwise-distinct B!=0 chart: modular units at 101,103,107")
        print("OPEN exact-Q promotion of the pairwise-distinct B!=0 chart")
    print("PASS pairwise-distinct B=0 chart: exact Q unit through moment 6")
    print("PASS 13 collision partitions: exact Q units through moment 6")
    print("PASS proper-hypergeometric moment formula recorded")


if __name__ == "__main__":
    main()
