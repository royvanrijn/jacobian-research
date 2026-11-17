#!/usr/bin/env python3
"""Exact-slice scout for the degree-ten five-channel coherent profile.

The full pairwise-distinct chart has cross-ratios ``lam`` and ``mu``.  This
script fixes one rational configuration before compiling the moments, so the
Wick Gram matrix is numeric and higher moments remain small.  It studies the
exceptional boundary of the best first-moment pivot: if

    mu_3 = A + B*a4,

then nonzero ``a4`` splits the chart into the main localization ``A*B != 0``
and the boundary ``A=B=0``.  Modular quotient saturation is discovery only;
a chart is promoted only when msolve returns the literal basis ``[1]`` over
the rationals.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

from research_gvc3_many_coherent_channels import (
    exact_unit,
    modular_saturation_cutoff,
    moment,
    primitive_polynomial,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "gvc3_degree10_five_channel_a4_boundary_slice.json"
)
DEGREES = (2, 4, 6, 8, 10)
GROUPS = tuple((index,) for index in range(5))
PRIMES = (101, 103, 107)


def sha256(expression: sp.Expr) -> str:
    return hashlib.sha256(str(expression).encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lam", type=sp.Rational, default=sp.Rational(2))
    parser.add_argument("--mu", type=sp.Rational, default=sp.Rational(3))
    parser.add_argument("--max-order", type=int, default=8)
    parser.add_argument("--primes", nargs="+", type=int, default=list(PRIMES))
    parser.add_argument("--singular", default="Singular")
    parser.add_argument("--modular-timeout", type=int, default=180)
    parser.add_argument("--exact-timeout", type=int, default=300)
    parser.add_argument("--msolve-threads", type=int, default=4)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    if arguments.lam in (0, 1) or arguments.mu in (0, 1, arguments.lam):
        raise SystemExit("lam and mu must define five pairwise-distinct directions")
    if arguments.max_order < 3:
        raise SystemExit("--max-order must be at least three")

    coefficients = sp.symbols("a0:5")
    a0, a1, a2, a3, a4 = coefficients
    variables = (a1, a2, a3, a4)
    parameters = (arguments.lam, arguments.mu)
    moments: dict[int, sp.Expr] = {}
    term_counts: dict[str, int] = {}
    moment_hashes: dict[str, str] = {}
    for order in range(2, arguments.max_order + 1):
        polynomial = primitive_polynomial(
            moment(DEGREES, GROUPS, order, coefficients, parameters).subs(a0, 1),
            variables,
        )
        term_counts[str(order)] = (
            0 if polynomial == 0 else len(sp.Poly(polynomial, *variables).terms())
        )
        moment_hashes[str(order)] = sha256(polynomial)
        if polynomial != 0:
            moments[order] = polynomial
        print(
            f"COMPILED moment {order} terms={term_counts[str(order)]}",
            flush=True,
        )

    first = sp.Poly(moments[3], a4)
    if first.degree() != 1:
        raise RuntimeError("the declared a4 pivot is not linear in moment three")
    coefficient = primitive_polynomial(first.coeff_monomial(a4), variables)
    remainder = primitive_polynomial(first.coeff_monomial(1), variables)
    boundary_equations = {
        2: remainder,
        3: coefficient,
        **{order: equation for order, equation in moments.items() if order > 3},
    }
    saturation = sp.prod(variables)

    attempts_by_prime: dict[str, list[dict[str, object]]] = {}
    unit_cutoffs: dict[str, int] = {}
    discovery_cutoff = 3
    for prime_index, prime in enumerate(arguments.primes):
        attempts: list[dict[str, object]] = []
        start = discovery_cutoff if prime_index else 3
        for cutoff in range(start, arguments.max_order + 1):
            selected = {
                order: equation
                for order, equation in boundary_equations.items()
                if order <= cutoff
            }
            result = modular_saturation_cutoff(
                selected,
                saturation,
                variables,
                prime,
                arguments.singular,
                arguments.modular_timeout,
                report_quotient_dimension=True,
            )
            result["attempted_cutoff"] = cutoff
            attempts.append(result)
            print(
                f"MODULAR p={prime} cutoff={cutoff} "
                f"status={result.get('status')} unit={result.get('unit')} "
                f"dim={result.get('dimension')} "
                f"qdim={result.get('quotient_dimension')}",
                flush=True,
            )
            if result.get("status") != "completed":
                break
            if result.get("unit") == 1:
                unit_cutoffs[str(prime)] = cutoff
                if prime_index == 0:
                    discovery_cutoff = cutoff
                break
        attempts_by_prime[str(prime)] = attempts

    exact_result = None
    if len(unit_cutoffs) == len(arguments.primes):
        exact_cutoff = max(unit_cutoffs.values())
        exact_result = exact_unit(
            boundary_equations,
            saturation,
            variables,
            exact_cutoff,
            arguments.exact_timeout,
            arguments.msolve_threads,
        )
        print(
            f"EXACT cutoff={exact_cutoff} unit={exact_result.get('unit')}",
            flush=True,
        )
        if exact_result.get("unit") != 1:
            raise RuntimeError("modular unit did not promote to an exact Q basis [1]")

    if exact_result and exact_result.get("unit") == 1:
        status = "exact characteristic-zero exclusion of one rational slice boundary"
        conclusion = (
            "the a4-pivot boundary is empty on the declared rational "
            "cross-ratio slice"
        )
    else:
        status = "bounded rational-slice discovery; not a generic-chart theorem"
        conclusion = (
            "the declared rational slice was not excluded through the compiled cutoff"
        )
    artifact = {
        "format": "gvc3-degree10-five-channel-a4-boundary-slice-v1",
        "status": status,
        "balanced_degree": 10,
        "laplacian_power": 5,
        "harmonic_degrees": list(DEGREES),
        "directions": "infinity, zero, one, lam, mu",
        "cross_ratio_slice": {"lam": str(arguments.lam), "mu": str(arguments.mu)},
        "coefficient_normalization": "a0=1",
        "coefficient_saturation": "a1*a2*a3*a4",
        "pivot": "mu3=A+B*a4",
        "pivot_boundary": "A=B=0",
        "A_terms": len(sp.Poly(remainder, *variables).terms()),
        "B_terms": len(sp.Poly(coefficient, *variables).terms()),
        "A_sha256": sha256(remainder),
        "B_sha256": sha256(coefficient),
        "max_order": arguments.max_order,
        "moment_term_counts": term_counts,
        "moment_sha256": moment_hashes,
        "modular_attempts": attempts_by_prime,
        "modular_unit_cutoffs": unit_cutoffs,
        "exact_Q_result": exact_result,
        "conclusion": conclusion,
        "not_in_scope": [
            "the generic two-cross-ratio boundary",
            "the main A*B!=0 pivot localization",
            "other rational or algebraic cross-ratio fibers",
            "noncoherent or repeated harmonic profiles",
        ],
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
