#!/usr/bin/env python3
"""Stage D: exact searches on thin Newton-directed Laurent supports."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import sympy as sp

from jcsearch.laurent import laurent_dict
from jcsearch.msolve import input_text, run
from jcsearch.newton import (laurent_support, leading_bracket_exponents,
                             stage_d_families)
from jcsearch.translated import (TranslatedTwoDivisorChart,
                                 cancellation_system, kernel_expressions,
                                 u, v, x, y)


PRIMES = (1000003, 1000033, 1000037)


def coordinates_in_basis(target, basis):
    """Exact coordinates of a Laurent expression, or None if absent."""
    coefficients = sp.symbols(f"k0:{len(basis)}")
    residual = sp.expand(sum(c*f for c, f in zip(coefficients, basis)) - target)
    equations = list(laurent_dict(residual, (u, v)).values())
    solution = sp.linsolve(equations, coefficients)
    if solution is sp.EmptySet:
        return None
    point = next(iter(solution))
    free = set().union(*(value.free_symbols for value in point)) - set(coefficients)
    substitutions = {symbol: 0 for symbol in free}
    return tuple(sp.simplify(value.subs(substitutions)) for value in point)


def normalized_equations(P, Q, A, B):
    residual = sp.expand(sp.Matrix([A, B]).jacobian((u, v)).det() + u)
    equations = list(laurent_dict(residual, (u, v)).values())
    equations += [
        P.subs({x: 0, y: 0}), Q.subs({x: 0, y: 0}),
        P.subs({x: 1, y: 0}), Q.subs({x: 1, y: 0}),
        sp.diff(P, x).subs({x: 0, y: 0}) - 1,
        sp.diff(P, y).subs({x: 0, y: 0}),
        sp.diff(Q, x).subs({x: 0, y: 0}),
        sp.diff(Q, y).subs({x: 0, y: 0}) - 1,
    ]
    return [sp.expand(e) for e in equations if sp.expand(e) != 0]


def solve_family(name, family, timeout):
    chart = TranslatedTwoDivisorChart()
    Psupport = laurent_support(family["P_polynomial_support"])
    Qsupport = laurent_support(family["Q_polynomial_support"])
    Pcancel = cancellation_system(Psupport, "pc")
    Qcancel = cancellation_system(Qsupport, "qc")
    Pbasis, Qbasis = kernel_expressions(Pcancel), kernel_expressions(Qcancel)

    aa = sp.symbols(f"a0:{len(Pbasis)}")
    bb = sp.symbols(f"b0:{len(Qbasis)}")
    A = sp.expand(sum(c*f for c, f in zip(aa, Pbasis)))
    B = sp.expand(sum(c*f for c, f in zip(bb, Qbasis)))
    P, Q = chart.pullback(A), chart.pullback(B)
    assert sp.denom(P) == 1 and sp.denom(Q) == 1

    equations = normalized_equations(P, Q, A, B)
    variables = aa + bb
    reachable = (1, 0) in leading_bracket_exponents(Psupport, Qsupport)
    identity_P = coordinates_in_basis(u*v, Pbasis)
    identity_Q = coordinates_in_basis(u - 1/(u*v), Qbasis)
    identity_control = identity_P is not None and identity_Q is not None
    record = {
        "family": name,
        "description": family["description"],
        "status": "prototype; not a faithful encoding of the (72,108) chain",
        "P_polynomial_support": family["P_polynomial_support"],
        "Q_polynomial_support": family["Q_polynomial_support"],
        "P_laurent_support_size": len(Psupport),
        "Q_laurent_support_size": len(Qsupport),
        "P_kernel_dimension": len(Pbasis),
        "Q_kernel_dimension": len(Qbasis),
        "P_pole_rank": Pcancel["certificate"].rank,
        "Q_pole_rank": Qcancel["certificate"].rank,
        "target_bracket_exponent_reachable": reachable,
        "identity_keller_control_in_kernels": identity_control,
        "equations": len(equations),
        "variables": len(variables),
        "runs": [],
    }

    for characteristic in PRIMES + (0,):
        try:
            result = run(equations, variables, prime=characteristic,
                         threads=4, timeout=timeout)
            item = {
                "prime": characteristic,
                "empty": result.empty,
                "positive_dimensional": result.positive_dimensional,
                "returncode": result.returncode,
            }
            if characteristic == 0:
                base = ROOT / "results" / "certificates" / f"stage_d_{name}_Q"
                base.with_suffix(".input").write_text(
                    input_text(equations, variables, 0))
                base.with_suffix(".msolve").write_text(result.output)
        except subprocess.TimeoutExpired:
            item = {"prime": characteristic, "timeout": timeout,
                    "empty": None}
        record["runs"].append(item)

        # A timeout on an easier modular run is already a useful scaling
        # certificate.  Do not spend the same bound on two more primes and Q.
        if item.get("timeout"):
            break

    exact = next((item for item in record["runs"] if item["prime"] == 0), None)
    if exact and exact.get("empty"):
        record["failure_class"] = (
            "collision-normalized nonlinear unit ideal; identity Keller control remains"
            if identity_control else
            "nonlinear unit ideal after exact pole cancellation"
        )
    elif not reachable:
        record["failure_class"] = "target Jacobian Laurent exponent absent"
    elif any(item.get("timeout") for item in record["runs"]):
        record["failure_class"] = "unclassified: F4 scaling timeout"
        record["exact_run"] = "skipped after earlier modular timeout"
    else:
        record["failure_class"] = "survivor requiring saturation/filtering"
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", action="append",
                        help="run only this named family (repeatable)")
    parser.add_argument("--timeout", type=int, default=180,
                        help="seconds per msolve run")
    args = parser.parse_args()
    families = stage_d_families()
    selected = args.family or list(families)
    unknown = set(selected) - set(families)
    if unknown:
        raise SystemExit(f"unknown families: {sorted(unknown)}")

    records = []
    for name in selected:
        print(f"building {name}", flush=True)
        record = solve_family(name, families[name], args.timeout)
        records.append(record)
        print(json.dumps(record, default=sp.sstr, indent=2), flush=True)
    destination = ROOT / "results" / "stage_d_newton_bands.json"
    destination.write_text(json.dumps(records, default=sp.sstr, indent=2))
    print(f"wrote {destination}")


if __name__ == "__main__":
    main()
