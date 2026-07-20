#!/usr/bin/env python3
"""Build the exact broad polynomial systems from Proposition 4.3 `(8,28)`.

The reduced pairs satisfy [P,Q]=x^2.  This script intentionally does not add
the affine collision normalization: the rational polygon transformations in
the proposition do not preserve those chosen affine points.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import sympy as sp

from jcsearch.msolve import input_text, run
from jcsearch.newton import (convex_lattice_points,
                             leading_bracket_exponents,
                             reduced_72_108_polygons)

x, y = sp.symbols("x y")


def edge_power_substitution(name, specification, Psupport, Qsupport, aa, bb):
    """Common square/cube boundary pattern suggested by the reduction proof."""
    pindex = {point: i for i, point in enumerate(Psupport)}
    qindex = {point: i for i, point in enumerate(Qsupport)}
    pscale, qscale, alpha, beta = sp.symbols("pscale qscale alpha beta")
    substitutions = {}

    if name == "case_with_vertical_vertices":
        pedge = [((k, 8 + k), sp.binomial(8, k)*(-alpha)**(8-k))
                 for k in range(9)]
        qedge = [((k, 12 + k), sp.binomial(12, k)*(-alpha)**(12-k))
                 for k in range(13)]
    else:
        pedge = [((k, 2*k), sp.binomial(8, k)*(-alpha)**(8-k))
                 for k in range(9)]
        qedge = [((k, 2*k), sp.binomial(12, k)*(-alpha)**(12-k))
                 for k in range(13)]
    pedge += [((8, 14+k), sp.binomial(2, k)*(-beta)**(2-k))
              for k in range(3)]
    qedge += [((12, 21+k), sp.binomial(3, k)*(-beta)**(3-k))
              for k in range(4)]
    for point, coefficient in pedge:
        substitutions[aa[pindex[point]]] = pscale*coefficient
    for point, coefficient in qedge:
        substitutions[bb[qindex[point]]] = qscale*coefficient
    remaining = tuple(c for c in aa + bb if c not in substitutions)
    return substitutions, remaining + (pscale, qscale, alpha, beta)


def build_case(name, specification, edge_powers=False):
    Psupport = convex_lattice_points(specification["P_vertices"])
    Qsupport = convex_lattice_points(specification["Q_vertices"])
    aa = sp.symbols(f"p0:{len(Psupport)}")
    bb = sp.symbols(f"q0:{len(Qsupport)}")
    P = sp.expand(sum(c*x**a*y**b for c, (a, b) in zip(aa, Psupport)))
    Q = sp.expand(sum(c*x**a*y**b for c, (a, b) in zip(bb, Qsupport)))
    substitutions = {}
    variables = aa + bb
    if edge_powers:
        substitutions, variables = edge_power_substitution(
            name, specification, Psupport, Qsupport, aa, bb)
        P = sp.expand(P.subs(substitutions, simultaneous=True))
        Q = sp.expand(Q.subs(substitutions, simultaneous=True))
    residual = sp.Poly(sp.expand(sp.diff(P, x)*sp.diff(Q, y)
                                 - sp.diff(P, y)*sp.diff(Q, x) - x**2), x, y)
    equations = [sp.expand(coefficient) for _, coefficient in residual.terms()]
    pindex = {point: i for i, point in enumerate(Psupport)}
    qindex = {point: i for i, point in enumerate(Qsupport)}
    vertex_coefficients = ([aa[pindex[p]] for p in specification["P_vertices"]]
                           + [bb[qindex[p]] for p in specification["Q_vertices"]])
    saturator = sp.expand(sp.prod(vertex_coefficients).subs(
        substitutions, simultaneous=True))
    return {
        "name": name,
        "P_support": Psupport,
        "Q_support": Qsupport,
        "P": P,
        "Q": Q,
        "equations": equations,
        "variables": variables,
        "saturator": saturator,
        "metrics": {
            "P_lattice_points": len(Psupport),
            "Q_lattice_points": len(Qsupport),
            "variables": len(variables),
            "jacobian_coefficient_equations": len(equations),
            "required_nonzero_vertex_coefficients": len(vertex_coefficients),
            "target_x2_exponent_reachable":
                (2, 0) in leading_bracket_exponents(Psupport, Qsupport),
            "edge_power_refinement": edge_powers,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", default=[],
                        help="case to run with modular F4SAT (repeatable)")
    parser.add_argument("--prime", type=int, default=1000003)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--edge-powers", action="store_true",
                        help="substitute common square/cube forms on two boundary edges")
    args = parser.parse_args()
    specifications = reduced_72_108_polygons()
    unknown = set(args.run) - set(specifications)
    if unknown:
        raise SystemExit(f"unknown cases: {sorted(unknown)}")

    output = []
    systems_dir = ROOT / "results" / "systems"
    systems_dir.mkdir(parents=True, exist_ok=True)
    for name, specification in specifications.items():
        system = build_case(name, specification, edge_powers=args.edge_powers)
        record = {
            "case": name,
            "source_scope": "Proposition 4.3 reduced Newton polygons",
            "refinement_scope": ("experimental common square/cube boundary pattern; cited leading-form implications still require independent audit"
                                 if args.edge_powers else "all lattice coefficients in the reduced polygon"),
            "bracket": "x^2",
            "collision_normalization": "not imposed; not preserved by reduction",
            "P_vertices": specification["P_vertices"],
            "Q_vertices": specification["Q_vertices"],
            **system["metrics"],
            "run": None,
        }
        suffix = "edgepowers" if args.edge_powers else "broad"
        source = systems_dir / f"reduced_72_108_{name}_{suffix}_modp.input"
        source.write_text(input_text(system["equations"] + [system["saturator"]],
                                     system["variables"], args.prime))
        record["msolve_input"] = str(source.relative_to(ROOT))
        record["msolve_mode"] = "last polynomial is saturator; invoke with -S"
        if name in args.run:
            try:
                result = run(system["equations"], system["variables"],
                             prime=args.prime, saturate=system["saturator"],
                             timeout=args.timeout)
                record["run"] = {"prime": args.prime,
                                 "empty_after_vertex_saturation": result.empty,
                                 "positive_dimensional": result.positive_dimensional,
                                 "returncode": result.returncode}
            except subprocess.TimeoutExpired:
                record["run"] = {"prime": args.prime, "timeout": args.timeout,
                                 "classification": "complexity datum, not emptiness"}
        output.append(record)
        print(json.dumps(record, indent=2))

    suffix = "edgepowers" if args.edge_powers else "broad"
    destination = ROOT / "results" / f"reduced_72_108_systems_{suffix}.json"
    destination.write_text(json.dumps(output, indent=2))
    print(f"wrote {destination}")


if __name__ == "__main__":
    main()
