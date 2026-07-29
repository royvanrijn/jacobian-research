#!/usr/bin/env python3
"""Finite moment-nullcone cutoff for quartic symbols with at most two roots.

After a linear change every such binary quartic is u^r*v^(4-r).  For
r=0,4 the first moment cuts out the one-sided hyperplane.  For r=1,2,3
this checker proves that the first four moment equations have radical equal
to the expected union of the two one-sided linear loci.
"""

from __future__ import annotations

from math import factorial
import json
from pathlib import Path
import shutil
import subprocess

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_variable_quartic_two_root_finite.json"
)

x, y, u, v = sp.symbols("x y u v")
a, b, c, d, e = sp.symbols("a b c d e")
P = a * x**4 + b * x**3 * y + c * x**2 * y**2 + d * x * y**3 + e * y**4


def moments(root_multiplicity: int) -> list[sp.Expr]:
    symbol = u**root_multiplicity * v ** (4 - root_multiplicity)
    values: list[sp.Expr] = []
    for order in range(1, 5):
        symbol_power = sp.Poly(sp.expand(symbol**order), u, v)
        polynomial_power = sp.Poly(sp.expand(P**order), x, y)
        value = sp.expand(
            sum(
                coefficient
                * polynomial_power.coeff_monomial(
                    x**x_order * y**y_order
                )
                * factorial(x_order)
                * factorial(y_order)
                for (x_order, y_order), coefficient
                in symbol_power.terms()
            )
        )
        values.append(
            sp.Poly(value, a, b, c, d, e).primitive()[1].as_expr()
        )
    return values


def singular_expression(expression: sp.Expr) -> str:
    return sp.sstr(expression).replace("**", "^")


def main() -> None:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"

    # The entries are the two prime one-sided ideals.  At the endpoints
    # there is only one component.
    expected_components = {
        0: (("e",),),
        1: (("a", "b", "c", "d"), ("d", "e")),
        2: (("a", "b", "c"), ("c", "d", "e")),
        3: (("a", "b"), ("b", "c", "d", "e")),
        4: (("a",),),
    }
    cutoffs: dict[str, int] = {}
    fiber_data: dict[str, dict[str, int]] = {}

    for root_multiplicity in range(5):
        values = moments(root_multiplicity)
        components = expected_components[root_multiplicity]
        if root_multiplicity in (0, 4):
            cutoff = 1
        else:
            cutoff = 4
        moment_generators = ",".join(
            singular_expression(value) for value in values[:cutoff]
        )
        program = f"""
ring r=0,(a,b,c,d,e),dp;
option(redSB);
ideal I={moment_generators};
ideal GI=std(I);
"""
        for index, component in enumerate(components, start=1):
            program += f"ideal J{index}={','.join(component)};\n"
        if len(components) == 1:
            program += "ideal J=J1;\n"
        else:
            program += "ideal J=intersect(J1,J2);\n"
        program += """
ideal GJ=std(J);
if (dim(GI)!=dim(GJ) || mult(GI)!=mult(GJ))
{
  print("BAD_DIMENSION_OR_DEGREE");
  exit(1);
}
if (size(reduce(I,GJ))!=0)
{
  print("BAD_CONTAINMENT");
  exit(1);
}
int generator;
for (generator=1;generator<=size(GJ);generator++)
{
  if (reduce(GJ[generator]^8,GI)!=0)
  {
    print("BAD_POWER_CERTIFICATE");
    exit(1);
  }
}
print("TWO_ROOT_FIBER");
print(dim(GI));
print(mult(GI));
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=True,
            timeout=60,
        )
        lines = completed.stdout.splitlines()
        assert lines[0] == "TWO_ROOT_FIBER", completed.stdout
        cutoffs[str(root_multiplicity)] = cutoff
        fiber_data[str(root_multiplicity)] = {
            "affine_cone_dimension": int(lines[1]),
            "projective_degree": int(lines[2]),
            "component_count": len(components),
            "radical_power_certificate": 8,
        }

    artifact = {
        "format": "two-variable-quartic-two-root-finite-v1",
        "field": "characteristic zero",
        "normal_forms": "u^r*v^(4-r), 0<=r<=4",
        "moment_cutoffs": cutoffs,
        "fibers": fiber_data,
        "consequence": (
            "the first four moments put every quartic rank-one tensor "
            "whose symbol has at most two roots in the one-sided nullcone"
        ),
        "written_source": (
            "extended-geometry/"
            "TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(
        "PASS quartic symbols with at most two roots: first four moments "
        "cut out the one-sided nullcone"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
