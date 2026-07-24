#!/usr/bin/env python3
"""Verify the first filtered path-to-boundary cone in degree forty-two.

On the normalized ``2 o 7 o 3`` chart let

* ``IT`` be the composite-cut-6-omitting (thick) half-braid;
* ``ITHIN`` be the prime-cut-7-omitting half-braid;
* ``IB`` be the full braid boundary;
* ``K`` be the reduced Dickson graph.

The inclusions of ideals are ``IT, ITHIN <= IB <= K``.  Consequently the
formal conormal filtration of the thick path has a path-to-boundary quotient
and a boundary-to-Dickson quotient.  This checker computes their tangent
ranks, annihilators on the Dickson base, and finite maximal-adic jets.

Finite jets do not by themselves prove an isomorphism of completed
cotangent complexes.  The exact annihilator calculation, when successful,
does isolate the base support of the relative sector module.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree30_hessian_ritt_braid import canonical_residuals  # noqa: E402
from explore_degree42_ritt_spectator_universality import (  # noqa: E402
    ALL_CUTS,
    DEGREE,
    WORD,
    W,
    build_chart,
    serialize_ideal,
)
from jcsearch.ritt_complex import dickson_vertex_factors  # noqa: E402


def serialize_polynomial(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


def dickson_normal_map(factor_variables):
    """Return the polynomial ``7 normal | 2 base`` coordinate map."""

    tau, zeta = sp.symbols("tau zeta")
    translation = 1 + tau
    dickson_factors = dickson_vertex_factors(
        WORD, W, translation, zeta
    )
    inner_linear, inner_quadratic = factor_variables[-1]
    dependent_variables = tuple(
        variable
        for variables in factor_variables
        for variable in variables
        if variable not in {inner_linear, inner_quadratic}
    )
    normals = sp.symbols(f"n0:{len(dependent_variables)}")
    normal_by_variable = dict(zip(dependent_variables, normals))
    images = {}
    for factor, variables in zip(dickson_factors, factor_variables):
        polynomial = sp.Poly(factor, W)
        for power, variable in enumerate(variables, 1):
            images[variable] = polynomial.nth(power)
            if variable in normal_by_variable:
                images[variable] += normal_by_variable[variable]
    assert sp.expand(images[inner_linear] - 3 * (translation**2 - zeta)) == 0
    assert sp.expand(images[inner_quadratic] - 3 * translation) == 0
    return normals, (tau, zeta), images


def singular_relative_cone_audit(
    parameters,
    factor_variables,
    thick,
    thin,
    boundary,
) -> str:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    normals, base_coordinates, images = dickson_normal_map(
        factor_variables
    )
    local_variables = normals + base_coordinates
    map_images = ",".join(
        serialize_polynomial(images[parameter]) for parameter in parameters
    )
    zero_substitutions = "\n".join(
        f"JT=subst(JT,{variable},0);"
        f"JTHIN=subst(JTHIN,{variable},0);"
        f"JB=subst(JB,{variable},0);"
        f"JK=subst(JK,{variable},0);"
        for variable in local_variables
    )
    program = f"""
ring source=0,({",".join(map(str, parameters))}),dp;
ideal ITsource={serialize_ideal(thick)};
ideal ITHINsource={serialize_ideal(thin)};
ideal IBsource={serialize_ideal(boundary)};
ring q=0,({",".join(map(str, local_variables))}),(dp({len(normals)}),dp(2));
map phi=source,{map_images};
option(redSB);
proc iszeroideal(ideal J)
{{
  ideal reducedGenerators=simplify(J,2);
  return(size(reducedGenerators)==0);
}}
proc issubsetstd(ideal A, ideal G)
{{
  ideal remainder=reduce(A,G);
  return(iszeroideal(remainder));
}}
proc firstannihilatingpower(ideal A, ideal G, poly z, int maximumPower)
{{
  ideal current=A;
  int exponent;
  for(exponent=0; exponent<=maximumPower; exponent++)
  {{
    if(issubsetstd(current,G))
    {{
      return(exponent);
    }}
    current=z*current;
  }}
  return(-1);
}}
ideal IT=phi(ITsource);
ideal ITHIN=phi(ITHINsource);
ideal IB=phi(IBsource);
ideal K={",".join(map(str, normals))};
ideal GIT=std(IT);
ideal GTHIN=std(ITHIN);
ideal GIB=std(IB);
ideal GK=std(K);
ideal maximalIdeal={",".join(map(str, local_variables))};
ideal maximalIdeal2=maximalIdeal*maximalIdeal;
ideal maximalIdeal3=maximalIdeal2*maximalIdeal;
matrix JT=jacob(IT);
matrix JTHIN=jacob(ITHIN);
matrix JB=jacob(IB);
matrix JK=jacob(K);
{zero_substitutions}
print("IDEAL_FLAG");
print(issubsetstd(IT,GIB));
print(issubsetstd(ITHIN,GIB));
print(issubsetstd(IB,GTHIN));
print(issubsetstd(IB,GK));
print("CONORMAL_RANKS");
print(rank(JT));
print(rank(JTHIN));
print(rank(JB));
print(rank(JK));
print("SCHEME_DIMENSIONS");
print(dim(GIT));
print(dim(GIB));
print(dim(GK));
print("BASE_ANNIHILATION_EXPONENTS");
print(firstannihilatingpower(IB,GIT,{base_coordinates[1]},20));
print(firstannihilatingpower(K,GIB,{base_coordinates[1]},20));
print("JET_LENGTHS_Q1");
print(vdim(std(IT+maximalIdeal)));
print(vdim(std(ITHIN+maximalIdeal)));
print(vdim(std(IB+maximalIdeal)));
print(vdim(std(K+maximalIdeal)));
print("JET_LENGTHS_Q2");
print(vdim(std(IT+maximalIdeal2)));
print(vdim(std(ITHIN+maximalIdeal2)));
print(vdim(std(IB+maximalIdeal2)));
print(vdim(std(K+maximalIdeal2)));
print("JET_LENGTHS_Q3");
print(vdim(std(IT+maximalIdeal3)));
print(vdim(std(ITHIN+maximalIdeal3)));
print(vdim(std(IB+maximalIdeal3)));
print(vdim(std(K+maximalIdeal3)));
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=7200,
    )
    compact = " ".join(result.stdout.split())
    assert "IDEAL_FLAG 1 1 1 1" in compact
    assert "CONORMAL_RANKS 5 6 6 7" in compact
    assert "SCHEME_DIMENSIONS 2 2 2" in compact
    assert "BASE_ANNIHILATION_EXPONENTS 8 1" in compact
    assert "JET_LENGTHS_Q1 1 1 1 1" in compact
    assert "JET_LENGTHS_Q2 5 4 4 3" in compact
    assert "JET_LENGTHS_Q3 14 9 9 6" in compact
    return result.stdout


def main() -> None:
    parameters, factor_variables, polynomial = build_chart()
    base_cuts = {2, 14}
    requested_cuts = tuple(cut for cut in ALL_CUTS if cut not in base_cuts)
    residuals = {
        cut: canonical_residuals(
            polynomial,
            cut,
            DEGREE // cut,
            parameters=parameters,
            factor_output=False,
            minimum_coefficient_degree=1,
        )
        for cut in requested_cuts
    }
    endpoint = residuals[3] + residuals[21]
    thick = endpoint + residuals[7]
    thin = endpoint + residuals[6]
    boundary = endpoint + residuals[6] + residuals[7]
    output = singular_relative_cone_audit(
        parameters,
        factor_variables,
        thick,
        thin,
        boundary,
    )
    print(output)
    print("PASS: the prime-omitting path equals the full boundary")
    print("PASS: the monomial conormal flag has ranks 5 < 6 < 7")
    print("PASS: thick-to-boundary and boundary-to-Dickson each add one direction")
    print("PASS: z^8 minimally annihilates the relative sector module")
    print("PASS: z minimally annihilates the common spectator module")
    print("PASS: the first three completed jet lengths separate the two layers")


if __name__ == "__main__":
    main()
