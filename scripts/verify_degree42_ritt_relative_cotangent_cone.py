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


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree30_hessian_ritt_braid import canonical_residuals  # noqa: E402
from explore_degree42_ritt_spectator_universality import (  # noqa: E402
    ALL_CUTS,
    DEGREE,
    build_chart,
    coefficient_point,
    dickson_graph,
    serialize_ideal,
)


def singular_relative_cone_audit(
    parameters,
    thick,
    thin,
    boundary,
    graph,
    parameter_numerator,
    inner_quadratic,
    monomial_point,
) -> str:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    point_substitutions = "\n".join(
        f"JT=subst(JT,{variable},{value});"
        f"JTHIN=subst(JTHIN,{variable},{value});"
        f"JB=subst(JB,{variable},{value});"
        f"JK=subst(JK,{variable},{value});"
        for variable, value in monomial_point.items()
    )
    parameter_text = str(parameter_numerator).replace("**", "^")
    program = f"""
ring q=0,({",".join(map(str, parameters))}),dp;
option(redSB);
proc iszeroideal(ideal J)
{{
  ideal reducedGenerators=simplify(J,2);
  return(size(reducedGenerators)==0);
}}
proc issubset(ideal A, ideal B)
{{
  ideal remainder=reduce(A,std(B));
  return(iszeroideal(remainder));
}}
ideal IT={serialize_ideal(thick)};
ideal ITHIN={serialize_ideal(thin)};
ideal IB={serialize_ideal(boundary)};
ideal K={serialize_ideal(graph)};
ideal zIdeal=K,({parameter_text});
ideal z2Ideal=K,({parameter_text})^2;
ideal z3Ideal=K,({parameter_text})^3;
ideal z4Ideal=K,({parameter_text})^4;
ideal z5Ideal=K,({parameter_text})^5;
ideal sectorAnnihilator=quotient(IT,IB);
ideal thinBoundaryAnnihilator=quotient(ITHIN,IB);
ideal boundarySpectatorAnnihilator=quotient(IB,K);
ideal thinSpectatorAnnihilator=quotient(ITHIN,K);
ideal maximalIdeal=K,({parameter_text}),{inner_quadratic}-3;
ideal maximalIdeal2=maximalIdeal*maximalIdeal;
ideal maximalIdeal3=maximalIdeal2*maximalIdeal;
ideal maximalIdeal4=maximalIdeal3*maximalIdeal;
ideal maximalIdeal5=maximalIdeal4*maximalIdeal;
matrix JT=jacob(IT);
matrix JTHIN=jacob(ITHIN);
matrix JB=jacob(IB);
matrix JK=jacob(K);
{point_substitutions}
print("IDEAL_FLAG");
print(issubset(IT,IB));
print(issubset(ITHIN,IB));
print(issubset(IB,K));
print("CONORMAL_RANKS");
print(rank(JT));
print(rank(JTHIN));
print(rank(JB));
print(rank(JK));
print("SECTOR_ANNIHILATOR_Z_POWERS");
print(issubset(sectorAnnihilator+K,zIdeal));
print(issubset(zIdeal,sectorAnnihilator+K));
print(issubset(sectorAnnihilator+K,z2Ideal));
print(issubset(z2Ideal,sectorAnnihilator+K));
print(issubset(sectorAnnihilator+K,z3Ideal));
print(issubset(z3Ideal,sectorAnnihilator+K));
print(issubset(sectorAnnihilator+K,z4Ideal));
print(issubset(z4Ideal,sectorAnnihilator+K));
print(issubset(sectorAnnihilator+K,z5Ideal));
print(issubset(z5Ideal,sectorAnnihilator+K));
print("SPECTATOR_ANNIHILATORS_AGREE_ON_K");
print(issubset(boundarySpectatorAnnihilator+K,thinSpectatorAnnihilator+K));
print(issubset(thinSpectatorAnnihilator+K,boundarySpectatorAnnihilator+K));
print("THIN_BOUNDARY_ANNIHILATOR_IS_UNIT");
print(vdim(std(thinBoundaryAnnihilator)));
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
print("JET_LENGTHS_Q4");
print(vdim(std(IT+maximalIdeal4)));
print(vdim(std(ITHIN+maximalIdeal4)));
print(vdim(std(IB+maximalIdeal4)));
print(vdim(std(K+maximalIdeal4)));
print("JET_LENGTHS_Q5");
print(vdim(std(IT+maximalIdeal5)));
print(vdim(std(ITHIN+maximalIdeal5)));
print(vdim(std(IB+maximalIdeal5)));
print(vdim(std(K+maximalIdeal5)));
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
    assert "IDEAL_FLAG 1 1 1" in compact
    assert "CONORMAL_RANKS 5 6 6 7" in compact
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
    graph, parameter_numerator, inner_quadratic = dickson_graph(
        factor_variables
    )
    output = singular_relative_cone_audit(
        parameters,
        thick,
        thin,
        boundary,
        graph,
        parameter_numerator,
        inner_quadratic,
        coefficient_point(factor_variables, 1, 0),
    )
    print(output)
    print("PASS: the monomial conormal flag has ranks 5 < 6 < 7")
    print("PASS: thick-to-boundary and boundary-to-Dickson each add one direction")


if __name__ == "__main__":
    main()
