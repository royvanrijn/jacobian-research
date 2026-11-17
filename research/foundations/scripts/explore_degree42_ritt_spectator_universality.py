#!/usr/bin/env python3
"""Explore spectator independence for the degree-42 Dickson braid.

The endpoint word is ``2 o 7 o 3``.  Its two half-braid paths to
``3 o 7 o 2`` omit respectively the composite cut 6 and the prime cut 7.
The degree-30 cut-6 sector predicts that the composite-omitting path is a
nilpotent thickening with:

* nilpotence index 4;
* defect annihilator ``(z^4)`` on the Dickson normalization;
* conductor slice ``Q[u,v]/(u^4,v^2)``;
* K-adic lengths ``(4,7,8,8)``.

This is the first direct test of whether these data depend on the outer pair
``{2,3}`` but not on the middle prime (5 in degree 30, 7 here).  Unlike the
degree-30 calculation, the raw degree-42 path ideals have additional global
scheme structure.  The script records that obstruction and tangent probes;
it deliberately does not call the spectator-independence statement proved.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from math import comb
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree30_hessian_ritt_braid import canonical_residuals  # noqa: E402
from jcsearch.ritt_complex import (  # noqa: E402
    compose_factors,
    dickson_vertex_factors,
)


W, Z = sp.symbols("W Z")
WORD = (2, 7, 3)
DEGREE = 42
ALL_CUTS = (2, 3, 6, 7, 14, 21)


def serialize_ideal(equations: list[sp.Expr]) -> str:
    return ",".join(str(equation).replace("**", "^") for equation in equations)


def build_chart() -> tuple[
    tuple[sp.Symbol, ...],
    tuple[tuple[sp.Symbol, ...], ...],
    sp.Expr,
]:
    factor_variables: list[tuple[sp.Symbol, ...]] = []
    factors: list[sp.Expr] = []
    for position, degree in enumerate(WORD):
        variables = sp.symbols(f"x{position}_1:{degree}")
        factor_variables.append(variables)
        factors.append(
            W**degree
            + sum(
                variables[power - 1] * W**power
                for power in range(1, degree)
            )
        )
    parameters = tuple(
        variable for variables in factor_variables for variable in variables
    )
    return parameters, tuple(factor_variables), compose_factors(tuple(factors), W)


def dickson_graph(
    factor_variables: tuple[tuple[sp.Symbol, ...], ...],
) -> tuple[list[sp.Expr], sp.Expr, sp.Symbol]:
    inner_linear, inner_quadratic = factor_variables[-1]
    translation = inner_quadratic / 3
    parameter = translation**2 - inner_linear / 3
    factors = dickson_vertex_factors(WORD, W, translation, parameter)
    coefficient_map: dict[sp.Symbol, sp.Expr] = {}
    for factor, variables in zip(factors, factor_variables):
        polynomial = sp.Poly(factor, W)
        coefficient_map.update(
            {
                variable: polynomial.nth(power)
                for power, variable in enumerate(variables, 1)
            }
        )
    free = {inner_linear, inner_quadratic}
    graph = [
        sp.together(variable - value).as_numer_denom()[0]
        for variable, value in coefficient_map.items()
        if variable not in free
    ]
    return graph, sp.together(parameter).as_numer_denom()[0], inner_quadratic


def singular_audit(
    parameters: tuple[sp.Symbol, ...],
    thick: list[sp.Expr],
    thin: list[sp.Expr],
    boundary: list[sp.Expr],
    graph: list[sp.Expr],
    parameter_numerator: sp.Expr,
    inner_quadratic: sp.Symbol,
    dickson_points: tuple[dict[sp.Symbol, sp.Expr], ...],
) -> str:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    tangent_commands = []
    for point_index, point in enumerate(dickson_points):
        for ideal_name in ("IT", "ITHIN", "IB", "K"):
            matrix_name = f"J{point_index}{ideal_name}"
            tangent_commands.append(f"matrix {matrix_name}=jacob({ideal_name});")
            tangent_commands.extend(
                f"{matrix_name}=subst({matrix_name},{variable},{value});"
                for variable, value in point.items()
            )
            tangent_commands.append(
                f'print("{point_index}_{ideal_name}");'
                f"print({len(parameters)}-rank({matrix_name}));"
            )
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
ideal K2=K*K;
ideal K3=K2*K;
ideal K4=K3*K;
ideal conductor=quotient(IT,K);
ideal conductorScheme2=K,({str(parameter_numerator).replace("**", "^")})^2;
ideal conductorScheme3=K,({str(parameter_numerator).replace("**", "^")})^3;
ideal conductorScheme4=K,({str(parameter_numerator).replace("**", "^")})^4;
ideal fixedTranslation={inner_quadratic}-3;
ideal conductorSlice=IT+conductor+fixedTranslation;
ideal conductorSliceK=conductorSlice+K;
ideal conductorSliceK2=conductorSlice+K2;
ideal conductorSliceK3=conductorSlice+K3;
ideal conductorSliceK4=conductorSlice+K4;
ideal maximalIdeal=K,({str(parameter_numerator).replace("**", "^")}),fixedTranslation;
ideal maximalIdeal2=maximalIdeal*maximalIdeal;
ideal maximalIdeal3=maximalIdeal2*maximalIdeal;
ideal maximalIdeal4=maximalIdeal3*maximalIdeal;
ideal maximalIdeal5=maximalIdeal4*maximalIdeal;
ideal socleIdeal=quotient(conductorSlice,maximalIdeal);
matrix sliceJacobian=jacob(conductorSlice);
print("THICK");
print(issubset(IT,K));
print(issubset(K,IT));
print("THIN");
print(issubset(ITHIN,K));
print(issubset(K,ITHIN));
print("BOUNDARY");
print(issubset(IB,K));
print(issubset(K,IB));
print("NILPOTENCE");
print(issubset(K2,IT));
print(issubset(K3,IT));
print(issubset(K4,IT));
print("CONDUCTOR_POWERS");
print(issubset(conductor+K,conductorScheme2));
print(issubset(conductorScheme2,conductor+K));
print(issubset(conductor+K,conductorScheme3));
print(issubset(conductorScheme3,conductor+K));
print(issubset(conductor+K,conductorScheme4));
print(issubset(conductorScheme4,conductor+K));
print("CONDUCTOR_SLICE");
print(dim(std(conductorSlice)));
print(vdim(std(conductorSlice)));
print({len(parameters)}-rank(sliceJacobian));
print("K_ADIC_LENGTHS");
print(vdim(std(conductorSliceK)));
print(vdim(std(conductorSliceK2)));
print(vdim(std(conductorSliceK3)));
print(vdim(std(conductorSliceK4)));
print("MAXIMAL_ADIC_LENGTHS");
print(vdim(std(conductorSlice+maximalIdeal)));
print(vdim(std(conductorSlice+maximalIdeal2)));
print(vdim(std(conductorSlice+maximalIdeal3)));
print(vdim(std(conductorSlice+maximalIdeal4)));
print(vdim(std(conductorSlice+maximalIdeal5)));
print("SOCLE_DIMENSION");
print(vdim(std(conductorSlice))-vdim(std(socleIdeal)));
{"".join(tangent_commands)}
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
    structural_obstruction = (
        "THICK 1 0 "
        "THIN 1 0 "
        "BOUNDARY 1 0 "
        "NILPOTENCE 0 0 0 "
    )
    assert structural_obstruction in compact, result.stdout + result.stderr
    return result.stdout


def coefficient_point(
    factor_variables: tuple[tuple[sp.Symbol, ...], ...],
    translation: int,
    parameter: int,
) -> dict[sp.Symbol, sp.Expr]:
    factors = dickson_vertex_factors(
        WORD,
        W,
        sp.Integer(translation),
        sp.Integer(parameter),
    )
    return {
        variable: sp.Poly(factor, W).nth(power)
        for factor, variables in zip(factors, factor_variables)
        for power, variable in enumerate(variables, 1)
    }


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
    assert all(residuals.values())

    # Both paths contain the endpoint cuts 3 and 21.  The thick path adds
    # cut 7 and hence omits composite cut 6; the thin path does the reverse.
    endpoint = residuals[3] + residuals[21]
    thick = endpoint + residuals[7]
    thin = endpoint + residuals[6]
    boundary = endpoint + residuals[6] + residuals[7]
    graph, parameter_numerator, inner_quadratic = dickson_graph(factor_variables)
    output = singular_audit(
        parameters,
        thick,
        thin,
        boundary,
        graph,
        parameter_numerator,
        inner_quadratic,
        (
            coefficient_point(factor_variables, 1, 2),
            coefficient_point(factor_variables, 1, 0),
        ),
    )
    print(output)
    print("RESULT: ordinary polynomial residuals were used, including degree one")
    print("RESULT: all three path/boundary ideals contain the Dickson graph")
    print("RESULT: neither half-path nor the boundary equals that graph globally")
    print("RESULT: naive global spectator independence fails before component isolation")


if __name__ == "__main__":
    main()
