#!/usr/bin/env python3
"""Verify the degree-thirty Ritt braid as a coefficient-decorated 2-complex.

Besides filling the six-vertex Ritt graph by its braid two-cell, this script
compares the two path incidence ideals on the normalized ``2 o 3 o 5`` chart.
It proves that their reductions and normalizations agree, but that one path
retains a nilpotent direction supported on the monomial divisor ``z=0``.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import re
from math import comb
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.ritt_complex import (
    MoveType,
    compose_factors,
    degree_thirty_braid_decorations,
    dickson,
    dickson_vertex_factors,
    permutation_ritt_complex,
    symmetric_braid_complex,
)
from explore_degree30_hessian_ritt_braid import (
    PARAMETERS,
    W,
    Z,
    canonical_residuals,
    p1,
    q1,
    q2,
    r1,
    r2,
    r3,
    r4,
)


def serialize_ideal(equations: list[sp.Expr]) -> str:
    return ",".join(str(equation).replace("**", "^") for equation in equations)


def graph_ideal() -> list[sp.Expr]:
    """Return the graph ideal of the Dickson ``A^2`` in the base chart."""

    translation = r4 / 5
    parameter = (2 * r4**2 - 5 * r3) / 25
    c5 = dickson(5, translation, parameter)
    c15 = dickson(15, translation, parameter)
    substitutions = {
        r2: 10 * translation**3 - 15 * parameter * translation,
        r1: 5 * translation**4
        - 15 * parameter * translation**2
        + 5 * parameter**2,
        q2: 3 * c5,
        q1: 3 * (c5**2 - parameter**5),
        p1: 2 * c15,
    }
    return [
        sp.together(variable - value).as_numer_denom()[0]
        for variable, value in substitutions.items()
    ]


def composition_tangent_matrix(
    outer_degree: int,
    translation: int,
    parameter: int,
) -> sp.Matrix:
    """Return the degree-30 normalized-composition tangent image."""

    inner_degree = 30 // outer_degree
    outer_parameters = sp.symbols(
        f"ta{outer_degree}_{inner_degree}_1:{outer_degree}"
    )
    inner_parameters = sp.symbols(
        f"tb{outer_degree}_{inner_degree}_1:{inner_degree}"
    )
    outer = Z**outer_degree + sum(
        outer_parameters[power - 1] * Z**power
        for power in range(1, outer_degree)
    )
    inner = W**inner_degree + sum(
        inner_parameters[power - 1] * W**power
        for power in range(1, inner_degree)
    )
    composition = sp.Poly(sp.expand(outer.subs(Z, inner)), W)
    factor_parameters = outer_parameters + inner_parameters
    jacobian = sp.Matrix(
        [
            [
                sp.diff(composition.nth(coefficient_degree), variable)
                for variable in factor_parameters
            ]
            for coefficient_degree in range(1, 30)
        ]
    )
    dickson_outer, dickson_inner = dickson_vertex_factors(
        (outer_degree, inner_degree),
        W,
        sp.Integer(translation),
        sp.Integer(parameter),
    )
    specialization = {
        outer_parameters[power - 1]: sp.Poly(dickson_outer, W).nth(power)
        for power in range(1, outer_degree)
    }
    specialization.update(
        {
            inner_parameters[power - 1]: sp.Poly(dickson_inner, W).nth(power)
            for power in range(1, inner_degree)
        }
    )
    return jacobian.subs(specialization)


def exact_rank(matrix: sp.Matrix) -> int:
    return DomainMatrix.from_Matrix(matrix).rank()


def tangent_intersection_dimension(matrices: list[sp.Matrix]) -> int:
    """Intersect injectively parametrized tangent images by a block kernel."""

    ambient_dimension = matrices[0].rows
    widths = [matrix.cols for matrix in matrices]
    assert [exact_rank(matrix) for matrix in matrices] == widths
    offsets = []
    total_width = 0
    for width in widths:
        offsets.append(total_width)
        total_width += width
    equations = sp.zeros(ambient_dimension * (len(matrices) - 1), total_width)
    for index in range(1, len(matrices)):
        row = (index - 1) * ambient_dimension
        equations[row : row + ambient_dimension, 0 : widths[0]] = matrices[0]
        column = offsets[index]
        equations[
            row : row + ambient_dimension,
            column : column + widths[index],
        ] = -matrices[index]
    return total_width - exact_rank(equations)


def ambient_tangent_audit() -> None:
    """Check that the path asymmetry survives before choosing the base chart."""

    cuts = (2, 3, 5, 6, 10, 15)
    path_cuts = {
        "endpoint": (2, 5, 6, 15),
        "first": (2, 3, 5, 6, 15),
        "second": (2, 5, 6, 10, 15),
        "boundary": cuts,
    }
    expected = {
        (1, 2): {"endpoint": 2, "first": 2, "second": 2, "boundary": 2},
        (1, 0): {"endpoint": 3, "first": 3, "second": 2, "boundary": 2},
    }
    expected_omissions = {
        (1, 2): {cut: 2 for cut in cuts},
        (1, 0): {2: 2, 3: 2, 5: 2, 6: 3, 10: 3, 15: 3},
    }
    for point, expected_dimensions in expected.items():
        polynomial_tangents = {
            cut: composition_tangent_matrix(cut, *point) for cut in cuts
        }
        for forget_linear_term in (False, True):
            tangents = {
                cut: (
                    matrix[1:, :] if forget_linear_term else matrix
                )
                for cut, matrix in polynomial_tangents.items()
            }
            dimensions = {
                name: tangent_intersection_dimension(
                    [tangents[cut] for cut in selected_cuts]
                )
                for name, selected_cuts in path_cuts.items()
            }
            assert dimensions == expected_dimensions
            omission_dimensions = {
                omitted_cut: tangent_intersection_dimension(
                    [
                        tangents[cut]
                        for cut in cuts
                        if cut != omitted_cut
                    ]
                )
                for omitted_cut in cuts
            }
            assert omission_dimensions == expected_omissions[point]


def singular_scheme_audit(
    endpoint: list[sp.Expr],
    first_path: list[sp.Expr],
    second_path: list[sp.Expr],
    boundary: list[sp.Expr],
    polynomial_endpoint: list[sp.Expr],
    polynomial_first_path: list[sp.Expr],
    polynomial_second_path: list[sp.Expr],
    polynomial_boundary: list[sp.Expr],
    dickson_graph: list[sp.Expr],
) -> str:
    """Run exact ideal, conductor, and tangent comparisons in Singular."""

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for the scheme audit"
    ideals = {
        "IE": endpoint,
        "IL": first_path,
        "IR": second_path,
        "IB": boundary,
        "IPE": polynomial_endpoint,
        "IPL": polynomial_first_path,
        "IPR": polynomial_second_path,
        "IPB": polynomial_boundary,
        "K": dickson_graph,
    }
    declarations = "\n".join(
        f"ideal {name}={serialize_ideal(equations)};"
        for name, equations in ideals.items()
    )
    program = f"""
ring q=0,({",".join(map(str, PARAMETERS))}),dp;
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
proc evaluatedrank(ideal I, int vp1, int vq1, int vq2, int vr1,
                   int vr2, int vr3, int vr4)
{{
  matrix J=jacob(I);
  J=subst(J,p1,vp1);
  J=subst(J,q1,vq1);
  J=subst(J,q2,vq2);
  J=subst(J,r1,vr1);
  J=subst(J,r2,vr2);
  J=subst(J,r3,vr3);
  J=subst(J,r4,vr4);
  return(rank(J));
}}
{declarations}
ideal Kgraph=K;
ideal quotientLeft=IL;
ideal quotientRight=K;
ideal K2=Kgraph*Kgraph;
ideal K3=K2*Kgraph;
ideal K4=K3*Kgraph;
ideal conductor=quotient(quotientLeft,quotientRight);
poly h=-2*r4^2+5*r3;
ideal conductorScheme=Kgraph,h^2;
ideal transverseSlice=IL,h,r4-5;
ideal fixedTranslation=r4-5;
ideal conductorSlice=IL+conductor+fixedTranslation;
ideal conductorSliceK=conductorSlice+Kgraph;
ideal conductorSliceK2=conductorSlice+K2;
ideal conductorSliceK3=conductorSlice+K3;
ideal conductorSliceK4=conductorSlice+K4;
ideal reducedPoint=Kgraph+h+fixedTranslation;

print("ENDPOINT_EQUALS_FIRST");
print(issubset(IE,IL));
print(issubset(IL,IE));
print("SECOND_EQUALS_BOUNDARY");
print(issubset(IR,IB));
print(issubset(IB,IR));
print("POLYNOMIAL_EQUALS_HESSIAN_PATHS");
print(issubset(IPE,IE));
print(issubset(IE,IPE));
print(issubset(IPL,IL));
print(issubset(IL,IPL));
print(issubset(IPR,IR));
print(issubset(IR,IPR));
print(issubset(IPB,IB));
print(issubset(IB,IPB));
print("SECOND_EQUALS_DICKSON");
print(issubset(IR,Kgraph));
print(issubset(Kgraph,IR));
print("FIRST_LIES_INSIDE_DICKSON");
print(issubset(IL,Kgraph));
print("DICKSON_LIES_INSIDE_FIRST");
print(issubset(Kgraph,IL));
print("SQUARE_ZERO");
print(issubset(K2,IL));
print("CUBE_ZERO");
print(issubset(K3,IL));
print("FOURTH_POWER_ZERO");
print(issubset(K4,IL));
print("CONDUCTOR_IS_DOUBLE_MONOMIAL_DIVISOR");
print(issubset(conductor+Kgraph,conductorScheme));
print(issubset(conductorScheme,conductor+Kgraph));
print("GENERIC_TANGENT_DIMENSIONS");
print(7-evaluatedrank(IL,550,267,33,-5,-20,0,5));
print(7-evaluatedrank(IR,550,267,33,-5,-20,0,5));
print("MONOMIAL_TANGENT_DIMENSIONS");
print(7-evaluatedrank(IL,2,3,3,5,10,10,5));
print(7-evaluatedrank(IR,2,3,3,5,10,10,5));
print("TRANSVERSE_SLICE");
print(dim(std(transverseSlice)));
print(vdim(std(transverseSlice)));
print(7-evaluatedrank(transverseSlice,2,3,3,5,10,10,5));
print("CONDUCTOR_SLICE");
print(dim(std(conductorSlice)));
print(vdim(std(conductorSlice)));
print(7-evaluatedrank(conductorSlice,2,3,3,5,10,10,5));
print("CONDUCTOR_SLICE_K_ADIC_LENGTHS");
print(vdim(std(conductorSliceK)));
print(vdim(std(conductorSliceK2)));
print(vdim(std(conductorSliceK3)));
print(vdim(std(conductorSliceK4)));
print("REDUCED_POINT_LENGTH");
print(vdim(std(reducedPoint)));
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=3600,
    )
    compact = " ".join(result.stdout.split())
    expected = (
        "ENDPOINT_EQUALS_FIRST 1 1 "
        "SECOND_EQUALS_BOUNDARY 1 1 "
        "POLYNOMIAL_EQUALS_HESSIAN_PATHS 1 1 1 1 1 1 1 1 "
        "SECOND_EQUALS_DICKSON 1 1 "
        "FIRST_LIES_INSIDE_DICKSON 1 "
        "DICKSON_LIES_INSIDE_FIRST 0 "
        "SQUARE_ZERO 0 "
        "CUBE_ZERO 0 "
        "FOURTH_POWER_ZERO 1 "
        "CONDUCTOR_IS_DOUBLE_MONOMIAL_DIVISOR 1 1 "
        "GENERIC_TANGENT_DIMENSIONS 2 2 "
        "MONOMIAL_TANGENT_DIMENSIONS 3 2 "
        "TRANSVERSE_SLICE 0 2 1 "
        "CONDUCTOR_SLICE 0 5 1 "
        "CONDUCTOR_SLICE_K_ADIC_LENGTHS 2 4 5 5 "
        "REDUCED_POINT_LENGTH 1"
    )
    assert expected in compact, result.stdout + result.stderr
    return result.stdout


def opposite_chart_scheme_audit() -> str:
    """Repeat the exact path comparison on the opposite ``5 o 3 o 2`` chart."""

    o1, o2, o3, o4, m1, m2, n1 = sp.symbols(
        "op1 op2 op3 op4 mid1 mid2 inn1"
    )
    parameters = (o1, o2, o3, o4, m1, m2, n1)
    outer = Z**5 + o4 * Z**4 + o3 * Z**3 + o2 * Z**2 + o1 * Z
    middle = Z**3 + m2 * Z**2 + m1 * Z
    inner = W**2 + n1 * W
    polynomial = sp.expand(outer.subs(Z, middle.subs(Z, inner)))
    polynomial_residuals = {
        cut: canonical_residuals(
            polynomial,
            cut,
            30 // cut,
            parameters=parameters,
            factor_output=False,
            minimum_coefficient_degree=1,
        )
        for cut in (2, 3, 6, 10)
    }
    assert all(equations[0] != 0 for equations in polynomial_residuals.values())
    residuals = {
        cut: equations[1:] for cut, equations in polynomial_residuals.items()
    }
    first_path = residuals[2] + residuals[3] + residuals[6]
    second_path = residuals[2] + residuals[6] + residuals[10]
    boundary = first_path + residuals[10]
    polynomial_first_path = (
        polynomial_residuals[2]
        + polynomial_residuals[3]
        + polynomial_residuals[6]
    )
    polynomial_second_path = (
        polynomial_residuals[2]
        + polynomial_residuals[6]
        + polynomial_residuals[10]
    )
    polynomial_boundary = polynomial_first_path + polynomial_residuals[10]

    translation = n1 / 2
    parameter = (translation**2 - m2 / 3) / 2
    factors = dickson_vertex_factors(
        (5, 3, 2), W, translation, parameter
    )
    chart_variables = ((o1, o2, o3, o4), (m1, m2), (n1,))
    coefficient_map = {}
    for factor, variables in zip(factors, chart_variables):
        factor_polynomial = sp.Poly(factor, W)
        coefficient_map.update(
            {
                variable: factor_polynomial.nth(power)
                for power, variable in enumerate(variables, 1)
            }
        )
    graph = [
        sp.together(variable - value).as_numer_denom()[0]
        for variable, value in coefficient_map.items()
        if variable not in {m2, n1}
    ]
    monomial_parameter = 3 * n1**2 - 4 * m2

    singular = shutil.which("Singular")
    assert singular is not None
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
ideal IL={serialize_ideal(first_path)};
ideal IR={serialize_ideal(second_path)};
ideal IB={serialize_ideal(boundary)};
ideal IPL={serialize_ideal(polynomial_first_path)};
ideal IPR={serialize_ideal(polynomial_second_path)};
ideal IPB={serialize_ideal(polynomial_boundary)};
ideal K={serialize_ideal(graph)};
ideal K2=K*K;
ideal K3=K2*K;
ideal K4=K3*K;
ideal quotientLeft=IL;
ideal quotientRight=K;
ideal conductor=quotient(quotientLeft,quotientRight);
ideal conductorScheme=K,({str(monomial_parameter).replace("**", "^")})^2;
print("OPPOSITE_LEFT_IN_DICKSON");
print(issubset(IL,K));
print("OPPOSITE_DICKSON_IN_LEFT");
print(issubset(K,IL));
print("OPPOSITE_RIGHT_EQUALS_DICKSON");
print(issubset(IR,K));
print(issubset(K,IR));
print("OPPOSITE_BOUNDARY_EQUALS_DICKSON");
print(issubset(IB,K));
print(issubset(K,IB));
print("OPPOSITE_POLYNOMIAL_EQUALS_HESSIAN");
print(issubset(IPL,IL));
print(issubset(IL,IPL));
print(issubset(IPR,IR));
print(issubset(IR,IPR));
print(issubset(IPB,IB));
print(issubset(IB,IPB));
print("OPPOSITE_NILPOTENCE");
print(issubset(K2,IL));
print(issubset(K3,IL));
print(issubset(K4,IL));
print("OPPOSITE_CONDUCTOR");
print(issubset(conductor+K,conductorScheme));
print(issubset(conductorScheme,conductor+K));
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=3600,
    )
    compact = " ".join(result.stdout.split())
    expected = (
        "OPPOSITE_LEFT_IN_DICKSON 1 "
        "OPPOSITE_DICKSON_IN_LEFT 0 "
        "OPPOSITE_RIGHT_EQUALS_DICKSON 1 1 "
        "OPPOSITE_BOUNDARY_EQUALS_DICKSON 1 1 "
        "OPPOSITE_POLYNOMIAL_EQUALS_HESSIAN 1 1 1 1 1 1 "
        "OPPOSITE_NILPOTENCE 0 0 1 "
        "OPPOSITE_CONDUCTOR 1 1"
    )
    assert expected in compact, result.stdout + result.stderr
    return result.stdout


def rotated_chart_scheme_audit(
    word: tuple[int, int, int],
    complex_,
    *,
    show_slice_basis: bool = False,
) -> str:
    """Audit one of the four remaining vertices of the braid hexagon."""

    word_label = "".join(map(str, word))
    factor_variables = []
    generic_factors = []
    for position, degree in enumerate(word):
        variables = sp.symbols(
            f"x{word_label}_{position}_1:{degree}"
        )
        factor_variables.append(variables)
        generic_factors.append(
            W**degree
            + sum(
                variables[power - 1] * W**power
                for power in range(1, degree)
            )
        )
    parameters = tuple(
        variable for variables in factor_variables for variable in variables
    )
    polynomial = compose_factors(tuple(generic_factors), W)

    vertex_by_word = {vertex.word: vertex for vertex in complex_.vertices}
    adjacency = {vertex.word: [] for vertex in complex_.vertices}
    for edge in complex_.edges:
        left, right = edge.endpoints
        adjacency[left].append(right)
        adjacency[right].append(left)
    endpoint = tuple(reversed(word))
    paths = []

    def extend(path: tuple[tuple[int, int, int], ...]) -> None:
        if len(path) == 4:
            if path[-1] == endpoint:
                paths.append(path)
            return
        for neighbor in adjacency[path[-1]]:
            if neighbor not in path:
                extend(path + (neighbor,))

    extend((word,))
    assert len(paths) == 2
    all_cuts = {2, 3, 5, 6, 10, 15}
    base_cuts = set(vertex_by_word[word].cuts)
    path_data = []
    for path in paths:
        path_cuts = set().union(
            *(set(vertex_by_word[vertex].cuts) for vertex in path)
        )
        omitted = all_cuts - path_cuts
        assert len(omitted) == 1
        path_data.append((path, path_cuts, omitted.pop()))
    thick_data = next(data for data in path_data if data[2] in {6, 10, 15})
    thin_data = next(data for data in path_data if data[2] in {2, 3, 5})
    decoration = next(
        item
        for item in degree_thirty_braid_decorations()
        if item.composite_omission == thick_data[2]
    )
    assert decoration.prime_omission == thin_data[2]

    requested_cuts = tuple(sorted(all_cuts - base_cuts))
    polynomial_residuals = {
        cut: canonical_residuals(
            polynomial,
            cut,
            30 // cut,
            parameters=parameters,
            factor_output=False,
            minimum_coefficient_degree=1,
        )
        for cut in requested_cuts
    }
    assert all(equations[0] != 0 for equations in polynomial_residuals.values())
    hessian_residuals = {
        cut: equations[1:] for cut, equations in polynomial_residuals.items()
    }

    def equations_for(
        path_cuts: set[int],
        residuals_by_cut: dict[int, list[sp.Expr]],
    ) -> list[sp.Expr]:
        return [
            equation
            for cut in sorted(path_cuts - base_cuts)
            for equation in residuals_by_cut[cut]
        ]

    thick_hessian = equations_for(thick_data[1], hessian_residuals)
    thin_hessian = equations_for(thin_data[1], hessian_residuals)
    boundary_hessian = equations_for(all_cuts, hessian_residuals)
    thick_polynomial = equations_for(thick_data[1], polynomial_residuals)
    thin_polynomial = equations_for(thin_data[1], polynomial_residuals)
    boundary_polynomial = equations_for(all_cuts, polynomial_residuals)

    inner_degree = word[-1]
    inner_top = factor_variables[-1][-1]
    translation = inner_top / inner_degree
    if inner_degree >= 3:
        inner_next = factor_variables[-1][-2]
        parameter = (
            comb(inner_degree, 2) * translation**2 - inner_next
        ) / inner_degree
        free_variables = {inner_top, inner_next}
    else:
        preceding_degree = word[-2]
        preceding_top = factor_variables[-2][-1]
        parameter = (
            translation**2 - preceding_top / preceding_degree
        ) / 2
        free_variables = {inner_top, preceding_top}
    dickson_factors = dickson_vertex_factors(
        word, W, translation, parameter
    )
    coefficient_map = {}
    for factor, variables in zip(dickson_factors, factor_variables):
        factor_polynomial = sp.Poly(factor, W)
        coefficient_map.update(
            {
                variable: factor_polynomial.nth(power)
                for power, variable in enumerate(variables, 1)
            }
        )
    graph = [
        sp.together(variable - value).as_numer_denom()[0]
        for variable, value in coefficient_map.items()
        if variable not in free_variables
    ]
    monomial_parameter = sp.together(parameter).as_numer_denom()[0]
    monomial_point = {
        variable: comb(degree, power)
        for degree, variables in zip(word, factor_variables)
        for power, variable in enumerate(variables, 1)
    }
    tangent_substitutions = "\n".join(
        f"sliceJacobian=subst(sliceJacobian,{variable},{value});"
        for variable, value in monomial_point.items()
    )
    slice_presentation_commands = ""
    slice_presentation_name = None
    if word == (3, 2, 5):
        first_coordinate = factor_variables[-1][-2]
        second_coordinate = factor_variables[-1][-3]
        slice_coordinates = (first_coordinate, second_coordinate)
        local_first = f"({first_coordinate}-10)"
        local_second = f"({second_coordinate}-10)"
        presentation_generators = (
            f"{local_first}^2,"
            f"{local_second}^2-6*{local_first}*{local_second}"
        )
        slice_presentation_name = "Q[u,v]/(u^2,v^2)"
    elif word == (3, 5, 2):
        first_coordinate = factor_variables[-2][-1]
        second_coordinate = factor_variables[-2][-3]
        slice_coordinates = (first_coordinate, second_coordinate)
        local_first = f"({first_coordinate}-5)"
        local_second = (
            f"({second_coordinate}-10-6*({first_coordinate}-5))"
        )
        presentation_generators = (
            f"{local_first}^4,"
            f"500*{local_second}^2"
            f"-1050*{local_first}^2*{local_second}"
            f"-53*{local_first}^3*{local_second}"
        )
        slice_presentation_name = "Q[u,v]/(u^4,v^2)"
    if slice_presentation_name is not None:
        eliminated_variables = [
            variable for variable in parameters if variable not in slice_coordinates
        ]
        elimination_monomial = "*".join(map(str, eliminated_variables))
        slice_presentation_commands = (
            "ideal sliceElimination=eliminate(conductorSlice,"
            f"{elimination_monomial}); "
            f"ideal slicePresentation={presentation_generators}; "
            'print("ROTATED_SLICE_PRESENTATION"); '
            "print(issubset(sliceElimination,slicePresentation)); "
            "print(issubset(slicePresentation,sliceElimination));"
        )
        if show_slice_basis:
            slice_presentation_commands += (
                ' print("ROTATED_CONDUCTOR_SLICE_ELIMINATION"); '
                "print(std(sliceElimination));"
            )
    elif show_slice_basis:
        raise ValueError(
            "explicit slice coordinates are only selected for 325 and 352"
        )

    singular = shutil.which("Singular")
    assert singular is not None
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
ideal IT={serialize_ideal(thick_hessian)};
ideal ITHIN={serialize_ideal(thin_hessian)};
ideal IB={serialize_ideal(boundary_hessian)};
ideal IPT={serialize_ideal(thick_polynomial)};
ideal IPTHIN={serialize_ideal(thin_polynomial)};
ideal IPB={serialize_ideal(boundary_polynomial)};
ideal K={serialize_ideal(graph)};
ideal K2=K*K;
ideal K3=K2*K;
ideal K4=K3*K;
ideal quotientThick=IT;
ideal quotientGraph=K;
ideal conductor=quotient(quotientThick,quotientGraph);
ideal conductorScheme2=K,({str(monomial_parameter).replace("**", "^")})^2;
ideal conductorScheme3=K,({str(monomial_parameter).replace("**", "^")})^3;
ideal conductorScheme4=K,({str(monomial_parameter).replace("**", "^")})^4;
ideal fixedTranslation={inner_top}-{inner_degree};
ideal conductorSlice=IT+conductor+fixedTranslation;
ideal conductorSliceK=conductorSlice+K;
ideal conductorSliceK2=conductorSlice+K2;
ideal conductorSliceK3=conductorSlice+K3;
ideal conductorSliceK4=conductorSlice+K4;
ideal maximalIdeal=K,({str(monomial_parameter).replace("**", "^")}),fixedTranslation;
ideal maximalIdeal2=maximalIdeal*maximalIdeal;
ideal maximalIdeal3=maximalIdeal2*maximalIdeal;
ideal maximalIdeal4=maximalIdeal3*maximalIdeal;
ideal maximalIdeal5=maximalIdeal4*maximalIdeal;
ideal socleLeft=conductorSlice;
ideal socleRight=maximalIdeal;
ideal socleIdeal=quotient(socleLeft,socleRight);
print("ROTATED_THICK");
print(issubset(IT,K));
print(issubset(K,IT));
print("ROTATED_THIN");
print(issubset(ITHIN,K));
print(issubset(K,ITHIN));
print("ROTATED_BOUNDARY");
print(issubset(IB,K));
print(issubset(K,IB));
print("ROTATED_SYNCHRONIZATION");
print(issubset(IPT,IT));
print(issubset(IT,IPT));
print(issubset(IPTHIN,ITHIN));
print(issubset(ITHIN,IPTHIN));
print(issubset(IPB,IB));
print(issubset(IB,IPB));
print("ROTATED_NILPOTENCE");
print(issubset(K2,IT));
print(issubset(K3,IT));
print(issubset(K4,IT));
print("ROTATED_CONDUCTOR_POWERS");
print(issubset(conductor+K,conductorScheme2));
print(issubset(conductorScheme2,conductor+K));
print(issubset(conductor+K,conductorScheme3));
print(issubset(conductorScheme3,conductor+K));
print(issubset(conductor+K,conductorScheme4));
print(issubset(conductorScheme4,conductor+K));
print("ROTATED_CONDUCTOR_SLICE");
print(dim(std(conductorSlice)));
print(vdim(std(conductorSlice)));
matrix sliceJacobian=jacob(conductorSlice);
{tangent_substitutions}
print({len(parameters)}-rank(sliceJacobian));
print("ROTATED_K_ADIC_LENGTHS");
print(vdim(std(conductorSliceK)));
print(vdim(std(conductorSliceK2)));
print(vdim(std(conductorSliceK3)));
print(vdim(std(conductorSliceK4)));
print("ROTATED_MAXIMAL_ADIC_LENGTHS");
print(vdim(std(conductorSlice+maximalIdeal)));
print(vdim(std(conductorSlice+maximalIdeal2)));
print(vdim(std(conductorSlice+maximalIdeal3)));
print(vdim(std(conductorSlice+maximalIdeal4)));
print(vdim(std(conductorSlice+maximalIdeal5)));
print("ROTATED_SOCLE_DIMENSION");
print(vdim(std(conductorSlice))-vdim(std(socleIdeal)));
{slice_presentation_commands}
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=3600,
    )
    compact = " ".join(result.stdout.split())
    expected_without_nilpotence = (
        "ROTATED_THICK 1 0 "
        "ROTATED_THIN 1 1 "
        "ROTATED_BOUNDARY 1 1 "
        "ROTATED_SYNCHRONIZATION 1 1 1 1 1 1 "
    )
    diagnostic = (
        f"word={word}, thick_omission={thick_data[2]}, "
        f"thin_omission={thin_data[2]}\n"
        + result.stdout
        + result.stderr
    )
    assert expected_without_nilpotence in compact, diagnostic
    nilpotence_match = re.search(
        r"ROTATED_NILPOTENCE ([01]) ([01]) ([01])",
        compact,
    )
    assert nilpotence_match is not None, compact
    nilpotence = tuple(map(int, nilpotence_match.groups()))
    assert nilpotence[0] == 0 and nilpotence[-1] == 1, diagnostic
    nilpotence_index = 3 if nilpotence[1] else 4
    conductor_match = re.search(
        r"ROTATED_CONDUCTOR_POWERS "
        r"([01]) ([01]) ([01]) ([01]) ([01]) ([01])",
        compact,
    )
    assert conductor_match is not None, diagnostic
    conductor_tests = tuple(map(int, conductor_match.groups()))
    conductor_exponent = next(
        (
            exponent
            for exponent, pair in zip(
                (2, 3, 4),
                (
                    conductor_tests[0:2],
                    conductor_tests[2:4],
                    conductor_tests[4:6],
                ),
            )
            if pair == (1, 1)
        ),
        None,
    )
    assert conductor_exponent is not None, diagnostic
    slice_match = re.search(
        r"ROTATED_CONDUCTOR_SLICE 0 ([0-9]+) ([0-9]+) "
        r"ROTATED_K_ADIC_LENGTHS "
        r"([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) "
        r"ROTATED_MAXIMAL_ADIC_LENGTHS "
        r"([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) "
        r"ROTATED_SOCLE_DIMENSION ([0-9]+)",
        compact,
    )
    assert slice_match is not None, diagnostic
    parsed_slice = tuple(map(
        int, slice_match.groups()
    ))
    slice_length, embedding_dimension = parsed_slice[:2]
    k_adic_lengths = parsed_slice[2:6]
    maximal_adic_lengths = parsed_slice[6:11]
    socle_dimension = parsed_slice[11]
    assert socle_dimension == 1, diagnostic
    observed_hilbert_vector = (
        maximal_adic_lengths[0],
        *(
            right - left
            for left, right in zip(
                maximal_adic_lengths,
                maximal_adic_lengths[1:],
            )
            if right > left
        ),
    )
    assert nilpotence_index == decoration.nilpotence_index, diagnostic
    assert conductor_exponent == decoration.conductor_power, diagnostic
    assert slice_length == decoration.transverse_slice.length, diagnostic
    assert (
        embedding_dimension
        == decoration.transverse_slice.embedding_dimension
    ), diagnostic
    assert (
        observed_hilbert_vector
        == decoration.transverse_slice.hilbert_vector
    ), diagnostic
    if slice_presentation_name is not None:
        assert "ROTATED_SLICE_PRESENTATION 1 1" in compact, diagnostic
    summary = (
        f"word {word_label}: composite omission {thick_data[2]} is thick; "
        f"prime omission {thin_data[2]} is reduced; "
        f"nilpotence index {nilpotence_index}; "
        f"annihilator (z^{conductor_exponent}); "
        f"slice length/embedding {slice_length}/{embedding_dimension}; "
        f"K-adic {tuple(k_adic_lengths)}; "
        f"maximal-adic {tuple(maximal_adic_lengths)}; "
        f"socle {socle_dimension}"
    )
    if slice_presentation_name is not None:
        summary += f"; presentation {slice_presentation_name}"
    if show_slice_basis:
        return summary + "\n" + result.stdout
    return summary


def main() -> None:
    rank_four = permutation_ritt_complex((2, 3, 5, 7), MoveType.CHEBYSHEV)
    relation_counts = {
        relation: sum(cell.relation == relation for cell in rank_four.two_cells)
        for relation in ("commuting", "braid")
    }
    assert (len(rank_four.vertices), len(rank_four.edges)) == (24, 36)
    assert relation_counts == {"commuting": 6, "braid": 8}
    assert rank_four.euler_characteristic == 2

    complex_ = symmetric_braid_complex((2, 3, 5), MoveType.CHEBYSHEV)
    assert len(complex_.vertices) == 6
    assert len(complex_.edges) == 6
    assert len(complex_.two_cells) == 1
    assert complex_.euler_characteristic == 1
    assert all(vertex.coefficients.dimension == 7 for vertex in complex_.vertices)
    assert all(edge.correspondence.dimension == 2 for edge in complex_.edges)
    decorations = degree_thirty_braid_decorations()
    assert tuple(
        decoration.transverse_slice.hilbert_vector
        for decoration in decorations
    ) == (
        (1, 1, 1, 1, 1),
        (1, 2, 1),
        (1, 2, 2, 2, 1),
    )
    assert tuple(
        decoration.transverse_slice.residue_field_tor_ranks
        for decoration in decorations
    ) == ((1, 1), (1, 2, 1), (1, 2, 1))

    translation, parameter = sp.symbols("t z")
    target = sp.expand(
        dickson(30, W + translation, parameter)
        - dickson(30, translation, parameter)
    )
    vertex_compositions = {}
    for word in complex_.words:
        factors = dickson_vertex_factors(word, W, translation, parameter)
        composition = compose_factors(factors, W)
        assert sp.expand(composition - target) == 0
        vertex_compositions[word] = composition
    for edge in complex_.edges:
        left, right = edge.endpoints
        assert sp.expand(vertex_compositions[left] - vertex_compositions[right]) == 0
    ambient_tangent_audit()

    factor_2 = Z**2 + p1 * Z
    factor_3 = Z**3 + q2 * Z**2 + q1 * Z
    factor_5 = W**5 + r4 * W**4 + r3 * W**3 + r2 * W**2 + r1 * W
    polynomial = sp.expand(factor_2.subs(Z, factor_3.subs(Z, factor_5)))
    polynomial_residuals = {
        cut: canonical_residuals(
            polynomial,
            cut,
            30 // cut,
            factor_output=False,
            minimum_coefficient_degree=1,
        )
        for cut in (3, 5, 10, 15)
    }
    assert all(equations[0] != 0 for equations in polynomial_residuals.values())
    residuals = {
        cut: equations[1:] for cut, equations in polynomial_residuals.items()
    }

    cell = complex_.two_cells[0]
    first_path, second_path = cell.paths
    assert first_path == ((2, 3, 5), (3, 2, 5), (3, 5, 2), (5, 3, 2))
    assert second_path == ((2, 3, 5), (2, 5, 3), (5, 2, 3), (5, 3, 2))

    # The base word already supplies cuts 2 and 6.  The path ideals are the
    # sums of the canonical residual ideals for the remaining vertex cuts.
    endpoint = residuals[5] + residuals[15]
    first_ideal = residuals[3] + endpoint
    second_ideal = endpoint + residuals[10]
    boundary = residuals[3] + endpoint + residuals[10]
    polynomial_endpoint = (
        polynomial_residuals[5] + polynomial_residuals[15]
    )
    polynomial_first = polynomial_residuals[3] + polynomial_endpoint
    polynomial_second = polynomial_endpoint + polynomial_residuals[10]
    polynomial_boundary = (
        polynomial_residuals[3]
        + polynomial_endpoint
        + polynomial_residuals[10]
    )
    output = singular_scheme_audit(
        endpoint,
        first_ideal,
        second_ideal,
        boundary,
        polynomial_endpoint,
        polynomial_first,
        polynomial_second,
        polynomial_boundary,
        graph_ideal(),
    )
    opposite_output = opposite_chart_scheme_audit()
    rotated_outputs = [
        rotated_chart_scheme_audit(word, complex_)
        for word in ((3, 2, 5), (5, 2, 3), (3, 5, 2), (2, 5, 3))
    ]
    print(output)
    print(opposite_output)
    print("\n".join(rotated_outputs))
    print("PASS: the six degree-30 decompositions form a filled S3 braid complex")
    print("PASS: every vertex and edge carries the same Dickson A^2 coefficient map")
    print("PASS: the two braid paths have the same smooth A^2 normalization")
    print("PASS: the first path has nilpotence index four and one excess tangent on z=0")
    print("PASS: the tangent asymmetry survives in polynomial and Hessian coefficient space")
    print("PASS: no unsynchronized first-order Hessian direction occurs at either test point")
    print("PASS: polynomial and Hessian path ideals agree on all six charts")
    print("PASS: exactly the three composite cuts kill the monomial braid excess")
    print("PASS: coordinate and annihilator slices have lengths two and five")
    print("PASS: its defect annihilator restricts scheme-theoretically to (z^2)")
    print("PASS: the same scheme defect occurs on the opposite 5 o 3 o 2 chart")
    print("PASS: the sector-dependent prime/composite rule holds on all six charts")
    print("PASS: the cut-15 slice is Q[u,v]/(u^2,v^2)")
    print("PASS: the cut-6 slice is Q[u,v]/(u^4,v^2)")
    print("PASS: the transverse conormal ranks and Koszul Tor ranks are recorded")


if __name__ == "__main__":
    main()
