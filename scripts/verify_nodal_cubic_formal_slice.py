#!/usr/bin/env python3
"""Exact first formal slice for the universal nodal cubic tensor.

For the nodal ternary-cubic symbol, this checker proves

    ker(C) / im(G_nod) = Q[y,z](-3),

with generator given by the tensor attached to Z^3.  It then decomposes
the fixed 24-dimensional quartic-kernel basis into the direct sum of the
22-dimensional quartic gauge image and a two-dimensional slice.  The
coordinate slice is spanned by the first two fixed basis directions; the
existing full-support sum/alternating-sum plane is another transverse
slice.  It computes the complete degree-five quadratic normal remainder,
quotients its dependence on the five-dimensional kernel of the quartic
gauge lift, and classifies the resulting pure-curvature zero scheme.  This
is a degree-five theorem, not an all-order formal normal form for the
universal nodal family.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import research_universal_cubic_quartic_kernel_saturation as frontier  # noqa: E402
import verify_cubic_formal_gauge_cokernel_atlas as atlas  # noqa: E402
import verify_cubic_symbol_dense_quartic_plane_saturation as dense  # noqa: E402
import verify_cubic_symbol_double_saturation as cubic  # noqa: E402
import verify_universal_cubic_cotangent_saturation as smooth  # noqa: E402


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "nodal_cubic_formal_slice.json"
)


def coefficient_column(vector: sp.Matrix, degree: int) -> sp.Matrix:
    variables = cubic.BASE_VARIABLES
    monomials = cubic.homogeneous_monomials(degree)
    return sp.Matrix(
        [
            sp.Poly(
                sp.expand(vector[row]), *variables
            ).coeff_monomial(monomial)
            for row in range(vector.rows)
            for monomial in monomials
        ]
    )


def rational_solution(
    matrix: sp.Matrix, targets: sp.Matrix
) -> sp.Matrix:
    solution, _parameters = matrix.gauss_jordan_solve(targets)
    free_symbols = set().union(
        *(entry.free_symbols for entry in solution)
    )
    return solution.subs({symbol: 0 for symbol in free_symbols})


def polynomial_lift(solution: sp.Matrix) -> sp.Matrix:
    x, y, z = cubic.BASE_VARIABLES
    assert solution.rows == 27
    return sp.Matrix(
        9,
        solution.cols,
        lambda row, column: sp.expand(
            solution[3 * row, column] * x
            + solution[3 * row + 1, column] * y
            + solution[3 * row + 2, column] * z
        ),
    )


def homogeneous_polynomial_lift(
    solution: sp.Matrix, degree: int
) -> sp.Matrix:
    """Turn coefficient blocks into a 3-by-3 homogeneous matrix."""

    monomials = cubic.homogeneous_monomials(degree)
    assert solution.rows == 9 * len(monomials)
    return sp.Matrix(
        3,
        3,
        lambda row, column: sp.expand(
            sum(
                solution[
                    (3 * row + column) * len(monomials) + index
                ]
                * monomial
                for index, monomial in enumerate(monomials)
            )
        ),
    )


def rational_solution_preserving(
    matrix: sp.Matrix,
    targets: sp.Matrix,
    parameters: tuple[sp.Symbol, ...],
) -> sp.Matrix:
    """Choose zero slack variables without specializing base parameters."""

    solution, _slack = matrix.gauss_jordan_solve(targets)
    free_symbols = set().union(
        *(entry.free_symbols for entry in solution)
    ).difference(parameters)
    return solution.subs({symbol: 0 for symbol in free_symbols})


def matrix_record(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [sp.sstr(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def gauge_velocity(derivation: sp.Matrix) -> sp.Matrix:
    relation_change = derivation * cubic.RELATION
    return sp.Matrix(
        (
            relation_change[2],
            -relation_change[1],
            relation_change[0],
        )
    )


def coordinate_action(
    tensor: sp.Matrix, velocity: sp.Matrix
) -> sp.Matrix:
    variables = cubic.BASE_VARIABLES
    return sp.Matrix(
        [
            sum(
                velocity[index]
                * sp.diff(tensor[component], variables[index])
                for index in range(3)
            )
            for component in range(tensor.rows)
        ]
    )


def slot_action(
    tensor: sp.Matrix, derivation: sp.Matrix
) -> sp.Matrix:
    values: list[sp.Expr] = []
    for triple in smooth.TRIPLES:
        value = sp.Integer(0)
        for position in range(3):
            for replacement in range(3):
                changed = list(triple)
                changed[position] = replacement
                changed_index = smooth.TRIPLES.index(
                    tuple(sorted(changed))
                )
                value += (
                    derivation[replacement, triple[position]]
                    * tensor[changed_index]
                )
        values.append(value)
    return sp.Matrix(values)


def two_slot_action(
    tensor: sp.Matrix,
    first_derivation: sp.Matrix,
    second_derivation: sp.Matrix,
) -> sp.Matrix:
    values: list[sp.Expr] = []
    for triple in smooth.TRIPLES:
        value = sp.Integer(0)
        for first_position in range(3):
            for second_position in range(first_position + 1, 3):
                for first_replacement in range(3):
                    for second_replacement in range(3):
                        changed = list(triple)
                        changed[first_position] = first_replacement
                        changed[second_position] = second_replacement
                        changed_index = smooth.TRIPLES.index(
                            tuple(sorted(changed))
                        )
                        value += (
                            first_derivation[
                                first_replacement,
                                triple[first_position],
                            ]
                            * second_derivation[
                                second_replacement,
                                triple[second_position],
                            ]
                            + second_derivation[
                                first_replacement,
                                triple[first_position],
                            ]
                            * first_derivation[
                                second_replacement,
                                triple[second_position],
                            ]
                        ) * tensor[changed_index]
        values.append(value)
    return sp.Matrix(values)


def second_action_polarization(
    tensor: sp.Matrix,
    first_derivation: sp.Matrix,
    second_derivation: sp.Matrix,
) -> sp.Matrix:
    """Coefficient of ab in the finite action of I+aD+bE."""

    variables = cubic.BASE_VARIABLES
    first_velocity = gauge_velocity(first_derivation)
    second_velocity = gauge_velocity(second_derivation)
    first_trace = sp.trace(first_derivation)
    second_trace = sp.trace(second_derivation)
    coordinate_cross = sp.Matrix(
        [
            sum(
                first_velocity[first_index]
                * second_velocity[second_index]
                * sp.diff(
                    tensor[component],
                    variables[first_index],
                    variables[second_index],
                )
                for first_index in range(3)
                for second_index in range(3)
            )
            for component in range(tensor.rows)
        ]
    )
    first_raw_linear = (
        coordinate_action(tensor, first_velocity)
        + slot_action(tensor, first_derivation)
    )
    second_raw_linear = (
        coordinate_action(tensor, second_velocity)
        + slot_action(tensor, second_derivation)
    )
    result = (
        coordinate_cross
        + slot_action(
            coordinate_action(tensor, second_velocity),
            first_derivation,
        )
        + slot_action(
            coordinate_action(tensor, first_velocity),
            second_derivation,
        )
        + two_slot_action(
            tensor, first_derivation, second_derivation
        )
        - first_trace * second_raw_linear
        - second_trace * first_raw_linear
        + (
            first_trace * second_trace
            + sp.trace(first_derivation * second_derivation)
        )
        * tensor
    )
    return result.applyfunc(sp.expand)


def finite_inverse_gauge_coefficients(
    base_tensor: sp.Matrix,
    perturbation: sp.Matrix,
    derivation: sp.Matrix,
    order: int,
) -> tuple[sp.Matrix, ...]:
    """Expand A_{-tD}(base+t*perturbation) through the given order."""

    scale = sp.Symbol("formal_scale")
    variables = cubic.BASE_VARIABLES
    basis_change = sp.eye(3) - scale * derivation
    velocity = gauge_velocity(derivation)
    transformed_variables = (
        sp.Matrix(variables) - scale * velocity
    )
    substitution = dict(zip(variables, transformed_variables))
    determinant = sp.Poly(
        sp.expand(basis_change.det()), scale
    )
    determinant_coefficients = [
        determinant.coeff_monomial(scale**index)
        for index in range(order + 1)
    ]
    inverse_coefficients = [sp.Integer(1)]
    for current_order in range(1, order + 1):
        inverse_coefficients.append(
            sp.expand(
                -sum(
                    determinant_coefficients[index]
                    * inverse_coefficients[current_order - index]
                    for index in range(
                        1, min(current_order, 3) + 1
                    )
                )
            )
        )

    coefficients = [
        sp.zeros(base_tensor.rows, 1)
        for _ in range(order + 1)
    ]
    for row, triple in enumerate(smooth.TRIPLES):
        numerator = sp.Integer(0)
        for first in range(3):
            for second in range(3):
                for third in range(3):
                    changed_index = smooth.TRIPLES.index(
                        tuple(sorted((first, second, third)))
                    )
                    numerator += (
                        basis_change[first, triple[0]]
                        * basis_change[second, triple[1]]
                        * basis_change[third, triple[2]]
                        * (
                            base_tensor[changed_index]
                            + scale * perturbation[changed_index]
                        ).subs(substitution, simultaneous=True)
                    )
        numerator_polynomial = sp.Poly(
            sp.expand(numerator), scale
        )
        numerator_coefficients = [
            numerator_polynomial.coeff_monomial(scale**index)
            for index in range(order + 1)
        ]
        for current_order in range(order + 1):
            coefficients[current_order][row] = sp.expand(
                sum(
                    numerator_coefficients[index]
                    * inverse_coefficients[current_order - index]
                    for index in range(current_order + 1)
                )
            )
    return tuple(coefficients)


def degree_five_curvature(
    nodal_tensor: sp.Matrix,
    gauge: sp.Matrix,
    generator: sp.Matrix,
    directions: sp.Matrix,
    gauge_lift: sp.Matrix,
    parameters: tuple[sp.Symbol, ...],
) -> tuple[sp.Matrix, sp.Matrix, tuple[int, ...], int]:
    """Project the deterministic first gauge's degree-five remainder."""

    y = cubic.y
    z = cubic.z
    quadratic_monomials = cubic.homogeneous_monomials(2)
    action_degree_five = sp.Matrix.hstack(
        *[
            coefficient_column(
                gauge[:, column] * monomial, 5
            )
            for column in range(gauge.cols)
            for monomial in quadratic_monomials
        ]
    )
    action_pivots = action_degree_five.rref()[1]
    quotient_basis = sp.Matrix.hstack(
        *[
            coefficient_column(monomial * generator, 5)
            for monomial in (y**2, y * z, z**2)
        ]
    )
    reduced_basis = sp.Matrix.hstack(
        *[
            action_degree_five[:, column]
            for column in action_pivots
        ],
        quotient_basis[:, 0],
        quotient_basis[:, 1],
        quotient_basis[:, 2],
    )
    assert len(action_pivots) == 39
    assert reduced_basis.rank() == 42
    coefficient_rows = reduced_basis.T.rref()[1]
    assert len(coefficient_rows) == 42
    quotient_projection = (
        reduced_basis[list(coefficient_rows), :].inv()[-3:, :]
    )

    def project(vector: sp.Matrix) -> sp.Matrix:
        coefficients = coefficient_column(vector, 5)
        return quotient_projection * coefficients[
            list(coefficient_rows), :
        ]

    derivations = [
        sp.Matrix(3, 3, list(gauge_lift[:, column]))
        for column in range(gauge_lift.cols)
    ]
    direction_gauges = [
        smooth.gauge_matrix(directions[:, column])
        for column in range(directions.cols)
    ]
    forms = [sp.Integer(0), sp.Integer(0), sp.Integer(0)]
    nonzero_cross_terms = 0
    for first in range(len(parameters)):
        diagonal = (
            second_action_polarization(
                nodal_tensor,
                derivations[first],
                derivations[first],
            )
            / 2
            - direction_gauges[first]
            * sp.Matrix(list(derivations[first]))
        )
        diagonal_projection = project(diagonal)
        for component in range(3):
            forms[component] += (
                diagonal_projection[component] * parameters[first] ** 2
            )
        for second in range(first + 1, len(parameters)):
            cross = (
                second_action_polarization(
                    nodal_tensor,
                    derivations[first],
                    derivations[second],
                )
                - direction_gauges[first]
                * sp.Matrix(list(derivations[second]))
                - direction_gauges[second]
                * sp.Matrix(list(derivations[first]))
            )
            cross_projection = project(cross)
            if any(cross_projection):
                nonzero_cross_terms += 1
            for component in range(3):
                forms[component] += (
                    cross_projection[component]
                    * parameters[first]
                    * parameters[second]
                )
    return (
        sp.Matrix([sp.factor(form) for form in forms]),
        quotient_projection,
        tuple(coefficient_rows),
        nonzero_cross_terms,
    )


def singular_module_certificate(
    compatibility: sp.Matrix,
    gauge: sp.Matrix,
    generator: sp.Matrix,
) -> dict[str, Any]:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    program = f"""
ring nodal_slice=0,(x,y,z),dp;
module C={atlas.singular_module(compatibility)};
module G={atlas.singular_module(gauge)};
module ETA={atlas.singular_module(generator)};
module K=syz(C);
module H=G,ETA;
H=std(H);
module DIFFERENCE=simplify(reduce(K,H),2);
module Q=std(modulo(K,G));
ideal ANN=std(quotient(G,K));
print("@@KERNEL_GENERATORS="+string(size(K)));
print("@@GENERATION_DIFFERENCE="+string(size(DIFFERENCE)));
print("@@ANNIHILATOR="+string(ANN));
print("@@HILBERT_NUMERATOR");
hilb(Q,1);
print("@@COMPLETE=1");
quit;
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "? error occurred" not in completed.stdout
    lines = completed.stdout.splitlines()

    def marked_integer(name: str) -> int:
        prefix = f"@@{name}="
        values = [
            int(line[len(prefix) :])
            for line in lines
            if line.startswith(prefix)
        ]
        assert len(values) == 1, (name, completed.stdout)
        return values[0]

    annihilator_prefix = "@@ANNIHILATOR="
    annihilators = [
        line[len(annihilator_prefix) :]
        for line in lines
        if line.startswith(annihilator_prefix)
    ]
    assert len(annihilators) == 1
    marker = lines.index("@@HILBERT_NUMERATOR")
    numerator = atlas.parse_numerator(lines[marker + 1])
    assert "@@COMPLETE=1" in lines
    return {
        "kernel_generator_count": marked_integer("KERNEL_GENERATORS"),
        "generation_difference_size": marked_integer(
            "GENERATION_DIFFERENCE"
        ),
        "annihilator_singular_syntax": annihilators[0],
        "hilbert_numerator_over_one_minus_t_cubed": numerator,
    }


def singular_intrinsic_curvature_certificate() -> dict[str, Any]:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    program = """
LIB "primdec.lib";
ring curvature=0,(a,b,c,d),dp;
ideal I=
  (a-3*c)^2-9*d^2,
  a*b+a*d-3*b*c,
  a^2-9*b^2;
I=std(I);
ideal Pplus=a-3*b,a-3*c+3*d;
ideal Pminus=a+3*b,a-3*c-3*d;
Pplus=std(Pplus);
Pminus=std(Pminus);
ideal REDUCED=std(intersect(Pplus,Pminus));
ideal RADICAL=std(radical(I));
ideal RADICAL_MINUS_PLANES=simplify(reduce(RADICAL,REDUCED),2);
ideal PLANES_MINUS_RADICAL=simplify(reduce(REDUCED,RADICAL),2);
poly embedded= a*c-3*c^2+3*b*d+3*d^2;
poly embedded_remainder=reduce(embedded,I);
ideal WITH_EMBEDDED=I,embedded;
WITH_EMBEDDED=std(WITH_EMBEDDED);
ideal RADICAL_MINUS_GENERATED=simplify(
  reduce(RADICAL,WITH_EMBEDDED),2
);
ideal GENERATED_MINUS_RADICAL=simplify(
  reduce(WITH_EMBEDDED,RADICAL),2
);
ideal maximal_action=
  reduce(a*embedded,I),
  reduce(b*embedded,I),
  reduce(c*embedded,I),
  reduce(d*embedded,I);
maximal_action=simplify(maximal_action,2);
print("@@RADICAL_MINUS_PLANES="+string(size(RADICAL_MINUS_PLANES)));
print("@@PLANES_MINUS_RADICAL="+string(size(PLANES_MINUS_RADICAL)));
print("@@EMBEDDED_REMAINDER="+string(embedded_remainder));
print("@@RADICAL_MINUS_GENERATED="+string(size(RADICAL_MINUS_GENERATED)));
print("@@GENERATED_MINUS_RADICAL="+string(size(GENERATED_MINUS_RADICAL)));
print("@@MAXIMAL_ACTION="+string(size(maximal_action)));
print("@@COMPLETE=1");
quit;
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "? error occurred" not in completed.stdout
    lines = completed.stdout.splitlines()

    def marked(name: str) -> str:
        prefix = f"@@{name}="
        values = [
            line[len(prefix) :]
            for line in lines
            if line.startswith(prefix)
        ]
        assert len(values) == 1, (name, completed.stdout)
        return values[0]

    result = {
        "radical_minus_two_planes_size": int(
            marked("RADICAL_MINUS_PLANES")
        ),
        "two_planes_minus_radical_size": int(
            marked("PLANES_MINUS_RADICAL")
        ),
        "embedded_generator_remainder": marked(
            "EMBEDDED_REMAINDER"
        ),
        "radical_minus_ideal_plus_generator_size": int(
            marked("RADICAL_MINUS_GENERATED")
        ),
        "ideal_plus_generator_minus_radical_size": int(
            marked("GENERATED_MINUS_RADICAL")
        ),
        "maximal_ideal_action_size": int(marked("MAXIMAL_ACTION")),
    }
    assert "@@COMPLETE=1" in lines
    assert result == {
        "radical_minus_two_planes_size": 0,
        "two_planes_minus_radical_size": 0,
        "embedded_generator_remainder": "ac-3c2+3bd+3d2",
        "radical_minus_ideal_plus_generator_size": 0,
        "ideal_plus_generator_minus_radical_size": 0,
        "maximal_ideal_action_size": 0,
    }
    return result


def main() -> None:
    cubic.FACTOR_SINGULAR_EXPRESSIONS = False
    x, y, z = cubic.BASE_VARIABLES
    compatibility = smooth.compatibility_matrix()
    nodal_tensor = atlas.symbol_tensor(cubic.CUBIC_STRATA["nodal"])
    gauge = smooth.gauge_matrix(nodal_tensor)
    generator = atlas.symbol_tensor(cubic.Z**3)

    assert (
        compatibility * nodal_tensor
    ).applyfunc(sp.expand) == sp.zeros(6, 1)
    assert (
        compatibility * gauge
    ).applyfunc(sp.expand) == sp.zeros(6, 9)
    assert (
        compatibility * generator
    ).applyfunc(sp.expand) == sp.zeros(6, 1)
    smooth.verify_dual_number_gauge_action(nodal_tensor, gauge)

    module_certificate = singular_module_certificate(
        compatibility, gauge, generator
    )
    assert module_certificate == {
        "kernel_generator_count": 10,
        "generation_difference_size": 0,
        "annihilator_singular_syntax": "x",
        "hilbert_numerator_over_one_minus_t_cubed": [1, -1],
    }

    action_degree_four = sp.Matrix.hstack(
        *[
            coefficient_column(gauge[:, column] * variable, 4)
            for column in range(gauge.cols)
            for variable in cubic.BASE_VARIABLES
        ]
    )
    quotient_basis = sp.Matrix.hstack(
        coefficient_column(y * generator, 4),
        coefficient_column(z * generator, 4),
    )
    slice_matrix = action_degree_four.row_join(quotient_basis)
    assert action_degree_four.rank() == 22
    assert slice_matrix.rank() == 24

    parameters, universal = frontier.universal_tensor()
    directions = sp.Matrix.hstack(
        *[
            sp.Matrix(
                [
                    sp.expand(universal[triple]).coeff(parameter)
                    for triple in smooth.TRIPLES
                ]
            )
            for parameter in parameters
        ]
    )
    direction_coefficients = sp.Matrix.hstack(
        *[
            coefficient_column(directions[:, column], 4)
            for column in range(directions.cols)
        ]
    )
    assert direction_coefficients.rank() == 24

    decomposition = rational_solution(
        slice_matrix, direction_coefficients
    )
    gauge_lift = polynomial_lift(decomposition[:27, :])
    slice_projection = decomposition[27:29, :]
    expected_projection = sp.zeros(2, 24)
    expected_projection[0, 1] = 1
    expected_projection[1, 0] = 1
    assert slice_projection == expected_projection
    slice_terms = sp.Matrix.hstack(
        *[
            slice_projection[0, column] * y * generator
            + slice_projection[1, column] * z * generator
            for column in range(24)
        ]
    )
    assert (
        gauge * gauge_lift + slice_terms - directions
    ).applyfunc(sp.expand) == sp.zeros(10, 24)

    x_generator_solution = rational_solution(
        action_degree_four,
        coefficient_column(x * generator, 4),
    )
    x_generator_lift = polynomial_lift(x_generator_solution)
    assert (
        gauge * x_generator_lift - x * generator
    ).applyfunc(sp.expand) == sp.zeros(10, 1)

    (
        curvature,
        curvature_projection,
        curvature_coefficient_rows,
        nonzero_cross_terms,
    ) = (
        degree_five_curvature(
            nodal_tensor,
            gauge,
            generator,
            directions,
            gauge_lift,
            parameters,
        )
    )
    (
        u1,
        u2,
        u3,
        u4,
        u5,
        u6,
        u7,
        u8,
        u9,
        u10,
        u11,
        u12,
        _u13,
        _u14,
        _u15,
        u16,
        u17,
        u18,
        u19,
        u20,
        _u21,
        _u22,
        _u23,
        _u24,
    ) = parameters
    expected_curvature = sp.Matrix(
        (
            -(
                24 * u1 * u11
                + 4 * u1 * u17
                + 4 * u1 * u19
                + 6 * u1 * u3
                - 6 * u1 * u8
                + 6 * u1 * u9
                + 36 * u10 * u2
                - 54 * u12 * u2
                - 81 * u2 * u6
                - 36 * u2 * u7
                + 3 * u3**2
                - 18 * u3 * u5
                + 27 * u5**2
                - 27 * u6**2
            )
            / 12,
            (
                -60 * u1 * u10
                + 24 * u1 * u12
                - 4 * u1 * u16
                - 12 * u1 * u18
                - 4 * u1 * u20
                + 36 * u1 * u6
                + 42 * u1 * u7
                - 42 * u11 * u2
                - 4 * u17 * u2
                - 4 * u19 * u2
                + 21 * u2 * u3
                + 42 * u2 * u8
                - 42 * u2 * u9
                + 18 * u3 * u4
                + 18 * u3 * u6
                - 54 * u4 * u5
            )
            / 12,
            (
                -36 * u1 * u11
                + 18 * u1 * u3
                + 36 * u1 * u8
                - 54 * u1 * u9
                - 6 * u10 * u2
                - 12 * u12 * u2
                - 4 * u16 * u2
                - 12 * u18 * u2
                - 4 * u2 * u20
                - 18 * u2 * u6
                + 6 * u2 * u7
                + 3 * u3**2
                - 27 * u4**2
            )
            / 12,
        )
    )
    assert (
        curvature - expected_curvature
    ).applyfunc(sp.expand) == sp.zeros(3, 1)
    assert nonzero_cross_terms == 30
    coordinate_curvature = curvature.subs(
        {parameter: 0 for parameter in parameters[2:]}
    )
    assert coordinate_curvature == sp.zeros(3, 1)
    dense_parameters = sp.symbols("dense_plus dense_minus")
    dense_substitution = {
        parameter: dense_parameters[0]
        + (-1) ** index * dense_parameters[1]
        for index, parameter in enumerate(parameters)
    }
    dense_curvature = curvature.subs(dense_substitution).applyfunc(
        sp.factor
    )
    expected_dense_curvature = sp.Matrix(
        (
            7
            * (
                4 * dense_parameters[0] ** 2
                - 13 * dense_parameters[0] * dense_parameters[1]
                + dense_parameters[1] ** 2
            )
            / 3,
            -(
                25 * dense_parameters[0] ** 2
                - 193 * dense_parameters[1] ** 2
            )
            / 12,
            -(
                55 * dense_parameters[0] ** 2
                - 14 * dense_parameters[0] * dense_parameters[1]
                + 97 * dense_parameters[1] ** 2
            )
            / 6,
        )
    )
    assert (
        dense_curvature - expected_dense_curvature
    ).applyfunc(sp.expand) == sp.zeros(3, 1)

    def project_degree_five(vector: sp.Matrix) -> sp.Matrix:
        coefficients = coefficient_column(vector, 5)
        return curvature_projection * coefficients[
            list(curvature_coefficient_rows), :
        ]

    lift_kernel_vectors = action_degree_four.nullspace()
    assert len(lift_kernel_vectors) == 5
    lift_kernel = [
        sp.Matrix(
            3,
            3,
            list(polynomial_lift(vector)[:, 0]),
        )
        for vector in lift_kernel_vectors
    ]
    for stabilizer in lift_kernel:
        assert (
            gauge * sp.Matrix(list(stabilizer))
        ).applyfunc(sp.expand) == sp.zeros(10, 1)

    # Stabilizer--stabilizer terms are degree-five gauge, as are all
    # stabilizer changes paired with one of the 22 removable directions.
    for first in range(5):
        for second in range(first, 5):
            assert project_degree_five(
                second_action_polarization(
                    nodal_tensor,
                    lift_kernel[first],
                    lift_kernel[second],
                )
            ) == sp.zeros(3, 1)
    direction_gauges = [
        smooth.gauge_matrix(directions[:, column])
        for column in range(24)
    ]
    lift_derivations = [
        sp.Matrix(3, 3, list(gauge_lift[:, column]))
        for column in range(24)
    ]
    for direction in range(2, 24):
        for stabilizer in lift_kernel:
            mixed = (
                second_action_polarization(
                    nodal_tensor,
                    lift_derivations[direction],
                    stabilizer,
                )
                - direction_gauges[direction]
                * sp.Matrix(list(stabilizer))
            )
            assert project_degree_five(mixed) == sp.zeros(3, 1)

    # Changes of lift act on the two slice--gauge cross coefficients by
    # this single rank-four map W -> Q_5 direct_sum Q_5.
    lift_kernel_action = sp.Matrix.vstack(
        *[
            sp.Matrix.hstack(
                *[
                    project_degree_five(
                        second_action_polarization(
                            nodal_tensor,
                            lift_derivations[slice_direction],
                            stabilizer,
                        )
                        - direction_gauges[slice_direction]
                        * sp.Matrix(list(stabilizer))
                    )
                    for stabilizer in lift_kernel
                ]
            )
            for slice_direction in (0, 1)
        ]
    )
    expected_lift_kernel_action = sp.Matrix(
        (
            (sp.Rational(1, 3), 0, 0, 0, sp.Rational(1, 2)),
            (0, sp.Rational(-7, 3), 0, 4, 0),
            (2, 0, 0, 0, sp.Rational(7, 2)),
            (0, -2, 0, sp.Rational(7, 2), 0),
            (sp.Rational(7, 3), 0, 0, 0, 4),
            (0, sp.Rational(-1, 3), 0, sp.Rational(1, 2), 0),
        )
    )
    assert lift_kernel_action == expected_lift_kernel_action
    assert lift_kernel_action.rank() == 4
    lift_invariant_functionals = sp.Matrix(
        (
            (-1, 0, -1, 0, 1, 0),
            (0, -1, 0, 1, 0, 1),
        )
    )
    assert (
        lift_invariant_functionals * lift_kernel_action
    ) == sp.zeros(2, 5)
    assert lift_invariant_functionals.rank() == 2

    curvature_polynomials = [
        sp.Poly(component, *parameters) for component in curvature
    ]
    slice_gauge_cross = sp.zeros(6, 22)
    for offset, direction in enumerate(range(2, 24)):
        for component in range(3):
            slice_gauge_cross[component, offset] = (
                curvature_polynomials[component].coeff_monomial(
                    parameters[0] * parameters[direction]
                )
            )
            slice_gauge_cross[3 + component, offset] = (
                curvature_polynomials[component].coeff_monomial(
                    parameters[1] * parameters[direction]
                )
            )
    intrinsic_cross_matrix = (
        lift_invariant_functionals * slice_gauge_cross
    )
    intrinsic_cross_forms = sp.Matrix(
        [
            sum(
                intrinsic_cross_matrix[row, offset]
                * parameters[offset + 2]
                for offset in range(22)
            )
            for row in range(2)
        ]
    ).applyfunc(sp.factor)
    expected_intrinsic_cross_forms = sp.Matrix(
        (
            sp.Rational(3, 4)
            * (u3 + 2 * u9 + 2 * u11),
            sp.Rational(3, 4)
            * (3 * u6 + 2 * u10 + 2 * u12),
        )
    )
    assert (
        intrinsic_cross_forms - expected_intrinsic_cross_forms
    ).applyfunc(sp.expand) == sp.zeros(2, 1)

    pure_gauge_curvature = curvature.subs(
        {parameters[0]: 0, parameters[1]: 0}
    ).applyfunc(sp.factor)
    expected_pure_gauge_curvature = sp.Matrix(
        (
            -(
                (u3 - 3 * u5) ** 2 - 9 * u6**2
            )
            / 4,
            sp.Rational(3, 2)
            * (u3 * u4 + u3 * u6 - 3 * u4 * u5),
            (u3**2 - 9 * u4**2) / 4,
        )
    )
    assert (
        pure_gauge_curvature - expected_pure_gauge_curvature
    ).applyfunc(sp.expand) == sp.zeros(3, 1)
    intrinsic_curvature_scheme = (
        singular_intrinsic_curvature_certificate()
    )

    # Continue the two reduced pure-curvature planes through collision
    # degree six, relative to the declared row-reduced quartic lift.
    branch_p, branch_q = sp.symbols("branch_p branch_q")
    degree_five_monomials = cubic.homogeneous_monomials(2)
    action_degree_five = sp.Matrix.hstack(
        *[
            coefficient_column(
                gauge[:, column] * monomial, 5
            )
            for column in range(gauge.cols)
            for monomial in degree_five_monomials
        ]
    )
    degree_five_lift_kernel = action_degree_five.nullspace()
    assert action_degree_five.rank() == 39
    assert len(degree_five_lift_kernel) == 15

    degree_six_monomials = cubic.homogeneous_monomials(3)
    action_degree_six = sp.Matrix.hstack(
        *[
            coefficient_column(
                gauge[:, column] * monomial, 6
            )
            for column in range(gauge.cols)
            for monomial in degree_six_monomials
        ]
    )
    quotient_basis_degree_six = sp.Matrix.hstack(
        *[
            coefficient_column(monomial * generator, 6)
            for monomial in (
                y**3,
                y**2 * z,
                y * z**2,
                z**3,
            )
        ]
    )
    action_degree_six_pivots = action_degree_six.rref()[1]
    reduced_basis_degree_six = sp.Matrix.hstack(
        *[
            action_degree_six[:, column]
            for column in action_degree_six_pivots
        ],
        quotient_basis_degree_six,
    )
    degree_six_coefficient_rows = (
        reduced_basis_degree_six.T.rref()[1]
    )
    assert action_degree_six.rank() == 60
    assert len(action_degree_six_pivots) == 60
    assert reduced_basis_degree_six.rank() == 64
    assert len(degree_six_coefficient_rows) == 64
    degree_six_projection = (
        reduced_basis_degree_six[
            list(degree_six_coefficient_rows), :
        ].inv()[-4:, :]
    )

    def project_degree_six(vector: sp.Matrix) -> sp.Matrix:
        coefficients = coefficient_column(vector, 6)
        return (
            degree_six_projection
            * coefficients[list(degree_six_coefficient_rows), :]
        ).applyfunc(sp.factor)

    branch_records: dict[str, dict[str, Any]] = {}
    for branch_name, branch_sign in (("plus", 1), ("minus", -1)):
        branch_substitution = {
            parameter: sp.Integer(0) for parameter in parameters
        }
        branch_substitution.update(
            {
                parameters[2]: 3 * branch_sign * branch_p,
                parameters[3]: branch_p,
                parameters[4]: branch_sign
                * (branch_p + branch_q),
                parameters[5]: branch_q,
            }
        )
        branch_parameter_vector = sp.Matrix(
            [
                branch_substitution[parameter]
                for parameter in parameters
            ]
        )
        branch_perturbation = directions * branch_parameter_vector
        branch_derivation = sp.Matrix(
            3,
            3,
            list(gauge_lift * branch_parameter_vector),
        )
        expansion = finite_inverse_gauge_coefficients(
            nodal_tensor,
            branch_perturbation,
            branch_derivation,
            3,
        )
        assert expansion[0] == nodal_tensor
        assert expansion[1] == sp.zeros(10, 1)

        degree_five_coefficients = coefficient_column(
            expansion[2], 5
        )
        assert (
            curvature_projection
            * degree_five_coefficients[
                list(curvature_coefficient_rows), :
            ]
        ) == sp.zeros(3, 1)
        correction_solution = rational_solution_preserving(
            action_degree_five,
            degree_five_coefficients,
            (branch_p, branch_q),
        )
        degree_five_correction = homogeneous_polynomial_lift(
            correction_solution, 2
        )
        assert (
            gauge * sp.Matrix(list(degree_five_correction))
            - expansion[2]
        ).applyfunc(sp.expand) == sp.zeros(10, 1)

        # Every alternative quadratic correction differs by a vector in
        # this 15-dimensional kernel.  Its action on Q_6 is zero.
        for kernel_vector in degree_five_lift_kernel:
            kernel_derivation = homogeneous_polynomial_lift(
                kernel_vector, 2
            )
            assert project_degree_six(
                smooth.gauge_matrix(branch_perturbation)
                * sp.Matrix(list(kernel_derivation))
            ) == sp.zeros(4, 1)

        degree_six_remainder = (
            expansion[3]
            - smooth.gauge_matrix(branch_perturbation)
            * sp.Matrix(list(degree_five_correction))
        ).applyfunc(sp.expand)
        degree_six_obstruction = project_degree_six(
            degree_six_remainder
        )
        expected_degree_six_obstruction = (
            sp.Rational(27, 8)
            * sp.Matrix(
                (
                    branch_q**3,
                    3
                    * branch_sign
                    * branch_p
                    * branch_q**2,
                    3 * branch_p**2 * branch_q,
                    branch_sign * branch_p**3,
                )
            )
        )
        assert (
            degree_six_obstruction
            - expected_degree_six_obstruction
        ).applyfunc(sp.expand) == sp.zeros(4, 1)
        branch_records[branch_name] = {
            "parameterization": {
                "u3": sp.sstr(3 * branch_sign * branch_p),
                "u4": sp.sstr(branch_p),
                "u5": sp.sstr(
                    branch_sign * (branch_p + branch_q)
                ),
                "u6": sp.sstr(branch_q),
            },
            "degree_five_correction": matrix_record(
                degree_five_correction
            ),
            "degree_six_obstruction": matrix_record(
                degree_six_obstruction
            ),
            "compact_form": (
                "27/8*(branch_q*y"
                + ("+" if branch_sign == 1 else "-")
                + "branch_p*z)^3*eta"
            ),
        }

    dense_parameter_matrix = sp.Matrix(
        24,
        2,
        lambda row, column: (
            1 if column == 0 else (-1) ** row
        ),
    )
    dense_slice_coordinates = (
        slice_projection * dense_parameter_matrix
    )
    assert dense_slice_coordinates == sp.Matrix(((1, -1), (1, 1)))
    assert dense_slice_coordinates.det() == 2

    basis = cubic.quartic_kernel_basis_tensors()
    coordinate_cotangent_result = cubic.run_singular_subspace(
        cubic.CUBIC_STRATA["nodal"],
        (basis[0], basis[1]),
        timeout=600,
    )
    dense_cotangent_result = cubic.run_singular_subspace(
        cubic.CUBIC_STRATA["nodal"],
        dense.dense_directions(),
        timeout=600,
    )
    assert coordinate_cotangent_result == (0, 6, 0, 0, 3)
    assert dense_cotangent_result == (0, 6, 0, 0, 3)

    exact_data = {
        "compatibility_matrix": matrix_record(compatibility),
        "nodal_gauge_matrix": matrix_record(gauge),
        "cyclic_generator_Z_cubed": matrix_record(generator),
        "x_generator_gauge_lift": matrix_record(x_generator_lift),
        "quartic_gauge_lift": matrix_record(gauge_lift),
        "slice_projection_y_eta_z_eta": matrix_record(
            slice_projection
        ),
        "dense_slice_coordinates": matrix_record(
            dense_slice_coordinates
        ),
        "degree_five_curvature_basis_y2_eta_yz_eta_z2_eta": (
            matrix_record(curvature)
        ),
        "degree_five_dense_slice_curvature": matrix_record(
            dense_curvature
        ),
        "quartic_lift_kernel_basis": matrix_record(
            sp.Matrix.hstack(
                *[
                    sp.Matrix(list(stabilizer))
                    for stabilizer in lift_kernel
                ]
            )
        ),
        "quartic_lift_kernel_action_on_slice_cross_terms": (
            matrix_record(lift_kernel_action)
        ),
        "lift_invariant_functionals": matrix_record(
            lift_invariant_functionals
        ),
        "intrinsic_slice_gauge_cross_forms": matrix_record(
            intrinsic_cross_forms
        ),
        "intrinsic_pure_gauge_curvature": matrix_record(
            pure_gauge_curvature
        ),
        "intrinsic_curvature_reduced_planes": [
            ["u3-3*u4", "u3-3*u5+3*u6"],
            ["u3+3*u4", "u3-3*u5-3*u6"],
        ],
        "intrinsic_curvature_embedded_generator": (
            "u3*u5-3*u5^2+3*u4*u6+3*u6^2"
        ),
        "degree_six_reduced_plane_continuation": branch_records,
    }
    exact_sha256 = hashlib.sha256(
        json.dumps(
            exact_data,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()
    artifact = {
        "schema": "nodal-cubic-formal-slice.v3",
        "mathematical_status": (
            "exact first-stage slice, intrinsic degree-five curvature, "
            "and relative degree-six continuation theorem"
        ),
        "basis_conventions": {
            "collision_variables": ["x", "y", "z"],
            "tensor_component_order": [
                list(triple) for triple in smooth.TRIPLES
            ],
            "quartic_parameter_order": [
                str(parameter) for parameter in parameters
            ],
            "slice_basis": ["y*eta", "z*eta"],
            "eta_source_cubic": "Z^3",
        },
        "module_certificate": module_certificate,
        "quartic_dimensions": {
            "compatible_space": 24,
            "gauge_image": 22,
            "slice": 2,
        },
        "coordinate_slice_basis_indices_zero_based": [0, 1],
        "coordinate_slice_cotangent_result": list(
            coordinate_cotangent_result
        ),
        "dense_slice_cotangent_result": list(
            dense_cotangent_result
        ),
        "degree_five_curvature": {
            "gauge_lift_convention": (
                "free variables in the exact rational quartic lift "
                "are set to zero"
            ),
            "quotient_basis": ["y^2*eta", "y*z*eta", "z^2*eta"],
            "quadratic_term_counts": [
                len(sp.Poly(entry, *parameters).terms())
                for entry in curvature
            ],
            "nonzero_cross_parameter_pairs": nonzero_cross_terms,
            "coordinate_slice_is_zero": True,
            "dense_slice_is_nonzero": True,
            "quartic_lift_kernel_dimension": 5,
            "lift_kernel_action_rank": 4,
            "lift_kernel_action_cokernel_dimension": 2,
            "intrinsic_cross_forms": [
                sp.sstr(entry) for entry in intrinsic_cross_forms
            ],
            "intrinsic_pure_gauge_curvature": [
                sp.sstr(entry) for entry in pure_gauge_curvature
            ],
            "intrinsic_curvature_scheme": (
                intrinsic_curvature_scheme
            ),
        },
        "degree_six_continuation": {
            "scope": (
                "the two reduced pure-curvature planes for the stored "
                "row-reduced quartic gauge lift"
            ),
            "degree_five_gauge_action_rank": 39,
            "degree_five_correction_kernel_dimension": 15,
            "degree_five_correction_kernel_action_on_Q6": "zero",
            "degree_six_gauge_action_rank": 60,
            "degree_six_quotient_dimension": 4,
            "branches": branch_records,
            "common_zero_locus_on_each_reduced_plane": (
                "branch_p=branch_q=0"
            ),
        },
        "exact_data": exact_data,
        "exact_data_sha256": exact_sha256,
        "proved": [
            "ker(C)/im(G_nodal) is Q[y,z](-3), generated by eta(Z^3)",
            (
                "the universal quartic space is the direct sum of the "
                "22-dimensional gauge image and the first-two-direction "
                "coordinate slice"
            ),
            (
                "the full-support sum/alternating-sum plane is transverse "
                "to the gauge image with determinant two"
            ),
            (
                "both transverse nodal slices have saturated cotangent "
                "presentation and the central length-six Ext block"
            ),
            (
                "for the stored deterministic quartic gauge lift, the "
                "complete degree-five normal remainder is the recorded "
                "three-component quadratic curvature"
            ),
            (
                "the five-dimensional quartic lift kernel acts with "
                "rank four on slice--gauge cross curvature, has the "
                "recorded two-dimensional cokernel, and acts trivially "
                "on pure gauge curvature"
            ),
            (
                "the displayed two cross forms and three pure-gauge "
                "quadrics are invariant under every quartic lift choice"
            ),
            (
                "the reduced pure-curvature zero scheme is the union of "
                "the two recorded planes, and the unreduced ideal has "
                "one embedded degree-two socle generator"
            ),
            (
                "for the stored row-reduced quartic lift, both reduced "
                "pure-curvature planes admit exact degree-five gauge "
                "corrections and have degree-six classes "
                "27/8*(q*y+p*z)^3*eta and "
                "27/8*(q*y-p*z)^3*eta"
            ),
            (
                "the 15-dimensional ambiguity in the degree-five "
                "correction acts trivially on the degree-six quotient"
            ),
        ],
        "not_proved": [
            (
                "independence of the degree-six class from the earlier "
                "five-dimensional quartic gauge-lift ambiguity"
            ),
            (
                "degree-six continuation over the embedded quadratic "
                "socle or the full slice--gauge curvature locus"
            ),
            "an all-order formal slice for the universal nodal family",
            "universal 24-parameter nodal cotangent saturation",
            "normality or Keller-open compatibility",
        ],
        "reproduce": (
            ".venv/bin/python "
            "scripts/verify_nodal_cubic_formal_slice.py"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("PASS: nodal gauge cokernel is Q[y,z](-3)")
    print("PASS: universal quartic space splits as gauge 22 plus slice 2")
    print("PASS: exact degree-five curvature has 30 nonzero cross pairs")
    print("PASS: lift-kernel quotient leaves two cross and three pure invariants")
    print("PASS: intrinsic zero scheme is two planes plus one embedded socle")
    print("PASS: reduced curvature planes have exact degree-six Veronese classes")
    print("PASS: degree-five correction ambiguity acts trivially on Q6")
    print("PASS: coordinate and dense transverse slices are cotangent-saturated")
    print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
