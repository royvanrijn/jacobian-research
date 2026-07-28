#!/usr/bin/env python3
"""Exact stratum audit for the reduced cubic Koszul defect.

For each projective ternary-cubic orbit representative, construct the
homogeneous generalized triple-cover algebra on

    M = coker(A --(z,-y,x)^T--> A^3),  A=Q[x,y,z].

The multiplication is recovered invariantly from its cubic tensor.  The
script then asks Singular for:

1. the (x,y,z)-saturation of an A-presentation of Omega_{B/A};
2. Ext_A^2(T,A), where T=B/Ann_B(Omega_{B/A}).

This is a leading-model computation.  It does not assert normality of the
homogeneous algebra, existence of a Keller open, or invariance under
higher-order tensor perturbations.
"""

from __future__ import annotations

import itertools
import math
import shutil
import subprocess
from functools import cache

import sympy as sp


x, y, z = sp.symbols("x y z")
X, Y, Z = sp.symbols("X Y Z")
BASE_VARIABLES = (x, y, z)
RELATION = sp.Matrix((z, -y, x))
STANDARD_BASIS = tuple(sp.eye(3).col(index) for index in range(3))
POLARIZATION_VARIABLES = sp.symbols("polarization_0:3")


def homogeneous_monomials(degree: int) -> list[sp.Expr]:
    return [
        x**x_degree * y**y_degree * z ** (degree - x_degree - y_degree)
        for x_degree in range(degree + 1)
        for y_degree in range(degree - x_degree + 1)
    ]


def polarized_value(
    cubic: sp.Expr,
    first: sp.Matrix,
    second: sp.Matrix,
    third: sp.Matrix,
) -> sp.Expr:
    substitution = {
        variable: (
            POLARIZATION_VARIABLES[0] * first[index]
            + POLARIZATION_VARIABLES[1] * second[index]
            + POLARIZATION_VARIABLES[2] * third[index]
        )
        for index, variable in enumerate((X, Y, Z))
    }
    polynomial = sp.Poly(
        sp.expand(cubic.subs(substitution)),
        *POLARIZATION_VARIABLES,
    )
    return sp.expand(
        polynomial.coeff_monomial(sp.prod(POLARIZATION_VARIABLES)) / 6
    )


def solve_cross_product(target: sp.Matrix) -> sp.Matrix:
    """Find a polynomial vector m with RELATION cross m = target."""

    if target == sp.zeros(3, 1):
        return sp.zeros(3, 1)
    target_degree = max(
        sp.Poly(entry, *BASE_VARIABLES).total_degree()
        for entry in target
        if entry != 0
    )
    candidate_monomials = [
        monomial
        for degree in range(2, target_degree)
        for monomial in homogeneous_monomials(degree)
    ]
    output_monomials = [
        monomial
        for degree in range(3, target_degree + 1)
        for monomial in homogeneous_monomials(degree)
    ]
    coefficient_count = 3 * len(candidate_monomials)
    coefficients = sp.symbols(
        f"cross_coefficient_0:{coefficient_count}"
    )
    candidate = sp.Matrix(
        [
            sum(
                coefficients[len(candidate_monomials) * row + column]
                * monomial
                for column, monomial in enumerate(candidate_monomials)
            )
            for row in range(3)
        ]
    )
    equations = []
    for entry in RELATION.cross(candidate) - target:
        polynomial = sp.Poly(entry, *BASE_VARIABLES)
        equations.extend(
            polynomial.coeff_monomial(monomial)
            for monomial in output_monomials
        )
    solution = next(iter(sp.linsolve(equations, coefficients)))
    target_parameters = set().union(
        *(entry.free_symbols for entry in target)
    ) - set(BASE_VARIABLES)
    free_parameters = (
        set().union(*(entry.free_symbols for entry in solution))
        - set(BASE_VARIABLES)
        - target_parameters
    )
    result = candidate.subs(dict(zip(coefficients, solution)))
    result = result.subs(
        {parameter: 0 for parameter in free_parameters}
    ).applyfunc(sp.expand)
    assert (
        RELATION.cross(result) - target
    ).applyfunc(sp.expand) == sp.zeros(3, 1)
    return result


def multiplication_table(
    cubic: sp.Expr,
    higher_tensor: dict[tuple[int, int, int], sp.Expr] | None = None,
) -> tuple[dict[tuple[int, int], sp.Matrix], dict[tuple[int, int], sp.Expr]]:
    """Recover trace-free and scalar multiplication from the cubic tensor."""

    trace_free_products: dict[tuple[int, int], sp.Matrix] = {}
    tensor_values = {
        triple: polarized_value(
            cubic,
            *(RELATION.cross(STANDARD_BASIS[index]) for index in triple),
        )
        + (higher_tensor or {}).get(triple, 0)
        for triple in itertools.combinations_with_replacement(range(3), 3)
    }
    for first in range(3):
        for second in range(first, 3):
            target = sp.Matrix(
                [
                    3 * tensor_values[tuple(sorted((first, second, third)))]
                    for third in range(3)
                ]
            )
            assert sp.expand(RELATION.dot(target)) == 0
            trace_free_products[first, second] = solve_cross_product(target)

    multiplication_lifts = [
        sp.Matrix.hstack(
            *[
                trace_free_products[tuple(sorted((first, second)))]
                for second in range(3)
            ]
        )
        for first in range(3)
    ]

    relation_eigenvalues = []
    for first, lift in enumerate(multiplication_lifts):
        if lift == sp.zeros(3, 3):
            relation_eigenvalues.append(sp.Integer(0))
            continue
        lift_degree = max(
            sp.Poly(entry, *BASE_VARIABLES).total_degree()
            for entry in lift
            if entry != 0
        )
        eigenvalue_monomials = [
            monomial
            for degree in range(2, lift_degree + 1)
            for monomial in homogeneous_monomials(degree)
        ]
        coefficients = sp.symbols(
            f"relation_coefficient_{first}_0:{len(eigenvalue_monomials)}"
        )
        eigenvalue = sum(
            coefficient * monomial
            for coefficient, monomial in zip(
                coefficients, eigenvalue_monomials
            )
        )
        equation_monomials = [
            monomial
            for degree in range(3, lift_degree + 2)
            for monomial in homogeneous_monomials(degree)
        ]
        equations = []
        for entry in lift * RELATION - eigenvalue * RELATION:
            polynomial = sp.Poly(entry, *BASE_VARIABLES)
            equations.extend(
                polynomial.coeff_monomial(monomial)
                for monomial in equation_monomials
            )
        solution = next(iter(sp.linsolve(equations, coefficients)))
        eigenvalue = sp.expand(
            eigenvalue.subs(dict(zip(coefficients, solution)))
        )
        assert (
            lift * RELATION - eigenvalue * RELATION
        ).applyfunc(sp.expand) == sp.zeros(3, 1)
        relation_eigenvalues.append(eigenvalue)

    scalar_products = {
        (first, second): sp.expand(
            sp.trace(
                multiplication_lifts[first]
                * multiplication_lifts[second]
            )
            - relation_eigenvalues[first] * relation_eigenvalues[second]
        )
        for first in range(3)
        for second in range(first, 3)
    }

    # The trace formula for the scalar part gives an associative algebra.
    def multiply_generator_by_vector(
        generator: int, vector: sp.Matrix
    ) -> tuple[sp.Expr, sp.Matrix]:
        scalar = sum(
            vector[index]
            * scalar_products[tuple(sorted((generator, index)))]
            for index in range(3)
        )
        return (
            sp.expand(scalar),
            (multiplication_lifts[generator] * vector).applyfunc(sp.expand),
        )

    for first, second, third in itertools.product(range(3), repeat=3):
        second_third = trace_free_products[
            tuple(sorted((second, third)))
        ]
        first_second = trace_free_products[
            tuple(sorted((first, second)))
        ]
        left_scalar, left_vector = multiply_generator_by_vector(
            first, second_third
        )
        right_scalar, right_vector = multiply_generator_by_vector(
            third, first_second
        )
        left_vector += (
            scalar_products[tuple(sorted((second, third)))]
            * STANDARD_BASIS[first]
        )
        right_vector += (
            scalar_products[tuple(sorted((first, second)))]
            * STANDARD_BASIS[third]
        )
        assert sp.expand(left_scalar - right_scalar) == 0
        difference = (left_vector - right_vector).applyfunc(sp.expand)
        assert RELATION.cross(difference).applyfunc(sp.expand) == sp.zeros(
            3, 1
        )

    return trace_free_products, scalar_products


def differential_relations(
    cubic: sp.Expr,
    higher_tensor: dict[tuple[int, int, int], sp.Expr] | None = None,
) -> list[list[sp.Expr]]:
    """An A-presentation of Omega_{B/A} on the 12 generators b_i de_j."""

    trace_free_products, scalar_products = multiplication_table(
        cubic, higher_tensor
    )

    def algebra_product(first: int, second: int) -> list[sp.Expr]:
        if first == 0:
            return [
                sp.Integer(index == second) for index in range(4)
            ]
        if second == 0:
            return [
                sp.Integer(index == first) for index in range(4)
            ]
        pair = tuple(sorted((first - 1, second - 1)))
        return [scalar_products[pair], *trace_free_products[pair]]

    def index(algebra_generator: int, differential_generator: int) -> int:
        return 4 * differential_generator + algebra_generator

    relations: list[list[sp.Expr]] = []

    # The coefficient-module relation in B, once for each differential.
    for differential in range(3):
        relation = [sp.Integer(0)] * 12
        relation[index(1, differential)] += z
        relation[index(2, differential)] -= y
        relation[index(3, differential)] += x
        relations.append(relation)

    # The differential of z e_1-y e_2+x e_3, multiplied by a B-generator.
    for algebra_generator in range(4):
        relation = [sp.Integer(0)] * 12
        relation[index(algebra_generator, 0)] += z
        relation[index(algebra_generator, 1)] -= y
        relation[index(algebra_generator, 2)] += x
        relations.append(relation)

    # Differentials of the six multiplication relations and their B-multiples.
    for first in range(3):
        for second in range(first, 3):
            pair = (first, second)
            for algebra_generator in range(4):
                relation = [sp.Integer(0)] * 12
                for coefficient_generator, coefficient in enumerate(
                    algebra_product(algebra_generator, first + 1)
                ):
                    relation[
                        index(coefficient_generator, second)
                    ] += coefficient
                for coefficient_generator, coefficient in enumerate(
                    algebra_product(algebra_generator, second + 1)
                ):
                    relation[
                        index(coefficient_generator, first)
                    ] += coefficient
                for differential in range(3):
                    relation[
                        index(algebra_generator, differential)
                    ] -= trace_free_products[pair][differential]
                relations.append(
                    [sp.expand(entry) for entry in relation]
                )

    assert len(relations) == 31
    return relations


def singular_polynomial(expression: sp.Expr) -> str:
    return sp.sstr(sp.factor(expression)).replace("**", "^")


def singular_vector(vector: list[sp.Expr]) -> str:
    return "[" + ",".join(map(singular_polynomial, vector)) + "]"


def singular_program(
    cubic: sp.Expr,
    higher_tensor: dict[tuple[int, int, int], sp.Expr] | None = None,
) -> str:
    relations = differential_relations(cubic, higher_tensor)

    # To compute Ann_B(Q), take the kernel of
    # B -> Q^3, b |-> (b de_1,b de_2,b de_3).  The first four syzygy
    # coordinates give its preimage in A^4, hence a presentation of T.
    annihilator_kernel_relations: list[list[sp.Expr]] = []
    for algebra_generator in range(4):
        relation = [sp.Integer(0)] * 36
        for differential in range(3):
            relation[
                12 * differential
                + 4 * differential
                + algebra_generator
            ] = 1
        annihilator_kernel_relations.append(relation)
    for block in range(3):
        for source_relation in relations:
            relation = [sp.Integer(0)] * 36
            relation[12 * block : 12 * (block + 1)] = source_relation
            annihilator_kernel_relations.append(relation)

    return f"""
LIB "primdec.lib";
LIB "homolog.lib";
ring coefficient_ring=0,(x,y,z),dp;
ideal collision=x,y,z;
module differential_presentation=
{",".join(map(singular_vector, relations))};
differential_presentation=std(differential_presentation);
list saturation_data=sat(differential_presentation,collision);
module saturated_presentation=saturation_data[1];
module saturation_quotient=simplify(
  reduce(saturated_presentation,differential_presentation),2
);
print("SATURATION_GENERATORS="+string(size(saturation_quotient)));

module annihilator_kernel_problem=
{",".join(map(singular_vector, annihilator_kernel_relations))};
module annihilator_syzygies=syz(annihilator_kernel_problem);
matrix support_presentation_matrix[4][ncols(annihilator_syzygies)];
int row_index,column_index;
for(column_index=1;
    column_index<=ncols(annihilator_syzygies);
    column_index++)
{{
  for(row_index=1;row_index<=4;row_index++)
  {{
    support_presentation_matrix[row_index,column_index]
      =annihilator_syzygies[column_index][row_index];
  }}
}}
module support_presentation=std(module(support_presentation_matrix));
print("SUPPORT_DIMENSION="+string(dim(support_presentation)));
module support_ext3=std(Ext_R(3,support_presentation));
print("EXT3_VECTOR_DIMENSION="+string(vdim(support_ext3)));
module support_ext2=Ext_R(2,support_presentation);
support_ext2=std(support_ext2);
print("EXT2_DIMENSION="+string(dim(support_ext2)));
print("EXT2_VECTOR_DIMENSION="+string(vdim(support_ext2)));
module maximal_times_free=
  [x,0,0],[y,0,0],[z,0,0],
  [0,x,0],[0,y,0],[0,z,0],
  [0,0,x],[0,0,y],[0,0,z];
module ext2_top_presentation=std(
  support_ext2+maximal_times_free
);
print("EXT2_TOP_DIMENSION="+string(vdim(ext2_top_presentation)));
module maximal_square_times_free=
  [x2,0,0],[xy,0,0],[xz,0,0],
  [y2,0,0],[yz,0,0],[z2,0,0],
  [0,x2,0],[0,xy,0],[0,xz,0],
  [0,y2,0],[0,yz,0],[0,z2,0],
  [0,0,x2],[0,0,xy],[0,0,xz],
  [0,0,y2],[0,0,yz],[0,0,z2];
module ext2_square_action=simplify(
  reduce(maximal_square_times_free,support_ext2),2
);
print("EXT2_SQUARE_ACTION_GENERATORS="+string(size(ext2_square_action)));
quit;
"""


def run_singular(
    cubic: sp.Expr,
    higher_tensor: dict[tuple[int, int, int], sp.Expr] | None = None,
) -> tuple[int, int, int, int, int, int, int]:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for this checker"
    result = subprocess.run(
        [singular, "-q"],
        input=singular_program(cubic, higher_tensor),
        text=True,
        capture_output=True,
        check=True,
    )
    values: dict[str, int] = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key in {
                "SATURATION_GENERATORS",
                "SUPPORT_DIMENSION",
                "EXT3_VECTOR_DIMENSION",
                "EXT2_DIMENSION",
                "EXT2_VECTOR_DIMENSION",
                "EXT2_TOP_DIMENSION",
                "EXT2_SQUARE_ACTION_GENERATORS",
            }:
                values[key] = int(value)
    assert set(values) == {
        "SATURATION_GENERATORS",
        "SUPPORT_DIMENSION",
        "EXT3_VECTOR_DIMENSION",
        "EXT2_DIMENSION",
        "EXT2_VECTOR_DIMENSION",
        "EXT2_TOP_DIMENSION",
        "EXT2_SQUARE_ACTION_GENERATORS",
    }, result.stdout + result.stderr
    return (
        values["SATURATION_GENERATORS"],
        values["SUPPORT_DIMENSION"],
        values["EXT3_VECTOR_DIMENSION"],
        values["EXT2_DIMENSION"],
        values["EXT2_VECTOR_DIMENSION"],
        values["EXT2_TOP_DIMENSION"],
        values["EXT2_SQUARE_ACTION_GENERATORS"],
    )


def singular_family_program(cubic: sp.Expr) -> str:
    """The exact family phi_h+t*psi_4 over Q[t,x,y,z]."""

    deformation_parameter = sp.symbols("t")
    family_tensor = {
        triple: deformation_parameter * value
        for triple, value in generic_quartic_tensor().items()
    }
    program = singular_program(cubic, family_tensor).replace(
        "ring coefficient_ring=0,(x,y,z),dp;",
        "ring coefficient_ring=0,(t,x,y,z),dp;",
    )
    return program.replace(
        'print("EXT2_VECTOR_DIMENSION="'
        "+string(vdim(support_ext2)));",
        'print("EXT2_VECTOR_DIMENSION="'
        "+string(vdim(support_ext2)));"
        'print("EXT2_MULTIPLICITY="+string(mult(support_ext2)));'
        "list parameter_saturation=sat(support_ext2,ideal(t));"
        "module parameter_torsion=simplify("
        "reduce(parameter_saturation[1],support_ext2),2);"
        'print("PARAMETER_TORSION_GENERATORS="'
        "+string(size(parameter_torsion)));"
        "ideal ext2_fitting=fitting(support_ext2,0);"
        "ideal ext2_support=std(radical(ext2_fitting));"
        "ideal collision_axis=std(ideal(x,y,z));"
        "ideal first_support_difference=simplify("
        "reduce(ext2_support,collision_axis),2);"
        "ideal second_support_difference=simplify("
        "reduce(collision_axis,ext2_support),2);"
        'print("COLLISION_AXIS_DIFFERENCE="'
        "+string(size(first_support_difference)"
        "+size(second_support_difference)));"
        "module central_ext2_presentation=std("
        "subst(support_ext2,t,0));"
        "module first_presentation_difference=simplify("
        "reduce(support_ext2,central_ext2_presentation),2);"
        "module second_presentation_difference=simplify("
        "reduce(central_ext2_presentation,support_ext2),2);"
        'print("PRESENTATION_PARAMETER_DIFFERENCE="'
        "+string(size(first_presentation_difference)"
        "+size(second_presentation_difference)));",
    )


def run_singular_family(cubic: sp.Expr) -> tuple[int, int, int, int, int]:
    """Return uniform cotangent, support, multiplicity, and t-torsion data."""

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for this checker"
    result = subprocess.run(
        [singular, "-q"],
        input=singular_family_program(cubic),
        text=True,
        capture_output=True,
        check=True,
    )
    wanted = {
        "SATURATION_GENERATORS",
        "EXT2_MULTIPLICITY",
        "PARAMETER_TORSION_GENERATORS",
        "COLLISION_AXIS_DIFFERENCE",
        "PRESENTATION_PARAMETER_DIFFERENCE",
    }
    values: dict[str, int] = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key in wanted:
                values[key] = int(value)
    assert set(values) == wanted, result.stdout + result.stderr
    return (
        values["SATURATION_GENERATORS"],
        values["EXT2_MULTIPLICITY"],
        values["PARAMETER_TORSION_GENERATORS"],
        values["COLLISION_AXIS_DIFFERENCE"],
        values["PRESENTATION_PARAMETER_DIFFERENCE"],
    )


@cache
def generic_quartic_tensor() -> dict[tuple[int, int, int], sp.Expr]:
    """A deterministic generic-looking element of the order-four kernel."""

    triples = list(itertools.combinations_with_replacement(range(3), 3))
    pairs = list(itertools.combinations_with_replacement(range(3), 2))
    input_monomials = homogeneous_monomials(4)
    output_monomials = homogeneous_monomials(5)
    columns = {
        (triple, monomial): index
        for index, (triple, monomial) in enumerate(
            itertools.product(triples, input_monomials)
        )
    }
    rows = {
        (pair, monomial): index
        for index, (pair, monomial) in enumerate(
            itertools.product(pairs, output_monomials)
        )
    }
    constraint = sp.MutableSparseMatrix(len(rows), len(columns), {})
    relation_terms = ((z, 1), (y, -1), (x, 1))
    for pair in pairs:
        for first, (relation_variable, sign) in enumerate(relation_terms):
            triple = tuple(sorted((first, *pair)))
            for monomial in input_monomials:
                constraint[
                    rows[(pair, sp.expand(relation_variable * monomial))],
                    columns[(triple, monomial)],
                ] += sign
    kernel_basis = constraint.nullspace()
    assert len(kernel_basis) == 24
    kernel_vector = sum(
        (
            (index + 1) * vector
            for index, vector in enumerate(kernel_basis)
        ),
        sp.zeros(constraint.cols, 1),
    )
    denominators = [
        int(sp.denom(entry)) for entry in kernel_vector if entry != 0
    ]
    common_denominator = math.lcm(*denominators)
    kernel_vector *= common_denominator
    result = {
        triple: sp.expand(
            sum(
                kernel_vector[columns[(triple, monomial)]] * monomial
                for monomial in input_monomials
            )
        )
        for triple in triples
    }
    for pair in pairs:
        assert sp.expand(
            z * result[tuple(sorted((0, *pair)))]
            - y * result[tuple(sorted((1, *pair)))]
            + x * result[tuple(sorted((2, *pair)))]
        ) == 0
    return result


CUBIC_STRATA = {
    "smooth": X**3 + Y**3 + Z**3,
    "nodal": Y**2 * Z - X**2 * (X + Z),
    "cuspidal": Y**2 * Z - X**3,
    "line-transverse-conic": Z * (X * Y - Z**2),
    "line-tangent-conic": Z * (Y * Z - X**2),
    "triangle": X * Y * Z,
    "concurrent-lines": X * Y * (X - Y),
    "double-line": X**2 * Y,
    "triple-line": X**3,
    "zero": sp.Integer(0),
}

SQUAREFREE_STRATA = {
    "smooth",
    "nodal",
    "cuspidal",
    "line-transverse-conic",
    "line-tangent-conic",
    "triangle",
    "concurrent-lines",
}


def main() -> None:
    for stratum, representative in CUBIC_STRATA.items():
        (
            saturation_generators,
            support_dimension,
            ext3_length,
            ext2_dimension,
            ext2_length,
            ext2_top_dimension,
            ext2_square_action_generators,
        ) = run_singular(representative)
        assert saturation_generators == 0
        assert ext3_length == 0
        if stratum in SQUAREFREE_STRATA:
            assert support_dimension == 2
            assert (ext2_dimension, ext2_length) == (0, 6)
            assert (ext2_top_dimension, ext2_square_action_generators) == (
                3,
                0,
            )
        elif stratum == "zero":
            assert support_dimension == 3
            assert (ext2_dimension, ext2_length) == (-1, 0)
        else:
            assert support_dimension == 3
            assert ext2_dimension == 1
            assert ext2_length == -1
        print(
            f"PASS: homogeneous {stratum}: cotangent saturation=0, "
            f"support dimension={support_dimension}, "
            f"Ext2 dimension={ext2_dimension}, "
            f"Ext2 length={ext2_length}, "
            f"Ext2 top dimension={ext2_top_dimension}, "
            f"Ext2 m^2-action generators={ext2_square_action_generators}"
        )

    print(
        "PASS: all squarefree homogeneous symbols put the finite defect in C/T"
    )
    print(
        "PASS: double- and triple-line homogeneous symbols fail support purity"
    )
    print(
        "PASS: the zero homogeneous tensor has no support-bidual defect "
        "but is nowhere generically etale"
    )

    quartic_tensor = generic_quartic_tensor()
    for stratum, representative in CUBIC_STRATA.items():
        result = run_singular(representative, quartic_tensor)
        assert result == (0, 2, 0, 0, 6, 3, 0)
        print(
            f"PASS: quartic-lifted {stratum}: cotangent saturation=0, "
            "Ext2 length=6"
        )
    print(
        "PASS: one exact order-four lift puts every cubic-symbol stratum "
        "in the same finite support-defect row"
    )


if __name__ == "__main__":
    main()
