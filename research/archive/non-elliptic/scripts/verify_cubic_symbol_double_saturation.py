#!/usr/bin/env python3
"""Exact stratum audit for the reduced cubic Koszul defect.

For each projective ternary-cubic orbit representative, construct the
homogeneous generalized triple-cover algebra on

    M = coker(A --(z,-y,x)^T--> A^3),  A=Q[x,y,z].

The multiplication is recovered invariantly from its cubic tensor.  The
script then asks Singular for:

1. the (x,y,z)-saturation of an A-presentation of Omega_{B/A};
2. Ext_A^2(T,A), where T=B/Ann_B(Omega_{B/A}).

On polynomial quartic complement families it also computes the relative
Nakayama module J/nJ of the Kähler different J=Ann_B(Omega) at the collision.

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
FACTOR_SINGULAR_EXPRESSIONS = True


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

    # A full quartic complement can carry eight coefficient parameters.
    # Solving the cross-product systems over that symbolic field is much
    # slower than using the intrinsic degree bounds: the trace-free table is
    # linear and the scalar table is quadratic in the tensor.  Interpolate
    # those two exact polynomials from parameter-free evaluations.
    parameter_symbols = tuple(
        sorted(
            set().union(
                *(
                    sp.sympify(value).free_symbols
                    for value in (higher_tensor or {}).values()
                )
            ).difference(BASE_VARIABLES),
            key=str,
        )
    )
    if parameter_symbols:
        for value in (higher_tensor or {}).values():
            parameter_polynomial = sp.Poly(
                sp.expand(value), *parameter_symbols
            )
            assert parameter_polynomial.total_degree() <= 1

        zero_substitution = {
            parameter: sp.Integer(0) for parameter in parameter_symbols
        }

        def evaluated_table(
            values: dict[sp.Symbol, int]
        ) -> tuple[
            dict[tuple[int, int], sp.Matrix],
            dict[tuple[int, int], sp.Expr],
        ]:
            substitution = {**zero_substitution, **values}
            evaluated_tensor = {
                triple: sp.expand(sp.sympify(value).subs(substitution))
                for triple, value in (higher_tensor or {}).items()
            }
            return multiplication_table(cubic, evaluated_tensor)

        central_trace_free, central_scalar = evaluated_table({})
        positive_tables = {
            parameter: evaluated_table({parameter: 1})
            for parameter in parameter_symbols
        }
        negative_tables = {
            parameter: evaluated_table({parameter: -1})
            for parameter in parameter_symbols
        }

        trace_free_linear: dict[
            sp.Symbol, dict[tuple[int, int], sp.Matrix]
        ] = {}
        scalar_linear: dict[
            sp.Symbol, dict[tuple[int, int], sp.Expr]
        ] = {}
        scalar_diagonal: dict[
            sp.Symbol, dict[tuple[int, int], sp.Expr]
        ] = {}
        for parameter in parameter_symbols:
            positive_trace_free, positive_scalar = positive_tables[
                parameter
            ]
            negative_trace_free, negative_scalar = negative_tables[
                parameter
            ]
            trace_free_linear[parameter] = {
                pair: (
                    (positive_trace_free[pair] - negative_trace_free[pair])
                    / 2
                ).applyfunc(sp.expand)
                for pair in central_trace_free
            }
            scalar_linear[parameter] = {
                pair: sp.expand(
                    (positive_scalar[pair] - negative_scalar[pair]) / 2
                )
                for pair in central_scalar
            }
            scalar_diagonal[parameter] = {
                pair: sp.expand(
                    (positive_scalar[pair] + negative_scalar[pair]) / 2
                    - central_scalar[pair]
                )
                for pair in central_scalar
            }

        trace_free_products = {
            pair: (
                central_trace_free[pair]
                + sum(
                    (
                        parameter * trace_free_linear[parameter][pair]
                        for parameter in parameter_symbols
                    ),
                    sp.zeros(3, 1),
                )
            ).applyfunc(sp.expand)
            for pair in central_trace_free
        }
        scalar_products = {
            pair: sp.expand(
                central_scalar[pair]
                + sum(
                    parameter * scalar_linear[parameter][pair]
                    + parameter**2
                    * scalar_diagonal[parameter][pair]
                    for parameter in parameter_symbols
                )
            )
            for pair in central_scalar
        }
        for first_index, first_parameter in enumerate(parameter_symbols):
            for second_parameter in parameter_symbols[first_index + 1 :]:
                _pair_trace_free, pair_scalar = evaluated_table(
                    {first_parameter: 1, second_parameter: 1}
                )
                for pair in scalar_products:
                    cross_coefficient = sp.expand(
                        pair_scalar[pair]
                        - central_scalar[pair]
                        - scalar_linear[first_parameter][pair]
                        - scalar_linear[second_parameter][pair]
                        - scalar_diagonal[first_parameter][pair]
                        - scalar_diagonal[second_parameter][pair]
                    )
                    scalar_products[pair] = sp.expand(
                        scalar_products[pair]
                        + first_parameter
                        * second_parameter
                        * cross_coefficient
                    )
        return trace_free_products, scalar_products

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


def trace_free_generator_action_matrices(
    cubic: sp.Expr,
    higher_tensor: dict[tuple[int, int, int], sp.Expr] | None = None,
) -> tuple[sp.Matrix, ...]:
    """Matrices for multiplication by e_1,e_2,e_3 on B=A plus M.

    The coefficient module is presented on ``1,e_1,e_2,e_3`` with the
    single relation ``z e_1-y e_2+x e_3``.  These matrices let the Singular
    certificate form ``n Ann_B(Omega)`` inside that fixed presentation.
    """

    trace_free_products, scalar_products = multiplication_table(
        cubic, higher_tensor
    )
    actions = []
    for generator in range(3):
        columns = [
            sp.Matrix(
                [
                    0,
                    *(sp.Integer(index == generator) for index in range(3)),
                ]
            )
        ]
        for second in range(3):
            pair = tuple(sorted((generator, second)))
            columns.append(
                sp.Matrix(
                    [
                        scalar_products[pair],
                        *trace_free_products[pair],
                    ]
                )
            )
        actions.append(
            sp.Matrix.hstack(*columns).applyfunc(sp.expand)
        )
    return tuple(actions)


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


@cache
def _cached_differential_relations(
    cubic: sp.Expr,
    higher_items: tuple[
        tuple[tuple[int, int, int], sp.Expr], ...
    ],
) -> tuple[tuple[sp.Expr, ...], ...]:
    higher_tensor = dict(higher_items) if higher_items else None
    return tuple(
        tuple(relation)
        for relation in differential_relations(cubic, higher_tensor)
    )


def cached_differential_relations(
    cubic: sp.Expr,
    higher_tensor: dict[tuple[int, int, int], sp.Expr] | None = None,
) -> list[list[sp.Expr]]:
    """Memoize the expensive symbolic family presentation in one replay."""

    higher_items = tuple(
        sorted(
            (
                (triple, sp.expand(value))
                for triple, value in (higher_tensor or {}).items()
            ),
            key=lambda item: item[0],
        )
    )
    return [
        list(relation)
        for relation in _cached_differential_relations(
            cubic, higher_items
        )
    ]


def singular_polynomial(expression: sp.Expr) -> str:
    normalized = (
        sp.factor(expression)
        if FACTOR_SINGULAR_EXPRESSIONS
        else sp.expand(expression)
    )
    return sp.sstr(normalized).replace("**", "^")


def singular_vector(vector: list[sp.Expr]) -> str:
    return "[" + ",".join(map(singular_polynomial, vector)) + "]"


def filtered_rees_vector(
    vector: list[sp.Expr],
    component_weights: tuple[int, ...],
    parameters: tuple[sp.Symbol, ...],
    rees_parameter: sp.Symbol,
    parameter_collision_weights: tuple[int, ...] | None = None,
) -> list[sp.Expr]:
    """Return the weighted Rees transform of one filtered relation.

    The collision variables ``x,y,z`` have weight one.  Deformation
    parameters have weight zero by default, but callers may give them
    positive collision weights.  The latter realizes a universal filtered
    coefficient such as the nodal normal-form variable ``u`` of weight one.
    The returned vector is normalized by its least filtered degree, so
    specialization at the Rees parameter zero is its initial relation and
    specialization at one is the original relation.
    """

    assert len(vector) == len(component_weights)
    parameter_collision_weights = parameter_collision_weights or (
        0,
    ) * len(parameters)
    assert len(parameters) == len(parameter_collision_weights)
    least_degree: int | None = None
    for raw_entry, component_weight in zip(vector, component_weights):
        entry = sp.sympify(raw_entry)
        if entry == 0:
            continue
        polynomial = sp.Poly(
            sp.expand(entry), *BASE_VARIABLES, *parameters
        )
        for monomial, _coefficient in polynomial.terms():
            filtered_degree = (
                sum(monomial[:3])
                + sum(
                    exponent * weight
                    for exponent, weight in zip(
                        monomial[3:], parameter_collision_weights
                    )
                )
                + component_weight
            )
            least_degree = (
                filtered_degree
                if least_degree is None
                else min(least_degree, filtered_degree)
            )
    assert least_degree is not None

    substitution = {
        variable: rees_parameter * variable for variable in BASE_VARIABLES
    }
    substitution.update(
        {
            parameter: rees_parameter**weight * parameter
            for parameter, weight in zip(
                parameters, parameter_collision_weights
            )
            if weight
        }
    )
    return [
        sp.expand(
            sp.sympify(raw_entry).subs(
                substitution, simultaneous=True
            )
            * rees_parameter ** (component_weight - least_degree)
        )
        for raw_entry, component_weight in zip(
            vector, component_weights
        )
    ]


def singular_rees_base_change_program(
    cubic: sp.Expr,
    tensors: tuple[dict[tuple[int, int, int], sp.Expr], ...],
    parameter_collision_weights: tuple[int, ...] | None = None,
) -> tuple[str, tuple[sp.Symbol, ...]]:
    """Build the strict-Rees certificate for cotangent and annihilator base change.

    The second module is the cokernel of ``B -> Omega^3``.  Its flatness,
    together with flatness of ``Omega``, forces the annihilator kernel and
    hence ``T=B/Ann(Omega)`` to commute with every parameter base change.
    """

    assert tensors
    parameters = sp.symbols(f"p0:{len(tensors)}")
    parameter_collision_weights = parameter_collision_weights or (
        0,
    ) * len(parameters)
    assert len(parameters) == len(parameter_collision_weights)
    rees_parameter = sp.Symbol("rees_t")
    family_tensor = {
        triple: sp.expand(
            sum(
                parameter * tensor[triple]
                for parameter, tensor in zip(parameters, tensors)
            )
        )
        for triple in tensors[0]
    }
    family_relations = cached_differential_relations(cubic, family_tensor)
    central_relations = cached_differential_relations(cubic)

    # deg(de_i)=2 and deg(e_j de_i)=4 in the central graded algebra.
    cotangent_weights = tuple(
        2 if component % 4 == 0 else 4 for component in range(12)
    )
    cotangent_rees = tuple(
        filtered_rees_vector(
            relation,
            cotangent_weights,
            parameters,
            rees_parameter,
            parameter_collision_weights,
        )
        for relation in family_relations
    )

    annihilator_cokernel_relations: list[list[sp.Expr]] = []
    central_annihilator_cokernel_relations: list[list[sp.Expr]] = []
    for block in range(3):
        for family_relation, central_relation in zip(
            family_relations, central_relations
        ):
            family_block = [sp.Integer(0)] * 36
            central_block = [sp.Integer(0)] * 36
            family_block[12 * block : 12 * (block + 1)] = (
                family_relation
            )
            central_block[12 * block : 12 * (block + 1)] = (
                central_relation
            )
            annihilator_cokernel_relations.append(family_block)
            central_annihilator_cokernel_relations.append(central_block)
    for algebra_generator in range(4):
        action = [sp.Integer(0)] * 36
        for block in range(3):
            action[
                12 * block + 4 * block + algebra_generator
            ] = 1
        annihilator_cokernel_relations.append(action)
        central_annihilator_cokernel_relations.append(action)

    annihilator_cokernel_weights = cotangent_weights * 3
    annihilator_cokernel_rees = tuple(
        filtered_rees_vector(
            relation,
            annihilator_cokernel_weights,
            parameters,
            rees_parameter,
            parameter_collision_weights,
        )
        for relation in annihilator_cokernel_relations
    )

    variable_names = ",".join(
        map(str, (*parameters, *BASE_VARIABLES, rees_parameter))
    )
    program = f'''\
LIB "primdec.lib";
ring rees_ring=0,({variable_names}),dp;
'''
    packets = (
        (
            "COTANGENT",
            cotangent_rees,
            central_relations,
        ),
        (
            "ANNIHILATOR_COKERNEL",
            annihilator_cokernel_rees,
            central_annihilator_cokernel_relations,
        ),
    )
    for label, rees_relations, central in packets:
        program += (
            f"module {label}_REES="
            + ",".join(map(singular_vector, rees_relations))
            + ";\n"
            + f"module {label}_CENTRAL="
            + ",".join(map(singular_vector, central))
            + ";\n"
            + f"{label}_REES=std({label}_REES);\n"
            + f"{label}_CENTRAL=std({label}_CENTRAL);\n"
            + f"list {label}_SATURATION_DATA="
            + f"sat({label}_REES,ideal({rees_parameter}));\n"
            + f"module {label}_SATURATED="
            + f"{label}_SATURATION_DATA[1];\n"
            + f"module {label}_SATURATION_QUOTIENT=simplify("
            + f"reduce({label}_SATURATED,{label}_REES),2);\n"
            + f"module {label}_INITIAL=std(subst("
            + f"{label}_SATURATED,{rees_parameter},0));\n"
            + f"module {label}_INITIAL_MINUS_CENTRAL=simplify("
            + f"reduce({label}_INITIAL,{label}_CENTRAL),2);\n"
            + f"module {label}_CENTRAL_MINUS_INITIAL=simplify("
            + f"reduce({label}_CENTRAL,{label}_INITIAL),2);\n"
            + f'print("{label}_REES_TORSION_GENERATORS="'
            + f"+string(size({label}_SATURATION_QUOTIENT)));\n"
            + f'print("{label}_INITIAL_PRESENTATION_DIFFERENCE="'
            + f"+string(size({label}_INITIAL_MINUS_CENTRAL)"
            + f"+size({label}_CENTRAL_MINUS_INITIAL)));\n"
        )
    program += "quit;\n"
    return program, parameters


def run_singular_rees_base_change_certificate(
    cubic: sp.Expr,
    tensors: tuple[dict[tuple[int, int, int], sp.Expr], ...],
    timeout: int | None = None,
    parameter_collision_weights: tuple[int, ...] | None = None,
) -> dict[str, int]:
    """Verify strict filtered base change on a full quartic subspace."""

    program, parameters = singular_rees_base_change_program(
        cubic,
        tensors,
        parameter_collision_weights,
    )
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for this checker"
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    wanted = {
        "COTANGENT_REES_TORSION_GENERATORS",
        "COTANGENT_INITIAL_PRESENTATION_DIFFERENCE",
        "ANNIHILATOR_COKERNEL_REES_TORSION_GENERATORS",
        "ANNIHILATOR_COKERNEL_INITIAL_PRESENTATION_DIFFERENCE",
    }
    values: dict[str, int] = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key in wanted:
                values[key] = int(value)
    assert set(values) == wanted, result.stdout + result.stderr
    return {
        "parameter_count": len(parameters),
        "cotangent_rees_torsion_generators": values[
            "COTANGENT_REES_TORSION_GENERATORS"
        ],
        "cotangent_initial_presentation_difference": values[
            "COTANGENT_INITIAL_PRESENTATION_DIFFERENCE"
        ],
        "annihilator_cokernel_rees_torsion_generators": values[
            "ANNIHILATOR_COKERNEL_REES_TORSION_GENERATORS"
        ],
        "annihilator_cokernel_initial_presentation_difference": values[
            "ANNIHILATOR_COKERNEL_INITIAL_PRESENTATION_DIFFERENCE"
        ],
    }


def singular_program(
    cubic: sp.Expr,
    higher_tensor: dict[tuple[int, int, int], sp.Expr] | None = None,
) -> str:
    relations = cached_differential_relations(cubic, higher_tensor)
    generator_actions = trace_free_generator_action_matrices(
        cubic, higher_tensor
    )

    # To compute Ann_B(Q), take the kernel of
    # B -> Q^3, b |-> (b de_1,b de_2,b de_3).  Singular's modulo(H1,H2)
    # returns the kernel of the map generated by H1 into coker(H2), so it
    # computes the required preimage in A^4 directly.  This avoids forming
    # a much larger syzygy module and then extracting its first four rows.
    annihilator_action_generators: list[list[sp.Expr]] = []
    for algebra_generator in range(4):
        relation = [sp.Integer(0)] * 36
        for differential in range(3):
            relation[
                12 * differential
                + 4 * differential
                + algebra_generator
            ] = 1
        annihilator_action_generators.append(relation)
    annihilator_quotient_relations: list[list[sp.Expr]] = []
    for block in range(3):
        for source_relation in relations:
            relation = [sp.Integer(0)] * 36
            relation[12 * block : 12 * (block + 1)] = source_relation
            annihilator_quotient_relations.append(relation)

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

module annihilator_action=
{",".join(map(singular_vector, annihilator_action_generators))};
module annihilator_relations=
{",".join(map(singular_vector, annihilator_quotient_relations))};
module support_presentation=std(
  modulo(annihilator_action,annihilator_relations)
);
print("SUPPORT_DIMENSION="+string(dim(support_presentation)));
module coefficient_relation=[0,z,-y,x];
matrix e1_action[4][4]=
{",".join(singular_polynomial(entry) for entry in generator_actions[0])};
matrix e2_action[4][4]=
{",".join(singular_polynomial(entry) for entry in generator_actions[1])};
matrix e3_action[4][4]=
{",".join(singular_polynomial(entry) for entry in generator_actions[2])};
module maximal_times_annihilator=
  x*support_presentation+y*support_presentation+z*support_presentation;
module trace_free_times_annihilator=
  e1_action*support_presentation+
  e2_action*support_presentation+
  e3_action*support_presentation;
module different_nakayama_denominator=std(
  coefficient_relation+maximal_times_annihilator+
  trace_free_times_annihilator
);
module different_generator_presentation=std(
  modulo(support_presentation,different_nakayama_denominator)
);
print("DIFFERENT_GENERATOR_DIMENSION="+
  string(dim(different_generator_presentation)));
print("DIFFERENT_GENERATOR_MULTIPLICITY="+
  string(mult(different_generator_presentation)));
print("DIFFERENT_MINIMAL_GENERATORS="+
  string(vdim(different_generator_presentation)));
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


def singular_family_program(
    cubic: sp.Expr,
    quartic_tensor: dict[tuple[int, int, int], sp.Expr] | None = None,
) -> str:
    """The exact family phi_h+t*psi_4 over Q[t,x,y,z]."""

    deformation_parameter = sp.symbols("t")
    family_tensor = {
        triple: deformation_parameter * value
        for triple, value in (
            quartic_tensor or generic_quartic_tensor()
        ).items()
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


def run_singular_family(
    cubic: sp.Expr,
    quartic_tensor: dict[tuple[int, int, int], sp.Expr] | None = None,
    timeout: int | None = None,
) -> tuple[int, int, int, int, int]:
    """Return uniform cotangent, support, multiplicity, and t-torsion data."""

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for this checker"
    result = subprocess.run(
        [singular, "-q"],
        input=singular_family_program(cubic, quartic_tensor),
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
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


def singular_subspace_program(
    cubic: sp.Expr,
    tensors: tuple[dict[tuple[int, int, int], sp.Expr], ...],
) -> tuple[str, tuple[sp.Symbol, ...]]:
    """Build the exact programme for a polynomial quartic subspace."""

    assert tensors
    parameters = sp.symbols(f"p0:{len(tensors)}")
    subspace_tensor = {
        triple: sp.expand(
            sum(
                parameter * tensor[triple]
                for parameter, tensor in zip(parameters, tensors)
            )
        )
        for triple in tensors[0]
    }
    parameter_names = ",".join(map(str, parameters))
    program = singular_program(cubic, subspace_tensor).replace(
        "ring coefficient_ring=0,(x,y,z),dp;",
        f"ring coefficient_ring=0,({parameter_names},x,y,z),dp;",
    )
    diagnostic_anchor = (
        'print("EXT2_VECTOR_DIMENSION="'
        "+string(vdim(support_ext2)));"
    )
    program = program.replace(
        diagnostic_anchor,
        diagnostic_anchor
        + 'print("EXT2_MULTIPLICITY="+string(mult(support_ext2)));'
        + "ideal ext2_fitting=fitting(support_ext2,0);"
        + "ideal ext2_support=std(radical(ext2_fitting));"
        + "ideal parameter_space=std(ideal(x,y,z));"
        + "ideal first_support_difference=simplify("
        + "reduce(ext2_support,parameter_space),2);"
        + "ideal second_support_difference=simplify("
        + "reduce(parameter_space,ext2_support),2);"
        + 'print("PARAMETER_SPACE_DIFFERENCE="'
        + "+string(size(first_support_difference)"
        + "+size(second_support_difference)));"
        + "module pruned_ext2_presentation=std(prune(support_ext2));"
        + "module central_ext2_presentation=std(prune("
        + "".join(f"subst(" for _ in parameters)
        + "support_ext2"
        + "".join(f",{parameter},0)" for parameter in parameters)
        + "));"
        + "module first_presentation_difference=simplify("
        + "reduce(pruned_ext2_presentation,"
        + "central_ext2_presentation),2);"
        + "module second_presentation_difference=simplify("
        + "reduce(central_ext2_presentation,"
        + "pruned_ext2_presentation),2);"
        + 'print("PRUNED_PRESENTATION_DIFFERENCE="'
        + "+string(size(first_presentation_difference)"
        + "+size(second_presentation_difference)));"
        + 'print("PRUNED_PRESENTATION_RANK="'
        + "+string(nrows(pruned_ext2_presentation)));",
    )
    different_anchor = (
        "module support_ext3=std(Ext_R(3,support_presentation));"
    )
    program = program.replace(
        different_anchor,
        "module pruned_different_generator_presentation=std("
        + "prune(different_generator_presentation));"
        + "ideal different_generator_fitting=fitting("
        + "pruned_different_generator_presentation,0);"
        + "ideal different_generator_support=std(radical("
        + "different_generator_fitting));"
        + "ideal different_parameter_space=std(ideal(x,y,z));"
        + "ideal first_different_support_difference=simplify("
        + "reduce(different_generator_support,"
        + "different_parameter_space),2);"
        + "ideal second_different_support_difference=simplify("
        + "reduce(different_parameter_space,"
        + "different_generator_support),2);"
        + 'print("DIFFERENT_PARAMETER_SPACE_DIFFERENCE="'
        + "+string(size(first_different_support_difference)"
        + "+size(second_different_support_difference)));"
        + "module central_different_generator_presentation=std(prune("
        + "".join(f"subst(" for _ in parameters)
        + "different_generator_presentation"
        + "".join(f",{parameter},0)" for parameter in parameters)
        + "));"
        + "module first_different_presentation_difference=simplify("
        + "reduce(pruned_different_generator_presentation,"
        + "central_different_generator_presentation),2);"
        + "module second_different_presentation_difference=simplify("
        + "reduce(central_different_generator_presentation,"
        + "pruned_different_generator_presentation),2);"
        + 'print("DIFFERENT_PRUNED_PRESENTATION_DIFFERENCE="'
        + "+string(size(first_different_presentation_difference)"
        + "+size(second_different_presentation_difference)));"
        + 'print("DIFFERENT_PRUNED_PRESENTATION_RANK="'
        + "+string(nrows(pruned_different_generator_presentation)));"
        + different_anchor,
    )
    return program, parameters


def run_singular_subspace_certificate(
    cubic: sp.Expr,
    tensors: tuple[dict[tuple[int, int, int], sp.Expr], ...],
    timeout: int | None = None,
) -> dict[str, int]:
    """Audit both saturation layers on a full quartic parameter space."""

    program, parameters = singular_subspace_program(cubic, tensors)
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for this checker"
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    wanted = {
        "SATURATION_GENERATORS",
        "EXT2_MULTIPLICITY",
        "PARAMETER_SPACE_DIFFERENCE",
        "PRUNED_PRESENTATION_DIFFERENCE",
        "PRUNED_PRESENTATION_RANK",
        "SUPPORT_DIMENSION",
        "EXT3_VECTOR_DIMENSION",
        "EXT2_DIMENSION",
        "EXT2_SQUARE_ACTION_GENERATORS",
        "DIFFERENT_GENERATOR_DIMENSION",
        "DIFFERENT_GENERATOR_MULTIPLICITY",
        "DIFFERENT_PARAMETER_SPACE_DIFFERENCE",
        "DIFFERENT_PRUNED_PRESENTATION_DIFFERENCE",
        "DIFFERENT_PRUNED_PRESENTATION_RANK",
    }
    values: dict[str, int] = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key in wanted:
                values[key] = int(value)
    assert set(values) == wanted, result.stdout + result.stderr
    return {
        "parameter_count": len(parameters),
        "cotangent_saturation_generators": values[
            "SATURATION_GENERATORS"
        ],
        "support_module_dimension": values["SUPPORT_DIMENSION"],
        "support_ext3_vector_dimension": values[
            "EXT3_VECTOR_DIMENSION"
        ],
        "support_ext2_dimension": values["EXT2_DIMENSION"],
        "support_ext2_multiplicity": values["EXT2_MULTIPLICITY"],
        "support_ext2_parameter_axis_radical_difference": values[
            "PARAMETER_SPACE_DIFFERENCE"
        ],
        "support_ext2_central_pruned_presentation_difference": values[
            "PRUNED_PRESENTATION_DIFFERENCE"
        ],
        "support_ext2_pruned_presentation_rank": values[
            "PRUNED_PRESENTATION_RANK"
        ],
        "support_ext2_collision_square_action_generators": values[
            "EXT2_SQUARE_ACTION_GENERATORS"
        ],
        "different_generator_module_dimension": values[
            "DIFFERENT_GENERATOR_DIMENSION"
        ],
        "different_generator_module_multiplicity": values[
            "DIFFERENT_GENERATOR_MULTIPLICITY"
        ],
        "different_generator_parameter_axis_radical_difference": values[
            "DIFFERENT_PARAMETER_SPACE_DIFFERENCE"
        ],
        "different_generator_central_pruned_presentation_difference": values[
            "DIFFERENT_PRUNED_PRESENTATION_DIFFERENCE"
        ],
        "different_generator_pruned_presentation_rank": values[
            "DIFFERENT_PRUNED_PRESENTATION_RANK"
        ],
    }


def run_singular_subspace(
    cubic: sp.Expr,
    tensors: tuple[dict[tuple[int, int, int], sp.Expr], ...],
    timeout: int | None = None,
) -> tuple[int, int, int, int, int]:
    """Audit a full polynomial parameter subspace of quartic tensors."""

    certificate = run_singular_subspace_certificate(
        cubic,
        tensors,
        timeout,
    )
    return (
        certificate["cotangent_saturation_generators"],
        certificate["support_ext2_multiplicity"],
        certificate["support_ext2_parameter_axis_radical_difference"],
        certificate[
            "support_ext2_central_pruned_presentation_difference"
        ],
        certificate["support_ext2_pruned_presentation_rank"],
    )


def run_singular_plane(
    cubic: sp.Expr,
    first_tensor: dict[tuple[int, int, int], sp.Expr],
    second_tensor: dict[tuple[int, int, int], sp.Expr],
    timeout: int | None = None,
) -> tuple[int, int, int, int, int]:
    """Audit a two-dimensional polynomial quartic parameter space."""

    return run_singular_subspace(
        cubic,
        (first_tensor, second_tensor),
        timeout,
    )


def run_plane_parameter_fitting(
    cubic: sp.Expr,
    first_tensor: dict[tuple[int, int, int], sp.Expr],
    second_tensor: dict[tuple[int, int, int], sp.Expr],
    timeout: int | None = None,
) -> tuple[int, int, int]:
    """Test S-flat rank six after m^2 truncation, S=Q[p0,p1]."""

    first_parameter, second_parameter = sp.symbols("p0 p1")
    plane_tensor = {
        triple: sp.expand(
            first_parameter * first_tensor[triple]
            + second_parameter * second_tensor[triple]
        )
        for triple in first_tensor
    }
    program = singular_program(cubic, plane_tensor).replace(
        "ring coefficient_ring=0,(x,y,z),dp;",
        "ring coefficient_ring=0,(p0,p1,x,y,z),dp;",
    )
    program = program.replace(
        "quit;",
        "module pruned_ext2=std(prune(support_ext2));"
        "int fitting_row,fitting_column;"
        'print("FITTING_ROWS="+string(nrows(pruned_ext2)));'
        'print("FITTING_COLUMNS="+string(ncols(pruned_ext2)));'
        "module maximal_square_free="
        "(x2+xy+xz+y2+yz+z2)*freemodule(nrows(pruned_ext2));"
        "module maximal_square_action=simplify("
        "reduce(maximal_square_free,pruned_ext2),2);"
        'print("FITTING_M2_ACTION="'
        "+string(size(maximal_square_action)));"
        "for(fitting_column=1;"
        "fitting_column<=ncols(pruned_ext2);"
        "fitting_column++)"
        "{for(fitting_row=1;"
        "fitting_row<=nrows(pruned_ext2);"
        "fitting_row++)"
        "{print("
        '"FITTING_ENTRY_"+string(fitting_row)+"_"'
        "+string(fitting_column)+"
        '"="+string(pruned_ext2[fitting_column][fitting_row]));}}'
        "quit;",
    )
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for this checker"
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    values: dict[str, str] = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key.startswith("FITTING_"):
                values[key] = value
    row_count = int(values["FITTING_ROWS"])
    column_count = int(values["FITTING_COLUMNS"])
    maximal_square_action = int(values["FITTING_M2_ACTION"])
    assert maximal_square_action == 0

    parameter_locals = {
        "p0": first_parameter,
        "p1": second_parameter,
        "x": x,
        "y": y,
        "z": z,
    }
    entries = {
        (row, column): sp.sympify(
            values[f"FITTING_ENTRY_{row}_{column}"].replace("^", "**"),
            locals=parameter_locals,
        )
        for column in range(1, column_count + 1)
        for row in range(1, row_count + 1)
    }
    truncated_basis = (sp.Integer(1), x, y, z)
    parameter_relations: list[list[sp.Expr]] = []
    for column in range(1, column_count + 1):
        for multiplier in truncated_basis:
            relation = [sp.Integer(0)] * (4 * row_count)
            for row in range(1, row_count + 1):
                polynomial = sp.Poly(
                    sp.expand(multiplier * entries[row, column]),
                    x,
                    y,
                    z,
                )
                for monomial_index, monomial in enumerate(
                    truncated_basis
                ):
                    relation[
                        4 * (row - 1) + monomial_index
                    ] = sp.expand(polynomial.coeff_monomial(monomial))
            if any(relation):
                parameter_relations.append(relation)

    parameter_program = f"""
LIB "homolog.lib";
ring parameter_ring=0,(p0,p1),dp;
module parameter_presentation=
{",".join(map(singular_vector, parameter_relations))};
parameter_presentation=std(parameter_presentation);
ideal fitting_six=std(fitting(parameter_presentation,6));
ideal fitting_five=std(fitting(parameter_presentation,5));
ideal fitting_six_difference=simplify(
  reduce(ideal(1),fitting_six),2
);
print("FITTING_SIX_UNIT_DIFFERENCE="
  +string(size(fitting_six_difference)));
print("FITTING_FIVE_GENERATORS="+string(size(fitting_five)));
quit;
"""
    parameter_result = subprocess.run(
        [singular, "-q"],
        input=parameter_program,
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    parameter_values: dict[str, int] = {}
    for line in parameter_result.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key in {
                "FITTING_SIX_UNIT_DIFFERENCE",
                "FITTING_FIVE_GENERATORS",
            }:
                parameter_values[key] = int(value)
    assert set(parameter_values) == {
        "FITTING_SIX_UNIT_DIFFERENCE",
        "FITTING_FIVE_GENERATORS",
    }, parameter_result.stdout + parameter_result.stderr
    return (
        maximal_square_action,
        parameter_values["FITTING_SIX_UNIT_DIFFERENCE"],
        parameter_values["FITTING_FIVE_GENERATORS"],
    )


@cache
def quartic_constraint_data() -> tuple[
    list[tuple[int, int, int]],
    list[sp.Expr],
    dict[tuple[tuple[int, int, int], sp.Expr], int],
    list[sp.Matrix],
]:
    """Return monomial coordinates and a basis of the order-four kernel."""

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
    return triples, input_monomials, columns, kernel_basis


@cache
def quartic_kernel_basis_tensors() -> tuple[
    dict[tuple[int, int, int], sp.Expr], ...
]:
    """A primitive integral basis of the 24-dimensional quartic kernel."""

    (
        triples,
        input_monomials,
        columns,
        kernel_basis,
    ) = quartic_constraint_data()
    result = []
    for vector in kernel_basis:
        denominators = [
            int(sp.denom(entry)) for entry in vector if entry != 0
        ]
        integral_vector = math.lcm(*denominators) * vector
        tensor = {
            triple: sp.expand(
                sum(
                    integral_vector[columns[(triple, monomial)]] * monomial
                    for monomial in input_monomials
                )
            )
            for triple in triples
        }
        result.append(tensor)
    return tuple(result)


@cache
def generic_quartic_tensor() -> dict[tuple[int, int, int], sp.Expr]:
    """A deterministic generic-looking element of the order-four kernel."""

    (
        triples,
        input_monomials,
        columns,
        kernel_basis,
    ) = quartic_constraint_data()
    kernel_vector = sum(
        (
            (index + 1) * vector
            for index, vector in enumerate(kernel_basis)
        ),
        sp.zeros(kernel_basis[0].rows, 1),
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
    for pair in itertools.combinations_with_replacement(range(3), 2):
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
