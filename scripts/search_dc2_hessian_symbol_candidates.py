#!/usr/bin/env python3
"""Search a small Hessian-to-symplectic bridge for DC_2 symbols.

This is an exact finite experiment, not a DC_2 obstruction theorem.  It has
two stages.

First, for a polynomial potential ``A`` in four Darboux variables, form the
Hamiltonian Hessian matrix ``N = Pi*Hess(A)``.  When ``N^4=0``, its pointwise
symplectic Cayley transform is polynomial.  We retain it only when its rows
are closed, so that it is the Jacobian of a polynomial map.  The search
enumerates every one- and two-monomial cubic/quartic support over its natural
affine coefficient base.

Second, the surviving one-pencil maps are composed in noncommuting pairs.
These two-step words are known polynomial automorphism controls.  They are
nevertheless useful for calibrating a symbol-hostility score: their raw
classical coordinates can have nonzero Moyal defects even though a
factor-by-factor Weyl lift exists.  The best preliminary rows are quantized
through the native-support hbar^5 strong-span test.

All ranks and polynomial identities are over Q.  Scores are meaningful only
inside the explicitly declared native-support lattice.
"""

from __future__ import annotations

import argparse
import itertools
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]

Q1, Q2, P1, P2 = VARIABLES = sp.symbols("q1 q2 p1 p2")
A_PARAMETER, B_PARAMETER = PARAMETERS = sp.symbols("a b")

# Bracket convention {p_i,q_j}=delta_ij, as in the DC_2 notes.
POISSON = sp.Matrix(
    (
        (0, 0, -1, 0),
        (0, 0, 0, -1),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
    )
)
PAIRS = tuple(itertools.combinations(range(4), 2))


Exponent = tuple[int, int, int, int]
TaggedExponent = tuple[int, Exponent]
SparseColumn = dict[TaggedExponent, sp.Expr]


@dataclass(frozen=True)
class PencilCount:
    degree_pattern: str
    templates: int
    nilpotent: int
    cayley_integrable: int
    square_zero: int
    moyal_flat: int


@dataclass(frozen=True)
class PreliminaryWord:
    key: str
    first_potential: str
    second_potential: str
    maximum_degree: int
    map_terms: int
    moyal_defect_terms: int
    first_nonzero_moyal_order: int | None
    localization_pole_width: int
    parameter_base_complexity: int
    preliminary_score: int


def monomials_of_degree(degree: int) -> tuple[sp.Expr, ...]:
    return tuple(
        sp.prod(variable**exponent for variable, exponent in zip(VARIABLES, exponents))
        for exponents in itertools.product(range(degree + 1), repeat=4)
        if sum(exponents) == degree
    )


def expression_key(expression: sp.Expr) -> str:
    return str(expression).replace("**", "^").replace("*", "_")


def polynomial_coefficients(expression: sp.Expr) -> dict[Exponent, sp.Expr]:
    if expression == 0:
        return {}
    return {
        monomial: coefficient
        for monomial, coefficient in sp.Poly(expression, VARIABLES).terms()
        if coefficient
    }


def poisson(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    return sp.expand(
        sum(
            POISSON[i, j]
            * sp.diff(left, VARIABLES[i])
            * sp.diff(right, VARIABLES[j])
            for i in range(4)
            for j in range(4)
            if POISSON[i, j]
        )
    )


def pi_power(left: sp.Expr, right: sp.Expr, power: int) -> sp.Expr:
    """Apply the constant Poisson bidifferential ``power`` times."""

    terms = tuple(
        (i, j, POISSON[i, j])
        for i in range(4)
        for j in range(4)
        if POISSON[i, j]
    )
    zero = (0, 0, 0, 0)
    states: dict[tuple[Exponent, Exponent], sp.Expr] = {(zero, zero): sp.Integer(1)}
    for _ in range(power):
        next_states: dict[tuple[Exponent, Exponent], sp.Expr] = {}
        for (left_orders, right_orders), coefficient in states.items():
            for left_index, right_index, sign in terms:
                new_left = list(left_orders)
                new_right = list(right_orders)
                new_left[left_index] += 1
                new_right[right_index] += 1
                key = (tuple(new_left), tuple(new_right))
                next_states[key] = next_states.get(key, 0) + sign * coefficient
        states = next_states

    value = 0
    for (left_orders, right_orders), coefficient in states.items():
        left_derivative = left
        right_derivative = right
        for variable, order in zip(VARIABLES, left_orders):
            if order:
                left_derivative = sp.diff(left_derivative, variable, order)
        for variable, order in zip(VARIABLES, right_orders):
            if order:
                right_derivative = sp.diff(right_derivative, variable, order)
        value += coefficient * left_derivative * right_derivative
    return sp.expand(value)


def hamiltonian_matrix(potential: sp.Expr) -> sp.Matrix:
    return POISSON * sp.hessian(potential, VARIABLES)


def is_zero_matrix(matrix: sp.Matrix) -> bool:
    return all(sp.expand(entry) == 0 for entry in matrix)


def cayley_jacobian(nilpotent: sp.Matrix) -> sp.Matrix:
    # For N^4=0, (I-N/2)^(-1)(I+N/2)
    # = I+N+N^2/2+N^3/4.
    square = nilpotent * nilpotent
    return sp.eye(4) + nilpotent + square / 2 + square * nilpotent / 4


def rows_are_closed(jacobian: sp.Matrix) -> bool:
    return all(
        sp.expand(
            sp.diff(jacobian[row, left], VARIABLES[right])
            - sp.diff(jacobian[row, right], VARIABLES[left])
        )
        == 0
        for row in range(4)
        for left in range(4)
        for right in range(left + 1, 4)
    )


def integrate_closed_rows(jacobian: sp.Matrix) -> tuple[sp.Expr, ...]:
    """Integrate closed polynomial rows by the radial homotopy."""

    outputs = []
    for row in range(4):
        one_form_contraction = sp.expand(
            sum(jacobian[row, column] * VARIABLES[column] for column in range(4))
        )
        value = 0
        coefficients = polynomial_coefficients(one_form_contraction)
        for exponents, coefficient in coefficients.items():
            monomial = sp.prod(
                variable**exponent
                for variable, exponent in zip(VARIABLES, exponents)
            )
            # The contraction already contributes one radial coordinate:
            # integral_0^1 C(t*x) x dt divides a contraction monomial of
            # total degree d by d.
            value += coefficient * monomial / sum(exponents)
        outputs.append(sp.expand(value))
    result = tuple(outputs)
    assert sp.Matrix(result).jacobian(VARIABLES) == jacobian
    return result


def euler_shear(potential: sp.Expr) -> tuple[sp.Expr, ...]:
    gradient = sp.Matrix([sp.diff(potential, variable) for variable in VARIABLES])
    return tuple(
        sp.expand(VARIABLES[index] + (POISSON * gradient)[index])
        for index in range(4)
    )


def compose(
    outer: tuple[sp.Expr, ...],
    inner: tuple[sp.Expr, ...],
) -> tuple[sp.Expr, ...]:
    substitution = dict(zip(VARIABLES, inner, strict=True))
    return tuple(
        sp.expand(expression.subs(substitution, simultaneous=True))
        for expression in outer
    )


def verify_symplectic_map(outputs: tuple[sp.Expr, ...]) -> None:
    jacobian = sp.Matrix(outputs).jacobian(VARIABLES)
    assert sp.expand(jacobian.det()) == 1
    assert (jacobian * POISSON * jacobian.T).applyfunc(sp.expand) == POISSON


def moyal_defects(
    outputs: tuple[sp.Expr, ...],
    powers: tuple[int, ...] = (3, 5),
) -> dict[tuple[int, int, int, Exponent], sp.Rational]:
    defects: dict[tuple[int, int, int, Exponent], sp.Rational] = {}
    for power in powers:
        for left, right in PAIRS:
            value = pi_power(outputs[left], outputs[right], power)
            for monomial, coefficient in polynomial_coefficients(value).items():
                defects[(power, left, right, monomial)] = coefficient
    return defects


def pencil_supports() -> tuple[tuple[str, tuple[sp.Expr, ...]], ...]:
    cubic = monomials_of_degree(3)
    quartic = monomials_of_degree(4)
    return (
        (
            "3",
            tuple((monomial,) for monomial in cubic)
            + tuple(itertools.combinations(cubic, 2)),
        ),
        (
            "3+4",
            tuple((left, right) for left in cubic for right in quartic),
        ),
        (
            "4",
            tuple((monomial,) for monomial in quartic)
            + tuple(itertools.combinations(quartic, 2)),
        ),
    )


def potential_from_support(support: tuple[sp.Expr, ...]) -> sp.Expr:
    if len(support) == 1:
        return A_PARAMETER * support[0]
    return A_PARAMETER * support[0] + B_PARAMETER * support[1]


def search_pencils() -> tuple[tuple[PencilCount, ...], tuple[sp.Expr, ...]]:
    counts = []
    primitive_monomials: set[sp.Expr] = set()
    for degree_pattern, supports in pencil_supports():
        nilpotent_count = 0
        integrable_count = 0
        square_zero_count = 0
        moyal_flat_count = 0
        for support in supports:
            potential = potential_from_support(support)
            nilpotent = hamiltonian_matrix(potential)
            square = nilpotent * nilpotent
            if not is_zero_matrix(square * square):
                continue
            nilpotent_count += 1
            jacobian = cayley_jacobian(nilpotent)
            if not rows_are_closed(jacobian):
                continue
            integrable_count += 1
            outputs = integrate_closed_rows(jacobian)
            verify_symplectic_map(outputs)
            if is_zero_matrix(square):
                square_zero_count += 1
                assert outputs == euler_shear(potential)
            if not moyal_defects(outputs):
                moyal_flat_count += 1
            if len(support) == 1:
                primitive_monomials.add(support[0])
        counts.append(
            PencilCount(
                degree_pattern=degree_pattern,
                templates=len(supports),
                nilpotent=nilpotent_count,
                cayley_integrable=integrable_count,
                square_zero=square_zero_count,
                moyal_flat=moyal_flat_count,
            )
        )
    return tuple(counts), tuple(sorted(primitive_monomials, key=str))


def word_map(first: sp.Expr, second: sp.Expr) -> tuple[sp.Expr, ...]:
    return compose(
        euler_shear(B_PARAMETER * second),
        euler_shear(A_PARAMETER * first),
    )


def specialize_word(outputs: tuple[sp.Expr, ...]) -> tuple[sp.Expr, ...]:
    return tuple(
        sp.expand(expression.subs({A_PARAMETER: 1, B_PARAMETER: 1}))
        for expression in outputs
    )


def preliminary_words(
    primitive_monomials: tuple[sp.Expr, ...],
) -> tuple[tuple[PreliminaryWord, tuple[sp.Expr, ...]], ...]:
    rows = []
    word_index = 0
    for first, second in itertools.combinations(primitive_monomials, 2):
        if poisson(first, second) == 0:
            continue
        word_index += 1
        outputs = word_map(first, second)
        verify_symplectic_map(outputs)
        sample = specialize_word(outputs)
        defects = moyal_defects(sample)
        defect_terms = len(defects)
        nonzero_orders = sorted({key[0] for key in defects})
        maximum_degree = max(
            sp.Poly(expression, VARIABLES).total_degree() for expression in sample
        )
        map_terms = sum(
            len(polynomial_coefficients(expression)) for expression in sample
        )
        # The base is A^2, and these maps need no localization.  Complexity
        # is penalized mildly; nonzero higher Moyal support dominates.
        preliminary_score = 12 * defect_terms - map_terms - 2 * 2
        row = PreliminaryWord(
            key=f"HPW{word_index:03d}",
            first_potential=str(first),
            second_potential=str(second),
            maximum_degree=maximum_degree,
            map_terms=map_terms,
            moyal_defect_terms=defect_terms,
            first_nonzero_moyal_order=(nonzero_orders[0] if nonzero_orders else None),
            localization_pole_width=0,
            parameter_base_complexity=2,
            preliminary_score=preliminary_score,
        )
        rows.append((row, sample))
    return tuple(rows)


def tuple_as_column(values: tuple[sp.Expr, ...]) -> SparseColumn:
    column: SparseColumn = {}
    for block, value in enumerate(values):
        for monomial, coefficient in polynomial_coefficients(value).items():
            column[(block, monomial)] = coefficient
    return column


def columns_matrix(
    columns: list[SparseColumn],
    constant: SparseColumn | None = None,
) -> tuple[list[TaggedExponent], sp.Matrix, sp.Matrix]:
    row_keys: set[TaggedExponent] = set(constant or {})
    for column in columns:
        row_keys.update(column)
    rows = sorted(row_keys)
    row_index = {row: index for index, row in enumerate(rows)}
    matrix = sp.zeros(len(rows), len(columns))
    for column_index, column in enumerate(columns):
        for row, coefficient in column.items():
            matrix[row_index[row], column_index] = coefficient
    constant_vector = sp.Matrix(
        [(constant or {}).get(row, 0) for row in rows]
    )
    return rows, matrix, constant_vector


def native_correction_data(
    outputs: tuple[sp.Expr, ...],
) -> tuple[list[tuple[int, Exponent]], list[SparseColumn]]:
    correction_basis = [
        (component, monomial)
        for component, output in enumerate(outputs)
        for monomial in sorted(polynomial_coefficients(output))
    ]
    columns = []
    for component, exponent in correction_basis:
        monomial = sp.prod(
            variable**degree
            for variable, degree in zip(VARIABLES, exponent)
        )
        values = []
        for left, right in PAIRS:
            value = 0
            if component == left:
                value += poisson(monomial, outputs[right])
            if component == right:
                value += poisson(outputs[left], monomial)
            values.append(sp.expand(value))
        columns.append(tuple_as_column(tuple(values)))
    return correction_basis, columns


def unpack_correction(
    vector: sp.Matrix,
    correction_basis: list[tuple[int, Exponent]],
) -> tuple[sp.Expr, ...]:
    values = [sp.Integer(0)] * 4
    for coefficient, (component, exponent) in zip(vector, correction_basis):
        monomial = sp.prod(
            variable**degree
            for variable, degree in zip(VARIABLES, exponent)
        )
        values[component] += coefficient * monomial
    return tuple(sp.expand(value) for value in values)


def order_three_defect(outputs: tuple[sp.Expr, ...]) -> SparseColumn:
    return tuple_as_column(
        tuple(
            sp.expand(pi_power(outputs[left], outputs[right], 3) / 24)
            for left, right in PAIRS
        )
    )


def order_five_defect(
    outputs: tuple[sp.Expr, ...],
    lower: tuple[sp.Expr, ...],
) -> SparseColumn:
    values = []
    for left, right in PAIRS:
        values.append(
            sp.expand(
                poisson(lower[left], lower[right])
                + pi_power(lower[left], outputs[right], 3) / 24
                + pi_power(outputs[left], lower[right], 3) / 24
                + pi_power(outputs[left], outputs[right], 5) / 1920
            )
        )
    return tuple_as_column(tuple(values))


def moyal_derivation(
    series: dict[int, sp.Expr],
    potential: sp.Expr,
    maximum_hbar_order: int = 4,
) -> dict[int, sp.Expr]:
    """Apply ``{-,potential}_M`` to a truncated even-hbar series."""

    output: dict[int, sp.Expr] = {}
    for hbar_order, expression in series.items():
        for added_order, power, denominator in (
            (0, 1, 1),
            (2, 3, 24),
            (4, 5, 1920),
        ):
            target_order = hbar_order + added_order
            if target_order > maximum_hbar_order:
                continue
            value = pi_power(expression, potential, power) / denominator
            if value:
                output[target_order] = sp.expand(
                    output.get(target_order, 0) + value
                )
    return {order: value for order, value in output.items() if value}


def exponential_moyal_pullback(
    expression: sp.Expr,
    potential: sp.Expr,
    maximum_hbar_order: int = 4,
) -> dict[int, sp.Expr]:
    """Compute ``exp({-,potential}_M)(expression)`` exactly."""

    total = {0: expression}
    current = {0: expression}
    factorial = 1
    for iteration in range(1, 33):
        current = moyal_derivation(
            current,
            potential,
            maximum_hbar_order,
        )
        factorial *= iteration
        if not current:
            return {
                order: sp.expand(value)
                for order, value in total.items()
                if value
            }
        for order, value in current.items():
            total[order] = sp.expand(total.get(order, 0) + value / factorial)
    raise AssertionError("factorwise Moyal exponential did not terminate")


def factorwise_quantization(
    row: PreliminaryWord,
    outputs: tuple[sp.Expr, ...],
) -> dict[str, object]:
    """Replay the known two-shear Weyl lift through hbar^5."""

    local_symbols = {str(variable): variable for variable in VARIABLES}
    first = sp.sympify(row.first_potential, locals=local_symbols)
    second = sp.sympify(row.second_potential, locals=local_symbols)
    outer_outputs = euler_shear(second)
    series = tuple(
        exponential_moyal_pullback(expression, first)
        for expression in outer_outputs
    )
    classical = tuple(item.get(0, 0) for item in series)
    correction_two = tuple(item.get(2, 0) for item in series)
    correction_four = tuple(item.get(4, 0) for item in series)
    assert classical == outputs

    for left, right in PAIRS:
        order_three_residual = sp.expand(
            poisson(correction_two[left], outputs[right])
            + poisson(outputs[left], correction_two[right])
            + pi_power(outputs[left], outputs[right], 3) / 24
        )
        assert order_three_residual == 0
        order_five_residual = sp.expand(
            poisson(correction_four[left], outputs[right])
            + poisson(outputs[left], correction_four[right])
            + poisson(correction_two[left], correction_two[right])
            + pi_power(correction_two[left], outputs[right], 3) / 24
            + pi_power(outputs[left], correction_two[right], 3) / 24
            + pi_power(outputs[left], outputs[right], 5) / 1920
        )
        assert order_five_residual == 0

    native_support = [
        set(polynomial_coefficients(output)) for output in outputs
    ]

    def correction_record(correction: tuple[sp.Expr, ...]) -> dict[str, object]:
        term_count = 0
        outside_count = 0
        components = []
        for component, value in enumerate(correction):
            support = set(polynomial_coefficients(value))
            term_count += len(support)
            outside_count += len(support - native_support[component])
            components.append(str(sp.factor(value)))
        return {
            "components": components,
            "term_count": term_count,
            "terms_outside_native_component_support": outside_count,
        }

    return {
        "method": (
            "factorwise Moyal exponential for the two exact shear "
            "automorphisms, specialized at (a,b)=(1,1)"
        ),
        "hbar2_correction": correction_record(correction_two),
        "hbar4_correction": correction_record(correction_four),
        "verified_relations": "Moyal canonical relations through hbar^5",
        "all_order_status": (
            "known exact Weyl automorphism by factorwise composition"
        ),
    }


def one_particular_solution(matrix: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
    solution_set = sp.linsolve((matrix, right))
    solution = next(iter(solution_set))
    free_parameters = set().union(
        *(entry.free_symbols for entry in solution)
    ) - set(VARIABLES)
    return sp.Matrix(
        [
            entry.subs({parameter: 0 for parameter in free_parameters})
            for entry in solution
        ]
    )


def admissible_gauge_rank(
    outputs: tuple[sp.Expr, ...],
    correction_basis: list[tuple[int, Exponent]],
    maximum_hamiltonian_degree: int = 4,
) -> tuple[int, int]:
    """Rank of bounded Hamiltonian gauges preserving native support."""

    hamiltonians = tuple(
        exponent
        for degree in range(1, maximum_hamiltonian_degree + 1)
        for exponent in itertools.product(range(degree + 1), repeat=4)
        if sum(exponent) == degree
    )
    correction_index = {
        basis_element: index for index, basis_element in enumerate(correction_basis)
    }
    raw_columns = []
    outside_rows: set[tuple[int, Exponent]] = set()
    for exponent in hamiltonians:
        hamiltonian = sp.prod(
            variable**degree
            for variable, degree in zip(VARIABLES, exponent)
        )
        inside: dict[tuple[int, Exponent], sp.Expr] = {}
        outside: dict[tuple[int, Exponent], sp.Expr] = {}
        for component, output in enumerate(outputs):
            variation = poisson(hamiltonian, output)
            for monomial, coefficient in polynomial_coefficients(variation).items():
                key = (component, monomial)
                if key in correction_index:
                    inside[key] = coefficient
                else:
                    outside[key] = coefficient
                    outside_rows.add(key)
        raw_columns.append((inside, outside))

    outside_rows_sorted = sorted(outside_rows)
    outside_index = {
        row: index for index, row in enumerate(outside_rows_sorted)
    }
    outside_matrix = sp.zeros(len(outside_rows_sorted), len(hamiltonians))
    inside_matrix = sp.zeros(len(correction_basis), len(hamiltonians))
    for column_index, (inside, outside) in enumerate(raw_columns):
        for row, coefficient in inside.items():
            inside_matrix[correction_index[row], column_index] = coefficient
        for row, coefficient in outside.items():
            outside_matrix[outside_index[row], column_index] = coefficient
    nullspace = outside_matrix.nullspace()
    admissible = (
        sp.Matrix.hstack(*nullspace)
        if nullspace
        else sp.zeros(len(hamiltonians), 0)
    )
    gauge_image = inside_matrix * admissible
    return admissible.cols, gauge_image.rank()


def quantization_profile(
    row: PreliminaryWord,
    outputs: tuple[sp.Expr, ...],
) -> dict[str, object]:
    correction_basis, current_columns = native_correction_data(outputs)
    defect_three = order_three_defect(outputs)
    rows_three, matrix_three, constant_three = columns_matrix(
        current_columns,
        defect_three,
    )
    rank_three = matrix_three.rank()
    augmented_three = matrix_three.row_join(constant_three).rank()
    factorwise = factorwise_quantization(row, outputs)
    if augmented_three != rank_three:
        admissible_gauges, gauge_rank = admissible_gauge_rank(
            outputs,
            correction_basis,
        )
        cokernel_three = len(rows_three) - rank_three
        final_score = (
            int(100 * sp.Rational(cokernel_three, len(rows_three)))
            + 10 * row.moyal_defect_terms
            + 50 * (augmented_three - rank_three)
            - 5 * gauge_rank
            - 10 * row.localization_pole_width
            - 3 * row.parameter_base_complexity
        )
        return {
            **asdict(row),
            "native_order_three": {
                "correction_dimension": len(correction_basis),
                "defect_dimension": len(rows_three),
                "rank": rank_three,
                "cokernel_dimension": len(rows_three) - rank_three,
                "augmented_rank": augmented_three,
                "section_rank_jump": augmented_three - rank_three,
            },
            "native_order_five": None,
            "native_gauge": {
                "hamiltonian_degree_bound": 4,
                "admissible_hamiltonian_dimension": admissible_gauges,
                "image_rank": gauge_rank,
            },
            "score": final_score,
            "quantization_status": "native-support obstruction at hbar^3",
            "factorwise_quantization": factorwise,
        }

    particular = one_particular_solution(matrix_three, -constant_three)
    kernel = matrix_three.nullspace()
    kernel_matrix = (
        sp.Matrix.hstack(*kernel)
        if kernel
        else sp.zeros(matrix_three.cols, 0)
    )
    base_lower = unpack_correction(particular, correction_basis)
    kernel_lowers = [
        unpack_correction(kernel_matrix[:, index], correction_basis)
        for index in range(kernel_matrix.cols)
    ]

    constant_five = order_five_defect(outputs, base_lower)
    lower_columns: list[SparseColumn] = []
    for variation in kernel_lowers:
        diagonal = tuple_as_column(
            tuple(
                poisson(variation[left], variation[right])
                for left, right in PAIRS
            )
        )
        shifted = order_five_defect(
            outputs,
            tuple(
                base_lower[index] + variation[index]
                for index in range(4)
            ),
        )
        keys = set(shifted) | set(constant_five) | set(diagonal)
        linear = {
            key: shifted.get(key, 0)
            - constant_five.get(key, 0)
            - diagonal.get(key, 0)
            for key in keys
        }
        lower_columns.extend(
            (
                {key: value for key, value in linear.items() if value},
                diagonal,
            )
        )
    for first, second in itertools.combinations(range(len(kernel_lowers)), 2):
        left_variation = kernel_lowers[first]
        right_variation = kernel_lowers[second]
        lower_columns.append(
            tuple_as_column(
                tuple(
                    poisson(left_variation[left], right_variation[right])
                    + poisson(right_variation[left], left_variation[right])
                    for left, right in PAIRS
                )
            )
        )

    strong_columns = current_columns + lower_columns
    rows_five, matrix_five, constant_five_vector = columns_matrix(
        strong_columns,
        constant_five,
    )
    rank_five = matrix_five.rank()
    augmented_five = matrix_five.row_join(constant_five_vector).rank()
    admissible_gauges, gauge_rank = admissible_gauge_rank(
        outputs,
        correction_basis,
    )

    cokernel_three = len(rows_three) - rank_three
    cokernel_five = len(rows_five) - rank_five
    average_cokernel_ratio = (
        sp.Rational(cokernel_three, len(rows_three))
        + sp.Rational(cokernel_five, len(rows_five))
    ) / 2
    # Integer presentation score: favor a large relative cokernel and raw
    # Moyal support, penalize gauge, poles, and coefficient-base complexity.
    final_score = (
        int(100 * average_cokernel_ratio)
        + 10 * row.moyal_defect_terms
        + 50 * (augmented_five - rank_five)
        - 5 * gauge_rank
        - 10 * row.localization_pole_width
        - 3 * row.parameter_base_complexity
    )
    return {
        **asdict(row),
        "native_order_three": {
            "correction_dimension": len(correction_basis),
            "defect_dimension": len(rows_three),
            "rank": rank_three,
            "cokernel_dimension": cokernel_three,
            "lift_kernel_dimension": kernel_matrix.cols,
            "augmented_rank": augmented_three,
            "section_rank_jump": augmented_three - rank_three,
        },
        "native_order_five": {
            "current_correction_columns": len(current_columns),
            "lower_lift_span_columns": len(lower_columns),
            "strong_columns": len(strong_columns),
            "defect_dimension": len(rows_five),
            "rank": rank_five,
            "cokernel_dimension": cokernel_five,
            "augmented_rank": augmented_five,
            "section_rank_jump": augmented_five - rank_five,
        },
        "native_gauge": {
            "hamiltonian_degree_bound": 4,
            "admissible_hamiltonian_dimension": admissible_gauges,
            "image_rank": gauge_rank,
        },
        "score": final_score,
        "quantization_status": (
            "native-support lift through hbar^5"
            if augmented_five == rank_five
            else "native-support obstruction at hbar^5"
        ),
        "factorwise_quantization": factorwise,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--shortlist",
        type=int,
        default=8,
        help="number of preliminary two-pencil words to quantize",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="optional deterministic JSON certificate path",
    )
    args = parser.parse_args()
    if args.shortlist < 1:
        raise SystemExit("--shortlist must be positive")

    pencil_counts, primitive_monomials = search_pencils()
    assert pencil_counts == (
        PencilCount("3", 210, 56, 40, 40, 40),
        PencilCount("3+4", 700, 124, 84, 84, 84),
        PencilCount("4", 630, 84, 60, 60, 60),
    )
    assert len(primitive_monomials) == 28

    word_rows = preliminary_words(primitive_monomials)
    assert len(word_rows) == 238
    nonflat_rows = [item for item in word_rows if item[0].moyal_defect_terms]
    assert len(nonflat_rows) == 24
    ordered = sorted(
        nonflat_rows,
        key=lambda item: (
            item[0].preliminary_score,
            item[0].moyal_defect_terms,
            -item[0].map_terms,
            item[0].key,
        ),
        reverse=True,
    )
    shortlisted = ordered[: min(args.shortlist, len(ordered))]
    profiles = [quantization_profile(row, outputs) for row, outputs in shortlisted]
    profiles.sort(
        key=lambda profile: (profile.get("score", -10**9), profile["key"]),
        reverse=True,
    )

    certificate = {
        "format": "dc2-hessian-symbol-optimization-v1",
        "scope": (
            "exact coefficient-uniform cubic/quartic sparse-family search "
            "and native-support Moyal calibration at (a,b)=(1,1); "
            "exceptional coefficient subloci in rejected supports are not "
            "classified, and two-step rows are known polynomial "
            "automorphism controls, not DC_2 candidates"
        ),
        "bridge": {
            "variables": [str(variable) for variable in VARIABLES],
            "poisson_convention": "{p_i,q_j}=delta_ij",
            "nilpotent_matrix": "N=Pi*Hess(A)",
            "nilpotence_gate": "N^4=0",
            "cayley_jacobian": "I+N+N^2/2+N^3/4",
            "integrability_gate": "each Cayley row is a closed polynomial one-form",
        },
        "one_pencil_search": {
            "counts": [asdict(count) for count in pencil_counts],
            "total_templates": sum(count.templates for count in pencil_counts),
            "total_cayley_integrable": sum(
                count.cayley_integrable for count in pencil_counts
            ),
            "higher_nilpotence_integrable": 0,
            "support_policy": (
                "identities must hold over the full natural Q[a] or "
                "Q[a,b] coefficient base; exceptional subloci of rejected "
                "supports are outside this search"
            ),
            "conclusion": (
                "every Cayley-integrable template is square-zero and has "
                "zero odd Moyal defect through the degree cutoff"
            ),
        },
        "two_pencil_search": {
            "primitive_monomial_pencils": len(primitive_monomials),
            "canonical_noncommuting_words": len(word_rows),
            "words_with_nonzero_raw_moyal_defect": len(nonflat_rows),
            "shortlist_size": len(profiles),
            "orientation_rule": (
                "one order per unordered pair; reversal is represented by "
                "inversion together with parameter sign change"
            ),
        },
        "score": {
            "hard_gates": [
                "polynomial four-tuple",
                "exact symplectic identity",
                "Jacobian determinant one",
            ],
            "native_correction_lattice": (
                "componentwise monomial support at (a,b)=(1,1)"
            ),
            "higher_is_better": [
                "relative native obstruction-cokernel size",
                "raw odd Moyal defect support",
                "a nonzero native-support section rank jump",
            ],
            "lower_is_better": [
                "bounded Hamiltonian gauge-image rank",
                "localization pole width",
                "parameter-base complexity",
            ],
            "warning": (
                "the scalar score compares only rows in this common native "
                "lattice; use the component vector, not the scalar, across "
                "other DC_2 filtrations"
            ),
        },
        "quantized_shortlist": profiles,
        "selection_conclusion": (
            "the direct sparse Hessian bridge yields no new nonautomorphic "
            "symbol: one-pencil rows are Moyal-flat, while the hostile "
            "two-pencil rows are known factorwise-quantizable controls"
        ),
        "next_family_gate": (
            "seek a genuinely higher-nilpotence Cayley-integrable family or "
            "complete polynomial rank-two admission for the reciprocal R21 "
            "packet before spending on higher PBW orders"
        ),
    }

    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)
    else:
        print(rendered, end="")

    print(
        "PASS: 1,540 sparse Hessian pencils were screened exactly; "
        "184 are Cayley-integrable square-zero controls"
    )
    print(
        "PASS: 238 canonical noncommuting two-pencil words were scored; "
        "24 have nonzero raw Moyal defect"
    )
    print(
        f"PASS: quantized the top {len(profiles)} native-support rows "
        "through hbar^5"
    )
    print(
        "EXPERIMENT: the smallest sparse Hessian bridge supplies controls, "
        "not a new nonautomorphic DC_2 symbol"
    )


if __name__ == "__main__":
    main()
