#!/usr/bin/env python3
"""Exact Ore/Weyl obstruction at a degree-five rank-two sample.

Run with the repository Python::

    .venv/bin/python scripts/explore_degree_five_a2_subprincipal.py

The calculation specializes to the admissible seed (kappa,tau)=(0,1), whose
quadratic completing shear is 69/7.  In the R-central Ore algebra the fiber
Poisson bracket is

    {f,g} = f_Z*delta(g) - delta(f)*g_Z,
    delta = 3*X^2*d_X + (2-6*X*Q)*d_Q.

Fiber Weyl ordering removes the even normal-order defects.  Its first genuine
defect is Pi^3(S,T)/24 at hbar^3.  This script tests whether that defect lies
in the bounded filtered image

    {s_2,T} + {S,t_2},

where s_2 has differential order at most three and Bernstein degree at most
25, while t_2 has order at most two and degree at most 21.  It then retains
the complete 42-dimensional hbar^3 solution space and proves over Q that the
hbar^5 cokernel equations contain the constant equation 1=0.  Thus the
standard parity-preserving filtered Weyl-symbol ansatz is obstructed at this
sample.  Odd-in-hbar symbols and wider nonstandard filtrations are outside
the certificate's scope.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from math import comb
from random import Random

import sympy as sp
from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref


Monomial = tuple[int, int, int]
SparsePoly = dict[Monomial, Fraction]


def add_term(target: dict[Monomial, Fraction], monomial: Monomial, value) -> None:
    if not value:
        return
    target[monomial] += value
    if not target[monomial]:
        del target[monomial]


def add(left: SparsePoly, right: SparsePoly, scale=1) -> SparsePoly:
    result = defaultdict(lambda: 0, left)
    for monomial, coefficient in right.items():
        add_term(result, monomial, scale * coefficient)
    return dict(result)


def scale(poly: SparsePoly, coefficient) -> SparsePoly:
    return {
        monomial: coefficient * value
        for monomial, value in poly.items()
        if coefficient * value
    }


def multiply(left: SparsePoly, right: SparsePoly) -> SparsePoly:
    result: dict[Monomial, Fraction] = defaultdict(lambda: 0)
    for (i, j, k), a in left.items():
        for (ii, jj, kk), b in right.items():
            add_term(result, (i + ii, j + jj, k + kk), a * b)
    return dict(result)


def d_z(poly: SparsePoly, times: int = 1) -> SparsePoly:
    result = dict(poly)
    for _ in range(times):
        differentiated: dict[Monomial, Fraction] = defaultdict(lambda: 0)
        for (i, j, k), coefficient in result.items():
            if k:
                add_term(differentiated, (i, j, k - 1), coefficient * k)
        result = dict(differentiated)
    return result


def delta(poly: SparsePoly, times: int = 1) -> SparsePoly:
    result = dict(poly)
    for _ in range(times):
        differentiated: dict[Monomial, Fraction] = defaultdict(lambda: 0)
        for (i, j, k), coefficient in result.items():
            # 3*X^2*d_X - 6*X*Q*d_Q
            if i or j:
                add_term(
                    differentiated,
                    (i + 1, j, k),
                    coefficient * (3 * i - 6 * j),
                )
            # 2*d_Q
            if j:
                add_term(differentiated, (i, j - 1, k), coefficient * 2 * j)
        result = dict(differentiated)
    return result


def poisson(left: SparsePoly, right: SparsePoly) -> SparsePoly:
    return add(multiply(d_z(left), delta(right)), multiply(delta(left), d_z(right)), -1)


def pi_power(left: SparsePoly, right: SparsePoly, power: int) -> SparsePoly:
    result: SparsePoly = {}
    for j in range(power + 1):
        term = multiply(
            delta(d_z(left, power - j), j),
            d_z(delta(right, power - j), j),
        )
        result = add(result, term, (-1) ** j * comb(power, j))
    return result


def sympy_poly_dict(expression, variables) -> SparsePoly:
    poly = sp.Poly(sp.expand(expression), *variables)
    result: SparsePoly = {}
    for monomial, coefficient in poly.terms():
        rational = sp.Rational(coefficient)
        result[monomial] = Fraction(int(rational.p), int(rational.q))
    return result


def degree_five_sample() -> tuple[SparsePoly, SparsePoly]:
    w = sp.symbols("w")
    X, Q, Z = sp.symbols("X Q Z")
    a, tau, s2 = sp.symbols("a tau s2")

    kappa = -(1 + 2 * a) / (1 + a)
    H = sp.factor(
        w**2
        * (w - 1)
        * (
            tau * w**2
            + (kappa / 2 - 2 * tau + 2) * w
            - kappa / 2
            + tau
            - 3
        )
    )
    p = sp.diff(H, w)
    q = sp.expand(w * p - H)
    W = Z + s2 * Q**2
    Y = Q - X * W / 3
    source_v = -3 * X * Y / (2 * a)
    source_gamma = 1 - 3 * X * Q / 2
    source_u = 1 + source_v
    marked = sp.expand(source_u * source_gamma)
    S = sp.cancel(
        -2
        * a
        * ((source_u + q.subs(w, marked) / source_gamma**2) / X**2)
        / 3
    )
    T = sp.cancel((1 + p.subs(w, marked) / source_gamma) / X)
    shear = sp.factor(
        (
            -12 * a * (a + 1) ** 2 * tau**2
            + 18 * a * (a + 1) * (4 * a + 5) * tau
            + 216 * a**3
            + 648 * a**2
            + 738 * a
            + 315
        )
        / (28 * (a + 1) ** 2)
    )
    sample = {
        a: sp.Rational(-1, 2),
        tau: 1,
        s2: shear.subs({a: sp.Rational(-1, 2), tau: 1}),
    }
    assert sample[s2] == sp.Rational(69, 7)
    return (
        sympy_poly_dict(S.subs(sample), (X, Q, Z)),
        sympy_poly_dict(T.subs(sample), (X, Q, Z)),
    )


def filtered_monomials(max_degree: int, max_z_order: int) -> list[Monomial]:
    """Monomials with deg(X)=deg(Q)=1 and leading deg(Z)=3."""

    return [
        (i, j, k)
        for k in range(max_z_order + 1)
        for i in range(max_degree - 3 * k + 1)
        for j in range(max_degree - 3 * k - i + 1)
    ]


def solve_bounded(
    S: SparsePoly,
    T: SparsePoly,
    rhs: SparsePoly,
    s_monomials: list[Monomial],
    t_monomials: list[Monomial],
    label: str,
) -> tuple[SparsePoly, SparsePoly] | None:
    """Solve {s,T}+{S,t}=rhs in a declared finite filtered space."""

    columns: list[SparsePoly] = []
    for monomial in s_monomials:
        columns.append(poisson({monomial: Fraction(1)}, T))
    for monomial in t_monomials:
        columns.append(poisson(S, {monomial: Fraction(1)}))

    row_monomials = sorted(
        set(rhs).union(*(set(column) for column in columns))
    )
    row_index = {monomial: index for index, monomial in enumerate(row_monomials)}

    prime = 32003
    field = GF(prime)

    def reduce_mod_prime(value: Fraction):
        return field(value.numerator) / field(value.denominator)

    augmented_rows = {index: {} for index in range(len(row_monomials))}
    for column_index, column in enumerate(columns):
        for monomial, coefficient in column.items():
            augmented_rows[row_index[monomial]][column_index] = reduce_mod_prime(
                coefficient
            )
    rhs_column = len(columns)
    for monomial, coefficient in rhs.items():
        augmented_rows[row_index[monomial]][rhs_column] = reduce_mod_prime(-coefficient)
    augmented_rows = {
        row: entries for row, entries in augmented_rows.items() if entries
    }

    print(f"rows={len(row_monomials)} columns={len(columns)} prime={prime}")
    print(
        f"{label}: S monomials={len(s_monomials)} "
        f"T monomials={len(t_monomials)} rhs monomials={len(rhs)}"
    )
    reduced, pivots, _ = sdm_irref(augmented_rows)
    if rhs_column in pivots:
        print(f"OBSTRUCTION: {label} is outside the bounded image modulo {prime}")
        return None
    print(f"PASS: bounded {label} correction exists modulo {prime}")
    print(f"rank={len(pivots)} nullity={len(columns) - len(pivots)}")
    particular = {}
    for reduced_row, pivot in enumerate(pivots):
        if pivot == rhs_column:
            continue
        rhs_value = reduced.get(reduced_row, {}).get(rhs_column, field.zero)
        if rhs_value:
            # The augmented row represents M*x-rhs=0.
            particular[pivot] = -rhs_value
    print(f"free-zero modular solution support={len(particular)}")

    # Lift the sparse modular support over Q.
    selected_columns = sorted(particular)
    selected_index = {
        original_column: index
        for index, original_column in enumerate(selected_columns)
    }
    exact_rows = {}
    for row, monomial in enumerate(row_monomials):
        entries = {}
        for original_column in selected_columns:
            coefficient = columns[original_column].get(monomial)
            if coefficient:
                entries[selected_index[original_column]] = QQ(
                    coefficient.numerator, coefficient.denominator
                )
        rhs_value = rhs.get(monomial)
        if rhs_value:
            entries[len(selected_columns)] = QQ(
                -rhs_value.numerator, rhs_value.denominator
            )
        if entries:
            exact_rows[row] = entries

    exact_reduced, exact_pivots, _ = sdm_irref(exact_rows)
    exact_rhs_column = len(selected_columns)
    if exact_rhs_column in exact_pivots:
        print(f"SCOPE: the modular sparse support for {label} does not lift over Q")
        return None

    exact_solution: dict[int, Fraction] = {}
    for reduced_row, pivot in enumerate(exact_pivots):
        if pivot == exact_rhs_column:
            continue
        rhs_value = exact_reduced.get(reduced_row, {}).get(exact_rhs_column, QQ.zero)
        if rhs_value:
            exact_solution[selected_columns[pivot]] = Fraction(
                -int(rhs_value.numerator), int(rhs_value.denominator)
            )

    s_correction: SparsePoly = {}
    t_correction: SparsePoly = {}
    for column, coefficient in exact_solution.items():
        if column < len(s_monomials):
            s_correction[s_monomials[column]] = coefficient
        else:
            t_correction[t_monomials[column - len(s_monomials)]] = coefficient
    exact_image = add(poisson(s_correction, T), poisson(S, t_correction))
    assert exact_image == rhs
    print(f"PASS: sparse bounded {label} correction lifts exactly over Q")
    print(
        f"exact support: S2={len(s_correction)} T2={len(t_correction)} "
        f"total={len(exact_solution)}"
    )
    return s_correction, t_correction


def modular_gauge_search(S: SparsePoly, T: SparsePoly, trials: int = 64) -> bool:
    """Test whether hbar^3 gauge freedom can remove the hbar^5 obstruction."""

    prime = 32003
    field = GF(prime)

    def mod_value(value):
        value = Fraction(value)
        return field(value.numerator) / field(value.denominator)

    def mod_poly(poly):
        return {monomial: mod_value(value) for monomial, value in poly.items()}

    Sm, Tm = mod_poly(S), mod_poly(T)
    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    columns3 = [poisson({m: field.one}, Tm) for m in s2_monomials]
    columns3 += [poisson(Sm, {m: field.one}) for m in t2_monomials]
    rhs3 = scale(pi_power(Sm, Tm, 3), -field.one / field(24))
    rows3 = sorted(set(rhs3).union(*(set(column) for column in columns3)))
    index3 = {monomial: index for index, monomial in enumerate(rows3)}
    augmented3 = {index: {} for index in range(len(rows3))}
    for column_index, column in enumerate(columns3):
        for monomial, coefficient in column.items():
            augmented3[index3[monomial]][column_index] = coefficient
    rhs3_column = len(columns3)
    for monomial, coefficient in rhs3.items():
        augmented3[index3[monomial]][rhs3_column] = -coefficient
    augmented3 = {row: entries for row, entries in augmented3.items() if entries}
    reduced3, pivots3, _ = sdm_irref(augmented3)
    assert rhs3_column not in pivots3
    free3 = [column for column in range(len(columns3)) if column not in pivots3]
    assert len(free3) == 42

    s4_monomials = filtered_monomials(21, 1)
    t4_monomials = filtered_monomials(17, 0)
    columns5 = [poisson({m: field.one}, Tm) for m in s4_monomials]
    columns5 += [poisson(Sm, {m: field.one}) for m in t4_monomials]
    base_rows5: dict[Monomial, dict[int, object]] = defaultdict(dict)
    for column_index, column in enumerate(columns5):
        for monomial, coefficient in column.items():
            base_rows5[monomial][column_index] = coefficient
    rhs5_column = len(columns5)

    def hbar3_solution(free_values):
        solution = dict(free_values)
        for reduced_row, pivot in enumerate(pivots3):
            if pivot == rhs3_column:
                continue
            row = reduced3.get(reduced_row, {})
            value = row.get(rhs3_column, field.zero)
            for free_column, free_value in free_values.items():
                value += row.get(free_column, field.zero) * free_value
            solution[pivot] = -value
        s2 = {
            s2_monomials[column]: value
            for column, value in solution.items()
            if column < len(s2_monomials) and value
        }
        t2 = {
            t2_monomials[column - len(s2_monomials)]: value
            for column, value in solution.items()
            if column >= len(s2_monomials) and value
        }
        return s2, t2

    def hbar5_rhs(s2, t2):
        defect = poisson(s2, t2)
        defect = add(defect, pi_power(s2, Tm, 3), field.one / field(24))
        defect = add(defect, pi_power(Sm, t2, 3), field.one / field(24))
        defect = add(defect, pi_power(Sm, Tm, 5), field.one / field(1920))
        return scale(defect, -1)

    def is_consistent(rhs):
        row_monomials = sorted(set(base_rows5).union(rhs))
        augmented = {}
        for row_index, monomial in enumerate(row_monomials):
            entries = dict(base_rows5.get(monomial, {}))
            if monomial in rhs:
                entries[rhs5_column] = -rhs[monomial]
            if entries:
                augmented[row_index] = entries
        _, pivots, _ = sdm_irref(augmented)
        return rhs5_column not in pivots

    assignments = [{}]
    assignments.extend([{free: field.one} for free in free3])
    random = Random(20260722)
    for _ in range(trials):
        assignments.append(
            {free: field(random.randrange(prime)) for free in free3}
        )
    for index, assignment in enumerate(assignments):
        s2, t2 = hbar3_solution(assignment)
        if is_consistent(hbar5_rhs(s2, t2)):
            print(f"PASS: hbar^3 gauge choice {index} removes hbar^5 obstruction")
            return True
    print(
        f"SCOPE: no hbar^5-compatible gauge found among {len(assignments)} "
        f"choices modulo {prime}"
    )

    # Measure the finite obstruction space actually seen by the quadratic
    # gauge family.  This prepares an exact finite-field elimination rather
    # than treating the random search as a certificate.
    base_s2, base_t2 = hbar3_solution({})
    basis_pairs = []
    for free in free3:
        shifted_s2, shifted_t2 = hbar3_solution({free: field.one})
        basis_pairs.append(
            (add(shifted_s2, base_s2, -1), add(shifted_t2, base_t2, -1))
        )
    constant = hbar5_rhs(base_s2, base_t2)
    zero_parameter_monomial = (0,) * len(basis_pairs)
    coefficient_vectors = [constant]
    parameter_monomials = [zero_parameter_monomial]
    single_values = []
    diagonal_quadratics = []
    for basis_s2, basis_t2 in basis_pairs:
        diagonal = scale(poisson(basis_s2, basis_t2), -1)
        value = hbar5_rhs(add(base_s2, basis_s2), add(base_t2, basis_t2))
        linear = add(add(value, constant, -1), diagonal, -1)
        single_values.append(value)
        diagonal_quadratics.append(diagonal)
        coefficient_vectors.extend((linear, diagonal))
        linear_monomial = [0] * len(basis_pairs)
        linear_monomial[len(single_values) - 1] = 1
        quadratic_monomial = linear_monomial.copy()
        quadratic_monomial[len(single_values) - 1] = 2
        parameter_monomials.extend(
            (tuple(linear_monomial), tuple(quadratic_monomial))
        )
    for i in range(len(basis_pairs)):
        for j in range(i + 1, len(basis_pairs)):
            cross = add(
                poisson(basis_pairs[i][0], basis_pairs[j][1]),
                poisson(basis_pairs[j][0], basis_pairs[i][1]),
            )
            coefficient_vectors.append(scale(cross, -1))
            cross_monomial = [0] * len(basis_pairs)
            cross_monomial[i] = cross_monomial[j] = 1
            parameter_monomials.append(tuple(cross_monomial))

    combined_rows: dict[Monomial, dict[int, object]] = defaultdict(dict)
    for column_index, column in enumerate(columns5):
        for monomial, coefficient in column.items():
            combined_rows[monomial][column_index] = coefficient
    offset = len(columns5)
    for coefficient_index, coefficient_vector in enumerate(coefficient_vectors):
        for monomial, coefficient in coefficient_vector.items():
            combined_rows[monomial][offset + coefficient_index] = coefficient
    indexed_combined = {
        row: combined_rows[monomial]
        for row, monomial in enumerate(sorted(combined_rows))
    }
    _, combined_pivots, _ = sdm_irref(indexed_combined)
    operator_rank = sum(pivot < offset for pivot in combined_pivots)
    obstruction_rank = len(combined_pivots) - operator_rank
    print(
        f"quadratic gauge coefficient vectors={len(coefficient_vectors)} "
        f"cokernel span rank={obstruction_rank}"
    )

    # Project the quadratic family to a basis of the left cokernel.  Row
    # reduction then leaves precisely the independent obstruction equations
    # in the 42 gauge parameters.
    output_monomials = sorted(
        set(base_rows5).union(*(set(vector_) for vector_ in coefficient_vectors))
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    transpose_rows = {}
    for column_index, column in enumerate(columns5):
        row = {
            output_index[monomial]: coefficient
            for monomial, coefficient in column.items()
        }
        if row:
            transpose_rows[column_index] = row
    transpose_rref, transpose_pivots, transpose_nonzero = sdm_irref(
        transpose_rows
    )
    left_kernel, _ = sdm_nullspace_from_rref(
        transpose_rref,
        field.one,
        len(output_monomials),
        transpose_pivots,
        transpose_nonzero,
    )
    incidence: dict[int, list[tuple[int, object]]] = defaultdict(list)
    for functional_index, functional in enumerate(left_kernel):
        for output_row, coefficient in functional.items():
            incidence[output_row].append((functional_index, coefficient))
    projected = [defaultdict(lambda: field.zero) for _ in left_kernel]
    for coefficient_index, coefficient_vector in enumerate(coefficient_vectors):
        for monomial, value in coefficient_vector.items():
            for functional_index, functional_value in incidence[
                output_index[monomial]
            ]:
                projected[functional_index][coefficient_index] += (
                    functional_value * value
                )
    projected_rows = {
        row: {column: value for column, value in equation.items() if value}
        for row, equation in enumerate(projected)
        if any(equation.values())
    }
    obstruction_rref, obstruction_pivots, _ = sdm_irref(projected_rows)
    assert len(obstruction_pivots) == obstruction_rank
    term_counts = sorted(len(equation) for equation in obstruction_rref.values())
    print(
        f"independent quadratic obstruction equations={len(obstruction_pivots)} "
        f"term-count range={term_counts[0]}..{term_counts[-1]}"
    )
    singleton_monomials = []
    for equation in obstruction_rref.values():
        if len(equation) == 1:
            coefficient_index = next(iter(equation))
            singleton_monomials.append(parameter_monomials[coefficient_index])
    print(
        "singleton obstruction monomial degrees="
        + str([sum(monomial) for monomial in singleton_monomials])
    )
    if zero_parameter_monomial in singleton_monomials:
        print(
            f"OBSTRUCTION CERTIFICATE: the hbar^5 ideal contains 1 modulo {prime}"
        )
    return False


def exact_gauge_obstruction(S: SparsePoly, T: SparsePoly) -> bool:
    """Prove or disprove gauge-independent hbar^5 obstruction over Q."""

    def qq(value):
        value = Fraction(value)
        return QQ(value.numerator, value.denominator)

    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    columns3 = [poisson({m: Fraction(1)}, T) for m in s2_monomials]
    columns3 += [poisson(S, {m: Fraction(1)}) for m in t2_monomials]
    rhs3 = scale(pi_power(S, T, 3), Fraction(-1, 24))
    rows3 = sorted(set(rhs3).union(*(set(column) for column in columns3)))
    index3 = {monomial: index for index, monomial in enumerate(rows3)}
    augmented3 = {index: {} for index in range(len(rows3))}
    for column_index, column in enumerate(columns3):
        for monomial, coefficient in column.items():
            augmented3[index3[monomial]][column_index] = qq(coefficient)
    rhs3_column = len(columns3)
    for monomial, coefficient in rhs3.items():
        augmented3[index3[monomial]][rhs3_column] = -qq(coefficient)
    augmented3 = {row: entries for row, entries in augmented3.items() if entries}
    reduced3, pivots3, _ = sdm_irref(augmented3)
    assert rhs3_column not in pivots3
    free3 = [column for column in range(len(columns3)) if column not in pivots3]
    assert len(free3) == 42

    def hbar3_solution(free_values):
        solution = {column: qq(value) for column, value in free_values.items()}
        for reduced_row, pivot in enumerate(pivots3):
            if pivot == rhs3_column:
                continue
            row = reduced3.get(reduced_row, {})
            value = row.get(rhs3_column, QQ.zero)
            for free_column, free_value in solution.items():
                if free_column in free3:
                    value += row.get(free_column, QQ.zero) * free_value
            solution[pivot] = -value

        def fraction(value):
            return Fraction(int(value.numerator), int(value.denominator))

        s2 = {
            s2_monomials[column]: fraction(value)
            for column, value in solution.items()
            if column < len(s2_monomials) and value
        }
        t2 = {
            t2_monomials[column - len(s2_monomials)]: fraction(value)
            for column, value in solution.items()
            if column >= len(s2_monomials) and value
        }
        return s2, t2

    base_s2, base_t2 = hbar3_solution({})
    basis_pairs = []
    for free in free3:
        shifted_s2, shifted_t2 = hbar3_solution({free: 1})
        basis_pairs.append(
            (add(shifted_s2, base_s2, -1), add(shifted_t2, base_t2, -1))
        )

    def hbar5_rhs(s2, t2):
        defect = poisson(s2, t2)
        defect = add(defect, pi_power(s2, T, 3), Fraction(1, 24))
        defect = add(defect, pi_power(S, t2, 3), Fraction(1, 24))
        defect = add(defect, pi_power(S, T, 5), Fraction(1, 1920))
        return scale(defect, -1)

    constant = hbar5_rhs(base_s2, base_t2)
    coefficient_vectors = [constant]
    parameter_monomials = [(0,) * len(basis_pairs)]
    for i, (basis_s2, basis_t2) in enumerate(basis_pairs):
        diagonal = scale(poisson(basis_s2, basis_t2), -1)
        value = hbar5_rhs(add(base_s2, basis_s2), add(base_t2, basis_t2))
        linear = add(add(value, constant, -1), diagonal, -1)
        coefficient_vectors.extend((linear, diagonal))
        linear_monomial = [0] * len(basis_pairs)
        linear_monomial[i] = 1
        diagonal_monomial = linear_monomial.copy()
        diagonal_monomial[i] = 2
        parameter_monomials.extend(
            (tuple(linear_monomial), tuple(diagonal_monomial))
        )
    for i in range(len(basis_pairs)):
        for j in range(i + 1, len(basis_pairs)):
            cross = add(
                poisson(basis_pairs[i][0], basis_pairs[j][1]),
                poisson(basis_pairs[j][0], basis_pairs[i][1]),
            )
            coefficient_vectors.append(scale(cross, -1))
            cross_monomial = [0] * len(basis_pairs)
            cross_monomial[i] = cross_monomial[j] = 1
            parameter_monomials.append(tuple(cross_monomial))

    s4_monomials = filtered_monomials(21, 1)
    t4_monomials = filtered_monomials(17, 0)
    columns5 = [poisson({m: Fraction(1)}, T) for m in s4_monomials]
    columns5 += [poisson(S, {m: Fraction(1)}) for m in t4_monomials]
    output_monomials = sorted(
        set().union(
            *(set(column) for column in columns5),
            *(set(vector_) for vector_ in coefficient_vectors),
        )
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    transpose_rows = {}
    for column_index, column in enumerate(columns5):
        row = {
            output_index[monomial]: qq(coefficient)
            for monomial, coefficient in column.items()
        }
        if row:
            transpose_rows[column_index] = row
    transpose_rref, transpose_pivots, transpose_nonzero = sdm_irref(
        transpose_rows
    )
    left_kernel, _ = sdm_nullspace_from_rref(
        transpose_rref,
        QQ.one,
        len(output_monomials),
        transpose_pivots,
        transpose_nonzero,
    )
    incidence: dict[int, list[tuple[int, object]]] = defaultdict(list)
    for functional_index, functional in enumerate(left_kernel):
        for output_row, coefficient in functional.items():
            incidence[output_row].append((functional_index, coefficient))
    projected = [defaultdict(lambda: QQ.zero) for _ in left_kernel]
    for coefficient_index, coefficient_vector in enumerate(coefficient_vectors):
        for monomial, value in coefficient_vector.items():
            for functional_index, functional_value in incidence[
                output_index[monomial]
            ]:
                projected[functional_index][coefficient_index] += (
                    functional_value * qq(value)
                )
    projected_rows = {
        row: {column: value for column, value in equation.items() if value}
        for row, equation in enumerate(projected)
        if any(equation.values())
    }
    obstruction_rref, obstruction_pivots, _ = sdm_irref(projected_rows)
    singleton_monomials = []
    for equation in obstruction_rref.values():
        if len(equation) == 1:
            coefficient_index = next(iter(equation))
            singleton_monomials.append(parameter_monomials[coefficient_index])
    constant_monomial = (0,) * len(basis_pairs)
    if constant_monomial in singleton_monomials:
        print(
            "PASS: exact hbar^5 cokernel projection contains the constant "
            "equation 1=0 over Q"
        )
        print(
            f"exact obstruction rank={len(obstruction_pivots)} "
            f"singleton degrees={[sum(m) for m in singleton_monomials]}"
        )
        return True
    print("SCOPE: no exact gauge-independent hbar^5 constant was found")
    return False


def audit_odd_subprincipal_kernel(S: SparsePoly, T: SparsePoly) -> None:
    """Record the larger branch excluded by parity-preserving symbols."""

    s1_monomials = filtered_monomials(27, 4)
    t1_monomials = filtered_monomials(23, 3)
    columns = [poisson({m: Fraction(1)}, T) for m in s1_monomials]
    columns += [poisson(S, {m: Fraction(1)}) for m in t1_monomials]
    row_monomials = sorted(set().union(*(set(column) for column in columns)))
    row_index = {monomial: index for index, monomial in enumerate(row_monomials)}
    rows = {index: {} for index in range(len(row_monomials))}
    for column_index, column in enumerate(columns):
        for monomial, coefficient in column.items():
            rows[row_index[monomial]][column_index] = QQ(
                coefficient.numerator, coefficient.denominator
            )
    rows = {row: entries for row, entries in rows.items() if entries}
    _, pivots, _ = sdm_irref(rows)
    assert len(columns) == 2132
    assert len(pivots) == 2075
    print(
        "SCOPE: allowing odd hbar symbols introduces an exact "
        "57-dimensional first-order kernel"
    )


def main() -> None:
    S, T = degree_five_sample()
    assert poisson(S, T) == {(0, 0, 0): Fraction(1)}
    assert pi_power(S, T, 1) == poisson(S, T)
    assert max(i + j + 3 * k for i, j, k in S) == 29
    assert max(i + j + 3 * k for i, j, k in T) == 25

    # The hbar^3 coefficient of the fiber Moyal commutator.
    rhs3 = scale(pi_power(S, T, 3), Fraction(-1, 24))
    correction2 = solve_bounded(
        S,
        T,
        rhs3,
        filtered_monomials(25, 3),
        filtered_monomials(21, 2),
        "hbar^3",
    )
    if correction2 is None:
        return
    S2, T2 = correction2

    # At hbar^5 the already chosen hbar^2 symbols contribute through both
    # the Poisson and cubic Moyal terms.
    defect5 = poisson(S2, T2)
    defect5 = add(defect5, pi_power(S2, T, 3), Fraction(1, 24))
    defect5 = add(defect5, pi_power(S, T2, 3), Fraction(1, 24))
    defect5 = add(defect5, pi_power(S, T, 5), Fraction(1, 1920))
    correction4 = solve_bounded(
        S,
        T,
        scale(defect5, -1),
        filtered_monomials(21, 1),
        filtered_monomials(17, 0),
        "hbar^5",
    )
    if correction4 is None:
        for relaxation in (2, 4, 6, 10, 20):
            print(f"testing hbar^5 degree relaxation +{relaxation}")
            correction4 = solve_bounded(
                S,
                T,
                scale(defect5, -1),
                filtered_monomials(21 + relaxation, 1),
                filtered_monomials(17 + relaxation, 0),
                f"hbar^5 relaxed +{relaxation}",
            )
            if correction4 is not None:
                break
    if correction4 is None:
        modular_gauge_search(S, T)
        assert exact_gauge_obstruction(S, T)
        audit_odd_subprincipal_kernel(S, T)
        return
    S4, T4 = correction4

    # There is no room for hbar^6 symbols of negative differential order.
    # Hence the remaining hbar^7 and hbar^9 coefficients must vanish exactly.
    defect7: SparsePoly = {}
    defect7 = add(defect7, pi_power(S, T, 7), Fraction(1, 322560))
    defect7 = add(defect7, pi_power(S2, T, 5), Fraction(1, 1920))
    defect7 = add(defect7, pi_power(S, T2, 5), Fraction(1, 1920))
    defect7 = add(defect7, pi_power(S4, T, 3), Fraction(1, 24))
    defect7 = add(defect7, pi_power(S2, T2, 3), Fraction(1, 24))
    defect7 = add(defect7, pi_power(S, T4, 3), Fraction(1, 24))
    defect7 = add(defect7, poisson(S2, T4))
    defect7 = add(defect7, poisson(S4, T2))

    defect9: SparsePoly = {}
    defect9 = add(defect9, pi_power(S, T, 9), Fraction(1, 92897280))
    defect9 = add(defect9, pi_power(S2, T, 7), Fraction(1, 322560))
    defect9 = add(defect9, pi_power(S, T2, 7), Fraction(1, 322560))
    defect9 = add(defect9, pi_power(S4, T, 5), Fraction(1, 1920))
    defect9 = add(defect9, pi_power(S2, T2, 5), Fraction(1, 1920))
    defect9 = add(defect9, pi_power(S, T4, 5), Fraction(1, 1920))
    defect9 = add(defect9, pi_power(S4, T2, 3), Fraction(1, 24))
    defect9 = add(defect9, pi_power(S2, T4, 3), Fraction(1, 24))
    defect9 = add(defect9, poisson(S4, T4))
    print(f"remaining hbar^7 monomials={len(defect7)}")
    print(f"remaining hbar^9 monomials={len(defect9)}")


if __name__ == "__main__":
    main()
