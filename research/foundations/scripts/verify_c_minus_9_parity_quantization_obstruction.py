#!/usr/bin/env python3
"""Exact bounded A_2-quantization obstruction for the c=-9 completion.

The calculation takes the actual c=-9 pair (S,T) from
``QUADRATIC_LADDER_AND_POISSON_AUDIT.md``.  It works in the R-central Ore
algebra

    [Z,a] = hbar*delta(a),
    delta = 3*X^2*d_X + (2-6*X*Q)*d_Q.

Since delta commutes with d_Z, Weyl symbols use the abelian Moyal operator
Pi = d_Z tensor delta - delta tensor d_Z.

For the standard parity-preserving, order-lowering filtered ansatz

    S_h = Weyl(S + hbar^2*S2),
    T_h = Weyl(T + hbar^2*T2),

the natural complete correction spaces are

    ord_Z(S2) <= 1, deg_B(S2) <= 11,
    ord_Z(T2) <= 0, deg_B(T2) <= 7,

where deg_B(X)=deg_B(Q)=1 and deg_B(Z)=3.  This script proves over Q that:

* the hbar^3 equation has rank 149 and a 10-dimensional solution torsor;
* after retaining all ten parameters, the coefficient of X^12 in the
  hbar^5 commutator defect is always -49.

Thus coefficient extraction at X^12 is an explicit nonzero obstruction
functional in this bounded class.  No Groebner exhaustion is used.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from math import comb

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices.sdm import sdm_irref


Monomial = tuple[int, int, int]
SparsePoly = dict[Monomial, object]


def add_term(target: dict[Monomial, object], monomial: Monomial, value) -> None:
    if not value:
        return
    target[monomial] += value
    if not target[monomial]:
        del target[monomial]


def add(left: SparsePoly, right: SparsePoly, factor=1) -> SparsePoly:
    result = defaultdict(lambda: 0, left)
    for monomial, coefficient in right.items():
        add_term(result, monomial, factor * coefficient)
    return dict(result)


def scale(poly: SparsePoly, factor) -> SparsePoly:
    return {
        monomial: sp.expand(factor * coefficient)
        for monomial, coefficient in poly.items()
        if sp.expand(factor * coefficient) != 0
    }


def multiply(left: SparsePoly, right: SparsePoly) -> SparsePoly:
    result: dict[Monomial, object] = defaultdict(lambda: 0)
    for (i, j, k), a in left.items():
        for (ii, jj, kk), b in right.items():
            add_term(result, (i + ii, j + jj, k + kk), a * b)
    return dict(result)


def d_z(poly: SparsePoly, times: int = 1) -> SparsePoly:
    result = dict(poly)
    for _ in range(times):
        differentiated: dict[Monomial, object] = defaultdict(lambda: 0)
        for (i, j, k), coefficient in result.items():
            if k:
                add_term(differentiated, (i, j, k - 1), k * coefficient)
        result = dict(differentiated)
    return result


def delta(poly: SparsePoly, times: int = 1) -> SparsePoly:
    result = dict(poly)
    for _ in range(times):
        differentiated: dict[Monomial, object] = defaultdict(lambda: 0)
        for (i, j, k), coefficient in result.items():
            if i or j:
                add_term(
                    differentiated,
                    (i + 1, j, k),
                    (3 * i - 6 * j) * coefficient,
                )
            if j:
                add_term(differentiated, (i, j - 1, k), 2 * j * coefficient)
        result = dict(differentiated)
    return result


def poisson(left: SparsePoly, right: SparsePoly) -> SparsePoly:
    return add(
        multiply(d_z(left), delta(right)),
        multiply(delta(left), d_z(right)),
        -1,
    )


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
    result: SparsePoly = {}
    for monomial, coefficient in sp.Poly(sp.expand(expression), *variables).terms():
        value = sp.Rational(coefficient)
        result[monomial] = Fraction(int(value.p), int(value.q))
    return result


def filtered_monomials(max_degree: int, max_z_order: int) -> list[Monomial]:
    return [
        (i, j, k)
        for k in range(max_z_order + 1)
        for i in range(max_degree - 3 * k + 1)
        for j in range(max_degree - 3 * k - i + 1)
    ]


def qq(value):
    value = Fraction(value)
    return QQ(value.numerator, value.denominator)


def rational(value):
    return sp.Rational(int(value.numerator), int(value.denominator))


def c_minus_9_pair() -> tuple[SparsePoly, SparsePoly]:
    X, Q, Z = sp.symbols("X Q Z")
    W = Z - 9 * Q**2
    Y = Q - X * W / 3
    U = 1 + X * Y
    S = sp.expand((U**3 * W + Y**2 * U * (4 + 3 * X * Y)) / 2)
    T = sp.expand(Y + 3 * X * U**2 * W + 3 * X * Y**2 * (4 + 3 * X * Y))
    return sympy_poly_dict(S, (X, Q, Z)), sympy_poly_dict(T, (X, Q, Z))


def complete_hbar3_solution(
    S: SparsePoly, T: SparsePoly
) -> tuple[SparsePoly, SparsePoly, tuple[sp.Symbol, ...]]:
    s_monomials = filtered_monomials(11, 1)
    t_monomials = filtered_monomials(7, 0)
    columns = [
        poisson({monomial: Fraction(1)}, T) for monomial in s_monomials
    ]
    columns += [
        poisson(S, {monomial: Fraction(1)}) for monomial in t_monomials
    ]
    rhs = scale(pi_power(S, T, 3), Fraction(-1, 24))

    output_monomials = sorted(
        set(rhs).union(*(set(column) for column in columns))
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    rhs_column = len(columns)
    rows = {index: {} for index in range(len(output_monomials))}
    for column_index, column in enumerate(columns):
        for monomial, coefficient in column.items():
            rows[output_index[monomial]][column_index] = qq(coefficient)
    for monomial, coefficient in rhs.items():
        rows[output_index[monomial]][rhs_column] = -qq(coefficient)
    rows = {row: entries for row, entries in rows.items() if entries}

    reduced, pivots, _ = sdm_irref(rows)
    assert rhs_column not in pivots
    free_columns = [
        column for column in range(len(columns)) if column not in pivots
    ]
    assert len(columns) == 159
    assert len(pivots) == 149
    assert len(free_columns) == 10

    parameters = sp.symbols(f"a0:{len(free_columns)}")
    solution = {
        column: parameter
        for column, parameter in zip(free_columns, parameters, strict=True)
    }
    for reduced_row, pivot in enumerate(pivots):
        row = reduced.get(reduced_row, {})
        value = rational(row.get(rhs_column, QQ.zero))
        for free_column, parameter in zip(
            free_columns, parameters, strict=True
        ):
            coefficient = row.get(free_column, QQ.zero)
            if coefficient:
                value += rational(coefficient) * parameter
        solution[pivot] = -sp.expand(value)

    S2 = {
        monomial: solution[column]
        for column, monomial in enumerate(s_monomials)
        if solution.get(column, 0) != 0
    }
    offset = len(s_monomials)
    T2 = {
        monomial: solution[offset + column]
        for column, monomial in enumerate(t_monomials)
        if solution.get(offset + column, 0) != 0
    }
    reconstructed = add(poisson(S2, T), poisson(S, T2))
    assert reconstructed == rhs
    return S2, T2, parameters


def main() -> None:
    S, T = c_minus_9_pair()
    assert poisson(S, T) == {(0, 0, 0): Fraction(1)}
    assert max(i + j + 3 * k for i, j, k in S) == 15
    assert max(i + j + 3 * k for i, j, k in T) == 11
    assert max(k for _, _, k in S) == 3
    assert max(k for _, _, k in T) == 2

    S2, T2, parameters = complete_hbar3_solution(S, T)

    # Coefficient of hbar^5 in the unnormalized Weyl commutator, equivalently
    # hbar^4 in [S_h,T_h]/hbar.
    defect5 = poisson(S2, T2)
    defect5 = add(defect5, pi_power(S2, T, 3), sp.Rational(1, 24))
    defect5 = add(defect5, pi_power(S, T2, 3), sp.Rational(1, 24))
    defect5 = add(defect5, pi_power(S, T, 5), sp.Rational(1, 1920))

    obstruction = sp.expand(defect5[(12, 0, 0)])
    assert obstruction == -49
    assert not (set(parameters) & obstruction.free_symbols)

    # The free-zero representative displays the full compact residue.
    free_zero = {parameter: 0 for parameter in parameters}
    specialized = {
        monomial: sp.factor(coefficient.subs(free_zero))
        for monomial, coefficient in defect5.items()
        if coefficient.subs(free_zero) != 0
    }
    assert specialized == {
        (12, 0, 0): sp.Integer(-49),
        (13, 1, 0): sp.Integer(147),
        (14, 2, 0): sp.Rational(-441, 4),
    }

    print("PASS: the exact c=-9 Ore symbols satisfy {S,T}=1")
    print("PASS: the complete hbar^3 operator has rank 149 and nullity 10 over Q")
    print("PASS: all ten hbar^3 correction parameters were retained")
    print("PASS: [X^12] of the hbar^5 defect is identically -49")
    print(
        "RESIDUE (free-zero representative): "
        "-49*X^12*(1-3*X*Q/2)^2"
    )
    print(
        "CONCLUSION: the natural parity-preserving filtered lift of the "
        "c=-9 pair is obstructed"
    )
    print("SCOPE: odd-hbar and filtration-enlarging corrections are not excluded")


if __name__ == "__main__":
    main()
