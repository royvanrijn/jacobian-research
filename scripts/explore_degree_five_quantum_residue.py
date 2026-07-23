#!/usr/bin/env python3
"""Audit the proposed quantum-residue obstruction on the degree-five family.

There are three different coefficient categories in the rank-two problem:

``poly``
    The bounded polynomial symbol spaces used by the filtered A_2 test.

``X-Laurent(N)``
    The same symbol bounds after allowing a pole of order at most ``N`` at
    ``X=0``.  A monomial ``X^i Q^j Z^k`` is retained when ``i >= -N`` and
    ``i+j+w_Z*k`` satisfies the chart's Bernstein bound
    (``w_Z=3`` generically and ``w_Z=4`` at ``kappa=-1``).

``formal-etale``
    A completion at a branch in the etale coordinates ``(S,T)``.  The
    formal Poincare lemma contracts the correction complex here, but its
    primitives are generally infinite series and need not be X-Laurent.

The last distinction is essential: formal-etale solvability does *not*
imply X-Laurent solvability.  This script makes the finite comparison which
is actually justified.  It:

* constructs ``S,T`` directly over an arbitrary exact coefficient field for
  the full degree-five seed surface, including the ``kappa=-1`` replacement
  chart;
* computes the affine hbar^3 correction space;
* computes the hbar^5 obstruction family over every hbar^3 lift, not merely
  one selected lift; and
* repeats the obstruction-span test in finite X-pole bands.

The hbar^5 span test is a rigorous *nonexistence* certificate when it says
``constant outside span``.  If the constant enters the linear span, the
test is inconclusive because the coefficients of the span vectors must also
come from one common point of the quadratic parameter family.

Heavy calculations default to the exact finite field GF(32003).  The family
formulas themselves are field-independent.  The existing
``explore_degree_five_a2_subprincipal.py`` supplies the independent exact-Q
certificate at ``(kappa,tau)=(0,1)``.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations_with_replacement
from math import comb
from typing import Callable

from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref

from explore_degree_five_a2_subprincipal import (
    add,
    d_z,
    multiply,
    pi_power,
    poisson,
    scale,
)


Monomial = tuple[int, int, int]
SparsePoly = dict[Monomial, object]
PRIME = 32003


def field_fraction(field, numerator: int, denominator: int = 1):
    """Convert a rational integer pair to an element of ``field``."""

    return field(numerator) / field(denominator)


def poly_power(poly: SparsePoly, exponent: int, one) -> SparsePoly:
    result = {(0, 0, 0): one}
    for _ in range(exponent):
        result = multiply(result, poly)
    return result


def evaluate_univariate(
    coefficients: list[object], argument: SparsePoly, one
) -> SparsePoly:
    """Evaluate low-degree coefficients by Horner's rule."""

    result: SparsePoly = {}
    scalar_one = {(0, 0, 0): one}
    for coefficient in reversed(coefficients):
        result = add(multiply(result, argument), scale(scalar_one, coefficient))
    return result


def degree_five_family(
    field,
    a,
    tau,
    shear=None,
    *,
    verify_canonical: bool = True,
) -> tuple[SparsePoly, SparsePoly]:
    """Return the completed generic-chart fiber symbols ``(S,T)``.

    The normalized parameter is

        kappa = -(1+2*a)/(1+a).

    ``shear`` defaults to the unique classical completing value.  The
    construction avoids symbolic cancellation: the double zero of H at the
    origin is divided out before substitution, so all cancellations happen
    coefficientwise in the sparse dictionary.
    """

    one = field.one
    c = lambda numerator, denominator=1: field_fraction(
        field, numerator, denominator
    )
    scalar_one = {(0, 0, 0): one}
    X = {(1, 0, 0): one}
    Q = {(0, 1, 0): one}
    Z = {(0, 0, 1): one}

    if not a or not a + one:
        raise ValueError("generic a-chart requires a != 0,-1")

    if shear is None:
        shear = (
            -12 * a * (a + one) ** 2 * tau**2
            + 18 * a * (a + one) * (4 * a + 5) * tau
            + 216 * a**3
            + 648 * a**2
            + 738 * a
            + 315
        ) / (28 * (a + one) ** 2)

    W = add(Z, scale(poly_power(Q, 2, one), shear))
    Y = add(Q, scale(multiply(X, W), -c(1, 3)))
    source_u = add(
        scalar_one,
        scale(multiply(X, Y), -c(3, 2) / a),
    )
    source_gamma = add(
        scalar_one,
        scale(multiply(X, Q), -c(3, 2)),
    )
    marked = multiply(source_u, source_gamma)

    # H=w^2(w-1)(tau*w^2+A*w+B), with kappa eliminated in favour of a.
    kappa = -(one + 2 * a) / (one + a)
    linear = kappa / 2 - 2 * tau + 2
    constant = -kappa / 2 + tau - 3
    h_coefficients = [field.zero for _ in range(6)]
    for index, coefficient in enumerate((constant, linear, tau)):
        h_coefficients[index + 3] += coefficient
        h_coefficients[index + 2] -= coefficient

    # p(w)=H'(w)=w*p1(w), q(w)=w*H'(w)-H(w)=w^2*q2(w).
    p1_coefficients = [
        (index + 2) * h_coefficients[index + 2] for index in range(4)
    ]
    q2_coefficients = [
        (index + 1) * h_coefficients[index + 2] for index in range(4)
    ]
    p1_marked = evaluate_univariate(p1_coefficients, marked, one)
    q2_marked = evaluate_univariate(q2_coefficients, marked, one)

    T_numerator = add(scalar_one, multiply(source_u, p1_marked))
    T = {
        (x_degree - 1, q_degree, z_degree): coefficient
        for (x_degree, q_degree, z_degree), coefficient in T_numerator.items()
        if coefficient
    }
    S_numerator = scale(
        add(
            source_u,
            multiply(multiply(source_u, source_u), q2_marked),
        ),
        -2 * a / 3,
    )
    S = {
        (x_degree - 2, q_degree, z_degree): coefficient
        for (x_degree, q_degree, z_degree), coefficient in S_numerator.items()
        if coefficient
    }

    if min(x_degree for x_degree, _, _ in S) < 0:
        raise AssertionError("the S numerator did not cancel its X^2 factor")
    if min(x_degree for x_degree, _, _ in T) < 0:
        raise AssertionError("the T numerator did not cancel its X factor")
    if verify_canonical:
        bracket = poisson(S, T)
        if bracket != {(0, 0, 0): one}:
            raise AssertionError(
                "the generated family is not canonical: "
                f"{len(bracket)} nonzero bracket terms"
            )
    return S, T


def exceptional_delta(poly: SparsePoly, times: int = 1) -> SparsePoly:
    """Ore derivation on the replacement chart ``kappa=-1``."""

    result = dict(poly)
    for _ in range(times):
        differentiated: dict[Monomial, object] = defaultdict(lambda: 0)
        for (x_degree, q_degree, z_degree), coefficient in result.items():
            # -2*X^3*d_X + 6*X^2*Q*d_Q
            if x_degree or q_degree:
                differentiated[(x_degree + 2, q_degree, z_degree)] += (
                    coefficient * (-2 * x_degree + 6 * q_degree)
                )
            # 2*d_Q
            if q_degree:
                differentiated[(x_degree, q_degree - 1, z_degree)] += (
                    coefficient * 2 * q_degree
                )
        result = {
            monomial: value
            for monomial, value in differentiated.items()
            if value
        }
    return result


def exceptional_poisson(left: SparsePoly, right: SparsePoly) -> SparsePoly:
    return add(
        multiply(d_z(left), exceptional_delta(right)),
        multiply(exceptional_delta(left), d_z(right)),
        -1,
    )


def exceptional_pi_power(
    left: SparsePoly, right: SparsePoly, power: int
) -> SparsePoly:
    result: SparsePoly = {}
    for index in range(power + 1):
        term = multiply(
            exceptional_delta(d_z(left, power - index), index),
            d_z(exceptional_delta(right, power - index), index),
        )
        result = add(result, term, (-1) ** index * comb(power, index))
    return result


def degree_five_exceptional_family(
    field,
    tau,
    shear=None,
    *,
    verify_canonical: bool = True,
) -> tuple[SparsePoly, SparsePoly]:
    """Return ``(S,T)`` on the replacement chart ``kappa=-1``."""

    one = field.one
    c = lambda numerator, denominator=1: field_fraction(
        field, numerator, denominator
    )
    scalar_one = {(0, 0, 0): one}
    X = {(1, 0, 0): one}
    Q = {(0, 1, 0): one}
    Z = {(0, 0, 1): one}
    if shear is None:
        shear = 2 * (2 * tau**2 - 15 * tau - 18) / 105

    W = add(Z, scale(multiply(X, Q), shear))
    source_u = add(scalar_one, multiply(X, W))
    source_gamma = add(
        scalar_one,
        multiply(poly_power(X, 2, one), Q),
    )
    marked = multiply(source_u, source_gamma)

    linear = c(3, 2) - 2 * tau
    constant = tau - c(5, 2)
    h_coefficients = [field.zero for _ in range(6)]
    for index, coefficient in enumerate((constant, linear, tau)):
        h_coefficients[index + 3] += coefficient
        h_coefficients[index + 2] -= coefficient
    p1_coefficients = [
        (index + 2) * h_coefficients[index + 2] for index in range(4)
    ]
    q2_coefficients = [
        (index + 1) * h_coefficients[index + 2] for index in range(4)
    ]
    p1_marked = evaluate_univariate(p1_coefficients, marked, one)
    q2_marked = evaluate_univariate(q2_coefficients, marked, one)

    T_numerator = add(scalar_one, multiply(source_u, p1_marked))
    T = {
        (x_degree - 1, q_degree, z_degree): coefficient
        for (x_degree, q_degree, z_degree), coefficient in T_numerator.items()
        if coefficient
    }
    S_numerator = scale(
        add(
            source_u,
            multiply(multiply(source_u, source_u), q2_marked),
        ),
        c(1, 2),
    )
    S = {
        (x_degree - 2, q_degree, z_degree): coefficient
        for (x_degree, q_degree, z_degree), coefficient in S_numerator.items()
        if coefficient
    }
    if min(x_degree for x_degree, _, _ in S) < 0:
        raise AssertionError("exceptional S did not cancel its X^2 factor")
    if min(x_degree for x_degree, _, _ in T) < 0:
        raise AssertionError("exceptional T did not cancel its X factor")
    if (
        verify_canonical
        and exceptional_poisson(S, T) != {(0, 0, 0): one}
    ):
        raise AssertionError("the exceptional family is not canonical")
    return S, T


@dataclass(frozen=True)
class CorrectionProfile:
    poisson: Callable[[SparsePoly, SparsePoly], SparsePoly]
    pi_power: Callable[[SparsePoly, SparsePoly, int], SparsePoly]
    z_weight: int
    s2_degree: int
    t2_degree: int
    s4_degree: int
    t4_degree: int
    s1_degree: int
    t1_degree: int
    r_degree: int
    R: tuple[tuple[Monomial, int], ...]


GENERIC_PROFILE = CorrectionProfile(
    poisson,
    pi_power,
    3,
    25,
    21,
    21,
    17,
    27,
    23,
    3,
    (((1, 0, 0), 2), ((2, 1, 0), -3)),
)
EXCEPTIONAL_PROFILE = CorrectionProfile(
    exceptional_poisson,
    exceptional_pi_power,
    4,
    28,
    24,
    24,
    20,
    30,
    26,
    4,
    (((1, 0, 0), 2), ((3, 1, 0), 2)),
)


def laurent_monomials(
    max_degree: int,
    max_z_order: int,
    pole_order: int,
    z_weight: int = 3,
) -> list[Monomial]:
    """Finite numerator model for a pole of order at most ``pole_order``."""

    return [
        (x_degree, q_degree, z_degree)
        for z_degree in range(max_z_order + 1)
        for x_degree in range(
            -pole_order,
            max_degree - z_weight * z_degree + 1,
        )
        for q_degree in range(
            max_degree - z_weight * z_degree - x_degree + 1
        )
    ]


def solve_affine(
    columns: list[SparsePoly],
    rhs: SparsePoly,
    field,
) -> tuple[dict[int, object], list[dict[int, object]], int]:
    """Solve ``columns*x=rhs`` and return a point and a kernel basis."""

    output_monomials = sorted(
        set(rhs).union(*(set(column) for column in columns))
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    rhs_column = len(columns)
    rows: dict[int, dict[int, object]] = {}
    for column_index, column in enumerate(columns):
        for monomial, coefficient in column.items():
            rows.setdefault(output_index[monomial], {})[
                column_index
            ] = coefficient
    for monomial, coefficient in rhs.items():
        rows.setdefault(output_index[monomial], {})[rhs_column] = -coefficient
    reduced, pivots, nonzero = sdm_irref(rows)
    if rhs_column in pivots:
        raise ValueError("affine correction equation is inconsistent")

    free_columns = [
        column for column in range(len(columns)) if column not in pivots
    ]
    particular: dict[int, object] = {}
    for reduced_row, pivot in enumerate(pivots):
        rhs_value = reduced.get(reduced_row, {}).get(rhs_column, field.zero)
        if rhs_value:
            particular[pivot] = -rhs_value

    # Remove the augmented column before asking SymPy for the homogeneous
    # nullspace.  Its basis is normalized at the ordered free columns.
    homogeneous_rows = {
        row: {
            column: value
            for column, value in entries.items()
            if column != rhs_column
        }
        for row, entries in rows.items()
    }
    homogeneous_reduced, homogeneous_pivots, homogeneous_nonzero = sdm_irref(
        homogeneous_rows
    )
    kernel, _ = sdm_nullspace_from_rref(
        homogeneous_reduced,
        field.one,
        len(columns),
        homogeneous_pivots,
        homogeneous_nonzero,
    )
    if len(kernel) != len(free_columns):
        raise AssertionError("unexpected affine-kernel dimension")
    return particular, kernel, len(pivots)


def split_pair(
    vector: dict[int, object],
    s_monomials: list[Monomial],
    t_monomials: list[Monomial],
) -> tuple[SparsePoly, SparsePoly]:
    split = len(s_monomials)
    return (
        {
            s_monomials[column]: value
            for column, value in vector.items()
            if column < split and value
        },
        {
            t_monomials[column - split]: value
            for column, value in vector.items()
            if column >= split and value
        },
    )


@dataclass
class ThirdOrderFamily:
    base: tuple[SparsePoly, SparsePoly]
    kernel: list[tuple[SparsePoly, SparsePoly]]
    operator_rank: int
    column_count: int


@dataclass
class FirstOrderAudit:
    column_count: int
    operator_rank: int
    kernel_dimension: int
    gauge_rank: int

    @property
    def quotient_dimension(self) -> int:
        return self.kernel_dimension - self.gauge_rank


def first_order_audit(
    S: SparsePoly,
    T: SparsePoly,
    field,
    profile: CorrectionProfile = GENERIC_PROFILE,
) -> FirstOrderAudit:
    """Compute the odd first-correction kernel and Hamiltonian gauge rank."""

    s_monomials = laurent_monomials(
        profile.s1_degree,
        4,
        0,
        profile.z_weight,
    )
    t_monomials = laurent_monomials(
        profile.t1_degree,
        3,
        0,
        profile.z_weight,
    )
    columns = [
        profile.poisson({monomial: field.one}, T)
        for monomial in s_monomials
    ]
    columns += [
        profile.poisson(S, {monomial: field.one})
        for monomial in t_monomials
    ]
    operator_rank = column_rank(columns)
    kernel_dimension = len(columns) - operator_rank

    R = {
        monomial: field(coefficient)
        for monomial, coefficient in profile.R
    }
    gauge_pairs: list[tuple[SparsePoly, SparsePoly]] = []
    for exponent in range(profile.s1_degree // profile.r_degree + 1):
        gauge_pairs.append(
            (poly_power(R, exponent, field.one), {})
        )
    gauge_pairs.append((T, {}))
    for exponent in range(profile.t1_degree // profile.r_degree + 1):
        gauge_pairs.append(
            ({}, poly_power(R, exponent, field.one))
        )
    for s_part, t_part in gauge_pairs:
        if add(
            profile.poisson(s_part, T),
            profile.poisson(S, t_part),
        ):
            raise AssertionError("declared target-Hamiltonian gauge is not closed")
    # The R-powers are independent and T contains Z; the two components are
    # disjoint.  Verify this dimension directly in the ambient correction
    # coordinates to guard against chart-specific degree mistakes.
    ambient_vectors: list[SparsePoly] = []
    for s_part, t_part in gauge_pairs:
        tagged = {
            (x_degree, q_degree, z_degree): value
            for (x_degree, q_degree, z_degree), value in s_part.items()
        }
        tagged.update(
            {
                (x_degree, q_degree, z_degree + 100): value
                for (x_degree, q_degree, z_degree), value in t_part.items()
            }
        )
        ambient_vectors.append(tagged)
    gauge_rank = column_rank(ambient_vectors)
    return FirstOrderAudit(
        len(columns),
        operator_rank,
        kernel_dimension,
        gauge_rank,
    )


def third_order_family(
    S: SparsePoly,
    T: SparsePoly,
    field,
    pole_order: int = 0,
    profile: CorrectionProfile = GENERIC_PROFILE,
) -> ThirdOrderFamily:
    s_monomials = laurent_monomials(
        profile.s2_degree,
        3,
        pole_order,
        profile.z_weight,
    )
    t_monomials = laurent_monomials(
        profile.t2_degree,
        2,
        pole_order,
        profile.z_weight,
    )
    columns = [
        profile.poisson({monomial: field.one}, T)
        for monomial in s_monomials
    ]
    columns += [
        profile.poisson(S, {monomial: field.one})
        for monomial in t_monomials
    ]
    rhs = scale(profile.pi_power(S, T, 3), -field.one / field(24))
    particular, kernel, rank = solve_affine(columns, rhs, field)
    return ThirdOrderFamily(
        split_pair(particular, s_monomials, t_monomials),
        [
            split_pair(vector, s_monomials, t_monomials)
            for vector in kernel
        ],
        rank,
        len(columns),
    )


def fifth_order_coefficients(
    S: SparsePoly,
    T: SparsePoly,
    family: ThirdOrderFamily,
    field,
    profile: CorrectionProfile = GENERIC_PROFILE,
) -> tuple[SparsePoly, list[SparsePoly]]:
    """Return the constant and all nonconstant coefficient vectors of O_5."""

    base_s, base_t = family.base

    def defect(s2, t2):
        value = profile.poisson(s2, t2)
        value = add(
            value,
            profile.pi_power(s2, T, 3),
            field.one / field(24),
        )
        value = add(
            value,
            profile.pi_power(S, t2, 3),
            field.one / field(24),
        )
        value = add(
            value,
            profile.pi_power(S, T, 5),
            field.one / field(1920),
        )
        return value

    constant = defect(base_s, base_t)
    nonconstant: list[SparsePoly] = []
    for basis_s, basis_t in family.kernel:
        diagonal = profile.poisson(basis_s, basis_t)
        shifted = defect(
            add(base_s, basis_s),
            add(base_t, basis_t),
        )
        linear = add(add(shifted, constant, -1), diagonal, -1)
        nonconstant.extend((linear, diagonal))
    for left, right in combinations_with_replacement(
        range(len(family.kernel)), 2
    ):
        if left == right:
            continue
        left_s, left_t = family.kernel[left]
        right_s, right_t = family.kernel[right]
        nonconstant.append(
            add(
                profile.poisson(left_s, right_t),
                profile.poisson(right_s, left_t),
            )
        )
    return constant, nonconstant


def column_rank(columns: list[SparsePoly]) -> int:
    output_monomials = sorted(set().union(*(set(column) for column in columns)))
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    rows: dict[int, dict[int, object]] = {}
    for column_index, column in enumerate(columns):
        for monomial, coefficient in column.items():
            rows.setdefault(output_index[monomial], {})[
                column_index
            ] = coefficient
    _, pivots, _ = sdm_irref(rows)
    return len(pivots)


@dataclass
class FifthOrderAudit:
    correction_columns: int
    correction_rank: int
    obstruction_coefficients: int
    span_rank: int
    augmented_rank: int

    @property
    def constant_outside_span(self) -> bool:
        return self.augmented_rank > self.span_rank


def fifth_order_audit(
    S: SparsePoly,
    T: SparsePoly,
    third: ThirdOrderFamily,
    field,
    correction_pole_order: int = 0,
    profile: CorrectionProfile = GENERIC_PROFILE,
) -> FifthOrderAudit:
    constant, nonconstant = fifth_order_coefficients(
        S,
        T,
        third,
        field,
        profile,
    )
    if correction_pole_order:
        s4_monomials = laurent_monomials(
            profile.s4_degree,
            1,
            correction_pole_order,
            profile.z_weight,
        )
        t4_monomials = laurent_monomials(
            profile.t4_degree,
            0,
            correction_pole_order,
            profile.z_weight,
        )
    else:
        s4_monomials = laurent_monomials(
            profile.s4_degree,
            1,
            0,
            profile.z_weight,
        )
        t4_monomials = laurent_monomials(
            profile.t4_degree,
            0,
            0,
            profile.z_weight,
        )
    correction_columns = [
        profile.poisson({monomial: field.one}, T)
        for monomial in s4_monomials
    ]
    correction_columns += [
        profile.poisson(S, {monomial: field.one})
        for monomial in t4_monomials
    ]
    correction_rank = column_rank(correction_columns)
    span_columns = correction_columns + nonconstant
    span_rank = column_rank(span_columns)
    augmented_rank = column_rank(span_columns + [constant])
    return FifthOrderAudit(
        len(correction_columns),
        correction_rank,
        len(nonconstant),
        span_rank,
        augmented_rank,
    )


def parse_rational(text: str, field):
    value = Fraction(text)
    return field(value.numerator) / field(value.denominator)


def run_seed(
    kappa_text: str,
    tau_text: str,
    max_pole: int,
    prime: int,
    exact: bool = False,
) -> None:
    field = QQ if exact else GF(prime)
    kappa = parse_rational(kappa_text, field)
    tau = parse_rational(tau_text, field)
    if kappa + field(2) == field.zero:
        raise ValueError("kappa=-2 is outside the normalized seed chart")
    if kappa + field.one == field.zero:
        S, T = degree_five_exceptional_family(field, tau)
        profile = EXCEPTIONAL_PROFILE
        chart = "kappa=-1 replacement chart"
    else:
        a = -(field.one + kappa) / (field(2) + kappa)
        S, T = degree_five_family(field, a, tau)
        profile = GENERIC_PROFILE
        chart = "generic a-chart"
    print(
        f"seed (kappa,tau)=({kappa_text},{tau_text}) over "
        f"{'Q' if exact else f'GF({prime})'}; "
        f"{chart}; support S={len(S)} T={len(T)}"
    )
    first = first_order_audit(S, T, field, profile)
    print(
        f"full first order: columns={first.column_count} "
        f"rank={first.operator_rank} kernel={first.kernel_dimension}; "
        f"Hamiltonian gauge={first.gauge_rank}; "
        f"quotient={first.quotient_dimension}"
    )

    for pole_order in range(max_pole + 1):
        third = third_order_family(
            S,
            T,
            field,
            pole_order,
            profile,
        )
        # Products of two hbar^2 corrections can have twice their pole.
        fifth = fifth_order_audit(
            S,
            T,
            third,
            field,
            correction_pole_order=2 * pole_order,
            profile=profile,
        )
        category = "poly" if pole_order == 0 else f"X-Laurent({pole_order})"
        verdict = (
            "OBSTRUCTED: constant outside span"
            if fifth.constant_outside_span
            else "INCONCLUSIVE: constant enters coefficient span"
        )
        print(
            f"{category}: h3 columns={third.column_count} "
            f"rank={third.operator_rank} nullity={len(third.kernel)}; "
            f"h5 correction columns={fifth.correction_columns} "
            f"rank={fifth.correction_rank}; "
            f"quadratic coefficients={fifth.obstruction_coefficients}; "
            f"span ranks={fifth.span_rank}->{fifth.augmented_rank}; {verdict}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kappa", default="0")
    parser.add_argument("--tau", default="1")
    parser.add_argument("--max-pole", type=int, default=2)
    parser.add_argument("--prime", type=int, default=PRIME)
    parser.add_argument(
        "--exact",
        action="store_true",
        help="use exact rational arithmetic instead of a finite field",
    )
    args = parser.parse_args()
    run_seed(args.kappa, args.tau, args.max_pole, args.prime, args.exact)


if __name__ == "__main__":
    main()
