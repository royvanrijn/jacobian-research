#!/usr/bin/env python3
"""Generate and eliminate the reduced F2 modified Laurent systems.

This is a formal-algebra front end, not yet a proof that the F2 corner chain
forces the displayed normal form for every member of the family.  It assumes

    P = C^r,
    Q = C^(2r-1) + lambda*C^(-1) + F,

with leading x-degree d for C.  Under that hypothesis it reconstructs the
finite coefficient system exactly, triangularly removes the coefficients
forced by polynomiality of P, changes to the polynomial coordinates of
P=C^r, and records the remaining residue as a direct Artinian/Fitting
presentation.  A separate Toeplitz determinant records ramification of the
fixed-parameter coefficient map.

The number of visible F rows is *not* silently fixed here.  If the bracket
in the chosen x-chart has x-degree h, the leading F row is forced at
x^(h+1-r*d), so there are h+1 visible rows.  The published r=2 calculation
uses h=2.  The certified F2 j=1 corner-chain chart instead has bracket X^4;
identifying that X-chart with either d=2 or d=3 modified chart is an open
normal-form step.

The r=2, d=2 and d=3 outputs regress the two published (50,75) modified
systems.  The r=3 outputs are the 14- and 22-equation candidate systems for
the first F2 frontier pair (75,125).  A uniform congruence-support test uses
the certified X^4 bracket, and a separate cubic-invariant test classifies a
monomial h=2 bracket without assuming a support mask for F_-6.  The full
systems remain conditional: neither d, the common-power chart, nor the
lower y-support ledger is forced by the currently certified corner chain.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import sympy as sp


EXPECTED_DIGESTS = {
    (2, 2): (
        "d7a2eedcd2dd58792d57954deec2ea778d01ffc66c6d961d46cac8447acf676a",
        "039cfb69ed65e1d847eeeb69012553074380a7a3f385956a43a763fb631059f1",
    ),
    (2, 3): (
        "1f34dd953bcf9043b03702fe8f6fd13ad74055e08787c2028d17a032b358c87c",
        "c32d7f6e29b4e99fb34411399b0062409ff122383f8b502c430ad33c675a7c1f",
    ),
    (3, 2): (
        "8fbe23edcb88dc4505714dc520d5ce46dabc58a7eba3860d37165e3cff835dbe",
        "c2648af67ff81484aa66775bfb2e78f3959520dba7e35d4c8c1f7311172963b0",
    ),
    (3, 3): (
        "d3101adb202ef35b0fa011300c48bb1963cac254612b834ece5256c2cc203ff4",
        "533af5cd52bbd4724ed1591dd20bc21c7fefcd4b75172a96803ad4c6a20ae070",
    ),
}


EXPECTED_COMPACT_DIGESTS = {
    (2, 2, 2): (
        "f9d9247df0ec0cd7c1335bf624c14d67a5b6789863b2171ac7f50eb74a6f0607",
        "c348ecc23f3b70fb6f6f23b17da0c5c1295bbee54067f9d03cb3f198fdb82b2c",
        "a0d6ca9ab0f8e6fc0d6919cf40e91bba1980ac1b78009471e942d49b78e70919",
    ),
    (2, 3, 2): (
        "272c5002d23d2a8e5409e6e8f78ac42ff5b9f2254c692ca2a3f7a0d782bcbc41",
        "04a808c84327f63737383a5fe17e6b7939aedb9413524cbe9c517fe8213ea1c0",
        "bee51ad99c27a167a4a11fb6ed9478e9981cbf0a085568a7b051cf80226e69c3",
    ),
    (3, 2, 2): (
        "7a76c4b866d61d25fa8d4dd90cf9f98ac225a741a003fea4ccae83b0df95f8bd",
        "6161ed5a2c083adbadb6fe176c57831f8b10cd844e704b2392e3b521f569c06a",
        "33ccfd3156b4343f8885d026e809c46648c8f0c4e17d828425e383d672c97742",
    ),
    (3, 3, 2): (
        "8eeb9fdbdf5b12a9fa0d6e0cb31d13b911f7f8e7ed554953407fa4e016710f4f",
        "26977c069dcbca353db5c43f5cb2dda1ab9332159c825e0346ca4d67d760b33e",
        "118005def4e72395e42399a38202f9fa85a246ddb8277868a8187ba33664a9f1",
    ),
    (3, 2, 4): (
        "7a76c4b866d61d25fa8d4dd90cf9f98ac225a741a003fea4ccae83b0df95f8bd",
        "97e160234e5fcdfbd6358fa8668ca5d962aaa381b1aaa8138fc19724c123920f",
        "33ccfd3156b4343f8885d026e809c46648c8f0c4e17d828425e383d672c97742",
    ),
    (3, 3, 4): (
        "8eeb9fdbdf5b12a9fa0d6e0cb31d13b911f7f8e7ed554953407fa4e016710f4f",
        "cea160eef9cd7e47b297f9e8d9ca34037b7b4e382ebc0ae8c84d1ea838b70eb6",
        "118005def4e72395e42399a38202f9fa85a246ddb8277868a8187ba33664a9f1",
    ),
}


EXPECTED_ESSENTIAL_RESIDUE_DIGESTS = {
    (2, 2, 2): "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    (2, 3, 2): "32f3343f4af9265666ced7ac50e4642b2763216a8cde3989c199f333b82548c4",
    (3, 2, 2): "9baf016bfda0f7049e276b7460a60ca0bff2112e59e6df9fcb9b47bc2d73fb29",
    (3, 3, 2): "00fae1edf284221b2dfba8604a0c973d34bf5cb52f59ccbb32a1d0b7222c5865",
    (3, 2, 4): "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    (3, 3, 4): "bcf6aa11956bda1261504c9abb74ff5136f9a8516e6ef13d5385a017795b6b53",
}


EXPECTED_ARTINIAN_REMAINDER_DIGESTS = {
    (2, 2, 2): "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    (2, 3, 2): "32f3343f4af9265666ced7ac50e4642b2763216a8cde3989c199f333b82548c4",
    (3, 2, 2): "9baf016bfda0f7049e276b7460a60ca0bff2112e59e6df9fcb9b47bc2d73fb29",
    (3, 3, 2): "46a369e8b4127cf18107f4f4050f2431dfc1c5d310bf3a5cb08bb0347b53459b",
    (3, 2, 4): "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945",
    (3, 3, 4): "bcf6aa11956bda1261504c9abb74ff5136f9a8516e6ef13d5385a017795b6b53",
}


def _z_name(exponent: int) -> str:
    if exponent < 0:
        return f"zm{-exponent}"
    return f"z{exponent}"


def _mul(left: list[sp.Expr], right: list[sp.Expr], bound: int) -> list[sp.Expr]:
    result = [sp.Integer(0)] * (bound + 1)
    for i, left_coefficient in enumerate(left):
        if left_coefficient == 0:
            continue
        for j, right_coefficient in enumerate(right[: bound + 1 - i]):
            if right_coefficient != 0:
                result[i + j] += left_coefficient * right_coefficient
    return [sp.expand(coefficient) for coefficient in result]


def _power(series: list[sp.Expr], exponent: int, bound: int) -> list[sp.Expr]:
    result = [sp.Integer(0)] * (bound + 1)
    result[0] = 1
    base = list(series)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = _mul(result, base, bound)
        remaining //= 2
        if remaining:
            base = _mul(base, base, bound)
    return result


def _inverse(series: list[sp.Expr], bound: int) -> list[sp.Expr]:
    assert series[0] == 1
    result = [sp.Integer(0)] * (bound + 1)
    result[0] = 1
    for degree in range(1, bound + 1):
        result[degree] = -sum(
            series[index] * result[degree - index]
            for index in range(1, min(degree, len(series) - 1) + 1)
        )
        result[degree] = sp.expand(result[degree])
    return result


def _inverse_rth_root(
    polynomial_series: list[sp.Expr], exponent: int, bound: int
) -> list[sp.Expr]:
    """Return B^(-1/exponent) modulo t^(bound+1), with B(0)=1.

    The recurrence is the coefficient form of

        exponent * B * S' + B' * S = 0,  S(0)=1.

    It avoids symbolic fractional-power expansion and is valid over every
    characteristic-zero coefficient ring.
    """

    assert exponent >= 2
    assert polynomial_series[0] == 1
    result = [sp.Integer(0)] * (bound + 1)
    result[0] = 1
    polynomial_degree = len(polynomial_series) - 1
    for degree in range(1, bound + 1):
        result[degree] = -sum(
            (exponent * (degree - index) + index)
            * polynomial_series[index]
            * result[degree - index]
            for index in range(1, min(degree, polynomial_degree) + 1)
        ) / (exponent * degree)
        result[degree] = sp.expand(result[degree])
    return result


def _canonical_polynomial(expression: sp.Expr, generators: list[sp.Symbol]) -> list:
    polynomial = sp.Poly(sp.expand(expression), *generators, domain=sp.QQ)
    return [
        [list(monomial), [int(coefficient.p), int(coefficient.q)]]
        for monomial, coefficient in polynomial.terms()
    ]


def _digest(polynomials: Iterable[list]) -> str:
    payload = json.dumps(list(polynomials), separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _canonical_matrix(
    matrix: sp.Matrix, generators: list[sp.Symbol]
) -> list[list[list]]:
    return [
        [_canonical_polynomial(matrix[row, column], generators) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _weighted_degree(
    monomial: tuple[int, ...], generator_weights: list[int]
) -> int:
    return sum(exponent * weight for exponent, weight in zip(monomial, generator_weights))


def _assert_weighted_homogeneous(
    expression: sp.Expr,
    generators: list[sp.Symbol],
    generator_weights: list[int],
    expected_weight: int,
) -> None:
    polynomial = sp.Poly(sp.expand(expression), *generators, domain=sp.QQ)
    weights = {
        _weighted_degree(monomial, generator_weights)
        for monomial, _ in polynomial.terms()
    }
    assert weights == {expected_weight}, (expression, weights, expected_weight)


@dataclass(frozen=True)
class SystemSummary:
    r: int
    d: int
    family_parameter_j: int
    reduced_pair: tuple[int, int]
    degree_pair: tuple[int, int]
    degree_pair_provenance: str
    first_block_equations: int
    residual_equations: int
    total_equations: int
    coefficient_variables: int
    free_coefficient_variables: int
    eliminated_coefficient_variables: int
    last_inspected_z_exponent: int
    residual_z_exponents: list[int]
    residual_weights: list[int]
    lambda_weight: int
    lambda_elimination_equation: int
    zero_residual_equations: list[int]
    visible_f_equations: list[int]
    residual_partition_assumption: str
    first_block_digest: str
    residual_digest: str
    first_block_term_counts: list[int]
    residual_term_counts: list[int]


@dataclass
class GeneratedSystem:
    summary: SystemSummary
    symbols_by_q: dict[int, sp.Symbol]
    lambda_symbol: sp.Symbol
    first_block: list[sp.Expr]
    raw_residual_block: list[sp.Expr]
    triangular_substitution: dict[sp.Symbol, sp.Expr]
    residual_block: list[sp.Expr]


@dataclass(frozen=True)
class CompactResidueSummary:
    r: int
    d: int
    jacobian_x_degree: int
    polynomial_degree: int
    truncation_degree: int
    polynomial_variables: list[str]
    zero_rows: list[int]
    visible_f_rows: list[int]
    visible_f_x_exponents: list[int]
    f_tail_continues_below_window: bool
    lambda_pivot_row: int
    lambda_is_independent_zero_row_pivot: bool
    lambda_eliminated_rows: list[int]
    lambda_eliminated_equations: int
    lambda_eliminated_variables: int
    essential_residue_rows: list[int]
    essential_residue_equations: int
    essential_residue_digest: str
    artinian_module_rank: int
    artinian_base_columns: int
    artinian_augmented_columns: int
    artinian_cokernel_rank: int
    artinian_maximal_minor_size: int
    artinian_fitting_index: int | None
    artinian_remainder_digest: str
    core_digest: str
    lambda_eliminated_digest: str
    toeplitz_matrix_digest: str
    toeplitz_shape: tuple[int, int]
    toeplitz_nonzero_witness: dict[str, object]
    coefficient_map_square_fitting_index: int


@dataclass
class CompactResidue:
    summary: CompactResidueSummary
    polynomial_symbols_by_q: dict[int, sp.Symbol]
    f_symbols_by_k: dict[int, sp.Symbol]
    core_block: list[sp.Expr]
    fiber_block: list[sp.Expr]
    lambda_substitution: sp.Expr
    lambda_eliminated_block: list[sp.Expr]
    essential_residue_block: list[sp.Expr]
    artinian_remainder_block: list[sp.Expr]
    toeplitz_matrix: sp.Matrix


def generate_system(r: int, d: int) -> GeneratedSystem:
    if r < 2:
        raise ValueError("r must be at least 2")
    if d < 2:
        raise ValueError("d must be at least 2")

    m = 2 * r - 1
    maximum_q = d * (3 * r - 1) - 1
    symbols_by_q = {
        q: sp.Symbol(_z_name(d - q)) for q in range(2, maximum_q + 1)
    }
    normalized_series = [sp.Integer(0)] * (maximum_q + 1)
    normalized_series[0] = 1
    for q, symbol in symbols_by_q.items():
        normalized_series[q] = symbol

    rth_power = _power(normalized_series, r, maximum_q)
    first_block: list[sp.Expr] = []
    substitution: dict[sp.Symbol, sp.Expr] = {}
    for k in range(1, d * m):
        coefficient_degree = d * r + k
        raw_equation = sp.expand(rth_power[coefficient_degree])
        first_block.append(raw_equation)

        equation = sp.expand(raw_equation.subs(substitution))
        deepest = symbols_by_q[coefficient_degree]
        assert equation.coeff(deepest) == r
        substitution[deepest] = sp.expand(-(equation - r * deepest) / r)

    lambda_symbol = sp.Symbol("lam")
    mth_power = _power(normalized_series, m, maximum_q)
    inverse = _inverse(normalized_series, max(0, d * r - 1 - d))
    raw_residual_block: list[sp.Expr] = []
    residual_block: list[sp.Expr] = []
    for k in range(1, d * r):
        equation = mth_power[d * m + k]
        if k >= d:
            equation += lambda_symbol * inverse[k - d]
        raw_residual_block.append(sp.expand(equation))
        residual_block.append(sp.factor(sp.expand(equation.subs(substitution))))

    all_z = [symbols_by_q[q] for q in sorted(symbols_by_q)]
    residual_z = [symbols_by_q[q] for q in range(2, d * r + 1)]
    first_generators = all_z + [lambda_symbol]
    residual_generators = residual_z + [lambda_symbol]
    first_canonical = [
        _canonical_polynomial(equation, first_generators) for equation in first_block
    ]
    residual_canonical = [
        _canonical_polynomial(equation, residual_generators)
        for equation in residual_block
    ]

    first_digest = _digest(first_canonical)
    residual_digest = _digest(residual_canonical)
    expected = EXPECTED_DIGESTS.get((r, d))
    if expected is not None:
        assert (first_digest, residual_digest) == expected

    all_weights = list(range(2, maximum_q + 1)) + [2 * d * r]
    residual_weights = list(range(2, d * r + 1))
    residual_generator_weights = residual_weights + [2 * d * r]
    for k, equation in enumerate(first_block, start=1):
        _assert_weighted_homogeneous(
            equation,
            first_generators,
            all_weights,
            d * r + k,
        )
    for k, equation in enumerate(residual_block, start=1):
        _assert_weighted_homogeneous(
            equation,
            residual_generators,
            residual_generator_weights,
            d * m + k,
        )

    for k, equation in enumerate(raw_residual_block, start=1):
        lambda_coefficient = sp.expand(equation).coeff(lambda_symbol)
        if k < d:
            assert lambda_coefficient == 0
        elif k == d:
            assert lambda_coefficient == 1

    first_count = d * m - 1
    residual_count = d * r - 1
    total_count = d * (3 * r - 1) - 2
    assert len(first_block) == first_count
    assert len(residual_block) == residual_count
    assert first_count + residual_count == total_count
    assert len(symbols_by_q) == total_count
    assert set(substitution) == {
        symbols_by_q[q] for q in range(d * r + 1, maximum_q + 1)
    }

    summary = SystemSummary(
        r=r,
        d=d,
        family_parameter_j=r - 2,
        reduced_pair=(r, m),
        degree_pair=(25 * r, 25 * m),
        degree_pair_provenance=(
            "F2 table: (m,n)=(j+2,2j+3), j=r-2, "
            "and v_11(A0)=v_11(5,20)=25"
        ),
        first_block_equations=first_count,
        residual_equations=residual_count,
        total_equations=total_count,
        coefficient_variables=total_count,
        free_coefficient_variables=residual_count,
        eliminated_coefficient_variables=first_count,
        last_inspected_z_exponent=1 - d * (3 * r - 2),
        residual_z_exponents=[d - q for q in range(2, d * r + 1)],
        residual_weights=residual_weights,
        lambda_weight=2 * d * r,
        lambda_elimination_equation=d,
        zero_residual_equations=list(range(1, max(1, d * r - 3))),
        visible_f_equations=list(range(max(1, d * r - 3), d * r)),
        residual_partition_assumption=(
            "historical modified chart with bracket x-degree h=2; "
            "use compact_residues for an explicit h-dependent partition"
        ),
        first_block_digest=first_digest,
        residual_digest=residual_digest,
        first_block_term_counts=[len(polynomial) for polynomial in first_canonical],
        residual_term_counts=[len(polynomial) for polynomial in residual_canonical],
    )
    return GeneratedSystem(
        summary=summary,
        symbols_by_q=symbols_by_q,
        lambda_symbol=lambda_symbol,
        first_block=first_block,
        raw_residual_block=raw_residual_block,
        triangular_substitution=substitution,
        residual_block=residual_block,
    )


def generate_compact_residue(
    system: GeneratedSystem, jacobian_x_degree: int
) -> CompactResidue:
    """Change to P-polynomial coordinates and eliminate lambda exactly.

    Put t=x^-1 and C=x^d A(t).  Once the negative coefficients of A^r
    have been removed, write

        B(t) = A(t)^r = 1 + p_2 t^2 + ... + p_(rd) t^(rd),
        S(t) = B(t)^(-1/r).

    Since 2r-1 is one less than 2r, the Q core is encoded by

        B(t)^2 S(t) + lambda t^(2rd) S(t).

    The returned block consists of its coefficients in degrees
    d(2r-1)+1 through d(2r-1)+rd-1.  Its derivative in the p variables
    is the Toeplitz matrix recorded below.  Separately, the residue itself
    is presented as an Artinian membership/Fitting condition after the
    freely adjustable visible F rows and the lambda pivot are removed.
    """

    r = system.summary.r
    d = system.summary.d
    m = 2 * r - 1
    polynomial_degree = d * r
    leading_q_degree = d * m
    truncation_degree = d * (3 * r - 1) - 1
    if not 0 <= jacobian_x_degree <= polynomial_degree - 2:
        raise ValueError(
            "jacobian x-degree must leave at least one negative F coefficient"
        )

    polynomial_symbols_by_q = {
        q: sp.Symbol(f"p{q}") for q in range(2, polynomial_degree + 1)
    }
    polynomial_series = [sp.Integer(0)] * (truncation_degree + 1)
    polynomial_series[0] = 1
    for q, symbol in polynomial_symbols_by_q.items():
        polynomial_series[q] = symbol

    root_inverse = _inverse_rth_root(polynomial_series, r, truncation_degree)
    polynomial_square = _power(polynomial_series, 2, truncation_degree)
    main_core = _mul(polynomial_square, root_inverse, truncation_degree)
    lambda_symbol = system.lambda_symbol
    core_block = []
    for k in range(1, polynomial_degree):
        equation = main_core[leading_q_degree + k]
        if k >= d:
            equation += lambda_symbol * root_inverse[k - d]
        core_block.append(sp.expand(equation))

    # The p-coordinates are a triangular change from the surviving root
    # coefficients, with diagonal r.  Verify every compact row against the
    # independently generated root-coordinate row.
    root_series = [sp.Integer(0)] * (polynomial_degree + 1)
    root_series[0] = 1
    for q in range(2, polynomial_degree + 1):
        root_series[q] = system.symbols_by_q[q]
    root_power = _power(root_series, r, polynomial_degree)
    polynomial_to_root = {
        polynomial_symbols_by_q[q]: root_power[q]
        for q in range(2, polynomial_degree + 1)
    }
    for compact, root_coordinate in zip(core_block, system.residual_block):
        assert sp.expand(compact.subs(polynomial_to_root) - root_coordinate) == 0

    tail_start = polynomial_degree - jacobian_x_degree - 1
    f_symbols_by_k = {
        k: sp.Symbol(f"f{k}") for k in range(tail_start, polynomial_degree)
    }
    fiber_block = [
        sp.expand(
            equation + f_symbols_by_k.get(k, sp.Integer(0))
        )
        for k, equation in enumerate(core_block, start=1)
    ]

    pivot = fiber_block[d - 1]
    assert sp.expand(pivot).coeff(lambda_symbol) == 1
    lambda_substitution = sp.expand(-(pivot - lambda_symbol))
    lambda_eliminated_rows = [
        k for k in range(1, polynomial_degree) if k != d
    ]
    lambda_eliminated_block = [
        sp.factor(sp.expand(fiber_block[k - 1].subs(lambda_symbol, lambda_substitution)))
        for k in lambda_eliminated_rows
    ]

    zero_rows = list(range(1, tail_start))
    visible_f_rows = list(range(tail_start, polynomial_degree))
    lambda_is_independent_zero_row_pivot = d < tail_start
    if lambda_is_independent_zero_row_pivot:
        essential_residue_rows = [k for k in zero_rows if k != d]
        essential_residue_block = [
            sp.factor(
                sp.expand(
                    core_block[k - 1].subs(
                        lambda_symbol,
                        lambda_substitution,
                    )
                )
            )
            for k in essential_residue_rows
        ]
    else:
        # Here the lambda row is already a freely adjustable F row.  Every
        # zero row precedes lambda and is independent of it.
        essential_residue_rows = zero_rows
        essential_residue_block = [
            sp.factor(core_block[k - 1]) for k in essential_residue_rows
        ]
        assert all(
            equation.coeff(lambda_symbol) == 0
            for equation in essential_residue_block
        )

    polynomial_generators = [
        polynomial_symbols_by_q[q] for q in range(2, polynomial_degree + 1)
    ]
    f_generators = [f_symbols_by_k[k] for k in sorted(f_symbols_by_k)]
    core_generators = polynomial_generators + [lambda_symbol]
    residue_generators = polynomial_generators + f_generators
    essential_generators = polynomial_generators

    polynomial_weights = list(range(2, polynomial_degree + 1))
    f_weights = [leading_q_degree + k for k in sorted(f_symbols_by_k)]
    for k, equation in enumerate(fiber_block, start=1):
        _assert_weighted_homogeneous(
            equation,
            core_generators + f_generators,
            polynomial_weights + [2 * polynomial_degree] + f_weights,
            leading_q_degree + k,
        )
    for k, equation in zip(lambda_eliminated_rows, lambda_eliminated_block):
        _assert_weighted_homogeneous(
            equation,
            residue_generators,
            polynomial_weights + f_weights,
            leading_q_degree + k,
        )
    for k, equation in zip(
        essential_residue_rows,
        essential_residue_block,
        strict=True,
    ):
        _assert_weighted_homogeneous(
            equation,
            essential_generators,
            polynomial_weights,
            leading_q_degree + k,
        )

    # Differentiate the compact coefficient operator.  If
    # H=B^(2-1/r)+lambda*t^(2rd)B^(-1/r), then
    # dH=G*dB, with G the series below.  Coefficient extraction makes the
    # residual Jacobian a Toeplitz matrix.
    polynomial_inverse = _inverse(polynomial_series, truncation_degree)
    inverse_power = _mul(root_inverse, polynomial_inverse, truncation_degree)
    polynomial_times_root_inverse = _mul(
        polynomial_series, root_inverse, truncation_degree
    )
    derivative_multiplier = [
        sp.Rational(m, r) * coefficient
        for coefficient in polynomial_times_root_inverse
    ]
    for degree in range(2 * polynomial_degree, truncation_degree + 1):
        derivative_multiplier[degree] -= (
            sp.Rational(1, r)
            * lambda_symbol
            * inverse_power[degree - 2 * polynomial_degree]
        )
        derivative_multiplier[degree] = sp.expand(derivative_multiplier[degree])

    toeplitz_matrix = sp.Matrix(
        polynomial_degree - 1,
        polynomial_degree - 1,
        lambda row, column: derivative_multiplier[
            leading_q_degree + (row + 1) - (column + 2)
        ],
    )
    direct_jacobian = sp.Matrix(core_block).jacobian(polynomial_generators)
    assert toeplitz_matrix == direct_jacobian

    matrix_canonical = _canonical_matrix(toeplitz_matrix, core_generators)
    matrix_digest = hashlib.sha256(
        json.dumps(matrix_canonical, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    witness_substitution = {
        symbol: ((index + 2) * (index + 3) + index * index) % 11 - 5
        for index, symbol in enumerate(polynomial_generators)
    }
    witness_substitution[lambda_symbol] = 1
    witness_determinant = sp.factor(
        toeplitz_matrix.subs(witness_substitution).det()
    )
    assert witness_determinant != 0

    core_canonical = [
        _canonical_polynomial(equation, core_generators) for equation in core_block
    ]
    residue_canonical = [
        _canonical_polynomial(equation, residue_generators)
        for equation in lambda_eliminated_block
    ]
    essential_canonical = [
        _canonical_polynomial(equation, essential_generators)
        for equation in essential_residue_block
    ]
    core_digest = _digest(core_canonical)
    residue_digest = _digest(residue_canonical)
    essential_digest = _digest(essential_canonical)
    expected_digests = EXPECTED_COMPACT_DIGESTS.get(
        (r, d, jacobian_x_degree)
    )
    if expected_digests is not None:
        assert (core_digest, residue_digest, matrix_digest) == expected_digests
    expected_essential_digest = EXPECTED_ESSENTIAL_RESIDUE_DIGESTS.get(
        (r, d, jacobian_x_degree)
    )
    if expected_essential_digest is not None:
        assert essential_digest == expected_essential_digest
    # Direct residue presentation.  In J=R[t]/(t^(L+1)), put
    # Y=B^(1/r), U=<1,...,t^(d(2r-1))>, and let T_h be the h+1 terminal
    # monomials.  The Q condition is
    #
    #     B^2 in Y*U + T_h + R*t^(2rd).
    #
    # The lambda column is independent exactly when its row precedes T_h.
    # Appending B^2 gives Fitt_(q-1) of the augmented cokernel.  Because
    # all base columns have distinct unit leading terms, its maximal-minor
    # ideal is exactly the ideal generated by essential_residue_block.
    artinian_module_rank = truncation_degree + 1
    lambda_column_epsilon = int(lambda_is_independent_zero_row_pivot)
    artinian_base_columns = (
        leading_q_degree
        + jacobian_x_degree
        + 2
        + lambda_column_epsilon
    )
    artinian_augmented_columns = artinian_base_columns + 1
    artinian_cokernel_rank = artinian_module_rank - artinian_base_columns
    assert artinian_cokernel_rank == len(essential_residue_block)
    assert artinian_augmented_columns <= artinian_module_rank + 1
    artinian_fitting_index = (
        artinian_cokernel_rank - 1
        if artinian_cokernel_rank > 0
        else None
    )

    # Build and row-reduce the presentation itself.  The Y*U columns have
    # successive unit pivots 0,...,d(2r-1); T_h and the independent lambda
    # column are monomials at their displayed degrees.  Reducing the appended
    # B^2 column therefore computes its maximal minors without expanding a
    # large determinant.
    root_series_in_polynomial_coordinates = _inverse(
        root_inverse,
        truncation_degree,
    )
    unit_check = _mul(
        root_series_in_polynomial_coordinates,
        root_inverse,
        truncation_degree,
    )
    assert unit_check == [sp.Integer(1)] + [
        sp.Integer(0)
    ] * truncation_degree
    artinian_remainder = list(polynomial_square)
    for pivot_degree in range(leading_q_degree + 1):
        pivot_coefficient = artinian_remainder[pivot_degree]
        for degree in range(pivot_degree, truncation_degree + 1):
            artinian_remainder[degree] = sp.expand(
                artinian_remainder[degree]
                - pivot_coefficient
                * root_series_in_polynomial_coordinates[
                    degree - pivot_degree
                ]
            )
    assert all(
        coefficient == 0
        for coefficient in artinian_remainder[: leading_q_degree + 1]
    )

    if lambda_is_independent_zero_row_pivot:
        artinian_remainder[2 * polynomial_degree] = sp.Integer(0)
    terminal_monomial_start = truncation_degree - jacobian_x_degree
    for degree in range(terminal_monomial_start, truncation_degree + 1):
        artinian_remainder[degree] = sp.Integer(0)

    artinian_residual_degrees = [
        leading_q_degree + k for k in essential_residue_rows
    ]
    artinian_residual_degree_set = set(artinian_residual_degrees)
    assert all(
        coefficient == 0 or degree in artinian_residual_degree_set
        for degree, coefficient in enumerate(artinian_remainder)
    )
    artinian_remainder_block = [
        sp.factor(artinian_remainder[degree])
        for degree in artinian_residual_degrees
    ]

    # Multiplication by the unit Y=B^(1/r) gives the explicit lower-unitriangular
    # change from the compact zero rows to these maximal-minor generators.
    lambda_on_residue = (
        lambda_substitution
        if lambda_is_independent_zero_row_pivot
        else sp.Integer(0)
    )
    high_q_series = [sp.Integer(0)] * (truncation_degree + 1)
    for k, equation in enumerate(core_block, start=1):
        high_q_series[leading_q_degree + k] = sp.expand(
            equation.subs(lambda_symbol, lambda_on_residue)
        )
    transformed_high_q = _mul(
        root_series_in_polynomial_coordinates,
        high_q_series,
        truncation_degree,
    )
    assert all(
        sp.expand(transformed_high_q[degree] - remainder) == 0
        for degree, remainder in zip(
            artinian_residual_degrees,
            artinian_remainder_block,
            strict=True,
        )
    )
    assert all(
        sp.expand(
            high_q_series[leading_q_degree + k] - equation
        )
        == 0
        for k, equation in zip(
            essential_residue_rows,
            essential_residue_block,
            strict=True,
        )
    )
    artinian_remainder_canonical = [
        _canonical_polynomial(equation, polynomial_generators)
        for equation in artinian_remainder_block
    ]
    artinian_remainder_digest = _digest(artinian_remainder_canonical)
    expected_artinian_remainder_digest = (
        EXPECTED_ARTINIAN_REMAINDER_DIGESTS.get(
            (r, d, jacobian_x_degree)
        )
    )
    if expected_artinian_remainder_digest is not None:
        assert (
            artinian_remainder_digest
            == expected_artinian_remainder_digest
        )
    summary = CompactResidueSummary(
        r=r,
        d=d,
        jacobian_x_degree=jacobian_x_degree,
        polynomial_degree=polynomial_degree,
        truncation_degree=truncation_degree,
        polynomial_variables=[str(symbol) for symbol in polynomial_generators],
        zero_rows=zero_rows,
        visible_f_rows=visible_f_rows,
        visible_f_x_exponents=[-k for k in visible_f_rows],
        f_tail_continues_below_window=True,
        lambda_pivot_row=d,
        lambda_is_independent_zero_row_pivot=lambda_is_independent_zero_row_pivot,
        lambda_eliminated_rows=lambda_eliminated_rows,
        lambda_eliminated_equations=polynomial_degree - 2,
        lambda_eliminated_variables=polynomial_degree - 1,
        essential_residue_rows=essential_residue_rows,
        essential_residue_equations=len(essential_residue_block),
        essential_residue_digest=essential_digest,
        artinian_module_rank=artinian_module_rank,
        artinian_base_columns=artinian_base_columns,
        artinian_augmented_columns=artinian_augmented_columns,
        artinian_cokernel_rank=artinian_cokernel_rank,
        artinian_maximal_minor_size=artinian_augmented_columns,
        artinian_fitting_index=artinian_fitting_index,
        artinian_remainder_digest=artinian_remainder_digest,
        core_digest=core_digest,
        lambda_eliminated_digest=residue_digest,
        toeplitz_matrix_digest=matrix_digest,
        toeplitz_shape=toeplitz_matrix.shape,
        toeplitz_nonzero_witness={
            "p_values": [
                witness_substitution[symbol] for symbol in polynomial_generators
            ],
            "lambda": witness_substitution[lambda_symbol],
            "determinant": str(witness_determinant),
        },
        coefficient_map_square_fitting_index=0,
    )
    return CompactResidue(
        summary=summary,
        polynomial_symbols_by_q=polynomial_symbols_by_q,
        f_symbols_by_k=f_symbols_by_k,
        core_block=core_block,
        fiber_block=fiber_block,
        lambda_substitution=lambda_substitution,
        lambda_eliminated_block=lambda_eliminated_block,
        essential_residue_block=essential_residue_block,
        artinian_remainder_block=artinian_remainder_block,
        toeplitz_matrix=toeplitz_matrix,
    )


def family_formula_audit() -> dict[str, object]:
    """Return the all-r formulas and their exact claim boundary."""

    r, d, h = sp.symbols("r d h", integer=True, positive=True)
    m = 2 * r - 1
    polynomial_degree = d * r
    truncation_degree = d * (3 * r - 1) - 1
    total_variables = d * (3 * r - 1) - 2
    first_rows = d * m - 1
    residual_rows = polynomial_degree - 1
    tail_start = polynomial_degree - h - 1
    assert sp.expand(first_rows + residual_rows - total_variables) == 0
    assert sp.expand(truncation_degree - 1 - total_variables) == 0
    assert sp.expand(polynomial_degree - tail_start - h - 1) == 0
    return {
        "family_table": {
            "parameter_change": "j=r-2",
            "reduced_pair": ["r", "2*r-1"],
            "initial_corner": [5, 20],
            "corner_total_degree": 25,
            "degree_pair": ["25*r", "25*(2*r-1)"],
        },
        "corner_chain_support_envelope": {
            "forced_P_vertices": [
                ["25*r", "20*r"],
                ["7*r", "2*r"],
                [4, 1],
            ],
            "forced_Q_vertices": [
                ["25*(2*r-1)", "20*(2*r-1)"],
                ["7*(2*r-1)", "2*(2*r-1)"],
                [1, 0],
            ],
            "source_monomial_change": "ell=5*i-j",
            "source_degree_envelope": (
                "0<=i, 0<=5*i-ell, 6*i-ell<=D; "
                "D=25*r for P and D=25*(2*r-1) for Q"
            ),
            "terminal_halfspaces": [
                "(2*r-1)*a-(7*r-4)*b<=r for P",
                "(2*r-1)*a-(7*r-4)*b<=2*r-1 for Q",
            ],
            "scope": (
                "exact B0 jet envelope and actual forced vertices, "
                "not an exhaustive lower Laurent support mask"
            ),
        },
        "normalized_series": (
            "t=x^-1; C=x^d*A(t); "
            "A(t)=1+sum_(q>=2) z_(d-q)*t^q"
        ),
        "finite_coefficient_window": {
            "last_inspected_t_index": "d*(3*r-1)-1",
            "last_inspected_C_x_exponent": "1-d*(3*r-2)",
            "inspected_coefficient_functions": "d*(3*r-1)-2",
            "power_rows_and_eliminated_variables": "d*(2*r-1)-1",
            "surviving_variables_and_Q_rows": "d*r-1",
            "scope": (
                "A(t) is an infinite formal series; these are finite-window "
                "bounds, not a lower support bound for C"
            ),
        },
        "triangular_elimination": {
            "polynomial_coordinate": (
                "B(t)=trunc_(<=d*r)(A(t)^r)="
                "1+sum_(q=2)^(d*r) p_q*t^q"
            ),
            "inverse_root": "S(t)=B(t)^(-1/r)",
            "Q_core": (
                "[t^(d*(2*r-1)+k)] "
                "(B(t)^2*S(t)+lambda*t^(2*d*r)*S(t))"
            ),
            "lambda_pivot": "k=d",
            "post_lambda_counts": {
                "equations": "d*r-2",
                "variables": "d*r-1",
            },
        },
        "bracket_tail_formula": {
            "assumption": "the chosen modified chart has [P,Q] of x-degree h",
            "leading_F_x_exponent": "h+1-d*r",
            "first_F_row": "d*r-h-1",
            "visible_F_rows": "h+1",
            "scope": "the formal F tail continues below the visible window",
            "historical_r2_chart": "h=2 gives three rows",
            "certified_F2_j1_X_chart": "h=4 gives five rows",
        },
        "direct_residue_fitting_formula": {
            "artinian_ring": (
                "R=Q[1/r,p_2,...,p_(d*r)], "
                "J=R[t]/(t^(d*(3*r-1)))"
            ),
            "membership": (
                "B^2 in B^(1/r)*<1,...,t^(d*(2*r-1))> "
                "+ <t^(L-h),...,t^L> + R*t^(2*d*r)"
            ),
            "lambda_column_epsilon": "1 if d<d*r-h-1, otherwise 0",
            "base_columns": "d*(2*r-1)+h+2+epsilon",
            "cokernel_rank": "q=d*r-h-2-epsilon",
            "residue_ideal": (
                "for q>0, Fitt_(q-1) of the cokernel after appending B^2; "
                "equivalently its maximal minors"
            ),
        },
        "coefficient_map_toeplitz_fitting_formula": {
            "multiplier": (
                "G(t)=((2*r-1)/r)*B(t)*S(t) "
                "-(lambda/r)*t^(2*d*r)*S(t)/B(t)"
            ),
            "entry": (
                "T_(k,q)=[t^(d*(2*r-1)+k-q)]G(t), "
                "1<=k<d*r, 2<=q<=d*r"
            ),
            "full_jacobian": (
                "with columns ordered eliminated-then-free, "
                "det(J_full)=r^(d*(3*r-1)-2)*det(T)"
            ),
            "kahler_statement": (
                "on a fixed lambda/F fiber, Fitt_0 of the coefficient-map "
                "Kahler module is generated by det(T); this is a generic "
                "ramification invariant, not the residue obstruction ideal"
            ),
        },
        "not_derived_from_corner_chain": [
            "the common-power equality beyond the forced leading bands",
            "the d=2 and d=3 chart list and its exhaustiveness",
            "the Laurent-y supports and endpoint nonvanishing of every C coefficient",
            "the F-tail supports in either modified chart",
            "a coordinate change carrying the certified X^4 bracket to historical h=2",
        ],
    }


def verify_published_r2_regression() -> None:
    gamma3 = generate_system(2, 2)
    z = gamma3.symbols_by_q
    lam = gamma3.lambda_symbol
    assert gamma3.first_block[0] == 2 * z[2] * z[3] + 2 * z[5]
    assert gamma3.raw_residual_block[0] == (
        3 * z[2] ** 2 * z[3] + 6 * z[2] * z[5] + 6 * z[3] * z[4] + 3 * z[7]
    )
    assert gamma3.raw_residual_block[1].coeff(lam) == 1
    assert gamma3.summary.total_equations == 8
    assert gamma3.summary.coefficient_variables == 8
    assert gamma3.summary.last_inspected_z_exponent == -7
    expected_gamma3_residue = [
        3 * z[3] * z[4],
        lam - sp.Rational(3, 2) * z[2] * z[3] ** 2
        + sp.Rational(3, 2) * z[4] ** 2,
        -z[3] * (z[3] ** 2 + 6 * z[2] * z[4]) / 2,
    ]
    assert all(
        sp.expand(actual - expected) == 0
        for actual, expected in zip(
            gamma3.residual_block, expected_gamma3_residue, strict=True
        )
    )

    gamma2 = generate_system(2, 3)
    z = gamma2.symbols_by_q
    lam = gamma2.lambda_symbol
    assert gamma2.first_block[0] == 2 * z[3] * z[4] + 2 * z[2] * z[5] + 2 * z[7]
    assert gamma2.first_block[-1].coeff(z[14]) == 2
    assert gamma2.raw_residual_block[2].coeff(lam) == 1
    assert gamma2.summary.total_equations == 13
    assert gamma2.summary.coefficient_variables == 13
    assert gamma2.summary.last_inspected_z_exponent == -11
    print("F2_50_75_PUBLISHED_REGRESSION_PASS")


def verify_r3_residue_presentations(
    compact_residues: list[CompactResidue],
) -> None:
    by_case = {
        (
            residue.summary.r,
            residue.summary.d,
            residue.summary.jacobian_x_degree,
        ): residue
        for residue in compact_residues
    }

    d2_h2 = by_case[(3, 2, 2)]
    assert d2_h2.summary.essential_residue_rows == [1]
    assert d2_h2.summary.artinian_module_rank == 16
    assert d2_h2.summary.artinian_augmented_columns == 16
    assert d2_h2.summary.artinian_fitting_index == 0
    assert sp.expand(
        d2_h2.essential_residue_block[0] - d2_h2.core_block[0]
    ) == 0
    assert sp.expand(
        d2_h2.artinian_remainder_block[0]
        - d2_h2.essential_residue_block[0]
    ) == 0

    d3_h2 = by_case[(3, 3, 2)]
    assert d3_h2.summary.essential_residue_rows == [1, 2, 4, 5]
    assert d3_h2.summary.artinian_module_rank == 24
    assert d3_h2.summary.artinian_augmented_columns == 21
    assert d3_h2.summary.artinian_fitting_index == 3
    assert len(d3_h2.artinian_remainder_block) == 4

    d2_h4 = by_case[(3, 2, 4)]
    assert d2_h4.summary.essential_residue_equations == 0
    assert d2_h4.summary.artinian_fitting_index is None

    d3_h4 = by_case[(3, 3, 4)]
    assert d3_h4.summary.essential_residue_rows == [1, 2]
    assert d3_h4.summary.artinian_fitting_index == 1
    print("F2_MODIFIED_ARTINIAN_RESIDUE_FITTING_PASS")


def endpoint_binomial_branch_audit(
    compact_residues: list[CompactResidue],
) -> dict[str, object]:
    """Verify the universal endpoint-binomial section of every residue.

    If ``n=d*r`` and ``B=1+a*t^n``, both ``B^(2-1/r)`` and
    ``B^(-1/r)`` are supported in degrees divisible by ``n``.  In the
    compiled row window ``1<=k<n``, the only possible core coefficient is
    therefore the lambda pivot ``k=d``.  The closed lambda value below kills
    that pivot.  This proves that every finite residue/Fitting ideal in this
    compiler is proper; it is not a claim that the point satisfies a
    Laurent-y support mask or the nonzero-bracket normalization.
    """

    endpoint_parameter = sp.Symbol("a_endpoint", nonzero=True)
    checked_cases: list[dict[str, object]] = []
    for residue in compact_residues:
        r = residue.summary.r
        d = residue.summary.d
        n = residue.summary.polynomial_degree
        lambda_value = -sp.Rational(
            (2 * r - 1) * (r - 1), 2 * r * r
        ) * endpoint_parameter**2
        substitution = {
            symbol: sp.Integer(0)
            for symbol in residue.polynomial_symbols_by_q.values()
        }
        substitution[residue.polynomial_symbols_by_q[n]] = endpoint_parameter
        substitution[sp.Symbol("lam")] = lambda_value
        substitution.update(
            {symbol: sp.Integer(0) for symbol in residue.f_symbols_by_k.values()}
        )
        assert all(
            sp.expand(equation.subs(substitution)) == 0
            for equation in residue.core_block
        )
        assert all(
            sp.expand(equation.subs(substitution)) == 0
            for equation in residue.fiber_block
        )
        checked_cases.append(
            {
                "r": r,
                "d": d,
                "jacobian_x_degree": residue.summary.jacobian_x_degree,
                "polynomial_degree_n": n,
                "core_rows_killed": n - 1,
            }
        )

    r_symbol = sp.Symbol("r", nonzero=True)
    a_symbol = sp.Symbol("a", nonzero=True)
    alpha = 2 - 1 / r_symbol
    beta = -1 / r_symbol
    second_coefficient = sp.factor(alpha * (alpha - 1) / 2)
    lambda_formula = -second_coefficient * a_symbol**2
    third_coefficient = sp.factor(
        alpha * (alpha - 1) * (alpha - 2) * a_symbol**3 / 6
        + lambda_formula * beta * a_symbol
    )
    expected_third = sp.factor(
        (2 * r_symbol - 1)
        * (r_symbol - 1)
        * a_symbol**3
        / (3 * r_symbol**3)
    )
    assert sp.factor(third_coefficient - expected_third) == 0

    return {
        "assumptions": "r>=2, d>=2, a!=0, n=d*r",
        "section": "B(t)=1+a*t^n",
        "lambda": "-((2*r-1)*(r-1)/(2*r^2))*a^2",
        "compiled_row_statement": (
            "R_k=0 for every 1<=k<n, with every visible F coefficient set "
            "to zero; hence every compiled residue/Fitting ideal is proper"
        ),
        "support_reason": (
            "d*(2*r-1)+k=2*n-d+k is divisible by n only at k=d, "
            "and 0<=k-d<n has B^(-1/r) support only at k=d"
        ),
        "first_nonzero_lower_core_row": {
            "row": "k=n+d",
            "Q_x_exponent": "-(n+d)",
            "coefficient": "((2*r-1)*(r-1)/(3*r^3))*a^3",
        },
        "bracket_warning": (
            "this endpoint section sets the first visible F coefficient to "
            "zero unless its row is the lambda pivot; a nonzero bracket can "
            "therefore exclude this special section even though it proves the "
            "residue ideal is not the unit ideal"
        ),
        "checked_compiler_cases": checked_cases,
    }


def congruence_support_bracket_audit(
    compact_residues: list[CompactResidue],
) -> dict[str, object]:
    """Test ``B(t) in K[t^d]`` against the first bracket-visible row.

    If the bracket has x-degree ``h``, its first visible row is
    ``k0=d*r-h-1``.  On the coefficient congruence section, both terms in
    that row have series indices congruent to ``-(h+1)`` modulo ``d``.
    Thus ``d`` not dividing ``h+1`` forces the leading F coefficient to
    vanish and contradicts a nonzero bracket of degree ``h``.
    """

    r_symbol, d_symbol, h_symbol = sp.symbols(
        "r d h", integer=True, positive=True
    )
    n_symbol = d_symbol * r_symbol
    k0_symbol = n_symbol - h_symbol - 1
    main_index = sp.expand(d_symbol * (2 * r_symbol - 1) + k0_symbol)
    lambda_index = sp.expand(k0_symbol - d_symbol)
    assert sp.expand(
        main_index - (d_symbol * (3 * r_symbol - 1) - h_symbol - 1)
    ) == 0
    assert sp.expand(
        lambda_index - (d_symbol * (r_symbol - 1) - h_symbol - 1)
    ) == 0
    assert sp.expand(main_index - lambda_index) == 2 * d_symbol * r_symbol

    checked_cases: list[dict[str, object]] = []
    for residue in compact_residues:
        r = residue.summary.r
        d = residue.summary.d
        h = residue.summary.jacobian_x_degree
        n = residue.summary.polynomial_degree
        k0 = n - h - 1
        assert 1 <= k0 < n
        congruence_parameters = {
            q: sp.Symbol(f"c{q // d}") for q in range(d, n + 1, d)
        }
        substitution = {
            symbol: congruence_parameters.get(q, sp.Integer(0))
            for q, symbol in residue.polynomial_symbols_by_q.items()
        }

        # Every compiled zero-row residue vanishes on the congruence section.
        assert all(
            sp.expand(equation.subs(substitution)) == 0
            for equation in residue.essential_residue_block
        )

        first_visible_core = sp.factor(
            residue.core_block[k0 - 1].subs(substitution)
        )
        excluded_by_support = (h + 1) % d != 0
        if excluded_by_support:
            assert first_visible_core == 0
        else:
            # Divisibility only says that this row may survive; the two
            # compiled h=2,d=3 cases verify that it is not identically zero.
            assert first_visible_core != 0
        checked_cases.append(
            {
                "r": r,
                "d": d,
                "jacobian_x_degree_h": h,
                "first_visible_row_k0": k0,
                "main_series_index": d * (3 * r - 1) - h - 1,
                "lambda_series_index": d * (r - 1) - h - 1,
                "d_divides_h_plus_1": not excluded_by_support,
                "first_visible_core_on_congruence_section": str(
                    first_visible_core
                ),
                "nonzero_bracket_verdict": (
                    "excluded" if excluded_by_support else "not decided"
                ),
            }
        )

    certified_x4_cases = [
        case
        for case in checked_cases
        if case["r"] == 3 and case["jacobian_x_degree_h"] == 4
    ]
    assert len(certified_x4_cases) == 2
    assert all(
        case["nonzero_bracket_verdict"] == "excluded"
        for case in certified_x4_cases
    )
    assert sp.divisors(5) == [1, 5]
    return {
        "assumptions": (
            "r>=2, d>=2, 0<=h<=d*r-2, B(t) belongs to "
            "K[y,y^(-1)][t^d], and the required bracket has nonzero "
            "x-degree-h leading coefficient"
        ),
        "first_visible_row": "k0=d*r-h-1",
        "series_indices": {
            "main": "d*(3*r-1)-(h+1)",
            "lambda": "d*(r-1)-(h+1)",
        },
        "theorem": (
            "if d does not divide h+1, the first visible F coefficient "
            "vanishes identically on the congruence section, so that "
            "section cannot realize a nonzero bracket of x-degree h"
        ),
        "certified_r3_X4_test": (
            "h=4 gives h+1=5; both d=2 and d=3 congruence sections are "
            "excluded because neither 2 nor 3 divides 5"
        ),
        "all_r_X4_consequence": (
            "for any r and integer d>=2 with d*r>=6, an X^4 congruence "
            "section can escape this gate only at d=5"
        ),
        "scope": (
            "this eliminates the congruence sections inside the conditional "
            "common-power charts, not the full d=2 or d=3 residues"
        ),
        "checked_compiler_cases": checked_cases,
    }


def verify_r3_d2_residue_geometry(
    residue: CompactResidue,
) -> dict[str, object]:
    """Exhibit the rational fourfold and a smooth torus point for d=2."""

    lam = sp.Symbol("lam")
    p2 = residue.polynomial_symbols_by_q[2]
    p3 = residue.polynomial_symbols_by_q[3]
    p4 = residue.polynomial_symbols_by_q[4]
    p5 = residue.polynomial_symbols_by_q[5]
    p6 = residue.polynomial_symbols_by_q[6]
    numerator = (
        7 * p2**4 * p3
        - 12 * p2**3 * p5
        - 36 * p2**2 * p3 * p4
        - 12 * p2 * p3**3
        + 54 * p2 * p4 * p5
        + 27 * p3**2 * p5
        + 27 * p3 * p4**2
    )
    p6_parameterization = sp.factor(
        -numerator / (54 * (p2 * p3 - 3 * p5))
    )
    assert sp.cancel(
        residue.essential_residue_block[0].subs(p6, p6_parameterization)
    ) == 0

    point = {
        p2: sp.Integer(1),
        p3: sp.Integer(1),
        p4: sp.Integer(1),
        p5: sp.Integer(1),
        p6: sp.Rational(55, 108),
    }
    assert point[p6] == p6_parameterization.subs(point)
    assert residue.essential_residue_block[0].subs(point) == 0
    lambda_value = sp.factor(residue.lambda_substitution.subs(point))
    f_values = {
        k: sp.factor(
            -residue.core_block[k - 1].subs(lam, lambda_value).subs(point)
        )
        for k in residue.f_symbols_by_k
    }
    variables = [residue.polynomial_symbols_by_q[q] for q in range(2, 7)]
    jacobian = sp.Matrix(residue.essential_residue_block).jacobian(
        variables
    ).subs(point)
    assert jacobian.rank() == 1
    assert lambda_value == sp.Rational(7525, 34992)
    assert f_values == {
        3: sp.Rational(25, 243),
        4: -sp.Rational(85, 2916),
        5: -sp.Rational(701, 8748),
    }
    print("F2_R3_D2_RESIDUE_RATIONAL_FOURFOLD_PASS")
    return {
        "geometry": (
            "on p2*p3-3*p5!=0 the one residue equation solves "
            "rationally for p6, so the residue is a rational fourfold"
        ),
        "p6_parameterization": str(p6_parameterization),
        "smooth_torus_point": {
            "p2_to_p6": ["1", "1", "1", "1", "55/108"],
            "lambda": str(lambda_value),
            "f3_to_f5": [str(f_values[k]) for k in (3, 4, 5)],
            "jacobian_rank": jacobian.rank(),
            "local_dimension": 4,
        },
    }


def verify_r3_d3_laurent_bracket_gate() -> dict[str, object]:
    """Classify a monomial bracket on the cubic congruence stratum.

    The reduced bracket determines ``F_-6'``, rather than ``F_-6`` itself.
    This audit keeps the resulting integration constant and proves that it
    cannot occur.  Thus the earlier single-monomial support assumption on
    ``F_-6`` is unnecessary: it is enough that the nonzero bracket
    coefficient is a Laurent monomial.
    """

    n, p, q = sp.symbols("N p q", integer=True)
    alpha, beta, z = sp.symbols("alpha beta z", nonzero=True)

    # If f=alpha*y^N+beta with N,alpha,beta nonzero, then the Euler
    # derivative gives a Bezout certificate for squarefreeness:
    #
    #     N*f-y*f' = N*beta.
    #
    # Here z abbreviates y^N.  Hence A^2|f forces A to be a Laurent unit.
    squarefree_certificate = sp.expand(
        n * (alpha * z + beta) - n * alpha * z
    )
    assert squarefree_certificate == n * beta

    # Write A=a*y^p.  If beta is nonzero, then
    # D=b*y^(N-2p)+c*y^(-2p), with b,c nonzero.  The three terms of D^2
    # have the following pairwise-distinct exponents when N!=0.
    d_square_exponents = [2 * n - 4 * p, n - 4 * p, -4 * p]
    assert [
        sp.expand(d_square_exponents[i] - d_square_exponents[j])
        for i, j in ((0, 1), (1, 2), (0, 2))
    ] == [n, n, 2 * n]

    # No two of those three exponents can both be zero unless N=p=0.
    # The determinants are for their coefficient vectors in (N,p).
    exponent_vectors = [(2, -4), (1, -4), (0, -4)]
    zero_pair_determinants = [
        left[0] * right[1] - left[1] * right[0]
        for left, right in (
            (exponent_vectors[0], exponent_vectors[1]),
            (exponent_vectors[0], exponent_vectors[2]),
            (exponent_vectors[1], exponent_vectors[2]),
        )
    ]
    assert zero_pair_determinants == [-4, -8, -4]

    # The single A^3 exponent 3p can cancel at most one of the three
    # distinct D^2 exponents.  Constancy would force the other two to be
    # zero, contradicting the determinant check.  Therefore beta=0.

    # With beta=0, A^2*D is a Laurent unit, so A=a*y^p and D=b*y^q.
    # Constancy in the nonconstant case forces equal exponents and
    # coefficient cancellation.  Euclid's lemma gives the unique lattice
    # ray (p,q,N)=(2*l,3*l,7*l).
    ell = sp.Symbol("ell", integer=True, nonzero=True)
    assert sp.gcd(2, 3) == 1
    assert sp.solve(
        [3 * p - 2 * q, 2 * p + q - n], (p, q), dict=True
    ) == [{p: 2 * n / 7, q: 3 * n / 7}]
    assert sp.expand(3 * (2 * ell) - 2 * (3 * ell)) == 0
    assert sp.expand(2 * (2 * ell) + 3 * ell - 7 * ell) == 0

    # The exceptional ray is sharp at the invariant level.  Taking
    # u=0, v=y^(2*l), w=y^(3*l)/3 gives A=-3*y^(2*l), D=9*y^(3*l),
    # lambda=0 and a nonconstant F_-6=(5/81)*y^(7*l).
    a_witness = -3 * z**2
    d_witness = 9 * z**3
    lambda_witness = sp.factor(
        -sp.Rational(5, 6561) * (3 * a_witness**3 + d_witness**2)
    )
    f6_witness = sp.factor(
        sp.Rational(5, 6561) * a_witness**2 * d_witness
    )
    assert lambda_witness == 0
    assert f6_witness == sp.Rational(5, 81) * z**7
    assert (2 + 1) % 7 == 3  # historical d=3 bracket weight y^2
    assert (0 + 1) % 7 == 1  # constant leading bracket coefficient
    assert (6 + 1) % 7 == 0  # first exceptional weight

    print("F2_R3_D3_LAURENT_BRACKET_GATE_PASS")
    return {
        "hypothesis": (
            "on the cubic congruence stratum, lambda is constant and the "
            "nonzero bracket is mu*y^s*x^2"
        ),
        "antiderivative": (
            "s=-1 is impossible in K[y,y^(-1)]; otherwise "
            "F_-6=c*y^(s+1)+c0 with c!=0"
        ),
        "integration_constant": {
            "verdict": "c0!=0 is impossible",
            "certificate": "N*f-y*f'=N*c0, so f is squarefree",
            "exponent_obstruction": (
                "after A is forced to be a Laurent unit, the three "
                "pairwise-distinct D^2 exponents cannot be made constant "
                "by the single A^3 exponent"
            ),
            "zero_pair_determinants": zero_pair_determinants,
        },
        "classification": (
            "every survivor has s+1=7*ell with ell!=0, lambda=0, "
            "A=a*y^(2*ell), D=b*y^(3*ell), and 3*a^3+b^2=0"
        ),
        "sharp_invariant_witness": {
            "u": "0",
            "v": "y^(2*ell)",
            "w": "y^(3*ell)/3",
            "A": "-3*y^(2*ell)",
            "D": "9*y^(3*ell)",
            "lambda": "0",
            "f6": "5*y^(7*ell)/81",
            "warning": (
                "this is a sharp witness for the finite invariant gate, "
                "not a reconstructed polynomial Keller map"
            ),
        },
        "consequence": (
            "the single-monomial support assumption on F_-6 is removed; "
            "a monomial bracket coefficient excludes the stratum unless "
            "s=7*ell-1 and lambda=0"
        ),
        "tested_weights": {
            "s=2_historical_d3_weight": "excluded because s+1=3",
            "s=0_constant_leading_coefficient": "excluded because s+1=1",
            "s=6_first_exceptional_weight": (
                "not excluded by this gate; requires lambda=0"
            ),
        },
    }


def verify_r3_d3_residue_geometry(
    residue: CompactResidue,
) -> dict[str, object]:
    """Reduce the congruence stratum and certify a smooth torus branch."""

    lam = sp.Symbol("lam")
    u, v, w = sp.symbols("u v w")
    congruence = {
        symbol: sp.Integer(0)
        for symbol in residue.polynomial_symbols_by_q.values()
    }
    congruence.update(
        {
            residue.polynomial_symbols_by_q[3]: u,
            residue.polynomial_symbols_by_q[6]: v,
            residue.polynomial_symbols_by_q[9]: w,
        }
    )
    assert all(
        sp.expand(equation.subs(congruence)) == 0
        for equation in residue.essential_residue_block
    )
    cubic_a = u**2 - 3 * v
    cubic_b = 2 * u**3 - 9 * u * v + 27 * w
    expected_lambda = sp.factor(
        -sp.Rational(5, 6561) * (3 * cubic_a**3 + cubic_b**2)
    )
    actual_lambda = sp.factor(residue.lambda_substitution.subs(congruence))
    assert sp.expand(actual_lambda - expected_lambda) == 0
    f_values = {
        k: sp.factor(
            -residue.core_block[k - 1]
            .subs(lam, actual_lambda)
            .subs(congruence)
        )
        for k in residue.f_symbols_by_k
    }
    expected_f6 = sp.factor(sp.Rational(5, 6561) * cubic_a**2 * cubic_b)
    assert sp.expand(f_values[6] - expected_f6) == 0
    assert f_values[7] == 0
    assert f_values[8] == 0

    # Compute the polynomial part of b(X)^(5/3) at X=infinity and verify its
    # exact Jacobian identity.  Here X=x^3 in the d=3 chart.
    maximum_reduced_degree = 8
    reduced_b = [sp.Integer(0)] * (maximum_reduced_degree + 1)
    reduced_b[0] = 1
    reduced_b[1] = u
    reduced_b[2] = v
    reduced_b[3] = w
    reduced_inverse_root = _inverse_rth_root(
        reduced_b, 3, maximum_reduced_degree
    )
    reduced_m = _mul(
        _power(reduced_b, 2, maximum_reduced_degree),
        reduced_inverse_root,
        maximum_reduced_degree,
    )
    X = sp.Symbol("X")
    reduced_p = X**3 + u * X**2 + v * X + w
    reduced_q = sum(
        reduced_m[index] * X ** (5 - index) for index in range(6)
    )
    u_prime, v_prime, w_prime = sp.symbols(
        "u_prime v_prime w_prime"
    )

    def y_derivative(expression: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(expression, u) * u_prime
            + sp.diff(expression, v) * v_prime
            + sp.diff(expression, w) * w_prime
        )

    reduced_bracket = sp.expand(
        sp.diff(reduced_p, X) * y_derivative(reduced_q)
        - y_derivative(reduced_p) * sp.diff(reduced_q, X)
    )
    assert sp.expand(
        reduced_bracket
        - (3 * X + u) * y_derivative(expected_lambda)
        - 3 * y_derivative(expected_f6)
    ) == 0

    # The rational point u=v=w=1 is smooth of dimension four, has p9!=0,
    # and keeps lambda and the first bracket-visible F coefficient nonzero.
    base_point = {
        symbol: sp.Integer(0)
        for symbol in residue.polynomial_symbols_by_q.values()
    }
    base_point.update(
        {
            residue.polynomial_symbols_by_q[3]: 1,
            residue.polynomial_symbols_by_q[6]: 1,
            residue.polynomial_symbols_by_q[9]: 1,
        }
    )
    variables = [residue.polynomial_symbols_by_q[q] for q in range(2, 10)]
    jacobian = sp.Matrix(residue.essential_residue_block).jacobian(
        variables
    ).subs(base_point)
    pivot_q = [2, 4, 5, 7]
    pivot_columns = [q - 2 for q in pivot_q]
    pivot_minor = sp.factor(jacobian[:, pivot_columns].det())
    assert pivot_minor != 0
    assert jacobian.rank() == 4
    assert actual_lambda.subs({u: 1, v: 1, w: 1}) == -sp.Rational(
        1880, 6561
    )
    assert expected_f6.subs({u: 1, v: 1, w: 1}) == sp.Rational(
        400, 6561
    )

    # Fix p3=p6=p9=1 and p8=eps.  The pivot minor gives a unique formal
    # solution for p2,p4,p5,p7.  Its exact two-jet puts every coordinate in
    # the nonzero torus over Q((eps)).
    eps = sp.Symbol("eps")
    formal_two_jet = {
        residue.polynomial_symbols_by_q[2]: 3 * eps,
        residue.polynomial_symbols_by_q[3]: 1,
        residue.polynomial_symbols_by_q[4]: 3 * eps**2,
        residue.polynomial_symbols_by_q[5]: 2 * eps,
        residue.polynomial_symbols_by_q[6]: 1,
        residue.polynomial_symbols_by_q[7]: eps**2,
        residue.polynomial_symbols_by_q[8]: eps,
        residue.polynomial_symbols_by_q[9]: 1,
    }
    for equation in residue.essential_residue_block:
        expanded = sp.expand(equation.subs(formal_two_jet))
        assert all(expanded.coeff(eps, order) == 0 for order in range(3))

    laurent_bracket_gate = verify_r3_d3_laurent_bracket_gate()

    # The first exceptional exponent (ell=1) is an exact solution of every
    # row in the compiled h=2 window, not merely of the eliminated invariant
    # relation.  It is still only a finite residual candidate because the
    # lower F tail and the bridge back to a polynomial Keller map are absent.
    y = sp.Symbol("y", nonzero=True)
    sharp_candidate = {
        symbol: sp.Integer(0)
        for symbol in residue.polynomial_symbols_by_q.values()
    }
    sharp_candidate.update(
        {
            residue.polynomial_symbols_by_q[6]: y**2,
            residue.polynomial_symbols_by_q[9]: y**3 / 3,
            lam: sp.Integer(0),
            residue.f_symbols_by_k[6]: sp.Rational(5, 81) * y**7,
            residue.f_symbols_by_k[7]: sp.Integer(0),
            residue.f_symbols_by_k[8]: sp.Integer(0),
        }
    )
    assert all(
        sp.expand(equation.subs(sharp_candidate)) == 0
        for equation in residue.fiber_block
    )
    sharp_leading_bracket = sp.factor(
        9 * sp.diff(sharp_candidate[residue.f_symbols_by_k[6]], y)
    )
    assert sharp_leading_bracket == sp.Rational(35, 9) * y**6

    print("F2_R3_D3_RESIDUE_SMOOTH_TORUS_BRANCH_PASS")
    print("F2_R3_D3_CUBIC_INVARIANT_REDUCTION_PASS")
    return {
        "congruence_stratum": "B(t)=1+u*t^3+v*t^6+w*t^9",
        "cubic_invariants": {
            "A": "u^2-3*v",
            "D": "2*u^3-9*u*v+27*w",
            "lambda": "-5*(3*A^3+D^2)/6561",
            "f6": "5*A^2*D/6561",
            "f7": "0",
            "f8": "0",
        },
        "reduced_bracket": (
            "for Pbar=X^3+u*X^2+v*X+w and Qbar the polynomial part "
            "of Pbar^(5/3), [Pbar,Qbar]=(3*X+u)*lambda' + 3*f6'"
        ),
        "laurent_bracket_gate": laurent_bracket_gate,
        "sharp_finite_residual_candidate": {
            "ell": 1,
            "B": "1+y^2*t^6+(y^3/3)*t^9",
            "lambda": "0",
            "visible_F": {
                "f6": "5*y^7/81",
                "f7": "0",
                "f8": "0",
            },
            "leading_bracket": "(35/9)*y^6*x^2",
            "verified_scope": "all eight compiled h=2 fiber rows",
            "warning": (
                "the lower F tail, exact Laurent support ledger, and "
                "polynomial-coordinate bridge are not reconstructed"
            ),
        },
        "smooth_base_point": {
            "nonzero_coordinates": {"p3": 1, "p6": 1, "p9": 1},
            "lambda": "-1880/6561",
            "f6": "400/6561",
            "pivot_variables": [f"p{q}" for q in pivot_q],
            "pivot_minor": str(pivot_minor),
            "jacobian_rank": jacobian.rank(),
            "local_dimension": 4,
        },
        "formal_coefficient_torus_branch": {
            "free_coordinates": "p3=p6=p9=1, p8=eps",
            "two_jet": (
                "p2=3*eps+O(eps^3), p4=3*eps^2+O(eps^3), "
                "p5=2*eps+O(eps^3), p7=eps^2+O(eps^3)"
            ),
            "field": "Q((eps))",
            "verdict": (
                "all p2,...,p9, lambda, and the first visible F "
                "coefficient are nonzero on the formal branch"
            ),
        },
    }


def verify_r3_residue_geometry(
    compact_residues: list[CompactResidue],
) -> dict[str, object]:
    by_case = {
        (
            residue.summary.r,
            residue.summary.d,
            residue.summary.jacobian_x_degree,
        ): residue
        for residue in compact_residues
    }
    return {
        "d2_h2": verify_r3_d2_residue_geometry(by_case[(3, 2, 2)]),
        "d3_h2": verify_r3_d3_residue_geometry(by_case[(3, 3, 2)]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--include-equations",
        action="store_true",
        help="include canonical residual equations in the JSON output",
    )
    args = parser.parse_args()

    verify_published_r2_regression()
    systems = [generate_system(r, d) for r in (2, 3) for d in (2, 3)]
    compact_residues = [
        generate_compact_residue(system, jacobian_x_degree=2)
        for system in systems
    ]
    compact_residues.extend(
        generate_compact_residue(system, jacobian_x_degree=4)
        for system in systems
        if system.summary.r == 3
    )
    verify_r3_residue_presentations(compact_residues)
    endpoint_binomial_branch = endpoint_binomial_branch_audit(
        compact_residues
    )
    congruence_support_bracket = congruence_support_bracket_audit(
        compact_residues
    )
    r3_residue_geometry = verify_r3_residue_geometry(compact_residues)
    payload: dict[str, object] = {
        "schema": "plane-jc.f2-modified-laurent-family.v3",
        "claim_boundary": (
            "exact formal systems conditional on the modified common-power chart; "
            "the F2 table derives the degree family, but the corner chain does not "
            "derive d=2/3, the lower y-supports, or a three-visible-row r=3 F-tail"
        ),
        "family_formulas": family_formula_audit(),
        "endpoint_binomial_branch": endpoint_binomial_branch,
        "congruence_support_bracket": congruence_support_bracket,
        "r3_residue_geometry": r3_residue_geometry,
        "systems": [asdict(system.summary) for system in systems],
        "compact_residues": [
            asdict(residue.summary) for residue in compact_residues
        ],
        "reproduction_command": (
            ".venv/bin/python plane-jc/cas/generate_f2_modified_system.py "
            "--include-equations --output "
            "artifacts/generated-results/jc2_f2_modified_laurent_family.json"
        ),
        "software_assumptions": [
            ".python-version",
            "requirements.txt",
            "characteristic-zero exact SymPy arithmetic",
        ],
    }
    if args.include_equations:
        payload["residual_equations"] = {
            f"r{system.summary.r}_d{system.summary.d}": [
                str(sp.factor(equation)) for equation in system.residual_block
            ]
            for system in systems
        }
        payload["lambda_eliminated_residue_equations"] = {
            (
                f"r{residue.summary.r}_d{residue.summary.d}_"
                f"h{residue.summary.jacobian_x_degree}"
            ): [
                str(sp.factor(equation))
                for equation in residue.lambda_eliminated_block
            ]
            for residue in compact_residues
        }
        payload["essential_residue_fitting_generators"] = {
            (
                f"r{residue.summary.r}_d{residue.summary.d}_"
                f"h{residue.summary.jacobian_x_degree}"
            ): [
                str(sp.factor(equation))
                for equation in residue.essential_residue_block
            ]
            for residue in compact_residues
        }
        payload["artinian_maximal_minor_generators"] = {
            (
                f"r{residue.summary.r}_d{residue.summary.d}_"
                f"h{residue.summary.jacobian_x_degree}"
            ): [
                str(sp.factor(equation))
                for equation in residue.artinian_remainder_block
            ]
            for residue in compact_residues
        }
        payload["lambda_substitutions"] = {
            (
                f"r{residue.summary.r}_d{residue.summary.d}_"
                f"h{residue.summary.jacobian_x_degree}"
            ): str(sp.factor(residue.lambda_substitution))
            for residue in compact_residues
        }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
        print(f"WROTE {args.output}")
    else:
        print(encoded, end="")
    print("F2_75_125_D2_SYSTEM_PASS")
    print("F2_75_125_D3_SYSTEM_PASS")
    print("F2_MODIFIED_POWER_COORDINATE_PASS")
    print("F2_MODIFIED_TOEPLITZ_FITTING_PASS")
    print("F2_MODIFIED_ENDPOINT_BINOMIAL_BRANCH_PASS")
    print("F2_MODIFIED_CONGRUENCE_SUPPORT_BRACKET_PASS")
    print("F2_MODIFIED_CHAIN_GAP_AUDIT_PASS")
    print("F2_MODIFIED_SYSTEM_FRONTEND_PASS")


if __name__ == "__main__":
    main()
