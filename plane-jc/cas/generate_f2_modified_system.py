#!/usr/bin/env python3
"""Generate the reduced F2 modified Laurent systems.

This is a formal-algebra front end, not yet a proof that the F2 corner chain
forces the displayed normal form for every member of the family.  It assumes

    P = C^r,
    Q = C^(2r-1) + lambda*C^(-1) + F,

with leading x-degree d for C and the same three-term terminal F block as in
the published r=2 calculation.  Under that hypothesis it reconstructs the
finite coefficient system exactly, triangularly removes the coefficients
forced by polynomiality of P, and audits the resulting weighted residual
system.

The r=2, d=2 and d=3 outputs regress the two published (50,75) modified
systems.  The r=3 outputs are the 14- and 22-equation candidate systems for
the first F2 frontier pair (75,125).
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


def _canonical_polynomial(expression: sp.Expr, generators: list[sp.Symbol]) -> list:
    polynomial = sp.Poly(sp.expand(expression), *generators, domain=sp.QQ)
    return [
        [list(monomial), [int(coefficient.p), int(coefficient.q)]]
        for monomial, coefficient in polynomial.terms()
    ]


def _digest(polynomials: Iterable[list]) -> str:
    payload = json.dumps(list(polynomials), separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


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
    reduced_pair: tuple[int, int]
    degree_pair: tuple[int, int]
    first_block_equations: int
    residual_equations: int
    total_equations: int
    lowest_z_exponent: int
    residual_z_exponents: list[int]
    residual_weights: list[int]
    lambda_weight: int
    lambda_elimination_equation: int
    zero_residual_equations: list[int]
    terminal_f_equations: list[int]
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
        reduced_pair=(r, m),
        degree_pair=(25 * r, 25 * m),
        first_block_equations=first_count,
        residual_equations=residual_count,
        total_equations=total_count,
        lowest_z_exponent=1 - d * (3 * r - 2),
        residual_z_exponents=[d - q for q in range(2, d * r + 1)],
        residual_weights=residual_weights,
        lambda_weight=2 * d * r,
        lambda_elimination_equation=d,
        zero_residual_equations=list(range(1, max(1, d * r - 3))),
        terminal_f_equations=list(range(max(1, d * r - 3), d * r)),
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
    assert gamma3.summary.lowest_z_exponent == -7

    gamma2 = generate_system(2, 3)
    z = gamma2.symbols_by_q
    lam = gamma2.lambda_symbol
    assert gamma2.first_block[0] == 2 * z[3] * z[4] + 2 * z[2] * z[5] + 2 * z[7]
    assert gamma2.first_block[-1].coeff(z[14]) == 2
    assert gamma2.raw_residual_block[2].coeff(lam) == 1
    assert gamma2.summary.total_equations == 13
    assert gamma2.summary.lowest_z_exponent == -11
    print("F2_50_75_PUBLISHED_REGRESSION_PASS")


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
    payload: dict[str, object] = {
        "claim_boundary": (
            "formal systems conditional on the F2 corner-chain normalization; "
            "the r=3 normalization and F-tail ledger remain to be proved"
        ),
        "systems": [asdict(system.summary) for system in systems],
    }
    if args.include_equations:
        payload["residual_equations"] = {
            f"r{system.summary.r}_d{system.summary.d}": [
                str(sp.factor(equation)) for equation in system.residual_block
            ]
            for system in systems
        }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    print("F2_75_125_D2_SYSTEM_PASS")
    print("F2_75_125_D3_SYSTEM_PASS")
    print("F2_MODIFIED_SYSTEM_FRONTEND_PASS")


if __name__ == "__main__":
    main()
