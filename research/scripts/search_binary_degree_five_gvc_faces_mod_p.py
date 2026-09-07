#!/usr/bin/env python3
"""Exhaust the two residual cubic-leading quintic faces over finite fields.

This is a finite-field experiment supporting the exact characteristic-zero
face calculations in ``verify_binary_degree_five_gvc_frontier.py``.  It also
checks every squarefree quartic cross-ratio and every projective quintic top
form over seven fields.  Linear/triangular equations are eliminated before
enumeration when possible.

The search covers the residual ``(d_x^2 d_y, x y^4)`` and
``(d_x^2 d_y, y^5)`` faces modulo 101, 103, and 107.  It records counts after
each pure moment and quotients the final nonzero survivors by the residual
diagonal scaling.  A bounded modular search is evidence, not an all-order
characteristic-zero proof.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "binary_degree_five_gvc_face_search.json"
)
PRIMES = (101, 103, 107)
TOP_PRIMES = (11, 13, 17, 19, 23, 29, 31)

x, y = sp.symbols("x y")
A, B, H, J = sp.symbols("A B H J")


def apply_operator(
    polynomial: sp.Expr,
    terms: dict[tuple[int, int], sp.Expr],
) -> sp.Expr:
    return sp.expand(
        sum(
            coefficient * sp.diff(polynomial, x, i, y, j)
            for (i, j), coefficient in terms.items()
        )
    )


def moment(
    polynomial: sp.Expr,
    terms: dict[tuple[int, int], sp.Expr],
    order: int,
) -> sp.Poly:
    result = sp.expand(polynomial**order)
    for _ in range(order):
        result = apply_operator(result, terms)
    return sp.Poly(result, x, y)


def primitive_coefficients(
    polynomial: sp.Poly,
    substitution: dict[sp.Symbol, sp.Expr],
    parameters: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, ...]:
    result = []
    for coefficient in polynomial.coeffs():
        specialized = sp.expand(coefficient.subs(substitution))
        if specialized == 0:
            continue
        _content, primitive = sp.Poly(
            specialized, *parameters, domain=sp.QQ
        ).primitive()
        result.append(sp.expand(primitive.as_expr()))
    return tuple(result)


def xy4_equations() -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]]:
    operator = {(2, 1): 1, (1, 3): H, (0, 5): J}
    polynomial = x * y**4 + A * x**2 * y**2 + B * x**3
    first_branch = {A: -6 * H}
    second = primitive_coefficients(
        moment(polynomial, operator, 2),
        first_branch,
        (B, H, J),
    )
    third = primitive_coefficients(
        moment(polynomial, operator, 3),
        first_branch,
        (B, H, J),
    )
    assert set(map(sp.expand, second)) == set(map(sp.expand, {
        B - 2 * H**2 + 140 * J,
        -H * (B + 18 * H**2 - 100 * J),
    }))
    assert set(map(sp.expand, third)) == set(map(sp.expand, {
        7 * B**2
        + 12 * B * H**2
        + 840 * B * J
        - 6480 * H**4
        + 15120 * H**2 * J
        + 166320 * J**2,
        H * (B - 15 * H**2 + 504 * J),
    }))
    return second, third


def y5_equations() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    operator = {(2, 1): 1, (1, 3): H, (0, 5): J}
    polynomial = y**5 + A * x * y**3 + B * x**2 * y
    first_branch = {B: -3 * A * H - 60 * J}
    equations = []
    for order in (2, 3, 4):
        coefficients = primitive_coefficients(
            moment(polynomial, operator, order),
            first_branch,
            (A, H, J),
        )
        assert len(coefficients) == 1
        equations.append(coefficients[0])
    return tuple(equations)  # type: ignore[return-value]


def evaluate(expression: sp.Expr, values: tuple[int, ...], prime: int) -> int:
    parameters = tuple(sorted(expression.free_symbols, key=str))
    polynomial = sp.Poly(expression, *parameters)
    assignment = dict(zip(parameters, values, strict=True))
    return int(polynomial.eval(assignment)) % prime


def evaluate_bhj(expression: sp.Expr, b: int, h: int, j: int, p: int) -> int:
    return int(expression.subs({B: b, H: h, J: j})) % p


@lru_cache(maxsize=None)
def ahj_terms(
    expression: sp.Expr,
) -> tuple[tuple[tuple[int, int, int], int], ...]:
    return tuple(
        (degrees, int(coefficient))
        for degrees, coefficient in sp.Poly(expression, A, H, J).terms()
    )


def evaluate_ahj(expression: sp.Expr, a: int, h: int, j: int, p: int) -> int:
    # Expanded direct evaluation is substantially faster than SymPy
    # substitution in the three-million-point loop.
    total = 0
    for (a_degree, h_degree, j_degree), coefficient in ahj_terms(expression):
        total += (
            coefficient
            * pow(a, a_degree, p)
            * pow(h, h_degree, p)
            * pow(j, j_degree, p)
        )
    return total % p


def compiled_terms(
    expressions: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> tuple[tuple[tuple[tuple[int, ...], int], ...], ...]:
    return tuple(
        tuple(
            (degrees, int(coefficient))
            for degrees, coefficient in sp.Poly(
                expression, *variables
            ).terms()
        )
        for expression in expressions
    )


def evaluate_compiled(
    polynomials: tuple[tuple[tuple[int, ...], int], ...],
    values: tuple[int, ...],
    prime: int,
) -> bool:
    powers = [
        [pow(value, exponent, prime) for exponent in range(7)]
        for value in values
    ]
    for polynomial in polynomials:
        total = 0
        for degrees, coefficient in polynomial:
            term = coefficient
            for index, degree in enumerate(degrees):
                term *= powers[index][degree]
            total += term
        if total % prime:
            return False
    return True


def quartic_top_equations() -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]]:
    lam = sp.symbols("lam")
    top_coefficients = sp.symbols("q0:6")
    polynomial = sum(
        top_coefficients[index] * x ** (5 - index) * y**index
        for index in range(6)
    )
    operator = {
        (3, 1): 1,
        (2, 2): -(lam + 1),
        (1, 3): lam,
    }
    first_branch = {
        top_coefficients[1]: (
            (lam + 1) * top_coefficients[2]
            - lam * top_coefficients[3]
        )
        / 2,
        top_coefficients[4]: (
            (lam + 1) * top_coefficients[3]
            - top_coefficients[2]
        )
        / (2 * lam),
    }
    by_order = []
    parameters = (
        lam,
        top_coefficients[0],
        top_coefficients[2],
        top_coefficients[3],
        top_coefficients[5],
    )
    for order in (2, 3):
        contracted = sp.Poly(
            moment(polynomial.subs(first_branch), operator, order),
            x,
            y,
        )
        equations = []
        for coefficient in contracted.coeffs():
            numerator = sp.together(coefficient).as_numer_denom()[0]
            _content, primitive = sp.Poly(
                numerator, *parameters, domain=sp.QQ
            ).primitive()
            equations.append(sp.expand(primitive.as_expr()))
        by_order.append(tuple(equations))
    assert tuple(map(len, by_order)) == (3, 4)
    return tuple(by_order)  # type: ignore[return-value]


def projective_four_tuples(prime: int):
    for a2 in range(prime):
        for a3 in range(prime):
            for a5 in range(prime):
                yield (1, a2, a3, a5)
    for a3 in range(prime):
        for a5 in range(prime):
            yield (0, 1, a3, a5)
    for a5 in range(prime):
        yield (0, 0, 1, a5)
    yield (0, 0, 0, 1)


def normalize_projective(
    coordinates: tuple[int, ...],
    prime: int,
) -> tuple[int, ...]:
    first = next(value for value in coordinates if value % prime)
    inverse = pow(first, -1, prime)
    return tuple(value * inverse % prime for value in coordinates)


def search_quartic_top(
    prime: int,
    second: tuple[sp.Expr, ...],
    third: tuple[sp.Expr, ...],
) -> dict[str, object]:
    lam = sp.symbols("lam")
    q0, _q1, q2, q3, _q4, q5 = sp.symbols("q0:6")
    variables = (lam, q0, q2, q3, q5)
    compiled_second = compiled_terms(second, variables)
    compiled_third = compiled_terms(third, variables)
    total_projective = prime**3 + prime**2 + prime + 1
    after_second = 0
    after_third = 0
    for lam_value in range(2, prime):
        expected = {
            normalize_projective((1, 0, 0, 0), prime),
            normalize_projective((0, 0, 0, 1), prime),
            normalize_projective((1, 10, 10, 1), prime),
            normalize_projective(
                (
                    pow(lam_value, 5, prime),
                    10 * pow(lam_value, 3, prime),
                    10 * pow(lam_value, 2, prime),
                    1,
                ),
                prime,
            ),
        }
        survivors = set()
        for coordinates in projective_four_tuples(prime):
            values = (lam_value, *coordinates)
            if not evaluate_compiled(compiled_second, values, prime):
                continue
            after_second += 1
            if evaluate_compiled(compiled_third, values, prime):
                after_third += 1
                survivors.add(coordinates)
        assert survivors == expected
    return {
        "prime": prime,
        "squarefree_cross_ratios": prime - 2,
        "projective_forms_per_cross_ratio": total_projective,
        "raw_projective_tuples": (prime - 2) * total_projective,
        "after_moment_2": after_second,
        "after_moment_3": after_third,
        "survivors_per_cross_ratio": 4,
        "survivor_description": (
            "x^5, y^5, (x+y)^5, (lambda*x+y)^5"
        ),
    }


def search_xy4(
    prime: int,
    second: tuple[sp.Expr, ...],
    third: tuple[sp.Expr, ...],
) -> dict[str, object]:
    # The first second-moment equation gives B=2H^2-140J.  Counting this
    # triangular chart is exactly the same as checking all p^3 triples.
    survivors_second: list[tuple[int, int, int]] = []
    for h in range(prime):
        for j in range(prime):
            b = (2 * h * h - 140 * j) % prime
            if all(
                evaluate_bhj(equation, b, h, j, prime) == 0
                for equation in second
            ):
                survivors_second.append((b, h, j))
    survivors_third = [
        candidate
        for candidate in survivors_second
        if all(
            evaluate_bhj(equation, *candidate, prime) == 0
            for equation in third
        )
    ]
    assert survivors_third == [(0, 0, 0)]
    return {
        "prime": prime,
        "raw_triples": prime**3,
        "after_moment_2": len(survivors_second),
        "after_moment_3": len(survivors_third),
        "final_survivors": [list(candidate) for candidate in survivors_third],
        "nonzero_stabilizer_orbits": 0,
    }


def search_y5(
    prime: int,
    equations: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> dict[str, object]:
    second, third, fourth = equations
    survivors_second = []
    for a in range(prime):
        for h in range(prime):
            for j in range(prime):
                if evaluate_ahj(second, a, h, j, prime) == 0:
                    survivors_second.append((a, h, j))
    survivors_third = [
        candidate
        for candidate in survivors_second
        if evaluate_ahj(third, *candidate, prime) == 0
    ]
    survivors_fourth = [
        candidate
        for candidate in survivors_third
        if evaluate_ahj(fourth, *candidate, prime) == 0
    ]
    expected = {
        (a, 0, 0) for a in range(prime)
    } | {
        (0, h, 0) for h in range(prime)
    }
    assert set(survivors_fourth) == expected
    return {
        "prime": prime,
        "raw_triples": prime**3,
        "after_moment_2": len(survivors_second),
        "after_moment_3": len(survivors_third),
        "after_moment_4": len(survivors_fourth),
        "final_survivor_description": (
            "{(A,0,0)} union {(0,H,0)}"
        ),
        "final_survivor_count": len(survivors_fourth),
        "nonzero_stabilizer_orbits": 2,
        "nullcone_screen": "two one-sided coordinate lines",
    }


def main() -> None:
    xy4_second, xy4_third = xy4_equations()
    y5_pure_equations = y5_equations()
    quartic_second, quartic_third = quartic_top_equations()
    quartic_records = [
        search_quartic_top(prime, quartic_second, quartic_third)
        for prime in TOP_PRIMES
    ]
    xy4_records = [
        search_xy4(prime, xy4_second, xy4_third) for prime in PRIMES
    ]
    y5_records = [
        search_y5(prime, y5_pure_equations) for prime in PRIMES
    ]
    artifact = {
        "format": "binary-degree-five-gvc-face-search-v1",
        "status": "experiment",
        "primes": list(PRIMES),
        "quartic_top_primes": list(TOP_PRIMES),
        "parameter_order": {
            "xy4": ["B", "H", "J"],
            "y5": ["A", "H", "J"],
        },
        "xy4": xy4_records,
        "y5": y5_records,
        "quartic_squarefree_top_forms": quartic_records,
        "total_raw_triples": sum(
            record["raw_triples"]
            for record in xy4_records + y5_records
        ),
        "total_raw_projective_top_forms": sum(
            record["raw_projective_tuples"] for record in quartic_records
        ),
        "scope": (
            "exhaustive on the two displayed weighted faces and on every "
            "squarefree quartic cross-ratio/projective top form over the "
            "listed finite fields; exact identities are checked first, but "
            "the modular counts are not an all-order proof"
        ),
    }
    assert artifact["total_raw_triples"] == 6_696_142
    assert artifact["total_raw_projective_top_forms"] == 2_082_612
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS xy4 face: only the origin survives moment three")
    print("PASS y5 face: exactly two one-sided lines survive moment four")
    print("PASS enumerated 6,696,142 triples over three primes")
    print("PASS enumerated 2,082,612 quartic top forms over seven primes")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
