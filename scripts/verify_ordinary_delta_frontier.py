#!/usr/bin/env python3
"""Exact symbolic certificates for low-rank ordinary-Laplacian GVC routes.

Requires SymPy >= 1.12.  All calculations are over Q.
"""
from __future__ import annotations

from collections.abc import Callable
import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring


def repeated_moment(op: Callable[[sp.Expr], sp.Expr], f: sp.Expr, m: int) -> sp.Expr:
    out = sp.expand(f**m)
    for _ in range(m):
        out = sp.expand(op(out))
    return sp.expand(out)


def assert_zero(expr: sp.Expr, label: str) -> None:
    value = sp.expand(expr)
    if value != 0:
        raise AssertionError(f"{label} failed: {sp.factor(value)}")


def verify_meng_yang_trace_obstruction() -> None:
    x1, x2, y1, y2, y3, alpha, beta = sp.symbols(
        "x1 x2 y1 y2 y3 alpha beta"
    )
    variables = (x1, x2, y1, y2, y3)
    u = 1 + x1 * x2
    A = y1 * u**3 + 3 * x1 * y2 * u**2 - x1**3 * y3
    B = (
        y1 * x2**2 * u * (4 + 3 * x1 * x2)
        + y2 * (x2 + 3 * x1 * x2**2 * (4 + 3 * x1 * x2))
        + y3 * (2 * x1 - 3 * x1**2 * x2)
    )
    psi = sp.expand(A**2 + alpha * A + beta * B)
    hessian = [[sp.diff(psi, variables[i], variables[j]) for j in range(5)] for i in range(5)]

    k_symbols: list[sp.Symbol] = []
    K: dict[tuple[int, int], sp.Symbol] = {}
    for i in range(5):
        for j in range(i, 5):
            kij = sp.Symbol(f"k{i + 1}{j + 1}")
            k_symbols.append(kij)
            K[i, j] = kij
            K[j, i] = kij

    trace_pairing = sp.expand(
        sum(K[i, j] * hessian[i][j] for i in range(5) for j in range(5))
    )
    poly = sp.Poly(trace_pairing, *variables)

    selected_exponents = [
        (6, 6, 0, 0, 0),
        (6, 5, 1, 0, 0),
        (6, 5, 0, 0, 0),
        (6, 4, 2, 0, 0),
        (6, 4, 1, 0, 0),
        (6, 4, 0, 0, 0),
        (6, 3, 0, 0, 0),
        (6, 2, 1, 0, 0),
        (6, 2, 0, 0, 0),
        (6, 0, 0, 0, 0),
        (5, 6, 1, 0, 0),
        (5, 5, 2, 0, 0),
        (5, 5, 1, 0, 0),
        (5, 3, 1, 0, 0),
        (4, 6, 2, 0, 0),
    ]
    coefficient_rows = []
    for exponent in selected_exponents:
        coeff = poly.coeff_monomial(exponent)
        coefficient_rows.append([sp.diff(coeff, k) for k in k_symbols])
    rank = sp.Matrix(coefficient_rows).rank()
    if rank != 15:
        raise AssertionError(f"Meng-Yang trace coefficient rank is {rank}, expected 15")

    expected_constant = 2 * (
        3 * alpha * K[0, 3] + 2 * beta * K[0, 4] + beta * K[1, 3] + K[2, 2]
    )
    assert_zero(
        poly.coeff_monomial((0, 0, 0, 0, 0)) - expected_constant,
        "Meng-Yang constant coefficient",
    )
    print("[ok] Meng-Yang Schur family: no nonzero constant Hessian trace pairing")


def verify_seven_variable_square_block() -> None:
    a, b, c, d, t, s, u = sp.symbols("a b c d t s u")
    q0, q1, q2, q3, q4 = sp.symbols("q0 q1 q2 q3 q4")
    p0 = (t + c) * (a * d + b * t)
    f = sp.expand(
        p0
        + q0 * t * u**2
        + q1 * s * t**2
        + q2 * c * u**2
        + q3 * c * s * t
        + q4 * c**2 * s
    )

    def delta(expr: sp.Expr) -> sp.Expr:
        return (
            sp.diff(expr, a, d)
            - sp.diff(expr, b, c)
            + sp.diff(expr, t, s)
            + sp.diff(expr, u, 2)
        )

    m1 = repeated_moment(delta, f, 1)
    assert_zero(
        m1 - (2 * (q0 + q1) * t + (2 * q2 + q3 + 1) * c),
        "7-variable first moment",
    )
    first_substitution = {q0: -q1, q3: -2 * q2 - 1}

    m2 = sp.expand(repeated_moment(delta, f, 2).subs(first_substitution))
    m3 = sp.expand(repeated_moment(delta, f, 3).subs(first_substitution))
    E1 = 6 * q2**2 + 2 * q2 - 2 * q4 + 1
    E2 = 6 * q1 * q2 + q1 - q2 + 2 * q4 - 1
    E3 = 6 * q1**2 + 4 * q2 + 3
    assert_zero(
        m2 - 4 * (E1 * c**2 - 2 * E2 * c * t + E3 * t**2),
        "7-variable second moment",
    )

    gb = sp.groebner([E1, E2, E3], q1, q2, q4, order="lex")
    g = 72 * q4**3 + 66 * q4**2 + 20 * q4 - 25
    if gb.reduce(g)[1] != 0:
        raise AssertionError("7-variable elimination cubic is not in the moment-2 ideal")

    c3 = sp.Poly(m3 / 36, c, t).coeff_monomial(c**3)
    h = (2 * q4 - 1) * (4 * q4 + 5)
    remainder = sp.factor(gb.reduce(c3)[1])
    if sp.expand(remainder + sp.Rational(2, 3) * h) != 0:
        raise AssertionError(f"unexpected moment-3 remainder: {remainder}")

    assert_zero(
        (32 * q4 + 34) * g - (288 * q4**2 + 354 * q4 + 275) * h - 525,
        "7-variable Bezout identity",
    )
    print("[ok] 7-variable complete square-block cubic lift dies at moment 3")


def verify_eight_variable_split_pair() -> None:
    coeff_ring, x, r, q, h, z, ell, w = ring("x,r,q,h,z,ell,w", QQ)
    params_ring = (x, r, q, h, z, ell, w)
    nstate = 8
    zero_exp = (0,) * nstate
    f: dict[tuple[int, ...], object] = {}

    def add(exp: tuple[int, ...], coefficient: object) -> None:
        f[exp] = f.get(exp, coeff_ring.zero) + coefficient

    def ex(**powers: int) -> tuple[int, ...]:
        indexes = {"a": 0, "b": 1, "c": 2, "d": 3, "t": 4, "s": 5, "u": 6, "v": 7}
        out = [0] * nstate
        for name, power in powers.items():
            out[indexes[name]] = power
        return tuple(out)

    one = coeff_ring.one
    for powers in (
        {"a": 1, "d": 1, "t": 1},
        {"a": 1, "c": 1, "d": 1},
        {"b": 1, "t": 2},
        {"b": 1, "c": 1, "t": 1},
    ):
        add(ex(**powers), one)
    add(ex(t=1, u=1, v=1), -2 * x)
    add(ex(s=1, t=2), x)
    add(ex(d=1, t=1, u=1), r)
    add(ex(c=1, u=1, v=1), -q - one)
    add(ex(c=1, s=1, t=1), q)
    add(ex(c=1, d=1, u=1), h)
    add(ex(c=2, s=1), z)
    add(ex(a=1, t=1, v=1), ell)
    add(ex(a=1, c=1, v=1), w)

    operator_terms = ((0, 3, 1), (1, 2, -1), (4, 5, 1), (6, 7, 1))

    def multiply(left: dict, right: dict) -> dict:
        out: dict[tuple[int, ...], object] = {}
        for le, lc in left.items():
            for re, rc in right.items():
                exponent = tuple(le[i] + re[i] for i in range(nstate))
                out[exponent] = out.get(exponent, coeff_ring.zero) + lc * rc
        return {e: co for e, co in out.items() if co}

    def power(base: dict, exponent: int) -> dict:
        out = {zero_exp: one}
        while exponent:
            if exponent & 1:
                out = multiply(out, base)
            exponent //= 2
            if exponent:
                base = multiply(base, base)
        return out

    def apply_operator(poly: dict) -> dict:
        out: dict[tuple[int, ...], object] = {}
        for exponent, coefficient in poly.items():
            for i, j, sign in operator_terms:
                if exponent[i] == 0 or exponent[j] == 0:
                    continue
                new_exp = list(exponent)
                factor = exponent[i] * exponent[j]
                new_exp[i] -= 1
                new_exp[j] -= 1
                key = tuple(new_exp)
                out[key] = out.get(key, coeff_ring.zero) + sign * factor * coefficient
        return {e: co for e, co in out.items() if co}

    moments: dict[int, dict[tuple[int, int], sp.Expr]] = {}
    for m in (2, 3, 4):
        poly = power(f, m)
        for _ in range(m):
            poly = apply_operator(poly)
        extracted: dict[tuple[int, int], sp.Expr] = {}
        for exponent, coefficient in poly.items():
            if any(exponent[i] for i in range(nstate) if i not in (2, 4)):
                raise AssertionError("split-pair moment retained an unexpected variable")
            extracted[(exponent[2], exponent[4])] = coefficient.as_expr()
        moments[m] = extracted

    xs, rs, qs, hs, zs, ls, ws = [entry.as_expr() for entry in params_ring]
    E1 = qs**2 + qs + hs * ws - 2 * zs + 1
    E2 = 4 * xs * qs + 2 * xs + rs * ws - qs + hs * ls - 4 * zs + 1
    E3 = 4 * xs**2 + rs * ls - 2 * qs + 1
    E4 = 4 * xs * qs * zs + 2 * xs * zs + rs * zs * ws - qs**2 - qs * hs * ws - 3 * qs * zs - qs + hs * zs * ls - 4 * zs**2 + zs
    E5 = 8 * xs**2 * zs + 4 * xs * qs**2 - 2 * xs * qs - 2 * xs * hs * ws - 4 * xs * zs - 2 * xs + 2 * rs * zs * ls - 4 * qs**2 - 12 * qs * zs - qs - hs * ws + 8 * zs - 1
    E6 = 12 * xs**2 * qs - 2 * xs**2 - xs * rs * ws - 7 * xs * qs - xs * hs * ls - 12 * xs * zs - xs + rs * qs * ls - 8 * qs**2 + 2 * qs - 2 * hs * ws + 12 * zs - 2
    E7 = 8 * xs**3 - 12 * xs * qs + rs * ls - rs * ws + 3 * qs - hs * ls + 4 * zs - 1

    assert_zero(moments[2][(2, 0)] - 4 * E1, "8-variable m2 c2")
    assert_zero(moments[2][(1, 1)] - 4 * E2, "8-variable m2 ct")
    assert_zero(moments[2][(0, 2)] - 4 * E3, "8-variable m2 t2")
    assert_zero(moments[3][(3, 0)] - 36 * E4, "8-variable m3 c3")
    assert_zero(moments[3][(2, 1)] - 36 * E5, "8-variable m3 c2t")
    assert_zero(moments[3][(1, 2)] - 36 * E6, "8-variable m3 ct2")
    assert_zero(moments[3][(0, 3)] - 36 * E7, "8-variable m3 t3")

    T = 36 * xs**4 + 7 * xs**2 * rs * ls - 58 * xs**2 * qs + 7 * xs**2 + xs * rs * ls - xs * rs * ws + 19 * xs * qs - xs * hs * ls + 20 * xs * zs - xs + rs**2 * ls**2 - 3 * rs * qs * ls + 2 * rs * ls + 11 * qs**2 - 7 * qs + hs * ws - 10 * zs + 2
    assert_zero(moments[4][(0, 4)] - 576 * T, "8-variable fourth-moment t4 coefficient")

    multipliers = [
        8 * xs**2 + 4 * xs - 12 * qs - 1,
        -2 * xs**2 - 2 * xs * qs + 3 * qs + 8 * zs + 1,
        3 * xs**2 + 2 * xs * qs - 8 * xs * zs - 2 * xs + rs * ls + 2 * qs**2 - 2 * qs - 4 * zs,
        -8,
        2 * (2 * xs + 1),
        -2 * (xs + qs + 1),
        3 * xs + 3 * qs + 1,
    ]
    assert_zero(
        T - 1 - sum(mult * eq for mult, eq in zip(multipliers, [E1, E2, E3, E4, E5, E6, E7])),
        "8-variable unit certificate",
    )
    print("[ok] 8-variable complete split-pair cubic lift dies at moment 4")


def verify_eight_variable_coarse_square_block() -> None:
    coeff_ring, A11, A12, A22, B11, B12, B22, gamma = ring(
        "A11,A12,A22,B11,B12,B22,gamma", QQ
    )
    params_ring = (A11, A12, A22, B11, B12, B22, gamma)
    nstate = 8
    zero_exp = (0,) * nstate
    f: dict[tuple[int, ...], object] = {}

    def add(exp: tuple[int, ...], coefficient: object) -> None:
        f[exp] = f.get(exp, coeff_ring.zero) + coefficient

    def ex(**powers: int) -> tuple[int, ...]:
        indexes = {"a": 0, "b": 1, "c": 2, "d": 3, "t": 4, "s": 5, "z1": 6, "z2": 7}
        out = [0] * nstate
        for name, power in powers.items():
            out[indexes[name]] = power
        return tuple(out)

    one = coeff_ring.one
    for powers in (
        {"a": 1, "d": 1, "t": 1},
        {"a": 1, "c": 1, "d": 1},
        {"b": 1, "t": 2},
        {"b": 1, "c": 1, "t": 1},
    ):
        add(ex(**powers), one)
    add(ex(s=1, t=2), -(A11 + A22))
    add(ex(s=1, c=1, t=1), -(one + 2 * (B11 + B22)))
    add(ex(s=1, c=2), gamma)
    add(ex(t=1, z1=2), A11)
    add(ex(t=1, z1=1, z2=1), 2 * A12)
    add(ex(t=1, z2=2), A22)
    add(ex(c=1, z1=2), B11)
    add(ex(c=1, z1=1, z2=1), 2 * B12)
    add(ex(c=1, z2=2), B22)

    operator_terms = ((0, 3, 1), (1, 2, -1), (4, 5, 1), (6, 6, 1), (7, 7, 1))

    def multiply(left: dict, right: dict) -> dict:
        out: dict[tuple[int, ...], object] = {}
        for le, lc in left.items():
            for re, rc in right.items():
                exponent = tuple(le[i] + re[i] for i in range(nstate))
                out[exponent] = out.get(exponent, coeff_ring.zero) + lc * rc
        return {e: co for e, co in out.items() if co}

    def power(base: dict, exponent: int) -> dict:
        out = {zero_exp: one}
        while exponent:
            if exponent & 1:
                out = multiply(out, base)
            exponent //= 2
            if exponent:
                base = multiply(base, base)
        return out

    def apply_operator(poly: dict) -> dict:
        out: dict[tuple[int, ...], object] = {}
        for exponent, coefficient in poly.items():
            for i, j, sign in operator_terms:
                if i == j:
                    if exponent[i] < 2:
                        continue
                    new_exp = list(exponent)
                    factor = exponent[i] * (exponent[i] - 1)
                    new_exp[i] -= 2
                else:
                    if exponent[i] == 0 or exponent[j] == 0:
                        continue
                    new_exp = list(exponent)
                    factor = exponent[i] * exponent[j]
                    new_exp[i] -= 1
                    new_exp[j] -= 1
                key = tuple(new_exp)
                out[key] = out.get(key, coeff_ring.zero) + sign * factor * coefficient
        return {e: co for e, co in out.items() if co}

    moments: dict[int, list[sp.Expr]] = {}
    for m in (2, 3, 4):
        poly = power(f, m)
        for _ in range(m):
            poly = apply_operator(poly)
        terms = []
        for exponent, coefficient in sorted(poly.items(), key=lambda item: (-item[0][2], -item[0][4])):
            if any(exponent[i] for i in range(nstate) if i not in (2, 4)):
                raise AssertionError("coarse-block moment retained an unexpected variable")
            terms.append(coefficient.as_expr())
        moments[m] = terms

    params = [entry.as_expr() for entry in params_ring]

    def primitive(expr: sp.Expr) -> sp.Expr:
        return sp.Poly(expr, *params, domain=QQ).primitive()[1].as_expr()

    e23 = [primitive(expr) for m in (2, 3) for expr in moments[m]]
    gb = sp.groebner(e23, *params, order="grevlex")
    g = 8 * params[-1]**3 + 12 * params[-1]**2 + 10 * params[-1] + 5
    if gb.reduce(g)[1] != 0:
        raise AssertionError("coarse 8-variable elimination cubic missing")

    normalized_m4 = [primitive(expr) for expr in moments[4]]
    remainders = [sp.factor(gb.reduce(expr)[1]) for expr in normalized_m4]
    h = 4 * params[-1]**2 + 10 * params[-1] + 15
    expected = [h, h, 3 * h, h, h]
    for index, (actual, wanted) in enumerate(zip(remainders, expected)):
        if sp.expand(actual - wanted) != 0 and sp.expand(actual + wanted) != 0:
            raise AssertionError(f"coarse m4 remainder {index}: {actual}")

    assert_zero(g + (2 - 2 * params[-1]) * h - 35, "coarse 8-variable Bezout identity")
    print("[ok] 8-variable complete coarse two-square cubic lift dies at moment 4")


def main() -> None:
    verify_meng_yang_trace_obstruction()
    verify_seven_variable_square_block()
    verify_eight_variable_split_pair()
    verify_eight_variable_coarse_square_block()
    print("All ordinary-Delta frontier certificates passed.")


if __name__ == "__main__":
    main()
