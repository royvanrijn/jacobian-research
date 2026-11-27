#!/usr/bin/env python3
"""Replay the support bounds and centered bracket identities for the triple line."""

from sympy import Matrix, Poly, Symbol, expand, eye, symbols


def check_support():
    x, y, v = symbols("x y v")
    terms = [
        x**i*y**j*(v - x*y)**k
        for i in range(5)
        for j in range(5 - i)
        for k in range(5 - i - j)
    ]
    terms += [x**i*y**(5 - i) for i in range(6)]
    y_polys = [Poly(term, y) for term in terms]
    allowed = {
        1: (4, 2, 0),
        2: (4, 2),
        3: (4, 0),
        4: (4,),
        5: (0,),
    }
    for m, bounds in allowed.items():
        high = [
            Poly(sum(poly.nth(j)*y**j for j in range(m + 1, 6)), x, y, v)
            for poly in y_polys
        ]
        keys = sorted({key for poly in high for key in poly.monoms() if key != (0, 0, 0)})
        matrix = Matrix([[poly.coeff_monomial(key) for poly in high] for key in keys])
        kernel = matrix.nullspace() if keys else eye(len(terms)).columnspace()
        for vector in kernel:
            coefficient = Poly(sum(poly.nth(m)*c for poly, c in zip(y_polys, vector)),
                               x, v)
            assert all(
                j < len(bounds) and i <= bounds[j]
                for (i, j), value in coefficient.terms() if value != 0
            )


def check_brackets():
    y = Symbol("y")
    q, k, f, g, dq, dk, df, dg = symbols("q k f g dq dk df dg")
    Q_y = 2*q*y
    DQ = dq*y**2 + dk
    P_y = f + 3*g*y**2
    DP = df*y + dg*y**3
    quadratic = Poly(expand(DP*Q_y - P_y*DQ), y)
    assert quadratic.nth(0) == -f*dk
    assert quadratic.nth(2) == expand(2*q*df - f*dq - 3*g*dk)

    s, c, d, alpha, beta, gamma, delta = symbols("s c d alpha beta gamma delta")
    dc, dd, da, dbt, dgm, ddt = symbols("dc dd da dbt dgm ddt")
    b, e = symbols("b e")
    Q_s = 3*s**2 + c
    DQ = dc*s + dd
    P_s = 4*s**3 + 3*alpha*s**2 + 2*beta*s + gamma
    DP = da*s**3 + dbt*s**2 + dgm*s + ddt
    cubic = Poly(expand(DP*Q_s - P_s*DQ), s)
    substitute = {
        da: 0,
        beta: 4*c/3 + b,
        dbt: 4*dc/3,
        gamma: 4*d/3 + alpha*c + e,
        dgm: 4*dd/3 + alpha*dc,
        ddt: alpha*dd + 4*c*dc/9 + 2*b*dc/3,
    }
    for power in (2, 3, 4, 5):
        assert expand(cubic.nth(power).subs(substitute, simultaneous=True)) == 0
    C, E = c + 3*b/2, d + 3*e/4
    assert expand(cubic.nth(1).subs(substitute, simultaneous=True)
                  + 4*(C*dd + E*dc)/3) == 0
    assert expand(cubic.nth(0).subs(substitute, simultaneous=True)
                  - 4*c*C*dc/9 + 4*E*dd/3) == 0


if __name__ == "__main__":
    check_support()
    check_brackets()
    print("Binary triple-line scalar: exact support and bracket identities verified.")
