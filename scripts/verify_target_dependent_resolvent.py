"""Exact regressions for arbitrary polynomial resolvent derivatives."""

from __future__ import annotations

import sympy as sp


x, y, z, T, Pvar, Qvar, Avar = sp.symbols("x y z T P Q A")


def direct_jacobian_formula(a: int, b: int, h_polynomial: sp.Expr) -> None:
    """Check det(P,Q,R)=-A^n H(s,P,Q) on a nonlinear seed."""

    f = y**2 + y + 1
    correction = y * Avar + Avar**2
    A = 1 + x * f
    B = A**b * z + correction.subs(Avar, A)
    P = A**a * B
    Q = y + x * A ** (a - 1) * B
    s = x / A
    jac_spq = sp.det(sp.Matrix([s, P, Q]).jacobian([x, y, z]))
    direct = sp.cancel(
        jac_spq * h_polynomial.subs({T: s, Pvar: P, Qvar: Q})
    )
    predicted = sp.cancel(
        -A ** (a + b - 2) * h_polynomial.subs({T: s, Pvar: P, Qvar: Q})
    )
    assert sp.cancel(direct - predicted) == 0


def bounded_generic_classification(n: int) -> tuple[int, int]:
    """Solve a generic bounded H and recover the unique D^n line."""

    # Use f(Y)=Y.  Then D=1-T(Q-PT) has total degree three, so D^n
    # lies in the full space of polynomials of total degree at most 3n.
    degree = 3 * n
    monomials = [
        Pvar**i * Qvar**j * T**k
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
        for k in range(degree + 1 - i - j)
    ]
    coefficients = sp.symbols(f"c0:{len(monomials)}")
    lam = sp.symbols("lam")
    generic_h = sum(c * monomial for c, monomial in zip(coefficients, monomials))
    d_power = sp.expand((1 - T * (Qvar - Pvar * T)) ** n)
    difference = sp.Poly(sp.expand(generic_h - lam * d_power), Pvar, Qvar, T)
    solved = sp.linsolve(difference.coeffs(), (*coefficients, lam))
    solution = next(iter(solved))

    # Exactly one free parameter remains, namely lambda, and every
    # coefficient is its corresponding coefficient in D^n.
    free = set().union(*(entry.free_symbols for entry in solution))
    assert free == {lam}
    for entry, monomial in zip(solution[:-1], monomials):
        expected = lam * sp.Poly(d_power, Pvar, Qvar, T).coeff_monomial(monomial)
        assert sp.expand(entry - expected) == 0
    assert solution[-1] == lam
    return len(monomials), len(difference.coeffs())


for a, b, h_polynomial in [
    (1, 2, 1 + T + Pvar * T**2),
    (2, 3, Qvar**2 - Pvar * T + T**4),
    (1, 4, (1 - T * (Qvar - Pvar * T)) ** 3),
]:
    direct_jacobian_formula(a, b, h_polynomial)

sizes = [bounded_generic_classification(n) for n in (1, 2)]
print(f"PASS: target-dependent derivative rigidity in bounded spaces {sizes}")
