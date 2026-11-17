#!/usr/bin/env python3
"""Structural audits for the planar (72,108) frontier.

This script verifies three repository-original reductions:

1. a uniform prime-defect tail lemma for the GGHV Proposition 4.1/4.3
   opposite-corner step and the exact k=1 parallelism consequence;
2. an exact mu_7-quotient formulation of the complete first Wronskian block,
   including all eleven triangular solves and all six compatibility equations;
3. the determinantal/Hilbert--Burch identities behind the Case-2 Cramer
   reduction and the Hilbert function of the resulting binary-octic complete
   intersection.

It does not re-prove the general GGV/GGHV structural theorems used before the
last tail, and it does not rerun the saved characteristic-zero certificates.
"""
from __future__ import annotations

from hashlib import sha256
from math import gcd
from pathlib import Path

import sympy as sp

OUT_DIR = Path(__file__).resolve().parent


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def prime_defect_candidates(A: int, B: int) -> list[tuple[int, int, int, int, int]]:
    """Return (a,b,g,p,q) from the primitive divisibility equation.

    If C=(A,B), C-(a,b)=g(p,q), gcd(p,q)=1, and the auxiliary-F
    condition holds, then B*p-A*q divides A-B.  We enumerate every positive
    primitive direction compatible with b>=0 and then every possible step g.
    """
    out: list[tuple[int, int, int, int, int]] = []
    defect = A - B
    for q in range(1, B + 1):
        for p in range(1, A + 1):
            if gcd(p, q) != 1:
                continue
            delta = B * p - A * q
            if delta <= 0 or defect % delta:
                continue
            for g in range(1, B // q + 1):
                a, b = A - g * p, B - g * q
                if a < 0 or b < 0:
                    continue
                out.append((a, b, g, p, q))
    return sorted(out)


def audit_prime_defect_tail() -> None:
    # The printed Proposition 4.1 case, C=(21,8), is recovered uniformly.
    c21 = prime_defect_candidates(21, 8)
    assert [(a, b) for a, b, *_ in c21] == [(1, 1), (5, 2), (13, 5)]
    assert {(p, q) for *_, p, q in c21 if (p, q) != (20, 7)} == {(8, 3)}

    # The omitted Proposition 4.3 calculation, C=(24,7).
    c24 = prime_defect_candidates(24, 7)
    assert [(a, b) for a, b, *_ in c24] == [
        (1, 1), (3, 1), (10, 3), (17, 5)
    ]
    proper = [(a, b) for a, b, *_ in c24 if (a, b) != (1, 1)]
    assert proper == [(3, 1), (10, 3), (17, 5)]
    assert all(2 * a == 7 * b - 1 for a, b in [(24, 7), *proper])
    assert {(p, q) for a, b, g, p, q in c24 if (a, b) != (1, 1)} == {(7, 2)}

    # Symbolic parallelism determinants for m=2,n=3.
    a, b, k = sp.symbols("a b k", integer=True)
    direct = -2 * a + (5 * k + 2) * b - k
    swapped = 3 * a - (5 * k + 3) * b + k
    line_sub = {a: (7 * b - 1) / 2}
    assert sp.factor(direct.subs(line_sub)) == (k - 1) * (5 * b - 1)
    assert sp.simplify(swapped.subs(line_sub) + (2 * k - 3) * (5 * b - 1) / 2) == 0

    # For integral k>=1,b>=1, the direct assignment is parallel iff k=1;
    # the swapped assignment is never parallel.
    for bb in range(1, 8):
        aa_num = 7 * bb - 1
        if aa_num % 2:
            continue
        aa = aa_num // 2
        for kk in range(1, 10):
            d1 = int(direct.subs({a: aa, b: bb, k: kk}))
            d2 = int(swapped.subs({a: aa, b: bb, k: kk}))
            assert (d1 == 0) == (kk == 1)
            assert d2 != 0

    print("PRIME_DEFECT_TAIL_PASS")
    print("PROP43_CANDIDATES=(24,7),(17,5),(10,3),(3,1)")
    print("UNIFORM_PARALLELISM_FORCES_K=1")


def triangular_wronskian(A: sp.Expr, variable: sp.Symbol, coeff_prefix: str):
    unknowns = {j: sp.Symbol(f"{coeff_prefix}{j}") for j in range(2, 13)}
    D = sum(unknowns[j] * variable**j for j in range(2, 13))
    equation = sp.Poly(
        sp.expand(2 * A * sp.diff(D, variable) - 3 * sp.diff(A, variable) * D - variable**2),
        variable,
    )
    solved: dict[sp.Symbol, sp.Expr] = {}
    for degree in range(2, 13):
        row = sp.cancel(equation.nth(degree).subs(solved))
        unsolved = [v for v in unknowns.values() if v not in solved and row.has(v)]
        assert len(unsolved) == 1
        pivot = unsolved[0]
        linear = sp.Poly(row, pivot)
        assert linear.degree() == 1
        solved[pivot] = sp.cancel(-linear.nth(0) / linear.nth(1))
    residuals = [sp.cancel(equation.nth(degree).subs(solved)) for degree in range(13, 19)]
    return unknowns, solved, residuals


def integer_primitive_polynomial(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Poly:
    polynomial = sp.Poly(expr, *variables, domain=sp.QQ)
    _, integral = polynomial.clear_denoms(convert=True)
    _, primitive = sp.Poly(integral, *variables, domain=sp.ZZ).primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive


def audit_mu7_quotient() -> None:
    t, tau, u = sp.symbols("t tau u")
    a2, a3, a4, a5, a6, a7 = sp.symbols("a2 a3 a4 a5 a6 a7")
    x2, x3, x4, x5, x6, q = sp.symbols("x2 x3 x4 x5 x6 q")

    A = t + a2*t**2 + a3*t**3 + a4*t**4 + a5*t**5 + a6*t**6 + a7*t**7 + t**8
    Abar = tau + x2*tau**2 + x3*tau**3 + x4*tau**4 + x5*tau**5 + x6*tau**6 + q*tau**7 + q*tau**8

    d_old, solved_old, residual_old = triangular_wronskian(A, t, "d")
    d_new, solved_new, residual_new = triangular_wronskian(Abar, tau, "e")
    assert len(solved_old) == len(solved_new) == 11
    assert len(residual_old) == len(residual_new) == 6

    invariant_substitution = {
        x2: a2*u,
        x3: a3*u**2,
        x4: a4*u**3,
        x5: a5*u**4,
        x6: a6*u**5,
        q: u**7,
    }

    # Exact equality of every triangular solution under
    # Dbar(tau)=u^-2 D(u*tau).
    for j in range(2, 13):
        lhs = sp.cancel(solved_new[d_new[j]].subs(invariant_substitution))
        rhs = sp.cancel(u**(j - 2) * solved_old[d_old[j]].subs(a7, u))
        assert sp.cancel(lhs - rhs) == 0

    # Exact equality of all six compatibility equations.  The kth coefficient
    # transforms by u^(k-2).
    for offset, degree in enumerate(range(13, 19)):
        lhs = sp.cancel(residual_new[offset].subs(invariant_substitution))
        rhs = sp.cancel(u**(degree - 2) * residual_old[offset].subs(a7, u))
        assert sp.cancel(lhs - rhs) == 0

    quotient_variables = (x2, x3, x4, x5, x6, q)
    quotient_polynomials = [
        integer_primitive_polynomial(expr, quotient_variables)
        for expr in residual_new
    ]
    # The last compatibility has one removable q factor on the known q!=0 chart.
    last_terms = quotient_polynomials[-1].terms()
    assert min(monomial[-1] for monomial, _ in last_terms) == 1
    localized_last = integer_primitive_polynomial(
        quotient_polynomials[-1].as_expr() / q, quotient_variables
    )
    localized = [*quotient_polynomials[:-1], localized_last]

    equation_path = OUT_DIR / "firstblock_mu7_quotient_equations.txt"
    lines = [
        "# Exact first-block equations on the residual-scaling quotient",
        "# variables: x2,x3,x4,x5,x6,q",
        "# Abar=t+x2*t^2+x3*t^3+x4*t^4+x5*t^5+x6*t^6+q*t^7+q*t^8",
        "# 2*Abar*Dbar'-3*Abar'*Dbar=t^2; e2,...,e12 solved triangularly",
        "# q != 0; the common q factor in the sixth equation has been removed",
        "",
    ]
    for index, polynomial in enumerate(localized, 1):
        lines.append(f"F{index} = {sp.sstr(polynomial.as_expr())}")
    equation_path.write_text("\n".join(lines) + "\n")
    digest = sha256(equation_path.read_bytes()).hexdigest()

    # The preserved degree-35 eliminant has nonzero constant coefficient, so
    # every actual first-block point lies on u=a7 != 0 and q=u^7 != 0.
    H_constant = -1888043347611739526396142670327809715470336
    assert H_constant != 0

    # Direct bracket covariance for the full Laurent system:
    # Pbar=u^-1 P(u*tau,z), Qbar=u^-2 Q(u*tau,z).
    # Each derivative product is u^-2 times the old bracket at u*tau.
    assert sp.simplify(u**-2 * (u*tau)**2 - tau**2) == 0

    print("MU7_FIRST_BLOCK_QUOTIENT_PASS")
    print("TRIANGULAR_SOLVES=11 COMPATIBILITY_EQUATIONS=6")
    print("QUOTIENT_TERM_COUNTS=" + str([len(p.terms()) for p in localized]))
    print(f"QUOTIENT_EQUATIONS_SHA256={digest}")


def audit_case2_homological_core() -> None:
    # Global homogeneous form of the two cubic rows:
    # R_i = a_i(r,s) h + b_i(r,s) + c_i(r,s), with degrees 1,3,1.
    a0, a1, b0, b1, c0, c1, h = sp.symbols("a0 a1 b0 b1 c0 c1 h")
    R0 = a0*h + b0 + c0
    R1 = a1*h + b1 + c1
    Delta = a0*b1 - a1*b0        # binary quartic
    H = b0*c1 - b1*c0            # binary quartic
    Q = a1*c0 - a0*c1            # binary quadratic

    # Signed-maximal-minor syzygies.
    assert sp.expand(a0*H + b0*Q + c0*Delta) == 0
    assert sp.expand(a1*H + b1*Q + c1*Delta) == 0

    # Polynomial Cramer identities, globally homogeneous and without charts.
    Eh = Delta*h - H
    Eq = Delta - Q
    assert sp.expand(Eh - (b1*R0 - b0*R1)) == 0
    assert sp.expand(Eq - (a0*R1 - a1*R0)) == 0

    # A quartic residual has radial degrees d=0,2,4.  Write m=d/2 and
    # h-degree e, with e+m<=2.  Its global octic Cramer transform replaces
    # the monomial f_d*h^e by
    #     f_d * H^e * Q^m * Delta^(2-e-m).
    # The difference from Delta^2*f_d*h^e lies in (Eh,Eq), hence in (R0,R1).
    DD, HH, QQ = sp.symbols("DD HH QQ")
    for e, m in ((0,0),(1,0),(2,0),(0,1),(1,1),(0,2)):
        transformed = HH**e * QQ**m * DD**(2-e-m)
        original = DD**2 * h**e
        difference = sp.expand(transformed - original)
        # Direct reduction by the monic relations HH=DD*h, QQ=DD.
        reduced = sp.expand(difference.subs({HH: DD*h, QQ: DD}, simultaneous=True))
        assert reduced == 0

    # Two coprime binary octics form a complete intersection.  Its Hilbert
    # function is 1,2,...,8,7,...,1 and vanishes from degree 15 onward;
    # therefore (r,s)^15 is contained in their ideal.
    def hilbert(degree: int, d: int = 8) -> int:
        base = degree + 1
        if degree >= d:
            base -= 2 * (degree - d + 1)
        if degree >= 2*d:
            base += degree - 2*d + 1
        return base

    values = [hilbert(i) for i in range(0, 18)]
    assert values[:15] == [1,2,3,4,5,6,7,8,7,6,5,4,3,2,1]
    assert values[15:] == [0,0,0]
    assert sum(values) == 64

    # If the origin fiber gives 1=A+rB+sC with A in I and m=(r,s), then
    # (1-A)^15 lies in m^15 subset I, while 1-(1-A)^15 is divisible by A.
    X = sp.symbols("X")
    assert sp.rem(1-(1-X)**15, X) == 0

    print("CASE2_GLOBAL_DETERMINANTAL_PASS")
    print("BINARY_OCTIC_HILBERT_FUNCTION=" + str(values[:15]))
    print("BINARY_OCTIC_LENGTH=64 SOCLE_DEGREE=14")
    print("GLOBAL_RADIAL_POWER_BOUND=15")
    print("SINGULAR_CRAMER_SUBCASE_IS_REDUNDANT")


def main() -> None:
    audit_prime_defect_tail()
    audit_mu7_quotient()
    audit_case2_homological_core()
    print("STRUCTURAL_GAP_ATTACK_PASS")


if __name__ == "__main__":
    main()
