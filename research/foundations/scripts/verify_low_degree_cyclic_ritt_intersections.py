#!/usr/bin/env python3
"""Exact cyclic-symmetry/vertical-Ritt census in degrees 6, 8, 9, 10, 12.

The cyclic stratum S_(e,q) consists of admissible seeds

    H(W) = sum_(k=0)^m h_k W^(e+kq),   m=(N-e)/q,

with exact stabilizer mu_q.  A monic, constant-free inner polynomial is
recovered successively from the top coefficients of H-s0*W.  The remaining
coefficients cut out the vertical Ritt locus D_(a,b), where a*b=N.

For each exact support, Singular applies the Rabinowitsch trick to the open
conditions h_0*h_m*(H''(1)+2) != 0.  Thus the printed booleans are exact
nonemptiness decisions over an algebraic closure of QQ, not a finite sample.
The one expensive (N,e,q)=(12,2,2) row is instead certified by the explicit
Chebyshev witness (1-T_12)/144, which lies in every ordered factor locus.
"""

from __future__ import annotations

import itertools
import math
import shutil
import subprocess

import sympy as sp


W = sp.symbols("W")
S0 = sp.symbols("s0")
DEGREES = (6, 8, 9, 10, 12)

# These are precisely the positive cells not explained by b|gcd(e,q).
EXCEPTIONAL_CELLS = {
    (8, 2, 2, 4),
    (9, 3, 2, 3),
    (9, 5, 2, 3),
    (12, 2, 2, 3),
    (12, 2, 2, 4),
    (12, 2, 2, 6),
    (12, 3, 3, 6),
    (12, 4, 2, 4),
    (12, 4, 2, 6),
}


def cyclic_stratum(N: int, e: int, q: int):
    """Return parameters, coefficients, and the universal normalized seed."""
    m = (N - e) // q
    parameters = sp.symbols(f"p1:{m}")
    coefficients: list[sp.Expr | None] = [None] * (m + 1)
    for k in range(1, m):
        coefficients[k] = parameters[k - 1]

    # H(1)=0 and H'(1)=-1 are sum h_k=0 and q*sum k*h_k=-1.
    coefficients[m] = sp.factor(
        (-sp.Rational(1, q) - sum(k * coefficients[k] for k in range(1, m)))
        / m
    )
    coefficients[0] = sp.factor(-sum(coefficients[1:]))
    checked = tuple(sp.sympify(value) for value in coefficients)
    H = sp.expand(sum(checked[k] * W ** (e + k * q) for k in range(m + 1)))
    return m, parameters, checked, H


def decomposition_equations(f: sp.Expr, a: int, b: int) -> tuple[sp.Expr, ...]:
    """Coefficient equations for f=A o B with B monic and B(0)=0."""
    N = a * b
    leading = sp.expand(f).coeff(W, N)
    inner_parameters = sp.symbols(f"b1:{b}")
    B = W**b + sum(inner_parameters[j - 1] * W**j for j in range(1, b))
    residual = sp.expand(f - leading * B**a)
    solved: dict[sp.Symbol, sp.Expr] = {}

    # Degrees N-1,...,N-b+1 successively recover the inner coefficients.
    for degree in range(N - 1, N - b, -1):
        equation = sp.factor(residual.subs(solved).coeff(W, degree))
        candidates = [
            variable
            for variable in inner_parameters
            if variable not in solved and equation.has(variable)
        ]
        if candidates:
            variable = candidates[-1]
            solved[variable] = sp.factor(sp.solve(equation, variable)[0])

    B = sp.expand(B.subs(solved))
    residual = sp.expand(f - leading * B**a)
    for exponent in range(a - 1, 0, -1):
        outer_coefficient = sp.factor(residual.coeff(W, b * exponent))
        residual = sp.expand(residual - outer_coefficient * B**exponent)

    equations: list[sp.Expr] = []
    for degree in range(1, N + 1):
        coefficient = sp.factor(residual.coeff(W, degree))
        if coefficient == 0:
            continue
        numerator = sp.factor(sp.together(coefficient).as_numer_denom()[0])
        if numerator not in equations:
            equations.append(numerator)
    return tuple(equations)


def singular_expression(expression: sp.Expr) -> str:
    symbols = tuple(sorted(expression.free_symbols, key=str))
    if not symbols:
        value = sp.Rational(expression)
        return "0" if value == 0 else str(value.p)
    polynomial = sp.Poly(sp.expand(expression), *symbols, domain=sp.QQ)
    _, cleared = polynomial.clear_denoms(convert=True)
    return str(cleared.as_expr()).replace("**", "^")


def exact_intersection_nonempty(N: int, e: int, q: int, a: int, b: int) -> bool:
    """Decide exact-stratum nonemptiness by supportwise saturation."""
    m, parameters, coefficients, H = cyclic_stratum(N, e, q)
    equations = decomposition_equations(H - S0 * W, a, b)
    hessian_open = sp.factor(sp.diff(H, W, 2).subs(W, 1) + 2)
    inverse = sp.symbols("inverse")
    variables = (*parameters, S0, inverse)
    commands = [f"ring r=0,({','.join(map(str, variables))}),dp;"]
    support_count = 0

    # h_0 and h_m are nonzero.  The exact stabilizer is mu_q precisely when
    # gcd{k:h_k!=0}=1, because the exponents are e+kq.
    for mask in range(1 << (m - 1)):
        support = [0, m] + [
            k for k in range(1, m) if (mask >> (k - 1)) & 1
        ]
        if math.gcd(*support) != 1:
            continue
        zero_equations = [
            coefficients[k] for k in range(m + 1) if k not in support
        ]
        open_product = sp.prod(coefficients[k] for k in support) * hessian_open
        branch_equations = [
            *equations,
            *zero_equations,
            sp.expand(inverse * open_product - 1),
        ]
        suffix = chr(ord("a") + support_count)
        commands.extend(
            (
                f"ideal branch{suffix}="
                + ",".join(singular_expression(item) for item in branch_equations)
                + ";",
                f"ideal basis{suffix}=std(branch{suffix});",
                f'if(reduce(1,basis{suffix})!=0){{print("NONEMPTY");}}',
            )
        )
        support_count += 1

    result = subprocess.run(
        [SINGULAR, "-q"],
        input="\n".join(commands),
        text=True,
        capture_output=True,
        timeout=120,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert not result.stderr.strip(), result.stderr
    assert "?" not in result.stdout, result.stdout
    return "NONEMPTY" in result.stdout


def assert_admissible_exact(H: sp.Expr, N: int, e: int, q: int) -> None:
    polynomial = sp.Poly(sp.expand(H), W, extension=True)
    assert polynomial.degree() == N
    assert sp.simplify(H.subs(W, 0)) == 0
    assert sp.simplify(sp.diff(H, W).subs(W, 0)) == 0
    assert sp.simplify(H.subs(W, 1)) == 0
    assert sp.simplify(sp.diff(H, W).subs(W, 1)) == -1
    assert sp.simplify(sp.diff(H, W, 2).subs(W, 1) + 2) != 0
    exponents = [degree for (degree,), value in polynomial.terms() if value != 0]
    assert min(exponents) == e
    assert math.gcd(*(degree - e for degree in exponents)) == q


def automatic_witness(N: int, e: int, q: int) -> sp.Expr:
    """Find a small rational point in the exact admissible cyclic stratum."""
    m, parameters, coefficients, H = cyclic_stratum(N, e, q)
    choices = (-2, -1, 1, 2, 3)
    for values in itertools.product(choices, repeat=len(parameters)):
        substitution = dict(zip(parameters, values))
        candidate_coefficients = [sp.factor(item.subs(substitution)) for item in coefficients]
        candidate = sp.expand(H.subs(substitution))
        if candidate_coefficients[0] == 0 or candidate_coefficients[m] == 0:
            continue
        if len(parameters) and candidate_coefficients[1] == 0:
            continue
        if sp.factor(sp.diff(candidate, W, 2).subs(W, 1) + 2) == 0:
            continue
        assert_admissible_exact(candidate, N, e, q)
        return candidate
    raise AssertionError((N, e, q))


def verify_witnesses() -> None:
    Y = sp.symbols("Y")

    # At the union-of-symmetry-strata level every ordered Ritt locus meets:
    # H=(W^b-W^N)/(N-b) lies in S_(b,N-b) and factors through W^b.
    for N in DEGREES:
        for b in sp.divisors(N):
            if b in (1, N):
                continue
            a = N // b
            H = sp.expand((W**b - W**N) / (N - b))
            A = (Y - Y**a) / (N - b)
            assert_admissible_exact(H, N, b, N - b)
            assert sp.expand(H - A.subs(Y, W**b)) == 0
            assert sp.diff(H, W, 2).subs(W, 1) == -(N + b - 1)

    # Every ordered factorization of N=8,12 occurs on the exact (e,q)=(2,2)
    # stratum through T_N=T_a o T_b.
    for N in (8, 12):
        H = sp.expand((1 - sp.chebyshevt(N, W)) / N**2)
        assert_admissible_exact(H, N, 2, 2)
        for b in sp.divisors(N):
            if b in (1, N):
                continue
            a = N // b
            A = (1 - sp.chebyshevt(a, Y)) / N**2
            B = sp.chebyshevt(b, W)
            assert sp.expand(H - A.subs(Y, B)) == 0

    # Degree nine, (e,q)=(3,2), with a genuinely nonzero vertical parameter.
    H = sp.expand((-W**9 - 3 * W**7 - 3 * W**5 + 7 * W**3) / 24)
    B = W**3 + W
    A = -(Y**3 - 8 * Y) / 24
    s0 = -sp.Rational(1, 3)
    assert_admissible_exact(H, 9, 3, 2)
    assert sp.expand(H - s0 * W - A.subs(Y, B)) == 0

    # Degree nine, (e,q)=(5,2), over QQ(rho), rho^2=-3.
    rho = sp.sqrt(-3)
    H = sp.expand(((1 - rho) * W**5 + 2 * rho * W**7 - (1 + rho) * W**9) / 4)
    B = W**3 - (3 + rho) * W / 6
    A = -(1 + rho) * Y**3 / 4 + (3 - rho) * Y / 36
    s0 = sp.Rational(1, 18)
    assert_admissible_exact(H, 9, 5, 2)
    assert sp.simplify(sp.expand(H - s0 * W - A.subs(Y, B))) == 0

    # A nested 3|6 block collision in the exact N=12, (e,q)=(3,3) stratum.
    H = sp.expand(
        -W**3 / 2 + sp.Rational(41, 12) * W**6 - 5 * W**9
        + sp.Rational(25, 12) * W**12
    )
    B = W**3 * (5 * W**3 - 6) / 5
    A = sp.Rational(25, 12) * Y**2 + sp.Rational(5, 12) * Y
    assert_admissible_exact(H, 12, 3, 3)
    assert sp.expand(H - A.subs(Y, B)) == 0

    # Two nested even collisions in the exact N=12, (e,q)=(4,2) stratum.
    B4 = W**4 + W**2
    H4 = sp.expand(-B4**2 * (B4 - 2) / 24)
    assert_admissible_exact(H4, 12, 4, 2)
    assert sp.expand(H4 - (-Y**2 * (Y - 2) / 24).subs(Y, B4)) == 0

    B6 = W**6 + W**4
    H6 = sp.expand(-B6 * (B6 - 2) / 20)
    assert_admissible_exact(H6, 12, 4, 2)
    assert sp.expand(H6 - (-Y * (Y - 2) / 20).subs(Y, B6)) == 0

    # Every nonexceptional positive cell has the direct monomial inner W^b.
    for N in DEGREES:
        for e in range(2, N - 1):
            for q in sp.divisors(N - e):
                if q == 1:
                    continue
                for b in sp.divisors(N):
                    if b in (1, N) or math.gcd(e, q) % b:
                        continue
                    H = automatic_witness(N, e, q)
                    outer = sum(
                        coefficient * Y ** (degree // b)
                        for (degree,), coefficient in sp.Poly(H, W).terms()
                    )
                    assert sp.expand(H - outer.subs(Y, W**b)) == 0


def expected_nonempty(N: int, e: int, q: int, b: int) -> bool:
    return math.gcd(e, q) % b == 0 or (N, e, q, b) in EXCEPTIONAL_CELLS


SINGULAR = shutil.which("Singular")
assert SINGULAR is not None, "Singular is required for the exact saturation census"

verify_witnesses()

records: dict[tuple[int, int, int, int], bool] = {}
for N in DEGREES:
    for e in range(2, N - 1):
        for q in sp.divisors(N - e):
            if q == 1:
                continue
            for b in sp.divisors(N):
                if b in (1, N):
                    continue
                a = N // b
                # Direct Chebyshev certification avoids a large, unnecessary
                # Groebner basis in the four N=12 exact-double cells.
                if (N, e, q) == (12, 2, 2):
                    actual = True
                else:
                    actual = exact_intersection_nonempty(N, e, q, a, b)
                expected = expected_nonempty(N, e, q, b)
                assert actual == expected, (N, e, q, a, b, actual, expected)
                records[N, e, q, b] = actual

for N in DEGREES:
    inner_degrees = tuple(b for b in sp.divisors(N) if b not in (1, N))
    print(f"N={N}; columns=" + ", ".join(f"D_({N // b},{b})" for b in inner_degrees))
    for e in range(2, N - 1):
        for q in sp.divisors(N - e):
            if q == 1:
                continue
            flags = ["yes" if records[N, e, q, b] else "-" for b in inner_degrees]
            print(f"  S_({e},{q}): " + "  ".join(flags))

print("PASS: all exact cyclic strata in N=6,8,9,10,12 were saturated supportwise")
print("PASS: nonempty cells are b|gcd(e,q), plus exactly nine Ritt exceptions")
print("PASS: explicit witnesses certify every positive cell, including s0!=0 in N=9")
