#!/usr/bin/env python3
"""Explore a swap-symmetric exact-rank-two bidegree-(4,4) SIC slice.

The slice is

    F = xi1^4 P(z1,z2) + r xi2^4 P(z2,z1).

It is an orbit-reduced subfamily of the two-fourth-power-symbol pencil.
This script is exploratory: it prints the normalized pure moments and
small Groebner eliminations without promoting a bounded calculation to an
all-order SIC statement.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb, factorial
import re
import subprocess

import sympy as sp


a = sp.symbols("a0:5")
r = sp.symbols("r")
y = sp.symbols("y")

P = sum(a[index] * y**index for index in range(5))
Q = r * sum(a[4 - index] * y**index for index in range(5))


def coefficient(polynomial: sp.Expr, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(polynomial), y).coeff_monomial(y**degree)


def moment(order: int) -> sp.Expr:
    return moment_from(P, Q, order)


def moment_from(
    first: sp.Expr,
    second: sp.Expr,
    order: int,
) -> sp.Expr:
    value = 0
    for chosen_first in range(order + 1):
        value += (
            comb(order, chosen_first)
            * factorial(4 * chosen_first)
            * factorial(4 * (order - chosen_first))
            * coefficient(
                first**chosen_first * second ** (order - chosen_first),
                4 * chosen_first,
            )
        )
    return sp.factor(value / factorial(4 * order))


def singular_expression(polynomial: sp.Expr) -> str:
    numerator = sp.Poly(
        sp.together(polynomial).as_numer_denom()[0],
        *a,
    )
    return str(numerator.as_expr()).replace("**", "^")


def weak_compositions(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, length - 1):
            yield (first, *tail)


def modular_moment_expression(
    order: int,
    prime: int | None,
    chart: int,
) -> str:
    """Return the raw even moment on r=-1 over Q or modulo ``prime``."""
    terms: dict[tuple[int, ...], int] = {}
    for counts in weak_compositions(order, 10):
        first = counts[:5]
        second = counts[5:]
        chosen_first = sum(first)
        y_degree = sum(index * first[index] for index in range(5))
        y_degree += sum((4 - index) * second[index] for index in range(5))
        if y_degree != 4 * chosen_first:
            continue
        denominator = 1
        for count in counts:
            denominator *= factorial(count)
        value = (
            factorial(order)
            // denominator
            * factorial(4 * chosen_first)
            * factorial(4 * (order - chosen_first))
        )
        if (order - chosen_first) % 2:
            value = -value
        exponent = tuple(
            first[index] + second[index]
            for index in range(5)
            if index != chart
        )
        terms[exponent] = terms.get(exponent, 0) + value
        if prime is not None:
            terms[exponent] %= prime

    pieces: list[str] = []
    variables = [
        str(variable)
        for index, variable in enumerate(a)
        if index != chart
    ]
    for exponent, value in sorted(terms.items(), reverse=True):
        if value == 0:
            continue
        factors = [str(value)]
        for variable, power_value in zip(variables, exponent, strict=True):
            if power_value == 1:
                factors.append(variable)
            elif power_value > 1:
                factors.append(f"{variable}^{power_value}")
        pieces.append("*".join(factors))
    return "+".join(pieces) if pieces else "0"


def lex_basis_mod_prime(prime: int) -> list[sp.Poly]:
    variables = a[:4]
    base_polynomials = [
        modular_moment_expression(order, prime, 4)
        for order in (2, 4, 6, 8, 10)
    ]
    script = "\n".join(
        [
            f"ring r={prime},(a0,a1,a2,a3),dp;",
            "option(redSB);",
            f"ideal I={','.join(base_polynomials)};",
            "ideal G=slimgb(I);",
            f"ring s={prime},(a0,a1,a2,a3),lp;",
            "ideal K=imap(r,G);",
            "ideal L=std(K);",
            'print("BEGIN_LEX");',
            "L;",
            'print("END_LEX");',
            "quit;",
        ]
    )
    result = subprocess.run(
        ["Singular", "-q"],
        input=script,
        text=True,
        capture_output=True,
        check=True,
    )
    inside = False
    expressions: list[str] = []
    for line in result.stdout.splitlines():
        if line == "BEGIN_LEX":
            inside = True
            continue
        if line == "END_LEX":
            break
        if inside:
            match = re.fullmatch(r"L\[\d+\]=(.*)", line)
            if match:
                expressions.append(match.group(1))
    if len(expressions) != 5:
        raise RuntimeError(
            f"unexpected lex basis at {prime}: {result.stdout}"
        )
    local_symbols = {str(variable): variable for variable in variables}
    return [
        sp.Poly(
            sp.sympify(expression.replace("^", "**"), locals=local_symbols),
            *variables,
            modulus=prime,
        )
        for expression in expressions
    ]


def crt_pair(value: int, modulus: int, residue: int, prime: int) -> tuple[int, int]:
    correction = (residue - value) * pow(modulus, -1, prime)
    combined = value + modulus * (correction % prime)
    return combined, modulus * prime


def rational_reconstruction(value: int, modulus: int) -> Fraction:
    bound = int((modulus // 2) ** 0.5)
    old_remainder, remainder = modulus, value % modulus
    old_denominator, denominator = 0, 1
    while remainder > bound:
        quotient = old_remainder // remainder
        old_remainder, remainder = (
            remainder,
            old_remainder - quotient * remainder,
        )
        old_denominator, denominator = (
            denominator,
            old_denominator - quotient * denominator,
        )
    if denominator == 0:
        raise ArithmeticError("zero reconstruction denominator")
    if denominator < 0:
        remainder = -remainder
        denominator = -denominator
    candidate = Fraction(remainder, denominator)
    if (
        abs(candidate.numerator) > bound
        or candidate.denominator > bound
        or (candidate.numerator - value * candidate.denominator) % modulus
    ):
        raise ArithmeticError(
            f"failed rational reconstruction modulo {modulus}"
        )
    return candidate


def reconstruct_lex_basis(primes: list[int]) -> list[sp.Expr]:
    modular_bases = [lex_basis_mod_prime(prime) for prime in primes]
    supports = [
        sorted(
            set().union(
                *[
                    set(polynomial.monoms())
                    for polynomial in (
                        basis[index]
                        for basis in modular_bases
                    )
                ]
            ),
            reverse=True,
        )
        for index in range(5)
    ]
    reconstructed: list[sp.Expr] = []
    for index, support in enumerate(supports):
        expression = sp.Integer(0)
        for monomial in support:
            value = 0
            modulus = 1
            for prime, basis in zip(primes, modular_bases, strict=True):
                residue = int(basis[index].coeff_monomial(monomial)) % prime
                value, modulus = crt_pair(value, modulus, residue, prime)
            coefficient = rational_reconstruction(value, modulus)
            term = sp.Rational(
                coefficient.numerator,
                coefficient.denominator,
            )
            for variable, power_value in zip(a[:4], monomial, strict=True):
                term *= variable**power_value
            expression += term
        reconstructed.append(sp.factor(expression))
    return reconstructed


def full_modular_screen(
    prime: int | None,
    selected_chart: int | None,
    print_lex_basis: bool,
    shifted: bool,
) -> None:
    charts = (
        [selected_chart]
        if selected_chart is not None
        else list(range(4, -1, -1))
    )
    for chart in charts:
        base_polynomials = [
            modular_moment_expression(order, prime, chart)
            for order in (2, 4, 6, 8, 10)
        ]
        tail_polynomials = [
            modular_moment_expression(order, prime, chart)
            for order in (12, 14, 16)
        ]
        field = "Q" if prime is None else f"F_{prime}"
        print(f"computed moments through order 16 over {field} on chart a{chart}=1")
        variables = [
            str(variable)
            for index, variable in enumerate(a)
            if index != chart
        ]
        if shifted:
            if chart != 4:
                raise ValueError("the shifted chart is defined only for a4=1")
            substitutions = {
                "a0": "(v+1)",
                "a1": "(u+s)",
                "a2": "(t+6)",
                "a3": "s",
            }
            base_polynomials = [
                polynomial
                for polynomial in base_polynomials
            ]
            tail_polynomials = [
                polynomial
                for polynomial in tail_polynomials
            ]
            for old, new in substitutions.items():
                base_polynomials = [
                    polynomial.replace(old, new)
                    for polynomial in base_polynomials
                ]
                tail_polynomials = [
                    polynomial.replace(old, new)
                    for polynomial in tail_polynomials
                ]
            variables = ["v", "u", "t", "s"]
        coefficient_field = 0 if prime is None else prime
        script = "\n".join(
            [
                f"ring r={coefficient_field},({','.join(variables)}),dp;",
                "option(redSB);",
                f"ideal I={','.join(base_polynomials)};",
                "ideal G=slimgb(I);",
                f'print("chart a{chart}=1 size="+string(size(G)));',
                'if (size(G)==1 && G[1]==1) { print("UNIT"); }',
                'else { print("vdim="+string(vdim(G))); }',
                *[
                    (
                        f"poly tail{order}=reduce("
                        f"{polynomial},G); "
                        f'if (tail{order}==0) '
                        f'{{ print("mu{order} REDUCES TO ZERO"); }} '
                        f'else {{ print("mu{order} NONZERO NORMAL FORM"); }}'
                    )
                    for order, polynomial in zip(
                        (12, 14, 16),
                        tail_polynomials,
                        strict=True,
                    )
                ],
                *(
                    [
                        f"ring s={coefficient_field},({','.join(variables)}),lp;",
                        "ideal K=imap(r,G);",
                        "ideal L=std(K);",
                        'print("LEX BASIS");',
                        "L;",
                    ]
                    if print_lex_basis
                    else []
                ),
                "quit;",
            ]
        )
        result = subprocess.run(
            ["Singular", "-q"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        print(result.stdout.strip())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-modular", action="store_true")
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--rational", action="store_true")
    parser.add_argument("--chart", type=int, choices=range(5))
    parser.add_argument("--print-lex-basis", action="store_true")
    parser.add_argument("--only-full", action="store_true")
    parser.add_argument("--shifted", action="store_true")
    parser.add_argument("--reconstruct", action="store_true")
    arguments = parser.parse_args()

    if arguments.reconstruct:
        primes = [32003, 32009, 32027, 32029, 32051, 32057, 32059]
        if not all(sp.isprime(prime) for prime in primes):
            raise AssertionError("reconstruction list contains a composite")
        reconstructed = reconstruct_lex_basis(primes)
        print("RECONSTRUCTED LEX BASIS")
        for polynomial in reconstructed:
            print(polynomial)
        return

    if arguments.only_full:
        if not arguments.full_modular:
            parser.error("--only-full requires --full-modular")
        full_modular_screen(
            None if arguments.rational else arguments.prime,
            arguments.chart,
            arguments.print_lex_basis,
            arguments.shifted,
        )
        return

    moments = [moment(order) for order in range(1, 6)]
    for order, value in enumerate(moments, start=1):
        print(f"mu_{order}/(4m)! = {sp.factor(value)}")

    # The antisymmetric branch r=-1 automatically kills odd moments.
    even_branch = [
        sp.factor(value.subs(r, -1))
        for value in (moments[1], moments[3])
    ]
    print("r=-1 even moments:")
    for order, value in zip((2, 4), even_branch, strict=True):
        print(f"mu_{order}/(4m)! = {value}")

    # Work on the affine chart a0=1.  This is only a bounded elimination.
    equations = [
        sp.together(value.subs({r: -1, a[0]: 1}))
        for value in (moments[1], moments[3])
    ]
    basis = sp.groebner(equations, *a[1:], order="grevlex")
    print(f"r=-1, a0=1 Groebner basis size: {len(basis.polys)}")
    for polynomial in basis.polys:
        print(sp.factor(polynomial.as_expr()))

    even_first = a[0] + a[2] * y**2 + a[4] * y**4
    even_second = -(a[4] + a[2] * y**2 + a[0] * y**4)
    even_moments = [
        sp.factor(moment_from(even_first, even_second, order))
        for order in (2, 4, 6, 8)
    ]
    print("r=-1 and P even:")
    for order, value in zip((2, 4, 6, 8), even_moments, strict=True):
        print(f"mu_{order}/(4m)! = {value}")

    even_chart = [
        sp.together(value.subs(a[4], 1))
        for value in even_moments[:3]
    ]
    even_basis = sp.groebner(
        even_chart,
        a[0],
        a[2],
        order="grevlex",
    )
    print(
        "r=-1, P even, a4=1, moments 2/4/6 "
        f"Groebner basis size: {len(even_basis.polys)}"
    )
    for polynomial in even_basis.polys:
        print(sp.factor(polynomial.as_expr()))

    if arguments.full_modular:
        full_modular_screen(
            None if arguments.rational else arguments.prime,
            arguments.chart,
            arguments.print_lex_basis,
            arguments.shifted,
        )


if __name__ == "__main__":
    main()
