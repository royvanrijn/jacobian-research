#!/usr/bin/env python3
"""Generic finite quotients on the reduced s0 common boundary.

An exact sparse recursion constructs the moments over ZZ.  Singular then
computes the principal quotient and the changed fiber bases on L=0, Q=0,
and L=Q=0 over characteristic-zero rational-function fields.  Two modular
replays independently check the same quotient shapes.  This is not by
itself a unit certificate for the full common boundary.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from math import factorial
from pathlib import Path
import re
import shutil
import subprocess

import sympy as sp

from explore_two_pair_sic_bidegree33_full_anchor import (
    PARAMETERS,
    Q_POLYNOMIALS,
    QUADRATIC_Q,
    WEIGHTS,
    chart_expression,
    moment_terms,
    prepare_s0_branch_for_msolve,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_boundary_generic_quotient.json"
)
PRIMES = (47, 101)
J_SPLIT_ROOTS = {47: (6, 41), 101: (48, 53)}
T2_NUMERATOR_ADAPTED = (
    "-4623840*L^3*t0-3107052*L^2*Q*s1+973440*L^2*s1^3"
    "-12654720*L^2*s1*t0^2-973440*L^2*s3+1155960*L^2"
    "-1770552*L*Q^2*t0+93366*L*Q*s1^2*t0-5811272*L*Q*t0^3"
    "-3651800*L*t0^5-1291059*Q^3*s1+430353*Q^2*s1^3"
    "-9089064*Q^2*s1*t0^2-430353*Q^2*s3+474012*Q^2"
    "+1164825*Q*s1^3*t0^2-17449125*Q*s1*t0^4"
    "-1164825*Q*s3*t0^2+1325250*Q*t0^2+768800*s1^3*t0^4"
    "-9994400*s1*t0^6-768800*s3*t0^4+912950*t0^4"
)
J_POINT_MINPOLY = (
    "441554190069069*bb^4+15795130399581456*bb^3"
    "+193851580108553334*bb^2+319468919863825776*bb"
    "+1067521643767708429"
)
J_POINT_ALPHA = (
    "3*(1026265600730531007*bb^3+41799868694363859156*bb^2"
    "+506411570533205547441*bb+545569851002913527492)"
    "/18525795986003750110"
)


def exact_convolve(
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> tuple[int, ...]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            answer[left_index + right_index] += left_coefficient * right_coefficient
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def exact_polynomial_powers(
    polynomial: tuple[int, ...],
    maximum: int,
) -> tuple[tuple[int, ...], ...]:
    powers = [(1,)]
    for _ in range(maximum):
        powers.append(exact_convolve(powers[-1], polynomial))
    return tuple(powers)


def exact_moment_terms(order: int) -> dict[tuple[int, ...], int]:
    """Construct one sparse moment over ZZ by the torus-weight recursion."""

    basis_powers = tuple(
        exact_polynomial_powers(polynomial, order)
        for polynomial in Q_POLYNOMIALS
    )
    quadratic_powers = exact_polynomial_powers(QUADRATIC_Q, order)
    parameter_order = (0, 6, 1, 5, 7, 11, 2, 4, 8, 10, 3, 9)
    exponents = [0] * len(PARAMETERS)
    answer: dict[tuple[int, ...], Fraction] = defaultdict(Fraction)

    @lru_cache(maxsize=None)
    def remaining_weight_bounds(position: int, degree_left: int) -> tuple[int, int]:
        remaining_weights = [WEIGHTS[index] for index in parameter_order[position:]]
        if not remaining_weights or degree_left == 0:
            return 0, 0
        return (
            min(0, degree_left * min(remaining_weights)),
            max(0, degree_left * max(remaining_weights)),
        )

    def visit(
        position: int,
        used_degree: int,
        weight: int,
        shift: int,
        denominator: int,
        q_polynomial: tuple[int, ...],
    ) -> None:
        if position == len(parameter_order):
            if weight != 0:
                return
            quadratic_exponent = order - used_degree
            product = exact_convolve(
                q_polynomial,
                quadratic_powers[quadratic_exponent],
            )
            contraction = sum(
                coefficient
                * factorial(3 * order - shift - q_degree)
                * factorial(shift + q_degree)
                for q_degree, coefficient in enumerate(product)
                if 0 <= shift + q_degree <= 3 * order
            )
            coefficient = Fraction(
                factorial(order) * contraction,
                denominator * factorial(quadratic_exponent),
            )
            if coefficient:
                answer[tuple(exponents)] += coefficient
            return

        parameter_index = parameter_order[position]
        parameter_weight = WEIGHTS[parameter_index]
        available = order - used_degree
        for exponent in range(available + 1):
            new_weight = weight + exponent * parameter_weight
            degree_left = available - exponent
            minimum, maximum = remaining_weight_bounds(position + 1, degree_left)
            if not minimum <= -new_weight <= maximum:
                continue
            exponents[parameter_index] = exponent
            visit(
                position + 1,
                used_degree + exponent,
                new_weight,
                shift + max(parameter_weight, 0) * exponent,
                denominator * factorial(exponent),
                exact_convolve(
                    q_polynomial,
                    basis_powers[parameter_index][exponent],
                ),
            )
        exponents[parameter_index] = 0

    visit(0, 0, 0, 0, 1, (1,))
    assert all(coefficient.denominator == 1 for coefficient in answer.values())
    return {
        exponents: coefficient.numerator
        for exponents, coefficient in answer.items()
        if coefficient
    }


def exact_chart_expression(terms: dict[tuple[int, ...], int]) -> str:
    serialized: list[str] = []
    for exponents, coefficient in sorted(terms.items()):
        factors: list[str] = []
        for variable, exponent in zip(PARAMETERS[1:], exponents[1:]):
            if exponent == 1:
                factors.append(variable)
            elif exponent > 1:
                factors.append(f"{variable}^{exponent}")
        monomial = "*".join(factors)
        if not monomial:
            serialized.append(str(coefficient))
        elif coefficient == 1:
            serialized.append(monomial)
        elif coefficient == -1:
            serialized.append(f"-{monomial}")
        else:
            serialized.append(f"{coefficient}*{monomial}")
    return "+".join(serialized).replace("+-", "-") or "0"


def normalized_factor(
    expression: sp.Expr,
    generators: tuple[sp.Symbol, ...],
) -> sp.Expr:
    return sp.Poly(expression, *generators, domain=sp.QQ).monic().as_expr()


def exact_principal_certificate(
    singular: str,
    polynomials: list[str],
) -> dict[str, object]:
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring r=(0,s1,s2,s3,t0,t1,t2),(s5,t4),dp;
poly p4={polynomials[1]};
poly p5={polynomials[2]};
poly p6={polynomials[3]};
poly p7={polynomials[4]};
ideal G=std(p4,p5);
poly r6=reduce(p6,G);
poly r7=reduce(p7,G);
print("GBSIZE "+string(size(G)));
print("VDIM "+string(vdim(G)));
print("R6SIZE "+string(size(r6)));
print("R7SIZE "+string(size(r7)));
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEADEXP "+string(leadexp(G[basisIndex])));
  print("GLEAD "+string(leadcoef(G[basisIndex])));
}}
poly z;
number c;
poly m;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  z=G[basisIndex];
  while(z!=0)
  {{
    c=leadcoef(z); m=leadmonom(z);
    print("GDEN "+string(denominator(c)));
    z=z-c*m;
  }}
}}
z=r6;
while(z!=0)
{{
  c=leadcoef(z); m=leadmonom(z);
  print("R6DEN "+string(denominator(c)));
  z=z-c*m;
}}
z=r7;
while(z!=0)
{{
  c=leadcoef(z); m=leadmonom(z);
  print("R7DEN "+string(denominator(c)));
  z=z-c*m;
}}
z=r6; c=leadcoef(z);
print("FIRSTTERM "+string(leadmonom(z)));
print("FIRSTNUM "+string(numerator(c)));
print("FIRSTDEN "+string(denominator(c)));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=180,
    )
    summaries = [
        re.search(rf"(?m)^{label} (\d+)$", completed.stdout)
        for label in ("GBSIZE", "VDIM", "R6SIZE", "R7SIZE")
    ]
    assert all(marker is not None for marker in summaries), completed.stdout[:1000]
    assert tuple(int(marker.group(1)) for marker in summaries if marker) == (3, 6, 6, 6)
    leading_exponents = tuple(
        tuple(int(value) for value in marker.split(","))
        for marker in re.findall(r"(?m)^LEADEXP ([0-9,]+)$", completed.stdout)
    )
    assert leading_exponents == ((2, 0), (1, 2), (0, 4))

    s1, s2, s3, t0, t1, t2 = sp.symbols("s1 s2 s3 t0 t1 t2")
    generators = (s1, s2, s3, t0, t1, t2)
    environment = {str(symbol): symbol for symbol in generators}
    linear = s1 * t0 - t1
    quadratic = s1**2 - s2 - sp.Rational(13, 3) * t0**2
    seventh_divisor = (
        9801 * s1**4
        - 19602 * s1**2 * s2
        - 23832 * s1**2 * t0**2
        - 60840 * s1 * t0 * t1
        + 9801 * s2**2
        + 54252 * s2 * t0**2
        + 75076 * t0**4
        + 30420 * t1**2
    )
    named_factors = {
        "L": normalized_factor(linear, generators),
        "Q": normalized_factor(quadratic, generators),
        "J": normalized_factor(seventh_divisor, generators),
    }
    leading_factor_names: list[str] = []
    for value in re.findall(r"(?m)^GLEAD (.*)$", completed.stdout):
        coefficient = sp.sympify(value.replace("^", "**"), locals=environment)
        _, factors = sp.factor_list(coefficient, *generators)
        assert len(factors) == 1 and factors[0][1] == 1
        normalized = normalized_factor(factors[0][0], generators)
        matches = [
            name for name, expected in named_factors.items()
            if normalized == expected
        ]
        assert len(matches) == 1
        leading_factor_names.append(matches[0])
    assert leading_factor_names == ["Q", "L", "J"]

    def denominator_support(tag: str) -> set[str]:
        observed: set[str] = set()
        for value in re.findall(rf"(?m)^{tag} (.*)$", completed.stdout):
            denominator = sp.sympify(value.replace("^", "**"), locals=environment)
            _, factors = sp.factor_list(denominator, *generators)
            for factor, _ in factors:
                normalized = normalized_factor(factor, generators)
                matches = [
                    name for name, expected in named_factors.items()
                    if normalized == expected
                ]
                assert len(matches) == 1
                observed.add(matches[0])
        return observed

    assert denominator_support("GDEN") == set()
    assert denominator_support("R6DEN") == {"L", "Q"}
    assert denominator_support("R7DEN") == {"L", "Q", "J"}

    first_term = re.search(r"(?m)^FIRSTTERM (.*)$", completed.stdout)
    first_numerator = re.search(r"(?m)^FIRSTNUM \((.*)\)$", completed.stdout)
    first_denominator = re.search(r"(?m)^FIRSTDEN \((.*)\)$", completed.stdout)
    assert first_term is not None and first_term.group(1) == "t4^3"
    assert first_numerator is not None and first_denominator is not None
    numerator = sp.Poly(
        sp.sympify(first_numerator.group(1).replace("^", "**"), locals=environment),
        *generators,
        domain=sp.QQ,
    )
    assert len(numerator.terms()) == 42 and numerator.degree(s3) == 1
    derivative = sp.Poly(sp.diff(numerator.as_expr(), s3), *generators, domain=sp.QQ)
    _, primitive_derivative = derivative.primitive()
    pivot = (
        430353 * s1**4
        - 860706 * s1**2 * s2
        - 1591461 * s1**2 * t0**2
        - 1946880 * s1 * t0 * t1
        + 430353 * s2**2
        + 2564901 * s2 * t0**2
        + 3802298 * t0**4
        + 973440 * t1**2
    )
    assert sp.expand(primitive_derivative.as_expr() + pivot) == 0
    t2_derivative = sp.diff(numerator.as_expr(), t2)
    assert numerator.degree(t2) == 1
    assert sp.diff(numerator.as_expr(), s3, t2) == 0
    assert sp.expand(
        t2_derivative + 3903051350016000 * linear * quadratic
    ) == 0
    denominator = sp.sympify(
        first_denominator.group(1).replace("^", "**"),
        locals=environment,
    )
    assert sp.expand(denominator - 39 * linear * quadratic) == 0
    adapted_linear, adapted_quadratic = sp.symbols("L Q")
    adapted_substitution = {
        t1: s1 * t0 - adapted_linear,
        s2: s1**2 - sp.Rational(13, 3) * t0**2 - adapted_quadratic,
    }
    adapted_auxiliary = 99 * adapted_quadratic + 155 * t0**2
    adapted_seventh = (
        adapted_auxiliary**2 + 30420 * adapted_linear**2
    )
    adapted_pivot = (
        32 * adapted_seventh
        + 1179 * adapted_quadratic * adapted_auxiliary
    )
    assert sp.expand(
        seventh_divisor.subs(adapted_substitution) - adapted_seventh
    ) == 0
    assert sp.expand(pivot.subs(adapted_substitution) - adapted_pivot) == 0
    # If J=H=0 and LQ is invertible, the second identity gives W=0;
    # the first then gives L^2=0, a contradiction.
    return {
        "field": "QQ(s1,s2,s3,t0,t1,t2)",
        "groebner_basis_size": 3,
        "leading_exponents_s5_t4": [list(values) for values in leading_exponents],
        "quotient_length": 6,
        "standard_basis": ["1", "t4", "t4^2", "t4^3", "s5", "s5*t4"],
        "groebner_leading_coefficient_factors": ["Q", "L", "J"],
        "mu6_remainder_terms": 6,
        "mu7_remainder_terms": 6,
        "mu6_denominator_factors": ["L", "Q"],
        "mu7_denominator_factors": ["L", "Q", "J"],
        "mu7_additional_divisor": str(seventh_divisor),
        "mu6_t4_cubed_numerator_terms": 42,
        "mu6_t4_cubed_denominator": "39*L*Q",
        "mu6_t4_cubed_t2_derivative": "-100078239744000",
        "s3_pivot": str(pivot),
        "adapted_coordinate_identities": {
            "W": "99*Q+155*t0^2",
            "J": "W^2+30420*L^2",
            "H": "32*J+1179*Q*W",
            "consequence": "(J,H):(L*Q)^infinity=(1)",
        },
    }


def principal_certificate(
    singular: str,
    prime: int,
    polynomials: list[str],
) -> dict[str, object]:
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring r=({prime},s1,s2,s3,t0,t1,t2),(s5,t4),dp;
poly p4={polynomials[1]};
poly p5={polynomials[2]};
poly p6={polynomials[3]};
poly p7={polynomials[4]};
ideal G=std(p4,p5);
poly r6=reduce(p6,G);
poly r7=reduce(p7,G);
print("GBSIZE "+string(size(G)));
print("VDIM "+string(vdim(G)));
print("R6SIZE "+string(size(r6)));
print("R7SIZE "+string(size(r7)));
print("R6 "+string(r6));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    summaries = [
        re.search(rf"(?m)^{label} (\d+)$", completed.stdout)
        for label in ("GBSIZE", "VDIM", "R6SIZE", "R7SIZE")
    ]
    assert all(marker is not None for marker in summaries), completed.stdout[:1000]
    assert tuple(int(marker.group(1)) for marker in summaries if marker) == (3, 6, 6, 6)

    denominators: list[str] = []
    for denominator in re.findall(r"/\(([^()]*)\)\*", completed.stdout):
        if denominator not in denominators:
            denominators.append(denominator)
    assert len(denominators) == 2

    s1, s2, s3, t0, t1, t2 = sp.symbols("s1 s2 s3 t0 t1 t2")
    environment = {
        str(symbol): symbol for symbol in (s1, s2, s3, t0, t1, t2)
    }
    parsed = [
        sp.Poly(
            sp.sympify(value.replace("^", "**"), locals=environment),
            s1,
            s2,
            s3,
            t0,
            t1,
            t2,
            modulus=prime,
        )
        for value in denominators
    ]
    linear = s1 * t0 - t1
    quadratic = s1**2 - s2 - (13 * pow(3, -1, prime) % prime) * t0**2
    expected = [
        sp.Poly(linear * quadratic, *parsed[0].gens, modulus=prime),
        sp.Poly(linear * quadratic**2, *parsed[0].gens, modulus=prime),
    ]
    assert parsed == expected
    return {
        "prime": prime,
        "groebner_basis_size": 3,
        "quotient_length": 6,
        "mu6_remainder_terms": 6,
        "mu7_remainder_terms": 6,
        "mu6_denominator_factors": [
            "s1*t0-t1",
            "s1^2-s2-(13/3)*t0^2",
        ],
        "mu6_denominator_products": ["L*Q", "L*Q^2"],
    }


def substitute(polynomial: str, replacements: tuple[tuple[str, str], ...]) -> str:
    for variable, value in replacements:
        polynomial = re.sub(rf"\b{variable}\b", value, polynomial)
    return polynomial


def exact_base_resultant_certificate(
    singular: str,
    polynomials: list[str],
) -> dict[str, object]:
    """Eliminate t2, then certify the first residual base divisor."""

    replacements = (
        ("t1", "(s1*t0-L)"),
        ("s2", "(s1^2-(13/3)*t0^2-Q)"),
        ("t2", "tt"),
    )
    adapted = [
        substitute(polynomial, replacements)
        for polynomial in polynomials
    ]
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring r=(0,s1,s3,t0,L,Q),(s5,t4),dp;
number tt=({T2_NUMERATOR_ADAPTED})/(93366*L*Q);
poly p3={adapted[0]};
poly p4={adapted[1]};
poly p5={adapted[2]};
poly p6={adapted[3]};
ideal G=std(p4,p5);
poly r6=reduce(p6,G);
poly z=r6;
number c;
poly m;
int coefficientIndex=0;
print("P "+string(numerator(leadcoef(p3))));
while(z!=0 && coefficientIndex<2)
{{
  coefficientIndex++;
  c=leadcoef(z); m=leadmonom(z);
  print("FIBERMONOMIAL "+string(m));
  print("C "+string(numerator(c)));
  z=z-c*m;
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=180,
    )
    fiber_monomials = re.findall(
        r"(?m)^FIBERMONOMIAL (.*)$",
        completed.stdout,
    )
    p_marker = re.search(r"(?m)^P \((.*)\)$", completed.stdout)
    c_markers = re.findall(r"(?m)^C \((.*)\)$", completed.stdout)
    assert fiber_monomials == ["s5*t4", "t4^2"]
    assert p_marker is not None and len(c_markers) == 2
    cubic = p_marker.group(1)
    first_coefficient, second_coefficient = c_markers
    s1, s3, t0, adapted_linear, adapted_quadratic = sp.symbols(
        "s1 s3 t0 L Q"
    )
    adapted_environment = {
        str(symbol): symbol
        for symbol in (s1, s3, t0, adapted_linear, adapted_quadratic)
    }
    cubic_polynomial = sp.Poly(
        sp.sympify(cubic.replace("^", "**"), locals=adapted_environment),
        s1,
        s3,
        t0,
        adapted_linear,
        adapted_quadratic,
        domain=sp.QQ,
    )
    assert len(cubic_polynomial.terms()) == 642
    assert cubic_polynomial.degree(s3) == 3
    assert cubic_polynomial.total_degree() == 21

    resultant_run = subprocess.run(
        [singular, "-q"],
        input=f"""
ring u=(0,s1,t0,L,Q),(s3),dp;
poly P={cubic};
poly C1={first_coefficient};
poly C2={second_coefficient};
number a=leadcoef(P);
poly p2=P-a*s3^3;
number b=leadcoef(p2);
poly p1=p2-b*s3^2;
number c=leadcoef(p1);
number d=leadcoef(p1-c*s3);
number e=leadcoef(C1);
poly q2=C1-e*s3^3;
number f=leadcoef(q2);
poly q1=q2-f*s3^2;
number cg=leadcoef(q1);
number h=leadcoef(q1-cg*s3);
number A=e*b-a*f;
number B=e*c-a*cg;
number C0=e*d-a*h;
number V1=A*(A*c-a*C0)-(A*b-a*B)*B;
number V0=A^2*d-(A*b-a*B)*C0;
poly gcdPC1=gcd(P,C1);
poly R1=resultant(P,C1,s3);
poly R2=resultant(P,C2,s3);
print("DEGREES "+string(deg(P))+" "+string(deg(C1))+" "
      +string(deg(C2))+" "+string(deg(gcdPC1)));
print("V1 "+string(numerator(V1)));
print("V0 "+string(numerator(V0)));
print("R1 "+string(R1));
print("R2 "+string(R2));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=180,
    )
    degree_marker = re.search(
        r"(?m)^DEGREES (\d+) (\d+) (\d+) (\d+)$",
        resultant_run.stdout,
    )
    first_resultant_marker = re.search(
        r"(?m)^R1 \((.*)\)$",
        resultant_run.stdout,
    )
    second_resultant_marker = re.search(
        r"(?m)^R2 \((.*)\)$",
        resultant_run.stdout,
    )
    v1_marker = re.search(r"(?m)^V1 \((.*)\)$", resultant_run.stdout)
    v0_marker = re.search(r"(?m)^V0 \((.*)\)$", resultant_run.stdout)
    assert degree_marker is not None
    assert degree_marker.groups() == ("3", "3", "3", "0")
    assert first_resultant_marker is not None and second_resultant_marker is not None
    assert v1_marker is not None and v0_marker is not None
    first_resultant = first_resultant_marker.group(1)
    second_resultant = second_resultant_marker.group(1)
    linear_pivot_coefficient = v1_marker.group(1)
    linear_pivot_constant = v0_marker.group(1)

    factor_run = subprocess.run(
        [singular, "-q"],
        input=f"""
ring b=0,(s1,t0,L,Q),dp;
poly R1={first_resultant};
poly R2={second_resultant};
poly V1={linear_pivot_coefficient};
poly V0={linear_pivot_constant};
poly D1=L^6*Q^6;
poly D2=L^9*Q^6;
poly S1=R1/D1;
poly S2=R2/D2;
print("EXACT1 "+string(size(R1))+" "+string(deg(R1))+" "
      +string(size(S1))+" "+string(deg(S1))+" "+string(R1-D1*S1==0));
print("EXACT2 "+string(size(R2))+" "+string(deg(R2))+" "
      +string(size(S2))+" "+string(deg(S2))+" "+string(R2-D2*S2==0));
print("VEXACT "+string(size(V1))+" "+string(deg(V1))+" "
      +string(size(V0))+" "+string(deg(V0)));
ring m47=47,(s1,t0,L,Q),dp;
poly S147=imap(b,S1);
poly S247=imap(b,S2);
poly V147=imap(b,V1);
poly V047=imap(b,V0);
poly G47=gcd(S147,S247);
poly G147=gcd(S147,V147);
poly G047=gcd(S147,V047);
list F47=factorize(S147);
list E47=factorize(S247);
print("MOD47 "+string(size(S147))+" "+string(deg(S147))+" "
      +string(size(F47[1]))+" "+string(F47[2][2])+" "
      +string(deg(F47[1][2]))+" "+string(size(S247))+" "
      +string(deg(S247))+" "+string(size(E47[1]))+" "
      +string(E47[2][2])+" "+string(deg(E47[1][2]))+" "
      +string(deg(G47)));
print("VP47 "+string(size(V147))+" "+string(deg(V147))+" "
      +string(size(V047))+" "+string(deg(V047))+" "
      +string(deg(G147))+" "+string(deg(G047)));
ring m101=101,(s1,t0,L,Q),dp;
poly S1101=imap(b,S1);
poly S2101=imap(b,S2);
poly V1101=imap(b,V1);
poly V0101=imap(b,V0);
poly G101=gcd(S1101,S2101);
poly G1101=gcd(S1101,V1101);
poly G0101=gcd(S1101,V0101);
list F101=factorize(S1101);
list E101=factorize(S2101);
print("MOD101 "+string(size(S1101))+" "+string(deg(S1101))+" "
      +string(size(F101[1]))+" "+string(F101[2][2])+" "
      +string(deg(F101[1][2]))+" "+string(size(S2101))+" "
      +string(deg(S2101))+" "+string(size(E101[1]))+" "
      +string(E101[2][2])+" "+string(deg(E101[1][2]))+" "
      +string(deg(G101)));
print("VP101 "+string(size(V1101))+" "+string(deg(V1101))+" "
      +string(size(V0101))+" "+string(deg(V0101))+" "
      +string(deg(G1101))+" "+string(deg(G0101)));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=180,
    )
    exact_markers = [
        re.search(
            rf"(?m)^{label} (\d+) (\d+) (\d+) (\d+) ([01])$",
            factor_run.stdout,
        )
        for label in ("EXACT1", "EXACT2")
    ]
    modular_markers = [
        re.search(
            rf"(?m)^{label} "
            r"(\d+) (\d+) (\d+) (\d+) (\d+) "
            r"(\d+) (\d+) (\d+) (\d+) (\d+) (\d+)$",
            factor_run.stdout,
        )
        for label in ("MOD47", "MOD101")
    ]
    exact_pivot_marker = re.search(
        r"(?m)^VEXACT (\d+) (\d+) (\d+) (\d+)$",
        factor_run.stdout,
    )
    modular_pivot_markers = [
        re.search(
            rf"(?m)^{label} (\d+) (\d+) (\d+) (\d+) (\d+) (\d+)$",
            factor_run.stdout,
        )
        for label in ("VP47", "VP101")
    ]
    assert all(marker is not None for marker in exact_markers)
    assert exact_markers[0].groups() == ("6702", "75", "6702", "63", "1")
    assert exact_markers[1].group(2) == "81"
    assert exact_markers[1].group(4) == "66"
    assert exact_markers[1].group(5) == "1"
    assert all(marker is not None for marker in modular_markers)
    assert modular_markers[0].groups() == (
        "6565", "63", "2", "1", "63",
        "6951", "66", "2", "1", "66", "0",
    )
    assert modular_markers[1].groups() == (
        "6633", "63", "2", "1", "63",
        "7038", "66", "2", "1", "66", "0",
    )
    assert exact_pivot_marker is not None
    assert exact_pivot_marker.group(2) == "60"
    assert exact_pivot_marker.group(4) == "63"
    assert all(marker is not None for marker in modular_pivot_markers)
    assert modular_pivot_markers[0].groups() == (
        "2074", "60", "5083", "63", "0", "0",
    )
    assert modular_pivot_markers[1].groups() == (
        "2090", "60", "5127", "63", "0", "0",
    )

    return {
        "adapted_base_parameters": ["s1", "t0", "L", "Q"],
        "eliminated_variable": "t2",
        "t2_denominator": "93366*L*Q",
        "mu3_polynomial_degree_in_s3": 3,
        "mu3_polynomial_total_degree": 21,
        "mu3_polynomial_terms": 642,
        "first_remaining_mu6_coefficients": ["s5*t4", "t4^2"],
        "first_remaining_mu6_degrees_in_s3": [3, 3],
        "generic_gcd_degree": 0,
        "first_resultant_factorization": "L^6*Q^6*R63",
        "first_resultant_degree": 75,
        "first_residual_factor_degree": 63,
        "first_residual_factor_terms_over_QQ": 6702,
        "residual_factor_modular_terms": {"47": 6565, "101": 6633},
        "residual_factor_irreducible_modulo": [47, 101],
        "second_resultant_factor": "L^9*Q^6*T66",
        "second_resultant_degree": 81,
        "second_residual_factor_degree": 66,
        "second_residual_factor_terms_over_QQ": int(
            exact_markers[1].group(3)
        ),
        "second_residual_factor_modular_terms": {"47": 6951, "101": 7038},
        "second_residual_factor_irreducible_modulo": [47, 101],
        "residual_resultants_gcd": "1",
        "linear_subresultant_identity": "V1*s3+V0=0",
        "linear_subresultant_coefficient_degree": 60,
        "linear_subresultant_constant_degree": 63,
        "linear_subresultant_terms_over_QQ": {
            "V1": int(exact_pivot_marker.group(1)),
            "V0": int(exact_pivot_marker.group(3)),
        },
        "linear_subresultant_modular_terms": {
            "47": {"V1": 2074, "V0": 5083},
            "101": {"V1": 2090, "V0": 5127},
        },
        "linear_subresultant_gcd_with_R63": {
            "V1": "1",
            "V0": "1",
        },
        "linear_subresultant_pivot": "s3=-V0/V1 on V1!=0",
        "first_resultant_sha256": hashlib.sha256(
            first_resultant.encode()
        ).hexdigest(),
        "second_resultant_sha256": hashlib.sha256(
            second_resultant.encode()
        ).hexdigest(),
        "linear_subresultant_coefficient_sha256": hashlib.sha256(
            linear_pivot_coefficient.encode()
        ).hexdigest(),
        "linear_subresultant_constant_sha256": hashlib.sha256(
            linear_pivot_constant.encode()
        ).hexdigest(),
        "scope": (
            "exact QQ resultant with degree-preserving irreducible reductions; "
            "the first two residual resultants are coprime and the first "
            "resultant divisor has a dense linear s3 pivot, but the "
            "codimension-two intersection is not excluded"
        ),
    }


def exact_j_divisor_point_certificate(
    singular: str,
    polynomials: list[str],
) -> dict[str, object]:
    """Exhibit a characteristic-zero mu3=J=0 rank-five point."""

    beta = sp.symbols("beta")
    minimal_polynomial = sp.Poly(
        sp.sympify(
            J_POINT_MINPOLY.replace("bb", "beta").replace("^", "**"),
            locals={"beta": beta},
        ),
        beta,
        domain=sp.QQ,
    )
    assert minimal_polynomial.degree() == 4
    assert minimal_polynomial.is_irreducible
    replacements = (
        ("t1", "0"),
        (
            "s2",
            f"(1-(13/3)-(({J_POINT_ALPHA}-155)/99))",
        ),
        ("s1", "1"),
        ("s3", "bb"),
        ("t0", "1"),
        ("t2", "0"),
    )
    adapted = [
        substitute(polynomial, replacements)
        for polynomial in polynomials
    ]
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring r=(0,bb),(s5,t4),dp;
minpoly={J_POINT_MINPOLY};
number aa={J_POINT_ALPHA};
number Qvalue=(aa-155)/99;
number Hvalue=1179*Qvalue*(99*Qvalue+155);
poly p3={adapted[0]};
poly p4={adapted[1]};
poly p5={adapted[2]};
ideal G=std(p4,p5);
print("FIELD "+string(aa^2+30420==0));
print("JZERO "+string((99*Qvalue+155)^2+30420==0));
print("LQHNONZERO "+string(Qvalue!=0)+" "+string(Hvalue!=0));
print("MU3ZERO "+string(p3==0));
print("QUOTIENT "+string(vdim(G))+" "+string(size(G)));
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEADEXP "+string(leadexp(G[basisIndex])));
  print("LEADNONZERO "+string(leadcoef(G[basisIndex])!=0));
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=180,
    )
    assert re.search(r"(?m)^FIELD 1$", completed.stdout)
    assert re.search(r"(?m)^JZERO 1$", completed.stdout)
    assert re.search(r"(?m)^LQHNONZERO 1 1$", completed.stdout)
    assert re.search(r"(?m)^MU3ZERO 1$", completed.stdout)
    quotient_marker = re.search(
        r"(?m)^QUOTIENT (\d+) (\d+)$",
        completed.stdout,
    )
    assert quotient_marker is not None
    assert quotient_marker.groups() == ("5", "3")
    leading_exponents = re.findall(
        r"(?m)^LEADEXP (\d+),(\d+)$",
        completed.stdout,
    )
    assert leading_exponents == [("2", "0"), ("0", "3"), ("1", "2")]
    assert re.findall(
        r"(?m)^LEADNONZERO ([01])$",
        completed.stdout,
    ) == ["1", "1", "1"]
    return {
        "field_generator": "beta",
        "minimal_polynomial": J_POINT_MINPOLY.replace("bb", "beta"),
        "minimal_polynomial_degree": 4,
        "minimal_polynomial_irreducible_over_QQ": True,
        "sqrt_minus_30420_representation": J_POINT_ALPHA.replace("bb", "beta"),
        "base_point": {
            "s1": "1",
            "s3": "beta",
            "t0": "1",
            "L": "1",
            "t2": "0",
            "Q": "(sqrt(-30420)-155)/99",
        },
        "mu3_zero": True,
        "J_zero": True,
        "L_Q_H_nonzero": True,
        "quotient_length": 5,
        "standard_basis": ["1", "t4", "t4^2", "s5", "s5*t4"],
        "initial_ideal": ["s5^2", "t4^3", "s5*t4^2"],
        "scope": (
            "one exact quartic-number-field special fiber on mu3=J=0; "
            "this proves nonemptiness of the characteristic-zero "
            "rank-five locus, not generic rank five on J=0"
        ),
    }


def exact_j_divisor_generic_certificate(
    polynomials: list[str],
) -> dict[str, object]:
    """Fraction-free Groebner certificate on the generic J divisor."""

    s1, s2, s3, s5, t0, t1, t2, t4, linear, quadratic, alpha = sp.symbols(
        "s1 s2 s3 s5 t0 t1 t2 t4 L Q alpha"
    )
    environment = {
        str(symbol): symbol
        for symbol in (
            s1,
            s2,
            s3,
            s5,
            t0,
            t1,
            t2,
            t4,
            linear,
            quadratic,
            alpha,
        )
    }
    parsed = [
        sp.sympify(
            polynomial.replace("^", "**"),
            locals=environment,
        )
        for polynomial in polynomials
    ]
    split_quadratic = (alpha * linear - 155 * t0**2) / 99
    replacements = {
        t1: s1 * t0 - linear,
        s2: s1**2 - sp.Rational(13, 3) * t0**2 - split_quadratic,
    }
    fiber_polynomials = [
        sp.Poly(
            sp.expand(parsed[index].subs(replacements)),
            s5,
            t4,
        )
        for index in (1, 2)
    ]
    assert [len(polynomial.terms()) for polynomial in fiber_polynomials] == [
        6,
        8,
    ]

    base_ring = sp.QQ.poly_ring(s1, s3, t0, linear, t2)
    zero = base_ring.zero
    one = base_ring.one
    quadratic_constant = 30420

    def coefficient_add(left, right):
        return left[0] + right[0], left[1] + right[1]

    def coefficient_negate(value):
        return -value[0], -value[1]

    def coefficient_multiply(left, right):
        return (
            left[0] * right[0]
            - quadratic_constant * left[1] * right[1],
            left[0] * right[1] + left[1] * right[0],
        )

    def coefficient_is_zero(value) -> bool:
        return value[0] == zero and value[1] == zero

    def coefficient_from_expression(expression):
        polynomial = sp.Poly(expression, alpha, domain=sp.EX)
        even = zero
        odd = zero
        for (exponent,), coefficient in polynomial.terms():
            reduced = base_ring.from_sympy(coefficient)
            reduced *= (-quadratic_constant) ** (exponent // 2)
            if exponent % 2:
                odd += reduced
            else:
                even += reduced
        return even, odd

    def clean(polynomial, is_zero=coefficient_is_zero):
        return {
            monomial: coefficient
            for monomial, coefficient in polynomial.items()
            if not is_zero(coefficient)
        }

    def fiber_polynomial(polynomial):
        return clean(
            {
                monomial: coefficient_from_expression(coefficient)
                for monomial, coefficient in polynomial.terms()
            }
        )

    def monomial_key(monomial):
        return (
            monomial[0] + monomial[1],
            -monomial[1],
            -monomial[0],
        )

    def leading_monomial(polynomial):
        return max(polynomial, key=monomial_key)

    def polynomial_add(left, right, is_zero=coefficient_is_zero):
        answer = dict(left)
        for monomial, coefficient in right.items():
            answer[monomial] = coefficient_add(
                answer.get(monomial, (zero, zero)),
                coefficient,
            )
        return clean(answer, is_zero)

    def polynomial_scale(polynomial, coefficient):
        return clean(
            {
                monomial: coefficient_multiply(value, coefficient)
                for monomial, value in polynomial.items()
            }
        )

    def polynomial_monomial_multiply(polynomial, monomial):
        return {
            (
                exponent[0] + monomial[0],
                exponent[1] + monomial[1],
            ): coefficient
            for exponent, coefficient in polynomial.items()
        }

    def monomial_divides(left, right) -> bool:
        return left[0] <= right[0] and left[1] <= right[1]

    def pseudo_reduce(polynomial, basis):
        remainder = {}
        polynomial = dict(polynomial)
        steps = 0
        while polynomial:
            monomial = leading_monomial(polynomial)
            coefficient = polynomial[monomial]
            for divisor in basis:
                divisor_monomial = leading_monomial(divisor)
                if not monomial_divides(divisor_monomial, monomial):
                    continue
                divisor_coefficient = divisor[divisor_monomial]
                multiplier = (
                    monomial[0] - divisor_monomial[0],
                    monomial[1] - divisor_monomial[1],
                )
                polynomial = polynomial_add(
                    polynomial_scale(polynomial, divisor_coefficient),
                    polynomial_scale(
                        polynomial_monomial_multiply(divisor, multiplier),
                        coefficient_negate(coefficient),
                    ),
                )
                remainder = polynomial_scale(remainder, divisor_coefficient)
                steps += 1
                break
            else:
                remainder[monomial] = coefficient
                del polynomial[monomial]
        return clean(remainder), steps

    def s_polynomial(left, right, add=polynomial_add, scale=polynomial_scale):
        left_monomial = leading_monomial(left)
        right_monomial = leading_monomial(right)
        least_common_multiple = (
            max(left_monomial[0], right_monomial[0]),
            max(left_monomial[1], right_monomial[1]),
        )
        left_multiplier = (
            least_common_multiple[0] - left_monomial[0],
            least_common_multiple[1] - left_monomial[1],
        )
        right_multiplier = (
            least_common_multiple[0] - right_monomial[0],
            least_common_multiple[1] - right_monomial[1],
        )
        return add(
            scale(
                polynomial_monomial_multiply(left, left_multiplier),
                right[right_monomial],
            ),
            scale(
                polynomial_monomial_multiply(right, right_multiplier),
                coefficient_negate(left[left_monomial]),
            ),
        )

    first, second = [
        fiber_polynomial(polynomial)
        for polynomial in fiber_polynomials
    ]
    middle, first_steps = pseudo_reduce(second, [first])
    last, second_steps = pseudo_reduce(
        s_polynomial(first, middle),
        [first, middle],
    )
    basis = [first, middle, last]
    assert [leading_monomial(polynomial) for polynomial in basis] == [
        (2, 0),
        (1, 2),
        (0, 3),
    ]
    assert [len(polynomial) for polynomial in basis] == [6, 7, 6]
    assert (first_steps, second_steps) == (2, 4)

    later_fiber_polynomials = [
        fiber_polynomial(
            sp.Poly(
                sp.expand(parsed[index].subs(replacements)),
                s5,
                t4,
            )
        )
        for index in (3, 4)
    ]
    later_reductions = [
        pseudo_reduce(polynomial, basis)
        for polynomial in later_fiber_polynomials
    ]
    standard_monomials = {
        (0, 0),
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 1),
    }
    assert all(
        set(remainder) <= standard_monomials
        for remainder, _steps in later_reductions
    )

    fraction_field = sp.QQ.frac_field(s1, s3, t0, linear, t2)
    fraction_zero = fraction_field.zero

    def lift_coefficient(value):
        return fraction_field.convert(value[0]), fraction_field.convert(value[1])

    lifted_basis = [
        {
            monomial: lift_coefficient(coefficient)
            for monomial, coefficient in polynomial.items()
        }
        for polynomial in basis
    ]

    def fraction_is_zero(value) -> bool:
        return value[0] == fraction_zero and value[1] == fraction_zero

    def fraction_add(left, right):
        return left[0] + right[0], left[1] + right[1]

    def fraction_negate(value):
        return -value[0], -value[1]

    def fraction_multiply(left, right):
        return (
            left[0] * right[0]
            - quadratic_constant * left[1] * right[1],
            left[0] * right[1] + left[1] * right[0],
        )

    def fraction_inverse(value):
        denominator = (
            value[0] * value[0]
            + quadratic_constant * value[1] * value[1]
        )
        return value[0] / denominator, -value[1] / denominator

    def fraction_divide(left, right):
        return fraction_multiply(left, fraction_inverse(right))

    def fraction_polynomial_add(left, right):
        answer = dict(left)
        for monomial, coefficient in right.items():
            answer[monomial] = fraction_add(
                answer.get(monomial, (fraction_zero, fraction_zero)),
                coefficient,
            )
        return clean(answer, fraction_is_zero)

    def fraction_polynomial_scale(polynomial, coefficient):
        return clean(
            {
                monomial: fraction_multiply(value, coefficient)
                for monomial, value in polynomial.items()
            },
            fraction_is_zero,
        )

    def fraction_s_polynomial(left, right):
        return s_polynomial(
            left,
            right,
            add=fraction_polynomial_add,
            scale=fraction_polynomial_scale,
        )

    def fraction_reduce(polynomial, divisors):
        remainder = {}
        polynomial = dict(polynomial)
        steps = 0
        while polynomial:
            monomial = leading_monomial(polynomial)
            coefficient = polynomial[monomial]
            for divisor in divisors:
                divisor_monomial = leading_monomial(divisor)
                if not monomial_divides(divisor_monomial, monomial):
                    continue
                multiplier = (
                    monomial[0] - divisor_monomial[0],
                    monomial[1] - divisor_monomial[1],
                )
                factor = fraction_divide(
                    coefficient,
                    divisor[divisor_monomial],
                )
                polynomial = fraction_polynomial_add(
                    polynomial,
                    fraction_polynomial_scale(
                        polynomial_monomial_multiply(divisor, multiplier),
                        fraction_negate(factor),
                    ),
                )
                steps += 1
                break
            else:
                remainder[monomial] = coefficient
                del polynomial[monomial]
        return clean(remainder, fraction_is_zero), steps

    final_remainder, final_steps = fraction_reduce(
        fraction_s_polynomial(lifted_basis[1], lifted_basis[2]),
        lifted_basis,
    )
    assert final_remainder == {}
    assert final_steps == 5


    leading_coefficient_metadata = []
    for polynomial in basis:
        coefficient = polynomial[leading_monomial(polynomial)]
        serialized = f"{coefficient[0]}|{coefficient[1]}"
        leading_coefficient_metadata.append(
            {
                "even_terms": len(coefficient[0].terms()),
                "odd_terms": len(coefficient[1].terms()),
                "sha256": hashlib.sha256(serialized.encode()).hexdigest(),
            }
        )
    return {
        "coefficient_field": (
            "QQ(alpha)(s1,s3,t0,L,t2), alpha^2=-30420"
        ),
        "split_divisor_equation": "99*Q+155*t0^2=alpha*L",
        "input_fiber_support_sizes": [6, 8],
        "groebner_basis_support_sizes": [6, 7, 6],
        "leading_monomials": ["s5^2", "s5*t4^2", "t4^3"],
        "standard_basis": ["1", "t4", "t4^2", "s5", "s5*t4"],
        "quotient_length": 5,
        "pseudo_reduction_steps": [first_steps, second_steps],
        "final_s_pair_reduction_steps": final_steps,
        "coprime_s_pair_skipped_by_product_criterion": [
            "s5^2",
            "t4^3",
        ],
        "leading_coefficients": leading_coefficient_metadata,
        "later_moment_normal_forms": {
            str(order): {
                "support_size": len(remainder),
                "pseudo_reduction_steps": steps,
                "coefficient_term_counts_even_odd": [
                    [
                        len(coefficient[0].terms()),
                        len(coefficient[1].terms()),
                    ]
                    for _monomial, coefficient in sorted(remainder.items())
                ],
            }
            for order, (remainder, steps) in zip(
                (6, 7),
                later_reductions,
            )
        },
        "scope": (
            "exact fraction-free characteristic-zero Groebner certificate "
            "and mu6,mu7 pseudo-normal forms over the generic J-divisor "
            "function field"
        ),
    }


def modular_j_divisor_replay(
    singular: str,
    polynomials: list[str],
    prime: int,
) -> dict[str, object]:
    """Replay the generic split-J quotient over one finite field."""

    replay: dict[str, object] = {}
    for root in J_SPLIT_ROOTS[prime]:
        q_value = f"(({root}*L-155*t0^2)/99)"
        replacements = (
            ("t1", "(s1*t0-L)"),
            (
                "s2",
                f"(s1^2-(13/3)*t0^2-{q_value})",
            ),
        )
        adapted = [
            substitute(polynomial, replacements)
            for polynomial in polynomials
        ]
        completed = subprocess.run(
            [singular, "-q"],
            input=f"""
ring r=({prime},s1,s3,t0,L,t2),(s5,t4),dp;
poly p4={adapted[1]};
poly p5={adapted[2]};
ideal G=std(p4,p5);
print("QUOTIENT "+string(vdim(G))+" "+string(size(G)));
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEADEXP "+string(leadexp(G[basisIndex])));
}}
""",
            text=True,
            capture_output=True,
            check=True,
            timeout=180,
        )
        quotient_marker = re.search(
            r"(?m)^QUOTIENT (\d+) (\d+)$",
            completed.stdout,
        )
        assert quotient_marker is not None
        assert quotient_marker.groups() == ("5", "3")
        leading_exponents = re.findall(
            r"(?m)^LEADEXP (\d+),(\d+)$",
            completed.stdout,
        )
        assert leading_exponents == [("2", "0"), ("0", "3"), ("1", "2")]
        replay[str(root)] = {
            "quotient_length": 5,
            "standard_basis": ["1", "t4", "t4^2", "s5", "s5*t4"],
            "initial_ideal": ["s5^2", "t4^3", "s5*t4^2"],
        }
    return {
        "prime": prime,
        "roots_of_minus_30420": list(J_SPLIT_ROOTS[prime]),
        "split_component_replays": replay,
        "scope": (
            "generic rational-function-field reconstruction evidence "
            "on both split J components; not a characteristic-zero theorem"
        ),
    }


def divisor_certificate(
    singular: str,
    prime: int,
    polynomials: list[str],
    name: str,
    parameters: str,
    replacements: tuple[tuple[str, str], ...],
    expected_exponents: tuple[tuple[int, int], ...],
    expected_length: int,
    standard_basis: tuple[str, ...],
) -> dict[str, object]:
    restricted = [
        substitute(polynomial, replacements)
        for polynomial in polynomials[1:5]
    ]
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring r=({prime},{parameters}),(s5,t4),dp;
poly p4={restricted[0]};
poly p5={restricted[1]};
poly p6={restricted[2]};
poly p7={restricted[3]};
ideal G=std(p4,p5);
poly r6=reduce(p6,G);
poly r7=reduce(p7,G);
print("GBSIZE "+string(size(G)));
print("VDIM "+string(vdim(G)));
print("R6SIZE "+string(size(r6)));
print("R7SIZE "+string(size(r7)));
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEADEXP "+string(leadexp(G[basisIndex])));
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    summaries = [
        re.search(rf"(?m)^{label} (\d+)$", completed.stdout)
        for label in ("GBSIZE", "VDIM", "R6SIZE", "R7SIZE")
    ]
    assert all(marker is not None for marker in summaries), completed.stdout[:1000]
    values = tuple(int(marker.group(1)) for marker in summaries if marker)
    assert values == (
        len(expected_exponents),
        expected_length,
        expected_length,
        expected_length,
    )
    leading_exponents = tuple(
        tuple(int(value) for value in marker.split(","))
        for marker in re.findall(r"(?m)^LEADEXP ([0-9,]+)$", completed.stdout)
    )
    assert leading_exponents == expected_exponents
    return {
        "stratum": name,
        "groebner_basis_size": len(expected_exponents),
        "leading_exponents_s5_t4": [list(values) for values in expected_exponents],
        "quotient_length": expected_length,
        "standard_basis": list(standard_basis),
        "mu6_remainder_terms": expected_length,
        "mu7_remainder_terms": expected_length,
    }


def certificate(prime: int) -> dict[str, object]:
    singular = shutil.which("Singular")
    assert singular is not None
    expressions = [
        chart_expression(moment_terms(order, prime), 0, prime)
        for order in range(2, 8)
    ]
    variables, polynomials = prepare_s0_branch_for_msolve(
        singular,
        expressions,
        prime,
        "s0-boundary",
        120,
    )
    assert variables == ("s1", "s2", "s3", "s5", "t0", "t1", "t2", "t4")
    inverse_three = pow(3, -1, prime)
    q_value = f"(s1^2-{13 * inverse_three % prime}*t0^2)"
    return {
        "prime": prime,
        "principal_open": principal_certificate(singular, prime, polynomials),
        "special_divisors": [
            divisor_certificate(
                singular,
                prime,
                polynomials,
                "L=0",
                "s1,s2,s3,t0,t2",
                (("t1", "(s1*t0)"),),
                ((2, 0), (0, 3)),
                6,
                ("1", "t4", "t4^2", "s5", "s5*t4", "s5*t4^2"),
            ),
            divisor_certificate(
                singular,
                prime,
                polynomials,
                "Q=0",
                "s1,s3,t0,t1,t2",
                (("s2", q_value),),
                ((1, 1), (0, 3), (3, 0)),
                5,
                ("1", "t4", "t4^2", "s5", "s5^2"),
            ),
            divisor_certificate(
                singular,
                prime,
                polynomials,
                "L=Q=0",
                "s1,s3,t0,t2",
                (("t1", "(s1*t0)"), ("s2", q_value)),
                ((0, 2), (2, 1), (3, 0)),
                5,
                ("1", "t4", "s5", "s5*t4", "s5^2"),
            ),
        ],
    }


def exact_certificate() -> dict[str, object]:
    singular = shutil.which("Singular")
    assert singular is not None
    expressions = [
        exact_chart_expression(exact_moment_terms(order))
        for order in range(2, 8)
    ]
    variables, polynomials = prepare_s0_branch_for_msolve(
        singular,
        expressions,
        0,
        "s0-boundary",
        120,
    )
    assert variables == ("s1", "s2", "s3", "s5", "t0", "t1", "t2", "t4")
    return {
        "principal_open": exact_principal_certificate(singular, polynomials),
        "principal_base_resultant": exact_base_resultant_certificate(
            singular,
            polynomials,
        ),
        "j_divisor_generic_quotient": exact_j_divisor_generic_certificate(
            polynomials,
        ),
        "j_divisor_mu3_point": exact_j_divisor_point_certificate(
            singular,
            polynomials,
        ),
        "j_divisor_modular_replays": [
            modular_j_divisor_replay(singular, polynomials, prime)
            for prime in PRIMES
        ],
        "special_divisors": [
            divisor_certificate(
                singular,
                0,
                polynomials,
                "L=0",
                "s1,s2,s3,t0,t2",
                (("t1", "(s1*t0)"),),
                ((2, 0), (0, 3)),
                6,
                ("1", "t4", "t4^2", "s5", "s5*t4", "s5*t4^2"),
            ),
            divisor_certificate(
                singular,
                0,
                polynomials,
                "Q=0",
                "s1,s3,t0,t1,t2",
                (("s2", "(s1^2-(13/3)*t0^2)"),),
                ((1, 1), (0, 3), (3, 0)),
                5,
                ("1", "t4", "t4^2", "s5", "s5^2"),
            ),
            divisor_certificate(
                singular,
                0,
                polynomials,
                "L=Q=0",
                "s1,s3,t0,t2",
                (
                    ("t1", "(s1*t0)"),
                    ("s2", "(s1^2-(13/3)*t0^2)"),
                ),
                ((0, 2), (2, 1), (3, 0)),
                5,
                ("1", "t4", "s5", "s5*t4", "s5^2"),
            ),
        ],
    }


def main() -> None:
    payload = {
        "characteristic_zero_certificate": exact_certificate(),
        "modular_replays": [certificate(prime) for prime in PRIMES],
        "scope": (
            "exact characteristic-zero rational-function-field quotient "
            "certificates with two modular replays; not a full boundary unit "
            "certificate"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("PASS over QQ: generic (mu_4,mu_5) quotient has length six")
    print("PASS over QQ: mu_6 support is L*Q and mu_7 adds divisor J")
    print("PASS over QQ: the t4^3 coefficient of mu_6 has a constant t2 pivot")
    print("PASS over QQ: the same coefficient has the alternate s3 pivot H")
    print("PASS over QQ: the residual base resultant is L^6*Q^6*R63")
    print("PASS R63 is irreducible by its degree-preserving reduction modulo 47")
    print("PASS the first two residual resultants are coprime off L*Q")
    print("PASS R63 has a dense rational linear s3 subresultant pivot")
    print("PASS over QQ: the generic J-divisor quotient has length five")
    print("PASS an exact mu_3=J=0 point has quotient length five")
    print("PASS both split J components have rank-five modular replays")
    print("PASS L=0, Q=0, and L=Q=0 have quotient lengths 6, 5, and 5")
    print("PASS modular replays agree in characteristics 47 and 101")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
