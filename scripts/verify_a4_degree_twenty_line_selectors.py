#!/usr/bin/env python3
"""Exact affine-line search in the degree-twenty A4 selector stratum.

This checker studies the six-dimensional root-linear valuation space

    <q0,q1,q2,q3,q4,q5>,  q5=a^2*T,

after pullback to the rational root chart.  It exactly classifies the two
fixed quadratic-factor strata and every rational affine-line component in
this full space.  A selector has a graph line ``U=m*V+n`` as a component
exactly when the coefficient matrix obtained by restricting the six strict
pullbacks to that line drops rank; vertical lines are treated separately.

An exact resultant-gcd certificate proves that the only rational affine line
is U=0.  Its selector kernel is two-dimensional: the known
[103:-16:0:8:0:0] direction and
q5-4*(3*q1+q2)=T*(a^2-4*rho), the old rational conic multiplied by T.  With
``--census-bound H`` the checker also performs a bounded rational-factor
census; ``--include-q5`` enlarges that census to the full six-dimensional
space.  With ``--conic-sieve`` it exhausts every projective quadratic over
``F_5`` and proves that the one nonstructural rank drop has no lift modulo
``5^4``.  With ``--cubic-sieve`` it factors all 3906 projective selector
members over ``F_5`` and classifies every irreducible cubic factor.  The
finite-field sieves are exact, but only constrain characteristic-zero
components whose reductions retain the tested degree and irreducibility.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import product
from math import gcd
import os
import re
import shutil
import subprocess
import sympy as sp


U, V, m, n = sp.symbols("U V m n")

H = (
    8 * U**3
    - 6 * U * V**2
    - 18 * U * V
    - 54 * U
    - 2 * V**3
    - 9 * V**2
    - 27 * V
    - 27
)
K = 4 * U**2 + 4 * U * V + 6 * U + V**2 + 3 * V + 9
M = U**2 + 2 * V**2 + 6 * V + 18
L = (
    U**3
    - 3 * U * V**2
    - 9 * U * V
    - 27 * U
    + 2 * V**3
    + 9 * V**2
    + 27 * V
    + 27
)
N1 = sp.expand(M * K)
N2 = (
    8 * U**3 * V
    + 12 * U**2 * V**2
    + 36 * U**2 * V
    + 108 * U**2
    + 6 * U * V**3
    + 36 * U * V**2
    + 108 * U * V
    + 162 * U
    + V**4
    + 9 * V**3
    + 27 * V**2
    + 54 * V
)
source_A = U**3 - V**3 - 9 * V**2 - 27 * V - 54
conic_chart_cubic = (
    U**3
    - 12 * U * V**2
    - 36 * U * V
    - 108 * U
    - 16 * V**3
    - 72 * V**2
    - 216 * V
    - 216
)
root_numerator = 3 * source_A * K**3 * L


def exact_quotient(numerator: sp.Expr, denominator: sp.Expr) -> sp.Expr:
    """Return an asserted exact polynomial quotient."""

    quotient = sp.cancel(numerator / denominator)
    assert sp.denom(quotient) == 1
    assert sp.expand(numerator - denominator * quotient) == 0
    return sp.expand(quotient)


# We clear H^5, the common denominator for q0,...,q4.  Integral scalar
# normalizations are Q1=4*q1, Q2=4*q2, Q3=3*q3, and Q4=q4.
pullback_0 = N1**3 * H**2
pullback_1 = H * (
    4 * N2 * root_numerator
    + (
        81 * N1 * N2**2
        + 243 * N1 * N2 * H
        + 729 * N1 * H**2
        - 72 * N2**3
        - 324 * N2**2 * H
        - 972 * N2 * H**2
        - 972 * H**3
    )
    * H
)
pullback_2 = (
    4 * N2**2 * root_numerator
    + 36 * root_numerator * H**2
    + (
        -243 * N1 * N2**2
        - 729 * N1 * N2 * H
        - 2187 * N1 * H**2
        + 216 * N2**3
        + 972 * N2**2 * H
        + 2916 * N2 * H**2
        + 2916 * H**3
    )
    * H**2
)
pullback_3 = H * (
    (3 * N1 + 4 * H) * root_numerator
    + (
        -54 * N1 * N2**2
        - 162 * N1 * N2 * H
        - 486 * N1 * H**2
        + 12 * N2**3
        - 324 * H**3
    )
    * H
)
pullback_4 = (
    N1 * N2 * root_numerator
    - 8 * root_numerator * H**2
    + (
        27 * N1 * N2**2
        + 81 * N1 * N2 * H
        + 243 * N1 * H**2
        - 24 * N2**3
        + 648 * H**3
    )
    * H**2
)

strict = tuple(
    exact_quotient(polynomial, K**3)
    for polynomial in (
        pullback_0,
        pullback_1,
        pullback_2,
        pullback_3,
        pullback_4,
    )
)

assert tuple(sp.Poly(f, U, V).total_degree() for f in strict) == (
    12,
    13,
    14,
    13,
    14,
)
assert sp.gcd(strict) == 1
strict_q5 = exact_quotient(N1**2 * root_numerator, K**3)
assert sp.Poly(strict_q5, U, V).total_degree() == 14
assert sp.expand(strict_q5 - 3 * source_A * M**2 * K**2 * L) == 0
full_strict = strict + (strict_q5,)
assert sp.gcd(full_strict) == 1
print("PASS: constructed the six strict degree-(12,13,14,13,14,14) pullbacks")


def divisibility_matrix(
    factor: sp.Expr,
    division_variable: sp.Symbol,
    polynomials: tuple[sp.Expr, ...] = strict,
) -> sp.Matrix:
    """Return the coefficient matrix for divisibility by a fixed factor."""

    remainders = tuple(
        sp.Poly(sp.rem(polynomial, factor, division_variable), U, V)
        for polynomial in polynomials
    )
    monomials = tuple(
        sorted(set().union(*(set(polynomial.monoms()) for polynomial in remainders)))
    )
    return sp.Matrix(
        [
            [polynomial.coeff_monomial(monomial) for polynomial in remainders]
            for monomial in monomials
        ]
    )


def assert_kernel(matrix: sp.Matrix, generators: tuple[sp.Matrix, ...]) -> None:
    """Assert that the displayed vectors span the complete kernel."""

    expected = sp.Matrix.hstack(*generators)
    assert matrix * expected == sp.zeros(matrix.rows, expected.cols)
    assert expected.rank() == len(generators)
    assert matrix.rank() + expected.rank() == matrix.cols


# The two quadratic factors which occur in the small reducible members have
# complete three-dimensional selector kernels in the full space.  Coordinates
# here use the integral basis (q0,4*q1,4*q2,3*q3,q4,q5).
rho_chart = V**2 + 3 * V + 9
selector_A_integral = sp.Matrix([0, 3, 1, 0, 0, 0])
selector_B_integral = sp.Matrix([0, 0, 0, 2, 1, 0])
selector_q0_integral = sp.Matrix([1, 0, 0, 0, 0, 0])
selector_q5_integral = sp.Matrix([0, 0, 0, 0, 0, 1])
assert_kernel(
    divisibility_matrix(K, U, full_strict),
    (selector_A_integral, selector_B_integral, selector_q5_integral),
)
assert_kernel(
    divisibility_matrix(M, U, full_strict),
    (selector_q0_integral, selector_B_integral, selector_q5_integral),
)
assert_kernel(
    divisibility_matrix(rho_chart, V, full_strict),
    (selector_A_integral,),
)

selector_A_strict = sp.expand(3 * strict[1] + strict[2])
assert sp.expand(
    selector_A_strict - 12 * rho_chart * source_A * K**3 * L
) == 0
selector_B_strict = sp.expand(2 * strict[3] + strict[4])
selector_B_residual = exact_quotient(selector_B_strict, 3 * M * K)
assert sp.Poly(selector_B_residual, U, V).total_degree() == 10
assert sp.gcd(M * K, selector_B_residual) == 1
print("PASS: the complete three-dimensional K/M divisibility kernels are exact")
print("PASS: 3*q1+q2 and 6*q3+q4 expose only the fixed quadratic factors")

restricted = tuple(
    sp.Poly(sp.expand(polynomial.subs(U, m * V + n)), V)
    for polynomial in strict
)
coefficient_matrix = sp.Matrix(
    [
        [polynomial.nth(exponent) for polynomial in restricted]
        for exponent in range(15)
    ]
)
assert coefficient_matrix.shape == (15, 5)

vertical_restricted = tuple(
    sp.Poly(sp.expand(polynomial.subs(V, n)), U)
    for polynomial in strict
)
vertical_matrix = sp.Matrix(
    [
        [polynomial.nth(exponent) for polynomial in vertical_restricted]
        for exponent in range(15)
    ]
)
assert vertical_matrix.shape == (15, 5)


def singular(expression: sp.Expr) -> str:
    """Serialize an integral SymPy polynomial for Singular."""

    return str(sp.expand(expression)).replace("**", "^")


# These thirty maximal minors are a compact exact certificate.  Their
# elimination ideal already contains every possible rational graph line.
# The complete maximal-minor ideals are then used on the two one-variable
# branches, where their standard bases are inexpensive.
certificate_rows = (
    (0, 1, 2, 3, 4),
    (0, 3, 6, 9, 12),
    (1, 4, 7, 10, 13),
    (2, 5, 8, 11, 14),
    (0, 5, 8, 11, 14),
    (1, 5, 9, 12, 14),
    (2, 4, 7, 11, 13),
    (3, 6, 9, 12, 14),
    (5, 7, 9, 11, 13),
    (6, 8, 10, 12, 14),
    (8, 9, 10, 11, 12),
    (10, 11, 12, 13, 14),
    (0, 3, 5, 11, 12),
    (0, 9, 10, 11, 13),
    (2, 7, 12, 13, 14),
    (4, 5, 7, 10, 14),
    (1, 4, 8, 9, 10),
    (0, 1, 8, 11, 12),
    (1, 3, 4, 6, 12),
    (4, 7, 8, 10, 11),
    (1, 5, 6, 12, 13),
    (5, 6, 8, 10, 13),
    (1, 3, 8, 9, 11),
    (6, 7, 10, 12, 13),
    (1, 2, 6, 7, 11),
    (0, 4, 7, 12, 13),
    (0, 2, 3, 4, 10),
    (0, 3, 6, 7, 9),
    (0, 5, 7, 9, 14),
    (3, 5, 11, 12, 14),
)
assert len(certificate_rows) == 30

determinant_program = []
for index, rows in enumerate(certificate_rows):
    submatrix_entries = ",\n".join(
        singular(coefficient_matrix[row, column])
        for row in rows
        for column in range(5)
    )
    determinant_program.append(
        f"matrix A{index}[5][5]={submatrix_entries};\n"
        f"poly f{index}=det(A{index});"
    )
determinant_program_text = "\n".join(determinant_program)
determinant_names = ",".join(
    f"f{index}" for index in range(len(certificate_rows))
)

n_zero_entries = ",\n".join(
    singular(coefficient_matrix[row, column].subs(n, 0))
    for row in range(15)
    for column in range(5)
)
vertical_entries = ",\n".join(
    singular(vertical_matrix[row, column])
    for row in range(15)
    for column in range(5)
)

singular_program = f"""
LIB "modstd.lib";
option(redSB);
ring r=0,(m,n),dp;
{determinant_program_text}
ideal I={determinant_names};
ideal G=modStd(I,1);
ideal E=eliminate(G,m);
poly h=n^2*(4*n^2+6*n+9)^8;
ideal EH=h;
if ((size(E)==1) && (reduce(h,E)==0) && (reduce(E[1],EH)==0))
{{
    print("PASS: exact graph-line elimination certificate");
}}
else
{{
    print("FAIL: exact graph-line elimination certificate");
}}

matrix C0[15][5]={n_zero_entries};
ideal I0=minor(C0,5);
ideal G0=std(I0);
ideal M2=m^2;
if ((size(G0)==1) && (reduce(m^2,G0)==0) && (reduce(G0[1],M2)==0))
{{
    print("PASS: n=0 graph-line ideal is (m^2)");
}}
else
{{
    print("FAIL: n=0 graph-line ideal");
}}

matrix CV[15][5]={vertical_entries};
ideal IV=minor(CV,5);
ideal GV=std(IV);
poly rho=n^2+3*n+9;
ideal RHO3=rho^3;
if ((size(GV)==1) && (reduce(rho^3,GV)==0) && (reduce(GV[1],RHO3)==0))
{{
    print("PASS: vertical-line ideal is ((n^2+3*n+9)^3)");
}}
else
{{
    print("FAIL: vertical-line ideal");
}}
quit;
"""

if shutil.which("Singular") is None:
    raise RuntimeError(
        "Singular is required for the exact determinantal certificate"
    )

# The stronger six-dimensional resultant certificate below supersedes this
# older q5=0 standard-basis program.  The text is retained as a compact
# independently inspectable certificate, but is not replayed a second time.

origin_matrix = coefficient_matrix.subs({m: 0, n: 0})
origin_kernel = sp.Matrix([309, -12, 0, 8, 0])
assert origin_matrix.rank() == 4
assert origin_matrix * origin_kernel == sp.zeros(15, 1)

# In the unscaled basis (q0,q1,q2,q3,q4), the vector is
# (103,-16,0,8,0).  Its common-H^5 model contains the artificial
# denominator-clearing factor H and the genuine line U=0.
selected_strict = sp.expand(
    309 * strict[0] - 12 * strict[1] + 8 * strict[3]
)
selected_residual = exact_quotient(selected_strict, 9 * U * H)
assert sp.Poly(selected_residual, U, V).total_degree() == 9
assert sp.gcd(U, selected_residual) == 1
assert sp.gcd(H, selected_residual) == 1

assert sp.discriminant(4 * n**2 + 6 * n + 9, n) == -108
assert sp.discriminant(n**2 + 3 * n + 9, n) == -27
print("PASS: [103:-16:0:8:0] gives U=0 times the known degree-nine residual")


def verify_full_affine_line_theorem() -> None:
    """Prove the rational affine-line theorem in all six directions."""

    full_restricted = tuple(
        sp.Poly(sp.expand(polynomial.subs(U, m * V + n)), V)
        for polynomial in full_strict
    )
    full_coefficient_matrix = sp.Matrix(
        [
            [polynomial.nth(exponent) for polynomial in full_restricted]
            for exponent in range(15)
        ]
    )
    assert full_coefficient_matrix.shape == (15, 6)

    full_vertical_restricted = tuple(
        sp.Poly(sp.expand(polynomial.subs(V, n)), U)
        for polynomial in full_strict
    )
    full_vertical_matrix = sp.Matrix(
        [
            [polynomial.nth(exponent) for polynomial in full_vertical_restricted]
            for exponent in range(15)
        ]
    )
    assert full_vertical_matrix.shape == (15, 6)

    # Four consecutive maximal minors suffice.  A common graph-line point
    # makes all three resultants with the first minor vanish.  Their exact
    # gcd has only n=0 and one anisotropic quadratic as support.  On n=0,
    # and for vertical lines, the gcds of the same four minors finish the
    # rational classification.
    full_certificate_rows = tuple(
        tuple(range(start, start + 6))
        for start in range(4)
    )
    full_determinant_program: list[str] = []
    for index, rows in enumerate(full_certificate_rows):
        graph_entries = ",\n".join(
            singular(full_coefficient_matrix[row, column])
            for row in rows
            for column in range(6)
        )
        vertical_entries = ",\n".join(
            singular(full_vertical_matrix[row, column])
            for row in rows
            for column in range(6)
        )
        full_determinant_program.extend((
            f"matrix F{index}[6][6]={graph_entries};",
            f"poly f{index}=det(F{index});",
            f"matrix V{index}[6][6]={vertical_entries};",
            f"poly v{index}=det(V{index});",
        ))

    full_singular_program = "\n".join((
        'option(redSB);',
        'ring r=0,(m,n),dp;',
        *full_determinant_program,
        'poly e1=resultant(f0,f1,m);',
        'poly e2=resultant(f0,f2,m);',
        'poly e3=resultant(f0,f3,m);',
        'poly egcd=gcd(gcd(e1,e2),e3);',
        'poly eh=n^16*(4*n^2+6*n+9)^81;',
        'ideal EI=egcd; ideal EH=eh;',
        'if ((reduce(egcd,EH)==0) && (reduce(eh,EI)==0))',
        '{ print("PASS: full graph-line resultant gcd"); }',
        'else { print("FAIL: full graph-line resultant gcd"); }',
        'poly g0=subst(f0,n,0); poly g1=subst(f1,n,0);',
        'poly g2=subst(f2,n,0); poly g3=subst(f3,n,0);',
        'poly ggcd=gcd(gcd(g0,g1),gcd(g2,g3));',
        'ideal GI=ggcd; ideal GH=m^4;',
        'if ((reduce(ggcd,GH)==0) && (reduce(m^4,GI)==0))',
        '{ print("PASS: full n=0 graph-line gcd is m^4"); }',
        'else { print("FAIL: full n=0 graph-line gcd"); }',
        'poly vgcd=gcd(gcd(v0,v1),gcd(v2,v3));',
        'poly vh=(n^2+3*n+9)^15;',
        'ideal VI=vgcd; ideal VH=vh;',
        'if ((reduce(vgcd,VH)==0) && (reduce(vh,VI)==0))',
        '{ print("PASS: full vertical-line gcd"); }',
        'else { print("FAIL: full vertical-line gcd"); }',
        'quit;',
    ))
    completed_full = subprocess.run(
        ["Singular", "-q"],
        input=full_singular_program,
        text=True,
        capture_output=True,
        timeout=600,
        check=True,
    )
    print(completed_full.stdout)
    assert "FAIL:" not in completed_full.stdout
    assert completed_full.stdout.count("PASS:") == 3

    full_origin_matrix = full_coefficient_matrix.subs({m: 0, n: 0})
    old_line_direction = sp.Matrix([309, -12, 0, 8, 0, 0])
    conic_line_direction = sp.Matrix([0, -3, -1, 0, 0, 1])
    assert_kernel(
        full_origin_matrix,
        (old_line_direction, conic_line_direction),
    )

    conic_line_strict = sp.expand(
        -3 * strict[1] - strict[2] + strict_q5
    )
    assert sp.expand(
        conic_line_strict
        - 3 * U * source_A * K**2 * conic_chart_cubic * L
    ) == 0
    assert sp.discriminant(4 * n**2 + 6 * n + 9, n) == -108
    assert sp.discriminant(n**2 + 3 * n + 9, n) == -27
    print("PASS: the complete U=0 selector kernel has dimension two")
    print("PASS: its q5 direction is T*(a^2-4*rho), the old conic times T")
    print("THEOREM: U=0 is the only rational affine line in the full space")


verify_full_affine_line_theorem()


def projective_parameters(
    bound: int,
    include_q5: bool,
) -> tuple[tuple[int, ...], ...]:
    """Primitive projective q-basis parameters in the symmetric height box."""

    parameters: list[tuple[int, ...]] = []
    values = range(-bound, bound + 1)
    dimension = 6 if include_q5 else 5
    for parameter in product(values, repeat=dimension):
        # Pure q0 is not an exact selector and is omitted.
        if parameter[1:] == (0,) * (dimension - 1):
            continue
        content = 0
        for value in parameter:
            content = gcd(content, abs(value))
        if content != 1:
            continue
        if next(value for value in parameter if value) < 0:
            continue
        parameters.append(parameter)
    return tuple(parameters)


def expected_reducible(parameter: tuple[int, ...]) -> bool:
    """Membership in the exact K- or M-factor selector planes."""

    first, second, third, fourth, fifth = parameter[:5]
    in_k_plane = (
        first == 0
        and second == 3 * third
        and fourth == 6 * fifth
    )
    in_m_plane = (
        second == 0
        and third == 0
        and fourth == 6 * fifth
    )
    return in_k_plane or in_m_plane


CONIC_MONOMIALS = (
    (2, 0),
    (1, 1),
    (0, 2),
    (1, 0),
    (0, 1),
    (0, 0),
)


def projective_degree_two_points(prime: int):
    """Yield normalized points of P^5 whose quadratic part is nonzero."""

    for first_nonzero in range(6):
        for tail in product(range(prime), repeat=5 - first_nonzero):
            point = (0,) * first_nonzero + (1,) + tail
            if any(point[:3]):
                yield point


def leading_monomial(polynomial: dict[tuple[int, int], int]):
    """Return the leading monomial for lexicographic U > V order."""

    return max(polynomial, key=lambda exponent: (exponent[0], exponent[1]))


def modular_remainder(
    polynomial: dict[tuple[int, int], int],
    factor: dict[tuple[int, int], int],
    prime: int,
) -> dict[tuple[int, int], int]:
    """Reduce a bivariate polynomial by one factor over a prime field."""

    dividend = {
        monomial: coefficient % prime
        for monomial, coefficient in polynomial.items()
        if coefficient % prime
    }
    factor_lead = leading_monomial(factor)
    inverse_lead = pow(factor[factor_lead], -1, prime)
    remainder: dict[tuple[int, int], int] = {}
    while dividend:
        dividend_lead = leading_monomial(dividend)
        dividend_coefficient = dividend[dividend_lead]
        if (
            dividend_lead[0] >= factor_lead[0]
            and dividend_lead[1] >= factor_lead[1]
        ):
            shift = (
                dividend_lead[0] - factor_lead[0],
                dividend_lead[1] - factor_lead[1],
            )
            scale = dividend_coefficient * inverse_lead % prime
            for monomial, coefficient in factor.items():
                shifted = (monomial[0] + shift[0], monomial[1] + shift[1])
                new_coefficient = (
                    dividend.get(shifted, 0) - scale * coefficient
                ) % prime
                if new_coefficient:
                    dividend[shifted] = new_coefficient
                else:
                    dividend.pop(shifted, None)
        else:
            remainder[dividend_lead] = dividend_coefficient
            del dividend[dividend_lead]
    return remainder


def modular_column_kernel(
    columns: list[dict[tuple[int, int], int]],
    prime: int,
) -> tuple[int, tuple[tuple[int, ...], ...]]:
    """Return rank and a deterministic kernel basis for sparse columns."""

    basis: dict[tuple[int, int], dict[tuple[int, int], int]] = {}
    basis_relations: dict[tuple[int, int], list[int]] = {}
    kernel: list[tuple[int, ...]] = []
    for column_index, column in enumerate(columns):
        working = dict(column)
        relation = [0] * len(columns)
        relation[column_index] = 1
        while working:
            pivot = leading_monomial(working)
            pivot_coefficient = working[pivot]
            if pivot not in basis:
                inverse = pow(pivot_coefficient, -1, prime)
                basis[pivot] = {
                    monomial: coefficient * inverse % prime
                    for monomial, coefficient in working.items()
                }
                basis_relations[pivot] = [
                    coefficient * inverse % prime for coefficient in relation
                ]
                break
            for monomial, coefficient in basis[pivot].items():
                new_coefficient = (
                    working.get(monomial, 0)
                    - pivot_coefficient * coefficient
                ) % prime
                if new_coefficient:
                    working[monomial] = new_coefficient
                else:
                    working.pop(monomial, None)
            relation = [
                (coefficient - pivot_coefficient * basis_coefficient) % prime
                for coefficient, basis_coefficient in zip(
                    relation, basis_relations[pivot]
                )
            ]
        else:
            kernel.append(tuple(relation))
    return len(basis), tuple(kernel)


def trim_univariate(coefficients: list[int]) -> list[int]:
    """Remove zero high coefficients while retaining the zero polynomial."""

    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def multiply_univariate(left: list[int], right: list[int]) -> list[int]:
    """Multiply dense integer polynomials stored in ascending order."""

    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        if left_coefficient:
            for right_index, right_coefficient in enumerate(right):
                answer[left_index + right_index] += (
                    left_coefficient * right_coefficient
                )
    return trim_univariate(answer)


def subtract_univariate(left: list[int], right: list[int]) -> list[int]:
    """Subtract dense integer polynomials stored in ascending order."""

    answer = [
        (left[index] if index < len(left) else 0)
        - (right[index] if index < len(right) else 0)
        for index in range(max(len(left), len(right)))
    ]
    return trim_univariate(answer)


def add_shifted_univariate(
    target: list[int], source: list[int], shift: int, scale: int
) -> None:
    """Add scale*V^shift*source to a dense polynomial in place."""

    if len(target) < len(source) + shift:
        target.extend([0] * (len(source) + shift - len(target)))
    for index, coefficient in enumerate(source):
        target[index + shift] += scale * coefficient


def integer_monic_conic_remainders(
    conic_coordinates: tuple[int, ...],
    polynomial_terms: tuple[tuple[tuple[int, int, int], ...], ...],
    width: int,
) -> tuple[tuple[int, ...], ...]:
    """Reduce all strict pullbacks by a numerical monic-U^2 conic."""

    a, b, c, d, e = conic_coordinates
    u_coefficient = [c, a]
    constant_coefficient = [e, d, b]
    maximum_u_degree = max(
        u_degree
        for polynomial in polynomial_terms
        for _, u_degree, _ in polynomial
    )
    constant_parts = [[1], [0]]
    u_parts = [[0], [1]]
    for exponent in range(1, maximum_u_degree):
        constant_parts.append(
            [
                -coefficient
                for coefficient in multiply_univariate(
                    constant_coefficient, u_parts[exponent]
                )
            ]
        )
        u_parts.append(
            subtract_univariate(
                constant_parts[exponent],
                multiply_univariate(u_coefficient, u_parts[exponent]),
            )
        )

    remainders = []
    for polynomial in polynomial_terms:
        constant_remainder = [0]
        u_remainder = [0]
        for coefficient, u_degree, v_degree in polynomial:
            add_shifted_univariate(
                constant_remainder,
                constant_parts[u_degree],
                v_degree,
                coefficient,
            )
            add_shifted_univariate(
                u_remainder,
                u_parts[u_degree],
                v_degree,
                coefficient,
            )
        assert len(constant_remainder) <= width
        assert len(u_remainder) <= width
        remainders.append(
            tuple(constant_remainder + [0] * (width - len(constant_remainder)))
            + tuple(u_remainder + [0] * (width - len(u_remainder)))
        )
    return tuple(remainders)


def integer_monic_cubic_remainders(
    cubic_coordinates: tuple[int, ...],
    polynomial_terms: tuple[tuple[tuple[int, int, int], ...], ...],
    width: int,
) -> tuple[tuple[int, ...], ...]:
    """Reduce all strict pullbacks by a numerical monic-U^3 cubic."""

    (
        u2v,
        uv2,
        v3,
        u2,
        uv,
        v2,
        u,
        v,
        constant,
    ) = cubic_coordinates
    u2_coefficient = [u2, u2v]
    u_coefficient = [u, uv, uv2]
    constant_coefficient = [constant, v, v2, v3]
    maximum_u_degree = max(
        u_degree
        for polynomial in polynomial_terms
        for _, u_degree, _ in polynomial
    )
    power_remainders = [
        ([1], [0], [0]),
        ([0], [1], [0]),
        ([0], [0], [1]),
    ]
    for exponent in range(3, maximum_u_degree + 1):
        components = []
        for component_index in range(3):
            component = [0]
            for coefficient, earlier_exponent in (
                (u2_coefficient, exponent - 1),
                (u_coefficient, exponent - 2),
                (constant_coefficient, exponent - 3),
            ):
                add_shifted_univariate(
                    component,
                    multiply_univariate(
                        coefficient,
                        power_remainders[earlier_exponent][component_index],
                    ),
                    0,
                    -1,
                )
            components.append(trim_univariate(component))
        power_remainders.append(tuple(components))

    remainders = []
    for polynomial in polynomial_terms:
        components = [[0], [0], [0]]
        for coefficient, u_degree, v_degree in polynomial:
            for component_index in range(3):
                add_shifted_univariate(
                    components[component_index],
                    power_remainders[u_degree][component_index],
                    v_degree,
                    coefficient,
                )
        assert all(len(component) <= width for component in components)
        remainders.append(tuple(
            coefficient
            for component in components
            for coefficient in component + [0] * (width - len(component))
        ))
    return tuple(remainders)


def modular_affine_solutions(
    matrix: list[list[int]], right_hand_side: list[int], prime: int
) -> tuple[tuple[int, ...], ...]:
    """Solve the exceptional lift system, enumerating its one free digit."""

    rows = [
        [coefficient % prime for coefficient in row]
        + [right_hand_side[index] % prime]
        for index, row in enumerate(matrix)
    ]
    row_index = 0
    pivots: list[int] = []
    for column_index in range(len(matrix[0])):
        pivot_index = next(
            (
                index
                for index in range(row_index, len(rows))
                if rows[index][column_index]
            ),
            None,
        )
        if pivot_index is None:
            continue
        rows[row_index], rows[pivot_index] = rows[pivot_index], rows[row_index]
        inverse = pow(rows[row_index][column_index], -1, prime)
        rows[row_index] = [
            coefficient * inverse % prime for coefficient in rows[row_index]
        ]
        for index in range(len(rows)):
            scale = rows[index][column_index]
            if index != row_index and scale:
                rows[index] = [
                    (coefficient - scale * pivot_coefficient) % prime
                    for coefficient, pivot_coefficient in zip(
                        rows[index], rows[row_index]
                    )
                ]
        pivots.append(column_index)
        row_index += 1

    variable_count = len(matrix[0])
    if any(
        not any(row[:variable_count]) and row[variable_count]
        for row in rows
    ):
        return ()
    free_variables = [
        index for index in range(variable_count) if index not in pivots
    ]
    assert free_variables == [8]
    solutions = []
    for free_value in range(prime):
        solution = [0] * variable_count
        solution[8] = free_value
        for pivot_row, pivot_column in reversed(list(enumerate(pivots))):
            solution[pivot_column] = (
                rows[pivot_row][variable_count]
                - sum(
                    rows[pivot_row][free_column] * solution[free_column]
                    for free_column in free_variables
                )
            ) % prime
        solutions.append(tuple(solution))
    return tuple(solutions)


def modular_unique_solution(
    matrix: list[list[int]], right_hand_side: list[int], prime: int
) -> tuple[int, ...] | None:
    """Solve a full-column-rank modular affine system, if consistent."""

    variable_count = len(matrix[0])
    rows = [
        [coefficient % prime for coefficient in row]
        + [right_hand_side[index] % prime]
        for index, row in enumerate(matrix)
    ]
    pivot_row = 0
    pivot_columns = []
    for column_index in range(variable_count):
        pivot_index = next(
            (
                index
                for index in range(pivot_row, len(rows))
                if rows[index][column_index]
            ),
            None,
        )
        assert pivot_index is not None
        rows[pivot_row], rows[pivot_index] = rows[pivot_index], rows[pivot_row]
        inverse = pow(rows[pivot_row][column_index], -1, prime)
        rows[pivot_row] = [
            coefficient * inverse % prime for coefficient in rows[pivot_row]
        ]
        for index in range(len(rows)):
            scale = rows[index][column_index]
            if index != pivot_row and scale:
                rows[index] = [
                    (coefficient - scale * pivot_coefficient) % prime
                    for coefficient, pivot_coefficient in zip(
                        rows[index], rows[pivot_row]
                    )
                ]
        pivot_columns.append(column_index)
        pivot_row += 1
    assert pivot_columns == list(range(variable_count))
    if any(
        not any(row[:variable_count]) and row[variable_count]
        for row in rows
    ):
        return None
    return tuple(rows[index][variable_count] for index in range(variable_count))


def modular_affine_space(
    matrix: list[list[int]], right_hand_side: list[int], prime: int
) -> tuple[tuple[int, ...], tuple[tuple[int, ...], ...]] | None:
    """Return a particular solution and nullspace basis over F_p."""

    variable_count = len(matrix[0])
    rows = [
        [coefficient % prime for coefficient in row]
        + [right_hand_side[index] % prime]
        for index, row in enumerate(matrix)
    ]
    pivot_row = 0
    pivots = []
    for column_index in range(variable_count):
        pivot_index = next(
            (
                index
                for index in range(pivot_row, len(rows))
                if rows[index][column_index]
            ),
            None,
        )
        if pivot_index is None:
            continue
        rows[pivot_row], rows[pivot_index] = rows[pivot_index], rows[pivot_row]
        inverse = pow(rows[pivot_row][column_index], -1, prime)
        rows[pivot_row] = [
            coefficient * inverse % prime for coefficient in rows[pivot_row]
        ]
        for index in range(len(rows)):
            scale = rows[index][column_index]
            if index != pivot_row and scale:
                rows[index] = [
                    (coefficient - scale * pivot_coefficient) % prime
                    for coefficient, pivot_coefficient in zip(
                        rows[index], rows[pivot_row]
                    )
                ]
        pivots.append(column_index)
        pivot_row += 1
    if any(
        not any(row[:variable_count]) and row[variable_count]
        for row in rows
    ):
        return None
    free_variables = [
        index for index in range(variable_count) if index not in pivots
    ]
    particular = [0] * variable_count
    for row_index, pivot in enumerate(pivots):
        particular[pivot] = rows[row_index][variable_count]
    basis = []
    for free_variable in free_variables:
        vector = [0] * variable_count
        vector[free_variable] = 1
        for row_index, pivot in enumerate(pivots):
            vector[pivot] = -rows[row_index][free_variable] % prime
        basis.append(tuple(vector))
    return tuple(particular), tuple(basis)


def assert_exceptional_conic_nonlift(
    polynomial_terms: tuple[tuple[tuple[int, int, int], ...], ...],
) -> None:
    """Prove that the lone nonstructural F_5 point dies modulo 5^4."""

    prime = 5
    width = 1 + max(sp.Poly(polynomial, U, V).total_degree() for polynomial in full_strict)
    base_conic_coordinates = (0, 1, 0, 3, 4)
    base_selector = (1, 4, 1, 0, 0, 0)
    base_conic = U**2 + V**2 + 3 * V + 4
    base_member = sp.expand(
        sum(
            coefficient * polynomial
            for coefficient, polynomial in zip(base_selector, full_strict)
        )
    )
    base_quotient = sp.quo(base_member, base_conic, U)

    def flat_remainder(polynomial: sp.Expr) -> list[int]:
        remainder = sp.Poly(sp.rem(polynomial, base_conic, U), U, V)
        return [
            int(remainder.coeff_monomial(U**u_degree * V**v_degree)) % prime
            for u_degree in (0, 1)
            for v_degree in range(width)
        ]

    jacobian_columns = [
        flat_remainder(-direction * base_quotient)
        for direction in (U * V, V**2, U, V, 1)
    ]
    jacobian_columns.extend(
        flat_remainder(full_strict[index]) for index in range(1, 6)
    )
    jacobian = [
        [column[row] for column in jacobian_columns]
        for row in range(2 * width)
    ]

    nodes = {base_conic_coordinates + base_selector[1:]}
    modulus = prime
    node_counts = []
    for _ in range(3):
        lifted_nodes = set()
        for node in nodes:
            remainders = integer_monic_conic_remainders(
                node[:5], polynomial_terms, width
            )
            selector = (1,) + node[5:]
            member_remainder = [
                sum(
                    selector[index] * remainders[index][coefficient_index]
                    for index in range(6)
                )
                for coefficient_index in range(2 * width)
            ]
            assert all(coefficient % modulus == 0 for coefficient in member_remainder)
            corrections = modular_affine_solutions(
                jacobian,
                [
                    (-coefficient // modulus) % prime
                    for coefficient in member_remainder
                ],
                prime,
            )
            for correction in corrections:
                lifted_nodes.add(
                    tuple(
                        (coordinate + modulus * digit) % (modulus * prime)
                        for coordinate, digit in zip(node, correction)
                    )
                )
        modulus *= prime
        nodes = lifted_nodes
        node_counts.append(len(nodes))
    assert node_counts == [5, 25, 0]
    print(
        "PASS: the exceptional F_5 conic has lift counts 5, 25, 0"
        " modulo 25, 125, 625"
    )


def finite_conic_sieve() -> None:
    """Exhaust the full projective degree-two factor chart over F_5."""

    prime = 5
    polynomial_dictionaries = tuple(
        {
            monomial: int(coefficient)
            for monomial, coefficient in sp.Poly(polynomial, U, V).terms()
        }
        for polynomial in full_strict
    )
    polynomial_terms = tuple(
        tuple(
            (coefficient, monomial[0], monomial[1])
            for monomial, coefficient in polynomial.items()
        )
        for polynomial in polynomial_dictionaries
    )
    hits = []
    point_count = 0
    for point in projective_degree_two_points(prime):
        point_count += 1
        factor = {
            monomial: coefficient
            for monomial, coefficient in zip(CONIC_MONOMIALS, point)
            if coefficient
        }
        rank, kernel = modular_column_kernel(
            [
                modular_remainder(polynomial, factor, prime)
                for polynomial in polynomial_dictionaries
            ],
            prime,
        )
        if rank < 6:
            hits.append((point, rank, kernel))
    assert point_count == 3875
    assert hits == [
        ((1, 0, 1, 0, 3, 4), 5, ((1, 4, 1, 0, 0, 0),)),
        (
            (1, 0, 2, 0, 1, 3),
            3,
            (
                (1, 0, 0, 0, 0, 0),
                (0, 0, 0, 2, 1, 0),
                (0, 0, 0, 0, 0, 1),
            ),
        ),
        (
            (1, 1, 4, 4, 2, 1),
            3,
            (
                (0, 3, 1, 0, 0, 0),
                (0, 0, 0, 2, 1, 0),
                (0, 0, 0, 0, 0, 1),
            ),
        ),
        ((0, 0, 1, 0, 3, 4), 5, ((0, 3, 1, 0, 0, 0),)),
    ]
    print(
        "PASS: all 3875 projective degree-two forms over F_5 were exhausted"
    )
    print(
        "PASS: the only rank drops are the reductions of M, K, rho_V,"
        " and one exceptional point"
    )
    assert_exceptional_conic_nonlift(polynomial_terms)
    print(
        "CERTIFIED SIEVE: a new characteristic-zero conic must reduce with"
        " lower degree or through a nonreduced K/M/rho_V neighborhood"
    )


CUBIC_MONOMIALS = (
    (3, 0),
    (2, 1),
    (1, 2),
    (0, 3),
    (2, 0),
    (1, 1),
    (0, 2),
    (1, 0),
    (0, 1),
    (0, 0),
)


def projective_points(coordinate_count: int, prime: int):
    """Yield normalized points of projective coordinate space over F_p."""

    for first_nonzero in range(coordinate_count):
        for tail in product(
            range(prime), repeat=coordinate_count - first_nonzero - 1
        ):
            yield (0,) * first_nonzero + (1,) + tail


def run_modular_cubic_factor_chunk(
    item: tuple[int, tuple[tuple[int, ...], ...]],
) -> tuple[int, tuple[str, ...]]:
    """Factor one chunk of projective selector members over F_5."""

    chunk_index, parameters = item
    lines = [
        "ring r=5,(U,V),dp;",
        *(
            f"poly R{index}={singular(polynomial)};"
            for index, polynomial in enumerate(full_strict)
        ),
        "poly F;",
        "ideal factors;",
        "int factor_index;",
    ]
    for parameter in parameters:
        expression = "+".join(
            f"({coefficient})*R{index}"
            for index, coefficient in enumerate(parameter)
            if coefficient
        )
        lines.extend((
            f"F={expression};",
            "factors=factorize(F,1);",
            "for (factor_index=1; factor_index<=size(factors);"
            " factor_index++)",
            "{",
            "  if (deg(factors[factor_index])==3)",
            "  {",
            '    print("CUBIC:"+string(factors[factor_index]));',
            "  }",
            "}",
        ))
    lines.append("quit;")
    completed = subprocess.run(
        ["Singular", "-q"],
        input="\n".join(lines),
        text=True,
        capture_output=True,
        timeout=600,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(
            "finite-field cubic factorization failed:\n"
            + completed.stdout[-2000:]
            + completed.stderr[-2000:]
        )
    factors = tuple(
        line.removeprefix("CUBIC:")
        for line in completed.stdout.splitlines()
        if line.startswith("CUBIC:")
    )
    return chunk_index, factors


def parse_singular_bivariate(expression: str) -> sp.Expr:
    """Parse Singular's compact U3+U2V notation into a SymPy expression."""

    expanded_exponents = re.sub(r"U([0-9]+)", r"U**\1", expression)
    expanded_exponents = re.sub(
        r"V([0-9]+)", r"V**\1", expanded_exponents
    )
    explicit_products = re.sub(
        r"(?<=[0-9UV])(?=[UV])", "*", expanded_exponents
    )
    return sp.expand(sp.sympify(explicit_products, locals={"U": U, "V": V}))


def normalized_modular_coordinates(
    expression: sp.Expr, monomials: tuple[tuple[int, int], ...], prime: int
) -> tuple[int, ...]:
    """Return the first-nonzero-one projective coefficient vector."""

    polynomial = sp.Poly(expression, U, V)
    coordinates = tuple(
        int(polynomial.coeff_monomial(U**u_degree * V**v_degree)) % prime
        for u_degree, v_degree in monomials
    )
    first = next(coordinate for coordinate in coordinates if coordinate)
    inverse = pow(first, -1, prime)
    return tuple(coordinate * inverse % prime for coordinate in coordinates)


def exceptional_cubic_lift_data(
    polynomial_terms: tuple[tuple[tuple[int, int, int], ...], ...],
    levels: int = 8,
) -> tuple[tuple[int, ...], tuple[int, ...], int]:
    """Lift the one nonstructural cubic residue through fixed 5-adic levels."""

    prime = 5
    width = 1 + max(
        sp.Poly(polynomial, U, V).total_degree() for polynomial in full_strict
    )
    base_factor_coordinates = (1, 0, 3, 4, 0, 4, 3, 2, 4, 4)
    base_cubic_coordinates = base_factor_coordinates[1:]
    base_selector = (0, 2, 4, 0, 0, 1)
    base_factor = sp.expand(sum(
        coefficient * U**u_degree * V**v_degree
        for coefficient, (u_degree, v_degree) in zip(
            base_factor_coordinates, CUBIC_MONOMIALS
        )
    ))
    base_member = sp.expand(sum(
        coefficient * polynomial
        for coefficient, polynomial in zip(base_selector, full_strict)
    ))
    base_quotient, base_remainder = sp.div(
        sp.Poly(base_member, U, V, modulus=prime),
        sp.Poly(base_factor, U, V, modulus=prime),
    )
    assert base_remainder.is_zero

    factor_dictionary = {
        monomial: int(coefficient) % prime
        for monomial, coefficient in sp.Poly(base_factor, U, V).terms()
    }

    def flat_remainder(expression: sp.Expr) -> list[int]:
        dictionary = modular_remainder(
            {
                monomial: int(coefficient) % prime
                for monomial, coefficient in sp.Poly(
                    expression, U, V
                ).terms()
            },
            factor_dictionary,
            prime,
        )
        return [
            dictionary.get((u_degree, v_degree), 0)
            for u_degree in range(3)
            for v_degree in range(width)
        ]

    factor_directions = tuple(
        U**u_degree * V**v_degree for u_degree, v_degree in CUBIC_MONOMIALS[1:]
    )
    columns = [
        flat_remainder(-direction * base_quotient.as_expr())
        for direction in factor_directions
    ]
    columns.extend(flat_remainder(polynomial) for polynomial in full_strict[:5])
    jacobian = [
        [column[row_index] for column in columns]
        for row_index in range(3 * width)
    ]
    assert modular_column_kernel(
        [
            {
                (row_index, 0): coefficient
                for row_index, coefficient in enumerate(column)
                if coefficient
            }
            for column in columns
        ],
        prime,
    ) == (14, ())

    node = base_cubic_coordinates + base_selector[:5]
    modulus = prime
    counts = []
    for _ in range(levels):
        remainders = integer_monic_cubic_remainders(
            node[:9], polynomial_terms, width
        )
        selector = node[9:] + (1,)
        member_remainder = [
            sum(
                selector[index] * remainders[index][coefficient_index]
                for index in range(6)
            )
            for coefficient_index in range(3 * width)
        ]
        assert all(coefficient % modulus == 0 for coefficient in member_remainder)
        correction = modular_unique_solution(
            jacobian,
            [
                (-coefficient // modulus) % prime
                for coefficient in member_remainder
            ],
            prime,
        )
        modulus *= prime
        if correction is None:
            counts.append(0)
            break
        node = tuple(
            (coordinate + modulus // prime * digit) % modulus
            for coordinate, digit in zip(node, correction)
        )
        counts.append(1)
    return tuple(counts), node, modulus


def h_transverse_slice_lift_counts(
    polynomial_terms: tuple[tuple[tuple[int, int, int], ...], ...],
    levels: int = 8,
) -> tuple[int, ...]:
    """Lift cubic factors in R3+x*R4+y*R5 near the artificial H factor."""

    prime = 5
    width = 1 + max(
        sp.Poly(polynomial, U, V).total_degree() for polynomial in full_strict
    )
    factor_coordinates = normalized_modular_coordinates(
        H, CUBIC_MONOMIALS, prime
    )
    factor_dictionary = {
        monomial: coefficient
        for monomial, coefficient in zip(CUBIC_MONOMIALS, factor_coordinates)
        if coefficient
    }
    factor_expression = sum(
        coefficient * U**u_degree * V**v_degree
        for (u_degree, v_degree), coefficient in factor_dictionary.items()
    )
    quotient, remainder = sp.div(
        sp.Poly(full_strict[3], U, V, modulus=prime),
        sp.Poly(factor_expression, U, V, modulus=prime),
    )
    assert remainder.is_zero

    def flat_remainder(expression: sp.Expr) -> list[int]:
        dictionary = modular_remainder(
            {
                monomial: int(coefficient) % prime
                for monomial, coefficient in sp.Poly(
                    expression, U, V
                ).terms()
            },
            factor_dictionary,
            prime,
        )
        return [
            dictionary.get((u_degree, v_degree), 0)
            for u_degree in range(3)
            for v_degree in range(width)
        ]

    columns = [
        flat_remainder(-U**u_degree * V**v_degree * quotient.as_expr())
        for u_degree, v_degree in CUBIC_MONOMIALS[1:]
    ]
    columns.extend(flat_remainder(full_strict[index]) for index in (4, 5))
    jacobian = [
        [column[row_index] for column in columns]
        for row_index in range(3 * width)
    ]
    base_space = modular_affine_space(
        jacobian, [0] * (3 * width), prime
    )
    assert base_space is not None
    assert len(base_space[1]) == 1
    assert base_space[1][0][-2:] == (4, 1)

    nodes = {factor_coordinates[1:] + (0, 0)}
    modulus = prime
    counts = []
    slice_terms = polynomial_terms[3:6]
    for _ in range(levels):
        lifted_nodes = set()
        for node in nodes:
            remainders = integer_monic_cubic_remainders(
                node[:9], slice_terms, width
            )
            x_coordinate, y_coordinate = node[9:]
            member_remainder = [
                remainders[0][coefficient_index]
                + x_coordinate * remainders[1][coefficient_index]
                + y_coordinate * remainders[2][coefficient_index]
                for coefficient_index in range(3 * width)
            ]
            assert all(
                coefficient % modulus == 0 for coefficient in member_remainder
            )
            affine_space = modular_affine_space(
                jacobian,
                [
                    (-coefficient // modulus) % prime
                    for coefficient in member_remainder
                ],
                prime,
            )
            if affine_space is None:
                continue
            particular, basis = affine_space
            assert len(basis) == 1
            for free_value in range(prime):
                correction = tuple(
                    (
                        particular[index]
                        + free_value * basis[0][index]
                    )
                    % prime
                    for index in range(len(node))
                )
                lifted_nodes.add(tuple(
                    (
                        coordinate
                        + modulus * correction[index]
                    )
                    % (modulus * prime)
                    for index, coordinate in enumerate(node)
                ))
        modulus *= prime
        nodes = lifted_nodes
        counts.append(len(nodes))
    return tuple(counts)


def assert_fixed_cubic_factor_tangents(
    hits: list[
        tuple[
            tuple[int, ...],
            int,
            tuple[tuple[int, ...], ...],
        ]
    ],
    polynomial_dictionaries: tuple[dict[tuple[int, int], int], ...],
) -> tuple[int, int]:
    """Audit rigid factors and isolate the artificial H-factor tangent locus."""

    prime = 5
    rigid_points = 0
    h_points = 0
    h_coordinates = normalized_modular_coordinates(
        H, CUBIC_MONOMIALS, prime
    )
    for coordinates, _, factor_kernel in hits:
        factor = {
            monomial: coefficient
            for monomial, coefficient in zip(CUBIC_MONOMIALS, coordinates)
            if coefficient
        }
        factor_expression = sum(
            coefficient * U**u_degree * V**v_degree
            for (u_degree, v_degree), coefficient in factor.items()
        )
        selector_points = set()
        for kernel_coordinates in projective_points(len(factor_kernel), prime):
            selector = tuple(
                sum(
                    kernel_coordinates[index] * factor_kernel[index][coordinate]
                    for index in range(len(factor_kernel))
                )
                % prime
                for coordinate in range(6)
            )
            first = next(coefficient for coefficient in selector if coefficient)
            inverse = pow(first, -1, prime)
            selector_points.add(tuple(
                coefficient * inverse % prime for coefficient in selector
            ))
        assert len(selector_points) == (
            prime ** len(factor_kernel) - 1
        ) // (prime - 1)

        for selector in sorted(selector_points):
            fixed_selector_index = next(
                index for index, coefficient in enumerate(selector) if coefficient
            )
            member_expression = sp.expand(sum(
                coefficient * polynomial
                for coefficient, polynomial in zip(selector, full_strict)
            ))
            quotient, remainder = sp.div(
                sp.Poly(member_expression, U, V, modulus=prime),
                sp.Poly(factor_expression, U, V, modulus=prime),
            )
            assert remainder.is_zero
            factor_columns = [
                modular_remainder(
                    {
                        monomial: int(coefficient) % prime
                        for monomial, coefficient in sp.Poly(
                            -U**u_degree * V**v_degree * quotient.as_expr(),
                            U,
                            V,
                        ).terms()
                    },
                    factor,
                    prime,
                )
                for u_degree, v_degree in CUBIC_MONOMIALS[1:]
            ]
            selector_columns = [
                modular_remainder(
                    polynomial_dictionaries[index], factor, prime
                )
                for index in range(6)
                if index != fixed_selector_index
            ]
            _, tangent_kernel = modular_column_kernel(
                factor_columns + selector_columns, prime
            )
            if coordinates == h_coordinates:
                assert len(tangent_kernel) >= len(factor_kernel) - 1
                h_points += 1
            else:
                assert len(tangent_kernel) == len(factor_kernel) - 1
                assert all(
                    all(coefficient == 0 for coefficient in relation[:9])
                    for relation in tangent_kernel
                )
                rigid_points += 1
    return rigid_points, h_points


def finite_irreducible_cubic_sieve() -> None:
    """Exhaust irreducible cubic factors of the selector plane over F_5."""

    prime = 5
    parameters = tuple(projective_points(6, prime))
    assert len(parameters) == (prime**6 - 1) // (prime - 1) == 3906
    chunk_size = 256
    chunks = tuple(
        parameters[index : index + chunk_size]
        for index in range(0, len(parameters), chunk_size)
    )
    # Keep this bounded audit polite on shared research machines.  Singular
    # already parallelizes across the two isolated processes here.
    workers = min(2, len(chunks), os.cpu_count() or 1)
    outputs: list[tuple[int, tuple[str, ...]]] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                run_modular_cubic_factor_chunk,
                (chunk_index, chunk),
            )
            for chunk_index, chunk in enumerate(chunks)
        ]
        for future in as_completed(futures):
            outputs.append(future.result())

    factor_coordinates = {
        normalized_modular_coordinates(
            parse_singular_bivariate(factor), CUBIC_MONOMIALS, prime
        )
        for _, factors in sorted(outputs)
        for factor in factors
    }
    polynomial_dictionaries = tuple(
        {
            monomial: int(coefficient)
            for monomial, coefficient in sp.Poly(polynomial, U, V).terms()
        }
        for polynomial in full_strict
    )
    polynomial_terms = tuple(
        tuple(
            (coefficient, monomial[0], monomial[1])
            for monomial, coefficient in polynomial.items()
        )
        for polynomial in polynomial_dictionaries
    )
    hits = []
    for coordinates in sorted(factor_coordinates):
        factor = {
            monomial: coefficient
            for monomial, coefficient in zip(CUBIC_MONOMIALS, coordinates)
            if coefficient
        }
        rank, kernel = modular_column_kernel(
            [
                modular_remainder(polynomial, factor, prime)
                for polynomial in polynomial_dictionaries
            ],
            prime,
        )
        assert rank < 6
        hits.append((coordinates, rank, kernel))

    expected_hits = {
        normalized_modular_coordinates(source_A, CUBIC_MONOMIALS, prime): (
            (0, 3, 1, 0, 0, 0),
            (0, 0, 0, 0, 0, 1),
        ),
        normalized_modular_coordinates(L, CUBIC_MONOMIALS, prime): (
            (0, 3, 1, 0, 0, 0),
            (3, 2, 0, 2, 1, 0),
            (0, 0, 0, 0, 0, 1),
        ),
        normalized_modular_coordinates(H, CUBIC_MONOMIALS, prime): (
            (1, 0, 0, 0, 0, 0),
            (0, 1, 0, 0, 0, 0),
            (0, 0, 0, 1, 0, 0),
        ),
        normalized_modular_coordinates(
            conic_chart_cubic, CUBIC_MONOMIALS, prime
        ): ((0, 2, 4, 0, 0, 1),),
    }
    assert hits == [
        (coordinates, 6 - len(kernel), kernel)
        for coordinates, kernel in sorted(expected_hits.items())
    ]
    rigid_tangent_points, h_tangent_points = assert_fixed_cubic_factor_tangents(
        hits, polynomial_dictionaries
    )
    assert (rigid_tangent_points, h_tangent_points) == (38, 31)

    print(
        "PASS: all 3906 projective selector members over F_5 were factored"
    )
    print(
        "PASS: the only irreducible cubic factors are A, L, H, and the"
        " known conic-chart cubic"
    )
    print(
        "PASS: all 38 non-H incidence points have no moving-factor tangent"
    )
    print(
        "PASS: every excess cubic-factor tangent lies on the artificial"
        " H-factor plane"
    )
    lift_counts, lift_node, lift_modulus = exceptional_cubic_lift_data(
        polynomial_terms
    )
    exact_lift_node = (
        0,
        -12,
        -16,
        0,
        -36,
        -72,
        -108,
        -216,
        -216,
        0,
        -3,
        -1,
        0,
        0,
    )
    balanced_lift_node = tuple(
        coordinate
        if coordinate <= lift_modulus // 2
        else coordinate - lift_modulus
        for coordinate in lift_node
    )
    assert lift_counts == (1,) * 8
    assert balanced_lift_node == exact_lift_node
    print(
        "PASS: the isolated cubic lift uniquely reconstructs the exact"
        " conic-chart factor"
    )
    h_slice_counts = h_transverse_slice_lift_counts(polynomial_terms)
    assert h_slice_counts == (5, 5, 25, 25, 125, 125, 625, 625)
    print(f"EXPERIMENT: H-transverse lift counts are {h_slice_counts}")
    print(
        "CERTIFIED SIEVE: a new characteristic-zero cubic must lose degree,"
        " become reducible, or enter the nontransverse H-adic neighborhood"
    )



def run_factor_chunk(
    item: tuple[int, tuple[tuple[int, ...], ...], bool],
) -> tuple[int, str]:
    """Factor one bounded census chunk in an isolated Singular process."""

    chunk_index, parameters, include_q5 = item
    lines = [
        "ring r=0,(U,V),dp;",
        f"poly H={singular(H)};",
    ]
    census_polynomials = full_strict if include_q5 else strict
    lines.extend(
        f"poly R{index}={singular(polynomial)};"
        for index, polynomial in enumerate(census_polynomials)
    )
    for local_index, parameter in enumerate(parameters):
        first, second, third, fourth, fifth = parameter[:5]
        coefficients: tuple[int, ...] = (
            12 * first,
            3 * second,
            3 * third,
            4 * fourth,
            12 * fifth,
        )
        if include_q5:
            coefficients += (12 * parameter[5],)
        expression = "+".join(
            f"({coefficient})*R{index}"
            for index, coefficient in enumerate(coefficients)
        )
        # The sharp plane uses only q0,q1,q3, whose common-H^5 model has
        # one artificial H factor coming from denominator clearing.
        if third == 0 and fifth == 0 and (
            not include_q5 or parameter[5] == 0
        ):
            expression = f"({expression})/H"
        expected = int(expected_reducible(parameter))
        lines.append(f"poly F{local_index}={expression};")
        lines.append(
            f"ideal I{local_index}=factorize(F{local_index},1);"
        )
        lines.append(
            f"int A{local_index}=size(I{local_index})>1;"
            f" if (A{local_index}!={expected})"
            f' {{ print("MISMATCH {parameter}"); }}'
        )
    lines.append("quit;")
    completed = subprocess.run(
        ["Singular", "-q"],
        input="\n".join(lines),
        text=True,
        capture_output=True,
        timeout=300,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(
            "bounded Singular factorization failed:\n"
            + completed.stdout[-2000:]
            + completed.stderr[-2000:]
        )
    return chunk_index, completed.stdout.strip()


def bounded_factor_census(bound: int, include_q5: bool) -> None:
    """Run an explicitly bounded rational reducibility experiment."""

    parameters = projective_parameters(bound, include_q5)
    chunk_size = 4000
    chunks = tuple(
        parameters[index : index + chunk_size]
        for index in range(0, len(parameters), chunk_size)
    )
    workers = min(8, len(chunks), os.cpu_count() or 1)
    outputs: list[tuple[int, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                run_factor_chunk,
                (chunk_index, chunk, include_q5),
            )
            for chunk_index, chunk in enumerate(chunks)
        ]
        for future in as_completed(futures):
            outputs.append(future.result())
    mismatches = "\n".join(
        output for _, output in sorted(outputs) if output
    )
    assert not mismatches, mismatches
    print(
        "EXPERIMENT: all"
        f" {len(parameters)} primitive selectors of q-height <= {bound}"
        " are rationally reducible exactly on the projectivized K/M kernels"
    )


argument_parser = argparse.ArgumentParser()
argument_parser.add_argument(
    "--census-bound",
    type=int,
    default=0,
    help="also run the bounded rational-factor census through this height",
)
argument_parser.add_argument(
    "--include-q5",
    action="store_true",
    help="include q5=a^2*T in the bounded census (not the line theorem)",
)
argument_parser.add_argument(
    "--conic-sieve",
    action="store_true",
    help="exhaust projective conic factors over F_5 and the exceptional lift",
)
argument_parser.add_argument(
    "--cubic-sieve",
    action="store_true",
    help="exhaust irreducible cubic factors over F_5 across all selectors",
)
argument_parser.add_argument(
    "--cubic-lift-only",
    action="store_true",
    help="continue the isolated exceptional cubic lift without refactoring",
)
arguments = argument_parser.parse_args()
if arguments.census_bound:
    if arguments.census_bound < 1:
        raise ValueError("--census-bound must be positive")
    bounded_factor_census(arguments.census_bound, arguments.include_q5)
if arguments.conic_sieve:
    finite_conic_sieve()
if arguments.cubic_sieve:
    finite_irreducible_cubic_sieve()
if arguments.cubic_lift_only:
    lift_polynomial_terms = tuple(
        tuple(
            (int(coefficient), monomial[0], monomial[1])
            for monomial, coefficient in sp.Poly(polynomial, U, V).terms()
        )
        for polynomial in full_strict
    )
    lift_counts, lift_node, lift_modulus = exceptional_cubic_lift_data(
        lift_polynomial_terms
    )
    print(f"EXCEPTIONAL CUBIC LIFT COUNTS: {lift_counts}")
    print(f"EXCEPTIONAL CUBIC LIFT MODULUS: {lift_modulus}")
    print(f"EXCEPTIONAL CUBIC LIFT NODE: {lift_node}")
