#!/usr/bin/env python3
"""Exact affine-line search in the degree-twenty A4 selector stratum.

This checker studies the complete root-linear valuation space

    <q0,q1,q2,q3,q4>

after pullback to the rational root chart.  A selector has an affine line
``U=m*V+n`` as a component exactly when the coefficient matrix obtained by
restricting the five strict pullbacks to that line drops rank.

Exact determinantal certificates prove that the only rational affine-line
component is U=0, for the already known selector [103:-16:0:8:0].  Thus
enlarging the sharp plane by q2 and q4 produces no new rational line.  The
same checker identifies the complete K- and M-factor planes.  With
``--census-bound H`` it also performs a bounded rational-factor census;
that census is an experiment, not a classification of nonlinear rational
components.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import gcd
import os
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
print("PASS: constructed the five strict degree-(12,13,14,13,14) pullbacks")


def divisibility_matrix(
    factor: sp.Expr,
    division_variable: sp.Symbol,
) -> sp.Matrix:
    """Return the coefficient matrix for divisibility by a fixed factor."""

    remainders = tuple(
        sp.Poly(sp.rem(polynomial, factor, division_variable), U, V)
        for polynomial in strict
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
# complete two-dimensional selector kernels.  Coordinates here use the
# integral basis (q0,4*q1,4*q2,3*q3,q4).
rho_chart = V**2 + 3 * V + 9
selector_A_integral = sp.Matrix([0, 3, 1, 0, 0])
selector_B_integral = sp.Matrix([0, 0, 0, 2, 1])
selector_q0_integral = sp.Matrix([1, 0, 0, 0, 0])
assert_kernel(
    divisibility_matrix(K, U),
    (selector_A_integral, selector_B_integral),
)
assert_kernel(
    divisibility_matrix(M, U),
    (selector_q0_integral, selector_B_integral),
)
assert_kernel(
    divisibility_matrix(rho_chart, V),
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
print("PASS: the complete K- and M-divisibility selector planes are exact")
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

completed = subprocess.run(
    ["Singular", "-q"],
    input=singular_program,
    text=True,
    capture_output=True,
    check=True,
)
print(completed.stdout)
assert "FAIL:" not in completed.stdout
assert completed.stdout.count("PASS:") == 3

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
print("PASS: the unique rational line has selector [103:-16:0:8:0]")
print("PASS: it is U=0 times the known degree-nine residual root-chart curve")
print("THEOREM: q2 and q4 add no rational affine-line selector")


def projective_parameters(bound: int) -> tuple[tuple[int, ...], ...]:
    """Primitive projective q-basis parameters in the symmetric height box."""

    parameters: list[tuple[int, ...]] = []
    values = range(-bound, bound + 1)
    for first in values:
        for second in values:
            for third in values:
                for fourth in values:
                    for fifth in values:
                        parameter = (first, second, third, fourth, fifth)
                        # Pure q0 is not an exact selector and is omitted.
                        if parameter[1:] == (0, 0, 0, 0):
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

    first, second, third, fourth, fifth = parameter
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


def run_factor_chunk(
    item: tuple[int, tuple[tuple[int, ...], ...]],
) -> tuple[int, str]:
    """Factor one bounded census chunk in an isolated Singular process."""

    chunk_index, parameters = item
    lines = [
        "ring r=0,(U,V),dp;",
        f"poly H={singular(H)};",
    ]
    lines.extend(
        f"poly R{index}={singular(polynomial)};"
        for index, polynomial in enumerate(strict)
    )
    for local_index, parameter in enumerate(parameters):
        first, second, third, fourth, fifth = parameter
        coefficients = (
            12 * first,
            3 * second,
            3 * third,
            4 * fourth,
            12 * fifth,
        )
        expression = "+".join(
            f"({coefficient})*R{index}"
            for index, coefficient in enumerate(coefficients)
        )
        # The sharp plane uses only q0,q1,q3, whose common-H^5 model has
        # one artificial H factor coming from denominator clearing.
        if third == 0 and fifth == 0:
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


def bounded_factor_census(bound: int) -> None:
    """Run an explicitly bounded rational reducibility experiment."""

    parameters = projective_parameters(bound)
    chunk_size = 4000
    chunks = tuple(
        parameters[index : index + chunk_size]
        for index in range(0, len(parameters), chunk_size)
    )
    workers = min(8, len(chunks), os.cpu_count() or 1)
    outputs: list[tuple[int, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(run_factor_chunk, item)
            for item in enumerate(chunks)
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
        " are rationally reducible exactly on the K/M planes"
    )


argument_parser = argparse.ArgumentParser()
argument_parser.add_argument(
    "--census-bound",
    type=int,
    default=0,
    help="also run the bounded rational-factor census through this height",
)
arguments = argument_parser.parse_args()
if arguments.census_bound:
    if arguments.census_bound < 1:
        raise ValueError("--census-bound must be positive")
    bounded_factor_census(arguments.census_bound)
