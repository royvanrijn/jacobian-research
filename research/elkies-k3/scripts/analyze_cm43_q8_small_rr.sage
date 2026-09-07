#!/usr/bin/env sage
"""Analyze the four-function RR space behind the CM-only marked q=8 step.

This is deliberately a small boundary calculation.  If R=P1-2*P2, then
the marked q=8 fiber has old-fiber degree four and

    D8 = 3*O + R - V,

where the effective vertical E7 divisor V has coefficients
``(2;3,5,6,9,7,6,4)`` on the affine and seven simple components.  On a
smooth old fiber, ``L(3*O+R)`` has basis ``1,x,y,z_R`` with
``z_R=(y+y(R))/(x-x(R))``.  Since R has pole six, clearing the finite base
poles only requires ``h6*z_R``.

The script first imposes the exact infinity regularity and the necessary
identity-component gate (the cleared numerator is divisible by T^2).  It
reports the resulting small kernels; exceptional-component valuations must
still be imposed before any kernel is called the q=8 pencil.
"""

import argparse
from pathlib import Path

from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument("--p", type=int, default=101)
parser.add_argument("--max-extra", type=int, default=20)
args = parser.parse_args()

load(str(Path(__file__).resolve().parent / "verify_cm43_q8_short_section.sage"))

p = ZZ(args.p)
assert p.is_prime()
field = GF(p)
RTp = PolynomialRing(field, "T")
Tp = RTp.gen()


def reduce_polynomial(poly):
    poly = RT(poly)
    assert all(QQ(value).denominator() % p for value in poly)
    return RTp([
        field(QQ(value).numerator())/field(QQ(value).denominator())
        for value in poly
    ])


Nx, Ny, h = map(
    reduce_polynomial,
    (q8_x.numerator(), q8_y.numerator(), h6),
)
assert (Nx.degree(), Ny.degree(), h.degree()) == (16, 24, 6)
h2 = h**2
h3 = h**3


def weight(kind):
    return {"one": 0, "x": 4, "y": 6}[kind]


def analyze(extra):
    columns = []

    def add(terms):
        column = {}
        for kind, polynomial in terms:
            for degree, value in enumerate(polynomial):
                if value:
                    key = (kind, degree)
                    column[key] = column.get(key, field(0))+value
        columns.append(column)

    # c*h*z_R contributes c*(h^3*y+Ny)/(h^2*x-Nx).
    for degree in range(extra+1):
        add((("one", Tp**degree*Ny), ("y", Tp**degree*h3)))
    for degree in range(extra+9):
        add((("x", Tp**degree*h2), ("one", -Tp**degree*Nx)))
    for degree in range(extra+5):
        add((("x2", Tp**degree*h2), ("x", -Tp**degree*Nx)))
    for degree in range(extra+3):
        add((("yx", Tp**degree*h2), ("y", -Tp**degree*Nx)))

    # The last two loops are b*x and d*y multiplied by h^2*x-Nx.
    # Extend the weights to their cleared-numerator monomials.
    cleared_weights = {
        "one": 0, "x": 4, "x2": 8, "y": 6, "yx": 10,
    }
    row_keys = set()
    for column in columns:
        for kind, degree in column:
            if degree+cleared_weights[kind] > 16 or degree <= 1:
                row_keys.add((kind, degree))
    row_keys = sorted(row_keys)
    row_index = {key: index for index, key in enumerate(row_keys)}
    entries = {}
    for column_index, column in enumerate(columns):
        for key, value in column.items():
            if key in row_index and value:
                entries[(row_index[key], column_index)] = value
    system = matrix(field, len(row_keys), len(columns), entries, sparse=True)
    rank = system.rank()
    nullity = system.ncols()-rank
    c_rank = None
    if nullity:
        kernel = system.right_kernel_matrix()
        c_rank = kernel[:, :extra+1].rank()
    print(
        f"CM43Q8SMALLRR|p={p}|M={extra}|rows={system.nrows()}"
        f"|columns={system.ncols()}|rank={rank}|nullity={nullity}"
        f"|kernel_c_rank={c_rank}",
        flush=True,
    )


for extra in range(args.max_extra+1):
    analyze(extra)

print("CM43Q8SMALLRR|status=NECESSARY_GATES_ONLY", flush=True)
