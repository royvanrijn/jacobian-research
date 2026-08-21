#!/usr/bin/env sage
"""Rejected intermediate-chamber Laurent RR diagnostic at CM-43.

The raw transported isotropic class has a fixed O/E7 part.  After its exact
28-step O/E7 reduction (verified separately) the intermediate class is

    D60_intermediate = Q79 + 4*O - 43*F.

It is not nef in the full CM-43 chamber: it intersects ``P1-P2`` by ``-4``.
Subtracting ``4*(P1-P2)`` and then the fixed section ``P3-P2`` leaves the old
fiber ``F``.  This script is retained only to reproduce the earlier bounded
Laurent experiment; it cannot construct a new q=60 pencil.

The intermediate relative smooth-fiber space is

    L(4*O+Q79) = <1,x,y,x^2,z_Q>,
    z_Q=(y+y(Q79))/(x-x(Q79)).

At a zero of the pole polynomial h, z_Q has a vertical pole and h*z_Q is the
finite-place integral generator.  At infinity h*z_Q has weight 60, so the
global module also needs Laurent powers of T; restricting coefficients to
ordinary polynomials is not saturated at the E7 fiber.  This script allows a
bounded pole T^{-K}, clears the denominator h^2*x-Nx, and imposes:

* regularity at infinity (cleared-numerator weight <= 120), and
* vanishing to order 43 at T=0 (the -43F term).

Finite-field ranks in this file diagnose only the rejected intermediate
presentation.
"""

import argparse
from pathlib import Path

from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument("--p", type=int, default=101)
parser.add_argument("--min-pole", type=int, default=56)
parser.add_argument("--max-pole", type=int, default=64)
parser.add_argument(
    "--max-c-degree", type=int, default=0,
    help="largest Laurent exponent in the coefficient of h*z_Q",
)
args = parser.parse_args()

load(str(Path(__file__).resolve().parent / "verify_cm43_humbert8_anchor.sage"))

p = ZZ(args.p)
assert p.is_prime()
field = GF(p)
RTp = PolynomialRing(field, "T")
Tp = RTp.gen()

Q79 = 4*point1-5*point2+point3
Nxq = RT(Q79[0].numerator())
Nyq = RT(Q79[1].numerator())
hq = RT(Q79[0].denominator().sqrt())
assert Q79[0].denominator() == hq**2
assert Q79[1].denominator() == hq**3
assert (Nxq.degree(), Nyq.degree(), hq.degree()) == (120, 180, 58)


def reduce_polynomial(poly):
    assert all(QQ(value).denominator() % p for value in poly)
    return RTp([
        field(QQ(value).numerator())/field(QQ(value).denominator())
        for value in poly
    ])


Nx, Ny, h = map(reduce_polynomial, (Nxq, Nyq, hq))
h2 = h**2
h3 = h**3

weights = {"one": 0, "x": 4, "x2": 8, "y": 6}


def shifted_terms(polynomial, shift):
    return ((degree+shift, value) for degree, value in enumerate(polynomial))


def analyze(pole_bound):
    columns = []
    names = []

    def add(name, terms):
        column = {}
        for kind, polynomial, shift in terms:
            for degree, value in shifted_terms(polynomial, shift):
                if value:
                    key = (kind, degree)
                    column[key] = column.get(key, field(0))+value
        columns.append(column)
        names.append(name)

    M = args.max_c_degree
    # c*h*z_Q contributes c*(h^3*y+Ny)/(h^2*x-Nx).
    for degree in range(-pole_bound, M+1):
        add(
            f"c[{degree}]",
            (("one", Ny, degree), ("y", h3, degree)),
        )

    # The four regular basis elements 1,x,x^2,y.  Their upper Laurent
    # exponents are forced by the weight M+60 of c*h*z_Q.
    regular_data = (
        ("a0", "one", 0),
        ("a1", "x", 4),
        ("a2", "x2", 8),
        ("b0", "y", 6),
    )
    for prefix, kind, basis_weight in regular_data:
        upper = M+60-basis_weight
        for degree in range(-pole_bound, upper+1):
            if kind == "one":
                terms = (("x", h2, degree), ("one", -Nx, degree))
            elif kind == "x":
                terms = (("x2", h2, degree), ("x", -Nx, degree))
            elif kind == "x2":
                terms = (("x3", h2, degree), ("x2", -Nx, degree))
            else:
                terms = (("yx", h2, degree), ("y", -Nx, degree))
            add(f"{prefix}[{degree}]", terms)

    cleared_weights = {
        "one": 0, "x": 4, "x2": 8, "x3": 12, "y": 6, "yx": 10,
    }
    row_keys = set()
    for column in columns:
        for kind, degree in column:
            if degree+cleared_weights[kind] > 120 or degree < 43:
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
        c_rank = kernel[:, :pole_bound+M+1].rank()
    print(
        f"CM43Q60NEFRR|p={p}|K={pole_bound}|M={M}"
        f"|rows={system.nrows()}|columns={system.ncols()}|rank={rank}"
        f"|nullity={nullity}|kernel_c_rank={c_rank}",
        flush=True,
    )


assert 0 <= args.min_pole <= args.max_pole
for pole_bound in range(args.min_pole, args.max_pole+1):
    analyze(pole_bound)

print("CM43Q60NEFRR|status=RANK_ANALYSIS_COMPLETE", flush=True)
