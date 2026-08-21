#!/usr/bin/env sage
"""Rejected non-nef diagnostic for the marked CM-43 q=60 class.

This file preserves the original sparse Riemann--Roch experiment so its
full-rank output remains reproducible.  It is not the live construction route.
The divisor used below is not movable on the CM-43 surface: the exact full
chamber audit in ``verify_cm43_marked_divisor_transport.sage`` subtracts the
fixed sections ``4*(P1-P2)`` and ``P3-P2`` and leaves exactly the old fiber
``F``.  Consequently no degree enlargement or additional II* gate can recover
a new q=60 pencil from this ansatz.

In the explicit divisor basis the q=60 fiber satisfies

    D60 = 11*O + Q79 - 42*F - Theta0,

where Theta0 is the identity component of the III* (E7) fiber.  On a smooth
old Kumar fiber, ``L(11*O+Q79)`` has basis

    1,x,...,x^5, y,xy,...,x^4*y, z_Q,
    z_Q=(y+y(Q79))/(x-x(Q79)).

If ``x(Q79)=Nx/h^2`` and ``y(Q79)=Ny/h^3``, then ``z_Q`` has common
denominator ``h*(h^2*x-Nx)`` and numerator ``h^3*y+Ny``.  This script constructs
the exact sparse linear conditions that the cleared numerator vanish to
order at least 42 at the old fiber T=0 and once more on its identity
component, while remaining regular at infinity.  The optional finite-field mode is intended for fast rank
and degree-bound selection before solving the same bounded system over QQ.

This is a historical ansatz/rank analyzer.  Its output is evidence about the
rejected non-nef presentation only and must not be interpreted as progress
toward a generic q=60 fibration.
"""

import argparse

from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument("--p", type=int, default=101)
parser.add_argument(
    "--min-extra", type=int, default=0,
    help="smallest degree M for the coefficient multiplying z_Q",
)
parser.add_argument(
    "--max-extra", type=int, default=48,
    help="largest degree M for the coefficient multiplying z_Q",
)
parser.add_argument("--over-qq", action="store_true")
parser.add_argument(
    "--omit-identity-gate", action="store_true",
    help="diagnose the 42F condition before imposing the extra Theta0 gate",
)
parser.add_argument(
    "--omit-fiber-gate", action="store_true",
    help="diagnose only infinity regularity in the proposed integral basis",
)
args = parser.parse_args()

RTQ = PolynomialRing(QQ, "T")
TQ = RTQ.gen()
r = -QQ(1225)/722
s = -QQ(93312)/442225
a4q = RTQ(
    2*r*s**2*TQ**3
    - (9*r*s + 4*r**2 + 4*r + 1)/3*TQ**4
)
a6q = RTQ(
    r*s**2*(3*s + 8*r - 2)/3*TQ**5
    - (
        54*r**2*s + 81*r*s - 16*r**3
        - 24*r**2 - 12*r - 2
    )/27*TQ**6
    + r**2*TQ**7
)

# The two height-5/2 sections and the generic height-four section.
c0 = QQ(1194481)/442225
W = QQ(663613890625)/34828517376
section_data = (
    (
        QQ(684775)/93312,
        (
            QQ(914233879)/294079625,
            QQ(26371835)/1181952,
            QQ(557834834375)/11609505792,
            QQ(540596465650390625)/6499837226778624,
        ),
    ),
    (
        -QQ(1765225)/93312,
        (
            -QQ(1085766121)/294079625,
            QQ(57241835)/1181952,
            -QQ(1437996415625)/11609505792,
            QQ(540596465650390625)/6499837226778624,
        ),
    ),
)
KQ = RTQ.fraction_field()
curve = EllipticCurve(KQ, [0, 0, 0, a4q, a6q])
points = []
for c, d in section_data:
    x = RTQ(TQ**2*(c0+c*TQ+W*TQ**2))
    y = RTQ(TQ**3*sum(d[index]*TQ**index for index in range(4)))
    points.append(curve((x, y)))
quartic_b = -3*s*TQ**2
quartic_c = 3*s**2*TQ + (2*r+1)*TQ**2
quartic_d = -s**3 - 2*r*s*TQ + r*TQ**2 - 2*s*TQ
x_long = quartic_d**2/(4*s**2) - quartic_c
x3 = RTQ(x_long + quartic_c/3)
y3 = RTQ(-quartic_b*s - quartic_d*x_long/(2*s))
point3 = curve((x3, y3))
q79 = 4*points[0] - 5*points[1] + point3
Nxq = RTQ(q79[0].numerator())
Nyq = RTQ(q79[1].numerator())
hq = RTQ(q79[0].denominator().sqrt())
assert q79[0].denominator() == hq**2
assert q79[1].denominator() == hq**3
assert (Nxq.degree(), Nyq.degree(), hq.degree()) == (120, 180, 58)

if args.over_qq:
    field = QQ
    RT = RTQ
    T = TQ
    Nx, Ny, h = Nxq, Nyq, hq
    characteristic_label = "QQ"
else:
    p = ZZ(args.p)
    assert p.is_prime()
    denominators = lcm(
        [coefficient.denominator() for coefficient in Nxq]
        + [coefficient.denominator() for coefficient in Nyq]
        + [coefficient.denominator() for coefficient in hq]
    )
    assert denominators % p
    field = GF(p)
    RT = PolynomialRing(field, "T")
    T = RT.gen()

    def reduce_polynomial(poly):
        return RT([
            field(value.numerator())/field(value.denominator()) for value in poly
        ])

    Nx, Ny, h = map(reduce_polynomial, (Nxq, Nyq, hq))
    characteristic_label = str(p)

h2 = h**2
h3 = h**3


def monomial_weight(kind, exponent):
    return 4*exponent if kind == "x" else 6+4*exponent


def analyze(extra_degree):
    columns = []
    names = []

    # A column is stored sparsely as ((kind, exponent, T-degree), value).
    # The cleared numerator is
    #
    #   g*h*(h^2*x-Nx) + c*(h^3*y+Ny).
    def add_column(name, terms):
        column = {}
        for kind, exponent, polynomial in terms:
            for degree, value in enumerate(polynomial):
                if value:
                    key = (kind, exponent, degree)
                    column[key] = column.get(key, field(0)) + value
        columns.append(column)
        names.append(name)

    for degree in range(extra_degree+1):
        add_column(
            f"c[{degree}]",
            (
                ("x", 0, T**degree*Ny),
                ("y", 0, T**degree*h3),
            ),
        )

    for exponent in range(6):
        bound = extra_degree+2-monomial_weight("x", exponent)
        for degree in range(bound+1):
            add_column(
                f"a{exponent}[{degree}]",
                (
                    ("x", exponent+1, T**degree*h3),
                    ("x", exponent, -T**degree*h*Nx),
                ),
            )

    for exponent in range(5):
        bound = extra_degree+2-monomial_weight("y", exponent)
        for degree in range(bound+1):
            add_column(
                f"b{exponent}[{degree}]",
                (
                    ("y", exponent+1, T**degree*h3),
                    ("y", exponent, -T**degree*h*Nx),
                ),
            )

    # At infinity the denominator h*(h^2*x-Nx) has weight 178, so regularity
    # requires numerator weight at most 178.  Choosing the old fiber F in
    # D60=11O+Q79-42F-Theta0 to be T=0 requires divisibility by T^42.
    # The extra Theta0 condition is tested on the raw nonsingular cubic,
    # where the denominator is generically nonzero; it kills the degree-42
    # coefficient too.  The denominator vanishes at the cusp, so this last
    # condition does not incorrectly add the exceptional E7 components.
    row_keys = set()
    for column in columns:
        for kind, exponent, degree in column:
            zero_cutoff = 41 if args.omit_identity_gate else 42
            if (
                degree+monomial_weight(kind, exponent) > 178
                or (not args.omit_fiber_gate and degree <= zero_cutoff)
            ):
                row_keys.add((kind, exponent, degree))
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
    kernel_c_rank = None
    if nullity and not args.over_qq:
        kernel = system.right_kernel_matrix()
        kernel_c_rank = kernel[:, :extra_degree+1].rank()
    print(
        "CM43Q60RR|field={}|M={}|rows={}|columns={}|rank={}|nullity={}"
        "|kernel_c_rank={}".format(
            characteristic_label, extra_degree, system.nrows(), system.ncols(),
            rank, nullity, kernel_c_rank,
        ),
        flush=True,
    )
    return system, names


last = None
assert 0 <= args.min_extra <= args.max_extra
for extra_degree in range(args.min_extra, args.max_extra+1):
    last = analyze(extra_degree)

print("CM43Q60RR|status=RANK_ANALYSIS_COMPLETE", flush=True)
