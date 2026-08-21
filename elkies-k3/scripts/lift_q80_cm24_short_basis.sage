#!/usr/bin/env sage
"""Hensel-lift the fully marked short-basis q=80 CM24 seed modulo 7.

The three CM24 basis sections all have nonzero E6 component class, so they
are represented by quadratic x and quartic y polynomials.  D1 and D2 use the
D5 spinor chart x/T=1; D1 and D3 pass through both the I4 and residual-I2
nodes.  The remaining component orientations and pairwise disjointness are
open conditions already checked on the seed.
"""

from sage.all import *
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--exponent", type=int, default=20)
args = parser.parse_args()
prime = 7

surface_names = ("d", "p", "q", "r", "b1", "b2", "b3", "b4", "e", "rho", "w")
section_names = (
    ("c1", "u10", "u11", "u12")
    + ("c2", "u20", "u21", "u22")
    + ("x30", "x31", "x32", "y30", "y31", "y32", "y33", "y34")
)
names = surface_names+section_names
ring = PolynomialRing(ZZ, names=names)
v = ring.gens_dict()
polynomials = PolynomialRing(ring, "T")
T = polynomials.gen()

A = T**2*(-3+v["p"]*T+v["q"]*T**2+v["r"]*T**3)
B = T**3*(
    2+v["b1"]*T+v["b2"]*T**2+v["b3"]*T**3+v["b4"]*T**4+v["e"]*T**5
)
discriminant = 4*A**3+27*B**2
equations = [A(1)+3*v["d"]**2, B(1)-2*v["d"]**3]
equations += [discriminant.derivative(order)(1) for order in range(1, 4)]
equations += [
    v["w"]**3+A(v["rho"])*v["w"]+B(v["rho"]),
    3*v["w"]**2+A(v["rho"]),
    discriminant.derivative()(v["rho"]),
]

X1 = T*(1+v["c1"]*T)
Y1 = T**2*(v["u10"]+v["u11"]*T+v["u12"]*T**2)
X2 = T*(1+v["c2"]*T)
Y2 = T**2*(v["u20"]+v["u21"]*T+v["u22"]*T**2)
X3 = v["x30"]+v["x31"]*T+v["x32"]*T**2
Y3 = sum(v[f"y3{index}"]*T**index for index in range(5))
for X, Y in ((X1, Y1), (X2, Y2), (X3, Y3)):
    equations += (Y**2-X**3-A*X-B).list()
for X, Y in ((X1, Y1), (X3, Y3)):
    equations += [X(1)-v["d"], Y(1)]
    equations += [X(v["rho"])-v["w"], Y(v["rho"])]
equations = tuple(ring(equation) for equation in equations if equation)

seed_values = (
    (3, 4, 3, 4, 0, 2, 6, 0, 2, 4, 1)
    + (2, 2, 1, 4)
    + (1, 2, 3, 3)
    + (4, 3, 3, 1, 2, 6, 2, 3)
)
seed = vector(ZZ, seed_values)
assert len(seed) == len(names)
field = GF(prime)
jacobian = matrix(
    field,
    [
        [equation.derivative(variable)(*seed) for variable in ring.gens()]
        for equation in equations
    ],
)
rank = jacobian.rank()
kernel = jacobian.right_kernel().dimension()
print(
    f"Q80CM24LIFT|stage=jacobian|variables={len(names)}|equations={len(equations)}|"
    f"rank={rank}|kernel={kernel}",
    flush=True,
)

point = seed
for exponent in range(1, args.exponent):
    modulus = prime**exponent
    values = vector(ZZ, [equation(*point) for equation in equations])
    assert all(value % modulus == 0 for value in values)
    rhs = vector(field, [field(-(value//modulus)) for value in values])
    if jacobian.augment(rhs).rank() != rank:
        print(
            f"Q80CM24LIFT|stage=obstruction|from_exponent={exponent}|status=NO_LIFT",
            flush=True,
        )
        raise SystemExit(1)
    correction = jacobian.solve_right(rhs)
    point += modulus*vector(ZZ, map(ZZ, correction))
    next_modulus = modulus*prime
    assert all(equation(*point) % next_modulus == 0 for equation in equations)
    print(
        f"Q80CM24LIFT|stage=lift|exponent={exponent+1}|modulus_bits={next_modulus.nbits()}",
        flush=True,
    )

modulus = prime**args.exponent
reconstructed = []
for value in point:
    try:
        reconstructed.append(ZZ(value % modulus).rational_reconstruction(modulus))
    except (ArithmeticError, ValueError):
        reconstructed.append(None)
print(
    "Q80CM24LIFT|stage=reconstruct|values="
    + ",".join(
        f"{name}:{value if value is not None else '?'}"
        for name, value in zip(names, reconstructed)
    ),
    flush=True,
)
if all(value is not None for value in reconstructed):
    candidate = vector(QQ, reconstructed)
    exact = all(equation(*candidate) == 0 for equation in equations)
    print(f"Q80CM24LIFT|stage=exact_check|pass={ZZ(exact)}", flush=True)
    if exact:
        print("Q80CM24LIFT|status=PASS", flush=True)
        raise SystemExit(0)
print("Q80CM24LIFT|status=PADIC_ONLY", flush=True)
