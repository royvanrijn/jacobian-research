#!/usr/bin/env sage
"""Extract the p-eliminant from an msolve q=80 P2 grevlex basis.

This avoids an expensive block-order recomputation.  The supplied reduced
Groebner basis defines a degree-eight quotient.  We enumerate its standard
monomials, construct multiplication by p via exact normal forms, and print
the minimal and characteristic polynomials.
"""

import argparse
from pathlib import Path
from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--basis", required=True)
args = parser.parse_args()

names = (
    tuple(f"x{index}" for index in range(5))
    + tuple(f"y{index}" for index in range(7))
    + ("h", "p")
)
base = QQ if args.prime == 0 else GF(args.prime)
ring = PolynomialRing(base, names=names, order="degrevlex")

lines = []
inside = False
for raw_line in Path(args.basis).read_text().splitlines():
    line = raw_line.strip()
    if line.startswith("["):
        inside = True
        line = line[1:]
    if not inside or not line or line.startswith("#"):
        continue
    terminal = line.endswith("]:")
    if terminal:
        line = line[:-2]
    if line.endswith(","):
        line = line[:-1]
    if line:
        lines.append(line)
    if terminal:
        break
groebner = tuple(ring(line.replace("^", "**")) for line in lines)
assert groebner

leading_exponents = tuple(
    tuple(polynomial.lm().exponents()[0]) for polynomial in groebner
)
standard_exponents = []
empty_degrees = 0
degree = 0
while empty_degrees < 2:
    before = len(standard_exponents)
    for exponent in IntegerVectors(degree, len(names)):
        exponent = tuple(exponent)
        if not any(
            all(left >= right for left, right in zip(exponent, leading))
            for leading in leading_exponents
        ):
            standard_exponents.append(exponent)
    empty_degrees = empty_degrees+1 if len(standard_exponents) == before else 0
    degree += 1
assert len(standard_exponents) == 8, standard_exponents

standard_monomials = tuple(
    prod(generator**exponent for generator, exponent in zip(ring.gens(), exponents))
    for exponents in standard_exponents
)
index = {exponents: position for position, exponents in enumerate(standard_exponents)}
p = ring.gen(len(names)-1)
matrix_p = matrix(base, len(standard_monomials))
for column, monomial in enumerate(standard_monomials):
    remainder = (p*monomial).reduce(groebner)
    for exponents, coefficient in remainder.dict().items():
        matrix_p[index[tuple(exponents)], column] = coefficient

print(
    f"Q80P2ELIM|prime={args.prime}|basis={len(groebner)}|"
    f"quotient_degree={len(standard_monomials)}|standard={standard_monomials}",
    flush=True,
)
# Sage's optimized small-prime dense charpoly currently trips an OpenBLAS
# bug on some Apple CPUs.  For finite fields, an 8x8 integer charpoly followed
# by reduction is exact and avoids that backend entirely.  The rational FLINT
# backend is already exact and stable.
polynomial_ring = PolynomialRing(base, "z")
if args.prime:
    characteristic = polynomial_ring(matrix(ZZ, matrix_p).charpoly("z"))
else:
    characteristic = polynomial_ring(matrix_p.charpoly("z"))
squarefree = characteristic // gcd(characteristic, characteristic.derivative())
print(f"Q80P2ELIM|characteristic={characteristic.factor()}", flush=True)
print(f"Q80P2ELIM|squarefree={squarefree.factor()}", flush=True)
