#!/usr/bin/env sage
"""Certify non-liftability of the retracted standard-P1 q=80 P3 seeds.

Their P2 has the wrong residual-I2 component.  The ambient surface is
rational.  We include a variable s with s^2+6=0 and
lift the simple-pole section seed in the embedding s=10 mod 53.  Full column
rank of the overdetermined Jacobian gives a unique digit at each stage.  The
lifted coordinates are then recognized as (a+b*sqrt(-6))/c by 3D LLL and
verified against every exact equation over the number field.
"""

from itertools import product as itertools_product
import argparse
from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument("--prime", type=int, choices=(11, 29, 53), default=53)
parser.add_argument("--exponent", type=int, default=8)
args = parser.parse_args()
prime = args.prime
target_exponent = args.exponent
names = (
    ("s", "z")
    + tuple(f"x{index}" for index in range(7))
    + tuple(f"y{index}" for index in range(10))
)
ring_q = PolynomialRing(QQ, names=names)
variables_q = ring_q.gens_dict()
poly_q = PolynomialRing(ring_q, "T")
T = poly_q.gen()
s_q = variables_q["s"]
z_q = variables_q["z"]

p = QQ(-105)/8
q = 18-2*p
e = (p-42)**2/36
rho = QQ(1)/49
node = QQ(51)/2401
A = T**2*(-3+p*T+q*T**2+(p-42)*T**3)
B = (
    2*T**3+(2*p+e-45)*T**4+(-9*p-4*e+186)*T**5
    +(12*p+6*e-299)*T**6+(-5*p-4*e+210)*T**7+e*T**8
)
Z = T-z_q
X = sum(variables_q[f"x{index}"]*T**index for index in range(7))
Y = sum(variables_q[f"y{index}"]*T**index for index in range(10))
identity = Y**2-X**3-A*X*Z**4-B*Z**6
equations_q = [s_q**2+6]
equations_q += [ring_q(identity[index]) for index in range(19)]
equations_q += [
    ring_q(X(1)-3*(1-z_q)**2),
    ring_q(Y(1)),
    ring_q(X(rho)-node*(rho-z_q)**2),
    ring_q(Y(rho)),
]

ring_z = PolynomialRing(ZZ, names=names)
equations_z = []
for equation in equations_q:
    denominator = lcm(coefficient.denominator() for coefficient in equation.coefficients())
    equations_z.append(ring_z(denominator*equation))

seed_data = {
    11: (
        (4, 5)
        + (4, 7, 5, 10, 8, 0, 3)
        + (3, 1, 8, 10, 5, 3, 5, 5, 0, 4)
    ),
    29: (
        (9, 26)
        + (23, 6, 14, 14, 5, 19, 25)
        + (4, 23, 3, 21, 21, 20, 16, 16, 12, 9)
    ),
    53: (
        (10, 18)
        + (15, 14, 11, 44, 50, 2, 42)
        + (47, 34, 48, 43, 1, 7, 10, 34, 31, 10)
    ),
}
seed = vector(ZZ, seed_data[prime])
assert len(seed) == len(names)
field = GF(prime)
substitution_seed = dict(zip(ring_z.gens(), map(field, seed)))
jacobian = matrix(
    field,
    [[equation.derivative(variable).subs(substitution_seed) for variable in ring_z.gens()]
     for equation in equations_z],
)
rank = jacobian.rank()
print(
    f"Q80P3PADIC|prime={prime}|stage=jacobian|rows={jacobian.nrows()}|"
    f"columns={jacobian.ncols()}|rank={rank}",
    flush=True,
)
assert rank == len(names), "the marked p=53 seed is not an ordinary isolated point"


def exact_values(point):
    substitution = dict(zip(ring_z.gens(), point))
    return vector(ZZ, [equation.subs(substitution) for equation in equations_z])


current = seed
for exponent in range(1, target_exponent):
    modulus = prime**exponent
    values = exact_values(current)
    assert all(value % modulus == 0 for value in values)
    rhs = vector(field, [field(-(value//modulus)) for value in values])
    delta = jacobian.solve_right(rhs)
    current += modulus*vector(ZZ, [ZZ(value) for value in delta])
    assert all(value % (modulus*prime) == 0 for value in exact_values(current))
    print(
        f"Q80P3PADIC|prime={prime}|stage=lift|exponent={exponent+1}|"
        f"modulus={modulus*prime}",
        flush=True,
    )

modulus = prime**target_exponent
s_residue = ZZ(current[0] % modulus)
assert (s_residue**2+6) % modulus == 0
K.<sqrt_minus6> = QuadraticField(-6)


def recognized_values(residue, count=5, coefficient_bound=8):
    lattice = matrix(ZZ, [
        [modulus, 0, 0],
        [-s_residue, 1, 0],
        [ZZ(residue % modulus), 0, 1],
    ]).LLL()
    values = {}
    for coefficients in itertools_product(
        range(-coefficient_bound, coefficient_bound+1), repeat=3
    ):
        if coefficients == (0, 0, 0):
            continue
        a, b, c = map(ZZ, vector(ZZ, coefficients)*lattice)
        if c == 0 or c % prime == 0:
            continue
        common = gcd((a, b, c))
        a, b, c = a//common, b//common, c//common
        if c < 0:
            a, b, c = -a, -b, -c
        value = K(a+b*sqrt_minus6)/c
        score = a*a+6*b*b+c*c
        values[value] = min(score, values.get(value, infinity))
    return tuple(
        value for value, _ in sorted(values.items(), key=lambda item: item[1])[:count]
    )


choices = [recognized_values(current[index]) for index in range(1, len(names))]
print(
    f"Q80P3PADIC|prime={prime}|stage=recognition|candidate_counts="
    + ",".join(map(str, map(len, choices))),
    flush=True,
)

ring_k = PolynomialRing(K, names=names)
equations_k = [ring_k(equation) for equation in equations_q]
hits = []
# The large 53-adic modulus normally makes the shortest representative unique;
# retain five candidates per coordinate, but test the shortest joint point first.
shortest = (sqrt_minus6,) + tuple(values[0] for values in choices)
if all(equation.subs(dict(zip(ring_k.gens(), shortest))) == 0 for equation in equations_k):
    hits.append(shortest)

print(f"Q80P3PADIC|prime={prime}|stage=exact_verify|hits={len(hits)}", flush=True)
for hit in hits:
    print(f"Q80P3PADIC|prime={prime}|HIT|z={hit[1]}|X={hit[2:9]}|Y={hit[9:19]}", flush=True)
print(
    f"Q80P3PADIC|prime={prime}|SUMMARY|status="
    + ("EXACT_P3_OVER_QSQRT_MINUS6" if hits else "RECOGNITION_NEEDS_WIDER_JOINT_BOX"),
    flush=True,
)
