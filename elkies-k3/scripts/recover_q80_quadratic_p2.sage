#!/usr/bin/env sage
"""Recover the q=80 P2 sections on the quadratic p-component exactly.

Input is the rational grevlex basis produced by ``export_q80_p2_msolve.sage``
and msolve.  Its quotient has degree eight.  Multiplication by ``p`` splits it
into the rational component and the component

    p^2-144*p+7371=0,  p=72 +/- 27*sqrt(-3).

Restricting all multiplication matrices to the two-dimensional quadratic
fiber recovers both section orientations by linear algebra.  This avoids a
fresh Groebner calculation over a number field.
"""

from sage.all import *
import argparse
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--basis", required=True)
parser.add_argument("--sign", type=int, choices=(-1, 1), default=1)
args = parser.parse_args()

names = (
    tuple(f"x{index}" for index in range(5))
    + tuple(f"y{index}" for index in range(7))
    + ("h", "p")
)
ring = PolynomialRing(QQ, names=names, order="degrevlex")

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

leading_exponents = tuple(tuple(poly.lm().exponents()[0]) for poly in groebner)
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
assert len(standard_exponents) == 8
standard_monomials = tuple(
    prod(generator**exponent for generator, exponent in zip(ring.gens(), exponents))
    for exponents in standard_exponents
)
standard_index = {
    exponents: index for index, exponents in enumerate(standard_exponents)
}


def multiplication_matrix(generator):
    answer = matrix(QQ, len(standard_monomials))
    for column, monomial in enumerate(standard_monomials):
        remainder = (generator*monomial).reduce(groebner)
        for exponents, coefficient in remainder.dict().items():
            answer[standard_index[tuple(exponents)], column] = coefficient
    return answer


multiplication = tuple(multiplication_matrix(generator) for generator in ring.gens())
p_matrix = multiplication[-1]
characteristic = p_matrix.charpoly("z")
z = characteristic.parent().gen()
assert characteristic == (z+QQ(105)/8)**4*(z**2-144*z+7371)**2

quadratic.<sqrt_minus_three> = QuadraticField(-3)
p_value = quadratic(72 + args.sign*27*sqrt_minus_three)
quadratic_multiplication = tuple(matrix(quadratic, value) for value in multiplication)
eigenspace = (
    quadratic_multiplication[-1]-p_value*identity_matrix(quadratic, 8)
).right_kernel().basis_matrix().transpose()
assert eigenspace.ncols() == 2


def restriction(matrix_value):
    target = matrix_value*eigenspace
    columns = [eigenspace.solve_right(target.column(index)) for index in range(2)]
    answer = matrix(quadratic, 2, 2, lambda row, column: columns[column][row])
    assert eigenspace*answer == target
    return answer


restricted = tuple(restriction(value) for value in quadratic_multiplication)
separator_index = next(
    index for index in range(5, 12) if not restricted[index].is_scalar()
)
separator_characteristic = restricted[separator_index].charpoly("theta")
quartic.<theta> = quadratic.extension(separator_characteristic)
restricted_quartic = tuple(matrix(quartic, value) for value in restricted)
eigenvector = (
    restricted_quartic[separator_index]-theta*identity_matrix(quartic, 2)
).right_kernel().basis()[0]
pivot = next(index for index, value in enumerate(eigenvector) if value)
coordinates = []
for matrix_value in restricted_quartic:
    image = matrix_value*eigenvector
    coordinate = image[pivot]/eigenvector[pivot]
    assert image == coordinate*eigenvector
    coordinates.append(coordinate)
solution = tuple(coordinates)

polynomials = PolynomialRing(quadratic, "T")
T = polynomials.gen()
q_value = 18-2*p_value
e_value = (p_value-42)**2/36
A = T**2*(-3+p_value*T+q_value*T**2+(p_value-42)*T**3)
B = (
    2*T**3 + (2*p_value+e_value-45)*T**4
    + (-9*p_value-4*e_value+186)*T**5
    + (12*p_value+6*e_value-299)*T**6
    + (-5*p_value-4*e_value+210)*T**7 + e_value*T**8
)
raw_discriminant = 4*A**3+27*B**2
residual = raw_discriminant // (T**7*(T-1)**4)
double_factor = gcd(residual, residual.derivative()).monic()
assert double_factor.degree() == 1
rho = -double_factor[0]
cubic_ring = PolynomialRing(quadratic, "x")
x = cubic_ring.gen()
node_factor = gcd(x**3+A(rho)*x+B(rho), 3*x**2+A(rho)).monic()
assert node_factor.degree() == 1
node = -node_factor[0]

print(
    f"Q80P2QUAD|stage=quotient|basis={len(groebner)}|degree=8|"
    f"quadratic_fiber=2|sign={args.sign}|p={p_value}",
    flush=True,
)
print(f"Q80P2QUAD|rho={rho}|node={node}|residual={residual.factor()}", flush=True)
print(
    f"Q80P2QUAD|quartic_relative_polynomial={separator_characteristic}",
    flush=True,
)
x_coefficients = solution[:5]
y_coefficients = solution[5:12]
quartic_polynomials = PolynomialRing(quartic, "t")
t = quartic_polynomials.gen()
X = sum(x_coefficients[i]*t**i for i in range(5))
Y = sum(y_coefficients[i]*t**i for i in range(7))
assert Y**2 == X**3+quartic_polynomials(A)*X+quartic_polynomials(B)
print(f"Q80P2QUAD|solution=1|X={x_coefficients}|Y={y_coefficients}", flush=True)
print("Q80P2QUAD|solution=2|X=same|Y=negative", flush=True)

# Audit the target open conditions and the resolved P1.P2 intersection.  The
# standard P1 has X=T+2*T^2; recover its polynomial square root directly.
X1 = t+2*t**2
P1_square = X1**3+quartic_polynomials(A)*X1+quartic_polynomials(B)
P1_roots = P1_square.sqrt(all=True)
assert len(P1_roots) == 2
P1Y = P1_roots[0]
rho_quartic = quartic(rho)
node_quartic = quartic(node)
passes_i2_node = X(rho_quartic) == node_quartic and Y(rho_quartic) == 0
print(
    f"Q80P2QUAD|P2_residual_I2_node={ZZ(passes_i2_node)}|P1Y={P1Y}",
    flush=True,
)
for orientation, oriented_y in ((1, Y), (-1, -Y)):
    D = X1-X
    S = P1Y+oriented_y
    H = D
    N = S**2-D**2*(X1+X)
    first = gcd(H, N)
    cancellation = gcd(H, N//first).monic()
    intersection = H.degree()-cancellation.degree()
    print(
        f"Q80P2QUAD|P1P2_orientation={orientation}|"
        f"cancellation_degree={cancellation.degree()}|intersection={intersection}",
        flush=True,
    )
print("Q80P2QUAD|status=PASS", flush=True)
