#!/usr/bin/env sage
"""Export a six-variable square-recursive q=80 P3 system to msolve.

For the simple-pole section write Z=T-z and choose

    X(0)=k^2,  leading(X)=l^2,  Y(0)=sign*k^3.

The I4 and residual-I2 node incidences solve X_4 and X_5 linearly.  Four
square equations are solved upward from Y_0 and four downward from Y_9.
Only the nine middle square equations and the two node values of Y remain.
This replaces the dense 14-variable P3 system without large one-sided
recursive denominators.
"""

from sage.all import *
import argparse
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--branch", choices=("rational", "quadratic"), default="quadratic")
parser.add_argument("--root-index", type=int, choices=(0, 1), default=0)
parser.add_argument("--relative-sign", type=int, choices=(-1, 1), default=1)
parser.add_argument("--no-saturation", action="store_true")
parser.add_argument(
    "--subset",
    default=None,
    help="comma-separated indices in the 11 unsaturated equations",
)
parser.add_argument("--out", required=True)
args = parser.parse_args()

if not is_prime(args.prime) or args.prime in (2, 3, 7):
    raise SystemExit("prime must be a prime away from 2,3,7")
base = GF(args.prime)
names = ("z", "k", "l", "x1", "x2", "x3") + (() if args.no_saturation else ("h",))
ring = PolynomialRing(base, names=names, order="degrevlex")
v = ring.gens_dict()
fraction = ring.fraction_field()
polynomials = PolynomialRing(fraction, "T")
T = polynomials.gen()

p_variable = PolynomialRing(base, "p").gen()
if args.branch == "rational":
    p = base(-105)/8
else:
    roots = sorted(
        (p_variable**2-144*p_variable+7371).roots(multiplicities=False), key=ZZ
    )
    if len(roots) != 2:
        raise SystemExit("quadratic p-factor does not split at this prime")
    p = base(roots[args.root_index])
q = 18-2*p
e = (p-42)**2/36
A = T**2*(-3+p*T+q*T**2+(p-42)*T**3)
B = (
    2*T**3 + (2*p+e-45)*T**4 + (-9*p-4*e+186)*T**5
    + (12*p+6*e-299)*T**6 + (-5*p-4*e+210)*T**7 + e*T**8
)
discriminant = 4*A**3+27*B**2
residual_delta = discriminant // (T**7*(T-1)**4)
double_factor = gcd(residual_delta, residual_delta.derivative()).monic()
if double_factor.degree() != 1:
    raise SystemExit("surface does not have a unique residual I2")
rho = base(-double_factor[0])
xnode_ring = PolynomialRing(base, "xnode")
xnode = xnode_ring.gen()
A_rho, B_rho = base(A(rho)), base(B(rho))
node_factor = gcd(xnode**3+A_rho*xnode+B_rho, 3*xnode**2+A_rho).monic()
if node_factor.degree() != 1:
    raise SystemExit("could not recover residual cubic node")
node = base(-node_factor[0])

z, k, l = map(fraction, (v["z"], v["k"], v["l"]))
known = k**2+v["x1"]*T+v["x2"]*T**2+v["x3"]*T**3+l**2*T**6
rhs_one = 3*(1-z)**2-known(1)
rhs_rho = node*(rho-z)**2-known(rho)
top_matrix = matrix(fraction, [[1, 1], [rho**4, rho**5]])
x4, x5 = top_matrix.solve_right(vector(fraction, [rhs_one, rhs_rho]))
X = known+x4*T**4+x5*T**5
Z = T-z
square = X**3+A*X*Z**4+B*Z**6

y_coefficients = [fraction.zero() for _ in range(10)]
y_coefficients[0] = k**3
y_coefficients[9] = fraction(args.relative_sign)*l**3
for degree in range(1, 5):
    partial = sum(y_coefficients[index]*T**index for index in range(10))
    known_coefficient = (partial**2)[degree]
    y_coefficients[degree] = (
        square[degree]-known_coefficient
    )/(2*y_coefficients[0])
for degree in range(17, 13, -1):
    index = degree-9
    partial = sum(y_coefficients[j]*T**j for j in range(10))
    known_coefficient = (partial**2)[degree]
    y_coefficients[index] = (
        square[degree]-known_coefficient
    )/(2*y_coefficients[9])
Y = sum(y_coefficients[index]*T**index for index in range(10))
identity = Y**2-square
assert all(identity[index] == 0 for index in tuple(range(5))+tuple(range(14, 19)))

equations = [ring(identity[index].numerator()) for index in range(5, 14)]
equations += [ring(Y(1).numerator()), ring(Y(rho).numerator())]
if args.subset is not None:
    subset = tuple(map(int, args.subset.split(",")))
    if any(index < 0 or index >= len(equations) for index in subset):
        raise SystemExit("subset index outside 0,...,10")
    equations = [equations[index] for index in subset]
if not args.no_saturation:
    open_product = z*(z-1)*(z-rho)*k*l*X(z)*Y(z)
    equations.append(ring(v["h"]*open_product.numerator()-open_product.denominator()))
equations = tuple(equation for equation in equations if equation)

output = Path(args.out)
output.parent.mkdir(parents=True, exist_ok=True)
with output.open("w") as handle:
    handle.write(",".join(names)+"\n")
    handle.write(str(args.prime)+"\n")
    for index, equation in enumerate(equations):
        handle.write(str(equation).replace("**", "^"))
        handle.write(",\n" if index+1 < len(equations) else "\n")

print(
    f"Q80P3REDUCED|stage=export|prime={args.prime}|branch={args.branch}|"
    f"root_index={args.root_index}|relative_sign={args.relative_sign}|p={p}|rho={rho}|node={node}|"
    f"variables={len(names)}|equations={len(equations)}|"
    f"terms={sum(len(equation.dict()) for equation in equations)}|out={output}",
    flush=True,
)
