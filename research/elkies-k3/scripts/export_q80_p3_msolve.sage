#!/usr/bin/env sage
"""Export the simple-pole P3 scheme on the exact rational q=80 surface.

The two required node incidences are built into the ansatz.  Writing Z=T-z,

    x=X/Z^2,  y=Y/Z^3,

we solve the two linear conditions on X at T=1 and T=1/49 for its top two
coefficients, and write Y=(T-1)(T-1/49)W.  An inverse variable enforces the
same open conditions as the bounded C++ P3 search.
"""

import argparse
from pathlib import Path
from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument("--prime", type=int, default=53)
parser.add_argument("--branch", choices=("rational", "quadratic"), default="rational")
parser.add_argument("--root-index", type=int, choices=(0, 1), default=0)
parser.add_argument("--native-saturation", action="store_true")
parser.add_argument("--no-saturation", action="store_true")
parser.add_argument("--fixed-z", default=None)
parser.add_argument("--out", required=True)
args = parser.parse_args()
if args.prime and (not is_prime(args.prime) or args.prime in (2, 3, 7)):
    raise SystemExit("prime must be zero (QQ) or a prime away from 2,3,7")
base = QQ if args.prime == 0 else GF(args.prime)

names = (
    (() if args.fixed_z is not None else ("z",))
    + tuple(f"x{index}" for index in range(5))
    + tuple(f"w{index}" for index in range(8))
    + (() if args.native_saturation or args.no_saturation else ("h",))
)
ring = PolynomialRing(base, names=names, order="degrevlex")
variables = ring.gens_dict()
polynomials = PolynomialRing(ring, "T")
T = polynomials.gen()
z = base(QQ(args.fixed_z)) if args.fixed_z is not None else variables["z"]

p_variable = PolynomialRing(base, "p").gen()
if args.branch == "rational":
    p = base(-105)/8
else:
    if args.prime == 0:
        raise SystemExit("the quadratic branch currently requires a finite prime")
    p_roots = sorted(
        (p_variable**2-144*p_variable+7371).roots(multiplicities=False),
        key=ZZ,
    )
    if len(p_roots) != 2:
        raise SystemExit("quadratic p-factor does not split at this prime")
    p = base(p_roots[args.root_index])
q = 18-2*p
e = (p-42)**2/36
A = T**2*(-3+p*T+q*T**2+(p-42)*T**3)
B = (
    2*T**3 + (2*p+e-45)*T**4 + (-9*p-4*e+186)*T**5
    + (12*p+6*e-299)*T**6 + (-5*p-4*e+210)*T**7 + e*T**8
)
raw_discriminant = 4*A**3+27*B**2
residual = raw_discriminant // (T**7*(T-1)**4)
double_factor = gcd(residual, residual.derivative()).monic()
if double_factor.degree() != 1:
    raise SystemExit("surface does not have a unique residual I2 at this prime")
rho = base(-double_factor[0])
cubic_ring = PolynomialRing(base, "xnode")
xnode = cubic_ring.gen()
A_rho = base(A(rho))
B_rho = base(B(rho))
node_factor = gcd(xnode**3+A_rho*xnode+B_rho, 3*xnode**2+A_rho).monic()
if node_factor.degree() != 1:
    raise SystemExit("could not recover the residual cubic node")
node = -node_factor[0]
Z = T-z

known = sum(variables[f"x{index}"]*T**index for index in range(5))
right_one = 3*(1-z)**2-known(1)
right_rho = node*(rho-z)**2-known(rho)
top_matrix = matrix(base, [[1, 1], [rho**5, rho**6]])
x5, x6 = top_matrix.solve_right(vector(ring, [right_one, right_rho]))
X = known+x5*T**5+x6*T**6
W = sum(variables[f"w{index}"]*T**index for index in range(8))
Y = (T-1)*(T-rho)*W

identity = Y**2-X**3-A*X*Z**4-B*Z**6
equations = [ring(identity[index]) for index in range(19)]
open_product = (
    z*(z-1)*(z-rho)*variables["x0"]*x6*X(z)*W(z)
)
if args.no_saturation:
    pass
elif args.native_saturation:
    equations.append(open_product)
else:
    equations.append(variables["h"]*open_product-1)

output = Path(args.out)
output.parent.mkdir(parents=True, exist_ok=True)
with output.open("w") as handle:
    handle.write(",".join(names)+"\n")
    handle.write(str(args.prime)+"\n")
    for index, equation in enumerate(equations):
        handle.write(str(equation).replace("**", "^"))
        handle.write(",\n" if index+1 < len(equations) else "\n")

print(
    f"Q80P3MSOLVE|stage=export|prime={args.prime}|variables={len(names)}|"
    f"equations={len(equations)}|native_saturation={args.native_saturation}|"
    f"no_saturation={args.no_saturation}|fixed_z={args.fixed_z}|"
    f"branch={args.branch}|root_index={args.root_index}|p={p}|rho={rho}|node={node}|"
    f"out={output}",
    flush=True,
)
