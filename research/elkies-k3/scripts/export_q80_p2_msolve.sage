#!/usr/bin/env sage
"""Export the intrinsic one-parameter q=80 P2 section scheme to msolve.

On the exact P1 locus the normalized surface has

    q=18-2*p,  e=(p-42)^2/36.

The exported equations impose a quartic polynomial x-coordinate and sextic
y-coordinate.  An inverse variable saturates away x(0)=0, the vanishing
leading coefficient, and x(1)=3; those are precisely the open conditions used
by the bounded finite-field P2 gate.  This is an algebraic section scheme, not
yet an identification with the target CM24 marking.
"""

import argparse
from pathlib import Path
from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument("--prime", type=int, default=53)
parser.add_argument("--fixed-p", default=None)
parser.add_argument("--out", required=True)
args = parser.parse_args()

if args.prime and (not is_prime(args.prime) or args.prime in (2, 3)):
    raise SystemExit("prime must be zero (QQ) or a prime greater than 3")
base = QQ if args.prime == 0 else GF(args.prime)
names = (
    tuple(f"x{index}" for index in range(5))
    + tuple(f"y{index}" for index in range(7))
    + ("h", "p")
)
ring = PolynomialRing(base, names=names, order="degrevlex")
variables = ring.gens_dict()
polynomials = PolynomialRing(ring, "T")
T = polynomials.gen()
p = variables["p"]
q = 18-2*p
e = (p-42)**2/base(36)
A = T**2*(-3+p*T+q*T**2+(p-42)*T**3)
B = (
    2*T**3 + (2*p+e-45)*T**4 + (-9*p-4*e+186)*T**5
    + (12*p+6*e-299)*T**6 + (-5*p-4*e+210)*T**7 + e*T**8
)
X = sum(variables[f"x{index}"]*T**index for index in range(5))
Y = sum(variables[f"y{index}"]*T**index for index in range(7))
identity = Y**2-X**3-A*X-B
equations = [ring(identity[index]) for index in range(13)]
open_product = variables["x0"]*variables["x4"]*(X(1)-3)
equations.append(variables["h"]*open_product-1)
if args.fixed_p is not None:
    equations.append(p-base(QQ(args.fixed_p)))

output = Path(args.out)
output.parent.mkdir(parents=True, exist_ok=True)
with output.open("w") as handle:
    handle.write(",".join(names)+"\n")
    handle.write(str(args.prime)+"\n")
    for index, equation in enumerate(equations):
        handle.write(str(equation).replace("**", "^"))
        handle.write(",\n" if index+1 < len(equations) else "\n")

print(
    f"Q80P2MSOLVE|stage=export|prime={args.prime}|variables={len(names)}|"
    f"equations={len(equations)}|fixed_p={args.fixed_p}|out={output}",
    flush=True,
)
