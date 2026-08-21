#!/usr/bin/env sage
"""Export the CM24 third-child polynomial-section scheme over GF(73).

The exact Jacobian is the one certified by
``reconstruct_q80_third_q12_jacobian_gf73.sage``.  At CM24 its saturated MW
basis has three representatives with ``P.O=0``, so their coordinates satisfy
``deg(x)<=4`` and ``deg(y)<=6``.  This script exports the thirteen coefficient
equations of ``y^2=x^3+A*x+B`` in msolve format.

The export is an exact finite-field system.  A solver result must be replayed
against the identity and component marking before it is used as an
equation-level neighbor certificate.
"""

import argparse
from pathlib import Path

from sage.all import GF, PolynomialRing


parser = argparse.ArgumentParser()
parser.add_argument("--out", required=True)
arguments = parser.parse_args()

finite = GF(73)
names = tuple(f"x{index}" for index in range(5)) + tuple(
    f"y{index}" for index in range(7)
)
ring = PolynomialRing(finite, names=names, order="degrevlex")
variables = ring.gens_dict()
polynomials = PolynomialRing(ring, "V")
V = polynomials.gen()

A = (
    6*V**8 + 16*V**7 + 47*V**6 + 33*V**5 + 58*V**4
    + 2*V**3 + 63*V**2 + 17*V + 23
)
B = (
    33*V**12 + 64*V**10 + 61*V**9 + 45*V**8 + 14*V**7
    + 20*V**6 + 54*V**5 + 8*V**4 + 50*V**3 + 57*V**2
    + 47*V + 43
)
X = sum(variables[f"x{index}"]*V**index for index in range(5))
Y = sum(variables[f"y{index}"]*V**index for index in range(7))
identity = Y**2 - X**3 - A*X - B
equations = [ring(identity[index]) for index in range(13)]
assert all(equation != 0 for equation in equations)

output = Path(arguments.out)
output.parent.mkdir(parents=True, exist_ok=True)
with output.open("w") as handle:
    handle.write(",".join(names) + "\n73\n")
    for index, equation in enumerate(equations):
        handle.write(str(equation).replace("**", "^"))
        handle.write(",\n" if index + 1 < len(equations) else "\n")

print(
    "Q80THIRDPOLYSECTIONS|prime=73|x_degree=4|y_degree=6|"
    f"variables={len(names)}|equations={len(equations)}|out={output}|"
    "status=PASS_EXPORTED",
    flush=True,
)
