#!/usr/bin/env sage
"""Export a square-recursive search for the height-21/2 H21 section.

At the fixed non-CM H21/H92 point, the E7+E8 Weierstrass model is

    y^2 = x^3 + (A1*T^3+A*T^4)*x
                    + B1*T^5+B*T^6+B2*T^7.

A section in the nonidentity E7 component with P.O=4 has the form

    x = T^2*N/Z^2,  y = T^3*M/Z^3,

with deg(Z)=4, deg(N)=10, and deg(M)=15.  The T^5 coefficient forces
N(0)=-(B1/A1)Z(0)^2.  At infinity N_10=k^2 and M_15=k^3.  We recover
M_1,...,M_7 from the low end and M_14,...,M_8 from the high end, leaving
only the middle square equations.  This is a bounded modular discovery
system; any hit must subsequently be lifted and checked exactly.
"""

from pathlib import Path
import argparse
import runpy

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"


def msolve_string(polynomial):
    return str(polynomial).replace("**", "^")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=23)
parser.add_argument(
    "--direct",
    action="store_true",
    help="retain the 15 M coefficients instead of recursively expanding them",
)
parser.add_argument(
    "--n-recursive",
    action="store_true",
    help="eliminate N_1,...,N_9 linearly from the low square coefficients",
)
parser.add_argument("--out", required=True, type=Path)
arguments = parser.parse_args()

p = ZZ(arguments.prime)
if not p.is_prime() or p in (2, 3):
    raise SystemExit("prime must be an odd prime other than 3")

namespace = runpy.run_path(str(ANCHOR))
r_value, s_value = namespace["EXPECTED_H21"]
coefficients = namespace["h21_coefficients"](r_value, s_value)
base = GF(p)
try:
    A1, A, B1, B, B2 = map(base, coefficients)
except (ZeroDivisionError, TypeError, ValueError) as error:
    raise SystemExit(f"bad reduction at p={p}: {error}")
if not A1 or not B1 or not base(-52203427).is_square():
    raise SystemExit("prime must preserve the E7 leading terms and split the orientation field")

if arguments.direct and arguments.n_recursive:
    raise SystemExit("choose at most one of --direct and --n-recursive")

core_names = (
    "d0", "d1", "d2", "d3", "k",
    "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9",
)
if arguments.n_recursive:
    names = core_names[:5] + tuple(f"m{index}" for index in range(9)) + ("h",)
else:
    names = core_names + (
        tuple(f"m{index}" for index in range(15)) if arguments.direct else ("m0",)
    ) + ("h",)
ring = PolynomialRing(base, names=names, order="degrevlex")
variables = ring.gens_dict()
fraction = ring.fraction_field()
polynomials = PolynomialRing(fraction, "T")
T = polynomials.gen()

d0, d1, d2, d3, k = map(fraction, (variables[name] for name in core_names[:5]))
m0 = fraction(variables["m0"])
Z = d0 + d1*T + d2*T**2 + d3*T**3 + T**4
n0 = -(fraction(B1)/fraction(A1))*d0**2
a4 = fraction(A1)*T**3 + fraction(A)*T**4
a6 = fraction(B1)*T**5 + fraction(B)*T**6 + fraction(B2)*T**7


def section_square(numerator):
    x_numerator = T**2*numerator
    raw = x_numerator**3 + a4*x_numerator*Z**4 + a6*Z**6
    if raw.valuation() < 6:
        raise ArithmeticError("the forced E7 cancellation did not remove T^5")
    return raw // T**6


if arguments.n_recursive:
    N = n0 + k**2*T**10
    lower_m = [fraction(variables[f"m{index}"]) for index in range(9)]
    for index in range(1, 10):
        degree = index-1
        constant = section_square(N)[degree]
        shifted = section_square(N+T**index)[degree]
        slope = shifted-constant
        if not slope:
            raise ArithmeticError(f"N_{index} is not linear in square coefficient {degree}")
        target = (sum(lower_m[j]*T**j for j in range(9))**2)[degree]
        N += ((target-constant)/slope)*T**index
else:
    N = n0 + sum(fraction(variables[f"n{index}"])*T**index for index in range(1, 10)) + k**2*T**10
X = T**2*N
raw_square = X**3 + a4*X*Z**4 + a6*Z**6
if raw_square.valuation() < 6:
    raise ArithmeticError("the forced E7 cancellation did not remove T^5")
square = raw_square // T**6
if square.degree() != 30:
    raise ArithmeticError(f"unexpected square degree {square.degree()}")

if arguments.n_recursive:
    m_coefficients = lower_m + [fraction.zero() for _ in range(7)]
    m_coefficients[15] = k**3
    for degree in range(29, 23, -1):
        index = degree-15
        partial = sum(m_coefficients[j]*T**j for j in range(16))
        known = (partial**2)[degree]
        m_coefficients[index] = (square[degree]-known)/(2*k**3)
    M = sum(m_coefficients[index]*T**index for index in range(16))
    remaining = tuple(range(9, 24))
elif arguments.direct:
    M = sum(fraction(variables[f"m{index}"])*T**index for index in range(15)) + k**3*T**15
    remaining = tuple(range(30))
else:
    m_coefficients = [fraction.zero() for _ in range(16)]
    m_coefficients[0] = m0
    m_coefficients[15] = k**3
    for degree in range(1, 8):
        partial = sum(m_coefficients[index]*T**index for index in range(16))
        known = (partial**2)[degree]
        m_coefficients[degree] = (square[degree]-known)/(2*m0)
    for degree in range(29, 22, -1):
        index = degree-15
        partial = sum(m_coefficients[j]*T**j for j in range(16))
        known = (partial**2)[degree]
        m_coefficients[index] = (square[degree]-known)/(2*k**3)
    M = sum(m_coefficients[index]*T**index for index in range(16))
    automatic = set(range(1, 8)) | set(range(23, 31))
    remaining = tuple(index for index in range(31) if index not in automatic)
identity = M**2-square
equations = [ring(identity[index].numerator()) for index in remaining]
equations = [equation for equation in equations if equation]
open_product = d0*k*m0
equations.append(ring(variables["h"]*open_product.numerator()-open_product.denominator()))

arguments.out.parent.mkdir(parents=True, exist_ok=True)
with arguments.out.open("w") as handle:
    handle.write(",".join(names)+"\n")
    handle.write(str(p)+"\n")
    for index, equation in enumerate(equations):
        handle.write(msolve_string(equation))
        handle.write(",\n" if index+1 < len(equations) else "\n")

print(
    f"H21SECTION|stage=export|prime={p}|sqrt_orientation={base(-52203427).sqrt()}|"
    f"direct={int(arguments.direct)}|"
    f"n_recursive={int(arguments.n_recursive)}|"
    f"variables={len(names)}|equations={len(equations)}|"
    f"terms={sum(len(equation.dict()) for equation in equations)}|out={arguments.out}",
    flush=True,
)
