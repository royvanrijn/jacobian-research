#!/usr/bin/env sage
"""Recover and certify the first Q80 marked section on a quadratic cover.

The exact formal reconstruction gives the rational x-coordinate
``X1=T+x2(t)T^2``.  Substitution into the global Weierstrass equation makes
the right side a rank-one square over a quadratic extension of ``QQ(t)``.
This script derives that extension, verifies the complete polynomial section
identity, and records a squarefree hyperelliptic model of the marked cover.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--local-parameter",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-cm24-slope-8-87-qq-local-parameter.json",
)
parser.add_argument(
    "--marked-functions",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-slope-8-87-qqcm-marked-functions.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-slope-8-87-first-marked-cover-qq.json",
)
args = parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parameter_payload = json.loads(args.local_parameter.read_text())
marked_payload = json.loads(args.marked_functions.read_text())
if parameter_payload.get("schema") != "q80-cm24-formal-branch-parameter-v1":
    raise ValueError("unexpected local parameter schema")
if marked_payload.get("schema") != "q80-qqcm-marked-functions-v1":
    raise ValueError("unexpected marked-function schema")
if "X1_2" not in marked_payload["functions"]:
    raise ValueError("the rational first-section x-coordinate is unavailable")

parameter_ring = PolynomialRing(QQ, "t")
t = parameter_ring.gen()
parameter_field = parameter_ring.fraction_field()


def read_function(record, ring=parameter_ring, field=parameter_field):
    numerator = ring(record["numerator"])
    denominator = ring(record["denominator"])
    if denominator == 0 or numerator.gcd(denominator) != 1:
        raise ArithmeticError("stored rational function is not reduced")
    return field(numerator/denominator)


centers = {"D": QQ(-1)/2, "P": QQ(9)/4, "Q": QQ(-9)/4, "E": QQ(-27)/32}
d, p, q, e = tuple(
    read_function(parameter_payload["functions"][name]) + centers[name]
    for name in ("D", "P", "Q", "E")
)
x2 = read_function(marked_payload["functions"]["X1_2"])

base_ring = PolynomialRing(parameter_field, "T")
T = base_ring.gen()
r = -3*d**2+3-p-q
A = T**2*(-3+p*T+q*T**2+r*T**3)

# Derive the I4 branch jets at T=1 exactly as in the Q80 surface certificate.
s_ring = PolynomialRing(parameter_field, "s")
s = s_ring.gen()
A_one = s_ring(A(T=1+s))
correction = (A_one+3*d**2)/(-3*d**2)
branch = 2*d**3*(
    1+QQ(3)/2*correction+QQ(3)/8*correction**2-QQ(1)/16*correction**3
)
jet_matrix = Matrix(
    parameter_field, 4, 4,
    lambda row, column: s_ring((1+s)**(4+column))[row],
)
fixed = vector(
    parameter_field,
    [s_ring(2*(1+s)**3+e*(1+s)**8)[index] for index in range(4)],
)
b1, b2, b3, b4 = tuple(
    jet_matrix.solve_right(vector(parameter_field, [branch[i] for i in range(4)])-fixed)
)
B = T**3*(2+b1*T+b2*T**2+b3*T**3+b4*T**4+e*T**5)

X1 = T+x2*T**2
rhs = X1**3+A*X1+B
quotient, remainder = rhs.quo_rem(T**4)
if remainder or quotient.degree() != 4:
    raise ArithmeticError("first marked x-coordinate has the wrong component profile")

# If quotient=c0*(1+l1*T+l2*T^2)^2, adjoining v^2=c0 gives the section.
c0 = parameter_field(quotient[0])
if c0 == 0:
    raise ArithmeticError("quadratic cover coefficient vanished")
l1 = parameter_field(quotient[1]/(2*c0))
l2 = parameter_field((quotient[2]/c0-l1**2)/2)
square_polynomial = c0*(1+l1*T+l2*T**2)**2
if quotient != square_polynomial:
    raise ArithmeticError("first marked section square identity failed")

# For v^2=num/den, w=v*den gives w^2=num*den.  Remove all even powers to
# obtain the canonical squarefree branch polynomial up to a rational square.
c0_numerator = parameter_ring(c0.numerator())
c0_denominator = parameter_ring(c0.denominator())
raw_cover = c0_numerator*c0_denominator
factorization_object = raw_cover.factor()
factorization = tuple(factorization_object)
squarefree_cover = parameter_ring(factorization_object.unit())
square_part = parameter_ring.one()
for factor, exponent in factorization:
    square_part *= factor**(exponent//2)
    if exponent % 2:
        squarefree_cover *= factor
squarefree_cover = parameter_ring(squarefree_cover)
if squarefree_cover.gcd(squarefree_cover.derivative()) != 1:
    raise ArithmeticError("marked-cover model is not squarefree")
cover_degree = squarefree_cover.degree()
cover_genus = (cover_degree-1)//2

# The selected CM orientation has Y1_2=-9/8*sqrt(-6), so c0(0)=-243/32.
if c0(0) != -QQ(243)/32:
    raise ArithmeticError("first marked cover has the wrong CM24 orientation")

output = {
    "schema": "q80-first-marked-cover-qq-v1",
    "status": "PASS_EXACT_FIRST_MARKED_COVER",
    "parameter": "t with CM24 at t=0",
    "surface": "y^2=x^3+A(t,T)x+B(t,T)",
    "section": {
        "X": f"T+({x2})*T^2",
        "Y": f"T^2*v*(1+({l1})*T+({l2})*T^2)",
        "cover_relation": f"v^2={c0}",
        "identity": "Y^2=X^3+A*X+B",
    },
    "squarefree_cover": {
        "equation": f"w^2={squarefree_cover}",
        "polynomial": str(squarefree_cover),
        "degree": int(cover_degree),
        "genus": int(cover_genus),
        "squarefree": True,
        "cm24_c0": str(c0(0)),
    },
    "inputs": {
        "local_parameter": {
            "path": str(args.local_parameter.relative_to(ROOT)),
            "sha256": sha256(args.local_parameter),
        },
        "marked_functions": {
            "path": str(args.marked_functions.relative_to(ROOT)),
            "sha256": sha256(args.marked_functions),
        },
    },
    "claim_boundary": {
        "proved": [
            "exact first marked section on the displayed quadratic cover",
            "exact global section identity",
            "squarefree cover degree and genus",
        ],
        "not_proved": [
            "identification of this cover with the ICARM rank-19 source curve",
            "the remaining two generic Q80 sections",
            "the third q12 Riemann--Roch pencil",
        ],
    },
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80FIRSTMARKEDCOVER|degree={cover_degree}|genus={cover_genus}|"
    "section_identity=0|status=PASS_EXACT_FIRST_MARKED_COVER",
    flush=True,
)
