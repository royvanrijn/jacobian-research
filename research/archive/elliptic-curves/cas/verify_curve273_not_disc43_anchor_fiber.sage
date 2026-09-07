#!/usr/bin/env sage
"""Exclude curve 273 as a rational fiber of the explicit disc-43 anchor.

This is only a direct-fiber exclusion.  It does not exclude curve 273 from a
one-parameter K3 family containing the discriminant-43 surface as a special
member.
"""

from pathlib import Path
import sys

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "elliptic-curves" / "cas"))

from icarm_curve273 import GENERAL_WEIERSTRASS_COEFFICIENTS  # noqa: E402


curve273 = EllipticCurve(QQ, list(GENERAL_WEIERSTRASS_COEFFICIENTS))
j273 = curve273.j_invariant()

# The exact rational E6+D4+2A2+A1 model from
# elkies-k3/scripts/verify_mw2_e6d4a2a2a1_qq.sage.
R = PolynomialRing(QQ, "t")
t = R.gen()
A = R([
    0, 0, -QQ(32447500)/583443, -QQ(906250)/194481,
    QQ(31250000)/194481, -QQ(19531250)/194481,
])
B = R([
    0, 0, 0, QQ(300827000000)/2315685267,
    QQ(340001171875)/1029193452, -QQ(498857421875)/257298363,
    QQ(29541015625)/10501974, -QQ(152587890625)/85766121,
    QQ(152587890625)/343064484,
])
Delta = -16 * (4*A**3 + 27*B**2)

# Clear denominators in j(t)=j273.  Every rational smooth fiber is among the
# rational roots of this polynomial.  The t^6 factor is the singular I0*
# fiber; the remaining degree-ten factor is irreducible over QQ.
j_equation = R(
    j273.denominator() * 1728 * 4*A**3
    - j273.numerator() * (4*A**3 + 27*B**2)
)
assert j_equation.degree() == 16
assert j_equation.valuation(t) == 6
remaining = R(j_equation // t**6).monic()
assert remaining.degree() == 10 and remaining.is_irreducible()
assert j_equation.roots(QQ) == [(QQ(0), 6)]
assert Delta(0) == 0

print(
    "CURVE273DISC43|j_equation_degree=16|singular_factor=t^6"
    "|remaining_degree=10|remaining_irreducible=1|smooth_rational_roots=0",
    flush=True,
)
print("CURVE273DISC43|status=PASS|scope=direct_anchor_fibers_only", flush=True)
