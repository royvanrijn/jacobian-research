#!/usr/bin/env sage
"""Verify the short marked section underlying the optimal CM-43 q=8 step.

The complete marked q=8 orbit calculation gives horizontal MW coordinates
(1,-2,0) in the explicit basis (P1,P2,P3).  This script evaluates
R=P1-2*P2 by exact function-field group law on the CM-43 Kumar equation and
certifies its pole-six form.  It is the equation-level starting datum for the
q=8 neighbor, and is far smaller than the pole-58 level-79 section.
"""

from pathlib import Path

from sage.all import *


BASE = Path(__file__).resolve().parent
load(str(BASE / "verify_cm43_humbert8_anchor.sage"))

# Variables point1, point2 and height43 are defined and exactly verified by
# the anchor script.  Both P1 and P2 meet the nonidentity E7 component.
q8_coordinates = vector(ZZ, (1, -2, 0))
q8_height = q8_coordinates*height43*q8_coordinates
assert q8_height == QQ(29)/2
e7_local_correction = QQ(3)/2
q8_zero_intersection = (q8_height+e7_local_correction-4)/2
assert q8_zero_intersection == 6

q8_point = point1-2*point2
q8_x = q8_point[0]
q8_y = q8_point[1]
h6 = q8_x.denominator().sqrt()
assert h6.degree() == 6
assert q8_x.denominator() == h6**2
assert q8_y.denominator() == h6**3
assert q8_x.numerator().degree() == 16
assert q8_y.numerator().degree() == 24
assert q8_y**2 == q8_x**3+a4*q8_x+a6

print(
    "CM43Q8SECTION|R=P1-2*P2|height=29/2|E7_label=nonidentity"
    "|R.O=6|x=N16/h6^2|y=N24/h6^3",
    flush=True,
)
print(f"CM43Q8SECTION|h6={h6}", flush=True)
print("CM43Q8SECTION|status=PASS", flush=True)
