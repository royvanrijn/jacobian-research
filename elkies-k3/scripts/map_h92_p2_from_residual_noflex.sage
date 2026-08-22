#!/usr/bin/env sage -python
"""Exact H92 tangent map using the known non-flex origin directly.

This avoids ``EllipticCurve_from_cubic``'s preliminary flex-scheme search over
QQ(w): the degree-10 residual point explicitly supplies the non-flex case.
"""

from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]
RECOVERY = ROOT / "elkies-k3/scripts/recover_h92_p2_coordinates.sage"
transport_source = RECOVERY.read_text().split("# A finite-field spot check confirms", 1)[0]
globals()["__file__"] = str(RECOVERY)
exec(compile(transport_source, str(RECOVERY), "exec"))
from sage.schemes.curves.constructor import Curve

residual = json.loads((ROOT / "artifacts/generated-results/elkies-k3-h92-p2-residual-crt.json").read_text())
assert residual["complete"]


def polynomial(coefficients, variable):
    return sum(Kw(QQ(value)) * variable**index for index, value in enumerate(coefficients))


target_base = Kw.gen() / (2 * Kw(s))
residual_z = polynomial(residual["residual_z"]["numerator"], target_base) / polynomial(residual["residual_z"]["denominator"], target_base)
residual_x = polynomial(residual["residual_x"]["numerator"], target_base) / polynomial(residual["residual_x"]["denominator"], target_base)

affine = PolynomialRing(Kw, names=("z", "X"))
z, X = affine.gens()
cubic = affine(sum(Kw(coefficient) * Kw.gen()**powers[2] * z**powers[0] * X**powers[1] for powers, coefficient in three_neighbor_cubic.dict().items()))
# The CRT lift is checked at a fresh prime by the residual recovery script.
# Expanding this exact substitution here is substantially more expensive than
# the final exact H92 Weierstrass identity below, which is the certificate
# retained for the recovered section.
print("H92P2NOFLEX|stage=residual_point_loaded", flush=True)


def tangent_residual(point_z, point_x):
    partial_z = cubic.derivative(z)(point_z, point_x)
    partial_x = cubic.derivative(X)(point_z, point_x)
    assert partial_x
    tangent_x = point_x - partial_z / partial_x * (z - point_z)
    intersection = affine(cubic.subs({X: tangent_x}))
    univariate = PolynomialRing(Kw, "z")(intersection)
    assert univariate.degree() == 3
    # The tangent has a double intersection at ``point_z``.  Vieta gives its
    # third intersection without an expensive exact polynomial division.
    qz = -univariate[2] / univariate[3] - 2 * point_z
    qx = tangent_x.subs({z: qz})
    return qz, qx


q_z, q_x = tangent_residual(residual_z, residual_x)
p3_z, p3_x = tangent_residual(q_z, q_x)
print("H92P2NOFLEX|stage=tangent_frame", flush=True)

projective = PolynomialRing(Kw, names=("z", "X", "Z"))
x, y, zed = projective.gens()
form = projective(sum(coefficient * x**powers[0] * y**powers[1] * zed**(3 - powers[0] - powers[1]) for powers, coefficient in cubic.dict().items()))
C = Curve(form)
P = C((residual_z, residual_x, Kw(1)))
P2 = C((q_z, q_x, Kw(1)))
P3 = C((p3_z, p3_x, Kw(1)))

# This is exactly the no-rational-flex branch of Sage's
# EllipticCurve_from_cubic constructor, with P,P2,P3 already known.
M = Matrix(Kw, [list(P), list(P2), list(P3)]).transpose()
F2 = projective(M.act_on_polynomial(form))
xyzM = [M[index, 0] * x + M[index, 1] * y + M[index, 2] * zed for index in range(3)]
T1 = [x*x, y*zed, x*zed]
S1 = x**2 * zed
F3 = projective(F2(T1) // S1)
xyzC = [coordinate(T1) for coordinate in xyzM]
a = Kw(F3.coefficient(x**3))
b = Kw(F3.coefficient(y*y*zed))
ab = a * b
T2 = [-x, y/b, ab*zed]
F4 = projective(F3(T2) / a)
xyzW = [coordinate(T2) for coordinate in xyzC]
S2 = a * S1(T2)
E = EllipticCurve(F4([x, y, 1]))
print("H92P2NOFLEX|stage=weierstrass_transform", flush=True)

# ``xyzW`` is the cubic-to-Weierstrass map.  The formulas involving
# ``M.inverse()`` are its inverse, hence are not applicable to ``q``.
q_image = [coordinate(q_z, q_x, Kw(1)) for coordinate in xyzW]
candidate_on_E = E(q_image)
print("H92P2NOFLEX|stage=tangent_image", flush=True)

_, formulas = anchor.parse_h92(ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt")
A1, A, B1, B, B2 = tuple(QQ(value(r, s)) for value in formulas)
h92 = EllipticCurve([0, 0, 0, A1*target_base**3 + A*target_base**4, B1*target_base**5 + B*target_base**6 + B2*target_base**7])
candidate = E.isomorphism_to(h92)(candidate_on_E)
print("H92P2NOFLEX|stage=h92_isomorphism", flush=True)
candidate_x, candidate_y = candidate.xy()
assert candidate_y**2 == candidate_x**3 + (A1*target_base**3 + A*target_base**4)*candidate_x + (B1*target_base**5 + B*target_base**6 + B2*target_base**7)
print(f"H92P2NOFLEX|status=PASS|x_degrees={candidate_x.numerator().degree()},{candidate_x.denominator().degree()}|y_degrees={candidate_y.numerator().degree()},{candidate_y.denominator().degree()}", flush=True)
print("H92P2NOFLEX|x=" + str(candidate_x), flush=True)
print("H92P2NOFLEX|y=" + str(candidate_y), flush=True)
