#!/usr/bin/env sage -python
"""Recover P2 over QQ(t) by resolving its tangent base point locally."""

from pathlib import Path
import json

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
RECOVERY = ROOT / "elkies-k3/scripts/recover_h92_p2_coordinates.sage"
transport_source = RECOVERY.read_text().split("# A finite-field spot check confirms", 1)[0]
globals()["__file__"] = str(RECOVERY)
exec(compile(transport_source, str(RECOVERY), "exec"))
payload = json.loads((ROOT / "artifacts/generated-results/elkies-k3-h92-p2-residual-crt.json").read_text())
assert payload["complete"]


def polynomial(coefficients, variable):
    return sum(Kw(QQ(value)) * variable**index for index, value in enumerate(coefficients))


target_base = Kw.gen() / (2 * Kw(s))
rz = polynomial(payload["residual_z"]["numerator"], target_base) / polynomial(payload["residual_z"]["denominator"], target_base)
rX = polynomial(payload["residual_x"]["numerator"], target_base) / polynomial(payload["residual_x"]["denominator"], target_base)
affine = PolynomialRing(Kw, names=("z", "X"))
z, X = affine.gens()
cubic = affine(sum(Kw(c) * Kw.gen()**powers[2] * z**powers[0] * X**powers[1]
                    for powers, c in three_neighbor_cubic.dict().items()))


def tangent(point_z, point_x):
    dz, dX = cubic.derivative(z)(point_z, point_x), cubic.derivative(X)(point_z, point_x)
    line = point_x - dz / dX * (z - point_z)
    restricted = PolynomialRing(Kw, "z")(cubic.subs({X: line}))
    third_z = -restricted[2] / restricted[3] - 2 * point_z
    return third_z, line.subs({z: third_z})


qz, qX = tangent(rz, rX)
p3z, p3X = tangent(qz, qX)
print("H92P2LOCAL|stage=tangent_frame", flush=True)

projective = PolynomialRing(Kw, names=("z", "X", "Z"))
x, y, Z = projective.gens()
form = projective(sum(c * x**powers[0] * y**powers[1] * Z**(3-powers[0]-powers[1])
                      for powers, c in cubic.dict().items()))
M = Matrix(Kw, [[rz, rX, 1], [qz, qX, 1], [p3z, p3X, 1]]).transpose()
F2 = projective(M.act_on_polynomial(form))
xyzM = [M[i, 0]*x + M[i, 1]*y + M[i, 2]*Z for i in range(3)]
T1, S1 = [x*x, y*Z, x*Z], x*x*Z
F3 = projective(F2(T1) // S1)
a, b = Kw(F3.coefficient(x**3)), Kw(F3.coefficient(y*y*Z))
T2 = [-x, y/b, a*b*Z]
F4 = projective(F3(T2) / a)
E = EllipticCurve(F4([x, y, 1]))
MI = M.inverse()
xyzI = [MI[i, 0]*x + MI[i, 1]*y + MI[i, 2]*Z for i in range(3)]
xyzIC = [f([x*Z, x*y, Z*Z]) for f in xyzI]
fwd = [f([-x, b*y, Z/(a*b)]) for f in xyzIC]
print("H92P2LOCAL|stage=weierstrass_chart", flush=True)

# The cubic-to-Weierstrass formulas have a base point at q.  Parameterize
# the smooth branch at q by z=qz+e and solve F(z,X)=0 coefficientwise.
series = PowerSeriesRing(Kw, "e", default_prec=12)
e = series.gen()
z_series = series(qz) + e
Fx = cubic.derivative(X)(qz, qX)
X_series = series(qX) - cubic.derivative(z)(qz, qX) / Fx * e
for degree in range(2, 10):
    error = series(cubic(z_series, X_series))
    X_series -= error[degree] / Fx * e**degree
assert not series(cubic(z_series, X_series)).truncate(10)
images = [series(f(z_series, X_series, 1)) for f in fwd]
valuation = min(image.valuation() for image in images)
limit = [Kw((image / e**valuation)[0]) for image in images]
candidate_on_E = E(limit)
print(f"H92P2LOCAL|stage=base_resolved|valuation={valuation}", flush=True)

_, formulas = anchor.parse_h92(ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt")
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)
h92 = EllipticCurve([0, 0, 0, A1*target_base**3+A*target_base**4,
                      B1*target_base**5+B*target_base**6+B2*target_base**7])
candidate = E.isomorphism_to(h92)(candidate_on_E)
cx, cy = candidate.xy()
assert cy**2 == cx**3 + (A1*target_base**3+A*target_base**4)*cx + (B1*target_base**5+B*target_base**6+B2*target_base**7)
print(f"H92P2LOCAL|status=PASS|x_degrees={cx.numerator().degree()},{cx.denominator().degree()}|y_degrees={cy.numerator().degree()},{cy.denominator().degree()}", flush=True)
print("H92P2LOCAL|x=" + str(cx), flush=True)
print("H92P2LOCAL|y=" + str(cy), flush=True)
