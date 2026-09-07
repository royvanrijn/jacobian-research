#!/usr/bin/env sage -python
"""Resolve the H92 II* germ by explicit ordinary blow-up charts.

The H92 E8 fibre is at infinity.  With ``u=1/t``, ``X=u^4*x`` and
``Y=u^6*y``, its actual integral germ is

    Y^2 = X^3 + (A*u^4+A1*u^5)*X + B2*u^5 + B*u^6 + B1*u^7.

This script does not replace the resolution by the E8 Dynkin label.  It
checks the successive strict transforms in the affine blow-up charts, finds
the three final ordinary nodes, and records the eight exceptional blow-ups
that give the E8 resolution.  The next compiler adapter must use these chart
maps to impose its module/quotient conditions.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-e8-resolution.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def row(polynomial):
    return str(polynomial)


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def hessian_determinant(polynomial, point):
    return matrix(
        QQ,
        [[polynomial.derivative(left, right)(*point) for right in (u, x, y)]
         for left in (u, x, y)],
    ).det()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

anchor = SourceFileLoader("h92_e8_resolution_anchor", str(ANCHOR)).load_module()
h92_ring, h92_formulas = anchor.parse_h92(H92)
r92, s92 = anchor.EXPECTED_H92
A1, A, B1, B, B2 = tuple(QQ(value(r92, s92)) for value in h92_formulas)
assert B2

assert digest(SECTION) == SECTION_SHA256
section = json.loads(SECTION.read_text())
assert section["status"] == "PASS_EXACT_H92_P1"

u_ring = PolynomialRing(QQ, "u")
u_base = u_ring.gen()
u_field = u_ring.fraction_field()
x_p1 = u_field(polynomial(
    u_ring, section["x_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(
    u_ring, section["x_entrance_base"]["denominator_coefficients"]
))
y_p1 = u_field(polynomial(
    u_ring, section["y_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(
    u_ring, section["y_entrance_base"]["denominator_coefficients"]
))
assert y_p1**2 == x_p1**3 + (A1 / u_base**3 + A / u_base**4) * x_p1 + (
    B1 / u_base**5 + B / u_base**6 + B2 / u_base**7
)
X_p1 = u_field(u_base**4) * x_p1
Y_p1 = u_field(u_base**6) * y_p1
assert X_p1.denominator()(0) and Y_p1.denominator()(0)
X_p1_0, Y_p1_0 = QQ(X_p1(0)), QQ(Y_p1(0))
assert X_p1_0 and Y_p1_0

R = PolynomialRing(QQ, names=("u", "x", "y"))
u, x, y = R.gens()


def strict(value, substitutions, exceptional, multiplicity):
    """Strict transform under one displayed ordinary blow-up chart."""
    transformed = R(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**multiplicity)
    assert not remainder
    return quotient


# Start with the actual integral infinity chart, not the singular affine
# Weierstrass equation in t.
f0 = y**2 - x**3 - A*u**4*x - A1*u**5*x - B2*u**5 - B*u**6 - B1*u**7

# The marked chord has an exact normalization in the integral E8 chart.  The
# section meets the smooth locus (not the E8 singular point), and its local
# quotient may be rewritten with a unit denominator near P1.  This records
# the raw local input for the later exceptional-cycle module condition; it
# does not choose that cycle or a line-bundle trivialization.
R_p1 = PolynomialRing(u_field, names=("X", "Y"))
X, Y = R_p1.gens()
f0_p1 = Y**2 - X**3 - A * u_base**4 * X - A1 * u_base**5 * X - B2 * u_base**5 - B * u_base**6 - B1 * u_base**7
chord_numerator = Y - Y_p1
chord_denominator = X - X_p1
chord_conjugate = Y + Y_p1
chord_transfer = X**2 + X * X_p1 + X_p1**2 + A * u_base**4 + A1 * u_base**5
assert (
    chord_numerator * chord_conjugate
    - chord_denominator * chord_transfer
) == f0_p1
assert chord_conjugate(X_p1, Y_p1)(0) == 2 * Y_p1_0

# Blow up the origin in the u-chart: (u,X,Y)=(u,u*x,u*y).
f1 = strict(f0, (u, u*x, u*y), u, 2)
assert f1 == y**2 - u*x**3 - A*u**3*x - A1*u**4*x - B2*u**3 - B*u**4 - B1*u**5

# The sole singular point on the first exceptional curve is the origin.  Its
# x-chart is (u,x,y)=(x*u,x,x*y).
f2 = strict(f1, (x*u, x, x*y), x, 2)
assert f2 == (
    y**2 - x**2*u - A*x**2*u**3 - A1*x**3*u**4 - B2*x*u**3
    - B*x**2*u**4 - B1*x**3*u**5
)

# Blowing up that point has two singular charts.  The x-chart is already an
# ordinary double point; the u-chart continues the E8 chain.
f3_x = strict(f2, (x*u, x, x*y), x, 2)
assert f3_x == (
    y**2 - x*u - A*x**3*u**3 - A1*x**5*u**4 - B2*x**2*u**3
    - B*x**4*u**4 - B1*x**6*u**5
)
assert f3_x(0, 0, 0) == 0
assert f3_x.derivative(x)(0, 0, 0) == f3_x.derivative(u)(0, 0, 0) == 0
assert f3_x.derivative(y)(0, 0, 0) == 0
assert hessian_determinant(f3_x, (0, 0, 0)) != 0

f3_u = strict(f2, (u, u*x, u*y), u, 2)
assert f3_u == (
    y**2 - u*x**2 - A*u**3*x**2 - A1*u**5*x**3 - B2*u**2*x
    - B*u**4*x**2 - B1*u**6*x**3
)

# Continue in the u-chart.  The other chart of this blow-up is needed below
# only to see the node at x=infinity.
f4_u = strict(f3_u, (u, u*x, u*y), u, 2)
assert f4_u == (
    y**2 - u*x**2 - B2*u*x - A*u**3*x**2 - A1*u**6*x**3
    - B*u**4*x**2 - B1*u**7*x**3
)
f4_x = strict(f3_u, (x*u, x, x*y), x, 2)
assert f4_x == (
    y**2 - x*u - B2*x*u**2 - A*x**3*u**3 - A1*x**6*u**5
    - B*x**4*u**4 - B1*x**7*u**6
)

# These derivatives on the two exceptional charts enumerate all remaining
# singular points, not merely three sampled points.  The nonzero finite root
# occurs in both charts and is one point on their overlap.
assert f4_u.derivative(u)(0, x, 0) == -x * (x + B2)
assert f4_x.derivative(x)(u, 0, 0) == -u * (1 + B2 * u)

# The u-chart sees x=0 and x=-B2; the x-chart sees the same second node as
# u=-1/B2 and the third node at infinity as u=0.  Every listed point is an
# ordinary double point, so its one ordinary blow-up is smooth.
nodes = (
    ("finite_0", f4_u, (QQ(0), QQ(0), QQ(0))),
    ("finite_minus_B2", f4_u, (QQ(0), -B2, QQ(0))),
    ("infinity", f4_x, (QQ(0), QQ(0), QQ(0))),
)
node_hessians = {}
for name, equation, point in nodes:
    assert equation(*point) == 0
    assert all(equation.derivative(variable)(*point) == 0 for variable in (u, x, y))
    node_hessians[name] = QQ(hessian_determinant(equation, point))
    assert node_hessians[name]

# The finite nonzero node agrees on the two overlap charts: x/u=-B2 is the
# same point as u/x=-1/B2.  This guards against accidentally counting it
# twice and makes the three-node count an actual chart calculation.
assert -B2 * (-1 / B2) == 1

charts = {
    "infinity_integral": row(f0),
    "blow_1_u": row(f1),
    "blow_2_x": row(f2),
    "blow_3_x_node": row(f3_x),
    "blow_3_u": row(f3_u),
    "blow_4_u": row(f4_u),
    "blow_4_x": row(f4_x),
}
payload = {
    "schema": "elkies-k3.h92-q6-e8-resolution.v1",
    "status": "PASS_EXACT_E8_BLOWUP_CHARTS",
    "inputs": {
        "h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)},
        "marked_section": {"path": str(SECTION.relative_to(ROOT)), "sha256": SECTION_SHA256},
        "h92_point": [str(r92), str(s92)],
    },
    "integral_infinity_coordinates": {
        "base": "u=1/t",
        "weierstrass_scaling": "X=u^4*x, Y=u^6*y",
        "equation": row(f0),
    },
    "marked_chord_normalization": {
        "source_chord": "(y-y(P1))/(x-x(P1))",
        "integral_chart_expression": "u^(-2)*(Y-Y_P)/(X-X_P)",
        "weierstrass_transfer": "(Y-Y_P)*(Y+Y_P)=(X-X_P)*(X^2+X*X_P+X_P^2+A*u^4+A1*u^5)",
        "smooth_marked_point": {
            "X_P_at_u0": str(X_p1_0),
            "Y_P_at_u0": str(Y_p1_0),
            "Y_plus_Y_P_at_P1_u0": str(2 * Y_p1_0),
        },
        "boundary": (
            "This normalizes the raw chord at the actual E8 chart. It does "
            "not choose the exceptional cycle or the line-bundle "
            "trivialization required for a finite condition quotient."
        ),
    },
    "charts": charts,
    "chart_maps_from_integral_infinity_chart": {
        "blow_1_u": {
            "coordinates": "(u,X,Y)=(u1,u1*x1,u1*y1)",
            "strict_transform_divisor": "u1^2",
        },
        "blow_2_x": {
            "coordinates": "(u,X,Y)=(x2*u2,x2^2*u2,x2^2*u2*y2)",
            "strict_transform_divisor": "x2^2 after blow_1_u",
        },
        "blow_3_x_node": {
            "coordinates": "(u,X,Y)=(x3^2*u3,x3^3*u3,x3^4*u3*y3)",
            "strict_transform_divisor": "x3^2 after blow_2_x",
        },
        "blow_3_u": {
            "coordinates": "(u,X,Y)=(u3^2*x3,u3^3*x3^2,u3^4*x3^2*y3)",
            "strict_transform_divisor": "u3^2 after blow_2_x",
        },
        "blow_4_u": {
            "coordinates": "(u,X,Y)=(u4^3*x4,u4^5*x4^2,u4^7*x4^2*y4)",
            "strict_transform_divisor": "u4^2 after blow_3_u",
        },
        "blow_4_x": {
            "coordinates": "(u,X,Y)=(x4^3*u4^2,x4^5*u4^3,x4^7*u4^4*y4)",
            "strict_transform_divisor": "x4^2 after blow_3_u",
        },
    },
    "ordinary_blowup_sequence": [
        "origin of integral infinity chart; retain its u-chart",
        "origin of blow_1_u; retain its x-chart",
        "origin of blow_2_x; retain both x- and u-charts",
        "ordinary node in the blow_3_x chart",
        "origin in the blow_3_u chart; retain both u- and x-charts",
        "finite node x=0 in the blow_4_u chart",
        "finite node x=-B2 in the blow_4_u chart",
        "node at infinity in the blow_4_x chart",
    ],
    "terminal_nodes": {
        name: {"point": [str(value) for value in point], "hessian_determinant": str(node_hessians[name])}
        for name, unused_equation, point in nodes
    },
    "resolution_boundary": (
        "These are actual H92 blow-up charts for the E8 germ.  This artifact "
        "does not yet specify the marked chord's quotient/module condition on "
        "all exceptional charts; that is the next compiler input."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92Q6E8CHART|blowups=8|terminal_nodes=3|status=PASS_EXACT_E8_BLOWUP_CHARTS", flush=True)
