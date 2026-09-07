#!/usr/bin/env sage -python
"""Derive the complete local E8 module for the first H3 marked chord.

Use the actual integral II* chart ``u=1/t, X=u^4*x, Y=u^6*y``.  If

    m=(y-y_P)/(x-x_P)=u^-2 Q,
    Q=(Y-Y_P)/(X-X_P),

then ``Q`` is the integral chord: it has the allowed horizontal poles at
``O`` and ``-P``, is regular at the E8 singular point and on every exceptional
chart, and is nonconstant on the affine II* component.  Thus the local module
for ``O+(-P)`` is ``<1,Q>``.  For the fixed representative
``D=O+(-P)-F_infinity``, the complete E8 local module is

    u * <1,Q>.

Consequently a finite Laurent coefficient combination ``a(u)+b(u)m`` lies in
this module precisely when ``ord_u(a)>=1`` and ``ord_u(b)>=3``.  This is not
the scalar ideal ``(u)``: the latter would incorrectly discard the permitted
horizontal chord pole.  The affine component is included by the integral
chord ``Q``; the resolved exceptional charts verify that its denominator is a
unit at the singular locus.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-e8-p1-branch-module.json"

exec(compile(CORE.read_text(), str(CORE), "exec"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def e8_marked_chord_condition(ambient):
    """Condition block for a declared Laurent ambient of ``a+b*m``.

    An ambient label ``('a',i)`` means ``u^i`` and ``('b',j)`` means
    ``u^j*m=u^(j-2)*Q``.  The obstruction coordinates are the independent
    associated-graded ``1`` and ``Q`` residues below ``u*<1,Q>``.
    """
    labels = []
    for kind, exponent in ambient:
        exponent = int(exponent)
        if kind == "a" and exponent < 1:
            labels.append(("1", exponent))
        elif kind == "b" and exponent < 3:
            labels.append(("Q", exponent - 2))
        elif kind not in ("a", "b"):
            raise ValueError("unknown E8 coefficient kind {}".format(kind))
    quotient_basis = tuple(dict.fromkeys(labels))

    def evaluator(entry):
        kind, exponent = entry
        exponent = int(exponent)
        if kind == "a" and exponent < 1:
            label = ("1", exponent)
        elif kind == "b" and exponent < 3:
            label = ("Q", exponent - 2)
        else:
            return vector(QQ, len(quotient_basis))
        return vector(QQ, [int(label == target) for target in quotient_basis])

    return quotient_condition(
        "H92 E8 marked-chord module",
        ambient,
        evaluator,
        quotient_basis,
        "u*<1,Q> on the complete II* fibre; Q is checked in actual E8 charts",
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

assert digest(SECTION) == SECTION_SHA256
section = json.loads(SECTION.read_text())
assert section["status"] == "PASS_EXACT_H92_P1"
anchor = SourceFileLoader("h92_q6_e8_module_anchor", str(ANCHOR)).load_module()
h92_ring, h92_formulas = anchor.parse_h92(H92)
r92, s92 = anchor.EXPECTED_H92
A1, A, B1, B, B2 = tuple(QQ(value(r92, s92)) for value in h92_formulas)

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()
x_p = u_field(polynomial(
    u_ring, section["x_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(
    u_ring, section["x_entrance_base"]["denominator_coefficients"]
))
y_p = u_field(polynomial(
    u_ring, section["y_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(
    u_ring, section["y_entrance_base"]["denominator_coefficients"]
))
assert y_p**2 == x_p**3 + (A1 / u**3 + A / u**4) * x_p + (
    B1 / u**5 + B / u**6 + B2 / u**7
)
X_p = u_field(u**4) * x_p
Y_p = u_field(u**6) * y_p
assert X_p.denominator()(0) and Y_p.denominator()(0)
X0, Y0 = QQ(X_p(0)), QQ(Y_p(0))
assert X0 and Y0 and Y0**2 == X0**3

# The singular point is X=Y=u=0.  X-X_P has nonzero constant -X0, hence is a
# unit in that local ring and in each chart lying above it.  The displayed
# Weierstrass identity also verifies cancellation at P; the remaining pole is
# precisely the allowed one at -P.
S = PolynomialRing(u_field, names=("X", "Y"))
X, Y = S.gens()
transfer = X**2 + X * X_p + X_p**2 + A * u**4 + A1 * u**5
equation = Y**2 - X**3 - (A * u**4 + A1 * u**5) * X - B2 * u**5 - B * u**6 - B1 * u**7
assert (Y - Y_p) * (Y + Y_p) - (X - X_p) * transfer == equation

# Actual retained chart pullbacks from derive_h92_q6_e8_resolution.sage.  In
# all of them X and Y vanish on the exceptional locus while X_P is a unit, so
# the chord denominator remains a unit.  The three terminal node blow-ups
# preserve this property by pullback.
chart_maps = {
    "blow_1_u": ("u1", "u1*x1", "u1*y1"),
    "blow_2_x": ("x2*u2", "x2^2*u2", "x2^2*u2*y2"),
    "blow_3_x_node": ("x3^2*u3", "x3^3*u3", "x3^4*u3*y3"),
    "blow_3_u": ("u3^2*x3", "u3^3*x3^2", "u3^4*x3^2*y3"),
    "blow_4_u": ("u4^3*x4", "u4^5*x4^2", "u4^7*x4^2*y4"),
    "blow_4_x": ("x4^3*u4^2", "x4^5*u4^3", "x4^7*u4^4*y4"),
}

regression_ambient = (
    ("a", -1), ("a", 0), ("a", 1),
    ("b", 1), ("b", 2), ("b", 3),
)
block = e8_marked_chord_condition(regression_ambient)
assert block["matrix"].rank() == 4
assert block["matrix"].right_kernel().basis_matrix() == matrix(QQ, [
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1],
])

payload = {
    "schema": "elkies-k3.h92-q6-e8-p1-branch-module.v2",
    "status": "PASS_EXACT_E8_MARKED_CHORD_MODULE",
    "inputs": {
        "h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)},
        "marked_section": {"path": str(SECTION.relative_to(ROOT)), "sha256": SECTION_SHA256},
        "resolution_source": "elkies-k3/scripts/derive_h92_q6_e8_resolution.sage",
    },
    "integral_chord": {
        "m": "(y-y_P)/(x-x_P)",
        "Q": "u^2*m=(Y-Y_P)/(X-X_P)",
        "X_P_at_u0": str(X0),
        "Y_P_at_u0": str(Y0),
        "denominator_at_singular_point": str(-X0),
        "cancellation_identity": "(Y-Y_P)(Y+Y_P)=(X-X_P)(X^2+X*X_P+X_P^2+A*u^4+A1*u^5)",
        "chart_maps": {name: {"u_X_Y": list(values)} for name, values in chart_maps.items()},
    },
    "module": {
        "for_O_plus_minus_P": "<1,Q>",
        "for_D_with_F_infinity": "u*<1,Q>",
        "coefficient_rule": "a(u)+b(u)*m belongs iff ord_u(a)>=1 and ord_u(b)>=3",
        "affine_component": "included: Q has its allowed pole at -P on the affine component and is nonconstant there",
    },
    "regression_block": {
        "ambient": [[kind, exponent] for kind, exponent in regression_ambient],
        "quotient_basis": list(block["quotient_basis"]),
        "matrix": [[str(value) for value in row] for row in block["matrix"].rows()],
        "rank": int(block["matrix"].rank()),
        "kernel_dimension": int(block["matrix"].right_kernel().dimension()),
    },
    "boundary": (
        "This completes the local E8 module. The global Laurent ambient, the "
        "common fibre representative for the E7 module, and smooth P.O "
        "collision conditions still need assembly before h0(D)=2 is certified."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6E8MODULE|module=u<1,Q>|condition=a>=u,b>=u3|rows=4|kernel=2|"
    "status=PASS_EXACT_E8_MARKED_CHORD_MODULE",
    flush=True,
)
