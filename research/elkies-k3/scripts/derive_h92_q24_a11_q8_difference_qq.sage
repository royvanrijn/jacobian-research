#!/usr/bin/env sage -python
"""Construct D=P12-O_pinned exactly from M and the small residual R.

The inputs M, R=P12-M, and O_pinned are exact marked sections on the same A11
Weierstrass equation.  This script uses fraction-free Jacobian projective
addition, normalizing only by exact polynomial gcds after each addition:

    P12 = M + R,
    D   = P12 - O_pinned.

The final degree-(16,24,6) section is checked literally over QQ and reduced
coefficient-for-coefficient to the pinned modular D.  No Groebner basis or
p-adic lifting is used here.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
BRIDGE = LOCAL / "q24-a11-bridge-m-section-marked-qq.json"
RESIDUAL = LOCAL / "q24-a11-q8-residual-section-qq.json"
ZERO = LOCAL / "q24-a11-pinned-zero-section-qq.json"
MARKED = LOCAL / "q24-a11-q8-horizontal-points-mod100003.json"
OUTPUT = LOCAL / "q24-a11-q8-difference-section-qq.json"

for path in (MODEL, BRIDGE, RESIDUAL, ZERO, MARKED):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

started = time.monotonic()


def log(stage, **fields):
    suffix = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"A11Q8DIFFQQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


model = json.loads(MODEL.read_text())
bridge = json.loads(BRIDGE.read_text())
residual = json.loads(RESIDUAL.read_text())
zero = json.loads(ZERO.read_text())
marked = json.loads(MARKED.read_text())
assert model["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert bridge["status"] == "PASS_EXACT_Q24_A11_BRIDGE_M_SECTION_MARKED_QQ"
assert residual["status"] == "PASS_EXACT_Q24_A11_Q8_RESIDUAL_SECTION_QQ"
assert zero["status"] == "PASS_EXACT_MARKED_A11_PINNED_ZERO_SECTION_QQ"
assert marked["status"] == "PASS_Q24_A11_Q8_HORIZONTAL_POINTS_RECONSTRUCTION_MODP"

RQ = PolynomialRing(QQ, "T")
T = RQ.gen()
A = RQ([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B = RQ([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])


def load_section(record):
    return tuple(
        RQ([QQ(value) for value in record[key]])
        for key in (
            "X_coefficients_low_to_high",
            "Y_coefficients_low_to_high",
            "Z_coefficients_low_to_high",
        )
    )


M = load_section(bridge["section"])
R = load_section(residual["section"])
O = load_section(zero["section"])
log("LOAD", M_degrees=tuple(poly.degree() for poly in M), R_degrees=tuple(poly.degree() for poly in R))


def verify(section, label):
    X, Y, Z = section
    if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
        raise ArithmeticError(f"{label} misses exact A11")


for label, section in (("M", M), ("R", R), ("O_pinned", O)):
    verify(section, label)


def normalize(section, label):
    X, Y, Z = map(RQ, section)
    common = Z.gcd(X).gcd(Y).monic()
    if common.degree() >= 0 and common != 1:
        Xq, Xr = X.quo_rem(common**2)
        Yq, Yr = Y.quo_rem(common**3)
        Zq, Zr = Z.quo_rem(common)
        if Xr or Yr or Zr:
            raise ArithmeticError(f"{label} common projective factor has wrong weights")
        X, Y, Z = Xq, Yq, Zq
    leading = QQ(Z.leading_coefficient())
    X = RQ(X / leading**2)
    Y = RQ(Y / leading**3)
    Z = RQ(Z / leading)
    if Z.leading_coefficient() != 1:
        raise ArithmeticError(f"{label} failed monic-Z normalization")
    verify((X, Y, Z), label)
    log(
        "NORMALIZE",
        label=label,
        removed_degree=int(common.degree()) if common != 1 else 0,
        degrees=tuple(poly.degree() for poly in (X, Y, Z)),
    )
    return X, Y, Z


def add_jacobian(P, Q, label):
    X1, Y1, Z1 = P
    X2, Y2, Z2 = Q
    Z1Z1 = Z1**2
    Z2Z2 = Z2**2
    U1 = X1 * Z2Z2
    U2 = X2 * Z1Z1
    S1 = Y1 * Z2 * Z2Z2
    S2 = Y2 * Z1 * Z1Z1
    H = U2 - U1
    if not H:
        raise ArithmeticError(f"{label} needs tangent/inverse special handling")
    I = (2 * H) ** 2
    J = H * I
    rr = 2 * (S2 - S1)
    V = U1 * I
    X3 = rr**2 - J - 2 * V
    Y3 = rr * (V - X3) - 2 * S1 * J
    Z3 = ((Z1 + Z2) ** 2 - Z1Z1 - Z2Z2) * H
    log("RAW_ADD", label=label, degrees=tuple(poly.degree() for poly in (X3, Y3, Z3)))
    return normalize((X3, Y3, Z3), label)


P12 = add_jacobian(M, R, "P12=M+R")
if tuple(poly.degree() for poly in P12) != (26, 39, 11):
    raise ArithmeticError("exact P12 has unexpected degree profile")
D = add_jacobian(P12, (O[0], -O[1], O[2]), "D=P12-O_pinned")
if tuple(poly.degree() for poly in D) != (16, 24, 6):
    raise ArithmeticError("exact translated difference has unexpected degree profile")

p = ZZ(marked["prime"])
F = GF(p)
RTp = PolynomialRing(F, "T")


def reduce_poly(poly):
    values = []
    for value in poly.list():
        value = QQ(value)
        if value.denominator() % p == 0:
            raise ArithmeticError("bad denominator in marked reduction")
        values.append(F(value.numerator()) / F(value.denominator()))
    return RTp(values)


def marked_section(record):
    return tuple(
        RTp(record[key])
        for key in (
            "X_coefficients_low_to_high",
            "Y_coefficients_low_to_high",
            "Z_coefficients_low_to_high",
        )
    )


if tuple(map(reduce_poly, P12)) != marked_section(marked["q8_target"]["section"]):
    raise ArithmeticError("exact P12 does not reduce to marked target")
if tuple(map(reduce_poly, D)) != marked_section(marked["marked_difference_P12_minus_Opinned"]["section"]):
    raise ArithmeticError("exact D does not reduce to marked difference")

max_bits = max(
    max(abs(ZZ(value.numerator())).nbits(), abs(ZZ(value.denominator())).nbits())
    for poly in D
    for value in poly
)
input_hashes = {
    str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
    for path in (MODEL, BRIDGE, RESIDUAL, ZERO, MARKED)
}
payload = {
    "schema": "elkies-k3.h3-q24-a11-q8-difference-section-qq.v1",
    "status": "PASS_EXACT_Q24_A11_Q8_DIFFERENCE_SECTION_QQ",
    "inputs": {"paths": list(input_hashes), "sha256": input_hashes},
    "construction": {
        "relations": ["P12=M+(P12-M)", "D=P12-O_pinned"],
        "method": "fraction-free Jacobian projective additions with exact weighted gcd normalization",
    },
    "section": {
        "X_coefficients_low_to_high": [str(value) for value in D[0].list()],
        "Y_coefficients_low_to_high": [str(value) for value in D[1].list()],
        "Z_coefficients_low_to_high": [str(value) for value in D[2].list()],
        "degrees_X_Y_Z": [int(poly.degree()) for poly in D],
        "exact_weierstrass_identity": True,
    },
    "marking": {
        "relation": "P12-O_pinned",
        "pinned_A11_MW": [0, 0, -1, 0, 0, 1],
        "P_dot_equation_zero": 6,
        "I12_component_depth_up_to_negation": 6,
        "exact_reduction_matches_marked_seed": True,
        "marked_reduction_prime": int(p),
    },
    "max_rational_coefficient_bits": int(max_bits),
    "large_Groebner_required": False,
    "proof_boundary": (
        "The displayed characteristic-zero section is an exact group-law consequence of exact marked M, "
        "P12-M, and O_pinned, satisfies the A11 equation literally, and matches the pinned modular target. "
        "The exact resolved H0 plane and child equation remain separate gates."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
log("RESULT", max_bits=max_bits, status=payload["status"])
print(f"OUTPUT|{OUTPUT}", flush=True)
