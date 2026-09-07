#!/usr/bin/env sage-python
"""Prove formal smoothness of the determinant-714 A4+2A6 MW1 branch.

The marked system has 40 variables and 47 displayed equations.  Fibre orders
``I5,I7,I7`` and component depths ``(1,2,1)`` force the degree-at-most-18
section residual to be divisible by ``t^2*(t-1)^4`` and to have degree at
most 16.  It therefore has eleven free quotient coefficients.  The pinned
unit Jacobian minor retains exactly residual coefficients 2 through 12, so
the other eight displayed residual equations are consequences of the exact
node/discriminant identity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_HENSEL = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-marked-gf7-hensel-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-formal-smoothness-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--hensel", type=Path, default=DEFAULT_HENSEL)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
hensel_path = arguments.hensel.resolve()
output_path = arguments.output.resolve()
hensel = json.loads(hensel_path.read_text())
if hensel["status"] != (
    "PASS_ONE_DIMENSIONAL_MARKED_TANGENT_AND_EXPLICIT_Z7_LIFT_TO_REQUESTED_PRECISION"
):
    raise ValueError("formal certificate requires the smooth GF(7) seed")
if hensel["system"]["equation_block_sizes"] != {
    "normalization": 1,
    "fibre_at_zero": 5,
    "fibre_at_one": 7,
    "fibre_at_infinity": 7,
    "component_marking": 8,
    "pole_one_section": 19,
}:
    raise ValueError("marked system block structure changed")

expected_pivot_rows = list(range(28)) + list(range(30, 41))
pivot_rows = hensel["jacobian_certificate"]["pivot_row_indices"]
minor = int(hensel["jacobian_certificate"]["pivot_minor_determinant_mod_prime"])
if pivot_rows != expected_pivot_rows or minor % 7 != 2:
    raise ArithmeticError("unit minor no longer retains residual coefficients 2..12")

coefficient_ring = PolynomialRing(ZZ, names=("A", "B", "X", "C", "H"))
A, B, X, C, H = coefficient_ring.gens()
D = 4 * A**3 + 27 * B**2
F = X**3 + A * X * C**4 + B * C**6
identity_right = D * C**4 * (H - B * C**2) - 9 * B * H**2 * C**2 + H**3
identity_left = 8 * A**3 * F
if identity_left != identity_right.subs({H: 2 * A * X + 3 * B * C**2}):
    raise ArithmeticError("discriminant/node identity failed")

field = GF(7)
coordinates = hensel["seed"]["coordinates_mod_prime"]
names = hensel["seed"]["coordinate_names"]
seed = dict(zip(names, coordinates))
a = [field(seed[f"a{index}"]) for index in range(9)]
support_a_values = {
    "zero": int(a[0]),
    "one": int(sum(a, field.zero())),
    "infinity_scaled": int(a[8]),
}
if not all(support_a_values.values()):
    raise ArithmeticError("A is not a unit at every marked support")

payload = {
    "schema": "elkies-k3.k3-cf7f-a4-2a6-mw1-formal-smoothness.v1",
    "status": "PASS_ONE_DIMENSIONAL_FORMALLY_SMOOTH_Z7_MARKED_FAMILY",
    "prime": 7,
    "inputs": {relative(hensel_path): digest(hensel_path)},
    "identity": {
        "definitions": {
            "D": "4*A^3+27*B^2",
            "H": "2*A*X+3*B*C^2",
            "F": "X^3+A*X*C^4+B*C^6",
        },
        "formula": "8*A^3*F = D*C^4*(H-B*C^2) - 9*B*H^2*C^2 + H^3",
        "verified_symbolically_over_Z": True,
    },
    "order_argument": {
        "at_zero": {
            "inputs": ["ord(D)>=5", "ord(H)>=1", "ord(M)>=1", "A(0) unit"],
            "conclusion": "ord(M^2-F)>=2",
        },
        "at_one": {
            "inputs": ["ord(D)>=7", "ord(H)>=2", "ord(M)>=2", "A(1) unit"],
            "conclusion": "ord(M^2-F)>=4",
        },
        "at_infinity": {
            "inputs": [
                "ord_infinity(D)>=7 in weight 24",
                "ord_infinity(H)>=1 in weight 14",
                "ord_infinity(M)>=1 in weight 9",
                "A(infinity) unit",
            ],
            "conclusion": "ord_infinity(M^2-F)>=2 in weight 18",
        },
        "global_consequence": (
            "The degree-at-most-18 residual is t^2*(t-1)^4 times a polynomial "
            "of degree at most ten. Its coefficients 2 through 12 kill the eleven "
            "quotient coefficients triangularly."
        ),
    },
    "independent_system": {
        "variable_count": 40,
        "displayed_equation_count": 47,
        "independent_equation_count": 39,
        "independent_displayed_row_indices": expected_pivot_rows,
        "retained_section_residual_coefficients": list(range(2, 13)),
        "forced_section_residual_coefficients": [0, 1] + list(range(13, 19)),
        "unit_minor_mod_7": minor,
        "tangent_dimension": 1,
    },
    "localization": {
        "A_at_supports_mod_7": support_a_values,
        "all_units": True,
    },
    "proof_boundary": {
        "proved": (
            "In the localization where A is a unit at all three supports, the eight "
            "omitted section equations follow from the exact identity, fibre orders, "
            "component jets, and residual coefficients 2 through 12. The retained "
            "39 equations have a unit Jacobian minor, so the marked branch is formally "
            "smooth of relative dimension one over Z_7."
        ),
        "not_proved": (
            "Formal smoothness does not algebraize or rationally parameterize the "
            "branch over Q, prove the primitive determinant-714 specialization, or "
            "construct a neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_k3_cf7f_a4_2a6_mw1_formal_smoothness.sage"
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != serialized:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)
print("K3CF7FFORMAL|independent=39|variables=40|minor=2|dimension=1|status=PASS")
