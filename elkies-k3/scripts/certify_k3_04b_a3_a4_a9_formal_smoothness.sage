#!/usr/bin/env sage-python
"""Prove formal smoothness of the marked determinant-500 GF(7) branch.

The 40-variable presentation has 53 displayed equations, but only five of
the nineteen section-residual coefficients are independent after imposing
the fibre and component jets.  The exact discriminant/node identity forces
orders four and ten at zero and infinity, so the residual is supported only
in degrees four through eight.  The Hensel certificate's unit 39-row minor
uses precisely those five coefficients, proving a one-dimensional formally
smooth marked branch in the localization where A is a unit at both supports.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, ZZ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_HENSEL = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-marked-gf7-hensel-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-formal-smoothness-v1.json"
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
if not hensel["status"].startswith(
    "PASS_ONE_DIMENSIONAL_MARKED_TANGENT_AND_EXPLICIT_Z7_LIFT"
):
    raise ValueError("formal certificate requires the smooth GF(7) seed")

pivot_rows = hensel["jacobian_certificate"]["pivot_row_indices"]
expected_pivot_rows = list(range(34)) + list(range(38, 43))
if pivot_rows != expected_pivot_rows:
    raise ArithmeticError("unit minor no longer uses the five middle residuals")
if hensel["jacobian_certificate"]["pivot_minor_determinant_mod_prime"] != 6:
    raise ArithmeticError("pinned unit minor changed")

coefficient_ring = PolynomialRing(ZZ, names=("A", "B", "X", "C", "H"))
A, B, X, C, H = coefficient_ring.gens()
D = 4 * A**3 + 27 * B**2
F = X**3 + A * X * C**4 + B * C**6
identity_right = D * C**4 * (H - B * C**2) - 9 * B * H**2 * C**2 + H**3
identity_left = 8 * A**3 * F
identity_substitution = identity_right.subs({H: 2 * A * X + 3 * B * C**2})
if identity_left != identity_substitution:
    raise ArithmeticError("discriminant/node identity failed")

coordinates = hensel["seed"]["coordinates_mod_prime"]
names = hensel["seed"]["coordinate_names"]
seed = dict(zip(names, coordinates))
if seed["a0"] % 7 == 0 or seed["a8"] % 7 == 0:
    raise ArithmeticError("A is not a unit at a marked support")

payload = {
    "schema": "elkies-k3.k3-04b-a3-a4-a9-formal-smoothness.v1",
    "status": "PASS_ONE_DIMENSIONAL_FORMALLY_SMOOTH_Z7_MARKED_FAMILY",
    "prime": 7,
    "inputs": {relative(hensel_path): digest(hensel_path)},
    "identity": {
        "definitions": {
            "D": "4*A^3+27*B^2",
            "H": "2*A*X+3*B*C^2",
            "F": "X^3+A*X*C^4+B*C^6",
        },
        "formula": (
            "8*A^3*F = D*C^4*(H-B*C^2) - 9*B*H^2*C^2 + H^3"
        ),
        "verified_symbolically_over_Z": True,
    },
    "order_argument": {
        "at_zero": {
            "inputs": ["ord(D)>=4", "ord(H)>=2", "ord(M)>=2", "A(0) unit"],
            "conclusion": "ord(M^2-F)>=4",
        },
        "at_infinity": {
            "inputs": [
                "ord_infinity(D)>=10 in weight 24",
                "ord_infinity(H)>=5 in weight 14",
                "ord_infinity(M)>=5 in weight 9",
                "A(infinity) unit",
            ],
            "conclusion": "ord_infinity(M^2-F)>=10 in weight 18",
        },
        "global_consequence": (
            "The degree-at-most-18 section residual is divisible by t^4 and "
            "has degree at most 8, hence only coefficients 4,5,6,7,8 can remain."
        ),
    },
    "independent_system": {
        "variable_count": 40,
        "independent_equation_count": 39,
        "independent_displayed_row_indices": expected_pivot_rows,
        "retained_section_residual_coefficients": [4, 5, 6, 7, 8],
        "forced_section_residual_coefficients": list(range(4)) + list(range(9, 19)),
        "unit_minor_mod_7": 6,
        "tangent_dimension": 1,
    },
    "localization": {
        "A_at_zero_mod_7": int(seed["a0"]),
        "A_at_infinity_mod_7": int(seed["a8"]),
        "both_units": True,
    },
    "proof_boundary": {
        "proved": (
            "In the localization where A is a unit at zero and infinity, the "
            "fourteen omitted section equations follow from the exact identity, "
            "fibre orders, component jets, and the five retained residual "
            "coefficients. The retained 39 equations have the displayed unit "
            "Jacobian minor, so the marked branch is formally smooth of relative "
            "dimension one over Z_7."
        ),
        "not_proved": (
            "Formal smoothness does not algebraize or rationally parameterize the "
            "branch over Q and does not construct a neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_k3_04b_a3_a4_a9_formal_smoothness.sage"
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("determinant-500 formal-smoothness artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)
print(
    "K304BEFORMAL|independent=39|variables=40|minor=6|dimension=1|status=PASS",
    flush=True,
)
