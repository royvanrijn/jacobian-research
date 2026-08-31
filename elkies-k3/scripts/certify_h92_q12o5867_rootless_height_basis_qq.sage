#!/usr/bin/env sage -python
"""Certify the exact q12/o5867 rank-17 section height lattice over QQ.

status: ACTIVE_PROOF
claim: 17 exact independent sections with Gram integrally pinned to R17
inputs: exact rootless model, exact lifted sections, modular selection, lattices
outputs: q12o5867-rootless-height-basis-qq.json

The terminal model has no reducible fibres.  For distinct polynomial sections
P,Q, Shioda's formula is <P,Q>=2-P.Q.  Their intersection is computed exactly
as the finite gcd degree of x_P-x_Q and y_P-y_Q plus the minimum normalized
order at infinity.  No Groebner basis or surface elimination is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q12o5867-smooth-rr-qq.json"
SECTIONS = LOCAL / "q12o5867-rootless-selected-basis-qq.json"
SELECTION = LOCAL / "q12o5867-rootless-mod131-selected-basis.json"
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = LOCAL / "q12o5867-rootless-height-basis-qq.json"
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


model = json.loads(MODEL.read_text())
lifted = json.loads(SECTIONS.read_text())
selection = json.loads(SELECTION.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN"
assert model["child"]["root_rank"] == 0
assert len(model["child"]["finite_fibres"]) == 1
assert model["child"]["finite_fibres"][0]["factor_degree"] == 24
assert model["child"]["finite_fibres"][0]["kodaira"] == "I1"
assert model["child"]["infinity"]["kodaira"] == "smooth"
assert lifted["status"] == "PASS_EXACT_QQ_Q12O5867_ROOTLESS_17_SELECTED_SECTIONS"
assert selection["status"] == "PASS_MOD131_Q12O5867_ROOTLESS_REGULAR_SHORT_BASIS"

R = PolynomialRing(QQ, "u")
A = R([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B = R([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])
points = []
for record in lifted["sections"]:
    section = record["section"]
    X = R([QQ(value) for value in section["x_coefficients_low_to_high"]])
    Y = R([QQ(value) for value in section["y_coefficients_low_to_high"]])
    assert X.degree() <= 4 and Y.degree() <= 6
    assert Y**2 == X**3+A*X+B
    points.append((X, Y))
assert len(points) == 17


def order_at_infinity(poly, weight):
    if not poly:
        return 100
    return ZZ(weight-poly.degree())


def section_intersection(left, right):
    dx = left[0]-right[0]
    dy = left[1]-right[1]
    finite = dx.gcd(dy).degree()
    infinity = min(order_at_infinity(dx, 4), order_at_infinity(dy, 6))
    assert finite >= 0 and infinity >= 0
    return ZZ(finite+infinity), ZZ(finite), ZZ(infinity)


gram = matrix(ZZ, 17, 17)
pair_records = []
for left in range(17):
    gram[left, left] = 4
    for right in range(left):
        intersection, finite, infinity = section_intersection(points[left], points[right])
        pairing = ZZ(2-intersection)
        gram[left, right] = gram[right, left] = pairing
        pair_records.append({
            "left": left,
            "right": right,
            "finite_intersection": int(finite),
            "infinity_intersection": int(infinity),
            "height_pairing": int(pairing),
        })

expected = matrix(ZZ, selection["height_gram"])
assert gram == expected
assert gram.is_positive_definite()
assert gram.det() == 948
pinned = load_matrix(PINNED)
transport = matrix(ZZ, selection["short_basis_to_pinned_basis"])
assert transport.det() == -1
assert transport*pinned*transport.transpose() == gram

payload = {
    "schema": "elkies-k3.h92-q12o5867-rootless-height-basis-qq.v1",
    "status": "PASS_EXACT_QQ_Q12O5867_ROOTLESS_RANK17_HEIGHT_BASIS_PINNED",
    "surface": {
        "fibres": "geometrically 24I1",
        "root_rank": 0,
        "torsion": "trivial",
        "torsion_reason": (
            "With no reducible fibres, every nonzero section has height "
            "4+2(P.O)>0, whereas torsion has height zero."
        ),
    },
    "section_count": 17,
    "height_gram": [[int(value) for value in row] for row in gram.rows()],
    "height_gram_determinant": int(gram.det()),
    "pair_intersections": pair_records,
    "basis_to_pinned_rank17": [[int(value) for value in row] for row in transport.rows()],
    "basis_to_pinned_determinant": int(transport.det()),
    "basis_to_pinned_relation": "C * pinned_rank17_gram * C^t = exact QQ section height Gram",
    "rank_conclusion": "unconditional Mordell-Weil rank at least 17 over QQ(u)",
    "saturation_boundary": (
        "The exact section lattice is integrally isometric to pinned R17 and has determinant 948. "
        "Promotion to the full geometric Mordell-Weil group, hence an unconditional saturation "
        "claim, still requires the separate Picard-rank/discriminant upper-bound and source-identity gate."
    ),
    "method": {
        "height_formula": "diagonal 4; off diagonal 2 minus exact section intersection",
        "large_Groebner_required": False,
        "elimination_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MODEL, SECTIONS, SELECTION, PINNED)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (MODEL, SECTIONS, SELECTION, PINNED)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867ROOTLESSHEIGHT|sections=17|det=948|transport_det=-1|torsion=0|"
    "status={}|seconds={:.3f}|output={}".format(
        payload["status"], payload["method"]["runtime_seconds"], OUTPUT,
    ), flush=True,
)
