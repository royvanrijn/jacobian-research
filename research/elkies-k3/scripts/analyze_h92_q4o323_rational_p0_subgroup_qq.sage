#!/usr/bin/env sage-python
"""Compute the exact height Gram of the lifted rational q323 P.O=0 shell.

The self-height of a section is obtained from P.O and its I3+I3+I4 component
profile.  Pairings use h(P-Q), so no resolved section-intersection or
Groebner calculation is needed.
"""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o323-component2-pointing-qq.json"
LIFTS = LOCAL / "q4o323-regular-p0-shell-qq.json"
ANCHOR = LOCAL / "q4o323-p0-shell-anchor-domains-mod61.json"
OUTPUT = LOCAL / "q4o323-rational-p0-subgroup-qq.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


model = json.loads(MODEL.read_text())
lifts = json.loads(LIFTS.read_text())
anchor = json.loads(ANCHOR.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O323_OLD_A11_COMPONENT2_POINTING"
assert lifts["status"] == "EXPERIMENTAL_Q4O323_REGULAR_P0_QQ_LIFTS"
assert anchor["status"] == "PASS_MOD61_Q4O323_REGULAR_P0_SHELL_ANCHOR_DOMAINS"

R = PolynomialRing(QQ, "u")
K = R.fraction_field()
A = R(model["global_short_model"]["A_coefficients_low_to_high"])
B = R(model["global_short_model"]["B_coefficients_low_to_high"])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])

exact = {
    int(record["shell_index"]): record
    for record in lifts["results"]
    if record["status"] == "PASS_EXACT_QQ_SECTION"
}
profiles = {
    int(index): tuple(data["raw_component_profile_I3_I3_I4"])
    for index, data in anchor["shell"]["regular_point_anchor_data"].items()
}


def point(record):
    section = record["section"]
    return E(
        K(R(section["x_coefficients_low_to_high"])),
        K(R(section["y_coefficients_low_to_high"])),
    )


points = {index: point(record) for index, record in exact.items()}
unused = set(points)
pairs = []
while unused:
    left = min(unused)
    inverse = [right for right in unused if points[right] == -points[left]]
    assert len(inverse) == 1
    right = inverse[0]
    assert right != left
    pairs.append((left, right))
    unused.remove(left)
    unused.remove(right)

representatives = [left for left, unused_right in pairs]


def point_dot_zero(P):
    if P.is_zero():
        return -2
    x = K(P[0])
    numerator_degree = R(x.numerator()).degree()
    denominator_degree = R(x.denominator()).degree()
    infinity_excess = max(0, numerator_degree-denominator_degree-4)
    assert denominator_degree % 2 == 0 and infinity_excess % 2 == 0
    return ZZ((denominator_degree+infinity_excess)//2)


def correction(profile):
    return sum(
        QQ(index*(order-index))/order
        for index, order in zip(profile, (3, 3, 4))
    )


def section_height(P, profile):
    return QQ(4+2*point_dot_zero(P))-correction(profile)


def difference_profile(left, right):
    return tuple(
        (profiles[left][slot]-profiles[right][slot]) % order
        for slot, order in enumerate((3, 3, 4))
    )


heights = {
    index: section_height(points[index], profiles[index])
    for index in representatives
}
gram = matrix(QQ, len(representatives), len(representatives))
for i, left in enumerate(representatives):
    gram[i, i] = heights[left]
    for j in range(i):
        right = representatives[j]
        difference = points[left]-points[right]
        difference_height = section_height(
            difference, difference_profile(left, right),
        )
        value = (heights[left]+heights[right]-difference_height)/2
        gram[i, j] = gram[j, i] = value

rank = gram.rank()
independent_columns = list(map(int, gram.pivots()))
independent_gram = gram.matrix_from_rows_and_columns(
    independent_columns, independent_columns,
)
payload = {
    "schema": "elkies-k3.h92-q4o323-rational-p0-subgroup-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O323_RATIONAL_P0_SUBGROUP_HEIGHT_GRAM",
    "exact_signed_section_count": len(points),
    "inverse_pair_count": len(pairs),
    "inverse_pairs": [list(pair) for pair in pairs],
    "representative_shell_indices": representatives,
    "raw_component_profiles_I3_I3_I4": {
        str(index): list(profiles[index]) for index in representatives
    },
    "self_heights": [str(heights[index]) for index in representatives],
    "height_gram": [[str(value) for value in row] for row in gram.rows()],
    "height_gram_rank": int(rank),
    "independent_representative_positions": independent_columns,
    "independent_height_gram_determinant": str(independent_gram.det()),
    "method": {
        "pairing_identity": "<P,Q>=(h(P)+h(Q)-h(P-Q))/2",
        "large_Groebner_required": False,
        "resolved_section_intersections_required": False,
    },
    "proof_boundary": (
        "This proves the exact height Gram and rank of the displayed rational subgroup. "
        "Identifying its integral embedding in the marked q323 MW lattice and testing "
        "membership of q207 are separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MODEL, LIFTS, ANCHOR)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (MODEL, LIFTS, ANCHOR)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323P0GRAM|signed={}|pairs={}|rank={}|det={}|output={}".format(
        len(points), len(pairs), rank, independent_gram.det(), OUTPUT,
    ), flush=True,
)
