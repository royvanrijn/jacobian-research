#!/usr/bin/env sage -python
"""Verify the explicit rank-17 endpoint published in arXiv:2608.25406v1.

Noam Elkies, "An elliptic K3 surface X/Q(t) with Mordell-Weil rank 17, I",
arXiv:2608.25406v1 (26 Aug 2026), now publishes the exact rootless rank-17
Weierstrass model and a height Gram matrix for 17 integral sections.

This script treats the paper as an *independent endpoint oracle* and checks:

* all 17 published sections, reconstructed from the quadratic chord data, lie
  on the published Weierstrass equation;
* the published model is a rootless elliptic K3 (24 simple discriminant zeros);
* the published 17x17 height Gram is positive definite with determinant 948;
* it has 1311 +/- pairs of norm-4 vectors, as stated in the paper;
* it is integrally equivalent to the repository's pinned ``rank17_gram.txt``.

The last check is much stronger than the previous determinant/genus/fingerprint
identification: it gives an explicit integral basis transformation between the
newly published Mordell-Weil basis and our pinned R17 lattice.

The separate coordinate matcher identifies the q12 equation-side endpoint with
this compact model.
"""

from hashlib import sha256
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, QuadraticForm, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
MODEL_DATA = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTION_DATA = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def quadratic_form_from_height_gram(G):
    """Q(v)=<v,v>/2, matching existing rank17 short-vector scripts."""
    coeff = []
    for i in range(G.nrows()):
        for j in range(i, G.ncols()):
            coeff.append(G[i, i] // 2 if i == j else G[i, j])
    return QuadraticForm(ZZ, G.nrows(), coeff)


R = PolynomialRing(QQ, "t")
t = R.gen()

# Theorem 4, arXiv:2608.25406v1, section 2.1.
S = (
    ZZ(307516108335972163537936) * t**8
    + ZZ(10476571005172375234427296) * t**7
    + ZZ(234256046667228607566274912) * t**6
    + ZZ(2020678721371875903158954848) * t**5
    + ZZ(8789387383568632081365832240) * t**4
    + ZZ(21430123310022469548285709072) * t**3
    + ZZ(28607402618712438778345257832) * t**2
    + ZZ(17860826619093915900857289304) * t
    + ZZ(4201305425690184127251888481)
)
T = (
    ZZ(1050290276365892761266194577222156800) * t**12
    + ZZ(67802587761728815952013525763236564480) * t**11
    + ZZ(2392486076703808362288120169049836903680) * t**10
    + ZZ(38126035250980128714796491999580538771200) * t**9
    + ZZ(372202978476351718721663756748866085220800) * t**8
    + ZZ(2373760737463050257069464720014664373086080) * t**7
    + ZZ(9904246958414858348647761354992989326760320) * t**6
    + ZZ(26905633537996991160744810870319331164617600) * t**5
    + ZZ(47243082583908684509409509915652973906060800) * t**4
    + ZZ(52862444598312784274784438443066814490530880) * t**3
    + ZZ(36435013603665838306995466090052055171475872) * t**2
    + ZZ(13865015501478235534002649882546248548532768) * t
    + ZZ(2193201312876924214657300134273061462776968)
)
A = -27 * S
B = QQ(27) / 4 * T

x1 = (
    ZZ(419884536396) * t**4
    - ZZ(6900780974412) * t**3
    + ZZ(84146613883956) * t**2
    + ZZ(448019664127620) * t
    + ZZ(304456582100883)
)
y1 = (
    -ZZ(1917605876395727232) * t**6
    - ZZ(102352278854532258864) * t**5
    - ZZ(1140847719698231045748) * t**4
    - ZZ(4035207954948742785564) * t**3
    - ZZ(3519107150812739581680) * t**2
    + ZZ(3523393851784245137088) * t
    + ZZ(2913630401455186533120)
)
assert y1**2 == x1**3 + A * x1 + B
assert A.degree() == 8 and B.degree() == 12

# Reconstruct all published ordinates from the displayed quadratic chords.
model_bytes = MODEL_DATA.read_bytes()
section_bytes = SECTION_DATA.read_bytes()
model_data = json.loads(model_bytes)
section_data = json.loads(section_bytes)
assert model_data["status"] == "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL"
assert (
    section_data["status"]
    == "PASS_TRANSCRIBED_PUBLISHED_R17_SECTIONS_AND_CHORDS"
)
A_data = R([QQ(value) for value in model_data["A_coefficients_low_to_high"]])
B_data = R([QQ(value) for value in model_data["B_coefficients_low_to_high"]])
assert A_data == A and B_data == B
published_sections = []
for expected_index, record in enumerate(section_data["sections"]):
    assert record["basis_index"] == expected_index
    x_coordinate = R([QQ(value) for value in record["x_coefficients_low_to_high"]])
    if expected_index == 0:
        y_coordinate = R(
            [QQ(value) for value in record["y_coefficients_low_to_high"]]
        )
        assert x_coordinate == x1 and y_coordinate == y1
    else:
        chord = record["chord"]
        reference_index = int(chord["reference_basis_index"])
        assert 0 <= reference_index < expected_index
        reference_x, reference_y = published_sections[reference_index]
        slope = R(
            [QQ(value) for value in chord["slope_coefficients_low_to_high"]]
        )
        assert slope.degree() <= 2
        y_coordinate = reference_y + slope * (x_coordinate - reference_x)
    assert x_coordinate.degree() <= 4 and y_coordinate.degree() <= 6
    assert (1 if y_coordinate.leading_coefficient() > 0 else -1) == record[
        "leading_y_sign"
    ]
    assert y_coordinate**2 == x_coordinate**3 + A * x_coordinate + B
    published_sections.append((x_coordinate, y_coordinate))
assert len(published_sections) == 17

Delta = R(-16 * (4 * A**3 + 27 * B**2))
assert Delta.degree() == 24
assert Delta.gcd(Delta.derivative()).degree() == 0
assert A.gcd(Delta).degree() == 0
assert B.gcd(Delta).degree() == 0

# Published height Gram, Theorem 4. Every diagonal entry is 4 because the
# displayed basis consists of integral sections.
Gpub = matrix(
    ZZ,
    [
        [ 4,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2, 1],
        [-2, 4, 2, 1, 0, 1, 0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1],
        [-2, 2, 4, 2, 1, 0, 1, 2, 1, 1, 1, 0, 1, 1, 2, 1, 0],
        [-2, 1, 2, 4, 1, 0, 1, 1, 2, 1, 1, 1, 0, 1, 0, 1, 0],
        [-2, 0, 1, 1, 4, 2, 0, 1, 1, 0, 1, 1, 1, 2, 1, 0,-1],
        [-2, 1, 0, 0, 2, 4, 1, 0, 0, 1, 2, 1, 1, 1, 1, 0,-1],
        [-2, 0, 1, 1, 0, 1, 4, 0, 0, 1, 1, 0, 1, 1, 2, 2,-2],
        [-2, 1, 2, 1, 1, 0, 0, 4, 1, 2, 1, 0, 1, 1, 1, 1, 0],
        [-2, 1, 1, 2, 1, 0, 0, 1, 4, 2, 2, 1, 2, 2, 1, 1,-1],
        [-2, 1, 1, 1, 0, 1, 1, 2, 2, 4, 2, 0, 2, 1, 2, 1,-1],
        [-2, 1, 1, 1, 1, 2, 1, 1, 2, 2, 4, 0, 1, 2, 2, 1,-2],
        [-2, 2, 0, 1, 1, 1, 0, 0, 1, 0, 0, 4, 0, 1, 0, 1, 1],
        [-2, 1, 1, 0, 1, 1, 1, 1, 2, 2, 1, 0, 4, 2, 2, 1,-2],
        [-2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 1, 2, 4, 2, 1,-2],
        [-2, 1, 2, 0, 1, 1, 2, 1, 1, 2, 2, 0, 2, 2, 4, 1,-2],
        [-2, 1, 1, 1, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0],
        [ 1, 1, 0, 0,-1,-1,-2, 0,-1,-1,-2, 1,-2,-2,-2, 0, 4],
    ],
)
assert Gpub.is_symmetric()
assert Gpub.det() == 948
assert all(value == 4 for value in Gpub.diagonal())
assert Gpub.is_positive_definite()

Qpub = quadratic_form_from_height_gram(Gpub)
# Existing scripts use Q(v)=height(v)/2, so height 4 occurs at Q=2.
short_pub = Qpub.short_vector_list_up_to_length(3, True)
assert len(short_pub[1]) == 0  # no height-2 roots
assert len(short_pub[2]) == 1311

Gpinned = load_matrix(PINNED)
assert Gpinned.nrows() == Gpinned.ncols() == 17
assert Gpinned.det() == 948 and Gpinned.is_positive_definite()
Qpinned = quadratic_form_from_height_gram(Gpinned)

# Sage delegates this positive-definite integral isometry test to PARI qfisom.
M = Qpub.is_globally_equivalent_to(Qpinned, return_matrix=True)
assert M is not False
M = matrix(ZZ, M)
assert abs(M.det()) == 1
# QuadraticForm's basis-change convention has changed in presentation across
# Sage interfaces; accept only a literal Gram identity and record its orientation.
if M.transpose() * Gpub * M == Gpinned:
    orientation = "M^T*Gpub*M=Gpinned"
elif M * Gpub * M.transpose() == Gpinned:
    orientation = "M*Gpub*M^T=Gpinned"
else:
    raise ArithmeticError("qfisom matrix returned without a matching Gram identity")

payload = {
    "schema": "elkies-k3.elkies-2026-published-r17-target.v1",
    "status": "PASS_EXACT_PUBLISHED_R17_IS_PINNED_R17",
    "source": {
        "arxiv": "2608.25406v1",
        "title": "An elliptic K3 surface X/Q(t) with Mordell-Weil rank 17, I",
        "submitted": "2026-08-26",
    },
    "published_equation": {
        "form": "y^2=x^3-27*S(t)*x+(27/4)*T(t)",
        "degrees_A_B_Delta": [int(A.degree()), int(B.degree()), int(Delta.degree())],
        "rootless_semistable": True,
        "published_section_identities": len(published_sections),
        "quadratic_chords_replayed": 16,
    },
    "inputs": {
        str(MODEL_DATA): sha256(model_bytes).hexdigest(),
        str(SECTION_DATA): sha256(section_bytes).hexdigest(),
    },
    "published_height_lattice": {
        "rank": 17,
        "determinant": int(Gpub.det()),
        "height_2_vector_pairs": len(short_pub[1]),
        "height_4_vector_pairs": len(short_pub[2]),
    },
    "pinned_identification": {
        "integrally_isometric": True,
        "basis_change_matrix": [list(map(int, row)) for row in M.rows()],
        "basis_change_determinant": int(M.det()),
        "gram_identity_orientation": orientation,
    },
    "published_high_rank_fibre_parameters": {
        "rank_at_least_25": "-2/377",
        "rank_at_least_26": "-308/251",
        "rank_at_least_27": "2456/135",
        "rank_at_least_28": "-9529/5471",
    },
    "proof_boundary": (
        "This checks every transcribed published section and chord and identifies the "
        "published rootless MW17 height lattice exactly with the repository's pinned "
        "R17 lattice. The separate coordinate matcher identifies the q12 endpoint."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "ELKIES2026R17|det=948|norm4_pairs=1311|isometry_det={}|status={}|output={}".format(
        M.det(), payload["status"], OUTPUT
    ),
    flush=True,
)
