#!/usr/bin/env sage-python
"""Verify the recovered eighteenth section on Elkies's first conic cover."""

from hashlib import sha256
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
MODEL_PATH = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS_PATH = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
COVER_PATH = ROOT / "elkies-k3/data/fibrations/elkies_2026_rank18_first_cover.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-rank18-first-cover.json"


def load(path):
    payload = path.read_bytes()
    return json.loads(payload), sha256(payload).hexdigest()


model, model_sha = load(MODEL_PATH)
section_data, sections_sha = load(SECTIONS_PATH)
cover, cover_sha = load(COVER_PATH)
assert model["status"] == "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL"
assert section_data["status"] == "PASS_TRANSCRIBED_PUBLISHED_R17_SECTIONS_AND_CHORDS"
assert cover["status"] == "PASS_EXACT_ELKIES_2026_FIRST_RANK18_COVER_SECTION"

R = PolynomialRing(QQ, "t")
t = R.gen()
K = R.fraction_field()
A = R([QQ(value) for value in model["A_coefficients_low_to_high"]])
B = R([QQ(value) for value in model["B_coefficients_low_to_high"]])
q = R([QQ(value) for value in cover["cover"]["q_coefficients_low_to_high"]])
assert q == 4225 * t**2 + 38636 * t + 289444

sections = []
for index, record in enumerate(section_data["sections"]):
    x_coordinate = R([QQ(value) for value in record["x_coefficients_low_to_high"]])
    if index == 0:
        y_coordinate = R([QQ(value) for value in record["y_coefficients_low_to_high"]])
    else:
        chord = record["chord"]
        reference_x, reference_y = sections[int(chord["reference_basis_index"])]
        slope = R([QQ(value) for value in chord["slope_coefficients_low_to_high"]])
        y_coordinate = reference_y + slope * (x_coordinate - reference_x)
    assert y_coordinate**2 == x_coordinate**3 + A * x_coordinate + B
    sections.append((x_coordinate, y_coordinate))

E = EllipticCurve(K, [A, B])
basis = [E(K(x_coordinate), K(y_coordinate)) for x_coordinate, y_coordinate in sections]
trace_vector = tuple(cover["trace_recovery"]["published_basis_vector"])
tau = sum((coefficient * point for coefficient, point in zip(trace_vector, basis)), E(0))
X, Y = tau[0], tau[1]
h = R([QQ(value) for value in cover["trace_recovery"]["h_coefficients_low_to_high"]])
M = R(
    [
        QQ(value)
        for value in cover["trace_recovery"][
            "line_slope_numerator_M_coefficients_low_to_high"
        ]
    ]
)
assert X.denominator() == h**2
assert Y.denominator() == h**3
# In a rootless elliptic K3, h(tau)=4+2(s_tau.s_0).  The denominator h^2
# has degree six, hence s_tau.s_0=3 and h(tau)=10.
assert h.degree() == 3 and cover["trace_recovery"]["height"] == 10

# Direct rootless-bisection compiler: the integral residual chord is the
# unique degree-<6 lift of the triple branch M=-Ny/Nx modulo h^2.
Nx = R(X * h**2)
Ny = R(Y * h**3)
assert Nx.gcd(h) == 1
assert M == R((-Ny * Nx.inverse_mod(h**2)) % (h**2))
chord_discriminant_numerator = (
    M**4 - 6 * M**2 * Nx - 8 * M * Ny - 3 * Nx**2 - 4 * A * h**4
)
assert chord_discriminant_numerator % (h**6) == 0

m = K(M / h)
n = -Y - m * X
sum_x = m**2 - X
product_x = (n**2 - B) / X
discriminant_x = sum_x**2 - 4 * product_x
assert discriminant_x == (QQ(127170526080) * h) ** 2 * q

section = cover["eighteenth_section"]
x0 = R([QQ(value) for value in section["x0_coefficients_low_to_high"]])
x1 = R([QQ(value) for value in section["x1_coefficients_low_to_high"]])
y0 = R([QQ(value) for value in section["y0_coefficients_low_to_high"]])
y1 = R([QQ(value) for value in section["y1_coefficients_low_to_high"]])
assert 2 * x0 == sum_x
assert 2 * x1 == QQ(127170526080) * h
assert y0 == m * x0 + n
assert y1 == m * x1

# Coefficients of 1 and u in the Weierstrass identity modulo u^2=q.
assert y0**2 + y1**2 * q == x0**3 + 3 * x0 * x1**2 * q + A * x0 + B
assert 2 * y0 * y1 == 3 * x0**2 * x1 + x1**3 * q + A * x1

U = PolynomialRing(K, "z")
z = U.gen()
L = K.extension(z**2 - q, "u")
u = L.gen()
EL = E.change_ring(L)
P18 = EL(L(x0) + L(x1) * u, L(y0) + L(y1) * u)
P18_conjugate = EL(L(x0) - L(x1) * u, L(y0) - L(y1) * u)
assert P18 + P18_conjugate == EL(tau)
assert P18 != P18_conjugate

# Parameterization from the rational point at infinity u/t=65.
S = PolynomialRing(QQ, "r")
r = S.gen()
t_of_r = (QQ(289444) - r**2) / (QQ(130) * r - QQ(38636))
u_of_r = QQ(65) * t_of_r + r
assert u_of_r**2 == q(t_of_r)

payload = {
    "schema": "elkies-k3.elkies-2026-rank18-first-cover-certificate.v1",
    "status": "PASS_EXACT_ELKIES_2026_FIRST_RANK18_COVER_SECTION",
    "inputs": {
        str(MODEL_PATH): model_sha,
        str(SECTIONS_PATH): sections_sha,
        str(COVER_PATH): cover_sha,
    },
    "checks": {
        "trace_height": 10,
        "trace_denominator_pattern": "h^2,h^3",
        "linear_double_pole_chord_recovery": True,
        "chord_discriminant_identity": True,
        "weierstrass_identity_mod_cover": True,
        "galois_trace_identity": True,
        "galois_anti_invariant_nonzero": True,
        "rational_parameterization": True,
    },
    "rank_consequence": {
        "generic_rank_lower_bound_on_cover": 18,
        "mechanism": "Lemma 5 of arXiv:2608.25406v1",
    },
    "proof_boundary": (
        "Exact recovery and verification of the new cover section. The rank-18 "
        "consequence invokes the paper's proved Galois-invariance lemma; this is not "
        "an assertion that every rational specialization has rank at least 18."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "ELKIES2026R18|trace_height=10|parameterized=true|"
    f"status={payload['status']}|output={OUTPUT}"
)
