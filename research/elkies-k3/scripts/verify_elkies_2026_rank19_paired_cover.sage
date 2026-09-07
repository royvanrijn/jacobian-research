#!/usr/bin/env sage-python
"""Verify both new sections and the E0-to-t maps on the paired cover."""

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, sage_eval


ROOT = Path(__file__).resolve().parents[2]
MODEL_PATH = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS_PATH = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
FIRST_PATH = ROOT / "elkies-k3/data/fibrations/elkies_2026_rank18_first_cover.json"
PAIRED_PATH = ROOT / "elkies-k3/data/fibrations/elkies_2026_rank19_paired_cover.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-rank19-paired-cover.json"


def load(path):
    payload = path.read_bytes()
    return json.loads(payload), sha256(payload).hexdigest()


model, model_sha = load(MODEL_PATH)
section_data, sections_sha = load(SECTIONS_PATH)
first, first_sha = load(FIRST_PATH)
paired, paired_sha = load(PAIRED_PATH)

R = PolynomialRing(QQ, "t")
t = R.gen()
K = R.fraction_field()
A = R([QQ(value) for value in model["A_coefficients_low_to_high"]])
B = R([QQ(value) for value in model["B_coefficients_low_to_high"]])
q1 = R([QQ(value) for value in first["cover"]["q_coefficients_low_to_high"]])
q2 = R([QQ(value) for value in paired["second_cover"]["q_coefficients_low_to_high"]])

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
    sections.append((x_coordinate, y_coordinate))

E = EllipticCurve(K, [A, B])
basis = [E(x_coordinate, y_coordinate) for x_coordinate, y_coordinate in sections]
second = paired["second_cover"]
trace_vector = tuple(second["trace_published_basis_vector"])
tau = sum((coefficient * point for coefficient, point in zip(trace_vector, basis)), E(0))
X, Y = tau[0], tau[1]
h = R([QQ(value) for value in second["h_coefficients_low_to_high"]])
M = R([QQ(value) for value in second["line_slope_numerator_M_coefficients_low_to_high"]])
assert X.denominator() == h**2 and Y.denominator() == h**3
assert h.degree() == 3 and second["trace_height"] == 10
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
assert sum_x**2 - 4 * product_x == (QQ(445172834496) * h) ** 2 * q2

new_section = second["new_section"]
x0 = R([QQ(value) for value in new_section["x0_coefficients_low_to_high"]])
x1 = R([QQ(value) for value in new_section["x1_coefficients_low_to_high"]])
y0 = R([QQ(value) for value in new_section["y0_coefficients_low_to_high"]])
y1 = R([QQ(value) for value in new_section["y1_coefficients_low_to_high"]])
assert 2 * x0 == sum_x and 2 * x1 == QQ(445172834496) * h
assert y0 == m * x0 + n and y1 == m * x1
assert y0**2 + y1**2 * q2 == x0**3 + 3 * x0 * x1**2 * q2 + A * x0 + B
assert 2 * y0 * y1 == 3 * x0**2 * x1 + x1**3 * q2 + A * x1

# Replay the quartic obtained after parameterizing the first conic.
S = PolynomialRing(QQ, "r")
r = S.gen()
denominator = 130 * r - 38636
t_of_r = (QQ(289444) - r**2) / denominator
F = S([QQ(value) for value in paired["quartic_model"]["F_coefficients_low_to_high"]])
assert F == denominator**2 * q2(t_of_r)

# Verify the recovered E0 -> quartic -> paired-cover maps in Q(E0).
P = PolynomialRing(QQ, names=("x", "y"))
x, y = P.gens()
PF = P.fraction_field()
equation = y**2 - (x**3 + QQ(1029367969) * x**2 - QQ(42900734074705920) * x)
mapping = paired["E0"]["map_to_quartic"]
locals_map = {"x": x, "y": y}
D_r = P(sage_eval(mapping["D_r"].replace("^", "**"), locals=locals_map))
N_r = P(sage_eval(mapping["N_r"].replace("^", "**"), locals=locals_map))
N_v = P(sage_eval(mapping["N_v"].replace("^", "**"), locals=locals_map))
r_E0 = PF(N_r) / (6 * PF(D_r))
v_E0 = -QQ(13) * PF(N_v) / (6 * PF(D_r) ** 2)
t_E0 = (QQ(289444) - r_E0**2) / (QQ(130) * r_E0 - QQ(38636))
u1_E0 = QQ(65) * t_E0 + r_E0
u2_E0 = v_E0 / (QQ(130) * r_E0 - QQ(38636))


def zero_mod_e0(value):
    numerator = P(PF(value).numerator())
    return P.ideal(equation).reduce(numerator) == 0


assert zero_mod_e0(v_E0**2 - F(r_E0))
assert zero_mod_e0(u1_E0**2 - q1(t_E0))
assert zero_mod_e0(u2_E0**2 - q2(t_E0))

# Exact rank-at-least-four certificate for the displayed E0 generators.
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(ROOT / "elliptic-curves/cas")]
from ecsearch.q12o5867_specialization import short_certificate_model  # noqa: E402
from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    source_point_to_target,
    verify_finite_quotient_certificate,
)

e0_model = tuple(Fraction(value) for value in paired["E0"]["weierstrass_coefficients_a1_a2_a3_a4_a6"])
e0_points = tuple(
    (Fraction(point[0]), Fraction(point[1])) for point in paired["E0"]["generators"]
)
E0 = EllipticCurve(QQ, list(e0_model))
assert all(E0(QQ(px), QQ(py)) for px, py in e0_points)
short_model, short_change = short_certificate_model(e0_model)
short_points = tuple(source_point_to_target(point, short_change) for point in e0_points)
independence = build_finite_quotient_certificate(
    short_model, short_points, relation_prime=3, prime_bound=500
)
verify_finite_quotient_certificate(short_model, short_points, independence)
assert independence["certified_independent"]

payload = {
    "schema": "elkies-k3.elkies-2026-rank19-paired-cover-certificate.v1",
    "status": "PASS_EXACT_ELKIES_2026_RANK19_PAIRED_COVER",
    "inputs": {
        str(MODEL_PATH): model_sha,
        str(SECTIONS_PATH): sections_sha,
        str(FIRST_PATH): first_sha,
        str(PAIRED_PATH): paired_sha,
    },
    "checks": {
        "first_section_certificate_imported": True,
        "second_trace_height": 10,
        "second_linear_double_pole_chord_recovery": True,
        "second_chord_discriminant": True,
        "second_section_weierstrass_identity": True,
        "quartic_from_first_conic_parameterization": True,
        "E0_to_quartic_identity": True,
        "E0_to_both_conics": True,
        "E0_generator_count": 4,
        "E0_generator_independence": independence,
    },
    "rank_consequence": {
        "E0_rank_lower_bound": 4,
        "generic_paired_cover_rank_lower_bound": 19,
    },
    "proof_boundary": (
        "The four E0 points are unconditionally independent, proving rank at least 4; "
        "the paper's exact rank-4 statement is not independently upper-bounded here. "
        "Rank at least 19 uses the paper's proved two-character Galois argument."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "ELKIES2026R19|E0_rank_lower_bound=4|maps=true|"
    f"status={payload['status']}|output={OUTPUT}"
)
