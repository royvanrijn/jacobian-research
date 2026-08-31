#!/usr/bin/env sage -python
"""Match the exact q12/orbit5867 endpoint to Elkies's 2026 model.

The Mobius map was discovered from the two degree-eight A-forms.  This replay
pins that map and verifies the complete A/B identities over QQ, including the
rational Weierstrass scaling and absence of a quadratic twist.  It also maps
the four published rank-25--28 specialization parameters into the old q12
coordinate, quantifying why the compact published coordinate is the preferred
search chart.
"""

import argparse
from hashlib import sha256
import json
from math import gcd
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_Q12 = ROOT / "artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json"
DEFAULT_PUBLISHED = (
    ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
)
DEFAULT_PUBLISHED_CERTIFICATE = (
    ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-h3-q12o5867-elkies-2026-coordinate-match.json"
)


def load_json(path):
    payload = path.read_bytes()
    return json.loads(payload), sha256(payload).hexdigest()


def projective_height(value):
    return max(abs(int(value.numerator())), int(value.denominator()))


def bits(value):
    return [
        int(abs(value.numerator()).nbits()),
        int(value.denominator().nbits()),
    ]


parser = argparse.ArgumentParser()
parser.add_argument("--q12-model", type=Path, default=DEFAULT_Q12)
parser.add_argument("--published-model", type=Path, default=DEFAULT_PUBLISHED)
parser.add_argument(
    "--published-certificate", type=Path, default=DEFAULT_PUBLISHED_CERTIFICATE
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

q12, q12_sha = load_json(args.q12_model)
published, published_sha = load_json(args.published_model)
published_certificate, published_certificate_sha = load_json(
    args.published_certificate
)
assert q12["status"] == "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN"
assert published["status"] == "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL"
assert (
    published_certificate["status"]
    == "PASS_EXACT_PUBLISHED_R17_IS_PINNED_R17"
)

R = PolynomialRing(QQ, "u")
u = R.gen()
A_q12 = R([QQ(value) for value in q12["child"]["minimal_A_coefficients_low_to_high"]])
B_q12 = R([QQ(value) for value in q12["child"]["minimal_B_coefficients_low_to_high"]])
A_published = R([QQ(value) for value in published["A_coefficients_low_to_high"]])
B_published = R([QQ(value) for value in published["B_coefficients_low_to_high"]])
assert [A_q12.degree(), B_q12.degree()] == [8, 12]
assert [A_published.degree(), B_published.degree()] == [8, 12]

# t=(a*u+b)/(c*u+d).  The normalization b=1 makes this record deterministic.
a = QQ(100922443800937835145344826635192511679007247756107475520000) / QQ(
    5626619595700588449039205357834474018054670901587758584778886968441
)
b = QQ(1)
c = -QQ(38200675730597483803366403060419923323594370543452712960000) / QQ(
    5626619595700588449039205357834474018054670901587758584778886968441
)
d = -QQ(12738075412603139735678193811035268846376) / QQ(
    51038889085719625040529308383385780794237
)
assert a * d - b * c != 0

numerator = a * u + b
denominator = c * u + d
A_transformed = R(denominator**8 * A_published(numerator / denominator))
B_transformed = R(denominator**12 * B_published(numerator / denominator))
k4 = A_q12[0] / A_transformed[0]
k6 = B_q12[0] / B_transformed[0]
assert A_q12 == k4 * A_transformed
assert B_q12 == k6 * B_transformed
assert k6**2 == k4**3

# If x_q12=s^2*(c*u+d)^4*x_pub and
# y_q12=s^3*(c*u+d)^6*y_pub, then s^4=k4 and s^6=k6.
s_squared = k6 / k4
assert s_squared**2 == k4
assert s_squared**3 == k6
assert s_squared.is_square()
s = s_squared.sqrt()

specializations = {}
for label, text in published["published_high_rank_fibre_parameters"].items():
    published_t = QQ(text)
    q12_u = (b - d * published_t) / (c * published_t - a)
    specializations[label] = {
        "published_t": str(published_t),
        "published_projective_height": projective_height(published_t),
        "q12_u": str(q12_u),
        "q12_projective_height": projective_height(q12_u),
        "q12_numerator_denominator_bits": bits(q12_u),
    }

payload = {
    "schema": "elkies-k3.h92-q12o5867-elkies-2026-coordinate-match.v1",
    "status": "PASS_EXACT_Q12O5867_IS_ELKIES_2026_PUBLISHED_MODEL",
    "inputs": {
        str(args.q12_model): q12_sha,
        str(args.published_model): published_sha,
        str(args.published_certificate): published_certificate_sha,
    },
    "base_change": {
        "formula": "t=(a*u+b)/(c*u+d)",
        "normalization": "b=1",
        "a": str(a),
        "b": str(b),
        "c": str(c),
        "d": str(d),
        "determinant_nonzero": True,
    },
    "weierstrass_isomorphism": {
        "formula_x": "x_q12=s^2*(c*u+d)^4*x_published",
        "formula_y": "y_q12=s^3*(c*u+d)^6*y_published",
        "s": str(s),
        "A_identity": "A_q12=s^4*(c*u+d)^8*A_published(t)",
        "B_identity": "B_q12=s^6*(c*u+d)^12*B_published(t)",
        "quadratic_twist": "trivial",
    },
    "published_specializations_in_q12_coordinate": specializations,
    "search_consequence": {
        "preferred_coordinate": "elkies_2026_published_t",
        "maximum_published_seed_height": max(
            row["published_projective_height"] for row in specializations.values()
        ),
        "minimum_q12_seed_height_bits": min(
            row["q12_projective_height"].bit_length()
            for row in specializations.values()
        ),
        "reason": (
            "Bounded rational searches in the raw q12 coordinate miss every published "
            "rank-25--28 seed; the compact published coordinate contains all four below "
            "projective height 10000."
        ),
    },
    "proof_boundary": (
        "Exact QQ base-coordinate and short-Weierstrass identification. Published rank "
        "lower bounds are taken from arXiv:2608.25406v1 and are not independently "
        "reproved by this coordinate replay."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "ELKIES2026MATCH|twist=trivial|published_seed_height_max={}|"
    "q12_seed_height_bits_min={}|status={}|output={}".format(
        payload["search_consequence"]["maximum_published_seed_height"],
        payload["search_consequence"]["minimum_q12_seed_height_bits"],
        payload["status"],
        args.output,
    ),
    flush=True,
)
