#!/usr/bin/env sage -python
"""Export and verify the exact level-474 H3 E7+E8/MW2 source family.

The characteristic-zero normalization certificate gives rational functions
``t(x)``, ``a(x)`` and a multiplier ``m(x)`` such that the published curve

    Y^2 = -27*X^6 + 198*X^4 - 171*X^2 + 576

maps to the H21/H92 component.  If ``Y0=Y/m(x)``, then the invariant
coordinates used in that certificate recover the H92 chart by

    r=(a+Y0)/2,  s=2/(Y0-a).

This checker performs that composition in the quadratic function field,
proves the degree-21 H21/H92 component equation, composes the five pinned H92
short-Weierstrass coefficients, and replays the published non-CM
specialization.  It exports the construction formulas without expanding the
large composed coefficients.

The result is the defining E7+E8/MW2 *source* family.  It does not claim that
the downstream rootless MW17 equation or its specialization to curve 273 has
yet been constructed.
"""

from sage.all import *

import argparse
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
NORMALIZATION = (
    ROOT
    / "artifacts/generated-results/elkies-k3-h21-h92-level474-normalization.json"
)
FACTOR = (
    ROOT / "artifacts/generated-results/elkies-k3-h21-h92-level474-factor-qq.json"
)
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
ANCHOR = (
    ROOT / "artifacts/generated-results/elkies-k3-h3-noncm-q6-source-anchor.json"
)
H21_DESCENT = (
    ROOT / "artifacts/generated-results/elkies-k3-h21-q6-section-descent.json"
)
H92_DESCENT = (
    ROOT / "artifacts/generated-results/elkies-k3-h92-section-descent.json"
)
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-h3-level474-source-family.json"
)

EXPECTED_HASHES = {
    NORMALIZATION: "a2b40633f8cbcbaa94292f540b938f47e080029471f81038f42883bf476a6304",
    FACTOR: "0e9c18c256531c7cc3268de5e7fcdb5e7ef8c422584845cf49ac1ebd46ccfa0e",
    H92: "427559ecd4c2c19d0a4ed7df1019c8a351ed34f454691e9ef1080a8834e74ea1",
    ANCHOR: "0560b1921c87ad2d8db6c293ce070cb30aa75626315801e3e4a71cad59573ea5",
    H21_DESCENT: "9ccdfc7b7a1ca79d549c161e9922051e9f90d7b89ddc81057e17188eedc2a4d2",
    H92_DESCENT: "fe525f75fa87c31afb34755fe63fc778349d2843010eb5c9b17ce6d8b8712e40",
    FRAME: "ba09ec834a7229e11e4ca687d187f663b6368c3e2fac9b5133bb1570e7031599",
}
EXPECTED_POINT = (QQ(13) / 7, QQ(12048) / 343)
EXPECTED_H92_POINT = (QQ(-3621005) / 690947, QQ(158286) / 143585)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repository_path(path):
    return str(path.resolve().relative_to(ROOT))


def parse_rational_function(expression, field):
    return field(expression.replace("^", "**"))


def evaluate_rational_function(value, argument, target_field):
    numerator = sum(
        target_field(coefficient) * argument**degree
        for degree, coefficient in enumerate(value.numerator().list())
    )
    denominator = sum(
        target_field(coefficient) * argument**degree
        for degree, coefficient in enumerate(value.denominator().list())
    )
    return numerator / denominator


def parse_h92_coefficients(path):
    ring = PolynomialRing(QQ, names=("r", "s"))
    r, s = ring.gens()
    environment = {"r": r, "s": s}
    expressions = []
    values = []
    for name in ("A1", "A", "B1", "B", "B2"):
        match = re.search(rf"\b{name}\s*=\s*(.*?);", path.read_text(), re.S)
        assert match is not None
        expression = re.sub(r"\s+", " ", match.group(1)).strip()
        expressions.append(expression)
        values.append(
            ring(sage_eval(expression.replace("^", "**"), locals=environment))
        )
    return ring, expressions, values


def load_frame(path):
    return matrix(
        ZZ,
        [
            [ZZ(entry) for entry in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ],
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--normalization", type=Path, default=NORMALIZATION)
parser.add_argument("--factor", type=Path, default=FACTOR)
parser.add_argument("--h92", type=Path, default=H92)
parser.add_argument("--anchor", type=Path, default=ANCHOR)
parser.add_argument("--h21-descent", type=Path, default=H21_DESCENT)
parser.add_argument("--h92-descent", type=Path, default=H92_DESCENT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

for path, expected_hash in EXPECTED_HASHES.items():
    actual_path = {
        NORMALIZATION: args.normalization,
        FACTOR: args.factor,
        H92: args.h92,
        ANCHOR: args.anchor,
        H21_DESCENT: args.h21_descent,
        H92_DESCENT: args.h92_descent,
        FRAME: FRAME,
    }[path]
    assert digest(actual_path) == expected_hash

normalization = json.loads(args.normalization.read_text())
factor = json.loads(args.factor.read_text())
anchor = json.loads(args.anchor.read_text())
h21_descent = json.loads(args.h21_descent.read_text())
h92_descent = json.loads(args.h92_descent.read_text())
assert normalization["status"] == "PASS_LEVEL474_NORMALIZATION"
assert factor["status"] == "PASS_CHARACTERISTIC_ZERO_FACTOR"
assert anchor["status"] == "PASS_H3_NONCM_Q6_SOURCE_ANCHOR"
assert h21_descent["status"] == "PASS_H21_SECTION_Q_DEFINED"
assert h92_descent["status"] == "PASS_H92_SECTION_Q_DEFINED"
assert anchor["inputs"]["normalization"]["sha256"] == digest(args.normalization)
assert anchor["inputs"]["factor"]["sha256"] == digest(args.factor)
assert anchor["inputs"]["h92"]["sha256"] == digest(args.h92)

# Parse the normalization parameter functions in their internal rational
# coordinate x.
internal_ring = PolynomialRing(QQ, "x")
internal_x = internal_ring.gen()
internal_field = internal_ring.fraction_field()
published_x_map = parse_rational_function(
    normalization["published_level_474"]["x"], internal_field
)
t_internal = parse_rational_function(
    normalization["parameter"]["t"], internal_field
)
a_internal = parse_rational_function(
    normalization["parameter"]["a"], internal_field
)
published_y_formula = normalization["published_level_474"]["y"]
assert published_y_formula.endswith("*Y")
y_multiplier_internal = parse_rational_function(
    published_y_formula[:-2], internal_field
)

# The published x-map is linear fractional.  Invert it exactly.
published_x_numerator = published_x_map.numerator()
published_x_denominator = published_x_map.denominator()
assert published_x_numerator.degree() <= 1
assert published_x_denominator.degree() <= 1
n0, n1 = [published_x_numerator[index] for index in range(2)]
d0, d1 = [published_x_denominator[index] for index in range(2)]

published_base = FunctionField(QQ, "X")
X = published_base.gen()
published_curve_ring = PolynomialRing(published_base, "YY")
YY_polynomial = published_curve_ring.gen()
published_sextic_ring = PolynomialRing(QQ, "X")
X_polynomial = published_sextic_ring.gen()
published_equation = normalization["published_level_474"]["equation"]
assert published_equation.startswith("y^2=")
published_sextic = published_sextic_ring(
    sage_eval(
        published_equation.split("=", 1)[1].replace("^", "**"),
        locals={"x": X_polynomial},
    )
)
assert published_sextic == (
    -27 * X_polynomial**6
    + 198 * X_polynomial**4
    - 171 * X_polynomial**2
    + 576
)
published_curve = published_base.extension(
    YY_polynomial**2 - published_base(published_sextic), names=("Y",)
)
Y = published_curve.gen()
X_curve = published_curve(X)
internal_x_curve = (n0 - X_curve * d0) / (X_curve * d1 - n1)
assert evaluate_rational_function(
    published_x_map, internal_x_curve, published_curve
) == X_curve

t_curve = evaluate_rational_function(t_internal, internal_x_curve, published_curve)
a_curve = evaluate_rational_function(a_internal, internal_x_curve, published_curve)
y_multiplier_curve = evaluate_rational_function(
    y_multiplier_internal, internal_x_curve, published_curve
)
Y0_curve = Y / y_multiplier_curve
r_curve = (a_curve + Y0_curve) / 2
s_curve = 2 / (Y0_curve - a_curve)
assert r_curve / s_curve == t_curve
assert t_curve * (s_curve - 1 / r_curve) == a_curve
assert r_curve + 1 / s_curve == Y0_curve

# Prove that the composed (r,s) lies on the exact degree-21 H21/H92
# component, not only on the quotient double cover.
h92_ring, h92_expressions, h92_coefficients = parse_h92_coefficients(args.h92)
r_symbol, s_symbol = h92_ring.gens()
component = h92_ring.zero()
for term in factor["factor"]["coefficients"]:
    component += (
        QQ(term["coefficient"])
        * r_symbol ** ZZ(term["r"])
        * s_symbol ** ZZ(term["s"])
    )
assert component(r_curve, s_curve) == 0

# Compose the exact short H92 coefficients over the published genus-two
# function field.  Nonvanishing at one exact specialization proves that the
# resulting generic Weierstrass equation is nondegenerate.
composed_h92 = [value(r_curve, s_curve) for value in h92_coefficients]

# Replay the record source point without specializing the quadratic function
# field object, which avoids ambiguous choices of the square-root embedding.
published_x_value, published_y_value = EXPECTED_POINT
internal_roots = (published_x_map - published_x_value).numerator().roots(QQ)
assert len(internal_roots) == 1 and internal_roots[0][1] == 1
internal_x_value = internal_roots[0][0]
t_value = t_internal(internal_x_value)
a_value = a_internal(internal_x_value)
y_multiplier_value = y_multiplier_internal(internal_x_value)
Y0_value = published_y_value / y_multiplier_value
r_value = (a_value + Y0_value) / 2
s_value = 2 / (Y0_value - a_value)
assert (r_value, s_value) == EXPECTED_H92_POINT
assert component(r_value, s_value) == 0
specialized_h92 = [value(r_value, s_value) for value in h92_coefficients]
assert specialized_h92 == [QQ(value) for value in anchor["h92"]["coefficients"]]
A1_value, A_value, B1_value, B_value, B2_value = specialized_h92
assert A1_value and A_value and B1_value and B_value and B2_value

# The short model has III* (E7) at tau=0 and II* (E8) at infinity.  The
# pinned frame supplies the two rational MW directions on the H21/H92
# intersection.
frame = load_frame(FRAME)
assert frame.nrows() == 17 and frame.det() == 948
root_gram = frame[:15, :15]
height_gram = (
    frame[15:, 15:]
    - frame[15:, :15] * root_gram.inverse() * frame[:15, 15:]
)
assert height_gram == matrix(QQ, [[QQ(21) / 2, 3], [3, 46]])

payload = {
    "schema": "elkies-k3.h3-level474-source-family.v1",
    "status": "PASS_EXACT_H3_SOURCE_FAMILY",
    "inputs": {
        "normalization": {
            "path": repository_path(args.normalization),
            "sha256": digest(args.normalization),
        },
        "component": {
            "path": repository_path(args.factor),
            "sha256": digest(args.factor),
        },
        "h92_short_model": {
            "path": repository_path(args.h92),
            "sha256": digest(args.h92),
        },
        "h3_frame": {
            "path": repository_path(FRAME),
            "sha256": digest(FRAME),
        },
        "h21_section_descent": {
            "path": repository_path(args.h21_descent),
            "sha256": digest(args.h21_descent),
        },
        "h92_section_descent": {
            "path": repository_path(args.h92_descent),
            "sha256": digest(args.h92_descent),
        },
    },
    "base_curve": {
        "coordinates": ["X", "Y"],
        "equation": "Y^2=-27*X^6+198*X^4-171*X^2+576",
        "genus": 2,
    },
    "h92_chart_map": {
        "internal_x_from_X": str(internal_x_curve),
        "t_internal": str(t_internal),
        "a_internal": str(a_internal),
        "published_y_multiplier_internal": str(y_multiplier_internal),
        "Y0": "Y/published_y_multiplier_internal(internal_x_from_X)",
        "r": "(a_internal(internal_x_from_X)+Y0)/2",
        "s": "2/(Y0-a_internal(internal_x_from_X))",
        "identities": [
            "r/s=t_internal(internal_x_from_X)",
            "t*(s-1/r)=a",
            "r+1/s=Y0",
            "H21_H92_degree21_component(r,s)=0",
        ],
    },
    "elliptic_k3_source": {
        "base_parameter": "tau",
        "equation": (
            "v^2=u^3+(A1*tau^3+A*tau^4)*u+"
            "(B1*tau^5+B*tau^6+B2*tau^7)"
        ),
        "coefficient_order": ["A1", "A", "B1", "B", "B2"],
        "coefficient_formulas_in_r_s": h92_expressions,
        "fiber_roots": ["E7 at tau=0", "E8 at tau=infinity"],
        "mw_rank": 2,
        "mw_height_gram": [["21/2", "3"], ["3", "46"]],
        "mw_directions_individually_rational": True,
    },
    "record_source_specialization": {
        "published_point": [str(published_x_value), str(published_y_value)],
        "h92_point": [str(r_value), str(s_value)],
        "h92_coefficients": [str(value) for value in specialized_h92],
        "matches_source_anchor": True,
    },
    "proof_boundary": (
        "This is an exact equation and chart map for the genus-two H3 "
        "E7+E8/MW2 source family.  The first q6 equation and the remaining "
        "neighbor chain must still be transported before an explicit "
        "rootless MW17 family or a specialization certificate for curve 273 "
        "is obtained."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "H3SOURCE|base_genus=2|component_degree=21|fibers=E7+E8|"
    "mw_rank=2|height_det=474",
    flush=True,
)
print(
    "H3SOURCE|record_x=13/7|record_y=12048/343|"
    "h92_r=-3621005/690947|h92_s=158286/143585",
    flush=True,
)
print("H3SOURCE|status=PASS_EXACT_H3_SOURCE_FAMILY", flush=True)
