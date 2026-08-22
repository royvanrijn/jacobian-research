#!/usr/bin/env sage -python
"""Extract the rational half of the residual construction modulo a prime."""

from pathlib import Path
import json
import os

from sage.all import GF, PolynomialRing, matrix
from importlib.machinery import SourceFileLoader


ROOT = Path(__file__).resolve().parents[2]
input_path = Path(os.environ.get(
    "H92P2_CANDIDATE_INPUT",
    ROOT / "artifacts/generated-results/h92-p2-candidate-mod-100003-500.json",
))
payload = json.loads(input_path.read_text())
field = GF(payload["prime"])
values = payload["values"]
degree = int(os.environ.get("H92P2_DOUBLE_X_DEGREE", "184"))

rows = [
    [-field(value[0])**index for index in range(degree + 1)]
    + [field(value[1]) * field(value[0])**index for index in range(degree + 1)]
    for value in values
]
kernel = matrix(field, rows).right_kernel()
assert kernel.dimension() == 1
relation = kernel.basis()[0]
denominator_raw = list(relation[degree + 1:])
while not denominator_raw[-1]:
    denominator_raw.pop()
denominator_scale = denominator_raw[-1]
numerator = [field(value / denominator_scale) for value in relation[:degree + 1]]
denominator = [field(value / denominator_scale) for value in denominator_raw]

base_ring = PolynomialRing(field, "t")
t = base_ring.gen()
function_field = base_ring.fraction_field()
double_x = (
    sum(function_field(value) * t**index for index, value in enumerate(numerator))
    / sum(function_field(value) * t**index for index, value in enumerate(denominator))
)

# x(2P)=(u^4-2Au^2-8Bu+A^2)/(4(u^3+Au+B)) on y^2=u^3+Au+B.
# The unique linear factor over F_p(t) gives the rational half P2.
anchor = SourceFileLoader(
    "h92_p2_anchor", str(ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage")
).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt")
A1, A, B1, B, B2 = (field(value(r, s)) for value in formulas)
A_t = A1*t**3 + A*t**4
B_t = B1*t**5 + B*t**6 + B2*t**7
ring = PolynomialRing(function_field, "u")
u = ring.gen()
halving = (
    u**4 - 4*double_x*u**3 - 2*A_t*u**2
    - (4*A_t*double_x + 8*B_t)*u + A_t**2 - 4*B_t*double_x
)
linear = [factor.monic() for factor, multiplicity in halving.factor() if factor.degree() == 1]
assert len(linear) == 1
half_x = -linear[0][0]
half_y = (half_x**3 + A_t*half_x + B_t).sqrt()
print(
    "H92P2HALF|double_x_degrees={},{}|half_x_degrees={},{}|half_y_degrees={},{}|status=PASS".format(
        double_x.numerator().degree(), double_x.denominator().degree(),
        half_x.numerator().degree(), half_x.denominator().degree(),
        half_y.numerator().degree(), half_y.denominator().degree(),
    )
)
output = os.environ.get("H92P2_HALF_OUTPUT")
if output:
    Path(output).write_text(json.dumps({
        "prime": int(field.characteristic()),
        "x": {
            "numerator": [int(value) for value in half_x.numerator().list()],
            "denominator": [int(value) for value in half_x.denominator().list()],
        },
        "y": {
            "numerator": [int(value) for value in half_y.numerator().list()],
            "denominator": [int(value) for value in half_y.denominator().list()],
        },
    }, sort_keys=True) + "\n")
