#!/usr/bin/env sage -python
"""Rank the actual source-q8 smooth principal-part block modulo a good prime.

The source q8 smooth line-bundle lattice requires every negative h-principal
part in the regular q,X frame to vanish.  This is the finite-field companion
to the exact principal-part compiler: it is a fast rank screen only, and does not
replace the characteristic-zero condition matrix.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, binomial, gcd, matrix


ROOT = Path(__file__).resolve().parents[2]
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
FRAME = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-collision-frame.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43.json"


def coefficient(field, value):
    value = QQ(value)
    denominator = field(ZZ(value.denominator()))
    if not denominator:
        raise ValueError("the chosen prime divides an input coefficient denominator")
    return field(ZZ(value.numerator())) / denominator


def polynomial(ring, field, values):
    return ring([coefficient(field, value) for value in values])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument(
    "--extra-h-power", type=int, default=0,
    help="raise every endpoint denominator exponent by this nonnegative amount",
)
parser.add_argument(
    "--extra-e7-pole", type=int, default=0,
    help="allow this additional individual E7 pole order before resolved-chart cancellation",
)
parser.add_argument("--include-kernel", action="store_true")
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be an odd prime different from 3")
if args.extra_h_power < 0:
    raise ValueError("extra-h-power must be nonnegative")
if args.extra_e7_pole < 0:
    raise ValueError("extra-e7-pole must be nonnegative")

p1 = json.loads(P1.read_text())
ambient = json.loads(AMBIENT.read_text())
frame = json.loads(FRAME.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert ambient["status"] == "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT"
assert frame["status"] == "PASS_EXACT_Q8_SMOOTH_COLLISION_FRAME"
assert ambient["ambient_dimension"] == 54

# The endpoint construction permits a nested, endpoint-compatible enlargement:
# retain each E8 floor e while replacing k by k+extra_h and permitting a
# controlled individual E7 slack i<=4k+d+extra_e7.  The seed is the case
# where both extras vanish.  A nonzero E7 slack is diagnostic only: its
# excess pole must later be killed by an actual resolved E7 quotient block,
# not retained as an endpoint condition.
if args.extra_h_power or args.extra_e7_pole:
    ambient_basis = []
    for family in ambient["families"]:
        k = int(family["h_power"]) + args.extra_h_power
        e = int(family["e8_minimal_u_power"])
        d = int(family["e7_allowed_t_denominator_power"])
        for i in range(e, 4*k+d+args.extra_e7_pole+1):
            ambient_basis.append({
                **family["generic_basis"],
                "u_power": i,
                "h_power": k,
            })
else:
    ambient_basis = ambient["ambient_basis"]

finite = GF(args.prime)
ring = PolynomialRing(finite, "u")
u = ring.gen()
field = ring.fraction_field()
h = polynomial(ring, finite, p1["structured_denominator"]["Z4_coefficients"])
assert h.degree() == 4 and gcd(h, h.derivative()) == 1
# The corrected marked E7 frame has base h powers through six.  In the
# regular q,X frame the worst term is m^9/h^6, with h-order -15.
pole_bound = 15 + args.extra_h_power
modulus = h**pole_bound
x_p = field(polynomial(ring, finite, p1["x_entrance_base"]["numerator_coefficients"])) / field(
    polynomial(ring, finite, p1["x_entrance_base"]["denominator_coefficients"])
)
y_p = field(polynomial(ring, finite, p1["y_entrance_base"]["numerator_coefficients"])) / field(
    polynomial(ring, finite, p1["y_entrance_base"]["denominator_coefficients"])
)
rho = field(h)*y_p/x_p
assert gcd(h, ring(rho.numerator())) == 1 and gcd(h, ring(rho.denominator())) == 1
coordinates = tuple((int(item["X_power"]), int(item["q_power"])) for item in frame["regular_degree_18_frame"])
coordinate_index = {entry: index for index, entry in enumerate(coordinates)}
residue_dimension = modulus.degree()


def residue(value):
    value = field(value)
    numerator, denominator = ring(value.numerator()), ring(value.denominator())
    assert gcd(denominator, modulus) == 1
    return ring((numerator * denominator.inverse_mod(modulus)) % modulus)


def column(entry):
    a, b = int(entry["x_power"]), int(entry["m_power"])
    i, k = int(entry["u_power"]), int(entry["h_power"])
    result = [finite.zero()] * (len(coordinates) * residue_dimension)
    for j in range(b + 1):
        exponent = 2*j-b-k-2*a
        if exponent >= 0:
            continue
        value = finite(binomial(b, j)) * u**i * rho**(b-j) * h**(pole_bound+exponent)
        remainder = residue(value)
        offset = coordinate_index[(a, j)] * residue_dimension
        for degree, value in enumerate(remainder.list()):
            result[offset + degree] += value
    return result


columns = [column(entry) for entry in ambient_basis]
condition = matrix(finite, len(coordinates)*residue_dimension, len(columns), lambda row, col: columns[col][row])
kernel = condition.right_kernel()
payload = {
    "schema": "elkies-k3.h92-q8-smooth-principal-parts-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_SMOOTH_BLOCK_RANK",
    "prime": int(args.prime),
    "extra_h_power": args.extra_h_power,
    "extra_e7_pole": args.extra_e7_pole,
    "ambient_basis": ambient_basis,
    "dimensions": {"rows": condition.nrows(), "columns": condition.ncols(), "rank": int(condition.rank()), "kernel": int(kernel.dimension())},
    "boundary": "This is a good-reduction rank screen for the actual smooth block, not a characteristic-zero pencil certificate.",
}
if args.include_kernel:
    payload["ambient_basis"] = ambient_basis
    payload["kernel_basis_rows"] = [[int(value) for value in row] for row in kernel.basis_matrix().rows()]
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92Q8SMOOTHMOD|prime={}|extra_h={}|extra_e7={}|rows={}|columns={}|rank={}|kernel={}|status=EXPERIMENTAL_MODULAR_SMOOTH_BLOCK_RANK".format(args.prime, args.extra_h_power, args.extra_e7_pole, condition.nrows(), condition.ncols(), condition.rank(), kernel.dimension()), flush=True)
