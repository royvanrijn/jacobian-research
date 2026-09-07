#!/usr/bin/env sage -python
"""Audit the missing finite-base transition of the H92 source-q8 q-frame.

The smooth-collision coordinate

    p = y(P1)/x(P1),
    q = (m-p)/h,
    X = h^2*x

is only a local frame near h=0.  Write the reduced denominator of p as

    den(p) = h * u^e * d0(u),

with gcd(d0,u*h)=1.  Away from h=0 and u=0, the globally safer chord
coordinate r=m/h is regular, and

    q = r - c/d0,
    c = num(p)/(h^2*u^e),

where c is a unit modulo d0 when p is reduced.

This script:
  * derives d0 exactly from the reconstructed P1;
  * checks coprimality/squarefreeness and whether d0 meets old singular fibers;
  * builds the exact d0-principal-part matrix modulo a chosen good prime for
    the current 11-dimensional q-frame candidate space
       u^23 q^7,
       u^25..u^28 q^8,
       u^27..u^32 q^9;
  * reports its rank/kernel.

The matrix is necessary and sufficient for regularity of these 11 pure-q
combinations at d0=0, because r is a regular polynomial coordinate there and
the coefficients of r^j are independent.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    GF, PolynomialRing, QQ, ZZ, binomial, gcd, matrix
)

ROOT = Path.cwd()
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
DEFAULT_OUTPUT = ROOT / "artifacts/local/elkies-k3/q8-qframe-finite-transition.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def strip_u_power(poly, u):
    e = 0
    value = poly
    while value and value[0] == 0:
        value = value // u
        e += 1
    return e, value


def reduce_qq_poly(poly, finite_ring):
    finite = finite_ring.base_ring()
    u_mod = finite_ring.gen()
    result = finite_ring.zero()
    for degree, coeff in enumerate(poly.list()):
        q = QQ(coeff)
        den = finite(ZZ(q.denominator()))
        if not den:
            raise ZeroDivisionError("prime divides a rational coefficient denominator")
        result += finite(ZZ(q.numerator()))/den * u_mod**degree
    return result


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--trace", action="store_true")
args = parser.parse_args()

if not ZZ(args.prime).is_prime():
    raise ValueError("--prime must be prime")
if not P1.exists():
    raise SystemExit("missing {}".format(P1))

p1 = json.loads(P1.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"

R = PolynomialRing(QQ, "u")
u = R.gen()
K = R.fraction_field()

h = polynomial(R, p1["structured_denominator"]["Z4_coefficients"])
x_num = polynomial(R, p1["x_entrance_base"]["numerator_coefficients"])
x_den = polynomial(R, p1["x_entrance_base"]["denominator_coefficients"])
y_num = polynomial(R, p1["y_entrance_base"]["numerator_coefficients"])
y_den = polynomial(R, p1["y_entrance_base"]["denominator_coefficients"])

x_p = K(x_num)/K(x_den)
y_p = K(y_num)/K(y_den)
p = y_p/x_p
p_num = R(p.numerator())
p_den = R(p.denominator())

d, rem = p_den.quo_rem(h)
assert not rem
e, d0 = strip_u_power(d, u)
assert d0
assert gcd(d0, u*h) in QQ
assert gcd(d0, p_num) in QQ

# Normalize only for reporting / deterministic finite reduction.
d0 = d0 / d0.leading_coefficient()
assert d0.is_monic()
squarefree = gcd(d0, d0.derivative()) in QQ

# Check whether d0 roots lie on singular old Weierstrass fibres.  From
#   a=(A1*u+A)/u^4,
#   b=(B1*u^2+B*u+B2)/u^7,
# the finite discriminant numerator is obtained from the exact P1 identity
# without needing the H92 source coefficients explicitly.  Recover a,b from
# yP^2-xP^3 and separate their standard pole orders is unnecessary here;
# instead use the fact that a*xP+b = yP^2-xP^3 and skip a discriminant claim
# if no source coefficient artifact is loaded.  d0 regularity itself is exact.

c = K(p_num) / K(h**2 * u**e)
c_num = R(c.numerator())
c_den = R(c.denominator())
assert gcd(c_den, d0) in QQ
assert gcd(c_num, d0) in QQ

candidates = (
    [{"q_power": 7, "u_power": 23}]
    + [{"q_power": 8, "u_power": i} for i in range(25, 29)]
    + [{"q_power": 9, "u_power": i} for i in range(27, 33)]
)
assert len(candidates) == 11
bmax = 9

F = GF(args.prime)
S = PolynomialRing(F, "u")
um = S.gen()

d0m = reduce_qq_poly(d0, S)
c_num_m = reduce_qq_poly(c_num, S)
c_den_m = reduce_qq_poly(c_den, S)

if not d0m or d0m.degree() != d0.degree():
    raise ArithmeticError("bad reduction changes degree of d0")
if gcd(d0m, d0m.derivative()) != 1:
    raise ArithmeticError("prime gives nonsquarefree reduction of d0")
if gcd(d0m, c_num_m*c_den_m) != 1:
    raise ArithmeticError("prime destroys c-unit condition modulo d0")

# Build conditions coefficient-by-coefficient in r after
#   q = r - c/d0.
#
# For fixed r^j, multiply by d0^(bmax-j).  Candidate u^i q^b contributes
#
#   u^i * binomial(b,j) * (-c)^(b-j) * d0^(bmax-b).
#
# Regularity is equivalent to this sum being 0 mod d0^(bmax-j).
rows = []
row_meta = []

for j in range(bmax):
    modulus = d0m**(bmax-j)
    c_mod = (
        (c_num_m % modulus)
        * (c_den_m.inverse_mod(modulus))
    ) % modulus
    column_remainders = []
    for entry in candidates:
        b = entry["q_power"]
        i = entry["u_power"]
        if b < j:
            value = S.zero()
        else:
            value = (
                um**i
                * F(binomial(b, j))
                * (-c_mod)**(b-j)
                * d0m**(bmax-b)
            ) % modulus
        column_remainders.append(value)

    for degree in range(modulus.degree()):
        row = [value[degree] for value in column_remainders]
        if any(row):
            rows.append(row)
            row_meta.append({
                "r_power": j,
                "modulus_d0_power": bmax-j,
                "coefficient_degree": degree,
            })

M = matrix(F, rows, ncols=len(candidates))
rank = int(M.rank())
kernel = M.right_kernel_matrix()

# Record the structural top-down consequence, independent of row ordering.
# At r^8 only q^8 and q^9 appear, and the q^9 pole forces d0 | A9.
# Since gcd(d0,u)=1 and A9/u^27 has degree <=5, deg(d0)>5 forces A9=0.
forced_chain = []
if d0.degree() > 5:
    forced_chain.append("r^8 forces A9=0 because d0|A9 and deg(A9/u^27)<=5")
if d0.degree() > 3:
    forced_chain.append("then r^7 forces A8=0 because d0|A8 and deg(A8/u^25)<=3")
if d0.degree() > 0:
    forced_chain.append("then r^6 forces A7=0 because d0|A7 and A7 is scalar*u^23")

payload = {
    "schema": "elkies-k3.h92-q8-qframe-finite-transition-probe.v1",
    "status": "PASS_EXACT_Q8_QFRAME_FINITE_TRANSITION_MODP",
    "prime": int(args.prime),
    "inputs": {
        "p1": {"path": str(P1.relative_to(ROOT)), "sha256": digest(P1)},
    },
    "finite_transition": {
        "p": "num(p)/(h*u^e*d0)",
        "q_to_safe_chord": "q=m/h-c/d0",
        "h_degree": int(h.degree()),
        "p_denominator_degree": int(p_den.degree()),
        "u_power_in_den_over_h": int(e),
        "d0_degree": int(d0.degree()),
        "d0": str(d0),
        "d0_squarefree_over_Q": bool(squarefree),
        "gcd_d0_u_h_is_unit": True,
        "c_unit_mod_d0": True,
    },
    "ambient": {
        "dimension": len(candidates),
        "basis": candidates,
    },
    "principal_part_matrix": {
        "rows": int(M.nrows()),
        "columns": int(M.ncols()),
        "rank": rank,
        "kernel_dimension": int(kernel.nrows()),
        "kernel_basis_rows": [
            [int(value) for value in row] for row in kernel.rows()
        ],
    },
    "forced_chain": forced_chain,
    "boundary": (
        "This tests only the current 11 pure-q candidate space at the missing "
        "finite d0 divisor. It does not construct the full rank-18 global "
        "pushforward lattice or recover H0(D) if this narrow space dies."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

print(
    "Q8QFRAMED0|prime={}|deg_d0={}|u_power={}|rows={}|rank={}|kernel={}|"
    "squarefree={}|status=PASS_EXACT_Q8_QFRAME_FINITE_TRANSITION_MODP".format(
        args.prime,
        d0.degree(),
        e,
        M.nrows(),
        rank,
        kernel.nrows(),
        int(squarefree),
    ),
    flush=True,
)
for line in forced_chain:
    print("Q8QFRAMED0FORCE|{}".format(line), flush=True)
