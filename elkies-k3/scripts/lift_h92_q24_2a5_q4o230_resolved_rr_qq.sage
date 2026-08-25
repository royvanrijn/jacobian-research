#!/usr/bin/env sage -python
"""Compile exact 2A5 --q4/orbit230--> A1+A4+A5 equations over QQ.

status: EXACT
claim: fraction-free quartic, exact globally minimal Jacobian, fibre profile

The exact H0 basis is normalized to ``(Z^2,0),(C,1)``.  This exposes the
constant function and avoids fourth powers of the large canonical BB entries.
The chord radicand is divided by Z^6 in QQ[U][T], after which the binary
quartic invariants give the short Jacobian.  No Groebner basis is used.
"""

import hashlib
import json
import sys
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, lcm


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SURFACE = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
P230 = LOCAL / "q24-2a5-p230-scaled-x-qq.json"
ROUTE = ROOT / "artifacts/generated-results/elkies-k3-h3-a5a5-q4o230-q6o1315-promoted-route-certificate.json"
OUTPUT = LOCAL / "q24-2a5-to-a1a4a5-q4o230-resolved-rr-qq.json"
CHECKPOINT = LOCAL / "q24-2a5-q4o230-jacobian-checkpoint.json"
SUPPORTS = LOCAL / "q24-2a5-q4o230-repeated-supports.json"

started = time.monotonic()


def log(stage, **fields):
    suffix = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"A5A5Q4O230QQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


surface = json.loads(SURFACE.read_text())
p230 = json.loads(P230.read_text())
route = json.loads(ROUTE.read_text())
assert surface["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert p230["status"] == "PASS_EXACT_QQ_P230_SECTION_AND_H0"
assert route["status"] == "PASS_EXACT_PROMOTED_DOUBLE_ZERO_EQUATION_COST_ROUTE_TO_PINNED_R17"
assert route["a11_splice"]["nodes"][1] == "A1+A4+A5/MW7"

RT = PolynomialRing(QQ, "T")
T = RT.gen()
A = RT([QQ(value) for value in surface["child"]["minimal_A_coefficients_low_to_high"]])
B = RT([QQ(value) for value in surface["child"]["minimal_B_coefficients_low_to_high"]])
section = p230["P230"]
X = RT([QQ(value) for value in section["X_coefficients_low_to_high"]])
Y = RT([QQ(value) for value in section["Y_coefficients_low_to_high"]])
Z = RT([QQ(value) for value in section["Z_coefficients_low_to_high"]])
assert Y**2 == X**3 + A * X * Z**4 + B * Z**6

constant_basis = p230["H0"]["constant_function_basis"]
assert len(constant_basis) == 2
AA0 = RT([QQ(value) for value in constant_basis[0]["AA_coefficients_low_to_high"]])
BB0 = QQ(constant_basis[0]["BB"])
AA1 = RT([QQ(value) for value in constant_basis[1]["AA_coefficients_low_to_high"]])
BB1 = QQ(constant_basis[1]["BB"])
assert AA0 == Z**2 and not BB0 and BB1 == 1 and AA1.degree() <= 3
assert not (AA0 * X - BB0 * Y) % Z**2
assert not (AA1 * X - BB1 * Y) % Z**2
log("LOAD", section_degrees="8,12,2", H0=2, basis="constant")

# Nested ring QQ[U][T]: T is the old base/quartic variable and U the new base.
RU = PolynomialRing(QQ, "U")
U = RU.gen()
RUT = PolynomialRing(RU, "T")


def nested(poly):
    return RUT([RU(value) for value in RT(poly).list()])


AA0b, AA1b, Zb, Xb, Yb, Ab = map(nested, (AA0, AA1, Z, X, Y, A))
N = AA1b - RUT(U) * AA0b
Db = RUT(U * BB0 - BB1)
raw = (
    N**4 - 6 * Xb * N**2 * Db**2 - 8 * Yb * N * Db**3
    - 3 * Xb**2 * Db**4 - 4 * Ab * Zb**4 * Db**4
)
quartic, remainder = raw.quo_rem(Zb**6)
assert not remainder and quartic.degree() == 4
assert all(RU(coefficient).degree() <= 4 for coefficient in quartic.list())
log("QUARTIC", old_degree=quartic.degree(), new_degree=max(RU(c).degree() for c in quartic))

e, d, c, b, a = [RU(value) for value in quartic.list()]
I = 12 * a * e - 3 * b * d + c**2
J = (
    72 * a * c * e + 9 * b * c * d - 27 * a * d**2
    - 27 * b**2 * e - 2 * c**3
)
A_child = RU(-27 * I)
B_child = RU(-27 * J)
Delta_child = RU(-16 * (4 * A_child**3 + 27 * B_child**2))
degree_profile = (A_child.degree(), B_child.degree(), Delta_child.degree())
log("JACOBIAN_RAW", A_degree=degree_profile[0], B_degree=degree_profile[1],
    Delta_degree=degree_profile[2])
assert degree_profile[0] <= 8 and degree_profile[1] <= 12 and degree_profile[2] <= 24

# In the infinity chart s=1/U, a section of O(weight) acquires order
# ``weight-degree``.  Reading this from the exact degrees avoids copying the
# million-bit coefficients merely to reverse three short polynomials.
infinity_orders = (
    8 - A_child.degree(),
    12 - B_child.degree(),
    24 - Delta_child.degree(),
)
assert infinity_orders == (0, 0, 2)
log("JACOBIAN", A_degree=degree_profile[0], B_degree=degree_profile[1],
    Delta_degree=degree_profile[2], infinity="I2")

if "--jacobian-checkpoint" in sys.argv:
    CHECKPOINT.write_text(json.dumps({
        "A_coefficients_low_to_high": [str(value) for value in A_child.list()],
        "B_coefficients_low_to_high": [str(value) for value in B_child.list()],
        "Delta_coefficients_low_to_high": [str(value) for value in Delta_child.list()],
        "quartic_coefficients_by_old_degree_each_low_to_high_in_U": [
            [str(value) for value in RU(coefficient).list()]
            for coefficient in quartic.list()
        ],
        "degrees_A_B_Delta": [int(value) for value in degree_profile],
        "infinity_orders_A_B_Delta": [int(value) for value in infinity_orders],
    }, sort_keys=True) + "\n")
    log("CHECKPOINT", output=CHECKPOINT)
    raise SystemExit(0)

# Recovering a QQ polynomial gcd is needlessly costly at these coefficient
# heights.  The two rational repeated supports were instead certified by
# simple-root Hensel lifts of Delta^(4) and Delta^(5), followed by exact QQ
# derivative valuations.  Divide them out exactly here.  Squarefreeness and
# the c4 gate for the remaining degree-11 factor descend from one good prime.
supports_data = json.loads(SUPPORTS.read_text())
assert supports_data["status"] == "PASS_EXACT_Q24_2A5_Q4O230_REPEATED_SUPPORTS"
support5 = QQ(supports_data["supports"]["I5"]["U"])
support6 = QQ(supports_data["supports"]["I6"]["U"])
linear5 = RU(U - support5)
linear6 = RU(U - support6)
repeated_factor = linear5**5 * linear6**6
nodal_factor, remainder = Delta_child.quo_rem(repeated_factor)
assert not remainder and nodal_factor.degree() == 11

Fp = GF(supports_data["prime"])
Rp = PolynomialRing(Fp, "u")
A_mod_p = Rp([Fp(value) for value in A_child.list()])
Delta_mod_p = Rp([Fp(value) for value in Delta_child.list()])
nodal_mod_p = Rp([Fp(value) for value in nodal_factor.list()])
assert A_mod_p.degree() == 8 and Delta_mod_p.degree() == 22
assert nodal_mod_p.degree() == 11 and nodal_mod_p.gcd(nodal_mod_p.derivative()) == 1
assert A_mod_p.gcd(Delta_mod_p) == 1
finite_factors = [(linear5, 5), (linear6, 6), (nodal_factor.monic(), 1)]
profile = sorted((factor.degree(), multiplicity) for factor, multiplicity in finite_factors)
assert profile == [(1, 5), (1, 6), (11, 1)]
assert sum(factor.degree() * multiplicity for factor, multiplicity in finite_factors) == 22
log("FIBRES", profile="I2+I5+I6+11I1", ADE="A1+A4+A5", euler=24)


def height_profile(polynomials):
    values = [value for poly in polynomials for value in poly.list()]
    return {
        "maximum_numerator_bits": int(max(abs(value.numerator()).nbits() for value in values)),
        "maximum_denominator_bits": int(max(value.denominator().nbits() for value in values)),
        "maximum_rational_bits": int(max(
            max(abs(value.numerator()).nbits(), value.denominator().nbits())
            for value in values
        )),
    }


input_paths = (SURFACE, P230, ROUTE, SUPPORTS)
payload = {
    "schema": "elkies-k3.q24-2a5-to-a1a4a5-q4o230-resolved-rr-qq.v1",
    "status": "PASS_EXACT_Q24_2A5_Q4O230_A1A4A5_RESOLVED_RR",
    "software": "SageMath 10.9 (conda-forge pinned repository environment)",
    "edge": {
        "source": "2A5/MW7",
        "q": 4,
        "orbit": 230,
        "target": "A1+A4+A5/MW7",
        "old_fibre_degree": 2,
        "zero_for_return": "old_A11_component_10",
    },
    "resolved_RR": {
        "divisor": "O + P230",
        "h0": 2,
        "basis": constant_basis,
        "exact_collision_congruences": True,
        "large_Groebner_required": False,
    },
    "quartic": {
        "old_base_variable": "T",
        "new_base_variable": "U",
        "degree_in_old_base": 4,
        "coefficients_by_old_degree_each_low_to_high_in_U": [
            [str(value) for value in RU(coefficient).list()]
            for coefficient in quartic.list()
        ],
        "fraction_free_identity": "raw chord radicand = Z^6 * quartic",
        "height_profile": height_profile([RU(value) for value in quartic.list()]),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": [str(value) for value in A_child.list()],
        "minimal_B_coefficients_low_to_high": [str(value) for value in B_child.list()],
        "discriminant_coefficients_low_to_high": [str(value) for value in Delta_child.list()],
        "degrees_A_B_Delta": [int(value) for value in degree_profile],
        "finite_fibres": [
            {
                "multiplicity": multiplicity,
                "support_degree": int(factor.degree()),
                "monic_support_coefficients_low_to_high": [str(value) for value in factor.list()],
            }
            for factor, multiplicity in finite_factors
        ],
        "infinity_kind": "I2",
        "infinity_orders_A_B_Delta": [int(value) for value in infinity_orders],
        "fibre_profile": "I2+I5+I6+11I1",
        "euler_number": 24,
        "root_lattice": "A1+A4+A5",
        "root_rank": 10,
        "root_determinant": 60,
        "MW_rank_if_rho19": 7,
        "semistable_c4_c6_gates": True,
        "height_profile": height_profile((A_child, B_child)),
    },
    "verification": {
        "exact_parent_section_identity": True,
        "exact_H0_collision": True,
        "exact_Z6_quartic_division": True,
        "exact_binary_quartic_invariants": True,
        "exact_squarefree_fibre_profile": True,
        "fibre_method": (
            "exact division by the Hensel-recovered I5/I6 supports; good-prime "
            "squarefreeness and c4 gates for the degree-11 residual factor"
        ),
        "exact_infinity_chart_I2": True,
        "route_root_data_matches": route["a11_splice"]["root_data"][1] == [10, 52, 60],
    },
    "proof_boundary": (
        "This proves the exact q4 quartic, globally minimal short Jacobian, and "
        "A1+A4+A5/MW7 fibre data. Equation-level identification of the selected "
        "old_A11_component_10 zero and the q4 return marking is the next gate."
    ),
    "inputs": {
        "paths": [str(path) for path in input_paths],
        "sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in input_paths},
    },
    "elapsed_seconds": round(time.monotonic() - started, 6),
}
assert payload["verification"]["route_root_data_matches"]
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5A5Q4O230QQ|quartic=4|A=8|B=12|fibres=I2+I5+I6+11I1|"
    "ADE=A1+A4+A5|MW=7|status={}|output={}".format(payload["status"], OUTPUT),
    flush=True,
)
