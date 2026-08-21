#!/usr/bin/env sage
"""Pin the deforming fifth-q4 pair23 equation and generic q6 join.

The default pair23 projection is a genus-three gauge.  Two rational secant
sections have values

  mu_A=(68*T^2+56*T+50)/(T^2+29*T+37),
  mu_B=(34*T^2+55*T+13)/(T^2+29*T+37).

This verifier applies s=(U-mu_A)/(mu_B-mu_A) to the irreducible pair23 branch
factor, restores the factorization unit in the binary-quartic Jacobian, and
certifies fibers I0*+3I4+3I2 over GF(73).  It also pins the exact lattice join:
the degree-47 generic q4 witness specializes to this pair23 torsion class and
its already certified q6 child is rootless/MW17.

This is a finite-field equation certificate plus an exact lattice transport;
it is not yet a characteristic-zero family or a rank-31 specialization.
"""

import hashlib
import json
from pathlib import Path

from sage.all import GF, FunctionField, PolynomialRing, gcd


ROOT = Path(__file__).resolve().parents[2]
PAIR23 = (
    ROOT / "artifacts/generated-results/"
    "q80-fifth-q4-marked-projection-pair23-gf73.json"
)
WINDOW = ROOT / "artifacts/local/q80-q4-deforming-window-10000-20000.json"
VERTICAL = (
    ROOT / "artifacts/local/"
    "q80-deforming-fifth-q4-vertical-compensation.json"
)
Q6 = (
    ROOT / "artifacts/generated-results/"
    "q80-alternate-fifth-q6-rootless-transport.json"
)
KNOWN_HASHES = {
    PAIR23: "6e9952d0a8c4748a499cd881585a615ce8b6538687779b69b536783a9c242f2e",
    WINDOW: "0cb51a1593af0a12acbb743ac7f03177849d0e828c8295ef8ffb647101b1f943",
    VERTICAL: "290350d05603aa76fb05153871439bafce45a58261215fc49808553f233fd345",
    Q6: "48381d91e288b2cefb85b1484d351d659748f801ea57d190453bd2db0a56eaab",
}
payloads = {}
for path, expected in KNOWN_HASHES.items():
    content = path.read_bytes()
    assert hashlib.sha256(content).hexdigest() == expected
    payloads[path] = json.loads(content)

pair23 = payloads[PAIR23]
window = payloads[WINDOW]
vertical = payloads[VERTICAL]
q6 = payloads[Q6]

q4_v = [
    -9, 8, -11, 10, -4, 0, 5, 1, -6,
    6, 1, -2, -1, -1, 1, 2, 0,
]
q6_v = [
    0, -2, 4, 2, -1, 2, 1, -1, 1,
    0, 1, -1, 1, 0, 0, 0, 0,
]
q4_row = next(row for row in window["candidates"] if row["v"] == q4_v)
assert q4_row["generic_child_roots"] == [1, 2, 2]
assert q4_row["generic_child_mw_rank"] == 16
assert q4_row["cm24_child_roots"] == [16, 66, 2048]
assert q4_row["old_fiber_degree"] == 47
vertical_row = next(
    row for row in vertical["rows"] if row["old_fiber_degree"] == 47
)
assert vertical_row["explicit_section_pair_matches"] == [[2, 3]]
assert vertical_row["residual_free_fibers"] == 0
assert q6["q4"]["v"] == q4_v
assert q6["q6"]["v"] == q6_v
assert q6["q6"]["child_root_data"] == [0, 0, 1]
assert q6["q6"]["child_MW_rank"] == 17

finite = GF(73, impl="modn")
source_ring = PolynomialRing(finite, names=("T", "U"))
T, U = source_ring.gens()
integral_cover = sum(
    finite(coefficient)*T**t_degree*U**u_degree
    for t_degree, u_degree, coefficient
    in pair23["integral_double_cover_terms_T_U_coefficient"]
)
branch = next(
    factor for factor, exponent in integral_cover.factor()
    if factor.degree(T) == 8 and exponent == 1
)

gauge_ring = PolynomialRing(finite, names=("s", "tau"))
s_poly, tau_poly = gauge_ring.gens()
denominator = tau_poly**2+29*tau_poly+37
mu_A_numerator = 68*tau_poly**2+56*tau_poly+50
mu_B_numerator = 34*tau_poly**2+55*tau_poly+13
gauge_numerator = (
    mu_A_numerator*(1-s_poly)+mu_B_numerator*s_poly
)
transformed = sum(
    gauge_ring(coefficient)*tau_poly**t_degree
    *gauge_numerator**u_degree*denominator**(4-u_degree)
    for (t_degree, u_degree), coefficient in branch.dict().items()
)
factor_degrees = tuple(
    (int(factor.degree(tau_poly)), int(exponent))
    for factor, exponent in transformed.factor()
)
assert sorted(factor_degrees) == sorted(
    ((1, 2), (1, 2), (1, 4), (1, 4), (4, 1))
), factor_degrees

base = FunctionField(finite, "s")
s = base.gen()
cover_ring = PolynomialRing(base, "tau")
tau = cover_ring.gen()
cover = cover_ring(sum(
    base(coefficient)*tau**exponents[1]*s**exponents[0]
    for exponents, coefficient in transformed.dict().items()
))
factorization = cover.factor()
twist = base(factorization.unit())
odd_part = cover_ring(1)
for factor, exponent in factorization:
    if int(exponent) % 2:
        odd_part *= factor
quartic = odd_part.monic()
assert quartic.degree() == 4
coefficients = list(quartic.list())+[base(0)]*5
e, d, c, b, a = coefficients[:5]
invariant_I = 12*a*e-3*b*d+c**2
invariant_J = 72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
jacobian_A = twist**2*(-27*invariant_I)
jacobian_B = twist**3*(-27*invariant_J)
delta_core = twist**6*(4*invariant_I**3-invariant_J**2)
assert jacobian_A.denominator() == 1
assert jacobian_B.denominator() == 1
assert delta_core.denominator() == 1
A = jacobian_A.numerator()
B = jacobian_B.numerator()
Delta = delta_core.numerator()
finite_scalings = []
for factor, _ in gcd(A, B).factor():
    scale_order = min(A.valuation(factor)//4, B.valuation(factor)//6)
    if scale_order <= 0:
        continue
    A //= factor**(4*scale_order)
    B //= factor**(6*scale_order)
    Delta //= factor**(12*scale_order)
    finite_scalings.append((str(factor), int(scale_order)))
assert not finite_scalings
assert (A.degree(), B.degree(), Delta.degree()) == (4, 9, 18)

delta_factors = tuple(
    (str(factor), int(exponent)) for factor, exponent in Delta.factor()
)
assert sorted(exponent for _, exponent in delta_factors) == [2, 2, 2, 4, 4, 4]
assert (8-A.degree(), 12-B.degree(), 24-Delta.degree()) == (4, 3, 6)
finite_signature = (
    ("I2", 3),
    ("I4", 3),
)
infinity_signature = "I0*"
root_data = (16, 66, 2048)
euler_number = 3*2+3*4+6
assert euler_number == 24

artifact = {
    "schema": "q80-deforming-fifth-pair23-gf73-v1",
    "status": "exact_finite_field_equation_and_lattice_join_certificate",
    "prime": 73,
    "source_artifacts": [
        {
            "path": str(path.relative_to(ROOT)),
            "sha256": digest,
        }
        for path, digest in KNOWN_HASHES.items()
    ],
    "pair23_secant_values": {
        "mu_A_keys": ["L01", "L03"],
        "mu_A_numerator_coefficients_low_to_high": [50, 56, 68],
        "mu_B_keys": ["L12"],
        "mu_B_numerator_coefficients_low_to_high": [13, 55, 34],
        "common_denominator_coefficients_low_to_high": [37, 29, 1],
        "gauge": "s=(U-mu_A)/(mu_B-mu_A)",
        "recognition_samples": 27,
        "withheld_samples": 8,
    },
    "transformed_branch_factor_degrees_exponents": [
        list(row) for row in factor_degrees
    ],
    "quadratic_twist": str(twist),
    "A_coefficients_low_to_high": list(map(int, A.list())),
    "B_coefficients_low_to_high": list(map(int, B.list())),
    "Delta_coefficients_low_to_high": list(map(int, Delta.list())),
    "Delta_factorization": [list(row) for row in delta_factors],
    "finite_fibers": [list(row) for row in finite_signature],
    "infinity_fiber": infinity_signature,
    "root_data": list(root_data),
    "geometric_CM24_MW_rank": 2,
    "euler_number": euler_number,
    "generic_join": {
        "q4": {"q": 4, "a": 2, "b": 2, "v": q4_v, "MW_rank": 16},
        "q6": {"q": 6, "a": 2, "b": 3, "v": q6_v, "MW_rank": 17},
        "q6_rootless": True,
        "q6_full_nefness_verifier": (
            "elkies-k3/scripts/verify_q80_alternate_final_q6_nef.sage"
        ),
    },
    "rank_claim": None,
    "remaining_gate": "lift the generic pair23 gauge and q6 pencil to characteristic zero",
    "reproduce": (
        "sage elkies-k3/scripts/"
        "verify_q80_deforming_fifth_pair23_gf73.sage"
    ),
}
output = (
    ROOT / "artifacts/generated-results/"
    "q80-deforming-fifth-pair23-gf73.json"
)
encoded = json.dumps(artifact, indent=2, sort_keys=True, default=int)+"\n"
output.write_text(encoded)
digest = hashlib.sha256(encoded.encode()).hexdigest()
print(
    "Q80DEFORMINGPAIR23|"
    f"gauge=s=(U-mu_A)/(mu_B-mu_A)|factor_degrees={factor_degrees}|"
    f"fibers=I0*+3I4+3I2|root_data={root_data}|MW=2|euler={euler_number}|"
    f"q4_MW=16|q6_rootless_MW=17|status=PASS_EXACT_JOIN",
    flush=True,
)
print(
    "Q80DEFORMINGPAIR23|"
    f"artifact={output}|sha256={digest}|"
    "characteristic_zero_lift=open|rank_claim=none|status=PASS_ARTIFACT_WRITE",
    flush=True,
)
