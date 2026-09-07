#!/usr/bin/env sage
"""Recover a correctly twisted marked fifth-child Jacobian over GF(73).

The marked-projection artifact certifies a degree-four squarefree cover, but
its monic quartic alone drops the unit of the factorization over ``GF(73)(u)``.
That unit is the required quadratic twist.  Restoring it clears all apparent
denominators in the binary-quartic invariants and gives the degree-(8,12)
K3 Weierstrass model.  This checker also reads off the exact semistable fiber
signature from the polynomial discriminant, including infinity.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, FunctionField, PolynomialRing, gcd


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--mode", choices=("pair01", "pair14", "pair23-ref04"), default="pair14"
)
parser.add_argument("--write-artifact", action="store_true")
arguments = parser.parse_args()

SOURCE_NAMES = {
    "pair01": "q80-fifth-q4-marked-projection-pair01-gf73.json",
    "pair14": "q80-fifth-q4-marked-projection-pair14-gf73.json",
    "pair23-ref04": "q80-fifth-q4-marked-projection-pair23-ref04-gf73.json",
}
SOURCE = ROOT / "artifacts/generated-results" / SOURCE_NAMES[arguments.mode]
KNOWN_SOURCE_SHA256 = {
    "pair01": "f992b085aa5740f37649a69b83ff9ee6b2eafba620b60932a940b9a3bb3cb0a4",
    "pair14": "e46c9925c6870a6f9185f36994a5aef682382bba7a9bf8d2adc3d897420988fa",
    "pair23-ref04": "f131d957f2b2e876815b94fc60ad406f77f81152b95194fdca44cceda206d117",
}

source_bytes = SOURCE.read_bytes()
source_sha256 = hashlib.sha256(source_bytes).hexdigest()
assert source_sha256 == KNOWN_SOURCE_SHA256[arguments.mode]
source = json.loads(source_bytes)
assert source["mode"] == (
    "pair23" if arguments.mode == "pair23-ref04" else arguments.mode
)
if arguments.mode == "pair23-ref04":
    assert source["reference_pair_zero_one"] == [0, 4]
assert source["squarefree_cover_degree_T"] == 4

finite = GF(73, impl="modn")
new_base = FunctionField(finite, "u")
u = new_base.gen()
cover_ring = PolynomialRing(new_base, "tau")
tau = cover_ring.gen()
cover = cover_ring(sum(
    new_base(coefficient)*tau**tau_degree*u**u_degree
    for tau_degree, u_degree, coefficient
    in source["integral_double_cover_terms_T_U_coefficient"]
))
factorization = cover.factor()
factor_degrees_exponents = tuple(
    (int(factor.degree()), int(exponent))
    for factor, exponent in factorization
)
assert factor_degrees_exponents == tuple(
    tuple(map(int, row)) for row in source["cover_factor_degrees_exponents"]
)

twist = new_base(factorization.unit())
odd_part = cover_ring(1)
for factor, exponent in factorization:
    if int(exponent) % 2:
        odd_part *= factor
monic_quartic = odd_part.monic()
assert monic_quartic.degree() == 4

coefficients = list(monic_quartic.list())+[new_base(0)]*5
e, d, c, b, a = coefficients[:5]
invariant_I = 12*a*e-3*b*d+c**2
invariant_J = (
    72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
)

# Scaling the quartic by its omitted unit scales I,J by twist^2,twist^3.
jacobian_A = twist**2*(-27*invariant_I)
jacobian_B = twist**3*(-27*invariant_J)
delta_core = twist**6*(4*invariant_I**3-invariant_J**2)
assert jacobian_A.denominator() == 1
assert jacobian_B.denominator() == 1
assert delta_core.denominator() == 1
A = jacobian_A.numerator()
B = jacobian_B.numerator()
Delta = delta_core.numerator()
# A T-dependent normalization of the marked degree-two coordinate can leave
# nonminimal common powers in the integral binary-quartic model.  Strip the
# exact x=f^2*x', y=f^3*y' factors before reading the elliptic surface.
finite_scalings = []
for factor, _exponent in gcd(A, B).factor():
    scale_order = min(A.valuation(factor)//4, B.valuation(factor)//6)
    if scale_order <= 0:
        continue
    A //= factor**(4*scale_order)
    B //= factor**(6*scale_order)
    Delta //= factor**(12*scale_order)
    finite_scalings.append((str(factor), int(scale_order)))
model_degrees = (A.degree(), B.degree(), Delta.degree())
print(
    "Q80FIFTHMARKEDJAC|"
    f"mode={arguments.mode}|finite_scalings={tuple(finite_scalings)}|"
    f"minimalized_degrees={model_degrees}|stage=finite_minimalization",
    flush=True,
)
assert model_degrees[:2] == (8, 12) and model_degrees[2] <= 24

twist_factors = tuple(
    (str(factor), int(exponent))
    for factor, exponent in twist.numerator().factor()
)
delta_factorization = tuple(Delta.factor())
delta_factors = tuple(
    (str(factor), int(exponent)) for factor, exponent in delta_factorization
)

def kodaira_data(ord_a, ord_b, ord_delta):
    if ord_a == 0 or ord_b == 0:
        n = int(ord_delta)
        return (f"I{n}", n-1, n*(n-1), n, n)
    if ord_delta == 2:
        return ("II", 0, 0, 1, 2)
    if ord_delta == 3:
        return ("III", 1, 2, 2, 3)
    if ord_delta == 4:
        return ("IV", 2, 6, 3, 4)
    if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
        return ("I0*", 4, 24, 4, 6)
    if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
        n = int(ord_delta-6)
        rank = n+4
        return (f"I{n}*", rank, 2*rank*(rank-1), 4, n+6)
    if ord_delta == 8:
        return ("IV*", 6, 72, 3, 8)
    if ord_delta == 9:
        return ("III*", 7, 126, 2, 9)
    if ord_delta == 10:
        return ("II*", 8, 240, 1, 10)
    raise ArithmeticError((ord_a, ord_b, ord_delta))


finite_signature = []
root_rank = 0
root_count = 0
root_determinant = 1
euler_number = 0
for factor, exponent in delta_factorization:
    ord_a = int(A.valuation(factor))
    ord_b = int(B.valuation(factor))
    ord_delta = int(exponent)
    kind, rank, count, determinant, euler = kodaira_data(
        ord_a, ord_b, ord_delta
    )
    degree = int(factor.degree())
    finite_signature.append(
        (str(factor), degree, ord_a, ord_b, ord_delta, kind)
    )
    root_rank += degree*rank
    root_count += degree*count
    root_determinant *= determinant**degree
    euler_number += degree*euler
finite_signature = tuple(finite_signature)

infinity_orders = (8-A.degree(), 12-B.degree(), 24-Delta.degree())
infinity_data = kodaira_data(*infinity_orders)
infinity_kind, rank, count, determinant, euler = infinity_data
root_rank += rank
root_count += count
root_determinant *= determinant
euler_number += euler
assert euler_number == 24
if arguments.mode == "pair14":
    assert (root_rank, root_count, root_determinant) == (15, 82, 360)

print(
    "Q80FIFTHMARKEDJAC|"
    f"twist={twist}|twist_factors={twist_factors}|"
    f"finite_scalings={tuple(finite_scalings)}|"
    f"degrees_A_B_Delta={(A.degree(), B.degree(), Delta.degree())}|"
    f"delta_factors={delta_factors}",
    flush=True,
)
print(
    "Q80FIFTHMARKEDJAC|"
    f"finite_signature={finite_signature}|"
    f"infinity_orders={infinity_orders}|infinity={infinity_kind}|"
    f"root_data={(root_rank, root_count, root_determinant)}|"
    f"geometric_CM24_MW={18-root_rank}|"
    f"euler={euler_number}|status=PASS_TWISTED_K3_FIBER_SIGNATURE",
    flush=True,
)
if arguments.mode == "pair23-ref04":
    # This genus-one gauge is the old I6+I4+4I3+2I1 source fibration,
    # not the compensated deforming fifth target found by the secant gauge.
    assert (root_rank, root_count, root_determinant) == (16, 66, 1944)

if arguments.write_artifact:
    output = {
        "schema": f"q80-fifth-q4-marked-{arguments.mode}-jacobian-gf73-v1",
        "mode": arguments.mode,
        "prime": 73,
        "source_artifact": str(SOURCE.relative_to(ROOT)),
        "source_sha256": source_sha256,
        "factor_degrees_exponents": [list(row) for row in factor_degrees_exponents],
        "quadratic_twist": str(twist),
        "quadratic_twist_factorization": [list(row) for row in twist_factors],
        "A_coefficients_low_to_high": list(map(int, A.list())),
        "B_coefficients_low_to_high": list(map(int, B.list())),
        "Delta_coefficients_low_to_high": list(map(int, Delta.list())),
        "Delta_factorization": [list(row) for row in delta_factors],
        "finite_fiber_signature": [list(row) for row in finite_signature],
        "infinity_orders_A_B_Delta": list(map(int, infinity_orders)),
        "infinity_fiber": infinity_kind,
        "root_rank": root_rank,
        "root_count": root_count,
        "root_determinant": root_determinant,
        "geometric_CM24_MW_rank": 18-root_rank,
        "euler_number": euler_number,
        "reproduce": (
            "sage elkies-k3/scripts/"
            "analyze_q80_fifth_q4_marked_jacobian_gf73.sage "
            f"--mode {arguments.mode} --write-artifact"
        ),
    }
    output_path = (
        ROOT / "artifacts/generated-results/"
        f"q80-fifth-q4-marked-{arguments.mode}-jacobian-gf73.json"
    )
    encoded = json.dumps(output, indent=2, sort_keys=True, default=int)+"\n"
    output_path.write_text(encoded)
    digest = hashlib.sha256(encoded.encode()).hexdigest()
    print(
        "Q80FIFTHMARKEDJAC|"
        f"artifact={output_path}|sha256={digest}|status=PASS_ARTIFACT_WRITE",
        flush=True,
    )
