#!/usr/bin/env sage -python
"""Lift R=P12-M in its regular split-I12 component-3 chart.

The ordinary degree-(4,6,0) polynomial-section chart is singular of rank
10/12 at the marked point.  At the split I12 fibre, however, R meets component
3.  In the shifted coordinate s this gives

    Z = 1,
    X = center mod s^3 + s^3*(q0+q1*s),
    Y = s^3*(r0+r1*s+r2*s^2+r3*s^3).

There are six variables and seven coefficient equations; the marked seed has
rank six.  Only q0,q1 are rationally reconstructed, after which Y is the
unique seed-compatible polynomial square root.  No Groebner basis is used.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, QQ, ZZ, Zmod, lcm, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
SEED = LOCAL / "q24-a11-q8-horizontal-points-mod100003.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--max-exponent", type=int, default=65536)
parser.add_argument("--reconstruct-min-exponent", type=int, default=1024)
parser.add_argument(
    "--checkpoint",
    type=Path,
    default=LOCAL / "q24-a11-q8-residual-resolved-hensel-p100003.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-q8-residual-section-qq.json",
)
args = parser.parse_args()

for path in (MODEL, SEED):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

model = json.loads(MODEL.read_text())
seed_artifact = json.loads(SEED.read_text())
assert model["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert seed_artifact["status"] == "PASS_Q24_A11_Q8_HORIZONTAL_POINTS_RECONSTRUCTION_MODP"
seed_section = seed_artifact["residual_P12_minus_M"]["section"]
p = ZZ(seed_artifact["prime"])
assert p == 100003

RQ = PolynomialRing(QQ, "T")
T = RQ.gen()
A = RQ([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B = RQ([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])
i12_factor = RQ(
    next(
        row["factor"]
        for row in model["child"]["discriminant_factorization"]
        if int(row["multiplicity"]) == 12
    )
)
beta = -QQ(i12_factor[0]) / QQ(i12_factor[1])

SQ = PolynomialRing(QQ, "s")
sQ = SQ.gen()
A_shift = SQ(A(sQ + beta))
B_shift = SQ(B(sQ + beta))
PSQ = PowerSeriesRing(QQ, "s", default_prec=14)
A_series = PSQ(A_shift)
B_series = PSQ(B_shift)
center = PSQ(-3 * B_series[0] / (2 * A_series[0]))
for unused in range(7):
    center = (center + (-A_series / 3) / center) / 2
if (center**2 + A_series / 3).valuation() < 13:
    raise ArithmeticError("exact I12 centre did not converge")
if (center**3 + A_series * center + B_series).valuation() != 12:
    raise ArithmeticError("centred cubic does not start in order 12")
C_shift = SQ([QQ(center[index]) for index in range(3)])


def rational_reduction(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


F = GF(p)
RTp = PolynomialRing(F, "T")
Tp = RTp.gen()
Rsp = PolynomialRing(F, "s")
sp = Rsp.gen()
beta_p = rational_reduction(beta, F)
X_seed = Rsp(RTp(seed_section["X_coefficients_low_to_high"])(sp + beta_p))
Y_seed = Rsp(RTp(seed_section["Y_coefficients_low_to_high"])(sp + beta_p))
Z_seed = Rsp(RTp(seed_section["Z_coefficients_low_to_high"])(sp + beta_p))
if Z_seed != 1:
    raise ArithmeticError("residual seed is not polynomial in the resolved chart")
C_seed = Rsp([rational_reduction(value, F) for value in C_shift.list()])
Q_seed, remainder = (X_seed - C_seed).quo_rem(sp**3)
R_seed, remainder_y = Y_seed.quo_rem(sp**3)
if remainder or remainder_y or Q_seed.degree() != 1 or R_seed.degree() != 3:
    raise ArithmeticError("residual seed misses component-3 chart")

names = ["q0", "q1", "r0", "r1", "r2", "r3"]
equation_degrees = tuple(range(6, 13))
unknown_count = 6
equation_count = 7
seed_values = [ZZ(Q_seed[index]) for index in range(2)] + [ZZ(R_seed[index]) for index in range(4)]


def denominator_data(poly):
    values = [QQ(value) for value in poly.list()]
    denominator = lcm([ZZ(value.denominator()) for value in values])
    return [ZZ(value * denominator) for value in values], denominator


A_num, A_den = denominator_data(A_shift)
B_num, B_den = denominator_data(B_shift)
C_num, C_den = denominator_data(C_shift)
reduction_cache = {}


def mod_ring(modulus):
    modulus = ZZ(modulus)
    ring = Zmod(modulus, is_field=(modulus == p))
    if modulus != p:
        ring.is_field.set_cache(False)
    return ring


def reduce_stored(numerators, denominator, modulus, ring):
    inverse = (denominator % modulus).inverse_mod(modulus)
    return [ring((value % modulus) * inverse % modulus) for value in numerators]


def chart_data(modulus):
    modulus = ZZ(modulus)
    if modulus in reduction_cache:
        return reduction_cache[modulus]
    ring = mod_ring(modulus)
    answer = (
        reduce_stored(A_num, A_den, modulus, ring),
        reduce_stored(B_num, B_den, modulus, ring),
        reduce_stored(C_num, C_den, modulus, ring),
    )
    reduction_cache[modulus] = answer
    if len(reduction_cache) > 2:
        del reduction_cache[next(iter(reduction_cache))]
    return answer


def polynomials(values, modulus):
    ring = mod_ring(modulus)
    RS = PolynomialRing(ring, "s")
    s = RS.gen()
    q0, q1, r0, r1, r2, r3 = map(ring, values)
    A_values, B_values, C_values = chart_data(modulus)
    Ap, Bp, Cp = RS(A_values), RS(B_values), RS(C_values)
    X = Cp + s**3 * (q0 + q1 * s)
    Y = s**3 * (r0 + r1 * s + r2 * s**2 + r3 * s**3)
    return RS, s, X, Y, Ap, Bp


def residual(values, modulus):
    unused_RS, unused_s, X, Y, Ap, Bp = polynomials(values, modulus)
    identity = Y**2 - X**3 - Ap * X - Bp
    if any(identity[index] for index in range(6)):
        raise ArithmeticError("component-3 chart lost automatic vanishing")
    return [identity[index] for index in equation_degrees]


def jacobian(values, modulus):
    RS, s, X, Y, Ap, unused_Bp = polynomials(values, modulus)
    x_factor = -3 * X**2 - Ap
    columns = []
    for index in range(2):
        derivative = x_factor * s ** (3 + index)
        columns.append([derivative[degree] for degree in equation_degrees])
    for index in range(4):
        derivative = 2 * Y * s ** (3 + index)
        columns.append([derivative[degree] for degree in equation_degrees])
    return matrix(
        mod_ring(modulus), equation_count, unknown_count,
        lambda row, column: columns[column][row],
    )


if any(int(value) for value in residual(seed_values, p)):
    raise ArithmeticError("marked residual misses resolved equations")
Jp = jacobian(seed_values, p)
if Jp.rank() != unknown_count:
    raise ArithmeticError(f"resolved residual rank {Jp.rank()}, expected {unknown_count}")
pivot_rows = list(Jp.transpose().pivots())
print(
    f"A11Q8RESIDUALRESOLVED|prime={p}|variables=6|equations=7|rank=6|"
    "status=PASS_REGULAR_SEED",
    flush=True,
)

input_hashes = {
    str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
    for path in (MODEL, SEED)
}
checkpoint = args.checkpoint.resolve()
output = args.output.resolve()
if checkpoint.exists():
    saved = json.loads(checkpoint.read_text())
    if saved.get("schema") != "elkies-k3.h3-q24-a11-q8-residual-resolved-hensel.v1":
        raise ArithmeticError("residual checkpoint schema mismatch")
    if saved.get("input_sha256") != input_hashes or ZZ(saved.get("prime")) != p:
        raise ArithmeticError("residual checkpoint prerequisites changed")
    exponent = int(saved["exponent"])
    values = [ZZ(value) for value in saved["coefficient_residues"]]
    if saved["pivot_equation_rows"] != pivot_rows:
        raise ArithmeticError("residual checkpoint pivots changed")
else:
    exponent = 1
    values = [value % p for value in seed_values]


def write_checkpoint():
    modulus = p**exponent
    payload = {
        "schema": "elkies-k3.h3-q24-a11-q8-residual-resolved-hensel.v1",
        "status": "PARTIAL_Q24_A11_Q8_RESIDUAL_RESOLVED_HENSEL",
        "prime": int(p),
        "exponent": exponent,
        "modulus_bits": int(modulus.nbits()),
        "variables": names,
        "pivot_equation_rows": pivot_rows,
        "coefficient_residues": [str(value % modulus) for value in values],
        "input_sha256": input_hashes,
        "large_Groebner_required": False,
    }
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    checkpoint.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def reconstruct_scalar(value, modulus):
    try:
        return QQ((ZZ(value) % modulus).rational_reconstruction(modulus))
    except (ArithmeticError, ValueError, ZeroDivisionError):
        return None


def reduce_exact(value, modulus):
    value = QQ(value)
    return ZZ(value.numerator()) * ZZ(value.denominator()).inverse_mod(modulus) % modulus


def exact_candidate():
    modulus = p**exponent
    q_values = [reconstruct_scalar(values[index], modulus) for index in range(2)]
    print(
        f"A11Q8RESIDUALRESOLVEDRECONSTRUCT|exponent={exponent}|"
        f"Q={sum(value is not None for value in q_values)}/2",
        flush=True,
    )
    if any(value is None for value in q_values):
        return None
    Xs = C_shift + sQ**3 * (q_values[0] + q_values[1] * sQ)
    right = Xs**3 + A_shift * Xs + B_shift
    if right.degree() != 12 or not QQ(right[12]).is_square():
        return None
    leading = QQ(right[12]).sqrt()
    if reduce_exact(leading, p) != seed_values[-1] % p:
        leading = -leading
    coefficients = [QQ.zero()] * 7
    coefficients[6] = leading
    for target_degree in range(11, 5, -1):
        unknown_index = target_degree - 6
        known = sum(
            coefficients[left] * coefficients[target_degree - left]
            for left in range(7)
            if 0 <= target_degree - left < 7
            and left != unknown_index
            and target_degree - left != unknown_index
        )
        coefficients[unknown_index] = (QQ(right[target_degree]) - known) / (2 * leading)
    Ys = SQ(coefficients)
    if Ys**2 != right or Ys.valuation() != 3:
        return None
    XT = RQ(Xs(T - beta))
    YT = RQ(Ys(T - beta))
    if YT**2 != XT**3 + A * XT + B:
        raise ArithmeticError("unshifted exact residual misses A11")
    # Full marked-prime reduction, including the selected y sign.
    XT_p = RTp([rational_reduction(value, F) for value in XT.list()])
    YT_p = RTp([rational_reduction(value, F) for value in YT.list()])
    if XT_p != RTp(seed_section["X_coefficients_low_to_high"]):
        raise ArithmeticError("exact residual x does not reduce to marked seed")
    if YT_p != RTp(seed_section["Y_coefficients_low_to_high"]):
        raise ArithmeticError("exact residual y does not reduce to marked seed")
    return XT, YT


answer = None
while exponent < args.max_exponent:
    next_exponent = min(2 * exponent, args.max_exponent)
    current_modulus = p**exponent
    quotient_modulus = p ** (next_exponent - exponent)
    target_modulus = current_modulus * quotient_modulus
    target_residual = [ZZ(value) for value in residual(values, target_modulus)]
    if any(value % current_modulus for value in target_residual):
        raise ArithmeticError("residual lift lost current precision")
    rhs = [(-value // current_modulus) % quotient_modulus for value in target_residual]
    J = jacobian(values, quotient_modulus)
    square = J[pivot_rows, :]
    right = matrix(
        mod_ring(quotient_modulus), unknown_count, 1,
        [rhs[index] for index in pivot_rows],
    )
    correction = square.solve_right(right).column(0)
    values = [
        (values[index] + current_modulus * ZZ(correction[index])) % target_modulus
        for index in range(unknown_count)
    ]
    exponent = next_exponent
    if any(int(value) for value in residual(values, target_modulus)):
        raise ArithmeticError("unused residual equation or Newton update failed")
    write_checkpoint()
    print(
        f"A11Q8RESIDUALRESOLVEDHENSEL|exponent={exponent}|bits={target_modulus.nbits()}|"
        "all_equations=7|status=PASS_STEP",
        flush=True,
    )
    if exponent >= args.reconstruct_min_exponent:
        answer = exact_candidate()
        if answer is not None:
            break

if answer is None and exponent >= args.reconstruct_min_exponent:
    answer = exact_candidate()
if answer is None:
    print(
        f"A11Q8RESIDUALRESOLVED_RESULT|exponent={exponent}|bits={(p**exponent).nbits()}|"
        "status=NEED_MORE_PRECISION",
        flush=True,
    )
    print(f"CHECKPOINT|{checkpoint}", flush=True)
    raise SystemExit(0)

X, Y = answer
max_bits = max(
    max(abs(ZZ(value.numerator())).nbits(), abs(ZZ(value.denominator())).nbits())
    for poly in (X, Y)
    for value in poly
)
payload = {
    "schema": "elkies-k3.h3-q24-a11-q8-residual-section-qq.v1",
    "status": "PASS_EXACT_Q24_A11_Q8_RESIDUAL_SECTION_QQ",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MODEL, SEED)],
        "sha256": input_hashes,
    },
    "resolved_chart": {
        "shift": f"s=T-({beta})",
        "variables": names,
        "variable_count": 6,
        "equation_count": 7,
        "seed_rank": 6,
        "I12_component_depth": 3,
        "hensel_prime": int(p),
        "hensel_exponent": exponent,
        "independently_reconstructed_coefficients": 2,
    },
    "section": {
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "Z_coefficients_low_to_high": ["1"],
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), 0],
        "exact_weierstrass_identity": True,
    },
    "marking": {
        "relation": "P12-M",
        "pinned_A11_MW": [-1, 0, -1, 0, 0, 0],
        "exact_reduction_matches_marked_seed": True,
    },
    "max_rational_coefficient_bits": int(max_bits),
    "method": "regular split-I12 component-3 chart; reconstruct two x scalars and derive y by exact square root",
    "large_Groebner_required": False,
    "proof_boundary": (
        "This exactly constructs the marked degree-one residual P12-M on A11. "
        "The exact translated q8 horizontal, H0 plane, quartic, and 2A5 Jacobian remain separate gates."
    ),
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"A11Q8RESIDUALRESOLVED_RESULT|degrees={(X.degree(),Y.degree(),0)}|"
    f"exponent={exponent}|max_bits={max_bits}|status={payload['status']}",
    flush=True,
)
print(f"OUTPUT|{output}", flush=True)
