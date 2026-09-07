#!/usr/bin/env sage -python
"""Lift P12-O_pinned in the regular split-I12 resolved chart.

The ordinary pole-order-six section chart is singular (rank 43/48) at the
marked modular point.  Resolve it by imposing the central I12 component in a
local coordinate s=T-beta:

    Z monic of degree 6,
    X = (center*Z^2 mod s^6) + s^6 Q,  deg Q <= 10,
    Y = s^6 R,                         deg R <= 18.

Keeping the leading coefficients of Q and R free gives 36 variables and 37
equations.  The marked seed has Jacobian rank 36.  A square full-rank row
subsystem therefore supports ordinary multivariable Hensel/Newton lifting;
the unused equation is checked at every precision.  No Groebner basis or
multivariate characteristic-zero solve is used.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, QQ, ZZ, Zmod, lcm, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--max-exponent", type=int, default=131072)
parser.add_argument("--reconstruct-min-exponent", type=int, default=1024)
parser.add_argument(
    "--checkpoint",
    type=Path,
    default=LOCAL / "q24-a11-q8-difference-resolved-hensel-p100003.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-q8-difference-section-qq.json",
)
args = parser.parse_args()

MODEL = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
SEED = LOCAL / "q24-a11-q8-horizontal-points-mod100003.json"
INPUTS = (MODEL, SEED)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

model = json.loads(MODEL.read_text())
seed_artifact = json.loads(SEED.read_text())
assert model["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert seed_artifact["status"] == "PASS_Q24_A11_Q8_HORIZONTAL_POINTS_RECONSTRUCTION_MODP"
seed_section = seed_artifact["marked_difference_P12_minus_Opinned"]["section"]
p = ZZ(seed_artifact["prime"])
if p != 100003:
    raise ArithmeticError("resolved chart is pinned at p=100003")

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
if i12_factor.degree() != 1:
    raise ArithmeticError("A11 I12 factor is not linear")
beta = -QQ(i12_factor[0]) / QQ(i12_factor[1])

SQ = PolynomialRing(QQ, "s")
sQ = SQ.gen()
A_shift = SQ(A(sQ + beta))
B_shift = SQ(B(sQ + beta))
PSQ = PowerSeriesRing(QQ, "s", default_prec=14)
s_series = PSQ.gen()
A_series = PSQ(A_shift)
B_series = PSQ(B_shift)
center = PSQ(-3 * B_series[0] / (2 * A_series[0]))
for unused in range(7):
    center = (center + (-A_series / 3) / center) / 2
if (center**2 + A_series / 3).valuation() < 13:
    raise ArithmeticError("exact I12 center did not converge")
if (center**3 + A_series * center + B_series).valuation() != 12:
    raise ArithmeticError("exact centered cubic does not start in order 12")
C_shift = SQ([QQ(center[index]) for index in range(6)])

names = (
    [f"z{i}" for i in range(6)]
    + ["c_inf", "ell_inf"]
    + [f"q{i}" for i in range(10)]
    + [f"r{i}" for i in range(18)]
)
NZ = 6
NQ = 10
NR = 18
unknown_count = len(names)
equation_degrees = tuple(range(12, 49))
equation_count = len(equation_degrees)
assert unknown_count == 36 and equation_count == 37


def rational_reduction(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


# Convert the marked seed to the shifted resolved coordinates.
F = GF(p)
RTp = PolynomialRing(F, "T")
Tp = RTp.gen()
Rsp = PolynomialRing(F, "s")
sp = Rsp.gen()
beta_p = rational_reduction(beta, F)
X_seed_T = RTp(seed_section["X_coefficients_low_to_high"])
Y_seed_T = RTp(seed_section["Y_coefficients_low_to_high"])
Z_seed_T = RTp(seed_section["Z_coefficients_low_to_high"])
X_seed = Rsp(X_seed_T(sp + beta_p))
Y_seed = Rsp(Y_seed_T(sp + beta_p))
Z_seed = Rsp(Z_seed_T(sp + beta_p))
C_seed = Rsp([rational_reduction(value, F) for value in C_shift.list()])
local_seed = Rsp((C_seed * Z_seed**2) % sp**6)
Q_seed, remainder = (X_seed - local_seed).quo_rem(sp**6)
if remainder or Q_seed.degree() != 10:
    raise ArithmeticError("marked seed missed the resolved X chart")
R_seed, remainder = Y_seed.quo_rem(sp**6)
if remainder or R_seed.degree() != 18:
    raise ArithmeticError("marked seed missed the resolved Y chart")
if Z_seed.degree() != 6 or Z_seed[6] != 1:
    raise ArithmeticError("marked seed Z is not monic of degree six")

seed_values = (
    [ZZ(Z_seed[index]) for index in range(6)]
    + [ZZ(Q_seed[10]), ZZ(R_seed[18])]
    + [ZZ(Q_seed[index]) for index in range(10)]
    + [ZZ(R_seed[index]) for index in range(18)]
)
assert len(seed_values) == unknown_count


def common_denominator_data(poly):
    values = [QQ(value) for value in poly.list()]
    denominator = lcm([ZZ(value.denominator()) for value in values])
    return [ZZ(value * denominator) for value in values], denominator


A_numerators, A_denominator = common_denominator_data(A_shift)
B_numerators, B_denominator = common_denominator_data(B_shift)
C_numerators, C_denominator = common_denominator_data(C_shift)
reduction_cache = {}


def mod_ring(modulus):
    modulus = ZZ(modulus)
    S = Zmod(modulus, is_field=(modulus == p))
    if modulus != p:
        S.is_field.set_cache(False)
    return S


def reduce_stored(numerators, denominator, modulus, S):
    if denominator.gcd(modulus) != 1:
        raise ZeroDivisionError("exact chart denominator is not a modular unit")
    inverse = (denominator % modulus).inverse_mod(modulus)
    return [S((value % modulus) * inverse % modulus) for value in numerators]


def reduced_chart_data(modulus):
    modulus = ZZ(modulus)
    if modulus in reduction_cache:
        return reduction_cache[modulus]
    S = mod_ring(modulus)
    answer = (
        reduce_stored(A_numerators, A_denominator, modulus, S),
        reduce_stored(B_numerators, B_denominator, modulus, S),
        reduce_stored(C_numerators, C_denominator, modulus, S),
    )
    reduction_cache[modulus] = answer
    if len(reduction_cache) > 2:
        del reduction_cache[next(iter(reduction_cache))]
    return answer


def chart_polynomials(values, modulus):
    S = mod_ring(modulus)
    RS = PolynomialRing(S, "s")
    s = RS.gen()
    z_values = values[:6]
    c_inf, ell_inf = values[6:8]
    q_values = values[8:18]
    r_values = values[18:36]
    Z = s**6 + sum(S(z_values[index]) * s**index for index in range(6))
    Q = S(c_inf) * s**10 + sum(S(q_values[index]) * s**index for index in range(10))
    RR = S(ell_inf) * s**18 + sum(S(r_values[index]) * s**index for index in range(18))
    A_values, B_values, C_values = reduced_chart_data(modulus)
    Ap = RS(A_values)
    Bp = RS(B_values)
    Cp = RS(C_values)
    local = RS((Cp * Z**2) % s**6)
    X = local + s**6 * Q
    Y = s**6 * RR
    return RS, s, X, Y, Z, Ap, Bp, Cp


def residual(values, modulus):
    unused_RS, unused_s, X, Y, Z, Ap, Bp, unused_Cp = chart_polynomials(
        values, modulus
    )
    polynomial = Y**2 - X**3 - Ap * X * Z**4 - Bp * Z**6
    if any(polynomial[index] for index in range(12)):
        raise ArithmeticError("resolved chart lost its automatic low-order vanishing")
    return [polynomial[index] for index in equation_degrees]


def jacobian(values, modulus):
    RS, s, X, Y, Z, Ap, Bp, Cp = chart_polynomials(values, modulus)
    x_factor = -3 * X**2 - Ap * Z**4
    z_factor = -4 * Ap * X * Z**3 - 6 * Bp * Z**5
    columns = []
    for index in range(6):
        dZ = s**index
        dX = RS((2 * Cp * Z * dZ) % s**6)
        derivative = x_factor * dX + z_factor * dZ
        columns.append([derivative[degree] for degree in equation_degrees])
    columns.append(
        [(x_factor * s**16)[degree] for degree in equation_degrees]
    )
    columns.append(
        [(2 * Y * s**24)[degree] for degree in equation_degrees]
    )
    for index in range(10):
        derivative = x_factor * s ** (6 + index)
        columns.append([derivative[degree] for degree in equation_degrees])
    for index in range(18):
        derivative = 2 * Y * s ** (6 + index)
        columns.append([derivative[degree] for degree in equation_degrees])
    if len(columns) != unknown_count:
        raise AssertionError("resolved Jacobian column count mismatch")
    return matrix(
        mod_ring(modulus),
        equation_count,
        unknown_count,
        lambda row, column: columns[column][row],
    )


if any(int(value) for value in residual(seed_values, p)):
    raise ArithmeticError("marked seed misses the resolved equations")
Jp = jacobian(seed_values, p)
if Jp.rank() != unknown_count:
    raise ArithmeticError(f"resolved seed rank {Jp.rank()}, expected {unknown_count}")
pivot_rows = list(Jp.transpose().pivots())
if len(pivot_rows) != unknown_count:
    raise ArithmeticError("resolved seed did not supply a square pivot block")
print(
    f"A11Q8DIFFRESOLVED|prime={p}|variables={unknown_count}|equations={equation_count}|"
    f"rank={Jp.rank()}|status=PASS_REGULAR_SEED",
    flush=True,
)

input_hashes = {
    str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
    for path in INPUTS
}
checkpoint = args.checkpoint.resolve()
output = args.output.resolve()

if checkpoint.exists():
    saved = json.loads(checkpoint.read_text())
    if saved.get("schema") != "elkies-k3.h3-q24-a11-q8-difference-resolved-hensel.v1":
        raise ArithmeticError("resolved checkpoint schema mismatch")
    if saved.get("input_sha256") != input_hashes or ZZ(saved.get("prime")) != p:
        raise ArithmeticError("resolved checkpoint prerequisites changed")
    exponent = int(saved["exponent"])
    values = [ZZ(value) for value in saved["coefficient_residues"]]
    if saved["pivot_equation_rows"] != pivot_rows:
        raise ArithmeticError("resolved checkpoint pivots changed")
else:
    exponent = 1
    values = [value % p for value in seed_values]


def write_checkpoint():
    modulus = p**exponent
    payload = {
        "schema": "elkies-k3.h3-q24-a11-q8-difference-resolved-hensel.v1",
        "status": "PARTIAL_Q24_A11_Q8_DIFFERENCE_RESOLVED_HENSEL",
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


def exact_candidate():
    modulus = p**exponent
    reconstructed = [reconstruct_scalar(value, modulus) for value in values]
    print(
        f"A11Q8DIFFRESOLVEDRECONSTRUCT|exponent={exponent}|"
        f"reconstructed={sum(value is not None for value in reconstructed)}/{unknown_count}",
        flush=True,
    )
    if any(value is None for value in reconstructed):
        return None
    z_values = reconstructed[:6]
    c_inf, ell_inf = reconstructed[6:8]
    q_values = reconstructed[8:18]
    r_values = reconstructed[18:36]
    Zs = sQ**6 + sum(z_values[index] * sQ**index for index in range(6))
    Qs = c_inf * sQ**10 + sum(q_values[index] * sQ**index for index in range(10))
    Rs = ell_inf * sQ**18 + sum(r_values[index] * sQ**index for index in range(18))
    Xs = SQ((C_shift * Zs**2) % sQ**6) + sQ**6 * Qs
    Ys = sQ**6 * Rs
    if Ys**2 != Xs**3 + A_shift * Xs * Zs**4 + B_shift * Zs**6:
        return None
    XT = RQ(Xs(T - beta))
    YT = RQ(Ys(T - beta))
    ZT = RQ(Zs(T - beta))
    if YT**2 != XT**3 + A * XT * ZT**4 + B * ZT**6:
        raise ArithmeticError("unshifted resolved section missed A11")
    return XT, YT, ZT, reconstructed


answer = None
while exponent < args.max_exponent:
    next_exponent = min(2 * exponent, args.max_exponent)
    current_modulus = p**exponent
    quotient_modulus = p ** (next_exponent - exponent)
    target_modulus = current_modulus * quotient_modulus
    target_residual = [ZZ(value) for value in residual(values, target_modulus)]
    if any(value % current_modulus for value in target_residual):
        raise ArithmeticError("resolved point lost its current Hensel precision")
    rhs = [(-value // current_modulus) % quotient_modulus for value in target_residual]
    J = jacobian(values, quotient_modulus)
    square = J[pivot_rows, :]
    right = matrix(
        mod_ring(quotient_modulus),
        unknown_count,
        1,
        [rhs[index] for index in pivot_rows],
    )
    correction = square.solve_right(right).column(0)
    values = [
        (values[index] + current_modulus * ZZ(correction[index])) % target_modulus
        for index in range(unknown_count)
    ]
    exponent = next_exponent
    if any(int(value) for value in residual(values, target_modulus)):
        raise ArithmeticError("resolved unused equation or Newton update failed")
    write_checkpoint()
    print(
        f"A11Q8DIFFRESOLVEDHENSEL|exponent={exponent}|bits={target_modulus.nbits()}|"
        f"rank={unknown_count}|all_equations={equation_count}|status=PASS_STEP",
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
        f"A11Q8DIFFRESOLVED_RESULT|exponent={exponent}|bits={(p**exponent).nbits()}|"
        "status=NEED_MORE_RESOLVED_HENSEL_PRECISION",
        flush=True,
    )
    print(f"CHECKPOINT|{checkpoint}", flush=True)
    raise SystemExit(0)

X, Y, Z, reconstructed = answer
max_bits = max(
    max(abs(ZZ(value.numerator())).nbits(), abs(ZZ(value.denominator())).nbits())
    for poly in (X, Y, Z)
    for value in poly
)
payload = {
    "schema": "elkies-k3.h3-q24-a11-q8-difference-section-qq.v1",
    "status": "PASS_EXACT_Q24_A11_Q8_DIFFERENCE_SECTION_QQ",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": input_hashes,
    },
    "resolved_chart": {
        "shift": f"s=T-({beta})",
        "variables": names,
        "variable_count": unknown_count,
        "equation_count": equation_count,
        "seed_rank": unknown_count,
        "I12_component_depth": 6,
        "pole_order": 6,
        "hensel_prime": int(p),
        "hensel_exponent": exponent,
    },
    "section": {
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [str(value) for value in Z.list()],
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "exact_weierstrass_identity": True,
    },
    "marking": {
        "relation": "P12-O_pinned",
        "equation_P_dot_O": 6,
        "equation_I12_component_depth_up_to_negation": 6,
        "exact_reduction_matches_marked_seed": True,
    },
    "max_rational_coefficient_bits": int(max_bits),
    "method": "regular split-I12 resolved chart plus checkpointed multivariable Hensel",
    "large_Groebner_required": False,
    "proof_boundary": (
        "This exactly constructs the marked A11 difference P12-O_pinned. "
        "Combining it with the exact pinned zero gives the second q8 horizontal section."
    ),
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11Q8DIFFRESOLVED_RESULT|degrees={}|exponent={}|max_bits={}|status={}".format(
        tuple(payload["section"]["degrees_X_Y_Z"]),
        exponent,
        max_bits,
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{output}", flush=True)
